// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date February 2019
// 
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <minus/minus.h>
#include <minus/chicago14a-io.h>
#include <minus/chicago-default.h>

using namespace MiNuS;
#define Float double
typedef minus_core<chicago> M;
static constexpr Float tol = 1e-3;
typedef std::complex<Float> complex;
using namespace std::chrono;


#define  M_VERBOSE 1     // display verbose messages
// comment this out plus the printouts if you prefer no debug and no dependency
// on this header
#include <minus/debug-util.h>


//  exit code. Conventions:
//  0 (Zero)	Success
// Non-zero	Failure
// 2	Incorrect usage
// 127	Command Not found
// 126	Not an executable 
// 
// The return value of a command is its exit status, or 128
// + N if the command is terminated by signal N. 
//
// Exit status is used to check the
// result (success/failure) of the execution of the command. If the exit status
// is zero, then the command is success. If the command is failed the exit status
// will be non-zero.
//
// The codes below are between 3 and 126 inclusive
//
#define SOLVER_FAILURE 3

void
print_usage()
{
  std::cerr << "Usage: minus [input solutions]\n\n";
  std::cerr << "If no argument is given, 'input' is assumed stdin,\n\
  'solutions' will be output to stdout\n";
  std::cerr << "Example: \n"
               "  minus input_file solutions_file\n"
               "  minus <input_file >solutions_file\n"
               "  minus -g       # (or --profile) : performs a default solve for profiling\n"
               "  minus -i       # (or --image_data) : reads point-tangents from stdin\n"
               "  minus -h       # (or --help) : print this help message\n"
               // "  minus -r       # (or --real)  :  outputs only real solutions\n"
               "  minus -AB      # (or --two_problems) : continue between 2 given problems\n"
            <<
  R"(-i | --image_data usage:
 
  Input format (notation _view_points_coords. any number of spaces and newlines optional. can be in
  one row or one column as well). This input format assumes tangent data for
  all points, but you specify which one to use in id0 and id1 below. When
  --use_all_tangents is passed (TODO), will try to select the better conditioned / least degenerate tangents 
 
  p000 p001        # If continuing from a standard internal problem to a new problem A, this is problem A
  p010 p011        # If continuing from two problems from A to B (flag -AB), this is also problem A
  p020 p021
  
  p100 p101
  p110 p111
  p120 p121
  
  p100 p101
  p110 p111
  p120 p121
 
  t000 t001
  t010 t011
  t020 t021
  
  t100 t101
  t110 t111
  t120 t121
  
  t100 t101
  t110 t111
  t120 t121
  
  id0 id1           # id \in {0,1,2} of the point to consider the tangent
  
  K00 K01 K02       # intrinsic parameters: only these elements
   0  K11 K22

  r000 r001 r002    # GROUND TRUTH (optional) if -gt flag provided, pass the ground truth here:
  r010 r011 r012    # default camera format if synthcurves flag passed: 
  r020 r021 r022    # just like a 3x4 [R|T] but transposed to better fit row-major:
   c00  c01  c02    #         | R |
                    # P_4x3 = | - |
  r100 r101 r102    #         | C'|
  r110 r111 r112    # 
  r120 r121 r122    #  
   c10  c11  c12    #                                                                                                                   # If two problems A->B are provided (flag -AB), this is only for problem B below
                    #
  r200 r201 r202    # 
  r210 r211 r212    # 
  r220 r221 r222    #
   c20  c21  c22    # 

  p000 p001         # If two problems A->B are provided (flag -AB), this is problem B
  p010 p011
  p020 p021
  
  p100 p101
  p110 p111
  p120 p121
  
  p100 p101
  p110 p111
  p120 p121
 
  t000 t001
  t010 t011
  t020 t021
  
  t100 t101
  t110 t111
  t120 t121
  
  t100 t101
  t110 t111
  t120 t121
  
  id0 id1           # id \in {0,1,2} of the point to consider the tangent

  # One way to use this is 
  #     synthdata | minus-chicago -i
  # where synthdata is provided in minus/scripts)";
             
  exit(1);
}

bool stdio_ = true;  // by default read/write from stdio
bool ground_truth_ = false;
bool two_problems_given_ = false;
bool reading_first_point_ = true;
std::ifstream infp_;
bool image_data_ = false;
bool profile_ = false;   // run some default solves for profiling
const char *input_ = "stdin";
const char *output_ = "stdout";
M::track_settings settings_;

void
print_num_steps(M::solution solutions[M::nsols])
{
  LOG("solution id x num steps:");
  unsigned sum=0;
  for (unsigned s=0; s < M::nsols; ++s) {
    LOG(s << " " << solutions[s].num_steps);
    sum += solutions[s].num_steps;
  }
  LOG("total number of steps: " << sum);
}

// Output solutions in ASCII matlab format
//
// ---------------------------------------------------------
// If in the future our solver is really fast, we may need Binary IO:
// complex solutions[NSOLS*NVE];
// 
// To read this output file in matlab, do:
// fid = fopen(fname,'r');
// a_raw = fread(fid,'double');
// fclose(fid);
//
// Reshape a to have proper real and imaginary parts
// a = a_raw(1:2:end) + i*a_raw(2:2:end);
// 
template <typename F=double>
static bool
mwrite(const M::solution s[M::nsols], const char *fname)
{
  bool scilab=false;
  std::string imag("+i*");
  if (scilab) imag = std::string("+%i*");
    
  std::ofstream fsols;
  std::streambuf *buf;
  
  if (stdio_) {
    buf = std::cout.rdbuf();
    std::cout << std::setprecision(20);
  } else {
    fsols.open(fname,std::ios::out);
    if (!fsols) {
      std::cerr << "minus: error, unable to open file name" << std::endl;
      return false;
    }
    buf = fsols.rdbuf();
    fsols << std::setprecision(20);
  }
  
  std::ostream out(buf);
  out << std::setprecision(20);
  out << "[";
  for (unsigned i=0; i <M::nsols; ++i) {
    for (unsigned var=0; var < M::nve; ++var) {
      out << s[i].x[var].real() << imag << s[i].x[var].imag();
      if (i*var +1 < M::nve * M::nsols) 
        out << std::endl;
      // BINARY fsols.write((char *)(s[i].x[var]),2*sizeof(double));
    }
  }
  out << "]\n";
  
  if (!stdio_) fsols.close();
  return true;
}
// Try to read n elements, filling in p in row-major order.
template <typename F=double>
static bool
read_block(std::istream &in, F *p, unsigned n)
{
  LOG("reading");
  const F *end = p + n;
  while (!in.eof() && p != end) {
      try {
        in >> *p++;
//        std::cerr << *(p-1) << std::endl;
        if (in.eof()) {
          std::cerr << "I/O Error: Premature input termination\n";
          return false;
        }
      } catch (std::istream::failure &E) {
        std::cerr << "I/O Error: Invalid input conversion or other error\n";
        return false;
      }
  }
  if (p != end) {
    std::cerr << "I/O Premature input termination\n";
    return false;
  }
  return true;
}

static bool
init_input(const char *fname, std::istream *inp)
{
  if (!stdio_) {
    infp_.open(fname, std::ios::in);
    if (!infp_) {
      std::cerr << "I/O Error opening input " << fname << std::endl;
      return false;
    }
    inp = &infp_;
  }
  inp->exceptions(std::istream::failbit | std::istream::badbit);
  return true;
}
  
//
// Reads the format specified in the print_usage() for the -i flag
// 
// This is processed into the global params_start_target_
// 
template <typename F=double>
static bool
iread(std::istream &in)
{
  LOG("reading p_");
  if (!read_block(in, (F *)data::p_, io::pp::nviews*io::pp::npoints*io::ncoords2d))
    return false;
  LOG("reading tgt_");
  if (!read_block(in, (F *)data::tgt_, io::pp::nviews*io::pp::npoints*io::ncoords2d))
    return false;
  unsigned tgt_ids[2];
  LOG("reading tgt_ids");
  if (!read_block<unsigned>(in, tgt_ids, 2))
    return false;
  if (reading_first_point_) {
    LOG("reading K_");
    if (!read_block(in, (F *) data::K_, io::ncoords2d*io::ncoords2d_h))
      return false;
    LOG("reading ground truth cams");
    if (ground_truth_ && !read_block(in, (F *) data::cameras_gt_, io::pp::nviews*4*3))
      return false;
    io::point_tangents2params_img(data::p_, data::tgt_, tgt_ids[0], tgt_ids[1],
        data::K_, data::params_start_target_);
    reading_first_point_ = false;
  } else { // when reading second point B, do not gammify A again
    static constexpr bool gammify_target_problem = false;
    io::point_tangents2params_img(data::p_, data::tgt_, tgt_ids[0], tgt_ids[1],
        data::K_, data::params_start_target_, gammify_target_problem);
  }
  return true;
}

// reads into the global variable params_
// Format is just like P01 variable in solveChicago in chicago.m2
// and contains the concatenated parameters of the start system
// and of the target system, with some randomization to improve conditioning.
// But here there is no imaginary 'i' string:
//
// P01(0).real()  P01(0).imag()      // I mean P01(0) or P01#0
// P01(1).real()  P01(1).imag()
// P01(2).real()  P01(2).imag()
// P01(3).real()  P01(3).imag()
// ...
//
// The file can also be one line, listing the above in row-major order like so:
// 
// P01(0).real()  
// P01(0).imag() 
// P01(1).real()
// P01(1).imag()
// ...
// 
// It is up to the user to build this from an actual input for a target system,
// be it point-tangents as in Ric's format, be it a linecomplex as in Hongy's format
//
// This format is generic enough to be adapted to M2 or matlab
template <typename F=double>
static bool
mread(std::istream &in)
{
  F *dparams = (F *)data::params_;
  while (!in.eof() && dparams != (F *)data::params_+2*2*M::f::nparams) {
      try {
      in >> *dparams++;
      // std::cerr << "reading " <<  *(dparams-1) << std::endl;;
      if (in.eof()) {
        std::cerr << "I/O Error: Premature input termination\n";
        return false;
      }
      in >> *dparams++;
      } catch (std::istream::failure &E) {
        std::cerr << "I/O Error: Invalid input conversion or other error\n";
        return false;
      }
  }
  if (dparams != (F *)data::params_+2*2*M::f::nparams)
    std::cerr << "I/O Premature input termination\n";
//  for (unsigned i=0; i < 2*NPARAMS; ++i)
//    std::cerr << "D " << params_[i] << std::endl;
  return true;
}

void
print_settings(const M::track_settings &settings)
{
  #ifdef M_VERBOSE
  std::cerr << "Track settings ------------------------------------------------\n";
  const char *names[11] = {
    "init_dt_",
    "min_dt_",
    "end_zone_factor_",
    "epsilon_",
    "epsilon2_",
    "dt_increase_factor_",
    "dt_decrease_factor_",
    "infinity_threshold_",
    "infinity_threshold2_",
    "max_corr_steps_",
    "num_successes_before_increase_"
  };
  Float *ptr = (Float *) &settings;
  for (int i=0; i < 9; ++i)
    std::cerr << names[i] << " = " << *ptr++ << std::endl;
  std::cerr << names[9] << " = " << settings.max_corr_steps_ << std::endl;
  std::cerr << names[10] << " = " << settings.num_successes_before_increase_ << std::endl;
  std::cerr << "---------------------------------------------------------------\n";
  #endif 
}

void
process_args(int argc, char **argv)
{
  settings_ = M::DEFAULT;
  --argc; ++argv;
  // switches that can show up only in 1st position
  
  enum {INITIAL_ARGS, AFTER_INITIAL_ARGS, IMAGE_DATA, MAX_CORR_STEPS, EPSILON} argstate = INITIAL_ARGS;
  bool incomplete = false;
  std::string arg;
  if (argc) {
    arg = std::string(*argv);
    if (arg == "-h" || arg == "--help")
      print_usage();
    if (arg == "-g" || arg == "--profile") {
      profile_ = true;
      argstate = AFTER_INITIAL_ARGS;
      --argc; ++argv;
    } else if (arg == "-i" || arg == "--image_data") {
      image_data_ = true; 
      argstate = IMAGE_DATA;
      --argc; ++argv;
    } else if (arg[0] != '-') {
      if (argc == 2) {
          input_ = argv[1];
          output_ = argv[2];
          stdio_ = false;
      } else {
          std::cerr << "minus: \033[1;91m error\e[m\n";
          print_usage();
      }
    }
    
    while (argc) { // second and beyond: above switches must already be set
      arg = std::string(*argv);
      LOG("parsing arg " + arg);
      
      // argstate >= AFTER_INITIAL_ARGS ----------------------------------------
      if (argstate == IMAGE_DATA) {
        if (arg == "-gt") {
          ground_truth_ = true;
          --argc; ++argv;
          argstate = IMAGE_DATA;
          continue;
        }
        if (arg == "-AB") {
          two_problems_given_ = true;
          --argc; ++argv;
          argstate = IMAGE_DATA;
          continue;
        }
        argstate = AFTER_INITIAL_ARGS;
        continue;
      }
      
      if (argstate == MAX_CORR_STEPS) {
        settings_.max_corr_steps_ = std::stoi(arg);
        --argc; ++argv;
        argstate = AFTER_INITIAL_ARGS;
        incomplete = false;
        continue;
      }
      
      if (argstate == EPSILON) {
        settings_.epsilon_ = std::stod(arg);
        settings_.epsilon2_ = settings_.epsilon_*settings_.epsilon_;
        --argc; ++argv;
        argstate = AFTER_INITIAL_ARGS;
        incomplete = false;
        continue;
      }

      // argstate == AFTER_INITIAL_ARGS ----------------------------------------
      if (arg == "--max_corr_steps") {
        --argc; ++argv;
        argstate = MAX_CORR_STEPS;
        incomplete = true;
        continue;
      }
      if (arg == "--epsilon") {
        --argc; ++argv;
        argstate = EPSILON;
        incomplete = true;
        continue;
      }
      std::cerr << "minus: \033[1;91m error\e[m\n - unrecognized argument " << arg << std::endl;;
      print_usage();
    }

    if (incomplete) {
      std::cerr << "LOG \033[1;91merror: argument incomplete, (internal state: " << argstate << ")\n";
      print_usage();
    }
  }
}

// Simplest possible command to compute the Chicago problem
// for estimating calibrated trifocal geometry from points and lines at points
//
// This is to be kept very simple C with only minimal C++ with Templates.
// If you want to complicate this, please create another executable.
// 
int
main(int argc, char **argv)
{
  std::istream *inp = &std::cin;
  
  process_args(argc, argv);

  if (image_data_) {
    LOG("param: input is image pixel data");
    if (ground_truth_)
      LOG("param: reading ground truth appended to input pixel data");
  }
  if (profile_)
    LOG("Running default solve for profiling");
  if (stdio_)
    LOG("reading from stdio");
  else
    LOG("reading from " << input_ << " writing to " << output_);

  print_settings(settings_);

  if (!profile_) { // read files: either stdio or physical
    init_input(input_, inp);
    if (image_data_) {  // read image pixel-based I/O parameters
      if (!iread<Float>(*inp))
        return 1;
      data::params_ = data::params_start_target_;
    } else {  // read raw I/O homotopy parameters (to be used as engine)
      if (!mread<Float>(*inp))  // reads into global params_
        return 1;
    }
  }
  
  static M::solution solutions[M::nsols];
  {
    LOG("\033[0;33mUsing 4 threads by default\e[m\n");
    #ifdef M_VERBOSE
    if (two_problems_given_)
      std::cerr 
        << "LOG \033[0;33mContinuing between two problems A -> B by internal 0->A then A->B\e[m\n" 
        << "LOG \033[0;33mStarting path tracker from random initial solution to first problem 0->A\e[m\n" 
        << std::endl;
    else
      std::cerr << "LOG \033[0;33mStarting path tracker from random initial solution to given problem\e[m\n" << std::endl;
    #endif 
    std::thread t[4];
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    //  unsigned retval = 
    //  ptrack(&MINUS_DEFAULT, start_sols_, params_, solutions);
    {
      t[0] = std::thread(M::track, settings_, data::start_sols_, data::params_, solutions, 0, 78);
      t[1] = std::thread(M::track, settings_, data::start_sols_, data::params_, solutions, 78, 78*2);
      t[2] = std::thread(M::track, settings_, data::start_sols_, data::params_, solutions, 78*2, 78*3);
      t[3] = std::thread(M::track, settings_, data::start_sols_, data::params_, solutions, 78*3, 78*4);
      t[0].join(); t[1].join(); t[2].join(); t[3].join();
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t2 - t1).count();
    #ifdef M_VERBOSE
    print_num_steps(solutions);
    std::cerr << "LOG \033[1;32mTime of solver: " << duration << "ms\e[m" << std::endl;
    #endif
  }

  if (two_problems_given_) {
    // First, lets make sure problem A was properly solved
    if (!io::has_valid_solutions(solutions)) { // if no ground-truth is provided, it will return error
      LOG("\033[1;91mFAIL:\e[m  no valid solutions in problem A");
      return SOLVER_FAILURE;                    // if it can detect that the solver failed by generic tests
    }
    
    // 
    // Continue between two problems A and B by continuing from an internal
    // problem R to A (to discover all solutions of A), then from A to B.
    //
    // format solutions (of A) to be similar to data::start_sols_
    complex sols_A_matrix[M::nsols][M::nve];
    io::solutions_struct2vector(solutions, sols_A_matrix);
    const complex * const sols_A = (complex *) sols_A_matrix;
    
    // reset solutions
    static const M::solution s0;
    for (unsigned s=0; s < M::nsols; ++s)
      solutions[s] = s0;

    // generate homotopy params_ -----------------------------------------------
    //
    // At this point:
    // params_start_target_ = [ P0gammified, PAgammified]
    //    
    // We want
    // params_start_target_ = [ PAgammified, PBbammified]
    //
    // First we do params_start_target  = [ PAgammified, PAgammified]
    memcpy(data::params_start_target_, 
           data::params_start_target_+M::f::nparams, M::f::nparams*sizeof(complex));
    
    LOG("\033[0;33mReading second problem B\e[m\n");
    // Now read problem B & extract parameters into 2nd half of
    // params_start_target_
    if (image_data_) {  // read image pixel-based I/O parameters
      if (!iread<Float>(*inp))
        return 1;
      data::params_ = data::params_start_target_;
    } else {  // read raw I/O homotopy parameters (to be used as engine)
      std::cerr << "When continuing from A to B, non-pixel input not implemented\n";
      return 1;
    }
    
    // Homotopy-continue from A to B ---------------------------------------
    LOG("\033[0;33mUsing 4 threads by default\e[m\n");
    #ifdef M_VERBOSE
    std::cerr << "LOG \033[0;33mStarting path tracker from A to B\e[m\n" << std::endl;
    #endif 
    std::thread t[4];
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    //  unsigned retval = 
    //  ptrack(&MINUS_DEFAULT, start_sols_, params_, solutions);
    {
      t[0] = std::thread(M::track, settings_, sols_A, data::params_, solutions, 0, 78);
      t[1] = std::thread(M::track, settings_, sols_A, data::params_, solutions, 78, 78*2);
      t[2] = std::thread(M::track, settings_, sols_A, data::params_, solutions, 78*2, 78*3);
      t[3] = std::thread(M::track, settings_, sols_A, data::params_, solutions, 78*3, 78*4);
      t[0].join(); t[1].join(); t[2].join(); t[3].join();
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t2 - t1).count();
    #ifdef M_VERBOSE
    print_num_steps(solutions);
    std::cerr << "LOG \033[1;32mTime of solver A -> B: " << duration << "ms\e[m" << std::endl;
    #endif
  }
  
  if (profile_) {
    // compare solutions to certain values from M2
    // two random entries
    if (std::abs(solutions[1].x[1] - complex(-.25177177692982444e1, -.84845195030295639)) <= tol &&
        std::abs(solutions[M::nsols-2].x[2] - complex(.7318330016224166, .10129116603501138)) <= tol)
      std::cerr << "LOG solutions look OK\n";
    else  {
      std::cerr << "LOG \033[1;91merror:\e[m solutions dont match m2. Errors: ";
      std::cerr << std::abs(solutions[1].x[2] - complex(-.25177177692982444e1, -.84845195030295639)) << ", "
          << std::abs(solutions[M::nsols-2].x[2] - complex(.7318330016224166, .10129116603501138)) << std::endl;
    }
  }
  
  if (!mwrite<Float>(solutions, output_)) return 2;

  // ---------------------------------------------------------------------------
  // test_final_solve_against_ground_truth(solutions);
  // optional: filter solutions using positive depth, etc.
  if (ground_truth_ || profile_) {
    // TODO(juliana) should we has_valid_solutions here? 
    io::RC_to_QT_format(data::cameras_gt_, data::cameras_gt_quat_);
    unsigned sol_id;
    bool found = io::probe_all_solutions(solutions, data::cameras_gt_quat_, &sol_id);
    if (found) {
      LOG("found solution at index: " << sol_id);
      LOG("number of iterations of solution: " << solutions[sol_id].num_steps);
      if (solutions[sol_id].status != M::REGULAR)
        LOG("PROBLEM found ground truth but it is not REGULAR: " << sol_id);
    } else {
      LOG("\033[1;91mFAIL:\e[m  ground-truth not found among solutions");
      return SOLVER_FAILURE; 
      // you can detect solver failure by checking this exit code.
      // if you use shell, see:
      // https://www.thegeekstuff.com/2010/03/bash-shell-exit-status
    }
  } else if (!io::has_valid_solutions(solutions)) { // if no ground-truth is provided, it will return error
    LOG("\033[1;91mFAIL:\e[m  no valid solutions");
    return SOLVER_FAILURE;                    // if it can detect that the solver failed by generic tests
  }
  return 0;
}
