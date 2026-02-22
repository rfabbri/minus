#ifndef minus_chicago_h_
#define minus_chicago_h_
// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date February 2019-2026
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
#include "cmd-util.h"

// exit code. Conventions:
// 0 (Zero)	Success
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
 
  Input format (indexing goes _view_points_coords. any number of spaces and newlines optional. can be in
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

bool ground_truth_ = false;
bool two_problems_given_ = false;
bool reading_first_point_ = true;
bool image_data_ = false;
bool profile_ = false;   // run some default solves for profiling
M::track_settings settings_; // general homotopy settings
M::f::settings ssettings_;   // specific settings (formulation-specific)

// print specialized settings
void
print_ssettings(const M::f::settings &ssettings) {
  #ifdef M_VERBOSE
  std::cerr << "Formulation-specific settings ---------------------------------\n";
  {
  const char *names[10] = {
    "prefilter_degeneracy_",
    "prefilter_area_degeneracy_eps_",
    "prefilter_angle_degeneracy_eps_"
  };
  std::cerr << names[0] << " = " << ssettings_.prefilter_degeneracy_ << std::endl;
  Float *ptr = (Float *) &ssettings;
  for (int i=1; i < 3; ++i)
    std::cerr << names[i] << " = " << *ptr++ << std::endl;
  std::cerr << "---------------------------------------------------------------\n";
  }
  #endif 
}

void
print_all_settings(const M::track_settings &settings, const M::f::settings &ssettings)
{
  print_settings(settings);
  print_ssettings(ssettings);
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
    if (!io::point_tangents2params_img(ssettings_, data::p_, data::tgt_, tgt_ids[0], tgt_ids[1],
        data::K_, data::params_start_target_)) {
      LOG("Data configuration close to degenerate, discarding");
      return false;
    }
    reading_first_point_ = false;
  } else { // when reading second point B, do not gammify A again
    static constexpr bool gammify_target_problem = false;
    io::point_tangents2params_img(ssettings_, data::p_, data::tgt_, tgt_ids[0], tgt_ids[1],
        data::K_, data::params_start_target_, gammify_target_problem);
  }
  return true;
}

void
process_args(int argc, char **argv)
{
  settings_ = M::DEFAULT;
  ssettings_ = M::f::DEFAULT;
  --argc; ++argv;
  // switches that can show up only in 1st position
  
  enum {
    INITIAL_ARGS, AFTER_INITIAL_ARGS, IMAGE_DATA, MAX_CORR_STEPS, 
    EPSILON, FILTER_DEGENERACY
  } argstate = INITIAL_ARGS;
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
        double epsilon = std::stod(arg);
        settings_.epsilon2_ = epsilon*epsilon;
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
      if (arg == "--prefilter_degeneracy=yes") {
        --argc; ++argv;
        ssettings_.prefilter_degeneracy_ = true;
        argstate = AFTER_INITIAL_ARGS;
        incomplete = false;
        continue;
      }
      if (arg == "--prefilter_degeneracy=no") {
        --argc; ++argv;
        ssettings_.prefilter_degeneracy_ = false;
        argstate = AFTER_INITIAL_ARGS;
        incomplete = false;
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
#endif  // minus_chicago_h
