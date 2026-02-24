#ifndef minus_linecircle_h_
#define minus_linecircle_h_
// 
// \author Ricardo Fabbri
// \date February 2026
//
// Solves the line-circle intersection problem with the best formulation and
// solver so far in terms of speed and robustness.
// 
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <minus/minus.h>
#include <minus/linecircle2a-io.h>
#include <minus/linecircle-default.h>

using namespace MiNuS;
#define Float double
typedef minus_core<linecircle> M;
static constexpr Float tol = 1e-3;
typedef std::complex<Float> complex;
using namespace std::chrono;


#define  M_VERBOSE 1     // display verbose messages
// comment this out plus the printouts if you prefer no debug and no dependency
// on this header
#include <minus/debug-util.h>
#include "cmd-util.h"
typedef minus_cmd_io<Float> cmd;


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
               "  minus -h       # (or --help) : print this help message\n"
               // "  minus -r       # (or --real)  :  outputs only real solutions\n"
               // "  minus -AB      # (or --two_problems) : continue between 2 given problems\n"
            <<
  R"(
  Input format:

  areal aimag breal bimag creal cimag dreal dimag ereal eimag freal fimag   
      # Space and new-line separators allowed.
      # These are coeffcients for the target system you want to solve:
      #     a(x^2 + y^2) + bx + c = 0, 
      #               dx + ey + f = 0

  # GROUND TRUTH (optional) if -gt flag provided, pass the ground truth here
  # first solution (x0,y0)
  x0real x0imag y0real y0imag            
  # second solution (x1,y1)
  x1real x1imag y1real y1imag)";

  exit(1);
}

bool ground_truth_ = false;
bool input_data_ = false;
bool profile_ = false;   // run some default solves for profiling
M::track_settings settings_;
M::f::settings ssettings_;   // specific settings (formulation-specific)

// print specialized settings, if any
void
print_ssettings(const M::f::settings &ssettings) {
  #ifdef M_VERBOSE
  // print here
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
  LOG("reading parameters a b c d e f directly."); // for other problems you may
                                                   // have to read data and
                                                   // convert to params
  if (!read_block(in, (F *)(data::params_start_target_+M::f::nparams), 2*M::f::nparams /*complex numbers*/))
    return false;
  LOG("reading ground truth solutions");
  if (ground_truth_ && !read_block(in, (F *)data::gt_sols_, 2*data::n_gt_sols_ /*complex numbers */))
    return false;
}

void
process_args(minus_cmd_io<Float> &cmd, int argc, char **argv)
{
  settings_ = M::DEFAULT;
  --argc; ++argv;
  // switches that can show up only in 1st position
  
  enum {
    INITIAL_ARGS, AFTER_INITIAL_ARGS, 
    INPUT_DATA /* input is problem/user data representation */, 
    MAX_CORR_STEPS, EPSILON, FILTER_DEGENERACY
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
    } else if (arg == "-i" || arg == "--input_data") {
      // Input is just the data specifying the target system to be solved.
      // 
      // Without this flag, the input is the full homotopy parameters comprising
      // of the concatenaded start and end system parameters, possibly
      // randomized. This is used if you want to set or randomize the
      // parameters outside MINUS, possibly together with setting your own start
      // solutions, or if you want to debug MINUS bypassing parameter
      // construction from the user's input representation (e.g., image points).
      // 
      // See print_usage() 
      input_data_ = true; 
      argstate = INPUT_DATA;
      --argc; ++argv;
    } else if (arg[0] != '-') { // if not a flag then two file names
      if (argc == 2) {
          cmd.input_ = argv[1];
          cmd.output_ = argv[2];
          cmd.stdio_ = false;
      } else {
          std::cerr << "minus: \033[1;91m error\e[m\n";
          print_usage();
      }
    }
    
    while (argc) { // second and beyond: above switches must already be set
      arg = std::string(*argv);
      LOG("parsing arg " + arg);
      
      // argstate >= AFTER_INITIAL_ARGS ----------------------------------------
      if (argstate == INPUT_DATA) {
        if (arg == "-gt") {
          ground_truth_ = true;
          --argc; ++argv;
          argstate = INPUT_DATA;
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
        // Flag to discard degenerate data with cheap rules (collinear points, etc).
        // 
        // Use with caution, since homotopy continuation is excellent at solving
        // nearly degenerate systems. Often the concern is not that HC may fail,
        // but that it may be slow.
        // 
        std::cerr << "LOG \033[1;91merror: no prefilter degeneracy option implemented\n";
        print_usage();
        // ssettings_.prefilter_degeneracy_ = true;
        // argstate = AFTER_INITIAL_ARGS;
        // incomplete = false;
        // continue;
      }
      if (arg == "--prefilter_degeneracy=no") {
        print_usage();
        std::cerr << "LOG \033[1;91merror: no prefilter degeneracy option implemented\n";
        // --argc; ++argv;
        // ssettings_.prefilter_degeneracy_ = false;
        // argstate = AFTER_INITIAL_ARGS;
        // incomplete = false;
        // continue;
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
#endif  // minus_linecircle_h_
