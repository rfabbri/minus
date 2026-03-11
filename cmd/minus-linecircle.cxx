// UNDER DEVELOPMENT -----------------------------------------------------------
// 
// FINAL TESTS TAKINNG PLACE march 2026
// 
//------------------------------------------------------------------------------
//
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date February 2019-2026
// 
#include "minus-linecircle.h"

// returns command exit status, see minu-problem.h for generic status
int
find_ground_truth(M::solution solutions[M::nsols])
{
  typedef minus_array<M::f::nve,Float> v;
  unsigned sol_id;
  // searches for ground-truth among solutions, possibly with normalizations and specific rules
  Float real_gt_sol[M::nve];
  
  if (v::get_real(data::gt_sols_[0], real_gt_sol)) {
    // real ground-truth: we run a custom, faster and more complete solution matcher 
    // PRO: just enconde your specific ground-truth as real to begin with
    bool found = io::probe_all_solutions(solutions, real_gt_sol, &sol_id);
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
  } else {
    LOG("WARNING: ground-truth is _not_ real, this is not the intended use for MINUS, only for debugging");
    // in the non-real case we run a generic solution matcher
    // PRO: remove this and keep only the real-specific part
    if (ground_truth_) {
      if (io::probe_solutions(solutions, data::gt_sols_[0])) {
        LOG("found complex ground-truth solution, even though MINUS is intended for real ground-truth");
      } else {
        LOG("\033[1;91mFAIL:\e[m  complex ground-truth not found among solutions");
        return SOLVER_FAILURE; 
      }
    }
  }
  return 0;
}

void
run_solver(M::solution solutions[M::nsols])
{
  LOG("\033[0;33mUsing 4 threads by default\e[m\n");
  #ifdef M_VERBOSE
  std::cerr << "LOG \033[0;33mStarting path tracker from random initial solution to given problem\e[m\n" << std::endl;
  #endif 
//  std::thread t[4];
  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  M::track(settings_, data::start_sols_, data::params_, solutions, 0, M::nsols);
  //  unsigned retval = 
  //  ptrack(&MINUS_DEFAULT, start_sols_, params_, solutions);
  {
//    t[0] = std::thread(M::track, settings_, data::start_sols_, data::params_, solutions, 0, 78);
//    t[1] = std::thread(M::track, settings_, data::start_sols_, data::params_, solutions, 78, 78*2);
//    t[2] = std::thread(M::track, settings_, data::start_sols_, data::params_, solutions, 78*2, 78*3);
//    t[3] = std::thread(M::track, settings_, data::start_sols_, data::params_, solutions, 78*3, 78*4);
//    t[0].join(); t[1].join(); t[2].join(); t[3].join();
  }
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(t2 - t1).count();
  #ifdef M_VERBOSE
  print_num_steps(solutions);
  std::cerr << "LOG \033[1;32mTime of solver: " << duration << "ms\e[m" << std::endl;
  #endif
}

// exact means the exact numerical value, without normalizing, for trying to
// reproduce a specific numeric behavior
void
probe_solutions_with_gt_exact(M::solution solutions[M::nsols])
{
  if (io::probe_solutions(solutions, data::gt_sols_[0])) // find exactly without normalizing
    std::cerr << "LOG solutions look OK\n";
  else
    std::cerr << "LOG \033[1;91merror:\e[m solutions dont match hardcoded ground-truth numerically.\n";
}


// Simplest possible command to compute the linecircle problem
// of intersecting a line and a circle 
//
// This is to be kept very simple C with only minimal C++ with Templates.
// If you want to complicate this, please create another executable.
// 
int
main(int argc, char **argv)
{
  std::cout << " UNDER DEVELOPMENT -----------------------------------------------------------\n";
  std::istream *inp = &std::cin;
  cmd c;
  
  process_args(c, argc, argv);
  print_all_settings(settings_, ssettings_);

  if (!profile_) { // Read files: either stdio or physical
    c.init_input(c.input_, inp);
    if (input_data_) {          // Read target problem data, which is then converted to
      if (!iread<Float>(*inp))  // the internally used problem parameters
        return 1;
      data::params_ = data::params_start_target_;
    } else {  // Read raw start+target homotopy parameters, possibly randomized. (To be used as engine)
      if (!c.mread(*inp))  // Reads into global params_
        return 1;
    }
  } // else, profile: the homotopy data is already hardcoded in data::params_

  alignas(64) static M::solution solutions[M::nsols];

  run_solver(solutions);

  if (profile_) 
    probe_solutions_with_gt_exact(solutions);
  
  if (!c.mwrite(solutions, c.output_)) return 2;

  if (ground_truth_ || profile_)
    return find_ground_truth(solutions);
  else if (!io::has_valid_solutions(solutions)) {   // if no ground-truth is provided, it will return error if
    LOG("\033[1;91mFAIL:\e[m  no valid solutions"); // it can detect that the solver failed by generic tests
    return SOLVER_FAILURE;                          // without using ground-truth, e.g., no real roots
  }                                                 // or problem-specific inequalities
  std::cout << " UNDER DEVELOPMENT -----------------------------------------------------------\n";
  return 0;
}
