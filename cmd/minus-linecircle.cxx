#include "minus-linecircle.h"

// Simplest possible command to compute the linecircle problem
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

  if (param_data_) {
    LOG("param: input is parameters");
    if (ground_truth_)
      LOG("param: reading ground truth appended to input pixel data");
  }
  if (profile_)
    LOG("Running default solve for profiling");
  if (stdio_)
    LOG("reading from stdio");
  else
    LOG("reading from " << input_ << " writing to " << output_);

  print_all_settings(settings_, ssettings_);

  if (!profile_) { // read files: either stdio or physical
    init_input(input_, inp);
    if (param_data_) {  // read image pixel-based I/O parameters
      if (!iread<Float>(*inp))
        return 1;
      data::params_ = data::params_start_target_;
    } else {  // read raw I/O homotopy parameters (to be used as engine)
      if (!mread<Float>(*inp))  // reads into global params_
        return 1;
    }
  }
  
  alignas(64) static M::solution solutions[M::nsols];
  {
    LOG("\033[0;33mUsing 4 threads by default\e[m\n");
    #ifdef M_VERBOSE
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

  
  if (profile_) {
    // compare solutions to certain hardcoded values from M2
    // two random entries
    if (std::abs(solutions[1].x[1] - solutions_gt_[0]) <= tol &&
        std::abs(solutions[M::nsols-2].x[2] - solutions_gt_[1]) <= tol)
      std::cerr << "LOG solutions look OK\n";
    else  {
      std::cerr << "LOG \033[1;91merror:\e[m solutions dont match m2. Errors: ";
      std::cerr << std::abs(solutions[1].x[2] - solutions_gt_[0]) << ", "
          << std::abs(solutions[M::nsols-2].x[2] - solutions_gt_[1]) << std::endl;
    }
  }
  
  if (!mwrite<Float>(solutions, output_)) return 2;

  // ---------------------------------------------------------------------------
  // test_final_solve_against_ground_truth(solutions);
  // optional: filter solutions using positive depth, etc.
  if (ground_truth_ || profile_) {
    // TODO(juliana) should we has_valid_solutions here? 
    unsigned sol_id;
    bool found = io::probe_all_solutions(solutions, data::sols_gt_, &sol_id);
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
