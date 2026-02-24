// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date February 2019-2026
// 
#include "minus-linecircle.h"

// Simplest possible command to compute the linecircle problem
// of intersecting a line and a circle 
//
// This is to be kept very simple C with only minimal C++ with Templates.
// If you want to complicate this, please create another executable.
// 
int
main(int argc, char **argv)
{
  std::istream *inp = &std::cin;
  cmd c;
  
  process_args(c, argc, argv);

  if (input_data_) {
    LOG("input is problem data (image pixel data) of the problem to be solved (target problem)"); // as opposed to start/target parameters
    if (ground_truth_)
      LOG("reading ground truth appended to input target problem data");
  }
  if (profile_)
    LOG("Running default solve for profiling");
  else if (c.stdio_)
    LOG("reading from stdio");
  else
    LOG("reading from " << c.input_ << " writing to " << c.output_);

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

  if (profile_) { // data::compare_to_hardcoded_gt(sols);
    // compare solutions to certain hardcoded values from M2
    bool solutions_bad = false;
    for (unsigned s=0; s < data::n_gt_sols_; ++s)
      for (unsigned v=0; v < M::nve; ++v)
        if (std::abs(solutions[data::gt_sols_id_[s]].x[v] - data::gt_sols_[s][v]) > tol) {
          solutions_bad = true;
          goto not_ok;
        }
    std::cerr << "LOG solutions look OK\n";
    not_ok: 
    if (solutions_bad) {
      std::cerr << "LOG \033[1;91merror:\e[m solutions dont match hardcoded ground-truth exactly (could be a normalization issue). Errors: ";
      for (unsigned s=0; s < data::n_gt_sols_; ++s) {
        std::cerr << "Solution id " << s << ", errors as pairs (variable id, error): " << std::endl;
        for (unsigned v=0; v < M::nve; ++v) {
          std::cerr << "\t" << v << "\t" << std::abs(solutions[data::gt_sols_id_[s]].x[v] - data::gt_sols_[s][v]) << std::endl;
        }
      }
    }
  }
  
  if (!c.mwrite(solutions, output_)) return 2;

  // ---------------------------------------------------------------------------
  // test_final_solve_against_ground_truth(solutions);
  // optional: filter solutions using problem-specific inequalities and
  // additional information
  if (ground_truth_ || profile_) {
    unsigned sol_id;
    bool found = io::probe_all_solutions(solutions, data::gt_sols_, &sol_id);
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
