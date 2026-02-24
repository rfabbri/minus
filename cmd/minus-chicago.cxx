// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date February 2019-2026
// 
#include "minus-chicago.h"

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

  if (input_data_) {
    LOG("input is problem data (image pixel data) of the problem to be solved (target problem)"); // as opposed to start/target parameters
    if (ground_truth_)
      LOG("reading ground truth appended to input target problem data (in pixels)");
  }
  if (profile_)
    LOG("Running default solve for profiling");
  else if (stdio_) 
    LOG("reading from stdio");
  else
    LOG("reading from " << input_ << " writing to " << output_);

  print_all_settings(settings_, ssettings_);

  minus_cmd_io cmd;

  if (!profile_) { // Read files: either stdio or physical
    cmd.init_input(input_, inp);
    if (input_data_) {          // Read target problem data, which is then converted to
      if (!iread<Float>(*inp))  // the internally used problem parameters
        return 1;
      data::params_ = data::params_start_target_;
    } else {  // Read raw start+target homotopy parameters, possibly randomized. (To be used as engine)
      if (!cmd.mread<Float>(*inp))  // Reads into global params_
        return 1;
    }
  } // else, profile: the homotopy data is already hardcoded in data::params_
  
  alignas(64) static M::solution solutions[M::nsols];
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
    if (input_data_) {  // read image pixel-based I/O parameters
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
    // override default data::compare_to_hardcoded_gt(sols);
    // 
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
  
  if (!cmd.mwrite<Float>(solutions, output_)) return 2;

  // ---------------------------------------------------------------------------
  // test_final_solve_against_ground_truth(solutions);
  // optional: filter solutions using problem-specific inequalities and
  // additional information
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
