// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date February 2019
//
// Tests more comprehensive runs of minus using the public interface
// 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <testlib/testlib_test.h>
#include <minus/minus.h>

#define Float double
typedef minus<chicago14a> M;
typedef minus_io<chicago14a> io;
static constexpr Float tol = 1e-3;
typedef std::complex<Float> complex;
using namespace std::chrono;

// Start solutions hardcoded for efficiency.
// If you want to play with different start sols,
// write another program that accepts start sols in runtime,
// but keep this one lean & mean.
#include <minus/chicago14a-default.hxx> 
// We include it separately so they don't clutter this app,
// neither minus.h, and can be reused by other progs
// TODO(developer note): make this part of Minus' template as a specialization. 
// But for efficiency I chose to do it outside.
// Perhaps a minus class should be written that wraps the lean minus_core.
// And in _that_ one, we put these default vectors depending on template tag.

#define  M_VERBOSE 1     // display verbose messages


void
test_against_ground_truth(const M::solution solutions[])
{
  // compare solutions to certain values from Macaulay2
  // two random entries. Just a sanity check against the original code prototype.
  // Not a Full comparison to ground truth cameras!
  bool ok=false;
  if (std::abs(solutions[1].x[1] - complex(-.25177177692982444e1, -.84845195030295639)) <= tol &&
      std::abs(solutions[M::nsols-2].x[2] - complex(.7318330016224166, .10129116603501138)) <= tol) {
    std::cerr << "LOG solutions look OK\n"; ok = true;
  } else  {
    std::cerr << "LOG \033[1;91merror:\e[m solutions dont match original code. Errors: ";
    std::cerr << std::abs(solutions[1].x[2] - complex(-.25177177692982444e1, -.84845195030295639)) << ", "
        << std::abs(solutions[M::nsols-2].x[2] - complex(.7318330016224166, .10129116603501138)) << std::endl;
  }
  TEST("Solutions match original code", ok, true);
}

void
test_full_solve()
{
  std::cerr << "Starting path tracker" << std::endl;
  M::solution solutions[M::nsols];
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  M::track_all(M::DEFAULT, start_sols_, params_, solutions); // <<<<<<<<< MEAT
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(t2 - t1).count();
  std::cerr << "LOG \033[1;32mTime of solver: " << duration << "ms\e[m" << std::endl;
  TEST("Did it track first solution?", solutions[0].t > 0, true);
  TEST("Did it track the last solution?", solutions[M::nsols-1].t > 0, true);
  test_against_ground_truth(solutions);
}

// fill in internal format for cameras_gt_
// into camera_gt_quaternion
// - convert each rotation to quaternion in the right order
// - make the rotations all relative to the first camera
static void
minus_initialize_gt()
{
  //  R01 = R1 * inv(R0);
  //  T01 = R1 * (C0 - C1);
  //  R12 = R2 * inv(R1);
  //  T12 = R2 * (C1 - C2);

  // rotation-center format (used in synthcurves dataset)
  // relative to the world, to internal quaternion-translation format relative
  // to first camera
  io::RC_to_QT_format(cameras_gt_, cameras_gt_quat_);
}

void
test_end_user_interface()
{
  // static data for points and cams

  // M::solve(M::DEFAULT, start_sols_, points, cameras);
  {
  std::cerr << "Starting path tracker in test_end_user_interface" << std::endl;
  M::solution solutions[M::nsols];
  io::point_tangents2params_img(p_, tgt_, 0, 1, K_, params_start_target_);
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  M::track_all(M::DEFAULT, start_sols_, params_start_target_, solutions);
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(t2 - t1).count();
  std::cerr << "LOG \033[1;32mTime of solver: " << duration << "ms\e[m" << std::endl;
  TEST("IO: Did it track first solution?", solutions[0].t > 0, true);
  TEST("IO: Did it track the last solution?", solutions[M::nsols-1].t > 0, true);
  test_against_ground_truth(solutions);
  // test solutions are good before formatting them out,
  // just like in minus-chicago

  Float cameras[M::nsols][2/*2nd and 3rd cams relative to 1st*/][4][3] = {};
  unsigned nsols_final = 0;
  unsigned id_sols[M::nsols] = {};
  io::all_solutions2cams(solutions, cameras, id_sols, &nsols_final);
  // to output

  // ---------------------------------------------------------------------------
  // test_final_solve_against_ground_truth(solutions);
  // optional: filter solutions using positive depth, etc.
  unsigned sol_id;
  minus_initialize_gt();
  bool found = io::probe_solutions(solutions, cameras_gt_quat_, &sol_id);
  TEST("IO: Found GT solution? ", found, true);
  if (found)
    std::cout << "found solution at index: " << sol_id << std::endl;
  }
}

void
test_minus()
{
  test_full_solve();
  test_end_user_interface();
}

TESTMAIN(test_minus);
