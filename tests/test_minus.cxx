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
#include <minus.h>

#define Float double
typedef minus<chicago14a> M;
static constexpr Float tol = 1e-3;
typedef std::complex<Float> complex;
using namespace std::chrono;

// Start solutions hardcoded for efficiency.
// If you want to play with different start sols,
// write another program that accepts start sols in runtime,
// but keep this one lean & mean.
#include <chicago14a-default.hxx> 
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
  
  static M::solution solutions[M::nsols];
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  M::track_all(M::DEFAULT, start_sols_, params_, solutions);
  
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(t2 - t1).count();
  std::cerr << "LOG \033[1;32mTime of solver: " << duration << "ms\e[m" << std::endl;
  
  TEST("Did it track first solution?", solutions[0].t > 0, true);
  TEST("Did it track the last solution?", solutions[M::nsols-1].t > 0, true);

  test_against_ground_truth(solutions);
}

/*

void
solutions2cams()
{
  // for each solution
    get_real();
    // build cams by using quat2rotm
}
*/

void
test_end_user_interface()
{
  // static data for points and cams

  // M::solve(M::DEFAULT, start_sols_, points, cameras);
  {
  io::point_tangents2params(points, tangents, params);
  M::track_all(M::DEFAULT, start_sols_, params_, solutions);
  
  TEST("Did it track first solution?", solutions[0].t > 0, true);
  TEST("Did it track the last solution?", solutions[M::nsols-1].t > 0, true);
  test_against_ground_truth(solutions);

  // test solutions are good before formatting them out,
  // just like in minus-chicago

//  io::solutions2cams(solutions, cameras);
//  test_final_solve_against_ground_truth(solutions);
  }
}

void
test_minus()
{
  test_full_solve();
  test_end_user_interface();
}

TESTMAIN(test_minus);
