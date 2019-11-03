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
typedef minus_io<chicago14a> io;
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
  M::track_all(M::DEFAULT, start_sols_, params_, solutions); // <<<<<<<<< MEAT
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(t2 - t1).count();
  std::cerr << "LOG \033[1;32mTime of solver: " << duration << "ms\e[m" << std::endl;
  TEST("Did it track first solution?", solutions[0].t > 0, true);
  TEST("Did it track the last solution?", solutions[M::nsols-1].t > 0, true);
  test_against_ground_truth(solutions);
}

//
// returns cameras[0:nsols_final][2][4][3]
//
// where the camera matrix P^t = [R|T]^t is cameras[sol_number][view_id][:][:]
// where view_id is 0 or 1 for second and third camera relative to the first,
// resp.
//
// This design is for cache speed. Translation in the camera matrix is stored
// such that its coordinates are memory contiguous.
// 
// The cameras array is fixed in size to NSOLS which is the max
// number of solutions, which perfectly fits in memory. The caller must pass an
// array with that minimum.
void
solutions2cams(M::solution raw_solutions[M::NSOLS], double cameras[M::NSOLS][2][4][3], 
    unsigned id_sols[M::NSOLS], unsigned *nsols_final)
{
  *nsols_final = 0;
  for (unsigned sol=0; sol < M::NSOLS; ++sol)
    double real_solutions[M::NNN];
    if (get_real(raw_solutions[sol], real_solutions)) {
      id_sols[(*nsols_final)++] = sol;
      // build cams by using quat2rotm
      solution2cams(real_solutions, (double [2][4][3] ) (cameras + sol));
    }
}

void 
solution2cams(F rs[NNN], double cameras[2][4][3])
{
  // camera 0 (2nd camera relative to 1st)
  quat2rotm(rs, (double [3][3]) cameras[0]);
  cameras[0][3][0] = rs[8];
  cameras[0][3][1] = rs[9];
  cameras[0][3][2] = rs[10];
  
  // camera 1 (3rd camera relative to 1st)
  quat2rotm(rs, (double [3][3]) cameras[1]);
  cameras[1][3][0] = rs[8];
  cameras[1][3][1] = rs[9];
  cameras[1][3][2] = rs[10];

  // quat12 rs(0:3), quat12 rs(4:7)
  //  T12 = solutions(9:11);
  //  T13 = solutions(12:14);
  //  R12 = quat2rotm(transpose(quat12));
  //  R13 = quat2rotm(transpose(quat13));
}

void
test_end_user_interface()
{
  // static data for points and cams

  // M::solve(M::DEFAULT, start_sols_, points, cameras);
  {
  io::point_tangents2params(points_gt, tangents, params);
  M::track_all(M::DEFAULT, start_sols_, params_, solutions);
  
  TEST("Did it track first solution?", solutions[0].t > 0, true);
  TEST("Did it track the last solution?", solutions[M::nsols-1].t > 0, true);
  test_against_ground_truth(solutions);

  // test solutions are good before formatting them out,
  // just like in minus-chicago

//  io::solutions2cams(solutions, cameras, &nsols_final);
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
