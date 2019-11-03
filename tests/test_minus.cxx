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


//  This is exactly the case in the slides in big-notes/trifocal/trifocal.key
// Extracted from the original synthcurves dataset
// synthcurves-multiview-3d-dataset/spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object
//
// Frame files: frame_..42, 54, 62
//
// Points:
// 
// 620
// 3011  tangents
// 3389  tangents
// 0-based ids. +1 to get file line
// 
// Extracting this data from those files: scripts/getlines.sh

// NOTE: The point order seems to match Hongyi,
// though he might have done an off-by-1 mistake for point indexing (3012 vs
// 3011)
double p[io::nviews][io::npoints][io::ncoords] = {
  // points for frame 42
  // + sed -n '3012p;3390p;621p' frame_0042-pts-2d.txt
  {  
    {141.01103052308988595 270.45312297462106699},
    {239.89822517853363593 86.442049763307068133},
    {286.7673976130331539 217.06531260627261304}
  },
  // points for frame 54
  // + sed -n '3012p;3390p;621p;' frame_0054-pts-2d.txt
  {
    {241.41513314836856807 447.15662243793082098},
    {123.95973916849976604 213.90676875312345828},
    {257.04360648826406077 159.4404341695463927}
  },
  // points for frame 62
  // + sed -n '3012p;3390p;621p' frame_0062-pts-2d.txt
  {
    {375.60750199363729962 277.22372936832925916},
    {295.57132984990698787 147.80261937455236421},
    {240.78946527513195974 410.13737156824942076}
  }
}


double tgt[io::nviews][io::ntangents][io::ncoords] = {
  // tangents for frame 42
  // + sed -n '3012p;3390p' frame_0042-tgts-2d.txt
  {
    {0.9536809622336909209 -0.3008199166827579818},
    {0.0082601187924503903515 -0.99996588463683822035}
  }
  // tangents for frame 54
  // + sed -n '3012p;3390p' frame_0054-tgts-2d.txt
  {
    {0.18491347256048701331 -0.9827548054655455001},
    {-0.99542450475950383648 0.095551322985590561587}
  }
  // tangents for frame 62
  // + sed -n '3012p;3390p' frame_0062-tgts-2d.txt
  {
    {0.77931350598248894102 -0.62663423094599701724},
    {0.76492323888624347283 0.64412144709812224619}
  }
}


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
  io::point_tangents2params(points_gt, tangents, params);
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
