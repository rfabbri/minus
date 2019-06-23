// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date February 2019
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

void
points2params(points, params)
{
  struct xy {
    double x, y;
  }
  
  double pts[9] {
    {x, y},
  }
  pts[0].x;
  
  // xxx
  // tangents for the first two point tracks
  double tgts[6] {
  }

  pts[0][0]; tgts[0][0];
  pts[0][1];
  pts[0][2];
  pts[1][0]; tgts[0][1];
}

void
test_world_to_camera()
{
  //  This is exactly the case in the slides in big-notes/trifocal/trifocal.key

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
  // Extracting this data from those files:
  // scripts/getlines.sh

  // points for frame 42
  // + sed -n '621p;3012p;3390p' frame_0042-pts-2d.txt
  286.7673976130331539 217.06531260627261304
  141.01103052308988595 270.45312297462106699
  239.89822517853363593 86.442049763307068133
  
  // points for frame 54
  // + sed -n '621p;3012p;3390p' frame_0054-pts-2d.txt
  257.04360648826406077 159.4404341695463927
  241.41513314836856807 447.15662243793082098
  123.95973916849976604 213.90676875312345828
  
  // points for frame 62
  // + sed -n '621p;3012p;3390p' frame_0062-pts-2d.txt
  295.57132984990698787 147.80261937455236421
  240.78946527513195974 410.13737156824942076
  375.60750199363729962 277.22372936832925916

  // tangents for frame 42
  // + sed -n '3012p;3390p' frame_0042-tgts-2d.txt
  0.9536809622336909209 -0.3008199166827579818
  0.0082601187924503903515 -0.99996588463683822035

  // tangents for frame 54
  // + sed -n '3012p;3390p' frame_0054-tgts-2d.txt
  0.18491347256048701331 -0.9827548054655455001
  -0.99542450475950383648 0.095551322985590561587
  
  // tangents for frame 62
  // + sed -n '3012p;3390p' frame_0062-tgts-2d.txt
  0.77931350598248894102 -0.62663423094599701724
  0.76492323888624347283 0.64412144709812224619


  // intrinsic params equal to every frame
  2584.9325098195013197 0 249.77137587221417903
  0 2584.7918606057692159 278.31267937919352562
  0 0 1

  // extrinsic params for each frame

  // extrinsics for frame 42
  // + cat frame_0042.extrinsic
  -0.097305153950172085242 -0.22322794404612877894 -0.96989741313794208821
  0.96072075769186959793 0.23341709945525662695 -0.15010690664274928263
  0.25989869710080021337 -0.94640675329312473618 0.1917470327448986267

  -295.2090359311167731 1074.0075457376335635 -236.40439390871563319
  // extrinsics for frame 54
  // + cat frame_0054.extrinsic
  0.81820480546085894158 0.07511341824191355987 -0.56999900940332604016
  0.54313649506229122466 -0.42609616539484057585 0.72349485524588375007
  -0.18853022052765816552 -0.90155423144469115648 -0.38943076883056443327

  194.82952402681169701 1020.3676638972305 431.76461692675769655
  // extrinsics for frame 62
  // + cat frame_0062.extrinsic
  -0.61853492444688140672 -0.60388598633423984374 -0.50272881631015808868
  -0.025677143402306701336 -0.62392573537516660132 0.78106168837247091918
  -0.78533765448129877473 0.4960225723146879373 0.37041379052099593361

  887.07508137499985423 -562.68690102473453862 -415.57529638919055515
  


  
  //  This is iccv-chicago-src/chicagoTestDataSph file 5linesCase1.txt
  //  which was used to generate the default run. It already inverts point  coords
  //
  //  5 lines:
  //0.879008519574604 0.476806063840702 0.038623735269322 0.894813299896884 -0.446440542880738 0.032208045798003 0.704559167216148 0.709645249326512 -0.033704562238054
  //-0.707245886556431 0.706967648445817 -0.027867353589865 0.740915644124991 0.671598100273408 -0.041627953405298 0.978562374828962 -0.205950670232225 0.014188247805869
  //-0.995986793114319 0.089500323696923 0.002805786591495 -0.159509552710649 0.987196385018730 0.016889155865539 0.590904880455181 -0.806741236242606 -0.029119177786994
  //0.383774765752577 0.923426731891359 0.018967747764466 0.996132185054309 0.087867342619014 -0.002506238319557 0.594754053692778 0.803907715857989 -0.038939157023577
  //0.998695342916704 0.051064782741988 0.007599111258801 -0.045862067346128 -0.998947781807808 -0.027129197306103 -0.610567512482729 0.791964211754958 0.030062848861502
  //Camera Motions:
  //0.456457606228917 0.889159884956647 0.032266897893194 -0.659449197822556 0.313742799357977 0.683148747596164 0.597304954949274 -0.333106821957916 0.729566060037158
  //-16.067155449351297 -772.430415639545799 304.233386926931757
  //-0.264893260894465 -0.442357100009250 0.856826561448759 -0.513078745490601 0.817000282830920 0.263174350535889 -0.816444585540554 -0.369906385353920 -0.443381895002389
  //-953.788605479449529 -308.108851792374651 1642.743688891120200


  // Steps to transform this to params:  minus/src-iccv-chicago/t
  //
  // // p0 is initial parameters
  // // sols is initial solutions - already in chicago-default.hxx start_sols_
  // (p0,sols) = readStartSys "startSys";
  // 
  // (pLines, x) = parseFile: from the above, returns 
  // pLines is a 15x3 matrix of line coefs
  //
  // solveChicago(p0, pLines /*targetLines*/, sols /*startSols*/);
  //    P0 := if (numcols p0 == 1) then p0 else transpose p0;
  //    P1 := lines2Parameters targetLines;
  //    P01 := (gammify P0)||(gammify P1);
  //    H01 := specialize(PH, P01); -- XXXX
  // 
  

  // XXX doing: prototype conversions in pseudocode to see if simplifications
  // emerge
}


void
test_end_user_interface()
{
  // static data for points and cams

  

  // Placeholder for M::solve(M::DEFAULT, start_sols_, points, cameras);
//  M::points2params(points, params);
//  
//  M::track_all(M::DEFAULT, start_sols_, params_, solutions);
//  M::solutions2cams(solutions, cameras);
}

void
test_minus()
{
  test_full_solve();
  test_end_user_interface();
}

TESTMAIN(test_minus);
