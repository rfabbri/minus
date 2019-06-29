// 
// \author Ricardo Fabbri
// \date June 2019
//
// This tests minus internals
// 
// This tests using minus.hxx directly
// instead of a separately instantiated version in the Templates/ folder
// 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <testlib/testlib_test.h>
#include <minus.hxx>

#define Float double
typedef minus<chicago14a> M;
typedef minus_util<Float> util;
typedef minus_io<chicago14a> io;
typedef std::complex<Float> complex;
using namespace std::chrono;


void
test_rand()
{
  std::cout << "Random number generation \n";
  { // gaussian random --------------------------------------------------------
  static constexpr unsigned NRAND = 20000.;
  Float d[NRAND];
  for (unsigned i=0; i < NRAND; ++i)
    d[i] = util::gauss(util::rnd);
  double s = 0;
  for (unsigned i=0; i < NRAND; ++i)
    s+= d[i];
  #ifdef M_VERBOSE     // display verbose messages
  std::cout << "***** avg: " << s/NRAND << std::endl;
  #endif
  TEST("Is the random generator gaussian(0,1000) average not too far from 0?", s/NRAND < 100, true);
  }

  { // randc ------------------------------------------------------------------
  constexpr Float tol = 0.001;
  constexpr unsigned n = 1000;
  complex z[n] = {};
  // std::ofstream log("log");
  for (unsigned i=0; i < n; ++i) {
    util::randc(z+i);
    // log << z[i] << std::endl;
  }
  
  bool test = true;
  for (unsigned i=0; i < n; ++i) {
    Float mag = std::abs(z[i]);
    if (mag > 1 + tol || mag < 1 - tol) {
      TEST("Is randc unit magnitude", true, false);
      test = false;
    }
  }
  if (test == true) TEST("Is randc unit magnitude", true, true);
  }

  { // rand_sphere ------------------------------------------------------------
    constexpr unsigned n = 7;
    constexpr Float tol = 1e-6;
    complex p[7];
    util::rand_sphere(p,7);
    std::ofstream log("log");
    
    Float real_mag = 0;
    Float complex_sum = 0;
    for (unsigned i=0; i < n; ++i) {
      log << p[i] << std::endl;
      real_mag += p[i].real()*p[i].real();
      complex_sum += p[i].imag();
    }
    TEST_NEAR("rand_sphere complex sum", complex_sum, 0, tol);
    TEST_NEAR("rand_sphere total real magnitude", std::sqrt(real_mag), 1, tol);
  }

  { // dot and cross -----------------------------------------------------------
    constexpr Float tol = 1e-6;
    std::cout << "dot and cross products\n";
    complex a[3] = {{5,3}, {3,4}, {6,5}};
    complex b[3] = {{-2,4}, {2,5}, {3,-4}};
    // complex cross_gt = { { , }, { , }, { ,-4} };
    complex dot_gt = complex({5,3})*complex({-2,4})
       + complex({3,4})*complex({2,5})+complex({6,5})*complex({3,-4});
    TEST_NEAR("Dot product", std::abs(dot_gt - minus_array<M::nnn,Float>::dot(a,b)), 0 , tol);
    // std::cout << "Result: " << dot_gt  << "computed: " << minus_array<M::nnn,Float>::dot(a,b) << std::endl;


    complex c_gt[3] = {{38, -40}, {-59, 25}, {17, 27}};
    complex c[3];
    minus_array<M::nnn,Float>::cross(a,b,c);
    Float m = 0; 
    for (unsigned i=0; i < 3; ++i) {
      m += std::abs(c_gt[i] - c[i]);
    }
    TEST_NEAR("Cross product", m, 0 , tol);
  }
}
#if 0
// XXX
                     // we only use the first half of the outer
                     // 2*M::nparams array 
                     // after this fn, complex part zero, but we will use this space later
                     // to gammify/randomize
lines2params(double plines[][3], complex params[static M::nparams])
{
  //    params (P1) is pDouble||pTriple||pChart  //  Hongyi: [pDouble; tripleChart; XR'; XT1'; XT2'];
  //    size              27       12        17 = 56

  // pDouble ----------------------------------------
  // converts 1st 9 lines to complex (real part zero)
  // 
  // 9x3 out of the 15x3 of the pairsiwe lines, linearized as 27x1
  // Tim: pDouble is matrix(targetLines^{0..8},27,1);
  // Order: row-major
  const double *pl = (const double *)plines;
  for (unsigned i=0; i < 27; ++i) params[i] = pl[i];

  // pTriple ----------------------------------------
  
  unsigned triple_intersections[6][3] = 
    {{0,3,9},{0+1,3+1,9+1},{0+2,3+2,9+2},{0,6,12},{0+1,6+1,12+1},{0+2,6+2,12+2}};

  complex (*params_lines)[2] = (complex (*)[2]) (params+27);
  // express each of the 6 tangents in the basis of the other pairwise lines
  // intersecting at the same point
  for (unsigned l=0; l < 6; ++l) {
    const complex *l0 = plines[triple_intersections[l][0]];
    const complex *l1 = plines[triple_intersections[l][1]];
    const complex *l2 = plines[triple_intersections[l][2]];
    l0l0 = dot(l0,l0); l0l1 = dot(l0,l1); l1l1 = dot(l1,l1);
    l2l0 = dot(l2,l0); l2l1 = dot(l2,l1);
    // cross([l0l0 l1l0 l2l0], [l0l1 l1l1 l2l1], l2_l0l1);
    {
      complex v1[3], v2[3];
      v1[0] = l0l0; v1[1] = l1l0; v1[2] = l2l0;
      v2[0] = l0l1; v2[1] = l1l1; v2[2] = l2l1;
      cross(v1, v2, l2_l0l1);
    }
    
    params_lines[l][0] = l2_l0l1[0]/l2_l0l1[2]; // divide by the last coord (see cross prod formula, plug direct)
    params_lines[l][1] = l2_l0l1[1]/l2_l0l1[2];
  }
  //        
  //    pChart: just unit rands 17x1
  //        sphere(7,1)|sphere(5,1)|sphere(5,1)
  //
  rand_sphere(params+27+12,7);
  rand_sphere(params+27+12+7,5);
  rand_sphere(params+27+12+7+5,5);
}
#endif

void
test_gamma()
{
  std::ofstream log("log");
  // sanity check
  complex p_test[M::nparams];
  
  // arbitrary values
  for (unsigned i=0; i < M::nparams; ++i)
    p_test[i] = 100;
  
  io::gammify(p_test);
  for (unsigned i=0; i < M::nparams; ++i)
    log << p_test[i] << std::endl;
}

void
test_io_shaping()
{
  test_gamma();
}

void
test_internals()
{
   test_rand();
   test_io_shaping();
}

TESTMAIN(test_internals);
