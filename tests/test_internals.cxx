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
    TEST_NEAR("Dot product", std::abs(dot_gt - minus_3d<Float>::dot(a,b)), 0 , tol);
    // std::cout << "Result: " << dot_gt  << "computed: " << minus_array<M::nnn,Float>::dot(a,b) << std::endl;


    complex c_gt[3] = {{38, -40}, {-59, 25}, {17, 27}};
    complex c[3];
    minus_3d<Float>::cross(a,b,c);
    Float m = 0; 
    for (unsigned i=0; i < 3; ++i) {
      m += std::abs(c_gt[i] - c[i]);
    }
    TEST_NEAR("Cross product", m, 0 , tol);
  }
}

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
test_lines2params()
{
  Float plines[15][3] = {0};
  complex params[M::nparams] = {0};
  io::lines2params(plines, params);
  TEST("lines2params sanity check", params[0], complex(0));
}

template <typename F=double>
void
get_params_start_target(F plines[15][3], C<F> * __restrict__ params/*[static 2*M::nparams]*/)
{
  io::lines2params(plines, params);
  io::gammify(params);
  io::gammify(params+M::nparams);
}

void
test_get_params_start_target()
{
  std::ofstream log("log_test_get_params_start_target");
  
  {
  complex params[2*M::nparams]; // start-target param pairs, P01 in chicago.m2, like params_ 
  Float plines[15][3];
  get_params_start_target<Float>(plines, params);
  TEST("get_params_start_target sanity check", params[0], complex(0));
  }
  
  {
  complex params[2*M::nparams]; // start-target param pairs, P01 in chicago.m2, like params_ 
  Float plines[15][3];
  get_params_start_target<Float>(plines, params);
  TEST("get_params_start_target sanity check", params[0], complex(0));
  for (unsigned i=0; i < 2*M::nparams; ++i)
    log << params[i] << std::endl;;
  }
  
}

void
test_io_shaping()
{
  test_gamma();
  test_lines2params();
  test_get_params_start_target();
}

void
test_internals()
{
   test_rand();
   test_io_shaping();
}

TESTMAIN(test_internals);
