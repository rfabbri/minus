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
#include <minus/minus.hxx>

#define Float double
typedef minus<chicago14a> M;
typedef minus_util<Float> util;
typedef minus_io<chicago14a> io;
typedef std::complex<Float> complex;
using namespace std::chrono;

// Start solutions hardcoded for efficiency.
// If you want to play with different start sols,
// write another program that accepts start sols in runtime,
// but keep this one lean & mean.
#include <minus/chicago14a-default.hxx> 

#if 0
template <unsigned N, typename F /* no complex number, only float types */>
struct minus_array_util { // Not speed critical -----------------------------------------
  static inline void 
  multiply_scalar_to_self(F *__restrict__ a, F b)
  {
    for (unsigned i = 0; i < N; ++i, ++a) *a = *a * b;
  }

  static inline void
  negate_self(F * __restrict__ a)
  {
    for (unsigned i = 0; i < N; ++i, ++a) *a = -*a;
  }

  static inline void 
  multiply_self(F * __restrict__ a, const F * __restrict__ b)
  {
    for (unsigned int i=0; i < N; ++i,++a,++b) *a *= *b;
  }

  static inline void 
  add_to_self(F * __restrict__ a, const F * __restrict__ b)
  {
    for (unsigned int i=0; i < N; ++i,++a,++b) *a += *b;
  }

  static inline void 
  add_scalar_to_self(F * __restrict__ a, F b)
  {
    for (unsigned int i=0; i < N; ++i,++a) *a += b;
  }

  static inline void 
  copy(const F * __restrict__ a, F * __restrict__ b)
  {
    memcpy(b, a, N*sizeof(F));
  }

  static inline F
  norm2(const F *__restrict__ a)
  {
    F val = 0;
    F const* __restrict__ end = a+N;
    while (a != end) val += std::abs(*a++);
    return val;
  }
};

typedef minus_array_util<Float> mu; 
#endif


Float eps_ = 1e-3;

Float plines_m2_[io::pp::nvislines][io::ncoords2d_h] = {
       {.879009, .476806, .0386237},
       {.894813, -.446441, .032208},
       {.704559, .709645, -.0337046},
       {-.707246, .706968, -.0278674},
       {.740916, .671598, -.041628},
       {.978562, -.205951, .0141882},
       {-.995987, .0895003, .00280579},
       {-.15951, .987196, .0168892},
       {.590905, -.806741, -.0291192},
       {.383775, .923427, .0189677},
       {.996132, .0878673, -.00250624},
       {.594754, .803908, -.0389392},
       {.998695, .0510648, .00759911},
       {-.0458621, -.998948, -.0271292},
       {-.610568, .791964, .0300628}
};

complex params_target_m2_[M::f::nparams] = {
 {.879009},
 {.476806},
 {.0386237},
 {.894813},
 {-.446441},
 {.032208},
 {.704559},
 {.709645},
 {-.0337046},
 {-.707246},
 {.706968},
 {-.0278674},
 {.740916},
 {.671598},
 {-.041628},
 {.978562},
 {-.205951},
 {.0141882},
 {-.995987},
 {.0895003},
 {.00280579},
 {-.15951},
 {.987196},
 {.0168892},
 {.590905},
 {-.806741},
 {-.0291192},
 {-.964279},
 {-.655833},
 {-.648147},
 {-.561685},
 {-1.08294},
 {.171924},
 {-.253346},
 {.779129},
 {.251946},
 {1.12584},
 {.0249},
 {1.00359},
 {.0203933},
 {.205885},
 {.272273},
 {.049322},
 {-.613848},
 {-.121431},
 {-.699339},
 {.141343},
 {.538675},
 {-.477745},
 {-.606832},
 {-.305558},
 {.128332},
 {.424345},
 {-.608765},
 {.0000133112},
 {.657926}
};

void
print(const Float *v, unsigned n)
{
  for (unsigned i=0; i < n; ++i)
    std::cout << v[i] << "  ";
  std::cout << std::endl;
}

void
print(const Float *v, unsigned nrows, unsigned ncols)
{
  for (unsigned i=0; i < nrows; ++i)
    print(v + i*ncols, ncols);
}

// test if two vectors of the same size are equal up to tolerance
bool
same_vectors(const Float *v, const Float *w, unsigned n, Float tol=eps_)
{
  for (unsigned i=0; i < n; ++i)
    if (std::fabs(v[i] - w[i]) > eps_) {
      std::cout << "v: \n";
      print(v, n);
      std::cout << "w: \n";
      print(w, n);
      printf("offending element i, v[i], w[i] = [%d], %g, %g\n", i, v[i], w[i]);
      return false;
    }
  return true;
}

// test if two vectors of the same size are equal up to tolerance
//bool
//same_homog_vectors(const Float *v, const Float *w, unsigned n, Float tol=eps_)
//{
//  for (unsigned i=0; i < n; ++i)
//  for (unsigned i=0; i < n; ++i)
//    if (std::fabs(v[i] - w[i]) > eps_) {
//      std::cout << "v: \n";
//      print(v, n);
//      std::cout << "w: \n";
//      print(w, n);
//      printf("offending element i, v[i] = [%d], %g\n", i, v[i]);
//      return false;
//    }
//  return true;
//}

// test if two vectors of the same size are equal up to tolerance
bool
same_matrices(const Float *vp, const Float *wp, unsigned nrows, unsigned ncols, Float tol=eps_)
{
  const Float (*v)[ncols] = (Float (*)[ncols]) vp, (*w)[ncols] = (Float (*)[ncols]) wp;
  for (unsigned i=0; i < nrows; ++i)
    for (unsigned j=0; j < ncols; ++j)
      if (std::fabs(v[i][j] - w[i][j]) > eps_) {
        std::cout << "v: \n";
        print((Float *) v, nrows, ncols);
        std::cout << "w: \n";
        print((Float *) w, nrows, ncols); 
        printf("offending element [i][j] v[i][j], w[i,j] = [%d][%d], %g, %g\n", i, j, v[i][j], w[i][j]);
        return false;
      }
  return true;
}

// test if two vectors of the same size are equal up to tolerance
bool
same_matrices_up_to_row_scale(const Float *vp, const Float *wp, unsigned nrows, unsigned ncols, Float tol=eps_)
{
  const Float (*v)[ncols] = (Float (*)[ncols]) vp, (*w)[ncols] = (Float (*)[ncols]) wp;
  for (unsigned i=0; i < nrows; ++i) {
    bool match = true;
    for (unsigned j=0; j < ncols; ++j) 
      if (std::fabs(v[i][j] - w[i][j]) > eps_)  {
        match = false;
        break;
      }
    if (!match) {
      for (unsigned j=0; j < ncols; ++j) 
        if (std::fabs(v[i][j] + w[i][j]) > eps_ && std::fabs(v[i][j] - w[i][j]) > eps_)  {
            std::cout << "v: \n";
            print((Float *) v, nrows, ncols);
            std::cout << "w: \n";
            print((Float *) w, nrows, ncols); 
            printf("offending element [i][j] v[i][j], w[i,j] = [%d][%d], %g, %g\n", i, j, v[i][j], w[i][j]);
            return false;
        }
    }
  }
  return true;
}

void
test_cross2()
{
  Float v[3] = {1, 0, 1};
  Float s[3] = {-1, 0, 1};
  Float r[3] = {-1};
  // Float gt[3] = {0, 2, 0};
  minus_3d<Float>::cross2(v,s,r);
  TEST("Cross product, test 2", std::fabs(minus_3d<Float>::dot(r, v)) > eps_, false);
  TEST("Cross product, test 3", std::fabs(minus_3d<Float>::dot(r, s)) > eps_, false);

  minus_3d<Float>::cross(v,s,r);
  TEST("Cross product, test 4", std::fabs(minus_3d<Float>::dot(r, v)) > eps_, false);
  TEST("Cross product, test 5", std::fabs(minus_3d<Float>::dot(r, s)) > eps_, false);
}

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
    // std::cout << "Result: " << dot_gt  << "computed: " << minus_array<M::nve,Float>::dot(a,b) << std::endl;


    complex c_gt[3] = {{38, -40}, {-59, 25}, {17, 27}};
    complex c[3];
    minus_3d<Float>::cross(a,b,c);
    Float m = 0; 
    for (unsigned i=0; i < 3; ++i) {
      m += std::abs(c_gt[i] - c[i]);
    }
    TEST_NEAR("Cross product", m, 0 , tol);
    test_cross2();
  }
}

void
test_gamma()
{
  {
  std::ofstream log("log");
  // sanity check
  complex p_test[M::f::nparams];
  
  // arbitrary values
  for (unsigned i=0; i < M::f::nparams; ++i)
    p_test[i] = 100;
  
  io::gammify(p_test);
  for (unsigned i=0; i < M::f::nparams; ++i)
    log << p_test[i] << std::endl;
  }
}

void
test_lines2params()
{
  { // sanity checks
  Float plines[io::pp::nvislines][io::ncoords2d_h] = {};
  complex params[M::f::nparams] = {};
  io::lines2params(plines, params);
  TEST("lines2params sanity check", params[0], complex(0));
  }

  { // does it match default from macaulay?
    // - define plines just as in hongy's 5-line file to macaulay
    // - see if the ungammified parameters match the harcoded one

    complex p1[M::f::nparams] = {};
    io::lines2params(plines_m2_, p1);
    TEST("lines2params matches m2 default run except for random pChart", same_vectors((Float *) params_target_m2_, (Float *) p1, M::f::nparams*2 - (7+5+5)*2), true);
    //    pChart: just unit rands 17x1
    //        sphere(7,1)|sphere(5,1)|sphere(5,1)

    Float nrm = minus_array<7,Float>::norm2(params_target_m2_+M::f::nparams - (7+5+5));
    TEST_NEAR("lines2params random pChart unit", nrm, 1, eps_);
    nrm = minus_array<5,Float>::norm2(params_target_m2_+M::f::nparams - (5+5));
    TEST_NEAR("lines2params random pChart unit", nrm, 1, eps_);
    nrm = minus_array<5,Float>::norm2(params_target_m2_+M::f::nparams - (5));
    TEST_NEAR("lines2params random pChart unit", nrm, 1, eps_);
  }
}


void 
test_point_tangents2lines()
{
  { // sanity check
    Float plines[io::pp::nvislines][io::ncoords2d_h] = {};
    Float pn[io::pp::nviews][io::pp::npoints][io::ncoords2d];
    Float tn[io::pp::nviews][io::pp::npoints][io::ncoords2d];
    
    io::invert_intrinsics(K_, p_[0], pn[0], io::pp::npoints);
    io::invert_intrinsics(K_, p_[1], pn[1], io::pp::npoints);
    io::invert_intrinsics(K_, p_[2], pn[2], io::pp::npoints);

    std::cout << "first normalized point" << std::endl;
    print(pn[0][0], 2);
    std::cout << "secondnormalized point" << std::endl;
    print(pn[0][1], 2);
    
    // don't use all three, but just invert all anyways.
    io::invert_intrinsics_tgt(K_, tgt_[0], tn[0], io::pp::npoints);
    io::invert_intrinsics_tgt(K_, tgt_[1], tn[1], io::pp::npoints);
    io::invert_intrinsics_tgt(K_, tgt_[2], tn[2], io::pp::npoints);
      
    io::point_tangents2lines(pn, tn, 0, 1, plines);
    
    TEST("point_tangents2lines sanity check", plines[1][1] != complex(0), true);
    TEST("point_tangents2lines sanity check", plines[2][0] != complex(0) && plines[2][1] != complex(0), true);

    // does it match macaulay2? 
    // actually, we are superior. more than twice the precision.
    TEST("point_tangents2lines's lines match m2", same_matrices_up_to_row_scale((Float *) plines, (Float *) plines_m2_, io::pp::nvislines, io::ncoords2d_h), true);
  }

  { // hardcoded simple input points and desired output lines
  }
}

void
test_get_params_start_target()
{
  std::ofstream log("log_test_get_params_start_target");
  
  {
  complex params[2*M::f::nparams]; // start-target param pairs, P01 in chicago.m2, like params_ 
  Float plines[io::pp::nvislines][io::ncoords2d_h] = {};
  io::get_params_start_target(plines, params);
  TEST("get_params_start_target sanity check", params[0], complex(0));
  }
  
  {
  complex params[2*M::f::nparams]; // start-target param pairs, P01 in chicago.m2, like params_ 
  Float plines[io::pp::nvislines][io::ncoords2d_h] = {};
  io::get_params_start_target(plines, params);
  TEST("get_params_start_target sanity check", params[0], complex(0));
  for (unsigned i=0; i < 2*M::f::nparams; ++i)
    log << params[i] << std::endl;;
  }
}

void
test_io_shaping()
{
  test_gamma();
  test_point_tangents2lines();
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
