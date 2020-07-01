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
#include <minus/debug-common.h>
#include <minus/chicago14a-io.h>
#include <minus/chicago-default.h> 

using namespace MiNuS;

static Float plines_m2_[io::pp::nvislines][io::ncoords2d_h] = {
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

static complex params_target_m2_[M::f::nparams] = {
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

static void
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

static void 
test_quat()
{
  {
    Float p[4]  = { 1, 1, 1, 1 };
    Float q[4]  = { 1.0001, 1.0001, 1.0001, 0.9999 };
    Float r[4]  = { 1.0, -1.0, 1.0, 0.9};
    Float d[4];

    util::normalize_quat(p); util::normalize_quat(q);

    TEST("rotation error to itself is zero?", util::rotation_error(p,p) < eps_, true);
    TEST("rotation error p,q", util::rotation_error(p,q) < eps_, true);
    TEST("rotation error p,r", util::rotation_error(p,r) > eps_, true);
  }
  {
    Float q_id[4] = {1, 0, 0, 0};
    Float v[3] = {3, 5, 11};
    Float vt[3] = {};
    util::quat_transform(q_id, v, vt);
      
    Float d[3];
    d[0] = v[0] - vt[0];
    d[1] = v[1] - vt[1];
    d[2] = v[2] - vt[2];

    TEST("transform by identity quat", std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) < eps_, true);
    
    Float r[3][3] = {
      {1, 0, 0},
      {0, 1, 0},
      {0, 0, 1}
    };
    
    Float qr[4] = {-1, -1, -1, -1};
    util::rotm2quat((Float *)r, qr);
    print(qr, 4);
    
    TEST("identity rotm to quat", util::rotation_error(qr,q_id) < eps_, true);
  }
  
  {
    Float q[4] = {0, 1, 0, 0};
    Float v[3] = {3, 5, 11};
    Float v_gt[3] = {3, -5, -11};
    Float vt[3] = {};
    util::quat_transform(q, v, vt);
      
    Float d[3];
    d[0] = v_gt[0] - vt[0];
    d[1] = v_gt[1] - vt[1];
    d[2] = v_gt[2] - vt[2];

    TEST_NEAR("transform by 0 1 0 0 quat", std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]), 0,  eps_);
  }
  
  {
    Float q[4] = {0, 0, 1, 0};
    Float v[3] = {3, 5, 11};
    Float v_gt[3] = {-3, 5, -11};
    Float vt[3] = {};
    util::quat_transform(q, v, vt);
      
    Float d[3];
    d[0] = v_gt[0] - vt[0];
    d[1] = v_gt[1] - vt[1];
    d[2] = v_gt[2] - vt[2];

    TEST_NEAR("transform by 0 0 1 0 quat", std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]), 0,  eps_);
  }
  
  {
    Float q[4] = {0, 0, 0, 1};
    Float v[3] = {3, 5, 11};
    Float v_gt[3] = {-3, -5, 11};
    Float vt[3] = {};
    util::quat_transform(q, v, vt);
      
    Float d[3];
    d[0] = v_gt[0] - vt[0];
    d[1] = v_gt[1] - vt[1];
    d[2] = v_gt[2] - vt[2];

    TEST_NEAR("transform by 0 0 0 1 quat", std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]), 0,  eps_);
  }
  
  {
    Float q[4] = {-1, 1, -3, 2};
    Float v[3] = {3, 5, 11};
    Float vt_gt[3] = { 4.466666666666665, -7.666666666666665, -8.733333333333334};
    Float vt[3] = {};
    
    util::normalize_quat(q);
    util::quat_transform(q, v, vt);

    std::cout << "q norm" << std::endl;
    print(q, 4);
    
    std::cout << "vt" << std::endl;
    print(vt, 3);
      
    Float d[3];
    d[0] = vt_gt[0] - vt[0];
    d[1] = vt_gt[1] - vt[1];
    d[2] = vt_gt[2] - vt[2];

    TEST_NEAR("transform by -1 1 -3 2 quat", std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]), 0,  eps_);
  }

  {
    Float q_gt[4] = {0.790532489152731, -0.321383239619688, -0.178689068700286,  -0.489736065256140};
    Float r_gt[9] = {
         0.456457606228917,  0.889159884956647,  0.032266897893194,
        -0.659449197822556,  0.313742799357976,  0.683148747596164,
         0.597304954949274, -0.333106821957916,  0.729566060037158
    };
    
    Float q[4] = {};
    Float r[3][3] = {};

    util::quat2rotm(q_gt, (Float *)r);

    TEST("quat2rotm rot matches gt", same_matrices((Float *) r, (Float *) r_gt, 3, 3), true);
    
    util::rotm2quat((Float *)r_gt, q);
    TEST("rotm2quat quat matches gt", same_vectors(q, q_gt, 4), true);
  }

  {
    std::cout << "\n--------------------------------------------------------------\n";
    io::RC_to_QT_format(data::cameras_gt_, data::cameras_gt_quat_);
    
    Float t01_gt[3] = {
      -016.0671554493513,  -772.4304156395458,   304.2333869269318
    };

    Float r01_gt[9] = {
       0.456457606228917,   0.889159884956647,   0.032266897893194,
      -0.659449197822556,   0.313742799357977,   0.683148747596164,
       0.597304954949274,  -0.333106821957916,   0.729566060037158
    };
    
    Float r01_computed[9] = {};
    
    util::quat2rotm(data::cameras_gt_quat_, r01_computed);
    
    TEST("RC_to_QT: relative rotation", same_matrices(r01_computed, r01_gt, 3, 3), true);
    TEST("RC_to_QT: relative translation", same_vectors(data::cameras_gt_quat_ + 8, t01_gt, 3), true);
    
    std::cout << "gt_quat:" << std::endl;
    print(data::cameras_gt_quat_,M::nve);
    
    std::cout << "--------------------------------------------------------------\n";
  }

  {
    Float r0_computed[9] = {};
    Float q0[4];
    
    util::rotm2quat((Float *) data::cameras_gt_[0], q0);

    std::cout << "q0" << std::endl;
    print(q0,4);
    
    util::quat2rotm(q0, r0_computed);
    
    TEST("rotm2quat -> quat2rotm", same_matrices(r0_computed, (Float*) data::cameras_gt_[0], 3, 3), true);
  }
    
  {
    std::cout << "\n--------------------------------------------------------------\n";
    Float q0[4], q1[4];
    Float r01_gt[9] = {
      0.456457606228917,  0.889159884956647,  0.032266897893194,
     -0.659449197822556,  0.313742799357977,  0.683148747596164,
      0.597304954949274, -0.333106821957916,  0.729566060037158
    };
    Float r01_computed[9] = {};

    util::rotm2quat((Float *) data::cameras_gt_[0], q0);
    util::rotm2quat((Float *) data::cameras_gt_[1], q1);

    std::cout << "q0: " << std::endl;
    print(q0, 4);
    
    std::cout << "q1: " << std::endl;
    print(q1, 4);
    
    Float d[4] = {};
    util::dquat(q1, q0, d);
    TEST_NEAR("dquat d unit", std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2] + d[3]*d[3]), 1, eps_);
    TEST_NEAR("dquat q1 unit", std::sqrt(q1[0]*q1[0] + q1[1]*q1[1] + q1[2]*q1[2] + q1[3]*q1[3]), 1, eps_);
    TEST_NEAR("dquat q0 unit", std::sqrt(q0[0]*q0[0] + q0[1]*q0[1] + q0[2]*q0[2] + q0[3]*q0[3]), 1, eps_);
    
    std::cout << "d: " << std::endl;
    print(d, 4);
    
    util::quat2rotm(d, r01_computed);
    TEST("relative rotation through quaternion", same_matrices(r01_computed, r01_gt, 3, 3), true);
    std::cout << "--------------------------------------------------------------\n";
  }
}

static void
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

#if 0
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
#endif
}

static void
test_gamma()
{
  { // sanity check
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

  { // 1
    //   - take gammified params
    //   - divide by non-gammified params
    // 2
    //   - take params = {1}
    //   - see what comes

    Float plines[io::pp::nvislines][io::ncoords2d_h] = {};
    Float pn[io::pp::nviews][io::pp::npoints][io::ncoords2d];
    Float tn[io::pp::nviews][io::pp::npoints][io::ncoords2d];
    
    // see if uno minus  default_gammas_m2 is less than 1
    io::invert_intrinsics(data::K_, data::p_[0], pn[0], io::pp::npoints);
    io::invert_intrinsics(data::K_, data::p_[1], pn[1], io::pp::npoints);
    io::invert_intrinsics(data::K_, data::p_[2], pn[2], io::pp::npoints);
    // don't use all three, but just invert all anyways.
    io::invert_intrinsics_tgt(data::K_, data::tgt_[0], tn[0], io::pp::npoints);
    io::invert_intrinsics_tgt(data::K_, data::tgt_[1], tn[1], io::pp::npoints);
    io::invert_intrinsics_tgt(data::K_, data::tgt_[2], tn[2], io::pp::npoints);
    
    // 1 
    complex params[M::f::nparams] = {};
    io::point_tangents2lines(pn, tn, 0, 1, plines);
    io::lines2params(plines, params); // ungammified target params
    
    for (unsigned i=0; i < M::f::nparams; ++i)
      params[i] = data::default_params_start_target_gammified_[i+M::f::nparams] / params[i];

    std::cout << "Default test gammas: " << std::endl;
    print(params, M::f::nparams - 17/*pChart are random*/, true);

    // 2
    complex uno[M::f::nparams];
    for (unsigned i=0; i < M::f::nparams; ++i)
      uno[i] = 1;
    io::gammify(uno);
    std::cout << "Minus gammas: " << std::endl;
    print(uno, M::f::nparams - 17, true);

    //    since it is random, all we do is to see if the blocks of similar gamma
    //    match
    //    
    //    uno_m2 = map((CC_53)^56,(CC_53)^1,{{-.982889+.1842*ii}, {-.982889+.1842*ii}, {-.982889+.1842*ii}, {-.636398-.771361*ii}, {-.636398-.771361*ii},
    //       {-.636398-.771361*ii}, {.840622-.541623*ii}, {.840622-.541623*ii}, {.840622-.541623*ii}, {-.750486+.660886*ii}, {-.750486+.660886*ii},
    //       {-.750486+.660886*ii}, {.980272-.197651*ii}, {.980272-.197651*ii}, {.980272-.197651*ii}, {-.942395-.334503*ii}, {-.942395-.334503*ii},
    //       {-.942395-.334503*ii}, {.43596-.899966*ii}, {.43596-.899966*ii}, {.43596-.899966*ii}, {-.509073+.860723*ii}, {-.509073+.860723*ii},
    //       {-.509073+.860723*ii}, {.981808-.189876*ii}, {.981808-.189876*ii}, {.981808-.189876*ii}, {-.982889-.1842*ii}, {-.750486-.660886*ii},
    //       {-.636398+.771361*ii}, {.980272+.197651*ii}, {.840622+.541623*ii}, {-.942395+.334503*ii}, {-.982889-.1842*ii}, {.43596+.899966*ii},
    //       {-.636398+.771361*ii}, {-.509073-.860723*ii}, {.840622+.541623*ii}, {.981808+.189876*ii}, {.99768+.068084*ii}, {.99768+.068084*ii},
    //       {.99768+.068084*ii}, {.99768+.068084*ii}, {.99768+.068084*ii}, {.99768+.068084*ii}, {.99768+.068084*ii}, {-.643374-.765552*ii},
    //       {-.643374-.765552*ii}, {-.643374-.765552*ii}, {-.643374-.765552*ii}, {-.643374-.765552*ii}, {.606663-.794959*ii}, {.606663-.794959*ii},
    //       {.606663-.794959*ii}, {.606663-.794959*ii}, {.606663-.794959*ii}})
  }
}
    

static void
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

static void 
test_point_tangents2lines()
{
  { // sanity check
    Float plines[io::pp::nvislines][io::ncoords2d_h] = {};
    Float pn[io::pp::nviews][io::pp::npoints][io::ncoords2d];
    Float tn[io::pp::nviews][io::pp::npoints][io::ncoords2d];
    
    io::invert_intrinsics(data::K_, data::p_[0], pn[0], io::pp::npoints);
    io::invert_intrinsics(data::K_, data::p_[1], pn[1], io::pp::npoints);
    io::invert_intrinsics(data::K_, data::p_[2], pn[2], io::pp::npoints);

    std::cout << "first normalized point" << std::endl;
    print(pn[0][0], 2);
    std::cout << "secondnormalized point" << std::endl;
    print(pn[0][1], 2);
    
    // don't use all three, but just invert all anyways.
    io::invert_intrinsics_tgt(data::K_, data::tgt_[0], tn[0], io::pp::npoints);
    io::invert_intrinsics_tgt(data::K_, data::tgt_[1], tn[1], io::pp::npoints);
    io::invert_intrinsics_tgt(data::K_, data::tgt_[2], tn[2], io::pp::npoints);
      
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

static void
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

static void
test_io_shaping()
{
  test_point_tangents2lines();  // OK
  test_gamma();
  test_lines2params();  // OK
  test_get_params_start_target();
}

void
test_internals()
{
  test_quat();
  test_rand();
  test_io_shaping();
}

TESTMAIN(test_internals);
