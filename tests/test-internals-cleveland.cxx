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
#include <minus/cleveland14a-io.h>
#include <minus/cleveland-default.h> 

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

# if 0
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
#endif

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

#if 0
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
#endif

static void
test_io_shaping()
{
  // test_gamma();
  test_lines2params();
  // test_get_params_start_target();
}

void
test_internals_cleveland()
{
  test_io_shaping();
}

TESTMAIN(test_internals_cleveland);
