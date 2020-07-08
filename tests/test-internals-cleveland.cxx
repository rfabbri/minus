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

// from PLMP/common.m2 with setRandomSeed 0
// varible yInd from encodey
// these are line coeffs multiplied by a random unit complex number
static Float plines_m2_[io::pp::nvislines][io::ncoords2d_h] = {
  {-.182283123647345,-.883468673051352},
  {-.138328815643857,+.116151008216494},
  {-.09467158260928,+.380351600587356},
  {.419747470991594,+.283856800884654},
  {-.446457989755368,+.522779323449015},
  {.201697065677576,+.479512995905775},
  {.218021454646735,+.241987687178185},
  {-.401476292438847,-.199589118152257},
  {.0669722635375493,+.82970132638757},
  {-.3464644645551,-.454311657883248},
  {.768447038961226,+.154238339346338},
  {-.202504372612249,+.135110896571111},
  {.566209343965705,+.0376781607063618},
  {-.03370955060354,+.550642594645581},
  {-.0235863595958664,+.610808822638324},
  {.509548617164191,-.0521603439732403},
  {-.324241065859168,+.441098529221508},
  {.407618403815177,+.521331526826217},
  {-.0825121737059541,-.957631139083466},
  {-.249954524242536,+.0951451280198519},
  {.0627185096936635,+.0259011894437367},
  {.318201746334172,-.436891632118813},
  {.408942305461928,+.132235310473869},
  {.463504382361441,+.555263047044862},
  {.410028526138918,+.147027764216868},
  {-.400567478525395,+.209300424197819},
  {.586954173423017,+.511354348941571},
  {-.260053331987163,-.0284556152307416},
  {-.409551778992493,+.135216140977863},
  {.862724192009962,+.0354040470583808},
  {.430519629907852,-.300141828932049},
  {.386435351066271,+.468880952952147},
  {.499001077289531,+.32616564358956},
  {-.464826396731534,-.40183328869152},
  {.137138074578468,-.0747352726186499},
  {.691416914061394,+.346434506766115}
};

// from PLMP/common.m2 with setRandomSeed 0
// output from encodey
static complex params_target_m2_[M::f::nparams] = {
  {-.16694512367336373, -.48016315075890303},
  {-.99157162029301263e-1, .50854514976583019e-1},
  {.40843596380086816, .74988960076388267},
  {.24005297935107414, -.55443473600880988},
  {.35579994981937468, .52664172553206179},
  {.45254063984093085, .16198317363546658},
  {-.32055954389973557 , .51298310725469221},
  {-.63763173411863083, -.97961198086125065e-1},
  {.45836599366709263, .88430207791748927e-1},
  {-.12271448292586033 , .42988808989281441},
  {-.35930990600489321, -.77521384557464845},
  {.92489100865958379e-1, .24803837837252718},
  {.42076533103869362 , -.51853029410746887},
  {-.11429224120858139, .51483346887061188},
  {.43086361121246536, .30053818684461087},
  {-.34229950262108566 , .65347171645133684},
  {-.23677249838081607, .32390509764976749},
  {.52594599113298623, .13494759146291802},
  {-.51110727522732258 , -.56481567464392135},
  {.16963877363442756, .41338937790927549},
  {.20953804942632459, .4197360139902635},
  {.26584142103519359, -.28562403901307354},
  {-.77053474409592271, -.25316335287360076},
  {.35896834045246645, .24713053768297635},
  {-.68370058123056654, -.26703762748309573e-1},
  {-.41144705510141732, -.23233994211259809},
  {.14328774543387834, .53669220083398772},
  {-.23774353167966386, -.35727122509457282},
  {-.5955536027299726, .35567311473951663},
  {.45177968652197747, .36130726734574298},
  {-.13434834793504769, -.47772067747921959},
  {.59130719950792443 , .54659366111760543},
  {.28302081264090906, -.15882021522509815},
  {-.75204866344735066, .30864929309597677},
  {-.35380812990953886, .2024173661172079},
  {.32854941449227759, .25507003159774794},
  {-.13165457056593666, -.29981297822303271},
  {.55958260236666435, .38157927002972025},
  {.78690558802403687, .48254060358968576},
  {.18341280875391097, .98478064286640177},
  {.26951487092715457e-1, .89315760170527014},
  {.48721832559873601, .41301632713026465},
  {.14922318657387967e-2, .16723995217160845},
  {.10127706933307021e1, -.10201803437563384e-2},
  {.94428842481176067, .24747728809068684},
  {.14129049569475999, .30135579844710858},
  {.65440350141506998, .20825075617745584},
  {.35475286173551113 , .303276361712708},
  {.66210753804900224, -.87105375587924472},
  {.29978728636557772, .44971878331292992},
  {.3152899304825087, .68788990969775671e-1},
  {.73625798204387827, .33603459149689929e-1},
  {.18904012029959902, .27838470658931337e-1}
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

    for (unsigned i=0; i < io::pp::nvislines*io::ncoords2d_h; ++i) 
      TEST("lines2params complex numbers same norm as m2 default run"
          std::norm(p1[i]), std::norm(params_target_m2_[i]), eps_);

    //    can't do this since the angle of each complex number
    //    is randomized in PLMP m2
//    TEST("lines2params matches m2 default run except for random pChart", 
//        same_vectors((Float *) params_target_m2_, (Float *) p1, 
//          M::f::nparams*2 - (7+5+5)*2), true);
    //    pChart: just unit rands 17x1
    //        sphere(7,1)|sphere(5,1)|sphere(5,1)

    Float nrm = minus_array<7,Float>::norm2(p1+M::f::nparams - (7+5+5));
    TEST_NEAR("lines2params random pChart unit", nrm, 1, eps_);
    nrm = minus_array<5,Float>::norm2(p1+M::f::nparams - (5+5));
    TEST_NEAR("lines2params random pChart unit", nrm, 1, eps_);
    nrm = minus_array<5,Float>::norm2(p1+M::f::nparams - (5));
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
