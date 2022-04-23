#ifndef minus_h_
#define minus_h_
// 
// \brief MInimal problem NUmerical continuation package
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date Created: Fri Feb  8 17:42:49 EST 2019
//
// \verbatim
// Modifications
//    Leykin Feb82019: Initial sketch as simplified code from Macaulay e/NAG.*
//    Tim    Feb2019:  Chicago-specific prototype in Macaulay2
//    Fabbri Mar162019: fully templated and optimized code
// \endverbatim
//
// OPTIMIZATIONS
//  - see Fabbri's trifocal.key in bignotes for basic results
//  - see CMakeLists.txt and README.md

#include <complex>
#include <random>

template <typename F=double>
using C = typename std::complex<F>;

// The problem solvers that this solver template currently supports
enum problem {chicago14a, chicago6a, phoenix10a /*, standard*/};

// each problem specializes this in their specific .h
template <problem P>
struct formulation_parameters;

// each problem specializes this in their specific .h
template <problem P>
struct problem_parameters;

// problem specific definitions that must be available before anything, at compile time
#include <minus/parameters.h>

template <problem P, typename F=double>
class minus_core { // fully static, not to be instantiated - just used for templating
  public: // ----------- Data structures --------------------------------------
  
  // Tracker parameters
  // Default values obtained in M2 by doing 
  // eg DEFAULT#tStep
  // or peek DEFAULT
  // 
  // We use underscore in case we want to make setters/getters with same name,
  // or members of Tracker class if more complete C++ desired
  struct track_settings; 
  typedef formulation_parameters<P> f;
  /* General content of formulation_parameters (each problem may add to this):
  template <problem P, typename F>
  struct minus_core<P, F>::formulation_parameters {
    // Specific values defined for each problem:
    static constexpr unsigned nve = NVE;          // the size of the system (Number of Variables or Equations)
    static constexpr unsigned nsols = NSOLS;      // the number of solutions
    static constexpr unsigned nparams = NPARAMS;  // the number of parameters
  }
  */
  // shortcuts to the formulation parameters useful to most users
  static constexpr unsigned nsols = f::nsols;
  static constexpr unsigned nve = f::nve;
  
  enum solution_status {
    UNDETERMINED,       // 0
    PROCESSING,         // 1
    REGULAR,            // 2 OK. rest is error.
    SINGULAR,           // 3 unused
    INFINITY_FAILED,    // 4
    MIN_STEP_FAILED,    // 5
    ORIGIN_FAILED,      // 6 unused
    INCREASE_PRECISION, // 7 unused
    DECREASE_PRECISION  // 8 unused
  };
  
  struct solution
  {
    C<F> x[f::nve];    // array of n coordinates
    F t;            // last value of parameter t used
    solution_status status;
    //  unsigned num_steps;  // number of steps taken along the path
    solution() : status(UNDETERMINED) { }
  };

  static const track_settings DEFAULT;
  
  public: // ----------- Functions --------------------------------------------
  
  ///// THE MEAT /////
  static void track(const track_settings &s, const C<F> s_sols[f::nve*f::nsols], 
      const C<F> params[2*f::nparams], solution raw_solutions[f::nsols], unsigned sol_min, unsigned sol_max);
  
  // helper function: tracks all, no begin or end to specify
  static void track_all(const track_settings &s, const C<F> s_sols[f::nve*f::nsols], 
      const C<F> params[2*f::nparams], solution raw_solutions[f::nsols])
  { track(s, s_sols, params, raw_solutions, 0, f::nsols); }
  
  private: // -----------------------------------------------------------------
  static constexpr unsigned NVEPLUS1 = f::nve+1;
  static constexpr unsigned NVEPLUS2 = f::nve+2;
  static constexpr unsigned NVE2 = f::nve*f::nve;
  static void evaluate_Hxt(const C<F> * __restrict x /*x, t*/,    const C<F> * __restrict params, C<F> * __restrict y /*HxH*/);
  static void evaluate_HxH(const C<F> * __restrict x /*x and t*/, const C<F> * __restrict params, C<F> * __restrict y /*HxH*/);
};

// TODO: make these static
template <problem P, typename F>
struct minus_core<P, F>::track_settings {
  track_settings():
    init_dt_(0.05),   // m2 tStep, t_step, raw interface code initDt
    min_dt_(1e-7),        // m2 tStepMin, raw interface code minDt
    end_zone_factor_(0.05),
    epsilon_(.0000001), // m2 CorrectorTolerance (chicago.m2, track.m2), raw interface code epsilon (interface2.d, NAG.cpp::rawSwetParametersPT)
    epsilon2_(epsilon_ * epsilon_), 
    dt_increase_factor_(2.),  // m2 stepIncreaseFactor
    dt_decrease_factor_(1./dt_increase_factor_),  // m2 stepDecreaseFactor not existent in DEFAULT, using what is in track.m2:77 
    infinity_threshold_(1e7), // m2 InfinityThreshold
    infinity_threshold2_(infinity_threshold_ * infinity_threshold_),
    max_corr_steps_(3),  // m2 maxCorrSteps (track.m2 param of rawSetParametersPT corresp to max_corr_steps in NAG.cpp)
    num_successes_before_increase_(20) // m2 numberSuccessesBeforeIncrease
  { }
  
  F init_dt_;   // m2 tStep, t_step, raw interface code initDt
  F min_dt_;        // m2 tStepMin, raw interface code minDt
  F end_zone_factor_;
  F epsilon_; // m2 CorrectorTolerance (chicago.m2, track.m2), raw interface code epsilon (interface2.d, NAG.cpp:rawSwetParametersPT)
  F epsilon2_; 
  F dt_increase_factor_;  // m2 stepIncreaseFactor
  F dt_decrease_factor_;  // m2 stepDecreaseFactor not existent in DEFAULT, using what is in track.m2:77 
  F infinity_threshold_; // m2 InfinityThreshold
  F infinity_threshold2_;
  unsigned max_corr_steps_;  // m2 maxCorrSteps (track.m2 param of rawSetParametersPT corresp to max_corr_steps in NAG.cpp)
  unsigned num_successes_before_increase_; // m2 numberSuccessesBeforeIncrease
};
// Original settings from Tim: Fri Feb 22 12:00:06 -03 2019 Git 0ec3340
// o9 = MutableHashTable{AffinePatches => DynamicPatch     }
//                      Attempts => 5
//                      Bits => infinity
//                      CorrectorTolerance => .000001
//                      EndZoneFactor => .05
//                      ErrorTolerance => 1e-8
//                      Field => CC
//                      gamma => 1
//                      InfinityThreshold => 1e7
//                      Iterations => 30
//                      maxCorrSteps => 3
//                      maxNumberOfVariables => 50
//                      MultistepDegree => 3
//                      NoOutput => true
//                      Normalize => false
//                      numberSuccessesBeforeIncrease => 2
//                      Precision => 53
//                      Predictor => RungeKutta4
//                      Projectivize => false
//                      ResidualTolerance => .0001
//                      SingularConditionNumber => 100000
//                      SLP => false
//                      SLPcorrector => false
//                      SLPpredictor => false
//                      Software => M2engine
//                      stepIncreaseFactor => 2
//                      tDegree => 1
//                      Tolerance => .000001
//                      tStep => .05
//                      tStepMin => 1e-7

// The user specializes this to their problem inside problem.hxx
// Needed to create this class since functions do not always support partial
// specialization
template <problem P, typename F>
struct eval {
  static void Hxt(const C<F> * __restrict x /*x, t*/,    const C<F> * __restrict params, C<F> * __restrict y /*HxH*/);
  static void HxH(const C<F> * __restrict x /*x and t*/, const C<F> * __restrict params, C<F> * __restrict y /*HxH*/);
};

template <problem P, typename F>
void minus_core<P, F>::evaluate_Hxt(const C<F> * __restrict x /*x, t*/, const C<F> * __restrict params, C<F> * __restrict y /*HxH*/)
{
  eval<P,F>::Hxt(x, params, y);
}

template <problem P, typename F>
void minus_core<P, F>::evaluate_HxH(const C<F> * __restrict x /*x, t*/, const C<F> * __restrict params, C<F> * __restrict y /*HxH*/)
{
  eval<P,F>::HxH(x, params, y);
}

// IO shaping: not used in tracker, but only for shaping user data
// The user specializes this to their problem inside problem.hxx
// This is a base template class for zero and first-order problems (points,
// lines, tangents or lines at points).
//
// This is the default I/O class for problems involving 0 or more points, lines,
// tangents (or lines at points with orientation). For zero points, for example, 
// some useless functions will be generated by the current template, but 
// for now this won't affect code performance.
// 
// We suggest creating another class minus_io_shaping_second_order, etc, for
// problems that may not involve points, tangents nor lines, eg., problems
// involving only conics or higher-order features.
// 
// Note: Needed to create this class since functions do not always support partial
// specialization
// https://stackoverflow.com/questions/1501357/template-specialization-of-particular-members
// 
// TODO: leave this class empty.
// Just use a mold / inherictance to reuse basic structure
// useful for all problems. Right now, this is somewhat specific to certain
// types of problems (points and lines) and formulations (minors-based /
// visible-line based). Feel free to ignore this in the specialization to your
// problem, except the constants and typedefs.
template <problem P, typename F=double>
struct minus_io_shaping {
  typedef minus_core<P, F> M;
  typedef struct M::solution solution;

  // cast to this to interpret real M::solution::x order
  // internal note: this order is eg  in parser.m2 l 68
  struct solution_shape {
    F q01[4];
    F q02[4];
    F t01[3];
    F t02[3];
  };
  
  static constexpr unsigned ncoords2d = 2;  // just a documented name for the number of inhomog coordinates
  static constexpr unsigned ncoords2d_h = 3;// just a name for the usual number of homog coordinates in P^2
  static constexpr unsigned ncoords3d = 3;  // just a documented name for the number of inhomog 3D coordinates
  typedef problem_parameters<P> pp;
  
  // shortcuts to the problem parameters
  static constexpr unsigned  nviews = pp::nviews;
  static constexpr unsigned  npoints = pp::npoints;
#if 0
  { // The basic structure, defined at each problem-specific .hxx
  // unsigned NVIEWS, unsigned NPOINTS /* per view*/, unsigned NFREELINES, unsigned NTANGENTS, 
           
  static constexpr unsigned nviews = NVIEWS; 
  static constexpr unsigned npoints = NPOINTS;
  static constexpr unsigned nfreelines = NFREELINES;
  // even though Chicago needs only 2 tangents, api assumes 3 tangents are given,
  // out of which two are selected by indexing. This is the most common use
  // case, where all features naturally have tangents. If strictly 2 tangents
  // are to be passed, you can leave the unused one as zeros throughout the API.
  static constexpr unsigned ntangents = NTANGENTS;
  // number of lines connecting each pair of points plus going through points
  // plus the number of free lines in the first order problem.
  // Note that tangent orientation may help ruling out solutions; this is why
  // we call it tangent, and not general lines at points. There is more
  // information at tangents which can be used as part of the model for curves.
  // The tangent orientation can be constrained by running without orientation
  // mattering at first, and then propagating these to neighboring features
  // along curves
  // for formulations based on all lines -- not all formulations use this
  static constexpr unsigned nvislines = ( (npoints*(npoints-1) >> 1) + ntangents + nfreelines ) * nviews; 
  }
#endif
  
  // nvislines = 15 for Chicago.
  // INPUT ---------------------------------------------------------------------
  static void point_tangents2params(const F p[pp::nviews][pp::npoints][ncoords2d], const F tgt[pp::nviews][pp::npoints][ncoords2d], unsigned id_tgt0, unsigned id_tgt1, C<F> * __restrict params/*[static 2*M::nparams]*/);
  static void point_tangents2params_img(const F p[pp::nviews][pp::npoints][ncoords2d], const F tgt[pp::nviews][pp::npoints][ncoords2d], unsigned id_tgt0, unsigned id_tgt1, const F K[/*3 or 2*/][ncoords2d_h], C<F> * __restrict params/*[static 2*M::nparams]*/);
  // this function is the same for all problems
  static void get_params_start_target(F plines[/*15 for chicago*/][ncoords2d_h], C<F> * __restrict params/*[static 2*M::nparams]*/);
  static void gammify(C<F> * __restrict params/*[ chicago: M::nparams]*/);
  static void point_tangents2lines(F p[pp::nviews][pp::npoints][ncoords2d], F tgt[pp::nviews][pp::npoints][ncoords2d], unsigned id_tgt0, unsigned id_tgt1, F plines[pp::nvislines][ncoords2d_h]);
  static void lines2params(const F plines[pp::nvislines][ncoords2d_h], C<F> * __restrict params/*[static M::n//params]*/);
  static void invert_intrinsics(const F K[/*3 or 2 ignoring last line*/][ncoords2d_h], const double pix_coords[][ncoords2d], double normalized_coords[][ncoords2d], unsigned npts);
  static void invert_intrinsics_tgt(const F K[/*3 or 2 ignoring last line*/][ncoords2d_h], const double pix_tgt_coords[][ncoords2d], double normalized_tgt_coords[][ncoords2d], unsigned npts);
  static void normalize_line(F line[ncoords2d_h]);
  static void normalize_lines(F lines[][ncoords2d_h], unsigned nlines);
  static void initialize_gt();
  static void RC_to_QT_format(const F rc[M::nviews-1][4][3], F qt[M::nve]);
  static void rotation_error(const F p[4], const F q[4]);

  // OUTPUT --------------------------------------------------------------------
  static void all_solutions2cams(solution raw_solutions[M::nsols], F cameras[M::nsols][2][4][3], unsigned id_sols[M::nsols], unsigned *nsols_final);
  static void solution2cams(F rs[M::f::nve], F cameras[2][4][3]);
  static bool probe_solutions(const typename M::solution solutions[M::nsols], solution_shape *probe_cameras,
      unsigned *solution_index);
  static bool probe_solutions(const typename M::solution solutions[M::nsols], F probe_cameras[M::nve],
      unsigned *solution_index);
  static bool probe_all_solutions(const typename M::solution solutions[M::nsols], solution_shape *probe_cameras,
      unsigned *solution_index);
  static bool probe_all_solutions(const typename M::solution solutions[M::nsols], F probe_cameras[M::nve],
      unsigned *solution_index);
  static bool has_valid_solutions(const typename M::solution solutions[M::nsols]);
};

// Shortcuts and aliases -------------------------------------------------------
// type alias used to hide a template parameter 
template<problem P, typename F=double>
using minus = minus_core<P, F>;  
// can now use minus<chicago14a>
// no need to do this:
// typedef minus<double, 312, 14, 56> minus_chicago14a;

template<problem P, typename F=double>
using minus_io = minus_io_shaping<P, F>;
// was:
// using minus_io = minus_io_shaping<3, 3, 0, 2, 312, 14, 56, P, double>;  // we now set numbers conditional on P
// TODO: remove this comment and below
//template<problem P>
//using minus6 = minus_core<312, 6, 45, P, double>;
//template<problem P>
//using minusPhoenix10a = minus_core<312, 6, 45, P, double>;
#endif  // minus_h_
