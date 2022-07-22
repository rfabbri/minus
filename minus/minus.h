#ifndef minus_h_
#define minus_h_
// 
// \brief MInimal problem NUmerical continuation package
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date Created: Fri Feb  8 17:42:49 EST 2019
//
// \verbatim
// Modifications
//    Leykin Feb82019: Initial sketch as simplified code from Macaulay2/NAG.*
//    Tim    Feb2019:  Chicago-specific prototype in Macaulay2
//    Fabbri Mar162019: fully templated and optimized code
// \endverbatim
//
// OPTIMIZATIONS
//  - see Fabbri's trifocal.key in bignotes for basic results
//  - see CMakeLists.txt and README.md

#include <complex>

#if defined(_MSC_VER)
#define __attribute__(x) /* blank - should simply ignore thanks to C preprocessor */
#endif
#include "internal-util.h"

namespace MiNuS {
  
template <typename F>
using C = typename std::complex<F>;

// The problem solvers that this solver template currently supports
enum problem {chicago14a, chicago6a, cleveland14a, phoenix10a /*, standard*/};

// The current best formulations for each problem
constexpr problem chicago = problem::chicago14a;
constexpr problem cleveland = problem::cleveland14a;
// You can now use solver<chicago> to default to the best formulation

// Each problem specializes this in their specific .h
template <problem P>
struct formulation_parameters;

// Each problem specializes this in their specific .h
template <problem P>
struct problem_parameters;

// Problem specific definitions that must be available before anything, at compile time
#include "parameters.h"

// Lowlevel API ----------------------------------------------------------------
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
    DECREASE_PRECISION, // 8 unused
    MAX_NUM_STEPS_FAIL, // 9 failed to converge in less than solution::num_steps 
  };
  
  struct solution
  {
    C<F> x[f::nve];    // array of n coordinates
    F t;               // last value of parameter t used
    solution_status status;
    unsigned num_steps;  // number of steps taken along the path
    solution() : status(UNDETERMINED), num_steps(0) { }
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
    epsilon_(0.000001), // m2 CorrectorTolerance (chicago.m2, track.m2), raw interface code epsilon (interface2.d, NAG.cpp::rawSwetParametersPT)
    epsilon2_(epsilon_ * epsilon_), 
    dt_increase_factor_(2.),  // m2 stepIncreaseFactor
    dt_decrease_factor_(1./dt_increase_factor_),  // m2 stepDecreaseFactor not existent in DEFAULT, using what is in track.m2:77 
    infinity_threshold_(1e7), // m2 InfinityThreshold
    infinity_threshold2_(infinity_threshold_ * infinity_threshold_),
    max_corr_steps_(4),  // m2 maxCorrSteps (track.m2 param of rawSetParametersPT corresp to max_corr_steps in NAG.cpp)
    num_successes_before_increase_(20), // m2 numberSuccessesBeforeIncrease
    max_num_steps_(650)
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
  unsigned max_num_steps_; // maximum number of steps per track. 
                       // Each step takes roughly 1 microseconds (tops)
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
// Internal data ---------------------------------------------------------------
// Data every problem has to declare by specializing this template
template <problem P, typename F=double>
struct minus_data {
};

// I/O -------------------------------------------------------------------------
// Generic I/O routines and defs common to all problems
template <typename F=double>
struct minus_io_common {
  // Variables and types -------------------------------------------------------
  static constexpr unsigned ncoords2d = 2;  // a documented name for the number of inhomog coordinates
  static constexpr unsigned ncoords2d_h = 3;// a name for the usual number of homog coordinates in P^2
  static constexpr unsigned ncoords3d = 3;  // a documented name for the number of inhomog 3D coordinates

  // Functions -----------------------------------------------------------------
  static void invert_intrinsics(const F K[/*3 or 2 ignoring last line*/][ncoords2d_h], const double pix_coords[][ncoords2d], double normalized_coords[][ncoords2d], unsigned npts);
  static void invert_intrinsics_tgt(const F K[/*3 or 2 ignoring last line*/][ncoords2d_h], const double pix_tgt_coords[][ncoords2d], double normalized_tgt_coords[][ncoords2d], unsigned npts);
  static void normalize_line(F l[ncoords2d_h]) {
    const F nrm = std::hypot(l[0], l[1]);
    l[0] /= nrm; l[1] /= nrm; l[2] /= nrm;
  }
  static void normalize_lines(F lines[][ncoords2d_h], unsigned nlines);
};

// Basic I/O function common to formulations that use
// 14 variables = 2* (quaternion + translation)
// This is not specialized to a problem in the implementation,
// but contains common implementations to all problems using 14a formulation
template <problem P, typename F=double>
struct minus_io_14a : public minus_io_common<F> {
  typedef minus_core<P, F> M;
  typedef problem_parameters<P> pp;
  typedef minus_io_common<F> io;
  typedef struct M::solution solution;
  // cast to this to interpret real M::solution::x order
  // internal note: this order is eg  in parser.m2 l 68
  struct solution_shape {
    F q01[4];
    F q02[4];
    F t01[3];
    F t02[3];
  };
  // Output --------------------------------------------------------------------
  static void RC_to_QT_format(const F rc[pp::nviews][4][3], F qt[M::nve]);
  static void all_solutions2cams(solution raw_solutions[M::nsols], F cameras[M::nsols][2][4][3], unsigned id_sols[M::nsols], unsigned *nsols_final);
  static void solution2cams(F rs[M::f::nve], F cameras[2][4][3])
  {
    typedef minus_util<F> u;
    // camera 0 (2nd camera relative to 1st)
    u::quat2rotm(rs, (F *) cameras[0]);
    cameras[0][3][0] = rs[8];
    cameras[0][3][1] = rs[9];
    cameras[0][3][2] = rs[10];
    
    // camera 1 (3rd camera relative to 1st)
    u::quat2rotm(rs+4, (F *) cameras[1]);
    cameras[1][3][0] = rs[11];
    cameras[1][3][1] = rs[12];
    cameras[1][3][2] = rs[13];

    // quat12 rs(0:3), quat12 rs(4:7)
    //  T12 = solutions(9:11);
    //  T13 = solutions(12:14);
    //  R12 = quat2rotm(transpose(quat12));
    //  R13 = quat2rotm(transpose(quat13));
  }

  static bool probe_solutions(const typename M::solution solutions[M::nsols], solution_shape *probe_cameras,
      unsigned *solution_index);
  static bool probe_solutions(const typename M::solution solutions[M::nsols], F probe_cameras[M::nve],
      unsigned *solution_index);
  static bool probe_all_solutions(const typename M::solution solutions[M::nsols], solution_shape *probe_cameras,
      unsigned *solution_index);
  static bool probe_all_solutions(const typename M::solution solutions[M::nsols], F probe_cameras[M::nve],
      unsigned *solution_index);
  static bool probe_all_solutions_quat(const F solutions_cameras[M::nsols][M::nve], solution_shape *probe_cameras,
    unsigned nsols, unsigned *solution_index);
  static bool probe_all_solutions_quat(const F solutions_cameras[M::nsols][M::nve], F probe_cameras[M::nve],
    unsigned nsols, unsigned *solution_index);
  // TODO: move this to generic minus_io - useful for all problems
  static void solutions_struct2vector(const typename M::solution solutions[M::nsols], C<F> sols_v[M::nsols][M::nve])
  {
    for (unsigned s=0; s < M::nsols; ++s)
      for (unsigned var=0; var < M::nve; ++var)
        sols_v[s][var] = solutions[s].x[var];
  }
};

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
// We suggest creating another class minus_io_second_order, etc, for
// problems that may not involve points, tangents nor lines, eg., problems
// involving only conics or higher-order features.
// 
// Feel free to ignore anything in this generic template in the specialization
// to your problem.
template <problem P, typename F=double>
struct minus_io : public minus_io_common<F> {
  // template specialization defined in problem-internals.h
  typedef problem_parameters<P> pp;
  typedef minus_core<P, F> M;
  typedef minus_io_common<F> io;
  // shortcuts to the problem parameters
  static constexpr unsigned  nviews = pp::nviews;
  static constexpr unsigned  npoints = pp::npoints;
  // Input ---------------------------------------------------------------------
  static void gammify(C<F> * __restrict params/*[ chicago: M::nparams]*/);
  // Output --------------------------------------------------------------------
  static bool has_valid_solutions(const typename M::solution solutions[M::nsols]);
};

// Highlevel API ---------------------------------------------------------------
template <problem P, typename F=double>
struct minus {
  // all specializations provide a solve() function
  // each with its own I/O parameters
};

#include "problem-defs.h"

} // namespace minus
#endif  // minus_h_
