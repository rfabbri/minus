#ifndef expminus_h_
#define expminus_h_
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

// problem specific definitions that must be available before anything, at compile time
#include <minus/parameters.h>

template <problem P, typename F=double>
class expminus_core { // fully static, not to be instantiated - just used for templating
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
  
  struct solution : public minus_core<P,F>::solution
  {
    unsigned num_steps;  // number of steps taken along the path: don't mess with order!
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
struct expminus_core<P, F>::track_settings {
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

template <problem P, typename F>
void expminus_core<P, F>::evaluate_Hxt(const C<F> * __restrict x /*x, t*/, const C<F> * __restrict params, C<F> * __restrict y /*HxH*/)
{
  eval<P,F>::Hxt(x, params, y);
}

template <problem P, typename F>
void expminus_core<P, F>::evaluate_HxH(const C<F> * __restrict x /*x, t*/, const C<F> * __restrict params, C<F> * __restrict y /*HxH*/)
{
  eval<P,F>::HxH(x, params, y);
}

// Shortcuts and aliases -------------------------------------------------------
// type alias used to hide a template parameter 
template<problem P, typename F=double>
using expminus = expminus_core<P, F>;  

#endif  // expminus_h_
