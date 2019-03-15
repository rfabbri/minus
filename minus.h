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
// \endverbatim
//
// OPTIMIZATIONS
//  - see trifocal.key in bignotes for basic results
//  - see CMakeLists.txt and README.md

#define NSOLS 312  /* template these */
#define NNN 14 
#define NNNPLUS1 15 
#define NNNPLUS2 16
#define NNN2 196  /* NNN squared */
#define NPARAMS 56 /* Number of parameters in parameter homotopy - eg, coefficients, etc, to represent NNNxNNN sys */
#include <complex>

// typedef std::complex<double> complex;

template <typename F=double>
using C = std::complex<F>;

// Tracker parameters
// Default values obtained in M2 by doing 
// eg DEFAULT#tStep
// or peek DEFAULT
// 
// We use underscore in case we want to make setters/getters with same name,
// or members of Tracker class if more complete C++ desired
template <typename F=double>
struct TrackerSettings {
  TrackerSettings():
    init_dt_(0.05),   // m2 tStep, t_step, raw interface code initDt
    min_dt_(1e-7),        // m2 tStepMin, raw interface code minDt
    end_zone_factor_(0.03),
    epsilon_(.001), // m2 CorrectorTolerance (chicago.m2, track.m2), raw interface code epsilon (interface2.d, NAG.cpp:rawSwetParametersPT)
    epsilon2_(epsilon_ * epsilon_), 
    max_corr_steps_(5),  // m2 maxCorrSteps (track.m2 param of rawSetParametersPT corresp to max_corr_steps in NAG.cpp)
    dt_increase_factor_(2.),  // m2 stepIncreaseFactor
    dt_decrease_factor_(1./dt_increase_factor_),  // m2 stepDecreaseFactor not existent in DEFAULT, using what is in track.m2:77 
    num_successes_before_increase_(20), // m2 numberSuccessesBeforeIncrease
    infinity_threshold_(1e6), // m2 InfinityThreshold
    infinity_threshold2_(infinity_threshold_ * infinity_threshold_)
  { }
  
  F init_dt_;   // m2 tStep, t_step, raw interface code initDt
  F min_dt_;        // m2 tStepMin, raw interface code minDt
  F end_zone_factor_;
  F epsilon_; // m2 CorrectorTolerance (chicago.m2, track.m2), raw interface code epsilon (interface2.d, NAG.cpp:rawSwetParametersPT)
  F epsilon2_; 
  unsigned max_corr_steps_;  // m2 maxCorrSteps (track.m2 param of rawSetParametersPT corresp to max_corr_steps in NAG.cpp)
  F dt_increase_factor_;  // m2 stepIncreaseFactor
  F dt_decrease_factor_;  // m2 stepDecreaseFactor not existent in DEFAULT, using what is in track.m2:77 
  unsigned num_successes_before_increase_; // m2 numberSuccessesBeforeIncrease
  F infinity_threshold_; // m2 InfinityThreshold
  F infinity_threshold2_;
};

enum SolutionStatus {
  UNDETERMINED,
  PROCESSING,
  REGULAR,
  SINGULAR,
  INFINITY_FAILED,
  MIN_STEP_FAILED,
  ORIGIN_FAILED,
  INCREASE_PRECISION,
  DECREASE_PRECISION
};

template <typename F=double>
struct Solution
{
  C<F> x[NNN];    // array of n coordinates
  F t;          // last value of parameter t used
  SolutionStatus status;
  //  unsigned num_steps;  // number of steps taken along the path
  Solution() : status(UNDETERMINED) { }
};



template <typename F=double>
class Minus {
  public:
  static const TrackerSettings<F> DEFAULT;
  
  static unsigned track_all(const TrackerSettings<F> &s, const C<F> s_sols[NNN*NSOLS], const C<F> params[2*NPARAMS], Solution<F> raw_solutions[NSOLS])
  {
    track(s, s_sols, params, raw_solutions, 0, NNN); // TODO: template start and end
  }
  
  static unsigned track(const TrackerSettings<F> &s, const C<F> s_sols[NNN*NSOLS], const C<F> params[2*NPARAMS], Solution<F> raw_solutions[NSOLS], unsigned sol_min, unsigned sol_max);
};

#endif  // minus_h_
