#ifndef comps_h_
#define comps_h_
//:
// \file
// \brief Vision Numerical Algebraic Geometry package for multiview problems
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date Fri Feb  8 17:42:49 EST 2019
//
// \verbatim
// Modifications
//    Leykin Feb82019: Initial sketch as simplified code from Macaulay e/NAG.*
//    Tim    Feb2019:  Chicago-specific prototype in Macaulay2
// \endverbatim
//
//
// OPTIMIZATIONS
//  - see trifocal.key in bignotes for basic results


#include <complex>
#include <cstring>

typedef std::complex<double> complex;

#define NSOLS 312  /* template these */
#define NNN 14 
#define NNNPLUS1 15 
#define NNNPLUS2 16
#define NNN2 196  /* NNN squared */
#define NPARAMS 56 /* Number of parameters in parameter homotopy - eg, coefficients, etc, to represent NNNxNNN sys */


// Tracker parameters
// Default values obtained in M2 by doing 
// eg DEFAULT#tStep
// or peek DEFAULT
// 
// We use underscore in case we want to make setters/getters with same name,
// or members of Tracker class if more complete C++ desired
struct TrackerSettings {
  double init_dt_ = 0.05;   // m2 tStep, t_step, raw interface code initDt
  double min_dt_ = 1e-7;        // m2 tStepMin, raw interface code minDt
  double end_zone_factor_ = 0.05;
  double epsilon_ = .000001; // m2 CorrectorTolerance (chicago.m2, track.m2), raw interface code epsilon (interface2.d, NAG.cpp:rawSwetParametersPT)
  double epsilon2_ = epsilon_ * epsilon_; 
  unsigned max_corr_steps_ = 3;  // m2 maxCorrSteps (track.m2 param of rawSetParametersPT corresp to max_corr_steps in NAG.cpp)
  double dt_increase_factor_ = 2.;  // m2 stepIncreaseFactor
  double dt_decrease_factor_ = 1./dt_increase_factor_;  // m2 stepDecreaseFactor not existent in DEFAULT, using what is in track.m2:77 
  unsigned num_successes_before_increase_ = 2; // m2 numberSuccessesBeforeIncrease
  double infinity_threshold_ = 1e7; // m2 InfinityThreshold
  double infinity_threshold2_ = infinity_threshold_ * infinity_threshold_;
};
// Current settings from Tim: Fri Feb 22 12:00:06 -03 2019 Git 0ec3340
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

static const TrackerSettings comps_DEFAULT;

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

struct Solution
{
  complex x[NNN];        // array of n coordinates
  double t;          // last value of parameter t used
  complex start_x[NNN];  // start of the path that produced x
  SolutionStatus status;
  unsigned num_steps;  // number of steps taken along the path
  Solution() : status(UNDETERMINED) { }
  void make(const complex* s_s) { memcpy(start_x, s_s, NNN*sizeof(complex)); }
};


unsigned ptrack(const TrackerSettings *t, const complex s_sols[NNN*NSOLS], const complex params[NPARAMS], Solution raw_solutions[NSOLS]);


































// --- TESTING ONLY: REMOVE ---------------------------------------------------

// in .h just for testing purposes
bool linear(
    const complex* A,  // size-by-size matrix of complex #s
    const complex* b,  // bsize-by-size RHS of Ax=b
    complex* x   // solution
    );

bool linear_eigen(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    );

bool linear_eigen2(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    );

bool linear_eigen3(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    );
bool linear_eigen3(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    );


#endif  // comps_h_
