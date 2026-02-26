#ifndef linecircle2a_default_data_hxx_
#define linecircle2a_default_data_hxx_
//
// Default start solutions and gammified system parameters
// for inclusion in minus app for default
// 
// to be included after minus.h
//
//
// Developer Notes:
// Start solutions hardcoded for efficiency.
// If you want to play with different start sols,
// write another program that accepts start sols in runtime,
// but keep this one lean & mean.
// We include it separately so they don't clutter this app,
// neither minus.h, and can be reused by other progs
// TODO(developer note): make this part of Minus' template as a specialization. 
// But for efficiency I chose to do it outside.
// Perhaps a minus class should be written that wraps the lean minus_core.
// And in _that_ one, we put these default vectors depending on template tag.


#include "minus.h"
#include "linecircle2a-default-data.h"

namespace MiNuS {


// You can reinterpret as 2D matrix in C as start_sols_[M::nsols][M::nve]
//
// /minus/tutorial/linecircle/linecircle2a-startSys
template <typename F>
alignas(64) const std::complex<F> minus_data<linecircle2a,F>::
start_sols_[M::nve*M::nsols] = {
  // solution 0
  {-.17590111751758006, .79337030532105957},    // x
  {.52121743263206199, .6735216912555235e-1},   // y

  // solution 1
  {.11550940600414354e1, .20377392928276108},   // x
  {-.69566546964533349, .86633155025362207}     // y
};

// Non-gammified (non-randomized)
// Start - target system parameters
//
// p0 || trash
//
// Where p0 are the system parameters associated with the start solution
//
//                                  actually just the start M::nparams
//                                  are initialized here, but we use 
//                                  the other M::nparams later 
//                                  to store target system params
//                                  the latter M::nparams are trash
template <typename F>
alignas(64) std::complex<F> minus_data<linecircle2a,F>::
params_start_target_[2*M::f::nparams] = {
  {.88881367728739857, -.458268749825746}, // start parameters corresponding to start_sols_
  {-.90958124762797854, .41552611706549814},
  {.55992250117654629, .82854498530629017},
  {-.96824002271365794, -.25002251581698631},
  {-.91443211037109806, -.40473931798413143},
  {.80682944155872807e-1, .99673981686413049}
  // ---------------------------------------------------------------------------
  // plus nparams we dont statically initialize and fill later
  // ...j
};

// Example randomized homotopy parameters for a specific input for testing and
// profiling.
// 
// A basic test case is to compare directly against Macaulay2 prototype when
// coding the optimized solver. A more complete time test is the -g profiling
// option of MINUS commands.
// 
// This is
//  P01 = p0 || transpose p1;
//  toExternalString(point P01)
// e.g., in tutorial/linecircle-end.m2 example, where:
// 
// p0: parameters describing the start system, and is the same as
// params_start_target_ 1st half, possibly randomized in PRO-level solvers.
// 
// p1: parameters describing a default target system to run when
// profiling. We usually pick a hard/slow system here.
//
// PRO: These are usually randomized/gammified in pro-level solvers.
//
// TODO: remove suffix _gammified
// -----------------------------------------------------------------------------
// Specific documentation for linecircle2a
//
// p0 and p1 are each the parameters of the linecircle2a problem, see linecircle2a.h
// Each hold the six parameters giving the coefficients of the equations 
template <typename F>
alignas(64) std::complex<F> minus_data<linecircle2a,F>::
default_params_start_target_gammified_[2*M::f::nparams] = {
  {.88881367728739857, -.458268749825746},
  {-.90958124762797854, .41552611706549814},
  {.55992250117654629, .82854498530629017},
  {-.96824002271365794, -.25002251581698631},
  {-.91443211037109806, -.40473931798413143},
  {.80682944155872807e-1, .99673981686413049},
  {2, 0},
  {3, 0},
  {4, 0},
  {6, 0},
  {-2, 0},
  {5.1, 0}
};

template <typename F>
const std::complex<F> * minus_data<linecircle2a,F>::
// Pro: 
// params_= default_params_start_target_gammified_;
params_= params_start_target_;

// Target solutions corresponding to default_params_start_target_gammified_
//
// These should be exactly numerically the same solutions we expect from the homotopy,
// e.g. as prototyped in Macaulay2. Note that there can be multiple different
// representations to the same solutions of a problem, e.g. homogeneous coordinates
template <typename F>
alignas(64) std::complex<F> minus_data<linecircle2a,F>::
gt_sols_[n_gt_sols_][M::nve] = {
  { // solution 1
  {-.83999999999999997, -.38032880511473233},  // x
  {30000000000000006e-1, -.11409864153441969}  // y
  }
  // we dont include more sols now, just want it to find the one above
  // TODO: make it real for your problem, to be more realistic in profiling
};

template <typename F>
const unsigned minus_data<linecircle2a,F>::
gt_sols_id_[minus_data<linecircle2a,F>::n_gt_sols_] = { // Id of the sols to test
 0, 1                               // internal solve. May compare just a few 
                                    // when there may be too many to hardcode.
};

// Ground-truth solutions in your data format could be here.
// see cameras_gt_ in chicago.m2

} // namespace minus

#endif   // linecircle2a_default_data_hxx_
