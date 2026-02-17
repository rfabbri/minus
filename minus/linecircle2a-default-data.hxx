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
// P0 || trash
//
// Where P0 are the system parameters associated with the start solution
//
//                                  actually just the start M::nparams
//                                  are initialized here, but we use 
//                                  the other M::nparams later 
//                                  to store target system params
//                                  the latter M::nparams are trash
template <typename F>
alignas(64) std::complex<F> minus_data<linecircle2a,F>::
params_start_target_[2*M::f::nparams] = {
  {.88881367728739857, -.458268749825746},
  {-.90958124762797854, .41552611706549814},
  {.55992250117654629, .82854498530629017},
  {-.96824002271365794, -.25002251581698631},
  {-.91443211037109806, -.40473931798413143},
  {.80682944155872807e-1, .99673981686413049}
  // plus nparams we dont statically initialize and fill later
  // ...j
};

#if 0
// Pro -------------------------------------------------------------------------
// 
// Example randomized parameters for a specific given input for testing
// 
// Used for testing and comparing to M2
//
// This is the line 
//  P01 = p0 || transpose p1;
//
//  toExternalString(point P01)
//
// In tutorial/linecircle-end.m2 example 
//
template <typename F>
alignas(64) std::complex<F> minus_data<linecircle2a,F>::
default_params_start_target_gammified_[2*M::f::nparams] = {
  {-.54232002300332649, .84017200182443086},
  {-.50519558694885758, -.86300487769618039},
  {.78557292739942153, -.61876908110950668},
  {-.22606308570353237, -.9741126635467775},
  {.98518486450687703, .17149572223984591},
  {-.15952718229455548, .98719353629830842},
  // plus nparams we dont statically initialize and fill later
  // in this case a b c d e and f for the line-circle eq
  {2, 0},
  {3, 0},
  {4, 0},
  {6, 0},
  {-2, 0},
  {5.1, 0}
};

#endif

template <typename F>
const std::complex<F> * minus_data<linecircle2a,F>::
// Pro: 
// params_= default_params_start_target_gammified_;
params_= params_start_target_;


// Input parameters corresponding to the above gammified homotopy parameters
template <typename F>
C<F> minus_data<linecircle2a,F>::
solutions_gt_[M::nve] = {
  {-.83999999999999997, -.38032880511473233},
  {30000000000000006e-1, -.11409864153441969}
};

// Ground-truth solutions in your data format could be here.
// see cameras_gt_ in chicago.m2

} // namespace minus

#endif   // linecircle2a_default_data_hxx_
