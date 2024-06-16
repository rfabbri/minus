#ifndef linecircle_default_data_hxx_
#define linecircle_default_data_hxx_
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
template <typename F>
alignas(64) const std::complex<F> minus_data<linecircle2a,F>::
start_sols_[M::nve*M::nsols] = {
  // solution 0
  {-.17713365790154765e1, .18062520635371109e1},  // x
  {-.23659492745314794e1, -.19271488053270542e1}, // y

  // solution 1
  {.97369023998695536, -.2199662496162279}, // x
  {.5699492277349707, -.18897940219334727}, // y
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
  {-.54232002300332649, .84017200182443086},
  {-.50519558694885758, -.86300487769618039},
  {.78557292739942153, -.61876908110950668},
  {-.22606308570353237, -.9741126635467775},
  {.98518486450687703, .17149572223984591},
  {-.15952718229455548, .98719353629830842}
  // plus nparams we dont statically initialize and fill later
  // ...j
};

// Example randomized parameters for a specific given input for testing
// 
// Used for testing and comparing to M2
//
// This is the line 
//  P01 = p0 || transpose p1;
//
//  toExternalString(point P01)
//
// In tutorial/linecircle-start.m2 example 
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
  2,
  3,
  4,
  6,
  -2,
  5.1 
};

template <typename F>
const std::complex<F> * minus_data<linecircle2a,F>::
params_= default_params_start_target_gammified_;


// Input parameters corresponding to the above gammified homotopy parameters

// Input point correspondences for testing could be hardcoded here
// see p_ and p_correct in chicago-default-data.hxx

template <typename F>
C<F> minus_data<linecircle2a,F>::
solutions_gt_[M::nve] = {
  {-.83999999999999997, -.38032880511473233},
  {30000000000000006e-1, -.11409864153441969}
};

// Ground-truth solutions in your data format could be here.
// For instance, the camera matrices in matrix format,
// while internally it would use Cayley or Quaternions format
// see cameras_gt_ in chicago.m2


// The tgt_ array is the same size as the p_ array.
// At each solve only two are used, but since usually all three points have
// tangents, we ask them as input anyways.
// This is in pixel-based image coordinates. 
  
r // namespace minus

#endif   // linecircle_default_data_hxx_
