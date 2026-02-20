#ifndef internal_util_hxx_
#define internal_util_hxx_

#include "internal-util.h"

// not really necessary:
#define EIGEN_UNROLLING_LIMIT 1000 
// _really_ necessary:
#define EIGEN_STRONG_INLINE __attribute__((always_inline)) inline
//#include "Eigen-latest/Core"
#include "Eigen/Core"

namespace MiNuS {
  
template <typename F>
std::random_device minus_util<F>::rd;

template <typename F>
std::mt19937 minus_util<F>::rnd{rd()};

template <typename F>
std::normal_distribution<F> minus_util<F>::gauss{0.0,1000.0};  

template <problem P, typename F> const typename 
minus_core<P, F>::track_settings minus_core<P, F>::DEFAULT;

const formulation_parameters<chicago14a>::settings formulation_parameters<chicago14a>::DEFAULT;

using namespace Eigen; // only used for linear solve
// construct to enable partial template instantiation for each problem specific
// linear solve, etc.
template <problem P, typename F>
struct numeric_subroutines {
 static void lsolve(
    Map<Matrix<C<F>, minus_core<P,F>::f::nve, minus_core<P,F>::f::nve +1>,Aligned> & __restrict m, 
    C<F> __restrict *ux);
};


} // namespace minus

// branch prediction annotation
#define unlikely(expr) __builtin_expect(!!(expr),0)
#define likely(expr)   __builtin_expect(!!(expr),1)



#endif // internal_util_hxx_
