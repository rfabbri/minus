#ifndef internal_util_hxx_
#define internal_util_hxx_

#include "internal-util.h"

namespace MiNuS {
  
template <typename F>
std::random_device minus_util<F>::rd;

template <typename F>
std::mt19937 minus_util<F>::rnd{rd()};

template <typename F>
std::normal_distribution<F> minus_util<F>::gauss{0.0,1000.0};  

template <problem P, typename F> const typename 
minus_core<P, F>::track_settings minus_core<P, F>::DEFAULT;

} // namespace minus

#endif // internal_util_hxx_
