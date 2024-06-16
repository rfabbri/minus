#ifndef chicago14a_default_data_h_
#define chicago14a_default_data_h_
// see .hxx for documentation

#include "minus.h"

namespace MiNuS {

template <typename F>
struct minus_data<linecircle,F> {
  typedef minus_core<linecircle> M;
  typedef minus_io<linecircle> io;
  typedef minus_io_14a<linecircle> io14;
  typedef std::complex<F> complex;
  static const complex start_sols_[M::nve*M::nsols];
  alignas (64) static complex params_start_target_[2*M::f::nparams];
  alignas (64) static complex default_params_start_target_gammified_[2*M::f::nparams];
  static const complex *params_;
};

} // namespace minus

#endif   // chicago14a_default_data_h_
