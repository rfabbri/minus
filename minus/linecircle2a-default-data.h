#ifndef linecircle_default_data_h_
#define linecircle_default_data_h_
// see .hxx for documentation

#include "minus.h"

namespace MiNuS {

template <typename F>
struct minus_data<linecircle2a,F> {
  typedef minus_core<linecircle2a> M;
  typedef minus_io<linecircle2a> io;
  typedef minus_io_14a<linecircle2a> io14;
  typedef std::complex<F> complex;
  static const complex start_sols_[M::nve*M::nsols];
  alignas (64) static complex params_start_target_[2*M::f::nparams];
  alignas (64) static complex default_params_start_target_gammified_[2*M::f::nparams];
  static const complex *params_;
  static constexpr unsigned n_gt_sols_ = M::f::nsols; // How many ground-truth sols we hardcode. 
  // Use less than M::f::nsols if you want to deal with a large number of
  // ground-truth solutions but just want to test some                           
  static complex gt_sols_[n_gt_sols_];
  static const unsigned gt_sols_id_[n_gt_sols_]; // id of gt-sol among all NVE
};

} // namespace minus

#endif   // linecircle_default_data_h_
