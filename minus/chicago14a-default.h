#ifndef chicago14a_default_h_
#define chicago14a_default_h_
// see .hxx for documentation

#include "minus.h"

namespace MiNuS {

#define Float double // XXX TODO: make these as part of template, otherwise dup sym when using two pbms
typedef minus_core<chicago> M;
typedef minus_io<chicago> io;
typedef std::complex<Float> complex;
  
template <typename F>
struct minus_data<chicago14a,F> {
  static const complex start_sols_[M::nve*M::nsols];
  static complex params_start_target_[2*M::f::nparams];
  static const complex default_params_start_target_gammified_[2*M::f::nparams];
  static const complex *params_;
  static Float p_[io::pp::nviews][io::pp::npoints][io::ncoords2d];
  static const Float p_correct_[io::pp::nviews][io::pp::npoints][io::ncoords2d];
  static Float cameras_gt_[io::pp::nviews][4][3];;
  static Float cameras_gt_quat_[M::nve];
  static Float tgt_[io::pp::nviews][io::pp::npoints][io::ncoords2d];
  static const Float tgt_correct_[io::pp::nviews][io::pp::npoints][io::ncoords2d];
  static Float K_[io::ncoords2d][io::ncoords2d_h];
};

}

#endif   // chicago14a_default_h_
