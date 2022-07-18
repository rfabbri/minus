#ifndef chicago14a_default_data_h_
#define chicago14a_default_data_h_
// see .hxx for documentation

#include "minus.h"

namespace MiNuS {

template <typename F>
struct minus_data<chicago14a,F> {
  typedef minus_core<chicago14a> M;
  typedef minus_io<chicago14a> io;
  typedef minus_io_14a<chicago14a> io14;
  typedef std::complex<F> complex;
  static const complex start_sols_[M::nve*M::nsols];
  static complex params_start_target_[2*M::f::nparams];
  static complex default_params_start_target_gammified_[2*M::f::nparams];
  static const complex *params_;
  static F p_[io::pp::nviews][io::pp::npoints][io::ncoords2d];
  static const F p_correct_[io::pp::nviews][io::pp::npoints][io::ncoords2d];
  static F cameras_gt_[io::pp::nviews][4][3];
  static F cameras_gt_quat_[M::nve];
  static F tgt_[io::pp::nviews][io::pp::npoints][io::ncoords2d];
  static const F tgt_correct_[io::pp::nviews][io::pp::npoints][io::ncoords2d];
  static F K_[io::ncoords2d][io::ncoords2d_h];
  // fill in internal format for cameras_gt_
  // into camera_gt_quaternion
  // - convert each rotation to quaternion in the right order
  // - make the rotations all relative to the first camera
  static void initialize_gt() {
    //  R01 = R1 * inv(R0);
    //  T01 = R1 * (C0 - C1);
    //  R12 = R2 * inv(R1);
    //  T12 = R2 * (C1 - C2);

    // rotation-center format (used in synthcurves dataset)
    // relative to the world, to internal quaternion-translation format relative
    // to first camera
    io14::RC_to_QT_format(cameras_gt_, cameras_gt_quat_);
  }
};

} // namespace minus

#endif   // chicago14a_default_data_h_
