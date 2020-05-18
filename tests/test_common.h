#ifndef test_common_h
#define test_common_h

// fill in internal format for cameras_gt_
// into camera_gt_quaternion
// - convert each rotation to quaternion in the right order
// - make the rotations all relative to the first camera
static void
minus_initialize_gt()
{
  //  R01 = R1 * inv(R0);
  //  T01 = R1 * (C0 - C1);
  //  R12 = R2 * inv(R1);
  //  T12 = R2 * (C1 - C2);

  // rotation-center format (used in synthcurves dataset)
  // relative to the world, to internal quaternion-translation format relative
  // to first camera
  io::RC_to_QT_format(cameras_gt_, cameras_gt_quat_);
}

#endif test_common_h
