#ifndef common_14a_io_h_
#define common_14a_io_h_

// Basic I/O function common to formulations that use
// 14 variables = 2* (quaternion + translation)
// This is not specialized to a problem in the implementation,
// but contains common implementations to all problems using 14a formulation
template <problem P, typename F=double>
struct minus_io_14a : public minus_io_common<F> {
  typedef minus_core<P, F> M;
  typedef problem_parameters<P> pp;
  typedef minus_io_common<F> io;
  typedef struct M::solution solution;
  // cast to this to interpret real M::solution::x order
  // internal note: this order is eg  in parser.m2 l 68
  struct solution_shape {
    F q01[4];
    F q02[4];
    F t01[3];
    F t02[3];
  };
  // Output --------------------------------------------------------------------
  static void RC_to_QT_format(const F rc[pp::nviews][4][3], F qt[M::nve]);
  static void all_solutions2cams(solution raw_solutions[M::nsols], F cameras[M::nsols][2][4][3], unsigned id_sols[M::nsols], unsigned *nsols_final);
  static void solution2cams(F rs[M::f::nve], F cameras[2][4][3])
  {
    typedef minus_util<F> u;
    // camera 0 (2nd camera relative to 1st)
    u::quat2rotm(rs, (F *) cameras[0]);
    const F n = sqrt(rs[8]*rs[8] + rs[9]*rs[9] + rs[10]*rs[10]) + sqrt(rs[11]*rs[11] + rs[12]*rs[12] + rs[13]*rs[13]);
    cameras[0][3][0] = rs[8]/n;
    cameras[0][3][1] = rs[9]/n;
    cameras[0][3][2] = rs[10]/n;
    
    // camera 1 (3rd camera relative to 1st)
    u::quat2rotm(rs+4, (F *) cameras[1]);
    cameras[1][3][0] = rs[11]/n;
    cameras[1][3][1] = rs[12]/n;
    cameras[1][3][2] = rs[13]/n;

    // quat12 rs(0:3), quat12 rs(4:7)
    //  T12 = solutions(9:11);
    //  T13 = solutions(12:14);
    //  R12 = quat2rotm(transpose(quat12));
    //  R13 = quat2rotm(transpose(quat13));
  }

  static bool probe_solutions(const typename M::solution solutions[M::nsols], solution_shape *probe_cameras,
      unsigned *solution_index);
  static bool probe_solutions(const typename M::solution solutions[M::nsols], F probe_cameras[M::nve],
      unsigned *solution_index);
  static bool probe_all_solutions(const typename M::solution solutions[M::nsols], solution_shape *probe_cameras,
      unsigned *solution_index);
  static bool probe_all_solutions(const typename M::solution solutions[M::nsols], F probe_cameras[M::nve],
      unsigned *solution_index);
  static bool probe_all_solutions_quat(const F solutions_cameras[M::nsols][M::nve], solution_shape *probe_cameras,
    unsigned nsols, unsigned *solution_index);
  static bool probe_all_solutions_quat(const F solutions_cameras[M::nsols][M::nve], F probe_cameras[M::nve],
    unsigned nsols, unsigned *solution_index);
  // TODO: move this to generic minus_io - useful for all problems
  static void solutions_struct2vector(const typename M::solution solutions[M::nsols], C<F> sols_v[M::nsols][M::nve])
  {
    for (unsigned s=0; s < M::nsols; ++s)
      for (unsigned var=0; var < M::nve; ++var)
        sols_v[s][var] = solutions[s].x[var];
  }
};

#endif // common_14a_io_h_ 
