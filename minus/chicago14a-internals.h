#ifndef chicago14a_internals_h
#define chicago14a_internals_h

namespace MiNuS {

template <typename F>
struct minus_io_14a<chicago14a, F> : public minus_io_common<F> {
  typedef minus_core<chicago14a, F> M;
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
  static void RC_to_QT_format(const F rc[3][4][3], F qt[M::nve]);
  static void all_solutions2cams(solution raw_solutions[M::nsols], F cameras[M::nsols][2][4][3], unsigned id_sols[M::nsols], unsigned *nsols_final);
  static void solution2cams(F rs[M::f::nve], F cameras[2][4][3]);
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
};
  
  
template <typename F>
struct minus_io<chicago14a, F> : public minus_io_14a<chicago14a, F> {
  // template specialization defined in problem-internals.h
  typedef problem_parameters<chicago14a> pp;
  typedef minus_core<chicago14a, F> M;
  // shortcuts to the problem parameters
  static constexpr unsigned  nviews = pp::nviews;
  static constexpr unsigned  npoints = pp::npoints;
  // Input ---------------------------------------------------------------------
  static void gammify(C<F> * __restrict__ params/*[ chicago: M::nparams]*/);
  // Output --------------------------------------------------------------------
  static bool has_valid_solutions(const typename M::solution solutions[M::nsols]);
};

} // namespace minus

#endif //chicago14a_internals_h
