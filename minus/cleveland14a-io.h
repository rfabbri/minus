#ifndef cleveland14a_internals_h
#define cleveland14a_internals_h

namespace MiNuS {
  
template <typename F>
struct minus_io<cleveland14a, F> : public minus_io_14a<cleveland14a, F> {
  // template specialization defined in problem-internals.h
  typedef problem_parameters<cleveland14a> pp;
  typedef minus_core<cleveland14a, F> M;
  typedef minus_io_common<F> io;
  // shortcuts to the problem parameters
  static constexpr unsigned  nviews = pp::nviews;
  static constexpr unsigned  npoints = pp::npoints;
  // Input ---------------------------------------------------------------------
  static void gammify(C<F> * __restrict__ params/*[ chicago: M::nparams]*/);
  // Chicago-specific input 
  // this function is the same for all problems
  static void get_params_start_target(
      F plines[/*15 for chicago*/][io::ncoords2d_h], 
      C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  static void points_lines2lines(
      const F p[pp::nviews][pp::npoints][io::ncoords2d], 
      const F l[pp::nviews][io::ncoords2d_h], 
      F plines[pp::nvislines][io::ncoords2d_h]);
  static void lines2params(const F plines[pp::nvislines][io::ncoords2d_h], C<F> * __restrict__ params/*[static M::n//params]*/);
  static void points_lines2params(
      const F p[pp::nviews][pp::npoints][io::ncoords2d], 
      const F l[pp::nviews][io::ncoords2d_h], 
      C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  // Output --------------------------------------------------------------------
  static bool has_valid_solutions(const typename M::solution solutions[M::nsols]);
};

} // namespace minus

#endif //cleveland14a_internals_h
