#ifndef 2v5p9a_internals_h
#define 2v5p9a_internals_h
// TODO: rename to io

namespace MiNuS {

template <typename F>
struct minus_io<2v5p9a, F> : public minus_io_14a<2v5p9a, F> {
  // template specialization defined in problem-internals.h
  typedef problem_parameters<2v5p9a> pp;
  typedef minus_core<2v5p9a, F> M;
  typedef minus_io_common<F> io;
  // shortcuts to the problem parameters
  static constexpr unsigned nviews  = pp::nviews;
  static constexpr unsigned npoints = pp::npoints;
  // Input ---------------------------------------------------------------------
  static void gammify(C<F> * __restrict__ params/*[ chicago: M::nparams]*/);
  // Chicago-specific input 
  static void point2params(const F p[pp::nviews][pp::npoints][io::ncoords2d], C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  static void point2params_img(const F p[pp::nviews][pp::npoints][io::ncoords2d], const F K[/*3 or 2*/][io::ncoords2d_h], C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  // this function is the same for all problems
  static void get_params_start_target(F plines[/*15 for chicago*/][io::ncoords2d_h], C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  // XXX remove next function or write points2internal_representation
  static void point_tangents2lines(const F p[pp::nviews][pp::npoints][io::ncoords2d], const F tgt[pp::nviews][pp::npoints][io::ncoords2d], unsigned id_tgt0, unsigned id_tgt1, F plines[pp::nvislines][io::ncoords2d_h]);
  static void lines2params(const F plines[pp::nvislines][io::ncoords2d_h], C<F> * __restrict__ params/*[static M::n//params]*/);
  // Output --------------------------------------------------------------------
  static bool has_valid_solutions(const typename M::solution solutions[M::nsols]);
};

} // namespace minus

#endif //2v5p9a_internals_h
