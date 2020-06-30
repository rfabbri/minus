#ifndef chicago14a_internals_h
#define chicago14a_internals_h

namespace MiNuS {
  
template <typename F>
struct minus_io_shaping<chicago14a, F> : public minus_io_14a_base<F> {
  typedef minus_core<chicago14a, F> M;
  typedef struct M::solution solution;
  typedef problem_parameters<chicago14a> pp;
  
  // Input ---------------------------------------------------------------------
  static void point_tangents2params(const F p[pp::nviews][pp::npoints][ncoords2d], const F tgt[pp::nviews][pp::npoints][ncoords2d], unsigned id_tgt0, unsigned id_tgt1, C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  static void point_tangents2params_img(const F p[pp::nviews][pp::npoints][ncoords2d], const F tgt[pp::nviews][pp::npoints][ncoords2d], unsigned id_tgt0, unsigned id_tgt1, const F K[/*3 or 2*/][ncoords2d_h], C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  // this function is the same for all problems
  static void get_params_start_target(F plines[/*15 for chicago*/][ncoords2d_h], C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  static void gammify(C<F> * __restrict__ params/*[ chicago: M::nparams]*/);
  static void point_tangents2lines(const F p[pp::nviews][pp::npoints][ncoords2d], const F tgt[pp::nviews][pp::npoints][ncoords2d], unsigned id_tgt0, unsigned id_tgt1, F plines[pp::nvislines][ncoords2d_h]);
  static void lines2params(const F plines[pp::nvislines][ncoords2d_h], C<F> * __restrict__ params/*[static M::nparams]*/);
};

} // namespace minus

#endif //chicago14a_internals_h
