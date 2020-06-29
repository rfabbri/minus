#ifndef cleveland14a_internals_h
#define cleveland14a_internals_h

namespace MiNuS {
  
template <typename F>
struct minus_io_shaping<cleveland14a, F>  : public minus_io_base<F> {
  typedef minus_core<cleveland14a, F> M;
  typedef struct M::solution solution;

  // cast to this to interpret real M::solution::x order
  // internal note: this order is eg  in parser.m2 l 68
  struct solution_shape {
    F q01[4];
    F q02[4];
    F t01[3];
    F t02[3];
  };
  
  static constexpr unsigned ncoords2d = 2;  // just a documented name for the number of inhomog coordinates
  static constexpr unsigned ncoords2d_h = 3;// just a name for the usual number of homog coordinates in P^2
  static constexpr unsigned ncoords3d = 3;  // just a documented name for the number of inhomog 3D coordinates
  typedef problem_parameters<cleveland14a> pp;
#if 0
  { // The basic structure, defined at each problem-specific .hxx
  // unsigned NVIEWS, unsigned NPOINTS /* per view*/, unsigned NFREELINES, unsigned NTANGENTS, 
           
  static constexpr unsigned nviews = NVIEWS; 
  static constexpr unsigned npoints = NPOINTS;
  static constexpr unsigned nfreelines = NFREELINES;
  // even though cleveland needs only 2 tangents, api assumes 3 tangents are given,
  // out of which two are selected by indexing. This is the most common use
  // case, where all features naturally have tangents. If strictly 2 tangents
  // are to be passed, you can leave the unused one as zeros throughout the API.
  static constexpr unsigned ntangents = NTANGENTS;
  // number of lines connecting each pair of points plus going through points
  // plus the number of free lines in the first order problem.
  // Note that tangent orientation may help ruling out solutions; this is why
  // we call it tangent, and not general lines at points. There is more
  // information at tangents which can be used as part of the model for curves.
  // The tangent orientation can be constrained by running without orientation
  // mattering at first, and then propagating these to neighboring features
  // along curves
  // for formulations based on all lines -- not all formulations use this
  static constexpr unsigned nvislines = ( (npoints*(npoints-1) >> 1) + ntangents + nfreelines ) * nviews; 
  }
#endif
  
  // nvislines = 15 for cleveland.
  // INPUT ---------------------------------------------------------------------
  static void point_tangents2params(const F p[pp::nviews][pp::npoints][ncoords2d], const F tgt[pp::nviews][pp::npoints][ncoords2d], unsigned id_tgt0, unsigned id_tgt1, C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  static void point_tangents2params_img(const F p[pp::nviews][pp::npoints][ncoords2d], const F tgt[pp::nviews][pp::npoints][ncoords2d], unsigned id_tgt0, unsigned id_tgt1, const F K[/*3 or 2*/][ncoords2d_h], C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  // this function is the same for all problems
  static void get_params_start_target(F plines[/*15 for cleveland*/][ncoords2d_h], C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  static void gammify(C<F> * __restrict__ params/*[ cleveland: M::nparams]*/);
  static void point_tangents2lines(const F p[pp::nviews][pp::npoints][ncoords2d], const F tgt[pp::nviews][pp::npoints][ncoords2d], unsigned id_tgt0, unsigned id_tgt1, F plines[pp::nvislines][ncoords2d_h]);
  static void lines2params(const F plines[pp::nvislines][ncoords2d_h], C<F> * __restrict__ params/*[static M::nparams]*/);
  static void invert_intrinsics(const F K[/*3 or 2 ignoring last line*/][ncoords2d_h], const double pix_coords[][ncoords2d], double normalized_coords[][ncoords2d], unsigned npts);
  static void invert_intrinsics_tgt(const F K[/*3 or 2 ignoring last line*/][ncoords2d_h], const double pix_tgt_coords[][ncoords2d], double normalized_tgt_coords[][ncoords2d], unsigned npts);
  static void normalize_line(F l[ncoords2d_h]) {
    const F nrm = std::hypot(l[0], l[1]);
    l[0] /= nrm; l[1] /= nrm; l[2] /= nrm;
  }
  static void normalize_lines(F lines[][ncoords2d_h], unsigned nlines);
  static void RC_to_QT_format(const F rc[pp::nviews][4][3], F qt[M::nve]);
  static void rotation_error(const F p[4], const F q[4]);

  // OUTPUT --------------------------------------------------------------------
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
  static bool has_valid_solutions(const typename M::solution solutions[M::nsols]);
};

} // namespace minus

#endif //cleveland14a_internals_h
