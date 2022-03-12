#ifndef problem_defs_h_
#define problem_defs_h_
// to be included in minus.h
// These are just APIs / declarations. The implementations are in the
// problem-specific .hxx's

#define P chicago14a
template <typename F>
struct minus <P, F> {
  typedef minus_core<P, F> M;
  typedef minus_io<P, F> io;
  typedef problem_parameters<P> pp;

  // documentation in P.hxx
  static bool solve(
      const F p[pp::nviews][pp::npoints][io::ncoords2d], 
      const F tgt[pp::nviews][pp::npoints][io::ncoords2d], 
      F solutions_cams[M::nsols][pp::nviews-1][4][3],  // first camera is always [I | 0]^t
      unsigned id_sols[M::nsols],
      unsigned *nsols_final);
 
  // documentation in P.hxx
  static bool solve_img(
      const F K[/*3 or 2 ignoring last line*/][io::ncoords2d_h],
      const F p[pp::nviews][pp::npoints][io::ncoords2d], 
      const F tgt[pp::nviews][pp::npoints][io::ncoords2d], 
      F solutions_cams[M::nsols][pp::nviews-1][4][3],  // first camera is always [I | 0]^t
      unsigned id_sols[M::nsols],
      unsigned *nsols_final);
};
#undef P
#define P cleveland14a
template <typename F>
struct minus <P, F> {
  typedef minus_core<P, F> M;
  typedef minus_io<P, F> io;
  typedef problem_parameters<P> pp;

  // documentation in P.hxx
  static bool solve(
      const F p[pp::nviews][pp::npoints][io::ncoords2d], 
      const F tgt[pp::nviews][pp::npoints][io::ncoords2d], 
      F solutions_cams[M::nsols][pp::nviews-1][4][3],  // first camera is always [I | 0]^t
      unsigned id_sols[M::nsols],
      unsigned *nsols_final);
 
  // documentation in P.hxx
  static bool solve_img(
      const F K[/*3 or 2 ignoring last line*/][io::ncoords2d_h],
      const F p[pp::nviews][pp::npoints][io::ncoords2d], 
      const F tgt[pp::nviews][pp::npoints][io::ncoords2d], 
      F solutions_cams[M::nsols][pp::nviews-1][4][3],  // first camera is always [I | 0]^t
      unsigned id_sols[M::nsols],
      unsigned *nsols_final);
};
#undef P
#endif  // problem_defs_h_
