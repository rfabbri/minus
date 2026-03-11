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
      unsigned *nsols_final,
      unsigned nthreads=4);
 
  // documentation in P.hxx
  static bool solve_img(
      const F K[/*3 or 2 ignoring last line*/][io::ncoords2d_h],
      const F p[pp::nviews][pp::npoints][io::ncoords2d], 
      const F tgt[pp::nviews][pp::npoints][io::ncoords2d], 
      F solutions_cams[M::nsols][pp::nviews-1][4][3],  // first camera is always [I | 0]^t
      unsigned id_sols[M::nsols],
      unsigned *nsols_final,
      unsigned nthreads=4);
};
#undef P
#if 0
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
#endif
#define P linecircle2a
template <typename F>
struct minus <P, F> {
  typedef minus_core<P, F> M;
  typedef minus_io<P, F> io;
  typedef problem_parameters<P> pp;

  static bool solve(
      const C<F> params_target[M::f::nparams], // p1 in end-linecircle2a.m2 
      C<F> solutions_final[M::nsols][M::nve],  //:< the real solutions with other optional filters
      typename M::solution solutions_target[M::nsols], //:< all the solutions to the target
      unsigned id_sols[M::nsols],
      unsigned *nsols_final,
      unsigned nthreads
    );
};
#undef P
#endif  // problem_defs_h_
