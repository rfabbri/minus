#ifndef linecircle2a_io_h
#define linecircle2a_io_h

namespace MiNuS {


#if 0
// specialize your minus_io class to the problem, in case you want 
// specific implementations for each function

template <typename F>
struct minus_io<linecircle2a, F> : public minus_io<linecircle2a, F> {
  typedef problem_parameters<linecircle2a> pp;
  typedef minus_core<linecircle2a, F> M;
  typedef minus_io_common<linecircle2a,F> io;
  // Input ---------------------------------------------------------------------
  // static void gammify(C<F> * __restrict params/*[ M::nparams]*/);
  // problem-specific input
  // Any function to convert data to parameters would go here
  // see point_tangents2params
  // Output --------------------------------------------------------------------
  // Shadows base class version if you want to specialize
  static bool has_valid_solutions(const typename M::solution solutions[M::nsols]);
};
#endif

} // namespace minus

#endif //linecircle2a_io_h
