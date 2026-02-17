#ifndef linecircle2a_io_h
#define linecircle2a_io_h

namespace MiNuS {

template <typename F>
struct minus_io<linecircle2a, F> : public minus_io_common<F> {
  // template specialization defined in problem-internals.h
  typedef problem_parameters<linecircle2a> pp;
  typedef minus_core<linecircle2a, F> M;
  typedef minus_io_common<F> io;
  // Input ---------------------------------------------------------------------
  // static void gammify(C<F> * __restrict params/*[ M::nparams]*/);
  // problem-specific input
  // Any function to convert data to parameters would go here
  // see point_tangents2params
  // Output --------------------------------------------------------------------
};

} // namespace minus

#endif //linecircle2a_io_h
