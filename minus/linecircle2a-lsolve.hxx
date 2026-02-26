#include "Eigen/LU"

#define P linecircle2a
namespace MiNuS {
template <typename F>
struct numeric_subroutines<P, F> {
 static void lsolve(
    Map<Matrix<C<F>, minus_core<P,F>::f::nve, minus_core<P,F>::f::nve +1>,Aligned> & __restrict m, 
    C<F> __restrict *ux);
};
// Solves Ax = b
//
// That is:
//
// given input m = [A | b]
// produces output x such that Ax = b
//
template <typename F>
__attribute__((always_inline)) inline void
numeric_subroutines<P,F>::
lsolve(
    Map<Matrix<C<F>, minus_core<P,F>::f::nve, minus_core<P,F>::f::nve +1>,Aligned> & __restrict m, 
    C<F> __restrict *ux)
{
  C<F> * const x = reinterpret_cast<C<F> *> (__builtin_assume_aligned(ux,64));
  //asm("#------ Lsolve begin"); // there is too many vmovsd moving data. It is sub-vectorized, using only xmm no y or zmm
  typedef minus_core<P, F> M;

  // Noob ----------------------------------------------------------------------
  PartialPivLU<Matrix<C<F>, M::f::nve, M::f::nve> > lu;
  Map<Matrix<C<F>, M::f::nve, 1>, Aligned> x_eigen(x);
  x_eigen = lu.compute(m.template block<M::f::nve,M::f::nve>(0,0)).solve(m.col(M::f::nve));

  // Pro -----------------------------------------------------------------------
  // Fine-tune Eigen's LU to your problem size and structure, e.g., limit pivoting
  // check out linecircle-lsolve-pro.hxx --------------------------------------
}

} // namespace minus
#undef P
