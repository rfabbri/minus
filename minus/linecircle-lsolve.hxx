// Specific to Chicago
// Solves Ax = b
//
// That is:
//
// given input m = [A | b]
// produces output x such that Ax = b
//
template <problem P, typename F>
__attribute__((always_inline)) inline void
lsolve(
    Map<Matrix<C<F>, minus_core<P,F>::f::nve, minus_core<P,F>::f::nve +1>,Aligned> & __restrict m, 
    C<F> __restrict *ux)
{
  C<F> * const x= reinterpret_cast<C<F> *> (__builtin_assume_aligned(ux,64));
  //asm("#------ Lsolve begin"); // there is too many vmovsd moving data. It is sub-vectorized, using only xmm no y or zmm
  typedef minus_core<P, F> M;
  static constexpr unsigned char rows = M::f::nve;

}
