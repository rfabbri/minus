// 
// Generic documentatin: -------------------------------------------------------
//    See Hxt-doc.md
// 
// Specific documentation: ------------------------------------------------------
//    linecircle2a
//      For this implementation, the parameters are just coefficients. 
//        
template <typename F>
inline __attribute__((always_inline)) void 
eval<linecircle2a, F>::
Hxt(const C<F>* __restrict ux /*x and t*/, const C<F> * __restrict uparams, C<F>* __restrict uy /*Hxt*/) 
{
  const C<F> *params = reinterpret_cast<C<F> *> (__builtin_assume_aligned(uparams,64));
  const C<F> *x = reinterpret_cast<C<F> *> (__builtin_assume_aligned(ux,64));
  C<F> *y = reinterpret_cast<C<F> *> (__builtin_assume_aligned(uy,64));

  const C<F> &X0 = x[0]; // x
  const C<F> &X1 = x[1]; // y
  const C<F> &X2 = x[2]; // t

  const C<F> &X3  = params[0];
  const C<F> &X4  = params[1];
  const C<F> &X5  = params[2];
  const C<F> &X6  = params[3];
  const C<F> &X7  = params[4];
  const C<F> &X8  = params[5];
  const C<F> &X9  = params[6];
  const C<F> &X10 = params[7];
  const C<F> &X11 = params[8];
  const C<F> &X12 = params[9];
  const C<F> &X13 = params[10];
  const C<F> &X14 = params[11];

  static constexpr C<F> C0 = 1;
  static constexpr C<F> C1 = -1;
  const C<F> G0 = C1 * X2;
  const C<F> G1 = C0 + G0;
  const C<F> G2 = G1 * X3;
  const C<F> G3 = X2 * X9;
  const C<F> G4 = G2 + G3;
  const C<F> G5 = X0 + X0;
  const C<F> G6 = G4 * G5;
  const C<F> G7 = G1 * X4;
  const C<F> G8 = X2 * X10;
  const C<F> G9 = G7 + G8;
  const C<F> G10 = G6 + G9;
  const C<F> G11 = G1 * X6;
  const C<F> G12 = X2 * X12;
  const C<F> G13 = G11 + G12;
  const C<F> G14 = X1 + X1;
  const C<F> G15 = G4 * G14;
  const C<F> G16 = G1 * X7;
  const C<F> G17 = X2 * X13;
  const C<F> G18 = G16 + G17;
  const C<F> G19 = X0 * X0;
  const C<F> G20 = X1 * X1;
  const C<F> G21 = G19 + G20;
  const C<F> G22 = C1 * X3;
  const C<F> G23 = G22 + X9;
  const C<F> G24 = G21 * G23;
  const C<F> G25 = C1 * X4;
  const C<F> G26 = G25 + X10;
  const C<F> G27 = X0 * G26;
  const C<F> G28 = G24 + G27;
  const C<F> G29 = C1 * X5;
  const C<F> G30 = G29 + X11;
  const C<F> G31 = G28 + G30;
  const C<F> G32 = C1 * X6;
  const C<F> G33 = G32 + X12;
  const C<F> G34 = X0 * G33;
  const C<F> G35 = C1 * X7;
  const C<F> G36 = G35 + X13;
  const C<F> G37 = X1 * G36;
  const C<F> G38 = G34 + G37;
  const C<F> G39 = C1 * X8;
  const C<F> G40 = G39 + X14;
  const C<F> G41 = G38 + G40;
  y[0] = G10; // NVExNVEPLUS1 matrix [Hx|Ht] as a 1D vector
  y[1] = G13;
  y[2] = G15;
  y[3] = G18;
  y[4] = G31;
  y[5] = G41;
}

// Computes constant values to be memoized. These values do not really depend on
// x, only on params. They can be final matrix values or intermediate ones to
// cache.
//
// TODO: if you fill-in this function, if it is not too big, make sure to force
// inlining even further and guarante that the dissassembled shows no function
// call. You have to test with inlining and without
template <typename F>
inline __attribute__((always_inline)) void 
eval<linecircle2a, F>::
Hxt_constants(const C<F> * __restrict ux, const C<F> * __restrict uparams, C<F> * __restrict uy /*Hxt*/) 
{
  /* If there are separated constants, independent of ux, evaluate them here before the
   * time-dependent evaluator. TODO: work w/team to get extract t-independent evaluators */
}
