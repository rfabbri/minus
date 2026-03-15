// 
// Generic documentation: -------------------------------------------------------
//    See HxH-doc.md
// 
template <typename F>
inline __attribute__((always_inline)) void 
eval<linecircle2a, F>::
HxH(const C<F>* __restrict ux /*x and t*/, const C<F> * __restrict uparams, C<F>* __restrict uy /*HxH*/) 
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
  
  const C<F> C0 = 1;
  const C<F> C1 = -1;
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
  const C<F> G22 = G4 * G21;
  const C<F> G23 = G9 * X0;
  const C<F> G24 = G22 + G23;
  const C<F> G25 = G1 * X5;
  const C<F> G26 = X2 * X11;
  const C<F> G27 = G25 + G26;
  const C<F> G28 = G24 + G27;
  const C<F> G29 = G13 * X0;
  const C<F> G30 = G18 * X1;
  const C<F> G31 = G29 + G30;
  const C<F> G32 = G1 * X8;
  const C<F> G33 = X2 * X14;
  const C<F> G34 = G32 + G33;
  const C<F> G35 = G31 + G34;
  y[0] = G10; // NVExNVEPLUS1 matrix [Hx|H] as a 1D vector, col-major
  y[1] = G13;
  y[2] = G15;
  y[3] = G18;
  y[4] = -G28;
  y[5] = -G35;
}

template <typename F>
inline __attribute__((always_inline)) void 
eval<linecircle2a, F>::
HxH_constants(const C<F>* __restrict ux /*x and t*/, const C<F> * __restrict uparams, C<F>* __restrict uy /*HxH*/) 
{
}

template <typename F>
inline __attribute__((always_inline)) void 
eval<linecircle2a, F>::
HxH_constants_all_sols(const C<F>* __restrict ux /*x and t*/, const C<F> * __restrict uparams, C<F>* __restrict uy /*HxH*/) 
{
}

template <typename F>
__attribute__((always_inline)) inline void
eval<linecircle2a, F>::
HxH_memoize(C<F> __restrict *block/*, C<F> * __restrict memo*/ /* constants */)
{
//  C<F> *const y = reinterpret_cast<C<F> *> (__builtin_assume_aligned(block,64));
//  C<F> *const yc = reinterpret_cast<C<F> *> (__builtin_assume_aligned(memo,64));

//  y[11]=y[13]=y[25]=y[27]=y[39]=y[41]=y[53]=y[55]=y[67]=
//        y[68]=y[81]=y[82]=y[95]=y[96]=y[109]=y[110]=y[124]=
//        y[125]=y[138]=y[139]=y[152]=y[153]=y[166]=y[167]=
//        y[180]=y[181]=y[194]=y[195]=0;
//  y[26]  = yc[0];
//  y[40]  = yc[1];
//  y[54]  = yc[2];
//  y[69]  = yc[3];
//  y[83]  = yc[4];
//  y[97]  = yc[5];
//  y[111] = yc[6];
//  y[123] = yc[7];
//  y[137] = yc[8];
//  y[151] = yc[9];
//  y[165] = yc[10];
//  y[179] = yc[11];
//  y[193] = yc[12];
//  y[207] = yc[13];
//  y[208] = yc[14];
//  y[209] = yc[15];
}
