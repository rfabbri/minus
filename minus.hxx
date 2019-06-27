#ifndef minus_hxx_
#define minus_hxx_
// 
// \brief MInimal problem NUmerical continuation package
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date Created: Fri Feb  8 17:42:49 EST 2019
// 
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include "Eigen/LU"
#include "minus.h"
#include "chicago14a.hxx"
// #include "chicago6a.hxx"

template <unsigned NNN, typename F>
struct minus_array { // Speed critical -----------------------------------------
  static inline void 
  multiply_scalar_to_self(C<F> *__restrict__ a, C<F> b)
  {
    for (unsigned i = 0; i < NNN; ++i, ++a) *a = *a * b;
  }

  static inline void
  negate_self(C<F> * __restrict__ a)
  {
    for (unsigned i = 0; i < NNN; ++i, ++a) *a = -*a;
  }

  static inline void 
  multiply_self(C<F> * __restrict__ a, const C<F> * __restrict__ b)
  {
    for (unsigned int i=0; i < NNN; ++i,++a,++b) *a *= *b;
  }

  static inline void 
  add_to_self(C<F> * __restrict__ a, const C<F> * __restrict__ b)
  {
    for (unsigned int i=0; i < NNN; ++i,++a,++b) *a += *b;
  }

  static inline void 
  add_scalar_to_self(C<F> * __restrict__ a, C<F> b)
  {
    for (unsigned int i=0; i < NNN; ++i,++a) *a += b;
  }

  static inline void 
  copy(const C<F> * __restrict__ a, C<F> * __restrict__ b)
  {
    memcpy(b, a, NNN*sizeof(C<F>));
  }

  static inline F
  norm2(const C<F> *__restrict__ a)
  {
    F val = 0;
    C<F> const* __restrict__ end = a+NNN;
    while (a != end) val += std::norm(*a++);
    return val;
  }

  static inline void
  cross(const C<F> v1[3], const C<F> v2[3], C<F> r[3])
  {
    r[0] = v1[1] * v2[2] - v1[2] * v2[1];
    r[1] = v1[2] * v2[0] - v1[0] * v2[2];
    r[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }
  static inline C<F>
  dot(const C<F> v1[3], const C<F> v2[3]) { return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
};

// Not performance critical ---------------------------------------------------
template <unsigned NNN, typename F>  // TODO remove NNN
struct minus_util {
  // Random unit array v of dimension n
  // only on real coordinates, with 0 complex ones
  // we are guaranteeing unifom sampling on the sphere,
  // but simpler rand() on each dimension then normalization also works
  static inline void 
  rand_sphere(C<F> v[5 /*can be 7*/], unsigned n) {
    F m=0;
    for (unsigned i=0; i < n; ++i) {
      F r = gauss(rnd);
      v[i] = C<F>{r};
      m += r*r;
    }
    m = std::sqrt(m);
    for (unsigned i=0; i < n; ++i)
      v[i] /= m;
  }
  // random complex
  static void randc(C<F> *z) { *z = C<F>{gauss(rnd), gauss(rnd)}; *z /= std::abs(*z); }
  static std::random_device rd;
  static std::mt19937 rnd;
  static std::normal_distribution<F> gauss;
};

template <unsigned NNN, typename F>
std::random_device minus_util<NNN,F>::rd;

template <unsigned NNN, typename F>
std::mt19937 minus_util<NNN,F>::rnd{rd()};

template <unsigned NNN, typename F>
std::normal_distribution<F> minus_util<NNN,F>::gauss{0.0,1000.0};  

template <unsigned NSOLS, unsigned NNN, unsigned NPARAMS, problem P, typename F>
const 
typename minus_core<NSOLS, NNN, NPARAMS, P, F>::track_settings minus_core<NSOLS, NNN, NPARAMS, P, F>::DEFAULT;

// THE MEAT //////////////////////////////////////////////////////////////////////
// t: tracker settings
// s_sols: start sols      
// params: params of target as specialized homotopy params - P01 in SolveChicago
// compute solutions sol_min...sol_max-1 within NSOLS
// 
template <unsigned NSOLS, unsigned NNN, unsigned NPARAMS, problem P, typename F>   // only one is NNN
void minus_core<NSOLS, NNN, NPARAMS, P, F>::
track(const track_settings &s, const C<F> s_sols[NNN*NSOLS], const C<F> params[2*NPARAMS], solution raw_solutions[NSOLS], unsigned sol_min, unsigned sol_max)
{
  C<F> Hxt[NNNPLUS1 * NNN] __attribute__((aligned(16))); 
  C<F> x0t0xtblock[2*NNNPLUS1] __attribute__((aligned(16)));
  C<F> dxdt[NNNPLUS1] __attribute__((aligned(16)));
  C<F> dxi[NNN] __attribute__((aligned(16)));
  C<F> *x0t0 = x0t0xtblock;  // t = real running in [0,1]
  C<F> *x0 = x0t0;
  F    *t0 = (F *) (x0t0 + NNN);
  C<F> *xt = x0t0xtblock + NNNPLUS1; 
  C<F> *x1t1 = xt;      // reusing xt's space to represent x1t1
  C<F> *const HxH=Hxt;  // HxH is reusing Hxt
  C<F> *const dx = dxdt;
  const C<F> *const RHS = Hxt + NNN2;  // Hx or Ht, same storage //// UNUSED:  C<F> *const LHS = Hxt;
  C<F> *const dx4 = dx;   // reuse dx for dx4
  F    *const dt = (F *)(dxdt + NNN);
  const F &t_step = s.init_dt_;  // initial step
  using namespace Eigen; // only used for linear solve
  Map<Matrix<C<F>, NNN, 1>,Aligned> dxi_eigen(dxi);
  Map<Matrix<C<F>, NNN, 1>,Aligned> dx4_eigen(dx4);
  Map<Matrix<C<F>, NNN, 1>,Aligned> &dx_eigen = dx4_eigen;
  Map<const Matrix<C<F>, NNN, NNN>,Aligned> AA((C<F> *)Hxt,NNN,NNN);  // accessors for the data
  Map<const Matrix<C<F>, NNN, 1>, Aligned > bb(RHS);
  static constexpr F the_smallest_number = 1e-13;
  typedef minus_array<NNN,F> v; typedef minus_array<NNNPLUS1,F> vp;
  PartialPivLU<Matrix<C<F>, NNN, NNN> > lu;

  solution *t_s = raw_solutions + sol_min;  // current target solution
  const C<F>* __restrict__ s_s = s_sols + sol_min*NNN;    // current start solution
  for (unsigned sol_n = sol_min; sol_n < sol_max; ++sol_n) { // solution loop
    t_s->status = PROCESSING;
    bool end_zone = false;
    v::copy(s_s, x0);
    *t0 = 0; *dt = t_step;
    unsigned predictor_successes = 0;

    // track H(x,t) for t in [0,1]
    while (t_s->status == PROCESSING && 1 - *t0 > the_smallest_number) {
      if (!end_zone && 1 - *t0 <= s.end_zone_factor_ + the_smallest_number)
        end_zone = true; // TODO: see if this path coincides with any other path on entry to the end zone
      if (end_zone) {
          if (*dt > 1 - *t0) *dt = 1 - *t0;
      } else if (*dt > 1 - s.end_zone_factor_ - *t0) *dt = 1 - s.end_zone_factor_ - *t0;
      /// PREDICTOR /// in: x0t0,dt out: dx
      /*  top-level code for Runge-Kutta-4
          dx1 := solveHxTimesDXequalsminusHt(x0,t0);
          dx2 := solveHxTimesDXequalsminusHt(x0+(1/2)*dx1*dt,t0+(1/2)*dt);
          dx3 := solveHxTimesDXequalsminusHt(x0+(1/2)*dx2*dt,t0+(1/2)*dt);
          dx4 := solveHxTimesDXequalsminusHt(x0+dx3*dt,t0+dt);
          (1/6)*dt*(dx1+2*dx2+2*dx3+dx4) */
      vp::copy(x0t0, xt);

      // dx1
      evaluate_Hxt(xt, params, Hxt); // Outputs Hxt
      dx4_eigen = lu.compute(AA).solve(bb);
      
      // dx2
      const C<F> one_half_dt = *dt*0.5;
      v::multiply_scalar_to_self(dx4, one_half_dt);
      v::add_to_self(xt, dx4);
      v::multiply_scalar_to_self(dx4, 2.);
      xt[NNN] += one_half_dt;  // t0+.5dt
      evaluate_Hxt(xt, params, Hxt);
      dxi_eigen = lu.compute(AA).solve(bb);

      // dx3
      v::multiply_scalar_to_self(dxi, one_half_dt);
      v::copy(x0t0, xt);
      v::add_to_self(xt, dxi);
      v::multiply_scalar_to_self(dxi, 4);
      v::add_to_self(dx4, dxi);
      evaluate_Hxt(xt, params, Hxt);
      dxi_eigen = lu.compute(AA).solve(bb);

      // dx4
      v::multiply_scalar_to_self(dxi, *dt);
      vp::copy(x0t0, xt);
      v::add_to_self(xt, dxi);
      v::multiply_scalar_to_self(dxi, 2);
      v::add_to_self(dx4, dxi);
      xt[NNN] = *t0 + *dt;               // t0+dt
      evaluate_Hxt(xt, params, Hxt);
      dxi_eigen = lu.compute(AA).solve(bb);
      v::multiply_scalar_to_self(dxi, *dt);
      v::add_to_self(dx4, dxi);
      v::multiply_scalar_to_self(dx4, 1./6.);

      // "dx1" = .5*dx1*dt, "dx2" = .5*dx2*dt, "dx3" = dx3*dt. Eigen vectorizes this:
      // dx4_eigen = (dx4_eigen* *dt + dx1_eigen*2 + dx2_eigen*4 + dx3_eigen*2)*(1./6.);
      
      // make prediction
      vp::copy(x0t0, x1t1);
      vp::add_to_self(x1t1, dxdt);
      
      /// CORRECTOR ///
      unsigned n_corr_steps = 0;
      bool is_successful;
      do {
        ++n_corr_steps;
        evaluate_HxH(x1t1, params, HxH);
        dx_eigen = lu.compute(AA).solve(bb);
        v::add_to_self(x1t1, dx);
        is_successful = v::norm2(dx) < s.epsilon2_ * v::norm2(x1t1);
      } while (!is_successful && n_corr_steps < s.max_corr_steps_);
      
      if (!is_successful) { // predictor failure
        predictor_successes = 0;
        *dt *= s.dt_decrease_factor_;
        if (*dt < s.min_dt_) t_s->status = MIN_STEP_FAILED; // slight difference to SLP-imp.hpp:612
      } else { // predictor success
        ++predictor_successes;
        std::swap(x1t1,x0t0);
        x0 = x0t0; t0 = (F *) (x0t0 + NNN); xt = x1t1;
        if (predictor_successes >= s.num_successes_before_increase_) {
          predictor_successes = 0;
          *dt *= s.dt_increase_factor_;
        }
      }
      if (v::norm2(x0) > s.infinity_threshold2_)
        t_s->status = INFINITY_FAILED;
    } // while (t loop)
    v::copy(x0, t_s->x); // record the solution
    t_s->t = *t0;
    if (t_s->status == PROCESSING) t_s->status = REGULAR;
    ++t_s; s_s += NNN;
  } // outer solution loop
}

#endif // minus_hxx_
