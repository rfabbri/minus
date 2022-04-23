#ifndef expminus_hxx_
#define expminus_hxx_
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
#include <minus/minus.h>
#include <minus/expminus.h>
#include <minus/internal-util.hxx>

// THE MEAT //////////////////////////////////////////////////////////////////////
// t: tracker settings
// s_sols: start sols      
// params: params of target as specialized homotopy params - P01 in SolveChicago
// compute solutions sol_min...sol_max-1 within NSOLS
// 
template <problem P, typename F> void 
expminus_core<P, F>::
track(const track_settings &s, const C<F> s_sols[f::nve*f::nsols], const C<F> params[2*f::nparams], solution raw_solutions[f::nsols], unsigned sol_min, unsigned sol_max)
{
  assert(sol_min <= sol_max && sol_max <= f::nsols);
  C<F> Hxt[NVEPLUS1 * f::nve] __attribute__((aligned(16))); 
  C<F> x0t0xtblock[2*NVEPLUS1] __attribute__((aligned(16)));
  C<F> dxdt[NVEPLUS1] __attribute__((aligned(16)));
  C<F> dxi[f::nve] __attribute__((aligned(16)));
  C<F> *x0t0 = x0t0xtblock;  // t = real running in [0,1]
  C<F> *x0 = x0t0;
  F    *t0 = (F *) (x0t0 + f::nve);
  C<F> *xt = x0t0xtblock + NVEPLUS1; 
  C<F> *x1t1 = xt;      // reusing xt's space to represent x1t1
  C<F> *const HxH=Hxt;  // HxH is reusing Hxt
  C<F> *const dx = dxdt;
  const C<F> *const RHS = Hxt + NVE2;  // Hx or Ht, same storage //// UNUSED:  C<F> *const LHS = Hxt;
  C<F> *const dx4 = dx;   // reuse dx for dx4
  F    *const dt = (F *)(dxdt + f::nve);
  const F &t_step = s.init_dt_;  // initial step
  using namespace Eigen; // only used for linear solve
  Map<Matrix<C<F>, f::nve, 1>,Aligned> dxi_eigen(dxi);
  Map<Matrix<C<F>, f::nve, 1>,Aligned> dx4_eigen(dx4);
  Map<Matrix<C<F>, f::nve, 1>,Aligned> &dx_eigen = dx4_eigen;
  Map<const Matrix<C<F>, f::nve, f::nve>,Aligned> AA((C<F> *)Hxt,f::nve,f::nve);  // accessors for the data
  Map<const Matrix<C<F>, f::nve, 1>, Aligned > bb(RHS);
  static constexpr F the_smallest_number = 1e-13;
  typedef minus_array<f::nve,F> v; typedef minus_array<NVEPLUS1,F> vp;
  typedef minus<P,F> m;
  PartialPivLU<Matrix<C<F>, f::nve, f::nve> > lu;

  solution *t_s = raw_solutions + sol_min;  // current target solution
  const C<F>* __restrict s_s = s_sols + sol_min*f::nve;    // current start solution
  for (unsigned sol_n = sol_min; sol_n < sol_max; ++sol_n) { // solution loop
    t_s->status = minus_core<P,F>::PROCESSING;
    bool end_zone = false;
    v::copy(s_s, x0);
    *t0 = 0; *dt = t_step;
    unsigned predictor_successes = 0, count = 0;  // number of steps

    // track H(x,t) for t in [0,1]
    while (t_s->status == m::PROCESSING && 1 - *t0 > the_smallest_number) {
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
      xt[f::nve] += one_half_dt;  // t0+.5dt
      evaluate_Hxt(xt, params, Hxt);
      dxi_eigen = lu.compute(AA).solve(bb);  // TODO: reuse pivots

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
      xt[f::nve] = *t0 + *dt;               // t0+dt
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
        is_successful = v::norm2(dx) < s.epsilon2_ * v::norm2(x1t1); // |dx|^2/|x1|^2 < eps2
      } while (!is_successful && n_corr_steps < s.max_corr_steps_);
      
      if (!is_successful) { // predictor failure
        predictor_successes = 0;
        *dt *= s.dt_decrease_factor_;
        if (*dt < s.min_dt_) t_s->status = m::MIN_STEP_FAILED; // slight difference to SLP-imp.hpp:612
      } else { // predictor success
        ++predictor_successes;
        std::swap(x1t1,x0t0);
        x0 = x0t0; t0 = (F *) (x0t0 + f::nve); xt = x1t1;
        if (predictor_successes >= s.num_successes_before_increase_) {
          predictor_successes = 0;
          *dt *= s.dt_increase_factor_;
        }
      }
      if (v::norm2(x0) > s.infinity_threshold2_)
        t_s->status = m::INFINITY_FAILED;
      if (++count > 500)
        break;
    } // while (t loop)
    v::copy(x0, t_s->x); // record the solution
    t_s->t = *t0; // TODO try to include this in the previous memcpy
    if (t_s->status == m::PROCESSING) t_s->status = m::REGULAR;
    t_s->num_steps = count;
    ++t_s; s_s += f::nve;
  } // outer solution loop
}

template <problem P, typename F> const typename 
expminus_core<P, F>::track_settings expminus_core<P, F>::DEFAULT;

#endif // expminus_hxx_
