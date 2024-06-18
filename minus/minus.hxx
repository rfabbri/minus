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
#include "minus.h"
#include "internal-util.hxx"

// not really necessary:
#define EIGEN_UNROLLING_LIMIT 1000 
// _really_ necessary:
#define EIGEN_STRONG_INLINE __attribute__((always_inline)) inline
//#include "Eigen-latest/Core"
//#include "Eigen/Core" // enough for 14a
#include "Eigen/LU" // Noob

#define unlikely(expr) __builtin_expect(!!(expr),0)
#define likely(expr)   __builtin_expect(!!(expr),1)

namespace MiNuS {

using namespace Eigen; // only used for linear solve

// #include "chicago14a-lsolve.hxx" TODO XXX make lsolve a template by problem // name
#include "linecircle2a-lsolve.hxx"

template <problem P, typename F>
__attribute__((always_inline)) inline void
memoize_Hxt(C<F> __restrict *block/*, C<F> * __restrict memo*/ /* constants */)
{
  C<F> *const y = reinterpret_cast<C<F> *> (__builtin_assume_aligned(block,64));
//  C<F> *const yc = reinterpret_cast<C<F> *> (__builtin_assume_aligned(memo,64));

  y[11]=y[13]=y[25]=y[27]=y[39]=y[41]=y[53]=y[55]=y[67]=
        y[68]=y[81]=y[82]=y[95]=y[96]=y[109]=y[110]=y[124]=
        y[125]=y[138]=y[139]=y[152]=y[153]=y[166]=y[167]=
        y[180]=y[181]=y[194]=y[195]=0;
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

template <problem P, typename F>
__attribute__((always_inline)) inline void
memoize_HxH(C<F> __restrict *block/*, C<F> * __restrict memo*/ /* constants */)
{
  C<F> *const y = reinterpret_cast<C<F> *> (__builtin_assume_aligned(block,64));
//  C<F> *const yc = reinterpret_cast<C<F> *> (__builtin_assume_aligned(memo,64));

  y[11]=y[13]=y[25]=y[27]=y[39]=y[41]=y[53]=y[55]=y[67]=
        y[68]=y[81]=y[82]=y[95]=y[96]=y[109]=y[110]=y[124]=
        y[125]=y[138]=y[139]=y[152]=y[153]=y[166]=y[167]=
        y[180]=y[181]=y[194]=y[195]=0;
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

// THE MEAT //////////////////////////////////////////////////////////////////////
// t: tracker settings
// s_sols: start sols      
// params: params of target as specialized homotopy params - P01 in SolveChicago
// compute solutions sol_min...sol_max-1 within NSOLS
// 
template <problem P, typename F> void 
minus_core<P, F>::
track(const track_settings &s, const C<F> s_sols_u[f::nve*f::nsols], const C<F> params_u[2*f::nparams], solution raw_solutions_u[f::nsols], unsigned sol_min, unsigned sol_max)
{
  const C<F> *s_sols = reinterpret_cast<C<F> *> (__builtin_assume_aligned(s_sols_u,64));
  const C<F> *params = reinterpret_cast<C<F> *> (__builtin_assume_aligned(params_u,64));
  solution *raw_solutions = reinterpret_cast<solution *> (__builtin_assume_aligned(raw_solutions_u,64));
  assert(sol_min <= sol_max && sol_max <= f::nsols);
  alignas(64) C<F> Hxt[NVEPLUS1 * f::nve]; 
  alignas(64) F x0t0f[f::nve*2+1];
  alignas(64) F xtf[f::nve*2+1];
  alignas(64) F dxdtf[f::nve*2+1];
  alignas(64) C<F> dxi[f::nve];
  C<F> *const x0t0 = (C<F> *) x0t0f;
  C<F> *const xt = (C<F> *) xtf;
  C<F> *const dxdt = (C<F> *) dxdtf;
  C<F> *const x0 = x0t0;
  F    *const t0 = (F *) (x0t0 + f::nve);
  F    *const t  = (F *) (xt + f::nve);
  C<F> *const x1t1 = xt;
  C<F> *const dx = dxdt;
  C<F> *const dx4 = dx;
  F    *const dt = (F *)(dxdt + f::nve);
  C<F> *const HxH=Hxt;
  Map<Matrix<C<F>, f::nve, NVEPLUS1>,Aligned> AA((C<F> *)Hxt,f::nve,NVEPLUS1);
  static constexpr F the_smallest_number = 1e-13; // XXX BENCHMARK THIS
  typedef minus_array<f::nve,F> v;

  //  alignas(64) C<F> ycHxt[16]; 
  //  alignas(64) C<F> ycHxH[13];
  // memoization_init() : 

  // C<F> previous[13];
  const F &t_step = s.init_dt_;  // initial step
  solution *t_s = raw_solutions + sol_min;  // current target solution
  const C<F>* __restrict s_s = s_sols + sol_min*f::nve;    // current start solution
  for (unsigned sol_n = sol_min; sol_n < sol_max; ++sol_n) { // solution loop
    t_s->status = PROCESSING;
    bool end_zone = false;
    v::copy(s_s, x0);
    *t0 = 0; *dt = t_step;
    char predictor_successes = 0;

    // track H(x,t) for t in [0,1]
    while (likely(t_s->status == PROCESSING && 1. - *t0 > the_smallest_number)) {
      if (unlikely(t_s->num_steps == s.max_num_steps_)) {
        t_s->status = MAX_NUM_STEPS_FAIL; // failed to reach solution in the available step budget
        break;
      }
      
      if (unlikely(!end_zone && 1. - *t0 <= s.end_zone_factor_ + the_smallest_number))
        end_zone = true; // TODO: see if this path coincides with any other path on entry to the end zone
      if (unlikely(end_zone)) {
          if (unlikely(*dt > 1. - *t0)) *dt = 1 - *t0;
      } else if (unlikely(*dt > 1. - s.end_zone_factor_ - *t0)) *dt = 1. - s.end_zone_factor_ - *t0;
      /// PREDICTOR /// in: x0t0,dt out: dx
      /*  top-level code for Runge-Kutta-4
          dx1 := solveHxTimesDXequalsminusHt(x0,t0);
          dx2 := solveHxTimesDXequalsminusHt(x0+(1/2)*dx1*dt,t0+(1/2)*dt);
          dx3 := solveHxTimesDXequalsminusHt(x0+(1/2)*dx2*dt,t0+(1/2)*dt);
          dx4 := solveHxTimesDXequalsminusHt(x0+dx3*dt,t0+dt);
          (1/6)*dt*(dx1+2*dx2+2*dx3+dx4) */
      v::fcopy(x0t0, xt);

      // dx1
      // evaluate_Hxt_constants(xt, params, ycHxt);
      memoize_Hxt<P,F>(Hxt);/*, ycHxt);*/
      evaluate_Hxt(xt, params, Hxt); // Outputs Hxt
      // dx4_eigen = lu.compute(AA).solve(bb);
      lsolve<P,F>(AA, dx4);
      
      // dx2
      const F one_half_dt = *dt*0.5;

      v::multiply_scalar_to_self(dx4, one_half_dt);

      v::add_to_self(xt, dx4);
      v::multiply_scalar_to_self(dx4, 2.);
      *t += one_half_dt;  // t0+.5dt
      evaluate_Hxt(xt, params, Hxt);
      memoize_Hxt<P,F>(Hxt);/*, ycHxt);*/
      lsolve<P,F>(AA, dxi);

      // dx3
      v::multiply_scalar_to_self(dxi, one_half_dt);
      v::copy(x0t0, xt);
      v::add_to_self(xt, dxi);
      v::multiply_scalar_to_self(dxi, 4.);
      v::add_to_self(dx4, dxi);
      evaluate_Hxt(xt, params, Hxt);
      memoize_Hxt<P,F>(Hxt);/*, ycHxt);*/
      lsolve<P,F>(AA, dxi);

      // dx4
      v::multiply_scalar_to_self(dxi, *dt);
      v::fcopy(x0t0, xt);
      v::add_to_self(xt, dxi);
      v::multiply_scalar_to_self(dxi, 2.);
      v::add_to_self(dx4, dxi);
      *t = *t0 + *dt;               // t0+dt
      evaluate_Hxt(xt, params, Hxt);
      memoize_Hxt<P,F>(Hxt);/*, ycHxt);*/
      lsolve<P,F>(AA, dxi);
      v::multiply_scalar_to_self(dxi, *dt);
      v::add_to_self(dx4, dxi);
      v::multiply_scalar_to_self(dx4, 1./6.);

      // "dx1" = .5*dx1*dt, "dx2" = .5*dx2*dt, "dx3" = dx3*dt. Eigen vectorizes this:
      // dx4_eigen = (dx4_eigen* *dt + dx1_eigen*2 + dx2_eigen*4 + dx3_eigen*2)*(1./6.);
      
      // make prediction
      v::fcopy(x0t0, x1t1);
      v::fadd_to_self((double *)x1t1, (double *)dxdt);

      
      /// CORRECTOR ///
      char n_corr_steps = 0;
      bool is_successful;
      //if (t_s->num_steps ==0)
      //evaluate_HxH_constants_all_sols(x1t1, params, ycHxH);
      //evaluate_HxH_constants(x1t1, params, ycHxH);
      /*
      {
        static std::mutex lock;
        const std::lock_guard<std::mutex> guard(lock);
         
        if ( t_s->num_steps > 1)
          for (unsigned i=0; i < 13; ++i) {
            F err = std::norm(ycHxH[i]-previous[i]);
            if (err > 1e-15) {
              std::cerr << "Different " << i << " -------------------------------------------- " << err << std::endl;
              std::cerr << "\tnow: " << ycHxH[i] << " previous: " << previous[i] << std::endl;
            }
          }
        for (unsigned i=0; i < 13; ++i)
          previous[i] = ycHxH[i];
      }
      */

      do {
        ++n_corr_steps;
        evaluate_HxH(x1t1, params, HxH);
        memoize_HxH<P,F>(HxH);//, ycHxH);
        lsolve<P,F>(AA, dx);
        v::add_to_self(x1t1, dx);
        is_successful = v::norm2(dx) < s.epsilon2_ * v::norm2(x1t1); // |dx|^2/|x1|^2 < eps2
      } while (likely(!is_successful && n_corr_steps < s.max_corr_steps_));
      
      if (unlikely(!is_successful)) { // predictor failure
        predictor_successes = 0;
        *dt *= s.dt_decrease_factor_;
        if (unlikely(*dt < s.min_dt_)) t_s->status = MIN_STEP_FAILED; // slight difference to SLP-imp.hpp:612
      } else { // predictor success
        ++predictor_successes;
        // std::swap(x1t1,x0t0);
        // x0 = x0t0; t0 = (F *) (x0t0 + f::nve); xt = x1t1;
        v::fcopy(x1t1, x0t0);
        if (unlikely(predictor_successes >= s.num_successes_before_increase_)) {
          predictor_successes = 0;
          *dt *= s.dt_increase_factor_;
        }
      }
      if (unlikely(v::norm2(x0) > s.infinity_threshold2_))
        t_s->status = INFINITY_FAILED;
      ++t_s->num_steps;
    } // while (t loop)
    memcpy(t_s, x0t0, (f::nve*2+1)*sizeof(F));
    if (t_s->status == PROCESSING) t_s->status = REGULAR;
    ++t_s; s_s += f::nve;
  } // outer solution loop
}

// I/O Base functions ----------------------------------------------------------

// RC: same format as cameras_gt_ and synthcurves dataset
// QT: same format as solution_shape
template <problem P, typename F>
inline void 
minus_io_14a<P, F>::
RC_to_QT_format(const F rc[pp::nviews][4][3], F qt[M::nve])
{
  typedef minus_util<F> u;
  F q0[4], q1[4], q2[4];

  u::rotm2quat((F *) rc[0], q0);
  u::rotm2quat((F *) rc[1], q1);
  u::rotm2quat((F *) rc[2], q2);

  // gt = q1 * conj(q0);
  // gt + 4 = q2 * conj(q0);
  u::dquat(q1, q0, qt);
  u::dquat(q2, q0, qt + 4);

  // gt + 8 = q1*(c0-c1)*q1.conj();
  // gt + 8 = quat_transform(q1,c0-c1);
  // gt + 8 + 3 = quat_transform(q2,c0-c2);
  F dc[3];
  dc[0] = rc[0][3][0] - rc[1][3][0];
  dc[1] = rc[0][3][1] - rc[1][3][1];
  dc[2] = rc[0][3][2] - rc[1][3][2];
  u::quat_transform(q1,dc, qt + 8);

  dc[0] = rc[0][3][0] - rc[2][3][0];
  dc[1] = rc[0][3][1] - rc[2][3][1];
  dc[2] = rc[0][3][2] - rc[2][3][2];
  u::quat_transform(q2,dc, qt + 8 + 3);
}

// Returns all real solutions
// The real_solutions array is fixed in size to NSOLS which is the max
// number of solutions, which perfectly fits in memory. The caller must pass an
// array with that minimum.
template <problem P, typename F>
inline void 
minus_io<P, F>::
all_real_solutions(typename M::solution raw_solutions[M::nsols], F real_solutions[M::nsols][M::nve], 
                   unsigned id_sols[M::nsols], unsigned *nsols_real)
{
  typedef minus_array<M::nve,F> v;
  *nsols_real = 0;
  id_sols[*nsols_real] = 0;
  for (unsigned sol=0; sol < M::nsols; ++sol) {
    if (raw_solutions[sol].status == M::REGULAR && v::get_real(raw_solutions[sol].x, real_solutions[*id_sols[nsols_real]]))
      id_sols[(*nsols_real)++] = sol;
  }
}

template <problem P, typename F>
inline void 
minus_io<P, F>::
all_regular_solutions(typename M::solution raw_solutions[M::nsols], C<F> regular_solutions[M::nsols][M::nve], 
                   unsigned id_sols[M::nsols], unsigned *nsols_regular)
{
  typedef minus_array<M::nve,F> v;
  *nsols_regular = 0;
  id_sols[*nsols_regular] = 0;
  for (unsigned sol=0; sol < M::nsols; ++sol) {
    if (raw_solutions[sol].status == M::REGULAR) {
      minus_array<M::f::nve,F>::copy(raw_solutions[sol].x, regular_solutions[*id_sols[nsols_regular]]);
      id_sols[(*nsols_regular)++] = sol;
    }
  }
}


//
// returns cameras[0:nsols_final][2][4][3]
//
// where the camera matrix P^t = [R|T]^t is cameras[sol_number][view_id][:][:]
// where view_id is 0 or 1 for second and third camera relative to the first,
// resp.
//
// This design is for cache speed. Translation in the camera matrix is stored
// such that its coordinates are memory contiguous.
// 
// The cameras array is fixed in size to NSOLS which is the max
// number of solutions, which perfectly fits in memory. The caller must pass an
// array with that minimum.
template <problem P, typename F>
inline void 
minus_io_14a<P, F>::
all_solutions2cams(solution raw_solutions[M::nsols], F cameras[M::nsols][2][4][3], 
                   unsigned id_sols[M::nsols], unsigned *nsols_final)
{
  typedef minus_array<M::nve,F> v;
  *nsols_final = 0;
  for (unsigned sol=0; sol < M::nsols; ++sol) {
    F real_solution[M::nve];
    if (raw_solutions[sol].status == M::REGULAR && v::get_real(raw_solutions[sol].x, real_solution)) {
      id_sols[(*nsols_final)++] = sol;
      // build cams by using quat2rotm
      solution2cams(real_solution, (F (*)[4][3] ) (cameras + sol));
    }
  }
}

// The camera parameter is cameras[img] which is a [4][3] array,
// where the first 3x3 block is R, and the 4th row is T. img is img 0 or 1,
// for 2nd and 3rd cams relative to 1st, resp.
// 
template <problem P, typename F>
inline bool 
minus_io_14a<P, F>::
probe_solutions(const typename M::solution solutions[M::nsols], solution_shape *probe_cameras,
    unsigned *solution_index)
{
  typedef minus_array<M::nve,F> v; typedef minus_util<F> u;
  static constexpr F eps = 1e-3;
  unsigned &sol=*solution_index;
  F real_solution[M::nve];
  for (sol = 0; sol < M::nsols; ++sol) 
    if (v::get_real(solutions[sol].x, real_solution)) {
      u::normalize_quat(real_solution);
      if (u::rotation_error(real_solution, probe_cameras->q01) < eps)
        return true;
    }
  return false;
}

// #undef NDEBUG

#ifndef NDEBUG
#include "debug-util.h"
#endif

/*
template <problem P, typename F>
inline bool 
minus_io_14a<P, F>::
probe_all_solutions(const typename M::solution solutions[M::nsols], solution_shape *probe_cameras,
    unsigned *solution_index)
{
  complex vsolutions[M::nve][M::nsols];
  
  // populate standard format solutions_v from solutions (path info)
  solutions_struct2vector(solutions, vsolutions);
    
  return probe_all_solutions(vsolutions, probe_cameras, solution_index);
}
*/

// like probe_solutions but tests all M::nsols in case more than one is close to
// the probe. Use this for debugging / investigation
template <problem P, typename F>
inline bool 
minus_io_14a<P, F>::
probe_all_solutions(const typename M::solution solutions[M::nsols], solution_shape *probe_cameras,
    unsigned *solution_index)
{
  typedef minus_array<M::nve,F> v; typedef minus_util<F> u;
  static constexpr F eps = 1e-3;
  F real_solution[M::nve];
  bool found=false;
  F min_rerror;
  for (unsigned sol = 0; sol < M::nsols; ++sol)  {
    if (!v::get_real(solutions[sol].x, real_solution))
      continue;
    u::normalize_quat(real_solution);
    F rerror = u::rotation_error(real_solution, probe_cameras->q01);
    if (rerror < eps) {
      if (found == true) {
#ifndef NDEBUG
        std::cerr << "Found another similar solution at " << sol << std::endl;
        std::cerr << "Error: " << rerror << std::endl;
#endif
        if (rerror < min_rerror) {
          min_rerror = rerror;
          *solution_index = sol;
        }
        
      } else {
#ifndef NDEBUG
        std::cerr << "Found a solution at " << sol << std::endl;
        std::cerr << "Error: " << rerror << std::endl;
#endif
        min_rerror = rerror;
        *solution_index = sol;
      }
      found = true;
    } else {
#ifndef NDEBUG
        std::cerr << "Solution is real but not close (sol,isvalid)" << sol << "," << solutions[sol].status << std::endl;
#endif
    }
  }

  if (!found)
    return false;

  // check the remaining parts of the solutions also match, not just rot
  v::get_real(solutions[*solution_index].x, real_solution);
  u::normalize_quat(real_solution+4);
  F rerror = u::rotation_error(real_solution+4, probe_cameras->q02);
  if (rerror < eps) {
#ifndef NDEBUG
    std::cerr << "probe: Rotation 02 also match\n";
#endif
    found = true;
  } else
    found = false;
  
  if (!found)
    return false;

  solution_shape *s = (solution_shape *) real_solution;

  F scale = std::sqrt(minus_3d<F>::dot(s->t01, s->t01));
  F scale_probe = std::sqrt(minus_3d<F>::dot(probe_cameras->t01, probe_cameras->t01));
  
  { // t01
  s->t01[0] /= scale; s->t01[1] /= scale; s->t01[2] /= scale;
  F dt[3];
  dt[0] = s->t01[0] - probe_cameras->t01[0]/scale_probe;
  dt[1] = s->t01[1] - probe_cameras->t01[1]/scale_probe;
  dt[2] = s->t01[2] - probe_cameras->t01[2]/scale_probe;
  
  if (minus_3d<F>::dot(dt, dt) < eps*eps) {
#ifndef NDEBUG
    std::cerr << "probe: translation 01 also match\n";
#endif
    found = true;
  } else {
    dt[0] = s->t01[0] + probe_cameras->t01[0]/scale_probe;
    dt[1] = s->t01[1] + probe_cameras->t01[1]/scale_probe;
    dt[2] = s->t01[2] + probe_cameras->t01[2]/scale_probe;
    if (minus_3d<F>::dot(dt, dt) < eps*eps) {
#ifndef NDEBUG
      std::cerr << "probe: translation 01 also match\n";
#endif
      found = true;
    } else {
      found = false;
#ifndef NDEBUG
      std::cerr << "probe: translation 01 DO NOT match\n";
#endif
    }
  }
  }
  
  if (!found)
    return false;
  
  { // t02
  s->t02[0] /= scale; s->t02[1] /= scale; s->t02[2] /= scale;
  F dt[3];
  dt[0] = s->t02[0] - probe_cameras->t02[0]/scale_probe;
  dt[1] = s->t02[1] - probe_cameras->t02[1]/scale_probe;
  dt[2] = s->t02[2] - probe_cameras->t02[2]/scale_probe;
  
  if (minus_3d<F>::dot(dt, dt) < eps*eps) {
#ifndef NDEBUG
    std::cerr << "probe: translation 02 also match\n";
#endif
    found = true;
  } else {
    //    std::cerr << "dt fail atttempt 1 " << std::endl;
    //    print(dt,3);
    //    std::cerr << "t02 fail atttempt 1 " << std::endl;
    //    print(s->t02,3);
    //    std::cerr << "probe t02 fail atttempt 1 " << std::endl;
    // print(probe_cameras->t02,3);
    
    dt[0] = s->t02[0] + probe_cameras->t02[0]/scale_probe;
    dt[1] = s->t02[1] + probe_cameras->t02[1]/scale_probe;
    dt[2] = s->t02[2] + probe_cameras->t02[2]/scale_probe;
    if (minus_3d<F>::dot(dt, dt) < eps*eps) {
#ifndef NDEBUG
      std::cerr << "probe: translation 02 also match\n";
#endif
      found = true;
    } else {
      found = false;
#ifndef NDEBUG
      std::cerr << "probe: translation 02 DO NOT match\n";
      std::cerr << "dt" << std::endl;
      print(dt,3);
#endif
    }
  }
  }
  return found;
}

// like probe_all_solutions but both solutions and ground truth probe are in
// quaternion-translation format (solution_shape)
//
// This uses real cameras as input, such as in the output of solve_img.
// 
template <problem P, typename F>
inline bool 
minus_io_14a<P, F>::
probe_all_solutions_quat(const F solutions_cameras[M::nsols][M::nve], solution_shape *probe_cameras,
    unsigned nsols, unsigned *solution_index)
{
#ifndef NDEBUG
  std::cerr << "Test xxxxxxxxxxx" << std::endl;
  std::cerr << "Nsols" <<  nsols << std::endl;
#endif
  typedef minus_util<F> u;
  static constexpr F eps = 1e-3;
  F real_solution[M::nve];
  bool found=false;
  F min_rerror;
  for (unsigned sol = 0; sol < nsols; ++sol)  {
    memcpy(real_solution, solutions_cameras[sol], M::nve*sizeof(F));
    u::normalize_quat(real_solution);
    F rerror = u::rotation_error(real_solution, probe_cameras->q01);
    if (rerror < eps) {
      if (found == true) {
#ifndef NDEBUG
        std::cerr << "Found another similar solution at " << sol << std::endl;
        std::cerr << "Error: " << rerror << std::endl;
#endif
        if (rerror < min_rerror) {
          min_rerror = rerror;
          *solution_index = sol;
        }
        
      } else {
#ifndef NDEBUG
        std::cerr << "Found a solution at " << sol << std::endl;
        std::cerr << "Error: " << rerror << std::endl;
#endif
        min_rerror = rerror;
        *solution_index = sol;
      }
      found = true;
    } 
  }

  if (!found)
    return false;

  // check the remaining parts of the solutions also match, not just rot
  memcpy(real_solution, solutions_cameras[*solution_index], M::nve*sizeof(F));
  u::normalize_quat(real_solution+4);
  F rerror = u::rotation_error(real_solution+4, probe_cameras->q02);
  if (rerror < eps) {
#ifndef NDEBUG
    std::cerr << "probe: Rotation 02 also match\n";
#endif
    found = true;
  } else
    found = false;
  
  if (!found)
    return false;

  solution_shape *s = (solution_shape *) real_solution;

  F scale = std::sqrt(minus_3d<F>::dot(s->t01, s->t01));
  F scale_probe = std::sqrt(minus_3d<F>::dot(probe_cameras->t01, probe_cameras->t01));
  
  { // t01
  s->t01[0] /= scale; s->t01[1] /= scale; s->t01[2] /= scale;
  F dt[3];
  dt[0] = s->t01[0] - probe_cameras->t01[0]/scale_probe;
  dt[1] = s->t01[1] - probe_cameras->t01[1]/scale_probe;
  dt[2] = s->t01[2] - probe_cameras->t01[2]/scale_probe;
  
  if (minus_3d<F>::dot(dt, dt) < eps*eps) {
#ifndef NDEBUG
    std::cerr << "probe: translation 01 also match\n";
#endif
    found = true;
  } else {
    dt[0] = s->t01[0] + probe_cameras->t01[0]/scale_probe;
    dt[1] = s->t01[1] + probe_cameras->t01[1]/scale_probe;
    dt[2] = s->t01[2] + probe_cameras->t01[2]/scale_probe;
    if (minus_3d<F>::dot(dt, dt) < eps*eps) {
#ifndef NDEBUG
      std::cerr << "probe: translation 01 also match\n";
#endif
      found = true;
    } else {
      found = false;
#ifndef NDEBUG
      std::cerr << "probe: translation 01 DO NOT match\n";
#endif
    }
  }
  }
  
  if (!found)
    return false;
  
  { // t02
  s->t02[0] /= scale; s->t02[1] /= scale; s->t02[2] /= scale;
  F dt[3];
  dt[0] = s->t02[0] - probe_cameras->t02[0]/scale_probe;
  dt[1] = s->t02[1] - probe_cameras->t02[1]/scale_probe;
  dt[2] = s->t02[2] - probe_cameras->t02[2]/scale_probe;
  
  if (minus_3d<F>::dot(dt, dt) < eps*eps) {
#ifndef NDEBUG
    std::cerr << "probe: translation 02 also match\n";
#endif
    found = true;
  } else {
    //    std::cerr << "dt fail atttempt 1 " << std::endl;
    //    print(dt,3);
    //    std::cerr << "t02 fail atttempt 1 " << std::endl;
    //    print(s->t02,3);
    //    std::cerr << "probe t02 fail atttempt 1 " << std::endl;
    // print(probe_cameras->t02,3);
    
    dt[0] = s->t02[0] + probe_cameras->t02[0]/scale_probe;
    dt[1] = s->t02[1] + probe_cameras->t02[1]/scale_probe;
    dt[2] = s->t02[2] + probe_cameras->t02[2]/scale_probe;
    if (minus_3d<F>::dot(dt, dt) < eps*eps) {
#ifndef NDEBUG
      std::cerr << "probe: translation 02 also match\n";
#endif
      found = true;
    } else {
      found = false;
#ifndef NDEBUG
      std::cerr << "probe: translation 02 DO NOT match\n";
      // std::cerr << "dt" << std::endl;
      // print(dt,3);
#endif
    }
  }
  }
  return found;
}

template <problem P, typename F>
inline bool
minus_io_14a<P, F>::
probe_solutions(const typename M::solution solutions[M::nsols], F probe_cameras[M::nve],
    unsigned *solution_index)
{
  return probe_solutions(solutions, (solution_shape *) probe_cameras, solution_index);
}

template <problem P, typename F>
inline bool
minus_io_14a<P, F>::
probe_all_solutions(const typename M::solution solutions[M::nsols], F probe_cameras[M::nve],
    unsigned *solution_index)
{
  return probe_all_solutions(solutions, (solution_shape *) probe_cameras, solution_index);
}

template <problem P, typename F>
inline bool
minus_io_14a<P, F>::
probe_all_solutions_quat(const F solutions_cameras[M::nsols][M::nve], F probe_cameras[M::nve],
    unsigned nsols, unsigned *solution_index)
{
  return probe_all_solutions_quat(solutions_cameras, (solution_shape *) probe_cameras, nsols, solution_index);
}


// For speed, assumes input point implicitly has 3rd homog coordinate is 1
// 
template <typename F>
inline void 
minus_io_common<F>::
invert_intrinsics(const F K[/*3 or 2 ignoring last line*/][ncoords2d_h], const double pix_coords[][ncoords2d], double normalized_coords[][ncoords2d], unsigned npts)
{
  for (unsigned p=0; p < npts; ++p) {
    const F *px = pix_coords[p];
    F *nrm = normalized_coords[p];
    nrm[1] = (px[1]-K[1][2])/K[1][1];
    nrm[0] = (px[0] - K[0][1]*nrm[1] - K[0][2])/K[0][0];
  }
}

// For speed, assumes input point implicitly has 3rd homog coordinate is 1
// 
template <typename F>
inline void 
minus_io_common<F>::
invert_intrinsics_tgt(const F K[/*3 or 2 ignoring last line*/][ncoords2d_h], const double pix_tgt_coords[][ncoords2d], double normalized_tgt_coords[][ncoords2d], unsigned npts)
{
  for (unsigned p=0; p < npts; ++p) {
    const F *tp = pix_tgt_coords[p];
    F *t = normalized_tgt_coords[p];
    t[1] = tp[1]/K[1][1];
    t[0] = (tp[0] - K[0][1]*t[1])/K[0][0];
  }
}

// Not sure if really necessary.
// Seemed to be important for numerics / error scales at some point.
// Normalizes line normals to unit
template <typename F>
inline void 
minus_io_common<F>::
normalize_lines(F lines[][ncoords2d_h], unsigned nlines)
{
  for (unsigned l=0; l < nlines; ++l)
    normalize_line(lines[l]);
}

} // namespace minus

#include "chicago14a.hxx"      // specific implementation of chicago 14a formulation
#include "linecircle2a.hxx"      // specific implementation of chicago 14a formulation
//#include "cleveland14a.hxx"      // specific implementation of cleveland 14a formulation now in PLMP
// #include <minus/phoenix10a.hxx>      // specific implementation of chicago 14a formulation
// #include "chicago6a.hxx"

#endif // minus_hxx_
