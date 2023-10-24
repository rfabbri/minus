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

namespace MiNuS {

using namespace Eigen; // only used for linear solve

// TODO: parameters restrict
template <problem P, typename F>
__attribute__((always_inline)) void
minus_core<P, F>::
lsolve(
    Map<const Matrix<C<F>, f::nve, f::nve>,Aligned> &matrix, 
    Map<const Matrix<C<F>, f::nve, 1>, Aligned > &b,
    Map<Matrix<C<F>, f::nve, 1>,Aligned> &x) 
{
   typedef Matrix<C<F>, f::nve, f::nve>  MatrixType;
   typedef PermutationMatrix<f::nve, f::nve> PermutationType;
   typedef Transpositions<f::nve, f::nve> TranspositionType;

   MatrixType m(matrix); // matrix holding LU together TODO: in-place
   PermutationType m_p;

   TranspositionType m_rowsTranspositions;
   // XXX modified by Fabbri to suit Chicago problem
   // static __attribute__((always_inline)) void unblocked_lu(
   //     MatrixType &m, 
   typename TranspositionType::StorageIndex* row_transpositions = &m_rowsTranspositions.coeffRef(0);
   static constexpr Index rows = f::nve;
   //Index first_zero_pivot = -1;
   for(Index k = 0; k < 14; ++k) {
     Index rrows = rows-k-1;

     Index row_of_biggest_in_col(k);
     F biggest_in_corner = std::norm(m(k,k));// std::norm(m.coeff(k,k));
     for (unsigned j=rows-1; j != k; --j) {
       F tmp;
       if ((tmp = std::norm(m(j,k))) > biggest_in_corner*1000) {
           biggest_in_corner = tmp;
           row_of_biggest_in_col = j;
           break;
       }
     }

     row_transpositions[k] = typename TranspositionType::StorageIndex(row_of_biggest_in_col);

     //if (biggest_in_corner != Score(0)) {
     if (k != row_of_biggest_in_col) {
       m.row(k).swap(m.row(row_of_biggest_in_col));
     }

     m.col(k).tail(rrows) /= m(k,k);
     // } else if (first_zero_pivot==-1)
       // the pivot is exactly zero, we record the index of the first pivot which is exactly 0,
       // and continue the factorization such we still have A = PLU
     //  first_zero_pivot = k;

     if (k < rows-1)
       m.bottomRightCorner(rrows,rrows).noalias() -= m.col(k).tail(rrows) * m.row(k).tail(rrows);
   }

   m_p = m_rowsTranspositions;

   // Step 1
   x = m_p * b;

   // TODO: use block indexing and std::vector-std::vector multiplication
   x(1)  -= m(1,0)*x(0);
   x(2)  -= m(2,0)*x(0)+m(2,1)*x(1);
   x(3)  -= m(3,0)*x(0)+m(3,1)*x(1)+m(3,2)*x(2);
   x(4)  -= m(4,0)*x(0)+m(4,1)*x(1)+m(4,2)*x(2)+m(4,3)*x(3);
   x(5)  -= m(5,0)*x(0)+m(5,1)*x(1)+m(5,2)*x(2)+m(5,3)*x(3)+m(5,4)*x(4);
   x(6)  -= m(6,0)*x(0)+m(6,1)*x(1)+m(6,2)*x(2)+m(6,3)*x(3)+m(6,4)*x(4)+m(6,5)*x(5);
   x(7)  -= m(7,0)*x(0)+m(7,1)*x(1)+m(7,2)*x(2)+m(7,3)*x(3)+m(7,4)*x(4)+m(7,5)*x(5)+m(7,6)*x(6);
   x(8)  -= m(8,0)*x(0)+m(8,1)*x(1)+m(8,2)*x(2)+m(8,3)*x(3)+m(8,4)*x(4)+m(8,5)*x(5)+m(8,6)*x(6)+m(8,7)*x(7);
   x(9)  -= m(9,0)*x(0)+m(9,1)*x(1)+m(9,2)*x(2)+m(9,3)*x(3)+m(9,4)*x(4)+m(9,5)*x(5)+m(9,6)*x(6)+m(9,7)*x(7)+m(9,8)*x(8);
   x(10) -= m(10,0)*x(0)+m(10,1)*x(1)+m(10,2)*x(2)+m(10,3)*x(3)+m(10,4)*x(4)+m(10,5)*x(5)+m(10,6)*x(6)+m(10,7)*x(7)+m(10,8)*x(8)+m(10,9)*x(9);
   x(11) -= m(11,0)*x(0)+m(11,1)*x(1)+m(11,2)*x(2)+m(11,3)*x(3)+m(11,4)*x(4)+m(11,5)*x(5)+m(11,6)*x(6)+m(11,7)*x(7)+m(11,8)*x(8)+m(11,9)*x(9)+m(11,10)*x(10);
   x(12) -= m(12,0)*x(0)+m(12,1)*x(1)+m(12,2)*x(2)+m(12,3)*x(3)+m(12,4)*x(4)+m(12,5)*x(5)+m(12,6)*x(6)+m(12,7)*x(7)+m(12,8)*x(8)+m(12,9)*x(9)+m(12,10)*x(10)+m(12,11)*x(11);
   x(13) -= (m(13,0)*x(0)+m(13,1)*x(1)+m(13,2)*x(2)+m(13,3)*x(3)+m(13,4)*x(4)+m(13,5)*x(5)+m(13,6)*x(6)+m(13,7)*x(7)+m(13,8)*x(8)+m(13,9)*x(9)+m(13,10)*x(10)+m(13,11)*x(11)+m(13,12)*x(12));

   // Step 2
   //m.template triangularView<UnitLower>().solveInPlace(x);

   x(13) /= m(13,13);
   x(12) -= m(12,13)*x(13); x(12) /= m(12,12);
   x(11) -= (m(11,12)*x(12)+m(11,13)*x(13)); x(11) /= m(11,11);
   x(10) -= (m(10,11)*x(11)+m(10,12)*x(12)+m(10,13)*x(13)); x(10) /= m(10,10);
   x(9)  -= (m(9,10)*x(10)+m(9,11)*x(11)+m(9,12)*x(12)+m(9,13)*x(13)); x(9) /= m(9,9);
   x(8)  -= (m(8,9)*x(9)+m(8,10)*x(10)+m(8,11)*x(11)+m(8,12)*x(12)+m(8,13)*x(13)); x(8) /= m(8,8);
   x(7)  -= (m(7,8)*x(8)+m(7,9)*x(9)+m(7,10)*x(10)+m(7,11)*x(11)+m(7,12)*x(12)+m(7,13)*x(13)); x(7) /= m(7,7);
   x(6)  -= (m(6,7)*x(7)+m(6,8)*x(8)+m(6,9)*x(9)+m(6,10)*x(10)+m(6,11)*x(11)+m(6,12)*x(12)+m(6,13)*x(13)); x(6) /= m(6,6);
   x(5)  -= (m(5,6)*x(6)+m(5,7)*x(7)+m(5,8)*x(8)+m(5,9)*x(9)+m(5,10)*x(10)+m(5,11)*x(11)+m(5,12)*x(12)+m(5,13)*x(13)); x(5) /= m(5,5);
   x(4)  -= (m(4,5)*x(5)+m(4,6)*x(6)+m(4,7)*x(7)+m(4,8)*x(8)+m(4,9)*x(9)+m(4,10)*x(10)+m(4,11)*x(11)+m(4,12)*x(12)+m(4,13)*x(13)); x(4) /= m(4,4);
   x(3)  -= (m(3,4)*x(4)+m(3,5)*x(5)+m(3,6)*x(6)+m(3,7)*x(7)+m(3,8)*x(8)+m(3,9)*x(9)+m(3,10)*x(10)+m(3,11)*x(11)+m(3,12)*x(12)+m(3,13)*x(13)); x(3) /= m(3,3);
   x(2)  -= (m(2,3)*x(3)+m(2,4)*x(4)+m(2,5)*x(5)+m(2,6)*x(6)+m(2,7)*x(7)+m(2,8)*x(8)+m(2,9)*x(9)+m(2,10)*x(10)+m(2,11)*x(11)+m(2,12)*x(12)+m(2,13)*x(13)); x(2) /= m(2,2);
   x(1)  -= (m(1,2)*x(2)+m(1,3)*x(3)+m(1,4)*x(4)+m(1,5)*x(5)+m(1,6)*x(6)+m(1,7)*x(7)+m(1,8)*x(8)+m(1,9)*x(9)+m(1,10)*x(10)+m(1,11)*x(11)+m(1,12)*x(12)+m(1,13)*x(13)); x(1) /= m(1,1);
   x(0)  -= (m(0,1)*x(1)+m(0,2)*x(2)+m(0,3)*x(3)+m(0,4)*x(4)+m(0,5)*x(5)+m(0,6)*x(6)+m(0,7)*x(7)+m(0,8)*x(8)+m(0,9)*x(9)+m(0,10)*x(10)+m(0,11)*x(11)+m(0,12)*x(12)+m(0,13)*x(13)); x(0) /= m(0,0);
}

      

// THE MEAT //////////////////////////////////////////////////////////////////////
// t: tracker settings
// s_sols: start sols      
// params: params of target as specialized homotopy params - P01 in SolveChicago
// compute solutions sol_min...sol_max-1 within NSOLS
// 
template <problem P, typename F> void 
minus_core<P, F>::
track(const track_settings &s, const C<F> s_sols[f::nve*f::nsols], const C<F> params[2*f::nparams], solution raw_solutions[f::nsols], unsigned sol_min, unsigned sol_max)
{
  assert(sol_min <= sol_max && sol_max <= f::nsols);
  alignas(16) C<F> Hxt[NVEPLUS1 * f::nve]; 
  alignas(16) C<F> x0t0xtblock[2*NVEPLUS1];
  alignas(16) C<F> dxdt[NVEPLUS1];
  alignas(16) C<F> dxi[f::nve];
  C<F> *x0t0 = x0t0xtblock;  // t = real running in [0,1]
  C<F> *x0 = x0t0;
  F    *t0 = (F *) (x0t0 + f::nve);
  C<F> *xt = x0t0xtblock + NVEPLUS1; 
  C<F> *x1t1 = xt;      // reusing xt's space to represent x1t1
  C<F> *const HxH=Hxt;  // HxH is reusing Hxt
  C<F> *const dx = dxdt; const C<F> *const RHS = Hxt + NVE2;  // Hx or Ht, same storage //// UNUSED:  C<F> *const LHS = Hxt;
  C<F> *const dx4 = dx;   // reuse dx for dx4
  F    *const dt = (F *)(dxdt + f::nve);
  const F &t_step = s.init_dt_;  // initial step
  Map<Matrix<C<F>, f::nve, 1>,Aligned> dxi_eigen(dxi);
  Map<Matrix<C<F>, f::nve, 1>,Aligned> dx4_eigen(dx4);
  Map<Matrix<C<F>, f::nve, 1>,Aligned> &dx_eigen = dx4_eigen;
  Map<const Matrix<C<F>, f::nve, f::nve>,Aligned> AA((C<F> *)Hxt,f::nve,f::nve);  // accessors for the data
  Map<const Matrix<C<F>, f::nve, 1>, Aligned > bb(RHS);
  static constexpr F the_smallest_number = 1e-13; // XXX BENCHMARK THIS
  typedef minus_array<f::nve,F> v; typedef minus_array<NVEPLUS1,F> vp;
  PartialPivLU<Matrix<C<F>, f::nve, f::nve> > lu;

  solution *t_s = raw_solutions + sol_min;  // current target solution
  const C<F>* __restrict s_s = s_sols + sol_min*f::nve;    // current start solution
  for (unsigned sol_n = sol_min; sol_n < sol_max; ++sol_n) { // solution loop
    t_s->status = PROCESSING;
    bool end_zone = false;
    v::copy(s_s, x0);
    *t0 = 0; *dt = t_step;
    unsigned predictor_successes = 0;

    // track H(x,t) for t in [0,1]
    while (t_s->status == PROCESSING 
        && 1 - *t0 > the_smallest_number) {
      if (t_s->num_steps == s.max_num_steps_) {
        t_s->status = MAX_NUM_STEPS_FAIL; // failed to reach solution in the available step budget
        break;
      }
      
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
      // dx4_eigen = lu.compute(AA).solve(bb);
      lsolve(AA, bb, dx4_eigen);
      
      // dx2
      const C<F> one_half_dt = *dt*0.5;
      v::multiply_scalar_to_self(dx4, one_half_dt);
      v::add_to_self(xt, dx4);
      v::multiply_scalar_to_self(dx4, 2.);
      xt[f::nve] += one_half_dt;  // t0+.5dt
      evaluate_Hxt(xt, params, Hxt);
      lsolve(AA, bb, dxi_eigen);

      // dx3
      v::multiply_scalar_to_self(dxi, one_half_dt);
      v::copy(x0t0, xt);
      v::add_to_self(xt, dxi);
      v::multiply_scalar_to_self(dxi, 4);
      v::add_to_self(dx4, dxi);
      evaluate_Hxt(xt, params, Hxt);
      lsolve(AA, bb, dxi_eigen);

      // dx4
      v::multiply_scalar_to_self(dxi, *dt);
      vp::copy(x0t0, xt);
      v::add_to_self(xt, dxi);
      v::multiply_scalar_to_self(dxi, 2);
      v::add_to_self(dx4, dxi);
      xt[f::nve] = *t0 + *dt;               // t0+dt
      evaluate_Hxt(xt, params, Hxt);
      lsolve(AA, bb, dxi_eigen);
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
        lsolve(AA, bb, dx_eigen); // TODO: always same AA, do not redo LU
        v::add_to_self(x1t1, dx);
        is_successful = v::norm2(dx) < s.epsilon2_ * v::norm2(x1t1); // |dx|^2/|x1|^2 < eps2
      } while (!is_successful && n_corr_steps < s.max_corr_steps_);
      
      if (!is_successful) { // predictor failure
        predictor_successes = 0;
        *dt *= s.dt_decrease_factor_;
        if (*dt < s.min_dt_) t_s->status = MIN_STEP_FAILED; // slight difference to SLP-imp.hpp:612
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
        t_s->status = INFINITY_FAILED;
      ++t_s->num_steps;
    } // while (t loop)
    v::copy(x0, t_s->x); // record the solution
    t_s->t = *t0; // TODO try to include this in the previous memcpy
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
#include "cleveland14a.hxx"      // specific implementation of cleveland 14a formulation now in PLMP
// #include <minus/phoenix10a.hxx>      // specific implementation of chicago 14a formulation
// #include "chicago6a.hxx"

#endif // minus_hxx_
