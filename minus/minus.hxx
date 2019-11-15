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
#include <minus/minus.h>

template <unsigned N, typename F>
struct minus_array { // Speed critical -----------------------------------------
  static inline void 
  multiply_scalar_to_self(C<F> *__restrict__ a, C<F> b)
  {
    for (unsigned i = 0; i < N; ++i, ++a) *a = *a * b;
  }

  static inline void
  negate_self(C<F> * __restrict__ a)
  {
    for (unsigned i = 0; i < N; ++i, ++a) *a = -*a;
  }

  static inline void 
  multiply_self(C<F> * __restrict__ a, const C<F> * __restrict__ b)
  {
    for (unsigned int i=0; i < N; ++i,++a,++b) *a *= *b;
  }

  static inline void 
  add_to_self(C<F> * __restrict__ a, const C<F> * __restrict__ b)
  {
    for (unsigned int i=0; i < N; ++i,++a,++b) *a += *b;
  }

  static inline void 
  add_scalar_to_self(C<F> * __restrict__ a, C<F> b)
  {
    for (unsigned int i=0; i < N; ++i,++a) *a += b;
  }

  static inline void 
  copy(const C<F> * __restrict__ a, C<F> * __restrict__ b)
  {
    memcpy(b, a, N*sizeof(C<F>));
  }

  static inline F
  norm2(const C<F> *__restrict__ a)
  {
    F val = 0;
    C<F> const* __restrict__ end = a+N;
    while (a != end) val += std::norm(*a++);
    return val;
  }
  
  // Get the real part of the solution vector s (e.g., member x of struct
  // solution).
  //  
  // rs: real solution; holds solution R12, t12, R13, T13 row-major
  //
  // \returns true if the solution is nearly real, false otherwise
  //
  // Not speed critical.
  static inline bool
  get_real(C<F> s[N], F rs[N])
  {
    // Hongyi function realSolutions = parseSolutionString(output)
    // solutions = reshape(solutions,[14,length(solutions)/14]);
    static const double eps = 10e-6;
    
    // TODO[improvement]
    // Fancy way to coNrt to real is to check if the complex number is close to
    // horizontal then get absolute value.
    /*
    for (unsigned var = 0; var < N; ++var)  // differs from Hongyi criterion
      if (s->x[var].real() < eps && s->x[var].real() >= eps
          || std::abs(std::tan(std::arg(s->x[var].imag()))) >= eps)
        return false;
    
    F real_solution[N];
    for (unsigned var = 0; var < N; ++var) 
      real_solution[var] = ((s->x[var].real() >= 0) ? 1 : -1) * std::abs(s->x[var]);
    */
    unsigned var = 0;
    for (; var < N; ++var)
      if (std::abs(s[var].imag()) >= eps) return false;
    for (var = N-1; var != (unsigned)-1; --var) 
      rs[var] = s[var].real();

    // quat12 rs(0:3), quat12 rs(4:7)
    //  T12 = solutions(9:11);
    //  T13 = solutions(12:14);
    //  R12 = quat2rotm(transpose(quat12));
    //  R13 = quat2rotm(transpose(quat13));
    return true;
  }
};

// Functions over 3 dimensions
template <typename F>
struct minus_3d {
  static inline void
  cross(const C<F> v1[3], const C<F> v2[3], C<F> r[3])
  {
    r[0] = v1[1] * v2[2] - v1[2] * v2[1];
    r[1] = v1[2] * v2[0] - v1[0] * v2[2];
    r[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }
  
  static inline C<F>
  dot(const C<F> v1[3], const C<F> v2[3]) { return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
  
  static inline void
  cross(const F v1[3], const F v2[3], F r[3])
  {
    r[0] = v1[1] * v2[2] - v1[2] * v2[1];
    r[1] = v1[2] * v2[0] - v1[0] * v2[2];
    r[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }
  
  static inline F
  dot(const F v1[3], const F v2[3]) { return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
  
  // same as cross but assumes points given in inhomogoeneous coordinates
  // as if v[2] = 1
  static inline void
  cross2(const F v1[2], const F v2[2], F r[3])
  {
    r[0] = v1[1] - v2[1];
    r[1] = v2[0] - v1[0];
    r[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }

  // inhomogeneous points with associated tangents to homogeneous line coefficients
  static inline void 
  point_tangent2line(const F p[2], const F tgt[2], F r[3])
  {
    r[0] = -tgt[1]; // normal vector
    r[1] = tgt[0];
    r[2] = p[0]*tgt[1] - p[1]*tgt[0]; // constant term
  }
};

// Not performance critical ---------------------------------------------------
template <typename F>
struct minus_util {
  // Random unit array v of dimension n
  // only on real coordinates, with 0 complex ones
  // we are guaranteeing unifom sampling on the sphere,
  // but simpler rand() on each dimension then normalization also works
  static inline void 
  rand_sphere(C<F> *v/*[chicago14a: 5 minimum, can be 7]*/, unsigned n) {
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
  static void randc(C<F> * __restrict__ z) { *z = C<F>{gauss(rnd), gauss(rnd)}; *z /= std::abs(*z); }
  static std::random_device rd;
  static std::mt19937 rnd;
  static std::normal_distribution<F> gauss;
  
  static inline void normalize_quat(F q[4])
  {
    const F norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    q[0] /= norm; q[1] /= norm; q[2] /= norm; q[3] /= norm;
  }
  
  // The quaternion order used in the following functions.
  // this is for historical reasons the one used in chicago problem.
  // 
  // Eigen and VXL use {x,y,z,w}, but these function use {w, x, y z} 
  // Ceres also uses {w, x, y, z}, where w is the scalar part.
  // 
  // Used to cast a q[4] vector to interpret entries and use algorithms not
  // mattering the memory order
  struct quaternion_shape { F w; F x; F y; F z; };
  
  // arbitrary quaternion to rotation matrix.
  // will normalize the quaternion in-place.
  // based on VXL/VNL
  static inline void quat2rotm(F qq[4], F r[9])
  {
    normalize_quat(q);
    solution_shape *q = (solution_shape *) qq;
    const F 
      x2 = q->x * q->x,  xy = q->x * q->y,  wx = q->w * q->x,
      y2 = q->y * q->y,  yz = q->y * q->z,  wy = q->w * q->y,
      z2 = q->z * q->z,  zx = q->z * q->x,  wz = q->w * q->z,
      w2 = q->w * q->w;
      
    *r++ = w2 + x2 - y2 - z2;    //  rot(0,0) = r[0]
    *r++ = F(2) * (xy - wz);     //  rot(0,1) = r[1] 
    *r++ = F(2) * (zx + wy);     //  rot(0,2) = r[2] 
    *r++ = F(2) * (xy + wz);     //  rot(1,0) = r[3] 
    *r++ = w2 - x2 + y2 - z2;    //  rot(1,1) = r[4]
    *r++ = F(2) * (yz - wx);     //  rot(1,2) = r[5] 
    *r++ = F(2) * (zx - wy);     //  rot(2,0) = r[6] 
    *r++ = F(2) * (yz + wx);     //  rot(2,1) = r[7] 
    *r   = w2 - x2 - y2 + z2;    //  rot(2,2) = r[8]
  }
  
  static inline void rotm2quat(F q[4], F r[9])
  {
    // use a struct to reinterpret q

    // This algorithm comes from  "Quaternion Calculus and Fast Animation",
    // Ken Shoemake, 1987 SIGGRAPH course notes
    F t = r[] + r[] + r[]; // trace
    if (t > F(0))
    {
      t = std::sqrt(t + F(1.0));
      q.w() = F(0.5)*t;
      t = F(0.5)/t;
      q.x() = (mat.coeff(2,1) - mat.coeff(1,2)) * t;
      q.y() = (mat.coeff(0,2) - mat.coeff(2,0)) * t;
      q.z() = (mat.coeff(1,0) - mat.coeff(0,1)) * t;
    }
    else
    {
      Index i = 0;
      if (mat.coeff(1,1) > mat.coeff(0,0))
        i = 1;
      if (mat.coeff(2,2) > mat.coeff(i,i))
        i = 2;
      Index j = (i+1)%3;
      Index k = (j+1)%3;

      t = std::sqrt(mat.coeff(i,i)-mat.coeff(j,j)-mat.coeff(k,k) + F(1.0));
      q.coeffs().coeffRef(i) = F(0.5) * t;
      t = F(0.5)/t;
      q.w() = (mat.coeff(k,j)-mat.coeff(j,k))*t;
      q.coeffs().coeffRef(j) = (mat.coeff(j,i)+mat.coeff(i,j))*t;
      q.coeffs().coeffRef(k) = (mat.coeff(k,i)+mat.coeff(i,k))*t;
    }
  }
};

template <typename F>
std::random_device minus_util<F>::rd;

template <typename F>
std::mt19937 minus_util<F>::rnd{rd()};

template <typename F>
std::normal_distribution<F> minus_util<F>::gauss{0.0,1000.0};  

template <problem P, typename F> const typename 
minus_core<P, F>::track_settings minus_core<P, F>::DEFAULT;

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
  PartialPivLU<Matrix<C<F>, f::nve, f::nve> > lu;

  solution *t_s = raw_solutions + sol_min;  // current target solution
  const C<F>* __restrict__ s_s = s_sols + sol_min*f::nve;    // current start solution
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
        is_successful = v::norm2(dx) < s.epsilon2_ * v::norm2(x1t1);
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
    } // while (t loop)
    v::copy(x0, t_s->x); // record the solution
    t_s->t = *t0; // TODO try to include this in the previous memcpy
    if (t_s->status == PROCESSING) t_s->status = REGULAR;
    ++t_s; s_s += f::nve;
  } // outer solution loop
}

#include <minus/chicago14a.hxx>      // specific implementation of chicago 14a formulation
#include <minus/phoenix10a.hxx>      // specific implementation of chicago 14a formulation
// XXX #include "phoenix10a.hxx"      // specific implementation of phoenix 10a formulation
// #include "chicago6a.hxx"

#endif // minus_hxx_
