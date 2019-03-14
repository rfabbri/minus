// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date Created: Fri Feb  8 17:42:49 EST 2019
// 
#include "minus.h"
#include <cstdio>
#include <iostream>
#include <iomanip>
#include "Eigen/Core"
#include "Eigen/LU"
//#include "Eigen/QR"
#include "chicago.hxx"

/* Lapack .h */
#if 0
extern "C" {
int sgesv_(int *n,      // number of rows in A
           int *nrhs,   // number of right hand sides
           double *a,   // n by n matrix A, on exit L&U from A=PLU
           int *lda,    // n
           int *ipiv,   // indices defining permutation P
           double *b,   // right-hand-side
           int *ldb,    // n
           int *info);  // error info

};
#endif


// Place any specific-type functions here


#if 0
// Original code solve_via_lapack_without_transposition
//
// \returns  lapack info code. 
//      >0 -> matrix is singular
//      <0 -> illegal value of an argument passed to lapack
bool linear(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  int info;
  static int permutation[NNN]; // unused
  static int bsize = 1;
  static int size = NNN;
  
  double *copyA = (double*)A;
  
  // TODO try to eliminate this memcpy and trash the original b if possible
  // memcpy  b -> x          NNN elements
  std::memcpy(x, b, 2*NNN*sizeof(double));

  double *copyb = (double*)x;  // result is stored in copyb
  
  sgesv_(&size, &bsize, copyA, &size, permutation, copyb, &size, &info);

  return info == 0;
}
#endif

/*
static bool
linear_eigen(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  using namespace Eigen;
  
  Map<Matrix<complex, NNN, 1> > xx(x);
  Map<const Matrix<complex, NNN, NNN> > AA(A,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > bb(b);
  
  xx = AA.colPivHouseholderQr().solve(bb);
  return true;
}
*/
  

static bool 
linear_eigen2(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  using namespace Eigen;
  
  Map<Matrix<complex, NNN, 1> > xx(x);
  Map<const Matrix<complex, NNN, NNN> > AA(A,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > bb(b);
  
  xx = AA.partialPivLu().solve(bb);
  return true; // TODO: better error handling
}

/*
static bool 
linear_eigen3(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  using namespace Eigen;
  
  Map<Matrix<complex, NNN, 1> > xx(x);
  Map<const Matrix<complex, NNN, NNN> > AA(A,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > bb(b);
  
  xx = AA.fullPivLu().solve(bb);
  return true;
}
*/

/*
static bool 
linear_eigen4(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  using namespace Eigen;
  
  Map<Matrix<complex, NNN, 1> > xx(x);
  Map<const Matrix<complex, NNN, NNN> > AA(A,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > bb(b);
  
  xx = AA.householderQr().solve(bb);
  return true;
}
*/

// 20s
// Direct inversion - good for smaller matrices, 5x5 etc
static bool 
linear_eigen5(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  using namespace Eigen;
  
  Map<Matrix<complex, NNN, 1> > xx(x);
  Map<const Matrix<complex, NNN, NNN> > AA(A,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > bb(b);
  
  // xx = AA.partialPivLu().solve(bb);
  // 
  Matrix<complex, NNN, NNN> inv = AA.inverse();
  xx = inv * bb;
  // xx += inv*(bb-AA*xx); // correct
  return true; // TODO: better error handling
}
// from compute_Grabner_basis, 5point.c
#if 0
bool linear_bundler(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
    for (unsigned i = 0; i < 10; i++) {
        /* Make the leading coefficient of row i = 1 */
        double leading = A[20 * i + i];
        matrix_scale(20, 1, A + 20 * i, 1.0 / leading, A + 20 * i);

        /* Subtract from other rows */
        for (j = i+1; j < 10; j++) {
            double leading2 = A[20 * j + i];
            double scaled_row[20];
            matrix_scale(20, 1, A + 20 * i, leading2, scaled_row);
            matrix_diff(20, 1, 20, 1, A + 20 * j, scaled_row, A + 20 * j);
        }
    }

    /* Now, do the back substitution */
    for (i = 9; i >= 0; i--) {
        for (j = 0; j < i; j++) {
            double scale = A[20 * j + i];
            double scaled_row[20];
            matrix_scale(20, 1, A + 20 * i, scale, scaled_row);
            matrix_diff(20, 1, 20, 1, A + 20 * j, scaled_row, A + 20 * j);
        }
    }
    
    /* Copy out results */
    for (i = 0; i < 10; i++) {
        memcpy(Gbasis + i * 10, A + i * 20 + 10, sizeof(double) * 10);
    }
}
#endif

static const double the_smallest_number = 1e-13;
static const double dbgtol = 1e-2;

static void
array_print(std::string name, const complex *a)
{
  std::cerr << name + "[NNN]: ";
  std::cerr << a[0] << " , " << a[1] << " , ... " <<  a[NNN-1] << std::endl;
}

static void
array_print_NNNplus1(std::string name, const complex *a)
{
  std::cerr << name + "[NNN+1]: ";
  std::cerr << a[0] << " , " << a[1] << " , ... " <<  a[NNN] << std::endl;
}

static void
array_print_n(std::string name, const complex *a, unsigned n)
{
  assert(n > 3);
  std::cerr << name + "[" + std::to_string(n) + "]: ";
  std::cerr << a[0] << " , " << a[1] << " , ... " <<  a[n-1] << std::endl;
}

static void
array_print_H(std::string name, const complex *a)
{
  std::cerr << name + "[][]: ";
  std::cerr << a[0] << " , " << a[1] << " , ... " <<  a[NNN*NNN-1] << std::endl;
  std::cerr << "\tlast col: ";
  std::cerr << a[NNN*NNN] << " , " << a[NNN*NNN+1] << " , ... " <<  a[NNN*(NNN+1)-1] << std::endl;
//  std::cerr << "\tlast col line-based: ";
//  std::cerr << a[NNN] << " , " << a[2*NNN+1] << " , ... " <<  a[NNN*(NNN+1)-1] << std::endl;
}

static void
array_print_H_full(const complex *a)
{
  std::cerr << "Full Hx real[][]: " << std::setprecision(20) << std::endl;
  for (unsigned r=0; r < NNN; ++r) {
    for (unsigned c=0; c < NNN; ++c)
      std::cerr << a[r+c*NNN].real() << "  ";
    std::cerr << std::endl;
  }
  
  std::cerr << "Full Hx imag[][]: " << std::setprecision(20) << std::endl;
  for (unsigned r=0; r < NNN; ++r) {
    for (unsigned c=0; c < NNN; ++c)
      std::cerr << a[r+c*NNN].imag() << "  ";
    std::cerr << std::endl;
  }
  
  std::cerr << "Full Ht real[]:" << std::endl;
  for (unsigned r=0; r < NNN; ++r) {
      std::cerr << a[r+NNN*NNN].real() << "  ";
    std::cerr << std::endl;
  }
  
  std::cerr << "Full Ht imag[]:" << std::endl;
  for (unsigned r=0; r < NNN; ++r) {
      std::cerr << a[r+NNN*NNN].imag() << "  ";
    std::cerr << std::endl;
  }
}

static inline void 
array_multiply_scalar_to_self(complex *__restrict__ a, complex b)
{
  // TODO: optimize. Idea: use Eigen's parallizable structures or VNL SSE
  for (unsigned i = 0; i < NNN; ++i, ++a) *a = *a * b;
}

static inline void
array_negate_self(complex * __restrict__ a)
{
  // TODO: optimize. Idea: use Eigen's parallizable structures or VNL SSE
  for (unsigned i = 0; i < NNN; ++i, ++a) *a = -*a;
}


static inline void 
array_multiply_self(complex * __restrict__ a, const complex * __restrict__ b)
{
  for (unsigned int i=0; i < NNN; ++i,++a,++b) *a *= *b;
}

static inline void 
array_add_to_self(complex * __restrict__ a, const complex * __restrict__ b)
{
  for (unsigned int i=0; i < NNN; ++i,++a,++b) *a += *b;
}

static inline void 
array_add_to_self_NNNplus1(complex * __restrict__ a, const complex * __restrict__ b)
{
  for (unsigned int i=0; i < NNNPLUS1; ++i,++a,++b) *a += *b;
}

static inline void 
array_add_scalar_to_self(complex * __restrict__ a, complex b)
{
  for (unsigned int i=0; i < NNN; ++i,++a) *a += b;
}

static inline void 
array_copy(
  const complex * __restrict__ a,
  complex * __restrict__ b)
{
  // for (int i = 0; i < n; i++, a++, b++) *b = *a;
  memcpy(b, a, NNN*sizeof(complex));
}

static inline void 
array_copy_NNNplus1(
  const complex * __restrict__ a,
  complex * __restrict__ b)
{
  // for (int i = 0; i < n; i++, a++, b++) *b = *a;
  memcpy(b, a, NNNPLUS1*sizeof(complex));
}


static inline double
array_norm2(const complex *__restrict__ a)
{
  double val = 0;
  complex const* __restrict__ end = a+NNN;
  while (a != end)
    val += std::norm(*a++);
  return val;
}

#define linear linear_eigen2

// THE MEAT //////////////////////////////////////////////////////////////////
// t: tracker settings
// s_sols: start sols      
// params: params of target as specialized homotopy params - P01 in SolveChicago
// compute solutions sol_min...sol_max-1 within NSOLS
//
// TODO: template min, max
// 
unsigned    
ptrack_subset(const TrackerSettings<double> &s, const complex s_sols[NNN*NSOLS], const complex params[2*NPARAMS], Solution raw_solutions[NSOLS], unsigned sol_min, unsigned sol_max)
{
  complex Hxt[NNNPLUS1 * NNN]; 
  complex x0t0xtblock[2*NNNPLUS1];
  complex dxdt[NNNPLUS1];
  complex dxi[NNN];
  complex *x0t0 = x0t0xtblock;  // t = real running in [0,1]
  complex *x0 = x0t0;
  double  *t0 = (double *) (x0t0 + NNN);
  complex *xt = x0t0xtblock + NNNPLUS1; 
  complex *x1t1 = xt;  // reusing xt's space to represent x1t1
  complex *const HxH=Hxt;  // HxH is reusing Hxt
  complex *const dx = dxdt;
  const complex *const RHS = Hxt + NNN2;  // Hx or Ht, same storage
  complex *const LHS = Hxt;
  complex *const dx4 = dx;   // reuse dx for dx4
  double *const dt = (double *)(dxdt + NNN);
  const double &t_step = s.init_dt_;  // initial step
  using namespace Eigen; // only used for linear solve
  Map<Matrix<complex, NNN, 1> > dxi_eigen(dxi);
  Map<Matrix<complex, NNN, 1> > dx4_eigen(dx4);
  Map<Matrix<complex, NNN, 1> > &dx_eigen = dx4_eigen;
  Map<const Matrix<complex, NNN, NNN> > AA((complex *)Hxt,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > bb(RHS);

  Solution *t_s = raw_solutions + sol_min;  // current target solution
  const complex* __restrict__ s_s = s_sols + sol_min*NNN;    // current start solution
  // for (unsigned sol_n = 0; sol_n < NSOLS; ++sol_n) { // outer loop
  for (unsigned sol_n = sol_min; sol_n < sol_max; ++sol_n) { // solution loop
    t_s->status = PROCESSING;
    bool end_zone = false;
    array_copy(s_s, x0);
    *t0 = 0;
    *dt = t_step;
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
          dx1 := solveHxTimesDXequalsMinusHt(x0,t0);
          dx2 := solveHxTimesDXequalsMinusHt(x0+(1/2)*dx1*dt,t0+(1/2)*dt);
          dx3 := solveHxTimesDXequalsMinusHt(x0+(1/2)*dx2*dt,t0+(1/2)*dt);
          dx4 := solveHxTimesDXequalsMinusHt(x0+dx3*dt,t0+dt);
          (1/6)*dt*(dx1+2*dx2+2*dx3+dx4) */
      array_copy_NNNplus1(x0t0, xt);

      // dx1
      evaluate_Hxt(xt, params, Hxt); // Outputs Hxt
      PartialPivLU<Matrix<complex, NNN, NNN> > lu(AA);
      dx4_eigen = lu.solve(bb);
      
      // dx2
      const complex one_half_dt = *dt*0.5;
      array_multiply_scalar_to_self(dx4, one_half_dt);
      array_add_to_self(xt, dx4);
      array_multiply_scalar_to_self(dx4, 2);
      xt[NNN] += one_half_dt;  // t0+.5dt
      evaluate_Hxt(xt, params, Hxt);
      dxi_eigen = lu.compute(AA).solve(bb);

      // dx3
      array_multiply_scalar_to_self(dxi, one_half_dt);
      array_copy(x0t0, xt);
      array_add_to_self(xt, dxi);
      array_multiply_scalar_to_self(dxi, 4);
      array_add_to_self(dx4, dxi);
      evaluate_Hxt(xt, params, Hxt);
      dxi_eigen = lu.compute(AA).solve(bb);

      // dx4
      array_multiply_scalar_to_self(dxi, *dt);
      array_copy_NNNplus1(x0t0, xt);
      array_add_to_self(xt, dxi);
      array_multiply_scalar_to_self(dxi, 2);
      array_add_to_self(dx4, dxi);
      xt[NNN] = *t0 + *dt;               // t0+dt
      evaluate_Hxt(xt, params, Hxt);
      dxi_eigen = lu.compute(AA).solve(bb);
      array_multiply_scalar_to_self(dxi, *dt);
      array_add_to_self(dx4, dxi);
      array_multiply_scalar_to_self(dx4, 1./6.);

      // "dx1" = .5*dx1*dt, "dx2" = .5*dx2*dt, "dx3" = dx3*dt. Eigen vectorizes this:
      // dx4_eigen = (dx4_eigen* *dt + dx1_eigen*2 + dx2_eigen*4 + dx3_eigen*2)*(1./6.);
      
      // make prediction
      array_copy_NNNplus1(x0t0, x1t1);
      array_add_to_self_NNNplus1(x1t1, dxdt);
      
      /// CORRECTOR ///
      unsigned n_corr_steps = 0;
      bool is_successful;
      do {
        ++n_corr_steps;
        evaluate_HxH(x1t1, params, HxH);
        dx_eigen = lu.compute(AA).solve(bb);
        array_add_to_self(x1t1, dx);
        is_successful = array_norm2(dx) < s.epsilon2_ * array_norm2(x1t1);
      } while (!is_successful && n_corr_steps < s.max_corr_steps_);
      
      if (!is_successful) { // predictor failure
        predictor_successes = 0;
        *dt *= s.dt_decrease_factor_;
        if (*dt < s.min_dt_) t_s->status = MIN_STEP_FAILED; // slight difference to SLP-imp.hpp:612
      } else { // predictor success
        ++predictor_successes;
        std::swap(x1t1,x0t0);
        x0 = x0t0;
        t0 = (double *) (x0t0 + NNN);
        xt = x1t1;
        if (predictor_successes >= s.num_successes_before_increase_) {
          predictor_successes = 0;
          *dt *= s.dt_increase_factor_;
        }
      }
      if (array_norm2(x0) > s.infinity_threshold2_)
        t_s->status = INFINITY_FAILED;
    } // while (t loop)
    // record the solution
    array_copy(x0, t_s->x);
    t_s->t = *t0;
    if (t_s->status == PROCESSING) t_s->status = REGULAR;
    ++t_s;
    s_s += NNN;
  } // outer solution loop

  return 0;  // in the past, n_sols was returned, which now is NSOLS
}
