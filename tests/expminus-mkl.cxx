// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date Created: Fri Mar  1 13:46:21 -03 2019
// 
#include "expminus.h"
#include <cstdio>
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#define MKL_Complex16 complex
#include <mkl.h>
//#include "Eigen/Dense"
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


// Original code solve_via_lapack_without_transposition
//
// \returns  lapack info code. 
//      >0 -> matrix is singular
//      <0 -> illegal value of an argument passed to lapack
bool linear_mkl(
    complex* A,  // NNN-by-NNN matrix of complex #s
    complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  lapack_int info;
  lapack_int permutation[NNN]; // unused
  static const lapack_int bsize = 1;
  static const lapack_int size = NNN;
  
  // TODO try to eliminate this memcpy and trash the original b if possible
  // memcpy  b -> x          NNN elements
  std::memcpy(x, b, 2*NNN*sizeof(double));

  zgesv(&size, &bsize, A, &size, permutation, x, &size, &info);

  return info == 0;
}
#if 0
// Original code solve_via_lapack_without_transposition
//
// \returns  lapack info code. 
//      >0 -> matrix is singular
//      <0 -> illegal value of an argument passed to lapack
bool linear(
    complex* A,  // NNN-by-NNN matrix of complex #s
    complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  static lapack_int permutation[NNN]; // unused
  static const lapack_int size = NNN;
  
  // TODO try to eliminate this memcpy and trash the original b if possible
  // memcpy  b -> x          NNN elements
  std::memcpy(x, b, 2*NNN*sizeof(double));

  return LAPACKE_zgesv(LAPACK_COL_MAJOR,size, 1, A, size, permutation, x, size);
}
#endif

/*
bool
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
  

bool 
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

bool 
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
*/

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
array_multiply_scalar_to_self(complex *a, complex b)
{
  // TODO: optimize. Idea: use Eigen's parallizable structures or VNL SSE
  for (unsigned i = 0; i < NNN; ++i, ++a) *a = *a * b;
}

static inline void
array_negate_self(complex *a)
{
  // TODO: optimize. Idea: use Eigen's parallizable structures or VNL SSE
  for (unsigned i = 0; i < NNN; ++i, ++a) *a = -*a;
}


static inline void 
array_multiply_self(complex *a, const complex *b)
{
  for (unsigned int i=0; i < NNN; ++i,++a,++b) *a *= *b;
}

static inline void 
array_add_to_self(complex *a, const complex *b)
{
  for (unsigned int i=0; i < NNN; ++i,++a,++b) *a += *b;
}

static inline void 
array_add_to_self_NNNplus1(complex *a, const complex *b)
{
  for (unsigned int i=0; i < NNNPLUS1; ++i,++a,++b) *a += *b;
}

static inline void 
array_add_scalar_to_self(complex *a, complex b)
{
  for (unsigned int i=0; i < NNN; ++i,++a) *a += b;
}

static inline void 
array_copy(
  const complex *a,
  complex *b)
{
  // for (int i = 0; i < n; i++, a++, b++) *b = *a;
  memcpy(b, a, NNN*sizeof(complex));
}

static inline void 
array_copy_NNNplus1(
  const complex *a,
  complex *b)
{
  // for (int i = 0; i < n; i++, a++, b++) *b = *a;
  memcpy(b, a, NNNPLUS1*sizeof(complex));
}


static inline double
array_norm2(const complex *a)
{
  double val = 0;
  complex const* const end = a+NNN;
  while (a != end)
    val += std::norm(*a++);
  return val;
  // return Eigen::Map<const Eigen::Matrix<complex, NNN, 1> > (a).squaredNorm();
}

#define linear linear_mkl

//void 
//array_copy_n(size_t n, const complex *a, complex *b)
//{
//   for (int i = 0; i < n; i++, a++, b++) *b = *a;
//  memcpy(b, a, n*sizeof(complex));
//}

// THE MEAT //////////////////////////////////////////////////////////////////
// t: tracker settings
// s_sols: start sols      
// params: params of target as specialized homotopy params - P01 in SolveChicago
unsigned    
exp_ptrack(const TrackerSettings *s, const complex s_sols[NNN*NSOLS], const complex params[2*NPARAMS], SolutionExp raw_solutions[NSOLS])
{
  // TODO: test by making variables static for a second run, some of these arrays may have to be zeroed
  // One huge clear instruction will work as they are sequential in mem.
  const complex one_half(0.5, 0);
  const double t_step = s->init_dt_;  // initial step
  complex x0t0[NNNPLUS1];  // t = real running in [0,1]
  complex *x0 = x0t0; double *t0 = (double *) (x0t0 + NNN);
  //  complex* x1 =  x1t1;
  //  complex* t1 = x1t1+NNN;
  complex dxdt[NNNPLUS1], *dx = dxdt, *dt = dxdt + NNN;
  complex Hxt[NNNPLUS1 * NNN], *HxH=Hxt;  // HxH is reusing Hxt
  complex *RHS = Hxt + NNN2;  // Hx or Ht, same storage
  complex *LHS = Hxt;
  complex xt[NNNPLUS1];
  complex dx1[NNN], dx2[NNN], dx3[NNN], dx4[NNN];
  complex *x1t1 = xt;  // reusing xt's space to represent x1t1
  
  // using namespace Eigen;
  //Map<const Matrix<complex, NNN, NNN> > eHx(Hxt,NNN,NNN);  // accessors for the data

  SolutionExp* t_s = raw_solutions;  // current target solution
  const complex* s_s = s_sols;    // current start solution
  unsigned solve_linear_calls = 0;
  unsigned num_total_iter = 0;
  // complex *p = t_s->path;
  for (unsigned sol_n = 0; sol_n < NSOLS; ++sol_n) { // outer loop
//     std::ostringstream s_str;
//     s_str << std::setw(std::ceil(3)) << std::setfill('0') << sol_n;
//     std::ofstream fpaths("path-solution-" + s_str.str() ,std::ios::out);
//     fpaths << std::setprecision(20);
//     std::ofstream finfo("info-solution-" + s_str.str(), std::ios::out);
//     finfo << std::setprecision(20);
     
    #ifdef M_VERBOSE
    std::cerr << "Trying solution #" << sol_n << std::endl;
    #endif
    t_s->make(s_s);  // cook a Solution: copy s_s to start node of path
    t_s->status = PROCESSING;
    bool end_zone = false;
    array_copy(s_s, x0);
    *t0 = 0;
    *dt = t_step;
    unsigned predictor_successes = 0, count = 0;  // number of steps

    // print start solution

    // PASS array_print("s_s",s_s);
    
    // track H(x,t) for t in [0,1]
    while (t_s->status == PROCESSING && 1 - *t0 > the_smallest_number) {
      if (!end_zone && 1 - *t0 <= s->end_zone_factor_ + the_smallest_number)
        end_zone = true; // TODO: see if this path coincides with any other path on entry to the end zone
      if (end_zone) {
          if (dt->real() > 1 - *t0) *dt = 1 - *t0;
      } else if (dt->real() > 1 - s->end_zone_factor_ - *t0) *dt = 1 - s->end_zone_factor_ - *t0;
      // PREDICTOR in: x0t0,dt
      //           out: dx
      // Runge Kutta
      /*  top-level code for Runge-Kutta-4
          dx1 := solveHxTimesDXequalsMinusHt(x0,t0);
          dx2 := solveHxTimesDXequalsMinusHt(x0+(1/2)*dx1*dt,t0+(1/2)*dt);
          dx3 := solveHxTimesDXequalsMinusHt(x0+(1/2)*dx2*dt,t0+(1/2)*dt);
          dx4 := solveHxTimesDXequalsMinusHt(x0+dx3*dt,t0+dt);
          (1/6)*dt*(dx1+2*dx2+2*dx3+dx4)
      */
      
      bool Axb_success = true;
      array_copy_NNNplus1(x0t0, xt);

      // dx1
      evaluate_Hxt(xt, params, Hxt); // Outputs Hxt
      // was: solve_via_lapack_without_transposition(n, LHS, 1, RHS, dx1);
      Axb_success &= linear(LHS,RHS,dx1);
      solve_linear_calls++;
      
//      Matrix<complex, NNN, NNN> Hxinv = eHx.inverse();
//      {
//        Map<Matrix<complex, NNN, 1> > xx(dx1);
//        Map<const Matrix<complex, NNN, 1> > bb(RHS);
//        xx = Hxinv * bb;
//      }
      
      // dx2
      complex one_half_dt = one_half* *dt;
      array_multiply_scalar_to_self(dx1, one_half_dt);
      array_add_to_self(xt, dx1); // x0+.5dx1*dt
      xt[NNN] += one_half_dt;  // t0+.5dt
      //
      evaluate_Hxt(xt, params, Hxt);
      Axb_success &= linear(LHS,RHS,dx2);
      solve_linear_calls++;
      //Hxinv = eHx.inverse();
//      {
//        Map<Matrix<complex, NNN, 1> > xx(dx2);
//        Map<const Matrix<complex, NNN, NNN> > AA(LHS,NNN,NNN);  // accessors for the data
//        Map<const Matrix<complex, NNN, 1> > bb(RHS);
//        xx = Hxinv * bb;
//        xx += Hxinv * (bb - AA*xx);
//      }
      
      // dx3
      array_multiply_scalar_to_self(dx2, one_half_dt);
      array_copy(x0t0, xt);
      array_add_to_self(xt, dx2); // x0+.5dx2*dt
      // xt[n] += one_half*(*dt); // t0+.5dt (SAME)
      //
      evaluate_Hxt(xt, params, Hxt);
      Axb_success &= linear(LHS,RHS,dx3);
      solve_linear_calls++;
      //Hxinv = eHx.inverse();
//      {
//        Map<Matrix<complex, NNN, 1> > xx(dx3);
//        Map<const Matrix<complex, NNN, NNN> > AA(LHS,NNN,NNN);  // accessors for the data
//        Map<const Matrix<complex, NNN, 1> > bb(RHS);
//        xx = Hxinv * bb; // predict
//        xx += Hxinv * (bb - AA*xx); //correct
//      }

      // dx4
      array_multiply_scalar_to_self(dx3, *dt);
      array_copy_NNNplus1(x0t0, xt);
      array_add_to_self(xt, dx3); // x0+dx3*dt
      xt[NNN] += *dt;               // t0+dt
      //
      evaluate_Hxt(xt, params, Hxt);
      Axb_success &= linear(LHS,RHS,dx4);
      solve_linear_calls++;
      //Hxinv = eHx.inverse();
//      {
//        Map<Matrix<complex, NNN, 1> > xx(dx4);
//        Map<const Matrix<complex, NNN, NNN> > AA(LHS,NNN,NNN);  // accessors for the data
//        Map<const Matrix<complex, NNN, 1> > bb(RHS);
//        xx = Hxinv * bb;
//        xx += Hxinv * (bb - AA*xx);
//      }

      // "dx1" = .5*dx1*dt, "dx2" = .5*dx2*dt, "dx3" = dx3*dt
      // TODO: make this into a single function loop directly
      array_multiply_scalar_to_self(dx4, *dt);
      array_multiply_scalar_to_self(dx1, 2);
      array_multiply_scalar_to_self(dx2, 4);
      array_multiply_scalar_to_self(dx3, 2);
      array_add_to_self(dx4, dx1);
      array_add_to_self(dx4, dx2);
      array_add_to_self(dx4, dx3);
      array_multiply_scalar_to_self(dx4, 1./6.);
      array_copy(dx4, dx);

      // make prediction
      array_copy_NNNplus1(x0t0, x1t1);
      array_add_to_self_NNNplus1(x1t1, dxdt);
      
      // CORRECTOR
      unsigned n_corr_steps = 0;
      bool is_successful;
      do {
        ++n_corr_steps;
        ++num_total_iter;
        evaluate_HxH(x1t1, params, HxH);
        
        Axb_success &= linear(LHS,RHS,dx);
        solve_linear_calls++;
//            {
//              Map<Matrix<complex, NNN, 1> > xx(dx);
//              Map<const Matrix<complex, NNN, NNN> > AA(LHS,NNN,NNN);  // accessors for the data
//              Map<const Matrix<complex, NNN, 1> > bb(RHS);
//              xx = Hxinv * bb;
//              xx += Hxinv * (bb - AA*xx);
//            }
        array_add_to_self(x1t1, dx);
        is_successful = array_norm2(dx) < s->epsilon2_ * array_norm2(x1t1);
      } while (!is_successful && n_corr_steps < s->max_corr_steps_);
      
      if (!is_successful) { // predictor failure
        predictor_successes = 0;
        *dt *= s->dt_decrease_factor_;
        if (dt->real() < s->min_dt_) t_s->status = MIN_STEP_FAILED; // slight difference to SLP-imp.hpp:612
      } else { // predictor success
        ++predictor_successes; count++;
        array_copy_NNNplus1(x1t1, x0t0);
        if (predictor_successes >= s->num_successes_before_increase_) {
          predictor_successes = 0;
          *dt *= s->dt_increase_factor_;
        }
//        for (unsigned kk=0; kk < NNN; ++kk)
//          fpaths << x0[kk] << ' ';
//        fpaths << std::endl;
      }
      if (array_norm2(x0) > s->infinity_threshold2_)
        t_s->status = INFINITY_FAILED;
      if (!Axb_success) t_s->status = SINGULAR;
      //if (sol_n == 0) std::cerr << "dt: " << dt->real() << std::endl;
    } // while 
    // record the solution
//    for (unsigned kk=0; kk < NNN; ++kk)
//      finfo << x0[kk] << ' ';
//    finfo << " \% final solution" << std::endl;
    
    array_copy(x0, t_s->x);
    t_s->t = *t0;
//    finfo << t_s->t;
//    finfo << " \% final time (<= 1)" << std::endl;
    if (t_s->status == PROCESSING) t_s->status = REGULAR;
    t_s->num_steps = count;
//    finfo << count;
//    finfo << " \% number of steps in path" << std::endl;
//    finfo << t_s->status ;
//    finfo << " \% status: " << std::endl;
    ++t_s;
    s_s += NNN;
//    fpaths.close();
//    finfo.close();
  } // outer solution loop

  std::cout << "# Solve linear calls: " << solve_linear_calls << std::endl;
  std::cout << "# Total number of iterations: " << num_total_iter << std::endl;

  return 0;  // in the past, n_sols was returned, which now is NSOLS
}

#if 0
/*
//
// compute solutions sol_min...sol_max-1 within NSOLS
// 
void
exp_ptrack_subset(const TrackerSettings *s, const complex s_sols[NNN*NSOLS], const complex params[2*NPARAMS], SolutionExp raw_solutions[NSOLS], unsigned sol_min, unsigned sol_max)
{
  // TODO: test by making variables static for a second run, some of these arrays may have to be zeroed
  // One huge clear instruction will work as they are sequential in mem.
  const complex one_half(0.5, 0);
  const double t_step = s->init_dt_;  // initial step
  complex x0t0[NNNPLUS1];  // t = real running in [0,1]
  complex *x0 = x0t0; double *t0 = (double *) (x0t0 + NNN);
  //  complex* x1 =  x1t1;
  //  complex* t1 = x1t1+NNN;
  complex dxdt[NNNPLUS1], *dx = dxdt, *dt = dxdt + NNN;
  complex Hxt[NNNPLUS1 * NNN], *HxH=Hxt;  // HxH is reusing Hxt
  const complex *RHS = Hxt + NNN2;  // Hx or Ht, same storage
  const complex *LHS = Hxt;
  complex xt[NNNPLUS1];
  complex dx1[NNN], dx2[NNN], dx3[NNN], dx4[NNN];
  complex *x1t1 = xt;  // reusing xt's space to represent x1t1
  
//  using namespace Eigen;
//  Map<const Matrix<complex, NNN, NNN> > eHx(Hxt,NNN,NNN);  // accessors for the data

  SolutionExp* t_s = raw_solutions+sol_min;  // current target solution
  const complex* s_s = s_sols+sol_min*NNN;    // current start solution
  // complex *p = t_s->path;
  for (unsigned sol_n = sol_min; sol_n < sol_max; ++sol_n) { // outer loop
//     std::ostringstream s_str;
//     s_str << std::setw(std::ceil(3)) << std::setfill('0') << sol_n;
//     std::ofstream fpaths("path-solution-" + s_str.str() ,std::ios::out);
//     fpaths << std::setprecision(20);
//     std::ofstream finfo("info-solution-" + s_str.str(), std::ios::out);
//     finfo << std::setprecision(20);
     
    #ifdef M_VERBOSE
    std::cerr << "Trying solution #" << sol_n << std::endl;
    #endif
    t_s->make(s_s);  // cook a Solution: copy s_s to start node of path
    t_s->status = PROCESSING;
    bool end_zone = false;
    array_copy(s_s, x0);
    *t0 = 0;
    *dt = t_step;
    unsigned predictor_successes = 0, count = 0;  // number of steps

    // print start solution

    // PASS array_print("s_s",s_s);
    
    // track H(x,t) for t in [0,1]
    while (t_s->status == PROCESSING && 1 - *t0 > the_smallest_number) {
      if (!end_zone && 1 - *t0 <= s->end_zone_factor_ + the_smallest_number)
        end_zone = true; // TODO: see if this path coincides with any other path on entry to the end zone
      if (end_zone) {
          if (dt->real() > 1 - *t0) *dt = 1 - *t0;
      } else if (dt->real() > 1 - s->end_zone_factor_ - *t0) *dt = 1 - s->end_zone_factor_ - *t0;
      // PREDICTOR in: x0t0,dt
      //           out: dx
      // Runge Kutta
      /*  top-level code for Runge-Kutta-4
          dx1 := solveHxTimesDXequalsMinusHt(x0,t0);
          dx2 := solveHxTimesDXequalsMinusHt(x0+(1/2)*dx1*dt,t0+(1/2)*dt);
          dx3 := solveHxTimesDXequalsMinusHt(x0+(1/2)*dx2*dt,t0+(1/2)*dt);
          dx4 := solveHxTimesDXequalsMinusHt(x0+dx3*dt,t0+dt);
          (1/6)*dt*(dx1+2*dx2+2*dx3+dx4)
      */
      
      bool Axb_success = true;
      array_copy_NNNplus1(x0t0, xt);

      // dx1
      evaluate_Hxt(xt, params, Hxt); // Outputs Hxt
      // was: solve_via_lapack_without_transposition(n, LHS, 1, RHS, dx1);
      Axb_success &= linear(LHS,RHS,dx1);
      
//      Matrix<complex, NNN, NNN> Hxinv = eHx.inverse();
//      {
//        Map<Matrix<complex, NNN, 1> > xx(dx1);
//        Map<const Matrix<complex, NNN, 1> > bb(RHS);
//        xx = Hxinv * bb;
//      }
      
      // dx2
      complex one_half_dt = one_half* *dt;
      array_multiply_scalar_to_self(dx1, one_half_dt);
      array_add_to_self(xt, dx1); // x0+.5dx1*dt
      xt[NNN] += one_half_dt;  // t0+.5dt
      //
      evaluate_Hxt(xt, params, Hxt);
      Axb_success &= linear(LHS,RHS,dx2);
      //Hxinv = eHx.inverse();
//      {
//        Map<Matrix<complex, NNN, 1> > xx(dx2);
//        Map<const Matrix<complex, NNN, NNN> > AA(LHS,NNN,NNN);  // accessors for the data
//        Map<const Matrix<complex, NNN, 1> > bb(RHS);
//        xx = Hxinv * bb;
//        xx += Hxinv * (bb - AA*xx);
//      }
      
      // dx3
      array_multiply_scalar_to_self(dx2, one_half_dt);
      array_copy(x0t0, xt);
      array_add_to_self(xt, dx2); // x0+.5dx2*dt
      // xt[n] += one_half*(*dt); // t0+.5dt (SAME)
      //
      evaluate_Hxt(xt, params, Hxt);
      Axb_success &= linear(LHS,RHS,dx3);
      //Hxinv = eHx.inverse();
//      {
//        Map<Matrix<complex, NNN, 1> > xx(dx3);
//        Map<const Matrix<complex, NNN, NNN> > AA(LHS,NNN,NNN);  // accessors for the data
//        Map<const Matrix<complex, NNN, 1> > bb(RHS);
//        xx = Hxinv * bb; // predict
//        xx += Hxinv * (bb - AA*xx); //correct
//      }

      // dx4
      array_multiply_scalar_to_self(dx3, *dt);
      array_copy_NNNplus1(x0t0, xt);
      array_add_to_self(xt, dx3); // x0+dx3*dt
      xt[NNN] += *dt;               // t0+dt
      //
      evaluate_Hxt(xt, params, Hxt);
      Axb_success &= linear(LHS,RHS,dx4);
      //Hxinv = eHx.inverse();
//      {
//        Map<Matrix<complex, NNN, 1> > xx(dx4);
//        Map<const Matrix<complex, NNN, NNN> > AA(LHS,NNN,NNN);  // accessors for the data
//        Map<const Matrix<complex, NNN, 1> > bb(RHS);
//        xx = Hxinv * bb;
//        xx += Hxinv * (bb - AA*xx);
//      }

      // "dx1" = .5*dx1*dt, "dx2" = .5*dx2*dt, "dx3" = dx3*dt
      // TODO: make this into a single function loop directly
      array_multiply_scalar_to_self(dx4, *dt);
      array_multiply_scalar_to_self(dx1, 2);
      array_multiply_scalar_to_self(dx2, 4);
      array_multiply_scalar_to_self(dx3, 2);
      array_add_to_self(dx4, dx1);
      array_add_to_self(dx4, dx2);
      array_add_to_self(dx4, dx3);
      array_multiply_scalar_to_self(dx4, 1./6.);
      array_copy(dx4, dx);

      // make prediction
      array_copy_NNNplus1(x0t0, x1t1);
      array_add_to_self_NNNplus1(x1t1, dxdt);
      
      // CORRECTOR
      unsigned n_corr_steps = 0;
      bool is_successful;
      do {
        ++n_corr_steps;
        evaluate_HxH(x1t1, params, HxH);
        
        Axb_success &= linear(LHS,RHS,dx);
//            {
//              Map<Matrix<complex, NNN, 1> > xx(dx);
//              Map<const Matrix<complex, NNN, NNN> > AA(LHS,NNN,NNN);  // accessors for the data
//              Map<const Matrix<complex, NNN, 1> > bb(RHS);
//              xx = Hxinv * bb;
//              xx += Hxinv * (bb - AA*xx);
//            }
        array_add_to_self(x1t1, dx);
        is_successful = array_norm2(dx) < s->epsilon2_ * array_norm2(x1t1);
      } while (!is_successful && n_corr_steps < s->max_corr_steps_);
      
      if (!is_successful) { // predictor failure
        predictor_successes = 0;
        *dt *= s->dt_decrease_factor_;
        if (dt->real() < s->min_dt_) t_s->status = MIN_STEP_FAILED; // slight difference to SLP-imp.hpp:612
      } else { // predictor success
        ++predictor_successes; count++;
        array_copy_NNNplus1(x1t1, x0t0);
        if (predictor_successes >= s->num_successes_before_increase_) {
          predictor_successes = 0;
          *dt *= s->dt_increase_factor_;
        }
//        for (unsigned kk=0; kk < NNN; ++kk)
//          fpaths << x0[kk] << ' ';
//        fpaths << std::endl;
      }
      if (array_norm2(x0) > s->infinity_threshold2_)
        t_s->status = INFINITY_FAILED;
      if (!Axb_success) t_s->status = SINGULAR;
    } // while 
    // record the solution
//    for (unsigned kk=0; kk < NNN; ++kk)
//      finfo << x0[kk] << ' ';
//    finfo << " \% final solution" << std::endl;
    
    array_copy(x0, t_s->x);
    t_s->t = *t0;
//    finfo << t_s->t;
//    finfo << " \% final time (<= 1)" << std::endl;
    if (t_s->status == PROCESSING) t_s->status = REGULAR;
    t_s->num_steps = count;
//    finfo << count;
//    finfo << " \% number of steps in path" << std::endl;
//    finfo << t_s->status ;
//    finfo << " \% status: " << std::endl;
    ++t_s;
    s_s += NNN;
//    fpaths.close();
//    finfo.close();
  } // outer solution loop
}
#endif
  
#if 0
  /*
    printf(
        "epsilon2 = %e, t_step = %lf, dt_min_dbl = %lf, dt_increase_factor_dbl "
        "= %lf, dt_decrease_factor_dbl = %lf\n",
        epsilon2,
        t_step,
        dt_min_dbl,
        dt_increase_factor_dbl,
        dt_decrease_factor_dbl);
        */
    return 0;
}
//int
//linear_bundle
#endif
