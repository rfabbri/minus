// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date Created: Fri Mar  1 13:46:21 -03 2019
// 
#include "expminus.h"
#include <cstdio>
#include <iostream>
#include <iomanip>
#include "Eigen/Dense"
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
  complex const* end = a+NNN;
  while (a != end)
    val += std::norm(*a++);
  return val;
}

#define linear linear_eigen2

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
  complex *x0 = x0t0; complex *t0 = x0t0 + NNN;
  //  complex* x1 =  x1t1;
  //  complex* t1 = x1t1+NNN;
  complex dxdt[NNNPLUS1], *dx = dxdt, *dt = dxdt + NNN;
  complex Hxt[NNNPLUS1 * NNN], *HxH=Hxt;  // HxH is reusing Hxt
  complex *LHS, *RHS;
  complex xt[NNNPLUS1];
  complex dx1[NNN], dx2[NNN], dx3[NNN], dx4[NNN];
  complex *x1t1 = xt;  // reusing xt's space to represent x1t1

  SolutionExp* t_s = raw_solutions;  // current target solution
  const complex* s_s = s_sols;    // current start solution
  // complex *p = t_s->path;
  // #pragma parallel for 
  for (unsigned sol_n = 0; sol_n < NSOLS; ++sol_n) { // outer loop
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
    while (t_s->status == PROCESSING && 1 - t0->real() > the_smallest_number) {
      if (1 - t0->real() <= s->end_zone_factor_ + the_smallest_number && !end_zone)
        end_zone = true; // TODO: see if this path coincides with any other path on entry to the end zone
      if (end_zone) {
          if (dt->real() > 1 - t0->real()) *dt = 1 - t0->real();
      } else if (dt->real() > 1 - s->end_zone_factor_ - t0->real()) *dt = 1 - s->end_zone_factor_ - t0->real();
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
      #ifdef M_VERBOSE
      std::cerr << "\tEntered while loop" << std::endl;
      std::cerr << "t0 real" << t0->real() << std::endl;
      #endif
      
      bool Axb_success = true;
      array_copy_NNNplus1(x0t0, xt);

      #ifdef M_VERBOSE
      array_print_NNNplus1("x0t0", x0t0);
      array_print_NNNplus1("xt", xt);
      array_print_n("xt(params)", xt+NNNPLUS1,2*NPARAMS);
      #endif
      
      // dx1
      evaluate_Hxt(xt, params, Hxt); // Outputs Hxt
      #ifdef M_VERBOSE
      array_print_H("Hxt",Hxt);

      if (std::abs(Hxt[0] - complex(-.0142658,+.0314315)) <= dbgtol &&
          std::abs(Hxt[NNN*NNN-1] - complex(0,0)) <= dbgtol)
        std::cerr << "Hx from Hxt PASS\n";
      else
        std::cerr << "Hx from Hxt FAIL\n";
      
      if (std::abs(Hxt[NNN*NNN] - complex(.0901308,.0653859)) <= dbgtol &&
          std::abs(Hxt[NNN*NNN+NNN-1] - complex(-.0644011,-.470543)) <= dbgtol)
        std::cerr << "Ht from Hxt PASS \n";
      else {
        std::cerr << "Ht from Hxt FAIL\n";
        std::cerr << "\t" << Hxt[NNN*NNN] - complex(.0901308,.0653859)
          << Hxt[NNN*NNN+NNN-1] - complex(-.0644011,-.470543);
      }
      #endif
      
      LHS = Hxt;
      RHS = Hxt + NNN2;
      #ifdef M_VERBOSE
      array_print_n("LHS",LHS, NNN*NNN);
      array_print("minus RHS",RHS);
      #endif
      // was: solve_via_lapack_without_transposition(n, LHS, 1, RHS, dx1);
      Axb_success &= linear(LHS,RHS,dx1);
      #ifdef M_VERBOSE
      array_print("dx1",dx1);
      // TODO: once this is working, use eigen2 for LU partial pivots, faster.
      // this is QR full collumn pivoting, which is as fast as lapack LU
      
      if (std::abs(dx1[0] - complex(120.129,-61.1552)) <= dbgtol &&
          std::abs(dx1[NNN-1] - complex(25.3738,3.0885)) <= dbgtol)
        std::cerr << "dx1 PASS\n";
      else
        std::cerr << "dx1 FAIL\n";

      #endif
      
      // dx2
      complex one_half_dt = one_half* *dt;
      array_multiply_scalar_to_self(dx1, one_half_dt);
      array_add_to_self(xt, dx1); // x0+.5dx1*dt
      xt[NNN] += one_half_dt;  // t0+.5dt
      //
      evaluate_Hxt(xt, params, Hxt);

      
      LHS = Hxt;
      RHS = Hxt + NNN2;
      Axb_success &= linear(LHS,RHS,dx2);
      #ifdef M_VERBOSE
      std::cerr << "second eval ---------" << std::endl;
      array_print_NNNplus1("xt", xt);
      array_print_H("Hxt",Hxt);
      array_print("dx2", dx2);
      
      if (std::abs(dx2[0] - complex(414.419,-394.218)) <= dbgtol &&
          std::abs(dx2[NNN-1] - complex(-246.31,-120.222)) <= dbgtol)
        std::cerr << "dx2 PASS\n";
      else
        std::cerr << "dx2 FAIL\n";
      #endif

      // dx3
      array_multiply_scalar_to_self(dx2, one_half_dt);
      array_copy(x0t0, xt);
      array_add_to_self(xt, dx2); // x0+.5dx2*dt
      // xt[n] += one_half*(*dt); // t0+.5dt (SAME)
      //
      evaluate_Hxt(xt, params, Hxt);
      #ifdef M_VERBOSE
      array_print_H_full(Hxt);
      #endif
      LHS = Hxt;
      RHS = Hxt + NNN2;
      Axb_success &= linear(LHS,RHS,dx3);
      #ifdef M_VERBOSE
      std::cerr << "third eval ---------" << std::endl;
      array_print_NNNplus1("xt", xt);
      array_print_H("Hxt",Hxt);
      array_print("dx3", dx3);
      if (std::abs(dx3[0] - complex(-11197.9,+14804.2)) <= 10*dbgtol &&
          std::abs(dx3[NNN-1] - complex(-17653.8,-4747.93)) <= dbgtol)
        std::cerr << "dx3 PASS\n";
      else {
        std::cerr << "dx3 FAIL\n";
        std::cerr << "\t\t" << dx3[0] - complex(-11197.9,+14804.2)
          << dx3[NNN-1] - complex(-17653.8,-4747.93) << std::endl;
      }
      #endif

      // dx4
      array_multiply_scalar_to_self(dx3, *dt);
      array_copy_NNNplus1(x0t0, xt);
      array_add_to_self(xt, dx3); // x0+dx3*dt
      xt[NNN] += *dt;               // t0+dt
      //
      evaluate_Hxt(xt, params, Hxt);
      #ifdef M_VERBOSE
      std::cerr << "fourth eval ---------" << std::endl;
      array_print_NNNplus1("xt", xt);
        std::cerr << "\txt full: {" << std::setprecision(20);
      for (unsigned id=0; id < NNN; ++id) {
        std::cerr << "{" << xt[id] << "} ,";
      }
      std::cerr << "{" << xt[NNN]  << "}};" << std::endl;
      
      if (std::abs(Hxt[0] - complex(72802.2,-81551.9)) <= 10*dbgtol &&
          std::abs(Hxt[NNN*NNN-1] - complex(0,0)) <= dbgtol)
        std::cerr << "Hx from Hxt PASS\n";
      else
        std::cerr << "Hx from Hxt FAIL\n";
      
      if (std::abs(Hxt[NNN*NNN] - complex(-248323000,+163093000)) <= 1e3 &&
          std::abs(Hxt[NNN*NNN+NNN-1] - complex(-.221521,-.47423)) <= dbgtol)
        std::cerr << "Ht from Hxt PASS \n";
      else
        std::cerr << "Ht from Hxt FAIL\n";
      array_print_H("Hxt",Hxt);
      array_print_H_full(Hxt);
      #endif
      LHS = Hxt;
      RHS = Hxt + NNN2;
      Axb_success &= linear(LHS,RHS,dx4);
      #ifdef M_VERBOSE
      array_print_NNNplus1("xt", xt);
      array_print("dx4", dx4);
      
      if (std::abs(dx4[0] - complex(885846,811668)) <= 1000*dbgtol &&
          std::abs(dx4[NNN-1] - complex(823211,-2145640)) <= 1000*dbgtol)
        std::cerr << "dx4 PASS\n";
      else {
        std::cerr << "dx4 FAIL\n";
        std::cerr << "\t\t" << dx4[0] - complex(885846,811668)
          << dx4[NNN-1] -  complex(823211,-2145640) << std::endl;
      }
      #endif

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
        #ifdef M_VERBOSE
        std::cerr << "\tCorrection step " << n_corr_steps << std::endl;
        array_print_NNNplus1("x1t1",x1t1);
        
        if (n_corr_steps == 1) {
          if (std::abs(HxH[0] - complex(26113200,11590100)) <= dbgtol &&
              std::abs(HxH[NNN*NNN-1] - complex(0,0)) <= dbgtol)
            std::cerr << "Hx from HxH PASS\n";
          else
            std::cerr << "Hx from HxH FAIL\n";
          
          if (std::abs(HxH[NNN*NNN] - complex(+1.741e11,-15110700000)) <= dbgtol &&
              std::abs(HxH[NNN*NNN+NNN-1] - complex(-3.71869e-13,+4.6629e-15)) <= dbgtol)
            std::cerr << "H from HxH PASS \n";
          else
            std::cerr << "H from HxH FAIL \n";
        }
        array_print_H("HxH",HxH);
        #endif
        
        
        LHS = HxH;         // Hx
        RHS = HxH + NNN2;  // -H
        Axb_success &= linear(LHS,RHS,dx);
        #ifdef M_VERBOSE
        array_print_n("LHS corr",LHS,NNN2);
        array_print("RHS corr",RHS);
        array_print("dx",dx);
        if (n_corr_steps == 1) {
          if (std::abs(dx[0] - complex(-2206.68,-2378.66)) <= dbgtol &&
              std::abs(dx[NNN-1] - complex(-1731.72,+6593.05)) <= dbgtol)
            std::cerr << "dx from first correction PASS \n";
          else
            std::cerr << "dx from first correction FAIL \n";
        }
        #endif
        array_add_to_self(x1t1, dx);
        is_successful = array_norm2(dx) < s->epsilon2_ * array_norm2(x1t1);
        #ifdef M_VERBOSE
        // printf("c: |dx|^2 = %lf\n",
        // norm2_complex_array<ComplexField>(n,dx));
        #endif
      } while (!is_successful && n_corr_steps < s->max_corr_steps_);
      
      if (!is_successful) { // predictor failure
        predictor_successes = 0;
        *dt *= s->dt_decrease_factor_;
        if (dt->real() < s->min_dt_) t_s->status = MIN_STEP_FAILED; // slight difference to SLP-imp.hpp:612
        #ifdef M_VERBOSE
        std::cerr << "\tPred failure" << std::endl;
        std::cerr << "decreasing dt: " << *dt;
        #endif
      } else { // predictor success
        ++predictor_successes; count++;
        array_copy_NNNplus1(x1t1, x0t0);
        if (predictor_successes >= s->num_successes_before_increase_) {
          predictor_successes = 0;
          *dt *= s->dt_increase_factor_;
          #ifdef M_VERBOSE
          std::cerr << "\tPred success" << std::endl;
          std::cerr << "increasing dt: " << *dt;
          #endif
        }
      }
      if (array_norm2(x0) > s->infinity_threshold2_)
        t_s->status = INFINITY_FAILED;
      if (!Axb_success) t_s->status = SINGULAR;
      #ifdef M_VERBOSE
      std::cerr << "\tAxb_success" << Axb_success << std::endl;
      #endif
      if (count < MAX_NSAMPLES) {
//        if (count % 10 == 0) { 
         // array_copy(x0, p);
          // p += NNN;
        // }
      }
    } // while 
    // record the solution
    array_copy(x0, t_s->x);
    t_s->t = t0->real();
    if (t_s->status == PROCESSING) t_s->status = REGULAR;
    // evaluate_HxH(x0t0, HxH);
    // cond_number_via_svd(HxH /*Hx*/, t_s->cond);
    t_s->num_steps = count;
    ++t_s;
    s_s += NNN;
  } // outer solution loop

  return 0;  // in the past, n_sols was returned, which now is NSOLS
}
  
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
