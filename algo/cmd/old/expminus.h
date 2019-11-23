#ifndef expminus_h_
#define expminus_h_
//:
// \file
// \brief Experimental minus
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date Fri Mar  1 13:37:49 -03 2019
//

#include <complex>
#include <cstring>
#include <minus.h>

// debug number of samples to output
#define MAX_NSAMPLES 100


struct SolutionExp
{
  complex path[MAX_NSAMPLES*NNN]; 
  complex x[NNN];        // array of n coordinates. Should equal path[num_steps-1];
  double t;          // last value of parameter t used
  complex start_x[NNN];  // start of the path that produced x
  SolutionStatus status;
  unsigned num_steps;  // number of steps taken along the path
  SolutionExp() : status(UNDETERMINED) { }
  void make(const complex* s_s) { memcpy(start_x, s_s, NNN*sizeof(complex)); }
};


unsigned exp_ptrack(const TrackerSettings *t, const complex s_sols[NNN*NSOLS], const complex params[NPARAMS], SolutionExp raw_solutions[NSOLS]);

void
exp_ptrack_subset(const TrackerSettings *s, const complex s_sols[NNN*NSOLS], const complex params[2*NPARAMS], SolutionExp raw_solutions[NSOLS], unsigned sol_min, unsigned sol_max);

// --- TESTING ONLY: REMOVE ---------------------------------------------------

// in .h just for testing purposes
#if 0
bool linear(
    const complex* A,  // size-by-size matrix of complex #s
    const complex* b,  // bsize-by-size RHS of Ax=b
    complex* x   // solution
    );
    */

#endif
    
bool linear_eigen(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    );

bool linear_eigen2(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    );

bool linear_eigen3(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    );
bool linear_eigen3(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    );


#endif  // expminus_h_
