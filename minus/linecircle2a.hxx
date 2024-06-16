#ifndef linecircle2a_hxx_
#define linecircle2a_hxx_
// to be included at the end of minus.hxx

namespace MiNuS {
  
template <typename F>
struct eval<linecircle2a, F> {
  static void inline  __attribute__((always_inline)) Hxt(const C<F> * __restrict x /*x, t*/,    const C<F> * __restrict params, C<F> * __restrict y /*Hxt*/);
  static void inline  __attribute__((always_inline)) HxH(const C<F> * __restrict x /*x and t*/, const C<F> * __restrict params, C<F> * __restrict y /*HxH*/);
  static void inline  __attribute__((always_inline)) Hxt_constants(const C<F> * __restrict x /*x, t*/,    const C<F> * __restrict params, C<F> * __restrict y /*Hxt*/);
  static void inline  __attribute__((always_inline)) HxH_constants(const C<F> * __restrict x /*x, t*/,    const C<F> * __restrict params, C<F> * __restrict y /*HxH*/);
  static void inline  __attribute__((always_inline)) HxH_constants_all_sols(const C<F> * __restrict x /*x, t*/,    const C<F> * __restrict params, C<F> * __restrict y /*HxH*/);
};

#include "linecircle2a-Hxt.hxx"
#include "linecircle2a-HxH.hxx"

// Problem and Formulation Paramers --------------------------------------------

} // namespace minus

#include "linecircle2a-io.h"

namespace MiNuS {

// --- gammify -----------------------------------------------------------------
//
// 9 random complex numbers (rand x + i rand y), non unit, seemingly uniform
// Corresponding to the 9 pairwise lines. Seems unit is a better idea
// numerically.
//
// gamma1 .. gamma9
// 
// diag0 Generate a 3*9 = 27 entry thing by duplicationg gammas
// gamma1
// gamma1
// gamma1
// gamma2
// gamma2
// gamma2
// ...
// gamma9
// gamma9
// gamma9
//
//  tripleIntersections := {{0,3,9},{0+1,3+1,9+1},{0+2,3+2,9+2},
//  {0,6,12},{0+1,6+1,12+1},{0+2,6+2,12+2}};
//
//  for each triple intersection i
//    Get the first two (point-point) lines
//    
//    diag1(i) = conjugate(gammas(tripleIntersection(i)(0)))
//    diag1(i+1) = conjugate(gammas(tripleIntersection(i)(1)))
//    
//  diag2 := 7 times a fixed random(); -- t chart gamma
//  diag3 := 5 times a fixed random(); -- q chart, cam 2, gamma
//  diag4 := 5 times ...               -- q chart, cam 3, gamma
//  p' := (diag0|diag1|diag2|diag3|diag4).*p;
//  total    27   12    7      5    5 = 56
//
template <typename F>
inline void 
minus_io<linecircle2a, F>::
gammify(C<F> * __restrict params /*[ chicago: M::nparams]*/)
{
  // Noob ----------------------------------------------------------------------
  // No - op

  // Pro -----------------------------------------------------------------------
  // see chigago14a.hxx gammify
}

// gammify_start_params: set to false if your start parameters are already
// gammified. 
template <typename F>
inline void
minus_io<linecircle2a, F>::
get_params_start_target(
    F plines[/*15 for chicago*/][io::ncoords2d_h], 
    C<F> * __restrict params/*[static 2*M::nparams]*/,
    bool gammify_start_params)
{
  // the user provides the start params in the first half of params.
  // we fill the second half and gammify both.
  lines2params(plines, params+M::f::nparams);
  if (gammify_start_params)
    gammify(params);
  gammify(params+M::f::nparams);
}

} // namespace minus

// Highlevel solver interface - Class minus ------------------------------------

#include <thread>
#include "linecircle2a-default-data.h"

namespace MiNuS {

// Input: data (e.g., point correspondences)
// 
// Output: solutions
// Output: nsols, the number of solutions
// Output: id_sols
// a vector of the ids of the points that lead to each solution:
// So each solution is actually cameras[id_sols[i]][view_id][:][:], for i=1 to
// nsols.
//
// This design is for cache speed
// 
// The cameras array is fixed in size to NSOLS which is the max
// number of solutions, which perfectly fits in memory. The caller must pass an
// array with that minimum.
// 
// returns false in case of numerical failure to find valid real solutions
// 
template <typename F>
inline bool
minus<linecircle2a, F>::solve(
    const C<F> params_final, // p1 in linecircle2a-end.m2 
    F solutions[M::nsols],  // first camera is always [I | 0]
    unsigned id_sols[M::nsols],
    unsigned *nsols_final,
    unsigned nthreads
    )
{
  typedef minus_data<linecircle2a,F> data;
  alignas(64) C<F> params[2*M::f::nparams];
  memcpy(params, data::params_start_target_, M::f::nparams*sizeof(C<F>));
  
  // You may want to convert your data e.g. points into params here

  alignas(64) typename M::solution solutions[M::nsols];
  alignas(64) typename M::track_settings settings = M::DEFAULT;

  unsigned npaths_per_thread = M::nsols/nthreads;
  assert(M::nsols % nthreads == 0);
  

  // TODO: improve https://stackoverflow.com/questions/55908791/creating-100-threads-in-c
  std::vector<std::thread> t; 
  t.reserve(nthreads);
  { // TODO: smarter way to select start solutions
    for (unsigned i = 0; i < nthreads; ++i)
      t.emplace_back(M::track, settings, data::start_sols_, params, solutions, 
          npaths_per_thread*i, npaths_per_thread*(i+1));

     for (auto &thr : t)
          thr.join();
  }
  if (!io::has_valid_solutions(solutions))
    return false;
 
  // you may want to decode solutions into a standard format, eg convert quaternions to 3x4 cams
  io::all_real_solutions(solutions, solutions_real, id_sols, nsols_final);

  return true;
}

//
// Performs tests to see if there are potentially valid solutions,
// without making use of ground truth. 
// 
template <typename F>
inline bool 
minus_io<linecircle2a, F>::
has_valid_solutions(const typename M::solution solutions[M::nsols])
{
  typedef minus_array<M::nve,F> v;
  F real_solution[M::nve];
  for (unsigned sol = 0; sol < M::nsols; ++sol) 
    if (solutions[sol].status == M::REGULAR && v::get_real(solutions[sol].x, real_solution))
      return true;
  return false;
}

} // namespace minus
#endif // linecircle2a_hxx_
