#ifndef chicago14a_hxx_
#define chicago14a_hxx_
// to be included at the end of minus.hxx

namespace MiNuS {
  
template <typename F>
struct eval<chicago14a, F> {
  static void inline  __attribute__((always_inline)) Hxt(const C<F> * __restrict x /*x, t*/,    const C<F> * __restrict params, C<F> * __restrict y /*Hxt*/);
  static void inline  __attribute__((always_inline)) HxH(const C<F> * __restrict x /*x and t*/, const C<F> * __restrict params, C<F> * __restrict y /*HxH*/);
  static void inline  __attribute__((always_inline)) Hxt_constants(const C<F> * __restrict x /*x, t*/,    const C<F> * __restrict params, C<F> * __restrict y /*Hxt*/);
  static void inline  __attribute__((always_inline)) HxH_constants(const C<F> * __restrict x /*x, t*/,    const C<F> * __restrict params, C<F> * __restrict y /*HxH*/);
  static void inline  __attribute__((always_inline)) HxH_constants_all_sols(const C<F> * __restrict x /*x, t*/,    const C<F> * __restrict params, C<F> * __restrict y /*HxH*/);
};

#include "chicago14a-Hxt.hxx"
#include "chicago14a-HxH.hxx"

// Problem and Formulation Paramers --------------------------------------------

} // namespace minus

#include "chicago14a-io.h"

namespace MiNuS {

// we only use the first half of the outer
// 2*M::nparams array 
// after this fn, complex part zero, but we will use this space later
// to gammify/randomize
template <typename F>
inline void 
minus_io<chicago14a, F>::
lines2params(const F plines[pp::nvislines][io::ncoords2d_h], C<F> * __restrict params/*[static 2*M::nparams]*/)
{
  typedef minus_util<F> util;
  typedef minus_3d<F> vec;
  //    params (P1) is pF||pTriple||pChart  //  Hongyi: [pF; tripleChart; XR'; XT1'; XT2'];
  //    size           27     12      17 = 56

  // pF ----------------------------------------
  // converts 1st 9 lines to C<F> (imaginary part zero)
  // remembering: 1st 9 lines are the ones between the points (no tangents)
  // 
  // 9x3 out of the 15x3 of the pairsiwe lines, linearized as 27x1
  // Tim: pF is matrix(targetLines^{0..8},27,1);
  // Order: row-major
  const F *pl = (const F *)plines;
  for (unsigned i=0; i < 27; ++i) params[i] = pl[i];

  // pTriple ----------------------------------------
  // At points that have tangents, there are 3 lines (triple intersects)
  static unsigned constexpr triple_intersections[6][3] = 
    {{0,3,9},{0+1,3+1,9+1},{0+2,3+2,9+2},{0,6,12},{0+1,6+1,12+1},{0+2,6+2,12+2}};

  C<F> (*params_lines)[2] = (C<F> (*)[2]) (params+27);
  // express each of the 6 tangents in the basis of the other pairwise lines
  // intersecting at the same point
  for (unsigned l=0; l < 6; ++l) {
    const F *l0 = plines[triple_intersections[l][0]];
    const F *l1 = plines[triple_intersections[l][1]];
    const F *l2 = plines[triple_intersections[l][2]];
    double l0l1 = vec::dot(l0,l1), l2l0 = vec::dot(l2,l0), l2l1 = vec::dot(l2,l1);
    double l2_l0l1[3]; // cross([l0l0 l1l0 l2l0], [l0l1 l1l1 l2l1], l2_l0l1);
    {
      F v1[3], v2[3];
      v1[0] = 1;    v1[1] = l0l1; v1[2] = l2l0;
      v2[0] = l0l1; v2[1] = 1;    v2[2] = l2l1;
      vec::cross(v1, v2, l2_l0l1);
    }
    params_lines[l][0] = l2_l0l1[0]/l2_l0l1[2]; // divide by the last coord (see cross prod formula, plug direct)
    params_lines[l][1] = l2_l0l1[1]/l2_l0l1[2]; // this will be unit-normalized if the lines are so
  }
  //        
  //    pChart: just unit rands 17x1
  //        sphere(7,1)|sphere(5,1)|sphere(5,1)
  //
  util::rand_sphere(params+27+12,7);
  util::rand_sphere(params+27+12+7,5);
  util::rand_sphere(params+27+12+7+5,5);
//  F c1[7] = {0.356520517738511 ,  0.450534892837314 ,  0.497658671520414 ,  0.530494023592847 ,0.350361054584548 ,  0.040309061260114 ,  0.128240708712460};
//  memcpy(params+27+12,c1,7*sizeof(F));

//  F c2[5] =  { 0.608716490477115 ,  0.290014694962129 ,  0.462945690541627 ,  0.548557032724055 ,0.173557426764642};
//  memcpy(params+27+12+7,c2,5*sizeof(F));
//  memcpy(params+27+12+7+5,c2,5*sizeof(F));
}

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
minus_io<chicago14a, F>::
gammify(C<F> * __restrict params /*[ chicago: M::nparams]*/)
{
  typedef minus_util<F> util;
  //  params = (diag0|diag1|diag2|diag3|diag4).*params;
  // diag0 --> pF in params ----------------------------------------------------
  C<F> (*p)[3] = (C<F> (*)[3]) params;
  C<F> gammas[9]; 
  for (unsigned l=0; l < 9; ++l) {
    util::randc(gammas+l);
    const C<F> &g = gammas[l];
    p[l][0] *= g; p[l][1] *= g; p[l][2] *= g;
  }
  
  // ids of two point-point lines at tangents
  static unsigned constexpr triple_intersect[6][2] = {{0,3},{0+1,3+1},{0+2,3+2},{0,6},{0+1,6+1},{0+2,6+2}};

  // diag1 --> pTriple in params -----------------------------------------------
  unsigned i = 9*3; // TODO: move to unsigned char
  for (unsigned tl=0; tl < 6; ++tl) {  // for each tangent line
    params[i++] *= std::conj(gammas[triple_intersect[tl][0]]);
    params[i++] *= std::conj(gammas[triple_intersect[tl][1]]);
  }
  
  // pChart gammas -------------------------------------------------------------
  C<F> g;
  // diag2 -- tchart gamma
  util::randc(&g); for (unsigned k=0; k < 7; ++k) params[i++] *= g;
  // diag3 -- qchart, cam 2, gamma
  util::randc(&g); for (unsigned k=0; k < 5; ++k) params[i++] *= g;
  // diag4 -- qchart, cam 3, gamma
  util::randc(&g); for (unsigned k=0; k < 5; ++k) params[i++] *= g;
  //  p = (diag0|diag1|diag2|diag3|diag4).*p;
  //  total  27   12    7      5    5 = 56
  assert(i == 56);
}

// Generate "visible" line representation from input point-tangents
// 
// The points that have tangents are indicated by the indices i0  < i1 < 3
// 
// pLines is a 15x3 matrix of line coefs  (we use view-line-point index, this
// is inverted to match Hongyi)
//    1    -- l_1_1 --
//    2    -- l_1_2 --
//    3    -- l_1_3 --
//    4    -- l_2_1 --
//    5    -- l_2_2 --
//    6    -- l_2_3 --
//    7    -- l_3_1 --
//    8    -- l_3_2 --
//    9    -- l_3_3 --
//    10   -- l_4_1 --
//    11   -- l_4_2 --
//    12   -- l_4_3 --
//    13   -- l_5_1 --
//    14   -- l_5_2 --
//    15   -- l_5_3 --
//    
//    l_line_view
//    
//    These lines are:
//
//    l_1: Point 1&2  (A, B)
//    l_2: Point 1&3  (A, C)
//    l_3: Point 2&3  (B, C)
//    l_4: Tangent at Point 1 (A)
//    l_5: Tangent at Point 2 (B)
//
// NOTE: the input tangent vector will be used as scratch so make a copy
// if you intend to reuse it 
//
// Input points and tangents in normalized image coordinates.

template <typename F>
bool 
minus_io<chicago14a, F>::
point_tangents2lines(const F p[pp::nviews][pp::npoints][io::ncoords2d], const F t[pp::nviews][pp::npoints][io::ncoords2d], unsigned i0, unsigned i1, F plines[pp::nvislines][io::ncoords2d_h])
{
  typedef minus_3d<F> vec;
  typedef minus_array<M::nve,F> v;
  
  assert (i0 < i1 && i1 < 3);
  unsigned i2 = (i0 == 0) ? ((i1 == 1) ? 2 : 1) : 0;

  static constexpr double eps = 1e-4; // very important to tune this as it will
                                      // save a lot of time if trash is
                                      // early-detected
  if (v::area2(p[0][i0],p[0][i1],p[0][i2])  < eps ||  // retinal area.  Could be spherical area 
      v::area2(p[1][i0],p[1][i1],p[1][i2])  < eps || 
      v::area2(p[2][i0],p[2][i1],p[2][i2])  < eps) {
#ifndef NDEDBUG
          std::cerr << "MINUS: area error ------------------------\n";
          std::cerr << "Areas: " << 
            v::area2(p[0][i0],p[0][i1],p[0][i2]) << " "  << 
            v::area2(p[1][i0],p[1][i1],p[1][i2]) << " " << 
            v::area2(p[2][i0],p[2][i1],p[2][i2]) << std::endl;
#endif
    return false;
  }
  
  // TODO: require inpunt points be on the view-sphere (unit-normalized), not h-normalized
  vec::cross2(p[0][i0], p[0][i1], plines[0]);
  vec::cross2(p[1][i0], p[1][i1], plines[1]);
  vec::cross2(p[2][i0], p[2][i1], plines[2]);
  
  vec::cross2(p[0][i0], p[0][i2], plines[3]);
  vec::cross2(p[1][i0], p[1][i2], plines[4]);
  vec::cross2(p[2][i0], p[2][i2], plines[5]);
  
  vec::cross2(p[0][i1], p[0][i2], plines[6]);
  vec::cross2(p[1][i1], p[1][i2], plines[7]);
  vec::cross2(p[2][i1], p[2][i2], plines[8]);

  // tangent at point p[i0]
  minus_3d<F>::point_tangent2line(p[0][i0], t[0][i0], plines[9]);
  minus_3d<F>::point_tangent2line(p[1][i0], t[1][i0], plines[10]);
  minus_3d<F>::point_tangent2line(p[2][i0], t[2][i0], plines[11]);
 
  // tangent at point p[i1]
  minus_3d<F>::point_tangent2line(p[0][i1], t[0][i1], plines[12]);
  minus_3d<F>::point_tangent2line(p[1][i1], t[1][i1], plines[13]);
  minus_3d<F>::point_tangent2line(p[2][i1], t[2][i1], plines[14]);

  io::normalize_lines(plines, pp::nvislines);

  if (v::abs_angle_between_lines(plines[0], plines[9])  < eps || 
      v::abs_angle_between_lines(plines[1], plines[10]) < eps ||
      v::abs_angle_between_lines(plines[2], plines[11]) < eps ||

      v::abs_angle_between_lines(plines[3], plines[9])  < eps ||
      v::abs_angle_between_lines(plines[4], plines[10]) < eps ||
      v::abs_angle_between_lines(plines[5], plines[11]) < eps ||
      
      v::abs_angle_between_lines(plines[6], plines[12]) < eps ||
      v::abs_angle_between_lines(plines[7], plines[13]) < eps ||
      v::abs_angle_between_lines(plines[8], plines[14]) < eps ||
      
      v::abs_angle_between_lines(plines[0], plines[12]) < eps ||
      v::abs_angle_between_lines(plines[1], plines[13]) < eps ||
      v::abs_angle_between_lines(plines[2], plines[14]) < eps) {
#ifndef NDEDBUG
    std::cerr << "MINUS: angle error ------------------------\n";
    std::cerr << "Angles: " << 
    v::abs_angle_between_lines(plines[0], plines[9])  << " " << 
    v::abs_angle_between_lines(plines[1], plines[10]) << " " <<
    v::abs_angle_between_lines(plines[2], plines[11]) << " " <<

    v::abs_angle_between_lines(plines[3], plines[9]) << " " <<
    v::abs_angle_between_lines(plines[4], plines[10]) << " " <<
    v::abs_angle_between_lines(plines[5], plines[11]) << " " <<
          
    v::abs_angle_between_lines(plines[6], plines[12]) << " " <<
    v::abs_angle_between_lines(plines[7], plines[13]) << " " <<
    v::abs_angle_between_lines(plines[8], plines[14]) << " " <<
          
    v::abs_angle_between_lines(plines[0], plines[12]) << " " <<
    v::abs_angle_between_lines(plines[1], plines[13]) << " " <<
    v::abs_angle_between_lines(plines[2], plines[14]);
#endif
    return false;
  }
  
  return true;
}

// gammify_start_params: set to false if your start parameters are already
// gammified. 
template <typename F>
inline void
minus_io<chicago14a, F>::
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

// \param[in] tgts: three tangents, one at each point, in normalized coordinates
// (inverted intrinsics).  Only two tangents will actually be used. If one of
// the points in each image has no reliable or well-defined tangents, you can
// pass anything (zeros or unallocated memory); it will be ignored. 
// only tgt[view][id_tgt0][:] and tgt[view][id_tgt1][:] will be used.
//
// id_tgt0  < id_tgt0 < 3
// 
template <typename F>
bool 
minus_io<chicago14a, F>::
point_tangents2params(
    const F p[pp::nviews][pp::npoints][io::ncoords2d], 
    const F tgt[pp::nviews][pp::npoints][io::ncoords2d], 
    unsigned id_tgt0, unsigned id_tgt1, 
    C<F> * __restrict params/*[static 2*M::nparams]*/, 
    bool gammify_start_params)
{
  // the user provides the start params in the first half of params.
  // we fill the second half and gammify both.
  F plines[pp::nvislines][io::ncoords2d_h];
  if (!point_tangents2lines(p, tgt, id_tgt0, id_tgt1, plines))
    return false;
  get_params_start_target(plines, params, gammify_start_params);
  return true;
}

// Same but for pixel input
template <typename F>
inline bool
minus_io<chicago14a, F>::
point_tangents2params_img(
    const F p[pp::nviews][pp::npoints][io::ncoords2d], 
    const F tgt[pp::nviews][pp::npoints][io::ncoords2d], 
    unsigned id_tgt0, unsigned id_tgt1, 
    const F K[/*3 or 2*/][io::ncoords2d_h], 
    C<F> * __restrict params/*[static 2*M::nparams]*/,
    bool gammify_start_params)
{
  F pn[pp::nviews][pp::npoints][io::ncoords2d];
  F tn[pp::nviews][pp::npoints][io::ncoords2d];
  
  // see if uno minus  default_gammas_m2 is less than 1
  io::invert_intrinsics(K, p[0], pn[0], pp::npoints);
  io::invert_intrinsics(K, p[1], pn[1], pp::npoints);
  io::invert_intrinsics(K, p[2], pn[2], pp::npoints);
  // don't use all three, but just invert all anyways.
  io::invert_intrinsics_tgt(K, tgt[0], tn[0], pp::npoints);
  io::invert_intrinsics_tgt(K, tgt[1], tn[1], pp::npoints);
  io::invert_intrinsics_tgt(K, tgt[2], tn[2], pp::npoints);
  return point_tangents2params(pn, tn, id_tgt0, id_tgt1, params/*[static 2*M::nparams]*/, gammify_start_params);
}

} // namespace minus

// Highlevel solver interface - Class minus ------------------------------------

#include <thread>
#include "chicago14a-default-data.h"

namespace MiNuS {


// Intrinsics already inverted 
// (inside RANSAC one will alredy have pre-inverted K)
//
// Input: points in pp:nviews views
// Input: tangents in pp:nviews views (e.g., SIFT orientations)
// Input: how to pick the tangent. For now, for Chicago we only consider
// the tangents on the first two points on each view.
// 
// Output: solutions_cams
// where the camera matrix P^t = [R|T]^t is cameras[sol_number][view_id][:][:]
// where view_id is 0 or 1 for second and third camera relative to the first,
// resp.
//
// Output: nsols, the number of solutions
// Output: id_sols
// a vector of the ids of the points that lead to each solution:
// So each solution is actually cameras[id_sols[i]][view_id][:][:], for i=1 to
// nsols.
//
// This design is for cache speed. Translation in the camera matrix is stored
// such that its coordinates are memory contiguous.
// 
// The cameras array is fixed in size to NSOLS which is the max
// number of solutions, which perfectly fits in memory. The caller must pass an
// array with that minimum.
// 
// returns false in case of numerical failure to find valid real solutions
// 
template <typename F>
inline bool
minus<chicago14a, F>::solve(
    const F p[pp::nviews][pp::npoints][io::ncoords2d], 
    const F tgt[pp::nviews][pp::npoints][io::ncoords2d], 
    F solutions_cams[M::nsols][pp::nviews-1][4][3],  // first camera is always [I | 0]
    unsigned id_sols[M::nsols],
    unsigned *nsols_final,
    unsigned nthreads
    )
{
  typedef minus_data<chicago14a,F> data;
  alignas(64) C<F> params[2*M::f::nparams];
  memcpy(params, data::params_start_target_, M::f::nparams*sizeof(C<F>));
  
  constexpr int id_tgt0 = 0; constexpr int id_tgt1 = 1; // TODO: select the best / least degenerate directions
  if (!io::point_tangents2params(p, tgt, id_tgt0, id_tgt1, params))
    return false;

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
 
  // decode solutions into 3x4 cams (actually 4x3 in mem)
  io::all_solutions2cams(solutions, solutions_cams, id_sols, nsols_final);

  // filter solutions that have no positive >1 depth for all three views
  return true;
}

// 
// same as solve() but intrinsics not inverted (input is in actual pixel units)
// returns false in case of numerical failure to find valid real solutions
// 
template <typename F>
inline bool
minus<chicago14a, F>::solve_img(
    const F K[/*3 or 2 ignoring last line*/][io::ncoords2d_h],
    const F p[pp::nviews][pp::npoints][io::ncoords2d], 
    const F tgt[pp::nviews][pp::npoints][io::ncoords2d], 
    F solutions_cams[M::nsols][pp::nviews-1][4][3],  // first camera is always [I | 0]
    unsigned id_sols[M::nsols],
    unsigned *nsols_final,
    unsigned nthreads)
{
  F pn[pp::nviews][pp::npoints][io::ncoords2d];
  F tn[pp::nviews][pp::npoints][io::ncoords2d];
  
  // see if uno minus  default_gammas_m2 is less than 1
  io::invert_intrinsics(K, p[0], pn[0], pp::npoints);
  io::invert_intrinsics(K, p[1], pn[1], pp::npoints);
  io::invert_intrinsics(K, p[2], pn[2], pp::npoints);
  // don't use all three, but just invert all anyways.
  io::invert_intrinsics_tgt(K, tgt[0], tn[0], pp::npoints);
  io::invert_intrinsics_tgt(K, tgt[1], tn[1], pp::npoints);
  io::invert_intrinsics_tgt(K, tgt[2], tn[2], pp::npoints);

  return solve(pn, tn, solutions_cams, id_sols, nsols_final, nthreads);
}

//
// Performs tests to see if there are potentially valid solutions,
// without making use of ground truth. 
// 
template <typename F>
inline bool 
minus_io<chicago14a, F>::
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
#endif // chicago14a_hxx_
