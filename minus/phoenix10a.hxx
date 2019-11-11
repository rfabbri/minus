#ifndef phoenix10a_hxx_
#define phoenix10a_hxx_

template <typename F>
struct eval<phoenix10a, F> {
  static void Hxt(const C<F> * __restrict__ x /*x, t*/,    const C<F> * __restrict__ params, C<F> * __restrict__ y /*HxH*/);
  static void HxH(const C<F> * __restrict__ x /*x and t*/, const C<F> * __restrict__ params, C<F> * __restrict__ y /*HxH*/);
};

// Evaluates Hx and Ht at the same time, reusing expressions.
// 
// Map from a multivariate poly with x 127-dimensional to y NVExNVEPLUS1 dimensional
// Where 127 = 14 for x, 1 for t, 2*56 total parameters. Returns y = [Hx|Ht]
// 
// cCode(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"Ht",gateMatrix{cameraVars})
// (Ask Tim for the way to use cCode so that the input orders are like this.
template <typename F>
inline void 
eval<phoenix10a, F>::
Hxt(const C<F> * __restrict__ x /*x, t*/, const C<F> * __restrict__ params, C<F> * __restrict__ y /*HxH*/) 
{
}

// Evaluates Hx and H at the same time, reusing expressions.
// 
// Map from a multivariate poly with x 127-dimensional to y NVExNVEPLUS1 dimensional
// Where 127 = 14 for x, 1 for t, 2*56 total parameters. Returns where y = [Hx|H]
// 
// cCode(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"H",gateMatrix{cameraVars})
// (Ask Tim for the way to use cCode so that the input orders are like this.
template <typename F>
inline void 
eval<phoenix10a, F>::
HxH(const C<F>* __restrict__ x /*x and t*/, const C<F> * __restrict__ params, C<F>* __restrict__ y /*HxH*/) 
{
}

//------------------------------------------------------------------------------

// TODO: what parameters to put here: NPOINTS = 2, NTANGENTS = 5:
// NTANGENTS1 = 2, NTANGENTS2=3 --> how to define these?
// What really matters is to specialize the P parameter.
// The rest we deal with as a consequence.
//
// 
// A base template class for zero and first-order problems (points,
// lines, tangents or lines at points).
template <typename F>
struct minus_io_shaping<3/*NVIEWS*/, 3/*NPOINTS*/, 0/*NFREELINES*/, 2/*NTANGENTS*/, 312/*NSOLS*/, 14/*NVE*/, 56/*NPARAMS*/,  phoenix10a, F> {
  typedef minus_core<312, 14, 56, phoenix10a, F> M;
  typedef struct M::solution solution;
  static constexpr unsigned nviews = 3; 
  static constexpr unsigned npoints = 3;
  static constexpr unsigned nfreelines = 0;
  // even though Chicago needs only 2 tangents, api assumes 3 tangents are given,
  // out of which two are selected by indexing. This is the most common use
  // case, where all features naturally have tangents. If strictly 2 tangents
  // are to be passed, you can leave the unused one as zeros throughout the API.
  static constexpr unsigned ntangents = 2;
  static constexpr unsigned ncoords = 2;  // just a documented name for the number of inhomog coordinates
  static constexpr unsigned ncoords_h = 3;  // just a name for the usual number of homog coordinates in P^2
  static constexpr unsigned ncoords3d = 3;  // just a documented name for the number of inhomog 3D coordinates
  // number of lines connecting each pair of points plus going through points
  // plus the number of free lines in the first order problem.
  // Note that tangent orientation may help ruling out solutions; this is why
  // we call it tangent, and not general lines at points. There is more
  // information at tangents which can be used as part of the model for curves.
  // The tangent orientation can be constrained by running without orientation
  // mattering at first, and then propagating these to neighboring features
  // along curves.
  static constexpr unsigned nvislines = ( (npoints*(npoints-1) >> 1) + ntangents + nfreelines ) * nviews; 
  // nvislines = 15 for Chicago.
  // INPUT ---------------------------------------------------------------------
  static void point_tangents2params(F p[nviews][npoints][ncoords], F tgt[nviews][npoints][ncoords], unsigned id_tgt0, unsigned id_tgt1, C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  // this function is the same for all problems
  static void get_params_start_target(F plines[/*15 for chicago*/][ncoords_h], C<F> * __restrict__ params/*[static 2*M::nparams]*/);
  static void gammify(C<F> * __restrict__ params/*[ chicago: M::nparams]*/);
  static void point_tangents2lines(F p[nviews][npoints][ncoords], F tgt[nviews][npoints][ncoords], unsigned id_tgt0, unsigned id_tgt1, F plines[nvislines][ncoords_h]);
  static void lines2params(F plines[nvislines][ncoords_h], C<F> * __restrict__ params/*[static M::nparams]*/);

  // OUTPUT --------------------------------------------------------------------
  static void all_solutions2cams(solution raw_solutions[M::nsols], F cameras[M::nsols][2][4][3], unsigned id_sols[M::nsols], unsigned *nsols_final);
  static void solution2cams(F rs[M::nve], F cameras[2][4][3]);
};

// we only use the first half of the outer
// 2*M::nparams array 
// after this fn, complex part zero, but we will use this space later
// to gammify/randomize
template <typename F>
inline void 
minus_io_shaping<3/*NVIEWS*/, 3/*NPOINTS*/, 0/*NFREELINES*/, 2/*NTANGENTS*/, 312/*NSOLS*/, 14/*NVE*/, 56/*NPARAMS*/,  phoenix10a, F>::
lines2params(F plines[nvislines][ncoords_h], C<F> * __restrict__ params/*[static 2*M::nparams]*/)
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
  unsigned triple_intersections[6][3] = 
    {{0,3,9},{0+1,3+1,9+1},{0+2,3+2,9+2},{0,6,12},{0+1,6+1,12+1},{0+2,6+2,12+2}};

  C<F> (*params_lines)[2] = (C<F> (*)[2]) (params+27);
  // express each of the 6 tangents in the basis of the other pairwise lines
  // intersecting at the same point
  for (unsigned l=0; l < 6; ++l) {
    const F *l0 = plines[triple_intersections[l][0]];
    const F *l1 = plines[triple_intersections[l][1]];
    const F *l2 = plines[triple_intersections[l][2]];
    double l0l0 = vec::dot(l0,l0), l0l1 = vec::dot(l0,l1), l1l1 = vec::dot(l1,l1),
    l2l0 = vec::dot(l2,l0), l2l1 = vec::dot(l2,l1);
    // cross([l0l0 l1l0 l2l0], [l0l1 l1l1 l2l1], l2_l0l1);
    double l2_l0l1[3]; 
    {
      F v1[3], v2[3];
      v1[0] = l0l0; v1[1] = l0l1; v1[2] = l2l0;
      v2[0] = l0l1; v2[1] = l1l1; v2[2] = l2l1;
      vec::cross(v1, v2, l2_l0l1);
    }
    params_lines[l][0] = l2_l0l1[0]/l2_l0l1[2]; // divide by the last coord (see cross prod formula, plug direct)
    params_lines[l][1] = l2_l0l1[1]/l2_l0l1[2];
  }
  //        
  //    pChart: just unit rands 17x1
  //        sphere(7,1)|sphere(5,1)|sphere(5,1)
  //
  util::rand_sphere(params+27+12,7);
  util::rand_sphere(params+27+12+7,5);
  util::rand_sphere(params+27+12+7+5,5);
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
minus_io_shaping<3/*NVIEWS*/, 3/*NPOINTS*/, 0/*NFREELINES*/, 2/*NTANGENTS*/, 312/*NSOLS*/, 14/*NVE*/, 56/*NPARAMS*/,  phoenix10a, F>::
gammify(C<F> * __restrict__ params /*[ chicago: M::nparams]*/)
{
  typedef minus_util<F> util;
  //  params = (diag0|diag1|diag2|diag3|diag4).*params;
  // diag0
  C<F> (*p)[3] = (C<F> (*)[3]) params;
  C<F> gammas[9]; 
  for (unsigned l=0; l < 9; ++l) {
    util::randc(gammas+l);
    const C<F> &g = gammas[l];
    p[l][0] *= g; p[l][1] *= g; p[l][2] *= g;
  }
  
  // ids of two point-point lines at tangents
  unsigned triple_intersect[6][2] = {{0,3},{0+1,3+1},{0+2,3+2},{0,6},{0+1,6+1},{0+2,6+2}};

  // diag1
  unsigned i = 9*3;
  for (unsigned tl=0; tl < 6; ++tl) {  // for each tangent line
    params[i++] *= std::conj(gammas[triple_intersect[tl][0]]);
    params[i++] *= std::conj(gammas[triple_intersect[tl][1]]);
  }
  
  C<F> g;
  // diag2 -- tchart gamma
  util::randc(&g); for (unsigned k=0; k < 7; ++k) params[i++] *= g;
  // diag3 -- qchart, cam 2, gamma
  util::randc(&g); for (unsigned k=0; k < 5; ++k) params[i++] *= g;
  // diag4 -- qchart, cam 3, gammg
  util::randc(&g); for (unsigned k=0; k < 5; ++k) params[i++] *= g;
  //  p = (diag0|diag1|diag2|diag3|diag4).*p;
  //  total  27   12    7      5    5 = 56
  assert(i == 56);
}


// Generate "visible" line representation from input point-tangents
// 
// The points that have tangents are indicated by the indices id_tgt0  < id_tgt0 < 3
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
template <typename F>
inline void 
minus_io_shaping<3/*NVIEWS*/, 3/*NPOINTS*/, 0/*NFREELINES*/, 2/*NTANGENTS*/, 312/*NSOLS*/, 14/*NVE*/, 56/*NPARAMS*/,  phoenix10a, F>::
point_tangents2lines(F p[nviews][npoints][ncoords], F t[nviews][npoints][ncoords], unsigned i0, unsigned i1, F plines[nvislines][ncoords_h])
{
  typedef minus_3d<F> vec;
  
  assert (i0 < i1 && i1 < 3);
  unsigned i2 = (i0 == 0) ? ((i1 == 1) ? 2 : 1) : 0;
  
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
  // TODO: test normalize to unit vectors for better numerics
}

// \param[in] tgts: three tangents, one at each point.
// Only two tangents will actually be used. If one of the points
// in each image has no reliable or well-defined tangents,
// you can pass anything (zeros or unallocated memory); 
// it will be ignored. 
// only tgt[view][id_tgt0][:] and tgt[view][id_tgt1][:] will be used.
//
// id_tgt0  < id_tgt0 < 3
// 
template <typename F>
inline void 
minus_io_shaping<3/*NVIEWS*/, 3/*NPOINTS*/, 0/*NFREELINES*/, 2/*NTANGENTS*/, 312/*NSOLS*/, 14/*NVE*/, 56/*NPARAMS*/,  phoenix10a, F>::
point_tangents2params(F p[nviews][npoints][ncoords], F tgt[nviews][npoints][ncoords], unsigned id_tgt0, unsigned id_tgt1, C<F> * __restrict__ params/*[static 2*M::nparams]*/)
{
  F plines[nvislines][ncoords_h];
  point_tangents2lines(p, tgt, id_tgt0, id_tgt1, plines);
  lines2params(plines, params);
  gammify(params);
k gammify(params+M::nparams);
}

template <typename F>
inline void
minus_io_shaping<3/*NVIEWS*/, 3/*NPOINTS*/, 0/*NFREELINES*/, 2/*NTANGENTS*/, 312/*NSOLS*/, 14/*NVE*/, 56/*NPARAMS*/,  phoenix10a, F>::
get_params_start_target(F plines[/*15 for chicago*/][ncoords_h], C<F> * __restrict__ params/*[static 2*M::nparams]*/)
{
  lines2params(plines, params);
  gammify(params);
  gammify(params+M::nparams);
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
template <typename F>
inline void 
minus_io_shaping<3/*NVIEWS*/, 3/*NPOINTS*/, 0/*NFREELINES*/, 2/*NTANGENTS*/, 312/*NSOLS*/, 14/*NVE*/, 56/*NPARAMS*/,  phoenix10a, F>::
all_solutions2cams(solution raw_solutions[M::nsols], F cameras[M::nsols][2][4][3], 
               unsigned id_sols[M::nsols], unsigned *nsols_final)
{
  *nsols_final = 0;
  for (unsigned sol=0; sol < M::nsols; ++sol) {
    F real_solutions[M::nve];
    if (get_real(raw_solutions[sol], real_solutions)) {
      id_sols[(*nsols_final)++] = sol;
      // build cams by using quat2rotm
      solution2cams(real_solutions, (F (*)[4][3] ) (cameras + sol));
    }
  }
}

template <typename F>
inline void 
minus_io_shaping<3/*NVIEWS*/, 3/*NPOINTS*/, 0/*NFREELINES*/, 2/*NTANGENTS*/, 312/*NSOLS*/, 14/*NVE*/, 56/*NPARAMS*/,  phoenix10a, F>::
solution2cams(F rs[M::nve], F cameras[2/*2nd and 3rd cams relative to 1st*/][4][3])
{
  // camera 0 (2nd camera relative to 1st)
  quat2rotm(rs, (F [3][3]) cameras[0]);
  cameras[0][3][0] = rs[8];
  cameras[0][3][1] = rs[9];
  cameras[0][3][2] = rs[10];
  
  // camera 1 (3rd camera relative to 1st)
  quat2rotm(rs, (F [3][3]) cameras[1]);
  cameras[1][3][0] = rs[11];
  cameras[1][3][1] = rs[12];
  cameras[1][3][2] = rs[13];

  // quat12 rs(0:3), quat12 rs(4:7)
  //  T12 = solutions(9:11);
  //  T13 = solutions(12:14);
  //  R12 = quat2rotm(transpose(quat12));
  //  R13 = quat2rotm(transpose(quat13));
}
#endif // phoenix10a_hxx
