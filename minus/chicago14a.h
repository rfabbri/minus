#ifndef chicago14a_h_
#define chicago14a_h_

template <>
struct formulation_parameters<chicago14a> {
  static constexpr unsigned nsols = 312;   // number of solutions
  static constexpr unsigned nve = 14;      // size of the system (Number of Variables or Equations)
                                           // this formulation: two quaternions and 
                                           // two translation vectors
  static constexpr unsigned nparams = 56;  // number of parameters
  //    params  is pF||pTriple||pChart  //  Hongyi: [pF; tripleChart; XR'; XT1'; XT2'];
  //    size       27     12     17 = 56
  //    pF: lines between points: 3lines*3views*3coordinates = 27
  //    pTriple: tangent lines each with 2 numbers: 2tangents*3views*2coords = 12
  //    pChart: 2xquaternions + 2translations + 1 homg coord per quat + 1 homg
  //    coord per [t01 t02] pair = 2*4+2*3+2+1 = 17
  //    internal note: see inGates in chicago.m2
};

template <>
struct problem_parameters<chicago14a> {
  static constexpr unsigned nviews = 3; 
  static constexpr unsigned npoints = 3;
  static constexpr unsigned nfreelines = 0;
  // even though Chicago needs only 2 tangents, api assumes 3 tangents are given,
  // out of which two are selected by indexing. This is the most common use
  // case, where all features naturally have tangents. If strictly 2 tangents
  // are to be passed, you can leave the unused one as zeros throughout the API.
  static constexpr unsigned ntangents = 2;
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
  // unsigned NVIEWS, unsigned NPOINTS /* per view*/, unsigned NFREELINES, unsigned NTANGENTS, 
};

#endif // chicago14a_h_
