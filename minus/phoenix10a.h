#ifndef phoenix10a_h_
#define phoenix10a_h_

template <>
struct formulation_parameters<phoenix10a> {
  static constexpr unsigned nsols = 480;   // number of solutions
  static constexpr unsigned nve = 10;      // size of the system (Number of Variables or Equations)
  static constexpr unsigned nparams = 42;  // number of parameters
};

template <>
struct problem_parameters<phoenix10a> {
  static constexpr unsigned nviews = 3;
  static constexpr unsigned npoints = 2;
  static constexpr unsigned nfreelines = 0;
  // even though Chicago needs only 2 tangents, api assumes 3 tangents are given,
  // out of which two are selected by indexing. This is the most common use
  // case, where all features naturally have tangents. If strictly 2 tangents
  // are to be passed, you can leave the unused one as zeros throughout the API.
  static constexpr unsigned ntangents[2] = {2, 3};
  // number of lines connecting each pair of points plus going through points
  // plus the number of free lines in the first order problem.
  // Note that tangent orientation may help ruling out solutions; this is why
  // we call it tangent, and not general lines at points. There is more
  // information at tangents which can be used as part of the model for curves.
  // The tangent orientation can be constrained by running without orientation
  // mattering at first, and then propagating these to neighboring features
  // along curves.
  static constexpr unsigned nvislines = ( (npoints*(npoints-1) >> 1) + ntangents[0] + ntangents[1] + nfreelines ) * nviews; 
  // nvislines = 15 for Chicago.
  // unsigned NVIEWS, unsigned NPOINTS /* per view*/, unsigned NFREELINES, unsigned NTANGENTS, 
};

#endif // phoenix10a_h_
