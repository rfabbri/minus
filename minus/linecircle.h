#ifndef linecircle_h_
#define linecircle_h_

template <>
struct formulation_parameters<linecircle> {
  static constexpr unsigned nsols = 2;    // number of solutions
  static constexpr unsigned nve = 2;      // size of the system (Number of Variables or Equations)
                                          // this formulation: two quaternions and 
                                          // two translation vectors
  static constexpr unsigned nparams = 6;  // number of parameters
};

// setup problem (data) parameters and attributes
template <>
struct problem_parameters<linecircle> {
  // static constexpr unsigned nviews = 3; 
  // static constexpr unsigned npoints = 3;
};

#endif // linecircle_h_
