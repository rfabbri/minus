#ifndef linecircle_h_
#define linecircle_h_

template <>
struct formulation_parameters<linecircle2a> {
  static constexpr unsigned nsols = 2;    // number of solutions
  static constexpr unsigned nve = 2;      // size of the system (Number of Variables or Equations)
  static constexpr unsigned nparams = 6;  // number of parameters
  
  struct settings;
  static const settings DEFAULT;
};

// Settings specific to each formulation such as prefilters to reject bad input,
// or postfilters to reject bad solutions
struct formulation_parameters<linecircle2a>::settings {
  settings():
  // Prefilter Parameters ------------------------------------------------------
  // These parameters are used in precomputations, such as early detection of
  // degeneracy before homotopy continuation
  // 
  // Very important to tune this as it will save a lot of time if trash is
  // early-detected
  // 
  // TODO: currently not implemented, but other formulations such as chicago14a
  // provide this
  prefilter_degeneracy_(true)
  // Perhaps filter circles that are degenerating,
  // or normalize line equtions to unit, etc.
  // 
  // Postfilter Parameters ------------------------------------------------------
  // 
  // TODO Epsilon for converting complex to real solutions etc.
  { }
  bool prefilter_degeneracy_; // prefilter degeneracy?
  // print() as a separate function in minus/cmd/cmd-util.h
};

// setup problem (data) parameters and attributes
template <>
struct problem_parameters<linecircle2a> {
  // static constexpr unsigned nviews = 3; 
  // static constexpr unsigned npoints = 3;
};

#endif // linecircle_h_
