#ifndef linecircle_default_h_
#define linecircle_default_h_
// Include this last of all in user app

#include "minus.h"
#include "linecircle2a-default-data.h" // default start system etc

// Convenience types to use Minus and the linecircle problem
typedef MiNuS::minus_core<MiNuS::linecircle> M;
typedef MiNuS::minus_io<MiNuS::linecircle> io;
typedef MiNuS::minus_io_14a<MiNuS::linecircle> io14;
typedef MiNuS::minus_data<MiNuS::linecircle, double> data;
extern template struct MiNuS::minus_data<MiNuS::linecircle, double>;

#endif   // linecircle_default_h_
