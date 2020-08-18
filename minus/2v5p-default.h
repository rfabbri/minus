#ifndef 2v5p_default_h_
#define 2v5p_default_h_
// Include this last of all in user app

#include "minus.h"
#include "2v5p9a-default-data.h" // default start system etc

// Convenience types to use Minus and the 2v5p problem
typedef MiNuS::minus_core<MiNuS::2v5p> M;
typedef MiNuS::minus_io<MiNuS::2v5p> io;
typedef MiNuS::minus_io_14a<MiNuS::2v5p> io14;
typedef MiNuS::minus_data<MiNuS::2v5p, double> data;

#endif   // 2v5p_default_h_
