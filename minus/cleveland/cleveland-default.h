#ifndef cleveland_default_h_
#define cleveland_default_h_
// Include this last of all in user app

#include "minus.h"
#include "cleveland14a-default-data.h" // default start system etc

// Convenience types to use Minus and the Chicago problem
typedef MiNuS::minus_core<MiNuS::cleveland> M;
typedef MiNuS::minus_io<MiNuS::cleveland> io;
typedef MiNuS::minus_io_14a<MiNuS::cleveland> io14;
typedef MiNuS::minus_data<MiNuS::cleveland, double> data;

#endif   // cleveland_default_h_
