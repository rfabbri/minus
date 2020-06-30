#ifndef cleveland_default_h_
#define cleveland_default_h_
// see .hxx for documentation

#include "minus.h"
#include "cleveland14a-default-data.h" // default start system etc

// Convenience types to use Minus and the Chicago problem
typedef MiNuS::minus_core<cleveland> M;
typedef MiNuS::minus_io<cleveland> io;
typedef MiNuS::minus_io_14a<cleveland> io14;
typedef MiNuS::minus_data<cleveland, Float> data;

#endif   // cleveland_default_h_
