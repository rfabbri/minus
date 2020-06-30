#ifndef chicago_default_h_
#define chicago_default_h_

#include "minus.h"
#include "chicago14a-default-data.h" // default start system etc

// Convenience types to use Minus and the Chicago problem
typedef MiNuS::minus_core<chicago> M;
typedef MiNuS::minus_io<chicago> io;
typedef MiNuS::minus_io_14a<chicago> io14;
typedef MiNuS::minus_data<chicago, Float> data;

#endif   // chicago_default_h_
