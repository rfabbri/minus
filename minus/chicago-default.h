#ifndef chicago_default_h_
#define chicago_default_h_
// Include this last of all in user app

#include "minus.h"
#include "chicago14a-default-data.h" // default start system etc

// Convenience types to use Minus and the Chicago problem
typedef MiNuS::minus_core<MiNuS::chicago> M;
typedef MiNuS::minus_io<MiNuS::chicago> io;
typedef MiNuS::minus_io_14a<MiNuS::chicago> io14;
typedef MiNuS::minus_data<MiNuS::chicago, double> data;
extern template struct MiNuS::minus_data<MiNuS::chicago, double>;

#endif   // chicago_default_h_
