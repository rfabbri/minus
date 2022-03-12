// Helps instantiate and reuse code,
// when including the .hxx directly causes inefficiencies and slowdown
// 
// \author Ricardo Fabbri
// \date Created: Fri Feb  8 17:42:49 EST 2019
#include <minus/minus.hxx>
#include <minus/chicago14a-default-data.hxx>

namespace MiNuS {
  
template class minus_core<chicago14a, double>;
template struct minus_io_14a<chicago14a, double>;
template struct minus_io<chicago14a, double>;
template struct minus<chicago14a, double>;
template struct minus_data<chicago14a, double>;
  
}
