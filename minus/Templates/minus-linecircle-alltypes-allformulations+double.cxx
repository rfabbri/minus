// Helps instantiate and reuse code,
// when including the .hxx directly causes inefficiencies and slowdown
// 
// \author Ricardo Fabbri
// \date Created: Fri Feb  8 17:42:49 EST 2019
#include <minus/minus.hxx>
#include <minus/linecircle2a-default-data.hxx>

namespace MiNuS {
  
template class minus_core<linecircle2a, double>;
template struct minus_io<linecircle2a, double>;
template struct minus<linecircle2a, double>;
template struct minus_data<linecircle2a, double>;
  
}
