// Helps instantiate and reuse code,
// when including the .hxx directly causes inefficiencies and slowdown
// 
// \author Ricardo Fabbri
// \date Created: Fri Feb  8 17:42:49 EST 2019
#include <minus/minus.hxx>
#include <minus/2v5p9a-default-data.hxx>

namespace MiNuS {
  
template class minus_core<2v5p9a, double>;
template class minus_io_14a<2v5p9a, double>;
template class minus_io<2v5p9a, double>;
template class minus<2v5p9a, double>;
template class minus_data<2v5p9a, double>;
  
}
