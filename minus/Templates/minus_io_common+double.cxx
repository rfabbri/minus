// Helps instantiate and reuse code,
// when including the .hxx directly causes inefficiencies and slowdown
// 
// \author Ricardo Fabbri
// \date Created: Fri Feb  8 17:42:49 EST 2019
#include <minus/minus.hxx>

namespace MiNuS {
  
template struct minus_io_common<double>;

}
