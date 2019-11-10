// Helps instantiate and reuse code,
// when including the .hxx directly causes inefficiencies and slowdown
// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date Created: Fri Feb  8 17:42:49 EST 2019
#include <minus.hxx>
template class minus_core<chicago14a, double>;
template class minus_io_shaping<chicago14a, double>;
template class minus_util<double>;
