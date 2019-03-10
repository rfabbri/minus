// traverses the SLP gate graph constructed in dag.sce
// outputs vectorized form of the operations.


// 1- traverse graph and build all leaf nodes / without any incoming edges
//    insert these nodes in Q+ and Q- stacks

exec tsort_ini.sce;
// exec tsort_queue.sce;
