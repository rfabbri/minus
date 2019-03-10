// traverses the SLP gate graph constructed in dag.sce
// outputs vectorized form of the operations.


// 1- traverse graph and build all leaf nodes / without any incoming edges
//    insert these nodes in Q+ and Q- stacks


Qplus = list();
Qtimes = list();
for i=1:max_n_nodes
  if nodes_in_graph(i) == 1 & isempty(dag_from(i))
    select node_type(i)
    case '+' then
      Qplus($+1) = i;
    case '*' then
      Qtimes($+1) = i;
    else 
      Qplus($+1) = i;
    end
  end
end


//Q = Qplus
//while isempty(Q)
//end

