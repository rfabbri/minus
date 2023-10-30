// builds a direct expression form an SLP
// an actual straight line program
// Each output gate Gi will be of the form
// Gi = F_i(input gates)
// Letting the compiler optimize the evaluations

// For each output gate Gi
//    Compute the tree back to the input gates
//    Evaluate it

// pseudocode:
//function de = direct_expression(Gi)
//  if (leaf or not left or right)
//    de = Gname
//  end
//  de = direct_expression(left) + mult_type(Gi) + direct_expression(right)
//endfunction

//function de = direct_expression(Gi)
//  if (leaf or not left or right)
//    de = Gname
//  end
//  de = direct_expression(left) + mult_type(Gi) + direct_expression(right)
//endfunction

for i=1:max_n_nodes
  if nodes_in_graph(i) == 1
    if isempty(dag_from(i))
      disp nodename(i)
    end
  end
end
