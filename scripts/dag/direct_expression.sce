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


function ids = get_roots()
  ids = [];
  for i=1:max_n_nodes
    if nodes_in_graph(i) == 1
      if isempty(dag_to(i))
        ids($+1) = i;
      end
    end
  end
endfunction

function de = direct_expression(i)
  if(nodes_in_graph(i) == 0)
    de = ''
  else
    if isempty(dag_from(i)) 
      de = nodename(i);
    end
    de = '(' + direct_expression(dag_from(i)(1)) + ')' + node_type(i) + '(' + direct_expression(dag_from(i)(2)) + ')'
  end
endfunction

disp roots:
l = get_roots()

for i=1:size(l,'*')
  direct_expression(l(i))
end
