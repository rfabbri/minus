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
  global dag_from;
  global nodes_in_graph;
  if(nodes_in_graph(i) == 0)
    de = ''
  else
    if isempty(dag_from(i)) 
      de = nodename(i);
    else
      p1a = '', p1b = ''; p2a = '', p2b = ''; p3a = '', p3b = '';
      if node_type(i) == '*'
        if ~isempty(dag_from(dag_from(i)(1)))
          if (node_type(dag_from(i)(1)) == '+')
            p1a = '(', p1b = ')';
          end
        end
        if ~isempty(dag_from(dag_from(i)(2)))
          if (node_type(dag_from(i)(2)) == '+')
            p2a = '(', p2b = ')';
          end
        end
      end
      de = p1a + direct_expression(dag_from(i)(1)) + p1b + node_type(i) ..
         + p2a + direct_expression(dag_from(i)(2)) + p2b;
      // At least two nodes forma gate. The rest mst be the same type
      // of addition otherwise it is not supported.
      // Since we assume it is addition, we do not parenthesize
      for k=3:length(dag_from(i))
        // for now only '+' will be the case here
        de = de + node_type(i) + p3a + direct_expression(dag_from(i)(k)) + p3b 
      end
    end
  end
endfunction

global dag_to_copy;
// remove node from dag_to(i)
function remove_from_dag_to(i,node)
  global dag_to_copy;
  
  for k=1:length(dag_to_copy(i))
    if dag_to_copy(i)(k) == node
      dag_to_copy(i)(k) = null()
      break;
    end
  end
end

function de = direct_expression_remove_out(i)
  global dag_from;
  global dag_to;
  global nodes_in_graph;
  global dag_to_copy;
  if nodes_in_graph(i) == 0 | isempty(dag_to_copy(i)) & ~isempty(dag_to(i))
    de = ''
  else
    if isempty(dag_from(i))
      de = nodename(i);
    else
      p1a = '', p1b = ''; p2a = '', p2b = ''; p3a = '', p3b = '';
      if node_type(i) == '*'
        if ~isempty(dag_from(dag_from(i)(1)))
          if (node_type(dag_from(i)(1)) == '+')
            p1a = '(', p1b = ')';
          end
        end
        if ~isempty(dag_from(dag_from(i)(2)))
          if (node_type(dag_from(i)(2)) == '+')
            p2a = '(', p2b = ')';
          end
        end
      end
      de = p1a + direct_expression(dag_from(i)(1)) + p1b + node_type(i) ..
         + p2a + direct_expression(dag_from(i)(2)) + p2b;

      remove_from_dag_to(dag_from(i)(1),i);
      remove_from_dag_to(dag_from(i)(2),i);

      // At least two nodes forma gate. The rest mst be the same type
      // of addition otherwise it is not supported.
      // Since we assume it is addition, we do not parenthesize
      for k=3:length(dag_from(i))
        // for now only '+' will be the case here
        de = de + node_type(i) + p3a + direct_expression(dag_from(i)(k)) + p3b 
        remove_from_dag_to(dag_from(i)(k),i);
      end
    end
  end
endfunction


// gets explicit expression for node i and below,
// without outputing expressions if they were alredy used
function p = pedag(i)
  dag_to_copy = dag_to;
  p = ''; p(i) = '';
  for k=i:-1:1
     p(k) = direct_expression_remove_out(k);
  end
endfunction

disp roots:
l = get_roots()
nr = size(l,'*');

de = '';
de(max_n_nodes) = '';
for i=1:size(l,'*')
  de(l(i)) = nodename(l(i)) + ' = ' + direct_expression(l(i));
end

del = de(l);

// p = pedag()

// write('expr.c',del);

//  "(C0+(C1*X14))*X15"
