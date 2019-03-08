n_Cs = 2;
n_Xs = 127:
n_Gs = 3439
n_total + n_Cs + n_Xs + n_Gs;
max_n_nodes = n_total;
dag_from = list();
dag_to = list();
nodes_in_graph = zeros(n_total,1);
node_type = list();
node_expr = list();

node_type(max_n_nodes) = [];
node_expr(max_n_nodes) = [];
dag_from(max_n_nodes) = [];
dag_to(max_n_nodes) = [];

// use struct
dag = list();
dag(1) = dag_from;
dag(2) = dag_to;

function i = slpdag_id(variablestr)
  select letter,
  case C then
    i = letternum + 1;
  case X then
    i = n_Cs + letternum + 1;
  case G then
    i = n_Cs + n_Xs + letternum + 1;
  end
endfunction

function [lhsvar, rhsvars, type] = get_vars(strexpr)
  // regexp(str, '/G[0-9][0-9]/')
  ts = tokens(strexpr);
  lhsvar = ts(1);
  if ts(2) ~= '='
    error('second token must be equal')
  end
  type = ts(4);
  if type ~= '+' & type ~= '*'
    error('fourth token must be an operand')
  end
  rhsvars = [];
  for i=3:2:size(ts,'*')
    rhsvars($+1) = ts(i);
  end
endfunction

// spdag_add("G50 = G39 + G40)
function slpdag_add(strexpr)
  get_vars(strexp);
  
  dag_to(39) = 50;
  dag_to(40) = 50;
  dag_from(50) = 39;
  dag_from(50) = 40;
  node_type(50) = '+'
  node_expr(50) = strexpr;
  nodes_in_graph(39) = 1;
  nodes_in_graph(40) = 1;
  nodes_in_graph(50) = 1;
endfunction
