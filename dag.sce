clear;

global n_Cs;
global n_Xs;
global n_Gs; 
global n_total;
global max_n_nodes;
global dag_from;
global dag_to;
global nodes_in_graph;
global node_type;
global node_expr;
n_Cs = 2;
n_Xs = 3;
n_Gs = 4; 
//n_Xs = 127:
//n_Gs = 3439
n_total = n_Cs + n_Xs + n_Gs;
max_n_nodes = n_total;
dag_from = list();
dag_to = list();
nodes_in_graph = zeros(n_total,1);
node_type = list();
node_expr = list();
node_type(max_n_nodes) = [];
node_expr(max_n_nodes) = [];
dag_from(max_n_nodes) = list();
dag_to(max_n_nodes) = list();
for i=1:max_n_nodes
  dag_from(i) = list();
  dag_to(i) = list();
end

// use struct
//dag = list();
//dag(1) = dag_from;
//dag(2) = dag_to;

function i = get_id(variablestr)
  letternum = eval(part(variablestr,2:$));
  select part(variablestr,1),
  case 'C' then
    i = letternum + 1;
  case 'X' then
    i = n_Cs + letternum + 1;
  case 'G' then
    i = n_Cs + n_Xs + letternum + 1;
  end
endfunction

function [lhsvar, rhsvars, type] = get_vars(strexpr)
  // regexp(str, '/G[0-9][0-9]/')
  ts = tokens(strexpr);
  if (isempty(ts))
    lhsvar = []
    return
  end
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
  global dag_from;
  global dag_to;
  global nodes_in_graph;
  global node_type;
  global node_expr;
  
  [lhs, rhs, op] = get_vars(strexpr);
  node_left = get_id(lhs);
  node_type(node_left) = '+'
  node_expr(node_left) = strexpr;
  nodes_in_graph(node_left) = 1;

  for i=1:size(rhs,'*')
    node_i = get_id(rhs(i));
    dag_to(node_i) = node_left;
    dag_from(node_left)($+1) = node_i;
    nodes_in_graph(node_i) = 1;
  end
endfunction

function vn = nodename(id)
  if id <= n_Cs
    vn = 'C' + string(id-1);
  elseif id <= n_Xs+n_Cs
    vn = 'X' + string(id-1-n_Cs);
  else
    vn = 'G' + string(id-1-n_Xs-n_Cs);
  end
endfunction

function gviz = export_graphviz()
  gviz = 'digraph {';
  // go through each in dag_to 
  for i=1:max_n_nodes
    if nodes_in_graph(i) == 1
      if ~isempty(dag_to(i));
        vn = nodename(i);
        for k=1:size(dag_to(i),'*')
          gviz = [gviz; vn + ' -> ' + nodename(dag_to(i)(k))];
        end
      end
    end
  end
  gviz = [gviz; '}'];
  write('chicago.dot',gv)
endfunction 

//function print_vectorized_topo_order()
//endfunction


//txt=mgetl('chicago.bare');

txt = [
  'G2 = G1 * X0'
  'G3 = X2 * X1'];

for i=1:size(txt,1)
  if isempty(tokens(txt(i)))
    continue
  end
  slpdag_add(txt(i))
end


//print_vectorized_topo_order();
