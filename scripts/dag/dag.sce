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
n_Cs = 3;
//n_Xs = 3;
//n_Gs = 4; 
//n_Xs = 8;
//n_Gs = 47;
n_Xs = 127;
n_Gs = 3440;
n_total = n_Cs + n_Xs + n_Gs;
max_n_nodes = n_total;
dag_from = list();
dag_to = list();
nodes_in_graph = zeros(n_total,1);
node_type = '';
node_expr = '';
node_type(max_n_nodes) = '';
node_expr(max_n_nodes) = '';
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
  node_type(node_left) = op;
  node_expr(node_left) = strexpr;
  nodes_in_graph(node_left) = 1;

  for i=1:size(rhs,'*')
    node_i = get_id(rhs(i));
    dag_to(node_i)($+1) = node_left;
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

function vn = nodenames(ids)
  vn = []
  for i=1:max(size(ids))
    if nodes_in_graph(ids(i))
      vn = [vn; nodename(ids(i))];
    else
      vn = [vn; string(ids(i))+ ' not in graph']
    end
  end
endfunction

function gviz = export_graphviz()
    
  // go through each in dag_to 
  if sum(nodes_in_graph) < 40
    gviz = ['digraph {'
      'splines=""line""'
      'ranksep=""1""'
      'node [shape=circle,color=dimgray,height=""0.3""]'];
    nt_times = '{ node [shape=diamond, height=""0.3"", width=""0.3"", color=orangered, fontcolor=orangered] ';
  else
    gviz = ['digraph {'
      'splines=""line""'
      'ranksep=""1""'
      'node [shape=circle,fontcolor=white,color=dimgray,label="""",height=""0.3""]'];
    nt_times = '{ node [shape=diamond, height=""0.3"", width=""0.3"", style=filled, color=orangered, fontcolor=orangered, fillcolor=orangered,label=""""] ';
  end
  for i=1:max_n_nodes
    if nodes_in_graph(i) == 1
      if node_type(i) == '*'
        nt_times = nt_times + ' ' + nodename(i);
      end
    end
  end
  nt_times = [nt_times; '}'];
  
  if sum(nodes_in_graph) < 40
    nt_plus = '{ node [shape=square, height=""0.3"", width=""0.3"", color=steelblue4, fontcolor=steelblue4] ';
  else
    nt_plus = '{ node [shape=square, height=""0.3"", width=""0.3"", style=filled, color=steelblue4, fontcolor=steelblue4, fillcolor=steelblue4,label=""""] ';
  end
  
  for i=1:max_n_nodes
    if nodes_in_graph(i) == 1
      if node_type(i) == '+'
        nt_plus = nt_plus + ' ' + nodename(i);
      end
    end
  end
  nt_plus= [nt_plus; '}'];
  gviz = [gviz; nt_times];    
  gviz = [gviz; nt_plus];    
  for i=1:max_n_nodes
    if nodes_in_graph(i) == 1
      gviz = [gviz; nodename(i)+'[label = ""' + string(i) + nodename(i) + '""];'];
    end
  end
  for i=1:max_n_nodes
    if nodes_in_graph(i) == 1
      if ~isempty(dag_to(i));
        vn = nodename(i);
        for k=1:size(dag_to(i))
          nt = node_type(dag_to(i)(k));
          if nt == '*'
            colorstr = 'color=orangered, fontcolor=orangered'
          else 
            colorstr = 'color=steelblue4, fontcolor=steelblue4'
          end
//          gviz = [gviz; vn + ' -> ' + nodename(dag_to(i)(k)) + '[label=""' + nt + '"", fontsize=""24.0""' + colorstr + '];'];
          gviz = [gviz; vn + ' -> ' + nodename(dag_to(i)(k)) + '[' + colorstr + '];'];
        end
      end
    end
  end
  gviz = [gviz; '}'];
  unix('rm -f chicago-tmp.dot');
  write('chicago-tmp.dot',gviz);
endfunction 

function gviz = export_graphviz_with_ranks()
    
  // go through each in dag_to 
  if sum(nodes_in_graph) < 40
    gviz = ['digraph {'
      'splines=""line""'
      'ranksep=""1""'
      'node [shape=circle,color=dimgray,height=""0.3""]'];
    nt_times = '{ node [shape=diamond, height=""0.3"", width=""0.3"", color=orangered, fontcolor=orangered] ';
  else
    gviz = ['digraph {'
      'splines=""line""'
      'ranksep=""1""'
      'node [shape=circle,fontcolor=white,color=dimgray,label="""",height=""0.3""]'];
    nt_times = '{ node [shape=diamond, height=""0.3"", width=""0.3"", style=filled, color=orangered, fontcolor=orangered, fillcolor=orangered,label=""""] ';
  end
  for i=1:max_n_nodes
    if nodes_in_graph(i) == 1
      if node_type(i) == '*'
        nt_times = nt_times + ' ' + nodename(i);
      end
    end
  end
  nt_times = [nt_times; '}'];

  if sum(nodes_in_graph) < 40
    nt_plus = '{ node [shape=square, height=""0.3"", width=""0.3"", style=filled, color=steelblue4, fontcolor=steelblue4] ';
  else
  nt_plus = '{ node [shape=square, height=""0.3"", width=""0.3"", style=filled, color=steelblue4, fontcolor=steelblue4, fillcolor=steelblue4,label=""""] ';
  end
  for i=1:max_n_nodes
    if nodes_in_graph(i) == 1
      if node_type(i) == '+'
        nt_plus = nt_plus + ' ' + nodename(i);
      end
    end
  end
  nt_plus= [nt_plus; '}'];
  gviz = [gviz; nt_times];    
  gviz = [gviz; nt_plus];    

  for r=1:size(rank_by_name)
    ranks = '{ rank = same; '
    for i=1:size(rank_by_name(r))
      ranks = ranks + rank_by_name(r)(i) + ';';
    end
    ranks = ranks + ' }'
    gviz = [gviz; ranks];
  end
  
  for i=1:max_n_nodes
    if nodes_in_graph(i) == 1
      if ~isempty(dag_to(i));
        vn = nodename(i);
        for k=1:size(dag_to(i))
          nt = node_type(dag_to(i)(k));
          if nt == '*'
            colorstr = 'color=orangered, fontcolor=orangered'
          else 
            colorstr = 'color=steelblue4, fontcolor=steelblue4'
          end
//          gviz = [gviz; vn + ' -> ' + nodename(dag_to(i)(k)) + '[label=""' + nt + '"", fontsize=""24.0""' + colorstr + '];'];
          gviz = [gviz; vn + ' -> ' + nodename(dag_to(i)(k)) + '[' + colorstr + '];'];
        end
      end
    end
  end
  gviz = [gviz; '}'];
  unix('rm -f chicago-tmp-rank.dot');
  write('chicago-tmp-rank.dot',gviz);
endfunction 

function determine_size()
endfunction

//function print_vectorized_topo_order()
//endfunction


//txt=mgetl('chicago.bare');

exec chicago_bare.sce;
//exec dag_input.sce;

for i=1:size(txt,1)
  if isempty(tokens(txt(i)))
    continue
  end
  slpdag_add(txt(i))
end

gv = export_graphviz();

translation = [];
for i=1:n_total
  if nodes_in_graph(i) == 1
    txt = nodename(i);
  else
    txt = 'not in graph';
  end
  translation = [translation; txt];
end

//print_vectorized_topo_order();
exec tsort.sce;
gstat = [string((1:n_total)') string(nodes_in_graph) translation node_type];
