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

exec chicago_bare_new.sce;
//exec dag_input.sce;
//txt=mgetl('chicago.bare');
txtc = strcat(txt,"",'r');
n_Xs = evstr(unix_g("echo " + "''"+ txtc + "''" + "|sed ''s/;/;\n/g'' | egrep -o X[0-9]+ |sort|tail -n 1|grep -o [0-9]*")(:))+1
n_Cs = evstr(unix_g("echo " + "''"+ txtc + "''" + "|sed ''s/;/;\n/g'' | egrep -o C[0-9]+ |sort|tail -n 1|grep -o [0-9]*")(:))+1
n_Gs = evstr(unix_g("echo " + "''"+ txtc + "''" + "|sed ''s/;/;\n/g'' | egrep -o G[0-9]+ |sort|tail -n 1|grep -o [0-9]*")(:))+1
//n_Xs = evstr(unix_g("echo " + "''"+ txtc + "''" + "|sed ''s/;/;\n/g'' | egrep -o X[0-9]+ |sort|wc -l")(:));
//n_Cs = evstr(unix_g("echo " + "''"+ txtc + "''" + "|sed ''s/;/;\n/g'' | egrep -o C[0-9]+ |sort|wc -l")(:));
//n_Gs = evstr(unix_g("echo " + "''"+ txtc + "''" + "|sed ''s/;/;\n/g'' | egrep -o G[0-9]+ |sort|wc -l")(:));
// n_Cs = 3;
//n_Xs = 3;
//n_Gs = 4; 
//n_Xs = 8;
//n_Gs = 47;
//n_Xs = 127;
//n_Gs = 3440;
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

exec dag_incl.sce;
//dag = list();
//dag(1) = dag_from;
//dag(2) = dag_to;

// Populate DAG
for i=1:size(txt,1)
  if isempty(tokens(txt(i)))
    continue
  end
  slpdag_add(txt(i))
end

// gv = export_graphviz();

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
// exec tsort.sce;
gstat = [string((1:n_total)') string(nodes_in_graph) translation node_type];
