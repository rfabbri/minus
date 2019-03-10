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

