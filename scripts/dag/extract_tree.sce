global extracted;
extracted = zeros(n_total,1);

function nodes_dependent_on_i = extract_tree(i)
  global extracted;
  nodes_dependent_on_i = '';
  if ~extracted(i)
    extracted(i) = 1;
    nodes_dependent_on_i = nodename(i) + ';';
    for k=1:length(dag_to(i))
       nodes_dependent_on_i = nodes_dependent_on_i + extract_tree(dag_to(i)(k));
    end
  end
endfunction

et = list();
for i=0:14
  et(i+1) = extract_tree(get_id('X'+ string(i)));
end

ets = ''
for i=1:15
  ets = ets + et(i);
end

//write('xt-dep',ets);

fd=mopen('out/xt-dep',"wt");
mfprintf(fd,ets);
mclose(fd);

// s/;//g
// uniq xt-dep |sort -n|uniq >xt-dep-uniq
