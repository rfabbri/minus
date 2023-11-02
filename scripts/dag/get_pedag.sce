// -----------------------------------------------------------------------------
// This gets the direct expressions of all nodes from the given one down to the first
// removing the nodes that are not used by anyone after this expansion
// 
// p=pedag(get_id('G8')); [p,gstat(1:size(p,1),:)]

p = pedag(get_id('G3439'));
p_c = p(p ~= "")
pnames = txt_lhs(p ~= "")
p_c = "const C<F> " + pnames + " = " + p_c + ";"

fd=mopen('out/Hxt/partial_expr.c',"wt");
for i=1:size(p_c,1)
  mfprintf(fd,p_c(i));
  mfprintf(fd,"\n");
end
mclose(fd);
