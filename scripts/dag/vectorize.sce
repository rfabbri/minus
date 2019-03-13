vslp = [];
// regenerate the vector expressions by traversing ranks
for i=2:size(rank)
  // Vectorize plus
  vvec = [];
  wvec = [];
  extract = [];
  for k=1:size(rankplus(i))
    vid = k-1;
    p = rankplus(i)(k);
    [lhs, rhs, op] = get_vars(node_expr(p));
    vvec = [vvec; 'v[' + string(vid) + '] = ' + rhs(1)];
    remaining_rhs = tokens(node_expr(p));
    remaining_rhs = remaining_rhs(5:$);
    wvec = [wvec; 'w[' + string(vid) + '] = ' + strcat(remaining_rhs)];
    extract = [extract; lhs + ' = v[' + string(vid) + '];']; // v *= v
    size(extract)
  end

  vslp = [vslp; vvec];
  vslp = [vslp; wvec];
  if (~isempty(vvec))
    vslp = [vslp; 'v += w'];
  end
  vslp = [vslp; extract];
  
  // Vectorize *
  vvec = [];
  wvec = [];
  extract = [];
  for k=1:size(ranktimes(i))
    vid = k-1;
    p = ranktimes(i)(k);
    [lhs, rhs, op] = get_vars(node_expr(p));
    vvec = [vvec; 'v[' + string(vid) + '] = ' + rhs(1)];
    remaining_rhs = tokens(node_expr(p));
    remaining_rhs = remaining_rhs(5:$);
    wvec = [wvec; 'w[' + string(vid) + '] = ' + strcat(remaining_rhs)];
    extract = [extract; lhs + ' = v[' + string(vid) + '];']; // v *= v
    p
    size(extract)
  end

  vslp = [vslp; vvec];
  vslp = [vslp; wvec];
  if (~isempty(vvec))
    vslp = [vslp; 'v *= w'];
  end
  vslp = [vslp; extract];
end
  
