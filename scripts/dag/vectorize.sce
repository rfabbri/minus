poweroftwo = %t;
vslp = [];
// regenerate the vector expressions by traversing ranks
// TODO: minimum size
for i=2:size(rank)

  istring = string(i);
  if size(rankplus(i)) ~= 0
    if (poweroftwo)
      vsizen = 2^ceil(log2(size(rankplus(i))))
    else
      vsizen = size(rankplus(i))
    end
    vsize = string(vsizen);
    vslp = [vslp; 'complex v' + string(i) + '[' + vsize + '],' + 'w' + string(i) + '[' + vsize + '];'];
    // Vectorize plus
    vvec = [];
    wvec = [];
    extract = [];
    for k=1:size(rankplus(i))
      vid = k-1;
      p = rankplus(i)(k);
      [lhs, rhs, op] = get_vars(node_expr(p));
      vvec = [vvec; 'v'+string(i)+'[' + string(vid) + '] = ' + rhs(1)];
      remaining_rhs = tokens(node_expr(p));
      remaining_rhs = remaining_rhs(5:$);
      wvec = [wvec; 'w' + string(i)+ '[' + string(vid) + '] = ' + strcat(remaining_rhs)];
      extract = [extract; 'const complex &' + lhs + ' = v' + string(i) + '[' + string(vid) + '];'];
      size(extract)
    end
    vslp = [vslp; vvec];
    vslp = [vslp; wvec];
    vslp = [vslp; 've' + istring + ' += we' + istring];
    vslp = [vslp; extract];
  end

  
  // Vectorize *
  if size(ranktimes(i)) ~= 0
    if (poweroftwo)
      vsizen = 2^ceil(log2(size(ranktimes(i))))
    else
      vsizen = size(ranktimes(i))
    end
    vsize = string(vsizen);
    vslp = [vslp; 'complex v' + string(i) + '[' + vsize + '],' + 'w' + string(i) + '[' + vsize + '];'];
    vvec = [];
    wvec = [];
    extract = [];
    for k=1:size(ranktimes(i))
      vid = k-1;
      p = ranktimes(i)(k);
      [lhs, rhs, op] = get_vars(node_expr(p));
      vvec = [vvec; 'v' + string(i) + '[' + string(vid) + '] = ' + rhs(1)];
      remaining_rhs = tokens(node_expr(p));
      remaining_rhs = remaining_rhs(5:$);
      wvec = [wvec; 'w' + istring + '[' + string(vid) + '] = ' + strcat(remaining_rhs)];
      extract = [extract; 'const complex &' + lhs + ' = v' + string(i) + '[' + string(vid) + '];'];
    end
    vvec
    vslp = [vslp; vvec];
    vslp = [vslp; wvec];
    vslp = [vslp; 've' + istring + ' *= we' + istring];
    vslp = [vslp; extract];
  end
end
  
vslp
