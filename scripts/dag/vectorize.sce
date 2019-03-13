poweroftwo = %t;
min_ranksize = 1000;
vslp = [];
// regenerate the vector expressions by traversing ranks
// TODO: minimum size
for i=2:size(rank)

  istring = string(i);
  rsize = size(rankplus(i));
  if rsize ~= 0
    if (poweroftwo)
      vsizen = 2^ceil(log2(rsize));
    else
      vsizen = rsize;
    end

    if vsizen > min_ranksize
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
    else
      for k=1:rsize
        p = rankplus(i)(k);
        vslp = [vslp; 'const complex ' + node_expr(p)];
      end
    end
  end

  
  // Vectorize *
  rsize = size(ranktimes(i));
  if rsize ~= 0
    if (poweroftwo)
      vsizen = 2^ceil(log2(rsize));
    else
      vsizen = rsize;
    end
    
    if vsizen > min_ranksize
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
    else
      for k=1:rsize
        p = ranktimes(i)(k);
        vslp = [vslp; 'const complex ' + node_expr(p)];
      end
    end
  end
end
  
vslp
