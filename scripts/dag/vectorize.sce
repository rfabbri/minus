// desired_sizes = [3  4  6  18  20  32  42  44  52];
desired_sizes = 10000;
poweroftwo = %f;
//min_ranksize = 10;
//max_ranksize = min_ranksize+30;
//min_ranksize = 0
//max_ranksize = 1000
vslp = [];
vslp = [vslp; 'using namespace Eigen;'];

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

    if sum(vsizen == desired_sizes) == 1
//    if vsizen >= min_ranksize & vsizen <= max_ranksize
//    if i == desired_rank
      vsize = string(vsizen);
      vslp = [vslp; 'complex vp' + string(i) + '[' + vsize + '],' + 'wp' + string(i) + '[' + vsize + '];'];
      estr='Map<Array <complex, ' + vsize +', 1> > evp' + istring + '(vp'+istring+'), ewp' + istring + '(wp'+istring+');';
      vslp = [vslp; estr];
      // Vectorize plus
      vvec = [];
      wvec = [];
      extract = [];
      for k=1:size(rankplus(i))
        vid = k-1;
        p = rankplus(i)(k);
        [lhs, rhs, op] = get_vars(node_expr(p));
        vvec = [vvec; 'vp'+string(i)+'[' + string(vid) + '] = ' + rhs(1) + ';'];
        remaining_rhs = tokens(node_expr(p));
        remaining_rhs = remaining_rhs(5:$);
        wvec = [wvec; 'wp' + string(i)+ '[' + string(vid) + '] = ' + strcat(remaining_rhs)];
        extract = [extract; 'const complex &' + lhs + ' = vp' + string(i) + '[' + string(vid) + '];'];
        size(extract)
      end
      vslp = [vslp; vvec];
      vslp = [vslp; wvec];
      vslp = [vslp; 'evp' + istring + ' += ewp' + istring + ';'];
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
    
//    if vsizen >= min_ranksize & vsizen <= max_ranksize
    if sum(vsizen == desired_sizes) == 1
//    if i == desired_rank
      vsize = string(vsizen);
      vslp = [vslp; 'complex vm' + string(i) + '[' + vsize + '],' + 'wm' + string(i) + '[' + vsize + '];'];
      disp 'mvslp:'
      vslp
      size(vslp)
      estr='Map<Array<complex, ' + vsize +', 1> > evm' + istring + '(vm'+istring+'), ewm' + istring + '(wm'+istring+');';
      vslp = [vslp; estr];
      disp 'mvslp:'
      vslp
      vvec = [];
      wvec = [];
      extract = [];
      for k=1:size(ranktimes(i))
        vid = k-1;
        p = ranktimes(i)(k);
        [lhs, rhs, op] = get_vars(node_expr(p));
        vvec = [vvec; 'vm' + string(i) + '[' + string(vid) + '] = ' + rhs(1) + ';'];
        remaining_rhs = tokens(node_expr(p));
        remaining_rhs = remaining_rhs(5:$);
        wvec = [wvec; 'wm' + istring + '[' + string(vid) + '] = ' + strcat(remaining_rhs)];
        extract = [extract; 'const complex &' + lhs + ' = vm' + string(i) + '[' + string(vid) + '];'];
      end
      vvec
      vslp = [vslp; vvec];
      vslp = [vslp; wvec];
      vslp = [vslp; 'evm' + istring + ' *= ewm' + istring + ';'];
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
cxxfile = 'chicago-eval-sorted-v64.cxx';
unix('rm -f ' + cxxfile);
write(cxxfile,vslp);
