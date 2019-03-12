Q = Qplus;
Qplus = list();
Qtype = '0+';
Qplus1 = list();
Qtimes1 = list();
rank = list();
level = 1;
while ~isempty(Q)
  vvec = [];
  wvec = [];
  extract = [];
  vid = 0;
  cvec = [];
  rank(level) = list();
  while ~isempty(Q)
    // pop node p from Q
    p = Q($); Q($) = null();
    disp ('doing node' + string(p) + ' : ' + nodename(p))
    rank(level)($+1) = nodename(p);

    select node_type(p)
    case '+'
      disp '===== ADD ================='
      [lhs, rhs, op] = get_vars(node_expr(p));
      vvec = [vvec; 'v[' + string(vid) + '] = ' + rhs(1)];
      remaining_rhs = tokens(node_expr(p));
      remaining_rhs = remaining_rhs(5:$);
      wvec = [wvec; 'w[' + string(vid) + '] = ' + strcat(remaining_rhs)];
      extract = [extract; lhs + ' = v[' string(vid) '];']; // v *= v
    case '*'
      disp '===== TIMES ==============='
      [lhs, rhs, op] = get_vars(node_expr(p));
      vvec = [vvec; 'v[' + string(vid) + '] = ' + rhs(1)];
      remaining_rhs = tokens(node_expr(p));
      remaining_rhs = remaining_rhs(5:$);
      wvec = [wvec; 'w[' + string(vid) + '] = ' + strcat(remaining_rhs)];
      extract = [extract; lhs + ' = v[' string(vid) '];']; // v *= v
    else
      cvec = [cvec; nodename(p)];
    end


    // remove every outgoing link
    for k=size(dag_to(p)):-1:1
      // remove dag_to(p)(k)
      q = dag_to(p)(k); dag_to(p)(k) = null();
      // remove p form dag_from(q)
      for jj = size(dag_from(q)):-1:1
        // can have multiple edges, as in squaring op, treated in separately dag_to
        if dag_from(q)(jj) == p
          dag_from(q)(jj) = null();
          break;
        end
      end
      // q becomes a leaf, add it to queue
      if isempty(dag_from(q))
        select node_type(q)
        case '+' then
          Qplus1($+1) = q;
        case '*' then
          Qtimes1($+1) = q;
          disp('inserting ' + nodename(q) + ' Qtype: ' + Qtype);
        else 
          error('node type must be plus or minus')
        end
      end
    end
    vid = vid+1;
  end
  disp 'type after first loop'
  Qtype
  if part(Qtype,2) == '+'
    if part(Qtype,1) == '0' // 0+
      if ~isempty(Qtimes)
        Q = Qtimes;
        Qtimes = list();
        Qtype = '0*';
        disp '0+ to 0*'
      elseif ~isempty(Qplus1)
        Q = Qplus1;
        Qplus1 = list();
        Qtype = '1+';
        level = level + 1;
        disp '0+ to 1+'
      else
        Q = Qtimes1;
        Qtimes1 = list();
        Qtype = '1*';
        level = level + 1;
        disp '0+ to 1*'
      end
    else // 1+
      if ~isempty(Qtimes1)
        Q = Qtimes1;
        Qtimes1 = list();
        Qtype = '1*';
        disp '1+ to 1*'
      elseif ~isempty(Qplus)
        Q = Qplus;
        Qplus = list();
        Qtype = '0+';
        level = level + 1;
        disp '1+ to 0+'
      else
        Q = Qtimes;
        Qtimes = list();
        Qtype = '0*';
        level = level + 1;
        disp '1+ to 0*'
      end
    end
  else // *
    if part(Qtype,1) == '0' // 0*
      if ~isempty(Qplus1)
        Q = Qplus1;
        Qplus1 = list();
        Qtype = '1+';
        level = level + 1;
        disp '0* to 1+'
      else
        Q = Qtimes1;
        Qtimes1 = list();
        Qtype = '1*';
        level = level + 1;
        disp '0* to 1*'
      end
    else // 1*
      if ~isempty(Qplus)
        Q = Qplus;
        Qplus = list();
        Qtype = '0+';
        level = level + 1;
        disp '1* to 0+'
      else
        Q = Qtimes;
        Qtimes = list();
        Qtype = '0*';
        level = level + 1;
        disp '1* to 0*'
      end
    end
  end
  vvec
  wvec
  extract
  cvec
end

// TODO: sort X's to adjacent vector positions
