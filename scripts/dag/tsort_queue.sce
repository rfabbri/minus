Q = Qplus;
Qplus = list();
Qtype = '+';
Qplus1 = list();
Qtimes1 = list();
while ~isempty(Q)
  vvec = [];
  wvec = [];
  extract = [];
  vid = 0;
  cvec = [];
  while ~isempty(Q)
    // pop node p from Q
    p = Q($); Q($) = null();
    disp ('doing node' + string(p) + ' : ' + nodename(p))

    select node_type(p)
    case '+'
      disp '===================='
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
        // can have multiple edges, as in squaring op
        if dag_from(q)(jj) == p
          dag_from(q)(jj) = null();
        end
      end
      // q becomes a leaf, add it to queue
      if isempty(dag_from(q))
        select node_type(q)
        case '+' then
          Qplus1($+1) = q;
        case '*' then
          Qtimes1($+1) = q;
        else 
          error('node type must be plus or minus')
        end
      end
    end
    vid = vid+1;
  end
  if Qtype == '+'
    if ~isempty(Qtimes)
      Q = Qtimes;
      Qtimes = list();
      Qtype = '*';
    elseif ~isempty(Qplus1)
      Q = Qplus1;
      Qplus1 = list();
      Qtype = '+';
    else
      Q = Qtimes1;
      Qtimes1 = list();
      Qtype = '*';
    end
  elseif Qtype == '*' 
    Q = Qplus1;
    Qplus1 = list();
    Qtimes = Qtimes1();
    Qtimes1 = list();
    Qtype = '+'
  end
  vvec
  wvec
  extract
  cvec
end
