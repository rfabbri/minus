Q = Qplus;
Qplus = list();
Qtype = '+';
Qplus1 = list();
Qtimes1 = list();
rank = list();
level = 1;
rank(1) = list();
while ~isempty(Q)
  vvec = [];
  wvec = [];
  extract = [];
  vid = 0;
  cvec = [];
  while ~isempty(Q)
    // pop node p from Q
    p = Q($); Q($) = null();
//    disp ('doing node' + string(p) + ' : ' + nodename(p))
    rank(level)($+1) = p;

    select node_type(p)
    case '+'
//      disp '===== ADD ================='
      [lhs, rhs, op] = get_vars(node_expr(p));
      vvec = [vvec; 'v[' + string(vid) + '] = ' + rhs(1)];
      remaining_rhs = tokens(node_expr(p));
      remaining_rhs = remaining_rhs(5:$);
      wvec = [wvec; 'w[' + string(vid) + '] = ' + strcat(remaining_rhs)];
      extract = [extract; lhs + ' = v[' string(vid) '];']; // v *= v
    case '*'
//      disp '===== TIMES ==============='
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
//          disp('inserting ' + nodename(q) + ' Qtype: ' + Qtype);
        case '*' then
          Qtimes1($+1) = q;
//          disp('inserting ' + nodename(q) + ' Qtype: ' + Qtype);
        else 
          error('node type must be plus or minus')
        end
      end
    end
    vid = vid+1;
  end
//  disp 'type after first loop'
  Qtype
  if Qtype == '+'
    if ~isempty(Qtimes)
      Q = Qtimes;
      Qtimes = list();
      Qtype = '*';
//      disp '0+ to 0*'
    elseif ~isempty(Qplus1)
      Q = Qplus1;
      Qtimes = Qtimes1;
      Qplus1 = list();
      Qtimes1 = list();
      Qtype = '+';
      level = level + 1;
      if (~isempty(Q))
        rank(level) = list();
      end
//      disp '0+ to 1+, reset to 0+'
    else
      Q = Qtimes1;
      Qtimes1 = list();
      Qtype = '*';
      level = level + 1;
      if (~isempty(Q))
        rank(level) = list();
      end
//      disp '0+ to 1* -> 0*'
    end
  else // *
    level = level + 1;
    if ~isempty(Qplus1)
      Q = Qplus1;
      Qtimes = Qtimes1;
      Qplus1 = list();
      Qtimes1 = list();
      Qtype = '+';
//      disp '0* to 1+ -> 0+'
    else
      Q = Qtimes1;
      Qtimes1 = list();
      Qtype = '*';
//      disp '0* to 1* -> 0*'
    end
    if (~isempty(Q))
      rank(level) = list();
    end
  end
  vvec
  wvec
  extract
  cvec
end

// TODO: sort X's and nearby variable indices to adjacent vector positions

// draw graph with ranks

size(aplat(rank))
sum(nodes_in_graph)
rs = zeros(size(rank),1);
rsplus = zeros(size(rank),1);
rstimes = zeros(size(rank),1);
rsconst = zeros(size(rank),1);
rankplus=list();
ranktimes=list();
rank_by_name = rank;
for i=1:size(rank)
  rs(i) = size(rank(i))
  rankplus(i) = list();
  ranktimes(i) = list();
  for k=1:rs(i)
    rank_by_name(i)(k) = nodename(rank(i)(k));
    select node_type(rank(i)(k));
    case '+'
      rsplus(i) = rsplus(i) + 1;
      rankplus(i)($+1) = rank(i)(k);
    case '*'
      rstimes(i) = rstimes(i) + 1;
      ranktimes(i)($+1) = rank(i)(k);
    else
      rsconst(i) = rsconst(i) + 1;
    end
  end
end

//exec dag.sce;
//gv = export_graphviz_with_ranks();

exec vectorize.sce;
size(vslp)
//exec dag_input.sce;
//size(txt)
