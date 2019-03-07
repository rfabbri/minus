debug needsPackage "SLPexpressions"
needs "./chicago/common.m2"
inv2 = A -> idGate/(A_(0,0)*A_(1,1)-A_(1,0)*A_(0,1)) * gateMatrix{{A_(1,1),-A_(0,1)},{-A_(1,0),A_(0,0)}}
idGate = inputGate 1
-*
a=inputGate a
b=inputGate b
c=inputGate c
d=inputGate d
g=gateMatrix{{a,b},{c,d}}
inv2 g
*-

inv = A -> (
    (m, n) := size A;
    if (m =!= n) then error("NOT SQUARE");
    if (m == 1) then gateMatrix{{idGate/A_(0,0)}}
    else if (m == 2) then inv2 A 
    else (
	k := floor(m/2);
	(A11, A12, A21, A22) := (A_{0..k}^{0..k},A_{k+1..n-1}^{0..k},A_{0..k}^{k+1..n-1},A_{k+1..n-1}^{k+1..n-1});
	I := inv A11;
	II := A21 * I;
	III := I * A12;
	IV := A21 * III;
	V := IV - A22;
	VI := inv V;
	C12 := III * VI;
	C21 := VI * II;
	VII := III * C21;
	C11 := I - VII;
	C22 := -VI;
	(C11|C12)||(C21|C22)
	)
    )

adj3 = M -> (
    n := numcols M;
    assert(n == numrows M);
    gateMatrix table(n,n,(i,j)->((-1)^(i+j))*det2(M_(delete(j,toList(0..n-1)))^(delete(i,toList(0..n-1)))))
    )

adj4 = M -> (
    n := numcols M;
    assert(n == numrows M);
    gateMatrix table(n,n,(i,j)->((-1)^(i+j))*det3(M_(delete(j,toList(0..n-1)))^(delete(i,toList(0..n-1)))))
    )

adjInv = A -> (
    (m, n) := size A;
    if (m =!= n) then error("NOT SQUARE");
    if (m == 1) then gateMatrix{{idGate/A_(0,0)}}
    else if (m == 2) then inv2 A 
    else if (m == 3) then idGate/det3(A) * adj3 A
--    else if (m == 4) then idGate/det4(A) * adj4 A
    else (
	k := floor(m/2);
	(split1,split2)=({0..k-1},{k..n-1});
	(A11, A12, A21, A22) := (A_split1^split1,A_split2^split1,
	                        A_split1^split2,A_split2^split2);
	I := adjInv A11;
	II := A21 * I;
	III := I * A12;
	IV := A21 * III;
	V := IV - A22;
	VI := adjInv V;
	C12 := III * VI;
	C21 := VI * II;
	VII := III * C21;
	C11 := I - VII;
	C22 := -VI;
	(C11|C12)||(C21|C22)
	)
    )

adjInv2 = A -> (
    (m, n) := size A;
    if (m =!= n) then error("NOT SQUARE");
    if (m == 1) then gateMatrix{{idGate/A_(0,0)}}
    else if (m == 2) then inv2 A 
    else if (m == 3) then (
	<< "hello 3 " << endl;
	idGate/det3(A) * adj3 A
	)
    else if (m == 4) then (
	<< "hello 4!" << endl;
	idGate/det4(A) * adj4 A
	)
    else (
	k := floor(m/2);
	<< k << " is the breakpoint " << endl;
	(split1,split2)=({0..k-1},{k..n-1});
	(A11, A12, A21, A22) := (A_split1^split1,A_split2^split1,
	                        A_split1^split2,A_split2^split2);
	I := adjInv2 A11;
	II := A21 * I;
	III := I * A12;
	IV := A21 * III;
	V := IV - A22;
	VI := adjInv2 V;
	C12 := III * VI;
	C21 := VI * II;
	VII := III * C21;
	C11 := I - VII;
	C22 := -VI;
	(C11|C12)||(C21|C22)
	)
    )


n=14
scan(flatten for i from 1 to n list for j from 1 to n list x_(i,j),
    s->value(toString(s)|" = inputGate " | toString(s)))
gnn=gateMatrix for i from 1 to n list for j from 1 to n list x_(i,j)    	
M = mutableMatrix(CC,n,n)
end--

restart
needs "inverse.m2"
iG=inv gnn;
aG =adjInv gnn;
aG2 =adjInv2 gnn;
{iG,aG,aG2}/depth
recursionLimit = 500
h1=cCode(iG,gnn)--7678
h2=cCode(aG,gnn)--7822
h2=cCode(aG2,gnn)--7822
elapsedTime E' = makeEvaluator(inv gnn,gnn)
elapsedTime E = makeEvaluator(adjInv gnn,gnn)
elapsedTime E'' = makeEvaluator(adjInv2 gnn,gnn)

x=random(CC^n,CC^n)
elapsedTime evaluate(E',mutableMatrix(x),M);
elapsedTime evaluate(E,mutableMatrix(x),M);
elapsedTime evaluate(E'',mutableMatrix(x),M);
elapsedTime inverse x;


x
evaluate(gnn,random(CC^16,CC^1))
