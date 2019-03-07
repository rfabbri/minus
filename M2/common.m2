-- IMPORTS 
debug needsPackage "NumericalAlgebraicGeometry"
debug SLPexpressions

-- fold along rows
rfold = L -> if (#L ==0) then random(FF^0,FF^0) else fold(L, (a,b) -> a||b)

-- fold along cols
cfold = L -> fold(L, (a,b) -> a|b)

-- scalar triple product of three column vectors
stp = method()
stp (GateMatrix, GateMatrix, GateMatrix) := (a,b,c) -> (transpose c) * crossProduct(a,b)
stp GateMatrix := M -> stp(M_{0},M_{1},M_{2})

-- FUNCTIONS

-- module generators for the kernel of a 2x4 gateMatrix M
kerGens24 = M -> gateMatrix {{M_(0,2)*M_(1,1)-M_(0,1)*M_(1,2), 0, M_(0,3)*M_(1,2)-M_(0,2)*M_(1,3), M_(0,3)*M_(1,1)-M_(0,1)*M_(1,3)}, {-M_(0,2)*M_(1,0)+M_(0,0)*M_(1,2), M_(0,3)*M_(1,2)-M_(0,2)*M_(1,3), 0, -M_(0,3)*M_(1,0)+M_(0,0)*M_(1,3)}, {M_(0,1)*M_(1,0)-M_(0,0)*M_(1,1), -M_(0,3)*M_(1,1)+M_(0,1)*M_(1,3), -M_(0,3)*M_(1,0)+M_(0,0)*M_(1,3), 0}, {0, M_(0,2)*M_(1,1)-M_(0,1)*M_(1,2), M_(0,2)*M_(1,0)-M_(0,0)*M_(1,2), M_(0,1)*M_(1,0)-M_(0,0)*M_(1,1)}}

size GateMatrix := M -> (numrows M, numcols M)
size Matrix := M -> (numrows M, numcols M)

-- evaluate a gateMatrix G at a matrix x
evaluate (GateMatrix, Matrix) := (G, x) -> (
    M := mutableMatrix(FF,numrows G,numcols G);
    E' := makeEvaluator(G,matrix{cameraVars|dataParams});
    evaluate(E',mutableMatrix(x),M);
    matrix M
    )

--random diagonal matrix
randDiag = n -> diagonalMatrix for i from 1 to n list random CC

dehomogenize = method(Options=>{})
dehomogenize (Matrix, ZZ) := o -> (v, n) -> (
    --assumes column vector
    (1/v_(n, 0))*v^(toList splice(0..n-1,n+1..numrows v-1))
    )
dehomogenize Matrix := o -> v -> dehomogenize(v, numrows v -1)

summary = L -> (
    n := #L;
    H := sort L;
    Q1 := (1/2) * (H#(floor((n-1)/3))+H#(ceiling((n-1)/3)));
    med := (1/2) * (H#(floor((n-1)/2))+H#(ceiling((n-1)/2)));
    Q3 := (1/2) * (H#(floor(2*(n-1)/3))+H#(ceiling(2*(n-1)/3)));
    mean := (sum L)/n;
    var := sum(L/(x-> (x - mean)^2))/(n-1);
    << "Min: " << toString(min L) << endl;
    << "1Q: " << toString(Q1) << endl;
    << "Med: " << toString(med) << endl;
    << "Avg: " << toString(sub(mean,RR)) << endl;
    << "3Q: " << toString(Q3) << endl;
    << "Max: " << toString(max L) << endl;
    << "Std Dev: " << toString(sqrt(var)) << endl;
    )


    

-- random element in the kernel of M
randKernel = method(Options=>{Tolerance=>1e-4})
randKernel (Matrix, InexactFieldFamily) := o -> (M, FF) -> (
    K := numericalKernel(M, o.Tolerance);
    K*random(FF^(numcols K), FF^1)
    )
randKernel Matrix := o -> M -> randKernel(M, CC)

-- write starting parameters and solutions to file
writeStartSys = method(Options=>{Filename=>"startSys"})
writeStartSys (Matrix, List) := o -> (M, sols) -> writeStartSys(point M, sols, o)
writeStartSys (Point, List) := o -> (p, sols) -> (
   assert(instance(o.Filename,String));
   f := openOut o.Filename;
   f << "Parameter values: " << endl;
   f << toExternalString p << endl;
   f << "Solutions : " << endl;
   for s in sols do f << toExternalString s << endl;
   close f;
   )

readStartSys = filename -> (
    l := separate("\n", get filename);
    p0 := value l#1;
    sols := for i from 3 to #l-2 list value l#i;
    (p0, sols)
    )

adjugate = method()
adjugate Thing := M -> (
    n := numcols M;
    assert(n == numrows M);
    matrix table(n,n,(i,j)->((-1)^(i+j))*det submatrix(M,delete(j,toList(0..n-1)),delete(i,toList(0..n-1))))
    )

-- not printing to high precision -- deprecated?
sol2String = p -> replace("\\{|\\}","",toString p.Coordinates)

-- produces gates for "small" determinants"
det2 = M -> M_(0,0)*M_(1,1)-M_(1,0)*M_(0,1)
det3 = M -> M_(0,0)*det2(M_{1,2}^{1,2})-M_(0,1)*det2(M_{0,2}^{1,2})+M_(0,2)*det2(M_{0,1}^{1,2})
det4 = M -> M_(0,0)*det3(M_{1,2,3}^{1,2,3})-M_(0,1)*det3(M_{0,2,3}^{1,2,3})+M_(0,2)*det3(M_{0,1,3}^{1,2,3})-M_(0,3)*det3(M_{0,1,2}^{1,2,3})

laplaceDet = M -> (
    (m, n) := size M;
    if (m=!=n) then error("not square matrix")
    else if (m>5) then error("no Laplace for matrices larger than 4x4")
    else if (m==2) then det2 M
    else if (m==3) then det3 M
    else -* m==4 *- det4 M
    )

-- jacobian of GateMatrix wrt. a list of inputGates
jacobian (GateMatrix, List) := (F,inGates) -> fold(apply(inGates,g->diff(g,F)),(a,b)->a|b)

-- get rotation matrix from cayley parameters
cay2R = method(Options=>{Normalized=>false})
cay2R (Thing,Thing,Thing) := o -> (X,Y,Z) -> (
    if instance(X, RR) then x := X_CC else x = X;
    if instance(Y, RR) then y := Y_CC else y = Y;
    if instance(Z, RR) then z := Z_CC else z = Z;
    M := matrix{
    {1+x*x-(y*y+z*z), 2*(x*y-z), 2*(x*z+y)},
    {2*(x*y+z), 1+y^2-(x*x+z*z), 2*(y*z-x)},
    {2*(x*z-y), 2*(y*z+x), 1 +z*z -(x*x+y*y)}
	};
    if o.Normalized then (1/(1+x^2+y^2+z^2)) * M else M
    )
cay2R List := o -> L -> cay2R(L#0, L#1, L#2, o)
cay2R Matrix := o -> M -> cay2R flatten entries M

-- get Cayley parameters from rotation matrix
R2Cay = method(Options=>{UnNormalize=>false})
R2Cay Matrix := o -> R -> (
    assert(numcols R == 3);
    assert(numrows R == 3);
    S := (R-id_(CC^3))*(R+id_(CC^3))^-1;
    (S_(2,1), S_(0,2), S_(1,0))
    )

-*/// TEST
restart
needs "common.m2"
(x, y, z) = (random CC, random CC, random CC)
R = cay2R(x, y, z, Normalized=>true)
R2Cay R
(x,y,z)

ab=
ac=
ad=
a^2+( / a)^2 + .. = 1
a^4 - a^2 + sumsof blahs squared =9

///*-

-- get rotation matrix from quaternion parameters
Q2R = method(Options=>{Normalized=>false, FF=>CC})
Q2R (Thing,Thing,Thing, Thing) := o -> (W, X,Y,Z) -> (
    if instance(W, RR) then w := W_CC else w = W;
    if instance(X, RR) then x := X_CC else x = X;
    if instance(Y, RR) then y := Y_CC else y = Y;
    if instance(Z, RR) then z := Z_CC else z = Z;
    M := matrix{
    {w*w+x*x-(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)},
    {2*(x*y+w*z), w^2+y^2-(x*x+z*z), 2*(y*z-w*x)},
    {2*(x*z-w*y), 2*(y*z+w*x), w^2 +z*z -(x*x+y*y)}
	};
    if o.Normalized then (1/(w^2+x^2+y^2+z^2)) * M else M
    )
Q2R List := o -> L -> Q2R(L#0, L#1, L#2, L#3, o)
Q2R Matrix := o -> M -> Q2R flatten entries M

-- get Cayley parameters from rotation matrix
R2Q = method(Options=>{UnNormalize=>false,FF=>CC})
R2Q Matrix := o -> R -> (
    assert(numcols R == 3);
    assert(numrows R == 3);
    c := (R_(2,1) - R_(1,2));
    b := (R_(0,2) - R_(2,0));
    a := (R_(1,0) - R_(0,1));
    w := (1/2)*sqrt(R_(0,0)+R_(1,1)+R_(2,2)+1);
    x := 1/(4*w) * c;
    y := 1/(4*w) * b;
    z := 1/(4*w) * a;
--    << w^2+x^2+y^2+z^2 << endl;
    (w, x, y, z)
    )

-*/// TEST
R=CC[W]
netList solveSystem {W^4-W^2+1/16}

clean T
T=QQ[a..d]
R=Q2R gens T
S = (R-id_(((QQ)^3)))*adjugate(R+id_((QQ)^3));
S
((first x)/(first L))*L
1/sqrt(sum(x/(y->y^2)))*x
L
///*-


-- cross product of col vectors -- takes Matrice or GateMatrix pair
crossProduct = (y,q) -> matrix{{y_(1,0)*q_(2,0)-y_(2,0)*q_(1,0)},{y_(2,0)*q_(0,0)-y_(0,0)*q_(2,0)},{y_(0,0)*q_(1,0)-y_(1,0)*q_(0,0)}}

--
randomLineThroughPoints = (P, FF) -> ( 
    m := numrows P; -- m = dim + 1
    n := numcols P; -- n = number of points
    assert(n<=2 and m>=3 and m<=4);
    K := numericalKernel(transpose P,1e-6);
    assert(numcols K == m-n); -- true if points are distinct
    transpose(K * random(FF^(m-n),FF^(m-2)))
    )

-- constructs problem data given a PL diagram D (complete visibility)
fabricatePair = (D, FF, nparams) -> (    
    (nLines,nGhosts,intersections) := D;
    worldPoints := sub(random(FF^4,FF^(#intersections)),CC);
    pointsOnLineIndices := apply(nLines+nGhosts, l->positions(intersections,i->member(l,i)));
    helperPoints := apply(pointsOnLineIndices, pp->sub(random(FF^4,FF^(2-#pp)),CC));
    sampleCameraParameters := for i from 1 to nvars list sub(random FF,CC);
    subTable := apply(sampleCameraParameters, cameraVars, (a,b) -> b=>inputGate a);
    sampleC := apply(C,cam -> (
	    M := mutableMatrix(CC, 3, 4);
	    evaluate(cam, mutableMatrix{sampleCameraParameters}, M);
	    matrix M
	    )
	    );
    (
	sampleCameraParameters,
	apply(sampleC, cam->(
	    	P := cam * worldPoints;
	    	L := matrix apply(nLines+nGhosts, l->(
		      Hl := helperPoints#l;
		      Pl := P_(pointsOnLineIndices#l);
		      Pl = Pl|cam*Hl;
	    	      line := randomLineThroughPoints(Pl, FF);
		      {(1/norm(2,line))*line} 
		       ));
	       (P,L)
	       ))
	) 
    )
-- convenience functions for minors
minors (GateMatrix, ZZ, Sequence, Boolean) := o ->  (M, k, S, laplace) -> (
    (Sm, Sn) := (first S, last S);
    (m,n) := (numrows M, numcols M);
    assert(k<=min(m,n));
    assert(all(Sm,s->#s==k));
    assert(all(Sn,s->#s==k));
    flatten apply(Sm,sm->apply(Sn, sn -> 
	    if (laplace) then laplaceDet submatrix(M,sm,sn)
	    else det submatrix(M,sm,sn)
	    ))
    )

allMinors = method(Options=>{Laplace=>false})
allMinors (GateMatrix, ZZ) := o -> (M, k) -> (
    (m, n ) := (numrows M, numcols M);
    s := (subsets(0..m-1,k),subsets(0..n-1,k));
    minors(M, k, s, o.Laplace)
    )

maxMinors = method(Options=>{Laplace=>false})
maxMinors GateMatrix := o -> M -> allMinors(M,min(numrows M, numcols M), Laplace=>o.Laplace)

-- indexes subsystem giving Jacobian rank
rowSelector = J0 -> (
    (m , n) := size J;
    J0' := J0^{0}; -- should technically check this is not identically zero
    i := 1;
    inds := new MutableList from {0};
    k := numericalRank J0';
    while ( (i<m) and (k < n)) do (
	if (numericalRank(J0'||J0^{i}) > k) then (
	    inds#k = i;
	    k = k+1;
	    J0' = J0'||J0^{i};
	    );
	i = i + 1;
	);
    if k < n then  << "WARNING: Jacobian has rank" << toString k << endl;
    toList inds
    )

log10 = x -> log(x)/log(10)

argCC = z -> atan((imaginaryPart z)/(realPart z))

-- complex number whose real and imag parts are standard normal
gaussCC = () -> (
    (u1,u2):=(random RR,random RR);
    sqrt(-2*log(u1))*cos(2*pi*u2)+ii*sqrt(-2*log(u1))*sin(2*pi*u2)
    )

-- random sample drawn from normal distriution N(mu, var^2)
rNorm = (mu,var) -> mu+var*(realPart gaussCC())_CC

-- random sample from (n-1)-sphere with radius r
sphere = (n,r) -> (
    l:=apply(n,i->rNorm(0,1));
    matrix{r/norm(2,l)*l}
    )

-- assumes "u" of unit length
householder=method()
householder (Matrix,ZZ) := (u,n) -> (
    if (numrows u > 1) then error("householder takes a row vector");
    R:=ring u;
    k:=numcols u;
    id_(R^(n-k))++(id_(R^k)-2*(transpose u)*u)
    )
householder (List,ZZ) := (u,n) -> householder(matrix {u},n)

randomOn = n -> diagonalMatrix(toList((n-1):1_RR)|{(-1)^(random 2)}) * fold(reverse apply(2..n,i->householder(sphere(i,1),n)),(a,b)->a*b)

randomCameraNormalized = () -> (
    R := randomOn 3;
    t := matrix{{random CC},{random CC},{random CC}};
--    t := transpose matrix{sphere(3,1)};
    tnorm := (1 / t_(2,0))* t;
    (R|tnorm)
    )

randomCameraNormalizedCayley = () -> (
    R := cay2R(random CC, random CC, random CC,Normalized=>true);
    t := matrix{{random CC},{random CC},{random CC}};
--    t := transpose matrix{sphere(3,1)};
    tnorm := (1 / t_(2,0))* t;
    (R|tnorm)
    )


randomCamera = () -> (
    R := randomOn 3;
    t := transpose matrix{sphere(3,1)};
    (R|t)
    )
