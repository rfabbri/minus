-- IMPORTS
needs "common.m2"
needsPackage "MonodromySolver"

-- GLOBALS
FF = CC

-- SET TRACKER OPTIONS HERE
-- null indicates default value 
scan({CorrectorTolerance=>null,
	EndZoneFactor=>null,
	InfinityThreshold => null, 
	maxCorrSteps => null, 
	NoOutput => null,
	numberSuccessesBeforeIncrease => null,
	Precision => null,
	Predictor => null,
	stepIncreaseFactor => null,
	tStep => null,
	tStepMin => null
	}, 
    opt -> setDefault(opt))

-- setting up input gates
-- "y" and "t" indices go (view, feature, coordinate)
-- "cq" indices go view, corresponding quaternion parameter (0 for constant)
varInSymbs = {x_2,y_2,z_2,x_3,y_3,z_3}
nvars=6
inSymbs=varInSymbs|
(for i from 1 to 3 list for j from 1 to 3 list for k from 1 to 3 list y_(i,j,k))//flatten//flatten |
(for i from 1 to 3 list for j from 1 to 2 list for k from 1 to 3 list t_(i,j,k))//flatten//flatten 
--(for i from 2 to 3 list for j from 0 to 3 list cq_(i,j))//flatten
scan(inSymbs,s->value(toString(s)|" = inputGate " | toString(s)))
inGates = toList(inSymbs)/value
cameraVars = take(inGates,nvars)
dataParams = drop(inGates,nvars)

-*
-- still using quaternions, just not as variables
w_2 = (cq_(2,0) + cq_(2,1)*x_2) + (cq_(2,2)*y_2 + cq_(2,3)*z_2)
w_3 = (cq_(3,0) + cq_(3,1)*x_3) + (cq_(3,2)*y_3 + cq_(3,3)*z_3)
*-

-- rotation and projection matrices
R11= gateMatrix(id_(FF^3))
R12 = compress cay2R(x_2,y_2,z_2)
R13 = compress cay2R(x_3,y_3,z_3)
Rots={R11,R12,R13}
RotEVs = { 
    makeEvaluator(R11, matrix{cameraVars}),
    makeEvaluator(R12, matrix{cameraVars}),
    makeEvaluator(R13, matrix{cameraVars})
  }

-- image points and tangents
imgPoints = (for j from 1 to 3 list for i from 1 to 3 list (
	transpose gateMatrix{for k from 1 to 3 list y_(i,j,k)}))
imgTans = (for j from 1 to 2 list for i from 1 to 3 list (
	transpose gateMatrix{for k from 1 to 3 list t_(i,j,k)}))

-- transfer points to first image plane
ptTransfers1=for j from 0 to 2 list imgPoints#j#0
ptTransfers2=for j from 0 to 2 list transpose R12 * imgPoints#j#1
ptTransfers3=for j from 0 to 2 list transpose R13 * imgPoints#j#2

-- matrices annihilating T1 and T2 , T1-T2 in kernel
T1ann = cfold for i from 0 to 2 list crossProduct(ptTransfers1#i,ptTransfers2#i);
T2ann = cfold for i from 0 to 2 list crossProduct(ptTransfers1#i,ptTransfers3#i);
--T1minusT2ann = cfold for i from 0 to 2 list crossProduct(ptTransfers2#i,ptTransfers3#i);--probably redundant
Tanneqs = {T1ann,T2ann}/stp//rfold;

-- matrices 
-- transfer points to first image plane
S=subsets(3,2)
lineTransfers12 = stp cfold for i from 0 to 2 list (transpose Rots#i) * crossProduct(imgPoints#0#i,imgPoints#1#i);
lineTransfers13 = stp cfold for i from 0 to 2 list (transpose Rots#i) * crossProduct(imgPoints#0#i,imgPoints#2#i);
lineTransfers23 = stp cfold for i from 0 to 2 list (transpose Rots#i) * crossProduct(imgPoints#1#i,imgPoints#2#i);
lineTransferEqs = rfold(lineTransfers12,lineTransfers13,lineTransfers23);


-*
-- BAD
--minors for T1ann & T2ann
T122s=new HashTable from (
    for j from 1 to 4 list for k from j+1 to 5 list 
    {j,k} => crossProduct(T1ann_{j},T1ann_{k}))//flatten;
T133s=(for i from 0 to 3 list for j from i+1 to 4 list for k from j+1 to 5 list(
     (transpose T1ann_{i})*T122s#{j,k})
     )//flatten//flatten;
T222s=new HashTable from (
    for j from 1 to 4 list for k from j+1 to 5 list 
    {j,k} => crossProduct(T2ann_{j},T2ann_{k}))//flatten;
T233s=(for i from 0 to 3 list for j from i+1 to 4 list for k from j+1 to 5 list(
     (transpose T2ann_{i})*T222s#{j,k})
     )//flatten//flatten;
*-

-- common line constraints

--2 determinantal tangent constraints
tangentConstraintMats = for j from 0 to 1 list(
    t := imgTans#j;
    y := imgPoints#j;
    cfold {crossProduct(t#0,y#0),transpose R12 * crossProduct(t#1,y#1),transpose R13 * crossProduct(t#2,y#2)}
    );
tangentConstraints = tangentConstraintMats/stp//rfold;

depth Tanneqs--17
depth lineTransferEqs--14
depth tangentConstraints--14
F=lineTransferEqs||tangentConstraints||Tanneqs;

fabricateChicago = FF -> (
    x := random(FF^6,FF^1);
    cqs := random(FF^8,FF^1);
    ws := ((transpose cqs^{0..3})*(matrix{{1}}||x^{0..2})) ||
    ((transpose cqs^{4..7})*(matrix{{1}}||x^{3..5}));
    P1 := id_(FF^3)|matrix{{0_FF},{0},{0}};
    P2 := cay2R(x^{0..2},Normalized=>true)|random(FF^3,FF^1);
    P3 := cay2R(x^{3..5},Normalized=>true)|random(FF^3,FF^1);
    X := random(FF^3,FF^3)||matrix{{1,1,1}};    
    Y := random(FF^3,FF^2)||matrix{{0,0}};            
    ys := {P1,P2,P3}/(m->rfold for i from 0 to 2 list m*X_{i}) // rfold;
    ts := {P1,P2,P3}/(m->rfold for i from 0 to 1 list m*Y_{i}) // rfold;
    p := ys||ts;
    (p, x)
    )
filterEval = (p,x) -> (
    -- false iff residual small
    resid := norm evaluate(F,x||p);
--    << "residual: " << resid << endl;
    (resid > 1e-4)
    )

gammify = p -> (
    ytGammas := for i from 1 to 15 list random CC;
    diag0 := flatten(ytGammas/(g->{g,g,g}));
--    diag1 := toList(4:1_CC); -- q chart, cam 2, gamma
--    diag2 := toList(4:1_CC); -- q chart, cam 3, gamma
    p' := diagonalMatrix(diag0)*p;--|diag1|diag2)*p;
    p'
    )

-- tests and setup
(p,x)=fabricateChicago CC
norm evaluate(F,x||p)
elapsedTime J=diff(gateMatrix {cameraVars},F);
J0 = (evaluate(J,x||p))
pivs=rowSelector J0
S=first SVD J0^pivs
log10((max S)/(min S))--looks good!
F'=F^pivs;
elapsedTime PH = parametricSegmentHomotopy(F', cameraVars, dataParams);


end
restart
needs "chicago.m2"

--both false?
filterEval(p,x)
filterEval(gammify p,x)

-- ad initio ~ 35s
elapsedTime (V,np)= monodromySolve(PH, 
    point p, {point x},Verbose=>true,NumberOfNodes=>4,
    TargetSolutionCount=>312,SelectEdgeAndDirection=>selectBestEdgeAndDirection,
    Potential=>potentialE, Randomizer=>gammify)

p1=matrix V.BasePoint
startSols=points V.PartialSols;
p2=random(CC^1,CC^45)
P01 = specialize(PH,
    transpose(p1|p2))
elapsedTime trackHomotopy(P01,startSols);-- ~4s

-- quality check
sols = solutionsWithMultiplicity points V.PartialSols;
L = (sols/(x -> (
	o := (transpose matrix x) ||
	   (transpose matrix V.BasePoint);
	S := first SVD evaluate(J,o);
	log10((max S)/(min S))
	)));
summary L -- better condition number stats than previous formulation

-- write to file
writeStartSys(V.BasePoint, sols, Filename => "startSys")

-- package debugging
uninstallPackage "MonodromySolver"
restart
installPackage "MonodromySolver"

uninstallPackage "SLPexpressions"
restart
installPackage "SLPexpressions"
