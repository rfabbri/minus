if instance(RERUNMONODROMY,Symbol) then RERUNMONODROMY = false 

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

-- constraints from transferring points to first image plane
pts12=for j from 0 to 2 list R12 * imgPoints#j#0
pts13=for j from 0 to 2 list R13 * imgPoints#j#0

-- matrices annihilating T1 and T2 , T1-T2 in kernel
T1ann = for i from 0 to 2 list crossProduct(imgPoints#i#1,pts12#i);
T2ann = for i from 0 to 2 list crossProduct(imgPoints#i#2,pts13#i);
Tanneqs = {T1ann,T2ann}/cfold/stp//rfold;

--BEGIN crosstalks...
-- T12,T13 up to normalization
T12 = crossProduct(T1ann#0,T1ann#1)
T13 = crossProduct(T2ann#0,T2ann#1)

-- components of projective translations
T121=T12_(0,0)
T123=T12_(2,0)
T131=T13_(0,0)
T133=T13_(2,0)

-- expressions appearing in A,B,C
topLA = T131*(imgPoints#0#2)_(2,0)-T133*(imgPoints#0#2)_(0,0)
topLB = T131*(imgPoints#1#2)_(2,0)-T133*(imgPoints#1#2)_(0,0)
topLC = T131*(imgPoints#2#2)_(2,0)-T133*(imgPoints#2#2)_(0,0)

bottomLA = T121*(imgPoints#0#1)_(2,0)-T123*(imgPoints#0#1)_(0,0)
bottomLB = T121*(imgPoints#1#1)_(2,0)-T123*(imgPoints#1#1)_(0,0)
bottomLC = T121*(imgPoints#2#1)_(2,0)-T123*(imgPoints#2#1)_(0,0)

topRA = (pts12#0)_(2,0)*(imgPoints#0#1_(0,0))-(pts12#0)_(0,0)*(imgPoints#0#1_(2,0))
topRB = (pts12#1)_(2,0)*(imgPoints#1#1_(0,0))-(pts12#1)_(0,0)*(imgPoints#1#1_(2,0))
topRC = (pts12#2)_(2,0)*(imgPoints#2#1_(0,0))-(pts12#2)_(0,0)*(imgPoints#2#1_(2,0))

bottomRA = (pts13#0)_(2,0)*(imgPoints#0#2_(0,0))-(pts13#0)_(0,0)*(imgPoints#0#2_(2,0))
bottomRB = (pts13#1)_(2,0)*(imgPoints#1#2_(0,0))-(pts13#1)_(0,0)*(imgPoints#1#2_(2,0))
bottomRC = (pts13#2)_(2,0)*(imgPoints#2#2_(0,0))-(pts13#2)_(0,0)*(imgPoints#2#2_(2,0))

AB=topLA*topRA*bottomLB*bottomRB-topLB*topRB*bottomLA*bottomRA
AC=topLA*topRA*bottomLC*bottomRC-topLC*topRC*bottomLA*bottomRA
crossTalks = gateMatrix{{AB},{AC}}
--expressions for R_{1j}*gamma_k^

-- END crosstalks

-- 3 common line constraints
S=subsets(3,2)
lineTransfers12 = stp cfold for i from 0 to 2 list (transpose Rots#i) * crossProduct(imgPoints#0#i,imgPoints#1#i);
lineTransfers13 = stp cfold for i from 0 to 2 list (transpose Rots#i) * crossProduct(imgPoints#0#i,imgPoints#2#i);
lineTransfers23 = stp cfold for i from 0 to 2 list (transpose Rots#i) * crossProduct(imgPoints#1#i,imgPoints#2#i);
lineTransferEqs = rfold(lineTransfers12,lineTransfers13,lineTransfers23);

--2 determinantal tangent constraints
tangentConstraintMats = for j from 0 to 1 list(
    t := imgTans#j;
    y := imgPoints#j;
    cfold {crossProduct(t#0,y#0),transpose R12 * crossProduct(t#1,y#1),transpose R13 * crossProduct(t#2,y#2)}
    );
tangentConstraints = tangentConstraintMats/stp//rfold;

-- circuit depth checks
depth Tanneqs--14
depth lineTransferEqs--11
depth tangentConstraints--11
depth crossTalks -- 20

F=lineTransferEqs||tangentConstraints||Tanneqs||crossTalks;

-- FUNCTIONS
--todo: not all of these are chicago-specific

-- fabricate a parameter point and solution over FF=CC or RR
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

-- filter for detecting path jumps off dominant component
filterEval = (p,x) -> (
    -- false iff residual small
    resid := norm evaluate(F,x||p);
--    << "residual: " << resid << endl;
    (resid > 1e-5)
    )

-- parameter randomizer
gammify = p -> (
    ytGammas := for i from 1 to 15 list random CC;
    diag0 := flatten(ytGammas/(g->{g,g,g}));
--    diag1 := toList(4:1_CC); -- q chart, cam 2, gamma
--    diag2 := toList(4:1_CC); -- q chart, cam 3, gamma
    p' := diagonalMatrix(diag0)*p;--|diag1|diag2)*p;
    p'
    )

-- blackbox M2 solver
solveChicagoM2 = method(Options=>{Gammify=>true})
solveChicagoM2 (Matrix, Matrix, List) := o -> (p0, p1, startSols) -> (
    P0 := if (numcols p0 == 1) then p0 else transpose p0;
    P1 := if (numcols p0 == 1) then p0 else transpose p0;
    if (o.Gammify) then (P0,P1) = (gammify P0, gammify P1);
    P01 := (gammify P0)||(gammify P1);
    H01 := specialize(PH, P01);
    trackHomotopy(H01, startSols)
    )

-- TESTS and EQUATION / AB INITIO SETUP
(p,x)=fabricateChicago CC
end
restart
needs "chicago.m2"
evaluate(F,x||p)
elapsedTime J=diff(gateMatrix {cameraVars},F);
J0 = (evaluate(J,x||p))
pivs=rowSelector J0
S=first SVD J0^pivs
log10((max S)/(min S))--looks good!
F'=F^pivs;
elapsedTime PH = parametricSegmentHomotopy(F', cameraVars, dataParams);
evalInputMatrix = gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters}
--both false?
filterEval(p,x)
filterEval(gammify p,x)

-- ad initio ~ 35s
if RERUNMONODROMY then elapsedTime (V,np)= monodromySolve(PH, 
    point p, {point x},Verbose=>true,NumberOfNodes=>4,
    TargetSolutionCount=>312,SelectEdgeAndDirection=>selectBestEdgeAndDirection,
    Potential=>potentialE, Randomizer=>gammify, FilterCondition=>filterEval)

end--

restart
RERUNMONODROMY=false
setRandomSeed 53817
needs "chicago.m2"

-- quality check
sols = solutionsWithMultiplicity points V.PartialSols;
L = (sols/(x -> (
	o := (transpose matrix x) ||
	   (transpose matrix V.BasePoint);
	S := first SVD evaluate(J,o);
	log10((max S)/(min S))
	)));
summary L -- better condition number stats than previous formulation
#sols

-- sort by condition number
S=(sort apply(L,sols,(l,s)->(l,s)))/last;

-- sample run
p1=matrix V.BasePoint
startSols=points V.PartialSols;
p2=random(CC^1,CC^45)
elapsedTime targetSols  = solveChicagoM2(p1,p2,S);-- ~4s
tally targetSols/status

-- write to file
writeStartSys(V.BasePoint, S, Filename => "startSys")

-- package debugging
uninstallPackage "MonodromySolver"
restart
installPackage "MonodromySolver"

uninstallPackage "SLPexpressions"
restart
installPackage "SLPexpressions"
