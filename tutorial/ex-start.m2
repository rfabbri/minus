-- EX-START                        Homotopy Continuation Tutorial                             EX-START
--
-- NAME
--      Fast HC Code Inteligencer - exactly how to craft your fast HC
--      Module: Start solution generator (a.k.a. training stage, run once offline)
--
-- DESCRIPTION
--      This script shows how to craft a fast solver for your problem
--      It is broken into Generic and Pro blocks
--
--      The design-goal of the file is to closely match what will be in the fast C++
--      solver, so you can start with a generic and incrementally craft and test out the pro features
--      
--      The file is not based on generic tutorial scripts. Rather, it is rather technical
--      based on techniques that originally solved a very hard problem and for the
--      first time ever before, trifical pose from points and tangents Fabbri, etal, CVPR'10.
--      It is important to keep in mind that, while there are novelties in the -- associated m2 
--      scripts of PLMP, this version of them is the definitive account
--      of what matterscarefully for writing a fast C++ solver
--
-- LEGEND
--      - Ex and Pro marked bellow mean Example (simple) and Pro (fast)
--      - YOU shows where can you make it faster for your problem

-- OTHERS
--      - Think of this as a prompt precise enough for LLM to help generate a fast solver
--
-- TODO
--      - update with PLMP versions of common, etc, making _sure_ it
--        doesn't impact performance
--
-- AUTHORS 
--      Ricardo Fabbri (C++ and M2) and Timothy Duff (M2)
--

-- code for generating various evaluators 
restart
needsPackage "SLPexpressions"
needsPackage "MonodromySolver"
needs "MinusUtility.m2" -- put cCode extesion here
-- "gateSystem" exists only in M2 v 1.14

-- Ex---------------------------------------------------------------------------
-- Pro doesn't use declareVariable
variables = declareVariable \ {x,y}
params = declareVariable \ {a,b,c,d,e,f}

-- Ex and Pro ------------------------------------------------------------------
GS = gateSystem(
    matrix{params},
    matrix{variables},
    transpose matrix{
	{a*(x^2+y^2)+b*x+c,
	    d*x+e*y+f}
	}
    )
-------------------------------------------------------------------------------

-- Ex and Pro ------------------------------------------------------------------
cameraVars = flatten entries vars GS
PH = parametricSegmentHomotopy GS

-- Pro only --------------------------------------------------------------------

-- YOU set tracker options here
-- null indicates default value 
-- TODO: show how to get default values here
scan({CorrectorTolerance=>1e-4,
	EndZoneFactor=>2e-1,
	InfinityThreshold => 1e6, 
	maxCorrSteps => 3, 
	noOutput => true,
	numberSuccessesBeforeIncrease => 2,
	Precision => null,
	Predictor => RungeKutta4,
	stepIncreaseFactor => 2,
	tStep => 5e-2,
	tStepMin => 1e-5
	}, 
    opt -> setDefault(opt)) -- setDefault is in MonodromySolver

setDefault(CorrectorTolerance=>1e-8)

-- Pro Gammify-style Randomization ---------------------------------------------------


-- YOU: adapt for your problem
-- Randomization for parameter point p
-- 
-- gammify = p -> (
--     gammas := for i from 1 to 9 list random CC;
--     diag0 := flatten(gammas/(g->{g,g,g}));
--     tripleIntersections := {{0,3,9},{0+1,3+1,9+1},{0+2,3+2,9+2},
-- 	{0,6,12},{0+1,6+1,12+1},{0+2,6+2,12+2}};
--     diag1 := flatten(tripleIntersections/(ind -> {
-- 	    conjugate gammas#(ind#0),
-- 	    conjugate gammas#(ind#1)
-- 	    }
-- 	));
--     diag2 := toList(7:random(CC)); -- t chart gamma
--     diag3 := toList(5:random(CC)); -- q chart, cam 2, gamma
--     diag4 := toList(5:random(CC)); -- q chart, cam 3, gamma
--     p' := diagonalMatrix(diag0|diag1|diag2|diag3|diag4)*p;
--     p'
--     )

-------------------------------------------------------------------------------

-- HxHt
h=cCode(
    transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"Ht"),
    gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters}
    )

-- HxH
h=cCode(transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"H"),gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters})

-- Maybe useful
-- cCode PH

-- Ex 1 --------------------------------------------------------------------------
-- monodromy needs an initial pair of parameter, solution
-- the command below won't work for most use cases
-- ie) need 
(p0, x0) = createSeedPair GS   

-- Pro 1 -----------------------------------------------------
--
-- YOU
--    - write a specialized function to sample from the parameter space
--    - TODO: write an example on how to do this for Ex
--        - Always base your final solver on chicago.m2 technique, currently
--          the fastest and proven. Example:


---- Pro ----- fabricateChicago = F -> ( -- gold standard Pro example ----------------------
---- Pro -----     D := (5,0,{{0,1,3},{0,2,4},{1,2}});
---- Pro -----     (P, L ) := fabricatePair(D, F, nparams);
---- Pro -----     P1 := id_(CC^3)|matrix{{0},{0},{0}};
---- Pro -----     P2 := (Q2R take(P,{0,3})) | 
---- Pro -----          transpose matrix({take(P,{8,10})});
---- Pro -----     P3 := Q2R take(P,{4,8})|transpose matrix{take(P,{11,13})};
---- Pro -----     projs := {P1,P2,P3};
---- Pro -----     allLines := L/last;
---- Pro -----     independentLineIndices := fold((last D)/(p -> set {p#0,p#1}),(a,b)->a+b);
---- Pro -----     dependentLineIndices := set(0..D#0-1)-independentLineIndices;
---- Pro -----     p0 := flatten entries fold(
---- Pro ----- 	    allLines/(m->m^(toList independentLineIndices))
---- Pro ----- 	    ,(a,b)->a|b);    
---- Pro -----     p1 = flatten flatten(
---- Pro ----- 	allLines/(m -> (
---- Pro ----- 	    n :=numericalKernel(transpose m^{0,1,3}, kTol);
---- Pro ----- 	    entries((1/n_(2,0))*n^{0,1})))
---- Pro -----     );
---- Pro -----     p2 = flatten flatten(
---- Pro -----     	    allLines/(m -> (
---- Pro ----- 	    n :=numericalKernel(transpose m^{0,2,4}, kTol);
---- Pro ----- 	    entries((1/n_(2,0))*n^{0,1})))
---- Pro -----     );
---- Pro -----     x := transpose matrix{take(P,nvars)};
---- Pro -----     p3 := flatten entries randKernel(transpose(x^{8..13}||matrix{{1_CC}}), F);
---- Pro -----     p4 := flatten entries randKernel(transpose(x^{0..3}||matrix{{1_CC}}), F);
---- Pro -----     p5 := flatten entries randKernel(transpose(x^{4..7}||matrix{{1_CC}}), F);
---- Pro -----     p := transpose matrix{p0|p1|p2|p3|p4|p5};
---- Pro -----     (p, x)
---- Pro -----     )

-- (p0, x0) = fabricateChicago(CC)   -- CC just means Complexes

-- filter path jumps during monodromy
-- filterEval = (p,x) -> (
--     -- false iff residual small
--     resid := norm evaluate(F,x||p);
-- --    << "residual: " << resid << endl;
--     (resid > 1e-4)
--     )
-- filterEval(p0,x0)  ----------------------------------------------------------

-- Ex 2
(V,np) = monodromySolve(GS,p0,{x0},Verbose=>true,NumberOfNodes=>5)

-- Pro 2 -------------------------------------------------------------------
-- elapsedTime (V,np)= monodromySolve(PH, 
--     point p0, {point x0},Verbose=>true,
--     FilterCondition=>filterEval,TargetSolutionCount=>312,SelectEdgeAndDirection=>selectBestEdgeAndDirection,
--     Potential=>potentialE, Randomizer=>gammify)
-------------------------------------------------------------------------------


-- Ex 3 -------------------------------------------------------------------
-- parameter point
V.BasePoint
-- corresponding solutions
points V.PartialSols -- Pro 3 --------------------------------------------------------------
-- quality check
-- sols = solutionsWithMultiplicity points V.PartialSols;
-- L = (sols/(x -> (
-- 	o := (transpose matrix x) ||
-- 	   (transpose matrix V.BasePoint);
-- 	S := first SVD evaluate(J',o);
-- 	log10((max S)/(min S))
-- 	)));
-- summary L -------------------------------------------------------------------
-- writeStartSys(V.BasePoint, sols, Filename => "startSys" | (toString CLBlocks#0)|(toString CLBlocks#1))
