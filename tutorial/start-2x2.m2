-- START-2x2                    Homotopy Continuation tutorial                  START-2x2
--
-- NAME
--      Fast HC Code Tutorial - exactly how to craft your fast HC
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
--      - Noob and Pro marked bellow mean Example (simple) and Pro (fast)
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
restart -- only useful for debugging
needsPackage "SLPexpressions"
needsPackage "MonodromySolver"
needs "MinusUtility.m2"

-- Noob---------------------------------------------------------------------------
-- Pro doesn't use declareVariable
variables = declareVariable \ {x,y}
params = declareVariable \ {a,b,c,d,e,f}

-- Noob and Pro ------------------------------------------------------------------
-- gateSystem exists only in M2 v 1.14
GS = gateSystem(
    matrix{params},
    matrix{variables},
    transpose matrix{
	{a*(x^2+y^2)+b*x+c,
	 d*x+e*y+f}
	}
    )
-------------------------------------------------------------------------------

-- Pro 
-- random seeds for reproducible runs
-- setRandomSeed H#CLBlocks

-- Noob and Pro ------------------------------------------------------------------
PH = parametricSegmentHomotopy GS

-- Pro ---------------------------------------------------------------------------

-- Set Monodromy tracker options here
-- 
-- This is for both monodromy and online tracking. You may want to use different
-- null indicates default value 
-- 
-- TODO: show how to get default values here
-- scan({CorrectorTolerance=>1e-4,
-- 	EndZoneFactor=>2e-1,
-- 	InfinityThreshold => 1e6, 
-- 	maxCorrSteps => 3, 
-- 	noOutput => true,
-- 	numberSuccessesBeforeIncrease => 2,
-- 	Precision => null,
-- 	Predictor => RungeKutta4,
-- 	stepIncreaseFactor => 2,
-- 	tStep => 5e-2,
-- 	tStepMin => 1e-5
-- 	}, 
--     opt -> setDefault(opt)) -- setDefault is in MonodromySolver
-- 
-- setDefault(CorrectorTolerance=>1e-8)

-- Pro Gammify-style Randomization -------------------------------------------------


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
symbols = flatten entries vars GS
h=cCode(
    transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"Ht"),
    gateMatrix{symbols|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters}
    )

-- HxH
h=cCode(transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"H"),gateMatrix{symbols|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters})

-- Maybe useful
-- cCode PH

-- Noob 1 --------------------------------------------------------------------------
-- monodromy needs an initial pair of parameter, solution
-- the command below won't work for most use cases
-- except for e.g. linear in parameters
-- Here you can write a random generator from parameters 
(p0, x0) = createSeedPair GS    -- 
print(evaluate(GS,p0,x0));

-- Pro 1 -----------------------------------------------------
--
-- YOU
--    - write a specialized function to sample from the parameter space
--    - TODO: write an example on how to do this for Noob
--        - Always base your final solver on chicago.m2 technique, currently
--          the fastest and proven. Example:
--
-- gold standard Pro example:
--
-- fabricateChicago = F -> ( 
--     D := (5,0,{{0,1,3},{0,2,4},{1,2}});
--     (P, L ) := fabricatePair(D, F, nparams);
--     P1 := id_(CC^3)|matrix{{0},{0},{0}};
--     P2 := (Q2R take(P,{0,3})) | 
--          transpose matrix({take(P,{8,10})});
--     P3 := Q2R take(P,{4,8})|transpose matrix{take(P,{11,13})};
--     projs := {P1,P2,P3};
--     allLines := L/last;
--     independentLineIndices := fold((last D)/(p -> set {p#0,p#1}),(a,b)->a+b);
--     dependentLineIndices := set(0..D#0-1)-independentLineIndices;
--     p0 := flatten entries fold(
-- 	    allLines/(m->m^(toList independentLineIndices))
-- 	    ,(a,b)->a|b);    
--     p1 = flatten flatten(
-- 	allLines/(m -> (
-- 	    n :=numericalKernel(transpose m^{0,1,3}, kTol);
-- 	    entries((1/n_(2,0))*n^{0,1})))
--     );
--     p2 = flatten flatten(
--     	    allLines/(m -> (
-- 	    n :=numericalKernel(transpose m^{0,2,4}, kTol);
-- 	    entries((1/n_(2,0))*n^{0,1})))
--     );
--     x := transpose matrix{take(P,nvars)};
--     p3 := flatten entries randKernel(transpose(x^{8..13}||matrix{{1_CC}}), F);
--     p4 := flatten entries randKernel(transpose(x^{0..3}||matrix{{1_CC}}), F);
--     p5 := flatten entries randKernel(transpose(x^{4..7}||matrix{{1_CC}}), F);
--     p := transpose matrix{p0|p1|p2|p3|p4|p5};
--     (p, x)
--     )

-- (p0, x0) = fabricateChicago(CC)   -- CC just means Complexes

-- filter path jumps during monodromy
-- filterEval = (p,x) -> (
--     -- false iff residual small
--     resid := norm evaluate(F,x||p);
-- --    << "residual: " << resid << endl;
--     (resid > 1e-4)
--     )
-- filterEval(p0,x0)  ----------------------------------------------------------

-- Noob 2
(V,np) = monodromySolve(GS,p0,{x0},Verbose=>true,NumberOfNodes=>5) -- first try with default NumberOfNodes

-- Pro 2 -------------------------------------------------------------------
-- elapsedTime (V,np)= monodromySolve(PH, 
--     point p0, {point x0},Verbose=>true,
--     FilterCondition=>filterEval,TargetSolutionCount=>312,SelectEdgeAndDirection=>selectBestEdgeAndDirection,
--     Potential=>potentialE, Randomizer=>gammify)
-------------------------------------------------------------------------------


-- Noob 3 -------------------------------------------------------------------
-- Pro
-- May want to inspect basepoint used for the initial solutions
-- V.BasePoint;
-- corresponding solutions
points V.PartialSols 
sols = solutionsWithMultiplicity points V.PartialSols; -- solutionsWithMultiplicity only a safe option, can try remove
                                                         -- it might speed up if removed

-- Pro 3 --------------------------------------------------------------
-- quality check
-- L = (sols/(x -> (
-- 	o := (transpose matrix x) ||
-- 	   (transpose matrix V.BasePoint);
-- 	S := first SVD evaluate(J',o);
-- 	log10((max S)/(min S))
-- 	)));
-- summary L -------------------------------------------------------------------

writeStartSys(V.BasePoint, sols, Filename => "startSys")
