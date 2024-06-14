-- Ex and Pro marked bellow mean Example (simple) and Pro (fast)


-- code for generating various evaluators 
restart
needsPackage "SLPexpressions"
needsPackage "MonodromySolver"
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

-- SET TRACKER OPTIONS HERE
-- null indicates default value 
scan({CorrectorTolerance=>1e-4,
	EndZoneFactor=>2e-1,
	InfinityThreshold => 1e6, 
	maxCorrSteps => 3, 
	NoOutput => true,
	numberSuccessesBeforeIncrease => 2,
	Precision => null,
	Predictor => RungeKutta4,
	stepIncreaseFactor => 2,
	tStep => 5e-2,
	tStepMin => 1e-5
	}, 
    opt -> setDefault(opt))

setDefault(CorrectorTolerance=>1e-8)


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
-- filter path jumps during monodromy
-- filterEval = (p,x) -> (
--     -- false iff residual small
--     resid := norm evaluate(F,x||p);
-- --    << "residual: " << resid << endl;
--     (resid > 1e-4)
--     )
-- 
-- (p0, x0) = fabricateChicago(CC)
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
points V.PartialSols

-- Pro 3 --------------------------------------------------------------
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
