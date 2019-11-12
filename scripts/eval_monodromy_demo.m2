-- code for generating various evaluators 
restart
needsPackage "SLPexpressions"
needsPackage "MonodromySolver"
-- "gateSystem" exists only in M2 v 1.14
variables = declareVariable \ {x,y}
params = declareVariable \ {a,b,c,d,e,f}
GS = gateSystem(
    matrix{params},
    matrix{variables},
    transpose matrix{
	{a*(x^2+y^2)+b*x+c,
	    d*x+e*y+f}
	}
    )

cameraVars = flatten entries vars GS
PH = parametricSegmentHomotopy GS

-- HxHt
h=cCode(
    transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"Ht"),
    gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters}
    )

-- HxH
h=cCode(transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"H"),gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters})

-- monodromy needs an initial pair of parameter, solution
-- the command below won't work for most use cases
-- ie) need 
(p0, x0) = createSeedPair GS
(V,np) = monodromySolve(GS,p0,{x0},Verbose=>true,NumberOfNodes=>5)
-- parameter point
V.BasePoint
-- corresponding solutions
points V.PartialSols
