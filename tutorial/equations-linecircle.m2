variables = declareVariable \ {x,y}
params = declareVariable \ {a,b,c,d,e,f}

-- gateSystem exists only in M2 >= 1v 1.14
GS = gateSystem(
    matrix{params},
    matrix{variables},
    transpose matrix{
	{a*(x^2+y^2)+b*x+c,
	 d*x+e*y+f}
	}
    )
