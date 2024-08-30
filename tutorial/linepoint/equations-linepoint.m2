variables = declareVariable \ {a,b}
params = declareVariable \ {x,y}

-- gateSystem exists only in M2 >= 1v 1.14
GS = gateSystem(
    matrix{params},
    matrix{variables},
    transpose matrix{
	{a*x+b*y,        -- problem at: x=1e-n, y=1 --> a*1e-n + b = 0 and b - = 1
                         --  a*1e-n + 1 = 0 --> a = -1/1e-n = -1*10^n
                         -- Fazer n -> inf , a um certo ponto o algoritmo vai ou
                         -- ficar muito lento ou falhar
                         -- por varios motivos
                         -- 1) precisao numerica
                         -- 2) o caminho geometrico no espaco das solucoes fica
                         -- infinitamente gigante na medida em que n->inf
                         -- OBS: isso independentemente do Gammify! A nao ser
                         -- complexificacao do problema "contorne" a
                         -- singularidade
	 b-1}
	}
    )

-- Pro 
-- random seeds for reproducible runs
-- setRandomSeed H#CLBlocks

PH = parametricSegmentHomotopy GS

