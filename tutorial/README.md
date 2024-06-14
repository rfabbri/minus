# TODO
## 4Tim

- createSeedPair vs fabricateChicago exactly?

-- Chicago.m2:125 what is this, is it ever used?
C = { 
    makeSLProgram(matrix{cameraParams},P1),
    makeSLProgram(matrix{cameraParams},P2),
    makeSLProgram(matrix{cameraParams},P3)
  }

## Evaluate to check equation
- How can we evaluate the the equations on the system?

    - evaluate(F,x0||p0) only for F gatematrix, 
    - but we ony have gateSystem in Ex, not gatematrix. 
    - gateSystem is >= M2 1.14, whats the corresponding eval fn?

Usually one creates gateSystem from Gatematrix:
GS = gateSystem(GS)


(p,x)=fabricateChicago CC
evaluate(F,x||p)

(p0,x0) = fabricateChicago(CC)
evaluate(F,x0||p0)
    
