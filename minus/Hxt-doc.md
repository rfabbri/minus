# Generic eval Hxt documentation: ------------------------------------------------------

## SYNOPSIS
   Evaluates variable and time partial derivatives Hx and Ht of the Homotopy
   equations simultaneously, reusing expressions.

## DESCRIPTION
   Maps from a multivariate polynomial whose variables and t are stored in 
   `x = [x_1,..,x_NVE, t]` and parameters are `params`, to `y = [Hx|-Ht]`, which is
   the full Jacobian of the homotopy equations with respect to the variables
   and time (with the last column negated), where:
   
   Hx is the Jacobian matrix of the homotopy equations with respect to the x
   variables `[x_1,...,x_NVE]`, and 
   
   Ht is the 1D vector of temporal derivatives of each equation.

## INPUT 
   x, NVEPLUS1-dimensional: NVE for variables  `[x1,x2,...]`, 1 for t
   params, 2NPARAMS-dimensional: 
       NPARAMS parameters for start system
       NPARAMS for target system.
       
       These may encode 
       1) what is varied depending on t, such as polynomial coefficients
       themselves or a function of them that blends between start and end
       system as t varies, and 2) what is kept constant and is not varied,
       such as randomization of start and end systems. 

## OUTPUT 
   y: NVExNVEPLUS1 full Jacobian matrix `[Hx|-Ht]` as a col-major 1D vector
   
   The last column is negated since it is what is needed for the linear solves in MINUS.

## SEE ALSO
   see `tutorial/*/start-*.m2` or `scripts/eval_monodromy_demo.m2` to see how to
   generate these from equations in Macaulay2, e.g.:
```Macaulay2
   -- something like:
   cCode(PH.GateHomotopy#"Hx" | - PH.GateHomotopy#"Ht",gateMatrix{cameraVars})
   -- we may transpose to get the right C/C++ order Minus expects
   -- see the correct working code in tutorial/linecircle/linecircle-start.m2
```
