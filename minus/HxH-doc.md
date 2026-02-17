# Generic eval HxH documentation: ------------------------------------------------------

# SYNOPSIS
   Evaluates variable partial derivatives Hx and the homotopy equations H
   themselves at the same time, reusing expressions.

# DESCRIPTION
   Maps from a multivariate polynomial whose variables and t are stored in 
   `x = [x_1,..,x_NVE, t]` and parameters are `params`, to `y = [Hx|H]`, where Hx is the Jacobian
   matrix of the homotopy equations with respect to the x variables `[x_1,...,x_NVE]`, 
   and H evaluates the homotopy equations themselves with parameters params at x (and t),
   `H = [f_1(x,t,params),...,f_NVE(x,t,params)]`. Some just put everything in one
   x (because all parameters are really all like variables), then we write
   `H(x) = [f_1(x),...,f_NVE(x)]`.

# INPUT 
   x, params: see documentation in the corresponding Hxt file and Hxt-doc.md

# OUTPUT 
   y: NVExNVEPLUS1 matrix [Hx|H] as a col-major 1D vector

```Macaulay2
    cCode(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"H",gateMatrix{cameraVars})
    -- we may transpose to get the right C/C++ order Minus expects
```
