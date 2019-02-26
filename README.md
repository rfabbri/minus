<!--<img src="vxd-logo.gif" width="215" height="84" /> -->


## MINUS: MInimial problem NUmerical continuation Solver
This package originated in solving medium-sized  (degree > 100) square problems
(ie, exactly-determined, 0-dimensional), notably in computer vision where
trifocal minimal problems from points and lines are of importance (as in
curve-based structure from motion, where lines are tangents to curves).
As of this date, such problems are too high degree to be solved symbolically.
 
For more details, see the
[website](http://multiview-3d-drawings.sourceforge.net)

## Usage
For use in your program, we provide a header-only library.
Simply do:
```
#include <minus>
```

Examples are avaialble in the tests/ subfolder
In your program, you can then use

```
  ptrack<14>(&VNAG_DEFAULT, start_sols, params, solutions);

```
To solve a 14x14 precompiled trifocal problem from lines on points ("Chicago").
You do need to know the size of the system in advance, for efficiency reasons.
This is not dynamic code, so no allocations are performed.


## Authors

Adapted and improved over Anton Leykin's and Timothy Duff's codes from Macaulay2 by
Ricardo Fabbri. The core code grew out of Macaulay2/e/NAG.cpp.

## Acknowledgements
This grew as part of ICERM's 2018 Nonlinear Algebra Program (Computer Vision
Working group) and 2019 Algebraic Vision research cluster. Both co-organized by
the authors.
