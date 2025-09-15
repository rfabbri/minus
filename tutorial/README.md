# CVPR 204 fast Homotopy Continuation tutorial

A new polynomial systems can be solved in a number of ways within our framework,
ranging from an initial Macaulay prototype, to a more advanced Macaulay
prototype to a fast C++ version. From the initial prototype to the final C++
code, there are a number of intermediate steps to follow. 

In this tutorial ywe use a toy example codenamed `linecircle` which you can use
as a template for your own solver. You should codename your own
problem/formulation into a string, say, `problem1`, and follow what is done for
linecircle to your own problem.

## Running the example Macaulay template

The example is called linecircle and computes the intersection of a line and a
circle.

The equations are stored in: `equations-linecircle.m2`

### Solve the start solutions in Macaulay

In a terminal, start Macaulay2, typing:
```bash
m2
```

Then

```
load("start-linecircle.m2")
```

This will run monodromy to compute the solutions in your file.
It will write out the following files into the current folder:

- `startSys` : the start solutions and corresponding base point (parameters).
  This is a text file in Macaulay format. It can be used both for your fast C++ solver,
  and for prototyping the fast online solution in the next step.
  
- `HxHt.cxx`, `HxH.cxx`: C++ evaluators for writing your optimized solver later on


### Solver for any target system in Macaulay

In file `end-linecircle.m2` variable `p1` will hold the desired target parameters.
For real problems, p1 will be constructed from problem data, say, image
correspondences.

In Macaulay, type:
```
load("end-linecircle.m2")
```

The solutions will be output to the screen. 

For chicago trifocal problem, these solutions take about 1min to compute.
This Macaulay solver will run a generic C++ homotopy continuation code under the
hood. This solver is code-matched to Minus fast C++ solver, with the
continuation parameters in Minus the same as the ones in Macaulay2, with some
additional parameters in Minus for performance.

## Optimizing the example Macaulay template

Inside the files `*-linecircle.m2` there are commented sections "Pro" after each
Macaulay command detaling what can be done to optimize the code. Try to
uncommment some of these and experiment with your code

## Building a fast C++ solver

The Macaulay code can be used to generate a fast C++ solver. The `starSols`,
and `HxHt` and `HxH` files used in the `start-linecircle.m2` can be used to
write a fast solver within the Minus C++ framework. After the necessary steps,
the final solver executables are available as `cmd/minus-linecircle` in the
binary folder. Since these steps are more involved, we refer to the next
sections.


<!------------------------------------------------------------------------------>

## Solving your own system

Give your new problem a name, say `problem1`, you should first
copy the corresponding files in `minus/tutorial/*linecircle` by substituting the
string `linecircle` to `problem1` in the filenames and their contents.

This will get you running the basic scripted solver in Macaulay.
You can then follow the steps to produce a more optiized solver still in
Macaulay as documented above for `linecircle`. 

There are steps to now generate your own C++ optimized solver within the
Minus C++ framework. We are in the process of releasing the precise steps.
An idea can be had in the toploevel `minus/README.md` file with accompanying
videos. The final steps will be released soon.
