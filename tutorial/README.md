# CVPR 204 fast Homotopy Continuation tutorial

A new polynomial systems can be solved in a number of ways within our framework,
ranging from an initial Macaulay prototype, to a more advanced Macaulay
prototype to a fast C++ version. From the initial prototype to the final C++
code, there are a number of intermediate steps to follow. 

## Running the example Macaulay2 template

The example is called linecircle and computes the intersection of a line and a
circle.

The equations are stored in: `equations-linecircle.m2`

### Solve the start solutions

In a terminal, type:
```bash
m2
```

Then

```
load("start-linecircle.m2")
```

This will run monodromy to compute the solutions in your file.
It will write out the following files into the current folder:

- `startSys` : the start solution and corresponding base point (parameters).
  This is a text file in Macaulay format. It can be used both for your fast C++ solver,
  and for prototyping the fast online solution in the next step.
  
- `HxHt.cxx`, `HxH.cxx`: C++ evaluators for writing your optimized solver later on


## Solving your own system

Give your new problem a name, say `problem1`


### Type your equations

Copy the template over:

```
cp equations-linecircle.m2 equations-problem1.m2
```

Edit the file typing your variables, parameters and equations, foloowing the
template as an example.


### Compute start stystems

First, copy over the solution template:
```
load("start-linecircle.m2")
```

In a terminal, type:
```bash
m2
```

Then

```
load("start-linecircle.m2")
```


