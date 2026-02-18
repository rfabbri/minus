# CVPR 2024 fast Homotopy Continuation tutorial

A new polynomial system can be solved in a number of ways within our framework,
ranging from an initial Macaulay prototype, to a more advanced "Pro" Macaulay
prototype to a fast C++ version. From the initial prototype to the final C++
code, there are a number of intermediate steps to follow. 

In this tutorial we use a toy example codenamed `linecircle` which you can use
as a template for your own solver. You should codename your own
problem/formulation into a string, say, `problem1`, and follow what is done for
linecircle to your own problem.

## Running the example Macaulay template

The example is called linecircle, stored in its own folder linecircle/
It simply computes the intersection of a line and a circle.

The equations are stored in: `equations-linecircle.m2`

### Solve the start solutions in Macaulay

In a terminal, start Macaulay2, by typing:
```bash
m2
```
inside `linecircle/`. Then

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

In the script `end-linecircle.m2` variable `p1` will hold the desired target parameters.
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
uncommment some of these and experiment with your code.

## Building a fast C++ solver

The Macaulay code can be used to generate a fast C++ solver. The `starSols`,
and `HxHt` and `HxH` files used in the `start-linecircle.m2` can be used to
write a fast solver within the Minus C++ framework. After the necessary steps,
the final solver executables are available as `cmd/minus-linecircle` in the
binary folder. Since these steps are more involved, we refer to the next
sections.


<!------------------------------------------------------------------------------>

## Solving your own system

Give your new problem a name, say `problem1`. 

### 1 build generic Macaulay2 scrips for your polynomial system
You should first
copy the corresponding files in `minus/tutorial/*linecircle` by substituting the
string `linecircle` to `problem1` in the filenames and their contents.

This should get you running the basic scripted solver in Macaulay.

### 2 fine-tune your Macaulay2 scrips
You can then follow the steps to produce a more optimized solver still in
Macaulay as documented above for `linecircle`. 

### 3 C++ fast solver
There are steps to now generate your own C++ optimized solver within the
Minus C++ framework. We are in the process of releasing the precise steps.
An idea can be had in the toplevel `minus/README.md` file with accompanying
videos. The final steps will be released soon, but you can start with.

#### 3.1 define a tag ID for your polynomial system
The tag encodes a specific problem, and a specific formulation and optimization
Example: linecircle2a  is <problem><system size><formulation tag>, where

<problem> = a non-numeric short string naming your problem (eg, linecircle)
<system size> = number of variables or equations of the square system (eg, 2 for 2x2)
<formulation tag> = a letter (or any string) to identify the specific
                    formulation. A formulation is a special set of equations and
                    solver optimizations. Eg, linecircle2a means formulation
                    'a', which represents a circle as a quadratic cartesian equation,
                    and a line as a homogeneous equation of degree 1, with
                    specific code optimizations. You could technically break
                    this tag into more informative components if you are
                    experimenting with diverse equations, representations,
                    coordinates, and source code optimizations.

Examples: chicago14a means a trifocal problem codenamed chicago, 14x14,
formulation and specific solver 'a' which corresponds to a certain
Cayley/homogeneous quaternion formulation. 


#### 3.2 fill in minus.h header with the polynomial system ID tag

Everywhere the sample tag linecircle2a appears, adapt to your system ID tag.

##### 3.2.1 Fill parameters.h (included by minus.h) with the include for your polynomial system

##### 3.2.2 Fill problem-defs (included by minus.h) with a specialized struct for your problem


#### 3.3 copy linecircle2a.h to yor system_id_tag.h and adapt it

Fill in this file with the number of solutions, etc.

#### 3.4 copy linecircle2a.hxx to yor system_id_tag.h and adapt it

The file linecircle2a.hxx has Noob/Pro markings on it to help adapt this file

#### 3.5 copy Macaulay2-generated C++ evaluators to Minus

Copy linecircle-Hxt.hxx to system_id_tag-Hxt.hxx and adapt, using
Macaulay2-generated c++ code

Copy linecircle-HxH.hxx to system_id_tag-HxH.hxx and adapt, using
Macaulay2-generated c++ code

##### 3.5.1 Split x to separate params
Automatically-generated Macaulay2 evaluators encode the variables, parameters,
and homotopy parameter t all in one big vector x. For performance reasons,
MINUS splits the parameters out into a separate vector. For now,
you must rename this by hand, e.g.:

const C<F> &X0 = x[0]; // x
const C<F> &X1 = x[1]; // y
const C<F> &X2 = x[2]; // t
const C<F> &X3 = x[3]; // HERE: replace x[3] with params [0]
                       //       replace x[4] with params [1]
                       // ....
                       // and so on until the last x.

TODO We could provide our own version of cCode for that.

#### 3.6 write basic functions used in I/O
- Adapt linecircle2a-io.h to your problem
- Fill in default data in linecircle2a-default-data.h
    - Copy startsys to an auxiliary file, say startsys.c
    - Open in vim, run :so minus/scripts/parse-m2string2.vim
    - Copy these to your linecircle2a-default-data.h at the appropriate places
    - Copy ground-truth solutions from end-linecirle.m2
    

#### 3.7 write the function solve in linecircle2a.hxx
- You need to set the start and target parameters, optionally converting from
  user data (e.g., point correspondences) to these parameters used in the
  equatons.
