<!--<img src="vxd-logo.gif" width="215" height="84" /> -->


## MINUS: MInimial problem NUmerical continuation Solver [![Build Status](https://travis-ci.org/rfabbri/minus.svg?branch=master)](https://travis-ci.org/rfabbri/minus)

This package originated in solving medium-sized  (eg, degree > 100) square problems
(ie, exactly-determined, 0-dimensional), notably in computer vision where
trifocal minimal problems from points and lines are of importance (as in
curve-based structure from motion, where lines are tangents to curves).
As of this date, such problems are too high degree to be solved symbolically,
while being low enough degree to allow an efficient implementation of a global
technique, rather than the usual local Levenberg-Marquardt of structure from
motion that requires an initialization. The solutions found using this technique
for square problems can be used to initialize Levenberg-Marquardt for an
overconstrained problem.

Minus is split into three parts:
- An efficient library for use in your C++ programs, no dependencies beyond
  standard C++ and Eigen (included)
- Simple commandline program easy to interface with other programs (eg, Matlab and Macaulay2)
- Optional: extensive tests useful for tuning the algorithm to a machine architecture and
  compiler. These are disabled by default, and requires [VXL core library](https://vxl.github.io) (*optional* - most users don't need this).
 
## Paper
The theory and practice associated to Minus is described in

"Trifocal Relative Pose from Lines at Points and its Efficient Solution", CVPR
2020, Arxiv March 23 2019 ([pdf](http://rfabbri.github.io/stuff/fabbri-kimia-etal-arxiv2019-v2.pdf)). 
For datasets and curve-based SfM code, see the [website](http://multiview-3d-drawings.sourceforge.net).

## Usage in C++ programs
For use in your program, we provide a C++ header-only library.
Simply do:
```C
#include <minus.hxx>
using namespace MiNuS;
```

And use `minus<chicago>` like so:
```C
  minus<chicago>::solve(p, tgt, solutions_cameras, &nsols_final); 
```
to solve a trifocal pose problem from lines at points ("Chicago"), using the default
formulation.  See the full example in `cmd/minus-chicago.cxx` and `tests/test-minus.cxx`.  This is
efficient *static* code, so no allocations are performed.  

The size and key parameters of the minimal problem are hardcoded as template
parameters in advance, for efficiency. `minus<>` is a shorthand for a generic
template, so you have full control to add your own compiled formulations, or
change the floating point implementation.  You can easily specify the
formulation by using, e.g., `minus<chicago14a>` instead of `minus<chicago>`,
which will use a 14x14 formulation for the Chicago problem.  Other
instances are available, e.g. `minus<chicago6a>` to solve a 6x6
formulation for the same problem instead. 

To solve another problem, say the `Cleveland` trifocal pose problem from mixed
points and lines, simply use the problem tag:
```C
  minus<cleveland>::solve(p, tgt, solutions_cameras, &nsols_final);
```

If your measurements are in pixels, use `solve_img` and pass your calibration matrix `K`:
```C
  minus<cleveland>::solve_img(K, p, tgt, solutions_cameras, &nsols_final);
```

### Further control on template parametrs

If you want to have full control on templating, say, to change from `double` to
`float`, `minus<problem>` is just a shorthand for `minus<problem,double>`.
So you can use `minus<chicago, float>`. See the section Hacking for how
to add your own problem formulation to Minus.

## Commandline programs

You will need to compile the commandline program to take advantage of
processor-specific code optimization. 

### Compiling
This requires CMake configure and make. Developer
tools are required. GCC 5 to 7 is strongly recommended.
You will have to have a recent version of Cmake installed (try it with your
system cmake, and if it doesn't work we recommend installing Cmake from git, it is
straightforward).  

```bash
cd minus
ccmake .           # press 'c' repreatedly (configure), then 'g' (generate)
make
```

The executables are located in the `bin` folder (`minus/bin`).
Each executable is optimized for a different minimal problem.
```bash
minus-chicago      # annihilates the Chicago problem!
```

Do an initial test

```bash
cd cmd
./minus-chicago -g         # -g profiles a predefined worst case, to get a time

Output:
  LOG Time of solver: xxxms
```

### Running

We will use Chicago as the basic example of a minimal problem to be solved.

The usage is as follows
```bash
minus-chicago input output 
```

Where input and output are ASCII text files. 
`input` encodes points and lines, and `output` has the solutions

If you are communicating to/from another program (Matlab or Macaulay2),
you should not use physical files, but use pipes (standard input and output).
Without any arguments, `minus` will read from stdin, and write to stdout.
That way, your script can do this:
```
cat input | minus-chicago 
# or
minus-chicago < input > output
# In Matlab, something like this:
solutions = system( pipe input to minus-chicago )  
# solutions are output of the command that goes directly into Matlab
```

By default, minus reads `input` in a raw format with least overhead,
the purpose being to use only the minus core solver from other programs,
which is fastest and has the least potential for bugs, but requires the calling
program to do a lot of pre-processing.

#### Image pixel data as input
```bash
minus-chicago -i
```
will read input in image data measured in pixels (before inverting intrinsics
K). It will read from standard input by default, which can be used with other
programs with in-memory I/O, avoiding physical files:

```bash
synthdata | minus-chicago -i              # synthdata is in minus/scripts/synthdata
```
or
```bash
minus-chicago -i  < input > output
```
The input to `minus-chicago -i` is as follows:
Input format (notation is `_view_points_coords`. any number of spaces and newlines optional. Can be in
one row or one column as well). This input format assumes tangent data for
all points, but you specify which one to use in id0 and id1 below. When
`--use_all_tangents` is passed, will try to select the better conditioned / least degenerate tangents 
 
```bash
  p000 p001
  p010 p011
  p020 p021
  
  p100 p101
  p110 p111
  p120 p121
  
  p100 p101
  p110 p111
  p120 p121
 
  t000 t001
  t010 t011
  t020 t021
  
  t100 t101
  t110 t111
  t120 t121
  
  t100 t101
  t110 t111
  t120 t121
  
  id0 id1           # id \in {0,1,2} of the point to consider the tangent
  
  K00 K01 K02       # intrinsic parameters: only these elements
   0  K11 K22
                    # GROUND TRUTH (optional) if -gt flag provided, pass the ground truth here:
  r000 r001 r002    # default camera format if synthcurves flag passed: 
  r010 r011 r012    # just like a 3x4 [R|T] but transposed to better fit row-major:
  r020 r021 r022    #         | R |
   c00  c01  c02    # P_4x3 = | - |
                    #         | C'|
  r100 r101 r102
  r110 r111 r112
  r120 r121 r122
   c10  c11  c12 
  
  r200 r201 r202
  r210 r211 r212
  r220 r221 r222
   c20  c21  c22
```

Further information is provided by typing `minus-chicago --help`.

#### Distributed parallelism for RANSAC

For now, each time you run minus-chicago inside RANSAC, you will have to call
minus again. But this is OK since you can parallelize your RANSAC using GNU
Parallel. Example:

```bash
parallel minus-chicago {1} {2} ::: "$inputfiles" ::: "$outputfiles"
```

Will identify the number of cores and threads in your CPU and distribute a copy
of minus-chicago for each input file describing a point/line configuration.

Parallel can be easily configured to run across many machines in the lab, for
instance.

(I've tried running at the cluster at Brown but it seems each Xeon is far less
strong than a i7 for a single run, though for parallel tasks it may overcome
that. I recommend running on a couple of machines you know minus runs fast)

Their format are as follows
```
input:
    Encodes input points and lines in text form.
    Represented as start-target parameter pairs as P01 in solveChicago/chicago.m2 from Tim.
    The format is close to that needed by the homotopy continuation code, to
    minimize bug in brittle C++ code.  The user can write scripts to get
    user-friendly data into this format.
    
    It looks like
    
    .391195550619826 -.00262962533857666
    .310140709227333 +.169842562835882
    -.725705624433656 +.441901252816163
    .139236887010717 +.482571706362417
    -.244506857304198 +.606302490573926
    -.394166300679963 -.40618253480102
    -.195460311312153 +.426521133558775
    ....

    Which is:

    real 0 imag 0
    real 1 imag 1
    real 2 imag 2
    ...

    You will have 2*NPARAMS lines (NPARAMS = 56 for Chicago)
    

    Example in bin/P01-chicago-5lines-spherical-case1

    It can also be formatted row-wise into a 1D text array.
    
output: 
    312 solutions in Matlab text format

    It looks like 
    [1.0256116377780581939+i*0.95290270838548252197
     0.087212057713832114025+i*0.010110978306756317896
     0.048192069754345270849+i*0.03896717224674180885
     1.5403585146403313555+i*0.018243345455052351056
     ...
     -0.00074739537885242771+i*-0.0029906069387749742439]

    Notice the i multiplying the imaginary part. If you need another format,
    please let me know.

    This matrix is NSOLS by NNN, where NSOLS 312 and NNN is 14 (number of
    variables).
```

For developers: the start system is compiled and don't need to be input


## Hacking

### Template internals

The code:
```C
  #include <minus.hxx>
  ...
  minus_core<chicago14a>::track(..)
```
Is shorthand for 
```C
  #include <minus.hxx>
  ...
  minus_core<chicago14a, double>::track(...)
```
where chicago14a is an enum (int) template parameter.

### Adding a new minimal problem formulation to Minus

These (rather informal) videos provide a walkthrough on how to add a new minimal
problem to MINUS:

How to add a new problem - Part 1 
  https://youtu.be/Yp5n6me04-Y
[![Add a new problem to MINUS - Part 1](https://img.youtube.com/vi/Yp5n6me04-Y/hqdefault.jpg)](https://youtu.be/Yp5n6me04-Y)

How to add a new problem - Part 2 
  https://youtu.be/KlSxoWQ0SdA

[![Add a new problem to MINUS - Part 2](https://img.youtube.com/vi/KlSxoWQ0SdA/hqdefault.jpg)](https://youtu.be/KlSxoWQ0SdA)

Let's say you have a new minimal problem, Chicago 6a, that is, variant
`a` to the 6x6 formulation of the trifocal pose from points at lines problem.
You might want to do use this since a small number of variables can mean a
drastically more efficient linear solve during solution tracking.

Steps:
- Include a new name `chicago6a` for your problem in the `enum problem` in
the beginning of `minus.h`. This is the table of problem tags of Minus.
- Specify the __constant parameters__ for the problem and formulation:
    - In the file `parameters.h`, write another include line for your problem:
```C
#include<chicago6a.h>
```
    - Write the header `chicago6a.h` by copying chicago14a.h and filling in the
      numbers and substituting the string `chicago14a` to `chicago6a`. You can add any constants
      that you feel are needed for your particular problem. For example, if your
      minmimal problem involves conics, you might want to add the number of conics.
      If your minimal problem formulation has multiple start solutions, for speed,
      you might want to indicate that as well. These are the solver settings that are
      known at compile time and can't be changed at runtime.
- Place your __evaluation functions__ into a file called `chiago6a.hxx`.
  using the existing file `chicago14a.hxx` to see how the function should be defined.
  Basically, you need to copy what is in `chicago14a.hxx`, include your function
  bodies according to the format, and substitute chicago14a to chicago6a.
    - Negate y entries by hand
    - X are refs: `%s/C X/const C<F> \&X/g`
    - Consts C are constexprs: `%s/C \(C[0-9] =\)/static constexpr C<F> \1/gc`
    - Gates G are variables `%s/C G/const C<F> G/gc` for gates (no ampersand)
    - At the end of `minus.hxx`, copy and paste the last include line to your
      problem, eg: `#include <minus/chicago6a.hxx>`
- Optional: create your app in `cmd/` immitating `cmd/minus-chicago.cxx`. If you
  have just added a new formulation solver for the same problem, make the
  necessary additions to the existing app, e.g., `cmd/minus-chicago.cxx`. 
- Optional: Default start system and gammified parameters. You might want to
  define a default start system. you might also want to define default gammified
  parameters to test your solver independent of I/O functions to convert from
  user input to gammified homotopy parameters.
    - mimmick the file `chicago14a-default.hxx` to create your own, e.g.,
    `chicago6a-default.hxx`. 
    - if your start system comes from our Macaulay2 scripts, use the commands in
    `scripts/extract-start-sols.vim` to help you translate to C++
    vector initialization format
    - in your app, say `cmd/minus-chicago.cxx` for chicago problems,
      for now you have to selectively include this. Ongoing work will remove
      this need, being only needed to do solver<chicago6a>.
- Optional: place a define to simplify your solver name: 
    this is a `using` clause towards the end of `minus.h`.
- Optional: If you are using Minus header-only, you are done! But for faster compile times 
and for including your formulation on the the official Minus codebase, and
for smaller codes, you can use our libminus, with an explicit instantiation. In
`Templates/minus-chicago-alltypes-allformulations+double.cxx`, copy the line
starting with `template class` and place your own instantiation. If it is a
different problem than Chicago, you may want to place it in a similar file with
chicago in the name. 
    
### Test Suite
  Enable tests in cmake, then indicate where your vxl build folder is located.
  When building VXL, only build VNL.


Most advanced programming tools work best under Linux.

### Mem leak with AddressSanitizer from Google

[https://github.com/google/sanitizers/wiki/AddressSanitizer]

You can build release with it, no need for slow debug. This is very fast.
Highly recommended for developing efficient code using vectors, pointers and buffers

Add this to `MINUS_EXTRA_CMAKE_CXX_FLAGS`:
```-fsanitize=address -fno-omit-frame-pointer
```
Now recompile minus and simply run it.  If nothing happens, you're golden. In
the event of any memleak, there will be a colorful output showing where it came
from, specially under Linux.

### Profiling

The best way is with kcachegrind + valgrind, by far. 
See [https://www.blogger.com/comment.g?blogID=7395958&postID=116062684092668856&bpli=1&pli=1]

In any system without valgrind or kcachegrind (eg, Macs), the easiest way is with gprof

Expect your program to take very very long - so reduce the problem / iterations
before running.

### Compilers

This was extensively tested with GCC 5
Do not use GCC 4 or 8

Intel ICC compiler with the same optimization flags as usual in Minus provided
a 2x DECREASE in speed. TODO: try other ICC-specific optimization flags

Also no success with GCC8 - slow or breaks fastmath.

Some versions of GCC might not detect your processor correctly. 
We use `-march=native` but if your processor is newer relative to your GCC,
then you might have to confirm:
```bash
 gcc -march=native -Q --help=target | grep -- '-march=' | cut -f3
```
For gcc 5, this returns broadwell, when my arch is kabylake which is compatible
to skylake, so gcc should be detecting as skylake ideally. Using newer GCC
guarantees that for me. You can also try
```bash
gcc-8 -march=native -E -v - </dev/null 2>&1 | grep cc1
```
To inspect the flags that are implied. I got good perf with gcc5, but it was
using by default `-mtune=generic`. If I put `-mtune=skylake` it might speedup
further.

To see command-line flags implied by eg `-march=native`, use:

```bash
gcc -march=native -E -v - </dev/null 2>&1 | grep cc1
```
If you want to see the compiler/precompiler defines set by certain parameters, do this:
```bash
echo | gcc -dM -E - -march=native
```
#### Selecting a compiler

Simply set the CC and CXX flags on the cmake step. Example:
```bash
CC=gcc CXX=g++ ccmake .
```
or, e.g., 
```bash
CC=gcc-7 CXX=g++-7 ccmake .
```
and so on.

Press 't' in ccmake to make sure all compiler-related paths are to the desired
compiler. If you have too many compilers around, you might want to force PATH
like so:

```bash
CC=/usr/bin/cc CXX=/usr/bin/c++ PATH=/usr/bin:$PATH ccmake .
```
Assuming you want to use the compiler at `/usr/bin`.

### Intel ICC compiler + MKL
Some tests were carried out with Intel ICC, but the gains were not significant
enough to justify a proprietary compiler. 

At Brown's CCV cluster, I used the following flags:

```
 MINUS_EXTRA_CMAKE_CXX_FLAGS	  -I${MKLROOT}/include -no-prec-div -ansi-alias  -mskylake-avx512 -xSKYLAKE-AVX512 -axSKYLAKE-AVX512
 MINUS_EXTRA_CMAKE_EXE_LINKER_FLAGS   -L/gpfs/runtime/opt/intel/2018.1/mkl/lib/intel64_lin   -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
```

### Studying assembly output
We refer to Eigen's documentation http://eigen.tuxfamily.org/index.php?title=Developer%27s_Corner#Studying_assembly_output 

For SSE/VMX use, good examples include:
Eigen/src/LU/arch/Inverse_SSE.h
As well as VXL's VNL

FMA could potentially speedup the evaluators.

Place this between code:
```asm
   asm("#------ BEGIN!");
   // some code to analyze
   asm("#------ END!");
```

use/adapt the script `scripts/minus-disassemble`

## Authors

Source code adapted and improved over Anton Leykin's and Timothy Duff's codes
from Macaulay2 by Ricardo Fabbri. The core code grew out of Macaulay2/e/NAG.cpp.
These improvements and specializations were jointly developed with Tomas Pajdla,
Benjamin Kimia and Hongyi Fan, with intensive discussion and testing from
remaining Arxiv paper co-authors while at ICERM/Brown University: Jonathan
Hauenstein, Margaret Regan, Elias Tsigaridas, Charles Wrampler, and David da
Costa de Pinho. The Chicago problem was originally formulated as a differential
version of curve-based trifocal estimation by Ricardo Fabbri, Peter Giblin, and
Benjamin Kimia.


## Acknowledgements
Minus was born out of Brown University/ICERM's 2018 Nonlinear Algebra Program (Computer Vision
Working group) and the 2019 Algebraic Vision research cluster, the former co-organized by
the authors Leykin, Hauenstein and the latter by Fabbri. Minus received an [NSF research
highlight](https://mathinstitutes.org/highlights/algebraic-computer-vision-advances-the-3d-reconstruction-of-curves-and-surfaces-from-multiple-views/).
