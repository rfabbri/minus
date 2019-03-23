<!--<img src="vxd-logo.gif" width="215" height="84" /> -->


## MINUS: MInimial problem NUmerical continuation Solver
This package originated in solving medium-sized  (eg, degree > 100) square problems
(ie, exactly-determined, 0-dimensional), notably in computer vision where
trifocal minimal problems from points and lines are of importance (as in
curve-based structure from motion, where lines are tangents to curves).
As of this date, such problems are too high degree to be solved symbolically,
while being low enough degree to allow a global technique rather than the usual
local Levenberg-Marquardt of structure from motion that requires an
initialization. In fact, the solutions found using this technique for square
problems can be used to initialize Levenberg-Marquardt for an overconstrained problem.

Minus is split into three parts:
- An efficient library for use in your C++ programs, no dependencies beyond
  standard C++ and Eigen (included)
- A simple commandline program easy to interface with other programs (eg, Matlab and Macaulay2)
- Optional: extensive tests useful for tuning the algorithm to a machine architecture and
  compiler. These are disabled by default, and requires [VXL core library](https://vxl.github.io) (*optional* - most users don't need this).
 
The theory and practice associated to Minus is described in
"Trifocal Relative Pose from Lines at Points and its Efficient Solution"
(Arxiv). For datasets and curve-based SfM code, see the [website](http://multiview-3d-drawings.sourceforge.net).

## Usage in C++ programs
For use in your program, we provide a C++ header-only library 
Simply do:
```C
#include <minus.hxx>
```

And use `minus<chicago>` like so:
```C
  minus<chicago>::track(minus<chicago>::DEFAULT, start_sols, params, solutions);
```
to solve a trifocal pose problem from lines at points ("Chicago"), using the default
formulation.  See the full example in `cmd/minus-chicago.cxx`.  This is
efficient *static* code, so no allocations are performed.  

The size and key parameters of the minimal problem are hardcoded as template
parameters in advance, for efficiency. `minus<>` is a shorthand for a generic
template, so you have full control to add your own compiled formulations, or
change the floating point implementation.  You can easily specify the
formulation by using, e.g., `minus<chicago14a>` instead of `minus<chicago>`,
which will use a 14x14 formulation for the Chicago problem.  Other
available instances are available, e.g. `minus<chicago6a>` to solve a 6x6
formulation for the same problem instead. 

To solve another problem, say the `Cleveland` trifocal pose problem from mixed
points and lines, simply use the problem tag:
```C
  minus<cleveland>::track(minus<cleveland>::DEFAULT, start_sols, params, solutions);
```

### Further control on template parametrs

If you want to have full control on templating, say, to change from `double` to
`float`, `minus<problem>` is just a shorthand for `minus_core<problem,double>`.
So you can use `minus_core<chicago14a, float>`. See the section Hacking for how
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
./minus-chicago -g         # -g profiles with a predefined input, to get a time

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
that. I recommend running on a couple of machines you now minus runs fast)

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
  minus<chicago14a>::track(..)
```
Is shorthand for 
```C
  #include <minus.hxx>
  ...
  minus_core<312, 14, 56, chicago14a, double>::track(...)
```
where chicago14a is an enum (int) template parameter.

### Adding a new minimal problem formulation to Minus

Let's say you have a new minimal problem, Chicago 6a, that is, variant
`a` to the 6x6 formulation of the trifocal pose from points at lines problem.

- Include a new name `chicago6a` for your problem in the `enum problem` in
`minus.h`. This is the table of problem tags of Minus.
- Place your evaluation functions into a file called `chiago6a.hxx`.
  using the existing file `chicago14a.hxx` to see how the function should be defined.
  Basically, you need to copy what is in `chicago14a.hxx`, include your function
  bodies according to the format, and substitute chicago14a to chicago6a.
    - Negate y entries by hand
    - X are refs: `%s/C X/const C<F> \&X/g`
    - Consts C are constexprs: `%s/C \(C[0-9] =\)/static constexpr C<F> \1/gc`
    - Gates G are variables `%s/C G/const C<F> G/gc` for gates (no ampersand)
- Optional: Default start system and gammified parameters
    - mimmick the file `chicago14a-default.hxx` to create your own, eg,
    `chicago6a-default.hxx`. 
    - if your start system comes from our Macaulay2 scripts, use the commands in
    `scripts/extract-start-sols.vim` to help you translate to C++
    vector initialization format
    - in your app, say `cmd/minus-chicago.cxx` for chicago problems,
      for now you have to selectively include this. Ongoing work will remove
      this need, being only needed to do minus<chicago6a>.
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
```-fsanitize=address -fno-omit-frame-pointer```

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


## Authors

Source code adapted and improved over Anton Leykin's and Timothy Duff's codes
from Macaulay2 by Ricardo Fabbri. The core code grew out of Macaulay2/e/NAG.cpp.
These improvements and specializations were jointly developed with Tomas Pajdla,
Benjamin Kimia and Hongyi Fan, with intensive discussion and testing from
remaining Arxiv paper co-authors while at ICERM/Brown University: Jonathan
Hauenstein, Margaret Regan, Elias Tsigaridas, Charles Wrampler, and David da
Costa de Pinho.


## Acknowledgements
Minus was born out of Brown University/ICERM's 2018 Nonlinear Algebra Program (Computer Vision
Working group) and the 2019 Algebraic Vision research cluster, the former co-organized by
the authors Leykin, Hauenstein and the latter by Fabbri.
