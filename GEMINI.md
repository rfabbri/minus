# MINUS - Minimal Problem Numerical Continuation Solver Project

## Description

This package contains a number of solvers for square nxn systems of polynomial
equations / algebraic equations, using homotopy continuation, which is basically
integrating an ODE derived from a family of polynomial systems, which tracks
solutions from a known internal system to a new, unknown system. That family is
called a homotopy. The algorithm itself is just a predictor-corrector
integrator, using typically a 4th-order Runge-Kutta with Newton's method. 
The whole procedure could be done with just Newton's method, so therefore the
method is called a "Glorified Newton's method". The key in this technique is
that by using complex numbers this procedure is always guaranteed to find the
true solution (in practice, by limiting runtime and precision we get the
solution about 95% of the time). The numerics code spend about 40% in linear
solves alone. It is the only scalable option when the equations get too
nonlinear with high dimensions, and is very stable compared to straightforward
variable elimination and solving 1D polynomials.

The difference between this code and other predict-correctors to solve a
differential equation from implicit nonlinear constraints is that 1) the
constraints are usually algebraic; 2) complex numbers are used, and 3) the
evaluators and Jacobians are symbolically available through compiled straightline
programs. Other than that, it is the same type of code used in videogame
dynamics, continuous collision engines, and real-time special effects like
elasticity and deformations used in videogames. As such, many of the code-level
optimizations found in these types of software are directly applicable here, and
I expect the AI tool (you) to provide me suggestions mapping to these highly
performant graphics/physics engine optimizations, whenever possible.

We will be coding high-performance, fast numerical code.
This is like programming a small module to be used a bigger system,
where the small module must be highly efficient and may be called millions of
times per second inside a bigger project such as structure from motion
(OpenMVG/Colmap), videogames, real-time video processing, autoomous vehicles,
robotic manipulation, etc. The downstream application we target is real-time.

## Constraints

The application must run in hundreds of microseconds at most per solve.


## IMPORTANT: Objective function to be optimized by AI
The objective function is the runtime per solve.

You can run a default worst-case solve with the command `minus-chicago`,
whose source code is in `cmd/minus-chicago.cxx`
When compiled, this generates minus-chicago, which you can run

```bash
minus-chicago -g 
```

Which will print on stderr a line with the runtime of the solver:

```
LOG Time of solver: 645ms
```

You can get to that line more easily by redirecting the solution itself to /dev/null:

```bash
minus-chicago -g  > /dev/null
```

And then filter by the string 'Time of solver:' to get the timing.

There are other problems which similarly go by e.g. `minus-linecircle` or `minus-<problemname>`,
and their timing on a default hard solve is similarly obtained with -g 


## Programming Language
- All projects we ask in this folder are in low-level C/C++ with arrays and templates, except for small scripts.
All array sizes are statically typed unless absolutely necessary.


## Folder structure


### Main folders
minus/
    The minus subfolder is where the core solver code resides. It generates a
    library but can also be used header-only in user code. It is written in
    high-performance C++ code for the CPU.

cmd/ 
    The fast executables that actually run the solvers. Their interface is such that
    it is fast and easy to call from another program.

tutorial/
    Step-by-step examples on how to get a prototype solver to work on Macaulay2,
    then handcrafting a fast solver from that in C++. If you want to add a new
    problem to Minus, this is where to go. Start at `tutorial/README.md`
    


### Secondary folders

scripts/
    This is non-essential. It contains a number of scripts to generate synthetic
    data for some solvers, to disassemble code, to call MINUS from Matlab, 
    to analyze evaluator SLPs into layered DAG diagrams, etc.

tests/
    Automated tests for the sovers and other tests. To run these, you need to
    enable VXL in CmakeLists (after downloading and bulding VXL core C++
    libraries).  TODO: we might make it dependent on Google's test suite in the
    future, although this test suite is also a bit heavy. This is currently just
    meant for more advanced users.

doc/   
    Auxiliary documentation , such as the actual
    equations used in the Chicago problem, explicitly in scalar format, for reference.
       
M2/ 
    Macaulay2 scripts that aid in tasks such as generating evaluators. These
    are currently scripts for higher-caliber problems such as Chicago (Trifocal
    problem, a 14x14 system with 312-degres of nonlinearity when expressed in
    Cayley/homogeneous quaternions).
       
algo/  
    This folder is unessential for the end-user, but contains additional tests
    that can be run with Make Tests, more advanced telemetry / data
    collection, etc. To compile this folder you will have to activate USE_VXL

config/ 
    Auxiliary CMake scripts to detect things, etc.

julia/ 
    Some extra-official scripts peple tried to implement in Julia of the
    problems used in Minus, for reference. We don't use those since all our
    evaluators actually come from Macaulay2. Our current workflow does not
    include Julia.
    
other/
    There is some other non-essential material here but mostly just historical notes.

### Main files
    minus/minus.h
    minus/minus.hxx
    cmd/minus-chicago.cxx

## Configuring and compiling
See @README.md

## C Coding Style
- Follow K&R style closely, the one in ANSI C programming language
- IMPORTANT Keep it ANSI-C: Avoid more recent C-constructs unless necessary or clearly
  better, in which case you should point out a more recent feature is being used
- Do not use arrays with a variable number of elements. Use the classic malloc
  only if the number of elements is variable, BUT IMPORTANT: stick to static arrays
  with static counts unless variable elements are strictly necessary
- C++ files will be named .h, .cxx, .hxx
- Beyond K&R and above, GNU C style guide from FSF and coding conventions from Google software engineers, For instance:
    - STRICTLY function names start at column 1 of the file
    - STRICTLY the function type goes on the line above the function name
- Lastly, you may fall back to conventions from VXL (formerly Target/TargetJr),
  such as folder structure and CMake conventions
- Spacing between symbols in for and while, print, etc like i = 3 + 4, not i=3+4.
- Minimize blank lines while keeping code legible
- Do not change one-liners in my code if I do not want to. Even if they are
  long.
- Do not repeat couts
- Make use of block comments /* */ where appropriate, not // all the time for multi-line
- For commenting large blocks of code use #if 0 .. #endif if /* */ is already in
  the block
- IMPORTANT Indentation is always in increments of 2 spaces. Never tab.
- names should be lowercase and using underscores, except for filenames that
  should use dash - rather than underscore since it is faster to type
- NEVER use spaces in filenames
- variables with local scope (local to a block) should only be used if performance is not at all
  sacrificed
- *Do not* use braces on one-liner if-else statements, while loops, for and similar
  constructs such as while if they do not require braces. But make extra effort
  to see if this really will not instroduce a bug.
- The style is high-performance pure C but only with touches of C++ when needed. The files are
  .cxx not .c
- Prefer include with <> delimiters, like #include <header.h>, leave "" as in
  #include "header.h" only for files relative to the same folder as this
  file. Warn me when there is an include "" specially in a library .h

### C Coding Technique
- Use the fastest, most efficient solution known in the world
- Use the best algorithm with the best asymptoptic complexity and best constants
  known for the problem and problem size in question
- Avoid low-level optimizations such as in-line assembly or unnecessary
  adornments and compiler extensions, but do suggest them if you think they will
  be very useful
- C++, templates and STL-specific code may be used only sparingly
    - When using anything like classes, prefer just struct and stay as close as
      possible to C in terms of efficiency
    - IMPORTANT Use C-style arrays, NEVER vector unless just prototyping or experimenting
    - Any C++ technnique, if used, must also be C++17 or earlier. Later C++
      technique should be discouraged since this project is efficiency and pure
      C-oriented. You can suggest a newer style construct if it is going to
      clearly provide speedups.
- if a quantity is unsigned, use 'unsigned' instead of 'int'
- An example of high-performance, fast C code with C++ templates and minimal C++
  relevant to this project can be found at the current project's files in
  `minus/`, and a fast user-level example command is `minus/cmd/minus-chicago.cxx`
  you can use as references, and other existing files in that MINUS project.
  
#### Commands
- Commandline executables will by default lie in cmd/ with a minus- or -cmd suffix or prefix on the name. 
- It must work with just stdin/stdout, where stdin it reads a list of
  numbers separated by spaces or newlines, and then it output to sdout, with any
  error handling that is suitable to the user done in stderr. 
- Be strictly UNIX style here with I/O, where *no news is good news*.
- No funky banners, just I/O, unless you do that on sderr. Keep stdout/stdin
  clean and efficient, as if doing in-memory interprocess communication.

## Tooling
- USE CMake with Make

### Debuggers
- DDD (gdb backend)
- Always compile C code with ample debug flags by default, -ggdb3 beyond just -g.
  TODO: You can hardcode that in CMakeLists.txt when Debug is activated.
- Use fencing by default in CMakeLists.txt in Debug mode: add flags -fsanitize=address -fno-omit-frame-pointer

### Linting: lowlevel and highlevel
- I may ask you "lintme" and you will lint all the codebase per the guidelines
  in this file
- IMPORTANT: rate the algorithms themselves, stating optimizations such as
  inefficient use of memory, possibility of efficient in-place version, and flag
  high computational complexity
- Suggest corner cases that have not been tested such as: empty or one-size
  input for a function that accepts a vector, integer overflow, etc.
- Output a grade from 0-10 and justify

## Initial setup
When I ask you "do the initial setup"  or similar, do:
- Verify that core dumps are droppable/enabled to drop in the current folder

## Bad mistakes you must not do
- IMPORTANT DO NOT use identation differently than what I told you!
- IMPORTANT DO NOT edit any files other than exactly the ones I asked for in a prompt. You may suggest which will require a careful review from me, never automated.

## Hacking
See @README_INTERNAL.md for additional experiments used while developing the
system, and for internal developer information
