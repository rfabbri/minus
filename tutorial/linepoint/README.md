# Line - Point problem

Problem: find a line that passes to a given point.

The variables are the coefficients of the line

An auxiliary eq is necessary b = 1 (in Chicago this is a random linear constraint)


### 

Note: MonodromySolve has no real function here . It will be passed a single
parameter/solution pair, and will output a single problem,solution pair.
One would think it wouldn't alter that only solution we gave it,
but it does:


#### Before monodromy:

p0 = point {{1,1}}
-- sols0: b = 1, a = -y/x = -1
sols0 = point { {-1,1} }


#### After monodromy: 

param { {toCC(-.59048606028471173p53,-.80704783786925527p53),toCC(-.99980236903670372p53,.19880212991688566p53e-1)} }
or {-.590486-.807048*ii, -.999802+.019880*ii}

sol { {toCC(-.57432507904448216p53,.81862732887471412p53),toCC(.1p53e1,.0p53)} }

or {{-.574325+.818627*ii, 1}}

After, eval goes: 

evaluate(GS,V.BasePoint,sols#0)
    oo34 = | 4.85723e-17ii 0 |


# TO DO
When to use matrix vs point?
