-- END-LINECIRCLE                     Homotopy Continuation Tutorial                   
--
-- NAME
--      Fast HC Code Inteligencer - exactly how to craft your fast HC
--      Module: Final, fast solver. 
--
-- DESCRIPTION
--
--      This script shows how to craft a fast solver for your problem
--
--      It is broken down into Generic and Pro blocks, the pro blocks you can
--      and adapt if needed
--
--      The design-goal of the file is to closely match what will be in the fast C++
--      solver, so you can start with a generic and incrementally craft and test out the pro features
--      
--      The file is not based on generic tutorial scripts. Rather, it is rather technical
--      based on techniques that originally solved a very hard problem and for the
--      first time ever before, trifical pose from points and tangents Fabbri, etal, CVPR'10.
--      It is important to keep in mind that, while there are novelties in the -- associated m2 
--      scripts of PLMP, this version of them is the definitive account
--      of what matterscarefully for writing a fast C++ solver
--    
--      Run noob-start first once
--
-- LEGEND
--      - Noob and Pro marked bellow mean Example (simple) and Pro (fast)
--      - YOU shows where can you make it faster for your problem

-- OTHERS
--      - Think of this as a prompt precise enough for LLM to help generate a fast solver
-- 
-- TODO
--      - update with PLMP versions of common, etc, making _sure_ it
--        doesn't impact performance
--
-- AUTHORS 
--      Ricardo Fabbri (C++ and M2) and Timothy Duff (M2)
--

restart -- only useful for debugging
needs "MinusUtility.m2"
load "equations-linepoint.m2"

-- Pro -------------------------------------------------------------------------
--
-- Set online tracker options here
-- 
-- Usually the same as for the start system monodromy
-- 
-- 'null' indicates default value, use getDefault(parameterName) to see it
-- scan({CorrectorTolerance=>1e-4,
-- 	EndZoneFactor=>2e-1,
-- 	InfinityThreshold => 1e6, 
-- 	maxCorrSteps => 3, 
-- 	noOutput => true,
-- 	numberSuccessesBeforeIncrease => 2,
-- 	Precision => null,
-- 	Predictor => RungeKutta4,
-- 	stepIncreaseFactor => 2,
-- 	tStep => 5e-2,
-- 	tStepMin => 1e-5
-- 	}, 
--     opt -> setDefault(opt)) -- check if originally we used 'null' XXX
-- 
-- 
-- setDefault(CorrectorTolerance=>1e-8) -- XXX is it used in the tracker or monodromy?
--


-- Noob and Pro ----------------------------------------------------------------
-- Read start system
(p0,sols0) = readStartSys "startSys";

-- Pro -------------------------------------------------------------------------
-- verify options are set as desired
-- netList((keys options trackHomotopy)/(opt ->(opt, getDefault opt)))
-- NAGtrace 3
-- setRandomSeed 0
-- 
-- Reads from file
-- p1 = data2parameters("target_system_data.txt'); 

p1 = matrix{{1,-1 + 0*ii}};
sols1gt = matrix{{
        1 + 0*ii, 
        1 + 0*ii
}}; -- ground truth solutions of target system


evaluate(GS,p1,sols1gt) -- shoold be 0

-- p1 = {a, b, c, d, e, f};
-- where {a*(x^2+y^2)+b*x+c, d*x+e*y+f}
--
-- Pro -------------------------------------------------------------------------
-- Convert your physical data to the parameters here.
-- p1 := data2parameters targetData; -- XXX just put any abcdef here

P01 = p0 || transpose p1;

-- Pro -------------------------------------------------------------------------
-- P01 := (gammify p0)||(gammify p1); -- include Noob and Pro gammify
-- XXX suggest how to write gammify function 

-- Noob and Pro ----------------------------------------------------------------
H01 = specialize(PH, P01);
sols1 = trackHomotopy(H01, sols0)

-- Pro
-- Pass options as 3rd argument that will override what was set with setDefault
-- trackHomotopy(H01, sols0, NewOptions)

-- Evaluate check
-- 

print "Target found solutions should evaluate to 0:"
print(apply(sols1, x -> evaluate(GS, point p1, x)))

print "Start solutions should evaluate to 0:"
sols0 = point \ sols0 -- convert to list of list
print(apply(sols0, x -> evaluate(GS, point p0, x)))


-- Pro -------------------------------------------------------------------------
-- J0 = evaluate(J,sols0||p0); -- Evaluates Jacobian if desired

-- Pro -------------------------------------------------------------------------
-- Example of technique to improve speed and robustness
-- 
-- -*
-- (p,x)=fabricateChicago CC
-- evaluate(F,x||p)
-- 
-- (p0,sols0) = fabricateChicago(CC)
-- evaluate(F,sols0||p0)
-- --J0 = evaluate(J,sols0||p0);
-- --272, 279, 286
-- S = subsets(4,2)
-- for s in S do (
--     << s << endl;
--     Jpivots = flatten for i from 0 to 4 list {4*i+s#0,4*i+s#1};
--     F' = F^(Jpivots);
--     J'=diff(matrix{cameraVars},F');
--     J'0 = evaluate(J', sols0||p0);
--     << J'0 << endl;
--     S=first SVD J'0;
--     << toString((max S)/(min S)) << endl;
--     << S << endl;
--     )
-- 
-- norm evaluate(F',sols0||p0)
-- transpose sols0^{6..11} * p0^{39..44}
-- p0_(45,0)
-- 
-- netList rsort for k from 272 to 272+15-1 list (
--     --try 289?
--     Jpivots = {2, 3, 6, 7, 10, 11, 14, 15, 18, 19,k,numrows F-3,numrows F-2,numrows F-1};--rowSelector J0;
--     F' = F^(Jpivots);
--     J'=diff(matrix{cameraVars},F');
--     S = first SVD evaluate(J', sols0||p0);
--     c = (max S)/(min S);
--     e = evaluate(diff(matrix{cameraVars},F^{k}),x||p);
--     n = norm e;
--     (c, k, n, e)
--     )
-- *-
-- -- INDICES FOR SQUARE SUBSYTEM
-- Jpivots = flatten(for i from 0 to 4 list {4*i+CLBlocks#0,4*i+CLBlocks#1})  | {274,numrows F-3,numrows F-2,numrows F-1};
-- 
-- assert(# Jpivots==nvars)
-- F' = F^(Jpivots);
-- J'=diff(matrix{cameraVars},F');
-- elapsedTime GS = gateSystem(gateMatrix{dataParams},gateMatrix{cameraVars},F');
-- elapsedTime PH = parametricSegmentHomotopy GS;
-- setDefault(CorrectorTolerance=>1e-8)
-- 
-- -- filter path jumps during monodromy
-- filterEval = (p,x) -> (
--     -- false iff residual small
--     resid := norm evaluate(F,x||p);
-- --    << "residual: " << resid << endl;
--     (resid > 1e-4)
--     )
-- 
-- seeds = {16314, 62236, 77740, 69544, 11665, 1055}
-- H = hashTable apply(subsets(4,2),seeds, (ss,s) -> ss=>s)
-- 
