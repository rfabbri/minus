-- EX-RUN                          Homotopy Continuation Tutorial                              EX-RUN
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
--      Run ex-start first once
--
-- LEGEND
--      - Ex and Pro marked bellow mean Example (simple) and Pro (fast)
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

restart
needs "parser.m2"
elapsedTime needs "chicago.m2"

(p0,sols) = readStartSys "startSys";
-- verify options are set as desired
netList((keys options trackHomotopy)/(opt ->(opt, getDefault opt)))
NAGtrace 3
setRandomSeed 0


-- (pLines, x) = parseFile(1,Tests=>true);
-- elapsedTime K = solveChicago(p0, pLines, sols);
