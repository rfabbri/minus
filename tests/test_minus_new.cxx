void
test_world_to_camera()
{
  //  This is exactly the case in the slides in big-notes/trifocal/trifocal.key
  // Extracted from the original synthcurves dataset
  // synthcurves-multiview-3d-dataset/spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object
  //
  // Frame files: frame_..42, 54, 62
  //
  // Points:
  // 
  // 620
  // 3011  tangents
  // 3389  tangents
  // 0-based ids. +1 to get file line
  // 
  // Extracting this data from those files: scripts/getlines.sh

  // NOTE: The point order seems to match Hongyi,
  // though he might have done an off-by-1 mistake for point indexing (3012 vs
  // 3011)
  double p[view][point][coord];
  // points for frame 42
  // + sed -n '3012p;3390p;621p' frame_0042-pts-2d.txt
  141.01103052308988595 270.45312297462106699
  239.89822517853363593 86.442049763307068133
  286.7673976130331539 217.06531260627261304
  
  // points for frame 54
  // + sed -n '3012p;3390p;621p;' frame_0054-pts-2d.txt
  241.41513314836856807 447.15662243793082098
  123.95973916849976604 213.90676875312345828
  257.04360648826406077 159.4404341695463927
  
  // points for frame 62
  // + sed -n '3012p;3390p;621p' frame_0062-pts-2d.txt
  375.60750199363729962 277.22372936832925916
  295.57132984990698787 147.80261937455236421
  240.78946527513195974 410.13737156824942076

  double tgt[view][tgt][coord];
  // tangents for frame 42
  // + sed -n '3012p;3390p' frame_0042-tgts-2d.txt
  0.9536809622336909209 -0.3008199166827579818
  0.0082601187924503903515 -0.99996588463683822035

  // tangents for frame 54
  // + sed -n '3012p;3390p' frame_0054-tgts-2d.txt
  0.18491347256048701331 -0.9827548054655455001
  -0.99542450475950383648 0.095551322985590561587
  
  // tangents for frame 62
  // + sed -n '3012p;3390p' frame_0062-tgts-2d.txt
  0.77931350598248894102 -0.62663423094599701724
  0.76492323888624347283 0.64412144709812224619


  // intrinsic params equal to every frame
  double K[][];
  2584.9325098195013197 0 249.77137587221417903
  0 2584.7918606057692159 278.31267937919352562
  0 0 1

  // extrinsic params for each frame
  double R[][][];
  double C[][][];

  // extrinsics for frame 42
  // + cat frame_0042.extrinsic
  -0.097305153950172085242 -0.22322794404612877894 -0.96989741313794208821
  0.96072075769186959793 0.23341709945525662695 -0.15010690664274928263
  0.25989869710080021337 -0.94640675329312473618 0.1917470327448986267

  -295.2090359311167731 1074.0075457376335635 -236.40439390871563319
  // extrinsics for frame 54
  // + cat frame_0054.extrinsic
  0.81820480546085894158 0.07511341824191355987 -0.56999900940332604016
  0.54313649506229122466 -0.42609616539484057585 0.72349485524588375007
  -0.18853022052765816552 -0.90155423144469115648 -0.38943076883056443327

  194.82952402681169701 1020.3676638972305 431.76461692675769655
  // extrinsics for frame 62
  // + cat frame_0062.extrinsic
  -0.61853492444688140672 -0.60388598633423984374 -0.50272881631015808868
  -0.025677143402306701336 -0.62392573537516660132 0.78106168837247091918
  -0.78533765448129877473 0.4960225723146879373 0.37041379052099593361

  887.07508137499985423 -562.68690102473453862 -415.57529638919055515
  


  // Generate lines
  // pLines is a 15x3 matrix of line coefs  (we use view-line-point index, this
  // is inverted to match Hongyi)
  //    1    -- l_1_1 --
  //    2    -- l_1_2 --
  //    3    -- l_1_3 --
  //    4    -- l_2_1 --
  //    5    -- l_2_2 --
  //    6    -- l_2_3 --
  //    7    -- l_3_1 --
  //    8    -- l_3_2 --
  //    9    -- l_3_3 --
  //    10   -- l_4_1 --
  //    11   -- l_4_2 --
  //    12   -- l_4_3 --
  //    13   -- l_5_1 --
  //    14   -- l_5_2 --
  //    15   -- l_5_3 --
  //    
  //    l_line_view
  //    
  //    These lines are:
  //
  //    l_1: Point 1&2  (A, B)
  //    l_2: Point 1&3  (A, C)
  //    l_3: Point 2&3  (B, C)
  //    l_4: Tangent at Point 1 (A)
  //    l_5: Tangent at Point 2 (B)
  
  //    XXX

  
  plines[15][3] = {
  };

  unsigned id []  = {
    {0, 1},
    {0, 2},
    {1, 2} 
  };

  // plines 1-9 (between points): 
  unsigned i=0;
  for (unsigned l=0; l < 3; ++l)
    for (unsigned v=0; v < 3; ++v)
       cross(pts[v][id[v][0]], pts[v][id[v][1]], plines[i++]);

  // plines 10-15 (tangents): 
  for (unsigned p=0; p < 2; ++p)
    for (unsigned v=0; v < 3; ++v)
       cross(p[v][p], p[v][p]+tgt[v][p], plines[i++]);
    
  complex params[2*M::nparams]; // start-target param pairs, P01 in chicago.m2, like params_ 
  
  // triple
  
  //  This is iccv-chicago-src/chicagoTestDataSph file 5linesCase1.txt
  //  which was used to generate the default run. It already inverts point  coords
  //
  //  5 lines (Hongy matlab this is visibleLines, inverted:
  //0.879008519574604 0.476806063840702 0.038623735269322 0.894813299896884 -0.446440542880738 0.032208045798003 0.704559167216148 0.709645249326512 -0.033704562238054
  //-0.707245886556431 0.706967648445817 -0.027867353589865 0.740915644124991 0.671598100273408 -0.041627953405298 0.978562374828962 -0.205950670232225 0.014188247805869
  //-0.995986793114319 0.089500323696923 0.002805786591495 -0.159509552710649 0.987196385018730 0.016889155865539 0.590904880455181 -0.806741236242606 -0.029119177786994
  //0.383774765752577 0.923426731891359 0.018967747764466 0.996132185054309 0.087867342619014 -0.002506238319557 0.594754053692778 0.803907715857989 -0.038939157023577
  //0.998695342916704 0.051064782741988 0.007599111258801 -0.045862067346128 -0.998947781807808 -0.027129197306103 -0.610567512482729 0.791964211754958 0.030062848861502
  //Camera Motions:
  //0.456457606228917 0.889159884956647 0.032266897893194 -0.659449197822556 0.313742799357977 0.683148747596164 0.597304954949274 -0.333106821957916 0.729566060037158
  //-16.067155449351297 -772.430415639545799 304.233386926931757
  //-0.264893260894465 -0.442357100009250 0.856826561448759 -0.513078745490601 0.817000282830920 0.263174350535889 -0.816444585540554 -0.369906385353920 -0.443381895002389
  //-953.788605479449529 -308.108851792374651 1642.743688891120200


  // Steps to transform this to params:  minus/src-iccv-chicago/t
  //
  // p0 is initial parameters OK
  // params_start_ in chicago14a-default.hxx OK

  // sols is initial solutions  OK
  // start_sols_already in chicago-default.hxx  OK
  // 
  // (p0,sols) = readStartSys "startSys"; OK
  // 
  // ------------------------------------------------------------ DONE
  // pLines = parseFile: from the above, returns  OK
  // pLines is a 15x3 matrix of line coefs made by Hongyi
  //    1    -- l_1_1 --
  //    2    -- l_1_2 --
  //    3    -- l_1_3 --
  //    4    -- l_2_1 --
  //    5    -- l_2_2 --
  //    6    -- l_2_3 --
  //    7    -- l_3_1 --
  //    8    -- l_3_2 --
  //    9    -- l_3_3 --
  //    10   -- l_4_1 --
  //    11   -- l_4_2 --
  //    12   -- l_4_3 --
  //    13   -- l_5_1 --
  //    14   -- l_5_2 --
  //    15   -- l_5_3 --
  // 
  //    l_line_view
  //       
  //    These lines are:
  //
  //    l_1: Point 1&2  (A, B)
  //    l_2: Point 1&3  (A, C)
  //    l_3: Point 2&3  (B, C)
  //    l_4: Tangent at Point 1 (A)
  //    l_5: Tangent at Point 2 (B)
  // 
  // --------------------------------------------------------------- 
  // //////////////////////////////solveChicago(p0, pLines /*targetLines*/, sols /*startSols*/);
  //    P1 := lines2Parameters targetLines;
  //    P01 := (gammify P0)||(gammify P1);  // which is what we need..
  //    
  // Each one in turn:
  //    
  // ---P1 := lines2Parameters targetLines; ------------- XXX
  //    P1 is pDouble||pTriple||pChart  //  Hongyi: [pDouble; tripleChart; XR'; XT1'; XT2'];
  //    size    27       12        17 = 56
  //    
  //    pDouble: 9x3 out of the 15x3 of the pairsiwe lines, linearized as 27x1
  //        Tim: pDouble is matrix(targetLines^{0..8},27,1);
  //        Order: row-major
  //        Hongyi: same
  //
  //
  //    pTriple: the 6 tangent lines in the coordinates of the point-point
  //        lines, in the form  (x1,y1,x2,y2,x3,y3,...): 12x1
  //        
  //        tripleIntersections := {{0,3,9},{0+1,3+1,9+1},{0+2,3+2,9+2},{0,6,12},{0+1,6+1,12+1},{0+2,6+2,12+2}};
  //        for each ind in tripleIntersections
  //            subm = get all ind lines in pLines
  //            n = numericalKernel(subm', kTol);
  //            // find intersection point - if we already have it, just use it!
  //            // if our data input is in the form of 3 points and lines
  //            // through them, simply convert keeping the invariant that the
  //            // line really goes through the points
  //            // if the tangent line does not go exactly through the point,
  //            // need to use svd to make sure it goes through it.
  //            // unless they have been randomized. Even then,
  //            // we should not find an average intersection point between the
  //            // 3 lines, but should keep the point as the more reliable
  //            // measurement and only adjust the line to the point.
  //            // This will be left for later.
  //
  //            // intersection point as non-normalized coordinates
  //            ind = (1/n_(2,0))*n^{0,1})
  //                     n^{0,1} is n(0:1,:)   if n is 3x1 -> n(0:1)
  //                     n_(2,0) is not n_{2,0} but simply -> n(2)
  //
  //
  //            Explanation: 
  //            The tripleChart computation starts with line indices:
  //            ind = [0 3 9;1 4 10; 2 5 11; 0 6 12; 1 7 13; 2 8 14]+1;
  //            
  //            Lets take ind [0 3 9]+1, or [1 4 10].
  //            This is l_1_1, l_2_1, l_4_1.
  //            These are lines 1, 2, and 4 at point 1.
  //            All these lines intersect at point 1.
         
  //            Now both your code and Tim's will now generate a 3x3 matrix
  //            for svd:
         
  //            -- l_1_1 --
  //            -- l_2_1 --
  //            -- l_4_1 --
         
  //            When you compute the numerical kernel, it means you are finding an
  //            intersection point. The intersection point is the point p such that the
  //            matrix times this p is as close to zero as possible. This means p satisfies
  //            all equations as closely as possible.
         
  //            This code only makes sense if your input is only lines and you want to find
  //            an approximate intersection, since in general they will not intersect
  //            exactly at the same point in the presence of noise.
  //
  //
  //            Code without svd: /Users/rfabbri/cprg/vxlprg/lemsvpe/minus/scripts/test_triple_intersections.m
  //                
  //        
  //    pChart: just unit rands 17x1
  //        sphere(7,1)|sphere(5,1)|sphere(5,1)
  //
  //
  // --- gammify ---
  //
  //
  // 9 random complex numbers (rand x + i rand y), non unit, seemingly uniform
  // Corresponding to the 9 pairwise lines. Seems unit is a better idea
  // numerically.
  //
  // gamma1 .. gamma9
  // 
  // diag0 Generate a 3*9 = 27 entry thing by duplicationg gammas
  // gamma1
  // gamma1
  // gamma1
  // gamma2
  // gamma2
  // gamma2
  // ...
  // gamma9
  // gamma9
  // gamma9

  //  tripleIntersections := {{0,3,9},{0+1,3+1,9+1},{0+2,3+2,9+2},
  //  {0,6,12},{0+1,6+1,12+1},{0+2,6+2,12+2}};
  
  //  for each triple intersection i
  //    Get the first two (point-point) lines
  //    
  //    diag1(i) = conjugate(gammas(tripleIntersection(i)(0)))
  //    diag1(i+1) = conjugate(gammas(tripleIntersection(i)(1)))
  //    
  //  diag2 := 7 times a fixed random(); -- t chart gamma
  //  diag3 := 5 times a fixed random(); -- q chart, cam 2, gamma
  //  diag4 := 5 times ...               -- q chart, cam 3, gamma
  //  p' := (diag0|diag1|diag2|diag3|diag4).*p;
  //  total    27   12    7      5    5 = 56
    
  
  // XXX doing: prototype conversions in pseudocode to see if simplifications
  // emerge
}

