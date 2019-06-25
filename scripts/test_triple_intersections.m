% prototyping a non-svd solution

pLinesind = [.879009  .476806 .0386237
-.707246 .706968 -.0278674
.383775  .923427 .0189677];

l = pLinesind';

l1 = l(:,1);
l2 = l(:,2);
l3 = l(:,3);

l1l1 = l1'*l1;
l1l2 = l1'*l2;
l2l1 = l1l2;
l2l2 = l2'*l2;
l3l1 = l3'*l1;
l3l2 = l3'*l2;

a = cross([l1l1 l2l1 l3l1], [l1l2 l2l2 l3l2]);

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

