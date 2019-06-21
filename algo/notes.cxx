//
// https://hackernoon.com/programmers-guide-to-linear-systems-4581e40fe6ad
//
// Macaulay2 tracker: trackHomotopy(H01, startSols, options);
// This is from the basic m2 engine: ./M2/Macaulay2/packages/NumericalAlgebraicGeometry/track.m2
// It calls rawHomotopyTrack c function interfaced through d
// These raw* functions are interfaced at d/interface2.d
//
//
//
// rawSetParametersPT in NAg.cpp init_dt -> raw code initDt -> tim m2 tStep == NULL (def?) -> my tracker init_dt
//
// min_dt --> raw code minDt -> tStepMin 1e-7
// 
//
// The main continuation is being done at:
// bool HomotopyConcrete<RT, FixedPrecisionHomotopyAlgorithm>::track(
// file SLP-imp.hpp
//
//
//
//
//
// To evaluate specialized homotopy as in the code, the first evaluation for sol
// 0 would be:
// evaluateHt(H01,sols#0,0) or evaluateHx(H01,sols#0,0)
//
//
//
//
//
// 
//
            matrix_scale(20, 1, A + 20 * i, leading2, scaled_row);
            matrix_diff(20, 1, 20, 1, A + 20 * j, scaled_row, A + 20 * j);
        }
    }

    /* Now, do the back substitution */
    for (i = 9; i >= 0; i--) {
        for (j = 0; j < i; j++) {
            double scale = A[20 * j + i];
            double scaled_row[20];
            matrix_scale(20, 1, A + 20 * i, scale, scaled_row);
            matrix_diff(20, 1, 20, 1, A + 20 * j, scaled_row, A + 20 * j);
        }
    }
    
    /* Copy out results */
    for (i = 0; i < 10; i++) {
        memcpy(Gbasis + i * 10, A + i * 20 + 10, sizeof(double) * 10);
    }



rfabbri@cortex:~/..inus/tests$ st n_ff
N	min	max	sum	mean	stddev
312	42	19258	218792	701.256	2106.94
rfabbri@cortex:~/..inus/tests$ man st
rfabbri@cortex:~/..inus/tests$ st n
N	min	max	sum	mean	stddev
312	42	19255	218733	701.067	2105.56
