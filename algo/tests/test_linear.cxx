#include <complex>
#include <cstring>
#include <testlib/testlib_test.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_vector.h>
#include <vbl/vbl_array_1d.h>
#include <vbl/vbl_attributes.h>
#include "Eigen/Core"
#include "Eigen/LU"
#include <chrono>

typedef std::complex<double> complex;
using namespace std::chrono;
#define NNN 14


static const float eps = 1e-4; // M2 is giving Ax-b up to 1e-5 

/* Lapack .h */
#if 0
extern "C" {
int sgesv_(int *n,      // number of rows in A
           int *nrhs,   // number of right hand sides
           double *a,   // n by n matrix A, on exit L&U from A=PLU
           int *lda,    // n
           int *ipiv,   // indices defining permutation P
           double *b,   // right-hand-side
           int *ldb,    // n
           int *info);  // error info

};
#endif


// Place any specific-type functions here


#if 0
// Original code solve_via_lapack_without_transposition
//
// \returns  lapack info code. 
//      >0 -> matrix is singular
//      <0 -> illegal value of an argument passed to lapack
bool linear(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  int info;
  static int permutation[NNN]; // unused
  static int bsize = 1;
  static int size = NNN;
  
  double *copyA = (double*)A;
  
  // TODO try to eliminate this memcpy and trash the original b if possible
  // memcpy  b -> x          NNN elements
  std::memcpy(x, b, 2*NNN*sizeof(double));

  double *copyb = (double*)x;  // result is stored in copyb
  
  sgesv_(&size, &bsize, copyA, &size, permutation, copyb, &size, &info);

  return info == 0;
}
#endif

/*
static bool
linear_eigen(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  using namespace Eigen;
  
  Map<Matrix<complex, NNN, 1> > xx(x);
  Map<const Matrix<complex, NNN, NNN> > AA(A,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > bb(b);
  
  xx = AA.colPivHouseholderQr().solve(bb);
  return true;
}
*/
  

static bool 
linear_eigen2(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  using namespace Eigen;
  
  Map<Matrix<complex, NNN, 1> > xx(x);
  Map<const Matrix<complex, NNN, NNN> > AA(A,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > bb(b);
  
  xx = AA.partialPivLu().solve(bb);
  return true; // TODO: better error handling
}

/*
static bool 
linear_eigen3(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  using namespace Eigen;
  
  Map<Matrix<complex, NNN, 1> > xx(x);
  Map<const Matrix<complex, NNN, NNN> > AA(A,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > bb(b);
  
  xx = AA.fullPivLu().solve(bb);
  return true;
}
*/

/*
static bool 
linear_eigen4(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  using namespace Eigen;
  
  Map<Matrix<complex, NNN, 1> > xx(x);
  Map<const Matrix<complex, NNN, NNN> > AA(A,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > bb(b);
  
  xx = AA.householderQr().solve(bb);
  return true;
}
*/

// 20s
// Direct inversion - good for smaller matrices, 5x5 etc
static bool 
linear_eigen5(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
  using namespace Eigen;
  
  Map<Matrix<complex, NNN, 1> > xx(x);
  Map<const Matrix<complex, NNN, NNN> > AA(A,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > bb(b);
  
  // xx = AA.partialPivLu().solve(bb);
  // 
  Matrix<complex, NNN, NNN> inv = AA.inverse();
  xx = inv * bb;
  // xx += inv*(bb-AA*xx); // correct
  return true; // TODO: better error handling
}
// from compute_Grabner_basis, 5point.c
#if 0
bool linear_bundler(
    const complex* A,  // NNN-by-NNN matrix of complex #s
    const complex* b,  // 1-by-NNN RHS of Ax=b  (bsize-by-NNN)
    complex* x   // solution
    )
{
    for (unsigned i = 0; i < 10; i++) {
        /* Make the leading coefficient of row i = 1 */
        double leading = A[20 * i + i];
        matrix_scale(20, 1, A + 20 * i, 1.0 / leading, A + 20 * i);

        /* Subtract from other rows */
        for (j = i+1; j < 10; j++) {
            double leading2 = A[20 * j + i];
            double scaled_row[20];
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
}
#endif


static void
print(const complex A[][NNN])
{
  // for now, print as vector
  for (unsigned i=0; i < NNN*NNN; ++i)
    std::cout << ((const complex *)A)[i] << std::endl;
}

static void
print(const complex *v)
{
  for (unsigned i=0; i < NNN; ++i)
    std::cout << v[i] << std::endl;
}

static void 
deb(complex *v, char *name)
{
  std::cout << "------------------ " << name << std::endl;
  print(v);
  std::cout << "------------------ !" << name << std::endl;
}

static void 
deb(complex v[][NNN], char *name)
{
  std::cout << "------------------ " << name << std::endl;
  print(v);
  std::cout << "------------------ !" << name << std::endl;
}

static void
transp(complex A[][NNN])
{
  for (unsigned i=0; i < NNN; ++i) 
    for (unsigned j=i+1; j < NNN; ++j)
      std::swap<complex>(A[i][j],A[j][i]);
}

static const float
Areal[NNN][NNN] =
{{-0.000003434320000, -0.000000164597000, -0.000005592940000, -0.000000488444000,  0.000052660000000,  0.003819990000000,  0.002195250000000,  0.002308850000000,  0.000050317500000,  0.000153031000000,  0.000058020700000, -0.000000088868500, -0.000000096726500, -0.000000122857000},
{-0.000003571970000, -0.000001876910000, -0.000000466329000, -0.000002115880000,  0.000613616000000,  0.003232570000000, -0.000948781000000, -0.000525815000000, -0.000032534300000,  0.000004592400000,  0.000093622900000, -0.000000121289000, -0.000000049506900,  0.000000029452700},
{ 0.000009915430000,  0.000003785870000,  0.000005756290000,  0.000004552740000, -0.001236430000000, -0.009591180000000, -0.000065404000000, -0.000984967000000,  0.000039872200000, -0.000166235000000, -0.000097623900000,  0.000000067569200,  0.000000015031900,  0.000000190989000},
{ 0.000005205040000,  0.000006428810000, -0.000010911400000,  0.000006723480000, -0.002102670000000, -0.003158650000000,  0.008339810000000,  0.007063030000000,  0.000121826000000,  0.000206442000000, -0.000176032000000,  0.000000120243000,  0.000000065235800, -0.000000017983400},
{-0.000040517800000, -0.000022522300000, -0.000001859680000, -0.000025472900000,  0.007338540000000,  0.036454100000000, -0.012766900000000, -0.007756110000000, -0.000547424000000,  0.000074537800000,  0.000208741000000,  0.000000013548000, -0.000000104352000, -0.000000177002000},
{-0.000004113200000, -0.000020704000000,  0.000056822700000, -0.000020393300000,  0.006738060000000, -0.003756360000000, -0.035645000000000, -0.031856000000000, -0.000270604000000, -0.000867531000000,  0.000129716000000,  0.000000036557600, -0.000000271177000, -0.000000182567000},
{-0.000000013554700, -0.000000019499100,  0.000000036473100, -0.000000020206600,  0.000006344820000,  0.000007357850000, -0.000026604700000, -0.000022759500000, -0.000000806556000, -0.000000987456000,  0.000000096172900,  0.000000029885300, -0.000000005967880, -0.000000059308300},
{ 0.000000045007700,  0.000000014845200,  0.000000033327600,  0.000000018301300, -0.000004817800000, -0.000044443500000, -0.000004635790000, -0.000008384310000,  0.000000439514000, -0.000000394510000, -0.000000981575000, -0.000000046341500, -0.000000048766500, -0.000000000105200},
{-0.000015909900000, -0.000011716200000,  0.000007841860000, -0.000012734500000,  0.003801810000000,  0.013296400000000, -0.010205600000000, -0.007731590000000,  0.000026650000000,  0.000315266000000, -0.000413498000000,  0.000000008109600, -0.000000033840200, -0.000000005620310},
{ 0.000012249700000,  0.000000967960000,  0.000018725100000,  0.000002124760000, -0.000313292000000, -0.013450900000000, -0.007098840000000, -0.007572090000000,  0.000445904000000, -0.000088444500000,  0.000470787000000,  0.000000083250900, -0.000000032612700, -0.000000039276000},
{ 0.000000760636000,  0.000000086582100,  0.000001097190000,  0.000000150621000, -0.000027434100000, -0.000830601000000, -0.000398421000000, -0.000432309000000,  0.000006374160000, -0.000017734700000, -0.000001935760000,  0.000000001248630, -0.000000001530010,  0.000000020339100},
{                 0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,  0.311033000000000, -0.177385000000000, -0.532478000000000, -0.906429000000000, -0.318025000000000, -0.618106000000000},
{ 0.206367000000000,  0.073497400000000, -0.329428000000000, -0.024551000000000,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0},
{                 0,                  0,                  0,                  0, -0.548244000000000,  0.111977000000000,  0.418035000000000,  0.000205810000000,                  0,                  0,                  0,                  0,                  0,                  0}};


static const float
Aimag[NNN][NNN] = 
{{-0.000003843760000, -0.000003089170000,  0.000002711870000, -0.000003304130000,  0.001004120000000,  0.003092410000000, -0.002964660000000, -0.002320460000000, -0.000086179900000, -0.000087158000000,  0.000122200000000, -0.000000149811000, -0.000000023653200,  0.000000125966000},
{ 0.000000598062000, -0.000001025790000,  0.000004321290000, -0.000000958995000,  0.000337448000000, -0.001125670000000, -0.002387580000000, -0.002218350000000, -0.000061416600000, -0.000122567000000,  0.000010519300000, -0.000000002397900,  0.000000048564600,  0.000000118954000},
{ 0.000002158390000,  0.000004680300000, -0.000010747100000,  0.000004735800000, -0.001527040000000, -0.000501380000000,  0.007208520000000,  0.006320380000000,  0.000095837000000,  0.000216593000000, -0.000129206000000,  0.000000088052300,  0.000000052783900, -0.000000059778500},
{-0.000010923400000, -0.000003131870000, -0.000009582890000, -0.000003970650000,  0.001023130000000,  0.010985100000000,  0.002012910000000,  0.002840370000000, -0.000020473600000,  0.000250932000000,  0.000078462700000, -0.000000054676900, -0.000000003243600, -0.000000237492000},
{ 0.000009992080000, -0.000010566100000,  0.000050156700000, -0.000009196590000,  0.003428130000000, -0.015451100000000, -0.026836200000000, -0.025203600000000, -0.000053286600000, -0.000781798000000,  0.000043466600000,  0.000000027371600, -0.000000201749000, -0.000000100177000},
{ 0.000050370100000,  0.000021735900000,  0.000021377200000,  0.000025630300000, -0.007098940000000, -0.047696000000000,  0.004361060000000, -0.000757140000000,  0.000607036000000, -0.000384670000000, -0.000222612000000, -0.000000005051430,  0.000000042383100,  0.000000164529000},
{ 0.000000035524100,  0.000000011251200,  0.000000027384900,  0.000000014175300, -0.000003664350000, -0.000035140700000, -0.000004381110000, -0.000007260650000,  0.000000328319000, -0.000000339314000, -0.000000777690000, -0.000000036055100, -0.000000038917300, -0.000000001629630},
{ 0.000000016077200,  0.000000024088300, -0.000000046678900,  0.000000024853100, -0.000007860250000, -0.000008099390000,  0.000033608600000,  0.000028866700000,  0.000001003970000,  0.000001253460000, -0.000000095491000, -0.000000036416600,  0.000000008784600,  0.000000074669800},
{ 0.000011359000000, -0.000000758610000,  0.000022220400000,  0.000000457230000,  0.000237576000000, -0.013035400000000, -0.009548730000000, -0.009696420000000,  0.000504572000000, -0.000050403000000,  0.000464285000000,  0.000000094690000, -0.000000041846000, -0.000000044950800},
{ 0.000012645500000,  0.000010536500000, -0.000010012000000,  0.000011221900000, -0.003420280000000, -0.010047300000000,  0.010412000000000,  0.008227970000000, -0.000093420800000, -0.000273949000000,  0.000304323000000, -0.000000020301300,  0.000000035930800,  0.000000011215000},
{ 0.000000730840000,  0.000000626364000, -0.000000643487000,  0.000000662445000, -0.000203923000000, -0.000567178000000,  0.000641222000000,  0.000511700000000,  0.000012433100000,  0.000009136600000, -0.000005058610000,  0.000000011555900,  0.000000005540040,  0.000000005054130},
{                 0,                  0,                  0,                  0,                  0,                  0,                  0,                  0, -0.293356000000000, -0.196075000000000,  0.225655000000000, -0.215290000000000,  0.033068000000000, -0.165647000000000},
{ 0.324348000000000,  0.234973000000000,  0.231526000000000,  0.888259000000000,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0,                  0},
{                 0,                  0,                  0,                  0,  0.226529000000000, -0.417930000000000, -0.081486800000000,  0.046138500000000,                  0,                  0,                  0,                  0,                  0,                  0}};

static complex b[NNN] = 
{{-0.000000000525644,- 0.000000000537608},
{ 0.000000000050204,+ 0.000000000169378},
{ 0.000000000079474,+ 0.000000000181265},
{ 0.000000000081347,- 0.000000000014625},
{ 0.000000000124128,- 0.000000002097500},
{-0.000000001613360,- 0.000000001028840},
{-0.000000000021083,+ 0.000000000010437},
{ 0.000000000009054,+ 0.000000000001964},
{ 0.000000002246100,- 0.000000007085690},
{-0.000000005798740,+ 0.000000000926283},
{-0.000000000004547,+ 0.000000000066179},
{ 0.635836000000000,- 0.118893000000000},
{ 0.921263000000000,+ 0.569246000000000},
{ 0.407041000000000,- 0.046188200000000}};


// Computed in Scilab
static const complex sol[NNN] = 
{{1.22168667179699431,+ 0.0498174125645148},
{ 1.01414252525281023,- 1.40608883721011702},
{ 0.14962345713095399,- 1.22815487409473612},
{-0.46428611099584222,- 0.03818401625592915},
{-0.73884057243577261,- 0.23266419315468043},
{-0.00075556528585659,+ 0.02021728028845612},
{-0.11764292844922636,- 0.04119032983965841},
{-0.02931068033016018,+ 0.14201886987190931},
{-0.21475574124524766,- 0.02152961056102606},
{-0.20469706313000821,+ 0.05335418324361420},
{-0.02153035145179618,- 0.29181713120224195},
{-0.81726804075896198,+ 0.01959585459446891},
{-0.1665840762930203,+ 0.87178442245356957},
{ 0.37196986306686414,+ 0.27537320811448124}};



//{{ 1.438980000000000,- 0.984509000000000},
//{ 2.090270000000000,+ 0.805246000000000},
//{-0.405540000000000,- 1.457750000000000},
//{-0.705762000000000,+ 0.166468000000000},
//{-2.172750000000000,+ 2.521430000000000},
//{-4.739900000000000,+ 0.491511000000000},
//{ 1.265470000000000,- 0.278627000000000},
//{ 0.221350000000000,+ 8.539289999999999},
//{ 0.081477200000000,- 0.061649800000000},
//{ 0.052804000000000,- 0.087110200000000},
//{ 0.103201000000000,+ 0.093230800000000},
//{ 0.494108000000000,+ 0.461057000000000},
//{ 2.118020000000000,- 1.405650000000000},
//{-2.557060000000000,+ 0.761854000000000}};

static void
join(complex *v, const float *real, const float *imag, unsigned size)
{
  for (unsigned i=0; i < size; ++i)
    v[i] = complex(real[i],imag[i]);
}


static void
test_near_v(const complex *v, const complex *w, int size)
{
  std::cout << "diff\n";
  float err=0;
  for (unsigned i=0; i < size; ++i) {
    err += std::abs(v[i]-w[i]);
    std::cout << err << std::endl;;
  }
}

static void
multiply(const complex *A, const complex *b, complex *res)
{
  vnl_matrix_fixed<complex,NNN,NNN> A_vnl(A);
  vnl_vector_fixed<complex, NNN> b_vnl(b);
  vnl_vector_fixed<complex, NNN> res_vnl(res);
  res_vnl = A_vnl*b_vnl;
  res_vnl.copy_out(res);
}

// tests solution of linear system in complex numbers
void 
test_linear()
{
  complex A[NNN][NNN];
  complex A_orig[NNN][NNN];

  join((complex *)A, (const float *)Areal, (const float *)Aimag, NNN*NNN);
  join((complex *)A_orig, (const float *)Areal, (const float *)Aimag, NNN*NNN);
  
  complex x_true[NNN];
  complex x[NNN];
  using namespace Eigen;
  

  transp(A);
  std::cout << "==========\n";
  print(b);
  
  Map<Matrix<complex, NNN, 1> > xx(x);
  Map<const Matrix<complex, NNN, NNN> > AA((complex *)A,NNN,NNN);  // accessors for the data
  Map<const Matrix<complex, NNN, 1> > b_eigen((const complex *)b);
  
  std::cout << "Is ground truth A_orig times the ground truth solution sol equal to the ground truth b?\n";
  complex bb_mul[NNN];
  multiply((const complex *)A_orig,sol,bb_mul);
  test_near_v(b,bb_mul,NNN);
  std::cout << "!Gt\n";
  PartialPivLU<Matrix<complex, NNN, NNN> > lu(AA);
  
  /* TIM 1903444 solvelinear calls, 6859039865ns time of solvelinear calls (6s)
   *
   * 19044 unoptimized:
   *            optimized: orders of magnitude
   * 0 0.429    14 - lapack
   * 1 12.174s  14,10f Householder, colpiv
   * 2 5s       7d,8f  - partial pivots     
   * 3 6s       13,38bad - full piv
   * 4 10s      12,10f Householder/partial pivoting
   * 
   * With less precision to work with, more accurate algos do better in time, 
   * less accurate do worse
   * 
   * */
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  bool retval;
//  for (unsigned i=0; i < 1903444; ++i) {
  for (unsigned i=0; i < 1000; ++i) {
 //    A[1][1] = i/20;
    // retval = linear_eigen2((const complex *)A, (const complex *)b, x);
    // xx = AA.partialPivLu().solve(b_eigen);
    lu.compute(AA);
    xx = lu.solve(b_eigen);
    __sync_synchronize();
  }
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
  std::cerr << "\033[1;32mTime of solver: " << duration << "ms\e[m" << std::endl;
  
  // printf("Retval: %i\n", retval);

  std::cout << "Computed solution\n";
  print(x);
  std::cout << "Ground-truth\n";
  print(sol);
  test_near_v(x,sol,NNN);
 
  std::cout << "Is A times the computed x equal to the ground-truth b?\n";
  complex bb[NNN];
  transp(A);
  multiply((const complex *)A,x,bb);
  test_near_v(b,bb,NNN);
}

TESTMAIN(test_linear);
