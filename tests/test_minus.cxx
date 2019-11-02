// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date February 2019
//
// Tests more comprehensive runs of minus using the public interface
// 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <testlib/testlib_test.h>
#include <minus.h>

#define Float double
typedef minus<chicago14a> M;
static constexpr Float tol = 1e-3;
typedef std::complex<Float> complex;
using namespace std::chrono;

// Start solutions hardcoded for efficiency.
// If you want to play with different start sols,
// write another program that accepts start sols in runtime,
// but keep this one lean & mean.
#include <chicago14a-default.hxx> 
// We include it separately so they don't clutter this app,
// neither minus.h, and can be reused by other progs
// TODO(developer note): make this part of Minus' template as a specialization. 
// But for efficiency I chose to do it outside.
// Perhaps a minus class should be written that wraps the lean minus_core.
// And in _that_ one, we put these default vectors depending on template tag.

#define  M_VERBOSE 1     // display verbose messages

void
test_against_ground_truth(const M::solution solutions[])
{
  // compare solutions to certain values from Macaulay2
  // two random entries. Just a sanity check against the original code prototype.
  // Not a Full comparison to ground truth cameras!
  bool ok=false;
  if (std::abs(solutions[1].x[1] - complex(-.25177177692982444e1, -.84845195030295639)) <= tol &&
      std::abs(solutions[M::nsols-2].x[2] - complex(.7318330016224166, .10129116603501138)) <= tol) {
    std::cerr << "LOG solutions look OK\n"; ok = true;
  } else  {
    std::cerr << "LOG \033[1;91merror:\e[m solutions dont match original code. Errors: ";
    std::cerr << std::abs(solutions[1].x[2] - complex(-.25177177692982444e1, -.84845195030295639)) << ", "
        << std::abs(solutions[M::nsols-2].x[2] - complex(.7318330016224166, .10129116603501138)) << std::endl;
  }
  TEST("Solutions match original code", ok, true);
}

void
test_full_solve()
{
  std::cerr << "Starting path tracker" << std::endl;
  
  static M::solution solutions[M::nsols];
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  M::track_all(M::DEFAULT, start_sols_, params_, solutions);
  
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(t2 - t1).count();
  std::cerr << "LOG \033[1;32mTime of solver: " << duration << "ms\e[m" << std::endl;
  
  TEST("Did it track first solution?", solutions[0].t > 0, true);
  TEST("Did it track the last solution?", solutions[M::nsols-1].t > 0, true);

  test_against_ground_truth(solutions);
}

//TODO to io::
//
// \param[in] tgts: three tangents, one at each point.
// Only two tangents will actually be used. If one of the points
// in each image has no reliable or well-defined tangents,
// you can pass anything (zeros or unallocated memory); 
// it will be ignored. 
// only tgt[view][id_tgt0][:] and tgt[view][id_tgt1][:] will be used.
//
// id_tgt0  < id_tgt0 < 3
// 
void
point_tangents2params(F p[3][3][2], F tgt[3][3][2], unsigned id_tgt0, unsigned id_tgt1, C<F> * __restrict__ params/*[static 2*M::nparams]*/)
{
  Float plines[15][3];
  point_tangents2lines(p, tgt, plines);
  // get_params_start_target<Float>(plines, params);
  io::lines2params<Float>(plines, params);
}


// unit quaternion to rot matrix
void
quat2rotm(const F q[4], F r[9])
{
  const F 
    x2 = q[0] * q[0],  xy = q[0] * q[1],  rx = q[3] * q[0],
    y2 = q[1] * q[1],  yz = q[1] * q[2],  ry = q[3] * q[1],
    z2 = q[2] * q[2],  zx = q[2] * q[0],  rz = q[3] * q[2],
    r2 = q[3] * q[3];
    
  *r++ = r2 + x2 - y2 - z2;    //  rot(0,0) = r[0]
  *r++ = 2. * (xy - rz);       //  rot(0,1) = r[1] 
  *r++ = 2. * (zx + ry);       //  rot(0,2) = r[2] 
  *r++ = 2. * (xy + rz);       //  rot(1,0) = r[3] 
  *r++ = r2 - x2 + y2 - z2;    //  rot(1,1) = r[4]
  *r++ = 2. * (yz - rx);       //  rot(1,2) = r[5] 
  *r++ = 2. * (zx - ry);       //  rot(2,0) = r[6] 
  *r++ = 2. * (yz + rx);       //  rot(2,1) = r[7] 
  *r   = r2 - x2 - y2 + z2;    //  rot(2,2) = r[8]
}

// \returns true if the solution is real, false otherwise
// 
//  rs: real solution; f holds solution R12, t12, R13, T13 row-major
bool
get_real(const solution *s, F rs[NNN])
{
  // Hongyi function realSolutions = parseSolutionString(output)
  // solutions = reshape(solutions,[14,length(solutions)/14]);

  eps = 10e-6;
  
  imagSolutions = imag(solutions);
  
  // TODO
  // Fancy way to convert to real is to check if the complex number is close to
  // horizontal then get absolute value.
  /*
  for (unsigned var = 0; var < NNN; ++var)  // differs from Hongyi criterion
    if (s->x[var].real() < eps && s->x[var].real() >= eps
        || std::abs(std::tan(std::arg(s->x[var].imag()))) >= eps)
      return false;
  
  F real_solution[NNN];
  for (unsigned var = 0; var < NNN; ++var) 
    real_solution[var] = ((s->x[var].real() >= 0) ? 1 : -1) * std::abs(s->x[var]);
  */

  for (unsigned var = 0; var < NNN; ++var)
    if (std::abs(s->x[var].imag()) >= eps)
        return false;
  
  F rs[NNN]; // real solution
  for (unsigned var = 0; var < NNN; ++var) 
    rs[var] = s->x[var].real();

  // quat12 rs(0:3), quat12 rs(4:7)
  F *p = rs;
  F norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
  p[0] /= norm; p[1] /= norm; p[2] /= norm; p[3] /= norm;
  p += 4;
  F norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
  p[0] /= norm; p[1] /= norm; p[2] /= norm; p[3] /= norm;

  //  T12 = solutions(9:11);
  //  T13 = solutions(12:14);
  //  R12 = quat2rotm(transpose(quat12));
  //  R13 = quat2rotm(transpose(quat13));

  return true;
}

void
solutions2cams()
{
  // for each solution
    get_real();
    // build cams by using quat2rotm
}

//TODO to io::

void
test_end_user_interface()
{
  // static data for points and cams

  // M::solve(M::DEFAULT, start_sols_, points, cameras);
  {
  io::point_tangents2params(points, tangents, params);
  M::track_all(M::DEFAULT, start_sols_, params_, solutions);
  
  TEST("Did it track first solution?", solutions[0].t > 0, true);
  TEST("Did it track the last solution?", solutions[M::nsols-1].t > 0, true);
  test_against_ground_truth(solutions);

  // test solutions are good before formatting them out,
  // just like in minus-chicago

//  io::solutions2cams(solutions, cameras);
//  test_final_solve_against_ground_truth(solutions);
  }
}

void
test_minus()
{
  test_full_solve();
  test_end_user_interface();
}

TESTMAIN(test_minus);
