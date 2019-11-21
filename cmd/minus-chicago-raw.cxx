// 
// \author Ricardo Fabbri based on original code by Anton Leykin 
// \date February 2019
// 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <minus/minus.h>

#define Float double
typedef minus<chicago14a> M;
static constexpr Float tol = 1e-3;
typedef std::complex<Float> complex;
using namespace std::chrono;

// Start solutions hardcoded for efficiency.
// If you want to play with different start sols,
// write another program that accepts start sols in runtime,
// but keep this one lean & mean.
#include <minus/chicago14a-default.hxx> 
// We include it separately so they don't clutter this app,
// neither minus.h, and can be reused by other progs
// TODO(developer note): make this part of Minus' template as a specialization. 
// But for efficiency I chose to do it outside.
// Perhaps a minus class should be written that wraps the lean minus_core.
// And in _that_ one, we put these default vectors depending on template tag.

#define  M_VERBOSE 1     // display verbose messages
void
print_usage()
{
  std::cerr << "Usage: minus [input solutions]\n\n";
  std::cerr << "If no argument is given, 'input' is assumed stdin,\n\
  'solutions' will be output to stdout\n";
  exit(1);
}

bool stdio_=true;  // by default read/write from stdio

// Output solutions in ASCII matlab format
//
// ---------------------------------------------------------
// If in the future our solver is really fast, we may need Binary IO:
// complex solutions[NSOLS*NVE];
// 
// To read this output file in matlab, do:
// fid = fopen(fname,'r');
// a_raw = fread(fid,'double');
// fclose(fid);
//
// Reshape a to have proper real and imaginary parts
// a = a_raw(1:2:end) + i*a_raw(2:2:end);
// 
template <typename F=double>
static bool
mwrite(const M::solution s[M::nsols], const char *fname)
{
  bool scilab=false;
  std::string imag("+i*");
  if (scilab) imag = std::string("+%i*");
    
  std::ofstream fsols;
  std::streambuf *buf;
  
  if (stdio_) {
    buf = std::cout.rdbuf();
    std::cout << std::setprecision(20);
  } else {
    fsols.open(fname,std::ios::out);
    if (!fsols) {
      std::cerr << "minus: error, unable to open file name" << std::endl;
      return false;
    }
    buf = fsols.rdbuf();
    fsols << std::setprecision(20);
  }
  
  std::ostream out(buf);
  out << std::setprecision(20);
  out << "[";
  for (unsigned i=0; i <M::nsols; ++i) {
    for (unsigned var=0; var < M::nve; ++var) {
      out << s[i].x[var].real() << imag << s[i].x[var].imag();
      if (i*var +1 < M::nve * M::nsols) 
        out << std::endl;
      // BINARY fsols.write((char *)(s[i].x[var]),2*sizeof(double));
    }
  }
  out << "]\n";
  
  if (!stdio_) fsols.close();
  return true;
}

// reads into the global variable start_sols_
// Format is just like P01 variable in solveChicago in chicago.m2
// But here there is no imaginary 'i' string:
//
// P01(0).real()  P01(0).imag()      // I mean P01(0) or P01#0
// P01(1).real()  P01(1).imag()
// P01(2).real()  P01(2).imag()
// P01(3).real()  P01(3).imag()
// ...
//
// The file can also be one line, listing the above in row-major order like so:
// 
// P01(0).real()  
// P01(0).imag() 
// P01(1).real()
// P01(1).imag()
// ...
// 
// It is up to the user to build this from an actual input for a target system,
// be it point-tangents as in Ric's format, be it a linecomplex as in Hongy's format
//
// This format is generic enough to be adapted to M2 or matlab
template <typename F=double>
static bool
mread(const char *fname)
{
  std::ifstream infp;
  std::istream *inp = &std::cin;
  
  if (!stdio_) {
    infp.open(fname, std::ios::in);
    if (!infp) {
      std::cerr << "I/O Error opening input " << fname << std::endl;
      return false;
    }
    inp = &infp;
  }
  
  std::istream &in = *inp;
  in.exceptions(std::istream::failbit | std::istream::badbit);
  unsigned i=0;
  F *dparams = (F *)params_;
  while (!in.eof() && dparams != (F *)params_+2*2*M::f::nparams) {
      try {
      in >> *dparams++;
      // std::cerr << "reading " <<  *(dparams-1) << std::endl;;
      if (in.eof()) {
        std::cerr << "I/O Error: Premature input termination\n";
        return false;
      }
      in >> *dparams++;
      } catch (std::istream::failure &E) {
        std::cerr << "I/O Error: Invalid input conversion or other error\n";
        return false;
      }
  }
  if (dparams != (F *)params_+2*2*M::f::nparams)
    std::cerr << "I/O Premature input termination\n";
//  for (unsigned i=0; i < 2*NPARAMS; ++i)
//    std::cerr << "D " << params_[i] << std::endl;
  return true;
}

// Simplest possible command to compute the Chicago problem
// for estimating calibrated trifocal geometry from points and lines at points
//
// This is to be kept very simple C with only minimal C++ with Templates.
// If you want to complicate this, please create another executable.
// 
int
main(int argc, char **argv)
{
  const char *input="stdin";
  const char *output="stdout";
  --argc;
  bool profile = false;   // run some default solves for profiling
    
  if (argc == 1) {
    if (std::string (argv[1]) == "-g" || std::string (argv[1]) == "--profile")
      profile = true;
  } else if (argc == 2) {
    input = argv[1];
    output = argv[2];
    stdio_ = false;
  } else if (argc != 0) {
    std::cerr << "minus: \033[1;91m error\e[m\n";
    print_usage();
  }

  #ifdef M_VERBOSE
  if (!profile) {
    std::cerr << "LOG Input being read from " << input << std::endl;
    std::cerr << "LOG Output being written to " << output << std::endl;
  } else
    std::cerr << "LOG Running default solve for profiling\n";
  #endif 

  if (!profile && !mread<Float>(input)) return 1; // reads into global params_
  
  static M::solution solutions[M::nsols];
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  #ifdef M_VERBOSE
  std::cerr << "LOG \033[0;33mStarting path tracker\e[m\n" << std::endl;
  #endif 
  //  unsigned retval = 
  //  ptrack(&MINUS_DEFAULT, start_sols_, params_, solutions);
  {
    #ifdef M_VERBOSE
    std::cerr << "LOG \033[0;33mUsing 4 threads by default\e[m\n" << std::endl;
    #endif 
    std::thread t[4];
    t[0] = std::thread(M::track, M::DEFAULT, start_sols_, params_, solutions, 0, 78);
    t[1] = std::thread(M::track, M::DEFAULT, start_sols_, params_, solutions, 78, 78*2);
    t[2] = std::thread(M::track, M::DEFAULT, start_sols_, params_, solutions, 78*2, 78*3);
    t[3] = std::thread(M::track, M::DEFAULT, start_sols_, params_, solutions, 78*3, 78*4);
    t[0].join(); t[1].join(); t[2].join(); t[3].join();
  }
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(t2 - t1).count();
  if (profile) {
    // compare solutions to certain values from M2
    // two random entries
    if (std::abs(solutions[1].x[1] - complex(-.25177177692982444e1, -.84845195030295639)) <= tol &&
        std::abs(solutions[M::nsols-2].x[2] - complex(.7318330016224166, .10129116603501138)) <= tol)
      std::cerr << "LOG solutions look OK\n";
    else  {
      std::cerr << "LOG \033[1;91merror:\e[m solutions dont match m2. Errors: ";
      std::cerr << std::abs(solutions[1].x[2] - complex(-.25177177692982444e1, -.84845195030295639)) << ", "
          << std::abs(solutions[M::nsols-2].x[2] - complex(.7318330016224166, .10129116603501138)) << std::endl;
    }
  }
  if (!mwrite<Float>(solutions, output)) return 2;
  #ifdef M_VERBOSE
  std::cerr << "LOG \033[1;32mTime of solver: " << duration << "ms\e[m" << std::endl;
  #endif
   
  return 0;
}
