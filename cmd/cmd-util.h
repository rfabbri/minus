#ifndef cmd_util_h_
#define cmd_util_h_
// 
// Header to be included in minus-*.cxx commands, after the appropriate defines
// 

#include <minus/internal-util.h>

namespace MiNuS {
  
// class with the maximum that is completely universal to all commands/problems 
template <typename F>
struct minus_cmd_io {
  bool stdio_ = true;  // by default read/write from stdio
  std::ifstream infp_;
  const char *input_ = "stdin";
  const char *output_ = "stdout";
  
  // ----------------------------------------------------------------------------
  
  // Output solutions in ASCII matlab format
  //
  // ---------------------------------------------------------
  // If in the future our solver is really fast, we may need Binary IO:
  // complex solutions[NSOLS*NVE];
  // 
  // To read this output file in matlab, do:
  // 
  // fid = fopen(fname,'r');
  // a_raw = fread(fid,'double');
  // fclose(fid);
  // 
  // % Reshape a to have proper real and imaginary parts
  // a = a_raw(1:2:end) + i*a_raw(2:2:end);
  // 
  bool
  mwrite(const M::solution s[M::nsols], const char *fname) const
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
  
  // Try to read n elements, filling in p in row-major order.
  static bool
  read_block(std::istream &in, F *p, unsigned n)
  {
    LOG("reading");
    const F *end = p + n;
    while (!in.eof() && p != end) {
        try {
          in >> *p++;
  //        std::cerr << *(p-1) << std::endl;
          if (in.eof()) {
            std::cerr << "I/O Error: Premature input termination\n";
            return false;
          }
        } catch (std::istream::failure &E) {
          std::cerr << "I/O Error: Invalid input conversion or other error\n";
          return false;
        }
    }
    if (p != end) {
      std::cerr << "I/O Premature input termination\n";
      return false;
    }
    return true;
  }
  
  bool
  init_input(const char *fname, std::istream *inp)
  {
    if (!stdio_) {
      infp_.open(fname, std::ios::in);
      if (!infp_) {
        std::cerr << "I/O Error opening input " << fname << std::endl;
        return false;
      }
      inp = &infp_;
    }
    inp->exceptions(std::istream::failbit | std::istream::badbit);
    return true;
  }
  
  // Reads directly into the global variable params_
  // if you want to build it or randomize it directly in external software
  // 
  // The file format is the concatenated parameters of the start system
  // and of the target system, optionally with some randomization to improve conditioning.
  // This is just like P01 variable in, e.g., minus/tutorial/linecircle/end-linecircle.m2
  // The file format does not include an explicit imaginary 'i' string,
  // it is just the real and imaginary parts in sequence, separated by space or
  // newline at will.
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
  // be it point-tangents as in Ric's differential-geometric format, be it a
  // linecomplex as in the visible lines formulation
  //
  // This format is generic enough to be adapted to M2 or matlab
  //
  bool
  mread(std::istream &in) const
  {
    F *dparams = (F *)data::params_;
    while (!in.eof() && dparams != (F *)data::params_+2*2*M::f::nparams) {
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
    if (dparams != (F *)data::params_+2*2*M::f::nparams)
      std::cerr << "I/O Premature input termination\n";
  //  for (unsigned i=0; i < 2*NPARAMS; ++i)
  //    std::cerr << "D " << params_[i] << std::endl;
    return true;
  }
};

void
print_num_steps(M::solution solutions[M::nsols])
{
  LOG("solution id x num steps:");
  unsigned sum=0;
  for (unsigned s=0; s < M::nsols; ++s) {
    LOG(s << " " << solutions[s].num_steps);
    sum += solutions[s].num_steps;
  }
  LOG("total number of steps: " << sum);
}

void 
print_settings(const M::track_settings &settings) { // print to stderr
  #ifdef M_VERBOSE
  {
  std::cerr << "Track settings ------------------------------------------------\n";
  const char *names[10] = {
    "init_dt_",
    "min_dt_",
    "end_zone_factor_",
    "epsilon2_",
    "dt_increase_factor_",
    "dt_decrease_factor_",
    "infinity_threshold2_",
    "max_num_steps_",
    "num_successes_before_increase_",
    "max_corr_steps_"
  };
  Float *ptr = (Float *) &settings;
  for (int i=0; i < 7; ++i)
    std::cerr << names[i] << " = " << *ptr++ << std::endl;
  std::cerr << names[7] << " = " << settings.max_num_steps_ << std::endl;
  std::cerr << names[8] << " = " << (int)settings.num_successes_before_increase_ << std::endl;
  std::cerr << names[9] << " = " << (int)settings.max_corr_steps_ << std::endl;
  }
  #endif
}

} // namespace minus

#endif // cmd_util_h_
