#ifndef debug_util_h_
#define debug_util_h_

#define LOG(msg) do { } while (0)

#ifdef M_VERBOSE
#include <iostream>
#endif

template <typename F>
inline void
pprint(const F *v, unsigned n, bool newline=false)
{
#ifdef M_VERBOSE
  for (unsigned i=0; i < n; ++i)
    std::cout << v[i] << ((newline)? "\n" : " ");
  std::cout << std::endl;
#endif
}

template <typename F>
inline void
pprint(const F *v, unsigned nrows, unsigned ncols)
{
#ifdef M_VERBOSE
  for (unsigned i=0; i < nrows; ++i)
    print(v + i*ncols, ncols);
#endif
}

//#ifndef NDEBUG
//
// LOG is a logging useful for the user
//
// LLOG is Lowlevel Logging - usually for telemetry or debugging
// 
#ifdef M_VERBOSE
  #undef LOG
  #define LOG(msg) do { \
    std::cerr << "LOG " << msg << std::endl; \
  } while(0)
  
 #define LLOG(x) std::cerr << x // usage: LLOG("Message " << value);
#else
 #define LLOG(x)
#endif // M_VERBOSE
// #endif // NDEBUG


#endif // debug_util_h_
