#ifndef debug_util_h_
#define debug_util_h_

#include <mutex>

template <typename F>
inline void
print(const F *v, unsigned n, bool newline=false)
{
  static std::mutex lock;
  const std::lock_guard<std::mutex> guard(lock);

  for (unsigned i=0; i < n; ++i)
    std::cout << v[i] << ((newline)? "\n" : " ");
  std::cout << std::endl;
}

template <typename F>
inline void
print(const F *v, unsigned nrows, unsigned ncols)
{
  for (unsigned i=0; i < nrows; ++i)
    print(v + i*ncols, ncols);
}

#define LOG(msg) do { } while (0)

//#ifndef NDEBUG
#ifdef M_VERBOSE
#undef LOG
#define LOG(msg) do { \
  std::cerr << "LOG " << msg << std::endl; \
} while(0)
#endif // M_VERBOSE
// #endif // NDEBUG

#endif // debug_util_h_
