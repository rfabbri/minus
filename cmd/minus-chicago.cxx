#include <iostream>

#define M_VERBOSE 1

static void
print_usage()
{
  std::cerr << "Usage: minus-chicago [input solutions]\n\n";
  std::cerr << "If no argument is given, 'input' is assumed stdin,\n\
  'solutions' will be output to stdout\n";
  exit(1);
}

// Simplest possible command to compute the Chicago problem
// for estimating calibrated trifocal geometry from points and lines at points
//
// This is to be kept very simple with only standard C.
// If you want to complicate this, please create another executable.
// 
int
main(int argc, char **argv)
{
  const char *input="stdin";
  const char *output="stdout";
  --argc;
    
  if (argc == 2) {
    input = argv[1];
    output = argv[2];
  } else if (argc != 0) {
    std::cerr << "minus-chicago: \033[1;91m error\e[m\n";
    print_usage();
  }

#ifdef M_VERBOSE
  std::cerr << "LOG: Input being read from " << input << std::endl;
  std::cerr << "LOG: Output being written to " << output << std::endl;
#endif 
  
  return 0;
}
