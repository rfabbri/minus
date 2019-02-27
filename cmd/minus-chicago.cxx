#include <iostream>

#define M_VERBOSE 1

static void
print_usage()
{
  std::cerr << "Usage: minus-chicago [input solutions]\n";
  std::cerr << "If no arg is given, 'input' is assumed input.txt,\n\
  'solutions' will be output to solutions.txt\n";
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
  const char *input="input.txt";
  const char *output="solutions.txt";
  --argc;
    
  if (argc == 2) {
    input = argv[1];
    output = argv[2];
  } else if (argc != 0) {
    std::cout << "minus-chicago: \033[1;91m error\e[m\n";
    print_usage();
  }

#ifdef M_VERBOSE
  std::cout << "LOG: Input being read from " << input << std::endl;
  std::cout << "LOG: Output being written to " << output << std::endl;
#endif 
  
  return 0;
}
