#include <testlib/testlib_register.h>

DECLARE( test_minus );
DECLARE( test_internals );

void
register_tests()
{
  REGISTER( test_minus );
  REGISTER( test_internals );
}

DEFINE_MAIN;
