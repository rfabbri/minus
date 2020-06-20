#include <testlib/testlib_register.h>

DECLARE( test_minus );
DECLARE( test_minus_cleveland );
DECLARE( test_internals );
DECLARE( test_internals_cleveland );

void
register_tests()
{
  REGISTER( test_minus );
  REGISTER( test_minus_cleveland );
  REGISTER( test_internals );
  REGISTER( test_internals_cleveland );
}

DEFINE_MAIN;
