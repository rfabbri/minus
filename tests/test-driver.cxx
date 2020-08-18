#include <testlib/testlib_register.h>

DECLARE( test_minus );
DECLARE( test_minus_cleveland );
DECLARE( test_minus_2v5p );
DECLARE( test_internals );
DECLARE( test_internals_cleveland );
DECLARE( test_internals_2v5p );

void
register_tests()
{
  REGISTER( test_minus );
  REGISTER( test_minus_cleveland );
  REGISTER( test_minus_2v5p );
  REGISTER( test_internals );
  REGISTER( test_internals_cleveland );
  REGISTER( test_internals_2v5p );
}

DEFINE_MAIN;
