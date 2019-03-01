#include <testlib/testlib_register.h>

DECLARE( test_minus );
DECLARE( test_expminus );

void
register_tests()
{
   REGISTER( test_minus );
   REGISTER( test_expminus );
}

DEFINE_MAIN;
