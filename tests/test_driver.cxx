#include <testlib/testlib_register.h>

DECLARE( test_minus );
DECLARE( test_expminus );
DECLARE( test_linear );

void
register_tests()
{
   REGISTER( test_minus );
   REGISTER( test_expminus );
   REGISTER( test_linear );
}

DEFINE_MAIN;
