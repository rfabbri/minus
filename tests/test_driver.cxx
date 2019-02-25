#include <testlib/testlib_register.h>

DECLARE( test_eno_poly );
DECLARE( test_eno_interp);
DECLARE( test_eno_zerox );
DECLARE( test_eno_1d );
DECLARE( test_eno_shock_1d );
DECLARE( test_eno_image );
DECLARE( test_eno_third_order );
DECLARE( test_brent_root );
DECLARE( test_vnag );

void
register_tests()
{
   REGISTER( test_eno_poly );
   REGISTER( test_eno_interp );
   REGISTER( test_eno_zerox );
   REGISTER( test_eno_1d );
   REGISTER( test_eno_shock_1d );
   REGISTER( test_eno_image );
   REGISTER( test_eno_third_order );
   REGISTER( test_vnag );
}

DEFINE_MAIN;
