# include <stdio.h>
#include "../include/common.h"
#include "../include/nkb_btrix.h"
#include "../include/nkb_cfft2d.h"
#include "../include/nkb_cholsky.h"
#include "../include/nkb_emit.h"
#include "../include/nkb_gmtry.h"
#include "../include/nkb_vpenta.h"
#include "../include/nkb_mxm.h"

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for NAS.

  Discussion:

    This is a version of the NAS kernel benchmark program,
    whose original version was created by David Bailey,
    dated 17 December 1984.

    Each of the tests begins by filling arrays with pseudorandom values
    generated by the recursion:
      x(n+1) = 5^7 * x(n)  (mod 2^30)
    This recursion will generate 2^28 (approx. 268 million) numbers
    before repeating.  For this scheme to work properly, the hardware
    multiply operation must be correct to 47 bits of precision.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 July 2015

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
    double er = 0;
    double er_total;
    double fp = 0;
    double fp_total;
    int i;
    char pn[9];
    double rt;
    double tm = 0;
    double tm_total;

    er_total = 0.0;
    fp_total = 0.0;
    tm_total = 0.0;

    timestamp ( );
    printf ( "\n" );
    printf ( "NAS:\n" );
    printf ( "  C version\n" );
    printf ( "\n" );
    printf ( "                The NAS kernel benchmark program\n" );
    printf ( "\n" );
    printf ( " Program        Error          FP Ops" );
    printf ( "        Seconds     MFLOPS\n" );
    printf ( "\n" );

    for ( i = 1; i <= 7; i++ )
    {
        if ( i == 1 )
        {
            strcpy ( pn, "BTRIX   " );
            btrix_test ( &er, &fp, &tm );
        }
        else if ( i == 2 )
        {
            strcpy ( pn, "CFFT2D  " );
            cfft2d_test ( &er, &fp, &tm );
        }
        else if ( i == 3 )
        {
            strcpy ( pn, "CHOLSKY " );
            cholsky_test ( &er, &fp, &tm );
        }
        else if ( i == 4 )
        {
            strcpy ( pn, "EMIT    " );
            emit_test ( &er, &fp, &tm );
        }
        else if ( i == 5 )
        {
            strcpy ( pn, "GMTRY   " );
            gmtry_test ( &er, &fp, &tm );
        }
        else if ( i == 6 )
        {
            strcpy ( pn, "MXM     " );
            mxm_test ( &er, &fp, &tm );
        }
        else if ( i == 7 )
        {
            strcpy ( pn, "VPENTA  " );
            vpenta_test ( &er, &fp, &tm );
        }
        rt = 1.0E-06 * fp / tm;
        printf ( " %s  %13.4g  %13.4g  %10.4f  %10.2f\n",
                 pn, er, fp, tm, rt );

        er_total = er_total + er;
        fp_total = fp_total + fp;
        tm_total = tm_total + tm;
    }

    strcpy ( pn, "Total   " );
    rt = 1.0E-06 * fp_total / tm_total;
    printf ( "\n" );
    printf ( " %s  %13.4g  %13.4g  %10.4f  %10.2f\n",
             pn, er_total, fp_total, tm_total, rt );

    printf ( "\n" );
    printf ( "NAS:\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    timestamp ( );

    return 0;
}
