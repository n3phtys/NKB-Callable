//
// Created by nephtys on 21.01.18.
//

#ifndef SCARLET_STUDY_NKB_VPENTA_H
#define SCARLET_STUDY_NKB_VPENTA_H

#include "./common.h"

void vpenta_test ( double *er, double *fp, double *tm );
void vpenta ( int n, double a[],
              double b[], double c[], double d[],  double e[], double f[], double x[],
              double y[] );



/******************************************************************************/

void vpenta_test ( double *er, double *fp, double *tm )

/******************************************************************************/
/*
  Purpose:

    VPENA_TEST tests VPENTA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
# define JL 0
# define JU 127
# define KL 0
# define KU 127
# define NJA 128
# define NJB 128

    double a[NJA*NJB];
    double ans;
    double b[NJA*NJB];
    double c[NJA*NJB];
    double d[NJA*NJB];
    double e[NJA*NJB];
    double f[NJA*NJB*3];
    double f7;
    double fx[NJA*NJB*3];
    int i;
    int it;
    int j;
    int jl = JL;
    int ju = JU;
    int k;
    int kl = KL;
    int ku = KU;
    int lf;
    int nja = NJA;
    int njb = NJB;
    double t;
    double t30;
    double time1;
    double x[NJA*NJB];
    double y[NJA*NJB];

    it = 400;
    ans = -0.354649411858726;
    lf = nja * njb * 3;
/*
  Random initialization.
*/
    f7 = 78125.0;
    t30 = 1073741824.0;
    t = f7 / t30;

    for ( j = kl; j <= ku; j++ )
    {
        for ( i = jl; i <= ju; i++ )
        {
            t = fmod ( f7 * t, 1.0 );
            a[i+j*nja] = t;
            t = fmod ( f7 * t, 1.0 );
            b[i+j*nja] = t;
            t = fmod ( f7 * t, 1.0 );
            c[i+j*nja] = t;
            t = fmod ( f7 * t, 1.0 );
            d[i+j*nja] = t;
            t = fmod ( f7 * t, 1.0 );
            e[i+j*nja] = t;
            for ( k = 0; k < 3; k++ )
            {
                t = fmod ( f7 * t, 1.0 );
                fx[i+j*nja+k*nja*njb] = t;
            }
        }
    }
/*
  Timing.
*/
    time1 = wtime ( );

    for ( i = 1; i <= it; i++ )
    {
        r8vec_copy ( lf, fx, f );
        vpenta ( JU, a, b, c, d, e, f, x, y );
    }

    *tm = wtime ( ) - time1;
/*
  Results.
*/
    *er = r8_abs ( ( f[18+18*nja+0*nja*njb] - ans ) / ans );
    *fp = ( double ) ( it * ku * ( 40 * ku - 53 ) );

    return;
# undef JL
# undef JU
# undef KL
# undef KU
# undef NJA
# undef NJB
}

/******************************************************************************/

void  vpenta ( int n, double a[],
               double b[], double c[], double d[],  double e[], double f[], double x[],
               double y[] )

/******************************************************************************/
/*
  Purpose:

    VPENTA inverts 3 pentadiagonal systems simultaneously.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
    int ju = n;
    int ku = n;
    int jl = 0;
    int kl = 0;
    int nja = n + 1;
    int njb = n + 1;


    int j;
    int jx;
    int k;
    double rld;
    double rld1;
    double rld2;
    double rldi;
/*
  Start forward generation process and sweep.
*/
    j = jl;
    for ( k = kl; k <= ku; k++ )
    {
        rld = c[j+k*nja];
        rldi = 1.0 / rld;
        f[j+k*nja+0*nja*njb] = f[j+k*nja+0*nja*njb] * rldi;
        f[j+k*nja+1*nja*njb] = f[j+k*nja+1*nja*njb] * rldi;
        f[j+k*nja+2*nja*njb] = f[j+k*nja+2*nja*njb] * rldi;
        x[j+k*nja] = d[j+k*nja] * rldi;
        y[j+k*nja] = e[j+k*nja] * rldi;
    }

    j = jl + 1;
    for ( k = kl; k <= ku; k++ )
    {
        rld1 = b[j+k*nja];
        rld = c[j+k*nja] - rld1 * x[j-1+k*nja];
        rldi = 1.0 / rld;
        f[j+k*nja+0*nja*njb] = ( f[j+k*nja+0*nja*njb]
                                 - rld1 * f[j-1+k*nja+0*nja*njb] ) * rldi;
        f[j+k*nja+1*nja*njb] = ( f[j+k*nja+1*nja*njb]
                                 - rld1 * f[j-1+k*nja+1*nja*njb] ) * rldi;
        f[j+k*nja+2*nja*njb] = ( f[j+k*nja+2*nja*njb]
                                 - rld1 * f[j-1+k*nja+2*nja*njb] ) * rldi;
        x[j+k*nja] = ( d[j+k*nja] - rld1 * y[j-1+k*nja] ) * rldi;
        y[j+k*nja] = e[j+k*nja] * rldi;
    }

    for ( j = jl + 2; j <= ju - 2; j++ )
    {
        for ( k = kl; k <= ku; k++ )
        {
            rld2 = a[j+k*nja];
            rld1 = b[j+k*nja] - rld2 * x[j-2+k*nja];
            rld = c[j+k*nja] - ( rld2 * y[j-2+k*nja] + rld1 * x[j-1+k*nja] );
            rldi = 1.0 / rld;
            f[j+k*nja+0*nja*njb] = ( f[j+k*nja+0*nja*njb]
                                     - rld2 * f[j-2+k*nja+0*nja*njb] - rld1 * f[j-1+k*nja+0*nja*njb] ) * rldi;
            f[j+k*nja+1*nja*njb] = ( f[j+k*nja+1*nja*njb]
                                     - rld2 * f[j-2+k*nja+1*nja*njb] - rld1 * f[j-1+k*nja+1*nja*njb] ) * rldi;
            f[j+k*nja+2*nja*njb] = ( f[j+k*nja+2*nja*njb]
                                     - rld2 * f[j-2+k*nja+2*nja*njb] - rld1 * f[j-1+k*nja+2*nja*njb] ) * rldi;
            x[j+k*nja] = ( d[j+k*nja] - rld1 * y[j-1+k*nja] ) * rldi;
            y[j+k*nja] = e[j+k*nja] * rldi;
        }
    }

    j = ju - 1;
    for ( k = kl; k <= ku; k++ )
    {
        rld2 = a[j+k*nja];
        rld1 = b[j+k*nja] - rld2 * x[j-2+k*nja];
        rld = c[j+k*nja] - ( rld2 * y[j-2+k*nja] + rld1 * x[j-1+k*nja] );
        rldi = 1.0 / rld;;
        f[j+k*nja+0*nja*njb] = ( f[j+k*nja+0*nja*njb]
                                 - rld2 * f[j-2+k*nja+0*nja*njb] - rld1 * f[j-1+k*nja+0*nja*njb] ) * rldi;
        f[j+k*nja+1*nja*njb] = ( f[j+k*nja+1*nja*njb]
                                 - rld2 * f[j-2+k*nja+1*nja*njb] - rld1 * f[j-1+k*nja+1*nja*njb] ) * rldi;
        f[j+k*nja+2*nja*njb] = ( f[j+k*nja+2*nja*njb]
                                 - rld2 * f[j-2+k*nja+2*nja*njb] - rld1 * f[j-1+k*nja+2*nja*njb] ) * rldi;
        x[j+k*nja] = ( d[j+k*nja] - rld1 * y[j-1+k*nja] ) * rldi;
    }

    j = ju;
    for ( k = kl; k <= ku; k++ )
    {
        rld2 = a[j+k*nja];
        rld1 = b[j+k*nja] - rld2 * x[j-2+k*nja];
        rld = c[j+k*nja] - ( rld2 * y[j-2+k*nja] + rld1 * x[j-1+k*nja] );
        rldi = 1.0 / rld;
        f[j+k*nja+0*nja*njb] = ( f[j+k*nja+0*nja*njb]
                                 - rld2 * f[j-2+k*nja+0*nja*njb] - rld1 * f[j-1+k*nja+0*nja*njb] ) * rldi;
        f[j+k*nja+1*nja*njb] = ( f[j+k*nja+1*nja*njb]
                                 - rld2 * f[j-2+k*nja+1*nja*njb] - rld1 * f[j-1+k*nja+1*nja*njb] ) * rldi;
        f[j+k*nja+2*nja*njb] = ( f[j+k*nja+2*nja*njb]
                                 - rld2 * f[j-2+k*nja+2*nja*njb] - rld1 * f[j-1+k*nja+2*nja*njb] ) * rldi;
    }
/*
  Back sweep solution.
*/
    for ( k = kl; k <= ku; k++ )
    {
        f[ju+k*nja+0*nja*njb] = f[ju+k*nja+0*nja*njb];
        f[ju+k*nja+1*nja*njb] = f[ju+k*nja+1*nja*njb];
        f[ju+k*nja+2*nja*njb] = f[ju+k*nja+2*nja*njb];
        f[ju-1+k*nja+0*nja*njb] = f[ju-1+k*nja+0*nja*njb]
                                  - x[ju-1+k*nja] * f[ju+k*nja+0*nja*njb];
        f[ju-1+k*nja+1*nja*njb] = f[ju-1+k*nja+1*nja*njb]
                                  - x[ju-1+k*nja] * f[ju+k*nja+1*nja*njb];
        f[ju-1+k*nja+2*nja*njb] = f[ju-1+k*nja+2*nja*njb]
                                  - x[ju-1+k*nja] * f[ju+k*nja+2*nja*njb];
    }

    for ( j = 2; j <= ju - jl; j++ )
    {
        jx = ju - j;
        for ( k = kl; k <= ku; k++ )
        {
            f[jx+k*nja+0*nja*njb] = f[jx+k*nja+0*nja*njb]
                                    - x[jx+k*nja] * f[jx+1+k*nja+0*nja*njb]
                                    - y[jx+k*nja] * f[jx+2+k*nja+0*nja*njb];
            f[jx+k*nja+1*nja*njb] = f[jx+k*nja+1*nja*njb]
                                    - x[jx+k*nja] * f[jx+1+k*nja+1*nja*njb]
                                    - y[jx+k*nja] * f[jx+2+k*nja+1*nja*njb];
            f[jx+k*nja+2*nja*njb] = f[jx+k*nja+2*nja*njb]
                                    - x[jx+k*nja] * f[jx+1+k*nja+2*nja*njb]
                                    - y[jx+k*nja] * f[jx+2+k*nja+2*nja*njb];
        }
    }
    return;
}



#endif //SCARLET_STUDY_NKB_VPENTA_H
