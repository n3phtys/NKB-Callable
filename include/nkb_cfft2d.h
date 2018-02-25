//
// Created by nephtys on 21.01.18.
//

#ifndef SCARLET_STUDY_NKB_CFFT2D_H
#define SCARLET_STUDY_NKB_CFFT2D_H


#include "custom_complex.h"


void cfft2d_test ( double *er, double *fp, double *tm );
void cfft2d1 ( int is, int m, int m1, int n, complex_number_t x[],
               complex_number_t w[], int ip[] );
void cfft2d2 ( int is, int m, int m1, int n, complex_number_t x[],
               complex_number_t w[], int ip[] );




/******************************************************************************/

void cfft2d_test ( double *er, double *fp, double *tm )

/******************************************************************************/
/*
  Purpose:

    CFFT2D_TEST is the test program for CFFT2D1 and CFFTD2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
# define M 128
# define M1 128
# define N 256

    double ans;
    double f7;
    int i;
//
//  IP must be dimensions 2*MAX(M,N)
//
    int ip[2*N];
    int it;
    int j;
    int k;
    int m = M;
    int m1 = M1;
    int n = N;
    double rmn;
    double t1;
    double t2;
    double t30;
    double time1;
    complex_number_t w1[M];
    complex_number_t w2[N];
    complex_number_t x[M1*N];


    it = 100;
    ans = 0.894799941219277;
    rmn = 1.0 / ( double ) ( m * n );
/*
  Random initialization.
*/
    f7 = 78125.0;
    t30 = 1073741824.0;
    t2 = f7 / t30;

    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < m; i++ )
        {
            t1 = fmod ( f7 * t2, 1.0 );
            t2 = fmod ( f7 * t1, 1.0 );
            x[i+j*m1] = complexify(t1 , t2);
        }
    }
    cfft2d1 ( 0, m, m1, n, x, w1, ip );
    cfft2d2 ( 0, m, m1, n, x, w2, ip );
/*
  Timing.
*/
    time1 = wtime ( );

    for ( k = 1; k <= it; k++ )
    {
        for ( j = 0; j < n; j++ )
        {
            for ( i = 0; i < m; i++ )
            {
                x[i+j*m1] = cmul(x[i+j*m1], rmn);
            }
        }

        cfft2d1 ( 1, m, m1, n, x, w1, ip );
        cfft2d2 ( 1, m, m1, n, x, w2, ip );
        cfft2d2 ( -1, m, m1, n, x, w2, ip );
        cfft2d1 ( -1, m, m1, n, x, w1, ip );
    }

    *tm = wtime ( ) - time1;
/*
  Results.
*/
    *er = fabs ( ( complex_real ( x[18+18*m1] ) - ans ) / ans );
    *fp = ( double ) ( it * m * n ) * ( 2.0
                                        + 10.0 * log ( ( double ) ( m * n ) ) / log ( 2.0 ) );

    return;
# undef M
# undef M1
# undef N
}
/******************************************************************************/

void  cfft2d1 ( int is, int m, int m1, int n, complex_number_t x[],
                complex_number_t w[], int ip[] )

/******************************************************************************/
/*
  Purpose:

    CFFT2D1 performs _Complex radix 2 FFT''s on the first dimension.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
    complex_number_t ct;
    complex_number_t cx;
    int i;
    int i1;
    int i2;
    int ii;
    int im;
    int j;
    int k;
    int l;
    int m2;
    double pi;
    double t;

    pi = 3.141592653589793;
/*
  If IS = 0 then initialize only.
*/
    m2 = m / 2;
    if ( is == 0 )
    {
        for ( i = 0; i < m2; i++ )
        {
            t = 2.0 * pi * ( double ) ( i ) / ( double ) ( m );
            w[i] = complexify(cos ( t ) , sin ( t ) );
        }
        return;
    }
/*
  Perform forward or backward FFT''s according to IS = 1 or -1.
*/
    for ( i = 0; i < m; i++ )
    {
        ip[0+i*2] = i + 1;
    }

    l = 1;
    i1 = 1;

    for ( ; ; )
    {
        i2 = 3 - i1;
        for ( j = l; j <= m2; j = j + l )
        {
            cx = w[j-l+1-1];
            if ( is < 0 )
            {
                cx = complex_conj ( cx );
            }
            for ( i = j - l + 1; i <= j; i++ )
            {
                ii = ip[i1-1+(i-1)*2];
                ip[i2-1+(i+j-l-1)*2] = ii;
                im = ip[i1-1+(i+m2-1)*2];
                ip[i2-1+(i+j-1)*2] = im;
                for ( k = 1; k <= n; k++ )
                {
                    ct = csub(x[ii-1+(k-1)*m1] , x[im-1+(k-1)*m1]);
                    x[ii-1+(k-1)*m1] = cadd(x[ii-1+(k-1)*m1] , x[im-1+(k-1)*m1]);
                    x[im-1+(k-1)*m1] = cmulc(ct , cx);
                }
            }
        }

        l = 2 * l;
        i1 = i2;

        if ( m2 < l )
        {
            break;
        }
    }

    for ( i = 1; i <= m; i++ )
    {
        ii = ip[i1-1+(i-1)*2];
        if ( i < ii )
        {
            for ( k = 1; k <= n; k++ )
            {
                ct = x[i-1+(k-1)*m1];
                x[i-1+(k-1)*m1] = x[ii-1+(k-1)*m1];
                x[ii-1+(k-1)*m1] = ct;
            }
        }
    }

    return;
}
/******************************************************************************/

void  cfft2d2 ( int is, int m, int m1, int n, complex_number_t x[],
                complex_number_t w[], int ip[] )

/******************************************************************************/
/*
  Purpose:

    CFFT2D2 performs _Complex radix 2 FFT''s on the second dimension.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
    complex_number_t ct;
    complex_number_t cx;
    int i;
    int i1;
    int i2;
    int ii;
    int im;
    int j;
    int k;
    int l;
    int n2;
    double pi;
    double t;

    pi = 3.141592653589793;
/*
  If IS = 0, then initialize only.
*/
    n2 = n / 2;
    if ( is == 0 )
    {
        for ( i = 0; i < n2; i++ )
        {
            t = 2.0 * pi * ( double ) ( i ) / ( double ) ( n );
            w[i] = complexify(cos ( t ) , sin ( t ) );
        }
        return;
    }
/*
  Perform forward or backward FFT''s according to IS = 1 or -1.
*/
    for ( i = 0; i < n; i++ )
    {
        ip[0+i*2] = i + 1;
    }

    l = 1;
    i1 = 1;

    for ( ; ; )
    {
        i2 = 3 - i1;

        for ( j = l; j <= n2; j = j + l )
        {
            cx = w[j-l+1-1];
            if ( is < 0 )
            {
                cx = complex_conj ( cx );
            }

            for ( i = j - l + 1; i <= j; i++ )
            {
                ii = ip[i1-1+(i-1)*2];
                ip[i2-1+(i+j-l-1)*2] = ii;
                im = ip[i1-1+(i+n2-1)*2];
                ip[i2-1+(i+j-1)*2] = im;
                for ( k = 1; k <= m; k++ )
                {
                    ct = csub(x[k-1+(ii-1)*m1] , x[k-1+(im-1)*m1]);
                    x[k-1+(ii-1)*m1] = cadd(x[k-1+(ii-1)*m1] , x[k-1+(im-1)*m1]);
                    x[k-1+(im-1)*m1] = cmulc(ct , cx);
                }
            }
        }

        l = 2 * l;
        i1 = i2;

        if ( n2 < l )
        {
            break;
        }

    }
    for ( i = 1; i <= n; i++ )
    {
        ii = ip[i1-1+(i-1)*2];
        if ( i < ii )
        {
            for ( k = 1; k <= m; k++ )
            {
                ct = x[k-1+(i-1)*m1];
                x[k-1+(i-1)*m1] = x[k-1+(ii-1)*m1];
                x[k-1+(ii-1)*m1] = ct;
            }
        }
    }
    return;
}


#endif //SCARLET_STUDY_NKB_CFFT2D_H
