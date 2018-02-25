//
// Created by nephtys on 21.01.18.
//

#ifndef SCARLET_STUDY_NKB_MXM_H
#define SCARLET_STUDY_NKB_MXM_H

#include "./common.h"




void mxm_test ( double *er, double *fp, double *tm );
void mxm ( double a[], double b[], double c[], int l, int m, int n );

/******************************************************************************/

void  mxm_test ( double *er, double *fp, double *tm )

/******************************************************************************/
/*
  Purpose:

    MXM_TEST tests MXM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
# define L 256
# define M 128
# define N 64

    double a[L*M];
    double ans;
    double b[M*N];
    double c[L*N];
    double f7;
    int i;
    int ii;
    int it;
    int j;
    int l = L;
    int m = M;
    int n = N;
    double t;
    double t30;
    double time1;

    it = 100;
    ans = 35.2026179738722;
/*
  Random initialization.
*/
    f7 = 78125.0;
    t30 = 1073741824.0;
    t = f7 / t30;

    for ( j = 0; j < m; j++ )
    {
        for ( i = 0; i < l; i++ )
        {
            t = fmod ( f7 * t, 1.0 );
            a[i+j*l] = t;
        }
    }

    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < m; i++ )
        {
            t = fmod ( f7 * t, 1.0 );
            b[i+j*m] = t;
        }
    }
/*
  Timing.
*/
    time1 = wtime ( );

    for ( ii = 1; ii <= it; ii++ )
    {
        mxm ( a, b, c, l, m, n );
    }
    *tm = wtime ( ) - time1;
/*
  Results.
*/
    *er = r8_abs ( ( c[18+18*l] - ans ) / ans );
    *fp = ( double ) ( 2 * it * l * m * n );

    return;
# undef L
# undef M
# undef N
}
/******************************************************************************/

void  mxm ( double a[], double b[], double c[], int l, int m, int n )

/******************************************************************************/
/*
  Purpose:

    MXM computes the matrix product C = A * B.

  Discussion:

    The function uses 4-way unrolled loops to carry out matrix multiplication.

    M must be a multiple of 4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{

    for (int k = 0; k < n; k++ )
    {
        for (int i = 0; i < l; i++ )
        {
            c[i+(k*l)] = 0.0;
        }
    }

    for (int j = 0; j < m; j += 4 )
    {
        for (int k = 0; k < n; k++ )
        {
            for (int i = 0; i < l; i++ )
            {
                c[i+(k*l)] = c[i+(k*l)]
                             + a[i+(j    *l)] * b[j  +(k*m)]
                             + a[i+((j+1)*l)] * b[j+1+(k*m)]
                             + a[i+((j+2)*l)] * b[j+2+(k*m)]
                             + a[i+((j+3)*l)] * b[j+3+(k*m)];
            }
        }
    }
    return;
}


#endif //SCARLET_STUDY_NKB_MXM_H
