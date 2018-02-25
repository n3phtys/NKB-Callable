//
// Created by nephtys on 21.01.18.
//

#ifndef SCARLET_STUDY_NKB_CHOLSKY_H
#define SCARLET_STUDY_NKB_CHOLSKY_H

#include "./common.h"


void cholsky_test ( double *er, double *fp, double *tm );
void cholsky ( int ida, int nmat, int m, int n, double a[], int nrhs, int idb,
               double b[] );


/******************************************************************************/

void cholsky_test ( double *er, double *fp, double *tm )

/******************************************************************************/
/*
  Purpose:

    CHOLSKY_TEST tests CHOLSKY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 July 2011

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
# define IDA 250
# define M 4
# define N 40
# define NMAT 250
# define NRHS 3

    double a[(IDA+1)*(M+1)*(N+1)];
    double ans;
    double ax[(IDA+1)*(M+1)*(N+1)];
    double b[(NRHS+1)*(NMAT+1)*(N+1)];
    double bx[(NRHS+1)*(NMAT+1)*(N+1)];
    double f7;
    int i;
    int i1;
    int ida = IDA;
    int it;
    int j;
    int k;
    int la;
    int lb;
    int m = M;
    int n = N;
    int nmat = NMAT;
    int nrhs = NRHS;
    double t;
    double t30;
    double time1;

    it = 200;
    ans = 5177.88531774562;
    la = ( ida + 1 ) * ( m + 1 ) * ( n + 1 );
    lb = ( nrhs + 1 ) * ( nmat + 1 ) * ( n + 1 );
/*
  Random initialization.
*/
    f7 = 78125.0;
    t30 = 1073741824.0;
    t = f7 / t30;

    i1 = 0;

    for ( k = 0; k <= n; k++ )
    {
        for ( j = -m; j <= 0; j++ )
        {
            for ( i = 0; i <= ida; i++ )
            {
                t = fmod ( f7 * t, 1.0 );
                ax[i1] = t;
                i1 = i1 + 1;
            }
        }
    }

    i1 = 0;
    for ( k = 0; k <= n; k++ )
    {
        for ( j = 0; j <= nmat; j++ )
        {
            for ( i = 0; i <= nrhs; i++ )
            {
                t = fmod ( f7 * t, 1.0 );
                bx[i1] = t;
                i1 = i1 + 1;
            }
        }
    }
/*
  Timing.
*/
    time1 = wtime ( );

    for ( j = 1; j <= it; j++ )
    {
        r8vec_copy ( la, ax, a );
        r8vec_copy ( lb, bx, b );
        cholsky ( ida, nmat, m, n, a, nrhs, ida, b );
    }
    *tm = wtime ( ) - time1;
/*
  Results.
*/
    i1 = 1+19*(nrhs+1)+19*(nrhs+1)*(nmat+1);
    *er = r8_abs ( ( b[i1] - ans ) / ans );
    *fp = ( double ) ( it * ( nmat + 1 ) * 4403 );

    return;
# undef IDA
# undef M
# undef N
# undef NMAT
# undef NRHS
}
/******************************************************************************/

void  cholsky ( int ida, int nmat, int m, int n, double a[], int nrhs, int idb,
                double b[] )

/******************************************************************************/
/*
  Purpose:

    CHOLSKY carries out Cholesky decomposition and back substitution.

  Discussion:

    The Cholesky decomposition is performed on a set of input matrices
    which are provided as a single three-dimensional array.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
    double eps;
    double epss[(nmat+1)*(n+1)];
    int i;
    int i0;
    int i1;
    int i2;
    int i3;
    int j;
    int jj;
    int k;
    int l;

    eps = 1.0E-13;
/*
  Cholesky decomposition.
*/
    for ( j = 0; j <= n; j++ )
    {
        i0 = i4_max ( -m, -j );
/*
  Off diagonal elements.
*/
        for ( i = i0; i <= -1; i++ )
        {
            for ( jj = i0 - i; jj <= -1; jj++ )
            {
                for ( l = 0; l <= nmat; l++ )
                {
                    i1 = l + (i+m)    * (ida+1) + j     * (ida+1)*(m+1);
                    i2 = l + (jj+m)   * (ida+1) + (i+j) * (ida+1)*(m+1);
                    i3 = l + (i+jj+m) * (ida+1) + j     * (ida+1)*(m+1);
                    a[i1] = a[i1] - a[i2] * a[i3];
                }
            }
            for ( l = 0; l <= nmat; l++ )
            {
                i1 = l + (i+m) * (ida+1) +    j  * (ida+1)*(m+1);
                i2 = l +    m  * (ida+1) + (i+j) * (ida+1)*(m+1);
                a[i1] = a[i1] * a[i2];
            }
        }
/*
  Store inverse of diagonal elements.
*/
        for ( l = 0; l <= nmat; l++ )
        {
            i1 = l + + m * (ida+1) + j * (ida+1)*(m+1);
            epss[l] = eps * a[i1];
        }
        for ( jj = i0; jj <= -1; jj++ )
        {
            for ( l = 0; l <= nmat; l++ )
            {
                i1 = l +     m  * (ida+1) + j * (ida+1)*(m+1);
                i2 = l + (jj+m) * (ida+1) + j * (ida+1)*(m+1);
                a[i1] = a[i1] - pow ( a[i2], 2 );
            }
        }

        for ( l = 0; l <= nmat; l++ )
        {
            i1 = l + m * (ida+1) + j * (ida+1)*(m+1);
            a[i1] = 1.0 / sqrt ( r8_abs ( epss[l] + a[i1] ) );
        }
    }
/*
  Solution.
*/
    for ( i = 0; i <= nrhs; i++ )
    {
        for ( k = 0; k <= n; k++ )
        {
            for ( l = 0; l <= nmat; l++ )
            {
                i1 = i+l*(nrhs+1)+k*(nrhs+1)*(idb+1);
                i2 = l + m * (ida+1) + k * (ida+1)*(m+1);
                b[i1] = b[i1] * a[i2];
            }
            for ( jj = 1; jj <= i4_min ( m, n - k ); jj++ )
            {
                for ( l = 0; l <= nmat; l++ )
                {
                    i1 = i+l*(nrhs+1)+(k+jj)*(nrhs+1)*(idb+1);
                    i2 = l + (-jj+m) * (ida+1) + (k+jj) * (ida+1)*(m+1);
                    i3 = i+l*(nrhs+1)+k*(nrhs+1)*(idb+1);
                    b[i1] = b[i1] - a[i2] * b[i3];
                }
            }
        }

        for ( k = n; 0 <= k; k-- )
        {
            for ( l = 0; l <= nmat; l++ )
            {
                i1 = i + l*(nrhs+1) + k*(nrhs+1)*(idb+1);
                i2 = l + m * (ida+1) + k * (ida+1)*(m+1);
                b[i1] = b[i1] * a[i2];
            }
            for ( jj = 1; jj <= i4_min ( m, k ); jj++ )
            {
                for ( l = 0; l <= nmat; l++ )
                {
                    i1 = i+l*(nrhs+1)+(k-jj)*(nrhs+1)*(idb+1);
                    i2 = l + (-jj+m) * (ida+1) + k * (ida+1)*(m+1);
                    i3 = i+l*(nrhs+1)+k*(nrhs+1)*(idb+1);
                    b[i1] = b[i1] - a[i2] * b[i3];
                }
            }
        }
    }
    return;
}




#endif //SCARLET_STUDY_NKB_CHOLSKY_H
