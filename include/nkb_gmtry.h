//
// Created by nephtys on 21.01.18.
//

#ifndef SCARLET_STUDY_NKB_GMTRY_H
#define SCARLET_STUDY_NKB_GMTRY_H

#include "./common.h"
#include "custom_complex.h"


void gmtry_test ( double *er, double *fp, double *tm );
void gmtry ( int nb, int nw, int nwall[], complex_number_t proj[], double rmatrx[],
             complex_number_t wall[], double xmax[], complex_number_t zcr[] );


/******************************************************************************/

void gmtry_test ( double *er, double *fp, double *tm )

/******************************************************************************/
/*
  Purpose:

    GMTRY_TEST tests GMTRY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
# define NB 5
# define NW 100

    double ans;
    double f7;
    int i;
    int it;
    int j;
    //int lw;
    int nb = NB;
    int nw = NW;
    int nwall[NB];
    complex_number_t proj[NW*NB];
    double rmatrx[NW*NB*NW*NB];
    double t1;
    double t2;
    double t30;
    double time1;
    complex_number_t wall[NW*NB];
    double xmax[NB];
    //complex_number_t z1;
    complex_number_t zcr[NW*NB];
    //complex_number_t zi;
    //complex_number_t zz;

    it = 2;
    ans = -2.57754233214174;
    //lw = 2 * nw * nb;
/*
  Random initialization.
*/
    f7 = 78125.0;
    t30 = 1073741824.0;
    t2 = f7 / t30;

    for ( j = 1; j <= nb; j++ )
    {
        nwall[j-1] = nw;
    }

    for ( j = 1; j <= nb; j++ )
    {
        for ( i = 1; i <= nw; i++ )
        {
            t1 = fmod ( f7 * t2, 1.0 );
            t2 = fmod ( f7 * t1, 1.0 );
            wall[i-1+(j-1)*NW] = complexify(t1, t2);
        }
    }
/*
  Timing.
*/
    time1 = wtime ( );

    for ( i = 1; i <= it; i++ )
    {
        gmtry ( nb, nw, nwall, proj, rmatrx, wall, xmax, zcr );
    }

    *tm = wtime ( ) - time1;
/*
  Results.
*/
    *er = r8_abs ( ( rmatrx[18+18*nw*nb] - ans ) / ans );
    *fp = ( double ) ( it ) * ( ( double ) ( 120 * ( nb * nw * nb * nw ) )
                                +  0.666 * ( double ) ( nb * nw * nb * nw * nb * nw ) );

    return;
# undef NB
# undef NW
}
/******************************************************************************/

void  gmtry ( int nb, int nw, int nwall[], complex_number_t proj[], double rmatrx[],
              complex_number_t wall[], double xmax[], complex_number_t zcr[] )

/******************************************************************************/
/*
  Purpose:

    GMTRY computes solid-related arrays.

  Discussion:

    This function was extracted from a vortex method program.
    It sets up arrays needed for the computation, and performs
    Gauss elimination on the matrix of wall-influence coefficients.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
    double arcl;
    double dum;
    int i;
    int i0;
    int j;
    int j0;
    int k;
    int k1;
    int k2;
    int kp;
    int kron;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
    int ks;
#pragma clang diagnostic pop
    int l;
    int l1;
    int l2;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
    int ls;
#pragma clang diagnostic pop
    int matdim;
    double period;
    double pi;
    double pidp;
    double r0;
    double sig2;
    double sigma;
    double ylimit;
    double ymax;
    double ymin;
    complex_number_t z1;
    complex_number_t zi;
    complex_number_t zz;

    pi = 3.141592653589793;
    period = 3.0;
/*
  Compute arclength.
*/
    matdim = 0;
    arcl = 0.0;
    ymin = 1.0E+30;
    ymax = -1.0E+30;
    pidp = pi / period;

    for ( l = 1; l <= nb; l++ )
    {
        matdim = matdim + nwall[l-1];
        for ( k = 1; k <= nwall[l-1]; k++ )
        {
            k2 = 1 + ( k % nwall[l-1] );
            arcl = arcl + complex_abs ( csub(wall[k-1+(l-1)*nw] , wall[k2-1+(l-1)*nw] ));
        }
    }
/*
  Compute core radius.
*/
    r0 = 0.5 * arcl / ( double ) ( matdim );
    sigma = r0 / 2.0;
/*
  Define creation points.
*/
    for ( l = 1; l <= nb; l++ )
    {
        for ( k = 1; k <= nwall[l-1]; k++ )
        {
            k1 = 1 + ( k + nwall[l-1] - 2 ) % ( nwall[l-1] );
            k2 = 1 + k % nwall[l-1];
            zz = csub(wall[k1-1+(l-1)*nw] , wall[k2-1+(l-1)*nw]);
            zcr[k-1+(l-1)*nw] = cadd(wall[k-1+(l-1)*nw] , cmulc(cdiv(complexify(0.0, r0) , complex_abs ( zz )) , zz));
        }
/*
  Check that wall and creation points are not crossed due to
  too sharp a concave kink or an error in defining the body.
  Also find highest, lowest and right-most point.
*/
        xmax[l-1] = complex_real ( zcr[0+(l-1)*nw] );
#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnusedValue"
        ls = 0;
#pragma clang diagnostic pop

        for ( k = 1; k <= nwall[l-1]; k++ )
        {
            ymin = r8_min ( ymin, complex_imag ( zcr[k-1+(l-1)*nw] ) );
            ymax = r8_max ( ymax, complex_imag ( zcr[k-1+(l-1)*nw] ) );
            xmax[l-1] = r8_max ( xmax[l-1], complex_real ( zcr[k-1+(l-1)*nw] ) );
            kp = 1 + ( k % nwall[l-1] );

            if ( 0.0 < complex_real (
                    cmulc( csub(zcr[kp-1+(l-1)*nw] , zcr[k-1+(l-1)*nw]) , complex_conj( csub(wall[kp-1+(l-1)*nw] , wall[k-1+(l-1)*nw]) ))
            ) )
            {
#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnusedValue"
                ls = l;
#pragma clang diagnostic pop
#pragma clang diagnostic push
#pragma ide diagnostic ignored "UnusedValue"
                ks = k;
#pragma clang diagnostic pop
            }
        }
    }
/*
  The "main period" will be between ylimit and ylimit + period.
*/
    ylimit = ( ymin - period + ymax ) / 2.0;
/*
  Project creation points into main period.  This is technical.
*/
    for ( l = 1; l <= nb; l++ )
    {
        for ( k = 1; k <= nwall[l-1]; k++ )
        {
            proj[k-1+(l-1)*nw] = csub(zcr[k-1+(l-1)*nw] , ( complexify(0.0, period *
                                                     ( ( int ) ( 5.0 + ( complex_imag ( zcr[k-1+(l-1)*nw] ) - ylimit )
                                                                       / period ) - 5.0 ))));
        }
    }
/*
  Compute matrix.
*/
    sig2 = pow ( 2.0 * pidp * sigma, 2 );
    i0 = 0;

    for ( l1 = 1; l1 <= nb; l1++ )
    {
        j0 = 0;
        for ( l2 = 1; l2 <= nb; l2++ )
        {
            if ( l1 == l2 )
            {
                kron = 1;
            }
            else
            {
                kron = 0;
            }
            for ( j = 1; j <= nwall[l2-1]; j++ )
            {
                rmatrx[i0+(j0+j-1)*nw*nb] = kron;
                z1 = complex_exp ( ( cmul( csub( wall[0+(l1-1)*nw] , zcr[j-1+(l2-1)*nw] ) , pidp)) );
                z1 = csub(z1 , cinverse(z1));
                dum = sig2 + pow ( complex_real ( z1 ), 2 ) + pow ( complex_imag ( z1 ), 2 );
                for ( i = 2; i <= nwall[l1-1]; i++ )
                {
                    zi = complex_exp ( cmul ( csub(wall[i-1+(l1-1)*nw] , zcr[j-1+(l2-1)*nw])  , pidp) );
                    zz = csub(zi , cinverse(zi));
                    rmatrx[i0+i-1+(j0+j-1)*nw*nb] = -0.25 / pi * log ( dum /
                                                                       ( sig2 + pow ( complex_real ( zz ), 2 ) + pow ( complex_imag ( zz ), 2 ) ) );
                }
            }
            j0 = j0 + nwall[l2-1];
        }
        i0 = i0 + nwall[l1-1];
    }
/*
  Gauss elimination.
*/
    for ( i = 1; i <= matdim; i++ )
    {
        rmatrx[i-1+(i-1)*nw*nb] = 1.0 / rmatrx[i-1+(i-1)*nw*nb];
        for ( j = i + 1; j <= matdim; j++ )
        {
            rmatrx[j-1+(i-1)*nw*nb] = rmatrx[j-1+(i-1)*nw*nb] * rmatrx[i-1+(i-1)*nw*nb];
            for ( k = i + 1; k <= matdim; k++ )
            {
                rmatrx[j-1+(k-1)*nw*nb] = rmatrx[j-1+(k-1)*nw*nb]
                                          - rmatrx[j-1+(i-1)*nw*nb] * rmatrx[i-1+(k-1)*nw*nb];
            }
        }
    }
    return;
}


#endif //SCARLET_STUDY_NKB_GMTRY_H
