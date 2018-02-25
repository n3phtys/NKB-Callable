//
// Created by nephtys on 21.01.18.
//

#ifndef SCARLET_STUDY_NKB_BTRIX_H
#define SCARLET_STUDY_NKB_BTRIX_H

#include "./common.h"


void btrix_test ( double *er, double *fp, double *tm );
void btrix ( int js, int je, int ls, int le, int k, int jd, int kd, int ld,
             int md, double a[], double b[], double c[], double s[] );


/******************************************************************************/

void btrix_test ( double *er, double *fp, double *tm )

/******************************************************************************/
/*
  Purpose:

    BTRIX_TEST tests BTRIX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
# define JD 30
# define KD 30
# define LD 30
# define MD 30

    double* a = (double*) malloc(sizeof(double) * 5*5*MD*MD);
    double ans;
    double* b = (double*) malloc(sizeof(double) * 5*5*MD*MD);
    double* bx = (double*) malloc(sizeof(double) * 5*5*MD*MD);
    double* c = (double*) malloc(sizeof(double) * 5*5*MD*MD);
    double f7;
    int i;
    int ii;
    int it;
    int j;
    int jd = JD;
    int je;
    int js;
    int k;
    int kd = KD;
    int l;
    int ld = LD;
    int le;
    int ls;
    int md = MD;
    int nb;
    int ns;
    double* s = (double*) malloc(sizeof(double) * JD*KD*LD*5);
    double* sx = (double*) malloc(sizeof(double) * JD*KD*LD*5);
    double t;
    double t30;
    double time1;

    js = 1;
    je = 28;
    ls = 1;
    le = 28;
    it = 20;
    ans = -0.286282658663962;
    nb = 25 * md * md;
    ns = jd * kd * ld * 5;
/*
  Random initialization.
*/
    f7 = 78125.0;
    t30 = 1073741824.0;
    t = f7 / t30;

    for ( l = 0; l < md; l++ )
    {
        for ( k = 0; k < md; k++ )
        {
            for ( j = 0; j < 5; j++ )
            {
                for ( i = 0; i < 5; i++ )
                {
                    t = fmod ( f7 * t, 1.0 );
                    a[i+j*5+k*5*5+l*5*5*md] = t;
                    t = fmod ( f7 * t, 1.0 );
                    bx[i+j*5+k*5*5+l*5*5*md] = t;
                    t = fmod ( f7 * t, 1.0 );
                    c[i+j*5+k*5*5+l*5*5*md] = t;
                }
            }
        }
    }

    for ( l = 0; l < 5; l++ )
    {
        for ( k = 0; k < ld; k++ )
        {
            for ( j = 0; j < kd; j++ )
            {
                for ( i = 0; i < jd; i++ )
                {
                    t = fmod ( f7 * t, 1.0 );
                    sx[i+j*jd+k*jd*kd+l*jd*kd*ld] = t;
                }
            }
        }
    }
/*
  Timing.
*/
    time1 = wtime ( );

    for ( ii = 1; ii <= it; ii++ )
    {
        r8vec_copy ( ns, sx, s );
        for ( k = 1; k <= kd; k++ )
        {
            r8vec_copy ( nb, bx, b );
            btrix ( js, je, ls, le, k, jd, kd, ld, md, a, b, c, s );
        }
    }

    *tm = wtime ( ) - time1;
/*
  Results.
*/
    *er = r8_abs ( ( s[18+18*jd+18*jd*kd+0*jd*kd*ld] - ans ) / ans );
    *fp = ( double ) ( it * md * ( le - 1 ) * 19165 );


    free(a);
    free(b);
    free(bx);
    free(c);
    free(s);
    free(sx);


    return;
# undef JD
# undef KD
# undef LD
# undef MD
}
/******************************************************************************/

void  btrix ( int js, int je, int ls, int le, int k, int jd, int kd, int ld,
              const int MD, double a[], double b[], double c[], double s[] )

/******************************************************************************/
/*
  Purpose:

    BTRIX is a block tridiagonal solver in one direction.

  Discussion:

    The array has four dimensions.  The routine solves along the
    "J" index.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
    const int md = MD;
    double c1;
    double c2;
    double c3;
    double c4;
    double c5;
    double d1;
    double d2;
    double d3;
    double d4;
    double d5;
    int j;
    int jem1;
    int l;
    double* l11 =(double*) malloc(sizeof(double) * MD);
    double* l21 =(double*) malloc(sizeof(double) * MD);
    double* l31 =(double*) malloc(sizeof(double) * MD);
    double* l41 =(double*) malloc(sizeof(double) * MD);
    double* l51 =(double*) malloc(sizeof(double) * MD);
    double* l22 =(double*) malloc(sizeof(double) * MD);
    double* l32 =(double*) malloc(sizeof(double) * MD);
    double* l33 =(double*) malloc(sizeof(double) * MD);
    double* l42 =(double*) malloc(sizeof(double) * MD);
    double* l43 =(double*) malloc(sizeof(double) * MD);
    double* l44 =(double*) malloc(sizeof(double) * MD);
    double* l52 =(double*) malloc(sizeof(double) * MD);
    double* l53 =(double*) malloc(sizeof(double) * MD);
    double* l54 =(double*) malloc(sizeof(double) * MD);
    double* l55 =(double*) malloc(sizeof(double) * MD);
    int m;
    int n;
    double* u12 =(double*) malloc(sizeof(double) * MD);
    double* u13 =(double*) malloc(sizeof(double) * MD);
    double* u14 =(double*) malloc(sizeof(double) * MD);
    double* u15 =(double*) malloc(sizeof(double) * MD);
    double* u23 =(double*) malloc(sizeof(double) * MD);
    double* u24 =(double*) malloc(sizeof(double) * MD);
    double* u25 =(double*) malloc(sizeof(double) * MD);
    double* u34 =(double*) malloc(sizeof(double) * MD);
    double* u35 =(double*) malloc(sizeof(double) * MD);
    double* u45 =(double*) malloc(sizeof(double) * MD);
/*
  Part 1.  Forward block sweep.
*/
    for ( j = js; j <= je; j++ )
    {
/*
  Step 1.  Construct L(I) in B.
*/
        if ( j != js )
        {
            for ( m = 0; m < 5; m++ )
            {
                for ( n = 0; n < 5; n++ )
                {
                    for ( l = ls; l <= le; l++ )
                    {
                        b[m+n*5+j*5*5+l*5*5*md] = b[m+n*5+j*5*5+l*5*5*md]
                                                  - a[m+0*5+j*5*5+l*5*5*md] * b[0+n*5+(j-1)*5*5+l*5*5*md]
                                                  - a[m+1*5+j*5*5+l*5*5*md] * b[1+n*5+(j-1)*5*5+l*5*5*md]
                                                  - a[m+2*5+j*5*5+l*5*5*md] * b[2+n*5+(j-1)*5*5+l*5*5*md]
                                                  - a[m+3*5+j*5*5+l*5*5*md] * b[3+n*5+(j-1)*5*5+l*5*5*md]
                                                  - a[m+4*5+j*5*5+l*5*5*md] * b[4+n*5+(j-1)*5*5+l*5*5*md];
                    }
                }
            }
        }
/*
  Step 2.  Compute L inverse.

  A.  Decompose L(I) into L and U.
*/
        for ( l = ls; l <= le; l++ )
        {
            l11[l] = 1.0 / b[0+0*5+j*5*5+l*5*5*md];
            u12[l] = b[0+1*5+j*5*5+l*5*5*md] * l11[l];
            u13[l] = b[0+2*5+j*5*5+l*5*5*md] * l11[l];
            u14[l] = b[0+3*5+j*5*5+l*5*5*md] * l11[l];
            u15[l] = b[0+4*5+j*5*5+l*5*5*md] * l11[l];
            l21[l] = b[1+0*5+j*5*5+l*5*5*md];
            l22[l] = 1.0 / ( b[1+1*5+j*5*5+l*5*5*md] - l21[l] * u12[l] );
            u23[l] = ( b[1+2*5+j*5*5+l*5*5*md] - l21[l] * u13[l] ) * l22[l];
            u24[l] = ( b[1+3*5+j*5*5+l*5*5*md] - l21[l] * u14[l] ) * l22[l];
            u25[l] = ( b[1+4*5+j*5*5+l*5*5*md] - l21[l] * u15[l] ) * l22[l];
            l31[l] = b[2+0*5+j*5*5+l*5*5*md];
            l32[l] = b[2+1*5+j*5*5+l*5*5*md] - l31[l] * u12[l];
            l33[l] = 1.0 / ( b[2+2*5+j*5*5+l*5*5*md] - l31[l] * u13[l] - l32[l]
                                                                         * u23[l] );
            u34[l]  = ( b[2+3*5+j*5*5+l*5*5*md] - l31[l] * u14[l] - l32[l] * u24[l] )
                      * l33[l];
            u35[l]  = ( b[2+4*5+j*5*5+l*5*5*md] - l31[l] * u15[l] - l32[l] * u25[l] )
                      * l33[l];
        }

        for ( l = ls; l <= le; l++ )
        {
            l41[l] = b[3+0*5+j*5*5+l*5*5*md];
            l42[l] = b[3+1*5+j*5*5+l*5*5*md] - l41[l] * u12[l];
            l43[l] = b[3+2*5+j*5*5+l*5*5*md] - l41[l] * u13[l] - l42[l] * u23[l];
            l44[l] = 1.0 / ( b[3+3*5+j*5*5+l*5*5*md] - l41[l] * u14[l]
                             - l42[l] * u24[l] - l43[l] * u34[l] );
            u45[l] = ( b[3+4*5+j*5*5+l*5*5*md] - l41[l] * u15[l] - l42[l] * u25[l]
                       - l43[l] * u35[l] ) * l44[l];
            l51[l] = b[4+0*5+j*5*5+l*5*5*md];
            l52[l] = b[4+1*5+j*5*5+l*5*5*md] - l51[l] * u12[l];
            l53[l] = b[4+2*5+j*5*5+l*5*5*md] - l51[l] * u13[l] - l52[l] * u23[l];
            l54[l] = b[4+3*5+j*5*5+l*5*5*md] - l51[l] * u14[l] - l52[l] * u24[l]
                     - l53[l] * u34[l];
            l55[l] = 1.0 / ( b[4+4*5+j*5*5+l*5*5*md] - l51[l] * u15[l] - l52[l] * u25[l]
                             - l53[l] * u35[l] - l54[l] * u45[l] );
        }
/*
  Step 3.  Solve for intermediate vector.

  A.  Construct the right hand side.
*/
        if ( j != js )
        {
            for ( m = 0; m < 5; m++ )
            {
                for ( l = ls; l <= le; l++ )
                {
                    s[j+k*jd+l*jd*kd+m*jd*kd*ld] = s[j+k*jd+l*jd*kd+m*jd*kd*ld]
                                                   - a[m+0*5+j*5*5+l*5*5*md] * s[j-1+k*jd+l*jd*kd+0*jd*kd*ld]
                                                   - a[m+1*5+j*5*5+l*5*5*md] * s[j-1+k*jd+l*jd*kd+1*jd*kd*ld]
                                                   - a[m+2*5+j*5*5+l*5*5*md] * s[j-1+k*jd+l*jd*kd+2*jd*kd*ld]
                                                   - a[m+3*5+j*5*5+l*5*5*md] * s[j-1+k*jd+l*jd*kd+3*jd*kd*ld]
                                                   - a[m+4*5+j*5*5+l*5*5*md] * s[j-1+k*jd+l*jd*kd+4*jd*kd*ld];
                }
            }
        }
/*
  B. Intermediate vector.

  Forward substitution.
*/
        for ( l = ls; l <= le; l++ )
        {
            d1 =   s[j+k*jd+l*jd*kd+0*jd*kd*ld] * l11[l];
            d2 = ( s[j+k*jd+l*jd*kd+1*jd*kd*ld] - l21[l] * d1 ) * l22[l];
            d3 = ( s[j+k*jd+l*jd*kd+2*jd*kd*ld] - l31[l] * d1 - l32[l] * d2 ) * l33[l];
            d4 = ( s[j+k*jd+l*jd*kd+3*jd*kd*ld] - l41[l] * d1 - l42[l] * d2
                   - l43[l] * d3 ) * l44[l];
            d5 = ( s[j+k*jd+l*jd*kd+4*jd*kd*ld] - l51[l] * d1 - l52[l] * d2
                   - l53[l] * d3 - l54[l] * d4 ) * l55[l];
/*
  Backward substitution.
*/
            s[j+k*jd+l*jd*kd+4*jd*kd*ld] = d5;
            s[j+k*jd+l*jd*kd+3*jd*kd*ld] = d4 - u45[l] * d5;
            s[j+k*jd+l*jd*kd+2*jd*kd*ld] = d3 - u34[l] * s[j+k*jd+l*jd*kd+3*jd*kd*ld]
                                           - u35[l] * d5;
            s[j+k*jd+l*jd*kd+1*jd*kd*ld] = d2 - u23[l] * s[j+k*jd+l*jd*kd+2*jd*kd*ld]
                                           - u24[l] * s[j+k*jd+l*jd*kd+3*jd*kd*ld] - u25[l] * d5;
            s[j+k*jd+l*jd*kd+0*jd*kd*ld] = d1 - u12[l] * s[j+k*jd+l*jd*kd+1*jd*kd*ld]
                                           - u13[l] * s[j+k*jd+l*jd*kd+2*jd*kd*ld]
                                           - u14[l] * s[j+k*jd+l*jd*kd+3*jd*kd*ld] - u15[l] * d5;
        }
/*
  Step 4.  Construct U(I) = inverse(L(I))*C(I+1) by columns and store in B.
*/
        if ( j != je )
        {
            for ( n = 0; n < 5; n++ )
            {
                for ( l = ls; l <= le; l++ )
                {
/*
  Forward substitution.
*/
                    c1 =   c[0+n*5+j*5*5+l*5*5*md] * l11[l];
                    c2 = ( c[1+n*5+j*5*5+l*5*5*md] - l21[l] * c1 ) * l22[l];
                    c3 = ( c[2+n*5+j*5*5+l*5*5*md] - l31[l] * c1 - l32[l] * c2 ) * l33[l];
                    c4 = ( c[3+n*5+j*5*5+l*5*5*md] - l41[l] * c1 - l42[l] * c2
                           - l43[l] * c3 ) * l44[l];
                    c5 = ( c[4+n*5+j*5*5+l*5*5*md] - l51[l] * c1 - l52[l] * c2
                           - l53[l] * c3 - l54[l] * c4 ) * l55[l];
/*
  Backward substitution.
*/
                    b[4+n*5+j*5*5+l*5*5*md] = c5;
                    b[3+n*5+j*5*5+l*5*5*md] = c4 - u45[l] * c5;
                    b[2+n*5+j*5*5+l*5*5*md] = c3 - u34[l] * b[3+n*5+j*5*5+l*5*5*md]
                                              - u35[l] * c5;
                    b[1+n*5+j*5*5+l*5*5*md] = c2 - u23[l] * b[2+n*5+j*5*5+l*5*5*md]
                                              - u24[l] * b[3+n*5+j*5*5+l*5*5*md] - u25[l] * c5;
                    b[0+n*5+j*5*5+l*5*5*md] = c1 - u12[l] * b[1+n*5+j*5*5+l*5*5*md]
                                              - u13[l] * b[2+n*5+j*5*5+l*5*5*md]
                                              - u14[l] * b[3+n*5+j*5*5+l*5*5*md] - u15[l] * c5;
                }
            }
        }
    }
/*
  Part 2.  Backward block sweep.
*/
    jem1 = je - 1;

    for ( j = jem1; js <= j; j-- )
    {
        for ( m = 0; m < 5; m++ )
        {
            for ( l = ls; l <= le; l++ )
            {
                s[j+k*jd+l*jd*kd+m*jd*kd*ld] = s[j+k*jd+l*jd*kd+m*jd*kd*ld]
                                               - b[m+0*5+j*5*5+l*5*5*md] * s[j+1+k*jd+l*jd*kd+0*jd*kd*ld]
                                               - b[m+1*5+j*5*5+l*5*5*md] * s[j+1+k*jd+l*jd*kd+1*jd*kd*ld]
                                               - b[m+2*5+j*5*5+l*5*5*md] * s[j+1+k*jd+l*jd*kd+2*jd*kd*ld]
                                               - b[m+3*5+j*5*5+l*5*5*md] * s[j+1+k*jd+l*jd*kd+3*jd*kd*ld]
                                               - b[m+4*5+j*5*5+l*5*5*md] * s[j+1+k*jd+l*jd*kd+4*jd*kd*ld];
            }
        }
    }



    free(l11);
    free(l21);
    free(l31);
    free(l41);
    free(l51);
    free(l22);
    free(l32);
    free(l33);
    free(l42);
    free(l43);
    free(l44);
    free(l52);
    free(l53);
    free(l54);
    free(l55);
    free(u12);
    free(u13);
    free(u14);
    free(u15);
    free(u23);
    free(u24);
    free(u25);
    free(u34);
    free(u35);
    free(u45);

    return;
}


void btrixWrapped(int ns, double *sx, int nb, double *bx, const int js, const int je, const int ls, const int le, const int jd,
                  const int kd, const int ld, const int md, double *a, double *b, double *c, double *s) {
    for (int i = 0; i < 10; i++) {
        r8vec_copy(ns, sx, s);
        for (int k = 1; k <= kd; k++) {
            r8vec_copy(nb, bx, b);
            btrix(js, je, ls, le, k, jd, kd, ld, md, a, b, c, s);
        }
    }
    return;
}


#endif //SCARLET_STUDY_NKB_BTRIX_H
