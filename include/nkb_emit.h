//
// Created by nephtys on 21.01.18.
//

#ifndef SCARLET_STUDY_NKB_EMIT_H
#define SCARLET_STUDY_NKB_EMIT_H

#include "common.h"
#include "custom_complex.h"




void emit_test ( double *er, double *fp, double *tm );
void emit ( int nb, int nw, double cp[], double dpds[],
            complex_number_t expmz[], complex_number_t expz[], complex_number_t force[],
            double gamma[], int nwall[], double ps[], double psi[],
            complex_number_t refpt[], double rhs[], double rmatrx[], double rmom[],
            complex_number_t wall[], complex_number_t z[], complex_number_t zcr[] );



/******************************************************************************/

void emit_test ( double *er, double *fp, double *tm )

/******************************************************************************/
/*
  Purpose:

    EMIT_TEST tests EMIT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{
# define NB 5
# define NV 1000
# define NVM 1500
# define NW 100

    double ans;
    double cp[NW*NB];
    double dpds[NW*NB];
    complex_number_t expmz[NVM];
    complex_number_t expz[NVM];
    double f7;
    complex_number_t force[NB];
    double gamma[NVM];
    int i;
    int it;
    int j;
    int nb = NB;
    int nv = NV;
    int nvm = NVM;
    int nw = NW;
    int nwall[NB];
    double ps[NVM];
    double psi[NW];
    complex_number_t refpt[NB];
    double rhs[NW*NB];
    double rmatrx[NW*NB*NW*NB];
    double rmom[NB];
    double t1;
    double t2;
    double t30;
    double time1;
    complex_number_t wall[NW*NB];
    complex_number_t z[NVM];
    complex_number_t zcr[NW*NB];

    it = 10;
    ans = 6.0088546832072;
/*
  Random initialization.
*/
    f7 = 78125.0;
    t30 = 1073741824.0;
    t2 = f7 / t30;

    for ( j = 0; j < nb; j++ )
    {
        nwall[j] = nw;
        refpt[j] = complexify(0.0, 0.0);
        force[j] = complexify(0.0, 0.0);
        rmom[j] = 0.0;
        for ( i = 0; i < nw; i++ )
        {
            t1 = fmod ( f7 * t2, 1.0 );
            t2 = fmod ( f7 * t1, 1.0 );
            wall[i+j*nw] = complexify(t1, t2);
            t1 = fmod ( f7 * t2, 1.0 );
            t2 = fmod ( f7 * t1, 1.0 );
            zcr[i+j*nw] = complexify(t1, t2);
            dpds[i+j*nw] = 0.0;
        }
    }

    for ( j = 0; j < nw * nb; j++ )
    {
        rmatrx[j+j*nw*nb] = 1.0;
        for ( i = 0; i < j; i++ )
        {
            t2 = fmod ( f7 * t2, 1.0 );
            rmatrx[i+j*nw*nb] = 0.001 * t2;
            rmatrx[j+i*nw*nb] = 0.001 * t2;
        }
    }

    for ( i = 0; i < nvm; i++ )
    {
        t1 = fmod ( f7 * t2, 1.0 );
        t2 = fmod ( f7 * t1, 1.0 );
        z[i] = complexify(t1, t2);
        t2 = fmod ( f7 * t2, 1.0 );
        gamma[i] = t2;
    }
/*
  Timing.
*/
    time1 = wtime ( );

    for ( i = 1; i <= it; i++ )
    {
        emit ( nb, nw, cp, dpds, expmz, expz, force,
               gamma, nwall, ps, psi, refpt, rhs, rmatrx, rmom, wall, z, zcr );
    }

    *tm = wtime ( ) - time1;
/*
  Results.
*/
    *er = r8_abs ( ( rhs[18] - ans ) / ans );
    *fp = ( double ) ( it * ( 56 * nv + nb * nw *
                                        ( 97 + 44 * nv + 2 * nb * nw ) ) );

    return;
# undef NB
# undef NV
# undef NVM
# undef NW
}
/******************************************************************************/

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
void  emit ( int nb, int nw, double cp[], double dpds[],
#pragma clang diagnostic pop
             complex_number_t expmz[], complex_number_t expz[], complex_number_t force[],
             double gamma[], int nwall[], double ps[], double psi[],
             complex_number_t refpt[], double rhs[], double rmatrx[], double rmom[],
             complex_number_t wall[], complex_number_t z[], complex_number_t zcr[] )

/******************************************************************************/
/*
  Purpose:

    EMIT creates new vortices according to certain boundary conditions.

  Discussion:

    This function was extracted from a vortex method program.
    It emits new vortices to satisfy the boundary condition.
    It also finishes computing pressure, forces, and other quantities.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 November 2010

  Author:

    Original FORTRAN77 version by David Bailey.
    C version by John Burkardt.
*/
{

    int nv = nb * nw * 2;
    int nvm = nv + (nb * nw);

    double chord;
    double cpm;
    double cupst;
    double delt;
    complex_number_t dum3;
    complex_number_t expmwk;
    complex_number_t expwkl;
    int i;
    int i0;
    int j;
    int k;
    int k1;
    int k2;
    int l;
    int matdim;
    double period;
    double pi;
    double pidp;
    double sig2;
    double sps;
    double u0;
    complex_number_t uupstr;

    period = 3.0;
    sig2 = 3.0;
    u0 = 4.0;
    matdim = nw * nb;
    delt = 1.0;
    chord = 5.0;
    pi = 3.141592653589793;
    uupstr = complexify(3.0 , 4.0);
/*
  Store exp(z(i)) and exp(-z(i)) to reduce work in inner loop.

  Note that the NV used here is a variable, whereas the NV in the
  calling program is a constant.  They are separate quantities.
*/
    pidp = pi / period;

    for ( i = 0; i < nv; i++ )
    {
        expz[i] = complex_exp ( cmul(z[i] , pidp) );
        expmz[i] = cinverse(expz[i]);
    }

    i0 = 0;
    cupst = pow ( complex_real ( uupstr ), 2 ) + pow ( complex_imag ( uupstr ), 2 );

    for ( l = 0; l < nb; l++ )
    {
        for ( k = 0; k < nwall[l]; k++ )
        {
            expwkl = complex_exp ( cmul(wall[k+l*nw], pidp) );
            expmwk = cinverse(expwkl);

            sps = 0.0;
            for ( i = 0; i < nv; i++ )
            {
                dum3 = csub(cmulc(expz[i] , expmwk) , cmulc(expwkl , expmz[i]));
                ps[i] = gamma[i] * log ( pow ( complex_real ( dum3 ), 2 ) +
                                         pow ( complex_imag ( dum3 ), 2 ) + sig2 );
                sps = sps + ps[i];
            }
            psi[k] = complex_imag ( cmulc(wall[k+l*nw] , complex_conj ( cadd(uupstr , complexify(0.0, u0))) ))
                     - sps * 0.25 / pi;
        }
/*
  Compute the right-hand side.
*/
        for ( k = 0; k < nwall[l]; k++ )
        {
            rhs[i0+k] = psi[k] - psi[0];
        }
        i0 = i0 + nwall[l];
    }
/*
  Solve the system.
*/
    for ( i = 0; i < matdim; i++ )
    {
        for ( j = i + 1; j < matdim; j++ )
        {
            rhs[j] = rhs[j] - rmatrx[j+i*nw*nb] * rhs[i];
        }
    }

    for ( i = matdim - 1; 0 <= i; i-- )
    {
        rhs[i] = rmatrx[i+i*nw*nb] * rhs[i];
        for ( j = 0; j < i; j++ )
        {
            rhs[j] = rhs[j] - rmatrx[j+i*nw*nb] * rhs[i];
        }
    }
/*
  Create new vortices.
*/
    i0 = 0;

    for ( l = 0; l < nb; l++ )
    {
        for ( k = 0; k < nwall[l]; k++ )
        {
/*
  Put the new vortex at the end of the array.
*/
            z[nv] = zcr[k+l*nw];
            gamma[nv] = rhs[i0+k];
/*
  Record the gain of linear and angular momentum.
*/
            force[l] = cadd(force[l] , cmul(z[nv], gamma[nv]));
            rmom[l] = rmom[l] + gamma[nv] * (
                    pow ( complex_real ( csub(z[nv] , refpt[l]) ), 2 ) +
                    pow ( complex_imag ( csub(z[nv] , refpt[l]) ), 2 ) );
            dpds[k+l*nw] = dpds[k+l*nw] - gamma[nv];
            nv = nv + 1;
        }
/*
  Filter and integrate pressure gradient to get pressure.
*/
        cp[0+l*nw] = 0.0;
        cpm = -1.0E+30;

        for ( k = 1; k < nwall[l]; k++ )
        {
            k1 = k % nwall[l];
            k2 = ( k + nwall[l] - 3 ) % nwall[l];
            cp[k+l*nw] = cp[k-1+l*nw] + ( 3.0 * ( dpds[k+l*nw] + dpds[k-1+l*nw] )
                                          + dpds[k1+l*nw] + dpds[k2+l*nw] ) / ( 4.0 * delt * cupst );
            cpm = r8_max ( cpm, cp[k+l*nw] );
        }
/*
  Normalize the pressure.
*/
        for ( k = 0; k < nwall[l]; k++ )
        {
            cp[k+l*nw] = cp[k+l*nw] - cpm;
        }
/*
  Finish computing force and moment, as time rate of change of linear
  and angular momentum.
*/
        force[l] = cdiv(cmulc(force[l] , complexify(0.0, 2.0)) , ( delt * chord * cupst ));
        rmom[l] = rmom[l] * 2.0 / ( delt * chord * chord * cupst );
        i0 = i0 + nwall[l];
    }
    return;
}


#endif //SCARLET_STUDY_NKB_EMIT_H
