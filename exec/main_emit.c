#include <stdio.h>
#include <stdlib.h>

#include "../include/nkb_emit.h"

int main(int argc, char *argv[])
{
    srand(RANDOM_SEED);

    int NB = 5; //[0;9]
    int NW = 100; //100 only?

    if (argc == 3) {
        NB = (int) strtol(argv[1], NULL, 10);
        NW = (int) strtol(argv[2], NULL, 10);
        printf("using NB=%d NW=%d\n", NB, NW);
    }


    const int NV = NB * NW * 2;
    const int NVM = NV + (NB * NW);

    double *cp = malloc(sizeof(double) * NW * NB);
    double *dpds = malloc(sizeof(double) *NW * NB);
    complex_number_t *expmz = malloc(sizeof(complex_number_t) * NVM);
    complex_number_t *expz = malloc(sizeof(complex_number_t) * NVM);
    complex_number_t *force = malloc(sizeof(complex_number_t) *NB);
    double *gamma = malloc(sizeof(double) * NVM);
    int nb = NB;
    int nvm = NVM;
    int nw = NW;
    int *nwall = malloc(sizeof(int) *NB);
    double *ps = malloc(sizeof(double) * NVM);
    double *psi = malloc(sizeof(double)*NW);
    complex_number_t *refpt = malloc(sizeof(complex_number_t) * NB);
    double *rhs = malloc(sizeof(double) * NW * NB);
    double *rmatrx = malloc(sizeof(double) * NW * NB * NW * NB);
    double *rmom = malloc(sizeof(double) * NB);
    complex_number_t *wall = malloc(sizeof(complex_number_t) * NW * NB);
    complex_number_t *z = malloc(sizeof(complex_number_t) * NVM);
    complex_number_t *zcr = malloc(sizeof(complex_number_t) * NW * NB);

    /*
Random initialization.
*/

    for (int j = 0; j < nb; j++) {
        nwall[j] = nw;
        refpt[j] = complexify(0.0, 0.0);
        force[j] = complexify(0.0, 0.0);
        rmom[j] = 0.0;
        for (int i = 0; i < nw; i++) {
            wall[i + j * nw] = complexify(fRand(0.0, 0.999) , (fRand(0.0, 0.999)));
            zcr[i + j * nw] = complexify(fRand(0.0, 0.999) , (fRand(0.0, 0.999)));
            dpds[i + j * nw] = 0.0;
        }
    }

    for (int j = 0; j < nw * nb; j++) {
        rmatrx[j + j * nw * nb] = 1.0;
        rhs[j] = 0.0;
        for (int i = 0; i < j; i++) {
            const double t2 = fRand(0.0, 0.999);
            rmatrx[i + j * nw * nb] = 0.001 * t2;
            rmatrx[j + i * nw * nb] = 0.001 * t2;
        }
    }

    for (int i = 0; i < nvm; i++) {
        const double t1 = fRand(0.0, 0.999);
        const double t2 = fRand(0.0, 0.999);
        z[i] = complexify(t1, t2);
        gamma[i] = fRand(0.0, 0.999);
    }


        emit(nb, nw, cp, dpds, expmz, expz, force, gamma, nwall, ps, psi, refpt, rhs, rmatrx, rmom, wall,
             z, zcr);


    free(cp);
    free(dpds);
    free(expmz);
    free(expz);
    free(force);
    free(gamma);

    free(refpt);
    free(wall);
    free(z);
    free(zcr);

    free(rhs);
    free(rmatrx);
    free(rmom);
    free(psi);
    free(ps);
    free(nwall);

    return 0;
}
