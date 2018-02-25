#include <stdio.h>
#include <stdlib.h>

#include "../include/nkb_vpenta.h"

int main(int argc, char *argv[]) {
    srand(RANDOM_SEED);
    int U = 64; //[4;1024]

    if (argc == 2) {
        U = (int) strtol(argv[1], NULL, 10);
        printf("using U=%d\n", U);
    }


    const int JL = 0;
    const int JU = U;
    const int KL = 0;
    const int KU = U;
    const int NJA = U + 1;
    const int NJB = U + 1;

    double *a = malloc(sizeof(double) * NJA * NJB);
    double *b = malloc(sizeof(double) * NJA * NJB);
    double *c = malloc(sizeof(double) * NJA * NJB);
    double *d = malloc(sizeof(double) * NJA * NJB);
    double *e = malloc(sizeof(double) * NJA * NJB);
    double *f = malloc(sizeof(double) * NJA * NJB * 3);
    double *fx = malloc(sizeof(double) * NJA * NJB * 3);
    const int jl = JL;
    const int ju = JU;
    const int kl = KL;
    const int ku = KU;
    const int nja = NJA;
    const int njb = NJB;
    const int lf = nja * njb * 3;
    double *x = malloc(sizeof(double) * NJA * NJB);
    double *y = malloc(sizeof(double) * NJA * NJB);


    /*
Random initialization.
*/

    for (int j = kl; j <= ku; j++) {
        for (int i = jl; i <= ju; i++) {
            a[i + j * nja] = fRand(0.0, 0.999);
            b[i + j * nja] = fRand(0.0, 0.999);
            c[i + j * nja] = fRand(0.0, 0.999);
            d[i + j * nja] = fRand(0.0, 0.999);
            e[i + j * nja] = fRand(0.0, 0.999);
            for (int k = 0; k < 3; k++) {
                fx[i + j * nja + k * nja * njb] = fRand(0.0, 0.999);
            }
        }
    }


    r8vec_copy(lf, fx, f);
    vpenta(U, a, b, c, d, e, f, x, y);


    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(f);
    free(fx);
    free(x);
    free(y);


    return 0;
}
