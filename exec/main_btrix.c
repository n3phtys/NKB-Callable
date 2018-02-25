#include <stdio.h>
#include <stdlib.h>

#include "../include/nkb_btrix.h"

int main(int argc, char *argv[]) {
    srand(RANDOM_SEED);

    int JD = 50; //[4;64]
    int KD = 50; //[4;64]
    int LE = 28;

    if (argc == 3) {
        JD = (int) strtol(argv[1], NULL, 10);
        KD = (int) strtol(argv[2], NULL, 10);
        LE = (int) strtol(argv[3], NULL, 10);
        printf("using JD=%d KD=%d LE=%d\n", JD, KD, LE);
    }


//4 dimensions
    const int LD = 50;
    const int MD = 50;

    double *a = malloc(sizeof(double) * 5 * 5 * MD * MD);
    double *b = malloc(sizeof(double) * 5 * 5 * MD * MD);
    double *bx = malloc(sizeof(double) * 5 * 5 * MD * MD);
    double *c = malloc(sizeof(double) * 5 * 5 * MD * MD);
    double *s = malloc(sizeof(double) * JD * KD * LD * 5);
    double *sx = malloc(sizeof(double) * JD * KD * LD * 5);
    const int jd = JD;
    const int kd = KD;
    const int ld = LD;
    const int md = MD;

    const int js = 1;
    const int je = 28;
    const int ls = 1;
    const int le = LE;
    const int nb = 25 * md * md;
    const int ns = jd * kd * ld * 5;
/*
  Random initialization.
*/


    for (int l = 0; l < md; l++) {
        for (int k = 0; k < md; k++) {
            for (int j = 0; j < 5; j++) {
                for (int i = 0; i < 5; i++) {
                    a[i + j * 5 + k * 5 * 5 + l * 5 * 5 * md] = fRand(0.0, 0.999);
                    bx[i + j * 5 + k * 5 * 5 + l * 5 * 5 * md] = fRand(0.0, 0.999);
                    c[i + j * 5 + k * 5 * 5 + l * 5 * 5 * md] = fRand(0.0, 0.999);
                }
            }
        }
    }

    for (int l = 0; l < 5; l++) {
        for (int k = 0; k < ld; k++) {
            for (int j = 0; j < kd; j++) {
                for (int i = 0; i < jd; i++) {
                    sx[i + j * jd + k * jd * kd + l * jd * kd * ld] = fRand(0.0, 0.999);
                }
            }
        }
    }


    btrixWrapped(ns, sx, nb, bx, js, je, ls, le, jd, kd, ld, md, a, b, c, s);

    free(a);
    free(b);
    free(bx);
    free(c);
    free(s);
    free(sx);


    return 0;
}
