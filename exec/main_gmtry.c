#include <stdio.h>
#include <stdlib.h>
#include "../include/nkb_gmtry.h"


int main(int argc, char *argv[])
{
    srand(RANDOM_SEED);

    int NB = 5; //[0;10]
    int NW = 100; //100 only?

    if (argc == 3) {
        NB = (int) strtol(argv[1], NULL, 10);
        NW = (int) strtol(argv[2], NULL, 10);
        printf("using NB=%d NW=%d\n", NB, NW);
    }



    int nb = NB;
    int nw = NW;
    int nwall[NB];
    complex_number_t *proj = malloc((NW * NB) * sizeof(complex_number_t));
    double *rmatrx = malloc(sizeof(double) * NW * NB * NW * NB);
    complex_number_t *wall = malloc(sizeof(complex_number_t) * NW * NB);
    double *xmax = malloc(sizeof(double) * NB);
    complex_number_t *zcr = malloc(sizeof(complex_number_t) * NW * NB);

    /*
Random initialization.
*/

    for (int j = 1; j <= nb; j++) {
        nwall[j - 1] = nw;
    }

    for (int j = 1; j <= nb; j++) {
        for (int i = 1; i <= nw; i++) {
            const double t1 = fRand(0.0, +0.9999);
            const double t2 = fRand(0.0, +0.9999);
            wall[i - 1 + (j - 1) * NW] = complexify(t1 , t2 );
        }
    }


        gmtry(nb, nw, nwall, proj, rmatrx, wall, xmax, zcr);


    free(proj);
    free(rmatrx);
    free(wall);
    free(xmax);
    free(zcr);



    return 0;
}
