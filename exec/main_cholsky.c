#include <stdio.h>
#include <stdlib.h>

#include "../include/nkb_cholsky.h"

int main(int argc, char *argv[])
{
    srand(RANDOM_SEED);


    int N = 40; //[1;128]
    int NRHS = 3; //[1;13]

    if (argc == 3) {
        N = (int) strtol(argv[1], NULL, 10);
        NRHS = (int) strtol(argv[2], NULL, 10);
        printf("using N=%d NRHS=%d\n", N, NRHS);
    }


    const int IDA = 250;
    const int M = 4;
    const int n = N;
    const int NMAT = 250;

    double a[(IDA + 1) * (M + 1) * (N + 1)];
    double ax[(IDA + 1) * (M + 1) * (N + 1)];
    double b[(NRHS + 1) * (NMAT + 1) * (N + 1)];
    double bx[(NRHS + 1) * (NMAT + 1) * (N + 1)];
    int i1;
    const int ida = IDA;
    const int m = M;
    const int nmat = NMAT;
    const int nrhs = NRHS;
    const int la = (ida + 1) * (m + 1) * (n + 1);
    const int lb = (nrhs + 1) * (nmat + 1) * (n + 1);
/*
  Random initialization.
*/

    i1 = 0;

    for (int k = 0; k <= n; k++) {
        for (int j = -m; j <= 0; j++) {
            for (int i = 0; i <= ida; i++) {
                ax[i1] = fRand(0.0, 0.999);
                i1 = i1 + 1;
            }
        }
    }

    i1 = 0;
    for (int k = 0; k <= n; k++) {
        for (int j = 0; j <= nmat; j++) {
            for (int i = 0; i <= nrhs; i++) {
                bx[i1] = fRand(0.0, 0.999);
                i1 = i1 + 1;
            }
        }
    }

        r8vec_copy(la, ax, a);
        r8vec_copy(lb, bx, b);
        cholsky(ida, nmat, m, n, a, nrhs, ida, b);




    return 0;
}
