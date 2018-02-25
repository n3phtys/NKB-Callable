#include <stdio.h>
#include <stdlib.h>

#include "../include/nkb_mxm.h"





int main(int argc, char *argv[])
{
    srand(RANDOM_SEED);

    int A = 256;
    int B = 128;
    int C = 64;

    if (argc == 4) {
        A = (int) strtol(argv[1], NULL, 10);
        B = (int) strtol(argv[2], NULL, 10);
        C = (int) strtol(argv[3], NULL, 10);
        printf("using A=%d B=%d C=%d\n", A, B, C);
    }


    double *a = malloc((A * B) * sizeof(double));
    double *b = malloc((B * C) * sizeof(double));
    double *c = malloc((A * C) * sizeof(double));

    fillArrayWithRandomValues(a, A * B);
    fillArrayWithRandomValues(b, B * C);

        mxm(a, b, c, A, B, C);


    free(a);
    free(b);
    free(c);


    return 0;
}
