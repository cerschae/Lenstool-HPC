#include <math.h>
#include <stdlib.h>
#include <stdio.h>
/* Cholesky decomposition. Taken from nR 2.9
   Inputs:
       n, integer
       a, n x n matrix
   Returns:
       The Cholesky factor L is returned in the
       lower triangle of a, except for its diagonal
       elements which are returned in p.

       p must be allocated before entering the function.
 */
void cholesky(double **a, int n, double *p)
{
    int i, j, k;
    double sum;

    for (i = 0; i < n; i++)
        for (j = i; j < n; j++)
        {
            sum = a[i][j];
            for ( k = i - 1; k >= 0; k-- )
                sum -= a[i][k] * a[j][k];

            if ( i == j )
            {
                if ( sum <= 0. )
                {
                    fprintf(stderr, "ERROR in cholesky.c. a, with rounding errors, is not positive definie\n");
                    exit(1);
                }
                p[i] = sqrt(sum);
            }
            else
                a[j][i] = sum / p[i];

        }
}

/* End of cholesky.c */

