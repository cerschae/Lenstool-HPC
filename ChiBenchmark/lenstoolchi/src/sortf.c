#include <stdlib.h>

//#define null  0


struct arbre
{
    double N;
    struct arbre *FG;
    struct arbre *FD;
};

static struct   arbre   *addf(double M, struct arbre *A, int (*comp)(double, double));
static void readf(struct arbre *R, double *B, int *i);


/*
*  sort double array
*/

void    sortf(int n, double *A, int (*comp)(double, double))
{
    struct  arbre   *racine;
    int i;

    racine = NULL;
    for (i = 0; i < n; i++)
        racine = addf(A[i], racine, comp);
    i = 0;
    readf(racine, A, &i);

}

/***************************************************************/

static struct   arbre   *addf(double M, struct arbre *A, int (*comp)(double, double))
{
    if (A == NULL)
    {
        A = (struct arbre *)malloc(sizeof(struct arbre));
        A->N = M;
        A->FG = A->FD = NULL;
    }
    else
    {
        if ((*comp)(M, A->N))
            (A->FG) = addf(M, A->FG, comp);
        else
            (A->FD) = addf(M, A->FD, comp);
    };
    return(A);

}

/***************************************************************/

static void readf(struct arbre *R, double *B, int *i)
{
    if (R != NULL)
    {
        readf(R->FG, B, i);
        B[(*i)++] = R->N;
        readf(R->FD, B, i);
    }
}

/***************************************************************/

int comp_asc(double A, double B)
{
    return(A <= B);
}

/***************************************************************/

int comp_desc(double A, double B)
{
    return(A <= B);
}
