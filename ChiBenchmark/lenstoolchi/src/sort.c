#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

static struct arbre* mettre(struct galaxie *M, struct arbre *A,
                            int (*comp)(struct galaxie *, struct galaxie *) );
static void lecture(struct arbre *R, struct galaxie *B, long int *i);


void    sort(long int n, struct galaxie *A, int (*comp)(struct galaxie *, struct galaxie *))
{
    struct arbre   *racine;
    struct galaxie *B;
    long int i;

    racine = NULL;
    for (i = 0L; i < n; i++)
        racine = mettre(&A[i], racine, comp);
    i = 0L;
    B = (struct galaxie *)malloc((unsigned long int)n*sizeof(struct galaxie));
    lecture(racine, B, &i);
    
    for (i = 0L; i < n; i++)
        A[i] = B[i];

    free(B);

}

/***************************************************************/

static struct  arbre   *mettre(struct galaxie *M, struct arbre *A,
                        int (*comp)(struct galaxie *, struct galaxie *) )
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
            A->FG = mettre(M, A->FG, comp);
        else
            A->FD = mettre(M, A->FD, comp);
    }
    return(A);

}

/***************************************************************/

static void lecture(struct arbre *R, struct galaxie *B, long int *i)
{
    if (R != NULL)
    {
        lecture(R->FG, B, i);
        B[(*i)++] = *(R->N);
        lecture(R->FD, B, i);
    }
}

/***************************************************************/

int comparer_z(struct galaxie *A, struct galaxie *B)
{
    return(A->z <= B->z);
}

/***************************************************************/

int comparer_tau(struct galaxie *A, struct galaxie *B)
{
    return(A->tau  > B->tau);
}

int comparer_pos(struct galaxie *A, struct galaxie *B)
{
    return(A->C.x > B->C.x);
}
