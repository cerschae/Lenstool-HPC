#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        classer             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    classer(struct galaxie image[NFMAX][NIMAX],
                struct galaxie *cimage, long int *ni, long int *ncistart )
{
    const extern  struct  g_source    S;
    const extern  struct  g_large     L;

    int    i, j, k;
    char   ngiant[NFMAX][IDSIZE];

    k = 0;
    for ( i = 0 ; i < S.ns ; i++ )
        for ( j = 0 ; j < NIMAX && strcmp(image[i][j].n, "") ; j++ )
            cimage[k++] = image[i][j];

    *ni = k;

    // note tous les numeros des sources donnant des arcs geants
    for (i = 0; i < *ni && cimage[i].tau > L.dlarge ; i++);
        strcpy(ngiant[i], cimage[i].n); 

    *ncistart = i;

}
