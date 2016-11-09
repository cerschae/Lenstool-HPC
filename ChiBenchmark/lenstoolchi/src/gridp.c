#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        gridp               */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/


void    gridp()
{
    const extern struct g_mode          M;
    const extern struct g_grille  G;
    const extern struct g_frame       F;
    extern struct point gimage[NGGMAX][NGGMAX];

    double  I, J, N;
    register    int i, j;

    NPRINTF(stderr, "SET: Grid Polar %dx%d\n", G.ngrid, G.ngrid);

    for (i = 0; i < G.ngrid; i++)
        for (j = 0; j < G.ngrid; j++)
        {
            I = i;
            J = j;
            N = G.ngrid - 1;
            I = I / N + 0.0001;
            J = J / N;
            gimage[i][j].x = I * F.rmax * cos(J * 2.*PI);
            gimage[i][j].y = I * F.rmax * sin(J * 2.*PI);
        }
}
