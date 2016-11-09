#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        grid                */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse
 * Initialize the gimage global variable in a rectangular shape.
 * gimage is defined with G.ngrid^2 (number in grille section)
 * coordinates that map the field of study defined in the champ
 * section.
 *
 ***************************************************************/
void    grid()
{
    const extern  struct  g_mode          M;
    const extern  struct  g_grille    G;
    const extern  struct  g_frame     F;
    extern  struct  point   gimage[NGGMAX][NGGMAX];

    double  I, J, N;
    register    int i, j;

    NPRINTF(stderr, "SET: Grid XY %dx%d\n", G.ngrid, G.ngrid);

    for (i = 0; i < G.ngrid; i++)
        for (j = 0; j < G.ngrid; j++)
        {
            I = i;
            J = j;
            N = G.ngrid - 1;
            I = I / N;
            J = J / N;
            gimage[i][j].x = F.xmin + I * (F.xmax - F.xmin);
            gimage[i][j].y = F.ymin + J * (F.ymax - F.ymin);
        };

}
