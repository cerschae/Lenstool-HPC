#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

static void Tsup(const struct point grille[NGGMAX][NGGMAX], int i, int j, struct triplet *T);
static void Tinf(const struct point grille[NGGMAX][NGGMAX], int i, int j, struct triplet *T);

/****************************************************************/
/*      nom:        inverse             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse
 * Create a link sequence of chaine structure. Each node contain
 * a triangle that contain the source P, the corresponding triangle
 * in the image plane and a link to another arclet of the same
 * familly.
 *
 * Return the number of images found in the gsource 2D map for the
 * source at position P in the source plane.
 *
 * The fist arclet is link to the Null pointer (Tsol). At the end, Tsol
 * points to the last arclet.
 *
 * Parameters :
 * - gsource : see description in e_unlensgrid.c
 * - P : a point in the source plane.
 * - Tsol : a list of arclets for a family
 ****************************************************************/
int inverse(const struct point gsource[][NGGMAX], struct point *P, struct bitriplet Tsol[NIMAX])
{
    const extern  struct  g_grille    G;
    const extern  struct  point   gimage[NGGMAX][NGGMAX];
    //const extern    struct  bitriplet   Tsol[NIMAX];
    int     nimage;

    register    int i, j;

    struct  triplet A, B; // Triplets Tsup and Tinf in source plane
    //struct    chaine  *maillon;

    //maillon = NULL;
    nimage = 0;

    for (i = 0 ; i < G.ngrid - 1 && nimage < NIMAX ; i++)
        for (j = 0 ; j < G.ngrid - 1 && nimage < NIMAX ; j++)
        {
            Tsup(gsource, i, j, &A);
            Tinf(gsource, i, j, &B);
            if (insidebord(P, &A))
            {
                //maillon=(struct chaine *)malloc(sizeof(struct chaine));
                Tsup(gimage, i, j, &Tsol[nimage].i);
                Tsol[nimage].s = A;
                //maillon->S=A;
                //maillon->F=Tsol;
                //Tsol=maillon;
                nimage++;
            }
            else if (inside(P, &B))
            {

                //maillon=(struct chaine *)malloc(sizeof(struct chaine));
                Tinf(gimage, i, j, &Tsol[nimage].i);
                Tsol[nimage].s = B;
                //maillon->S=B;
                //maillon->F=Tsol;
                //Tsol=maillon;
                nimage++;
            };
        };

    /*  if (nimage>0)   //you can't delete it as it is referenced in Tsol
            free((struct chaine *) maillon);*/

    return(nimage);
}

/*********************************************************************/
/* Return the 3 coordinates {a,b,c} contrained in grille at the 3 corners
 * of the looking upwards triangle {(i,j),(i+1,j),(i,j+1)}
 * Parameters :
 * - grille : a square grid of points of size NGGMAX
 * - i,j : integer coordinates in the grid
 *  */
static void Tsup(const struct point grille[NGGMAX][NGGMAX], int i, int j, struct triplet *T)
{
    T->a = grille[i][j];
    T->b = grille[i+1][j];
    T->c = grille[i][j+1];
}

/*********************************************************************/
/* Return the 3 coordinates {a,b,c} contrained in grille at the 3 corners
 * of the looking downwards triangle {(i+1,j+1),(i+1,j),(i,j+1)}
 * Parameters :
 * - grille : a square grid of points of size NGGMAX
 * - i,j : integer coordinates in the grid
 *  */
static void Tinf(const struct point grille[NGGMAX][NGGMAX], int i, int j, struct triplet *U)
{
    U->a = grille[i+1][j+1];
    U->b = grille[i+1][j];
    U->c = grille[i][j+1];
}
/*********************************************************************/
