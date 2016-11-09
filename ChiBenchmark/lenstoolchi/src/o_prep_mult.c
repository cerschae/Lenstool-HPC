#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>


/****************************************************************/
/*      nom:        o_prep_mult         */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/*                              */
/* Modified :                           */
/*      EJ (31/08/2005)                 */
/****************************************************************/
#undef SWAPI
#undef SWAPD
#define SWAPI(a,b) {swapi = a; a = b; b = swapi;}
#define SWAPD(a,b) {swapd = a; a = b; b = swapd;}


void    o_prep_mult(int ntmult, struct galaxie *mult)
{

    /*************  declaration de common et locale ****************/

    extern struct g_mode  M;
    extern struct galaxie multi[NFMAX][NIMAX];
    extern struct g_image I;
    extern struct z_lim   zlim[];
    extern struct pot       lens[];
    extern struct cline     cl[];
    int    i, j, l, m;
    int  k;     // number of unknown redshifts
    int    nimages; // number of families in limages
    double swapd;
    int    swapi;
    char     str1[IDSIZE];
    char   limages[ZMBOUND][IDSIZE];
    int    matched;

    NPRINTF(stderr, "INFO: multiple images: %d\n", ntmult);

    /* classement par z croissant */

    sort(ntmult, mult, comparer_z);

    // make that 2 consecutives images don't have the same name
    if ( !strcmp(mult[1].n, mult[0].n) )
        sprintf( mult[0].n, "%s.%d", mult[0].n, 0 );

    j = 0;
    k = 0;
    multi[0][0] = mult[0];

    I.mult[0] = 1;
    NPRINTF(stderr, "%s %d %d %d\n", mult[0].n, 0, 0, I.mult[0]);
    /* k : index d'une source
       j : index de l'image de cette source
       multi[n_src][n_arc] : tableau 2D classe en source/images
       mult  : tableau 1D issu du fichier d'images*/

    // sort the images according to their identifier
    for ( i = 1; i < ntmult && k < NFMAX; i++ )
    {
        // make that 2 consecutives images don't have the same name
        sprintf( str1, "%s.%d", mult[i].n, 0 ); //create again the previous name
        if ( !strcmp(str1, multi[k][0].n) )
            sprintf( mult[i].n, "%s.%d", mult[i].n, j + 1 );

        if (!indexCmp(mult[i].n, multi[k][0].n))
        {
            multi[k][++j] = mult[i];
            I.mult[k] = j + 1;
        }
        else
        {
            j = 0;
            if ( !strcmp( mult[i+1].n, mult[i].n ) )
                sprintf( mult[i].n, "%s.%d", mult[i].n, j);

            multi[++k][0] = mult[i];
            I.mult[k] = 1;
        }
        NPRINTF(stderr, "%s %d %d %d\n", mult[i].n, k, j, I.mult[k]);
    }

    if ( k >= NFMAX && i < ntmult )
    {
        fprintf(stderr, "ERROR: too many systems in %s (maximum %d)\n",
                I.multfile, NFMAX);
        exit(-1);
    }

    // display a log of the images on stderr
    I.n_mult = k + 1;
    for (i = 0; i < I.n_mult; i++)
    {
        NPRINTF(stderr, "\tID%d:Mult %d\n", i, I.mult[i]);
        multi[i][I.mult[i]] = multi[i][0];
        for (j = 0; j <= I.mult[i]; j++)
        {
            NPRINTF(stderr, "\t%s (%.3lf,%.3lf) (%.3lf,%.3lf,%.3lf) %.3lf\n",
                    multi[i][j].n, multi[i][j].C.x, multi[i][j].C.y,
                    multi[i][j].E.a, multi[i][j].E.b, multi[i][j].E.theta,
                    multi[i][j].z);
        }

        // no single images without redshift
        if ( I.mult[i] == 1 && multi[i][0].z == 0 )
        {
            fprintf(stderr, "ERROR: no multiple images in file %s\n",
                    I.multfile);
            exit(-1);
        };
    };

    // display a log on stderr and count the unknown redshift images
    k = 0;
    for (i = 0; i < I.n_mult; i++)
    {
        if (multi[i][0].z < 0.01)
            k++;

        NPRINTF(stderr, "INFO: image %s z=%.3lf unknown=%d\n",
                multi[i][0].n, multi[i][0].z, k);
    }

    if ( k != 0 )
        NPRINTF(stderr, "INFO: %d multiple images with unknown redshift\n", k);

    // Check if there are some critical line redshift optimization
    int ncrit = 0;
    for ( i = 0; i < I.nzlim; i++ )
    {
        int res = splitzmlimit(zlim[i].n, limages);
        for( l = 0; l < res; l++ )
            if( !strncmp(limages[l], "cl", 2) )
                ncrit ++ ;
    }
    if ( ncrit != 0 )
        NPRINTF(stderr, "INFO: %d critical lines with unknown redshift\n", ncrit);

    // Check the number of unknown redshifts
    if ( getNzmlimit() != k + ncrit)
    {
        fprintf(stderr, "ERROR: %d images with unknown redshifts in %s, but only %d optimised with z_m_limit statements\n", k, I.multfile, getNzmlimit());
        exit(-1);
    };

    // Assign the z_m_limit index to the right image with unknown redshift
    // for all the redshift unknown families in multi[k][0] list
    for ( i = 0 ; i < k ; i++ )
    {
        // check that the names in multi[k][0] and zlim[I.nzlim] match together
        l = 0;  // lookup index in the zlim list
        matched = 0;
        while ( !matched && l < I.nzlim )
        {
            // look up in all families in zlim[l].n
            nimages = splitzmlimit( zlim[l].n, limages );
            m = 0;
            // !!!THE ACCESS TO limages[nimages] MAY PRODUCE A SEGMENTATION FAULT
            while ( m < nimages && indexCmp( limages[m], multi[i][0].n ) )
            {
                // test with .0 appended to limages[m] as we did in multi (cf l.68)
                sprintf( str1, "%s.0", limages[m] );
                if ( !indexCmp(multi[i][0].n, str1) )
                    strcpy( limages[m], str1 );
                else
                    m++;
            }
            // pack limages back to zlim[l]
            packzmlimit( limages, nimages, zlim[l].n );

            if ( m < nimages )
                matched = 1;    // we found matching families --> exit
            else
                l++;    // we search in the next zlim[l].n
        }

        if ( !matched && l == I.nzlim )
        {
            fprintf(stderr, "ERROR: Image name %s not found in the z_m_limit list\n",
                    multi[i][0].n);
            exit(-1);
        }

        /*
        // The matching families are in zlim[l] and in multi[i]
        // Copy zlim[l] to zlim[i] (ie. sort the zlim
        // list in the multi array order).
        if( l != i )
        {
            // swap the zlim names between index l and i
            strcpy( str1, zlim[l].n );
            strcpy( zlim[l].n, zlim[i].n );
            strcpy( zlim[i].n, str1 );

            // swap the zlim info
            SWAPI(zlim[l].opt, zlim[i].opt);
            SWAPI(zlim[l].bk, zlim[i].bk);
            SWAPD(zlim[l].min, zlim[i].min);
            SWAPD(zlim[l].max, zlim[i].max);
            SWAPD(zlim[l].dderr, zlim[i].dderr);
        }
        */
        NPRINTF(stderr, "INFO: image name:%s zl name:%s\n", multi[i][0].n, zlim[l].n);
    }

    // Append dratio information to each zlim and multi elements.
    for ( i = 0 ; i < I.nzlim ; i++ )
    {
        zlim[i].ddmin = dratio(lens[0].z, zlim[i].min);
        zlim[i].ddmax = dratio(lens[0].z, zlim[i].max);

        nimages = splitzmlimit( zlim[i].n, limages );
        for ( l = 0; l < nimages ; l++ )
        {
            // for critical lines
            if (!strncmp(limages[l], "cl", 2))
            {
                sscanf(limages[l], "cl%d", &m);
                if( zlim[i].bk == 1)
                {
                    cl[m].z = .5 * (zlim[i].min + zlim[i].max);
                    NPRINTF(stderr, "INFO: z_m limits cl%d : cur=%lf min=%lf max=%lf err=%lf\n", m, cl[m].z, zlim[i].min, zlim[i].max, zlim[i].dderr);
                }
                else if( zlim[i].bk == 3 || zlim[i].bk == 0. )
                {
                    cl[m].z = zlim[i].min;       // mu
                    NPRINTF(stderr, "INFO: z_m limits cl%d : cur=%lf mean=%lf sigma=%lf err=%lf\n", m, cl[m].z, zlim[i].min, zlim[i].max, zlim[i].dderr);
                }
            }
            else
            {
                // for multiple images
                m = 0;
                while ( indexCmp( limages[l], multi[m][0].n ) && strncmp(limages[l], "cl", 2) ) m++;
                if ( zlim[i].bk == 1 )
                {
                    // flat prior
                    for( j = 0; j < I.mult[m]; j++)
                        multi[m][j].z  = .5 * (zlim[i].min + zlim[i].max);
                    NPRINTF(stderr, "INFO: z_m limits %s : cur=%lf min=%lf max=%lf err=%lf\n", multi[m][0].n, multi[m][0].z, zlim[i].min, zlim[i].max, zlim[i].dderr);
                }
                else if (zlim[i].bk == 3 || zlim[i].bk == 0. )
                {
                    // gaussian prior
                    for( j = 0; j < I.mult[m]; j++)
                        multi[m][j].z  = zlim[i].min;       // mu
                    NPRINTF(stderr, "INFO: z_m limits %s : cur=%lf mean=%lf sigma=%lf err=%lf\n", multi[m][0].n, multi[m][0].z, zlim[i].min, zlim[i].max, zlim[i].dderr);
                }
            }
        }
    }

    // Compute dratio for the images
    for (i = 0; i < I.n_mult; i++)
        for(j = 0; j < I.mult[i]; j++)
        {
                if ( lens[0].z < multi[i][j].z )
                    multi[i][j].dl0s = distcosmo2(lens[0].z, multi[i][j].z);
                else
                    multi[i][j].dl0s = 0.;
        
                multi[i][j].dos = distcosmo1(multi[i][j].z);
                multi[i][j].dr = multi[i][j].dl0s / multi[i][j].dos;
        }    

    for (i = 0; i < I.n_mult; i++)
        NPRINTF(stderr, "INFO: image %s z=%lf dos=%lf dl0s=%lf dlsds=%lf\n",
                multi[i][0].n, multi[i][0].z, multi[i][0].dos, multi[i][0].dl0s,multi[i][0].dr);
}

/* In the case of bounded images for the z_m_limit redshift optimization,
 * split the N bounded names in a list of names.
 * Return the number of bounded families.
 */
int splitzmlimit(char name[IDSIZE], char limage[ZMBOUND][IDSIZE])
{
    unsigned int  i, j, k;

    j = k = 0;
    for ( i = 0; i < strlen(name); i++ )
    {
        if ( name[i] != ' ' )
            // append characters to limages
            limage[j][k++] = name[i];
        else
        {
            // terminate the string limages[j] with \0
            limage[j][k] = '\0';
            j++;
            k = 0; // pass to the next family in name

            // The test j > ZMBOUND has been done in the r_image.c file
        }
    }

    // Terminate the last name
    limage[j][k] = '\0';
    return j+1;
}

/* Pack an array of family names to a single string.
 */
void    packzmlimit(char limages[ZMBOUND][IDSIZE], int nimages, char name[IDSIZE])
{
    int  i;

    strcpy( name, limages[0] );

    for ( i = 1; i < nimages; i++ )
    {
        strcat(name, " ");
        strcat(name, limages[i]);
    }
}

/* Return the total number of families with a redshift to optimize
 */
int getNzmlimit()
{
    extern struct g_image I;
    extern struct z_lim   zlim[];

    char limages[ZMBOUND][IDSIZE];
    int i, res;

    res = 0;
    for ( i = 0; i < I.nzlim; i++ )
        res += splitzmlimit(zlim[i].n, limages);

    return res;
}
