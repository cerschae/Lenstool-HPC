#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        e_unlens            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * Return a catalog of sources from a catalog of arclets for the
 * s_source() function.
 */
void    e_unlens( long int na,
                  struct galaxie *arclet,
                  long int *ns,
                  struct galaxie *source )
{
//    const extern struct g_source S;
//    const extern struct g_image  I;
//    const extern struct pot      lens[];
//    const extern struct galaxie  multi[NFMAX][NIMAX];

    struct ellipse ampli;
    long int    i; //, j;

    // for each arclet
    for ( i = 0 ; i < na ; i++ )
    {
        // assign a redshift to each arclet
//        if (arclet[i].z == 0.)
//        {
//            // look for the multi that corresponds to this arclet
//            for ( j = 0; j < I.n_mult && indexCmp( multi[j][0].n, arclet[i].n ); j++);
//            if ( j < I.n_mult )
//                arclet[i].z = multi[j][0].z;
//            else
//                arclet[i].z = S.zs;
//        }

        // Compute its DLS/DS ratio
//        arclet[i].dr = dratio(lens[0].z, arclet[i].z);
        source[i+(*ns)].c = arclet[i].c;  // point like of extended image?
        ampli = e_unmag_gal(&arclet[i]);

        // if point like image then source shape parameters = 0
        if (arclet[i].c == 's')
            source[i+(*ns)].E.a =
                source[i+(*ns)].E.b =
                    source[i+(*ns)].E.theta = 0.;
        else
        {
            // else compute source shape parameters through amplification
            if (arclet[i].E.b != 0.)
                isoima(&arclet[i].E, &ampli, &source[i+(*ns)].E);
            else
            {
                source[i+(*ns)].E.a =
                    source[i+(*ns)].E.b =
                        arclet[i].E.a * fabs(ampli.a * ampli.b);
                source[i+(*ns)].E.theta = 0.;
            }
        }

        // compute source position
        e_dpl(&arclet[i].C, arclet[i].dr, &source[i+(*ns)].C);
        // ... name
        strcpy(source[i+(*ns)].n, arclet[i].n);
        // ... redshift
        source[i+(*ns)].z = arclet[i].z;
        // ... magnitude
        source[i+(*ns)].mag = arclet[i].mag - 2.5 * log10(fabs(ampli.a * ampli.b));
        // ... DLS/DS ratio
        source[i+(*ns)].dl0s = arclet[i].dl0s;
        source[i+(*ns)].dos = arclet[i].dos;
        source[i+(*ns)].dr = arclet[i].dr;
        // ... intensity
        source[i+(*ns)].I0 = arclet[i].I0;
    }

    (*ns) += na;
}
