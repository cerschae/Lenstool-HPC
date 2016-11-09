#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        critic              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/
/* Parameters
 * xmin, ymin, xmax, ymax define the size of the working image
 *
 * Global variables used :
 * - CL, radial, tangent, lens, nrline, ntline, flagr, flagt, G
 * - in e_zeroamp() : G, lens
 * - in next() : G, lens
 * - in e_dpl() : G, lens, lens_table
 * - in dratio() : C
 * - in chsigne() : G, lens
 * - in follow() : CL, radial, nrline, flagr, G, lens, lens_table
 * - in followi() : CL, ntline, flagt, tangent, G, lens, lens_table
 * */
#undef NPOINT
#define NPOINT 4

void    criticinv(double xmin, double ymin, double xmax, double ymax)
{
    const extern struct g_cline CL;
    extern struct biline  radial[], tangent[];
    const extern struct pot     lens[];
    const extern struct g_grille  G;
    extern int nrline, ntline, flagr, flagt;

    struct point A, B, DPL, N, O, OI, OS, P, Q;
    struct point LP[NLMAX+1];   /*path from the bottom left hand corner to the
                                farthest lens*/
    struct point LDPL[NLMAX+1]; /*distance in units of ppas in x and y directions
                                between 2 consecutive lenses*/
    int lnpas[NLMAX+1];
    int nu[NLMAX+1];    /*contains the lens order on the LP path*/
    int fsel[NLMAX+1];  /*keep in memory that nu[i] lens has been considered
                            in the LP path*/
    double ppas;    /* sampling size.(Diag/NPOINT)*/
    double dc;
    double Dmin;
    double J;
    double dl0s, dos, dlsds;  // distances between lens[0] and CL.nz[i], and ratio
    register int j, k, ii;
    register int i;
    int signe_flag;


    /* initialisation de quelques constantes */

    nrline = 0;
    flagr = 1;
    ntline = 0;
    flagt = 1;
    J = NPOINT;
    ppas = sqrt((xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin)) / NPOINT;

    /* liste des points a suivre */
    /* if only 1 lens*/
    if (G.nlens == 1)
    {
        LP[0] = lens[0].C;
        LDPL[0].x = Max(xmax - lens[0].C.x, xmin - lens[0].C.x) / J;
        LDPL[0].y = Max(ymax - lens[0].C.y, ymin - lens[0].C.y) / J;
        if ((lens[0].type == 5) || (lens[0].type == 7) ||
            (lens[0].type == 0) || (lens[0].type == 1))
        {
            LP[0].x += LDPL[0].x;
            LP[0].y += LDPL[0].y;
        };

        lnpas[0] = NPOINT;
    } /*end of if only 1 lens*/
    else
    {
        LP[0].x = xmin; /*LP is the path from the bottom left hand corner to the farthest lens*/
        LP[0].y = ymin;
        Dmin = 9999.;
        for (i = 0; i < G.nlens; i++)
        {
            nu[i] = 99;
            fsel[i] = 0;
        };

        /* look for the closest lens from the bottom left hand corner*/
        for (i = 0; i < G.nlens; i++)
        {
            dc = dist(LP[0], lens[i].C);
            if (dc < Dmin)
            {
                Dmin = dc; /* Dmin contains the minimum distance from the corner*/
                nu[0] = i; /* nu[0] contains the index for the 1st closest lens*/
            };
        };

        LP[1] = lens[nu[0]].C;
        fsel[nu[0]] = 1;

        /*for each lens*/
        for (i = 1; i < G.nlens; i++)
        {
            Dmin = 9999.;
            /*find the closest lens from the LP[i] point*/
            for (j = 0; j < G.nlens; j++)
            {
                if (fsel[j] != 1)
                {
                    dc = dist(LP[i], lens[j].C);
                    if (dc < Dmin)
                    {
                        Dmin = dc;
                        nu[i] = j; /*keep the lens order on the LP path*/
                    };
                };
            };

            LP[i+1] = lens[nu[i]].C;
            fsel[nu[i]] = 1;  /*keep in memory that nu[i] lens has been considered in the LP path*/
        }; /*end of for each lens*/

        /*for each lens*/
        for (i = 0; i < G.nlens; i++)
        {
            if ((lens[nu[i]].type == 5) || (lens[nu[i]].type == 7) ||
                (lens[nu[i]].type == 0) || (lens[nu[i]].type == 1))
            {
                LP[i+1].x += 0.01;
                LP[i+1].y += 0.015;
            };
            /*number of ppas in x and y between 2 consecutive lenses*/
            lnpas[i] = (int) (dist(LP[i], LP[i+1]) / 2. / ppas);
            LDPL[i].x = (LP[i+1].x - LP[i].x) / lnpas[i];
            LDPL[i].y = (LP[i+1].y - LP[i].y) / lnpas[i];
        };
    }; /*end of there are more than 1 lens*/

    /*for each critical line to draw*/
    for (k = 0; k < CL.nplan; k++)
    {
        dl0s = distcosmo2(lens[0].z, CL.cz[k]);
        dos = distcosmo1(CL.cz[k]);
        dlsds = dl0s / dos;

        /*for each lens*/
        for (ii = 0; ii < G.nlens; ii++)
        {
            A = LP[ii];
            DPL = LDPL[ii];
            /*for each point along the path*/
            for (i = 1; i < lnpas[ii]; A = B, i++)
            {
                B.x = A.x + DPL.x;
                B.y = A.y + DPL.y;
                signe_flag = chsigne(A, B, dl0s, dos, CL.cz[k]);

                if (signe_flag == 1)
                    /*search for the radial critical line*/
                {
                    radial[nrline].i = flagr;
                    OI = O = radial[nrline].I = e_zeroamp(A, B, dl0s, dos, CL.cz[k]);
                    e_dpl(&O, dlsds, &OI);
                    radial[nrline++].S = OI;
                    N.x = O.x - DPL.y / 2.;
                    N.y = O.y + DPL.x / 2.;
                    radial[nrline].i = flagr;
                    P = radial[nrline].I = next(N, O, CL.cpas / 2., dl0s, dos, CL.cz[k]);
                    e_dpl(&P, dlsds, &radial[nrline++].S);
                    radial[nrline].i = flagr;
                    Q = radial[nrline].I = next(O, P, CL.cpas / 2., dl0s, dos, CL.cz[k]);
                    e_dpl(&Q, dlsds, &radial[nrline++].S);
                    follow(P, Q, O, dl0s, dos, CL.cz[k]);
                    radial[nrline].i = flagr;
                    radial[nrline].I = OI;
                    radial[nrline++].S = OS;
                    flagr++;
                } /*end of if (signe_flag==1)*/

                else if (signe_flag == 2)
                {
                    tangent[ntline].i = flagt;
                    OI = O = tangent[ntline].I = e_zeroamp(A, B, dl0s, dos, CL.cz[k]);
                    e_dpl(&O, dlsds, &OS);
                    tangent[ntline++].S = OS;
                    N.x = O.x - DPL.y / 2.;
                    N.y = O.y + DPL.x / 2.;
                    tangent[ntline].i = flagt;
                    P = tangent[ntline].I = next(N, O, CL.cpas, dl0s, dos, CL.cz[k]);
                    e_dpl(&P, dlsds, &tangent[ntline++].S);
                    tangent[ntline].i = flagt;
                    Q = tangent[ntline].I = next(O, P, CL.cpas, dl0s, dos, CL.cz[k]);
                    e_dpl(&Q, dlsds, &tangent[ntline++].S);
                    followi(P, Q, O, dl0s, dos, CL.cz[k]);
                    tangent[ntline].i = flagt;
                    tangent[ntline].I = OI;
                    tangent[ntline++].S = OS;
                    flagt++;
                }; /*end of if (signe_flag==2) */
            }; /*end of for(i=1;i<lnpas[ii];A=B,i++)*/
        }; /*end of for(ii=0;ii<G.nlens;ii++)*/
    }; /*end of for(k=0;k<CL.nplan;k++)*/

}
