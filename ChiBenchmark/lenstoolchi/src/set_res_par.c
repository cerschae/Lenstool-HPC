#include<stdio.h>
#include<float.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include "lt.h"

/****************************************************************/
/*      nom:        set_res_par         */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************
 * For each potential, reset the potential type dependent parameters
 * from the optimised parameters b0 and epot.
 */

void    set_res_par()
{
    extern struct g_mode    M;
    extern struct g_source  S;
    extern struct g_grille  G;
    const extern struct g_cosmo   C;
    extern struct pot   lens[], lmin[], lmax[], prec[];
    extern int block[][NPAMAX];

    register int i;
    //register int ii,jj;
    double test;
    double  d1;

    double GG = 10.867;

//  extern double *v_xx;
//  extern double *v_yy;
//  extern double **map_p;
//  extern double **tmp_p;
//  extern double **map_axx;
//  extern double **map_ayy;

//  char    mode[20],nature[20],comment[1024],type[20];

    /*
    * update the program values for each matter clump
    * to match the entry values
    */

    /*
    * set each clump parameters
    */

    if ( lens[0].type != 10 )
    {
        for (i = 0; i < G.nlens; i++)
        {

            /* calcul du rapport Dls/Ds */
            lens[i].dlsds = dratio(lens[i].z, S.zs);

            d1 = d0 / C.h * distcosmo1(lens[i].z);

            /*
            *  dynamical parameters
            */

            // always rescale rcore, scale_radius, etc.
            lens[i].rckpc = lens[i].rc * d1;
            // and rcut ... if defined!
            if ( lens[i].rcut != DBL_MAX )
                lens[i].rcutkpc = lens[i].rcut * d1;

            switch (lens[i].type)
            {
                case(1):
                    lens[i].sigma = sqrt(lens[i].b0 / 4. / pia_c2);
                    lens[i].ct = lens[i].b0 * lens[i].dlsds;
                    lens[i].cr = 0.;
                    break;
                case(-1):
                    lens[i].sigma = sqrt(lens[i].b0 / 4. / pia_c2);
                    lens[i].ct = lens[i].b0 * lens[i].dlsds;
                    lens[i].cr = 0.;
                    break;
                case(7):
                    lens[i].masse = lens[i].b0 / (4.*RTA * GM_c2) * d1;
                    lens[i].ct = sqrt(lens[i].b0 * lens[i].dlsds);
                    lens[i].cr = 0.;
                    break;
                case(9):
                    lens[i].pmass = lens[i].b0 * cH2piG * C.h / distcosmo1(lens[i].z);
                    lens[i].ct = 0.;
                    lens[i].cr = 0.;
                    break;
                case(5):
                    lens[i].sigma = vol * sqrt(lens[i].b0 / 18. / RTA);
                    break;
                case(8):
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    break;
                case(81):
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    // total mass
                    lens[i].masse = 4 * M_PI / 3 * M_PI / GG * (lens[i].sigma / 1000) * (lens[i].sigma / 1000) *
                                    lens[i].rcut * d0 / C.h * distcosmo1(lens[i].z);

                    break;

                case(82):
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    break;
                case(83):
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    break;
                case(84):
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    break;
                case(85):
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    break;
                case(86):
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    break;
                case(87):
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    break;
                case(89):
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    break;
                case(88):
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    break;
                case(10):
                    break;
                case(12):
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    e_nfw_rs2c(lens[i].sigma, lens[i].rckpc, &lens[i].pmass, &lens[i].beta, &lens[i].masse, lens[i].z);
                    break;
                case(13):
                    lens[i].sigma = lens[i].b0 * 1e12 * cH0_4piG * C.h / distcosmo1(lens[i].z);
                    break;
                case(14):
                    break;
                default:
                    NPRINTF(stderr, "WARN: Clump %d with unknown type --> default: Pseudo-Elliptical Potential with Core Radius\n", i);
                    lens[i].sigma = sqrt(lens[i].b0 / 6. / pia_c2);
                    test = lens[i].dlsds * lens[i].dlsds * lens[i].b0 * lens[i].b0
                           - lens[i].rc * lens[i].rc;

                    if (test > 0.)
                    {
                        lens[i].ct = sqrt(test);
                        lens[i].cr = sqrt(pow(lens[i].b0 * lens[i].dlsds, 2. / 3.) *
                                          pow(lens[i].rc, 4. / 3.) - lens[i].rc * lens[i].rc);
                    }
                    else
                        lens[i].ct = lens[i].cr = 0.;


                    if (lens[i].type > 20)
                        updatecut(i);

                    break;

            } //end of switch lens.type

        } //end of for each potential


        // optimise limits as well
        for ( i = 0; i < G.nplens[0]; i++ )
        {
            d1 = d0 / C.h * distcosmo1(lens[i].z);
            if ( block[i][RC] != 0 )
            {
                lmin[i].rckpc = lmin[i].rc * d1;
                lmax[i].rckpc = lmax[i].rc * d1;
                prec[i].rckpc = prec[i].rc * d1;
            }

            if ( lens[i].rcut != DBL_MAX && block[i][RCUT] != 0 )
            {
                lmin[i].rcutkpc = lmin[i].rcut * d1;
                lmax[i].rcutkpc = lmax[i].rcut * d1;
                prec[i].rcutkpc = prec[i].rcut * d1;
            }
        }

    } // end of spline test

}
