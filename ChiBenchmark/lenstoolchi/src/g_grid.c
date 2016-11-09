#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        g_grid              */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/****************************************************************/

void    g_grid(int igrid, int ngrid, double zgrid)
{
    const extern  struct  g_mode          M;
    const extern  struct  pot     lens[];
    const extern  struct  g_grille    G;
    const extern  struct  g_frame     F;
    extern struct point gsource_global[NGGMAX][NGGMAX];
   
   //const extern    struct  chaine  *Tsol;
   
    
    struct bitriplet Tsol[NIMAX];
    register    int i, j, k;
    int nimage, kbest;


    double  xmin, xmax, ymin, ymax, I, J, N;
    double  alpha = 1.8;

    /*coefficient recducteur du lieu de tirage dans le plan source*/
    double  dlsds, dt, Dmin;
    FILE    *OUT;
    struct  point   **gi, **gs, gimage[NIMAX];

    NPRINTF(stderr, "COMP: grids in the Source and Image plane\n");

    dlsds = dratio(lens[0].z, zgrid);
    ngrid = Min(NGMAX, ngrid);

    gi = (struct point **)al_sq_point(ngrid, ngrid);
    gs = (struct point **)al_sq_point(ngrid, ngrid);

    xmin = F.xmin + .0001;
    ymin = F.ymin + .0001;
    xmax = F.xmax - .0001;
    ymax = F.ymax - .0001;

    /*ie. grid in the source plane to the image plane*/
    if (igrid == 1)
    {
        e_unlensgrid(gsource_global, dlsds);
        for (i = 0; i < ngrid; i++)
        {
            for (j = 0; j < ngrid; j++)
            {
                I = i;
                J = j;
                N = ngrid - 1;
                I = I / N;
                J = J / N;
                if (G.pol == 0)
                {
                    gs[i][j].x = (xmin + I * (xmax - xmin)) / alpha;
                    gs[i][j].y = (ymin + J * (ymax - ymin)) / alpha;
                }
                else
                {
                    I = I + 0.0001;
                    gs[i][j].x = I * F.rmax * cos(J * 2.*PI) / alpha;
                    gs[i][j].y = I * F.rmax * sin(J * 2.*PI) / alpha;
                };
                //Tsol=NULL;
                nimage = inverse((const struct point (*)[NGGMAX]) gsource_global, &gs[i][j], Tsol);
                if (nimage > 0)
                    nimage = e_test_P(Tsol, nimage, &gs[i][j], gimage, dlsds, 0.00001);
                Dmin = 100000.;
                kbest = 0;
                for (k = 0; k < nimage; k++)
                {
                    if (j == 0)
                        kbest = k;
                    else
                    {
                        dt = fabs(e_pot(gimage[k], dlsds) - e_pot(gi[i][j-1], dlsds));
                        if (dt < Dmin)
                        {
                            Dmin = dt;
                            kbest = k;
                        };
                    };
                };
                gi[i][j] = gimage[kbest];
            };
        };
    }
    /*if igrid == 2 ie. grid in the image plane to the source plane*/
    else
    {
        for (i = 0; i < ngrid; i++)
            for (j = 0; j < ngrid; j++)
            {
                I = i;
                J = j;
                N = ngrid - 1;
                I = I / N;
                J = J / N;
                if (G.pol == 0)
                {
                    gi[i][j].x = (xmin + I * (xmax - xmin));
                    gi[i][j].y = (ymin + J * (ymax - ymin));
                }
                else
                {
                    I = I + 0.0001;
                    gi[i][j].x = I * F.rmax * cos(J * 2.*PI);
                    gi[i][j].y = I * F.rmax * sin(J * 2.*PI);
                };

                e_dpl(&gi[i][j], dlsds, &gs[i][j]);
            };
    };

    OUT = fopen("gs1.dat", "w");
    if ( M.iref == 3 )
        fprintf(OUT, "#REFERENCE 3 %lf %lf\n", M.ref_ra, M.ref_dec);
    else
        fprintf(OUT, "#REFERENCE 0");

    for (i = 0; i < ngrid; i++)
    {
        for (j = 0; j < ngrid; j++)
            fprintf(OUT, "%d\t%lf\t%lf\n", j, gs[i][j].x, gs[i][j].y);


        if (++i < ngrid)
            for (j = ngrid - 1; j >= 0; j--)
                fprintf(OUT, "%d\t%lf\t%lf\n", j, gs[i][j].x, gs[i][j].y);

    };
    fclose(OUT);

    OUT = fopen("gs2.dat", "w");
    if ( M.iref == 3 )
        fprintf(OUT, "#REFERENCE 3 %lf %lf\n", M.ref_ra, M.ref_dec);
    else
        fprintf(OUT, "#REFERENCE 0");

    for (j = 0; j < ngrid; j++)
    {
        for (i = 0; i < ngrid; i++)
            fprintf(OUT, "%d\t%lf\t%lf\n", i, gs[i][j].x, gs[i][j].y);

        if (++j < ngrid)
            for (i = ngrid - 1; i >= 0; i--)
                fprintf(OUT, "%d\t%lf\t%lf\n", i, gs[i][j].x, gs[i][j].y);

    };
    fclose(OUT);

    OUT = fopen("gi1.dat", "w");
    if ( M.iref == 3 )
        fprintf(OUT, "#REFERENCE 3 %lf %lf\n", M.ref_ra, M.ref_dec);
    else
        fprintf(OUT, "#REFERENCE 0");

    for (i = 0; i < ngrid; i++)
    {
        for (j = 0; j < ngrid; j++)
            fprintf(OUT, "%d\t%lf\t%lf\n", j, gi[i][j].x, gi[i][j].y);

        if (++i < ngrid)
            for (j = ngrid - 1; j >= 0; j--)
                fprintf(OUT, "%d\t%lf\t%lf\n", j, gi[i][j].x, gi[i][j].y);
    };
    fclose(OUT);

    OUT = fopen("gi2.dat", "w");
    if ( M.iref == 3 )
        fprintf(OUT, "#REFERENCE 3 %lf %lf\n", M.ref_ra, M.ref_dec);
    else
        fprintf(OUT, "#REFERENCE 0");

    for (j = 0; j < ngrid; j++)
    {
        for (i = 0; i < ngrid; i++)
            fprintf(OUT, "%d\t%lf\t%lf\n", i, gi[i][j].x, gi[i][j].y);

        if (++j < ngrid)
            for (i = ngrid - 1; i >= 0; i--)
                fprintf(OUT, "%d\t%lf\t%lf\n", i, gi[i][j].x, gi[i][j].y);
    };
    fclose(OUT);
}
