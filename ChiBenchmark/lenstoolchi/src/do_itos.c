#include <math.h>
#include <constant.h>
#include <fonction.h>

/*
*       nom:        do_itos - lesn-tool
*       auteur:     Jean-Paul Kneib
*       date:       oct 94
*       place:      IoA cambridge
*/

/*
* transform a list of pixels of an image to the source plan
* quick version - PSF = dirac
*
* Global variables used :
* - M, imFrame, ps
* - in e_dpl() : G, lens, lens_table
* - in err_invim() : iline, radial, tangent, nrline, ntline
*/

void    do_itos(double **im, struct pixlist *pl, int npl, double dlsds,
                double **source, double **erreur,
                int **imult)
{
    extern  struct g_pixel  imFrame;
    register    int i, j, k, is, js;
    const extern  struct  g_mode  M;
    const extern  struct g_pixel  ps;
    struct  point   A, B;
    double  xx;

    double xpix, ypix, xw, yw;
    int     offscl;  // 0 if OK, 1 if offscale

    // TO FORCE TO USE THE FIRST LOOP
    // TODO : CORRECT THE SECOND LOOP (WCS)
    imFrame.wcsinfo = NULL;

    if (imFrame.wcsinfo == NULL)
    {
        /*According to the subsampling in the image plane
         *for each subpixel in the image plane*/
        for (k = 0; k < npl; k++)
            for (i = 0; i < imFrame.ech; i++)
                for (j = 0; j < imFrame.ech; j++)
                {
                    xx = (double)imFrame.ech;
                    A.y = (double)pl[k].i + imFrame.pixely * (-.5 + (double)i / xx);
                    A.x = (double)pl[k].j + imFrame.pixelx * (-.5 + (double)j / xx);

                    /*convert relative wc in arcsec from image plane to source plane*/
                    e_dpl(&A, dlsds, &B);

                    /*convert from relative wc in arcsec to pixel in the source image*/
                    is = (int) ((B.y - ps.ymin) / ps.pixely);
                    js = (int) ((B.x - ps.xmin) / ps.pixelx);

                    /*fill the multiplicity, source and error images*/
                    if (is >= 0 && js >= 0 && is < ps.ny && js < ps.nx)
                    {
                        imult[is][js] += 1;
                        source[is][js] = ((imult[is][js] - 1) * source[is][js] +
                                          pl[k].flux)
                                         / ((double)imult[is][js]);
                        erreur[is][js] = ((imult[is][js] - 1) * erreur[is][js] +
                                          pl[k].flux * pl[k].flux)
                                         / ((double)imult[is][js]);
                    }

                }   /*end of for each subpixel in the image plane*/

    } /*end of if(imFrame.wcsinfo==NULL)*/
    else
    {
        for (k = 0; k < npl; k++)
        {
            A.y = imFrame.ymin + imFrame.pixely * (((double) pl[k].i) - .5 * (imFrame.ech - 1.) / imFrame.ech);
            xx = imFrame.xmin + imFrame.pixelx * (((double) pl[k].j) - .5 * (imFrame.ech - 1.) / imFrame.ech);

            for (i = 0; i < imFrame.ech; i++)
            {
                A.x = xx;
                for (j = 0; j < imFrame.ech; j++)
                {
                    e_dpl(&A, dlsds, &B);
                    if (ps.wcsinfo != NULL)
                    {
                        xw = B.x / (-3600.) / cos(M.ref_dec * DTR) + M.ref_ra;
                        yw = B.y / 3600. + M.ref_dec;
                        wcs2pix(ps.wcsinfo, xw, yw, &xpix, &ypix, &offscl);
                        is = (int)xpix;
                        js = (int)ypix;
                    }
                    else
                    {
                        is = (int) ((B.y - ps.ymin) / ps.pixely);
                        js = (int) ((B.x - ps.xmin) / ps.pixelx);
                    }
                    if ((is >= 0) && (js >= 0) && (is < ps.ny) && (js < ps.nx))
                    {
                        imult[is][js] += 1;
                        source[is][js] = ((imult[is][js] - 1) * source[is][js] +
                                          pl[k].flux)
                                         / ((double)imult[is][js]);
                        erreur[is][js] = ((imult[is][js] - 1) * erreur[is][js] +
                                          pl[k].flux * pl[k].flux)
                                         / ((double)imult[is][js]);
                    };
                    A.x += imFrame.pixelx / imFrame.ech * ((double) j);
                };
                A.y += imFrame.pixely / imFrame.ech * ((double) i);
            };
        };
    } /*end of if imFrame.wcsinfo!=NULL*/

    for (is = 0; is < ps.ny; is++)
        for (js = 0; js < ps.nx; js++)
            erreur[is][js] -= source[is][js] * source[is][js];

}
