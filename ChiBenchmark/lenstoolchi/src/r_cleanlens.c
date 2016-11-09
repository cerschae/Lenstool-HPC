#include<stdio.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<lt.h>

/****************************************************************/
/*                 Program  : r_invim               */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/****************************************************************/

void r_cleanlens(FILE *IN, FILE *OUT)
{
    extern  struct  g_mode      M;
    extern  struct  g_pixel     ps, imFrame,wFrame, PSF;
    extern  struct  g_cube      cubeFrame;
    char    second[20], file[50], third[FILENAME_SIZE+10];
    double  reel;
    int i = 0;

    /* You can oversample your source plane image compared to the sampling
     * of your image plane image. (eg. 1 pixel out of 2 is 0 in source plane)
     */

    fmot(IN, second);
    while (strcmp(second, "end"))
    {
        flire(IN, third);
        CHECK_THIRD(FILENAME_SIZE+10)

        if (!strcmp(second, "cleanset"))
        {
            sscanf(third, "%d%lf", &M.iclean, &M.zclean);
            fprintf(OUT, "\t%s\t%d %lf\n", second, M.iclean, M.zclean);
        }
        else if (!strcmp(second, "imframe"))
        {
            sscanf(third, "%d%s", &imFrame.format, imFrame.pixfile);
            fprintf(OUT, "\t%s\t%d %s\n", second, imFrame.format, imFrame.pixfile);
        }
        else if (!strcmp(second, "wframe"))
       	{
            sscanf(third, "%d%s", &wFrame.format, wFrame.pixfile);
            fprintf(OUT, "\t%s\t%d %s\n", second, wFrame.format, wFrame.pixfile);
        }
        else if (!strcmp(second, "cubeframe"))
       	{
            sscanf(third, "%d%s", &cubeFrame.format, cubeFrame.pixfile);
            fprintf(OUT, "\t%s\t%d %s\n", second, cubeFrame.format, cubeFrame.pixfile);
        }
        else if (!strcmp(second, "psfframe"))
        {
            sscanf(third, "%d%s", &PSF.format, PSF.pixfile);
            fprintf(OUT, "\t%s\t%d %s\n", second, PSF.format, PSF.pixfile);
        }
        else if (!strcmp(second, "sframe"))
        {
            sscanf(third, "%s", ps.pixfile);
            fprintf(OUT, "\t%s\t %s\n", second, ps.pixfile);
        }
        else if (!strcmp(second, "c_image"))
        {
            sscanf(third, "%s", M.centerfile);
            fprintf(OUT, "\t%s\t %s\n", second, M.centerfile);
        }
        else if (!strcmp(second, "ncont"))
        {
            sscanf(third, "%d%s", &imFrame.ncont, imFrame.outfile);
            fprintf(OUT, "\t%s\t%d %s\n", second, imFrame.ncont, imFrame.outfile);
        }
        else if (!strcmp(second, "contour"))
        {
            sscanf(third, "%d%s", &i, file);
            if( i > 0 )
            {
                strcpy(imFrame.contfile[i-1], file);
                fprintf(OUT, "\t%s\t%d %s\n", second, i, imFrame.contfile[i-1]);
            }
        }
        else if (!strcmp(second, "column"))
        {
            sscanf(third, "%d", &imFrame.column);
            fprintf(OUT, "\t%s\t%d\n", second, imFrame.column);
        }
        else if (!strcmp(second, "echant"))
        {
            sscanf(third, "%d", &imFrame.ech);
            fprintf(OUT, "\t%s\t%d\n", second, imFrame.ech);
        }
        else if (!strcmp(second, "s_echant"))
        {
            sscanf(third, "%d", &ps.ech);
            fprintf(OUT, "\t%s\t%d\n", second, ps.ech);
        }
        else if (!strcmp(second, "s_n"))
        {
            sscanf(third, "%d", &ps.nx);
            fprintf(OUT, "\t%s\t%d\n", second, ps.nx);
            ps.ny = ps.nx;
        }
        else if (!strcmp(second, "header"))
        {
            sscanf(third, "%d", &imFrame.header);
            fprintf(OUT, "\t%s\t%d\n", second, imFrame.header);
        }
        else if (!strcmp(second, "pixel"))
        {
            sscanf(third, "%lf", &reel);
            imFrame.pixelx = imFrame.pixely = reel;
            fprintf(OUT, "\t%s\t%lf\n", second, imFrame.pixelx);
        }
        else if (!strcmp(second, "pixelx"))
        {
            sscanf(third, "%lf", &imFrame.pixelx);
            fprintf(OUT, "\t%s\t%lf\n", second, imFrame.pixelx);
        }
        else if (!strcmp(second, "pixely"))
        {
            sscanf(third, "%lf", &imFrame.pixely);
            fprintf(OUT, "\t%s\t%lf\n", second, imFrame.pixely);
        }
        else if (!strcmp(second, "xmin"))
        {
            sscanf(third, "%lf", &imFrame.xmin);
            fprintf(OUT, "\t%s\t%lf\n", second, imFrame.xmin);
        }
        else if (!strcmp(second, "ymin"))
        {
            sscanf(third, "%lf", &imFrame.ymin);
            fprintf(OUT, "\t%s\t%lf\n", second, imFrame.ymin);
        }
        else if (!strcmp(second, "s_xmin"))
        {
            sscanf(third, "%lf", &ps.xmin);
            fprintf(OUT, "\t%s\t%lf\n", second, ps.xmin);
        }
        else if (!strcmp(second, "s_ymin"))
        {
            sscanf(third, "%lf", &ps.ymin);
            fprintf(OUT, "\t%s\t%lf\n", second, ps.ymin);
        }
        else if (!strcmp(second, "s_xmax"))
        {
            sscanf(third, "%lf", &ps.xmax);
            fprintf(OUT, "\t%s\t%lf\n", second, ps.xmax);
        }
        else if (!strcmp(second, "s_ymax"))
        {
            sscanf(third, "%lf", &ps.ymax);
            fprintf(OUT, "\t%s\t%lf\n", second, ps.ymax);
        };
        // Read the next line
        fmot(IN, second);
    }

    fprintf(OUT, "\t%s\n", second);

}
