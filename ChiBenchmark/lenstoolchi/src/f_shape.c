#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*      nom:        f_shape             */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse            */
/* lecture de fichiers ellipses, source ou image        */
/****************************************************************
 *
 * flag=1 : ID RA DEC A B THETA Z MAG  (default)
 * flag=2 : ID RA DEC A B THETA Z MAG VARE1 VARE2 (for WL or source with Sersic index)
 * flag=3 : ID RA DEC E1 E2 Z VARE1 VARE2 (for WL)
 */
void    f_shape( long int *istart,
                 struct galaxie *liste,
                 char *name, int flag )
{
    const extern struct g_mode          M;
    const extern struct g_image         I; 
    FILE    *IN;
    long int i;
    char    line[256];

    int     iref;
    double  ra, dec;
    int     e_scan = 0;   // error during the line scanning
    int     n_scan;  // number of scanned arguments
    double     dummy;  // check that the line has been fully read
    double  e1, e2;

    i = (*istart);

    // default values
    iref = ra = dec = 0;

    // Read the input arclet file
    NPRINTF(stderr, "READ: %s\n", name);
    IN = fopen(name, "r");

    // Info on the file format
    switch (flag)
    {
        case(2):
            NPRINTF(stderr, "File format %d (ID RA DEC a b theta z mag var_e1 var_e2)\n", flag);
            break;
        case(3):
            NPRINTF(stderr, "File format %d (ID RA DEC e1 e2 z var_e1 var_e2)\n", flag);
            break;    
        default:
            NPRINTF(stderr, "File format %d (ID RA DEC a b theta z mag)\n", flag);
            break;
    }

    if ( IN == NULL || ferror(IN) )
    {
        fprintf( stderr, "ERROR: Error reading %s\n", name);
        exit(-1);
    }

    while ( fgets(line, 256, IN) != NULL && !feof(IN) && !ferror(IN) && !e_scan )
    {

        // Initialise variables
        liste[i].dl0s = liste[i].dos = liste[i].dr = -1;
        liste[i].var1 = liste[i].var2 = 0.;
        liste[i].grad2.a = liste[i].grad2.c = 1;
        liste[i].grad.x = liste[i].grad.y = 0;
        liste[i].np_grad = NULL;
        liste[i].np_grad2a = NULL;
        liste[i].np_grad2b = NULL;
        liste[i].np_grad2c = NULL;
        liste[i].type = 3;   // gaussian profile (default)

        if ( strstr(line, "#REFERENCE" ) != NULL )
        {
            NPRINTF( stderr, "%s", line );
            getRADEC( line, &iref, &ra, &dec );
            continue;
        }

        if ( line[0] == '#' )
            continue;

        if ( flag == 3 )
        {
	// Format: ID RA DEC e1 e2 z var_e1 var_e2
            n_scan = sscanf(line, "%s%lf%lf%lf%lf%lf%lf%lf%lf",
                        liste[i].n, &liste[i].C.x, &liste[i].C.y,
                        &e1, &e2, &liste[i].z, 
                        &liste[i].var1, &liste[i].var2, &dummy);
            if ( n_scan != 8 )  e_scan = 1;

            // Convert flag 3 mode e1, e2 into a, b ,theta lenstool
            liste[i].E.theta = 0.5 * atan2(e2, e1);
            e1 = sqrt(e1 * e1 + e2 * e2);
            liste[i].E.a = sqrt(1. + e1);
            liste[i].E.b = sqrt(1. - e1);
        }
        else if ( flag == 2 )
        {
	// Format: ID RA DEC a b theta z mag var_e1 var_e2
            n_scan = sscanf(line, "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                        liste[i].n, &liste[i].C.x, &liste[i].C.y,
                        &liste[i].E.a, &liste[i].E.b, &liste[i].E.theta,
                        &liste[i].z, &liste[i].mag, 
                        &liste[i].var1, &liste[i].var2, &dummy);
           if( n_scan != 10 ) e_scan = 1;
        }
        else
        {
	// Format: ID RA DEC a b theta z mag 
            n_scan = sscanf(line, "%s%lf%lf%lf%lf%lf%lf%lf%lf",
                        liste[i].n, &liste[i].C.x, &liste[i].C.y,
                        &liste[i].E.a, &liste[i].E.b,
                        &liste[i].E.theta, &liste[i].z, &liste[i].mag, &dummy);
           if( n_scan != 8 ) e_scan = 1;
        }

        if( e_scan == 1 )
        {
            int ncol = 8;  // valid for default and flag=3
            if ( flag == 2 )     ncol = 10;
            fprintf(stderr, "ERROR: reading catalog %s. %d/%d columns found for line %ld\n", name, n_scan, ncol, i + 1);
            exit(-1);
        }

        // convert input to absolute coordinates
        convertXY( &liste[i].C.x, &liste[i].C.y, iref, ra, dec);

        // convert to output relative coordinates
        if ( M.iref == 1 || M.iref == 3 )
        {
            // Relative coordinates between -180 and 180
            liste[i].C.x -= M.ref_ra;
            if ( liste[i].C.x > 180. ) liste[i].C.x -= 360.;
            if ( liste[i].C.x < -180. ) liste[i].C.x += 360.;
            liste[i].C.x *= -3600 * cos(M.ref_dec * DTR);
            liste[i].C.y -= M.ref_dec;
            liste[i].C.y *= 3600;
        }
        else if ( M.iref == 2 )
        {
            liste[i].C.x -= M.ref_ra;
            liste[i].C.x -= M.ref_dec;
        }

        liste[i].E.theta *= DTR;
        if ( liste[i].E.a == 0. || liste[i].E.b == 0. )
            liste[i].c = 's';
        else
            liste[i].c = 'g';

        if ( liste[i].I0 == 0. )
            liste[i].I0 = 50.;

        i++;
    }


    fclose(IN);

    (*istart) = i;


    NPRINTF(stderr, "%ld arclets read\n", i);
}

/* Return the ra and dec values in degrees and scaling factors
 * parameters :
 * - line : the #REFERENCE line
 */
void getRADEC(char *line,
              int *iref, double *ra, double *dec )
{
    double  ss0, tt0;
    int     hh0, mm0, dd0, nn0;

    char    ref1[20], ref2[20];

    *iref = 0;
    *ra = *dec = 0.;

    if ( sscanf(line, "%*s%d%s%s",
                iref, ref1, ref2) == 3 )
    {
        // Convert ra and dec strings to double
        if ( *iref == 1 )
        {
            sscanf(ref1, "%d:%d:%lf", &hh0, &mm0, &ss0);
            sscanf(ref2, "%d:%d:%lf", &dd0, &nn0, &tt0);

            *ra = ((double)(hh0) + ((double)(mm0)) / 60. + ss0 / 3600.) * 15.;

            if (dd0 < 0)
                *dec = ((double)(dd0)) - ((double)(nn0)) / 60 - tt0 / 3600;
            else
                *dec = ((double)(dd0)) + ((double)(nn0)) / 60 + tt0 / 3600;
        }
        else if ( *iref == 3 || *iref == 2 )
        {
            sscanf(ref1, "%lf", ra);
            sscanf(ref2, "%lf", dec);
        }
    }
}

/* Convert a relative coordinates to absolute coordinates
 */
void convertXY( double *x, double *y, int iref, double ref_ra, double ref_dec )
{
    // Convert the input values to absolute WCS coordinates
    if ( iref == 1 || iref == 3 )
    {
        *x /= -3600.*cos(ref_dec * DTR);
        *x += ref_ra;
        *y /= 3600.;
        *y += ref_dec;
    }
    else if ( iref == 2 ) // image coordinates
    {
        *x += ref_ra;
        *y += ref_dec;
    }
}

