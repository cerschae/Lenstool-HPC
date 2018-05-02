#include <iostream>
#include <cstring>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "convert_coordinates.hpp"
#include <iomanip>

/* Convert relative coordinate to absolute coordinates
 */
void convertXY_to_abs( double *x, double *y, int iref, double ref_ra, double ref_dec )
{
	double DTR=acos(-1.)/180.;
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

/* Convert absolute coordinate to relative coordinates
 */
void convertXY_to_rela( double *x, double *y, int iref, double ref_ra, double ref_dec )
{
	double DTR=acos(-1.)/180.;
	// convert to output relative coordinates
	if (  iref == 1 || iref == 3 )
	{
		*x -= ref_ra;
		*x *= -3600 * cos(ref_dec * DTR);
		*y -= ref_dec;
		*y *= 3600;
	}
	else if ( iref == 2 )
	{
		*x -= ref_ra;
		*y -= ref_dec;
	}
}
