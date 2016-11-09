#ifndef GRILLE_H
#define GRILLE_H


#include "wcs.h"
/*
*  grille.h
*  kneib jean-paul
*  septembre 94
*  IoA cambridge
*  lens tool V2
*/

/*
*  useful header
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
*  useful function definition
*/
#define	Min(A,B)	((A)<(B)?(A):(B))
#define	Max(A,B)	((A)>(B)?(A):(B))
#define FPRINTF         if (M.verbose > 0) fprintf
#define NPRINTF         if (M.verbose == 1) fprintf

/*
*  dimension header
*/
#include "dimension.h"

/*
*  constants header
*/
#include "constant.h"

/*
*  structures header
*/
#include "structure.h"


#endif // if GRILLE_H
