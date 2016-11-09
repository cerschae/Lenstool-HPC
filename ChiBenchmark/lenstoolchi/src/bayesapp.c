//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//            Bayesian Inference
//
// Filename:  bayesapp.c
//
// Purpose:   Test operation of BayeSys with
//            UserEmpty,UserTry1,UserTry2,UserInsert1,UserInsert2,UserDelete1.
//=============================================================================
//
// BayeSys3 can operate with arbitrary likelihood functions.
// In this example, an atom has 2 coordinates interpreted as
//            location in [0,1] = Cube[0]  (restricted to [0.0.75])
// and
//         flux in [0,infinity) = -10.0*log(Cube[1])
// to correspond with testing the operation of MassInf with FluxUnit0=10
//
// These are only example programs, involving inefficient calculations that
// repeatedly build an object from scratch with "UserBuild" instructions.
//
//=============================================================================
// History:   JS       2 Jan 2002, 4 Jan 2003, 10 Feb 2003, 20 Aug 2003
//            Copyright (c) 2002,2003 Maximum Entropy Data Consultants Ltd.
//-----------------------------------------------------------------------------
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "bayesys3.h"
#include "userstr.h"
#include "dimension.h"
#include "structure.h"
#include "fonction.h"
#include "constant.h"
#include "lt.h"

static int minusone; // count chi2 errors in image plane
static double nchi2;    // count number of computed chi2
static time_t start, end; // computation time
static int UserBuild(double*, CommonStr*, ObjectStr*, int, int);
#define POTENTIAL 0

struct matrix e_grad2_pot_ptr(const struct point *pi, const struct pot *ilens);

//=============================================================================
//  MANAGE THE THREAD SAFE ACCESS TO nchi2
//=============================================================================
#include <pthread.h>
static pthread_mutex_t mp = PTHREAD_MUTEX_INITIALIZER;

void init_mutex()
{
    int ret = pthread_mutex_init(&mp, NULL);
    if( ret != 0 )
    {
        fprintf(stderr, "ERROR initilializing mutex in bayesapp.c\n");
        exit(1);
    }
}

void destroy_mutex()
{
    int ret = pthread_mutex_destroy(&mp);
    if( ret != 0 )
    {
        fprintf(stderr, "ERROR destroying mutex in bayesapp.c\n");
        exit(1);
    }
	
}

void increase_nchi2()
{
    int ret = pthread_mutex_lock(&mp);
    nchi2 = nchi2 + 1;
    pthread_mutex_unlock(&mp);
}


//=============================================================================
//                           SIMPLE CODE
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserEmpty
//
// Purpose:   Set Lhood = logLikelihood(coords) for empty sample
//            with 0 atoms in Object, and initialise its other information.
//
//            Lhood := log L(empty)
//
//            I have already put the number 0 in Object->Natoms,
//            so you can use a simple "UserBuild" instruction.
//-----------------------------------------------------------------------------
//
int UserEmpty(        //   O  >=0, or -ve error code
			  double*    Lhood,     //   O  loglikelihood
			  CommonStr* Common,    // I O  general information
			  ObjectStr* Object)    // I O  sample object, output new Lhood
{
    UserBuild(Lhood, Common, Object, 0, 1);

    if( Common->Valency )
        *Lhood = 0.;  // Empty value is set afterwards by FluxEmpty with Mock[]=0

    return 1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserTry1
//
// Purpose:   d(logLikelihood(coords)) after supposedly adding one new atom
//            to Object, without any side effects on it.
//
//            If Valency = 0,
//              dLtry := d logL(...,x)
//            else mimic the above with                                cool
//              dLtry := (1/cool) log INTEGRAL dPrior(z) dlogL(...,x,z)
//
//            I have already put the new atom x after the last location
//            of the Atoms list, at Object->Cubes[ Object->Natoms ], so
//            if Valency=0 you can use simple "UserBuild" instructions.
//-----------------------------------------------------------------------------
//
// External definition
int rescaleParam(int id, int type, int ipx, double *pval);

int UserTry1(         //   O  +ve = OK, 0 = DO NOT USE, -ve = error
			 double*    dLtry,     //   O  trial d(logLikelihood) value
			 CommonStr* Common,    // I    general information
			 ObjectStr* Object)    // I    sample object (DO NOT UPDATE Lhood)
{
    int    Natoms = Object->Natoms;
    double Lhood  = Object->Lhood;                      // existing Lhood
    double Ltry;
    int    OK;
    extern struct g_grille  G;
	
    *dLtry = 0.0;
    OK = 1;
	
    OK = UserBuild(&Ltry, Common, Object, Natoms + 1, 0); // trial
    if ( OK > 0 )
        *dLtry = Ltry - Lhood;                          // increment
	
    return  OK;                                         // OK?
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserTry2
//
// Purpose:   d(logLikelihood(coords)) after supposedly adding one,
//            then two new atoms to Object, without any side effects on it.
//
//   If Valency = 0,
//     dLtry1 := d logL(...,x1)
//
//     dLtry2 := d logL(...,x1,x2)
//
//   else mimic the above with                                    cool
//     dLtry1 := (1/cool) log INTEGRAL dPrior(z1) dlogL(...,x1,z1)
//                                                                         cool
//     dLtry2 := (1/cool) log INTEGRAL dPrior(z1,z2) dlogL(...,x1,z1,x2,z2)
//
//            I have already put the new atoms x1,x2 after the last location
//            of the Atoms list, at Object->Cubes[ Object->Natoms ] and
//                               at Object->Cubes[ Object->Natoms + 1 ]
//            so if Valency=0 you can use simple "UserBuild" instructions.
//-----------------------------------------------------------------------------
//
int UserTry2(         //   O  +ve = OK, 0 = DO NOT USE, -ve = error
			 double*    dLtry1,    //   O  trial d(logLikelihood) value for 1st atom
			 double*    dLtry2,    //   O  trial d(logLikelihood) value for both atoms
			 CommonStr* Common,    // I    general information
			 ObjectStr* Object)    // I    sample object (DO NOT UPDATE Lhood)
{
    int    Natoms = Object->Natoms;
    double Lhood  = Object->Lhood;                           // existing Lhood
    double Ltry1, Ltry2;
    int    OK;
    extern struct g_grille  G;
	
    *dLtry1 = *dLtry2 = 0.0;
    OK = 1;
	
    OK = UserBuild(&Ltry1, Common, Object, Natoms + 1, 0);   //trial for 1 more
    if ( OK > 0 )
    {
        *dLtry1 = Ltry1 - Lhood;                             // increment for 1
        OK = UserBuild(&Ltry2, Common, Object, Natoms + 2, 0); //trial for 2 more
        if ( OK > 0 )
            *dLtry2 = Ltry2 - Lhood;                         // increment for 2
    }
	
    return  OK;                             // return OK only if both trials OK
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserInsert1
//
// Purpose:   Insert 1 new atom into Object, keeping it up-to-date, and
//            set d(loglikelihood(coords)).
//
//            If Valency = 0,
//              dL := d logL(...,x)
//
//            else
//                                                cool
//              sample z from  Prior(z) L(...,x,z)
//
//              and set    dL := d logL(...,x,z)  at the sampled z.
//
//            I have already put the new atom x at the last location of
//            the updated Atoms list, at Object->Cubes[ Object->Natoms - 1 ],
//-----------------------------------------------------------------------------
//
int UserInsert1(      //   O  >=0, or -ve error code
				double*    dL,        //   O  d(loglikelihood)
				CommonStr* Common,    // I    general information
				ObjectStr* Object)    // I O  sample object
{
    int    Natoms = Object->Natoms;              // new number
    double Lold   = Object->Lhood;               // not yet updated
    double Lnew;
    extern struct g_grille  G;
	
    UserBuild(&Lnew, Common, Object, Natoms, 1); // new updated state
    *dL = Lnew - Lold;                           // I will update Object->Lhood
	
    return 1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserInsert2
//
// Purpose:   Insert 2 new atoms into Object, keeping it up-to-date, and
//            set d(loglikelihood(fluxes)).
//
//            If Valency = 0,
//              dL := d logL(...,x1,x2)
//
//            else
//                                                                cool
//              sample z1,z2 from  Prior(z1,z2) L(...,x1,z1,x2,z2)
//
//              and set  dL := d logL(...,x1,z1,x2,z2)  at the sampled z1,z2.
//
//            I have already put the new atoms x1,x2 at the last location
//            of the Atoms list, at Object->Cubes[ Object->Natoms - 2 ] and
//                               at Object->Cubes[ Object->Natoms - 1 ]
//            so if Valency=0 you can use a simple "UserBuild" instruction.
//-----------------------------------------------------------------------------
//
int UserInsert2(      //   O  >=0, or -ve error code
				double*    dL,        //   O  d(loglikelihood)
				CommonStr* Common,    // I    general information
				ObjectStr* Object)    // I O  sample object
{
    int    Natoms = Object->Natoms;              // new number
    double Lold   = Object->Lhood;               // not yet updated
    double Lnew;
    extern struct g_grille  G;
	
    UserBuild(&Lnew, Common, Object, Natoms, 1); // new updated state
    *dL = Lnew - Lold;                           // I will update Object->Lhood
	
    return 1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserDelete1
//
// Purpose:   Delete 1 old atom from Object, keeping it up-to-date, and
//            set d(loglikelihood(fluxes)).
//
//            dL := d logL(...)
//
//            I have already put the old atom after the last location of
//            the updated Atoms list, at Object->Cubes[ Object->Natoms ],
//            so if Valency=0 you can use a simple "UserBuild" instruction.
//-----------------------------------------------------------------------------
//
int UserDelete1(      //   O  >=0, or -ve error code
				double*    dL,        //   O  d(loglikelihood)
				CommonStr* Common,    // I    general information
				ObjectStr* Object)    // I O  sample object
{
    int    Natoms = Object->Natoms;              // new number
    double Lold   = Object->Lhood;               // not yet updated
    double Lnew;
    extern struct g_grille  G;
	
    UserBuild(&Lnew, Common, Object, Natoms, 1); // new updated state
    *dL = Lnew - Lold;                           // I will update Object->Lhood
	
    return 1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserBuild
//
// Purpose:   Build Lhood, and any other data info such as Mock from scratch.
//            Can be called with any number of atoms
//
//            This example implementation has
//            Mock[0] = SUM[atoms] flux ,             flux = - 10*log(Cube[1])
//            Mock[1] = SUM[atoms] flux * location ,  location = Cube[0]
//-----------------------------------------------------------------------------
//
void  o_set_lens_ptr(struct pot *ilens, int ipx, double x);

static int UserBuild( //   O  +ve = OK, 0 = DO NOT USE, -ve = error
					 double*    Lhood,     //   O  loglikelihood
					 CommonStr* Common,    // I    General information
					 ObjectStr* Object,    // I(O) Cubes in, perhaps Mock out
					 int        Natoms,    // I    # atoms
					 int        output)    // I    +ve if Mock to be written out, -ve if not
{
	
    // MassInf optimization
    if( Common->Valency )
    {
        // Weak-Lensing is done in UserFoot() under the linear approximation
        *Lhood = Object->Lhood;
        return 1; 
		
        // Compute the likelihood for the SL 
        extern struct g_grille  G;
        double *np_b0 = (double *)calloc(G.nlens - G.nmsgrid, sizeof(double));
        long int ilens;
        int k, i, error;
        double **Cubes = Object->Cubes;
        double C, Lnorm;
        UserCommonStr *UserCommon = (UserCommonStr *) Common->UserCommon;
		
        struct point vB[NFMAX];
        extern struct g_image I;
        extern struct galaxie multi[NFMAX][NIMAX];
        int j;
        for( i = 0; i < I.n_mult; i++ )
        {
            vB[i].x = vB[i].y = 0.;
            for( j = 0; j < I.mult[i]; j++ )
            {
                vB[i].x += multi[i][j].C.x;
                vB[i].y += multi[i][j].C.y;
            }
        }
		
        for( k = 0; k < Natoms; k++ )
        {
         	   ilens = (int)((G.nlens - G.nmsgrid) * Cubes[k][Common->Ndim - 1]);
         	   np_b0[ilens] += Cubes[k][Common->Ndim];
               for( i = 0; i < I.n_mult; i++ )
                   for( j = 0; j < I.mult[i]; j++ )
                   {
                       vB[i].x -= Cubes[k][Common->Ndim] * multi[i][j].np_grad[ilens].x;
                       vB[i].y -= Cubes[k][Common->Ndim] * multi[i][j].np_grad[ilens].y;
                   }
        }

        double *Mock = Object->Mock;
        int nimage = 0;
        int nmult = Common->Ndata / 2;
        for( i = 0; i < I.n_mult; i++ )
            for( j = 0; j < I.mult[i]; j++ )
            {
                Mock[G.nlens - G.nmsgrid + nimage] = vB[i].x / I.mult[i];
                Mock[G.nlens - G.nmsgrid + nimage + nmult] = vB[i].y / I.mult[i];
                nimage++;
            }
        free(np_b0);
        increase_nchi2();
        return 1;

        // Compute Likelihood
        UserCommon->Nchi2++;
        C = 0.; Lnorm = 0.;
        error = o_chi_lhood0(&C, &Lnorm, np_b0);
        *Lhood = -0.5 * (C + Lnorm);
		
	    // Track NaN errors
	    if ( *Lhood != *Lhood )
	    {
    		fprintf(stderr, "ERROR: Lhood=NaN found in bayesapp.c:UserBuild(). Dump cube...\n");
    		for ( i = 0; i < Common->Ndim; i++ )
    			fprintf(stderr, "%lf ", np_b0[i]);
    		fprintf(stderr,"\n");
    		exit(2);
	    }
		
        free(np_b0);
        increase_nchi2();
        if ( error )
        {
            minusone++;
            return 0;
        }
        return 1;
    }  
	
	
    //----------------------------------------------------------------------
    // Case with WL only but with moving clumps (Common->Valency == 0)
    // Copy lens[G.nmsgrid], Natoms times
    //----------------------------------------------------------------------
    extern struct g_image I;
    extern struct g_grille G;
    double** Cubes  = Object->Cubes;  // I    Cubes in [0,1)  [Natoms][Ndim]
	
    if( I.stat == 8 && G.nmsgrid != G.nlens )
    {
        const double Sqrt2Pi = 2.50662827463100050240;
        extern struct pot lens[];
        extern int    block[][NPAMAX];
        extern long int narclet;
        extern struct galaxie arclet[];
        extern double*  tmpCube; 
        UserCommonStr* UserCommon = (UserCommonStr*)Common->UserCommon;
		
        int      Ndata  = Common->Ndata;  // I    # data
        double*  Data   = Common->Data;   // I    Data            [Ndata]
        double*  Acc    = Common->Acc;    // I    Accuracies      [Ndata]
        int      nDim   = Common->Ndim; 
        struct pot mylens = lens[0];      // rough copy, but sufficient
        struct matrix grad2;
		
        long int ilens, k;
        int ipx, i;
        int valid = 1;
        *Lhood = 0.;
		
        // Create Mock data
        double *Mock = (double *) calloc(Ndata, sizeof(double));
		
        for( i = 0; i < Natoms; i ++ )
        {
            int iparam = 0;
            memcpy(tmpCube, Cubes[i], (unsigned int)(Common->Ndim * sizeof(double)));
            for ( ipx = CX; ipx <= PMASS; ipx++ )
                if ( block[G.nmsgrid][ipx] != 0 )
                {
                    valid *= rescaleParam(G.nmsgrid, POTENTIAL, ipx, &tmpCube[iparam]);
                    o_set_lens_ptr(&mylens, ipx, tmpCube[iparam]);
                    iparam ++;
                }
			
            if( block[G.nmsgrid][B0] != 0 ) 
            {
                mylens.sigma = mylens.b0;
                mylens.b0 = 6.*pia_c2 * mylens.sigma * mylens.sigma;
            }
			
            for( k = 0; k < narclet; k++ )
            {
                grad2 = e_grad2_pot_ptr(&arclet[k].C, &mylens);
                Mock[k] += (grad2.a - grad2.c) * arclet[k].dr;
                Mock[narclet + k] += 2. * grad2.b * arclet[k].dr;
            }
        }
		
        // Accumulate logLikelihood from grid
        double   Lnorm  = 0.0;            // Likelihood normalisation
        double   C      = 0.0;            // chisquared
        for( k = 0; k < Ndata; k++ )
        {
            C += (Mock[k] - Data[k]) * Acc[k] * Acc[k] * (Mock[k] - Data[k]);
            Lnorm += log(Acc[k] / Sqrt2Pi);     // (fussy to include this)
        }
		
        *Lhood += Lnorm - C / 2.0;
        UserCommon->Nchi2++;
        if( output )
            for( k = 0; k < Ndata; k++ )
                Object->Mock[k] = Mock[k];
		
        free(Mock);
        return 1; 
    }
	
    //----------------------------------------------------------------------
    // Standard case where Common->Valency == 0  and 0 < Natoms < 1 (No MassInf optimization)
    //----------------------------------------------------------------------
    const double Sqrt2 = 1.414213562373095;
    extern double*  tmpCube; 
    int i, j;
    int error;
    UserCommonStr* UserCommon = (UserCommonStr*)Common->UserCommon;
    int valid = 0;
	
    *Lhood = 0.;
    if ( G.nmsgrid == G.nlens )
    {
        // Check for single atom 
        if ( Natoms != 1 )
            return valid;

        // Copy by value of the cube (cube is modified in rescaleCube_1Atom)
        for( j = 0; j < Common->Ndim; j++ )
	        tmpCube[j] = Cubes[0][j];
		
        /* Set the clumps parameters from the cube and eventually reject the whole object*/
        valid = rescaleCube_1Atom(tmpCube, Common->Ndim);
		
        /* Compute the likelihood for this object */
        if ( valid )
        {
            UserCommon->Nchi2++;
            *Lhood = o_lhood(&error);
			
            if ( ! error )
            {
        	    // Track NaN errors
        	    if ( *Lhood != *Lhood )
        	    {
            		fprintf(stderr, "ERROR: Lhood=NaN found in bayesapp.c:UserBuild(). Dump cube...\n");
            		for ( i = 0; i < Common->Ndim; i++ )
            			fprintf(stderr, "%lf ", tmpCube[i]);
            		fprintf(stderr,"\n");
            		exit(2);
        	    }
    			
                increase_nchi2();
            }
            else
                return 0;
			
#ifdef DEBUG
            FILE *dbg;
            dbg = fopen("debug.dat", "a");
            for( i = 0; i < Common->Ndim; i++ )
        	    fprintf(dbg, "%lf ",tmpCube[i]);
            fprintf(dbg, "\n");
            fclose(dbg);
#endif
			
        }
    }
    else
    {
        extern struct g_grille  G;
        double *np_b0 = (double *)calloc(G.nlens - G.nmsgrid, sizeof(double));

        for( i = 0; i < Natoms; i++ )
        {
            long int ilens = (int)((G.nlens - G.nmsgrid) * Cubes[i][0]);
            np_b0[ilens] += -10.0 * log(Cubes[i][1]);
        }
 
        // Compute Likelihood
        UserCommon->Nchi2++;
        double C = 0.; double Lnorm = 0.;
        error = o_chi_lhood0(&C, &Lnorm, np_b0);
        *Lhood = -0.5 * (C + Lnorm);

	    // Track NaN errors
	    if ( *Lhood != *Lhood )
	    {
    		fprintf(stderr, "ERROR: Lhood=NaN found in bayesapp.c:UserBuild(). Dump cube...\n");
    		for ( i = 0; i < Common->Ndim; i++ )
    			fprintf(stderr, "%lf ", np_b0[i]);
    		fprintf(stderr,"\n");
    		exit(2);
	    }

        free(np_b0);
        increase_nchi2();
        if ( error )
        {
            minusone++;
            return 0;
        }
        return 1;

    }

    return valid;
	
    return 1;
}

//=============================================================================
//       Dummy procedure to link to BayeSys3 when not running MassInf
//=============================================================================
int UserFoot(
			 double*    Cube,
			 CommonStr* Common,
			 int*       ibits,
			 double*    zbits,
			 int*       nbits)
{
    extern int *ibitstmp;
    extern double **zbitstmp1;
    extern struct g_grille G;
    extern struct g_image I;
	
    long int ilens;
    int nDim = Common->Ndim;
	
	ilens =  (int)((G.nlens - G.nmsgrid + I.n_mult) * Cube[nDim-1]);
	
    if( ilens < G.nlens - G.nmsgrid )
    {
        nbits[0] = Common->Ndata; nbits[1] = 0; nbits[2] = 0;
        memcpy(ibits, ibitstmp, (unsigned int)(Common->Ndata * sizeof(int)));;
        memcpy(zbits, zbitstmp1[ilens], (unsigned int)(Common->Ndata * sizeof(double)));
    }
    else
    {
        int j;
        int isrc = ilens - (G.nlens - G.nmsgrid);
	
        nbits[0] = 0; nbits[1] = I.mult[isrc]; nbits[2] = I.mult[isrc];

        int nimage = 0;
        for( j = 0; j < isrc; j ++ )
            nimage += I.mult[j];

        int nmult = nimage;
        for( j = isrc; j < I.n_mult; j++ )
            nmult += I.mult[j];
        
        for( j = 0; j < I.mult[isrc]; j++ )
        {
            ibits[j] = nimage;
            ibits[j+I.mult[isrc]] = nimage + nmult;
            zbits[j] = 1.; 
            zbits[j+I.mult[isrc]] = 1.;
            nimage++;
        }
    }

    increase_nchi2();
	
    return (Cube && Common && ibits && zbits && nbits) ? 1 : 1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function: getParName
//
// Purpose:  Return the name of the parameter corresponding to the ipx parameters
//=============================================================================
char* getParName(int ipx, char* name, int type)
{
    switch (ipx)
    {
        case(CX):
            strcpy( name, "x (arcsec)");
            break;
        case(CY):
            strcpy( name, "y (arcsec)");
            break;
        case(EPOT):
            strcpy( name, "epot");
            break;
        case(EMASS):
            if ( type == 14 )
                strcpy( name, "gamma" );
            else
                strcpy( name, "emass" );
            break;
        case(THETA):
            strcpy( name, "theta (deg)" );
            break;
        case(PHI):
            strcpy( name, "phi (deg)" );
            break;
        case(RC):
            if ( type == 13 )
                strcpy( name, "Re (arcsec)" );
            else if ( type == 12 || type == 16 )
                strcpy( name, "rs (arcsec)" );
            else
                strcpy( name, "rc (arcsec)" );
            break;
        case(B0):
            if ( type == 13 )
                strcpy( name, "Sigma_e (10^8 Msol/kpc^2)" );
            else
                strcpy( name, "sigma (km/s)" );
            break;
        case(ALPHA):
            if ( type == 13 )
                strcpy( name, "n" );
            else
                strcpy( name, "alpha" );
            break;
        case(BETA):
            if ( type == 12 )
                strcpy( name, "c" );
            else
                strcpy( name, "beta" );
            break;
        case(RCUT):
            if ( type == 12 )
                strcpy( name, "r200 (arcsec)" );
            else
                strcpy( name, "rcut (arcsec)" );
            break;
        case(MASSE):
            strcpy( name, "M200 (Msol)" );
            break;
        case(ZLENS):
            strcpy( name, "z_lens" );
            break;
        case(RCSLOPE):
            strcpy( name, "rc_slope" );
            break;
        case(PMASS):
            if ( type == 12 )
                strcpy( name, "rhos (Msol/Mpc^3)" );
            else
                strcpy( name, "pmass (g/cm2)" );
            break;
        case(OMEGAM):
            strcpy( name, "omegaM" );
            break;
        case(OMEGAX):
            strcpy( name, "omegaX" );
            break;
        case(WX):
            strcpy( name, "wX" );
            break;
        case(WA):
            strcpy( name, "wa" );
            break;
        case(SCX):
            strcpy( name, "x (arcsec)" );
            break;
        case(SCY):
            strcpy( name, "y (arcsec)" );
            break;
        case(SA):
            strcpy( name, "a (arcsec)" );
            break;
        case(SB):
            strcpy( name, "b (arcsec)" );
            break;
        case(SEPS):
            strcpy( name, "eps" );
            break;
        case(STHETA):
            strcpy( name, "theta (deg)" );
            break;
        case(SINDEX):
            strcpy( name, "index" );
            break;
        case(SFLUX):
            strcpy( name, "mag" );
            break;
        case(VFCX):
            strcpy( name, "velocity x (arcsec)");
            break;
        case(VFCY):
            strcpy( name, "velocity y (arcsec)");
            break;
        case(VFVT):
            strcpy( name, "vt (km/s)");
            break;
        case(VFRT):
            strcpy( name, "rt (kpc)");
            break;
        case(VFI):
            strcpy( name, "inclination (deg)");
            break;
        case(VFTHETA):
            strcpy( name, "velocity theta (deg)");
            break;
        case(VFLCENT):
            strcpy( name, "lambda (Angstroms)");
            break;
        case(VFSIGMA):
            strcpy( name, "velocity dispersion (km/s)");
            break;
    }
    return name;
}


/* Convert input Lenstool values to the corresponding values in <bayes.dat> file.
 * bayes2lt() function in readBayesModels.c
 */
static double lt2bayes(int i, int ipx, double val)
{
    extern struct pot       lens[];
    extern struct g_cosmo C;
	
    switch (ipx)
    {
        case(THETA):
            val *= RTD;
            break;
        case(B0):
            // copied from o_print_res()
            if ( lens[i].type <= 1 )
                val = sqrt(val / 4. / pia_c2);
            else if ( lens[i].type == 13 )
                val = val * 1e4 * (cH0_4piG * C.h / distcosmo1(lens[i].z)); // in 10^8 Msol/kpc^2
            else
                val = sqrt(val / 6. / pia_c2);
			
            break;
            /* Everything in solar masses
			 case(PMASS):
			 if( lens[i].type == 12 )
			 val /= 1e14; // rhos in 10^14 Msol -> Msol
			 break;
			 case(MASSE):
			 if( lens[i].type == 12 )
			 val /= 1e14; // masse in 10^14 Msol -> Msol
			 break;
			 */
        default:
            break;
    }
	
    return val;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function: getParVal
//
// Purpose:  Return the parameter value corresponding to the ilens and ipx parameters
//=============================================================================
static double getParVal(int i, int ipx)
{
    double x;
	
    x = o_get_lens(i, ipx);
	
    return lt2bayes(i, ipx, x);
}

static void addStat(double val, UserCommonStr *UserCommon, int param, long int samp)
{
    double var, var1, mean;
	
    UserCommon->err[param] += val * val;
    UserCommon->avg[param] += val;
    mean = UserCommon->avg[param] / samp;
    var1  = UserCommon->err[param] / samp - mean * mean;
	var = abs(var1);                           //TV
    UserCommon->sum += var / mean / mean;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserMonitor
//
// Purpose:   I have provided a new ensemble of objects.
//
//        1.  For probabilistic exploration, restrict Common->cool to <= 1.
//                    cool = 0 is prior, seen at the first call.
//                0 < cool < 1 is "burn-in"
//                    cool = 1 is evolve posterior,
//                    cool -> infinity is maximum likelihood
//            Return with a POSITIVE return code to exit BayeSys,
//            usually after the end of an evolution stage with cool = 1.
//            You can use Common->Nsystem, which counts calls to UserMonitor.
//
//        2.  You should reset any nuisance parameters x in UserCommon and/or
//            each of your UserObject structures, preferably by sampling from
//            their conditonal probabilities Pr(x|Objects) and/or Pr(x|Object).
//
//        3.  Collect your statistics and display your diagnostics.
//-----------------------------------------------------------------------------
//
int UserMonitor(            //   O  0 = continue, +ve = finish, -ve = abort
				CommonStr* Common,    // I    Full general information
				ObjectStr* Objects)   // I    Full ensemble of sample objects  [ENSEMBLE]
{
	
    // Track NaN errors
    if ( Common->cool != Common->cool )
    {
    	fprintf(stderr, "ERROR: Common->cool=NaN found in bayesapp.c:UserMonitor().\n");
    	exit(2);
    }
	
    // Variables used everywhere
    const double Sqrt2Pi = 2.50662827463100050240;
    extern struct g_mode    M;
    extern struct g_grille  G;
    UserCommonStr* UserCommon = (UserCommonStr*)Common->UserCommon;
    extern double *tmpCube;
    char  bayesName[30],burninName[30];   
    FILE  *bayes;
    int   k, ipx; // counters
    long int i, j;
    double lhood0, chi2;
	
    
    //.............................................................................
    // USER IMPOSES FINAL TEMPERATURE AND TERMINATION CONDITION !
    if ( Common->cool >= 1.0 && M.inverse == 3 )
        Common->cool = 1.0;
	
    //.............................................................................
    // Save the models in burnin.dat/bayes.dat
    strcpy(bayesName, "bayes");
    strcpy(burninName, "burnin");
#ifdef DEBUG
    sprintf(bayesName, "%s.dbg", bayesName);
    sprintf(burninName, "%s.dbg", burninName);
#endif
#ifdef MPI
    sprintf(bayesName, "%s.%d", bayesName, ltID);
    sprintf(burninName, "%s.%d", burninName, ltID);
#endif
    sprintf(bayesName, "%s.dat", bayesName);
    sprintf(burninName, "%s.dat", burninName);
	
    if ( Common->cool  < 1. )
        bayes = fopen( burninName, "a");
    else
        bayes = fopen( bayesName, "a");
	
    UserCommon->sum = 0;
    for ( k = 0; k < Common->ENSEMBLE; k++ )
    {
        double **Cubes = Objects[k].Cubes;
        extern double *np_b0;
		
        if( Common->Valency )
        {
            extern struct g_image I;
            extern double **zbitstmp1;
            int Ndata = Common->Ndata; 
            double *Data = Common->Data;
            double *Acc = Common->Acc;
            long int ilens;
			
            // Initialize arrays
            double *Mock = (double *) calloc(Ndata, sizeof(double));
         	memset(np_b0, 0, (unsigned long int)((G.nlens - G.nmsgrid) * sizeof(double)));
			
            // MassInf optimization (commented: 1 optimized pot.) 
      	    for( i = 0; i < Objects[k].Natoms; i++ )
            {
                //  Not implemented in UserBuild()
				//                for( j = 0; j < Common->Ndim - 2; j ++ )
				//                    tmpCube[j] += Cubes[i][j];  // TODO: rescale between 0..1 with Unif prior
				
                // Cubes[][Common->Ndim - 1] is coordinate
         	   ilens = (int)((G.nlens - G.nmsgrid + I.n_mult) * Cubes[i][Common->Ndim - 1]);
               if( ilens < G.nlens - G.nmsgrid )
               {
             	   np_b0[ilens] += Cubes[i][Common->Ndim];

                   for( j = 0; j < Ndata; j++ )
                       Mock[j] += Cubes[i][Common->Ndim] * zbitstmp1[ilens][j];
               }
               else
               {
                   int isrc = ilens - (G.nlens - G.nmsgrid);
                   int nimage = 0;
                   for( j = 0; j < isrc; j ++ )
                       nimage += I.mult[j];

                   int nmult = nimage;
                   for( j = isrc; j < I.n_mult; j++ )
                       nmult += I.mult[j];
                    
                   for( j = 0; j < I.mult[isrc]; j++ )
                   {
                       Mock[nimage] += Cubes[i][Common->Ndim + 1];
                       Mock[nimage + nmult] += Cubes[i][Common->Ndim + 2];
                       nimage++;
                   }
               }
            }
			
            // Compute Likelihood from SL
            chi2 = 0.; lhood0 = 0.;
			
			// Accumulate logLikelihood from grid
            for( j = 0; j < Ndata; j++ )
            {
                chi2 += (Mock[j] - Data[j]) * Acc[j] * Acc[j] * (Mock[j] - Data[j]);
                lhood0 += log(Sqrt2Pi / Acc[j]); 
            }  
 
            // Compute Likelihood (not to include otherwise there is no match 
            // when reading with readBayesModel.c)
            //int error = o_chi_lhood0(&chi2, &lhood0, np_b0);  // TODO: catch error
            free(Mock);
        }
        else if( G.nmsgrid != G.nlens )
        {
            // CONDITION NOT USED BECAUSE EVERYTHING IS DONE WITH MASSINF NOW
            long int ilens;

         	memset(np_b0, 0, (unsigned long int)((G.nlens - G.nmsgrid) * sizeof(double)));
      	    for( i = 0; i < Objects[k].Natoms; i++ )
            {
         	   ilens = (int)((G.nlens - G.nmsgrid) * Cubes[i][0]);
               np_b0[ilens] += -10.0 * log(Cubes[i][1]);
            }

            int error = o_chi_lhood0(&chi2, &lhood0, np_b0);  // TODO: catch error
        }
        else
        {
            // Reset the tmpCube array 
         	memset(tmpCube, 0, (unsigned long int)(UserCommon->Ndim * sizeof(double)));
			
            int iparam = 0;
			
            // Natoms should be 1, but still allow case with Natoms > 1  
      	    for( i = 0; i < Objects[k].Natoms; i++ )
				for( j = 0; j < Common->Ndim; j ++ )
					tmpCube[iparam++] = Cubes[i][j];
			
            // Normal case: UserCommon->Ndim is the number of free parameters
            // But eg, UserCommon->Ndim can be Natoms times the same pot
            rescaleCube_1Atom(tmpCube, UserCommon->Ndim);
			
            // Compute Likelihood
            int error = o_chi_lhood0(&chi2, &lhood0, NULL);  // TODO: catch error

            // Loop to the next iteration if error
            if( error )
                continue;
        }
		
        double o_lhood_rez = (-0.5*(chi2 + lhood0));
        fprintf( bayes, "%ld %lf ", UserCommon->Nsample + 1, o_lhood_rez);
		
        // SAVE: save lens parameters
        int iparam = 0;
        long int samp = UserCommon->Nsample * Common->ENSEMBLE + k + 1;
        extern int    block[][NPAMAX];
        double val = 0.;
		
        for ( i = 0; i < G.no_lens; i++ )
            for ( ipx = CX; ipx <= PMASS; ipx++ )
                if ( block[i][ipx] != 0 )
                {
                    val = getParVal(i, ipx);
                    fprintf(bayes, "%lf ", val );
                    addStat(val, UserCommon, iparam++, samp);
                }
		
        int nzeros= 0;
        for ( i = 0; i < G.nlens - G.nmsgrid; i++ )
        {
            val = np_b0[i];
			
            if ( fabs(val) > 1E-6 )
            {
                if( nzeros > 0 )  
                    fprintf(bayes, "%dx0.0 %lf ", nzeros, val);
                else
                    fprintf(bayes, "%lf ", val);
                nzeros = 0;
            }
            else
                nzeros++;
			
            addStat(val, UserCommon, iparam++, samp);
        }
        if ( nzeros > 0 )
            fprintf(bayes, "%dx0.0 ", nzeros);
		
        // SAVE: source positions for grid model (for DEBUG purposes for now)
        extern struct g_image I;
        if ( Common->Valency && I.n_mult != 0)
        {
            double *Ps = (double *)calloc(I.n_mult * 2, sizeof(double));
            long int ilens;
            for ( i = 0; i < Objects[k].Natoms; i++ )
            {
                ilens = (int)((G.nlens - G.nmsgrid + I.n_mult) * Cubes[i][Common->Ndim - 1]);
                if( ilens >= G.nlens - G.nmsgrid )
                {
                    int isrc = ilens - (G.nlens - G.nmsgrid);
                    Ps[isrc] += Cubes[i][Common->Ndim + 1];
                    Ps[isrc + I.n_mult] += Cubes[i][Common->Ndim + 2];
                }
            }

            for ( i = 0; i < I.n_mult * 2; i++ )
                fprintf(bayes, "%lf ", Ps[i]);
            free(Ps);
        }
			
        // SAVE: source parameters for M.iclean == 2
        if ( M.iclean == 2 )
        {
            extern int sblock[NFMAX][NPAMAX];
            extern struct g_source   S;
            extern struct galaxie    source[NFMAX];
			
            for( i = 0; i < S.ns; i++ )
                for( ipx =SCX; ipx <= SFLUX; ipx++ )
                    if( sblock[i][ipx] != 0 )
                    {
                        val = o_get_source(&source[i], ipx);
                        if( ipx == STHETA )    val *= RTD;
                        fprintf(bayes, "%f ", val);
                        addStat(val, UserCommon, iparam++, samp);
                    }
        }
		
        // SAVE: save cosmological parameters
        extern int cblock[NPAMAX];
        for ( ipx = OMEGAM; ipx <= WA; ipx++ )
            if ( cblock[ipx] != 0 )
            {
                val = getParVal(0, ipx);
                fprintf(bayes, "%lf ", val);
                addStat(val, UserCommon, iparam++, samp);
            }
		
        // SAVE: save zmlimit parameters
        extern struct g_image   I;
        extern struct z_lim     zlim[];
        extern struct galaxie   multi[NFMAX][NIMAX];
        char limages[ZMBOUND][IDSIZE];
        int nimages;
		
        for ( ipx = 0; ipx < I.nzlim; ipx++ )
            if ( zlim[ipx].bk != 0 )
            {
                // look for the images families corresponding to the zmlimit to optimize
                nimages = splitzmlimit(zlim[ipx].n, limages);
                i = 0;
                while ( indexCmp( multi[i][0].n, limages[0] ) ) i++;
                val = multi[i][0].z;
                fprintf( bayes, "%lf ", val );
                addStat(val, UserCommon, iparam++, samp);
                /*              i = 0;
				 while( indexCmp( multi[i][0].n, zlim[ipx].n ) ) i++;
				 for( k = 0; k < I.mult[i]; k++ )
				 fprintf( bayes, "%lf ", multi[i][k].z );
				 */
            }

        // SAVE: save z_a_limit parameter
        extern struct z_lim     zalim;
        if( zalim.bk != 0 )
        {
            val = I.zarclet;
            fprintf( bayes, "%lf ", val);
            addStat(val, UserCommon, iparam++, samp);
        }
		
        // SAVE: save velocity field parameters
        extern int vfblock[NPAMAX];
        for ( ipx = VFCX; ipx <= VFSIGMA; ipx++ )
            if ( vfblock[ipx] != 0 )
            {
                val = getParVal(0, ipx);
                if( ipx == VFTHETA )    val *= RTD;
                if( ipx == VFI )    val *= RTD;
                fprintf(bayes, "%lf ", val);
                addStat(val, UserCommon, iparam++, samp);
            }
		
        // SAVE : save pot parameters
        extern struct g_pot     P[NPOTFILE];
        for( i = 0; i < G.npot; i++ )
        {
            struct g_pot  *pot = &P[i];
			
            if ( pot->ftype != 0 )
            {
                if ( pot->ircut != 0 )
                {
                    fprintf( bayes, "%lf ", pot->cut );
                    addStat(pot->cut, UserCommon, iparam++, samp);
                }
                if ( pot->isigma != 0 )
                {
                    fprintf( bayes, "%lf ", pot->sigma );
                    addStat(pot->sigma, UserCommon, iparam++, samp);
                }
                if ( pot->islope != 0 )
                {
                    fprintf( bayes, "%lf ", pot->slope );
                    addStat(pot->slope, UserCommon, iparam++, samp);
                }
                if ( pot->ivdslope != 0 )
                {
                    fprintf( bayes, "%lf ", pot->vdslope );
                    addStat(pot->vdslope, UserCommon, iparam++, samp);
                }
                if ( pot->ivdscat != 0 )
                {
                    fprintf( bayes, "%lf ", pot->vdscat );
                    addStat(pot->vdscat, UserCommon, iparam++, samp);
                }
                if ( pot->ircutscat != 0 )
                {
                    fprintf( bayes, "%lf ", pot->rcutscat );
                    addStat(pot->rcutscat, UserCommon, iparam++, samp);
                }
                if ( pot->ia != 0 )
                {
                    fprintf( bayes, "%lf ", pot->a );
                    addStat(pot->rcutscat, UserCommon, iparam++, samp);
                }
                if ( pot->ib != 0 )
                {
                    fprintf( bayes, "%lf ", pot->b );
                    addStat(pot->rcutscat, UserCommon, iparam++, samp);
                }
            }
        }
		
        // SAVE: the nuisance parameters
        extern struct sigposStr sigposAs;
        if ( sigposAs.bk != 0 )
        {
            for ( i = 0; i < I.n_mult; i++ )
                for ( j = 0; j < I.mult[i]; j++ )
                {
                    val = sqrt(I.sig2pos[i][j]);
                    fprintf( bayes, "%lf ", val );
                    addStat(val, UserCommon, iparam++, samp);
                }
        }
		
        if ( I.dsigell != -1. )
        {
            val = sqrt(I.sig2ell);
            fprintf( bayes, "%lf ", val );
            addStat(val, UserCommon, iparam++, samp);
        }
		
        // append the evidence
        if ( Common->cool  < 1. ) fprintf( bayes, "%lf ", Common->Evidence);
        
        //append #Chi2 in last column in bayes.dat
        if ( Common->cool  >= 1. ) fprintf( bayes, "%lf ", chi2); 
		
        fprintf( bayes, "\n" );
    }
	
    fclose(bayes);
	
	//.............................................................................
	// Print the status line of STDOUT
    int Nrem;
    if( Common->cool < 1. )
    {
        double chi2rate = difftime(end, start);
        chi2rate = chi2rate > 0 ? chi2rate : -1;
        chi2rate = nchi2 / chi2rate;
		
        Nrem = Common->Nsystem;    // # to go, as (2 * Nsample++) catches up with Nsystem++
		
        printf("                                                                                \r");
        fprintf(stdout,                     // (stderr flushes output immediately)
                "Burn-in : %10.6lf %5d %12.3lf %10.3lf   %d/%.0lf  %.0lfchi2/s\r",
                Common->cool, Nrem, chi2, Common->Evidence, minusone, nchi2, chi2rate);
    }
    else
    {
        if ( M.itmax == 0 )
            M.itmax = Common->Nsystem;
		
        Nrem = M.itmax - UserCommon->Nsample;  // or any other setting
        UserCommon->sum /= Common->ENSEMBLE * Common->Ndim;
        fprintf(stdout, "Sampling : %10.6lf %7ld/%d %12.3lf %10.3lf[CTRL-C to interrupt]\r",
                Common->cool, UserCommon->Nsample + 1, M.itmax + 1, chi2, sqrt(UserCommon->sum));
    }
    fflush(stdout);
    minusone = 0;
    nchi2 = 0.;
    start = end;
    time(&end);
	
	//.............................................................................
	// Set any NUISANCE PARAMETERS you may have here
	// Could over-ride system's choices of FluxUnit (both Common and Objects) here
	// Could CALIBRATE data here (corrupting Evidence and Information)
	//.............................................................................
	// BURN-IN
    extern int optInterrupt; // Global variable modified in o_run_bayes.c/signalReset()
    if ( Common->cool < 1.0 && !optInterrupt) // final temperature not reached..
        return 0;             // ...normal return, keep going
	//.............................................................................
	// EVOLUTION: Accumulate any desired statistics
    for ( k = 0; k < Common->ENSEMBLE; k++ )
        UserCommon->atoms += Objects[k].Natoms; // # atoms
	
    // # samples of full ensemble
    UserCommon->Nsample ++;
	
    if ( Nrem > 0 && !optInterrupt)                                // not finished...
        return 0;                                  // ..normal return
	//.............................................................................
	// EXIT WITH POSITIVE FINISH CODE WHEN EVOLUTION IS DEEMED COMPLETE
    NPRINTF(stderr, "\n" );
    return 1;
}
