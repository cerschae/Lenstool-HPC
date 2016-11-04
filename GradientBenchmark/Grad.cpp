#include <Grad.h>



struct point module_potentialDerivatives_totalGradient(const int *Nlens, const struct point *pImage, PotentialSet *lens )
{
    struct point grad, clumpgrad;
	grad.x=0;
	grad.y=0;

	//This here could be done with function pointer to better acomodate future ass distributions functions
	// However I'm unsure of the time of function pointers -> ask gilles
	//for the moment lens and Nlens is organised the following way :  1. SIS, 2. PIEMD

	//SIS is the first
	for(int i=0; i<Nlens[0]; i++){
		clumpgrad=grad_halo_sis(pImage,i,&lens[0]);  //compute gradient for each clump separately

		if(clumpgrad.x == clumpgrad.x or clumpgrad.y == clumpgrad.y){ //nan check
		grad.x+=clumpgrad.x;
		grad.y+=clumpgrad.y;
		}  // add the gradients
	}

	//PIEMD is the second
	for(int i=0; i<Nlens[1]; i++){
		clumpgrad=grad_halo_piemd(pImage,i,&lens[1]);  //compute gradient for each clump separately

		if(clumpgrad.x == clumpgrad.x or clumpgrad.y == clumpgrad.y){ //nan check
		grad.x+=clumpgrad.x;
		grad.y+=clumpgrad.y;
		}  // add the gradients
	}

    return(grad);
}

/**@brief Return the gradient of the projected lens potential for one PIEMD clump. Uses SoA insteand of AoS lenses for speed
 *!!! You have to multiply by dlsds to obtain the true gradient for the expressions, see the papers :
 *JP Kneib & P Natarajan, Cluster Lenses, The Astronomy and Astrophysics Review (2011) for 1 and 2 and JP Kneib PhD (1993) for 3
*
* @param pImage 	point where the result is computed in the lens plane
* @param lens		mass distribution
*/

struct point grad_halo_piemd(const struct point *pImage, int iterator,PotentialSet *lens)
{
    struct point true_coord, true_coord_rotation, result;
    double R, angular_deviation;
    complex zis;

    result.x = result.y = 0.;

    /*positionning at the potential center*/
    true_coord.x = pImage->x - lens->x[iterator];  // Change the origin of the coordinate system to the center of the clump
    true_coord.y = pImage->y - lens->y[iterator];


	 /* PIEMD */
	/*rotation of the coordiante axes to match the potential axes*/
	true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[iterator]);
	/*Doing something....*/
	zis = piemd_1derivatives_ci05(true_coord_rotation.x, true_coord_rotation.y, lens->ellipticity_potential[iterator], lens->rcore[iterator]);

	result.x=lens->b0[iterator] * zis.re;
	result.y=lens->b0[iterator] * zis.im;

    return result;
}

/**@brief Return the gradient of the projected lens potential for one SIS clump. Uses SoA insteand of AoS lenses for speed
 *!!! You have to multiply by dlsds to obtain the true gradient for the expressions, see the papers :
 *JP Kneib & P Natarajan, Cluster Lenses, The Astronomy and Astrophysics Review (2011) for 1 and 2 and JP Kneib PhD (1993) for 3
*
* @param pImage 	point where the result is computed in the lens plane
* @param lens		mass distribution
*/

struct point grad_halo_sis(const struct point *pImage, int iterator,PotentialSet *lens)
{
    struct point true_coord, true_coord_rotation, result;
    double R, angular_deviation;
    complex zis;

    result.x = result.y = 0.;

    /*positionning at the potential center*/
    true_coord.x = pImage->x - lens->x[iterator];  // Change the origin of the coordinate system to the center of the clump
    true_coord.y = pImage->y - lens->y[iterator];

    /*Elliptical Isothermal Sphere*/
	/*rotation of the coordiante axes to match the potential axes*/
	true_coord_rotation = rotateCoordinateSystem(true_coord, lens->ellipticity_angle[iterator]);

	R=sqrt(true_coord_rotation.x*true_coord_rotation.x*(1-lens->ellipticity[iterator]/3.)+true_coord_rotation.y*true_coord_rotation.y*(1+lens->ellipticity[iterator]/3.));	//ellippot = ellipmass/3
	result.x=(1-lens->ellipticity[iterator]/3.)*lens->b0[iterator]*true_coord_rotation.x/(R);
	result.y=(1+lens->ellipticity[iterator]/3.)*lens->b0[iterator]*true_coord_rotation.y/(R);

    return result;
}


/**** usefull functions for PIEMD profile : see old lenstool ****/

/**  I*w,v=0.5 Kassiola & Kovner, 1993 PIEMD, paragraph 4.1
*
* Global variables used :
* - none
*/

complex piemd_1derivatives_ci05(double x, double y, double eps, double rc)
{
    double  sqe, cx1, cxro, cyro, rem2;
    complex zci, znum, zden, zis, zres;
    double norm;

    sqe = sqrt(eps);
    cx1 = (1. - eps) / (1. + eps);
    cxro = (1. + eps) * (1. + eps);
    cyro = (1. - eps) * (1. - eps);
    rem2 = x * x / cxro + y * y / cyro;
    /*zci=cpx(0.,-0.5*(1.-eps*eps)/sqe);
    znum=cpx(cx1*x,(2.*sqe*sqrt(rc*rc+rem2)-y/cx1));
    zden=cpx(x,(2.*rc*sqe-y));
    zis=pcpx(zci,lncpx(dcpx(znum,zden)));
    zres=pcpxflt(zis,b0);*/

    // --> optimized code
    zci.re = 0;
    zci.im = -0.5 * (1. - eps * eps) / sqe;
    znum.re = cx1 * x;
    znum.im = 2.*sqe * sqrt(rc * rc + rem2) - y / cx1;
    zden.re = x;
    zden.im = 2.*rc * sqe - y;
    norm = zden.re * zden.re + zden.im * zden.im;     // zis = znum/zden
    zis.re = (znum.re * zden.re + znum.im * zden.im) / norm;
    zis.im = (znum.im * zden.re - znum.re * zden.im) / norm;
    norm = zis.re;
    zis.re = log(sqrt(norm * norm + zis.im * zis.im));  // ln(zis) = ln(|zis|)+i.Arg(zis)
    zis.im = atan2(zis.im, norm);
//  norm = zis.re;
    zres.re = zci.re * zis.re - zci.im * zis.im;   // Re( zci*ln(zis) )
    zres.im = zci.im * zis.re + zis.im * zci.re;   // Im( zci*ln(zis) )
    //zres.re = zis.re*b0;
    //zres.im = zis.im*b0;

    return(zres);
}

/// Useful functions

// changes the coordinates of point P into a new basis (rotation of angle theta)
//	y'	y    x'
//	 *      |   /
//	   *    |  /  theta
//	     *  | /
//	       *|--------->x
struct point rotateCoordinateSystem(struct point P, double theta)
{
    struct  point   Q;

    Q.x = P.x*cos(theta) + P.y*sin(theta);
    Q.y = P.y*cos(theta) - P.x*sin(theta);

    return(Q);
}


/** @brief This module function calculates profile depended information like the impactparameter b0 and the potential ellipticity epot
 *
 * @param lens: mass distribution for which to calculate parameters
*/

void module_readParameters_calculatePotentialparameter(Potential *lens){

	switch (lens->type)
    {

	case(5): /*Elliptical Isothermal Sphere*/
		//impact parameter b0
		lens->b0 = 4* pi_c2 * lens->vdisp * lens->vdisp ;
		//ellipticity_potential
		lens->ellipticity_potential = lens->ellipticity/3 ;
	    break;

	case(8): /* PIEMD */
		//impact parameter b0
		lens->b0 = 6.*pi_c2 * lens->vdisp * lens->vdisp;
		//ellipticity_parameter
	    if ( lens->ellipticity == 0. && lens->ellipticity_potential != 0. ){
			// emass is (a2-b2)/(a2+b2)
			lens->ellipticity = 2.*lens->ellipticity_potential / (1. + lens->ellipticity_potential * lens->ellipticity_potential);
			//printf("1 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
		}
		else if ( lens->ellipticity == 0. && lens->ellipticity_potential == 0. ){
			lens->ellipticity_potential = 0.00001;
			//printf("2 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
		}
		else{
			// epot is (a-b)/(a+b)
			lens->ellipticity_potential = (1. - sqrt(1 - lens->ellipticity * lens->ellipticity)) / lens->ellipticity;
			//printf("3 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
		}
        break;

	default:
			std::cout << "ERROR: LENSPARA profil type of clump "<< lens->name << " unknown : "<< lens->type << std::endl;
            //printf( "ERROR: LENSPARA profil type of clump %s unknown : %d\n",lens->name, lens->type);
		break;
    };

}
