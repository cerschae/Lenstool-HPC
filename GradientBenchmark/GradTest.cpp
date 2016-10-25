#include <GradTest.h>
#include <structure.h>


complex piemd_1derivatives(double x, double y, double eps, double rc)
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

void gradtest()
{
	/** Testing for Gradient function **/

	//////Setting the problem
	int big(1);
	runmode_param runmode;
	runmode.nhalos = big;
	point image;
	image.x = image.y = 2;

	///////theoretical gradient calculation
	//Init
	point gradtheo;
	gradtheo.x = gradtheo.y = 0;
	double b0(0), vdisp(1.), ellipticity_potential(0), ellipticity(0.11);
	complex zis;
	//Computation
	ellipticity_potential = (1. - sqrt(1 - ellipticity * ellipticity)) / ellipticity;
	b0 = 6.*pi_c2 * vdisp * vdisp;


	//Doing something....
    zis = piemd_1derivatives(image.x, image.y, ellipticity_potential, 1);

    gradtheo.x=b0* zis.re;
    gradtheo.y=b0 * zis.im;



	///////Gradient calculation using functions

    //Init
	point grad;
	Potential  *ilens;
	Potential lens[big];

	grad.x = grad.y = 0;
	for (int i = 0; i <big; ++i){
		ilens = &lens[i];

	    ilens->position.x = ilens->position.y = 0.;
	    ilens->type = 8;
	    ilens->ellipticity = 0.11;
	    ilens->ellipticity_potential = 0.;
	    ilens->ellipticity_angle = 0.;
	    ilens->vdisp = 1.;
	    ilens->rcut = 5.;
	    ilens->rcore = 1;
	    ilens->weight = 0;
	    ilens->rscale = 0;
	    ilens->exponent = 0;
	    ilens->alpha = 0.;
	    ilens->einasto_kappacritic = 0;
	    ilens->z = 0.4;
	    module_readParameters_calculatePotentialparameter(ilens);

	}
	//Init SoA
	PotentialSet lenses;
	lenses.type = 	new int[big];
	lenses.x  = 	new double[big];
	lenses.y = 		new double[big];
	lenses.b0 = 	new double[big];
	lenses.ellipticity_angle = new double[big];
	lenses.ellipticity = new double[big];
	lenses.ellipticity_potential = new double[big];
	lenses.rcore = 	new double[big];
	lenses.rcut = 	new double[big];
	lenses.z = 		new double[big];

	for (int i = 0; i <big; ++i){
		lenses.type[i] = 	lens[i].type;
		lenses.x[i]  = 		lens[i].position.x;
		lenses.y[i] = 		lens[i].position.y;
		lenses.b0[i] = 		lens[i].b0;
		lenses.ellipticity_angle[i] = lens[i].ellipticity_angle;
		lenses.ellipticity[i] = lens[i].ellipticity;
		lenses.ellipticity_potential[i] = lens[i].ellipticity_potential;
		lenses.rcore[i] = 	lens[i].rcore;
		lenses.rcut[i] = 	lens[i].rcut;
		lenses.z[i] = 		lens[i].z;

	}
	// Computation
	grad = module_potentialDerivatives_totalGradient(&runmode,&image, &lenses);


	////////Checking

	std::cout << "Theoretically we should have "<< gradtheo.x << "  and  " << gradtheo.y << std::endl;
	std::cout << "and the function returns "<< grad.x << "  and  " << grad.y << std::endl;

}

