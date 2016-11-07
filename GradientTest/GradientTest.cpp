#include <iostream>
#include <string.h>
#include "structure.h"
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <../GradientBenchmark/gradient.hpp>

#define EPS 	0.000000000000001

complex piemd_1derivatives(double x, double y, double eps, double rc);

void rotatecoordinatestest();
void piemdtest();
void sistest();

int main()
{

	//Rotation of coordinates
	rotatecoordinatestest();
	//Tests for PIEMD potential
	piemdtest();
	//Tests for SIS potential
	sistest();
}

void rotatecoordinatestest(){
	//Init
    point P;
    P.x = P.y = 2;
    double theta = 0.3; //rad

    point coord = rotateCoordinateSystem(P,theta);
    if (coord.x != P.x*cos(theta) + P.y*sin(theta)){
		std::cerr << "Rotatecoordinatestest X: Theoretically we should have "
				<< P.x*cos(theta) + P.y*sin(theta) << "  and  "
				<< coord.x << std::endl;
		exit(0);
    }
    if (coord.y != P.x*cos(theta) - P.y*sin(theta)){
		std::cerr << "Rotatecoordinatestest Y: Theoretically we should have "
				<<  P.x*cos(theta) - P.y*sin(theta) << "  and  "
				<< coord.y << std::endl;
		exit(0);
    }
}

void piemdtest() {
	//////Setting the problem
	int big(1);
	int N[2];
	N[0] = 0;
	N[1] = 1;
	runmode_param runmode;
	point image;
	runmode.nhalos = big;
	image.x = image.y = 2;
	/////// Test for 1 PIEMD mass distribution:
	//1: module_readParameters_calculatePotentialparameter
	//2: piemd_1derivatives_ci05
	//3: gradient
	//1: module_readParameters_calculatePotentialparameter
	///////theoretical b0, ellipticity calculation
	//Init
	double b0(0), ellipticity_potential(0);
	double vdisp(1.), ellipticity(0.11);
	//Computation
	ellipticity_potential = (1. - sqrt(1 - ellipticity * ellipticity))
			/ ellipticity;
	b0 = 6. * pi_c2 * vdisp * vdisp;
	///////numerical b0, ellipticity_potential calculation
	//Init
	Potential* ilens;
	Potential lens[big];
	for (int i = 0; i < big; ++i) {
		ilens = &lens[i];
		ilens->position.x = ilens->position.y = -1.;
		ilens->type = 8;
		ilens->ellipticity = 0.11;
		ilens->ellipticity_potential = 0.;
		ilens->ellipticity_angle = 0.4;
		ilens->vdisp = 1.;
		ilens->rcut = 5.;
		ilens->rcore = 1;
		ilens->weight = 0;
		ilens->rscale = 0;
		ilens->exponent = 0;
		ilens->alpha = 0.;
		ilens->einasto_kappacritic = 0;
		ilens->z = 0.4;
		//Computation
		module_readParameters_calculatePotentialparameter(ilens);
	}
	if (b0 != lens[0].b0) {
		std::cerr << "PIEMD b0: Theoretically we should have " << b0 << "  and  "
				<< lens[0].b0 << std::endl;
		exit(0);
	}
	if (ellipticity_potential != lens[0].ellipticity_potential) {
		std::cerr << "PIEMD ellipticity_potential: Theoretically we should have "
				<< ellipticity_potential << "  and  "
				<< lens[0].ellipticity_potential << std::endl;
		exit(0);
	}
	//2: piemd_1derivatives_ci05
	//Init
	complex ztheo, zcomput;
	///////theoretical piemd
	ztheo = piemd_1derivatives(image.x, image.y, ellipticity_potential, 1);
	///////numerical piemd
	zcomput = piemd_1derivatives_ci05(image.x, image.y, ellipticity_potential,1);
	if (ztheo.re != zcomput.re) {
		std::cerr << "PIEMD z_piemd.re: Theoretically we should have " << ztheo.re
				<< "  and  " << zcomput.re << std::endl;
		exit(0);
	}
	if (ztheo.im != zcomput.im) {
		std::cerr << "PIEMD z_piemd.im: Theoretically we should have " << ztheo.im
				<< "  and  " << zcomput.im << std::endl;
		exit(0);
	}
	//3: gradient
	///////theoretical
	//Init
	point gradtheo, grad, true_coord, true_coord_rot;
	grad.x = grad.y = 0;
	//Computation
    	true_coord.x = image.x + 1;  // Change the origin of the coordinate system to the center of the clump
    	true_coord.y = image.y + 1;
    	true_coord_rot = rotateCoordinateSystem(true_coord,0.4);
	ztheo = piemd_1derivatives(true_coord_rot.x, true_coord_rot.y, ellipticity_potential, 1);
	gradtheo.x = b0 * ztheo.re;
	gradtheo.y = b0 * ztheo.im;
	///////Computation
	//Init SoA
	Potential_SOA L[2];
	Potential_SOA lenses;
	lenses.type = new int[big];
	lenses.position_x = new double[big];
	lenses.position_y = new double[big];
	lenses.b0 = new double[big];
	lenses.ellipticity_angle = new double[big];
	lenses.ellipticity = new double[big];
	lenses.ellipticity_potential = new double[big];
	lenses.rcore = new double[big];
	lenses.rcut = new double[big];
	lenses.z = new double[big];
	for (int i = 0; i < big; ++i) {
		lenses.type[i] = lens[i].type;
		lenses.position_x[i] = lens[i].position.x;
		lenses.position_y[i] = lens[i].position.y;
		lenses.b0[i] = lens[i].b0;
		lenses.ellipticity_angle[i] = lens[i].ellipticity_angle;
		lenses.ellipticity[i] = lens[i].ellipticity;
		lenses.ellipticity_potential[i] = lens[i].ellipticity_potential;
		lenses.rcore[i] = lens[i].rcore;
		lenses.rcut[i] = lens[i].rcut;
		lenses.z[i] = lens[i].z;
	}
	// Computation
	L[1]=lenses;
	grad = module_potentialDerivatives_totalGradient_SOA(N, &image, L);
	//std::cerr << big << std::endl;
	if (fabs(gradtheo.x - grad.x ) > EPS) {
		std::cerr << "PIEMD grad.x: Theoretically we should have " << gradtheo.x
				<< "  and  " << grad.x << std::endl;
		exit(0);
	}
	if (fabs(gradtheo.y - grad.y) > EPS) {
		std::cerr << "PIEMD grad.y: Theoretically we should have " << gradtheo.y
				<< "  and  " << grad.y << std::endl;
		exit(0);
	}
	
	grad = module_potentialDerivatives_totalGradient_SOA_AVX(N, &image, L);
	//std::cerr << big << std::endl;
	if (fabs(gradtheo.x - grad.x ) > EPS) {
		std::cerr << "PIEMD AVX grad.x: Theoretically we should have " << gradtheo.x
				<< "  and  " << grad.x << std::endl;
		exit(0);
	}
	if (fabs(gradtheo.y - grad.y) > EPS) {
		std::cerr << "PIEMD AVX grad.y: Theoretically we should have " << gradtheo.y
				<< "  and  " << grad.y << std::endl;
		exit(0);
	}
	std::cout << "PIEMD Everythings fine, No test failed " << std::endl;
}

void sistest() {
	//////Setting the problem
	int big(1);
	int N[2];
	N[0] = 1;
	N[1] = 0;
	runmode_param runmode;
	point image;
	runmode.nhalos = big;
	image.x = image.y = 2;
	/////// Test for 1 SIS mass distribution:
	//1: module_readParameters_calculatePotentialparameter
	//2: gradient
	//1: module_readParameters_calculatePotentialparameter
	///////theoretical b0, ellipticity calculation
	//Init
	double b0(0), ellipticity_potential(0);
	double vdisp(1.), ellipticity(0.11);

	//Computation
	ellipticity_potential = ellipticity/3;
	b0 = 4. * pi_c2 * vdisp * vdisp;

	///////numerical b0, ellipticity_potential calculation
	//Init
	Potential* ilens;
	Potential lens[big];
	for (int i = 0; i < big; ++i) {
		ilens = &lens[i];
		ilens->position.x = ilens->position.y = -1.;
		ilens->type = 5;
		ilens->ellipticity = 0.11;
		ilens->ellipticity_potential = 0.;
		ilens->ellipticity_angle = 0.4;
		ilens->vdisp = 1.;
		ilens->rcut = 5.;
		ilens->rcore = 1;
		ilens->weight = 0;
		ilens->rscale = 0;
		ilens->exponent = 0;
		ilens->alpha = 0.;
		ilens->einasto_kappacritic = 0;
		ilens->z = 0.4;
		//Computation
		module_readParameters_calculatePotentialparameter(ilens);
	}
	if (b0 != lens[0].b0) {
		std::cerr << "SIS b0: Theoretically we should have " << b0 << "  and  "
				<< lens[0].b0 << std::endl;
		exit(0);
	}
	if (ellipticity_potential != lens[0].ellipticity_potential) {
		std::cerr << "SIS ellipticity_potential: Theoretically we should have "
				<< ellipticity_potential << "  and  "
				<< lens[0].ellipticity_potential << std::endl;
		exit(0);
	}

	//2: gradient
	///////theoretical
	//Init
	double R;
	point gradtheo, grad, true_coord, true_coord_rot;
	grad.x = grad.y = 0;
	//Computation
    true_coord.x = image.x + 1;  // Change the origin of the coordinate system to the center of the clump
    true_coord.y = image.y + 1;
    true_coord_rot = rotateCoordinateSystem(true_coord,0.4);
	R=sqrt(true_coord_rot.x*true_coord_rot.x*(1-ellipticity/3.)+true_coord_rot.y*true_coord_rot.y*(1+ellipticity/3.));
	gradtheo.x=(1-ellipticity/3.)*b0*true_coord_rot.x/(R);
	gradtheo.y=(1+ellipticity/3.)*b0*true_coord_rot.y/(R);

	///////Computation
	//Init SoA
	Potential_SOA L[2];
	Potential_SOA lenses;
	lenses.type = new int[big];
	lenses.position_x = new double[big];
	lenses.position_y = new double[big];
	lenses.b0 = new double[big];
	lenses.ellipticity_angle = new double[big];
	lenses.ellipticity = new double[big];
	lenses.ellipticity_potential = new double[big];
	lenses.rcore = new double[big];
	lenses.rcut = new double[big];
	lenses.z = new double[big];
	for (int i = 0; i < big; ++i) {
		lenses.type[i] = lens[i].type;
		lenses.position_x[i] = lens[i].position.x;
		lenses.position_y[i] = lens[i].position.y;
		lenses.b0[i] = lens[i].b0;
		lenses.ellipticity_angle[i] = lens[i].ellipticity_angle;
		lenses.ellipticity[i] = lens[i].ellipticity;
		lenses.ellipticity_potential[i] = lens[i].ellipticity_potential;
		lenses.rcore[i] = lens[i].rcore;
		lenses.rcut[i] = lens[i].rcut;
		lenses.z[i] = lens[i].z;
	}
	L[0]=lenses;
	// Computation
	grad = module_potentialDerivatives_totalGradient_SOA(N, &image, L);
	if (fabs(gradtheo.x - grad.x) > EPS) {
		std::cerr << "SIS grad.x: Theoretically we should have " << gradtheo.x
				<< "  and  " << grad.x << std::endl;
		exit(0);
	}
	if (fabs(gradtheo.y - grad.y) > EPS) {
		std::cerr << "SIS grad.y: Theoretically we should have " << gradtheo.y
				<< "  and  " << grad.y << std::endl;
		exit(0);
	}
	std::cout << "SIS Everythings fine, No test failed " << std::endl;
}

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
