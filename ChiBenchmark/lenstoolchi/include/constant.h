/*
* constants definition
*/
#ifndef CONSTANT_H
#define CONSTANT_H

#define PI 	3.141592653589793238462643
#define LOG10   2.302585092994046
#define INVG    2.325968e8	// in Msol/(km/s)^2/Mpc 
#define pi_in_a	648000.		/* pi en arcsecond */
#define RTD 	57.2957795130823208768	/* 1 rad en deg  = 180/pi */
#define DTR 	0.01745329251994329577	/* 1 deg en rad  = pi/180 */
#define RTA 	206264.81	/* 1 rad en arecsecond = 648000/pi */
#define pia_c2	7.209970e-06	/* pi en arcsecond/ c^2 =  648000/vol/vol */

#define vol	299792.50	/* en km/s	    cf. Weinberg*/
#define GM_c2	1.475		/* 1.e12 G.M_sol/c2 en 1.e12 km cf. Weinberg*/
#define kpot	4.59e-7		/* 2(arcsec)/c^2 in arecsec/(km/s)^2 */
#define ikpot	2178644.6	/* c^2/2(arcsec) in (km/s)^2/arecsec */
#define h0      50.
#define MCRIT	7.36126993e11	/* c^3/4Gh0/RTA/RTA in M_sol/arcsec^2 (h0=50) */
#define MCRIT12	.2343165	/* c^3/4piGh0/RTA/RTA in 1e12 M_sol/arcsec^2 (h0=50) */
#define PC      30.8563		/* 1 pc en 1.e12 km */
#define Mpc     30856300.		/* 1 Mpc en 1.e12 km */
#define D0Mpc   5995.85         /* c/H0 en Mpc */
#define D0pc    5995850000.     /* c/H0 en pc */
#define d0	29.068701	/* vol/h0*1000/rta -- c/H0 en h-1.kpc/arcsecond  (h0=50)*/
#define D0	896952.55
		 /* vol/h0*Mpc/rta -- c/H0 en  1.e12 h-1.km/arcsecond  */
#define cH2piG  0.11585881	/* cH0/2piG en g/cm^2 (H0=50) */
#define cH4piG  0.057929405	/* cH0/4piG en g/cm^2 (H0=50) */
#define cH0_4piG  2.7730112e-4	/* cH0/4piG en 10^12 M_Sol/kpc^2 (H0=50) */
#define E2max	.125
#define tyaa	.45964488	/* (1/H0)*(pi/648000)^2 (h0=50) in year */
#define th_a2_day 167.8852919   /* (1/H0)*(pi/648000)^2 (h0=50) in days */

#define Msol_in_g	1.989e33
#define cm_in_arcsec	4.982444e15
#define inte	0.012537756 /* cm_in_arcsec^2/Msol_in_g  */

#endif // if CONSTANT_H
