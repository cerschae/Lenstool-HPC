#include<stdio.h>
#include<time.h>
#include<string.h>
#include<math.h>
#include<float.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>


/****************************************************************/
/*      nom:        set_default         */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/03/94            */
/*      place:      ESO LaSilla         */
/****************************************************************/

void    set_default()
{
    extern struct g_mode    M;
    extern struct g_source  S;
    extern struct g_image   I;
    extern struct g_grille  G;
    extern struct g_msgrid  H;
    extern struct g_frame   F;
    extern struct g_large   L;
    extern struct g_cosmo   C;
    extern struct g_cline   CL;
    extern struct g_observ  O;/* default flux dispersion^2 */
    extern struct g_pot     P[NPOTFILE];
    extern struct ipot      ip;
    extern struct sigposStr sigposAs;
    extern struct g_pixel   ps, imFrame;
    extern struct g_dyn     Dy;
    extern struct z_lim     zalim;
    int i;

    /* runmode */

    M.verbose = 0;
    M.seeing = 0.;
    M.sort = 0;
    M.minchi0 = 0.1;
    M.ichi2 = 0;
    M.source = 0;
    M.image = 0;
    M.centerfile[0] = 0;
    M.pixel = 0;
    M.iclean = 0;
    M.cube = 0;

    /* image */

    I.forme = 1;
    I.mult_abs = 0;
    I.Anfilt = 15;
    I.Anpixseeing = 6;
    I.sigell = 0.3;
    I.dsigell = -1.;
    I.sig2amp = 0.04; /* default flux dispersion^2 */
    I.Dmag = 0.1;
    I.zarclet = -1.;

    // z_a_limit
    zalim.bk = 0;
    zalim.min = -1;
    zalim.max = -1;

    // sigposArcsec position dispersion
    sigposAs.bk = 0; // not optimised
    sigposAs.min = 0.2;
    sigposAs.max = 4.;
    sigposAs.prec = 0.01;

    /* crtitical line */

    CL.nplan = 0;
    CL.zone = 0;
    CL.npzone = 100;
    strcpy(CL.zonefile, "zone.dat");
    CL.limitHigh = 10;  // (arcsec) limitHigh for marchingSquares
    CL.cpas = 1;        // (arcsec) limitLow for marchingSquares
    CL.dmax = 0.;       // field is computed from the champ section
    strcpy(CL.algorithm, "MARCHINGSQUARES" );

    /* large */

    strcpy(L.iname, "giant");
    L.dlarge = 2.;
    L.vitesse = -1;
    L.iso = 0;
    L.nmaxiso = 0;
    L.scale = 1;
    L.zonex = 0;
    L.zoney = 0;
    L.profil = 0;
    L.pt = 0;
    L.ncourbe = 0;
    L.npt = 0;

    /* grille */

    G.pol = 0;
    G.no_lens = 0;
    G.nlens = 0;
    G.ngrid = 30;
    G.nlens_crit = 0;
    G.npot = 0;
    G.nplens[0] = -1;
    G.nmsgrid = -1;  // default no multiscale grid

    // multi-scale grid
    H.threshold = 0.;
    H.levels = 3;
    H.param = 3;

    /* source */

    S.zs = 0.8;
    S.ns = 0;
    S.grid = 0;
    S.rand = ((int)time(NULL) % 100) - 50;  //-2;
    S.distz = 0;
    S.emax = 0.3;
    S.zsmin = 0.5;
    S.zsmax = 1.5;
    S.par1 = 0.;
    S.par2 = 0.;
    S.par3 = 0.;
    S.taille = 0.3;
    S.lfalpha = -1;
    S.lfm_star = -20;
    S.lfm_min = -24;
    S.lfm_max = -16;
    /* luminosity function (moved from r_source.c) */
    S.lfalpha = -1.11;
    S.lfm_star = -20.9;
    S.lfm_min = -25.;
    S.lfm_max = -16.9;


    /* frame */

    F.xmin = -20.;
    F.xmax = 20.;
    F.ymin = -20.;
    F.ymax = 20.;
    F.rmax = 25.;

    /* image plane */
    imFrame.xmin = -20.;
    imFrame.xmax = 20.;
    imFrame.ymin = -20.;
    imFrame.ymax = 20.;
    imFrame.pixelx = imFrame.pixely = 0.;
    imFrame.ech = 2.; // number of source plane pixels for 1 image plane pixel
    imFrame.nx = 40;
    imFrame.ny = 40;
    imFrame.column = 1;
    imFrame.pixfile[0] = 0;
    imFrame.ncont = 0;
    imFrame.outfile[0] = 0;
    for( i = 0; i < 10; i++ )
        imFrame.contfile[0][0] = 0;
    imFrame.header = 0;

    /* source plane */
    ps.xmin = -20.;
    ps.xmax = 20.;
    ps.ymin = -20.; 
    ps.ymax = 20.;
    ps.pixelx = ps.pixely = 1.;
    ps.ech = 1.; // pixel scale in the source plane (default no oversampling)
    ps.nx = 40;  // size of the source plane image
    ps.ny = 40;
    ps.pixfile[0] = 0;

    /* potfile */

    for( i = 0; i < NPOTFILE; i++ )
    {
        P[i].potid = i;
        P[i].select = 0;  // not selection by their Einstein radius
        P[i].mag0 = 17.;
        P[i].isigma = 0;
        P[i].ircut = 0;
        P[i].islope = 0;
        P[i].ivdslope = 0;
        P[i].ivdscat = 0;
        P[i].ircutscat = 0;
        P[i].type = 81;    // default PIEMD
        P[i].corekpc = P[i].core = -1.0;
        P[i].cutkpc1 = P[i].cut1 = P[i].cut = DBL_MAX;
        P[i].sigma1 = P[i].sigma = -1.0;
        P[i].slope1 = 4.;
        P[i].vdslope1 = 4.;
        P[i].vdscat = 0.;
        P[i].rcutscat = 0.;
        P[i].a1 = P[i].a = DBL_MAX;  // default not defined
        P[i].b1 = P[i].b = DBL_MAX;  
    }

    /* cosmology */
    
    C.model = 1;
    C.H0 = 50.;
    C.h = 1.;
    C.omegaM = 1.;
    C.omegaX = 0.;
    C.kcourb = 0.;
    C.wX = -1.;
    C.wa = 0.;

    /* observe */

    O.prec = 0.2;
    O.bruit = 0;
    O.idum = -1;
    O.setseeing = 0;
    O.setbin = 0;
    O.bin = 0;
    O.setseeing = 0;
    O.filtre = 0;
    O.seeing = 1.;
    O.r0st = 0;
    O.SKY = 0;
    O.gain = 1.;

    // Parameters
    ip.pmax = 11;
    ip.extend = 1;
    ip.map = 0;

    // Dynamics 
    Dy.dyntype = 0;
}

