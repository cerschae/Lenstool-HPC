extern double angdist(double z1, double z2, double omm, double oml, double lilh);
extern double radeigenj(double zbeta, double rsc, double rhosc, double Dl, double Dlsr, double Dsr, double lilh, double sigcrit, double mtol, double hiscale, double hr, double rhocrit, double x);
extern double radeigen(double zbeta, double rsc, double rhosc, double Dl, double Dlsr, double Dsr, double lilh, double sigcrit, double mtol, double hiscale, double hr, double rhocrit, double x);
extern double surfdens(double th);
extern double scmass(double x);
extern double taneigenj(double zbeta, double rsc, double rhosc, double Dl, double Dlst, double Dst, double lilh, double sigcrit, double rhocrit, double mtol, double hiscale, double hr, double x);
extern double taneigen(double zbeta, double rsc, double rhosc, double Dl, double Dlst, double Dst, double lilh, double sigcrit, double rhocrit, double mtol, double hiscale, double hr, double x);
extern double sigcrit(double Ds, double Dl, double Dls);
extern double rhocrit(double lilh, double omegam, double omegal, double z);
extern double sigjaffe(double mtol, double hiscale, double dmscale, double hmscale, double x);
extern double lumintj(double x);
extern double zbrentmodr(double (*fx)(double zbeta, double rsc, double rhosc, double Dl, double Dlsr, double Dsr, double lilh, double sigcrit, double mtol, double hiscale, double hr, double rhocrit, double x), double zbeta, double rsc, double rhosc, double Dl, double Dlsr, double Dsr, double lilh, double sigcrit, double mtol, double hiscale, double hr, double rhocrit, double x1, double x2, double tol);

extern double zbrentmodt(double (*func)(double zbeta, double rsc, double rhosc, double Dl, double Dlst, double Dst, double lilh, double sigcrit, double rhocrit, double mtol, double hiscale, double hr, double x), double zbeta, double rsc, double rhosc, double Dl, double Dlst, double Dst, double lilh, double sigcrit, double rhocrit, double mtol, double hiscale, double hr, double x1, double x2, double tol);

extern int zbracmodt(double (*func)(double zbeta, double rsc, double rhosc, double Dl, double Dlst, double Dst, double lilh, double sigcrit, double rhocrit, double mtol, double hiscale, double hr, double x), double zbeta, double rsc, double rhosc, double Dl, double Dlst, double Dst, double lilh, double sigcrit, double rhocrit, double mtol, double hiscale, double hr, double *x1, double *x2);

extern int zbracmodr(double (*fx)(double zbeta, double rsc, double rhosc, double Dl, double Dlsr, double Dsr, double lilh, double sigcrit, double mtol, double hiscale, double hr, double rhocrit, double x), double zbeta, double rsc, double rhosc, double Dl, double Dlsr, double Dsr, double lilh, double sigcrit, double mtol, double hiscale, double hr, double rhocrit, double *xrguessmin, double *xrguessmax);



extern double sighern(double mtol, double hiscale, double dmscale, double hmscale, double x);
extern double lumint(double x);
extern double qsimp(double (*func)(double), double a, double b);
extern double surfdens2(double x);
extern double einmass(double angle,double Ds, double Dl, double Dls);
