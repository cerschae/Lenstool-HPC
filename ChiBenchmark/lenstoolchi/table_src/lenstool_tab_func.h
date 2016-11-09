extern double scmass(double x);
extern double surfdens2(double x, double alpha);
extern double qromo(double (*func)(double), double a, double b);
extern double midpnt(double (*func)(double), double a, double b, int n);
extern void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
extern double   nfwg_dpl(double r,double rs,double kappas,double alpha);
extern double nfwg_kappa(double r,double rs, double kappas, double alpha);
extern int	read_bin(int tot_count);
