#include<math.h>
#include<fonction.h>

/* -------------------------------------------------------------------
generate a redshift for a source following a Smail et al. distribution
------------------------------------------------------------------*/

// Sampling of the Smail et al. function. Make sure this number is significantly larger 
// than your sample of objects for which you want to get a redshift, so
// that each return redshift is unique
#define KSMAIL 1000

// Return a pointer to a structure that help speedup the following calculations
// (see gsl: General-Discrete-Distributions.html)
gsl_ran_discrete_t* smailpreproc()
{
    extern struct g_source S;
    int i;
    double z, p[KSMAIL];

    // Compute the Smail et al. array in [zsmin, zsmax]
    for( i = 0; i < KSMAIL; i++ )
    {
        z = S.zsmin + i * (S.zsmax - S.zsmin) / (KSMAIL-1);
        p[i] = pow(z, S.par1) * exp(-pow(z/S.par2, S.par3));
    }

    return gsl_ran_discrete_preproc(KSMAIL, p);
}

/* Return a redshift validating the Smail et al. function defined with parameters
 * in g_source structure.
 * Input :
 *   - r : random seed
 *   - g : gsl structure obtained from smailpreproc()
 */
double d_rndzsmail(gsl_rng * r, gsl_ran_discrete_t * g)
{
    extern struct g_source S;
    size_t i;
    double zrnd;

    i = gsl_ran_discrete (r, g);
    zrnd = S.zsmin + i * (S.zsmax - S.zsmin) / (KSMAIL - 1);

    return zrnd;
}

