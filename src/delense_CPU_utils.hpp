#pragma once

#include <structure_hpc.hpp>
#include "gradient.hpp"

#define GET_INDEX2D(y, LDY, x, LDX)         (LDX*y + x)
#define GET_INDEX3D(y, LDY, x, LDX, z, LDZ) (LDX*LDZ*y + LDZ*x + z)

/** @brief Tranform a point from image to source plane. Result stored in sourcepoint argument
 *
 * Tranform a point from image to source plane using lensequation
 *
 * @param image_point    image position
 * @param dlsds          dls/ds
 * @param nhalos         number of halos
 * @param potential_param        gravitational potential information
 * @param source_point   address where source information will be stored
 *
 *
 */
static
void mychi_transformImageToSourcePlane(const runmode_param *runmode, const struct point *image_point, double dlsds, const struct Potential *lens, struct point *source_point)
{   // dlsds is the distance between lens and source divided by the distance observer-source
        struct point Grad;  // gradient

        Grad = module_potentialDerivatives_totalGradient(runmode->nhalos, image_point, lens);
        //Grad = module_potentialDerivatives_totalGradient_SOA(image_point, lens, runmode->Nlens);

        source_point->x = image_point->x - dlsds*Grad.x;
        source_point->y = image_point->y - dlsds*Grad.y;
        //printf("dlsds %f", dlsds);
}

static
void mychi_transformImageToSourcePlane_SOA(const int Nlens, const struct point *image_point, double dlsds, const struct Potential_SOA *lens, struct point *source_point)
{   // dlsds is the distance between lens and source divided by the distance observer-source
        struct point Grad;  // gradient
        Grad = module_potentialDerivatives_totalGradient_SOA(image_point, lens, Nlens);
        //
        source_point->x = image_point->x - dlsds*Grad.x;
        source_point->y = image_point->y - dlsds*Grad.y;
        //printf("dlsds %f", dlsds);
}
//
//
//
static
inline
void mychi_transformImageToSourcePlane_SOA_Packed( const struct point *image_point, double dlsds, struct point *source_point, double *grad_x, double * grad_y, int grad_id)
{

        source_point->x = image_point->x - dlsds*grad_x[grad_id];
        source_point->y = image_point->y - dlsds*grad_y[grad_id];
        //int world_rank;
        //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        //if (world_rank == 1)
        //printf("      %d: %f %f = %f %f - dlsds = %f grad id = %d grad = (%f %f)\n", world_rank, source_point->x, source_point->y,image_point->x, image_point->y, dlsds, grad_id, grad_x[grad_id], grad_y[grad_id]);
        //printf("dlsds %f", dlsds);
}
//
//
//
/** @brief Tranform a triangle from image to source plane. Result stored in S triangle argument
 *
 * Return a triplet of points in the source plane corresponding to the triplet
 * of images. dlsds is the lens efficiency at the source redshift.
 * I is the triangle in the image plane (input), S is the same triangle in the source plane (output)
 *
 * @param I      triangle in image plane
 * @param dlsds  dls/ds
 * @param nhalos         number of halos
 * @param potential_param        gravitational potential information
 * @param S      address where triangle source information will be stored
 *
 *
 */
static
inline
void mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_upper( struct triplet *I, double dlsds, struct triplet *S, double *grad_x, double * grad_y, int grad_id, int nbgridcell)
{
        mychi_transformImageToSourcePlane_SOA_Packed( &I->a, dlsds,   &S->a, grad_x, grad_y, grad_id             );
        mychi_transformImageToSourcePlane_SOA_Packed( &I->b, dlsds,   &S->b, grad_x, grad_y, grad_id + nbgridcell);
        mychi_transformImageToSourcePlane_SOA_Packed( &I->c, dlsds,   &S->c, grad_x, grad_y, grad_id +          1);
}
//
static
inline
void mychi_transformtriangleImageToSourcePlane_SOA_grid_gradient_lower( struct triplet *I, double dlsds, struct triplet *S, double*
grad_x, double * grad_y, int grad_id, int nbgridcell)
{
        mychi_transformImageToSourcePlane_SOA_Packed( &I->a, dlsds,   &S->a, grad_x, grad_y, grad_id + nbgridcell + 1);
        mychi_transformImageToSourcePlane_SOA_Packed( &I->b, dlsds,   &S->b, grad_x, grad_y, grad_id +              1);
        mychi_transformImageToSourcePlane_SOA_Packed( &I->c, dlsds,   &S->c, grad_x, grad_y, grad_id + nbgridcell    );
}



/** @brief Return the scalar triple product (a*b).c of the 3 vectors A[x,y,1], B[x,y,1], C[x,y,1].
 * If 2 of the 3 vectors are equal, colinear or form an orthogonal basis,
 * the triple product is 0.
 * This is also the determinant of the matrix
 *   | Ax  Bx  Cx |
 *   | Ay  By  Cy |
 *   |  1   1   1 |
 */
static
inline
double mychi_determinant(const struct point *A,
                const struct point *B,
                const struct point *C)
{
        return( B->x*C->y - B->y*C->x +
                        A->x*B->y - A->y*B->x +
                        A->y*C->x - A->x*C->y );
}
/** @brief Return 1 if P is inside the triangle T, 0 otherwise.
 * Return 1 if P is inside the triangle T, 0 otherwise.
 * @param P  a point
 * @param T  a triplet of points.
 *
 *
 */
static
inline
int mychi_inside(const struct point *P, struct triplet *T)
{
        double  s, s1, s2, d;

        d  = mychi_determinant(&T->a, &T->b, &T->c);
        s  = mychi_determinant(&T->a, &T->b, P)*d;
        if (s < 0.) return 0;
        s1 = mychi_determinant(&T->b, &T->c, P)*d;
        if (s1 < 0.) return 0;
        s2 = mychi_determinant(&T->c, &T->a, P)*d;
        if (s2 < 0.) return 0;
        return 1;

        //return((s > 0.) && (s1 > 0.) && (s2 > 0.));  // If all determinants are positive,
        // the point must be inside the triangle
}


/*
   int
   mychi_inside2(const struct point *A, const struct point *B, const struct point *C)
   {

// Compute vectors
v0 = C - A;
v1 = B - A;
v2 = P - A;

// Compute dot products
dot00 = dot(v0, v0);
dot01 = dot(v0, v1);
dot02 = dot(v0, v2);
dot11 = dot(v1, v1);
dot12 = dot(v1, v2);

// Compute barycentric coordinates
invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
u = (dot11 * dot02 - dot01 * dot12) * invDenom;
v = (dot00 * dot12 - dot01 * dot02) * invDenom;

// Check if point is in triangle
return (u >= 0) && (v >= 0) && (u + v < 1);
}
*/

/** @brief Return 1 if P is inside the triangle T or on its border, 0 otherwise.
 *
 * Return 1 if P is inside the triangle T or on its border, 0 otherwise.
 * @param  P  a point
 * @param  T  a triplet of points.
 *
 *
 */
static
inline
int mychi_insideborder(const struct point *P, struct triplet *T)
{
        double  s, s1, s2, d;

        d  = mychi_determinant(&T->a, &T->b, &T->c);
        s  = mychi_determinant(&T->a, &T->b, P)*d;
        if (s < 0.) return 0;
        s1 = mychi_determinant(&T->b, &T->c, P)*d;
        if (s1 < 0.) return 0;
        s2 = mychi_determinant(&T->c, &T->a, P)*d;
        if (s2 < 0.) return 0;
        return 1;
        //return((s >= 0.) && (s1 >= 0.) && (s2 >= 0.));  // If all determinants are positive or 0,
        // the point must be inside the triangle or on its border

}
/** @brief Barycentre of a triplet/triangle
 *
 * A is a structure triplet that contains 3 structures point a,b and c
 * Return value B is a point
 *
 *
 */
static
inline
struct  point   mychi_barycenter(struct triplet *A)
{
        struct  point   B;

        B.x = (A->a.x + A->b.x + A->c.x) / 3.;
        B.y = (A->a.y + A->b.y + A->c.y) / 3.;
        return(B);
}

/** @brief Euclidean distance between 2 points
 *
 * Euclidean distance between 2 points
 *
 */
static
inline
double mychi_dist(struct point A, struct point B)
{
        double  x, y;
        x = A.x - B.x;
        y = A.y - B.y;
        return(sqrt(x*x + y*y));
}






