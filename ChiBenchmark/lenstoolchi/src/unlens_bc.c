#include<stdio.h>
#include<math.h>
#include "fonction.h"
#include "constant.h"
#include"dimension.h"
#include "structure.h"

//#define DDEBUG
/****************************************************************/
/*      nom:        e_unlens            */
/*      auteur:     Jean-Paul Kneib         */
/*      date:       10/02/92            */
/*      place:      Toulouse
 *****************************************************************
 * This function performs a kind of iterative inversion process to get
 * the position of one arclet knowing the position of the source Bs.

 * Return in multib a list of arclet positions in the image plane
 * corresponding to the Bs source position. Return the number of
 * valid positions in multib.
 *
 * If it has been possible to find a triangle in the source plane that
 * contains Bs and is smaller than the minimum surface required (defined
 * by the DMIN constant), then multib[i] contains the barycenter
 * of this triangle in the image plane, otherwise it contains the observed
 * arclet position multi[i] and the warn global variable is incremented.
 *
 *
 * Parameters :
 * - n : number of arclets in multi[]
 * - multi[] : list of arclets for a source
 * - Ps[] : list of sources corresponding to the arclets of multi[]
 * - Bs : barycenter in the source plane of the list of arclets in multi[]
 *
 * Global varibles used :
 * - distmin (constant)
 * - in e_im_prec() : it, G, lens, lens_table
 * - in e_transform() : G, lens, lens_table
 * - in e_amp() : G, lens, lens_table
 */
void e_powell( struct point Bs, double dlsds,
               struct point *p, struct point *xi,
               double ftol, double *fret);
void amoeba(struct point Bs, double dlsds, struct point *p, double *y, double ftol );


//function for single arclet 
int  unlens_bc_single(struct point Psi,
		      struct point Bs,
		      struct galaxie *multii,
		      struct point *multibi,
		      int n_famille);

//det_stop: if det_stop than we should stop imediately 
//          It can happen in openmp mode then chi2_img 
//          have already decided return -1. 
int  unlens_bc(const struct point *Ps,
               struct point Bs,
               struct galaxie *multi,
               struct point *multib,
               int n,
	       int n_famille,
	       int *det_stop)
{
    const extern  double  distmin[NFMAX];

    int    nimages; // number of valid images in multib
    int    i, j;
    double d, dmax;


    nimages = 0;

    for (i = 0; i < n; i++)
    {
       if (*det_stop)
	 continue;
       if (unlens_bc_single(Ps[i], Bs, &(multi[i]), &(multib[i]), n_famille))
	 nimages++;
    }

    // Reorder the predicted images to make them closer to their observed counterparts
    for( i = 0; i < n; i++ )
    {
        int jclose = -1;
        dmax = 100000.; // in arcsec
        for( j = 0; j < n; j++ )
        {
            d = dist(multi[i].C, multib[j]);
            if( d < dmax )
            {
                dmax = d; 
                jclose = j;
            }
        }

        // swap the predicted points to match the closer observed images
        if( jclose != i )
        {
            struct point tmp;
            tmp.x = multib[i].x;
            tmp.y = multib[i].y;
            multib[i].x = multib[jclose].x;
            multib[i].y = multib[jclose].y;
            multib[jclose].x = tmp.x;
            multib[jclose].y = tmp.y;
        }
    }

    // Test if 2 images are at the same place
    multib[n] = multib[0];
    for( i = 0; i < n; i++ )
    {
        if( dist(multib[i], multib[i+1]) < 0.1 )
        {
            // compute the distances to the observed images, and keep the smallest one for stats in o_chires
            double d1 = dist(multib[i], multi[i].C); 
            double d2 = dist(multib[i+1], multi[i+1].C); 
            if( d1 > d2 )
            {
                multib[i].x = multi[i].C.x;
                multib[i].y = multi[i].C.y;
            }
            else
            {
                multib[i+1].x = multi[i+1].C.x;
                multib[i+1].y = multi[i+1].C.y;
            }
        }
    }
    
    return nimages; // number of valid images in multib
}

//function for single arclet 
int  unlens_bc_single(struct point Psi,
		      struct point Bs,
		      struct galaxie *multii,
		      struct point *multibi,
		      int n_famille)
{                                                                                        
   const extern  double  distmin[NFMAX];

   struct bitriplet TE;
   struct bitriplet Tfinal; // list of final triangles found. (source and image planes)
   // bi-triangle definition

   /*  TE.i.a=multii->C;
    TE.i.b.x=TE.i.a.x +2.*(Bs.x-Psi.x) -(Bs.y-Psi.y);
    TE.i.b.y=TE.i.a.y +2.*(Bs.y-Psi.y) -(Bs.x-Psi.x);
    TE.i.c.x=TE.i.a.x +2.*(Bs.x-Psi.x) +(Bs.y-Psi.y);
    TE.i.c.y=TE.i.a.y +2.*(Bs.y-Psi.y) +(Bs.x-Psi.x); */

   /*distance between the barycenter of the sources and a particular source
    * normalized by the amplification factor at the arclet position*/
   double amp = fabs(e_amp_gal(multii, NULL)); 
   double dd = 0.7 * dist(Bs, Psi) / amp;  // distance scaled back to image plane

   /*minimal distance between the images of a familly*/
   double d_min = 0.4 * distmin[n_famille];
   
   /*if the distance in the source plane is larger than the
    * smallest distance between 2 images*/
   if ( d_min < dd )
     dd = d_min;
   
   /* TE.i represents an equilateral triangle with Cx,Cy at its center*/
   TE.i.a.x = multii->C.x;
   TE.i.a.y = multii->C.y + 2.*dd;
   TE.i.b.x = multii->C.x + 1.7 * dd;
   TE.i.b.y = multii->C.y - dd;
   TE.i.c.x = multii->C.x - 1.7 * dd;
   TE.i.c.y = multii->C.y - dd;
   

   /* Compute a triplet of simulated sources for the corresponding
    * triplet of simulated images. We assume that Bs is in TE.s */
   e_transform(&TE.i, multii->dr, &TE.s);
   
   /*count the loops needed to reach the smallest source triangle*/
   int it = 0;
   
   /*return the smallest couple of triangles in which Bs is located
    * in the source triangle TE.s and the image triangle TE.i contains
    * or is close to the current arclet position multii->C in the image
    * plane. Tfinal.i is a subtriangle of TE.i and Tfinal.s its associated
    * triangle in the source plane.
    * It's a kind of dichotomic process to get the
    * position of one arclet and its source with precision. */
   
   e_im_prec(&TE, &Bs, multii->dr, &it, &Tfinal);
   
   
   //      Di=dist(barycentre(Tfinal.i),multii->C);
   //      CG.x=Tfinal.i.a.x;
   //      CG.y=Tfinal.i.a.y;
   // NPRINTF(stderr,"Ai=%.3lf Aig=%.3lf Di=%.3lf\n",1./e_amp(multii->C,dlsds),1./e_amp(CG,dlsds),Di);
   // 
   double D = dist(barycentre(&Tfinal.s), Bs) / amp; //,Di;
   
   if ( D < 0.1 )
     {
	*multibi = barycentre(&Tfinal.i);
	
	return 1; //nimages++
     }
   else
     //Bs is too far and cannot be reached by the e_im_prec() function
     {	
	*multibi = multii->C;
	//NPRINTF(stderr,"WARNING: could not find the searched image\n");
	// NPRINTF(stderr,"Ai=%.3lf Aig=%.3lf D:%lf i:%d j:%d it:%d\n",1./e_amp(multii->C,dlsds),1./e_amp(CG,dlsds),D,n_famille,i,it);
	// 
	return 0; //do not count the current image
     }
}
