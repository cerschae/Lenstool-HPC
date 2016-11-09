#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>

/****************************************************************/
/*                 Program  : grille				*/
/*                 Version  : Oct 2011			*/
/*                 Location :               */
/*                 Auteur   : Tomas Verdugo   */
/*		 This subprogram was added by TV         */
/*                                                  */
/****************************************************************/

void r_dynfile(FILE *IN,FILE *OUT)
{
	//extern struct	g_pot	P;
	extern struct	g_dyn	Dy;
	extern  struct  g_mode          M;
	extern struct pot	     lens[],lmin[];
	extern int             block[NLMAX][NPAMAX];
	char	second[20],third[120];
	
	fmot(IN,second);
	while(strcmp(second,"end"))
	{
		flire(IN,third);
		
		if (!strcmp(second,"dyntype"))
		{
	            sscanf(third,"%d",&Dy.dyntype); 
	            fprintf(OUT,"\t%s\t%d\n",second,Dy.dyntype);
		}
		else if (!strcmp(second,"dynnumber"))
		{
			sscanf(third,"%d",&Dy.dynnumber); 
			fprintf(OUT,"\t%s\t%d\n",second,Dy.dynnumber);
		}		
		else if (!strcmp(second,"velocity"))
		{
	            sscanf(third,"%lf",&Dy.dynvel); 
	            fprintf(OUT,"\t%s\t%lf\n",second,Dy.dynvel);
		}
		else if (!strcmp(second,"e_velocity"))
		{
	            sscanf(third,"%lf",&Dy.dynevel); 
	            fprintf(OUT,"\t%s\t%lf\n",second,Dy.dynevel);
		}
		else if (!strcmp(second,"indmass"))
		{
			sscanf(third,"%lf",&Dy.indmass); 
			fprintf(OUT,"\t%s\t%lf\n",second,Dy.indmass);
		}
		else if (!strcmp(second,"e_indmass"))
		{
			sscanf(third,"%lf",&Dy.indemass); 
			fprintf(OUT,"\t%s\t%lf\n",second,Dy.indemass);
		}
		else if (!strcmp(second,"refradius_kpc"))
		{
			sscanf(third,"%lf",&Dy.refradius); 
			fprintf(OUT,"\t%s\t%lf\n",second,Dy.refradius);
		}
		
		
		
		

		// Read the next line
		fmot(IN,second);
	}
	
	fprintf(OUT,"\t%s\n",second);
	
	
	///////	
	//The following are some errors so that the code doesn't run in vain	
	///////	
	
	
	
	if (Dy.dyntype != 12){
		fprintf(stderr," ERROR: Dyn Opt only works with type 12 (NFW), at least now \n");
		exit(-1);
	}
	
	
	
	//if( (D.dynnumber == 1 || D.dynnumber == 2) && (block[0][RC] == 0 ||  block[0][B0] == 0) ){
		//fprintf(stderr," ERROR: In order to run dyn_opt, the core_radius_kpc and v_disp need to be free parameters \n");
		//exit(-1);
	//}
	
	
	//if( Dy.refradius <= lmin[0].rckpc ){
		//fprintf(stderr," ERROR: reference radius is less than the inferior limit in rs\n");
		//NPRINTF(stderr,"Values refradius < lmin_rckpc: %lf %s %lf\n",Dy.refradius,"<",lmin[0].rckpc);
		//exit(-1);
	//}
	
	

}
