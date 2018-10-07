//
// This file is part of lenstoolhpc
// authors: gilles.fourestey@epfl.ch
// 
#include <math.h>
#include <cuda_runtime.h>
#include "module_cosmodistances.hpp"
#include "module_readParameters.hpp"

#include "allocation.hpp"

void PotentialSOAAllocation(Potential_SOA **lens_SOA, const int nhalos)
{
	
	Potential_SOA* p;
#if (defined __WITH_GPU)  && (defined __UNIFIED_MEM) 
#warning "Using unified memory"
	cudaError error = cudaMallocManaged(lens_SOA, sizeof(Potential_SOA));
	//if (error == 0) printf("Allocation error\n");
	p = *lens_SOA;
        cudaMallocManaged(&(p->type), nhalos*sizeof(int));
        cudaMallocManaged(&p->position_x         , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->position_y           , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->b0                   , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->ellipticity_angle    , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->ellipticity          , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->ellipticity_potential, nhalos*sizeof(size_t));
        cudaMallocManaged(&p->rcore                , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->rcut                 , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->vdisp                 , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->z                    , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->anglecos             , nhalos*sizeof(size_t));
        cudaMallocManaged(&p->anglesin             , nhalos*sizeof(size_t));
	cudaMallocManaged(&p->SOA_index            , nhalos*sizeof(int));
	//printf("Unified memory: %p\n", p->type); fflush(stdout);
#else
	//p = *lens_SOA;
	//printf("p = %p, %p\n", lens_SOA, p);
	*lens_SOA = (Potential_SOA*) malloc(sizeof(Potential_SOA));
	p = *lens_SOA;
	//printf("p = %p, %p\n", lens_SOA, p);
        p->type 		 = (int*) malloc(sizeof(int)*nhalos);
	p->position_x 		 = (double*) malloc(sizeof(double)*nhalos);
	p->position_y 		 = (double*) malloc(sizeof(double)*nhalos);
	p->b0 			 = (double*) malloc(sizeof(double)*nhalos);
	p->ellipticity_angle 	 = (double*) malloc(sizeof(double)*nhalos);
	p->ellipticity 		 = (double*) malloc(sizeof(double)*nhalos);
	p->ellipticity_potential = (double*) malloc(sizeof(double)*nhalos);
	p->rcore 		 = (double*) malloc(sizeof(double)*nhalos);
	p->rcut 		 = (double*) malloc(sizeof(double)*nhalos);
	p->vdisp		 = (double*) malloc(sizeof(double)*nhalos);
	p->z 	                 = (double*) malloc(sizeof(double)*nhalos);
	p->anglecos 		 = (double*) malloc(sizeof(double)*nhalos);
	p->anglesin 		 = (double*) malloc(sizeof(double)*nhalos);
	p->SOA_index 		 = (int*) malloc(sizeof(int)*nhalos);
        //lens_SOA->position_x            = new type_t[nhalos];
        //lens_SOA->position_y            = new type_t[nhalos];
        //lens_SOA->b0                    = new type_t[nhalos];
        //lens_SOA->ellipticity_angle     = new type_t[nhalos];
        //lens_SOA->ellipticity           = new type_t[nhalos];
        //lens_SOA->ellipticity_potential = new type_t[nhalos];
        //lens_SOA->rcore                 = new type_t[nhalos];
        //lens_SOA->rcut                  = new type_t[nhalos];
        //lens_SOA->z                     = new type_t[nhalos];
        //lens_SOA->anglecos              = new type_t[nhalos];
        //lens_SOA->anglesin              = new type_t[nhalos];
#endif
}



void PotentialSOADeallocation(Potential_SOA *lens_SOA)
{
#if (defined __WITH_GPU)  && (defined __UNIFIED_MEM) 
#warning "Using unified memory"
	cudaFree(lens_SOA->type);
	cudaFree(lens_SOA->position_x);
	cudaFree(lens_SOA->position_y);
	cudaFree(lens_SOA->b0);
	cudaFree(lens_SOA->ellipticity_angle);
	cudaFree(lens_SOA->ellipticity);
	cudaFree(lens_SOA->ellipticity_potential);
	cudaFree(lens_SOA->rcore);
	cudaFree(lens_SOA->rcut);
	cudaFree(lens_SOA->vdisp);
	cudaFree(lens_SOA->z);
	cudaFree(lens_SOA->anglecos);
	cudaFree(lens_SOA->anglesin);
	cudaFree(lens_SOA->SOA_index);
#else
	free(lens_SOA->type);
	free(lens_SOA->position_x);
	free(lens_SOA->position_y);
	free(lens_SOA->b0);
	free(lens_SOA->ellipticity_angle);
	free(lens_SOA->ellipticity);
	free(lens_SOA->ellipticity_potential);
	free(lens_SOA->rcore);
	free(lens_SOA->rcut);             
	free(lens_SOA->z);             
	free(lens_SOA->anglecos);
	free(lens_SOA->anglesin);       
	free(lens_SOA);
#endif
}


void read_potentialSOA_ntypes(std::string infile, int N_type[])
{
        int ind = 0;
        std::string first, second, third, line1, line2;

        std::ifstream IN(infile.c_str(), std::ios::in);
        //First sweep throught the runmode file to find N_type (number of types)
        if ( IN )
        {
                while(std::getline(IN,line1))
                {
                        std::istringstream read1(line1); // create a stream for the line
                        read1 >> first;
                        if (!strncmp(first.c_str(), "potent", 6))  // Read in potential
                        {
                                while(std::getline(IN,line2))
                                {
                                        std::istringstream read2(line2);
                                        read2 >> second >> third;
                                        if (strcmp(second.c_str(), "end") == 0)  // Move to next potential and initialize it
                                        {
                                                break; // Break while loop and move to next potential
                                        }

                                        if ( !strcmp(second.c_str(), "profil") ||  // Get profile
                                                        !strcmp(second.c_str(), "profile") )
                                        {
                                                if(!strcmp(third.c_str(), "PIEMD") ||
                                                                !strcmp(third.c_str(), "1") )
                                                {
                                                        ind=atoi(third.c_str());
                                                        N_type[ind] += 1;
                                                }

                                                else if(!strcmp(third.c_str(), "NFW") ||
                                                                !strcmp(third.c_str(), "2") )
                                                {
                                                        ind=atoi(third.c_str());
                                                        N_type[ind] += 1;
                                                }
                                                else if(!strcmp(third.c_str(), "SIES") ||
                                                                !strcmp(third.c_str(), "3") )
                                                {
                                                        ind=atoi(third.c_str());
                                                        N_type[ind] += 1;
                                                }
                                                else if(!strncmp(third.c_str(), "point", 5) ||
                                                                !strcmp(third.c_str(), "4") )
                                                {
                                                        ind=atoi(third.c_str());
                                                        N_type[ind] += 1;
                                                }
                                                else if(!strcmp(third.c_str(), "SIE") ||
                                                                !strcmp(third.c_str(), "5") )
                                                {
                                                        ind=atoi(third.c_str());
                                                        N_type[ind] += 1;
                                                }
                                                else if(!strcmp(third.c_str(), "8") )   //PIEMD
                                                {
                                                        ind=atoi(third.c_str());
                                                        N_type[ind] += 1;
                                                }
                                                else if(!strcmp(third.c_str(), "81") )  //PIEMD81
                                                {
                                                        ind=atoi(third.c_str());
                                                        N_type[ind] += 1;
                                                        //std::cerr << "Type First: " << ind << std::endl;
                                                }
                                                else if(!strcmp(third.c_str(), "14") )  //PIEMD81
                                                {
                                                        ind=atoi(third.c_str());
                                                        N_type[ind] += 1;
                                                        //std::cerr << "Type First: " << ind << std::endl;
                                                }
                                                else{
                                                        printf( "ERROR: Unknown Lensprofile, Emergency stop\n");
                                                        exit (EXIT_FAILURE);
                                                }
                                        }
                                }
                        }
                }
        }

        IN.close();
        IN.clear();
        IN.open(infile.c_str(), std::ios::in);

}



void module_readParameters_PotentialSOA_local(std::string infile, Potential_SOA *lens_SOA, int nhalos, int n_tot_halos, cosmo_param cosmology)
{
	//printf("lenses_SOA = %p\n", lens_SOA); fflush(stdout);
	//printf("lenses_SOA->type = %p\n", lens_SOA->type); fflush(stdout);
	//
        double DTR=acos(-1.)/180.;      /* 1 deg in rad  = pi/180 */
	//
        double core_radius_kpc = 0.;
        double cut_radius_kpc  = 0.;

        int N_type     [100];
        int Indice_type[100];
        int ind, initial_index;
        Potential lens_temp;
        //Used to store the initial index of lenses
        initial_index = 0;

        //Init of N_types and Indice_type (Number of lenses of a certain type)
        for(int i=0;i < 100; ++i){
                N_type[i] = 0;
                Indice_type[i] = 0;
        }
        //First sweep through the runmode file to find N_type (number of types)
        read_potentialSOA_ntypes(infile, N_type);

        //Calcuting starting points for each type in lens array
        for(int i=1;i < 100; ++i)
	{
                Indice_type[i] = N_type[i]+Indice_type[i-1];
                //printf("%d %d \n ",N_type[i], Indice_type[i]);
        }

        std::string first, second, third, line1, line2;
        std::ifstream IN(infile.c_str(), std::ios::in);
        if(IN){
                while(std::getline(IN,line1))
                {
                        first = "";
                        std::istringstream read1(line1); // create a stream for the line
                        read1 >> first;
                        //std::cerr << " 1: " << first << std::endl;
                        if (!strncmp(first.c_str(), "potent", 6))  // Read in potential
                        {

                                lens_temp.position.x = lens_temp.position.y = 0.;
                                lens_temp.ellipticity = 0;
                                lens_temp.ellipticity_potential = 0.;
                                lens_temp.ellipticity_angle = 0.;
                                lens_temp.vdisp = 0.;
                                lens_temp.rcut = 0.;
                                lens_temp.rcore = 0;
                                lens_temp.b0 = 0;
                                core_radius_kpc = 0.;
                                cut_radius_kpc = 0;
                                lens_temp.weight = 0;
                                lens_temp.rscale = 0;
                                lens_temp.exponent = 0;
                                lens_temp.alpha = 0.;
                                lens_temp.einasto_kappacritic = 0;
                                lens_temp.z = 0;
                                while(std::getline(IN,line2))
                                {
                                        //Init temp potential



                                        std::istringstream read2(line2);
                                        read2 >> second >> third;
                                        //std::cerr << line2 << std::endl;
                                        //std::cerr << " 2: " << second << std::endl;
                                        if (strcmp(second.c_str(), "end") == 0)  // Move to next potential and initialize it
                                        {
                                                if ( lens_temp.z == 0. )  // Check if redshift from current halo was initialized
                                                {
                                                        fprintf(stderr, "ERROR: No redshift defined for potential at position x: %f and y: %f \n", lens_temp.position.x , lens_temp.position.y);
                                                        exit(-1);
                                                }
                                                break; // Break while loop and move to next potential
                                        }

                                        //Find profile
                                        if ( !strcmp(second.c_str(), "profil") ||  // Get profile
                                                        !strcmp(second.c_str(), "profile") )
                                        {
                                                lens_temp.type=atoi(third.c_str());
                                                //std::cerr << lens_temp.type << std::endl;
                                        }

                                        else if (!strcmp(second.c_str(), "name"))    // Get name of lens
                                        {
                                                sscanf(third.c_str(),"%s",lens_temp.name);
                                        }

                                        else if (!strcmp(second.c_str(), "x_centre") ||  // Get x center
                                                        !strcmp(second.c_str(), "x_center") )
                                        {
                                                lens_temp.position.x=atof(third.c_str());
                                                //std::cout << "PositionX : " << std::setprecision(15) << lens_temp.position.x << std::endl;
                                        }
                                        else if (!strcmp(second.c_str(), "y_centre") ||  // Get y center
                                                        !strcmp(second.c_str(), "y_center") )
                                        {
                                                lens_temp.position.y=atof(third.c_str());
                                        }
                                        else if ( !strcmp(second.c_str(), "ellipticitymass") || !strcmp(second.c_str(), "ellipticity") || !strcmp(second.c_str(), "ellipticite") )  // Get ellipticity
                                        {
                                                lens_temp.ellipticity=atof(third.c_str());
                                                //lens_temp.ellipticity=lens_temp.ellipticity/3.;
                                        }
                                        else if (!strcmp(second.c_str(), "ellipticity_angle") || !strcmp(second.c_str(), "angle_pos"))  // Get ellipticity angle
                                        {
                                                lens_temp.ellipticity_angle=atof(third.c_str());
                                                lens_temp.ellipticity_angle *= DTR;
                                        }
                                        else if ( !strcmp(second.c_str(), "rcore") || !strcmp(second.c_str(), "core_radius"))  // Get core radius
                                        {
                                                lens_temp.rcore=atof(third.c_str());
                                        }
                                        else if (!strcmp(second.c_str(), "rcut") || !strcmp(second.c_str(), "cut_radius"))  // Get cut radius
                                        {
                                                lens_temp.rcut=atof(third.c_str());
                                        }
                                        else if ( !strcmp(second.c_str(), "core_radius_kpc"))  // Get core radius
                                        {
                                                core_radius_kpc=atof(third.c_str());
                                        }
                                        else if (!strcmp(second.c_str(), "cut_radius_kpc"))  // Get cut radius
                                        {
                                                cut_radius_kpc=atof(third.c_str());
                                        }
                                        else if (!strcmp(second.c_str(), "NFW_rs") ||  // Get scale radius of NFW
                                                        !strcmp(second.c_str(), "rscale"))
                                        {
                                                lens_temp.rscale=atof(third.c_str());
                                        }
                                        else if (!strcmp(second.c_str(), "exponent") )  // Get exponent
                                        {
                                                lens_temp.exponent=atof(third.c_str());
                                        }
                                        else if (!strcmp(second.c_str(), "alpha") )  // Get alpha
                                        {
                                                lens_temp.alpha=atof(third.c_str());
                                        }
                                        else if (!strcmp(second.c_str(), "einasto_kappacritic") ||  // Get critical kappa
                                                        !strcmp(second.c_str(), "kappacritic"))
                                        {
                                                lens_temp.einasto_kappacritic=atof(third.c_str());
                                        }
                                        else if (!strcmp(second.c_str(), "z_lens"))  // Get redshift
                                        {
                                                lens_temp.z=atof(third.c_str());
                                                //std::cerr << lens_temp.z << std::endl;
                                        }
                                        else if (!strcmp(second.c_str(), "v_disp"))  // Get Dispersion velocity
                                        {
                                                lens_temp.vdisp=atof(third.c_str());
                                                //std::cerr << "vdisp : "<< third << " " << lens_temp.vdisp << std::endl;
                                        }
                                        else if ( !strncmp(second.c_str(), "virial_mass", 6) ||  // Get virial mass
                                                        !strcmp(second.c_str(), "masse") ||
                                                        !strcmp(second.c_str(), "m200") ||
                                                        !strcmp(second.c_str(), "mass") )
                                        {
                                                lens_temp.weight=atof(third.c_str());
                                        }



                                } // closes inner while loop

                                // Converting distance in kpc to arcsec.
                                double d1 = d0/cosmology.h*module_cosmodistances_observerObject(lens_temp.z,cosmology);
                                //printf(" D1 HPC : %f %f %f %f\n",d1, d0,cosmology.h,lens_temp.z );
                                // Set rcore value in kpc or in arcsec.
                                if ( core_radius_kpc != 0. )
                                        lens_temp.rcore = core_radius_kpc / d1;
                                else
                                        core_radius_kpc = lens_temp.rcore * d1;

                                // Set rcut value in kpc or in arcsec.
                                if ( core_radius_kpc != 0. )
                                        lens_temp.rcore = core_radius_kpc / d1;
                                else
                                        core_radius_kpc = lens_temp.rcore * d1;

                                // Set rcut value in kpc or in arcsec.
                                if ( cut_radius_kpc != 0. )
                                {
                                        //std::cerr << "d1 " << d1 << std::endl;
                                        lens_temp.rcut = cut_radius_kpc / d1;}
                                else
                                        cut_radius_kpc = lens_temp.rcut * d1;

                                //Calculate parameters like b0, potential ellipticity and anyother parameter depending on the profile
                                module_readParameters_calculatePotentialparameter(&lens_temp);

                                //assign value to SOA
                                //std::cerr << "Type + indice :" << lens_temp.type << Indice_type[lens_temp.type-1] << std::endl;
                                //printf("%p %p\n", lens_SOA, lens_SOA->type[0]);
                                if(Indice_type[lens_temp.type-1] <nhalos)
				{
                                        ind = Indice_type[lens_temp.type-1];
                                        //std::cerr<< ind << std::endl;
                                        lens_SOA->type[ind]                  = lens_temp.type;
                                        lens_SOA->position_x[ind]            = lens_temp.position.x;
                                        lens_SOA->position_y[ind]            = lens_temp.position.y;
                                        lens_SOA->b0[ind]                    = lens_temp.b0;
                                        lens_SOA->vdisp[ind]                 = lens_temp.vdisp;
                                        lens_SOA->ellipticity_angle[ind]     = lens_temp.ellipticity_angle;
                                        lens_SOA->ellipticity[ind]           = lens_temp.ellipticity;
                                        lens_SOA->ellipticity_potential[ind] = lens_temp.ellipticity_potential;
					lens_SOA->vdisp[ind]                 =          lens_temp.vdisp;
                                        lens_SOA->rcore[ind]                 = lens_temp.rcore;
                                        lens_SOA->rcut[ind]                  = lens_temp.rcut;
                                        lens_SOA->z[ind]                     = lens_temp.z;
                                        lens_SOA->anglecos[ind]              = cos(lens_temp.ellipticity_angle);
                                        lens_SOA->anglesin[ind]              = sin(lens_temp.ellipticity_angle);
                                        //Store new index for bayes map purposes
                                        lens_SOA->SOA_index[initial_index] = ind;
					//
                                        initial_index += 1;
                                        Indice_type[lens_temp.type-1] += 1;
                                }
                        }  // closes if loop

                }  // closes while loop
        }
        IN.close();
}

