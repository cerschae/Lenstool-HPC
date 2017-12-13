/**
* @file   module_readParameters.cpp
* @Author Christoph Schaefer, EPFL (christophernstrerne.schaefer@epfl.ch)
* @date   July 2015
* @version 0,1
* @brief  read parameters, images, sources and any other relevant information
*
* read parameter file
* 
*/



// Include
//===========================================================================================================
#include <iostream>
#include <cstring>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "module_readParameters.hpp"
#include <iomanip>





//===========================================================================================================
// Function definitions

/** @brief This module function reads the cosmology parameters from the parameter file
* 
* This module function reads the cosmology parameters from the parameter file. The given pointer is 
* updated with the results.
* Read an infile: The structure of the file is as follows:
	|    .
	|    .
	|    .
	|cosmology
	|		model	x <---- this is a value
	|		 H0	x
	|		  .	.
	|		  .	.
	|		 end
	|    .
	|    .
	|  finish
* 	
* @param Cosmology cosmology
* @param infile path to file	 
*/



void module_readParameters_readCosmology(std::string infile, cosmo_param &Cosmology)
{

    std::string first, second, third, line1, line2;

    Cosmology.model =1;


    std::ifstream IN(infile.c_str(),std::ios::in); // open file
    if ( IN )
    {
	std::getline(IN,line1);
	std::istringstream read1(line1); // create a stream for the line
	read1 >> first;
	while(strncmp(first.c_str(), "fini",4) != 0 ) // read line by line until finish is reached
        {
            if ( strncmp(first.c_str(), "cosmolog", 8) == 0){ // if the line contains information about cosmology

		    std::getline(IN,line2);    // read the line word by word
		    std::istringstream read2(line2);
		    read2 >> second >> third;

		    while(strncmp(second.c_str(), "end",3) != 0)    // Go ahead until "end" is reached
        	    { 

        		if ( !strcmp(second.c_str(), "model") )  // set model of universe
        		{                                                          
            		Cosmology.model=atoi(third.c_str());
        		}
       			else if ( !strcmp(second.c_str(), "H0") ) // set Hubble constant
        		{
            		Cosmology.H0=atof(third.c_str());
            		Cosmology.h = Cosmology.H0 / 50.;
        		}
        		else if ( !strcmp(second.c_str(), "omegaM") || !strcmp(second.c_str(), "omega") ) // set density of matter
        		{
           		Cosmology.omegaM=atof(third.c_str());
        		}
        		else if ( !strcmp(second.c_str(), "omegaX") || !strcmp(second.c_str(), "lambda") ) // set cosmological constant
        		{
            		Cosmology.omegaX=atof(third.c_str());
        		}
        		else if ( !strcmp(second.c_str(), "wX") || !strcmp(second.c_str(), "q") || !strcmp(second.c_str(), "w0") )     // set  "q" for Model 2, "wX" for Model 3, "w0" for Model 4
        		{
            		Cosmology.wX=atof(third.c_str());
        		}
        		else if ( !strcmp(second.c_str(), "wa") || !strcmp(second.c_str(), "n") || !strcmp(second.c_str(), "delta") || !strcmp(second.c_str(), "w1") )     // set  "n" for Model 2, "delta" for model 3, "w1" for model 4  
        		{
            		Cosmology.wa=atof(third.c_str());
        		}
        		else if ( !strcmp(second.c_str(), "omegaK") ) // set universe curvature
        		{
            		Cosmology.curvature=atof(third.c_str());
        		}
			// read next line
			std::getline(IN,line2);
		        std::istringstream read2(line2);
		        read2 >> second >> third;

    		   }

    		// if a flat Universe
    		if ( Cosmology.curvature == 0. )
    		{
        	Cosmology.omegaX = 1 - Cosmology.omegaM;
    		}
    		else
        	Cosmology.curvature = Cosmology.omegaM + Cosmology.omegaX - 1;
		
	    }
	    // read next line
	    std::getline(IN,line1);
	    std::istringstream read1(line1);
	    read1 >> first;

        }

        IN.close();

    }
    else
    {
        fprintf(stderr, "ERROR: file %s not found\n", infile.c_str());  // Exit if file was not found
        exit(-1);
    }

}

/** @brief This module function reads the number of sources, arclets, etc. in the parameter file. We need to know this to allocate the memory
* 
* Function to read the number of multiple images and clumps
* Check if there is 1 limit for each clump set
* The program reads the number of potentials and limits defined, and checks whether there are the same number
* Then it opens the imfile to read the numbers of images and set of images
* Reads also the size of the multiple images area, the size of a pixel, and the runmode.
* 	
* @param infile path to file
* @param runmode runmode parameter
*/

void read_runmode(std::istream &IN, struct runmode_param *runmode){
	std::string  first, second, third, fourth, fifth, line1, line2;
	//sscanf variables
	double in1, in2;	//%lf in1, then (type_t)in1 ;
	std::getline(IN,line2);
					std::istringstream read2(line2);
					read2 >> second >> third >> fourth;    // Read in 4 words
		    		while (strncmp(second.c_str(), "end", 3))    // Read until we reach end
		    		{
						if ( !strcmp(second.c_str(), "nbgridcells") )
		        		{
							sscanf(line2.c_str(), " %*s %d ", &runmode->nbgridcells);
		        		}

		        		if ( !strcmp(second.c_str(), "source") )
		        		{
							char filename[FILENAME_SIZE];
		            		sscanf(line2.c_str(), " %*s %d %s ", &runmode->source, &filename);
		            		runmode->sourfile = filename;
		        		}

						if ( !strcmp(second.c_str(), "image") )
		        		{
							char filename[FILENAME_SIZE];
		            		sscanf(line2.c_str(), " %*s %d %s ", &runmode->image, &filename);
		            		runmode->imagefile = filename;
		        		}
		        		if ( !strcmp(second.c_str(), "inverse") )
		        		{
							sscanf(line2.c_str(), " %*s %d ", &runmode->inverse);
		        		}
		        		if ( !strcmp(second.c_str(), "mass") )
		        		{
							sscanf(line2.c_str(), " %*s %d %d %lf %lf", &runmode->mass, &runmode->mass_gridcells, &in1, &in2);
							runmode->z_mass = (type_t)in1;
							runmode->z_mass_s =(type_t)in2;
		        		}
		        		if ( !strcmp(second.c_str(), "poten") )
		        		{
							sscanf(line2.c_str(), " %*s %d %d %lf", &runmode->potential, &runmode->pot_gridcells, &in1);
							runmode->z_pot = (type_t)in1;
		        		}
		        		if ( !strcmp(second.c_str(), "dpl") )
		        		{
							sscanf(line2.c_str(), " %*s %d %d %lf", &runmode->dpl, &runmode->dpl_gridcells, &in1);
							runmode->z_dpl = (type_t)in1;
		        		}
		        		if ( !strcmp(second.c_str(), "grid") )
		        		{
							sscanf(line2.c_str(), " %*s %d %d %lf", &runmode->grid, &runmode->gridcells, &in1);
							runmode->zgrid = (type_t)in1;
		        		}
		        		if ( !strcmp(second.c_str(), "ampli") )
		        		{
							sscanf(line2.c_str(), " %*s %d %d %lf", &runmode->amplif, &runmode->amplif_gridcells, &in1);
							runmode->z_amplif = (type_t)in1;
							//std::cerr<<runmode_>ampli << << <<std::endl;
									        		}
						if ( !strcmp(second.c_str(), "arclets") )
		        		{
		            		runmode->arclet = 1;  // Not supported yet
		        		}
						if ( !strcmp(second.c_str(), "time") )
		        		{
							sscanf(line2.c_str(), " %*s %d ", &runmode->time);

		        		}
						if( !strncmp(second.c_str(), "Debug",5) )    // Read in if we are in debug mode
						{
						runmode->debug = 1;
						}

						// Read the next line
						std::getline(IN,line2);
						std::istringstream read2(line2);
						read2 >> second >> third;
		    		}

}

void read_runmode_potential(std::istream &IN, int &numberPotentials){
	std::string  first, second, third, fourth, fifth, line1, line2;

	numberPotentials += 1;
	std::getline(IN,line2);
	std::istringstream read2(line2);
	double z(0);
	read2 >> second >> third;
	//std::cout << second  << third << std::endl;
	while (strncmp(second.c_str(), "end", 3))    // Read until we reach end
	{
		if (!strcmp(second.c_str(), "z_lens"))  // Get redshift
		{
			z=atof(third.c_str());
		}
		// Read the next line
		std::getline(IN,line2);
		std::istringstream read2(line2);
		read2 >> second >> third;
	}
	// Check if redshift from current halo was initialized
	if ( z == 0. )
		{
			fprintf(stderr, "ERROR: A redshift is not defined for a potential \n");
			exit(-1);
		}
}

void read_runmode_image(std::istream &IN, struct runmode_param *runmode){
	std::string  first, second, third, fourth, fifth, line1, line2;

    std::getline(IN,line2);
 	std::istringstream read2(line2);
 	double z(0);
 	read2 >> second >> third;
 	while (strncmp(second.c_str(), "end", 3))    // Read until we reach end
		{
			if ( !strcmp(second.c_str(), "multfile") )
			{
				char filename[FILENAME_SIZE];
				sscanf(line2.c_str(), " %*s %d %s ", &runmode->multi, &filename);
				runmode->imagefile = filename;
			}
			//std::cout << runmode->multi << runmode->imagefile << std::endl;
			// Read the next line
			std::getline(IN,line2);
			std::istringstream read2(line2);
			read2 >> second >> third;
		}
}

void read_runmode_potfile(std::istream &IN, struct runmode_param *runmode){
	std::string  first, second, third, fourth, fifth, line1, line2;

    runmode->potfile = 1;

    std::getline(IN,line2);
	std::istringstream read2(line2);
	read2 >> second >> third >> fourth;    // Read in 4 words
	while (strncmp(second.c_str(), "end", 3))    // Read until we reach end
	{

	if ( !strcmp(second.c_str(), "filein") )
	{
		//std::cerr<< third << " " << fourth << " " << runmode->potfilename << std::endl;
		runmode->potfilename = fourth;
		//std::cerr<< line2 << runmode->potfilename;
		break;

	}
	// Read the next line
	std::getline(IN,line2);
	std::istringstream read2(line2);
	read2 >> second >> third;
	}
}

void read_runmode_countimages(struct runmode_param *runmode){
	std::string line2;
	int old_j = 0;
	int j = 0;
	int imageIndex = 0;

	if (runmode->image == 1 or runmode->inverse == 1 or runmode->time >= 1 or runmode->multi >= 1){

		std::string imageFile_name = runmode->imagefile;
		std::ifstream IM(imageFile_name.c_str(), std::ios::in);

		if ( IM )
		{
			while( std::getline(IM,line2) )  // Read every line
			{
			if ( strncmp(line2.c_str(), "#", 1) == 0){	//Skipping commented lines
				continue;
			}
			else{
				int j=atoi(line2.c_str());
				if(j != old_j){				//If a new set has been reached, increase the nset iterator and update old j
					runmode->nsets +=1;
					old_j = j;
				}
				imageIndex++;
				}
			}
		}
		else{

			fprintf(stderr, "ERROR: file %s not found\n", imageFile_name.c_str());
			exit(-1);
		}

	runmode->nimagestot=imageIndex;
	IM.close();
	}
}

void read_runmode_countsources(struct runmode_param *runmode){
	std::string  imageFile_name,line1, line2;

	if (runmode->source == 1 and runmode->image == 0 and runmode->multi == 0 ){
		imageFile_name = runmode->sourfile;
		//printf("Source to image mode activated\n");

	std::ifstream IM(imageFile_name.c_str(), std::ios::in);
	//printf(" Booo says hello again \n");
    	if ( IM ){
			int i = 0;
	        while( std::getline(IM,line1) ){    // Read until we reach the end
	            i++;
			}
			runmode->nsets = i ;
	    	}
		else{

			fprintf(stderr, "ERROR: file %s not found\n", imageFile_name.c_str());
			exit(-1);
		}

	runmode->nimagestot= 0;	// Important
	IM.close();

	}
}

void read_runmode_countpotfile(struct runmode_param *runmode){
	std::string  imageFile_name,line1, line2;

	if (runmode->potfile == 1){
			imageFile_name = runmode->potfilename;

			std::ifstream IM(imageFile_name.c_str(), std::ios::in);

	    	if ( IM ){
				int i = 0;
		        while( std::getline(IM,line1) ){    // Read until we reach the end
	                if ( line1[0] == '#' ){
	                    continue;
	                }
		            i++;
				}
				runmode->npotfile = i ;
		    	}
			else{

				fprintf(stderr, "ERROR: file %s not found\n", imageFile_name.c_str());
				exit(-1);
			}

		IM.close();

		}
}


void module_readParameters_readRunmode(std::string infile, struct runmode_param *runmode)
{

	std::string  first, second, third, fourth, fifth, line1, line2;
	int Nsis(0), Npiemd(0);
	/// Set default values

	runmode->nbgridcells = 1000;
	runmode->source = 0;
	runmode->image = 0;
	runmode->nsets = 0;
	runmode->nhalos = 0;
	runmode->multi = 0;
	runmode->mass = 0;
	runmode->mass_gridcells = 1000;
	runmode->z_mass = 0.4;
	runmode->z_mass_s = 0.8;
	runmode->potential = 0;
	runmode->pot_gridcells = 1000;
	runmode->potfile = 0; //weird seg fault due to this line, haven't figured out why
	runmode->npotfile = 0;
	runmode->z_pot = 0.8;
	runmode->dpl = 0;
	runmode->dpl_gridcells = 1000;
	runmode->z_dpl = 0.8;
	runmode->inverse = 0;
	runmode->arclet = 0;
	runmode->debug = 0;
	runmode->nimagestot = 0;
	runmode->nsets = 0;
	runmode->gridcells = 1000;
	std::cerr << sizeof(*runmode) << std::endl;
	runmode->cline = 0;
	runmode->time = 0;


	int j=0;
	std::string imageFile_name;
	int imageIndex=0;
	int numberPotentials=0, numberLimits=0;

/*************************** read nhalos, imfile_name, pixel_size, multipleimagesarea_size and runmode from configuration file *********************/

std::ifstream IN(infile.c_str(), std::ios::in);
    if ( IN )
    {

	std::getline(IN,line1);
	std::istringstream read1(line1); // create a stream for the line
	read1 >> first;  // Read the first word
	//std::cout<<first;
        while(  strncmp(first.c_str(), "fini",4) != 0  )    // Continue until we reach finish
        {
			if ( strncmp(first.c_str(), "runmode", 7) == 0){    // Read in runmode information
				read_runmode(IN,runmode);
		    }
			else if (!strncmp(first.c_str(), "potent", 6)) // each time we find a "potential" string in the configuration file, we add an halo
            {
				read_runmode_potential(IN,numberPotentials);
				//std::cerr<<numberPotentials << std::endl;
            }
            else if (!strncmp(first.c_str(), "image", 5)) // same for the limit of the potential
            {
            	read_runmode_image(IN,runmode);
            }
            else if (!strncmp(first.c_str(), "limit", 5)) // same for the limit of the potential
            {
                numberLimits++;
            }
            else if (!strncmp(first.c_str(), "cline", 5))
            {
                runmode->cline = 1;
            }

            else if (!strncmp(first.c_str(), "potfile", 7))
            {
            	read_runmode_potfile(IN,runmode);
            }
	    // read the next line
	    std::getline(IN,line1);
	    std::istringstream read1(line1);
	    read1 >> first;
	    //std::cout<<first;
        }

        IN.close();
	
        //if(numberLimits!=numberPotentials) printf("Warning : Number of clumps different than number of limits in %s\n", infile.c_str()); // must be limits for each halo
	runmode->nhalos=numberPotentials;

    }
    else
    {
        fprintf(stderr, "ERROR: file %s not found\n", infile.c_str());
        exit(-1);
    }
    //
    //getting nimage value (number of images), nset value (number of sources) and npotfile
    //if image or multi mode is activated get nimage and nset
    read_runmode_countimages(runmode);
	//if source mode is activated, get nset
    read_runmode_countsources(runmode);
    //if potfile mode, count number of potential in potfile
    read_runmode_countpotfile(runmode);


std::cerr <<"nsets: " <<runmode->nsets <<" nhalos: " << runmode->nhalos << " nimagestot: " << runmode->nimagestot << " npotfile: " << runmode->npotfile << " multi: " << runmode->multi<< std::endl;

}




/** @brief read the positions, redshifts and numbers of multiple images from the images file
* 
* This module function reads in the strong lensing images
* Output: image coordinates, image shapes (semi-minor, semi-major of ellipse and orientation angle), source redshifts, number of images per set
* 	
** the images file must contain
	1 	x_center	y_center	major axis	minor axis	ellipticity angle	redshift
	1	   .		   .		   .		    .			.		    .
	1	   .		   .		   .		    .			.		    .
	2	   .		   .		   .		    .			.		    .
	2	   .		   .		   .		    .			.		    .
	2	   .		   .		   .		    .			.		    .
	2	   .		   .		   .		    .			.		    .
	2	   .		   .		   .		    .			.		    .
	.	   .		   .		   .		    .			.		    .
	.	   .		   .		   .		    .			.		    .
	.	   .		   .		   .		    .			.		    .
   nsetofimages	   .		   .		   .		    .			.		    .

   At each line we store in j the index of the set of images to store the redshifts and the number of images of each set
   At each line we add an images in the position of images' array
* 
* @param runmode runmode parameter
* @param image	array where images will be stored
* @param nImagesSet array where the number of images per Set will be stored
*/


/// read the positions, redshifts and numbers of multiple images from the images file
void module_readParameters_readImages(const struct runmode_param *runmode, struct galaxy image[], int nImagesSet[])
{

	std::string second, line1;
	int imageIndex=0;
	int j=0;
	//cast variables
	double cast_x,cast_y,cast_a,cast_b,cast_theta,cast_z,cast_mag;
	



for(int i=0; i<runmode->nsets; i++){
	nImagesSet[i]=0;
}

/*********initialisation of nimages array*********/
for(int i=0; i<runmode->nimagestot; i++){
	
	/* Initialise here the variables of the images*/
	image[i].center.x = image[i].center.y = 0;
	image[i].shape.a = image[i].shape.b = image[i].shape.theta = 0;
	image[i].redshift = 0; 

    }

//printf("imagefile :%d \n", runmode->imagefile.c_str

// Read in images
    std::ifstream IM(runmode->imagefile.c_str(),std::ios::in);  
	if ( IM )
    	{
			int i = 0;
        	while( std::getline(IM,line1) )    // Read until we reach the end
        	{
                        // Read in the parameters, * means we read in a parameter but do not store it
            
			sscanf(line1.c_str(), "%d %lf %lf %lf %lf %lf %lf %lf", &j, &cast_x, &cast_y, &cast_a, &cast_b, &cast_theta, &cast_z, &cast_mag);
			//Casting
			image[i].center.x =(type_t)cast_x;
			image[i].center.y =(type_t)cast_y;
			image[i].shape.a =(type_t)cast_a;
			image[i].shape.b =(type_t)cast_b;
			image[i].shape.theta =(type_t)cast_theta;
			image[i].redshift =(type_t)cast_z;
			image[i].mag =(type_t)cast_mag;
			//Variables
			nImagesSet[j-1]++;  // Increase the counter of the number of system for system with number j-1 by 1
            imageIndex++;
            i++;
		}
	}
	else{
			std::cout << "Error : file "  << runmode->imagefile <<  " not found" << std::endl;
			exit(-1);
		}
	
	IM.close();
}

/** @brief read the positions, redshifts and numbers of multiple sources from the sources file
* 
* This module function reads in the strong lensing sources
* Output: sources coordinates, sources shapes (semi-minor, semi-major of ellipse and orientation angle), source redshifts
* 	
*  the images file must contain
	 	x_center	y_center	major axis	minor axis	ellipticity angle	redshift
	1	   .		   .		   .		    .			.		    .
	2	   .		   .		   .		    .			.		    .
	3	   .		   .		   .		    .			.		    .
   
* 
* @param runmode runmode parameter
* @param source	array where sources will be stored
*/


/// read the positions, redshifts and numbers of multiple images from the images file
void module_readParameters_readSources(struct runmode_param *runmode, struct galaxy source[])
{

	std::string second, line1;
	int j=0;
	//cast variables
	double cast_x,cast_y,cast_a,cast_b,cast_theta,cast_z,cast_mag;
	
	//Calculating nset
	std::ifstream IM(runmode->sourfile.c_str(),std::ios::in);
	if ( IM ){
		int i = 0;
        while( std::getline(IM,line1) ){    // Read until we reach the end
            i++;
		}
		runmode->nsets = i ;
	}
	else{
			std::cout << "Error : file "  << runmode->sourfile <<  " not found" << std::endl;
			exit(-1);
		}
		
	IM.close();

/*********initialisation of nimages array*********/
for(int i=0; i<runmode->nsets; i++){
	
	/* Initialise here the variables of the images*/
	source[i].center.x = source[i].center.y = 0;
	source[i].shape.a = source[i].shape.b = source[i].shape.theta = 0;
	source[i].redshift = 0; 
	source[i].mag = 0; 

    }

	std::ifstream IM2(runmode->sourfile.c_str(),std::ios::in);
// Read in images
	if ( IM2 )
    	{
			int i = 0;
        	while( std::getline(IM2,line1) )    // Read until we reach the end
        	{
                        // Read in the parameters, * means we read in a parameter but do not store it
            
			sscanf(line1.c_str(), "%d %lf %lf %lf %lf %lf %lf %lf", &j, &cast_x, &cast_y, &cast_a, &cast_b, &cast_theta, &cast_z, &cast_mag);
			//Casting
			source[i].center.x =(type_t)cast_x;
			source[i].center.y =(type_t)cast_y;
			source[i].shape.a =(type_t)cast_a;
			source[i].shape.b =(type_t)cast_b;
			source[i].shape.theta =(type_t)cast_theta;
			source[i].redshift =(type_t)cast_z;
			source[i].mag =(type_t)cast_mag;
			
            i++;
		}
	}
	else{
			std::cout << "Error : file "  << runmode->sourfile <<  " not found" << std::endl;
			exit(-1);
		}
	
	IM2.close();
}

/** @brief This module function reads the critical and caustic line information
 *@param infile path to file
* @param cline cline parameter variable
*/



#if 0
void module_readParameters_CriticCaustic(std::string infile, cline_param *cline){
	
    std::string first, second, third, line1, line2;
    //cast variables
    double cast_1;
    //Default value initialiasation
    cline->nplan = 0;
    for(int i =0; i < NPZMAX; ++i){
		cline->cz[i] = 0;
	    cline->dos[i] = 0;
	    cline->dls[i] = 0;
	    cline->dlsds[i] = 0;
	}

    cline->limitLow = 1;
    cline->dmax = 1;
    cline->limitHigh = 10;
    cline->nbgridcells = 1000;

	std::ifstream IN(infile.c_str(), std::ios::in);
    if ( IN ){
		while(std::getline(IN,line1)){
			std::istringstream read1(line1); // create a stream for the line
			read1 >> first;
			if ( strncmp(first.c_str(), "cline", 5) == 0){    // Read in runmode information
	            std::getline(IN,line2);
				std::istringstream read2(line2);
				read2 >> second >> third;    // Read in 4 words
	    		while (strncmp(second.c_str(), "end", 3))    // Read until we reach end
	    		{	
	        		if ( !strcmp(second.c_str(), "nplan") ){
						sscanf(third.c_str(), "%d", &cline->nplan);
				            if ( cline->nplan > NPZMAX ){
				                cline->nplan = NPZMAX;
							}
							int j = 0;
				            while ( read2 >> third )
				            {
				                sscanf(third.c_str(), "%lf", &cast_1);
				                cline->cz[j] =(type_t)cast_1;
				                //printf(" zf %f \n",cline->cz[j]);
				                j++;
				            }
						}
	            		
					if ( !strcmp(second.c_str(), "dmax") )
	        		{
	            		sscanf(third.c_str(), "%lf", &cast_1);
	            		cline->dmax =(type_t)cast_1;
	            		cline->xmax = cline->dmax;
	            		cline->xmin = -cline->dmax;
	            		cline->ymax = cline->dmax;
	            		cline->ymin = -cline->dmax;	            		
	        		}
	        		if ( !strcmp(second.c_str(), "pas") || !strcmp(second.c_str(), "step") || !strcmp(second.c_str(), "limitLow") )
	        		{
						sscanf(third.c_str(), "%lf", &cast_1);
						cline->limitLow =(type_t)cast_1;
	        		}
	        		if ( !strcmp(second.c_str(), "limitHigh") )
	        		{
						sscanf(third.c_str(), "%lf", &cast_1);
						cline->limitHigh =(type_t)cast_1;
	        		}
	        		if ( !strcmp(second.c_str(), "nbgridcells") )
	        		{
						sscanf(third.c_str(), "%d", &cline->nbgridcells);
	        		}
					// Read the next line
					std::getline(IN,line2);
					std::istringstream read2(line2);
					read2 >> second >> third;
	    		}
		    }  
		} // closes while loop    

}  // closes if(IN)



IN.close();
	
}
#endif

/** @brief This module function reads the potfile information
 *@param infile path to file
* @param potfile_param potfile parameter variable
*/

void setpotfiles_param(potfile_param *potfile, cosmo_param cosmology){


}

void module_readParameters_readpotfiles_param(std::string infile, potfile_param *potfile, cosmo_param cosmology){

    std::string first, second, third, line1, line2;
    //cast variables
    double cast_1, cast_2;
    //Default value potfile initialiasation

    potfile->potid = 0;
    potfile->ftype = 1;
    potfile->type = 0;
    potfile->zlens = 0;
    potfile->mag0 = 0;
    potfile->sigma = 0;
    potfile->core = -1.0;
    potfile->corekpc = -1;
    potfile->ircut = 0;
    potfile->cut1 = 0;
    potfile->cut2 = 0;
    potfile->cutkpc1 = 0;
    potfile->cutkpc2 = 0;
    potfile->islope = 0;
    potfile->slope1 = 0;
    potfile->slope2 = 0;
    potfile->npotfile = 0;

	std::ifstream IN(infile.c_str(), std::ios::in);
    if ( IN ){
		while(std::getline(IN,line1)){
			std::istringstream read1(line1); // create a stream for the line
			read1 >> first;
			if ( strncmp(first.c_str(), "potfile", 7) == 0){    // Read in potfile information
	            std::getline(IN,line2);
				std::istringstream read2(line2);
				read2 >> second >> third;    // Read in 4 words
	    		while (strncmp(second.c_str(), "end", 3))    // Read until we reach end
	    		{

					if ( !strcmp(second.c_str(), "filein") )
	        		{
						sscanf(line2.c_str(), " %*s %d %s ", &potfile->ftype, potfile->potfile);
	        		}
					else if ( !strcmp(second.c_str(), "type") )
	        		{
						sscanf(third.c_str(), "%d", &potfile->type);
	        		}
					else if ( !strcmp(second.c_str(), "zlens") || !strcmp(second.c_str(), "z_lens"))
	        		{
						sscanf(third.c_str(), "%lf", &cast_1);
						potfile->zlens =(type_t)cast_1;
	        		}

	        		else if (!strcmp(second.c_str(), "mag0") || !strcmp(second.c_str(), "r200"))
					{
						sscanf(third.c_str(), "%lf", &cast_1);
						potfile->mag0 =(type_t)cast_1;
					}
					else if (!strcmp(second.c_str(), "select"))
					{
						sscanf(third.c_str(), "%d", &potfile->select);
					}
					else if (!strcmp(second.c_str(), "core"))
					{
						sscanf(third.c_str(), "%lf", &cast_1);
						potfile->core =(type_t)cast_1;
					}
					else if (!strcmp(second.c_str(), "corekpc"))
					{
						sscanf(third.c_str(), "%lf", &cast_1);
						potfile->corekpc =(type_t)cast_1;
					}
					else if (!strcmp(second.c_str(), "cut"))
					{
						sscanf(line2.c_str(), " %*s %d%lf%lf", &potfile->ircut, &cast_1, &cast_2);
						potfile->cut1 =(type_t)cast_1;
						potfile->cut2 =(type_t)cast_2;
					}
					else if (!strcmp(second.c_str(), "cutkpc"))
					{
						sscanf(line2.c_str(), " %*s %d%lf%lf", &potfile->ircut, &cast_1, &cast_2);
						potfile->cutkpc1 =(type_t)cast_1;
						potfile->cutkpc2 =(type_t)cast_2;
					}
					else if (!strcmp(second.c_str(), "slope") || !strcmp(second.c_str(), "m200slope"))
					{
						sscanf(line2.c_str(), " %*s %d%lf%lf", &potfile->islope, &cast_1, &cast_2);
						potfile->slope1 =(type_t)cast_1;
						potfile->slope2 =(type_t)cast_2;
					}
					else if (!strcmp(second.c_str(), "sigma"))
					{
						sscanf(line2.c_str(), " %*s %d%lf%lf", &potfile->isigma, &cast_1, &cast_2);
						potfile->sigma1 =(type_t)cast_1;
						potfile->sigma2 =(type_t)cast_2;
					}
					else if (!strcmp(second.c_str(), "vdslope") || !strcmp(second.c_str(), "c200slope"))
					{
						sscanf(line2.c_str(), " %*s %d%lf%lf", &potfile->ivdslope, &cast_1, &cast_2);
						potfile->vdslope1 =(type_t)cast_1;
						potfile->vdslope2 =(type_t)cast_2;
					}
					else if (!strncmp(second.c_str(), "vdscat", 6))
					{
						sscanf(line2.c_str(), " %*s %d%lf%lf", &potfile->ivdscat, &cast_1, &cast_2);
						potfile->vdscat1 =(type_t)cast_1;
						potfile->vdscat2 =(type_t)cast_2;
					}
					else if (!strncmp(second.c_str(), "rcutscat", 8))
					{
						sscanf(line2.c_str(), " %*s %d%lf%lf", &potfile->ircutscat, &cast_1, &cast_2);
						potfile->rcutscat1 =(type_t)cast_1;
						potfile->rcutscat2 =(type_t)cast_2;
					}
					else if (!strcmp(second.c_str(), "a") || !strcmp(second.c_str(), "m200"))
					{
						sscanf(line2.c_str(), " %*s %d%lf%lf", &potfile->ia, &cast_1, &cast_2);
						potfile->a1 =(type_t)cast_1;
						potfile->a2 =(type_t)cast_2;
					}
					else if (!strcmp(second.c_str(), "b") || !strcmp(second.c_str(), "c200"))
					{
						sscanf(line2.c_str(), " %*s %d%lf%lf", &potfile->ib, &cast_1, &cast_2);
						potfile->b1 =(type_t)cast_1;
						potfile->b2 =(type_t)cast_2;
					}

					// Read the next line
					std::getline(IN,line2);
					std::istringstream read2(line2);
					read2 >> second >> third;
	    		}
		    }
		} // closes while loop

}  // closes if(IN)

IN.close();

//Calculate vdisp, rcut and rcore

//*********************************************************************
// Set the Potfile current and limiting values
//*********************************************************************
if ( potfile->ftype <= 4 )
{
    // Scale potfile SIGMA
	potfile->sigma = potfile->sigma1;

    // ... and potfile RCUT
    if ( potfile->cut1 == 0 && potfile->cutkpc1 != 0 )
    {
    	potfile->cut1 = potfile->cutkpc1 / (d0 / cosmology.h * module_cosmodistances_observerObject(potfile->zlens,cosmology));
    	potfile->cut2 = potfile->cutkpc2 / (d0 / cosmology.h * module_cosmodistances_observerObject(potfile->zlens,cosmology));
    }
    potfile->cut = potfile->cut1;

    // ... and potfile RCORE
    if ( potfile->core == -1.0 && potfile->corekpc != -1 )
    	potfile->core = potfile->corekpc / (d0 / cosmology.h * module_cosmodistances_observerObject(potfile->zlens,cosmology));

    // ... and potfile RCUT SLOPE
    potfile->slope = potfile->slope1;

    // ... and potfile VDSLOPE
    potfile->vdslope = potfile->vdslope1;
}

}

/** @brief This module function loads the potential from the potfile into a lens
 *@param infile path to file
* @param cline cline parameter variable
*/



void module_readParameters_readpotfiles(const runmode_param *runmode, potfile_param *potfile, Potential *lens){

	std::string second, line1;
	double aa,bb;
	double DTR=acos(-1.)/180.;	/* 1 deg in rad  = pi/180 */
	//cast variables
	double cast_x,cast_y,cast_theta,cast_lum,cast_mag;


// Read in potentials
    std::ifstream IM(potfile->potfile,std::ios::in);
	if ( IM )
    	{
			int i = runmode->nhalos;
        	while( std::getline(IM,line1) )    // Read until we reach the end
        	{

                // Skip commented lines with #
                if ( line1[0] == '#' )
                    continue;

                // Default initialisation of clump
                lens[i].type = potfile->type;
                lens[i].z = potfile->zlens;
                lens[i].ellipticity_potential = lens[i].ellipticity = 0.;
                lens[i].alpha  = 0.;
                lens[i].rcut = 0.;
                lens[i].rcore = 0.;
                lens[i].rscale = 0.;
                lens[i].mag = 0.;
                lens[i].lum = 0.;
                lens[i].vdisp = 0.;
                lens[i].position.x = lens[i].position.y = 0.;
                lens[i].ellipticity_angle = 0.;
                lens[i].weight = 0;
                lens[i].exponent = 0;
                lens[i].einasto_kappacritic = 0;


				// Read a line of the catalog
				if ( potfile->ftype == 1 || potfile->ftype == 3 )
				{
						sscanf( line1.c_str(), "%s%lf%lf%lf%lf%lf%lf%lf",
							&lens[i].name, &cast_x, &cast_y,
								 &aa, &bb, &cast_theta, &cast_mag, &cast_lum);
						lens[i].ellipticity = (type_t) (aa*aa-bb*bb)/(aa*aa+bb*bb);
						if ( lens[i].ellipticity < 0 )
						{
							fprintf( stderr, "ERROR: The potfile clump %s has a negative ellipticity.\n", lens[i].name );
							exit(-1);
						}
						//goto NEXT;

				}


				//Casting
				lens[i].position.x =(type_t)cast_x;
				lens[i].position.y =(type_t)cast_y;
				lens[i].theta =(type_t)cast_theta;
				lens[i].lum =(type_t)cast_lum;
				lens[i].mag =(type_t)cast_mag;
				//general parameters
				lens[i].vdisp = potfile->sigma;
				lens[i].rcore = potfile->core;
				lens[i].rcut = potfile->cut;
				lens[i].ellipticity_angle = lens[i].theta* DTR;

			    //Calculate parameters like b0, potential ellipticity and anyother parameter depending on the profile
			    module_readParameters_calculatePotentialparameter(&lens[i]);

				//Variables
				potfile->npotfile++;
				i++;
		}
	}
	else{
			std::cout << "Error : file "  << potfile->potfile <<  " not found" << std::endl;
			exit(-1);
		}

	IM.close();

}

void module_readParameters_readpotfiles_SOA(const runmode_param *runmode, potfile_param *potfile, Potential_SOA *lens){

	std::string second, line1;
	double aa,bb;
	double DTR=acos(-1.)/180.;	/* 1 deg in rad  = pi/180 */
	//cast variables
	double cast_x,cast_y,cast_theta,cast_lum,cast_mag, cast_name;


// Read in potentials
    std::ifstream IM(potfile->potfile,std::ios::in);
	if ( IM )
    	{
			int i = runmode->nhalos;
        	while( std::getline(IM,line1) )    // Read until we reach the end
        	{

                // Skip commented lines with #
                if ( line1[0] == '#' )
                    continue;

                std::cerr << "Turn " << i << std::endl;

                // Default initialisation of clump
                lens->type[i] = potfile->type;
                lens->z[i] = potfile->zlens;
                lens->ellipticity_potential[i] = lens->ellipticity[i] = 0.;
                lens->alpha[i]  = 0.;
                lens->rcut[i] = 0.;
                lens->rcore[i] = 0.;
                lens->rscale[i] = 0.;
                lens->mag[i] = 0.;
                lens->lum[i] = 0.;
                lens->vdisp[i] = 0.;
                lens->position_x[i] = lens->position_y[i] = 0.;
                lens->ellipticity_angle[i] = 0.;
                lens->weight[i] = 0;
                lens->exponent[i] = 0;
                lens->einasto_kappacritic[i] = 0;

                std::cerr << "Init finished "<< std::endl;


				// Read a line of the catalog
				if ( potfile->ftype == 1 || potfile->ftype == 3 )
				{
						sscanf( line1.c_str(), "%lf%lf%lf%lf%lf%lf%lf%lf",
							&cast_name, &cast_x, &cast_y,
								 &aa, &bb, &cast_theta, &cast_mag, &cast_lum);
						lens->ellipticity[i] = (type_t) (aa*aa-bb*bb)/(aa*aa+bb*bb);
						if ( lens->ellipticity[i] < 0 )
						{
							fprintf( stderr, "ERROR: The potfile clump %d has a negative ellipticity.\n", i );
							exit(-1);
						}
						//goto NEXT;

				}


				//Casting
				lens->name[i] =(type_t)cast_name;
				lens->position_x[i] =(type_t)cast_x;
				lens->position_y[i] =(type_t)cast_y;
				lens->theta[i] =(type_t)cast_theta;
				lens->lum[i] =(type_t)cast_lum;
				lens->mag[i] =(type_t)cast_mag;
				//general parameters
				lens->vdisp[i] = potfile->sigma;
				lens->rcore[i] = potfile->core;
				lens->rcut[i] = potfile->cut;
				lens->ellipticity_angle[i] = lens->theta[i]* DTR;
				lens->anglecos[i]	     = cos(lens->ellipticity_angle[i]);
				lens->anglesin[i] 	     = sin(lens->ellipticity_angle[i]);

			    //Calculate parameters like b0, potential ellipticity and anyother parameter depending on the profile
			    module_readParameters_calculatePotentialparameter_SOA(lens,i);

				//Variables
				potfile->npotfile++;
				i++;
		}









	}
	else{
			std::cout << "Error : file "  << potfile->potfile <<  " not found" << std::endl;
			exit(-1);
		}

	IM.close();

}



/** @brief read the information about arclets
* !Not used! Will be reworked
* This module function reads in the arclet images for weak lensing
* @param Arclet file
*/


void module_readParameters_arclets(std::string arclets_filename, point arclets_position[], ellipse arclets_shape[], double arclets_redshift[])
{
	std::string second, line1;
	int j=0;
	//cast variables
	double cast_x,cast_y,cast_a,cast_b,cast_theta,cast_z;


/******************** read the arclets file *******************/
/* the arclets file must contain
			id	x_center y_center   a   b   theta   redshift
line :	1          	 .	   .        .       .   .     .        .        
	2
	3
	.
	.
	narclets
*/

    std::ifstream IM(arclets_filename.c_str(),std::ios::in);
	if ( IM )
    	{
        	while( std::getline(IM,line1))    // Read every line
        	{       // Read in parameters, * means we read the parameter but don't store it (*s)   
			sscanf(line1.c_str(), "%*s %lf %lf %lf %lf %lf %lf", &cast_x, &cast_y, &cast_a, &cast_b, &cast_theta, &cast_z);
			//Casting
			arclets_position[j].x =(type_t)cast_x;
			arclets_position[j].y =(type_t)cast_y;
			arclets_shape[j].a =(type_t)cast_a;
			arclets_shape[j].b =(type_t)cast_b;
			arclets_shape[j].theta =(type_t)cast_theta;
			arclets_redshift[j] =(type_t)cast_z;

			j++;
		}
		IM.close();
	}
	else
	{
		printf("Error : file %s not found\n",arclets_filename.c_str());
		exit(-1);
	}
}


/** @brief This module function reads in if a parameter will be optimized by the MCMC or stay fixed. 
 * 
 * This module function reads in if a parameter will be optimized by the MCMC or stay fixed. 
 * If it will be optimized, it specifies its minimum and maximum allowed values. Unless declared otherwise by the user, input values are fixed and won't be optimized.
* 
* read an infile :
	|  .
	|  .
	|limit
	|	x_center x  x  x  x  <--- these values contains : block  min  max  accuracy
	|	y_center x  x  x  x		if block=1 this is a free parameter, otherwise not
	|	   .	 .  .  .  .		min and max define the extremal value of this parameter
	|	   .	 .  .  .  .		accuracy is a criterium of convergence for the MCMC
	|	  end
	|  .
	|  .
	|limit
	|	x_center x  x  x  x
	|	y_center x  x  x  x
	|	   .	 .  .  .  .
	|	   .	 .  .  .  .
	|	  end
	|  .
	|  .
	|finish
and fills the variables with these values
* @param infile path to input file
* @param host_potentialoptimization array where limits will be stored
* @param nhalos number of mass distributions

*/

void module_readParameters_limit(std::string infile, struct potentialoptimization host_potentialoptimization[], int nhalos )
{
    std::string first, second, line1, line2;
    int i=0;
	//cast variables
	double cast_min,cast_max,cast_sigma;

 type_t DTR=acos(-1.)/180.;	/* 1 deg in rad  = pi/180 */

/*** initialize the block variables to zero (= not to be optimized) ***/
for(int index=0; index<nhalos; index++)
{
	host_potentialoptimization[index].position.x.block=0;
	host_potentialoptimization[index].position.y.block=0;
	host_potentialoptimization[index].weight.block=0;
	host_potentialoptimization[index].ellipticity_angle.block=0;
	host_potentialoptimization[index].ellipticity.block=0;
	host_potentialoptimization[index].rcore.block=0;
	host_potentialoptimization[index].rcut.block=0;
	host_potentialoptimization[index].rscale.block=0;
	host_potentialoptimization[index].exponent.block=0;
	host_potentialoptimization[index].alpha.block=0;
	host_potentialoptimization[index].einasto_kappacritic.block=0;
	host_potentialoptimization[index].z.block=0;
}




// Read in 
std::ifstream IN(infile.c_str(), std::ios::in);
    if ( IN )
    {
	while(std::getline(IN,line1))
	{   
	    	std::istringstream read1(line1);
	    	read1 >> first;


		if (!strncmp(first.c_str(), "limit", 5))  // Read the limits
		{
                while(std::getline(IN,line2))
                {
                    std::istringstream read2(line2);
                    read2 >> second;    // Read in 1 word
                    if (strcmp(second.c_str(), "end") == 0) break;  // stop reading at "end" and move to next potential limit section
                

                    if (!strcmp(second.c_str(), "x_centre") ||    // Read in for x center
	                !strcmp(second.c_str(), "x_center") )
                    {
                        sscanf(line2.c_str(), "%*s %d %lf %lf %lf", &host_potentialoptimization[i].position.x.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].position.x.min =(type_t)cast_min;
                        host_potentialoptimization[i].position.x.max =(type_t)cast_max;
                        host_potentialoptimization[i].position.x.sigma =(type_t)cast_sigma;
                    }
                    else if (!strcmp(second.c_str(), "y_centre") ||  // Read in for y center
            		 !strcmp(second.c_str(), "y_center") )
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].position.y.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].position.y.min =(type_t)cast_min;
                        host_potentialoptimization[i].position.y.max =(type_t)cast_max;
                        host_potentialoptimization[i].position.y.sigma =(type_t)cast_sigma;
                    }
                    else if ( !strcmp(second.c_str(), "v_disp") )  // Read in for ellipticity
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].vdisp.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].vdisp.min =(type_t)cast_min;
                        host_potentialoptimization[i].vdisp.max =(type_t)cast_max;
                        host_potentialoptimization[i].vdisp.sigma =(type_t)cast_sigma;
                    }
                    else if ( !strcmp(second.c_str(), "ellip_pot") )  // Read in for ellipticity
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].ellipticity_potential.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].ellipticity_potential.min =(type_t)cast_min;
                        host_potentialoptimization[i].ellipticity_potential.max =(type_t)cast_max;
                        host_potentialoptimization[i].ellipticity_potential.sigma =(type_t)cast_sigma;
                    }
                    else if ( !strcmp(second.c_str(), "ellipticitymass") || !strcmp(second.c_str(), "ellipticity") || !strcmp(second.c_str(), "ellipticite") )  // Read in for ellipticity
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].ellipticity.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].ellipticity.min =(type_t)cast_min;
                        host_potentialoptimization[i].ellipticity.max =(type_t)cast_max;
                        host_potentialoptimization[i].ellipticity.sigma =(type_t)cast_sigma;
                    }
                    else if (!strcmp(second.c_str(), "ellipticity_angle") || !strcmp(second.c_str(), "angle_pos"))  // Read in for ellipticity angle
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].ellipticity_angle.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].ellipticity_angle.min =(type_t)cast_min;
                        host_potentialoptimization[i].ellipticity_angle.max =(type_t)cast_max;
                        host_potentialoptimization[i].ellipticity_angle.sigma =(type_t)cast_sigma;

                        host_potentialoptimization[i].ellipticity_angle.min *= DTR;
                        host_potentialoptimization[i].ellipticity_angle.max *= DTR;
                        host_potentialoptimization[i].ellipticity_angle.sigma *= DTR;
                    }
                    else if ( !strcmp(second.c_str(), "rcut") || !strcmp(second.c_str(), "cut_radius"))    // Read in for r cut
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].rcut.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].rcut.min =(type_t)cast_min;
                        host_potentialoptimization[i].rcut.max =(type_t)cast_max;
                        host_potentialoptimization[i].rcut.sigma =(type_t)cast_sigma;
                    }
                    else if ( !strcmp(second.c_str(), "rcore") || !strcmp(second.c_str(), "core_radius"))  // Read in for core radius
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].rcore.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].rcore.min =(type_t)cast_min;
                        host_potentialoptimization[i].rcore.max =(type_t)cast_max;
                        host_potentialoptimization[i].rcore.sigma =(type_t)cast_sigma;
                    }
                    else if ( !strcmp(second.c_str(), "NFW_rs") ||    // Read in for NFW scale radius
                              !strcmp(second.c_str(), "rscale") )
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].rscale.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].rscale.min =(type_t)cast_min;
                        host_potentialoptimization[i].rscale.max =(type_t)cast_max;
                        host_potentialoptimization[i].rscale.sigma =(type_t)cast_sigma;
                    }
            	else if ( !strcmp(second.c_str(), "exponent") )    // Read in for exponent
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].exponent.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].exponent.min =(type_t)cast_min;
                        host_potentialoptimization[i].exponent.max =(type_t)cast_max;
                        host_potentialoptimization[i].exponent.sigma =(type_t)cast_sigma;
                    }
            	else if ( !strcmp(second.c_str(), "alpha") )    // Read in for alpha
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].alpha.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].alpha.min =(type_t)cast_min;
                        host_potentialoptimization[i].alpha.max =(type_t)cast_max;
                        host_potentialoptimization[i].alpha.sigma =(type_t)cast_sigma;
                    }
            	else if ( !strcmp(second.c_str(), "einasto_kappacritic") )  // Read in for critical kappa
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].einasto_kappacritic.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].einasto_kappacritic.min =(type_t)cast_min;
                        host_potentialoptimization[i].einasto_kappacritic.max =(type_t)cast_max;
                        host_potentialoptimization[i].einasto_kappacritic.sigma =(type_t)cast_sigma;
                    }
                    else if (!strcmp(second.c_str(), "virial_mass") ||  // Read in for virial mass
                             !strcmp(second.c_str(), "masse") ||
                             !strcmp(second.c_str(), "m200") ||
                             !strcmp(second.c_str(), "mass"))
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].weight.block,
                        		&cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].weight.min =(type_t)cast_min;
                        host_potentialoptimization[i].weight.max =(type_t)cast_max;
                        host_potentialoptimization[i].weight.sigma =(type_t)cast_sigma;
                    }
                    else if ( !strcmp(second.c_str(), "z_lens") )    // Read in for redshift
                    {
                        sscanf(line2.c_str(), "%*s%d%lf%lf%lf", &host_potentialoptimization[i].z.block,
                               &cast_min, &cast_max, &cast_sigma);
                        host_potentialoptimization[i].z.min =(type_t)cast_min;
                        host_potentialoptimization[i].z.max =(type_t)cast_max;
                        host_potentialoptimization[i].z.sigma =(type_t)cast_sigma;
                    }
                }  // end of inner while loop
                i++; // Move to next potential
    		}
	}
    }
IN.close();
}


/** @brief This module function reads in the potential form and its parameters (e.g. NFW).
* 
* @param infile path to input file
* @param lens array where mass distributions will be stored
* @param nhalos number of mass distributions
*/


void module_readParameters_Potential(std::string infile, Potential lens[], int nhalos)
{

	double DTR=acos(-1.)/180.;	/* 1 deg in rad  = pi/180 */
    Potential  *ilens;
    std::string first, second, third, line1, line2;
    int i=0;

std::ifstream IN(infile.c_str(), std::ios::in);
    if ( IN )
    {
	while(std::getline(IN,line1))
        {
	    std::istringstream read1(line1); // create a stream for the line
	    read1 >> first;
            if (!strncmp(first.c_str(), "potent", 6))  // Read in potential
		{
/***********************************************/

    ilens = &lens[i];

    ilens->position.x  		 = ilens->position.y = 0.;
    ilens->ellipticity 		 = 0;
    ilens->ellipticity_potential = 0.;
    ilens->ellipticity_angle 	 = 0.;
    ilens->rcut 		 = 0.;
    ilens->rcore 		 = 0;
    ilens->weight 		 = 0;
    ilens->rscale 		 = 0;
    ilens->exponent 		 = 0;
    ilens->alpha 		 = 0.;
    ilens->einasto_kappacritic   = 0;
    ilens->z 			 = 0;


    while(std::getline(IN,line2))
    {
    	std::istringstream read2(line2);
    	read2 >> second >> third;
    	if (strcmp(second.c_str(), "end") == 0)  // Move to next potential and initialize it
        {
   	    if ( ilens->z == 0. )  // Check if redshift from current halo was initialized
            {
                fprintf(stderr, "ERROR: No redshift defined for potential %d\n", i);
                exit(-1);
            }
      
            break; // Break while loop and move to next potential

        }
    
 	// Read in values
 	
        if ( !strcmp(second.c_str(), "profil") ||  // Get profile
             !strcmp(second.c_str(), "profile") )
        {
		if(!strcmp(third.c_str(), "PIEMD") ||
		   !strcmp(third.c_str(), "1") )
		{
			ilens->type=1;
			strcpy(ilens->type_name,"PIEMD");//ilens->type_name="PIEMD";
		}
		if(!strcmp(third.c_str(), "NFW") ||
		   !strcmp(third.c_str(), "2") )
		{
			ilens->type=2;
			strcpy(ilens->type_name,"NFW");//ilens->type_name="NFW";
		}
		if(!strcmp(third.c_str(), "SIES") ||
		   !strcmp(third.c_str(), "3") )
		{
			ilens->type=3;
			strcpy(ilens->type_name,"SIES");//ilens->type_name="SIES";
		}
		if(!strncmp(third.c_str(), "point", 5) ||
		   !strcmp(third.c_str(), "4") )
		{
			ilens->type=4;
			strcpy(ilens->type_name,"point");//ilens->type_name="point";
		}
		if(!strcmp(third.c_str(), "SIE") ||
		   !strcmp(third.c_str(), "5") )
		{
			ilens->type=5;
			strcpy(ilens->type_name,"SIE");//ilens->type_name="point";
		}
		if(!strcmp(third.c_str(), "8") )
		{
			ilens->type=8;
			strcpy(ilens->type_name,"PIEMD1");//ilens->type_name="point";
		}
		if(!strcmp(third.c_str(), "81") )
		{
			ilens->type=81;
			strcpy(ilens->type_name,"PIEMD81");//ilens->type_name="point";
		}
		
        }
        
        else if (!strcmp(second.c_str(), "name"))    // Get name of lens
        {
            sscanf(third.c_str(),"%s",ilens->name);
        }
        
        else if (!strcmp(second.c_str(), "x_centre") ||  // Get x center
		 !strcmp(second.c_str(), "x_center") )
        {
            ilens->position.x=atof(third.c_str());
            //std::cout << "PositionX : " << std::setprecision(15) << ilens->position.x << std::endl;
        }
        else if (!strcmp(second.c_str(), "y_centre") ||  // Get y center
		 !strcmp(second.c_str(), "y_center") )
        {
            ilens->position.y=atof(third.c_str());
        }
        else if ( !strcmp(second.c_str(), "ellipticitymass") || !strcmp(second.c_str(), "ellipticity") )  // Get ellipticity
        {
            ilens->ellipticity=atof(third.c_str());
            //ilens->ellipticity=ilens->ellipticity/3.;
        }
        else if (!strcmp(second.c_str(), "ellipticity_angle") || !strcmp(second.c_str(), "angle_pos"))  // Get ellipticity angle
        {
            ilens->ellipticity_angle=atof(third.c_str());
            ilens->ellipticity_angle *= DTR;
        }
        else if ( !strcmp(second.c_str(), "rcore") || !strcmp(second.c_str(), "core_radius"))  // Get core radius
        {
            ilens->rcore=atof(third.c_str());
        }
        else if (!strcmp(second.c_str(), "rcut") || !strcmp(second.c_str(), "cut_radius"))  // Get cut radius
        {
            ilens->rcut=atof(third.c_str());
        }
		else if (!strcmp(second.c_str(), "NFW_rs") ||  // Get scale radius of NFW
			 !strcmp(second.c_str(), "rscale"))
	        {
	            ilens->rscale=atof(third.c_str());
	        }
		else if (!strcmp(second.c_str(), "exponent") )  // Get exponent
	        {
	            ilens->exponent=atof(third.c_str());
	        }
		else if (!strcmp(second.c_str(), "alpha") )  // Get alpha
	        {
	            ilens->alpha=atof(third.c_str());
	        }
		else if (!strcmp(second.c_str(), "einasto_kappacritic") ||  // Get critical kappa 
			 !strcmp(second.c_str(), "kappacritic"))
        {
            ilens->einasto_kappacritic=atof(third.c_str());
        }
        else if (!strcmp(second.c_str(), "z_lens"))  // Get redshift
        {
            ilens->z=atof(third.c_str());
        }
        else if (!strcmp(second.c_str(), "v_disp"))  // Get Dispersion velocity
        {
            ilens->vdisp=atof(third.c_str());
        }
        else if ( !strncmp(second.c_str(), "virial_mass", 6) ||  // Get virial mass
                  !strcmp(second.c_str(), "masse") ||
                  !strcmp(second.c_str(), "m200") ||
                  !strcmp(second.c_str(), "mass") )
        {
            ilens->weight=atof(third.c_str());
        }


    } // closes inner while loop
    
    //Calculate parameters like b0, potential ellipticity and anyother parameter depending on the profile
    module_readParameters_calculatePotentialparameter(ilens);
  
    i++; // Set counter to next potential




}  // closes if loop

}  // closes while loop

if(i==0){
	printf("Parameter potential not found in the file %s",infile.c_str());
	exit(-1);
	}
if(i>1){
	for(int j=1; j<i; j++)
		{
		if(lens[j].z!=lens[j-1].z) printf("Warning : Some halos have different redshifts ! They should be in the same plane (same redshift) !\n");
		}
	}

}  // closes if(IN)

IN.close();

}

/** @brief Finds the amount of potential of same types. Has to be updated when introducing a new type of potential.
*
*/


void read_potentialSOA_ntypes(std::string infile, int N_type[]){

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

					if(!strcmp(third.c_str(), "NFW") ||
					   !strcmp(third.c_str(), "2") )
					{
						ind=atoi(third.c_str());
						N_type[ind] += 1;
					}
					if(!strcmp(third.c_str(), "SIES") ||
					   !strcmp(third.c_str(), "3") )
					{
						ind=atoi(third.c_str());
						N_type[ind] += 1;
					}
					if(!strncmp(third.c_str(), "point", 5) ||
					   !strcmp(third.c_str(), "4") )
					{
						ind=atoi(third.c_str());
						N_type[ind] += 1;
					}
					if(!strcmp(third.c_str(), "SIE") ||
					   !strcmp(third.c_str(), "5") )
					{
						ind=atoi(third.c_str());
						N_type[ind] += 1;
					}
					if(!strcmp(third.c_str(), "8") )	//PIEMD
					{
						ind=atoi(third.c_str());
						N_type[ind] += 1;
					}
					if(!strcmp(third.c_str(), "81") )	//PIEMD81
					{
						ind=atoi(third.c_str());
						N_type[ind] += 1;
						//std::cerr << "Type First: " << ind << std::endl;
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

void module_readParameters_PotentialSOA_direct(std::string infile, Potential_SOA *lens_SOA, int nhalos, int npotfiles, cosmo_param cosmology){


	double DTR=acos(-1.)/180.;	/* 1 deg in rad  = pi/180 */

	double core_radius_kpc = 0.;
	double cut_radius_kpc = 0.;

	int N_type[100];
	int Indice_type[100];
	int ind;
	Potential lens_temp;
	//Init of lens_SOA
	lens_SOA->name = 	new type_t[nhalos+npotfiles];
	lens_SOA->type = 	new int[nhalos+npotfiles];
	lens_SOA->position_x  = 	new type_t[nhalos+npotfiles];
	lens_SOA->position_y = 		new type_t[nhalos+npotfiles];
	lens_SOA->b0 = 	new type_t[nhalos+npotfiles];
	lens_SOA->vdisp = 	new type_t[nhalos+npotfiles];
	lens_SOA->ellipticity_angle = new type_t[nhalos+npotfiles];
	lens_SOA->ellipticity = new type_t[nhalos+npotfiles];
	lens_SOA->ellipticity_potential = new type_t[nhalos+npotfiles];
	lens_SOA->rcore = 	new type_t[nhalos+npotfiles];
	lens_SOA->rcut = 	new type_t[nhalos+npotfiles];
	lens_SOA->z = 		new type_t[nhalos+npotfiles];
	lens_SOA->anglecos = 	new type_t[nhalos+npotfiles];
	lens_SOA->anglesin = 		new type_t[nhalos+npotfiles];

	lens_SOA->mag = 	new type_t[nhalos+npotfiles];
	lens_SOA->lum = 		new type_t[nhalos+npotfiles];
	lens_SOA->weight = 	new type_t[nhalos+npotfiles];
	lens_SOA->exponent = 		new type_t[nhalos+npotfiles];
	lens_SOA->einasto_kappacritic = 		new type_t[nhalos+npotfiles];
	lens_SOA->rscale = 		new type_t[nhalos+npotfiles];
	lens_SOA->alpha = 		new type_t[nhalos+npotfiles];
	lens_SOA->theta = 		new type_t[nhalos+npotfiles];




	//Init of N_types and Indice_type (Number of lenses of a certain type)
	for(int i=0;i < 100; ++i){
		N_type[i] = 0;
		Indice_type[i] = 0;
	}



    //First sweep through the runmode file to find N_type (number of types)
    read_potentialSOA_ntypes(infile,N_type);

    //Calcuting starting points for each type in lens array
	for(int i=1;i < 100; ++i){
		Indice_type[i] = N_type[i]+Indice_type[i-1];
		//printf("%d %d \n ",N_type[i], Indice_type[i]);
	}

    std::string first, second, third, line1, line2;
    std::ifstream IN(infile.c_str(), std::ios::in);
	if(IN){
	while(std::getline(IN,line1))
        {
	    std::istringstream read1(line1); // create a stream for the line
	    read1 >> first;
            if (!strncmp(first.c_str(), "potent", 6))  // Read in potential
		{
			lens_temp.position.x = lens_temp.position.y = 0.;
			lens_temp.ellipticity = 0;
			lens_temp.ellipticity_potential = 0.;
			lens_temp.ellipticity_angle = 0.;
			lens_temp.vdisp = 0.;
			lens_temp.rcut = 0.;
			lens_temp.rcore = 0;
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
	double d1 = d0 / cosmology.h * module_cosmodistances_observerObject(lens_temp.z,cosmology);

	// Set rcore value in kpc or in arcsec.
	if ( core_radius_kpc != 0. )
		lens_temp.rcore = core_radius_kpc / d1;
	else
		core_radius_kpc = lens_temp.rcore * d1;

	// Set rcut value in kpc or in arcsec.
	if ( cut_radius_kpc != 0. )
		lens_temp.rcut = cut_radius_kpc / d1;
	else if ( lens_temp.rcut != 0. )
		cut_radius_kpc = lens_temp.rcut * d1;

    //Calculate parameters like b0, potential ellipticity and anyother parameter depending on the profile
    module_readParameters_calculatePotentialparameter(&lens_temp);

    //assign value to SOA
    //std::cerr << "Type + indice :" << lens_temp.type << Indice_type[lens_temp.type-1] << std::endl;
	if(Indice_type[lens_temp.type-1] <nhalos){

		ind = Indice_type[lens_temp.type-1];
		//std::cerr<< ind << std::endl;
		lens_SOA->type[ind] 		     = lens_temp.type;
		lens_SOA->position_x[ind]  	     = lens_temp.position.x;
		lens_SOA->position_y[ind] 	     = lens_temp.position.y;
		lens_SOA->b0[ind] 		     = 		lens_temp.b0;
		lens_SOA->vdisp[ind] 		     = 		lens_temp.vdisp;
		lens_SOA->ellipticity_angle[ind]     = lens_temp.ellipticity_angle;
		lens_SOA->ellipticity[ind]           = lens_temp.ellipticity;
		lens_SOA->ellipticity_potential[ind] = lens_temp.ellipticity_potential;
		lens_SOA->rcore[ind] 		     = lens_temp.rcore;
		lens_SOA->rcut[ind] 		     = lens_temp.rcut;
		lens_SOA->z[ind] 		     = lens_temp.z;
		lens_SOA->anglecos[ind] 	     = cos(lens_temp.ellipticity_angle);
		lens_SOA->anglesin[ind] 	     = sin(lens_temp.ellipticity_angle);

		Indice_type[lens_temp.type-1] += 1;
	}
}  // closes if loop

}  // closes while loop
	}
	IN.close();

}

/** @brief This module function converts potentials to potentieal_SOA (AOS to SOA conversion)
*
*/

void module_readParameters_PotentialSOA(std::string infile, Potential *lens, Potential_SOA *lens_SOA, int nhalos){

	lens_SOA->type = 	new int[nhalos];
	lens_SOA->position_x  = 	new type_t[nhalos];
	lens_SOA->position_y = 		new type_t[nhalos];
	lens_SOA->b0 = 	new type_t[nhalos];
	lens_SOA->ellipticity_angle = new type_t[nhalos];
	lens_SOA->ellipticity = new type_t[nhalos];
	lens_SOA->ellipticity_potential = new type_t[nhalos];
	lens_SOA->rcore = 	new type_t[nhalos];
	lens_SOA->rcut = 	new type_t[nhalos];
	lens_SOA->z = 		new type_t[nhalos];
	lens_SOA->anglecos = 	new type_t[nhalos];
	lens_SOA->anglesin = 		new type_t[nhalos];

	int N_type[100];
	int Indice_type[100];
	int ind;

	for(int i=0;i < 100; ++i){
			N_type[i] = 0;
			Indice_type[i] = 0;
	}

	for (int i = 0; i < nhalos; ++i){
		N_type[lens[i].type] += 1;
	}
	for(int i=1;i < 100; ++i){
		Indice_type[i] = N_type[i]+Indice_type[i-1];
		//printf("%d %d \n ",N_type[i], Indice_type[i]);
	}


	for (int i = 0; i < nhalos; ++i){

		if(Indice_type[lens[i].type-1] <nhalos){

			ind = Indice_type[lens[i].type-1];
			//printf ("%d ", ind);
			lens_SOA->type[ind] 		     = lens[i].type;
			lens_SOA->position_x[ind]  	     = lens[i].position.x;
			//std::cerr << lens_SOA[1].position_x[*i_point] << " " << lens[i].position.x  << std::endl;
			lens_SOA->position_y[ind] 	     = lens[i].position.y;
			lens_SOA->b0[ind] 		     = 		lens[i].b0;
			lens_SOA->ellipticity_angle[ind]     = lens[i].ellipticity_angle;
			lens_SOA->ellipticity[ind]           = lens[i].ellipticity;
			lens_SOA->ellipticity_potential[ind] = lens[i].ellipticity_potential;
			lens_SOA->rcore[ind] 		     = lens[i].rcore;
			lens_SOA->rcut[ind] 		     = lens[i].rcut;
			lens_SOA->z[ind] 		     = lens[i].z;
			lens_SOA->anglecos[ind] 	     = cos(lens[i].ellipticity_angle);
			lens_SOA->anglesin[ind] 	     = sin(lens[i].ellipticity_angle);

			Indice_type[lens[i].type-1] += 1;
		}
	}

	//printf("Bla anglecos = %f\n", lens_SOA->anglecos[0]);

}

/** @brief This module function calculates profile depended information like the impactparameter b0 and the potential ellipticity epot
 * 
 * @param lens: mass distribution for which to calculate parameters
*/

void module_readParameters_calculatePotentialparameter(Potential *lens){

	


	switch (lens->type)
    {
	
	case(5): /*Elliptical Isothermal Sphere*/
		//impact parameter b0
		lens->b0 = 4* pi_c2 * lens->vdisp * lens->vdisp ;
		//ellipticity_potential 
		lens->ellipticity_potential = lens->ellipticity/3 ;
	    break;
	    
	case(8): /* PIEMD */
		//impact parameter b0
		lens->b0 = 6.*pi_c2 * lens->vdisp * lens->vdisp;
		//ellipticity_parameter
	    if ( lens->ellipticity == 0. && lens->ellipticity_potential != 0. ){
			// emass is (a2-b2)/(a2+b2)
			lens->ellipticity = 2.*lens->ellipticity_potential / (1. + lens->ellipticity_potential * lens->ellipticity_potential);
			//printf("1 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
		}
		else if ( lens->ellipticity == 0. && lens->ellipticity_potential == 0. ){
			lens->ellipticity_potential = 0.00001;
			//printf("2 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
		}
		else{
			// epot is (a-b)/(a+b)
			lens->ellipticity_potential = (1. - sqrt(1 - lens->ellipticity * lens->ellipticity)) / lens->ellipticity;
			//printf("3 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
		}
        break;

	case(81): /* PIEMD */
		//impact parameter b0
		lens->b0 = 6.*pi_c2 * lens->vdisp * lens->vdisp;
		//ellipticity_parameter
	    if ( lens->ellipticity == 0. && lens->ellipticity_potential != 0. ){
			// emass is (a2-b2)/(a2+b2)
			lens->ellipticity = 2.*lens->ellipticity_potential / (1. + lens->ellipticity_potential * lens->ellipticity_potential);
			printf("1 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
		}
		else if ( lens->ellipticity == 0. && lens->ellipticity_potential == 0. ){
			lens->ellipticity_potential = 0.00001;
			printf("2 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
		}
		else{
			// epot is (a-b)/(a+b)
			lens->ellipticity_potential = (1. - sqrt(1 - lens->ellipticity * lens->ellipticity)) / lens->ellipticity;
			printf("3 : %f %f \n",lens->ellipticity,lens->ellipticity_potential);
		}
        break;

	default:
            printf( "ERROR: LENSPARA profil type of clump %s unknown : %d\n",lens->name, lens->type);
		break;
    };
	
}

void module_readParameters_calculatePotentialparameter_SOA(Potential_SOA *lens, int ind){




	switch (lens->type[ind])
    {

	case(5): /*Elliptical Isothermal Sphere*/
		//impact parameter b0
		lens->b0[ind] = 4* pi_c2 * lens->vdisp[ind] * lens->vdisp[ind] ;
		//ellipticity_potential
		lens->ellipticity_potential[ind] = lens->ellipticity[ind]/3 ;
	    break;

	case(8): /* PIEMD */
		//impact parameter b0
		lens->b0[ind] = 6.*pi_c2 * lens->vdisp[ind] * lens->vdisp[ind];
		//ellipticity_parameter
	    if ( lens->ellipticity[ind] == 0. && lens->ellipticity_potential[ind] != 0. ){
			// emass is (a2-b2)/(a2+b2)
			lens->ellipticity[ind] = 2.*lens->ellipticity_potential[ind] / (1. + lens->ellipticity_potential[ind] * lens->ellipticity_potential[ind]);
			printf("1 : %f %f \n",lens->ellipticity[ind],lens->ellipticity_potential[ind]);
		}
		else if ( lens->ellipticity[ind] == 0. && lens->ellipticity_potential[ind] == 0. ){
			lens->ellipticity_potential[ind] = 0.00001;
			printf("2 : %f %f \n",lens->ellipticity[ind],lens->ellipticity_potential[ind]);
		}
		else{
			// epot is (a-b)/(a+b)
			lens->ellipticity_potential[ind] = (1. - sqrt(1 - lens->ellipticity[ind] * lens->ellipticity[ind])) / lens->ellipticity[ind];
			printf("3 : %f %f \n",lens->ellipticity[ind],lens->ellipticity_potential[ind]);
		}
        break;

	case(81): /* PIEMD */
		//impact parameter b0
		lens->b0[ind] = 6.*pi_c2 * lens->vdisp[ind] * lens->vdisp[ind];
		//ellipticity_parameter
	    if ( lens->ellipticity[ind] == 0. && lens->ellipticity_potential[ind] != 0. ){
			// emass is (a2-b2)/(a2+b2)
			lens->ellipticity[ind] = 2.*lens->ellipticity_potential[ind] / (1. + lens->ellipticity_potential[ind] * lens->ellipticity_potential[ind]);
			printf("1 : %f %f \n",lens->ellipticity[ind],lens->ellipticity_potential[ind]);
		}
		else if ( lens->ellipticity[ind] == 0. && lens->ellipticity_potential[ind] == 0. ){
			lens->ellipticity_potential[ind] = 0.00001;
			printf("2 : %f %f \n",lens->ellipticity[ind],lens->ellipticity_potential[ind]);
		}
		else{
			// epot is (a-b)/(a+b)
			lens->ellipticity_potential[ind] = (1. - sqrt(1 - lens->ellipticity[ind] * lens->ellipticity[ind])) / lens->ellipticity[ind];
			printf("3 : %f %f \n",lens->ellipticity[ind],lens->ellipticity_potential[ind]);
		}
        break;

	default:
            printf( "ERROR: LENSPARA profil type of clump %f unknown : %d N %d \n",lens->name[ind], lens->type[ind], ind);
		break;
    };

}



/** @brief This module function reads the grid information
 * 
* @param infile path to input file
* @param grid array where grid information will be stored
*/

void module_readParameters_Grid(std::string infile, grid_param *grid)
{

    std::string first, second, third, line1, line2;
    double cast_1;
    double dmax ;

std::ifstream IN(infile.c_str(), std::ios::in);
    if ( IN )
    {
	while(std::getline(IN,line1))
        {
	    std::istringstream read1(line1); // create a stream for the line
	    read1 >> first;
            if (!strncmp(first.c_str(), "grille", 7) || !strncmp(first.c_str(), "frame", 5) )  // Read in potential
		{

    while(std::getline(IN,line2))
    {
    	std::istringstream read2(line2);
    	read2 >> second >> third;
    	if (strcmp(second.c_str(), "end") == 0)  // Move to next potential and initialize it
        {
            break; // Break while loop and move to next potential
        }
    
 	// Read in values
 	
        if (!strcmp(second.c_str(), "dmax"))
        {
            sscanf(third.c_str(), "%lf", &dmax);
            //printf( "\t%s\t%lf\n", second.c_str(), dmax);
            //grid->dmax = (type_t) cast_1;
            grid->xmin = -dmax;
            grid->xmax = dmax;
            grid->ymin = -dmax;
            grid->ymax = dmax;
            grid->rmax = dmax * sqrt(2.);
        }
        else if (!strcmp(second.c_str(), "xmin"))
        {
            sscanf(third.c_str(), "%lf", &cast_1);
            grid->xmin = (type_t) cast_1;
            //printf("\t%s\t%lf\n", second.c_str(), grid->xmin);
        }
        else if (!strcmp(second.c_str(), "xmax"))
        {
            sscanf(third.c_str(), "%lf", &cast_1);
            grid->xmax = (type_t) cast_1;
            //printf("\t%s\t%lf\n", second.c_str(), grid->xmax);
        }
        else if (!strcmp(second.c_str(), "ymin"))
        {
            sscanf(third.c_str(), "%lf", &cast_1);
            grid->ymin = (type_t) cast_1;
            //printf( "\t%s\t%lf\n", second.c_str(), grid->ymin);
        }
        else if (!strcmp(second.c_str(), "ymax"))
        {
            sscanf(third.c_str(), "%lf", &cast_1);
            grid->ymax = (type_t) cast_1;
            //printf( "\t%s\t%lf\n", second.c_str(), grid->ymax);
        }
}

}  // closes if loop

}  // closes while loop

}  // closes if(IN)



IN.close();

}


/* @brief Get number of sets for cleanlens mode
 * !!Not used. Will be reworked!!
 * 
 * This module function reads in the "clean lens" mode sources: These sources are read in and lensed only once assuming a fixed mass model
// Then we can see the predicted location of the images and compare with their real positions
* 
* @param clean lens file, number of clean lens images
* @return coordinates for each image, shape of each image, redshifts, number of sets, number of images per set 
*/


// Get number of sets for cleanlens mode
void module_readParameters_SingleLensingSourcesNumberSets(std::string infile, int &nsetofimages_cleanlens )
{
    std::string line1;
    int setNumber=0;
    nsetofimages_cleanlens=0;


    // Get number of sets of images
    std::ifstream IM(infile.c_str(),std::ios::in);
	if ( IM )
    	{
        	while( std::getline(IM,line1))
        	{ 
                        // Read in parameters
			sscanf(line1.c_str(), "%d", &setNumber);
			if(setNumber > nsetofimages_cleanlens) nsetofimages_cleanlens = setNumber;
		}
	}
	IM.close();
}



/** @brief Prints out cosmology
*/

void module_readParameters_debug_cosmology(int DEBUG, cosmo_param cosmology)
{

if (DEBUG == 1)  // If we are in debug mode
{
	printf("\n\n####################### Summary of Input Parameters #######################\n");
    printf("Cosmology: Cosmology model = %d, omegaM = %lf, omegaX = %lf, curvature = %lf, wX = %lf, wa = %lf, H0 =  %lf, h = %lf\n", cosmology.model, cosmology.omegaM, cosmology.omegaX, cosmology.curvature, cosmology.wX, cosmology.wa, cosmology.H0, cosmology.h);
	
}

}

/** @brief Prints out runmode
*/

void module_readParameters_debug_runmode(int DEBUG, runmode_param runmode)
{
if (DEBUG == 1)  // If we are in debug mode
{

    printf("Runmode: nhalos = %d, nsets = %d, nimagestot = %d, source = %d, image = %d, arclet  = %d, cline = %d, inverse= %d, multi= %d DEBUG = %d\n", runmode.nhalos, runmode.nsets, runmode.nimagestot, runmode.source, runmode.image, runmode.arclet, runmode.cline, runmode.inverse, runmode.multi, runmode.debug);

}

}
/** @brief Prints out images
*/
void module_readParameters_debug_image(int DEBUG, galaxy image[], int nImagesSet[], int nsets)
{
if (DEBUG == 1)  // If we are in debug mode
{
	int index = 0;
	for ( int i = 0; i < nsets; ++i){
		for( int j = 0; j < nImagesSet[i]; ++j){
			printf("Image [%d]: x = %lf, y = %lf, shape: a = %f, b = %f, theta = %lf, redshift = %lf,  nImagesSet = %d,\n",index , image[index].center.x, image[index].center.y,  image[index].shape.a, image[index].shape.b, image[index].shape.theta, image[index].redshift,  nImagesSet[i]);
			index +=1;
			//printf( "%d \n", index );
		}
	}
}

}
/** @brief Prints out source
*/
void module_readParameters_debug_source(int DEBUG, galaxy source[], int nsets)
{
if (DEBUG == 1)  // If we are in debug mode
{
	for ( int i = 0; i < nsets; ++i){
		printf("Source[%d]: x = %lf, y = %lf, shape: a = %f, b = %f, theta = %lf, redshift = %lf. \n\n", i, source[i].center.x, source[i].center.y,  source[i].shape.a, source[i].shape.b, source[i].shape.theta, source[i].redshift);
	}
}

}
/** @brief Prints out massdistributions
*/
void module_readParameters_debug_potential(int DEBUG, Potential lenses[], int nhalos)
{
if (DEBUG == 1)  // If we are in debug mode
{
	for ( int i = 0; i < nhalos; ++i){
		printf("Potential[%d]: x = %f, y = %f, vdisp = %f, type = %d, type_name = %s, name = %s,\n \t ellipticity = %f, ellipticity_pot = %f, ellipticity angle (radians) = %f, rcore = %f, rcut = %f,\n \t rscale = %f, exponent = %f, alpha = %f, einasto kappa critical = %f, z = %f\n", i,lenses[i].position.x, lenses[i].position.y, lenses[i].vdisp, lenses[i].type, lenses[i].type_name, lenses[i].name, lenses[i].ellipticity, lenses[i].ellipticity_potential, lenses[i].ellipticity_angle, lenses[i].rcore, lenses[i].rcut, lenses[i].rscale, lenses[i].exponent, lenses[i].alpha, lenses[i].einasto_kappacritic, lenses[i].z);
	}
}

}

void module_readParameters_debug_potential_SOA(int DEBUG, Potential_SOA lenses, int nhalos)
{
if (DEBUG == 1)  // If we are in debug mode
{
	for ( int i = 0; i < nhalos; ++i){
		printf("Potential[%d]: x = %f, y = %f, vdisp = %f, type = %d, \n \t ellipticity = %f, ellipticity_pot = %f, ellipticity angle (radians) = %f, rcore = %f, rcut = %f,\n \t z = %f\n, angle cos %f, sin %f \n", i,lenses.position_x[i], lenses.position_y[i], lenses.vdisp[i], lenses.type[i], lenses.ellipticity[i], lenses.ellipticity_potential[i], lenses.ellipticity_angle[i], lenses.rcore[i], lenses.rcut[i], lenses.z[i], lenses.anglecos[i], lenses.anglesin[i]);
	}
}

}

/** @brief Prints out potfile_param
*/
void module_readParameters_debug_potfileparam(int DEBUG, potfile_param *potfile)
{
if (DEBUG == 1)  // If we are in debug mode
{

	printf("Potfile: potid = %d, filename = %s, ftype = %d, type = %d, zlens = %f, mag0 = %f, sigma = %f,\n \t core = %f, corekpc = %f, ircut = %d, cut1 = %f, cut2 = %f,\n \t cutkpc1 = %f, cutkpc2 = %f, islope = %d, slope1 = %f, slope2 = %f\n",    potfile->potid,potfile->potfile,potfile->ftype,potfile->type,potfile->zlens,potfile->mag0,potfile->sigma,potfile->core,potfile->corekpc,potfile->ircut,potfile->cut1,potfile->cut2,potfile->cutkpc1,potfile->cutkpc2,potfile->islope,potfile->slope1,potfile->slope2);
}

}
/** @brief Prints out cline
*/
void module_readParameters_debug_criticcaustic(int DEBUG, cline_param cline)
{
if (DEBUG == 1)  // If we are in debug mode
{

    printf("Cline: Number of planes = %d ", cline.nplan);
    for( int i=0; i < cline.nplan; ++i){
		printf(" z[%d] = %f, ",i ,cline.cz[i]);
		}  
	printf(" dmax = %f, Low = %f, High= %f, Nb of gridcells %d \n", cline.dmax , cline.limitLow ,cline.limitHigh, cline.nbgridcells);

}

}
/** @brief Prints out limits
*/
void module_readParameters_debug_limit(int DEBUG, struct potentialoptimization host_potentialoptimization)
{
if (DEBUG == 1)  // If we are in debug mode
{

    printf("DEBUG: Center.x B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.position.x.block, host_potentialoptimization.position.x.min, host_potentialoptimization.position.x.max, host_potentialoptimization.position.x.sigma);
    printf("DEBUG: Center.y B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.position.y.block, host_potentialoptimization.position.y.min, host_potentialoptimization.position.y.max, host_potentialoptimization.position.y.sigma);
    printf("DEBUG: weight.y B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.weight.block, host_potentialoptimization.weight.min, host_potentialoptimization.weight.max, host_potentialoptimization.weight.sigma);
    printf("DEBUG: ellipticity_angle B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.ellipticity_angle.block, host_potentialoptimization.ellipticity_angle.min, host_potentialoptimization.ellipticity_angle.max, host_potentialoptimization.ellipticity_angle.sigma);
    printf("DEBUG: ellipticity B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.ellipticity.block, host_potentialoptimization.ellipticity.min, host_potentialoptimization.ellipticity.max, host_potentialoptimization.ellipticity.sigma);
    printf("DEBUG: ellipticity_potential B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.ellipticity_potential.block, host_potentialoptimization.ellipticity_potential.min, host_potentialoptimization.ellipticity_potential.max, host_potentialoptimization.ellipticity_potential.sigma);
    printf("DEBUG: rcore B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.rcore.block, host_potentialoptimization.rcore.min, host_potentialoptimization.rcore.max, host_potentialoptimization.rcore.sigma);
    printf("DEBUG: rcut B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.rcut.block, host_potentialoptimization.rcut.min, host_potentialoptimization.rcut.max, host_potentialoptimization.rcut.sigma);
    printf("DEBUG: rscale B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.rscale.block, host_potentialoptimization.rscale.min, host_potentialoptimization.rscale.max, host_potentialoptimization.rscale.sigma);
    printf("DEBUG: exponent B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.exponent.block, host_potentialoptimization.exponent.min, host_potentialoptimization.exponent.max, host_potentialoptimization.exponent.sigma);
    printf("DEBUG: vdisp B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.vdisp.block, host_potentialoptimization.vdisp.min, host_potentialoptimization.vdisp.max, host_potentialoptimization.vdisp.sigma);
    printf("DEBUG: alpha B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.alpha.block, host_potentialoptimization.alpha.min, host_potentialoptimization.alpha.max, host_potentialoptimization.alpha.sigma);
    printf("DEBUG: z B = %d, min = %f, max = %f, sigma = %f \n", host_potentialoptimization.z.block, host_potentialoptimization.z.min, host_potentialoptimization.z.max, host_potentialoptimization.z.sigma);
            
}

}






