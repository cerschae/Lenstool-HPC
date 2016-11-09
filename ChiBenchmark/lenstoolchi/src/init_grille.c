#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "fonction.h"
#include "dimension.h"
#include "structure.h"
#include "lt.h"

/****************************************************************/
/*                 Program  : grille                */
/*                 Version  : 1 mai 1992            */
/*                 Location : Obs. Toulouse         */
/*                 Auteur   : jean-paul             */
/****************************************************************/

int init_grille(char *infile, int noedit)
{
    extern struct g_mode    M;
    extern struct g_grille  G;
    extern struct g_source  S;
    extern struct g_pot     P[NPOTFILE];
    extern struct pot       lens[];
	extern struct g_dyn	    Dy;      //   TV Oct2011

    char *editor;

    long int   ils, ilens;
    long int   isrc = 0;
    int   isMsgrid = 0;  // boolean 1 defined 0 not defined (default)
    char  first[150];
    int   i;
#ifndef DEBUG
    char  command[50];
#endif

    FILE  *IN, *OUT;

    set_default();

    if ( !noedit )
    {

#ifndef DEBUG
       editor=getenv("EDITOR");
       if(editor!=NULL)
       {
              sprintf(command, "%s %s", editor, infile);
       }
       else
       {
              sprintf(command, "vi %s", infile);
       }
       system(command);
#endif

        M.verbose = 1;
    }

    OUT = fopen("para.out", "w");

    if ( OUT == NULL )
    {
        fprintf(stderr, "ERROR: Opening para.out for writting\n");
        exit(-1);
    }

    IN = fopen(infile, "r");
    if ( IN != NULL )
    {
        ilens = ils = 0;
        do
        {
            flire(IN, first);
            strtok(first, " ");
            if ( strlen(first) > 0 ) fprintf(OUT, "%s\n", first);
            fflush(OUT);

            if (!strcmp(first, "runmode"))
                r_runmode(IN, OUT);
            else if (!strcmp(first, "cleanlens"))
                r_cleanlens(IN, OUT);
            else if (!strncmp(first, "image", 5))
                r_image(IN, OUT);
            else if (!strcmp(first, "source"))
                r_source(IN, OUT);
            else if (!strcmp(first, "shapemodel"))
            {
                r_shapemodel(IN, OUT, S.ns);
                S.ns ++;
                if( M.source == 0 )
                    M.source = -1;  // to distinguish with M.source == 0 ( tirer(): random sources)
            }
            else if (!strcmp(first, "shapelimit"))
            {
                r_shapelimit(IN, OUT, isrc);
                isrc ++;
            }
            else if (!strcmp(first, "cline"))
                r_cline(IN, OUT);
            else if (!strncmp(first, "grille", 6) || !strncmp(first, "grid", 4))
                r_grille(IN, OUT);
            else if (!strncmp(first, "potent", 6))
            {
                // set the name of the clump
                lens[ilens].n[0] = '\0';
                //TODO: better set the O prefix after we completely read the .par file
                if ( ilens < G.no_lens )
                    strcpy(lens[ilens].n, "O" );

                sprintf(lens[ilens].n, "%s%ld", lens[ilens].n, ilens + 1);

                r_potentiel(IN, OUT, ilens);
                ilens++;
            }
            else if (!strncmp(first, "limit", 5))
            {
                r_limit(IN, OUT, ils);
                ils++;
            }
            else if (!strncmp(first, "potfile", 7))
            {
                // check before adding more potfile
                if( G.npot == NPOTFILE )
                {
                    fprintf(stderr, "ERROR: Increase NPOTFILE to allow more potfile. Current value NPOTFILE=%d\n", NPOTFILE);
                    exit(1);
                }

                r_potfile(IN, OUT, &P[G.npot]);
                G.nplens[G.npot] = ilens;
                ilens += set_potfile(G.nplens[G.npot], &P[G.npot]);
                ils = ilens;  // synchronize limits with potentials
                G.npot++;

            }
            else if ( !strncmp(first, "multiscale", 5))
            {
                r_msgrid(IN, OUT);
                isMsgrid = 1;
            }
			else if (!strncmp(first,"dynfile",7))              //TV Oct2011 
			{   
				r_dynfile(IN,OUT); 
			} 
            else if (!strncmp(first, "grande", 6))
                r_large(IN, OUT);
            else if (!strncmp(first, "observ", 6))
                r_observ(IN, OUT);
            else if (!strcmp(first, "champ"))
                r_frame(IN, OUT);
            else if (!strcmp(first, "vfield"))
                r_vfield(IN, OUT);
            else if (!strcmp(first, "vfieldlimit"))
                r_vfieldlimit(IN, OUT);
            else if (!strncmp(first, "cosmolog", 8))
                r_cosmologie(IN, OUT);
            else if (!strcmp(first, "cosmolimit"))
                r_cosmolimit(IN, OUT);
        }
        while (strcmp(first, "fini") != 0 );

        fclose(IN);
        fflush(OUT);

        /* Potentials organisation in the lens[] list
         *
         *    Indiv. optim. pot. | Fixed pot. |Potfile 0 | Potfile 1 | Accelerated pot. |
         *                       ^            ^          ^           ^                  ^
         *                   G.no_lens    G.nplens[0]  G.nplens[1]  G.nmsgrid         G.nlens
         *
         */

        // Create multiscale grid of potentials with the correct field
        if ( isMsgrid )
        {
            G.nmsgrid = ilens;
            multiscale_grid(&ilens);
        }

        /*
         * Treatment of erroneous statements in the .par 'grille' section
         */
        // if there is no clumps defined
        if ( ilens == 0 )
        {
            fprintf(stderr, "ERROR: no potential defined\n");
            exit(-1);
        }
        // if there is less clumps than written in the .par file
        if ( ilens < G.nlens )
            G.nlens = ilens;

        // if we want to optimize more clumps than available
        if ( G.no_lens > G.nlens )
            G.no_lens = G.nlens;

        // if the multiscale section or nmsgrid keyword are not defined
        if ( G.nmsgrid == -1 )
            G.nmsgrid = G.nlens;

        // set the last potfile
        G.nplens[G.npot] = G.nmsgrid;

        set_lens_par(OUT);

        // Commented EJ111009 because no way to optimise ellipticity
        // of potfile lenses. 
        // If vdisp/rcut have to be optimised independently, uncomment OR 
        // move potfile lens out of potfile.
        // If scaling Relation of rcut BUT not on vdisp (or reverse), 
        // implement a test in o_scale_pot for block[B0] or block[RCUT]
        //
        // if we want to optimize galaxies from the potfile
        // (after set_lens_par() to initialise correctly optimised potfile gals)
        // Uncommented EG30310 for testing
        i = 0;
        while( G.no_lens > G.nplens[i] && i < G.npot )
        {
           G.nplens[i] = G.no_lens;
           i++;
        }
        int npot_to_remove = i > 0 ? i - 1 : 0;  // npot_to_remove always >= 0

        // reduce the number of potfiles
        for( i = 0; i < npot_to_remove; i++ )
            G.nplens[i] = G.nplens[i + 1]; 

        G.npot -= npot_to_remove; 

        // Accelerated clumps cannot be individually optimised
        if ( G.no_lens > G.nmsgrid )
            G.nmsgrid = G.no_lens;

        NPRINTF(stderr, "Total number of clumps :%ld  ", G.nlens);
        NPRINTF(stderr, "(no_lens:%ld, ", G.no_lens);
        for( i = 0; i < G.npot; i++ )
            NPRINTF(stderr, "nplens[%d]:%ld, ", i, G.nplens[i]);

        NPRINTF(stderr, "nmsgrid:%ld, nlens:%ld)\n", G.nmsgrid, G.nlens);

        if( G.nlens > NLMAX ) 
        {
            fprintf(stderr, "ERROR: Only %d clumps allowed. Change dimension.h to increase this number\n", NLMAX);
            exit(1);
        }

        fclose(OUT);
    }
    else
    {
        fprintf(stderr, "ERROR: file %s not found\n", infile);
        exit(-1);
    }

    return(1);
}

