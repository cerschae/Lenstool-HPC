#include<stdio.h>
#include<math.h>
#include<fonction.h>
#include<constant.h>
#include<dimension.h>
#include<structure.h>
#include<stdlib.h>
#include<string.h>

//einasto
#define PATH1 "mat_2nplus1.txt" // l'emplacement du fichier contenant les valeurs
#define PATH2 "mat_2n.txt"
#define PATH3 "mat_2nplus2.txt"

#define CMAX 20
#define LMAX 80
#define MAX_CHAR_PER_LINE 1650
#define SEP " "

extern float Tab1[LMAX][CMAX];
extern float Tab2[LMAX][CMAX];
extern float Tab3[LMAX][CMAX];

void read_table_einasto()
{

	FILE* m_File1;
	FILE* m_File2;
	FILE* m_File3;
	int i;
	int j;
	
	char szbuff1[MAX_CHAR_PER_LINE];
	char* token1;
	char szbuff2[MAX_CHAR_PER_LINE];
	char* token2;
	char szbuff3[MAX_CHAR_PER_LINE];
	char* token3;

	for(i=0;i<LMAX;i++)
	{
		for (j=0;j<CMAX;j++)
		{
			Tab1[i][j]=-1;
		}
	}

	i=0;
	m_File1=fopen(PATH1, "r");
	if (! m_File1) return;
	while (!feof(m_File1) && i<LMAX)
	{
		j=0;
		fgets(szbuff1,MAX_CHAR_PER_LINE,m_File1);
		token1=strtok(szbuff1,SEP);
		while(token1 !=NULL && j<CMAX)
		{
			Tab1[i][j]=atof(token1);
			token1=strtok(NULL,SEP);
			j++;
		}
		i++;
	}
	
	i=0;
	j=0;

	for(i=0;i<LMAX;i++)
	{
		for (j=0;j<CMAX;j++)
		{
			Tab2[i][j]=-1;
		}
	}

	i=0;
	m_File2=fopen(PATH2, "r");
	if (! m_File2) return;
	while (!feof(m_File2) && i<LMAX)
	{
		j=0;
		fgets(szbuff2,MAX_CHAR_PER_LINE,m_File2);
		token2=strtok(szbuff2,SEP);
		while(token2 !=NULL && j<CMAX)
		{
			Tab2[i][j]=atof(token2);
			token2=strtok(NULL,SEP);
			j++;
		}
		i++;
	}
	
	i=0;
	j=0;

	for(i=0;i<LMAX;i++)
	{
		for (j=0;j<CMAX;j++)
		{
			Tab3[i][j]=-1;
		}
	}

	i=0;
	m_File3=fopen(PATH3, "r");
	if (! m_File3) return;
	while (!feof(m_File3) && i<LMAX)
	{
		j=0;
		fgets(szbuff3,MAX_CHAR_PER_LINE,m_File3);
		token3=strtok(szbuff3,SEP);
		while(token3 !=NULL && j<CMAX)
		{
			Tab3[i][j]=atof(token3);
			token3=strtok(NULL,SEP);
			j++;
		}
		i++;
	}
	
	i=0;
	j=0;
	
	fclose(m_File1);
	fclose(m_File2);
	fclose(m_File3);
	
}
