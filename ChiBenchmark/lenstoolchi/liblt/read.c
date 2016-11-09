#include<stdio.h>
#include<ctype.h>
#include<stdlib.h>
#include<string.h>

/* Read a word in fichier and return it in mot as a
 * \0 terminated string.
 *
 * A word finishes when \n, ' ', \t or # is found.
 */
void fmot(FILE *fichier, char mot[])
{
    int i = 0;
    char c;

    c = getc(fichier);
    while ( c == ' ' || c == '\t' )
    {
        c = getc(fichier);
    }

    while ( c != '\n' && c != ' ' && c != '\t' && c != '#' && c != '\r' )
    {
        if ( c != '\r' )
            mot[i++] = c ;

        c = getc(fichier);
    }

    mot[i] = '\0';
}

/* Remove the final blank or TAB characters
 */
void chop(char * str)
{
    int i = strlen(str);
    if ( i == 0 ) return;
    while ( ( str[i-1] == ' ' || str[i-1] == '\t' ) && i > 0 )
    {
        i--;
    }
    str[i] = '\0';
}

/* Convert a NUL terminated string to upper case.
 * The argument is modified.
 * Return a pointer to the string.
 */
char* upcase(char * str)
{
    int i = 0;
    while ( str[i] != '\0' )
    {
        str[i] = islower(str[i]) ? toupper(str[i]) : str[i];
        i ++;
    }
    return str;
}

/* Read a line in fichier and return it in phrase as
 * a \0 terminated string.
 *
 * A line finishes when \n is found.
 * The returned line do not include comments nor trailing blank
 * characters (' ' or '\t' ).
 */
void flire(FILE *fichier, char phrase[])
{
    int i = 0;
    char c;
   
    c = getc(fichier);
    while ( !feof(fichier) && !ferror(fichier) &&
            c != '\n' && c != '#' )
    {
        if ( c != '\r' )  // ignore criage return of DOS newline
            phrase[i++] = c;

        c = getc(fichier);
    }

    // dump the rest of the line until the \n character if necessary
    if ( c == '#' )
        while ( !feof(fichier) && !ferror(fichier) && getc(fichier) != '\n' );

    phrase[i] = '\0';

    chop(phrase);
}

/* Increment by n*8 the position in the FILE fichier*/
void ftab(FILE *fichier, int n)
{
    int i;
    for (i = 0; i < n*8; i++) getc(fichier);
}

/* Increment by n the position in the FILE fichier*/
void fblanc(FILE *fichier, int n)
{
    int i;
    for (i = 0; i < n; i++) getc(fichier);
}

/* Compares str1 to str2 character by character
 * Return
 * - 0 if str1 and str2 have the same root (eg C2a and C2b).
 * - 1 otherwise (eg C23a and C2a)
 * */
int indexCmp(const char *str1, const char *str2)
{
    char a[50], b[50];  // temporary variables

    // if <str1> and <str2> are exactly the same
    if ( !strcmp(str1, str2) )
        return(0);

    strcpy(a, str1);
    strcpy(b, str2);

    if ( strlen(a) > 1 && strlen(b) > 1 )
    {
        /*Remove the last character of <str1>*/
        a[strlen(a)-1] = '\0';

        /*Remove the last character of <str2>*/
        b[strlen(b)-1] = '\0';

        /*Compare the 2 strings*/
        if (!strcmp(a, b))
            return(0);
        else
            return(1);
    }
    else
        // <str1> or <str2> have only 1 character
        return(1);
}

/*
 * Like the Unix wc tool, wc return the number of lines in the
 * file.
 *
 * It returns -1 if an error occured during the count and 0 for
 * an empty file.
 *
 * The file pointer is reset to the beginning of the file.
 * parameters :
 * - file : a pointer to a FILE
 */
long int wc(FILE * file)
{
    long int     count = 0; // newline counter
    fpos_t  pos;

    fgetpos(file, &pos);

    rewind(file);
    while ( !ferror(file) && !feof(file))
    {
        if ( fgetc(file) == '\n' )
            count++;
    }

    if ( ferror(file) )
        count = -1;
    else
        fsetpos(file, &pos);

    return count;
}

/* Return the number of words in a string
 * Words are separated by spaces
 * str string is not modified, but must be lower than 128 characters
 */
int getWords(char *str)
{
	char *pch;
	char tmp[128];
	int count;
	
	if( strlen(str) > 128 ) 
	{
		fprintf( stderr, "ERROR (getWords) : string size greater than 128\n");
		exit(1);
	}
	strcpy(tmp, str);
	count = 0; pch = strtok(tmp, " ");
	while( pch != NULL ) { pch = strtok(NULL, " "); count++; }
	return count;
}


