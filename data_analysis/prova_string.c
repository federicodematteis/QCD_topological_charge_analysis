#include <stdio.h>

int main()
{
    char* a = "a"; 
    char* extension = ".txt";
    char fileSpec[strlen(a)+strlen(extension)+1];
    FILE *out;

    snprintf( fileSpec, sizeof( fileSpec ), "%s%s", a, extension );

    out = fopen( fileSpec, "w" );
    fclose(out);
    return 0;
}