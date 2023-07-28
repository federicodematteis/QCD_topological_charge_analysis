/**********************************************************
 *
 * File hash_for_Q.c
 *
 *
 **********************************************************/

#define MAIN_PROGRAM

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

#define dscr 200
#define xmin -4.0
#define xmax 4.0

/*---------------------------------------------------------------*/

int main()
{
	int i, n;
	int index;
	double histo[dscr+1];
	FILE* doc;
	double Q;
	char a[4];
	char txt[5];
	char nome_file[15];
	char nome_histo[15];

	strcpy(nome_file, "log_Q_n"); 
	strcpy(nome_histo, "histo_n"); 
	strcpy(txt, ".txt");

	for(n=0; n<=400; n += 40)
	{ 
		sprintf(a, "%d", n);
		int counter=0;
		
		while(a[counter]!='\0')
		{
			counter++ ;
		}
		
		for(i=7; i<7+counter; i++)
		{
			nome_file[i]= a[i-7];
			nome_histo[i]= a[i-7];
		}

		for(i=7+counter; i<7+counter+4; i++)
		{
			nome_file[i]=txt[i-(7+counter)];
			nome_histo[i]=txt[i-(7+counter)];
		}

		nome_file[i]='\0';
		nome_histo[i]='\0';

		for(i=0; i<dscr; i++)
			histo[i]=0.;

			/*++++++++++ compute histo ++++++++++*/
		
		doc=fopen(nome_file, "rt");

		while(fscanf(doc, "%lf", &Q)==1)
		{ 
			index=(int)( 0.5 + (double)dscr*(Q-xmin)/(xmax-xmin) );
			
			if(index<dscr+1 && index>=0)
				histo[index]+=1;
		}
	
		fclose(doc);

		doc = fopen(nome_histo, "wt");
		
		for(i=0; i<dscr+1; i++)
		{
			fprintf(doc, "%f %f\n", xmin+(xmax-xmin)*(double)i/(double)dscr, histo[i] );
		}

		fclose(doc);
	}
	exit(EXIT_SUCCESS);	
}
