/**********************************************************
 *
 * File Q_square.c
 *
 *
 **********************************************************/

#define MAIN_PROGRAM

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

void inttochar(int n, char *a)
{
        int i,j;
        char v[3];
        for(i=0; i<3; i++)
        {
                a[i] = '0';
        }
        sprintf(v,"%d",n);      /* converts the integer n into a string of chars v */
        for(i=1; i<3; i++)
        {
                if(v[i]=='\0')
                {
                        break;
                }
        }
        for(j=3-i; j<3; j++)
        {
                a[j]=v[j-3+i];
        }
        return;
}


/*------------------------------------------------------------------------*/


int main(int argc , char* argv[])
{
        int i, n, maxnum, step;
        FILE* doc;
        double Q;
        char a[3];
	char nome_file[15];
	int n_sweep;
	double Q_square;
	double sigma_Q_square;
/*	double a=0.087;
	int d=16;
	double L=((double)d*a);
	double lattice_vol=pow(L,3);
*/
        strcpy(nome_file, "log_Q_n000.txt");

        maxnum=atoi(argv[1]);
        step=atoi(argv[2]);
    
        for(n=0; n<=maxnum; n += step)
        {
		inttochar(n,a);
		
		for(i=0;i<3;i++)
                {
                        nome_file[7+i]=a[i];
                }

                doc=fopen(nome_file, "rt");
		
		Q_square=0;
		n_sweep=0;

                while(fscanf(doc, "%lf", &Q)==1)
                {
                        Q_square += Q*Q;
			n_sweep += 1;
                }
		
		Q_square/=(double)n_sweep;

                fclose(doc);
                
		doc=fopen(nome_file, "rt");

		sigma_Q_square = 0;

		for(i=0; i<n_sweep; i++)
		{
			fscanf(doc, "%lf", &Q);

			sigma_Q_square += (Q*Q - Q_square)*(Q*Q - Q_square);
		}

		sigma_Q_square /= (double)n_sweep*(double)n_sweep;
		sigma_Q_square = sqrt(sigma_Q_square);


                fclose(doc);

		printf("Wilson t = %d\t,   Q^2 = %f +/- %f \n", n, Q_square, sigma_Q_square);
        }
        
	
	exit(EXIT_SUCCESS);

 
}
