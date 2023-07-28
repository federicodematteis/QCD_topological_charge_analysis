#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/*
#############################################
#
# Uso una funzione di hash per effettuare un 
# binnaggio del valore di Q estratto dal log 
# file.
# Associo a Q un valore intero usando una 
# normalizzazione tra Q_max e Q_min.
# Riempio il BIN dell'istogramma ed effettuo 
# il rescaling dei bin.
#
#############################################
*/




/* Funzione per invertire l'ordine degli elementi di un array*/
void mirror(int arr[], int n)
{
    for (int dim = 0, aug = n - 1; dim < aug; dim++, aug--)
    {
        int t = arr[dim];
        arr[dim] = arr[aug];
        arr[aug] = t;
    }
}

int main()
{
    double Q_max = 4.0;
    double Q_min = -4.0;
    int N_CNFG;
    int NBIN = 200;
    int N_bin = NBIN+1;
    int n_cnfg=0;
    double value_Q;
    int Q_histo[N_bin];
    int i;

    for (i = 0; i < N_bin; i++)
    {
        Q_histo[i] = 0;
    }
    int bin;

    FILE *data_Q;
    data_Q = fopen("data_Q/log_Q_n240.txt", "r");
    if (data_Q == NULL)
    {
        printf("Cannot open file for reading topological charge data.\n");
        exit(1);
    }

    FILE *hist_Q;
    hist_Q = fopen("hist_Q/hist_Q.txt", "wt");
    if (hist_Q == NULL)
    {
        printf("Cannot open file for writing histogram of the topological charge.\n");
        exit(1);
    }

    while (fscanf(data_Q, "%lf\n", &value_Q)==1)
    {   
        int BIN =  ((int) (NBIN*((Q_max-value_Q)/(double)(Q_max-Q_min))));

        if (BIN<NBIN && BIN>=0)
        {   
            Q_histo[BIN]=Q_histo[BIN]+1;

            n_cnfg = n_cnfg+1;    
        }
    }
    
    N_CNFG=n_cnfg;

    /*
    ######################################################
    # Nota: devo specchiare l'istogramma, a causa della 
    # scelta della conversione a (int) responsabile del cambio 
    # di segno di Q in -(int)Q durante il binnaggio. 
    ######################################################
    */

    int n = sizeof(Q_histo)/sizeof(Q_histo[0]);

    mirror(Q_histo,n); 

    for(bin=0; bin < NBIN+1; bin++){
       
        fprintf(hist_Q, "%lf ", Q_min + bin*((Q_max-Q_min)/(double)NBIN));
        fprintf(hist_Q, "%lf\n", ((double)NBIN/8.)*( (double)Q_histo[bin]));

    }



    fclose(data_Q);
    fclose(hist_Q);

    return 0;
}

