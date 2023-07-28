#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>



int main(){
    /*variables and files*/


        int N_CNFG=3981;
        double L = 12*12*12*12;
        double T = 2.7984;
        int t_wilson = 11;

        FILE *data_Q_1;
        FILE *data_Q_2;
        FILE *data_Q_3;
        FILE *data_Q_4;
        FILE *data_Q_5;
        FILE *data_Q_6;
        FILE *data_Q_7;
        FILE *data_Q_8;
        FILE *data_Q_9;
        FILE *data_Q_10;
        FILE *data_Q_11;
        
        FILE *autocorrelation_Q;
        FILE *autocorrelation_Q_binned;

        double variance_Q[t_wilson];
        double Q_square_mean[t_wilson]
        double variance_Q_independent[t_wilson];
        double Q_mean[t_wilson];

        double Q[N_CNFG][t_wilson];
        double Q_binned[N_CNFG/Dbin][t_wilson];


        /*double *Q_binned =(double*)calloc((n_cnfg/Dbin),sizeof(double));*/


        data_Q_1 = fopen("data_Q/log_Q_n000.txt", "r");
        data_Q_2 = fopen("data_Q/log_Q_n040.txt", "r");
        data_Q_3 = fopen("data_Q/log_Q_n080.txt", "r");
        data_Q_4 = fopen("data_Q/log_Q_n120.txt", "r");
        data_Q_5 = fopen("data_Q/log_Q_n160.txt", "r");
        data_Q_6 = fopen("data_Q/log_Q_n200.txt", "r");
        data_Q_7 = fopen("data_Q/log_Q_n240.txt", "r");
        data_Q_8 = fopen("data_Q/log_Q_n280.txt", "r");
        data_Q_9 = fopen("data_Q/log_Q_n320.txt", "r");
        data_Q_10 = fopen("data_Q/log_Q_n360.txt", "r");
        data_Q_11 = fopen("data_Q/log_Q_n400.txt", "r");
        
        autocorrelation_Q= fopen("autocorrelation/autocorrelation_Q.txt", "wt");
        autocorrelation_Q_binned= fopen("autocorrelation/autocorrelation_Q_binned.txt", "wt");

        if (autocorrelation_Q == NULL)
        {
            printf("Cannot open file for writing.\n");
            exit(1);
        }

        if (autocorrelation_Q_binned == NULL)
        {
            printf("Cannot open file for writing.\n");
            exit(1);
        }

        if (data_Q_1 == NULL)
        {
            printf("Cannot open file for reading.\n");
            exit(1);
        }
        
        if (data_Q_2 == NULL)
        {
            printf("Cannot open file for reading.\n");
            exit(1);
        }
        
        if (data_Q_3 == NULL)
        {
            printf("Cannot open file for reading.\n");
            exit(1);
        }
        
        if (data_Q_4 == NULL)
        {
            printf("Cannot open file for reading.\n");
            exit(1);
        }
        
        if (data_Q_5 == NULL)
        {
            printf("Cannot open file for reading.\n");
            exit(1);
        }
        
        if (data_Q_6 == NULL)
        {
            printf("Cannot open file for reading.\n");
            exit(1);
        }
        
        if (data_Q_7 == NULL)
        {
            printf("Cannot open file for reading.\n");
            exit(1);
        }
        
        if (data_Q_8 == NULL)
        {
            printf("Cannot open file for reading.\n");
            exit(1);
        }
        
        if (data_Q_9 == NULL)
        {
            printf("Cannot open file for reading.\n");
            exit(1);
        }
        
        if (data_Q_10 == NULL)
        {
            printf("Cannot open file for reading.\n");
            exit(1);
        }
        
        if (data_Q_11 == NULL)
        {
            printf("Cannot open file for reading.\n");
            exit(1);
        }

    for()
    {

        double value_Q;
        int n_cnfg=0;
        int i=0;    
        int k;
        int Dbin_max = 10;

        int j;
        int Dbin=2;
        double b=0;
        double a=0;

        
        /*uso la definizione generale di varianza per la variabile q^2*/
        while (fscanf(data_Q_1, "%lf\n", &value_Q)==1)
        {  

            Q[i][t_wilson]=value_Q;

            Q_square_mean = Q_square_mean + Q[i][t_wilson]*Q[i][t_wilson];
            Q_mean = Q_mean + Q[i];
            

            a = a + Q[i][t_wilson]*Q[i][t_wilson]*Q[i][t_wilson]*Q[i][t_wilson];
            b = b + Q[i][t_wilson]*Q[i][t_wilson];

            n_cnfg = n_cnfg+1;
            i = i+1;
        }
                while (fscanf(data_Q_1, "%lf\n", &value_Q)==1)
        {  

            Q[i][t_wilson]=value_Q;

            Q_square_mean = Q_square_mean + Q[i][t_wilson]*Q[i][t_wilson];
            Q_mean = Q_mean + Q[i];
            

            a = a + Q[i][t_wilson]*Q[i][t_wilson]*Q[i][t_wilson]*Q[i][t_wilson];
            b = b + Q[i][t_wilson]*Q[i][t_wilson];

            n_cnfg = n_cnfg+1;
            i = i+1;
        }


        
        /*calcolo la varianza di Q^2 ipotizzando le Q estratte indipendenti*/
        Q_square_mean[t_wilson]=Q_square_mean/(double)n_cnfg;

        variance_Q[t_wilson] = a/(double)n_cnfg-b*b/((double)n_cnfg*(double)n_cnfg);
        
        Q_mean[t_wilson] = Q_mean/(double)n_cnfg;

        for(j=0; j<n_cnfg; j++)
        {
            variance_Q_independent[t_wilson] = variance_Q_independent[t_wilson] 
            +(Q[j][t_wilson]*Q[j][t_wilson]-Q_square_mean)*(Q[j][t_wilson]*Q[j][t_wilson]-Q_square_mean);
        }

        variance_Q_independent[t_wilson] = variance_Q_independent[t_wilson]/((double)(n_cnfg)*(double)(n_cnfg));

        printf("Before the binning\n");
        printf("err(Q^2) more generally is %lf\n", sqrt(variance_Q[t_wilson]));
        printf("err(Q^2) for indepedent ipothesys %lf\n", sqrt(variance_Q_independent[t_wilson]));

        /*siccome le due varianze sono diverse ipotizzo che esista correlazione tra le Q a diversi sweep*/
        /*calcolo l'autocorrelazione di Q^2*/

        for(k=0; k<n_cnfg; k++){
            
            double gamma = 0;
            double gamma_0 = 0;

            for(i=0; i < n_cnfg-Dbin_max; i++)
            {
                gamma =  gamma + Q[i][t_wilson]* Q[(i+k)%n_cnfg][t_wilson];
                gamma_0 = gamma_0 + Q[i][t_wilson]*Q[i][t_wilson];
            }

            gamma = gamma/(double)(n_cnfg-Dbin_max) - (Q_mean*Q_mean);
            gamma_0 = gamma_0 / (double)(n_cnfg-Dbin_max) - (Q_mean*Q_mean);
            
            fprintf(autocorrelation_Q, "%d ", k);
            fprintf(autocorrelation_Q, "%lf\n", gamma/gamma_0);

        }

        /*
        ##########################################################
        esiste una piccola autocorrelazione, devo binnare le Q
        Noto che il tempo di autocorrelazione è ordine 1 in tau:
        Gamma(t_m)/Gamma(0)= exp(t_m/1).
        Nel binnaggio delle Q deciso di eliminare 1 sample ogni 2
        ##########################################################
        */

        for(i = 0; i<n_cnfg; i++){
            /*crea Q_binned di dimensione pari*/
            /*tengo solo le entrate di posizione pari del vettore delle Q*/

            if(i%Dbin==0){
                Q_binned[i/Dbin][t_wilson]=Q[i][t_wilson];
                /*printf("prova pari %d\n %d\n", i%2, i);*/
            }
        }
        
        printf(" dimensione di Q^2 binnato %d\n", n_cnfg/Dbin);


        /*ho binnato il vettore di Q, ora ricalcolo l'autocorrelazione
        e verifico che sia nulla*/

        for(i=0; i<n_cnfg/Dbin; i++){

            Q_mean_binned[t_wilson] = Q_mean_binned[t_wilson] + Q_binned[i][t_wilson];
        }

        Q_mean_binned[t_wilson]=Q_mean_binned[t_wilson]/(n_cnfg/Dbin);

        for(k=0; k<n_cnfg/Dbin; k++){
            
            double gamma = 0;
            double gamma_0 = 0;

            for(i=0; i < (n_cnfg/Dbin) - Dbin_max ; i++)
            {
                gamma =  gamma + Q_binned[i][t_wilson]*Q_binned[(i+k)%(n_cnfg/Dbin)][t_wilson];
                gamma_0 = gamma_0 + Q_binned[i][t_wilson]*Q_binned[i][t_wilson];
            }

            gamma = gamma/(double)((n_cnfg/Dbin)-Dbin_max) - (Q_mean_binned[t_wilson]*Q_mean_binned[t_wilson]);
            gamma_0 = gamma_0 /(double)((n_cnfg/Dbin)-Dbin_max) - (Q_mean_binned[t_wilson]*Q_mean_binned[t_wilson]);
        
            fprintf(autocorrelation_Q_binned, "%d ", k);
            fprintf(autocorrelation_Q_binned, "%lf\n", gamma/gamma_0);
        }
        
        /*ricalcolo l'errore usando la formula per variabili non necessariamente indipendenti*/

        a=b=0;
        
        printf("%d\n", n_cnfg/Dbin);
        for(i=0; i < n_cnfg/Dbin; i++)
        {
            a = a + Q_binned[i][t_wilson]*Q_binned[i][t_wilson]*Q_binned[i][t_wilson]*Q_binned[i][t_wilson];
            b = b + Q_binned[i][t_wilson]*Q_binned[i][t_wilson];
            
            /*implemento la media di Q^2 dopo il binnaggio*/
            
        }
        
        variance_Q_binned[t_wilson] = a/(n_cnfg/Dbin)- ( b/(n_cnfg/Dbin) )*(b/(n_cnfg/Dbin));
        
        
        Q_square_mean_binned[t_wilson] = b/(n_cnfg/Dbin);

        for(i=0; i < n_cnfg/Dbin ; i++)
        {
            variance_Q_binned_independent[t_wilson] = variance_Q_binned_independent[t_wilson] 
            + (Q_binned[i][t_wilson]*Q_binned[i][t_wilson]-Q_square_mean_binned[t_wilson])*(Q_binned[i][t_wilson]*Q_binned[i][t_wilson]-Q_square_mean_binned[t_wilson]);
        }

    
        variance_Q_binned_independent[t_wilson] = variance_Q_binned_independent[t_wilson]/(double)((n_cnfg/Dbin)*(n_cnfg/Dbin));

        printf("after the binning\n");
        printf("errore sulla stima di Q^2 dopo il binnaggio: %lf\n", sqrt(variance_Q_binned[t_wilson]));
        printf("errore sulla stima di Q^2 dopo il binnaggio: %lf\n", sqrt(variance_Q_binned_independent[t_wilson]));
        printf(" media di Q^2 prima del binnaggio %lf\n", Q_square_mean[t_wilson]);
        printf(" media di Q^2 dopo il binnaggio %lf\n", Q_square_mean_binned[t_wilson]);

        printf(" t_o^2*Chi =  %lf\n", (Q_square_mean[t_wilson]/(L))*(T*T));

        fclose(data_Q_1);
        fclose(data_Q_2);
        fclose(data_Q_3);
        fclose(data_Q_4);
        fclose(data_Q_5);
        fclose(data_Q_6);
        fclose(data_Q_7);
        fclose(data_Q_8);
        fclose(data_Q_9);
        fclose(data_Q_10);
        fclose(data_Q_11);

        fclose(autocorrelation_Q);
        fclose(autocorrelation_Q_binned);


    }


    /*manca solo da calcolare la media di Q^2 dopo il binnaggio*/

    return 0;
}