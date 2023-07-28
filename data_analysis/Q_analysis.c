#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>



int main(){
    /*variables and files*/

    double value_Q;
    int n_cnfg=0;
    double Q[981];
    double L = 14.*14.*14.*14.;
    double T = 3.7960;

    FILE *data_Q_wf0;
    FILE *data_Q_wf;
    FILE *autocorrelation_Q;
    FILE *autocorrelation_Q_binned;

    double variance_Q=0;
    double Q_square_mean=0;
    double b=0;
    double a=0;
    double variance_Q_independent=0;
    double Q_mean=0;

    int i=0;    
    int k;
    int Dbin_max = 30;

    int j;
    int Dbin=6;
    double Q_binned[n_cnfg/Dbin];
    /*double *Q_binned =(double*)calloc((n_cnfg/Dbin),sizeof(double));*/


    data_Q_wf0 = fopen("log_Q_n400.txt", "r");
    data_Q_wf= fopen("log_Q_n400.txt", "r");
    autocorrelation_Q= fopen("autocorrelation/autocorrelation_Q.txt", "wt");
    autocorrelation_Q_binned= fopen("autocorrelation/autocorrelation_Q_binned.txt", "wt");

    if (data_Q_wf0 == NULL)
    {
        printf("Cannot open file for reading.\n");
        exit(1);
    }

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

    if (data_Q_wf == NULL)
    {
        printf("Cannot open file for reading.\n");
        exit(1);
    }
    
    /*uso la definizione generale di varianza per la variabile q^2*/
    while (fscanf(data_Q_wf0, "%lf\n", &value_Q)==1 && n_cnfg<981)
    {  

        Q[i]=value_Q;

        Q_square_mean = Q_square_mean + Q[i]*Q[i];
        Q_mean = Q_mean + Q[i];
        

        a = a + Q[i]*Q[i]*Q[i]*Q[i];
        b = b + Q[i]*Q[i];

        n_cnfg = n_cnfg+1;
        i = i+1;
    }
    
    /*calcolo la varianza di Q^2 ipotizzando le Q estratte indipendenti*/
    Q_square_mean=(14./16.)*(14./16.)*(14./16.)*(14./16.)*(Q_square_mean/(double)n_cnfg);

    variance_Q = a/(double)n_cnfg-b*b/((double)n_cnfg*(double)n_cnfg);
    
    Q_mean = Q_mean/(double)n_cnfg;

    for(j=0; j<n_cnfg; j++)
    {
        variance_Q_independent = variance_Q_independent +(Q[j]*Q[j]-Q_square_mean)*(Q[j]*Q[j]-Q_square_mean);
    }

    variance_Q_independent = variance_Q_independent/((double)(n_cnfg)*(double)(n_cnfg));

    printf("Before the binning\n");
    printf("err(Q^2) more generally is %lf\n", sqrt(variance_Q));
    printf("err(Q^2) for indepedent ipothesys %lf\n", sqrt(variance_Q_independent));

    /*siccome le due varianze sono diverse ipotizzo che esista correlazione tra le Q a diversi sweep*/
    /*calcolo l'autocorrelazione di Q^2*/

    for(k=0; k<Dbin_max; k++){
        
        double gamma = 0;
        double gamma_0 = 0;

        for(i=0; i < n_cnfg-Dbin_max; i++)
        {
            gamma =  gamma + Q[i]*Q[(i+k)];
            gamma_0 = gamma_0 + Q[i]*Q[i];
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
            Q_binned[i/Dbin]=Q[i];
            /*printf("prova pari %d\n %d\n", i%2, i);*/
        }
    }
    
    printf(" dimensione di Q^2 binnato %d\n", n_cnfg/Dbin);


    /*ho binnato il vettore di Q, ora ricalcolo l'autocorrelazione
    e verifico che sia nulla*/

    double Q_mean_binned=0;

    for(i=0; i<n_cnfg/Dbin; i++){

        Q_mean_binned = Q_mean_binned + Q_binned[i];
    }

    Q_mean_binned=Q_mean_binned/(n_cnfg/Dbin);

    for(k=0; k<Dbin_max; k++){
        
        double gamma = 0;
        double gamma_0 = 0;

        for(i=0; i < (n_cnfg/Dbin) - Dbin_max ; i++)
        {
            gamma =  gamma + Q_binned[i]*Q_binned[(i+k)];
            gamma_0 = gamma_0 + Q_binned[i]*Q_binned[i];
        }

        gamma = gamma/(double)((n_cnfg/Dbin)-Dbin_max) - (Q_mean_binned*Q_mean_binned);
        gamma_0 = gamma_0 /(double)((n_cnfg/Dbin)-Dbin_max) - (Q_mean_binned*Q_mean_binned);
       
        fprintf(autocorrelation_Q_binned, "%d ", k);
        fprintf(autocorrelation_Q_binned, "%lf\n", gamma/gamma_0);
    }
    
    /*ricalcolo l'errore usando la formula per variabili non necessariamente indipendenti*/

    a=b=0;
    double variance_Q_binned_independent=0;
    printf("%d\n", n_cnfg/Dbin);
    for(i=0; i < n_cnfg/Dbin; i++)
    {
        a = a + Q_binned[i]*Q_binned[i]*Q_binned[i]*Q_binned[i];
        b = b + Q_binned[i]*Q_binned[i];
        
        /*implemento la media di Q^2 dopo il binnaggio*/
        
    }
    

    double variance_Q_binned = a/(n_cnfg/Dbin)- ( b/(n_cnfg/Dbin) )*(b/(n_cnfg/Dbin));
    
    double Q_square_mean_binned =(14./16.)*(14./16.)*(14./16.)*(14./16.)*(b/(n_cnfg/Dbin));

    for(i=0; i < n_cnfg/Dbin ; i++)
    {
        variance_Q_binned_independent = variance_Q_binned_independent + (Q_binned[i]*Q_binned[i]-Q_square_mean_binned)*(Q_binned[i]*Q_binned[i]-Q_square_mean_binned);
    }

   
    variance_Q_binned_independent = variance_Q_binned_independent/(double)((n_cnfg/Dbin)*(n_cnfg/Dbin));

    printf("after the binning\n");
    printf("errore sulla stima di Q^2 dopo il binnaggio: %lf\n", sqrt(variance_Q_binned));
    printf("errore sulla stima di Q^2 dopo il binnaggio: %lf\n", sqrt(variance_Q_binned_independent));
    printf(" media di Q^2 prima del binnaggio %lf\n", Q_square_mean);
    printf(" media di Q^2 dopo il binnaggio %lf\n", Q_square_mean_binned);

    printf(" t_o^2*Chi =  %lf\n", (Q_square_mean/(L))*(T*T));

    fclose(data_Q_wf);
    fclose(data_Q_wf0);
    fclose(autocorrelation_Q);
    fclose(autocorrelation_Q_binned);





    /*manca solo da calcolare la media di Q^2 dopo il binnaggio*/

    return 0;
}