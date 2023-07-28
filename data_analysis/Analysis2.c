#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*questo for mi permette di ciclare sui nomi dei file*/


double err_prop_div(double a, double b, double sigma_a, double sigma_b){
    
    double sigma_div = sqrt( (1./b*b)*sigma_a*sigma_a+((a*a)/(b*b*b*b))*sigma_b*sigma_b );
    return sigma_div;
}

double linear_extrapolation(double x1, double x2, double y1, double y2, double x)
{
    double y = ((y2-y1)/(x2-x1))*x - (x1/(x2-x1))*(y2-y1) + y1;
    return y;
}

double linear_extrapolation_error(double x1, double x2, double y1, double y2, double x, double errX, double errY1, double  errY2)
{
    double a,b,c,d,e;
    /*prepare errors*/
    a = sqrt(errY1*errY1+errY2*errY2);
    b = sqrt(x*x*a*a+(y2-y1)*(y2-y1)*errX*errX);
    c = (1./(x2-x1))*b;

    d = (x1/(x2-x1))*a;
    e = sqrt(c*c+d*d);
    /*compute error*/
    
    double err = sqrt(e*e+errY1*errY1);

    return err;
}

double linear_extrapolation_error_1(double x1, double x2, double y1, double y2, double x, double errX, double errY1, double  errY2)
{
    double a,b,c,d,e;
    /*prepare errors*/
    a = sqrt(errY1*errY1+errY2*errY2);
    b = x/(x2-x1)*a;

    c = (x1/(x2-x1))*a;
    d = sqrt(c*c+b*b);
    /*compute error*/
    
    double err = sqrt(d*d+errY1*errY1);

    return err;
}

double linear_extrapolation_error_2(double x1, double x2, double y1, double y2, double x, double errX, double errY1, double  errY2)
{
    double err = sqrt(((x-x1)/(x2-x1))*((x-x1)/(x2-x1))*errY1*errY1 + (1-(x-x1)/(x2-x1))*(1-(x-x1)/(x2-x1))*errY2*errY2);
    return err;
}

double linear_extrapolation_error_3(double x1, double x2, double y1, double y2, double x, double errX, double errY1, double  errY2)
{
    double sigma1 = sqrt(y1*y1*(errX/(x2-x1))*(errX/(x2-x1))+((x-x1)/(x2-x1))*((x-x1)/(x2-x1))*errY1*errY1);
    double sigma2 = sqrt(y2*y2*(errX/(x2-x1))*(errX/(x2-x1))+((x-x1)/(x2-x1))*((x-x1)/(x2-x1))*errY2*errY2);
    double err = sqrt(sigma1*sigma1+sigma2*sigma2+errY1*errY1);
    return err;
}


int main()
{
    /*
    Le informazioni che mi servono sulla simulazione sono:
    1)Il numero di step nel Wilson flow e il valore dello step.
    2)La costante adimensionale t_0/a^2 del reticolo utilizzato 
    3)Il numero di configurazioni di campo generate dalla simulazione.
    4)Il valore della scala tipica del reticolo (L/a)^4.

    Devo usare un file di setup per la simuazione per automatizzare 
    tutta l'analisi, in cui greppo dal log file tutte queste info
    per darle al C.
    */

    double VOL = (16*16*16*16)*((14./16.)*(14./16.)*(14./16.)*(14./16.));
    int step = 50;
    double t_0 = 3.7960;
    double sigma_t_0 = 0.0012;
    int Wilson_flow_max = 500;
    int N_CNFG = 4905;
    int bin = 2;

    int n_cnfg;
    int t_wilson=0;
    int T_Wilson=0; 
    double value_Q;
    char buf[16];

    double Q[11][N_CNFG];
    double Q_binned[11][N_CNFG/bin];

    while(t_wilson <= Wilson_flow_max)
    {  
        
        snprintf(buf, 16, "log_Q_n%d.txt", t_wilson); 
        FILE *data_Q = fopen(buf, "r");

        if (data_Q == NULL) {
            printf("error in buffer while processing files: %s\n", buf);
            return 1;
        }
        /*n_cnfg<np.min(N_CNFG, MAX_CNFG)*/
        for(n_cnfg=0; n_cnfg<N_CNFG; n_cnfg++)
        {
            fscanf(data_Q, "%lf\n", &value_Q);
            Q[T_Wilson][n_cnfg]=value_Q;

            if(n_cnfg%bin==0){
                Q_binned[T_Wilson][n_cnfg/bin] = Q[T_Wilson][n_cnfg];
            }
            
            /*n_cnfg = n_cnfg + 1;*/
            /*printf("%lf\n", Q[T_Wilson][n_cnfg]);*/
            
        }
        printf("%d\n", n_cnfg);

        fclose(data_Q);

        T_Wilson = T_Wilson + 1;
        t_wilson = t_wilson + step;
           
    }

    int j;
    double Q_SQUARE_MEAN[11];
    double CHI_tsquare[11];
    
    for(j=0; j<11; j++)
    {
        double Q_square_mean=0;

        for(n_cnfg=0; n_cnfg<N_CNFG; n_cnfg++){

            Q_square_mean = Q_square_mean+Q[j][n_cnfg]*Q[j][n_cnfg];
        }

        Q_square_mean = Q_square_mean/(double)N_CNFG;
        Q_SQUARE_MEAN[j] = Q_square_mean*((14./16.)*(14./16.)*(14./16.)*(14./16.));
        CHI_tsquare[j] = (Q_SQUARE_MEAN[j]/(VOL))*( t_0*t_0 );  /*note that: t_0 = t_0/a^2*/
        
        printf("<Q^2> = %lf ", Q_SQUARE_MEAN[j]);
        printf("K(t_w)*t0^2 = %lf\n", CHI_tsquare[j]);
    }
    
    double Q_SQUARE_MEAN_binned[11];
    double CHI_tsquare_binned[11];


    for(j=0; j<11; j++)
    {
        double Q_square_mean_binned=0;

        for(n_cnfg=0; n_cnfg<N_CNFG/bin; n_cnfg++){

            Q_square_mean_binned = Q_square_mean_binned+Q_binned[j][n_cnfg]*Q_binned[j][n_cnfg];
        }

        Q_square_mean_binned = Q_square_mean_binned/(double)(N_CNFG/bin);
        Q_SQUARE_MEAN_binned[j] = Q_square_mean_binned*((14./16.)*(14./16.)*(14./16.)*(14./16.));
        CHI_tsquare_binned[j] = (Q_SQUARE_MEAN_binned[j]/(VOL))*( t_0*t_0 );
        /*printf("%lf\n", CHI_tsquare_binned[j]);
        /*printf("%lf\n", Q_square_mean);*/
    }

    /*voglio estrapolare il valore di Chi*t_0^2 al tempo di riferimento 
    t_0/a^2 =2.789 
    Per farlo eseguo una media pesata delle due misure di Chi*t_0^2 
    ai tempi di Wilson vicini a t_0;
    equilvalentemente posso trovare la retta passante per i due punti
    misurati di Chi*t_0^2 ed estrapolare la misura che mi interessa.
    In entrambi i casi devo propagare gli errori sulle singole osservabili.
    */
    
    double SIGMA_Q_SQUARE[11];

    double SIGMA_Q_SQUARE_binned[11];

    /*printf("*********************\n");
    printf("Pre binning\n");
    printf("<Q^2>    ||  sigma(Q^2)\n");*/

    for(j=0; j<11; j++)
    {
        double Q_square_sigma = 0;

        for(n_cnfg=0; n_cnfg<N_CNFG; n_cnfg++)
        {
            Q_square_sigma = Q_square_sigma +  (Q[j][n_cnfg]*Q[j][n_cnfg]-Q_SQUARE_MEAN[j])*(Q[j][n_cnfg]*Q[j][n_cnfg]-Q_SQUARE_MEAN[j]);
        }

        Q_square_sigma = Q_square_sigma/((double)N_CNFG*(double)N_CNFG);
        Q_square_sigma = sqrt(Q_square_sigma);
        SIGMA_Q_SQUARE[j] = Q_square_sigma;
       
        /*printf("%lf     ", Q_SQUARE_MEAN[j]);
        printf("%lf\n", Q_square_sigma);*/
    }

    /*printf("*********************\n");
    printf("Post binning\n");
    printf("<Q^2>    ||  sigma(Q^2)\n");*/

    for(j=0; j<11; j++)
    {
        double Q_square_sigma_binned = 0;

        for(n_cnfg=0; n_cnfg<N_CNFG/bin; n_cnfg++)
        {
            Q_square_sigma_binned = Q_square_sigma_binned +  (Q_binned[j][n_cnfg]*Q_binned[j][n_cnfg]-Q_SQUARE_MEAN_binned[j])*(Q_binned[j][n_cnfg]*Q_binned[j][n_cnfg]-Q_SQUARE_MEAN_binned[j]);
        }

        Q_square_sigma_binned = Q_square_sigma_binned/((double)(N_CNFG/bin)*(double)(N_CNFG/bin));
        Q_square_sigma_binned = sqrt(Q_square_sigma_binned);
        SIGMA_Q_SQUARE_binned[j] = Q_square_sigma_binned;

        /*printf("%lf     ", Q_SQUARE_MEAN_binned[j]);
        printf("%lf\n", Q_square_sigma_binned);*/
    }



    double sigma_t_0_square = 2*sigma_t_0*t_0;
    double SIGMA_CHI_tsquare[11];
    double SIGMA_CHI_tsquare_binned[11];
    
    for(j=0; j<11; j++){

        /*SIGMA_CHI_tsquare[j] = 
            sqrt((CHI_tsquare[j]/(t_0*t_0))*(CHI_tsquare[j]/(t_0*t_0))*(sigma_t_0_square*sigma_t_0_square)
            +t_0*t_0 * (1./(double)N_CNFG)*SIGMA_Q_SQUARE[j]*(1./(double)N_CNFG)*SIGMA_Q_SQUARE[j]);

        SIGMA_CHI_tsquare_binned[j] = 
            sqrt((CHI_tsquare_binned[j]/(t_0*t_0))*(CHI_tsquare_binned[j]/(t_0*t_0))*(sigma_t_0_square*sigma_t_0_square)
            +t_0*t_0 * (1./(double)N_CNFG)*SIGMA_Q_SQUARE_binned[j]*(1./(double)N_CNFG)*SIGMA_Q_SQUARE_binned[j]);*/

        SIGMA_CHI_tsquare[j] = 
            sqrt((Q_SQUARE_MEAN[j]/(VOL))*(Q_SQUARE_MEAN[j]/(VOL))*(sigma_t_0_square*sigma_t_0_square)
            +t_0*t_0 * (1./VOL)*SIGMA_Q_SQUARE[j]*(1./VOL)*SIGMA_Q_SQUARE[j]);

        SIGMA_CHI_tsquare_binned[j] = 
            sqrt((Q_SQUARE_MEAN_binned[j]/(VOL))*(Q_SQUARE_MEAN_binned[j]/(VOL))*(sigma_t_0_square*sigma_t_0_square)
            +t_0*t_0 * (1./VOL)*SIGMA_Q_SQUARE_binned[j]*(1./VOL)*SIGMA_Q_SQUARE_binned[j]);
        
        SIGMA_CHI_tsquare[j] = 
            SIGMA_Q_SQUARE[j]*(t_0*t_0)*(1./VOL);
        
        SIGMA_CHI_tsquare_binned[j] = 
            SIGMA_Q_SQUARE_binned[j]*(t_0*t_0)*(1./VOL);

    }

    /*ora ho il vettore contenente gli errori sulla stima di Chi*t_w^2 per ogni wilson time
    lo devo usare per estrapolare Chi*t0^2

    double weight_6 = 1/((t_0-2.4));
    double weight_7 = 1/((2.8-t_0));
    double weight_variance_6 = 1/(SIGMA_CHI_tsquare[6]*SIGMA_CHI_tsquare[6]);
    double weight_variance_7 = 1/(SIGMA_CHI_tsquare[7]*SIGMA_CHI_tsquare[7]);

    double weight_variance_6_binned = 1/(SIGMA_CHI_tsquare_binned[6]*SIGMA_CHI_tsquare_binned[6]);
    double weight_variance_7_binned = 1/(SIGMA_CHI_tsquare_binned[7]*SIGMA_CHI_tsquare_binned[7]); 

    double weight_variance_Qsquare_6 = 1/(SIGMA_Q_SQUARE[6]*SIGMA_Q_SQUARE[6]);
    double weight_variance_Qsquare_7 = 1/(SIGMA_Q_SQUARE[7]*SIGMA_Q_SQUARE[7]);

    double weight_variance_Qsquare_6_binned = 1/(SIGMA_Q_SQUARE_binned[6]*SIGMA_Q_SQUARE_binned[6]);
    double weight_variance_Qsquare_7_binned = 1/(SIGMA_Q_SQUARE_binned[7]*SIGMA_Q_SQUARE_binned[7]);

    double Q_square_t0 = ( Q_SQUARE_MEAN[6]*weight_6*weight_variance_Qsquare_6+Q_SQUARE_MEAN[7]*weight_7*weight_variance_Qsquare_7 )/(weight_6*weight_variance_Qsquare_6+weight_7*weight_variance_Qsquare_6);
    double CHI_t0_square = (CHI_tsquare[6]*weight_6*weight_variance_6+CHI_tsquare[7]*weight_7*weight_variance_7)/(weight_6*weight_variance_6+weight_7*weight_variance_7);
    
    double Q_square_t0_binned = ( Q_SQUARE_MEAN_binned[6]*weight_6*weight_variance_Qsquare_6_binned+Q_SQUARE_MEAN_binned[7]*weight_7*weight_variance_Qsquare_7_binned )/(weight_6*weight_variance_Qsquare_6_binned+weight_7*weight_variance_Qsquare_7_binned);
    double CHI_t0_square_binned = (CHI_tsquare_binned[6]*weight_6*weight_variance_6_binned+CHI_tsquare_binned[7]*weight_7*weight_variance_7_binned)/ (weight_6*weight_variance_6_binned+weight_7*weight_variance_7_binned);
    
    
    printf("****************************************************\n");
    printf("Q^2[t_w] estrapolata t=t0 ");
    printf("%.9lf\n", Q_square_t0);
    printf("Chi*(t_w)^2 estrapolata t=t0 ");
    printf("%.9lf\n", CHI_t0_square);
    printf("rieseguo l'analisi Binnando le cariche topologiche\n");
    printf("Q^2[t_w] estrapolata t=t0 ");
    printf("%.9lf\n", Q_square_t0_binned);
    printf("Chi*(t_w)^2 estrapolata t=t0 \n");
    printf("%.9lf\n", CHI_t0_square_binned);
    */

    /*****************************************************************************************/
    /*rifaccio tutta l'estrapolazione usando un fit lineare per i due dati vicini a t_0

    printf("****************************************************\n");
    printf("estrapolazione prima del binnaggio delle Q\n");
    printf("****************************************************\n");

    double Q_square_t0_ex = linear_extrapolation(2.4, 2.8, Q_SQUARE_MEAN[6], Q_SQUARE_MEAN[7], t_0);
    double sigma_Q_square_t0 = linear_extrapolation_error_3(2.4, 2.8, Q_SQUARE_MEAN[6], Q_SQUARE_MEAN[7], 2.7984, 0.0009, SIGMA_Q_SQUARE[6], SIGMA_Q_SQUARE[7]);

    printf("%.9lf ", Q_square_t0_ex);
    printf("%.9lf\n", sigma_Q_square_t0);

    double Chi_t0square_ex = linear_extrapolation(2.4, 2.8, CHI_tsquare[6], CHI_tsquare[7], t_0);
    double sigma_Chi_t0square = linear_extrapolation_error_3(2.4, 2.8, CHI_tsquare[6], CHI_tsquare[7], 2.7984, 0.0009, SIGMA_CHI_tsquare[6], SIGMA_CHI_tsquare[7]);

    printf("%.9lf ", Chi_t0square_ex);
    printf("%.9lf\n", sigma_Chi_t0square);

    printf("****************************************************\n");
    printf("estrapolazione dopo il binnaggio delle Q\n");
    printf("****************************************************\n");
   

    double Q_square_t0_ex_bin = linear_extrapolation(2.4, 2.8, Q_SQUARE_MEAN_binned[6], Q_SQUARE_MEAN_binned[7], t_0);
    printf("%.9lf\n", Q_square_t0_ex_bin);
    double Chi_t0square_ex_bin = linear_extrapolation(2.4, 2.8, CHI_tsquare_binned[6], CHI_tsquare_binned[7], t_0);
    printf("%.9lf\n", Chi_t0square_ex_bin); */

    /*double Chi_t0square_ex = linear_extrapolation(2.4, 2.8, CHI_tsquare[6], CHI_tsquare[7], t_0);
    double sigma_Chi_t0square = linear_extrapolation_error_3(2.4, 2.8, CHI_tsquare[6], CHI_tsquare[7], t_0, sigma_t_0, SIGMA_CHI_tsquare[6], SIGMA_CHI_tsquare[7]);*/


    double Q_square_mean_t0[6];
    double CHi_t0_square[6];
    double Sigma_Q_square_t0[6];
    double Sigma_Chi_t0square[6];

    double Q_square_mean_t0_binned[6];
    double CHi_t0_square_binned[6];
    double Sigma_Q_square_t0_binned[6];
    double Sigma_Chi_t0square_binned[6];

    double fraction = 0;
    int sample = 0;
    
    while(fraction <= 1.25 )
    {
        double t_low; 
        double t_up;
        int index_low;
        int index_up;

        for(j=0; j<11; j++)
        {
            if(j*step <= 100*t_0*fraction && 100*t_0*fraction < (j+1)*step)
            {
                t_low = (double)(j*step)/100.;
                t_up = (double)((j+1)*step)/100.;
                index_up = j+1;
                index_low = j;
            }

        }
        

        Q_square_mean_t0[sample]=linear_extrapolation(t_low, t_up, Q_SQUARE_MEAN[index_low], Q_SQUARE_MEAN[index_up], fraction*t_0);
        Sigma_Q_square_t0[sample] = linear_extrapolation_error_2(t_low, t_up, Q_SQUARE_MEAN[index_low], Q_SQUARE_MEAN[index_up], fraction*t_0, fraction*sigma_t_0, SIGMA_Q_SQUARE[index_low], SIGMA_Q_SQUARE[index_up]);
        CHi_t0_square[sample] = linear_extrapolation(t_low, t_up, CHI_tsquare[index_low], CHI_tsquare[index_up], fraction*t_0);
        Sigma_Chi_t0square[sample] = linear_extrapolation_error_2(t_low, t_up, CHI_tsquare[index_low], CHI_tsquare[index_up], fraction*t_0, fraction*sigma_t_0, SIGMA_CHI_tsquare[index_low], SIGMA_CHI_tsquare[index_up]);
        
        Q_square_mean_t0_binned[sample]=linear_extrapolation(t_low, t_up, Q_SQUARE_MEAN_binned[index_low], Q_SQUARE_MEAN_binned[index_up], fraction*t_0);
        Sigma_Q_square_t0_binned[sample] = linear_extrapolation_error_2(t_low, t_up, Q_SQUARE_MEAN_binned[index_low], Q_SQUARE_MEAN_binned[index_up], fraction*t_0, fraction*sigma_t_0, SIGMA_Q_SQUARE_binned[index_low], SIGMA_Q_SQUARE_binned[index_up]);
        CHi_t0_square_binned[sample] = linear_extrapolation(t_low, t_up, CHI_tsquare_binned[index_low], CHI_tsquare_binned[index_up], fraction*t_0);
        Sigma_Chi_t0square_binned[sample] = linear_extrapolation_error_2(t_low, t_up, CHI_tsquare_binned[index_low], CHI_tsquare_binned[index_up], fraction*t_0, fraction*sigma_t_0, SIGMA_CHI_tsquare_binned[index_low], SIGMA_CHI_tsquare_binned[index_up]);
        
        fraction = fraction + 0.25;
        sample = sample + 1;
        
    }

    printf("\nInterpolo le osservabili a step di 0.25 nel Wilson flow\n \n");
    printf("  <Q^2>   ||   sgm(Q^2)   ||   K*t0^2   ||  sgm(K*t0^2)  ||  Ktw/Kt0  ||  sgm(Ktw/Kt0)\n");
    printf("======================================================================================\n");
    
    for(j=0; j<6; j++){
        printf("%lf       %lf       %lf       %lf       %lf       %lf\n", Q_square_mean_t0[j], Sigma_Q_square_t0[j], CHi_t0_square[j], Sigma_Chi_t0square[j], CHi_t0_square[j]/CHi_t0_square[4], err_prop_div(CHi_t0_square[j], CHi_t0_square[4], Sigma_Chi_t0square[j], Sigma_Chi_t0square[4] ));

    }

    printf("\nInterpolazione delle variabili dopo aver binnato la matrice delle cariche topologiche\n \n");
    printf("  <Q^2>   ||   sgm(Q^2)   ||   K*t0^2   ||  sgm(K*t0^2)  ||  Ktw/Kt0  ||  sgm(Ktw/Kt0)\n");
    printf("======================================================================================\n");

    for(j=0; j<6; j++){
        printf("%lf       %lf       %lf       %lf       %lf       %lf\n", Q_square_mean_t0_binned[j], Sigma_Q_square_t0_binned[j], CHi_t0_square_binned[j], Sigma_Chi_t0square_binned[j], CHi_t0_square_binned[j]/CHi_t0_square_binned[4], err_prop_div(CHi_t0_square_binned[j], CHi_t0_square_binned[4], Sigma_Chi_t0square_binned[j], Sigma_Chi_t0square_binned[4] ));

    }





    return 0;
}