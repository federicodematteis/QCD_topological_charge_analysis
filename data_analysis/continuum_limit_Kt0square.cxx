//#include "trdstyle.C"
#include<cmath>

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>

double err_prop_div(double a, double b, double sigma_a, double sigma_b){
    
    double sigma_div = sqrt( (1./b*b)*sigma_a*sigma_a+((a*a)/(b*b*b*b))*sigma_b*sigma_b );
    return sigma_div;
}

void continuum_limit_Kt0square(){
    gStyle->SetOptStat(0);
    //setTDRStyle();
    std::ifstream in_file;

    std::vector<double> kt0square;
    std::vector<double> kt0square_err;
    std::vector<double> sommer_scale_inverse;
    std::vector<double> sommer_scale_inverse_err;   
    
    double index[3];
    double KT0SQUARE[3];
    double KT0SQUARE_ERR[3];
    double SOMMER_SCALE_INV[3];
    double SOMMER_SCALE_INV_ERR[3];
    
    int lattice = 3;

    double t_0 = 0.176*0.176;
    double sigma_t_0 = sqrt(2*0.176*0.176*0.004*0.004);

    for ( int i = 0; i < lattice; i++)
    {
        index[i]=i;
        //creo un vettore degli errori sul lattice spacing nulli 
        SOMMER_SCALE_INV_ERR[i] = 0; 
    }
    

    double x,y,z,w;

    in_file.open("Chi_t0_square_050.txt");

    if (!in_file) {
        std::cout << "error1" << std::endl;
        exit(1);
    }
  
    while (in_file >> x >> y >> z >> w) {
    	sommer_scale_inverse.push_back(x);
        sommer_scale_inverse_err.push_back(y);
        kt0square.push_back(z);
        kt0square_err.push_back(w);
        
    }

    for (int i=0; i< kt0square.size(); i++){
        std::cout << kt0square[i] << std::endl;
    }

    for (int i = 0; i < lattice; i++)

    {
        KT0SQUARE[i]=kt0square[i];
        KT0SQUARE_ERR[i]=kt0square_err[i];
        SOMMER_SCALE_INV[i]=sommer_scale_inverse[i];
        
    }

//TGraph obj 
    //TGraph *E = new TGraph(dbin, SPACING, ENERGY);
    // provo TGraphErrors
    
    TGraphErrors *CHIT0SQUARE = new TGraphErrors(lattice, SOMMER_SCALE_INV, KT0SQUARE, SOMMER_SCALE_INV_ERR, KT0SQUARE_ERR);

    TMultiGraph *mg2 = new TMultiGraph();
    CHIT0SQUARE->SetMarkerSize(2);
    CHIT0SQUARE->SetMarkerStyle(20);
    CHIT0SQUARE->SetLineColor(kRed);
    CHIT0SQUARE->SetMarkerColor(kBlack);
    CHIT0SQUARE->SetLineWidth(2);
    mg2->Add(CHIT0SQUARE);
    mg2->Draw("AP");
    TLegend *p_legend = new TLegend();
    //p_legend->AddEntry(gamma,"Autocorrelation");
    //p_legend->Draw();
    

//Now we want to perform a fit and find the parameter theta

    TF1 *lin = new TF1("lin", "[0]+[1]*x", 0, 1);
    lin->SetParameter(0,1);
    lin->SetParameter(1,1);
    CHIT0SQUARE->Fit("lin", "EX0");

    double par1=lin->GetParameter(0);
    double par2=lin->GetParameter(1);
    double K = par1/(t_0*t_0);

    double sigma_par1 = lin->GetParError(0);
    double sigma_K = err_prop_div(par1, t_0*t_0, sigma_par1, sqrt(2*t_0*t_0*sigma_t_0*sigma_t_0));
    double sigma_K_MeV = (197.3*(1./4)*sqrt(sqrt(K))*(1./K)*sigma_K);

    double sigma_K_giusti = err_prop_div(par1, 0.1108*0.1108, sigma_par1, sqrt(2*0.0017*0.0017*0.1108*0.1108) );
    
    double K_MeV = sqrt(sqrt(K))*197.3;


    double Pion_decay_const = 130.2/sqrt(2); //MeV
    double sigma_Pion_dc = (1./sqrt(2))*0.8; //MeV 0.8 is the Pion decay Constant error f_pi (not capital)
    printf("%lf\n", sigma_Pion_dc);

    double Eta_mass = 547.862; //MeV
    double sigma_Eta_mass = 0.017; //MeV
    double k0_mass = 497.611; //MeV
    double sigma_k0_mass = 0.013; //MeV
    double kpm_mass = 493.677; //MeV
    double sigma_kpm_mass = 0.016; //MeV

    double c = (6.*K_MeV*K_MeV*K_MeV*K_MeV)/(Pion_decay_const*Pion_decay_const) - Eta_mass*Eta_mass + 2*kpm_mass*kpm_mass;
    double Eta_prime_mass = sqrt(c);

    double SIGMA_K_MeV = sigma_K_MeV*sigma_K_MeV*sigma_K_MeV*sigma_K_MeV;
    printf("%lf", K_MeV*K_MeV*K_MeV*K_MeV);

    //double sigma_c = sqrt( ((6.*6.)/(Pion_decay_const*Pion_decay_const*Pion_decay_const*Pion_decay_const))*SIGMA_K_MeV*SIGMA_K_MeV 
    //+ 2*Eta_mass*Eta_mass*sigma_Eta_mass*sigma_Eta_mass + 4.*2.*kpm_mass*kpm_mass*sigma_kpm_mass*sigma_kpm_mass  );
    //double sigma_c = (6./Pion_decay_const*Pion_decay_const)*SIGMA_K_MeV;
    //double sigma_Eta_prime_mass = 0.5*(1./sqrt(c))*sigma_c;
    //double sigma_Eta_prime_mass = (3./(Pion_decay_const*Pion_decay_const*Eta_prime_mass))*SIGMA_K_MeV;
    double temp = err_prop_div(K_MeV*K_MeV*K_MeV*K_MeV, Pion_decay_const*Pion_decay_const, SIGMA_K_MeV, sqrt(2)*Pion_decay_const*sigma_Pion_dc);
    double sigma_c = sqrt( (6*temp)*(6*temp) + 2*Eta_mass*Eta_mass*sigma_Eta_mass*sigma_Eta_mass + 4.*2.*kpm_mass*kpm_mass*sigma_kpm_mass*sigma_kpm_mass);
    double sigma_Eta_prime_mass = (0.5/sqrt(c))*sigma_c;
    printf("\n \nKt0^2    || sigma(Kt0^2) || K [fm^-4] || sgm(K) [fm^-4] ||   [K MeV]^4   || sgm [K MeV]^4 ||  Kr0^4  || sgm (Kr0^4) ||   M_n'  || sgm (M_n')\n");
    printf("=================================================================================================================================================\n \n");
    printf("%lf       %lf      %lf       %lf        %lf^4      %lf^4      %lf     %lf   %lf   %lf\n \n", par1, sigma_par1 , K, sigma_K , K_MeV , sigma_K_MeV, par1/(0.1108*0.1108), sigma_K_giusti, Eta_prime_mass, sigma_Eta_prime_mass);
    
    double M_n_prime_th = 957.78; //MeV
    double sigma_Mnprime_th = 0.06; //MeV

    double k_t0_square_th = 0.000667;
    double sigma_t0_square = 0.000007;
    printf(" t_value M_n' %lf\n", (M_n_prime_th-Eta_prime_mass)/(sqrt(sigma_Eta_prime_mass*sigma_Eta_prime_mass+sigma_Mnprime_th*sigma_Mnprime_th)));
    printf(" t_value Kt0^2' %lf\n", (k_t0_square_th-par1)/(sqrt(sigma_par1*sigma_par1+sigma_t0_square*sigma_t0_square)));

}
