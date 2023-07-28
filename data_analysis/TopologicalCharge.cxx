//#include "trdstyle.C"
#include<cmath>
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>


void TopologicalCharge(){

    gStyle->SetOptStat(0);
    //setTDRStyle();

    std::ifstream in_file;
    std::vector<float> Chi;
    float measured;
    int counts=0;
    int nbin = 200;
    float minChi = -4.;
    float maxChi = 4.;
    float delta_value = (maxChi-minChi)/nbin;
    float chi;
    // Costruisco l'histo per la carica topologica
    TCanvas* palette  = new TCanvas ("Topological Charge", "", 800, 600);
    TH1F *Topological_charge = new TH1F("h_correction", "", nbin, minChi, maxChi);


    in_file.open("log_Q_n320.txt");
    if (!in_file) {
        std::cout << "error1" << std::endl;
        exit(1);
    }
      
    while (in_file >> measured) {
        Chi.push_back(measured);
        counts +=1;
        chi = measured;
        Topological_charge->Fill(chi);
    }

    // for (int i=0; i < Chi.size(); i++)
    // {
    //     chi = Chi[i];
    //     Topological_charge->Fill(chi); 
    // }

    // Method for noralization of histogram bin by bin.
    

    // float SUM_scale;
    // for (int i=0; i < nbin; i++)
    // {
    //     double scale = Topological_charge->GetBinContent(i)/(double)counts;
        
    //     Topological_charge->SetBinContent(i, scale);
    //     SUM_scale = SUM_scale+scale;
    // }
    // printf("%f", SUM_scale);

    Topological_charge->SetLineWidth(1.);
    Topological_charge->SetLineColor(kBlack);
    Topological_charge->SetStats(0);
    Topological_charge->SetFillColor(kCyan);
    Topological_charge->SetFillStyle(3003);

    Topological_charge->Scale(1./Topological_charge->Integral(), "width");

    Topological_charge->Draw("same");
    Topological_charge->GetXaxis()->SetTitle("\\ Q^{t_w}");
    Topological_charge->GetYaxis()->SetTitle("Counts");

    TLegend *legend = new TLegend();
    legend->AddEntry(Topological_charge, "");
    legend->SetHeader("Topological charge distribution");
    legend->Draw();    

}
