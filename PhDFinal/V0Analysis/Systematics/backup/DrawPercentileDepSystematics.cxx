#include <iostream>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TPad.h"
#include "TPaveText.h"
#include <TLatex.h>

void DrawPercentileDepSystematics(
    TString fWhichEstimator = "V0M",
    TString lType = "Lambda",
    Double_t lLoMult = 0., 
	Double_t lHiMult = 100.,
	Double_t lLoEE = 0., 
	Double_t lHiEE = 100.){

    TString fWhichOtherEstimator = "SPDClusters";
    if(fWhichEstimator.Contains("SPDClusters")) fWhichOtherEstimator = "V0M";

    int mult[10] = {0,5,10,15,20,30,40,50,70,100};
    TFile* f[9];
    TH1D* h[9];
    TFile* fMB = TFile::Open(Form("SystematicsFinalResults-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root",lType.Data(),lLoMult,lHiMult,lLoEE,lHiEE));
    TH1D* hMB = (TH1D*)fMB->Get("hSystTot");

    for (int i = 0; i<9; i++){
        double spdlo = lLoMult, spdhi = lHiMult, v0lo = lLoEE, v0hi = lHiEE;
        if (fWhichEstimator.Contains("SPDClusters")) {
            spdlo = mult[i];
            spdhi = mult[i+1];
        }
        if (fWhichEstimator.Contains("V0M")) {
            v0lo = mult[i];
            v0hi = mult[i+1];
        }
        f[i] = TFile::Open(Form("SystematicsFinalResults-Lambda-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root",spdlo,spdhi,v0lo,v0hi));
        h[i] = (TH1D*)f[i]->Get("hSystTot");
        h[i]->SetTitle(Form("%i-%i ",mult[i],mult[i+1]));
    }

    TCanvas * c = new TCanvas("c", "", 1200,1200);
    hMB->SetFillColor(kOrange-5);
    c->SetFillColor(kWhite);
    c->SetLeftMargin(0.17);
    c->SetRightMargin(0.17);
    c->SetBottomMargin(0.17);

    hMB->SetTitle("");
    hMB->SetFillColor(18);
    hMB->SetMarkerStyle(20);
    hMB->SetLineStyle(0);
    hMB->GetYaxis()->SetTitle("Rel. Syst. Uncertainty");
    hMB->GetYaxis()->SetRangeUser(-0.0001,0.+0.25);
    hMB->SetMarkerStyle(0);
    hMB->Draw("E2");


    TLegend* l = new TLegend(0.5,0.7,0.7,0.85);
    l->SetTextSize(0.022);
    l->SetBorderSize(0);  
    for (int i = 0; i<9; i++){
        h[i]->SetLineColor(i+1);
        l->AddEntry(h[i],Form("%i-%i %s",mult[i],mult[i+1],"%"),"L");
    }
     l->AddEntry(hMB,"integr. mult.","F");
    l->SetHeader(Form("%s selection (percentiles):",fWhichEstimator.Data()));
    l->Draw("SAME");

    for (int i = 0; i<9; i++){
        h[i]->Draw("SAME");
        l->Draw("SAME");
    }

    TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
    xlabel-> SetNDC();
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.06);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    xlabel-> DrawLatex(0.35, 0.82, "#Lambda + #bar{#Lambda}");

    TLatex *xlabel2 = new TLatex();
    xlabel2->SetTextFont(42);
    xlabel2-> SetNDC();
    xlabel2-> SetTextColor(1);
    xlabel2-> SetTextSize(0.03);
    xlabel2-> SetTextAlign(22);
    xlabel2-> SetTextAngle(0);
    xlabel2-> DrawLatex(0.35, 0.75, Form("%s fixed [%.0f-%.0f]",fWhichOtherEstimator.Data(),lLoEE,lHiEE));

    c->SaveAs(Form("images/MultDepSystematics%sL.png",fWhichEstimator.Data()));
  
}