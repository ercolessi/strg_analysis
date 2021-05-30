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

void DrawSystContr(
    TString fWhichEstimator = "V0M",
    Double_t lLoMult = 0., 
	Double_t lHiMult = 100.,
	Double_t lLoEE = 0., 
	Double_t lHiEE = 100.){

    TString fWhichOtherEstimator = "ZDC";
    if(fWhichEstimator.Contains("ZDC")) fWhichOtherEstimator = "V0M";

    int mult[4] = 
    //{0,30,60,100};
    {0,10,30,100};
    TFile* f[3];
    TH1D* h[3];
    TFile* fMB = TFile::Open(Form("SystematicsFinalResults-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",lLoMult,lHiMult,lLoEE,lHiEE));
    TH1D* hMB = (TH1D*)fMB->Get("hSystTot");

    for (int i = 0; i<3; i++){
        double v0lo = lLoMult, v0hi = lHiMult, zdclo = lLoEE, zdchi = lHiEE;
        if (fWhichEstimator.Contains("V0M")) {
            v0lo = mult[i];
            v0hi = mult[i+1];
        }
        if (fWhichEstimator.Contains("ZDC")) {
            zdclo = mult[i];
            zdchi = mult[i+1];
        }
        f[i] = TFile::Open(Form("../CLEAN/sistematiche/CLEAN/SystematicsFinalResults-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",v0lo,v0hi,zdclo,zdchi));
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

    h[0]->SetLineColor(kRed);
    h[1]->SetLineColor(kBlue);
    h[2]->SetLineColor(kGreen);

    TLegend* l = new TLegend(0.5,0.7,0.7,0.85);
    l->SetTextSize(0.022);
    l->SetBorderSize(0);  
    l->AddEntry(h[0],Form("%i-%i %",mult[0],mult[1]),"L");
    l->AddEntry(h[1],Form("%i-%i %",mult[1],mult[2]),"L");
    l->AddEntry(h[2],Form("%i-%i %",mult[2],mult[3]),"L");
    l->AddEntry(hMB,"integr. mult.","F");
    l->SetHeader(Form("%s selection (percentiles):",fWhichEstimator.Data()));
    l->Draw("SAME");

    h[0]->Draw("SAME");
    h[1]->Draw("SAME");
    h[2]->Draw("SAME");

    TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
    xlabel-> SetNDC();
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.06);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    xlabel-> DrawLatex(0.35, 0.82, "#Xi");

    TLatex *xlabel2 = new TLatex();
    xlabel2->SetTextFont(42);
    xlabel2-> SetNDC();
    xlabel2-> SetTextColor(1);
    xlabel2-> SetTextSize(0.03);
    xlabel2-> SetTextAlign(22);
    xlabel2-> SetTextAngle(0);
    xlabel2-> DrawLatex(0.35, 0.75, Form("%s fixed [%.0f-%.0f]",fWhichOtherEstimator.Data(),lLoEE,lHiEE));

    c->SaveAs(Form("images/MultDepSystematics%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.png",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE));
  
}