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

void DrawSystYield(){

    TFile* fV00100 = TFile::Open("CorrSystematics-Yield-Xi-V0M-sel-V0M_000_100_ZDC_000_100.root");
    TFile* fV0030 = TFile::Open("CorrSystematics-Yield-Xi-ZDC-sel-V0M_000_030_ZDC_000_100.root");
    TFile* fV070100 = TFile::Open("CorrSystematics-Yield-Xi-ZDC-sel-V0M_070_100_ZDC_000_100.root");
    TFile* fZDC0100 = TFile::Open("CorrSystematics-Yield-Xi-ZDC-sel-V0M_000_100_ZDC_000_100.root");
    TFile* fZDC030 = TFile::Open("CorrSystematics-Yield-Xi-V0M-sel-V0M_000_100_ZDC_000_030.root");
    TFile* fZDC70100 = TFile::Open("CorrSystematics-Yield-Xi-V0M-sel-V0M_000_100_ZDC_070_100.root");


   TH1D* hV00100 = (TH1D*)fV00100->Get("hSystTot");
    TH1D* hV0030 = (TH1D*)fV0030->Get("hSystTot");
    TH1D* hV070100 = (TH1D*)fV070100->Get("hSystTot");
    TH1D* hZDC0100 = (TH1D*)fZDC0100->Get("hSystTot");
    TH1D* hZDC030 = (TH1D*)fZDC030->Get("hSystTot");
    TH1D* hZDC70100 = (TH1D*)fZDC70100->Get("hSystTot");

    hV00100->SetLineColor(kBlack);
    hZDC0100->SetLineColor(kBlue);
    hV0030->SetLineColor(kRed+1);
    hV070100->SetLineColor(kGreen+1);
    hZDC030->SetLineColor(kMagenta+1);
    hZDC70100->SetLineColor(kAzure+8);

    hV00100->SetLineWidth(3);
    hZDC0100->SetLineWidth(3);
    hV0030->SetLineWidth(3);
    hV070100->SetLineWidth(3);
    hZDC030->SetLineWidth(3);
    hZDC70100->SetLineWidth(3);

    TCanvas * c = new TCanvas("c", "", 1500,800);
    c->SetFillColor(kWhite);
    c->SetLeftMargin(0.17);
    c->SetRightMargin(0.17);
    c->SetBottomMargin(0.17);
    c->SetGridy();
    c->SetGridx();

    hV00100->SetTitle("");
    hV00100->SetMarkerStyle(20);
    hV00100->GetYaxis()->SetTitle("Rel. Syst. Uncertainty (%)");
    hV00100->GetXaxis()->SetTitle("percentile (%)");
    hV00100->GetYaxis()->SetTitleOffset(1.5);
    hV00100->GetYaxis()->SetTitleSize(0.04);
    hV00100->GetXaxis()->SetTitleOffset(1.);
    hV00100->GetXaxis()->SetTitleSize(0.04);
    hV00100->GetYaxis()->SetRangeUser(-0.0001,0.+0.1);   
     hZDC0100->Draw("HIST ");
    hV0030->Draw("HIST SAME");
    hV070100->Draw("HIST SAME");
    TLegend* f = new TLegend(0.4,0.7,0.6,0.85);
    f->SetTextSize(0.025);
    f->SetBorderSize(0);  
    f->AddEntry(hZDC0100,"V0M fixed [0,100]%","L");
    f->AddEntry(hV0030,"V0M fixed [0,30]%","L");
    f->AddEntry(hV070100,"V0M fixed [70,100]%","L");    
   // f->SetHeader(Form("%s percentile","ZDC"));
    f->Draw("SAME");
    c->Draw();
 
   

  
   


    TCanvas * c1 = new TCanvas("c", "", 1500,800);
    c1->SetFillColor(kWhite);
    c1->SetLeftMargin(0.17);
    c1->SetRightMargin(0.17);
    c1->SetBottomMargin(0.17);
    c1->SetGridy();
    c1->SetGridx();
    hV00100->Draw("HIST");
   
        hZDC030->Draw("HIST SAME");
    hZDC70100->Draw("HIST SAME");

      TLegend* l = new TLegend(0.6,0.7,0.8,0.85);
    l->SetTextSize(0.025);
    l->SetBorderSize(0);  
    l->AddEntry(hV00100,"ZDC fixed [0,100]%","L");
    l->AddEntry(hZDC030,"ZDC fixed [0,30]%","L");
    l->AddEntry(hZDC70100,"ZDC fixed [70,100]%","L");
   // l->SetHeader(Form("%s percentile","V0M"));
    l->Draw("SAME");
    c1->Draw();


    
    /*

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
    xlabel2-> DrawLatex(0.35, 0.75, Form("%s fixed [%.0f-%.0f]",fWhichOtherEstimator.Data(),lLoEE,lHiEE));*/

    //c->SaveAs(Form("images/MultDepSystematics%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.png",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE));
  
}
