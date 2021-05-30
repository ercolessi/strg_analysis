#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TH1.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TAttFill.h"

#endif
using namespace std;
#include <TString.h>


void DrawContributoYieldsSyst( 
    TString fWhichEstimator = "ZDC", 
    Double_t lLoMult = 0., 
    Double_t lHiMult = 100., 
    Double_t lLoEE = 0., 
    Double_t lHiEE = 100.){

    Double_t* systcounter;
       
    if (fWhichEstimator.Contains("V0M")) {
        if (lLoEE>0) {
            double syst[3] = {6,8,10};
            systcounter = syst;
        }
        if (lHiEE<100) {
            double syst[3] = {4,7,10};
            systcounter = syst;
        }
    }
    if (fWhichEstimator.Contains("ZDC")) {
        double syst[3] = {4,7,10};
            systcounter = syst;
        
        if (lLoMult>0) {
            double syst[3] = {3,6,9};
            systcounter = syst;
        }
    }

    TFile* fileOther = TFile::Open(Form("Test2SystTestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE),"READ");
    TFile* fileExtrap = TFile::Open(Form("sistematiche/Final-Extrap-Syst-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE));
    TFile* fileR = TFile::Open(Form("sistematiche/CorrSystematics-Yield-%s-%s-sel-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE));
        
    TH1D* hOther = (TH1D*)fileOther->Get("hSystContribOther");
    TH1D* hExtrap = (TH1D*)fileExtrap->Get("hSystUncorr");
    TH1D* hR = (TH1D*)fileR->Get("hSystTot");

    TH1D* hSum = (TH1D*)hOther->Clone("hSum"); 
    for (int bin = 1; bin <= hSum->GetNbinsX(); bin++){
        int binr = 1;
        if (bin >= systcounter[0] && bin < systcounter[1]) binr = 2;
        if (bin >= systcounter[1] && bin <= systcounter[2]) binr = 3;
        double sum = TMath::Sqrt( 
            hOther->GetBinContent(bin)*hOther->GetBinContent(bin) + 
            hExtrap->GetBinContent(bin)*hExtrap->GetBinContent(bin) +
            hR->GetBinContent(binr)*hR->GetBinContent(binr)         
         );
        hSum->SetBinContent(bin,sum);
    }

    hSum->SetLineColor(kBlack);
    hSum->SetLineWidth(2);
    hOther->SetLineColor(kBlue);
    hR->SetLineColor(kGreen+2);
    hExtrap->SetLineColor(kRed);
    hExtrap->SetLineWidth(1);

    TLegend* lnee = new TLegend (0.27029276,0.70894,0.556995,0.83088);
    lnee->SetBorderSize(0);
    lnee->AddEntry(hSum ,"Total","L");
    lnee->AddEntry(hExtrap ,"Extrapolation","L");
    lnee->AddEntry(hR ,"Cut Variations (R factor)","L");
    lnee->AddEntry(hOther ,"Other contribution","L");
    lnee->SetTextSize(0.032);
    lnee->SetTextFont(42);
    
   

    TCanvas* c = new TCanvas("c","",1100,1000);
    c->SetRightMargin(0.09);
    c->SetLeftMargin(0.25);
    c->SetBottomMargin(0.15);
    hSum->GetYaxis()->SetRangeUser(0.,0.15);
    hSum->GetYaxis()->SetTitle("Syst. Contribution %");
    if (fWhichEstimator.Contains("ZDC")) hSum->GetXaxis()->SetTitle(Form("%s percentile","(#sqrt{s} - ZDC)" ));
    if (fWhichEstimator.Contains("V0M")) hSum->GetXaxis()->SetTitle(Form("%s percentile","V0M" ));

    TLatex *xlaber = new TLatex();
    xlaber->SetTextFont(42);
    xlaber-> SetNDC();
    xlaber-> SetTextColor(1);
    xlaber-> SetTextSize(0.03);
    xlaber-> SetTextAlign(22);
    xlaber-> SetTextAngle(0);
    

    hSum->SetStats(0);
    hSum->Draw();
    hOther->Draw("SAME HIST");
    hR->Draw("SAME HIST");
    hExtrap->Draw("SAME HIST");
    lnee->Draw("SAME");
    if (fWhichEstimator.Contains("V0M")) xlaber-> DrawLatex(0.45, 0.6, Form("%s fixed [%.0f-%.0f]","(#sqrt{s} - ZDC)",lLoEE,lHiEE));
    if (fWhichEstimator.Contains("ZDC")) xlaber-> DrawLatex(0.45, 0.6, Form("%s fixed [%.0f-%.0f]","V0M",lLoMult,lHiMult));


    c->SaveAs(Form("ContributoYields_%s%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.png","Xi",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE));

}