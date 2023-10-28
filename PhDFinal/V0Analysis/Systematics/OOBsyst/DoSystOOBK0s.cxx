#include <iostream>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TPad.h"
#include "TPaveText.h"
#include <TLatex.h>

double Error(double errA, double errB);
void DivideAndComputeRogerBarlow( TH1D* h1, TH1D *h2 );
double PassRogerBarlowCriterion(int nsigmas, TH1D* h, int bin);

void DoSystOOBK0s(TString lCascType = "K0Short"){

    TH1D* hdev[5], *htot[5],*htotDef, *hResPart[5], *hResAntiPart[5];

    TString particle = "", antiparticle = "";
    if (lCascType.Contains("Xi")){
        particle = "XiMinus";
        antiparticle = "XiPlus";
    }
    if (lCascType.Contains("Omega")){
        particle = "OmegaMinus";
        antiparticle = "OmegaPlus";
    }
    if (lCascType.Contains("Lambda")){
        particle = "Lambda";
        antiparticle = "AntiLambda";
    }
    if(lCascType.Contains("K0Short"))
    {
        particle = "K0Short";
        antiparticle = "";
    }

    //
    // files
    // default
    TFile* fileDefPart  = TFile::Open(Form("../files/K0Short/Results-%s-13TeV-SPDClusters_000_100_V0M_000_100_FDNoFD.root",
        particle.Data()),"READ");
    TH1D* hDefPart = (TH1D*)fileDefPart->Get(Form("fHistPt%s",particle.Data()));
    // ITS All
    TFile* fileITSAllPtPart  = TFile::Open(Form("Results-%s-13TeV-SPDClusters_000_100_V0M_000_100_FDNoFD_ITSOOB10GeV.root",
        particle.Data()),"READ");
    hResPart[0] = (TH1D*)fileITSAllPtPart->Get(Form("fHistPt%s",particle.Data()));
    // ITS < 2 GeV/c
    TFile* fileITSlt2Part  = TFile::Open(Form("Results-%s-13TeV-SPDClusters_000_100_V0M_000_100_FDNoFD_ITSOOB2GeV.root",
        particle.Data()),"READ");
    hResPart[1] = (TH1D*)fileITSlt2Part->Get(Form("fHistPt%s",particle.Data()));
    // TOF > 2 GeV/c
    TFile* fileTOFgt2Part  = TFile::Open(Form("Results-%s-13TeV-SPDClusters_000_100_V0M_000_100_FDNoFD_TOFOOB.root",
        particle.Data()),"READ");
    hResPart[2] = (TH1D*)fileTOFgt2Part->Get(Form("fHistPt%s",particle.Data()));

    for (int i = 0; i < 3; i++){

        DivideAndComputeRogerBarlow(hResPart[i], hDefPart);
    }

    //Compute Max Deviation
    double binvalue[3];
    double maxvalue = 0.;
    TH1D *hMaxDev = (TH1D *)hDefPart->Clone("hMaxDev");
    hMaxDev->Reset();
    for (int i = 1; i<= hResPart[0]->GetNbinsX(); i++){
        for (int k = 0; k < 3; k++){
        double dev = PassRogerBarlowCriterion(1,hResPart[k],i);
        binvalue[k] = TMath::Abs(dev-1);
        }
        maxvalue = binvalue[0];
        int counter = 0;
        for (int k = 0; k < 3; k++){
        maxvalue = max(maxvalue,binvalue[k]);
        if (maxvalue == binvalue[k]) counter = k;
        }
        hMaxDev->SetBinError(i,hResPart[counter]->GetBinError(i));
        hMaxDev->SetBinContent(i,TMath::Abs(maxvalue));
        //if (TMath::Abs(maxvalue) == 0) hMaxDev->SetBinContent(i,0.0);
    }

    hResPart[0]->SetLineColor(kRed);
    hResPart[0]->SetMarkerColor(kRed);
    hResPart[0]->SetMarkerStyle(20);
    hResPart[0]->SetMarkerSize(1.8);
    hResPart[0]->SetLineWidth(1);
    hResPart[1]->SetLineColor(kOrange);
    hResPart[1]->SetMarkerColor(kOrange);
    hResPart[1]->SetMarkerStyle(20);
    hResPart[1]->SetMarkerSize(1.8);
    hResPart[1]->SetLineWidth(1);
    hResPart[2]->SetLineColor(kMagenta+2);
    hResPart[2]->SetMarkerColor(kMagenta+2);
    hResPart[2]->SetMarkerStyle(20);
    hResPart[2]->SetMarkerSize(1.8);
    hResPart[2]->SetLineWidth(1);


    TCanvas* c1 = new TCanvas("c1","c1",1600,800);
    c1->Divide(2);
    c1->cd(1);
    c1->cd(1)->SetRightMargin(0.09);
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetBottomMargin(0.15);

    TLegend* legend = new TLegend (0.3,.75,0.91,0.89);
    legend->SetBorderSize(0);
    legend->AddEntry(hResPart[0],"ITS match 1 daughter all p_{T}","LEP");
    legend->AddEntry(hResPart[1],"ITS match 1 daughter p_{T} < 2 GeV/c","LEP");
    legend->AddEntry(hResPart[2],"TOF match 1 daughter p_{T} > 2 GeV/c","LEP");

    hResPart[0]->GetYaxis()->SetRangeUser(0.85,1.15);
    hResPart[0]->SetStats(0);
    hResPart[0]->GetYaxis()->SetTitle("Y_{OOB cond. varied} / Y_{OOB cond. def.} (K0Short)");
    hResPart[0]->Draw();
    hResPart[0]->SetTitle("");
    hResPart[0]->GetXaxis()->SetRangeUser(0., 10.);

    for (int i = 1; i < 3; i++){
        hResPart[i]->Draw("SAME");
    }
    legend->Draw("SAME");

    c1->cd(2);
    c1->cd(2)->SetRightMargin(0.09);
    c1->cd(2)->SetLeftMargin(0.18);
    c1->cd(2)->SetBottomMargin(0.15);
    hMaxDev->SetStats(0);
    hMaxDev->GetYaxis()->SetTitle("max. rel. deviation (if > 1 #sigma_{RB})");
    hMaxDev->GetYaxis()->SetRangeUser(0.,.05);
    hMaxDev->GetXaxis()->SetRangeUser(0., 10.);
    hMaxDev->SetLineColor(kBlack);
    hMaxDev->GetYaxis()->SetTitleOffset(2.4);
    hMaxDev->SetLineWidth(2);
    hMaxDev->SetTitle("");
    TH1F *hMaxDev1 = (TH1F *)hMaxDev->Clone("hMaxDev1");
    hMaxDev1->Scale(1./2);
    hMaxDev1->GetYaxis()->SetTitle("#frac{max. rel. deviation (if > 1 #sigma_{RB})}{2}");
    hMaxDev1->Draw("HIST");

    c1->SaveAs(Form("images/OOBSyst%s.png",lCascType.Data()));
    c1->SaveAs(Form("images/OOBSyst%s.root",lCascType.Data()));

    //Output File
    TFile* Write = new TFile (Form("OOBSyst%s.root",lCascType.Data()), "RECREATE");
    hMaxDev->Write();

}

double Error(double errA, double errB){

    double errtot = TMath::Sqrt( TMath::Abs( errA*errA - errB*errB ) );
    return errtot;

};


//----------------------------------------------------------------------------------------------------
void DivideAndComputeRogerBarlow( TH1D* h1, TH1D *h2 ){
  //Use Roger Barlow "sigma_{delta}" as errors for ratios
  Double_t lh1NBins = h1->GetNbinsX();
  Double_t lh2NBins = h2->GetNbinsX();

  if( lh1NBins != lh2NBins ){
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return;
  }

  Double_t lSigmaDelta[100];
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){
    //Computation of roger barlow sigma_{delta}
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(h1->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
    //Computation of relationship to h2 for plotting in ratio plot
    if ( h2->GetBinContent(i) > 1e-12 ){
      lSigmaDelta[i] /= h2->GetBinContent(i);
    }else{
      lSigmaDelta[i] = 0;
    }
  }
  //Regular Division
  h1 -> Divide( h2 );
  //Replace Erorrs
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){
    h1->SetBinError(i, lSigmaDelta[i]);
  }
}

//----------------------------------------------------------------------------------------------------
double PassRogerBarlowCriterion(int nsigmas, TH1D* h, int bin){
  double dev = h->GetBinContent(bin);
  double RBsigma = h->GetBinError(bin);

  if (TMath::Abs(dev-1)>nsigmas*RBsigma) {return dev;}
  else {cout << "NO" << endl;return 1.;}
}