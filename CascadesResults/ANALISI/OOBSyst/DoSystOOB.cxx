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

void DoSystOOB(){

    TString fWhichParticle = "XiMinus";
    TString fWhichAntiParticle = "XiPlus";
    TH1D* hdev[5], *htot[5],*htotDef, *hResPart[5], *hResAntiPart[5];
    //
    // files
    // default
    TFile* fileDefPart  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100.root",
        fWhichParticle.Data()),"READ");
    TH1D* hDefPart = (TH1D*)fileDefPart->Get("fHistPtXiMinus");
    TFile* fileDefAntiPart  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100.root",
        fWhichAntiParticle.Data()),"READ");
    TH1D* hDefAntiPart = (TH1D*)fileDefAntiPart->Get("fHistPtXiPlus");
    // ITS All
    TFile* fileITSAllPtPart  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100-AllPtITSrefit.root",
        fWhichParticle.Data()),"READ");
    hResPart[0] = (TH1D*)fileITSAllPtPart->Get("fHistPtXiMinus");
    TFile* fileITSAllPtAntiPart  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100-AllPtITSrefit.root",
        fWhichAntiParticle.Data()),"READ");
    hResAntiPart[0] = (TH1D*)fileITSAllPtAntiPart->Get("fHistPtXiPlus");
    // ITS > 2 GeV/c
   /* TFile* fileITSgt2Part  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100-ITSgt2GeVc.root",
        fWhichParticle.Data()),"READ");
    hResPart[1] = (TH1D*)fileITSgt2Part->Get("fHistPtXiMinus");
    TFile* fileITSgt2AntiPart  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100-ITSgt2GeVc.root",
        fWhichAntiParticle.Data()),"READ");
    hResAntiPart[1] = (TH1D*)fileITSgt2AntiPart->Get("fHistPtXiPlus");*/
    // ITS < 2 GeV/c
    TFile* fileITSlt2Part  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100-ITSlt2GeVc.root",
        fWhichParticle.Data()),"READ");
    hResPart[1] = (TH1D*)fileITSlt2Part->Get("fHistPtXiMinus");
    TFile* fileITSlt2AntiPart  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100-ITSlt2GeVc.root",
        fWhichAntiParticle.Data()),"READ");
    hResAntiPart[1] = (TH1D*)fileITSlt2AntiPart->Get("fHistPtXiPlus");
    // TOF > 2 GeV/c
    TFile* fileTOFgt2Part  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100-TOFgt2GeVc.root",
        fWhichParticle.Data()),"READ");
    hResPart[2] = (TH1D*)fileTOFgt2Part->Get("fHistPtXiMinus");
    TFile* fileTOFgt2AntiPart  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100-TOFgt2GeVc.root",
        fWhichAntiParticle.Data()),"READ");
    hResAntiPart[2] = (TH1D*)fileTOFgt2AntiPart->Get("fHistPtXiPlus");
    // TOF < 2 GeV/c
    /*TFile* fileTOFlt2Part  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100-TOFlt2GeVc.root",
        fWhichParticle.Data()),"READ");
    hResPart[4] = (TH1D*)fileTOFlt2Part->Get("fHistPtXiMinus");
    TFile* fileTOFlt2AntiPart  = TFile::Open(Form("Results-%s-13TeV-V0M_000_100_ZDC_000_100-TOFlt2GeVc.root",
        fWhichAntiParticle.Data()),"READ");
    hResAntiPart[4] = (TH1D*)fileTOFlt2AntiPart->Get("fHistPtXiPlus");*/

    htotDef = (TH1D*)hDefPart -> Clone("htotDef");
    htotDef->Reset();
    for (int bin = 1; bin <= htotDef->GetNbinsX(); bin++){
        htotDef->SetBinContent(bin, hDefPart->GetBinContent(bin) + hDefAntiPart->GetBinContent(bin));
        htotDef->SetBinError(bin, TMath::Sqrt(hDefPart->GetBinError(bin)*hDefPart->GetBinError(bin)+hDefAntiPart->GetBinError(bin)*hDefAntiPart->GetBinError(bin)));
    }        

    for (int i = 0; i < 3; i++){

        htot[i] = (TH1D*)hResPart[i] -> Clone(Form("htot%i",i));
        htot[i]->Reset();
        
        for (int bin = 1; bin <= htot[i]->GetNbinsX(); bin++){
            htot[i]->SetBinContent(bin, hResPart[i]->GetBinContent(bin) + hResAntiPart[i]->GetBinContent(bin));
            htot[i]->SetBinError(bin, TMath::Sqrt(hResPart[i]->GetBinError(bin)*hResPart[i]->GetBinError(bin)+hResAntiPart[i]->GetBinError(bin)*hResAntiPart[i]->GetBinError(bin)));
        }

        hdev[i] = (TH1D*)htot[i] -> Clone(Form("hdev%i",i));
        //hdev[i]->Reset();

        DivideAndComputeRogerBarlow(hdev[i],htotDef);

        /*for (int bin = 1; bin <= htot[i]->GetNbinsX(); bin++){
            hdev[i]->SetBinContent(bin,TMath::Abs( (htot[i]->GetBinContent(bin)-htotDef->GetBinContent(bin) ) / htotDef->GetBinContent(bin)));
            hdev[i]->SetBinError(bin,(Error(htot[i]->GetBinError(bin),htotDef->GetBinError(bin))) / htotDef->GetBinContent(bin));
        }        */
    }

      //Compute Max Deviation
    double binvalue[3];
    double maxvalue = 0.;
    TH1D* hMaxDev = (TH1D*)htot[0] -> Clone("hMaxDev");
    hMaxDev->Reset();
    for (int i = 1; i<= hdev[0]->GetNbinsX(); i++){
        for (int k = 0; k < 3; k++){ 
        double dev = PassRogerBarlowCriterion(1,hdev[k],i);
        binvalue[k] = TMath::Abs(dev-1);
        cout << "bin " << i << binvalue[k] << endl;
        }
        maxvalue = binvalue[0];
        int counter = 0;
        for (int k = 0; k < 3; k++){ 
        maxvalue = max(maxvalue,binvalue[k]);
        if (maxvalue == binvalue[k]) counter = k;
        }   
        cout << "max val " << maxvalue << endl;
        hMaxDev->SetBinError(i,hdev[counter]->GetBinError(i));
        hMaxDev->SetBinContent(i,TMath::Abs(maxvalue));
        //if (TMath::Abs(maxvalue) == 0) hMaxDev->SetBinContent(i,0.0);
    }

     hdev[0]->SetLineColor(kRed);
    hdev[0]->SetMarkerColor(kRed);
    hdev[0]->SetMarkerStyle(20);
    hdev[0]->SetMarkerSize(1.8);
    hdev[0]->SetLineWidth(1);
     hdev[1]->SetLineColor(kOrange);
    hdev[1]->SetMarkerColor(kOrange);
    hdev[1]->SetMarkerStyle(20);
    hdev[1]->SetMarkerSize(1.8);
    hdev[1]->SetLineWidth(1);
     hdev[2]->SetLineColor(kMagenta+2);
    hdev[2]->SetMarkerColor(kMagenta+2);
    hdev[2]->SetMarkerStyle(20);
    hdev[2]->SetMarkerSize(1.8);
    hdev[2]->SetLineWidth(1);
     /*hdev[3]->SetLineColor(kGre;en+2);
    hdev[3]->SetMarkerColor(kGreen+2);
    hdev[3]->SetMarkerStyle(20);
    hdev[3]->SetMarkerSize(1.2);
    hdev[3]->SetLineWidth(1);
     hdev[4]->SetLineColor(kBlue);
    hdev[4]->SetMarkerColor(kBlue);
    hdev[4]->SetMarkerStyle(20);
    hdev[4]->SetMarkerSize(1.2);
    hdev[4]->SetLineWidth(1);*/ 


    TCanvas* c = new TCanvas("c","c",1600,1500);
    c->SetRightMargin(0.09);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);

    TLegend* legend = new TLegend (0.3,.75,0.91,0.89);  
    legend->SetBorderSize(0);
    legend->AddEntry(hdev[0],"ITS match 1 daughter all p_{T}","LEP");
    legend->AddEntry(hdev[1],"ITS match 1 daughter p_{T} < 2 GeV/c","LEP");
    legend->AddEntry(hdev[2],"TOF match 1 daughter p_{T} > 2 GeV/c","LEP");

    hdev[0]->GetYaxis()->SetRangeUser(0.9,1.1);
    hdev[0]->SetStats(0);
    hdev[0]->GetYaxis()->SetTitle("Y_{OOB cond. varied} / Y_{OOB cond. def.}");
    hdev[0]->Draw();
    hdev[0]->SetTitle("");

    for (int i = 1; i < 3; i++){
        hdev[i]->Draw("SAME");
    }
    legend->Draw("SAME");

    TCanvas* c1 = new TCanvas("c1","c1",1600,1500);
    c1->SetRightMargin(0.09);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.15);
    hMaxDev->SetStats(0);
    hMaxDev->GetYaxis()->SetTitle("max. rel. deviation (if > 1 #sigma_{RB})");
    hMaxDev->GetYaxis()->SetRangeUser(0.,.03);
    hMaxDev->SetLineColor(kBlack);
    hMaxDev->SetLineWidth(2);
    hMaxDev->SetTitle("");
    hMaxDev->Draw("HIST");

    //Output File
    TFile* Write = new TFile ("OOBSystMB.root", "RECREATE");
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