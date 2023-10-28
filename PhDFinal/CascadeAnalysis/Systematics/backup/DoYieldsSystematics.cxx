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
#include "TH3.h"
#include "TH2.h"
#include "TGraphErrors.h"


#endif
using namespace std;
#include <TString.h>




/* definition of the fields in the histogram returned */
enum EValue_t {
  kYield = 1,
  kYieldStat,
  kYieldSysHi,
  kYieldSysLo,
  kMean,
  kMeanStat,
  kMeanSysHi,
  kMeanSysLo, 
  kExtra,
  kChi

};

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) ;

TH1 *YieldMean(TH1 *hstat, TH1 *hsys, TF1 *f = NULL, Double_t minfit=0.4,Double_t maxfit=8., TString logfilename="log.root", Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q");
void YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, Double_t &integral, Double_t &mean, Double_t &extra, Bool_t printinfo=kFALSE);
TH1* YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, Double_t min, Double_t binwidth = 0.01);
TH1 * YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, Double_t max, Double_t binwidth = 0.1);
TH1 * YieldMean_ReturnRandom(TH1 *hin);
TH1 * YieldMean_ReturnCoherentRandom(TH1 *hin);
TH1 *YieldMean_ReturnExtremeHisto(TH1 *hin, Float_t sign = 1.);
TH1 *YieldMean_ReturnExtremeHardHisto(TH1 *hin);
TH1 *YieldMean_ReturnExtremeSoftHisto(TH1 *hin);
TH1 * YieldMean_ReturnExtremeLowHisto(TH1 *hin);
TH1 * YieldMean_ReturnExtremeHighHisto(TH1 *hin);

void beautifygraph(TGraphErrors* g);
TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);

void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 );
Bool_t PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin);
TH1F* ComputeYieldsforSystematics( 
    TString filename = "",
    TString lCascType = "Omega", 
    TString fWhichVarEstimator = "V0M", 
    TString fWhichFixedEstimator = "SPDClusters", 
    Double_t lFixedLo = 0., 
    Double_t lFixedHi = 100., 
    Double_t lFixedLo2 = 0., 
    Double_t lFixedHi2 = 100., 
    TString lWhichSystVar = "V0Radius",
    Float_t* percentile = 0x0,
    Int_t nbin = 0
);

//--------------------------------------------------------------------------------
//----------------------------- MAIN FUNCTION ------------------------------------
//--------------------------------------------------------------------------------
void DoYieldsSystematics(
  TString fWhichVarEstimator = "V0M", 
  TString fWhichFixedEstimator = "SPDClusters" ,
  TString fWhichParticle = "Omega", 
  Double_t lFixedLo2 = 10., 
  Double_t lFixedHi2 = 20.,
  Double_t lFixedLo = 0., 
  Double_t lFixedHi = 100., 
  Bool_t kDoMB = kTRUE  ){
    cout<<"  "<<endl;
    cout<<"-------------------------------------------------------------"<<endl;
    cout<<"                  Make systematics for yields                   "<<endl;
    cout<<"-------------------------------------------------------------"<<endl;
    cout<<endl;

    gROOT->SetBatch(kTRUE);

    Float_t * percentile;
    Long_t percbinnumb;
    Float_t p1[] = {0,5,10,20,30,40,50,100};
    // {0,10,30,50,100};
    //{0,10,20,30,40,50,60,70,100};
    //{0,5,10,20,30,40,50,100};
    //{0,10,20,30,40,50,60,70,100};
    //{0,20,30,40,50,60,70,100};
    //{0,5,10,20,30,40,50,100}; //7
    Long_t n1 = sizeof(p1)/sizeof(Float_t) - 1;
    Float_t pmb[] = {0,100}; //7
    Long_t nmb = sizeof(pmb)/sizeof(Float_t) - 1;
  
    if (kDoMB){
      percentile = pmb;
      percbinnumb = nmb;
    } else{
      percentile = p1;
      percbinnumb = n1;
    }

    TString outputfilename = Form("SystematicsYields-%s-13TeV-%s_%03.0f_%03.0f.root",fWhichParticle.Data(),fWhichFixedEstimator.Data(),lFixedLo2, lFixedHi2);
    if (kDoMB) outputfilename = Form("SystematicsYields-%s-13TeV-INELgt0.root",fWhichParticle.Data());

    // Selection classes
    TH1F* hV0Radius = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"V0Radius",percentile,percbinnumb);
    TH1F* hCascRadius = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"CascRadius",percentile,percbinnumb);
    TH1F* hNegToPV = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"DCANegToPV",percentile,percbinnumb);
    TH1F* hPosToPV = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"DCAPosToPV",percentile,percbinnumb);
    TH1F* hV0ToPV = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"DCAV0ToPV",percentile,percbinnumb);
    TH1F* hBachToPV = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"DCABachToPV",percentile,percbinnumb);
    TH1F* hV0Daught = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"DCAV0Daughters",percentile,percbinnumb);
    TH1F* hCascDaught = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"DCACascDaughters",percentile,percbinnumb);
    TH1F* hV0CosPA = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"V0CosPA",percentile,percbinnumb);
    TH1F* hCascCosPA = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"CascCosPA",percentile,percbinnumb);
    TH1F* hV0Mass = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"V0Mass",percentile,percbinnumb);
    TH1F* hCompetingSpecies = 0x0;
    TH1F* hPLT = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"ProperLifetime",percentile,percbinnumb);
    TH1F* hTPCNSigmas = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"TPCPIDNSigmas",percentile,percbinnumb);
    TH1F* hTPCNClusters = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"TPCNClusters",percentile,percbinnumb);
    TH1F* hSigExtBinCount = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"SigExtBinCount",percentile,percbinnumb);
    //TH1F* hGeantFluka = ComputeYieldsforSystematics(outputfilename,fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,lFixedLo2,lFixedHi2,"GeantFlukaCorrection",percentile,percbinnumb);
    
    // Write Max Dev on output File
    TFile* Write = new TFile (outputfilename, "UPDATE");
    hV0Radius->Write();
    hCascRadius->Write();
    hNegToPV->Write();
    hPosToPV->Write();
    hBachToPV->Write();
    hV0Daught->Write();
    hCascDaught->Write();
    hV0CosPA->Write();
    hCascCosPA->Write();
    hV0ToPV->Write();
    hPLT->Write();
    hV0Mass->Write();
    hTPCNClusters->Write();
    hTPCNSigmas->Write();
    hSigExtBinCount->Write();
    //hGeantFluka->Write();

  }

TH1F* ComputeYieldsforSystematics( 
    TString filename = "",
    TString lCascType = "Xi", 
    TString fWhichVarEstimator = "V0M", 
    TString fWhichFixedEstimator = "SPDClusters", 
    Double_t lFixedLo = 0., 
    Double_t lFixedHi = 100., 
    Double_t lFixedLo2 = 0., 
    Double_t lFixedHi2 = 100., 
    TString lWhichSystVar = "V0Radius",
    Float_t* percentile = 0x0,
    Int_t nbin = 0
){

    cout << "\n---------------------------------------------------------------" << endl;
    cout << endl << "\n---> Do Systematics for: "<< lWhichSystVar << endl;
    cout << "\n---------------------------------------------------------------" << endl;

    double minfit, maxfit;
    Double_t mass = 1.32171;

    TString fWhichParticle = "", fWhichAntiParticle = "";
    if (lCascType.Contains("Xi")){
        fWhichParticle = "XiMinus";
        fWhichAntiParticle = "XiPlus";
        minfit = .6;
        maxfit = 6.5;
        //mass is default
    }
    if (lCascType.Contains("Omega")){
        fWhichParticle = "OmegaMinus";
        fWhichAntiParticle = "OmegaMinus";
        minfit = 0.9;
        maxfit = 6.5;
        mass = 1.67245;
    }
    if (lCascType.Contains("Lambda")){
        fWhichParticle = "Lambda";
        fWhichAntiParticle = "AntiLambda";
        minfit = 0.4;
        maxfit = 8.;
        mass = 1.115683;
    }
     
    int const binnumber = nbin;

    Float_t perc1 [nbin];
    Float_t perc2 [nbin];

    for (int i = 0; i<(nbin); i++){
      perc1[i] = percentile[i];
      perc2[i] = percentile[i+1];
    }
   
    // Histos
    TH1F* lHistPt[binnumber][5];
    TH1F* lSystPt[binnumber][5];
    TH1F* lHistPtClone[binnumber][5];
    TH1F* lAntiHistPt[binnumber][5];
    TH1F* lAntiSystPt[binnumber][5];
    TH1F* hout[binnumber][5];
    Double_t Yield[binnumber][5];
    Double_t YieldStat[binnumber][5];

    TFile* filesyst = TFile::Open(Form("SystematicsFinalResults-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root",lCascType.Data(),0.,100.,0.,100.));
    TH1F* hsyst = (TH1F*)filesyst->Get("hSystTot");
  
    // Fit function
    TF1* LevyTsallisfunc = LevyTsallis("LevyTsallisfunc", mass);

    // Get spectra for this cut
    TFile* lResultsFile[binnumber][5];
    TFile* lAntiResultsFile[binnumber][5];
    for (int i = 0; i < binnumber; i++){ //loop over selection classes           
        double v0min = lFixedLo, v0max = lFixedHi, spdmin = lFixedLo, spdmax = lFixedHi;
        if (fWhichFixedEstimator.Contains("SPDClusters")) {
          v0min = perc1[i];
          v0max = perc2[i];
        }
        if (fWhichFixedEstimator.Contains("V0M")) {
          spdmin = perc1[i];
          spdmax = perc2[i];
        }   
        //
        TString lSystFile = "sistematiche/Results-Systematics";
        lSystFile.Append( Form( "-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f-", fWhichParticle.Data(), spdmin, spdmax, v0min, v0max) );
        TString lAntiSystFile = "sistematiche/Results-Systematics";
        lAntiSystFile.Append( Form( "-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f-", fWhichAntiParticle.Data(), spdmin, spdmax, v0min, v0max) );
        //
        lResultsFile[i][0] = TFile::Open(Form( "sistematiche/Results-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root", fWhichParticle.Data(), spdmin, spdmax, v0min, v0max));
        lAntiResultsFile[i][0] = TFile::Open(Form( "sistematiche/Results-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root", fWhichAntiParticle.Data(), spdmin, spdmax, v0min, v0max));
        
        //Var cut (numerator)
        for (int j = 1; j < 5; j++){
          //Particle
          TString lDataFilename = lSystFile + lWhichSystVar + Form("-%i.root",j);
          lResultsFile[i][j] = TFile::Open(lDataFilename.Data());
          lHistPt[i][j] = (TH1F*)lResultsFile[i][j]->Get(Form("fHistPt%s", fWhichParticle.Data()));
          //Anti Particle
          TString lAntiDataFilename = lAntiSystFile + lWhichSystVar + Form("-%i.root",j);
          lAntiResultsFile[i][j] = TFile::Open(lAntiDataFilename.Data());
          lAntiHistPt[i][j] = (TH1F*)lAntiResultsFile[i][j]->Get(Form("fHistPt%s", fWhichAntiParticle.Data()));
        }  
        
        //Def cut (denominator)
        lHistPt[i][0] = (TH1F*)lResultsFile[i][0]->Get(Form("fHistPt%s", fWhichParticle.Data()));       
        lAntiHistPt[i][0] = (TH1F*)lAntiResultsFile[i][0]->Get(Form("fHistPt%s", fWhichAntiParticle.Data()));
      
        // Add Spectra
        for (int j = 0; j < 5; j++){
          lHistPtClone[i][j] = (TH1F*)lHistPt[i][j]->Clone(Form("lHistPt%i%i",i,j));
          lHistPtClone[i][j]->Reset(); //be clean shall we?
          
          for (int bin = 1; bin <=  lHistPt[i][j]->GetNbinsX(); bin ++){
            lHistPtClone[i][j]->SetBinContent(bin, lHistPt[i][j]->GetBinContent(bin)+ lAntiHistPt[i][j]->GetBinContent(bin));
            lHistPtClone[i][j]->SetBinError(bin, TMath::Sqrt(lHistPt[i][j]->GetBinError(bin)*lHistPt[i][j]->GetBinError(bin)+
                                                      lAntiHistPt[i][j]->GetBinError(bin)*lAntiHistPt[i][j]->GetBinError(bin)));   
          }     
        }
        for (int j = 0; j < 5; j++){
          lSystPt[i][j] = (TH1F*)lHistPtClone[i][j]->Clone(Form("lSystPt%i%i",i,j));
          for (int bin = 1; bin <=lSystPt[i][j]->GetNbinsX(); bin ++){
            lSystPt[i][j]->SetBinError(bin, hsyst->GetBinContent(bin)*lHistPtClone[i][j]->GetBinContent(bin));   
          }      
        }
    }  

    // Get the yields
      for (int i = 0; i < binnumber; i++){ //loop over selection classes  
      for (int j = 0; j < 5; j++){
        hout[i][j] = (TH1F*)YieldMean(lHistPtClone[i][j],lSystPt[i][j],LevyTsallisfunc,minfit,maxfit);
        // Get yield
        Yield[i][j] = hout[i][j]->GetBinContent(1);
        YieldStat[i][j] = hout[i][j]->GetBinContent(2);
      }
    }//end loop over classes 

    TH1F* hYield[5];
    for (int j = 0; j < 5; j++){
      hYield[j] = new TH1F(Form("hYield%i-%s",j,fWhichVarEstimator.Data()), Form(";percentile %s (%);dN/dy",fWhichVarEstimator.Data()),binnumber,0.,(double)binnumber);
      for (int i = 0; i<(nbin); i++){
        hYield[j]->GetXaxis()->SetBinLabel(i+1,Form("%.0f-%.0f %s",perc1[i],perc2[i],"%"));
      }
      for (int nbin = 1; nbin <= hYield[j]->GetNbinsX(); nbin++){
        hYield[j]->SetBinContent(nbin, Yield[nbin-1][j]);
        hYield[j]->SetBinError(nbin, YieldStat[nbin-1][j]);
      }
    }

    // Divide following Roger-Barlow prescription<
    for (int j = 1; j < 5; j++){
      DivideAndComputeRogerBarlow(hYield[j],hYield[0]);
    }

    TH1F* hMaxDev = (TH1F*)hYield[0]->Clone("hMaxDev");
    hMaxDev->Reset();

    //Compute Max Deviation
    double binvalue[4];
    double maxvalue = 0.;
    for (int i = 1; i<= hYield[0]->GetNbinsX(); i++){
      for (int k = 1; k < 5; k++){ 
        binvalue[k-1] = 0;
        // Pass Roger Barlow at 1 sigmas
        if (PassRogerBarlowCriterion(1, hYield[k], i)) binvalue[k-1] = TMath::Abs(hYield[k]->GetBinContent(i)-1);
      }
      maxvalue = binvalue[0];
      int counter = 1;
      for (int k = 2; k < 5; k++){ 
        maxvalue = max(maxvalue,binvalue[k-1]);
        if (maxvalue == binvalue[k-1]) counter = k;
      }   
      hMaxDev->SetBinError(i,hYield[counter]->GetBinError(i));
      hMaxDev->SetBinContent(i,TMath::Abs(maxvalue));
    }
    hMaxDev->SetTitle(Form("Syst %s",lWhichSystVar.Data()));
    hMaxDev->SetName(Form("hMaxDev%s%s",lWhichSystVar.Data(),fWhichVarEstimator.Data()));

    TLatex *labbig = new TLatex();
    labbig->SetTextFont(42);
    labbig->SetNDC();
    labbig->SetTextColor(1);
    labbig->SetTextSize(0.04);
    labbig->SetTextAlign(22);
    labbig->SetTextAngle(0);

    // Prepare Canvases
    //Max Deviation
    hMaxDev->GetXaxis()->SetRangeUser(0.,(double)binnumber);
    hMaxDev->GetYaxis()->SetRangeUser(-0.0005,.2);
    hMaxDev->SetYTitle("max rel. dev.");
    hMaxDev->SetTitle(Form("%s",lWhichSystVar.Data()));
    hMaxDev->GetYaxis()->SetTitleSize(0.05);
    hMaxDev->GetYaxis()->SetTitleOffset(1.0);
    hMaxDev->GetXaxis()->SetTitleSize(0.05);
    hMaxDev->GetXaxis()->SetTitleOffset(1.);
    hMaxDev->SetMarkerStyle(20);
    hMaxDev->SetMarkerSize(1.5);
    hMaxDev->SetMarkerColor(kBlack);
    hMaxDev->SetLineColor(kBlack);
    TCanvas* maxdev = new TCanvas("maxdev"," ",1000,800);
    maxdev->SetRightMargin(0.15);
    maxdev->SetLeftMargin(0.15);
    maxdev->SetBottomMargin(0.15);
    maxdev->SetGridy();
    maxdev->SetGridx();
    hMaxDev->SetStats(kFALSE);  
    hMaxDev->GetXaxis()->SetLabelSize(0.05);
    hMaxDev->Draw();
    if (lFixedLo2 == 0 && lFixedHi2 == 100.) {
      labbig-> DrawLatex(0.45, 0.75, Form("MB (INEL>0)")); 
    } else labbig-> DrawLatex(0.45, 0.75, Form("applied to %s [%.0f-%.0f] points",fWhichFixedEstimator.Data(),lFixedLo2,lFixedHi2)); 
    maxdev->SaveAs(Form("imagesAN/yieldsyst/fixed%s/Yields-%s-MaxRelDev%s-%s_%.0f%0.f.png",fWhichFixedEstimator.Data(),fWhichParticle.Data(),lWhichSystVar.Data(),fWhichFixedEstimator.Data(),lFixedLo2,lFixedHi2));

    // Same yields raios in canvases
    //Cut Variation
    TCanvas* cutvar = new TCanvas("cutvar"," ",900,800);
    cutvar->SetRightMargin(0.08);
    cutvar->SetLeftMargin(0.15);
    cutvar->SetBottomMargin(0.15);
    cutvar->SetGridy();
    //
    TLegend* legend = new TLegend (0.7,.7,0.91,0.89);  
    legend->SetBorderSize(0);
    //
    hYield[1]->SetMarkerColor(kRed);
    hYield[2]->SetMarkerColor(kBlue+1);
    hYield[3]->SetMarkerColor(kGreen+3);
    hYield[4]->SetMarkerColor(kAzure+9);
    hYield[1]->SetLineColor(kRed);
    hYield[2]->SetLineColor(kBlue+1);
    hYield[3]->SetLineColor(kGreen+3);
    hYield[4]->SetLineColor(kAzure+9);    
    legend->AddEntry(hYield[1],"very loose","LEP");
    legend->AddEntry(hYield[2],"loose","LEP");
    legend->AddEntry(hYield[3],"tight","LEP");
    legend->AddEntry(hYield[4],"very tight","LEP");

    for (int k = 1; k < 5; k++){
      hYield[k]->SetName(Form("hYield-%s-%s-%i",lWhichSystVar.Data(),fWhichVarEstimator.Data(),k));
      hYield[k]->SetTitle(Form("hYield-%s-%i",lWhichSystVar.Data(),k));
      hYield[k]->GetXaxis()->SetRangeUser(0.,(double)binnumber);
      hYield[k]->GetYaxis()->SetRangeUser(0.8,1.2);
      hYield[k]->SetYTitle("#LT dN/dy #GT^{syst-cut} / #LT dN/dy #GT^{def-cut}");
      hYield[k]->SetTitle(Form("%s",lWhichSystVar.Data()));
      hYield[k]->GetYaxis()->SetTitleSize(0.05);
      hYield[k]->GetYaxis()->SetTitleOffset(1.3);
      hYield[k]->GetXaxis()->SetTitleSize(0.05);
      hYield[k]->GetXaxis()->SetTitleOffset(0.8);
      hYield[k]->GetXaxis()->SetLabelSize(0.05);
      hYield[k]->SetMarkerStyle(20);
      hYield[k]->SetMarkerSize(1.5);
    }
    hYield[1]->SetStats(kFALSE);
    //hYield[1]->Draw();
    for (int k = 1; k < 5; k++){
      hYield[k]->Draw("SAME");
    }
    legend->SetTextSize(0.035);
    legend->Draw("SAME");    
    labbig->DrawLatex(0.35, 0.86, Form("%s [%.0f-%0.f]%s",fWhichFixedEstimator.Data(),lFixedLo,lFixedHi,"%"));

    TLatex *xlabel = new TLatex();
    //xlabel->SetTextFont(42);
    xlabel-> SetNDC();
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.06);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    if (lFixedLo2 == 0 && lFixedHi2 == 100.) {
      labbig-> DrawLatex(0.45, 0.75, Form("MB (INEL>0)")); 
    } else labbig-> DrawLatex(0.45, 0.75, Form("applied to %s [%.0f-%.0f] points",fWhichFixedEstimator.Data(),lFixedLo2,lFixedHi2)); 
    xlabel-> DrawLatex(0.35, 0.65, Form("#%s",lCascType.Data()));  
     
    
    cutvar->SaveAs(Form("imagesAN/yieldsyst/fixed%s/Yields-%s-CutVar-%s-%s_%.0f%.0f.png",fWhichFixedEstimator.Data(),fWhichParticle.Data(),lWhichSystVar.Data(),fWhichFixedEstimator.Data(),lFixedLo2,lFixedHi2));

    // Write on file, update
    TFile* lOutputFile = TFile::Open(filename, "UPDATE");
    TDirectoryFile *lDirCutVar = new TDirectoryFile(Form("%sVariationCuts%s",fWhichVarEstimator.Data(),lWhichSystVar.Data()),Form("Variation Cuts %s",lWhichSystVar.Data()));
    lDirCutVar->cd();
    for (int k = 1; k < 5; k++){
      hYield[k]->Write();
    }
    //lOutputFile->cd();
    lOutputFile->Close();
    for (int i = 0; i < binnumber; i++){ 
      for (int j = 1; j < 5; j++){
        lResultsFile[i][j]->Close();
      }
    }

    return hMaxDev;
  }

//--------------------------------------------------------------------------------------------

void beautifygraph(TGraphErrors* g){
    g->GetYaxis()->SetTitleSize(0.05);
    g->GetYaxis()->SetTitleOffset(1.1);
    g->GetXaxis()->SetTitleSize(0.05);
    g->GetXaxis()->SetTitleOffset(1.1);
    g->SetTitle(""); 
}

TH1 *
YieldMean(TH1 *hstat, TH1 *hsys, TF1 *f = NULL, Double_t minfit=.4,Double_t maxfit=8., TString logfilename="log.root",Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q")
{
  if(maxfit>max)
    max=maxfit; 
  if(minfit<min)
    min=minfit; 


  /* set many iterations when fitting the data so we don't
     stop minimization with MAX_CALLS */
  TVirtualFitter::SetMaxIterations(1000000);

  /* create output histo */
  Double_t integral, mean, extra;
  TH1 *hout = new TH1D("hout", "", 9, 0, 9);
  TH1 *hlo, *hhi;
  
  /* create histo with stat+sys errors */
  TH1 *htot = (TH1 *)hstat->Clone(Form("%sfittedwith%s",hstat->GetName(),"LevyTsallis"));
  for (Int_t ibin = 0; ibin < htot->GetNbinsX(); ibin++) {
    htot->SetBinError(ibin + 1, TMath::Sqrt(hsys->GetBinError(ibin + 1) * hsys->GetBinError(ibin + 1) + hstat->GetBinError(ibin + 1) * hstat->GetBinError(ibin + 1)));
  }

  /*
   *   measure the central value 
   */
  Int_t fitres;
  Int_t trials = 0;
  trials = 0;
  do {
    fitres = htot->Fit(f, opt,"",minfit,maxfit);
    Printf("Trial: %d", trials++);
    if(trials > 20) {
      Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
      break;
    }
  }
  while (fitres != 0);
  TFile* filewithfits=TFile::Open(logfilename.Data(),"UPDATE");
  htot->Write();
  //f->SetName(Form("Levyfitto%s",hstat->GetName()));
  f->Write();
  filewithfits->Close();    
  delete filewithfits;   
  
  cout<<" Fit sys+stat for " <<f->GetName()<<endl;    
  cout<<"NDF="<<f->GetNDF()<<" Chi^2="<<f->GetChisquare()<<" Chi^2/NDF="<<f->GetChisquare()/f->GetNDF()<<endl;
  hout->SetBinContent(kChi, f->GetChisquare()/f->GetNDF());

  hlo = YieldMean_LowExtrapolationHisto(htot, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(htot, f, max, hiprecision);
  YieldMean_IntegralMean(htot, hlo, hhi, integral, mean, extra, kTRUE);
  hout->SetBinContent(kYield, integral);
  hout->SetBinContent(kMean, mean);
  hout->SetBinContent(kExtra, extra);

  /*
   * STATISTICS
   */
  
  TCanvas *cCanvasStat = new TCanvas("cCanvasStat");
  cCanvasStat->Divide(2, 1);
  
  /*
   * measure statistical error
   */

  /* fit with stat error */
  trials = 0;
  do {
    fitres = hstat->Fit(f, opt,"",minfit,maxfit);
    Printf("Trial: %d", trials++);
    if(trials > 10) {
      Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
      break;
    }
  }
  while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hstat, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hstat, f, max, hiprecision);
  
  /* random generation with integration (coarse) */
  TH1 *hIntegral_tmp = new TH1F("hIntegral_tmp", "", 1000, 0.75 * integral, 1.25 * integral);
  TH1 *hMean_tmp = new TH1F("hMean_tmp", "", 1000, 0.75 * mean, 1.25 * mean);
  for (Int_t irnd = 0; irnd < 100; irnd++) {
    /* get random histogram */
    TH1 *hrnd = YieldMean_ReturnRandom(hstat);
    /* fit */
    TH1 *hrndlo = YieldMean_ReturnCoherentRandom(hlo);
    TH1 *hrndhi = YieldMean_ReturnCoherentRandom(hhi);
    /* integrate */
    YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean, extra);
    hIntegral_tmp->Fill(integral);
    hMean_tmp->Fill(mean);
    delete hrnd;
    delete hrndlo;
    delete hrndhi;
   }
  /* random generation with integration (fine) */
  TH1 *hIntegral = new TH1F("hIntegral", "", 100, 
                            hIntegral_tmp->GetMean() - 10. * hIntegral_tmp->GetRMS(),
                            hIntegral_tmp->GetMean() + 10. * hIntegral_tmp->GetRMS());
  TH1 *hMean = new TH1F("hMean", "", 100,
                        hMean_tmp->GetMean() - 10. * hMean_tmp->GetRMS(),
                        hMean_tmp->GetMean() + 10. * hMean_tmp->GetRMS());
  for (Int_t irnd = 0; irnd < 1000; irnd++) {
    /* get random histogram */
    TH1 *hrnd = YieldMean_ReturnRandom(hstat);
    /* fit */
    TH1 *hrndlo = YieldMean_ReturnCoherentRandom(hlo);
    TH1 *hrndhi = YieldMean_ReturnCoherentRandom(hhi);
    /* integrate */
    YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean, extra);
    hIntegral->Fill(integral);
    hMean->Fill(mean);
    delete hrnd;
    delete hrndlo;
    delete hrndhi;
  }
  TF1 *gaus = (TF1 *)gROOT->GetFunction("gaus");
  
  cCanvasStat->cd(1);
  hIntegral->Fit(gaus, "q");
  integral = hout->GetBinContent(kYield) * gaus->GetParameter(2) / gaus->GetParameter(1);
  hout->SetBinContent(kYieldStat, integral);
  
  cCanvasStat->cd(2);
  hMean->Fit(gaus, "q");
  mean = hout->GetBinContent(kMean) * gaus->GetParameter(2) / gaus->GetParameter(1);
  hout->SetBinContent(kMeanStat, mean);
  
  /*
   * SYSTEMATICS
   */

  TCanvas *cCanvasSys = new TCanvas("cCanvasYieldSys");
  cCanvasSys->Divide(2, 1);
  cCanvasSys->cd(1)->DrawFrame(min, 1.e-3, max, 1.e3);
  hsys->SetMarkerStyle(20);
  hsys->SetMarkerColor(1);
  hsys->SetMarkerSize(1);
  hsys->Draw("same");
  cCanvasSys->cd(2)->DrawFrame(min, 1.e-3, max, 1.e3);
  hsys->Draw("same");
  
  /*
   * systematic error high
   */

  TH1 *hhigh = YieldMean_ReturnExtremeHighHisto(hsys);
  trials = 0;
  do {
    fitres = hhigh->Fit(f, opt,"",minfit,maxfit);
    Printf("Trial: %d", trials++);
    if(trials > 10) {
      Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
      break;
    }
  }
  while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hhigh, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hhigh, f, max, hiprecision);
  YieldMean_IntegralMean(hhigh, hlo, hhi, integral, mean, extra);
  integral = TMath::Abs(integral - hout->GetBinContent(kYield));
  hout->SetBinContent(kYieldSysHi, integral);

  cCanvasSys->cd(1);
  f->SetLineColor(2);
  f->DrawCopy("same");
  
  /*
   * systematic error hard
   */

  TH1 *hhard = YieldMean_ReturnExtremeHardHisto(hsys);
  trials = 0;
  do {
    fitres = hhard->Fit(f, opt,"",minfit,maxfit);
    Printf("Trial: %d", trials++);
    if(trials > 10) {
      Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
      break;
    }
  }
  while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hhard, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hhard, f, max, hiprecision);
  YieldMean_IntegralMean(hhard, hlo, hhi, integral, mean, extra);
  mean = TMath::Abs(mean - hout->GetBinContent(kMean));
  hout->SetBinContent(kMeanSysHi, mean);

  cCanvasSys->cd(2);
  f->SetLineColor(2);
  f->DrawCopy("same");
  
  /*
   * systematic error low
   */

  TH1 *hlow = YieldMean_ReturnExtremeLowHisto(hsys);
  trials = 0;
  do {
    fitres = hlow->Fit(f, opt,"",minfit,maxfit);
    Printf("Trial: %d", trials++);
    if(trials > 10) {
      Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
      break;
    }
  }
  while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hlow, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hlow, f, max, hiprecision);
  YieldMean_IntegralMean(hlow, hlo, hhi, integral, mean, extra);
  integral = TMath::Abs(integral - hout->GetBinContent(kYield));
  hout->SetBinContent(kYieldSysLo, integral);

  cCanvasSys->cd(1);
  f->SetLineColor(4);
  f->DrawCopy("same");

  /*
   * systematic error soft
   */

  TH1 *hsoft = YieldMean_ReturnExtremeSoftHisto(hsys);
  trials = 0;
  do {
    fitres = hsoft->Fit(f, opt,"",minfit,maxfit);
    Printf("Trial: %d", trials++);
    if(trials > 10) {
      Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
      break;
    }
  }
  while (fitres != 0);
  hlo = YieldMean_LowExtrapolationHisto(hsoft, f, min, loprecision);
  hhi = YieldMean_HighExtrapolationHisto(hsoft, f, max, hiprecision);
  YieldMean_IntegralMean(hsoft, hlo, hhi, integral, mean, extra);
  mean = TMath::Abs(mean - hout->GetBinContent(kMean));
  hout->SetBinContent(kMeanSysLo, mean);

  cCanvasSys->cd(2);
  f->SetLineColor(4);
  f->DrawCopy("same");

  return hout;
}

TH1 *
YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, Double_t min, Double_t binwidth)
{
  /* find lowest edge in histo */
  Int_t binlo;
  Double_t lo;
  for (Int_t ibin = 1; ibin < h->GetNbinsX() + 1; ibin++) {
    if (h->GetBinContent(ibin) != 0.) {
      binlo = ibin;
      lo = h->GetBinLowEdge(ibin);
      break;
    }
  }
  
  Int_t nbins = (lo - min) / binwidth;
  if(nbins<1)
  return 0x0;   
  TH1 *hlo = new TH1F("hlo", "", nbins, min, lo);
  
  /* integrate function in histogram bins */
  Double_t cont, err, width;
  for (Int_t ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
    width = hlo->GetBinWidth(ibin + 1);
    cont = f->Integral(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2));//, (Double_t *)0, 1.e-6);
    err = f->IntegralError(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), (Double_t *)0, (Double_t *)0, 1.e-6);
    hlo->SetBinContent(ibin + 1, cont / width);
    hlo->SetBinError(ibin + 1, err / width);
  }

  return hlo;
}

TH1 *
YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, Double_t max, Double_t binwidth)
{
  /* find highest edge in histo */
  Int_t binhi;
  Double_t hi;
  for (Int_t ibin = h->GetNbinsX(); ibin > 0; ibin--) {
    if (h->GetBinContent(ibin) != 0.) {
      binhi = ibin + 1;
      hi = h->GetBinLowEdge(ibin + 1);
      break;
    }
  }
  if(max<hi) {
  Printf("Warning! You should probably set a higher max value (Max = %f, hi = %f)", max, hi);
    return 0x0;
  }
  Int_t nbins = (max - hi) / binwidth;
 if(nbins<1)
  return 0x0;  
  TH1 *hhi = new TH1F("hhi", "", nbins, hi, max);
  
  /* integrate function in histogram bins */
  Double_t cont, err, width;
  for (Int_t ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
    width = hhi->GetBinWidth(ibin + 1);
    cont = f->Integral(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2));//, (Double_t *)0, 1.e-6);
    err = f->IntegralError(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), (Double_t *)0, (Double_t *)0, 1.e-6);
    hhi->SetBinContent(ibin + 1, cont / width);
    hhi->SetBinError(ibin + 1, err / width);
  }

  return hhi;
}

TH1 *
YieldMean_ReturnRandom(TH1 *hin)
{
  TH1 *hout = (TH1 *)hin->Clone("hout");
  hout->Reset();
  Double_t cont, err;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    cont = hin->GetBinContent(ibin + 1);
    err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, gRandom->Gaus(cont, err));
    hout->SetBinError(ibin + 1, err);
  }
  return hout;
}

TH1 *
YieldMean_ReturnCoherentRandom(TH1 *hin)
{
 if(!hin)
  return 0x0;   
  TH1 *hout = (TH1 *)hin->Clone("hout");
  hout->Reset();
  Double_t cont, err, cohe;
  cohe = gRandom->Gaus(0., 1.);
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    cont = hin->GetBinContent(ibin + 1);
    err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, cont + cohe * err);
    hout->SetBinError(ibin + 1, err);
  }
  return hout;
}

TH1 *
YieldMean_ReturnExtremeHighHisto(TH1 *hin)
{
  TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehigh", hin->GetName()));
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    Double_t val = hin->GetBinContent(ibin + 1);
    Double_t err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, val + err);
  }
  return hout;
}

TH1 *
YieldMean_ReturnExtremeLowHisto(TH1 *hin)
{
  TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremelow", hin->GetName()));
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    Double_t val = hin->GetBinContent(ibin + 1);
    Double_t err = hin->GetBinError(ibin + 1);
    hout->SetBinContent(ibin + 1, val - err);
  }
  return hout;
}

TH1 *
YieldMean_ReturnExtremeSoftHisto(TH1 *hin)
{
  return YieldMean_ReturnExtremeHisto(hin, -1.);
}

TH1 *
YieldMean_ReturnExtremeHardHisto(TH1 *hin)
{
  return YieldMean_ReturnExtremeHisto(hin, 1.);
}

TH1 *
YieldMean_ReturnExtremeHisto(TH1 *hin, Float_t sign)
{
  Double_t ptlow, pthigh;
  for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    ptlow = hin->GetBinLowEdge(ibin + 1);
    break;
  }
  for (Int_t ibin = hin->GetNbinsX(); ibin >= 0; ibin--) {
    if (hin->GetBinError(ibin + 1) <= 0.) continue;
    pthigh = hin->GetBinLowEdge(ibin + 2);
    break;
  }

  Double_t mean = hin->GetMean();
  Double_t maxdiff = 0.;
  TH1 *hmax = NULL;
  for (Int_t inode = 0; inode < hin->GetNbinsX(); inode++) {

    Double_t ptnode = hin->GetBinCenter(inode + 1);
    TH1 *hout = (TH1 *)hin->Clone(Form("%s_extremehard", hin->GetName()));
    
    for (Int_t ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      Double_t val = hin->GetBinContent(ibin + 1);
      Double_t err = hin->GetBinError(ibin + 1);
      Double_t cen = hin->GetBinCenter(ibin + 1);
      if (cen < ptnode)
        err *= -1. + (cen - ptlow) / (ptnode - ptlow);
      else
        err *= (cen - ptnode) / (pthigh - ptnode);

      hout->SetBinContent(ibin + 1, val + sign * err);
    }

    Double_t diff = TMath::Abs(mean - hout->GetMean());
    if (diff > maxdiff) {
      //      printf("found max at %f\n", ptnode);
      if (hmax) delete hmax;
      hmax = (TH1 *)hout->Clone("hmax");
      maxdiff = diff;
    }
    delete hout;
  }
  return hmax;
}

void YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, Double_t &integral, Double_t &mean, Double_t &extra, Bool_t printinfo)
{
  
  /*
   * compute integrals
   */
  
  Double_t cont, err, width, cent;
  Double_t I = 0., IX = 0., Ierr = 0., IXerr = 0., Ilerr = 0., IXlerr = 0.;
  Double_t M = 0., Merr = 0., Mlerr = 0., C;
  Double_t E = 0;
  Double_t dataonly=0.0;

  /* integrate the data */
  for (Int_t ibin = 0; ibin < hdata->GetNbinsX(); ibin++) {
    cent = hdata->GetBinCenter(ibin + 1);
    width = hdata->GetBinWidth(ibin + 1);
    cont = width * hdata->GetBinContent(ibin + 1);
    err = width * hdata->GetBinError(ibin + 1);
    if (err <= 0.) continue;
    I += cont;
    IX += cont * cent;
  }
  
  dataonly=I; 
  /* integrate low */
  if(hlo)
  { 
    for (Int_t ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
      cent = hlo->GetBinCenter(ibin + 1);
      width = hlo->GetBinWidth(ibin + 1);
      cont = width * hlo->GetBinContent(ibin + 1);
     err = width * hlo->GetBinError(ibin + 1);
     if (err <= 0.) continue;
     I += cont;
     IX += cont * cent;
     E += cont;
   }
 }
  /* integrate high */
 if(printinfo)  
    cout<<"low part data only = "<<dataonly<<" total = "<<I<<" ratio= "<<dataonly/I<<endl;  
 if(hhi)
 {
   for (Int_t ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
     cent = hhi->GetBinCenter(ibin + 1);
     width = hhi->GetBinWidth(ibin + 1);
     cont = width * hhi->GetBinContent(ibin + 1);
     err = width * hhi->GetBinError(ibin + 1);
     if (err <= 0.) continue;
     I += cont;
     IX += cont * cent;
     E += cont;
   }
}
  /* set values */
  integral = I;
  mean = IX / I;
  extra = E;
  if(printinfo) 
    cout<<"low+high data only = "<<dataonly<<" total = "<<I<<" ratio= "<<dataonly/I<<endl;  
}


//---------------------------------------------------------------------
Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{
  
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-5, 1.e5);
  return fLevyTsallis;
}  
//------------------------------------------------------

//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop + errorfrombottom) );
    }
    return 1.;
}

//----------------------------------------------------------------------------------------------------
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 ){ 
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
Bool_t PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin){
  double dev = h->GetBinContent(bin);
  double RBsigma = h->GetBinError(bin);

  if( TMath::Abs(dev-1) > (nsigmas*RBsigma) ) {return kTRUE;}
  else {return kFALSE;}
}