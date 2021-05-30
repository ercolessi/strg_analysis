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
  kExtra
};
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) ;

TH1 *YieldMean(TH1 *hstat, TF1 *f = NULL, Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q",TString logfilename="log.root",Double_t minfit=0.6,Double_t maxfit=6.5);
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

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 );
Bool_t PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin);
TH1F* ComputeYieldsforSystematics( TString filename, TString fWhichEstimator = "V0M",  TString fWhichParticle = "XiMinus",  Double_t lLoMult = 0.,  Double_t lHiMult = 100.,  Double_t lLoEE = 0., Double_t lHiEE = 100.,  TString lWhichSystVar = "V0Radius");

//--------------------------------------------------------------------------------
//----------------------------- MAIN FUNCTION ------------------------------------
//--------------------------------------------------------------------------------
void DoYieldsSystematicsFullXi(
  TString fWhichEstimator = "ZDC", 
  TString fWhichParticle = "XiMinus", 
  Double_t lLoMult = 0., 
  Double_t lHiMult = 100., 
  Double_t lLoEE = 0., 
  Double_t lHiEE = 100.){

  cout<<"  "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"                  Make systematics for yields                   "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<endl;

  Double_t percentileV0[] = {0,30,50,100};
 // if (fWhichEstimator.Contains("ZDC") && lLoMult == 70) Double_t percentileV0[4] = {0.,30,50,100};
  const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
  Double_t percentileZDC[] = {0,40,70,100};
  const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;
  Double_t * percentile;
  int tempbinnumber = 0;
  if (fWhichEstimator.Contains("V0M")) {
      tempbinnumber = nbinV0;
      percentile = percentileV0;
  } 
  else if (fWhichEstimator.Contains("ZDC")) {
      tempbinnumber = nbinZDC;
      percentile = percentileZDC;
  } 
 // else {cout << "No valid name for estimator... its V0M or ZDC" << endl; return;}
  const int binnumber = tempbinnumber;
  TString outputfilename = Form("Systematics-Yield-%s-%s-sel-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE);

  // Selection classes
  TH1F* hV0Radius = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"V0Radius");
  TH1F* hCascRadius = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"CascRadius");
  TH1F* hNegToPV = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"DCANegToPV");
  TH1F* hPosToPV = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"DCAPosToPV");
  TH1F* hV0ToPV = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"DCAV0ToPV");
  TH1F* hBachToPV = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"DCABachToPV");
  TH1F* hV0Daught = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"DCAV0Daughters");
  TH1F* hCascDaught = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"DCACascDaughters");
  TH1F* hV0CosPA = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"V0CosPA");
  TH1F* hCascCosPA = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"CascCosPA");
  TH1F* hV0Mass = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"V0Mass");
  TH1F* hCompetingSpecies = 0x0;
  TH1F* hPLT = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"ProperLifetime");
  TH1F* hTPCNSigmas = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"TPCPIDNSigmas");
  TH1F* hTPCNClusters = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"TPCNClusters");
  TH1F* hSigExtBinCount = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"SigExtBinCount");
  //TH1F* hGeantFluka = ComputeYieldsforSystematics(outputfilename,fWhichEstimator.Data(),fWhichParticle.Data(),lLoMult,lHiMult,lLoEE,lHiEE,"GeantFlukaCorrection");
 
  // Write Max Dev on output File
  TFile* Write = new TFile (outputfilename, "UPDATE");
 /* hV0Radius->Write();
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
  hSigExtBinCount->Write();*/
  //hGeantFluka->Write();

}

//--------------------------------------------------------------------------------
//----------------------------- OTHER FUNCTIONS ----------------------------------
//--------------------------------------------------------------------------------

TH1F* ComputeYieldsforSystematics(
  TString filename,
  TString fWhichEstimator = "ZDC", 
  TString fWhichParticle = "XiMinus", 
  Double_t lLoMult = 0., 
  Double_t lHiMult = 100., 
  Double_t lLoEE = 0., 
  Double_t lHiEE = 100.,
  TString lWhichSystVar = "V0Radius") {

    cout << "\n---------------------------------------------------------------" << endl;
    cout << endl << "\n---> Do Systematics for: "<< lWhichSystVar << endl;
    cout << "\n---------------------------------------------------------------" << endl;

    TString fWhichAntiParticle = "XiPlus";

    Double_t percentileV0[] = {0,30,50,100};
    //if (fWhichEstimator.Contains("ZDC") && lLoMult == 70) percentileV0 = {0.,30,50,100};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    Double_t percentileZDC[] = {0,40,70,100};
    const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;
    Double_t percentileMB[] = {0,100};
    const Long_t nbinMB = sizeof(percentileMB)/sizeof(Double_t) - 1;
    Double_t * percentile;
    int tempbinnumber = 0;
    TString estimator = "", selection = "",estimatorevt = "";
    TFile* filenormcorr;
    if (fWhichEstimator.Contains("V0M")) {
        tempbinnumber = nbinV0;
        percentile = percentileV0;
        estimator = "V0";
        selection = "multsel";
        filenormcorr = TFile::Open("../NormalizationCorrections/EventCountLoss_16d3.root");
        if (lLoEE > 0) {
          estimator = "V0FixLowEE";
          selection = "multsel_fixedlowEE";
          filenormcorr = TFile::Open("../NormalizationCorrections/EventCountLoss_16d3_ZDC70100.root");
        }
        if (lHiEE < 100) {
          estimator = "V0FixHighEE";
          selection = "multsel_fixedhighEE";   
          filenormcorr = TFile::Open("../NormalizationCorrections/EventCountLoss_16d3_ZDC030.root");
        }
    } 
    else if (fWhichEstimator.Contains("ZDC")) {
        tempbinnumber = nbinZDC;
        percentile = percentileZDC;
        estimator = "ZDC";
        selection = "EEsel";
        estimatorevt = estimator;
        filenormcorr = TFile::Open("../NormalizationCorrections/EventCountLoss_16d3_ZDC030.root");
        if (lLoMult > 0) {
          estimator = "ZDCFixLowmult";
          estimatorevt = "V0FixLowmult";
          selection = "EEsel_fixedlowmult";
          filenormcorr = TFile::Open("../NormalizationCorrections/EventCountLoss_16d3_ZDC030.root");
        }
        if (lHiMult < 100) {
          estimator = "ZDCFixHighmult";
          estimatorevt = "V0FixHighmult";
          selection = "EEsel_fixedhighmult";
          filenormcorr = TFile::Open("../NormalizationCorrections/EventCountLoss_16d3_ZDC030.root");
        }

    } 
    else if (fWhichEstimator.Contains("MB")) {
        tempbinnumber = nbinMB;
        percentile = percentileMB;
        lLoMult = 0.;
        lHiMult = 100.;
        lLoEE = 0.;
        lHiEE = 100.;
        estimator = "V0";
        selection = "multsel";
        filenormcorr = TFile::Open("../NormalizationCorrections/EventCountLoss_16d3_MB.root");
    } 
    else {cout << "No valid name for estimator... its V0M, ZDC or MB" << endl;}
    const int binnumber = tempbinnumber;

    TH1F* hEevt;
    if (fWhichEstimator.Contains("ZDC")) {hEevt = (TH1F*)filenormcorr->Get(Form("EventLoss/hevtloss%s",estimatorevt.Data()));}
    else {hEevt = (TH1F*)filenormcorr->Get(Form("EventLoss/hevtloss%s",estimator.Data()));};
    TH1F* hsgnloss[binnumber];
    for (int i =0; i<binnumber;i++){
        hsgnloss[i] = (TH1F*)filenormcorr->Get(Form("SgnLoss/%s/%s/fHistptSel%s_%.0f-%.0f_%s",
                                        "XiMinus",selection.Data(),estimator.Data(),percentile[i],percentile[i+1],"XiMinus"));
    }

    // Histos
    TH1D* lHistPt[binnumber][5];
    TH1D* lSystPt[binnumber][5];
    TH1D* lAntiHistPt[binnumber][5];
    TH1D* lAntiSystPt[binnumber][5];
    TH1D* hout[binnumber][5];
    Double_t Yield[binnumber][5];
    Double_t YieldStat[binnumber][5];
  
    // Fit function
    TF1* LevyTsallisfunc = LevyTsallis("LevyTsallisfunc", 1.321);

    // Get spectra for this cut
    TFile* lResultsFile[binnumber][5];
    TFile* lAntiResultsFile[binnumber][5];
    for (int i = 0; i < binnumber; i++){ //loop over selection classes           
      double zdcmin = lLoEE, zdcmax = lHiEE, v0min = lLoMult, v0max = lHiMult;
      if (fWhichEstimator.Contains("V0M")) {
        v0min = percentile[i];
        v0max = percentile[i+1];
      }
      if (fWhichEstimator.Contains("ZDC")) {
        zdcmin = percentile[i];
        zdcmax = percentile[i+1];
      }   
      //
      TString lSystFile = "Results/Results-Systematics";
      lSystFile.Append( Form( "-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f-", fWhichParticle.Data(), v0min, v0max, zdcmin, zdcmax) );
      TString lAntiSystFile = "Results/Results-Systematics";
      lAntiSystFile.Append( Form( "-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f-", fWhichAntiParticle.Data(), v0min, v0max, zdcmin, zdcmax) );
      //
      lResultsFile[i][0] = TFile::Open(Form( "Results/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichParticle.Data(), v0min, v0max, zdcmin, zdcmax));
      lAntiResultsFile[i][0] = TFile::Open(Form( "Results/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichAntiParticle.Data(), v0min, v0max, zdcmin, zdcmax));
      for (int j = 1; j < 5; j++){
        //Particle
        TString lDataFilename = lSystFile + lWhichSystVar + Form("-%i.root",j);
        lResultsFile[i][j] = TFile::Open(lDataFilename.Data());
        lHistPt[i][j] = (TH1D*)lResultsFile[i][j]->Get(Form("fHistPt%s", fWhichParticle.Data()));
        //Anti Particle
        TString lAntiDataFilename = lAntiSystFile + lWhichSystVar + Form("-%i.root",j);
        lAntiResultsFile[i][j] = TFile::Open(lAntiDataFilename.Data());
        lAntiHistPt[i][j] = (TH1D*)lAntiResultsFile[i][j]->Get(Form("fHistPt%s", fWhichAntiParticle.Data()));
      }  
      lHistPt[i][0] = (TH1D*)lResultsFile[i][0]->Get(Form("fHistPt%s", fWhichParticle.Data()));
      lAntiHistPt[i][0] = (TH1D*)lAntiResultsFile[i][0]->Get(Form("fHistPt%s", fWhichAntiParticle.Data()));
      // Add Spectra
      for (int j = 0; j < 5; j++){
        for (int bin = 1; bin <=  lHistPt[i][j]->GetNbinsX(); bin ++){
          lHistPt[i][j]->AddBinContent(bin, lAntiHistPt[i][j]->GetBinContent(bin));
          lHistPt[i][j]->SetBinError(bin, TMath::Sqrt(lHistPt[i][j]->GetBinError(bin)*lHistPt[i][j]->GetBinError(bin)+
                                                    lAntiHistPt[i][j]->GetBinError(bin)*lAntiHistPt[i][j]->GetBinError(bin)));
        }
      }
      //correct for sgn loss and event loss
      double sgn = 0, sgne = 0, evt = 0;
      for (int j = 0; j < 5; j++){
        for (int b = 1; b <= lHistPt[i][0]->GetNbinsX(); b++){

          sgn = hsgnloss[i]->GetBinContent(b+1);
          sgne = hsgnloss[i]->GetBinError(b+1);
          evt = hEevt->GetBinContent(i+1);
        
          if (sgn!=0){
            lHistPt[i][j]->SetBinContent(b,lHistPt[i][j]->GetBinContent(b)*evt/sgn);                    
            lHistPt[i][j]->SetBinError(b,ErrorInRatio(lHistPt[i][j]->GetBinContent(b),lHistPt[i][j]->GetBinError(b),
                                    sgn,sgne));
          }
        }
      }
    }//end loop over classes 

    //Correct for sgn loss and event loss


    
    // Get the yields
    for (int i = 0; i < binnumber; i++){ //loop over selection classes  
      for (int j = 0; j < 5; j++){
        hout[i][j] = (TH1D*)YieldMean(lHistPt[i][j],LevyTsallisfunc);
        // Get yield
        Yield[i][j] = hout[i][j]->GetBinContent(1);
        YieldStat[i][j] = hout[i][j]->GetBinContent(2);
      }
    }//end loop over classes 

    TH1F* hYield[5];
    for (int j = 0; j < 5; j++){
      hYield[j] = new TH1F(Form("hYield%i-%s",j,fWhichEstimator.Data()), Form(";percentile %s (%);dN/dy",fWhichEstimator.Data()),binnumber,percentile);
      for (int nbin = 1; nbin <= hYield[j]->GetNbinsX(); nbin++){
        hYield[j]->SetBinContent(nbin, Yield[nbin-1][j]);
        hYield[j]->SetBinError(nbin, YieldStat[nbin-1][j]);
      }
    }

    // Divide following Roger-Barlow prescription
    for (int j = 1; j < 5; j++){
      cout << hYield[j]->GetBinContent(3) << "   " << hYield[0]->GetBinContent(3) << endl;
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
    hMaxDev->SetName(Form("hMaxDev%s%s",lWhichSystVar.Data(),fWhichEstimator.Data()));

    // Prepare Canvases
    //Max Deviation
    hMaxDev->GetXaxis()->SetRangeUser(0.,100.);
    hMaxDev->GetYaxis()->SetRangeUser(-0.0005,.1);
    hMaxDev->SetYTitle("max rel. dev.");
    hMaxDev->SetTitle(Form("%s",lWhichSystVar.Data()));
    hMaxDev->GetYaxis()->SetTitleSize(0.05);
    hMaxDev->GetYaxis()->SetTitleOffset(1.0);
    hMaxDev->GetXaxis()->SetTitleSize(0.05);
    hMaxDev->GetXaxis()->SetTitleOffset(1.);
    hMaxDev->SetMarkerStyle(20);
    hMaxDev->SetMarkerSize(1.1);
    hMaxDev->SetMarkerColor(kBlack);
    hMaxDev->SetLineColor(kBlack);
    TCanvas* maxdev = new TCanvas("maxdev"," ",1000,800);
    maxdev->SetRightMargin(0.15);
    maxdev->SetLeftMargin(0.15);
    maxdev->SetBottomMargin(0.15);
    maxdev->SetGridy();
    maxdev->SetGridx();
    hMaxDev->SetStats(kFALSE);  
    hMaxDev->Draw();
    maxdev->SaveAs(Form("images/yields/%s/Yields-%s-MaxRelDev%s-%s.png",fWhichEstimator.Data(),fWhichParticle.Data(),lWhichSystVar.Data(),fWhichEstimator.Data()));

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
      hYield[k]->SetName(Form("hYield-%s-%s-%i",lWhichSystVar.Data(),fWhichEstimator.Data(),k));
      hYield[k]->SetTitle(Form("hYield-%s-%i",lWhichSystVar.Data(),k));
      hYield[k]->GetXaxis()->SetRangeUser(0.,100.);
      hYield[k]->GetYaxis()->SetRangeUser(0.9,1.1);
      hYield[k]->SetYTitle("Yield^{syst-cut} / Yield^{def-cut}");
      hYield[k]->SetTitle(Form("%s",lWhichSystVar.Data()));
      hYield[k]->GetYaxis()->SetTitleSize(0.05);
      hYield[k]->GetYaxis()->SetTitleOffset(1.1);
      hYield[k]->GetXaxis()->SetTitleSize(0.05);
      hYield[k]->GetXaxis()->SetTitleOffset(0.8);
      hYield[k]->SetMarkerStyle(20);
      hYield[k]->SetMarkerSize(1.1);
    }
    hYield[1]->SetStats(kFALSE);
    hYield[1]->Draw();
    for (int k = 2; k < 5; k++){
      hYield[k]->Draw("SAME");
    }
    legend->SetTextSize(0.035);
    legend->Draw("SAME");    
    cutvar->SaveAs(Form("images/yields/%s/Yields-%s-CutVar-%s-%s.png",fWhichEstimator.Data(),fWhichParticle.Data(),lWhichSystVar.Data(),fWhichEstimator.Data()));

    // Write on file, update
    TFile* lOutputFile = TFile::Open(filename, "UPDATE");
    TDirectoryFile *lDirCutVar = new TDirectoryFile(Form("%sVariationCuts%s",fWhichEstimator.Data(),lWhichSystVar.Data()),Form("Variation Cuts %s",lWhichSystVar.Data()));
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
TH1 *
YieldMean(TH1 *hstat, TF1 *f = NULL, Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q",TString logfilename="log.root",Double_t minfit=0.6,Double_t maxfit=6.5)
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
  
  /* create histo with to fit */
  TH1 *htot = (TH1 *)hstat->Clone(Form("%sfittedwith%s",hstat->GetName(),f->GetName()));

  /*
   *   measure the central value 
   */
  Int_t fitres;
  Int_t trials = 0;
  trials = 0;
  do {
    fitres = htot->Fit(f, opt,"",minfit,maxfit);
    Printf("Trial: %d", trials++);
    if(trials > 10) {
      Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
      break;
    }
  }
  while (fitres != 0);
  TFile* filewithfits=TFile::Open(logfilename.Data(),"UPDATE");
  htot->Write();
  filewithfits->Close();    
  delete filewithfits;   
  
  cout<<" Fit sys+stat for " <<f->GetName()<<endl;    
  cout<<"NDF="<<f->GetNDF()<<" Chi^2="<<f->GetChisquare()<<" Chi^2/NDF="<<f->GetChisquare()/f->GetNDF()<<endl;

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
  hIntegral->Fit(gaus, "");
  integral = hout->GetBinContent(kYield) * gaus->GetParameter(2) / gaus->GetParameter(1);
  hout->SetBinContent(kYieldStat, integral);
  
  cCanvasStat->cd(2);
  hMean->Fit(gaus, "");
  mean = hout->GetBinContent(kMean) * gaus->GetParameter(2) / gaus->GetParameter(1);
  hout->SetBinContent(kMeanStat, mean);

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
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e6);
  return fLevyTsallis;
}  
//------------------------------------------------------

//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop - errorfrombottom) );
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