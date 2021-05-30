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
  kExtra,
  kChi

};
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) ;

TH1 *YieldMean(TH1 *hstat, TH1 *hsys, TF1 *f = NULL, Double_t minfit=0.6,Double_t maxfit=6.5, Double_t min = 0., Double_t max = 100., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q",TString logfilename="log.root");
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

void DoMeanPt(
  TString fWhichEstimator = "ZDC", 
  Double_t lLoMult = 0., 
  Double_t lHiMult = 100., 
  Double_t lLoEE = 0., 
  Double_t lHiEE = 100.) {
 
    TString fWhichOtherEstimator = "V0M";
    if (fWhichEstimator.Contains("V0M")) fWhichOtherEstimator = "ZDC";
    //
    Double_t percentileV0[] = {0.,1.,5.,10.,15.,20.,30.,40.,50.,70.,100.};
    Double_t systpercentileV0[4] = {0.,15.,40.,100.};
    Int_t systcounterV0[3] = {5,8,10};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    //
    Double_t percentileZDC[] = {0.,20.,30.,40.,50.,60.,70.,80.,90.,100.};
    Double_t systpercentileZDC[4] = {0.,40.,70.,100.};
    Int_t systcounterZDC[3] = {4,7,9};
    const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;
    //
    Double_t * percentile, *systpercentile, *dNch, *PosSystdNch, *NegSystdNch, *StatdNch;
    Int_t *systcounter;
    double dNchMB = 6.89;
    
    Double_t dNchZDC[nbinZDC] = { 10.07, 8.95, 7.88, 7.00, 6.26, 5.70, 5.18, 4.61, 3.79 };        
    Double_t PosSystdNchZDC[nbinZDC] = {0.12,0.11,0.10,0.09,0.08,0.07,0.07,0.06,0.05};
    Double_t NegSystdNchZDC[nbinZDC] = {0.08,0.07,0.07,0.06,0.05,0.05,0.04,0.04,0.03};
    Double_t StatdNchZDC[nbinZDC] = {0.};
    //OFFICIAL FOR V0M 
    Double_t dNchV0[nbinV0] = {25.75, 19.83, 16.12, 13.76, 12.06, 10.11, 8.07, 6.48, 4.64, 2.52};     
    Double_t PosSystdNchV0[nbinV0] = { +0.35, +0.27, +0.22, +0.19, +0.17, +0.14, +0.11, +0.09, +0.07, +0.04 };
    Double_t NegSystdNchV0[nbinV0] = { 0.29, 0.22, 0.18, 0.16, 0.14, 0.11, 0.09, 0.07, 0.05, +0.03 };
    Double_t StatdNchV0[nbinV0] = {0.};
    //
    Double_t dNchZDC_V0030[nbinZDC] = { 15.84,14.55,13.80,13.20,12.70,12.24,11.78,11.29,10.72  };
    Double_t PosSystdNchZDC_V0030[nbinZDC] = { 0.19 , 0.17 ,0.17 ,0.16 ,0.15 ,0.15 ,0.14 ,0.13 ,0.12 };
    Double_t NegSystdNchZDC_V0030[nbinZDC] = { 0.19 , 0.17 ,0.17 ,0.16 ,0.15 ,0.15 ,0.14 ,0.13 ,0.12 };
     
    Double_t OffdNchMB = 6.89;
    Double_t OffSystdNchMB = 0.11;
    TString namesystsgnloss = ""; 
    TString nameMCVaried = "";

    // Choose scenario and initialize correct variables
    int tempbinnumber = 0;
    if (fWhichEstimator.Contains("V0M")) {
        tempbinnumber = nbinV0;
        dNch = dNchV0;
        PosSystdNch = PosSystdNchV0;
        NegSystdNch = NegSystdNchV0;
        StatdNch = StatdNchV0;
        nameMCVaried = "hV0MZDC030";
        percentile = percentileV0;
        systpercentile = systpercentileV0;
        systcounter = systcounterV0;
    } 
    else if (fWhichEstimator.Contains("ZDC")) {
        tempbinnumber = nbinZDC;
        dNch = dNchZDC;
        PosSystdNch = PosSystdNchZDC;
        NegSystdNch = NegSystdNchZDC;
        StatdNch = StatdNchZDC;
        namesystsgnloss = Form("%s-%s","EEsel","ZDC");
        nameMCVaried = "hZDC";
        if (lHiMult < 100) {
          dNch = dNchZDC_V0030;
          PosSystdNch = PosSystdNchZDC_V0030;
          NegSystdNch = NegSystdNchZDC_V0030;
          namesystsgnloss = Form("%s-%s","EEsel_fixedhighmult","ZDCFixHighmult");
          nameMCVaried = "hZDCV0M030";
        }
        percentile = percentileZDC;
        systpercentile = systpercentileZDC;
        systcounter = systcounterZDC;
    } 
    else {cout << "No valid name for estimator... its V0M or ZDC" << endl; return;}
    const int binnumber = tempbinnumber;
    //    
    double errpercentile[binnumber], centrpercentile[binnumber];
    for (Int_t n = 0; n < binnumber; n++) {
      centrpercentile[n] = (percentile[n] + percentile[n + 1]) / 2;
      errpercentile[n] = (percentile[n + 1] - percentile[n]) / 2;
    }
    //
    // Histos
    TH1D* lHistPt[binnumber];
    TH1D* lSystPt[binnumber];
    TH1D* hout[binnumber];
    TH1D* hSystSpectra[3];
    Double_t Mean[binnumber];
    Double_t MeanStat[binnumber];
    Double_t MeanSysHi[binnumber];
    Double_t MeanSysLo[binnumber];
    Double_t Chi[binnumber];
    //
    //Systematics Files
    TFile* fileSyst[3];
    TFile* fileSgnLossSyst = TFile::Open("NormalizationCorrections/SgnLossSyst.root");
    TH1D* hsgnloss[binnumber];
    TFile* fileYieldsSyst = TFile::Open(Form("sistematiche/CorrSystematics-Yield-%s-%s-sel-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE));	
    //MC Varied sys
    TFile* fileMCVaried = TFile::Open("MCVariedSyst.root");
    TH1D* hMCVariedUncorr = (TH1D*)fileMCVaried->Get(nameMCVaried.Data());
    TH1D* hMCVariedTot = (TH1D*)fileMCVaried->Get(Form("%sTot",nameMCVaried.Data()));

    //Take systematics
    TH1F* hSystNormYields = (TH1F*)fileYieldsSyst->Get("hSystTot");
    for (int i = 0; i < 3; i++){
      double zdcmin = lLoEE, zdcmax = lHiEE, v0min = lLoMult, v0max = lHiMult;
      if (fWhichEstimator.Contains("V0M")) {
        v0min = systpercentile[i];
        v0max = systpercentile[i+1];
      }
      if (fWhichEstimator.Contains("ZDC")) {
        zdcmin = systpercentile[i];
        zdcmax = systpercentile[i+1];
      }
      //files
	    fileSyst[i] = TFile::Open(Form("sistematiche/SystematicsFinalResults-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",0.,100.,0.,100.));//v0min,v0max,zdcmin,zdcmax));	
      // histo x systematics
      hSystSpectra[i] = (TH1D*)fileSyst[i]->Get("hSystTot");
    }
    //
    //
    // Fit function
    TF1* LevyTsallisfunc = LevyTsallis("LevyTsallisfunc", 1.321);
    //
    //Get Spectra
    Double_t lLo, lHi;
    if (fWhichEstimator.Contains("V0M")) {
       lLo = lLoEE;
       lHi = lHiEE;
     }
     if (fWhichEstimator.Contains("ZDC")) {
       lLo = lLoMult;
       lHi = lHiMult;
     }
    TFile* lResultsFile = TFile::Open(Form("StatSpectra%s-Xi-%s_%03.0f_%03.0f.root",fWhichEstimator.Data(),fWhichOtherEstimator.Data(),lLo,lHi));
    for (int i = 0; i < binnumber; i++){    
      //syst from sgn loss
      hsgnloss[i] = (TH1D*)fileSgnLossSyst->Get(Form("%s/hRatioClone%i",namesystsgnloss.Data(),i));
      //
      double zdcmin = lLoEE, zdcmax = lHiEE, v0min = lLoMult, v0max = lHiMult;
      if (fWhichEstimator.Contains("V0M")) {
        v0min = percentile[i];
        v0max = percentile[i+1];
      }
      if (fWhichEstimator.Contains("ZDC")) {
        zdcmin = percentile[i];
        zdcmax = percentile[i+1];
      }
      lHistPt[i] = (TH1D*)lResultsFile->Get(Form("XiSpectra_Stats_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",v0min,v0max,zdcmin,zdcmax)); 
      //
      lSystPt[i] = (TH1D*)lHistPt[i]->Clone(Form("XiSpectra_Syst_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",v0min,v0max,zdcmin,zdcmax));
      for (int jbin = 1; jbin <= lSystPt[i]->GetNbinsX(); jbin++){
      lSystPt[i]->SetBinError(jbin,(lSystPt[i]->GetBinContent(jbin)*
                                                    TMath::Sqrt( hSystSpectra[0]->GetBinContent(jbin)*hSystSpectra[0]->GetBinContent(jbin) 
                                                    + TMath::Abs(hsgnloss[i]->GetBinContent(jbin+1)-1)*TMath::Abs(hsgnloss[i]->GetBinContent(jbin+1)-1)
                                                    )
                ));
      }
    }

    // Compute Yields
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
        
        hout[i] = (TH1D*) YieldMean(lHistPt[i],lSystPt[i], LevyTsallisfunc);
        // Get <pT>
        Mean[i] = hout[i]->GetBinContent(kMean);
        MeanStat[i]= hout[i]->GetBinContent(kMeanStat);
        MeanSysHi[i]= hout[i]->GetBinContent(kMeanSysHi);
        MeanSysLo[i]= hout[i]->GetBinContent(kMeanSysLo);
    }//end loop over classes 
     
    // Write numbers
    std::ofstream outfile;
    outfile.open(Form("MeanpT_%s_%.0f-%.0f_%s_%.0f-%.0f.txt","V0M",lLoMult,lHiMult,"ZDC",lLoEE,lHiEE), std::ofstream::out | std::ofstream::trunc); // std::ios_base::app
    outfile << Form("\n %s [%.0f-%.0f] -  %s [%.0f-%.0f]","V0M",lLoMult,lHiMult,"ZDC",lLoEE,lHiEE) << endl;
    outfile << "\nMean pT \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << Mean[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Mean pT Stat \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << MeanStat[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Mean pT Syst High \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << MeanSysHi[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Mean pT Syst Low \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << MeanSysLo[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << " \n\n\n" << endl;     

    //
    TGraphErrors *MeanNchStat = new TGraphErrors(binnumber,dNch,Mean,StatdNch,MeanStat);
    TGraphErrors *MeanPercStat = new TGraphErrors(binnumber,centrpercentile,Mean,StatdNch,MeanStat);
    //
   
    TFile* Write = new TFile(Form("MeanpT-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE), "recreate");
    MeanNchStat->SetName("MeanpTvsNch_Stat");
    MeanPercStat->SetName("MeanpTvsPerc_Stat");
    MeanNchStat->Write();
    MeanPercStat->Write();

    //Write->Close();

   

  
}

//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
TH1 *
YieldMean(TH1 *hstat, TH1 *hsys, TF1 *f = NULL, Double_t minfit=.6,Double_t maxfit=6.5, Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q",TString logfilename="log.root")
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
  TH1 *htot = (TH1 *)hstat->Clone(Form("%sfittedwith%s",hstat->GetName(),f->GetName()));
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
  //f->SetName("");
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
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e5);
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
