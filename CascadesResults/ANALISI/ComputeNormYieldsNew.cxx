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
#include "TDirectory.h"
//#include "src/Utils.h"

#include <memory>

#endif

using namespace utils;
using namespace std;
/* definition of the fields in the histogram returned */

namespace yieldmeannew{

  enum EValue_t {
    kYield = 1,
    kYieldStat,
    kYieldSysTot,
    kYieldSysHiCorr,
    kYieldSysLoCorr,
    kMean,
    kMeanStat,
    kMeanSysTot,
    kMeanSysHardCorr,
    kMeanSysSoftCorr,
    kFitRes
  };

  void YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, double &integral, double &mean,Bool_t printinfo=kFALSE);
  TH1 * YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, double min, double binwidth = 0.01);
  TH1 * YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, double max, double binwidth = 0.1);
  TH1 * YieldMean_ReturnRandom(TH1 *hin);
  TH1 * YieldMean_ReturnCoherentRandom(TH1 *hin);
  TH1 * YieldMean_ReturnExtremeHisto(TH1 *hin, const char* hout_title, float sign = 1.);
  TH1 * YieldMean_ReturnExtremeHardHisto(TH1 *hin, const char* hout_title);
  TH1 * YieldMean_ReturnExtremeSoftHisto(TH1 *hin, const char* hout_title);
  TH1 * YieldMean_ReturnExtremeLowHisto(TH1 *hin, const char* hout_title);
  TH1 * YieldMean_ReturnExtremeHighHisto(TH1 *hin, const char* hout_title);
  int Fitter(TH1* histo, TF1* func, Option_t *opt);
  void RandomShifter(TH1* hin, TH1* hhi, TH1* hlo, TH1*& hIntegral, TH1*& hMean, float integral_limit, float mean_limit, int nRepetitions = 1000);
  void SaveToFile(const char* file_name, const char* dir_name, TCanvas& canvas);
  //TH1 * MergeHistograms(TH1 *hdata, TH1 *hlo, TH1 *hhi, const char* hout_name);

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) ;

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);

void ComputeNormYields(
  TString fWhichEstimator = "V0M", 
  Double_t lLoMult = 0., 
  Double_t lHiMult = 100., 
  Double_t lLoEE = 0., 
  Double_t lHiEE = 100.) {
 
    TString fWhichOtherEstimator = "V0M";
    if (fWhichEstimator.Contains("V0M")) fWhichOtherEstimator = "ZDC";
    //
    Double_t percentileV0[] = {0.,1.,5,10,15,20,30,40,50,70,100};
    Double_t systpercentileV0[4] = {0.,15.,40.,100.};
    Int_t systcounterV0[3] = {5,8,10};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    //
    Double_t percentileZDC[] = {0,20,30,40,50,60,70,80,90,100};
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
     
    Double_t dNchZDC_V070100[nbinZDC] = { 2.646, 2.829 ,2.827 ,2.786 ,2.732 ,2.693 ,2.659 ,2.582 ,2.364 }; 
    Double_t PosSystdNchZDC_V070100[nbinZDC] = { 0.034 , 0.037 ,0.037 ,0.036 ,0.036 ,0.035 ,0.035 ,0.034 ,0.031  };
    Double_t NegSystdNchZDC_V070100[nbinZDC] = { 0.034 , 0.037 ,0.037 ,0.036 ,0.036 ,0.035 ,0.035 ,0.034 ,0.031}; 
   
    Double_t dNchV0_ZDC030[nbinV0] = {26.35 , 20.51 , 16.86 , 14.54 , 12.83 , 10.91 , 8.835 , 7.165 , 5.234 ,2.701 };
    Double_t PosSystdNchV0_ZDC030[nbinV0] = { 0.32 , 0.25 ,  0.20 , 0.17 , 0.15 , 0.13 , 0.106 , 0.087 , 0.065 , 0.035 };
    Double_t NegSystdNchV0_ZDC030[nbinV0] = {0.32 , 0.25 , 0.20 , 0.17 , 0.15 , 0.13 , 0.106 , 0.087 , 0.065 , 0.035 };
    Double_t dNchV0_ZDC70100[nbinV0] = {21.30 , 17.79 , 14.39 , 12.60 , 11.02 , 9.324 , 7.568 , 6.202 , 4.580 , 2.494  };
    Double_t PosSystdNchV0_ZDC70100[nbinV0] = {0.37 , 0.23 ,0.17 ,0.15, 0.13, 0.113, 0.092, 0.078, 0.059 ,0.032};
    Double_t NegSystdNchV0_ZDC70100[nbinV0] = {0.37 , 0.23 ,0.17 ,0.15, 0.13, 0.113, 0.092, 0.078, 0.059 ,0.032};
    
       //
    Double_t OffdNchMB = 6.89;
    Double_t OffSystdNchMB = 0.11;
    TString namesystsgnloss = ""; 

    // Choose scenario and initialize correct variables
    int tempbinnumber = 0;
    if (fWhichEstimator.Contains("V0M")) {
        tempbinnumber = nbinV0;
        dNch = dNchV0;
        PosSystdNch = PosSystdNchV0;
        NegSystdNch = NegSystdNchV0;
        StatdNch = StatdNchV0;
        if (lLoEE > 0) {
          dNch = dNchV0_ZDC70100;
          PosSystdNch = PosSystdNchV0_ZDC70100;
          NegSystdNch = NegSystdNchV0_ZDC70100;
          namesystsgnloss = Form("%s-%s","multsel_fixedlowEE","V0FixLowEE");
          
          }
        if (lHiEE < 100) {
          dNch = dNchV0_ZDC030;
          PosSystdNch = PosSystdNchV0_ZDC030;
          NegSystdNch = NegSystdNchV0_ZDC030;
          namesystsgnloss = Form("%s-%s","multsel_fixedhighEE","V0FixHighEE");
          
        }
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
        if (lLoMult > 0) {
          dNch = dNchZDC_V070100;
          PosSystdNch = PosSystdNchZDC_V070100;
          NegSystdNch = NegSystdNchZDC_V070100;
          namesystsgnloss = Form("%s-%s","EEsel_fixedlowmult","ZDCFixLowmult");
        }
        if (lHiMult < 100) {
          dNch = dNchZDC_V0030;
          PosSystdNch = PosSystdNchZDC_V0030;
          NegSystdNch = NegSystdNchZDC_V0030;
          namesystsgnloss = Form("%s-%s","EEsel_fixedhighmult","ZDCFixHighmult");
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
    Double_t Yield[binnumber];
    Double_t YieldStat[binnumber];
    Double_t YieldSysHi[binnumber];
    Double_t YieldSysLo[binnumber];
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
    TFile* fileExtrapSyst = TFile::Open(Form("sistematiche/Final-Extrap-Syst-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,30.));//lHiEE));	
    
    //Take systematics
    TH1F* hSystNormYields = (TH1F*)fileYieldsSyst->Get("hSystTot");
    TH1F* hSystExtrapYields = (TH1F*)fileExtrapSyst->Get("hSystUncorr");
    TH1F* hSystExtrapYieldsCorr = (TH1F*)fileExtrapSyst->Get("hSystCorr");
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
                                                   // + TMath::Abs(hsgnloss[i]->GetBinContent(jbin+1)-1)*TMath::Abs(hsgnloss[i]->GetBinContent(jbin+1)-1)
                                                    )
                ));
      }
    }

    // Compute Yields
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
        
        hout[i] = (TH1D*) YieldMeanNew(lHistPt[i],lSystPt[i], 0x0,0x0,0x0,0x0,LevyTsallisfunc);
        // Get yield
        Yield[i] = hout[i]->GetBinContent(kYield);
        YieldStat[i] = hout[i]->GetBinContent(kYieldStat);
        YieldSysHi[i] = hout[i]->GetBinContent(kYieldSysTot);
        YieldSysLo[i] = 0;//hout[i]->GetBinContent(kYieldSysLo);
        // Get <pT>
       
    }//end loop over classes 

    //
    double YieldMB = 0.0273555395;//houtMB->GetBinContent(1);
    double YieldStatMB = 0.0001852769;//houtMB->GetBinContent(2);
    double YieldSysHiMB = 0.0019195267;//houtMB->GetBinContent(3);
    double YieldSysLoMB = 0.0019195267;//houtMB->GetBinContent(4);
    //
    //Compute Normalized Yields
    double NormYield[binnumber];
    double NormYieldStat[binnumber];
    double PosNormYieldSys[binnumber];
    double NegNormYieldSys[binnumber];
    double NormdNch[binnumber];
    double AvYield[binnumber];
    double AvYieldStat[binnumber];
    double PosAvYieldSys[binnumber];
    double NegAvYieldSys[binnumber];
    double PosAvYieldSysCorr[binnumber];
    double NegAvYieldSysCorr[binnumber];
    double AvdNch[binnumber];
    double TotUncorr[binnumber];
    //  
    for(int i = 0; i < binnumber ; i++){
      NormYield[i] = (Yield[i]/dNch[i])/(YieldMB/dNchMB);
      NormYieldStat[i] = ErrorInRatio(Yield[i],YieldStat[i],YieldMB, YieldStatMB)*dNchMB/dNch[i];
      AvYield[i] = (Yield[i]/dNch[i]);
      AvYieldStat[i] = YieldStat[i]/dNch[i];


      if (i < systcounter[0])   {
        PosNormYieldSys[i] = NormYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(1)*hSystNormYields->GetBinContent(1)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+PosSystdNch[i]/dNch[i]*PosSystdNch[i]/dNch[i]+OffSystdNchMB/OffdNchMB*OffSystdNchMB/OffdNchMB);
        NegNormYieldSys[i] = NormYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(1)*hSystNormYields->GetBinContent(1)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+NegSystdNch[i]/dNch[i]*NegSystdNch[i]/dNch[i]+OffSystdNchMB/OffdNchMB*OffSystdNchMB/OffdNchMB);
        PosAvYieldSys[i] = AvYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(1)*hSystNormYields->GetBinContent(1)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+PosSystdNch[i]/dNch[i]*PosSystdNch[i]/dNch[i]);
        NegAvYieldSys[i] = AvYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(1)*hSystNormYields->GetBinContent(1)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+NegSystdNch[i]/dNch[i]*NegSystdNch[i]/dNch[i]);    
        PosAvYieldSysCorr[i] = AvYield[i]*TMath::Sqrt(YieldSysHi[i]/Yield[i]*YieldSysHi[i]/Yield[i]+hSystExtrapYieldsCorr->GetBinContent(i+1)*hSystExtrapYieldsCorr->GetBinContent(i+1)+PosSystdNch[i]/dNch[i]*PosSystdNch[i]/dNch[i]);
        NegAvYieldSysCorr[i] = AvYield[i]*TMath::Sqrt(YieldSysLo[i]/Yield[i]*YieldSysLo[i]/Yield[i]+hSystExtrapYieldsCorr->GetBinContent(i+1)*hSystExtrapYieldsCorr->GetBinContent(i+1)+NegSystdNch[i]/dNch[i]*NegSystdNch[i]/dNch[i]);    
       // cout << (hSystNormYields->GetBinContent(1)*hSystNormYields->GetBinContent(1)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+PosSystdNch[i]/dNch[i]*PosSystdNch[i]/dNch[i]) - (YieldSysHi[i]/Yield[i]*YieldSysHi[i]/Yield[i]+hSystExtrapYieldsCorr->GetBinContent(i+1)*hSystExtrapYieldsCorr->GetBinContent(i+1)+PosSystdNch[i]/dNch[i]*PosSystdNch[i]/dNch[i]) <<endl;
        cout << " uncorr Y " << hSystNormYields->GetBinContent(1) << "     " << "tot Y " << YieldSysHi[i]/Yield[i] << endl;
          cout << " uncorr E " << hSystExtrapYields->GetBinContent(i+1) << "     " << "tot E " << hSystExtrapYieldsCorr->GetBinContent(i+1)  << endl;
  
      }
      if (i >= systcounter[0] && i < systcounter[1])  {
        PosNormYieldSys[i] = NormYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(2)*hSystNormYields->GetBinContent(2)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+PosSystdNch[i]/dNch[i]*PosSystdNch[i]/dNch[i]+OffSystdNchMB/OffdNchMB*OffSystdNchMB/OffdNchMB);
        NegNormYieldSys[i] = NormYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(2)*hSystNormYields->GetBinContent(2)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+NegSystdNch[i]/dNch[i]*NegSystdNch[i]/dNch[i]+OffSystdNchMB/OffdNchMB*OffSystdNchMB/OffdNchMB);
        PosAvYieldSys[i] = AvYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(2)*hSystNormYields->GetBinContent(2)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+PosSystdNch[i]/dNch[i]*PosSystdNch[i]/dNch[i]);
        NegAvYieldSys[i] = AvYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(2)*hSystNormYields->GetBinContent(2)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+NegSystdNch[i]/dNch[i]*NegSystdNch[i]/dNch[i]); 
        PosAvYieldSysCorr[i] = AvYield[i]*TMath::Sqrt(YieldSysHi[i]/Yield[i]*YieldSysHi[i]/Yield[i]+hSystExtrapYieldsCorr->GetBinContent(i+1)*hSystExtrapYieldsCorr->GetBinContent(i+1)+PosSystdNch[i]/dNch[i]*PosSystdNch[i]/dNch[i]);
        NegAvYieldSysCorr[i] = AvYield[i]*TMath::Sqrt(YieldSysLo[i]/Yield[i]*YieldSysLo[i]/Yield[i]+hSystExtrapYieldsCorr->GetBinContent(i+1)*hSystExtrapYieldsCorr->GetBinContent(i+1)+NegSystdNch[i]/dNch[i]*NegSystdNch[i]/dNch[i]);    
     cout << " uncorr Y " << hSystNormYields->GetBinContent(2) << "     " << "tot Y " << YieldSysHi[i]/Yield[i] << endl;
          cout << " uncorr E " << hSystExtrapYields->GetBinContent(i+1) << "     " << "tot E " << hSystExtrapYieldsCorr->GetBinContent(i+1)  << endl;
  
      }      
      if (i >= systcounter[1] && i < systcounter[2]) {
        PosNormYieldSys[i] = NormYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(3)*hSystNormYields->GetBinContent(3)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+PosSystdNch[i]/dNch[i]*PosSystdNch[i]/dNch[i]+OffSystdNchMB/OffdNchMB*OffSystdNchMB/OffdNchMB);
        NegNormYieldSys[i] = NormYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(3)*hSystNormYields->GetBinContent(3)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+NegSystdNch[i]/dNch[i]*NegSystdNch[i]/dNch[i]+OffSystdNchMB/OffdNchMB*OffSystdNchMB/OffdNchMB);
        PosAvYieldSys[i] = AvYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(3)*hSystNormYields->GetBinContent(3)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+PosSystdNch[i]/dNch[i]*PosSystdNch[i]/dNch[i]);
        NegAvYieldSys[i] = AvYield[i]*TMath::Sqrt(hSystNormYields->GetBinContent(3)*hSystNormYields->GetBinContent(3)+hSystExtrapYields->GetBinContent(i+1)*hSystExtrapYields->GetBinContent(i+1)+NegSystdNch[i]/dNch[i]*NegSystdNch[i]/dNch[i]);           
        PosAvYieldSysCorr[i] = AvYield[i]*TMath::Sqrt(YieldSysHi[i]/Yield[i]*YieldSysHi[i]/Yield[i]+hSystExtrapYieldsCorr->GetBinContent(i+1)*hSystExtrapYieldsCorr->GetBinContent(i+1)+PosSystdNch[i]/dNch[i]*PosSystdNch[i]/dNch[i]);
        NegAvYieldSysCorr[i] = AvYield[i]*TMath::Sqrt(YieldSysLo[i]/Yield[i]*YieldSysLo[i]/Yield[i]+hSystExtrapYieldsCorr->GetBinContent(i+1)*hSystExtrapYieldsCorr->GetBinContent(i+1)+NegSystdNch[i]/dNch[i]*NegSystdNch[i]/dNch[i]);    
           cout << " uncorr Y " << hSystNormYields->GetBinContent(3) << "     " << "tot Y " << YieldSysHi[i]/Yield[i] << endl;
          cout << " uncorr E " << hSystExtrapYields->GetBinContent(i+1) << "     " << "tot E " << hSystExtrapYieldsCorr->GetBinContent(i+1)  << endl;
  
      }   
    }
   
    // Write numbers
    std::ofstream outfile;
    outfile.open(Form("Yields_%s_%.0f-%.0f_%s_%.0f-%.0f.txt","V0M",lLoMult,lHiMult,"ZDC",lLoEE,lHiEE), std::ofstream::out | std::ofstream::trunc); // std::ios_base::app
    outfile << Form("\n %s [%.0f-%.0f] -  %s [%.0f-%.0f]","V0M",lLoMult,lHiMult,"ZDC",lLoEE,lHiEE) << endl;
    outfile << "\nYields \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << Yield[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Yields Stat \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << YieldStat[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Yields Syst High \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << YieldSysHi[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Yields Syst Low \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << YieldSysLo[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << " \n\n\n" << endl;
  //
    outfile << "\nYields Normalized\n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << NormYield[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Yields Stat \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << NormYieldStat[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Yields Syst High \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << PosNormYieldSys[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Yields Syst Low \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << NegNormYieldSys[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << " \n\n\n" << endl;


    // Do Plots
    TGraphErrors *YieldsNchStat = new TGraphErrors(binnumber,dNch,Yield,StatdNch,YieldStat);
    TGraphAsymmErrors *YieldsNchSyst = new TGraphAsymmErrors(binnumber,dNch,Yield,NegSystdNch,PosSystdNch,YieldSysLo,YieldSysHi);
    //
    TGraphErrors *MeanNchStat = new TGraphErrors(binnumber,dNch,Mean,StatdNch,MeanStat);
    TGraphAsymmErrors *MeanNchSyst = new TGraphAsymmErrors(binnumber,dNch,Mean,NegSystdNch,PosSystdNch,MeanSysLo,MeanSysHi);
    //
    TGraphErrors *NormYieldsNchStat = new TGraphErrors(binnumber,dNch,NormYield,StatdNch,NormYieldStat);
    TGraphAsymmErrors *NormYieldsNchSyst = new TGraphAsymmErrors(binnumber,dNch,NormYield,NegSystdNch,PosSystdNch,NegNormYieldSys,PosNormYieldSys);
    //
    TGraphErrors *NormYieldspercentile = new TGraphErrors(binnumber,centrpercentile,NormYield,StatdNch,NormYieldStat);
    TGraphAsymmErrors *NormYieldsSystpercentile = new TGraphAsymmErrors(binnumber,centrpercentile,NormYield,errpercentile,errpercentile,NegNormYieldSys,PosNormYieldSys);
    //
    TGraphErrors *AvYieldsNchStat = new TGraphErrors(binnumber,dNch,AvYield,StatdNch,AvYieldStat);
    TGraphAsymmErrors *AvYieldsNchSyst = new TGraphAsymmErrors(binnumber,dNch,AvYield,NegSystdNch,PosSystdNch,NegAvYieldSys,PosAvYieldSys);
    TGraphAsymmErrors *AvYieldsNchSystCorr = new TGraphAsymmErrors(binnumber,dNch,AvYield,NegSystdNch,PosSystdNch,NegAvYieldSysCorr,PosAvYieldSysCorr);
   
    //
    TGraphErrors *AvYieldspercentile = new TGraphErrors(binnumber,centrpercentile,AvYield,StatdNch,AvYieldStat);
    TGraphAsymmErrors *AvYieldsSystpercentile = new TGraphAsymmErrors(binnumber,centrpercentile,AvYield,errpercentile,errpercentile,NegAvYieldSys,PosAvYieldSys);
    TGraphAsymmErrors *AvYieldsSystCorrpercentile = new TGraphAsymmErrors(binnumber,centrpercentile,AvYield,errpercentile,errpercentile,NegAvYieldSysCorr,PosAvYieldSysCorr);
    AvYieldsNchStat->SetName("AvYieldsNchStat");
    AvYieldsNchSyst->SetName("AvYieldsNchSyst");
    AvYieldsNchSystCorr->SetName("AvYieldsNchSystCorr");
    AvYieldspercentile->SetName("AvYieldsStatpercentile");
    AvYieldsSystpercentile->SetName("AvYieldsSystpercentile");
    AvYieldsSystCorrpercentile->SetName("AvYieldsSystCorrpercentile");
    //
    TGraphErrors *ChiGraph = new TGraphErrors(binnumber,centrpercentile,Chi,errpercentile,StatdNch); 

    // Draw stuff


    // sample histo 
    TH1D* g = new TH1D("g", " ", 40, 0., 40.);
    g->SetBinContent(1,0.);
    g->SetBinContent(29,2.);
    g->SetLineColor(kWhite);
    TCanvas* c = new TCanvas("c","",800,800);
    c->SetRightMargin(0.09);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetGridy();
    c->SetGridx();
    g->GetXaxis()->SetRangeUser(0.,30);
    g->GetYaxis()->SetRangeUser(0.,0.14);
    g->GetYaxis()->SetTitle("#LT dN/dy #GT");
    g->GetYaxis()->SetTitleSize(0.05);
    g->GetYaxis()->SetTitleOffset(1.1);
    g->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}");
    g->GetXaxis()->SetTitleSize(0.05);
    g->GetXaxis()->SetTitleOffset(0.9);
    g->SetTitle("");
    g->SetStats(0);
    g->Draw();

   TPavesText *pst = new TPavesText(2.258112,0.03395189,5.67125,0.04311871,5,"br");
   pst->SetBorderSize(1);
   pst->SetLineWidth(0);
   pst->SetTextFont(42);
   pst->SetTextSize(0.03104787);
   TText *pst_LaTex = pst->AddText("stat.");
   pst_LaTex = pst->AddText("syst.");
   pst->Draw("SAME");
   TMarker *marker = new TMarker(2.777503,0.04075613,2);
   marker->SetMarkerStyle(2);
   marker->SetMarkerSize(2.9);
   marker->Draw();
   marker = new TMarker(2.777503,0.03650348,25);
   marker->SetMarkerStyle(25);
   marker->SetMarkerSize(2.8);
   marker->Draw("SAME");
   TLatex* tex = new TLatex(2.480708,0.04378024,"V0M selected");
   tex->SetTextSize(0.03622251);
   tex->SetLineWidth(2);
   tex->Draw("SAME");

    YieldsNchSyst->SetName("YieldsNchSyst");
    YieldsNchStat->SetName("YieldsNchStat");
    YieldsNchSyst->SetFillStyle(3000);
    YieldsNchSyst->SetFillColor(kRed);
    YieldsNchSyst->SetLineColor(kRed);
    YieldsNchSyst->Draw("SAME PE2");
    YieldsNchStat->SetMarkerStyle(20);
    YieldsNchStat->SetMarkerColor(kRed);
    YieldsNchStat->SetLineColor(kRed);
    YieldsNchStat->SetMarkerSize(1.8); 
    YieldsNchStat->Draw("SAME ep");
    //
   /* TH1D* g1 = new TH1D("g1", " ", 40, 0., 40.);
    g1->SetBinContent(1,0.);
    g1->SetBinContent(29,2.);
    g1->SetLineColor(kWhite);
    g1->SetStats(0);
    TCanvas* d = new TCanvas("d","",800,800);
    d->SetRightMargin(0.09);
    d->SetLeftMargin(0.15);
    d->SetBottomMargin(0.15);
    g1->GetXaxis()->SetRangeUser(0.,28.);
    g1->GetYaxis()->SetRangeUser(0.8,1.7);
    g1->GetYaxis()->SetTitle("<p_{T}>");
    g1->GetYaxis()->SetTitleSize(0.05);
    g1->GetYaxis()->SetTitleOffset(0.9);
    g1->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetXaxis()->SetTitleOffset(0.9);
    g1->SetTitle("");
    g1->Draw();
    MeanNchSyst->SetFillStyle(3004);
    MeanNchSyst->SetFillColor(kBlue);
    MeanNchSyst->SetLineColor(kBlue);
    MeanNchSyst->Draw("SAME PE2");
    MeanNchStat->SetMarkerColor(kBlue);
    MeanNchStat->SetMarkerStyle(20);
    MeanNchStat->SetMarkerSize(1.2); 
    MeanNchStat->Draw("SAME EP");*/
    //
    TH1D* g2 = new TH1D("g2", " ", 40, 0., 40.);
    g2->SetBinContent(1,0.);
    g2->SetBinContent(29,2.);
    g2->SetLineColor(kWhite);
    g2->SetStats(0);
    TCanvas* n = new TCanvas("n","",1000,900);
    n->SetRightMargin(0.09);
    n->SetLeftMargin(0.15);
    n->SetBottomMargin(0.15);
    n->SetGridy();
    g2->SetStats(0);
    g2->GetXaxis()->SetRangeUser(0.,28);
    g2->GetYaxis()->SetRangeUser(0.,2.3);
    g2->GetYaxis()->SetTitle("1/n_{ch} <dN/dy> / (1/n_{ch} <dN/dy>)_{MB}");
    g2->GetYaxis()->SetTitleSize(0.04);
    g2->GetYaxis()->SetTitleOffset(1.2);
    g2->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    g2->GetXaxis()->SetTitleSize(0.04);
    g2->GetXaxis()->SetTitleOffset(1.2);
    g2->Draw();
    NormYieldsNchSyst->SetFillColor(kBlue);
    NormYieldsNchSyst->SetFillStyle(3000);
    NormYieldsNchSyst->Draw("SAME E2");
    NormYieldsNchStat->SetLineColor(kBlue);
    NormYieldsNchStat->SetMarkerColor(kBlue);
    NormYieldsNchStat->SetTitle("");
    NormYieldsNchStat->SetMarkerStyle(22);
    NormYieldsNchStat->SetMarkerSize(1.2); 
    NormYieldsNchStat->Draw("SAME ep");
    NormYieldsNchStat->SetName("NormYieldsNchStat");
    NormYieldsNchSyst->SetName("NormYieldsNchSyst");
    //
    TH1D* g3 = new TH1D("g3", " ", 102, -1., 101.);
    g3->SetBinContent(1,0.);
    g3->SetBinContent(29,2.);
    g3->SetLineColor(kWhite);
    g3->SetStats(0);
    TCanvas* p = new TCanvas("p","",1000,900);
    p->SetRightMargin(0.09);
    p->SetLeftMargin(0.15);
    p->SetBottomMargin(0.15);
    p->SetGridy();
    g3->SetStats(0);
    g3->GetXaxis()->SetRangeUser(-1.,101);
    g3->GetYaxis()->SetRangeUser(0.,2.3);
    g3->GetYaxis()->SetTitle("1/n_{ch} <dN/dy> / (1/n_{ch} <dN/dy>)_{MB}");
    g3->GetYaxis()->SetTitleSize(0.04);
    g3->GetYaxis()->SetTitleOffset(1.2);
    g3->GetXaxis()->SetTitle(Form("percentile %s (%)",fWhichEstimator.Data()));
    g3->GetXaxis()->SetTitleSize(0.04);
    g3->GetXaxis()->SetTitleOffset(1.2);
    g3->Draw();

    TLegend* l = new TLegend (0.6,0.75,0.89,0.89);

    NormYieldspercentile->GetXaxis()->SetRangeUser(-1.,101);
    NormYieldspercentile->GetYaxis()->SetRangeUser(0.,1.5);
    NormYieldspercentile->GetYaxis()->SetTitle("1/n_{ch} <dN/dy> / (1/n_{ch} <dN/dy>)_{MB}");
    NormYieldspercentile->GetYaxis()->SetTitleSize(0.04);
    NormYieldspercentile->GetYaxis()->SetTitleOffset(1.2);
    NormYieldspercentile->GetXaxis()->SetTitle(Form("percentile %s (%)",fWhichEstimator.Data()));
    NormYieldspercentile->GetXaxis()->SetTitleSize(0.04);
    NormYieldspercentile->GetXaxis()->SetTitleOffset(1.2);
    NormYieldspercentile->SetLineColor(kBlue);
    NormYieldspercentile->SetMarkerColor(kBlue);
    NormYieldspercentile->SetTitle("");
    NormYieldspercentile->SetMarkerStyle(20);
    NormYieldspercentile->SetMarkerSize(1.2); 
    //
    NormYieldsSystpercentile->GetXaxis()->SetRangeUser(-1.,101);
    NormYieldsSystpercentile->GetYaxis()->SetRangeUser(0.,1.5);
    NormYieldsSystpercentile->GetYaxis()->SetTitle("1/n_{ch} <dN/dy> / (1/n_{ch} <dN/dy>)_{MB}");
    NormYieldsSystpercentile->GetYaxis()->SetTitleSize(0.04);
    NormYieldsSystpercentile->GetYaxis()->SetTitleOffset(1.2);
    NormYieldsSystpercentile->GetXaxis()->SetTitle(Form("percentile %s (%)",fWhichEstimator.Data()));
    NormYieldsSystpercentile->GetXaxis()->SetTitleSize(0.04);
    NormYieldsSystpercentile->GetXaxis()->SetTitleOffset(1.2);
    NormYieldsSystpercentile->SetLineColor(kBlue);
    NormYieldsSystpercentile->SetFillColor(kBlue-10);
    NormYieldsSystpercentile->SetFillStyle(1001);
    NormYieldsSystpercentile->Draw("SAME E2");
    NormYieldspercentile->Draw("SAME ep");
    NormYieldspercentile->SetName("NormYieldsvspercentile_Stat");
    NormYieldsSystpercentile->SetName("NormYieldsvspercentile_Syst");

    l->SetBorderSize(0);
    l->AddEntry(NormYieldspercentile,"stat","EP");
    l->AddEntry(NormYieldsSystpercentile,"syst","F");
    l->SetTextSize(0.028);
    
    l->Draw("SAME");


    //Chi graph
    TCanvas* chi = new TCanvas("chi","",1700,900);
    chi->SetRightMargin(0.09);
    chi->SetLeftMargin(0.15);
    chi->SetBottomMargin(0.15);
    chi->SetGridy();
    ChiGraph->SetMarkerStyle(20);
    ChiGraph->SetMarkerSize(1.6);
    ChiGraph->SetMarkerColor(kRed);
    ChiGraph->SetLineColor(kRed);
    ChiGraph->SetLineWidth(2);
    ChiGraph->SetName("ChiNDFGraph");
    TH1D* bkg = new TH1D("bkg", " ", 102, -10., 110.);
    bkg->SetBinContent(1,0.);
    bkg->SetBinContent(50,2.);
    bkg->SetLineColor(kWhite);
    bkg->SetStats(0);
    bkg->GetYaxis()->SetTitle("#Chi^{2} / ndf");
    bkg->GetXaxis()->SetTitle(Form("percentile %s (%)",fWhichEstimator.Data()));
    bkg->Draw();
    ChiGraph->Draw("SAME EP");


    TFile* Write = new TFile(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE), "recreate");
    NormYieldsNchSyst->Write();
    NormYieldsNchStat->Write();
    NormYieldspercentile->Write();
    NormYieldsSystpercentile->Write();
    ChiGraph->Write();
    YieldsNchStat->Write();
    YieldsNchSyst->Write();

    AvYieldsNchSyst->Write();
    AvYieldsNchSystCorr->Write();
    AvYieldsNchStat->Write();
    AvYieldspercentile->Write();
    AvYieldsSystpercentile->Write();
    AvYieldsSystCorrpercentile->Write();

    //Write->Close();

   

  
}


TH1 *
  YieldMeanNew(TH1 *hstat, TH1 *hsys, TH1 *hsys_pt_uncorr, TH1 *hsys_pt_corr, TH1 *hsys_mult_corr, TF1* &fout, TF1 *f = nullptr, double min = 0., double max = 10., double loprecision = 0.01, double hiprecision = 0.1, bool store_log = false, const char* logfilename="log.root", const char* path = "", Option_t *opt = "0q")
  {
    /* set many iterations when fitting the data so we don't
      stop minimization with MAX_CALLS */
    TVirtualFitter::SetMaxIterations(1000000);

    /* create output histo */
    double integral, mean, integral_ref, mean_ref;
    TH1 *hout = new TH1D(Form("hout_%s",hstat->GetName()), "", kFitRes, 0, kFitRes);
    TH1 *hlo = nullptr, *hhi = nullptr;

    /* create histo with stat+sys errors */
    TH1 *htot = (TH1 *)hstat->Clone(Form("%sfittedwith%s",hstat->GetName(),f->GetName()));
    for (int ibin = 0; ibin < htot->GetNbinsX(); ibin++) {
      htot->SetBinError(ibin + 1, TMath::Sqrt(hsys->GetBinError(ibin + 1) * hsys->GetBinError(ibin + 1) + hstat->GetBinError(ibin + 1) * hstat->GetBinError(ibin + 1)));
    }

    /*
    *   measure the central value
    */
    int fitres = Fitter(htot,f,opt);
    hout->SetBinContent(kFitRes,fitres);

    fout = new TF1(*f);
    fout->SetLineColor(kRed);
    std::cout<<"Fit sys+stat for " <<f->GetName()<<std::endl;
    std::cout<<"NDF="<<f->GetNDF()<<" Chi^2="<<f->GetChisquare()<<" Chi^2/NDF="<<f->GetChisquare()/f->GetNDF()<<std::endl;

    hlo = YieldMean_LowExtrapolationHisto(htot, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(htot, f, max, hiprecision);

    YieldMean_IntegralMean(htot, hlo, hhi, integral, mean, true);
    hout->SetBinContent(kYield, integral);
    hout->SetBinContent(kMean, mean);
    integral_ref = integral;
    mean_ref = mean;

    TCanvas cCanvasDefault("cCanvasDefault");
    cCanvasDefault.DrawFrame(min,0,max,1.3 * hstat->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    htot->Draw("pesame");
    //hlo->Draw("pesame");
    //hhi->Draw("pesame");
    fout->Draw("same");


    if(store_log){
      SaveToFile(logfilename,path,cCanvasDefault);
    }

    /*
    * STATISTICS
    */

    TCanvas cCanvasStat("cCanvasStat");
    cCanvasStat.Divide(2, 1);

    /*
    * measure statistical error
    */

    /* fit with stat error */
    Fitter(hstat,f,opt);
    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hstat, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hstat, f, max, hiprecision);

    TH1* hIntegral = nullptr;
    TH1* hMean = nullptr;

    RandomShifter(hstat, hhi, hlo, hIntegral, hMean, integral_ref, mean_ref);

    TF1 *gaus = (TF1 *)gROOT->GetFunction("gaus");
    
    cCanvasStat.cd(1);
    hIntegral->Fit(gaus, "q");
    integral = hout->GetBinContent(kYield) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kYieldStat, integral);

    cCanvasStat.cd(2);
    hMean->Fit(gaus, "q");
    mean = hout->GetBinContent(kMean) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kMeanStat, mean);
    
    if(store_log){
      SaveToFile(logfilename,path,cCanvasStat);
    }

    TCanvas cCanvasStatExtra("cCanvasStatExtra");
    cCanvasStatExtra.DrawFrame(min,0,max,1.3 * hstat->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hstat->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasStatExtra);
    }
    delete hMean;
    delete hIntegral;

    /*
    * PT-UNCORRELATED SYSTEMATICS
    */

    /* fit with pt-uncorrelated syst error */
    TCanvas cCanvasSystPtUncorr("cCanvasSystPtUncorr");
    cCanvasSystPtUncorr.Divide(2, 1);

    Fitter(hsys_pt_uncorr,f,opt);
    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hsys_pt_uncorr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hsys_pt_uncorr, f, max, hiprecision);

    TH1* hIntegralSyst = nullptr;
    TH1* hMeanSyst = nullptr;

    RandomShifter(hsys_pt_uncorr, hhi, hlo, hIntegralSyst, hMeanSyst, integral_ref, mean_ref);
    
    cCanvasSystPtUncorr.cd(1);
    hIntegralSyst->Fit(gaus, "q");
    integral = hout->GetBinContent(kYield) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kYieldSysTot, integral);

    cCanvasSystPtUncorr.cd(2);
    hMeanSyst->Fit(gaus, "q");
    mean = hout->GetBinContent(kMean) * gaus->GetParameter(2) / gaus->GetParameter(1);
    hout->SetBinContent(kMeanSysTot, mean);
    
    if(store_log){
      SaveToFile(logfilename,path,cCanvasSystPtUncorr);
    }

    TCanvas cCanvasSystPtUncorrExtra("cCanvasSystPtUncorrExtra");
    cCanvasSystPtUncorrExtra.DrawFrame(min,0,max,1.3 * hsys_pt_uncorr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsys_pt_uncorr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasSystPtUncorrExtra);
    }

    delete hMeanSyst;
    delete hIntegralSyst;

    /*
    * PT-CORRELATED SYSTEMATICS
    */

    TCanvas cCanvasSysPtCorr("cCanvasYieldSysPtCorr");
    cCanvasSysPtCorr.Divide(2, 1);
    cCanvasSysPtCorr.cd(1)->DrawFrame(min, 0.7 * hsys_pt_corr->GetMinimum(), max, 1.3 * hsys_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsys_pt_corr->SetMarkerStyle(20);
    hsys_pt_corr->SetMarkerColor(1);
    hsys_pt_corr->SetMarkerSize(1);
    hsys_pt_corr->Draw("same");
    cCanvasSysPtCorr.cd(2)->DrawFrame(min, 0.7 * hsys_pt_corr->GetMinimum() , max, 1.3 * hsys_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsys_pt_corr->Draw("same");

    /*
    * systematic error high
    */

    TH1 *hhigh_pt_corr = YieldMean_ReturnExtremeHighHisto(hsys_pt_corr,"pt_corr_extreme_high");
    Fitter(hhigh_pt_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hhigh_pt_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhigh_pt_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hhigh_pt_corr, hlo, hhi, integral, mean);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    integral = sqrt(Sq(integral)+Sq(hout->GetBinContent(kYieldSysTot)));
    hout->SetBinContent(kYieldSysTot, integral);

    TCanvas cCanvasHighPtCorr("cCanvasHighPtCorr");
    cCanvasHighPtCorr.DrawFrame(min,0,max,1.3 * hhigh_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hhigh_pt_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasHighPtCorr);
    }

    cCanvasSysPtCorr.cd(1);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
    * systematic error hard
    */

    TH1 *hhard_pt_corr = YieldMean_ReturnExtremeHardHisto(hsys_pt_corr,"pt_corr_extreme_hard");
    Fitter(hhard_pt_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hhard_pt_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhard_pt_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hhard_pt_corr, hlo, hhi, integral, mean);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    mean = sqrt(Sq(mean) + Sq(hout->GetBinContent(kMeanSysTot)));
    hout->SetBinContent(kMeanSysTot, mean);

    TCanvas cCanvasHardPtCorr("cCanvasHardPtCorr");
    cCanvasHardPtCorr.DrawFrame(min,0,max,1.3 * hhard_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hhard_pt_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasHardPtCorr);
    }

    cCanvasSysPtCorr.cd(2);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
    * systematic error low
    */

    TH1 *hlow_pt_corr = YieldMean_ReturnExtremeLowHisto(hsys_pt_corr, "pt_corr_extreme_low");
    Fitter(hlow_pt_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hlow_pt_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hlow_pt_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hlow_pt_corr, hlo, hhi, integral, mean);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    integral = sqrt(Sq(integral) + Sq(hout->GetBinContent(kYieldSysTot)));
    hout->SetBinContent(kYieldSysTot, integral);

    TCanvas cCanvasLowPtCorr("cCanvasLowPtCorr");
    cCanvasLowPtCorr.DrawFrame(min,0,max,1.3 * hlow_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hlow_pt_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasLowPtCorr);
    }

    cCanvasSysPtCorr.cd(1);
    f->SetLineColor(4);
    f->DrawCopy("same");

    /*
    * systematic error soft
    */

    TH1 *hsoft_pt_corr = YieldMean_ReturnExtremeSoftHisto(hsys_pt_corr, "pt_corr_extreme_soft");
    Fitter(hsoft_pt_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hsoft_pt_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hsoft_pt_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hsoft_pt_corr, hlo, hhi, integral, mean);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    mean = sqrt(Sq(mean) + Sq(hout->GetBinContent(kMeanSysTot)));
    hout->SetBinContent(kMeanSysTot, mean);

    TCanvas cCanvasSoftPtCorr("cCanvasSoftPtCorr");
    cCanvasSoftPtCorr.DrawFrame(min,0,max,1.3 * hsoft_pt_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsoft_pt_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasSoftPtCorr);
    }

    cCanvasSysPtCorr.cd(2);
    f->SetLineColor(4);
    f->DrawCopy("same");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasSysPtCorr);
    }

    /*
    * MULTIPLICITY CORRELATED SYSTEMATICS
    */

    TCanvas cCanvasSysMultCorr("cCanvasYieldSysMultCorr");
    cCanvasSysMultCorr.Divide(2, 1);
    cCanvasSysMultCorr.cd(1)->DrawFrame(min, 0.7 * hsys_mult_corr->GetMinimum(), max, 1.3 * hsys_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsys_mult_corr->SetMarkerStyle(20);
    hsys_mult_corr->SetMarkerColor(1);
    hsys_mult_corr->SetMarkerSize(1);
    hsys_mult_corr->Draw("same");
    cCanvasSysMultCorr.cd(2)->DrawFrame(min, 0.7 * hsys_mult_corr->GetMinimum() , max, 1.3 * hsys_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsys_mult_corr->Draw("same");

    /*
    * systematic error high
    */

    TH1 *hhigh_mult_corr = YieldMean_ReturnExtremeHighHisto(hsys_mult_corr,"mult_corr_extreme_high");
    Fitter(hhigh_mult_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hhigh_mult_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhigh_mult_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hhigh_mult_corr, hlo, hhi, integral, mean);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    hout->SetBinContent(kYieldSysHiCorr, integral);

    TCanvas cCanvasHighMultCorr("cCanvasHighMultCorr");
    cCanvasHighMultCorr.DrawFrame(min,0,max,1.3 * hhigh_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hhigh_mult_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasHighMultCorr);
    }

    cCanvasSysMultCorr.cd(1);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
    * systematic error hard
    */

    TH1 *hhard_mult_corr = YieldMean_ReturnExtremeHardHisto(hsys_mult_corr,"mult_corr_extreme_hard");
    Fitter(hhard_mult_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hhard_mult_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hhard_mult_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hhard_mult_corr, hlo, hhi, integral, mean);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    hout->SetBinContent(kMeanSysHardCorr, mean);

    TCanvas cCanvasHardMultCorr("cCanvasHardMultCorr");
    cCanvasHardMultCorr.DrawFrame(min,0,max,1.3 * hhard_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hhard_mult_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasHardMultCorr);
    }

    cCanvasSysMultCorr.cd(2);
    f->SetLineColor(2);
    f->DrawCopy("same");

    /*
    * systematic error low
    */

    TH1 *hlow_mult_corr = YieldMean_ReturnExtremeLowHisto(hsys_mult_corr, "mult_corr_extreme_low");
    Fitter(hlow_mult_corr,f,opt);

    delete hlo;
    delete hhi;
    hlo = YieldMean_LowExtrapolationHisto(hlow_mult_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hlow_mult_corr, f, max, hiprecision);
    YieldMean_IntegralMean(hlow_mult_corr, hlo, hhi, integral, mean);
    integral = TMath::Abs(integral - hout->GetBinContent(kYield));
    hout->SetBinContent(kYieldSysLoCorr, integral);

    TCanvas cCanvasLowMultCorr("cCanvasLowMultCorr");
    cCanvasLowMultCorr.DrawFrame(min,0,max,1.3 * hlow_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hlow_mult_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasLowMultCorr);
    }

    cCanvasSysMultCorr.cd(1);
    f->SetLineColor(4);
    f->DrawCopy("same");

    /*
    * systematic error soft
    */

    TH1 *hsoft_mult_corr = YieldMean_ReturnExtremeSoftHisto(hsys_mult_corr, "mult_corr_extreme_soft");
    Fitter(hsoft_mult_corr,f,opt);

    delete hlo;
    delete hhi;

    hlo = YieldMean_LowExtrapolationHisto(hsoft_mult_corr, f, min, loprecision);
    hhi = YieldMean_HighExtrapolationHisto(hsoft_mult_corr, f, max, hiprecision);
 
    YieldMean_IntegralMean(hsoft_mult_corr, hlo, hhi, integral, mean);
    mean = TMath::Abs(mean - hout->GetBinContent(kMean));
    hout->SetBinContent(kMeanSysSoftCorr, mean);

    TCanvas cCanvasSoftMultCorr("cCanvasSoftMultCorr");
    cCanvasSoftMultCorr.DrawFrame(min,0,max,1.3 * hsoft_mult_corr->GetMaximum(),";#it{p}_{T} (GeV/#it{c});");
    hsoft_mult_corr->Draw("pesame");
    hlo->Draw("pesame");
    hhi->Draw("pesame");    

    if(store_log){
      SaveToFile(logfilename,path,cCanvasSoftMultCorr);
    }

    cCanvasSysMultCorr.cd(2);
    f->SetLineColor(4);
    f->DrawCopy("same");

    if(store_log){
      SaveToFile(logfilename,path,cCanvasSysMultCorr);
    }

    std::cout << "dN / dy = " << hout->GetBinContent(kYield) << " +- " << hout->GetBinContent(kYieldStat);
    std::cout << " +- " << hout->GetBinContent(kYieldSysTot);
    std::cout << " ( + " << hout->GetBinContent(kYieldSysHiCorr) << " - " << hout->GetBinContent(kYieldSysLoCorr) << ") ";
    std::cout << "\n<pT> = " << hout->GetBinContent(kMean) << " +- " << hout->GetBinContent(kMeanStat);
    std::cout << " + " << hout->GetBinContent(kMeanSysTot);
    std::cout << " ( + " << hout->GetBinContent(kMeanSysHardCorr) << " - " << hout->GetBinContent(kMeanSysSoftCorr) << ") ";
    std::cout << std::endl;

    delete hlo;
    delete hhi;
    return hout;
  }

  TH1 *
  YieldMean_LowExtrapolationHisto(TH1 *h, TF1 *f, double min, double binwidth)
  {
    /* find lowest edge in histo */
    int binlo = -1;
    double lo = 0.;
    for (int ibin = 1; ibin < h->GetNbinsX() + 1; ibin++) {
      if (h->GetBinContent(ibin) != 0.) {
        binlo = ibin;
        lo = h->GetBinLowEdge(ibin);
        break;
      }
    }
    if (binlo == -1) ::Fatal("YieldMean_LowExtrapolationHisto","Lower bin is undefined !?");

    int nbins = (lo - min) / binwidth;
    TH1 *hlo = new TH1F(Form("%s_lo",h->GetName()), "", nbins, min, lo);

    /* integrate function in histogram bins */
    double cont, err, width;
    for (int ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
      width = hlo->GetBinWidth(ibin + 1);
      cont = f->Integral(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), 1.e-6);//(double *)0, 1.e-6);
      err = f->IntegralError(hlo->GetBinLowEdge(ibin + 1), hlo->GetBinLowEdge(ibin + 2), (double *)0, (double *)0, 1.e-6);
      hlo->SetBinContent(ibin + 1, cont / width);
      hlo->SetBinError(ibin + 1, err / width);
    }

    return hlo;
  }

  TH1 *
  YieldMean_HighExtrapolationHisto(TH1 *h, TF1 *f, double max, double binwidth)
  {
    /* find highest edge in histo */
    int binhi = -1;
    double hi = 0.;
    for (int ibin = h->GetNbinsX(); ibin > 0; ibin--) {
      if (h->GetBinContent(ibin) != 0.) {
        binhi = ibin + 1;
        hi = h->GetBinLowEdge(ibin + 1);
        break;
      }
    }
    if (binhi == -1) ::Fatal("YieldMean_HighExtrapolationHisto","Higher bin is undefined !?");
    if(max<hi) {
      Printf("Warning! You should probably set a higher max value (Max = %f, hi = %f)", max, hi);
    }
    int nbins = (max - hi) / binwidth;
    TH1 *hhi = new TH1F(Form("%s_hi",h->GetName()), "", nbins, hi, max);

    /* integrate function in histogram bins */
    double cont, err, width;
    for (int ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
      width = hhi->GetBinWidth(ibin + 1);
      cont = f->Integral(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), 1.e-6); //(double *)0, 1.e-6);
      err = f->IntegralError(hhi->GetBinLowEdge(ibin + 1), hhi->GetBinLowEdge(ibin + 2), (double *)0, (double *)0, 1.e-6);
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
    double cont, err;
    for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
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
    TH1 *hout = (TH1 *)hin->Clone("hout");
    hout->Reset();
    double cont, err, cohe;
    cohe = gRandom->Gaus(0., 1.);
    for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      cont = hin->GetBinContent(ibin + 1);
      err = hin->GetBinError(ibin + 1);
      hout->SetBinContent(ibin + 1, cont + cohe * err);
      hout->SetBinError(ibin + 1, err);
    }
    return hout;
  }

  TH1 *
  YieldMean_ReturnExtremeHighHisto(TH1 *hin, const char* hout_title)
  {
    TH1 *hout = (TH1 *)hin->Clone(hout_title);
    for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      double val = hin->GetBinContent(ibin + 1);
      double err = hin->GetBinError(ibin + 1);
      hout->SetBinContent(ibin + 1, val + err);
      hout->SetBinError(ibin + 1, err);
    }
    return hout;
  }

  TH1 *
  YieldMean_ReturnExtremeLowHisto(TH1 *hin, const char* hout_title)
  {
    TH1 *hout = (TH1 *)hin->Clone(hout_title);
    for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      double val = hin->GetBinContent(ibin + 1);
      double err = hin->GetBinError(ibin + 1);
      hout->SetBinContent(ibin + 1, val - err);
      hout->SetBinError(ibin + 1, err);
    }
    return hout;
  }

  TH1 *
  YieldMean_ReturnExtremeSoftHisto(TH1 *hin, const char* hout_title)
  {
    return YieldMean_ReturnExtremeHisto(hin, hout_title, -1.);
  }

  TH1 *
  YieldMean_ReturnExtremeHardHisto(TH1 *hin, const char* hout_title)
  {
    return YieldMean_ReturnExtremeHisto(hin, hout_title, 1.);
  }

  TH1 *
  YieldMean_ReturnExtremeHisto(TH1 *hin, const char* hout_title, float sign)
  {
    double ptlow = -1., pthigh = -1.;
    for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      ptlow = hin->GetBinLowEdge(ibin + 1);
      break;
    }
    for (int ibin = hin->GetNbinsX(); ibin >= 0; ibin--) {
      if (hin->GetBinError(ibin + 1) <= 0.) continue;
      pthigh = hin->GetBinLowEdge(ibin + 2);
      break;
    }
    if (ptlow < 0. || pthigh < 0.) ::Fatal("YieldMean_ReturnExtremeHisto","Problem in the determination of the pt region");

    double mean = hin->GetMean();
    double maxdiff = 0.;
    TH1 *hmax = nullptr;
    for (int inode = 0; inode < hin->GetNbinsX(); inode++) {

      double ptnode = hin->GetBinCenter(inode + 1);
      TH1 *hout = (TH1 *)hin->Clone("extreme_tmp");

      for (int ibin = 0; ibin < hin->GetNbinsX(); ibin++) {
        if (hin->GetBinError(ibin + 1) <= 0.) continue;
        double val = hin->GetBinContent(ibin + 1);
        double err = hin->GetBinError(ibin + 1);
        double cen = hin->GetBinCenter(ibin + 1);
        double new_err = 0.;
        if (cen < ptnode)
          new_err = err * (-1. + (cen - ptlow) / (ptnode - ptlow));
        else
          new_err *= err * ((cen - ptnode) / (pthigh - ptnode));

        hout->SetBinContent(ibin + 1, val + sign * new_err);
        hout->SetBinError(ibin + 1, err);
      }

      double diff = TMath::Abs(mean - hout->GetMean());
      if (diff > maxdiff) {
        //      printf("found max at %f\n", ptnode);
        if (hmax) delete hmax;
        hmax = (TH1 *)hout->Clone(Form("%s_hmax",hout_title));
        maxdiff = diff;
      }
      delete hout;
    }
    return hmax;
  }

  void YieldMean_IntegralMean(TH1 *hdata, TH1 *hlo, TH1 *hhi, double &integral, double &mean,Bool_t printinfo)
  {

    /*
    * compute integrals
    */

    double cont, err, width, cent;
    double I = 0., IX = 0., Ierr = 0., IXerr = 0., Ilerr = 0., IXlerr = 0.;
    double M = 0., Merr = 0., Mlerr = 0., C;
    double dataonly=0.0;

    /* integrate the data */
    for (int ibin = 0; ibin < hdata->GetNbinsX(); ibin++) {
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
    for (int ibin = 0; ibin < hlo->GetNbinsX(); ibin++) {
      cent = hlo->GetBinCenter(ibin + 1);
      width = hlo->GetBinWidth(ibin + 1);
      cont = width * hlo->GetBinContent(ibin + 1);
      err = width * hlo->GetBinError(ibin + 1);
      if (err <= 0.) continue;
      I += cont;
      IX += cont * cent;
    }
    /* integrate high */
    for (int ibin = 0; ibin < hhi->GetNbinsX(); ibin++) {
      cent = hhi->GetBinCenter(ibin + 1);
      width = hhi->GetBinWidth(ibin + 1);
      cont = width * hhi->GetBinContent(ibin + 1);
      err = width * hhi->GetBinError(ibin + 1);
      if (err <= 0.) continue;
      I += cont;
      IX += cont * cent;
    }

    /* set values */
    integral = I;
    mean = IX / I;
    if(printinfo)
      std::cout<<"data only = "<<dataonly<<" total = "<<I<<" ratio= "<<dataonly/I<<std::endl;
  }

  void RandomShifter(TH1* hin, TH1* hhi, TH1* hlo, TH1*& hIntegral, TH1*& hMean, float integral_limit, float mean_limit, int nRepetitions){

    TH1F hIntegral_tmp("hIntegral_tmp", "", 1000, 0.75 * integral_limit, 1.25 * integral_limit); //NOTE: only used for the range of the next histograms
    TH1F hMean_tmp("hMean_tmp", "", 1000, 0.75 * mean_limit, 1.25 * mean_limit);
    
    double integral, mean;
    
    for (int irnd = 0; irnd < 100; irnd++) {
      /* get random histogram */
      TH1 *hrnd = YieldMean_ReturnRandom(hin); //QUESTION: why not coherent?
      /* fit */
      TH1 *hrndlo = YieldMean_ReturnCoherentRandom(hlo);
      TH1 *hrndhi = YieldMean_ReturnCoherentRandom(hhi);
      /* integrate */
      YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean);
      hIntegral_tmp.Fill(integral);
      hMean_tmp.Fill(mean);
      delete hrnd;
      delete hrndlo;
      delete hrndhi;
    }
    if(hIntegral){
      delete hIntegral;
    }
    hIntegral  = new TH1F(Form("hIntegral_%s",hin->GetName()), ";d#it{N}/d#it{y};", 100, hIntegral_tmp.GetMean() - 10. * hIntegral_tmp.GetRMS(), hIntegral_tmp.GetMean() + 10. * hIntegral_tmp.GetRMS());
    if(hMean){
      delete hMean;
    }
    hMean = new TH1F(Form("hMean_%s",hin->GetName()), ";<#it{p}_{T}> (GeV/#it{c});", 100, hMean_tmp.GetMean() - 10. * hMean_tmp.GetRMS(), hMean_tmp.GetMean() + 10. * hMean_tmp.GetRMS());
    for (int irnd = 0; irnd < nRepetitions; irnd++) {
      /* get random histogram */
      TH1 *hrnd = YieldMean_ReturnRandom(hin); //QUESTION: why not coherent?
      /* fit */
      TH1 *hrndlo = YieldMean_ReturnCoherentRandom(hlo);
      TH1 *hrndhi = YieldMean_ReturnCoherentRandom(hhi);
      /* integrate */
      YieldMean_IntegralMean(hrnd, hrndlo, hrndhi, integral, mean);
      hIntegral->Fill(integral);
      hMean->Fill(mean);
      delete hrnd;
      delete hrndlo;
      delete hrndhi;
    }
  }

  int Fitter(TH1* histo, TF1* func, Option_t *opt){
    int fitres;
    int trials = 0;
    do {
      fitres = histo->Fit(func, opt);
      //Printf("Trial: %d", trials++);
      trials++;
      if(trials > 10) {
        Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
        break;
      }
    }
    while (fitres != 0);
    return fitres;
  }

  void SaveToFile(const char* file_name, const char* dir_name, TCanvas& canvas){
    TFile file(file_name,"UPDATE");
    bool isDir = file.cd(dir_name);
    if(!isDir){
      file.mkdir(dir_name);
      file.cd(dir_name);
    }
    canvas.Write();
    file.Close();
  }
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
        return TMath::Sqrt( TMath::Abs(errorfromtop - errorfrombottom) );
    }
    return 1.;
}