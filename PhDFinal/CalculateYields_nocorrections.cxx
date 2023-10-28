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

TH1 *YieldMean(TH1 *hstat, TH1 *hsys, TF1 *f = NULL, Double_t minfit=0.4,Double_t maxfit=8., TString logfilename="log.root", Double_t min = 0., Double_t max = 100., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q");
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
void beautifygraph(TGraphErrors* g, TString xtitle, TString ytitle);
TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);
void Init(Int_t lClassCode,
          TString lWhichEnergyEstimator,
          TString lWhichMultEstimator,
          Bool_t lDoMB,
          TString fWhichParticle,
          vector<Double_t> &percentileMult_low,
          vector<Double_t> &percentileMult_high,
          vector<Double_t> &percentileEnergy_low,
          vector<Double_t> &percentileEnergy_high);
void beautifygraph(TGraphErrors *g);

void CalculateYields_nocorrections(
    TString fWhichMultEstimator = "SPDClusters",
    TString fWhichEnergyEstimator = "V0M",
    TString fWhichParticle = "Omega",
    Int_t lClassCode = 4,
    Bool_t lDoMB = kFALSE,
    Bool_t lDoSyst = kFALSE)
{

  //
  TString fParticle = "";
  TString fAntiParticle = "";
  double minfit, maxfit;
  Double_t mass;
  TString folder = "CascadeAnalysis/results";
  Double_t YieldMB, ErrYieldMB;

  TString add = "";
  if (fWhichParticle.Contains("Xi"))
  {
    fParticle = "XiMinus";
    fAntiParticle = "XiPlus";
    minfit = 0.6;
    maxfit = 6.5;
    mass = 1.32171;
    YieldMB = 0.02701;
    ErrYieldMB = 0.000175954;
  }
  if (fWhichParticle.Contains("Omega"))
  {
    fParticle = "OmegaMinus";
    fAntiParticle = "OmegaPlus";
    minfit = .9;
    maxfit = 5.5;
    mass = 1.67245;
    YieldMB = 0.0024070792;
    ErrYieldMB = 0;
  }
  if (fWhichParticle.Contains("Lambda"))
  {
    fParticle = "Lambda";
    fAntiParticle = "AntiLambda";
    minfit = 0.4;
    maxfit = 8.;
    mass = 1.115683;
    folder = "V0Analysis/results";
    add = "_FDUseMCRatio";
    YieldMB = 0.249202;
    ErrYieldMB = 0.0004292;
  }
  if (fWhichParticle.Contains("K0Short"))
  {
    fParticle = "K0Short";
    fAntiParticle = "";
    minfit = 0.;
    maxfit = 10.;
    mass = 0.497611;
    folder = "V0Analysis/results";
    add = "_FDUseMCRatio";
    YieldMB = 0.422277;
    ErrYieldMB = 0.000361306;
  }
  // Init(fWhichEstimator, fWhichSelEstimator, lLoSel, lHiSel, lDoMB, fWhichParticle, percentile, dSystNch, Nch);
  // Double_t percentileMult_low[] = {50, 40, 30, 20, 0};
  // Double_t percentileMult_high[] = {70, 50, 40, 30, 25 };
  // Double_t percentileEnergy_low[] = {20, 30, 40, 50, 60 };
  // Double_t percentileEnergy_high[] = {30, 40, 50, 60, 100 };

  // Percentile
  vector<Double_t> percentileMult_low;
  vector<Double_t> percentileMult_high;
  vector<Double_t> percentileEnergy_low;
  vector<Double_t> percentileEnergy_high;

  Init(lClassCode, fWhichEnergyEstimator, fWhichMultEstimator, lDoMB, fWhichParticle, percentileMult_low, percentileMult_high, percentileEnergy_low, percentileEnergy_high);
  const int binnumber = percentileMult_low.size();

  // Open file
  TString outputfilename;
  if (!lDoMB)
  {
    outputfilename = Form("yields/%sYieldsRAW_%s%s_class%i.root",
                          fWhichParticle.Data(), fWhichMultEstimator.Data(), fWhichEnergyEstimator.Data(), lClassCode);
  }
  else
  {
    outputfilename = Form("yields/%sYieldsRAW_INELgt0.root", fWhichParticle.Data());
  }
  TFile *lResults[binnumber];
  TFile *lAntiResults[binnumber];

  //////////////////////////////////////////////
  ////////////// ESTIMATORS ////////////////////
  //////////////////////////////////////////////

  // centrality
  double errpercentile[binnumber], centrpercentile[binnumber];
  for (Int_t n = 0; n < binnumber; n++)
  {
    centrpercentile[n] = 0;
    errpercentile[n] = 0;
  }

  // Nch ZDC file
  TFile *NchZDCfile = TFile::Open("Selections/NchRawContainer_reference.root");

  Double_t Nch[binnumber];
  Double_t StatNch[binnumber];
  Double_t SystNch[binnumber];
  Double_t NchNormToMB[binnumber];
  Double_t SystNchNormToMB[binnumber];
  Double_t NchMB;
  Double_t SystNchMB;

  Double_t ZDCSum[binnumber];
  Double_t StatZDCSum[binnumber];
  Double_t SystZDCSum[binnumber];
  Double_t ZDCSumNormToMB[binnumber];
  Double_t SystZDCSumNormToMB[binnumber];
  Double_t ZDCSumMB;
  Double_t SystZDCSumMB;

  //
  TH3D *ZNHisto;
  TH3D *NchHisto;
  if (fWhichEnergyEstimator.Contains("ZDC")) {
     if (fWhichMultEstimator.Contains("Clusters")){
       ZNHisto = (TH3D *)NchZDCfile->Get("hznsum_spdzdc");
       NchHisto = (TH3D *)NchZDCfile->Get("hspd_spdzdc");
     }
  } else if (fWhichEnergyEstimator.Contains("V0M")){
    if (fWhichMultEstimator.Contains("Clusters")){
       ZNHisto = (TH3D *)NchZDCfile->Get("hznsum_spdv0m");
       NchHisto = (TH3D *)NchZDCfile->Get("hspd_spdv0m");
    }
  }
  //
  TH1D *dummyznMB = ZNHisto->ProjectionX("dummyznMB", ZNHisto->GetYaxis()->FindBin(0. + 1e-10), ZNHisto->GetYaxis()->FindBin(100.), ZNHisto->GetZaxis()->FindBin(0. + 1e-10), ZNHisto->GetZaxis()->FindBin(100.));
  ZDCSumMB = dummyznMB->GetMean();
  SystZDCSumMB = ZDCSumMB * 0.03;
  //
  TH1D *dummynchMB = NchHisto->ProjectionX("dummynchMB", NchHisto->GetYaxis()->FindBin(0. + 1e-10), NchHisto->GetYaxis()->FindBin(100.), NchHisto->GetZaxis()->FindBin(0. + 1e-10), NchHisto->GetZaxis()->FindBin(100.));
  NchMB = dummynchMB->GetMean();
  SystNchMB = NchMB * 0.03;

  cout << NchMB << endl;

  int miny = 1, maxy = 1, minz = 1, maxz = 1;
  for (int i = 0; i < binnumber; i++)
  {
    //
    StatZDCSum[i] = 0;
    StatNch[i] = 0;

    miny = ZNHisto->GetYaxis()->FindBin(percentileMult_low[i] + 1e-10);
    maxy = ZNHisto->GetYaxis()->FindBin(percentileMult_high[i]);
    minz = ZNHisto->GetZaxis()->FindBin(percentileEnergy_low[i] + 1e-10);
    maxz = ZNHisto->GetZaxis()->FindBin(percentileEnergy_high[i]);

    //
    TH1D *dummyzdc = ZNHisto->ProjectionX(Form("dummyzdc%i", i), miny, maxy, minz, maxz);
    ZDCSum[i] = dummyzdc->GetMean();
    SystZDCSum[i] = ZDCSum[i] * 0.03;
    ZDCSumNormToMB[i] = ZDCSum[i] / ZDCSumMB;
    SystZDCSumNormToMB[i] = ZDCSumNormToMB[i] * TMath::Sqrt(2) * 0.03;
    //
    TH1D *dummynch = NchHisto->ProjectionX(Form("dummynch%i", i), miny, maxy, minz, maxz);
    Nch[i] = dummynch->GetMean();
    SystNch[i] = Nch[i] * 0.03;
    NchNormToMB[i] = Nch[i] / NchMB;
    SystNchNormToMB[i] = NchNormToMB[i] * TMath::Sqrt(2) * 0.03;
  }


  //////////////////////////////////////////////////////////
  ////////////////////// YIELDS ////////////////////////////
  //////////////////////////////////////////////////////////

  TH1D *lHistPt[binnumber], *lSystPt[binnumber];
  TH1D *hout[binnumber];
  Double_t Yield[binnumber];
  Double_t YieldStat[binnumber];
  Double_t YieldSysHi[binnumber];
  Double_t YieldSysLo[binnumber];
  Double_t AvYield[binnumber];
  Double_t AvYieldStat[binnumber];
  Double_t AvYieldSystLo[binnumber];
  Double_t AvYieldSystHi[binnumber];
  Double_t NormYield[binnumber];
  Double_t NormYieldStat[binnumber];
  Double_t NormYieldSystHi[binnumber];
  Double_t NormYieldSystLo[binnumber];
  Double_t Mean[binnumber];
  Double_t MeanStat[binnumber];
  Double_t MeanSysHi[binnumber];
  Double_t MeanSysLo[binnumber];

  TH1D *lAntiHistPt[binnumber], *lAntiSystPt[binnumber];

  //TFile *f = TFile::Open("correctedspectra/CorrSpectra-Omega-13TeV_INELgt0.root");

  // Get spectra
  for (int i = 0; i < binnumber; i++)
  {
    double multmin, multmax, eemin, eemax;
    //
    multmin = percentileMult_low[i];
    multmax = percentileMult_high[i];
    eemin = percentileEnergy_low[i];
    eemax = percentileEnergy_high[i];

    lResults[i] = TFile::Open(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f%s.root",
                                   folder.Data(), fParticle.Data(), fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax, add.Data()));

    if (!fParticle.Contains("K0Short")){
      lAntiResults[i] = TFile::Open(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f%s.root",
                                          folder.Data(), fAntiParticle.Data(), fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax,add.Data()));
    }

    lHistPt[i] = (TH1D *)lResults[i]->Get(Form("fHistPt%s",fParticle.Data()));
    lHistPt[i]->SetName(Form("fHistPt%i", i));
    lSystPt[i] = (TH1D *)lHistPt[i]->Clone(Form("lSystPt%i", i));

    if (!fParticle.Contains("K0Short")){
      lAntiHistPt[i] = (TH1D *)lAntiResults[i]->Get(Form("fHistPt%s",fAntiParticle.Data()));
      lAntiHistPt[i]->SetName(Form("fAntiHistPt%i", i));
      lAntiSystPt[i] = (TH1D *)lAntiHistPt[i]->Clone(Form("lAntiSystPt%i", i));
    }

    // apply systematics
    for (int b = 1; b <= lSystPt[i]->GetNbinsX(); b++)
    {
      lSystPt[i]->SetBinError(b, 0.05 * lSystPt[i]->GetBinContent(b));
      if (!fParticle.Contains("K0Short")){
        lAntiSystPt[i]->SetBinError(b, 0.01 * lAntiSystPt[i]->GetBinContent(b));

        lHistPt[i]->SetBinContent(b, lHistPt[i]->GetBinContent(b) + lAntiHistPt[i]->GetBinContent(b));
        lHistPt[i]->SetBinError(b, TMath::Sqrt(lHistPt[i]->GetBinError(b) * lHistPt[i]->GetBinError(b) + lAntiHistPt[i]->GetBinError(b) * lAntiHistPt[i]->GetBinError(b)));
      }
    }

    /*lHistPt[i] = (TH1D *)f->Get(Form("FinalSpectra/FinalPtSpectrumStat_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
    lSystPt[i] = (TH1D *)lHistPt[i]->Clone(Form("lSystPt%i", i));
    for (int b = 1; b <= lSystPt[i]->GetNbinsX(); b++) {
      lSystPt[i]->SetBinError(b, 0.05 * lSystPt[i]->GetBinContent(b));
    }*/
  }

  // Fit functions
  TF1 *LevyTsallisfunc[binnumber];

  // Compute Yields
  for (int i = 0; i < binnumber; i++)
  { // loop over selection classes
    // Fit function
    LevyTsallisfunc[i] = LevyTsallis(Form("LevyTsallisfunc%i", i), mass);

    cout << lHistPt[i]->GetBinContent(1) << endl;

    hout[i] = (TH1D *)YieldMean(lHistPt[i], lSystPt[i], LevyTsallisfunc[i], minfit, maxfit, Form("yields/logfiles/logfile%s_%s%s_class%i.root", fWhichParticle.Data(), fWhichMultEstimator.Data(), fWhichEnergyEstimator.Data(), lClassCode));

    // Get yield
    Yield[i] = hout[i]->GetBinContent(kYield);
    YieldStat[i] = hout[i]->GetBinContent(kYieldStat);
    YieldSysHi[i] = hout[i]->GetBinContent(kYieldSysHi);
    YieldSysLo[i] = hout[i]->GetBinContent(kYieldSysLo);
    // Get <pT>
    Mean[i] = hout[i]->GetBinContent(kMean);
    MeanStat[i] = hout[i]->GetBinContent(kMeanStat);
    MeanSysHi[i] = hout[i]->GetBinContent(kMeanSysHi);
    MeanSysLo[i] = hout[i]->GetBinContent(kMeanSysLo);
    // Yields / Nch
    AvYield[i] = Yield[i] / Nch[i];
    AvYieldStat[i] = YieldStat[i] / Nch[i];
    AvYieldSystLo[i] = AvYield[i] * TMath::Sqrt(std::pow(YieldSysLo[i] / Yield[i], 2) + std::pow(SystNch[i] / Nch[i], 2));
    AvYieldSystHi[i] = AvYield[i] * TMath::Sqrt(std::pow(YieldSysHi[i] / Yield[i], 2) + std::pow(SystNch[i] / Nch[i], 2));
    // Yields / Nch (Self-normalised to INEL>0)
    NormYield[i] = (Yield[i] / Nch[i]) / (YieldMB / NchMB);
    NormYieldStat[i] = ErrorInRatio(Yield[i], YieldStat[i], YieldMB, ErrYieldMB) * NchMB / Nch[i];

  } // end loop over classes

  // Do Plots
  TGraphErrors *YieldsNchStat = new TGraphErrors(binnumber, Nch, Yield, StatNch, YieldStat);
  TGraphErrors *YieldsPercStat = new TGraphErrors(binnumber, centrpercentile, Yield, errpercentile, YieldStat);
  TGraphErrors *AvYieldsNchStat = new TGraphErrors(binnumber, Nch, AvYield, StatNch, AvYieldStat);
  TGraphErrors *AvYieldsZDCSumStat = new TGraphErrors(binnumber, ZDCSum, AvYield, StatNch, AvYieldStat);
  TGraphErrors *NormYieldsNchStat = new TGraphErrors(binnumber, Nch, NormYield, StatNch, NormYieldStat);
  TGraphErrors *NormYieldsZDCSumStat = new TGraphErrors(binnumber, ZDCSum, NormYield, StatNch, NormYieldStat);
  TGraphErrors *AvYieldsPercStat = new TGraphErrors(binnumber, centrpercentile, AvYield, errpercentile, AvYieldStat);
  TGraphErrors *NormYieldsNormNchStat = new TGraphErrors(binnumber, NchNormToMB, NormYield, StatNch, NormYieldStat);
  TGraphErrors *NormYieldsNormZDCSumStat = new TGraphErrors(binnumber, ZDCSumNormToMB, NormYield, StatNch, NormYieldStat);
  TGraphErrors *MeanptNormNchStat = new TGraphErrors(binnumber, NchNormToMB, Mean, StatNch, MeanStat);
  TGraphErrors *MeanptNormZDCSumStat = new TGraphErrors(binnumber, ZDCSumNormToMB, Mean, StatNch, MeanStat);

  //
  NormYieldsNormNchStat->SetName("NormYieldsNormNchStat");
  NormYieldsNormNchStat->GetYaxis()->SetTitle("(#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}})/(#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}})_{INEL>0}");
  NormYieldsNormNchStat->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");

  NormYieldsNormZDCSumStat->SetName("NormYieldsNormZDCSumStat");
  NormYieldsNormZDCSumStat->GetYaxis()->SetTitle("(#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}})/(#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}})_{INEL>0}");
  NormYieldsNormZDCSumStat->GetXaxis()->SetTitle("ZDCSum/ZDCSum^{MB}");

  //
  MeanptNormNchStat->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT");
  MeanptNormNchStat->GetXaxis()->SetTitle("#LT n_{ch} / n_{ch}^{MB}#GT");
  MeanptNormNchStat->SetName("MeanptNormNchStat");
  //
  MeanptNormZDCSumStat->SetName("MeanptNormZDCSumStat");
  MeanptNormZDCSumStat->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT");
  MeanptNormZDCSumStat->GetXaxis()->SetTitle("#LT ZDC #GT/#LT ZDC #GT_{MB}");
  //
  YieldsNchStat->GetYaxis()->SetTitle("#LT dN/dy #GT");
  YieldsNchStat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}");
  YieldsNchStat->SetName("YieldsNchStat");

  //
  YieldsPercStat->GetYaxis()->SetTitle("#LT dN/dy #GT");
  YieldsPercStat->GetXaxis()->SetTitle(Form("%s percentile", ""));
  YieldsPercStat->SetName("YieldsPercStat");
  beautifygraph(YieldsPercStat);
  //
  AvYieldsNchStat->SetName("AvYieldsNchStat");
  AvYieldsNchStat->GetYaxis()->SetTitle("#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}} ");
  AvYieldsNchStat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}");
  beautifygraph(AvYieldsNchStat);
  //
  //
  AvYieldsZDCSumStat->SetName("AvYieldsZDCSumStat");
  AvYieldsZDCSumStat->GetYaxis()->SetTitle("#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}} ");
  AvYieldsZDCSumStat->GetXaxis()->SetTitle("#LT ZDC Sum #GT^{RAW} (a.u.)");
  beautifygraph(AvYieldsZDCSumStat);
  //
  NormYieldsNchStat->SetName("NormYieldsNchStat");
  NormYieldsNchStat->GetYaxis()->SetTitle("(#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}})/(#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}})_{INEL>0}");
  NormYieldsNchStat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}");
  // beautifygraph(NormYieldsNchStat);
  NormYieldsZDCSumStat->SetName("NormYieldsZDCSumStat");
  NormYieldsZDCSumStat->GetYaxis()->SetTitle("(#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}})/(#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}})_{INEL>0}");
  NormYieldsZDCSumStat->GetXaxis()->SetTitle("#LT ZDC Sum #GT^{RAW} (a.u.)");
  beautifygraph(NormYieldsZDCSumStat);
  //
  AvYieldsPercStat->SetName("AvYieldsPercStat");
  AvYieldsPercStat->GetYaxis()->SetTitle("#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}} ");
  AvYieldsPercStat->GetXaxis()->SetTitle(Form("%s percentile", ""));
  beautifygraph(AvYieldsPercStat);

  TFile *Write = new TFile(outputfilename.Data(), "recreate");
  YieldsNchStat->Write();
  AvYieldsNchStat->Write();
  AvYieldsZDCSumStat->Write();
  NormYieldsNchStat->Write();
  NormYieldsZDCSumStat->Write();
  MeanptNormNchStat->Write();
  MeanptNormZDCSumStat->Write();
  NormYieldsNormNchStat->Write();
  NormYieldsNormZDCSumStat->Write();

  Write->cd();
}

void beautifygraph(TGraphErrors *g)
{
  g->GetYaxis()->SetTitleSize(0.05);
  g->GetYaxis()->SetTitleOffset(1.1);
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetTitleOffset(1.1);
  g->SetTitle("");
}

void Init(Int_t lClassCode,
		  TString lWhichEnergyEstimator,
		  TString lWhichMultEstimator,
		  Bool_t lDoMB,
		  TString fWhichParticle,
		  vector<Double_t> &percentileMult_low,
		  vector<Double_t> &percentileMult_high,
		  vector<Double_t> &percentileEnergy_low,
		  vector<Double_t> &percentileEnergy_high)
{
	// class 0 --> standalone
	Double_t percentileEnergy_low_0[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
	Double_t percentileEnergy_high_0[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
	Double_t percentileMult_low_0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t percentileMult_high_0[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
	Long_t n0 = sizeof(percentileEnergy_low_0) / sizeof(Double_t);

  // class 10 --> kStandalone for Omegas
  Double_t percentileEnergy_low_10[] = {0, 5, 15, 30, 50};
  Double_t percentileEnergy_high_10[] = {5, 15, 30, 50, 100};
  Double_t percentileMult_low_10[] = {0, 0, 0, 0, 0};
  Double_t percentileMult_high_10[] = {100, 100, 100, 100, 100};
  Long_t n10 = sizeof(percentileEnergy_low_10) / sizeof(Double_t);

  // class 11 --> kHighMult
  Double_t percentileMult_low_11[] = {10, 10, 10, 10};
  Double_t percentileMult_high_11[] = {20, 20, 20, 20};
  Double_t percentileEnergy_low_11[] = {0, 10, 30, 50};
  Double_t percentileEnergy_high_11[] = {10, 30, 50, 100};
  Long_t n11 = sizeof(percentileMult_low_11) / sizeof(Double_t);

  // class 12 --> kHighMult
  Double_t percentileMult_low_12[] = {40, 40, 40, 40};
  Double_t percentileMult_high_12[] = {50, 50, 50, 50};
  Double_t percentileEnergy_low_12[] = {0, 20, 40, 70};
  Double_t percentileEnergy_high_12[] = {20, 40, 70, 100};
  Long_t n12 = sizeof(percentileMult_low_12) / sizeof(Double_t);

  // class 1 --> kHighMult
  Double_t percentileMult_low_1[] = {10, 10, 10, 10, 10, 10, 10};
  Double_t percentileMult_high_1[] = {20, 20, 20, 20, 20, 20, 20};
  Double_t percentileEnergy_low_1[] = {0, 5, 10, 20, 30, 40, 50};
  Double_t percentileEnergy_high_1[] = {5, 10, 20, 30, 40, 50, 100};
  Long_t n1 = sizeof(percentileMult_low_1) / sizeof(Double_t);

  // class 2 --> kLowMult
  Double_t percentileMult_low_2[] = {40, 40, 40, 40, 40, 40, 40};
  Double_t percentileMult_high_2[] = {50, 50, 50, 50, 50, 50, 50};
  Double_t percentileEnergy_low_2[] = {0, 20, 30, 40, 50, 60, 70};
  Double_t percentileEnergy_high_2[] = {20, 30, 40, 50, 60, 70, 100};
  Long_t n2 = sizeof(percentileMult_low_2) / sizeof(Double_t);

  // class 3 --> fixed high ZN
  Double_t percentileMult_low_3[] = {10, 40, 60, 70, 80};
  Double_t percentileMult_high_3[] = {40, 60, 70, 80, 100};
  Double_t percentileEnergy_low_3[] = {70, 60, 40, 40, 40};
  Double_t percentileEnergy_high_3[] = {100, 100, 100, 80, 70};
  Long_t n3 = sizeof(percentileMult_low_3) / sizeof(Double_t);

  // class 4 --> fixed low ZN
  Double_t percentileMult_low_4[] = {0, 10, 20, 30, 50};
  Double_t percentileMult_high_4[] = {20, 30, 40, 50, 100};
  Double_t percentileEnergy_low_4[] = {40, 30, 30, 20, 0};
  Double_t percentileEnergy_high_4[] = {60, 70, 50, 50, 30};
  Long_t n4 = sizeof(percentileMult_low_4) / sizeof(Double_t);

  // class 5 --> fixed very low ZN
  Double_t percentileMult_low_5[] = {0, 10, 20, 30};
  Double_t percentileMult_high_5[] = {10, 20, 30, 50};
  Double_t percentileEnergy_low_5[] = {20, 10, 0, 0};
  Double_t percentileEnergy_high_5[] = {30, 30, 20, 10};
  Long_t n5 = sizeof(percentileMult_low_5) / sizeof(Double_t);

  if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 0) {
		if (!lDoMB) {
			for (Int_t i = 0; i < n0; i++) {
				percentileMult_low.push_back(percentileMult_low_0[i]);
				percentileMult_high.push_back(percentileMult_high_0[i]);
				percentileEnergy_low.push_back(percentileEnergy_low_0[i]);
				percentileEnergy_high.push_back(percentileEnergy_high_0[i]);
			}
		} else {
        for (Int_t i = 0; i < 2; i++){
				percentileMult_low.push_back(0.);
				percentileMult_high.push_back(100.);
				percentileEnergy_low.push_back(0.);
				percentileEnergy_high.push_back(100.);
		}
      }
	}

  if (lClassCode == 10)
  {
      for (Int_t i = 0; i < n10; i++)
      {
    percentileMult_low.push_back(percentileMult_low_10[i]);
    percentileMult_high.push_back(percentileMult_high_10[i]);
    percentileEnergy_low.push_back(percentileEnergy_low_10[i]);
    percentileEnergy_high.push_back(percentileEnergy_high_10[i]);
      }
      cout << "------------------------------------------" << endl;
      cout << " Initializing kStandalone class for Omegas..." << endl;
  }

  if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 1) {
		for (Int_t i = 0; i < n1; i++) {
			percentileMult_low.push_back(percentileMult_low_1[i]);
			percentileMult_high.push_back(percentileMult_high_1[i]);
			percentileEnergy_low.push_back(percentileEnergy_low_1[i]);
			percentileEnergy_high.push_back(percentileEnergy_high_1[i]);
		}
	}

	if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 2) {
		for (Int_t i = 0; i < n2; i++) {
			percentileMult_low.push_back(percentileMult_low_2[i]);
			percentileMult_high.push_back(percentileMult_high_2[i]);
			percentileEnergy_low.push_back(percentileEnergy_low_2[i]);
			percentileEnergy_high.push_back(percentileEnergy_high_2[i]);
		}
	}

  if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 3) {
		for (Int_t i = 0; i < n3; i++) {
			percentileMult_low.push_back(percentileMult_low_3[i]);
			percentileMult_high.push_back(percentileMult_high_3[i]);
			percentileEnergy_low.push_back(percentileEnergy_low_3[i]);
			percentileEnergy_high.push_back(percentileEnergy_high_3[i]);
		}
	}

  if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 4) {
		for (Int_t i = 0; i < n4; i++) {
			percentileMult_low.push_back(percentileMult_low_4[i]);
			percentileMult_high.push_back(percentileMult_high_4[i]);
			percentileEnergy_low.push_back(percentileEnergy_low_4[i]);
			percentileEnergy_high.push_back(percentileEnergy_high_4[i]);
		}
	}

  if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 5) {
		for (Int_t i = 0; i < n5; i++) {
			percentileMult_low.push_back(percentileMult_low_5[i]);
			percentileMult_high.push_back(percentileMult_high_5[i]);
			percentileEnergy_low.push_back(percentileEnergy_low_5[i]);
			percentileEnergy_high.push_back(percentileEnergy_high_5[i]);
		}
	}

  if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 11)
  {
    for (Int_t i = 0; i < n11; i++)
    {
      percentileMult_low.push_back(percentileMult_low_11[i]);
      percentileMult_high.push_back(percentileMult_high_11[i]);
      percentileEnergy_low.push_back(percentileEnergy_low_11[i]);
      percentileEnergy_high.push_back(percentileEnergy_high_11[i]);
    }
  }

  if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 12)
  {
    for (Int_t i = 0; i < n12; i++)
    {
      percentileMult_low.push_back(percentileMult_low_12[i]);
      percentileMult_high.push_back(percentileMult_high_12[i]);
      percentileEnergy_low.push_back(percentileEnergy_low_12[i]);
      percentileEnergy_high.push_back(percentileEnergy_high_12[i]);
    }
  }
}

void beautifygraph(TGraphErrors* g, TString xtitle, TString ytitle){
    g->GetYaxis()->SetTitleSize(0.05);
    g->GetYaxis()->SetTitleOffset(1.1);
    g->GetXaxis()->SetTitleSize(0.05);
    g->GetXaxis()->SetTitleOffset(1.1);
    g->SetTitle("");
    g->GetXaxis()->SetTitle(xtitle);
    g->GetYaxis()->SetTitle(ytitle);
}

TH1 *
YieldMean(TH1 *hstat, TH1 *hsys, TF1 *f = NULL, Double_t minfit=.4,Double_t maxfit=8., TString logfilename="log.root",Double_t min = 0., Double_t max = 100., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q")
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
  TH1 *htot = (TH1 *)hstat->Clone(Form("%sfittedwith%s",hstat->GetName(),"Func"));
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
    if(trials > 10) {
      Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
      break;
    }
  }
  while (fitres != 0);
  TFile* filewithfits=TFile::Open(logfilename.Data(),"UPDATE");
  htot->Write();
  f->SetName(Form("Levyfitto%s",hstat->GetName()));
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
  fLevyTsallis->SetParLimits(2, 1.e-5, 1.e5);//4
  fLevyTsallis->SetParLimits(3, 1.e-4, 1.e4);//5
  return fLevyTsallis;
}
//------------------------------------------------------

//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {

    double lSigmaDelta = TMath::Sqrt( TMath::Abs( TMath::Power(Aerr,2) - TMath::Power(Berr,2) ) );
    //Computation of relationship to h2 for plotting in ratio plot
    if ( B > 1e-12 ){
      lSigmaDelta /= B;
    }else{
      lSigmaDelta = 0;
    }
  return lSigmaDelta;
}
