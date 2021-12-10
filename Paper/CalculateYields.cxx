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

TH1 *YieldMean(TH1 *hstat, TH1 *hsys, TF1 *f = NULL, Double_t minfit=0.4,Double_t maxfit=8., TString logfilename="log.root", Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q");
TH1 *YieldMeanNoSyst(TH1 *hstat, TF1 *f = NULL, Double_t minfit = 0.4,Double_t maxfit = 8., Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "0q",TString logfilename="log.root");
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

void CalculateYields(
  TString fWhichEstimator = "ZDC", 
  TString fWhichSelEstimator = "SPDClusters", 
  Double_t lLoSel = 10., 
  Double_t lHiSel = 15.,
  TString period = "17j",
  TString fWhichParticle = "Lambda") { 

    TString fParticle = "";
    TString fAntiParticle = "";
    if (fWhichParticle.Contains("Xi")){
      fParticle = "XiMinus";
      fAntiParticle = "XiPlus";
    }
    if (fWhichParticle.Contains("Omega")){
      fParticle = "OmegaMinus";
      fAntiParticle = "OmegaPlus";
    }  
    if (fWhichParticle.Contains("Lambda")){
      fParticle = "Lambda";
      fAntiParticle = "Lambda";
    }  

    TString outputfilename = Form("yields/Yields%s_SelectedWith%s%.0f%.0f_%s.root",fWhichEstimator.Data(),fWhichSelEstimator.Data(),lLoSel,lHiSel,period.Data());


    Double_t percentile[] = 
    {20,30,40,50,60,80,100};  
   // {0,5,10,15,20,30,40};
    const Long_t binnumber = sizeof(percentile)/sizeof(Double_t) - 1;
    //    
    double errpercentile[binnumber], centrpercentile[binnumber];
    for (Int_t n = 0; n < binnumber; n++) {
      centrpercentile[n] = (percentile[n] + percentile[n + 1]) / 2;
      errpercentile[n] = (percentile[n + 1] - percentile[n]) / 2;
    }

    Double_t dNch[binnumber];
    Double_t StatdNch[binnumber];

    TFile* Nchfile = TFile::Open(Form("NchRawContainer_%s.root",period.Data()));
    //period.Data()));
    TH3D* NchHisto = 0x0;
    if (
    (fWhichEstimator.Contains("SPD") && fWhichSelEstimator.Contains("V0M")) ||
    (fWhichEstimator.Contains("V0M") && fWhichSelEstimator.Contains("SPD"))
    ) NchHisto = (TH3D*)Nchfile->Get("hspd_spdv0m");     
    if (
    (fWhichEstimator.Contains("SPD") && fWhichSelEstimator.Contains("ZDC")) ||
    (fWhichEstimator.Contains("ZDC") && fWhichSelEstimator.Contains("SPD"))
    ) NchHisto = (TH3D*)Nchfile->Get("hspd_spdzdc"); 
    if (
    (fWhichEstimator.Contains("V0M") && fWhichSelEstimator.Contains("ZDC")) ||
    (fWhichEstimator.Contains("ZDC") && fWhichSelEstimator.Contains("V0M"))
    ) NchHisto = (TH3D*)Nchfile->Get("hspd_v0mzdc");
    //

    int miny=1, maxy=1, minz=1, maxz=1;
    for (int i = 0; i < binnumber; i++){ 
      StatdNch[i] = 0.;
      //
      if(fWhichEstimator.Contains("SPD")) {
        miny = NchHisto->GetYaxis()->FindBin(percentile[i]); 
        maxy = NchHisto->GetYaxis()->FindBin(percentile[i+1]-1e-6);
        minz = NchHisto->GetZaxis()->FindBin(lLoSel);
        maxz = NchHisto->GetZaxis()->FindBin(lHiSel);
      }
      if(fWhichEstimator.Contains("ZDC")) {
        miny = NchHisto->GetYaxis()->FindBin(lLoSel);
        maxy = NchHisto->GetYaxis()->FindBin(lHiSel);
        minz = NchHisto->GetZaxis()->FindBin(percentile[i]); 
        maxz = NchHisto->GetZaxis()->FindBin(percentile[i+1]-1e-6);
      }
      if(fWhichEstimator.Contains("V0M") && fWhichSelEstimator.Contains("SPD")) {
        miny = NchHisto->GetYaxis()->FindBin(lLoSel);
        maxy = NchHisto->GetYaxis()->FindBin(lHiSel);
        minz = NchHisto->GetZaxis()->FindBin(percentile[i]); 
        maxz = NchHisto->GetZaxis()->FindBin(percentile[i+1]-1e-6);
      }
      if(fWhichEstimator.Contains("V0M") && fWhichSelEstimator.Contains("ZDC")) {
        miny = NchHisto->GetYaxis()->FindBin(percentile[i]); 
        maxy = NchHisto->GetYaxis()->FindBin(percentile[i+1]-1e-6);
        minz = NchHisto->GetZaxis()->FindBin(lLoSel);
        maxz = NchHisto->GetZaxis()->FindBin(lHiSel);
      }
      //
      TH1D* dummy = NchHisto->ProjectionX(Form("dummy%i",i),miny,maxy,minz,maxz);  
      dNch[i] = dummy->GetMean(); 
    }  

    // Histos
    TH1D* lHistPtPart[binnumber], * lHistPtAntiPart[binnumber],* lHistPtPartSyst[binnumber], * lHistPtAntiPartSyst[binnumber], *lHistPt[binnumber], *lSystPt[binnumber];
    TH1D* hout[binnumber], * houtpart[binnumber], * houtantipart[binnumber];
    Double_t YieldPart[binnumber];
    Double_t YieldPartStat[binnumber];
    Double_t YieldAntiPart[binnumber];
    Double_t YieldAntiPartStat[binnumber];
    Double_t Yield[binnumber];
    Double_t YieldStat[binnumber];
    Double_t AvYield[binnumber];
    Double_t AvYieldStat[binnumber];
    Double_t YieldSysHi[binnumber];
    Double_t YieldSysLo[binnumber];
    Double_t Mean[binnumber];
    Double_t MeanStat[binnumber];
    Double_t MeanSysHi[binnumber];
    Double_t MeanSysLo[binnumber];
    Double_t RatioPAP[binnumber];
    Double_t RatioPAPStat[binnumber];
    //
   
    // Fit function
    TF1* LevyTsallisfunc = LevyTsallis("LevyTsallisfunc", 1.321);
    //
    TFile* lResultsFilePart[binnumber], * lResultsFileAntiPart[binnumber];
    TString prefix = "V0Analysis/results";
    //Form("/home/fercoles/strg_analysis/FullStatistics/RawSpectra/%s",period.Data());

    for (int i = 0; i < binnumber; i++){ 

     if (
       (fWhichEstimator.Contains("SPD")) ||
       (fWhichEstimator.Contains("V0M") && fWhichSelEstimator.Contains("ZDC"))
       )
      {
        lResultsFilePart[i] = TFile::Open(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root",prefix.Data(),fParticle.Data(),fWhichEstimator.Data(),percentile[i],percentile[i+1],fWhichSelEstimator.Data(),lLoSel,lHiSel));
        lResultsFileAntiPart[i] = TFile::Open(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root",prefix.Data(),fAntiParticle.Data(),fWhichEstimator.Data(),percentile[i],percentile[i+1],fWhichSelEstimator.Data(),lLoSel,lHiSel));  
      }
      //
      if (
        (fWhichEstimator.Contains("ZDC")) ||
        (fWhichEstimator.Contains("V0M") && fWhichSelEstimator.Contains("SPD"))
        ) 
      {
        lResultsFilePart[i] = TFile::Open(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root",prefix.Data(),fParticle.Data(),fWhichSelEstimator.Data(),lLoSel,lHiSel,fWhichEstimator.Data(),percentile[i],percentile[i+1]));
        lResultsFileAntiPart[i] = TFile::Open(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root",prefix.Data(),fAntiParticle.Data(),fWhichSelEstimator.Data(),lLoSel,lHiSel,fWhichEstimator.Data(),percentile[i],percentile[i+1]));
      }

      lHistPtPart[i] = (TH1D*)lResultsFilePart[i]->Get(Form("fHistPt%s", fParticle.Data())); 
      lHistPtAntiPart[i] = (TH1D*)lResultsFileAntiPart[i]->Get(Form("fHistPt%s",fAntiParticle.Data())); 
      
      lHistPt[i] = (TH1D*)lHistPtPart[i]->Clone(Form("lHistPt%i",i));
      lHistPt[i]->Reset();
      lSystPt[i] = (TH1D*)lHistPtPart[i]->Clone(Form("lHistPt%i",i));
      lSystPt[i]->Reset();
      lHistPtPartSyst[i] = (TH1D*)lHistPtPart[i]->Clone(Form("lHistPtPartSyst%i",i));
      lHistPtPartSyst[i]->Reset();
      lHistPtAntiPartSyst[i] = (TH1D*)lHistPtPart[i]->Clone(Form("lHistPtAntiPartSyst%i",i));
      lHistPtAntiPartSyst[i]->Reset();

      for(int b = 1; b <= lHistPtPart[i]->GetNbinsX(); b++){
        lHistPt[i]->SetBinContent(b, (lHistPtPart[i]->GetBinContent(b)+lHistPtAntiPart[i]->GetBinContent(b)));
        lHistPt[i]->SetBinError(b, TMath::Sqrt(
            lHistPtPart[i]->GetBinError(b)*lHistPtPart[i]->GetBinError(b) +
            lHistPtAntiPart[i]->GetBinError(b)*lHistPtAntiPart[i]->GetBinError(b) 
            )
          );
        lSystPt[i]->SetBinContent(b, (lHistPtPart[i]->GetBinContent(b)+lHistPtAntiPart[i]->GetBinContent(b)));
        lSystPt[i]->SetBinError(b, lHistPt[i]->GetBinContent(b)*0.06
        /*TMath::Sqrt(
            lHistPtPart[i]->GetBinError(b)*lHistPtPart[i]->GetBinError(b) +
            lHistPtAntiPart[i]->GetBinError(b)*lHistPtAntiPart[i]->GetBinError(b) 
            )*/
          );
        lHistPtPartSyst[i]->SetBinContent(b, lHistPtPart[i]->GetBinContent(b));
        lHistPtPartSyst[i]->SetBinError(b, lHistPtPart[i]->GetBinContent(b)*0.06);
        lHistPtAntiPartSyst[i]->SetBinContent(b, lHistPtAntiPart[i]->GetBinContent(b));
        lHistPtAntiPartSyst[i]->SetBinError(b, lHistPtAntiPart[i]->GetBinContent(b)*0.06);        
      }
     
    }

    // Compute Yields
    for (int i = 0; i < binnumber; i++){ //loop over selection classes         
        hout[i] = (TH1D*) YieldMean(lHistPt[i],lSystPt[i],LevyTsallisfunc,0.4,8.,Form("logfile_%sSelectedWith%s%.0f%.0f_%s.root",fWhichEstimator.Data(),fWhichSelEstimator.Data(),lLoSel,lHiSel,period.Data()));
        // Get yield
        Yield[i] = hout[i]->GetBinContent(kYield);
        YieldStat[i] = hout[i]->GetBinContent(kYieldStat);
        YieldSysHi[i] = hout[i]->GetBinContent(kYieldSysHi);
        YieldSysLo[i] = hout[i]->GetBinContent(kYieldSysLo);
        // Get <pT>
        Mean[i] = hout[i]->GetBinContent(kMean);
        MeanStat[i]= hout[i]->GetBinContent(kMeanStat);
        MeanSysHi[i]= hout[i]->GetBinContent(kMeanSysHi);
        MeanSysLo[i]= hout[i]->GetBinContent(kMeanSysLo);
        //
        houtpart[i] = (TH1D*) YieldMean(lHistPtPart[i],lHistPtPartSyst[i],LevyTsallisfunc,0.4,8.);
        // Get yield part
        YieldPart[i] = houtpart[i]->GetBinContent(kYield);
        YieldPartStat[i] = houtpart[i]->GetBinContent(kYieldStat);
        houtantipart[i] = (TH1D*) YieldMean(lHistPtAntiPart[i],lHistPtAntiPartSyst[i],LevyTsallisfunc,0.4,8.);
        // Get yield anti-part
        YieldAntiPart[i] = houtantipart[i]->GetBinContent(kYield);
        YieldAntiPartStat[i] = houtantipart[i]->GetBinContent(kYieldStat);
        }//end loop over classes 

    for(int i = 0; i < binnumber ; i++){
      AvYield[i] = Yield[i]/dNch[i];
      AvYieldStat[i] = YieldStat[i]/dNch[i];
      //
      RatioPAP[i] = YieldPart[i]/YieldAntiPart[i];
      RatioPAPStat[i] = ErrorInRatio(YieldPart[i],YieldPartStat[i],YieldAntiPart[i],YieldAntiPartStat[i]);
    }
   
    // Write numbers
    std::ofstream outfile;
    outfile.open(Form("Yields%s_SelectedWith%s%.0f%.0f_%s.txt",fWhichEstimator.Data(),fWhichSelEstimator.Data(),lLoSel,lHiSel,period.Data()), std::ofstream::out | std::ofstream::trunc); // std::ios_base::app
    outfile << "\nYields Part + AntiPart\n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << Yield[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Yields Stat Part + AntiPart\n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << YieldStat[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "#LT p_{T} #GT (GeV/c) Part + AntiPart\n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << Mean[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "#LT p_{T} #GT (GeV/c) Stat Part + AntiPart\n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << MeanStat[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "\nYields Particle \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << YieldPart[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Yields Particle Stat \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << YieldPartStat[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "\nYields Anti-Particle \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << YieldPart[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << "Yields Anti-Particle Stat \n{ ";
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
         outfile << YieldPartStat[i] << ", ";
    }
    outfile << " }\n" << endl;
    outfile << " \n\n\n" << endl;

    // Do Plots
    TGraphErrors *YieldsNchStat = new TGraphErrors(binnumber,dNch,Yield,StatdNch,YieldStat);
    TGraphErrors *AvYieldsNchStat = new TGraphErrors(binnumber,dNch,AvYield,StatdNch,AvYieldStat);
    TGraphErrors *MeanptNchStat = new TGraphErrors(binnumber,dNch,Mean,StatdNch,MeanStat);
    TGraphErrors *MeanptPercStat = new TGraphErrors(binnumber,centrpercentile,Mean,errpercentile,MeanStat);
    TGraphErrors *AvYieldsPercStat = new TGraphErrors(binnumber,centrpercentile,AvYield,errpercentile,AvYieldStat);
    TGraphErrors *PartYieldsNchStat = new TGraphErrors(binnumber,dNch,YieldPart,StatdNch,YieldPartStat);
    TGraphErrors *PartYieldsPercStat = new TGraphErrors(binnumber,centrpercentile,YieldPart,errpercentile,YieldPartStat);
    TGraphErrors *AntiPartYieldsNchStat = new TGraphErrors(binnumber,dNch,YieldAntiPart,StatdNch,YieldAntiPartStat);
    TGraphErrors *AntiPartYieldsPercStat = new TGraphErrors(binnumber,centrpercentile,YieldAntiPart,errpercentile,YieldAntiPartStat);
    TGraphErrors *RatioPAPYieldsNchStat = new TGraphErrors(binnumber,dNch,RatioPAP,StatdNch,RatioPAPStat);
    TGraphErrors *RatioPAPYieldsPercStat = new TGraphErrors(binnumber,centrpercentile,RatioPAP,errpercentile,RatioPAPStat);
   

    PartYieldsNchStat->SetName("PartYieldsNchStat");
    PartYieldsNchStat->GetYaxis()->SetTitle("#LT dN/dy #GT (#Xi^{-})");
    PartYieldsNchStat ->GetYaxis()->SetTitleSize(0.05);
    PartYieldsNchStat->GetYaxis()->SetTitleOffset(1.1);
    PartYieldsNchStat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}");
    PartYieldsNchStat->GetXaxis()->SetTitleSize(0.05);
    PartYieldsNchStat->GetXaxis()->SetTitleOffset(0.9);
    PartYieldsPercStat->SetName("PartYieldsPercStat");
    PartYieldsPercStat->GetYaxis()->SetTitle("#LT dN/dy #GT (#Xi^{-}) ");
    PartYieldsPercStat->GetYaxis()->SetTitleSize(0.05);
    PartYieldsPercStat->GetYaxis()->SetTitleOffset(1.1);
    PartYieldsPercStat->GetXaxis()->SetTitle(Form("%s percentile",fWhichEstimator.Data()));
    PartYieldsPercStat->GetXaxis()->SetTitleSize(0.05);
    PartYieldsPercStat->GetXaxis()->SetTitleOffset(0.9);
    AntiPartYieldsNchStat->SetName("AntiPartYieldsNchStat");
    AntiPartYieldsNchStat->GetYaxis()->SetTitle("#LT dN/dy #GT (#bar{#Xi}^{+}})");
    AntiPartYieldsNchStat->GetYaxis()->SetTitleSize(0.05);
    AntiPartYieldsNchStat->GetYaxis()->SetTitleOffset(1.1);
    AntiPartYieldsNchStat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}");
    AntiPartYieldsNchStat->GetXaxis()->SetTitleSize(0.05);
    AntiPartYieldsNchStat->GetXaxis()->SetTitleOffset(0.9);
    AntiPartYieldsPercStat->SetName("AntiPartYieldsPercStat");
    AntiPartYieldsPercStat->GetYaxis()->SetTitle("#LT dN/dy #GT (#bar{#Xi}^{+}})");
    AntiPartYieldsPercStat ->GetYaxis()->SetTitleSize(0.05);
    AntiPartYieldsPercStat ->GetYaxis()->SetTitleOffset(1.1);
    AntiPartYieldsPercStat ->GetXaxis()->SetTitle(Form("%s percentile",fWhichEstimator.Data()));
    AntiPartYieldsPercStat ->GetXaxis()->SetTitleSize(0.05);
    AntiPartYieldsPercStat ->GetXaxis()->SetTitleOffset(0.9);
    RatioPAPYieldsNchStat->SetName("RatioPAPYieldsNchStat");
    RatioPAPYieldsNchStat->GetYaxis()->SetTitle("#frac{#LT dN/dy #GT (#Xi^{-})}{#LT dN/dy #GT (#bar{#Xi}^{+})}");
    RatioPAPYieldsNchStat->GetYaxis()->SetTitleSize(0.05);
    RatioPAPYieldsNchStat->GetYaxis()->SetTitleOffset(1.1);
    RatioPAPYieldsNchStat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}");
    RatioPAPYieldsNchStat->GetXaxis()->SetTitleSize(0.05);
    RatioPAPYieldsNchStat->GetXaxis()->SetTitleOffset(0.9);
    RatioPAPYieldsPercStat->SetName("RatioPAPYieldsPercStat");
    RatioPAPYieldsPercStat->GetYaxis()->SetTitle("#frac{#LT dN/dy #GT (#Xi^{-})}{#LT dN/dy #GT (#bar{#Xi}^{+})} ");
    RatioPAPYieldsPercStat->GetYaxis()->SetTitleSize(0.05);
    RatioPAPYieldsPercStat ->GetYaxis()->SetTitleOffset(1.1);
    RatioPAPYieldsPercStat ->GetXaxis()->SetTitle(Form("%s percentile",fWhichEstimator.Data()));
    RatioPAPYieldsPercStat ->GetXaxis()->SetTitleSize(0.05);
    RatioPAPYieldsPercStat ->GetXaxis()->SetTitleOffset(0.9);

    AvYieldsNchStat->SetName("AvYieldsNchStat");
    AvYieldsNchStat->GetYaxis()->SetTitle("#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}} ");
    AvYieldsNchStat->GetYaxis()->SetTitleSize(0.05);
    AvYieldsNchStat->GetYaxis()->SetTitleOffset(1.1);
    AvYieldsNchStat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}");
    AvYieldsNchStat->GetXaxis()->SetTitleSize(0.05);
    AvYieldsNchStat->GetXaxis()->SetTitleOffset(0.9);
    AvYieldsNchStat->SetTitle("");

    AvYieldsPercStat->SetName("AvYieldsPercStat");
    AvYieldsPercStat->GetYaxis()->SetTitle("#frac{#LT dN/dy #GT}{#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}} ");
    AvYieldsPercStat->GetYaxis()->SetTitleSize(0.05);
    AvYieldsPercStat->GetYaxis()->SetTitleOffset(1.1);
    AvYieldsPercStat->GetXaxis()->SetTitle(Form("%s percentile",fWhichEstimator.Data()));
    AvYieldsPercStat->GetXaxis()->SetTitleSize(0.05);
    AvYieldsPercStat->GetXaxis()->SetTitleOffset(0.9);
    AvYieldsPercStat->SetTitle("");
   
    YieldsNchStat->GetYaxis()->SetTitle("#LT dN/dy #GT");
    YieldsNchStat->GetYaxis()->SetTitleSize(0.05);
    YieldsNchStat->GetYaxis()->SetTitleOffset(1.1);
    YieldsNchStat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}");
    YieldsNchStat->GetXaxis()->SetTitleSize(0.05);
    YieldsNchStat->GetXaxis()->SetTitleOffset(0.9);
    YieldsNchStat->SetTitle("");
    YieldsNchStat->SetName("YieldsNchStat");
    YieldsNchStat->SetMarkerStyle(20);
    YieldsNchStat->SetMarkerColor(kRed);
    YieldsNchStat->SetLineColor(kRed);
    YieldsNchStat->SetMarkerSize(1.8); 
    YieldsNchStat->Draw("SAME ep");

    MeanptNchStat->GetXaxis()->SetRangeUser(0.,28.);
    MeanptNchStat->GetYaxis()->SetRangeUser(0.8,1.7);
    MeanptNchStat->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT");
    MeanptNchStat->GetYaxis()->SetTitleSize(0.05);
    MeanptNchStat->GetYaxis()->SetTitleOffset(0.9);
    MeanptNchStat->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT^{RAW}_{|#eta|<0.5}");
    MeanptNchStat->GetXaxis()->SetTitleSize(0.05);
    MeanptNchStat->GetXaxis()->SetTitleOffset(0.9);
    MeanptNchStat->SetTitle(""); 
    MeanptNchStat->SetMarkerColor(kBlue);
    MeanptNchStat->SetMarkerStyle(20);
    MeanptNchStat->SetMarkerSize(1.2); 
    MeanptNchStat->SetName("MeanptNchStat");
    MeanptNchStat->Draw("SAME EP");

    MeanptPercStat->SetName("MeanptPercStat");
    MeanptPercStat->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT");
    MeanptPercStat->GetYaxis()->SetRangeUser(0.8,1.7);
    MeanptPercStat->GetYaxis()->SetTitleSize(0.05);
    MeanptPercStat->GetYaxis()->SetTitleOffset(1.1);
    MeanptPercStat->GetXaxis()->SetTitle(Form("%s percentile",fWhichEstimator.Data()));
    MeanptPercStat->GetXaxis()->SetTitleSize(0.05);
    MeanptPercStat->GetXaxis()->SetTitleOffset(0.9);
    MeanptPercStat->SetTitle("");
    //

    cout << "\n\n Extracted raw dNch in this selection:\n\n { " ;
    for (int k = 0; k < binnumber;k++){       
        if (k==binnumber-1) {cout << dNch[k] << " }\n" << endl;}
        else cout << dNch[k] << ", ";

    }

    TFile* Write = new TFile(outputfilename.Data(), "recreate");
   
    YieldsNchStat->Write();
    AvYieldsNchStat->Write();
    AvYieldsPercStat->Write();
    MeanptNchStat->Write();
    MeanptPercStat->Write();
    PartYieldsNchStat->Write();
    PartYieldsPercStat->Write();
    AntiPartYieldsNchStat->Write();
    AntiPartYieldsPercStat->Write();
    RatioPAPYieldsNchStat->Write();
    RatioPAPYieldsPercStat->Write(); 
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
YieldMeanNoSyst(TH1 *hstat, TF1 *f = NULL, Double_t minfit=.4,Double_t maxfit=8., Double_t min = 0., Double_t max = 10., Double_t loprecision = 0.01, Double_t hiprecision = 0.1, Option_t *opt = "",TString logfilename="log.root")
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
  /*for (Int_t ibin = 0; ibin < htot->GetNbinsX(); ibin++) {
    htot->SetBinError(ibin + 1, hstat->GetBinError(ibin + 1));
  }*/

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

  return hlo;
}

//////////////////////////////////////////////////////////////////////
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
  fLevyTsallis->SetParLimits(2, 1.e-2, 1.e2);
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
