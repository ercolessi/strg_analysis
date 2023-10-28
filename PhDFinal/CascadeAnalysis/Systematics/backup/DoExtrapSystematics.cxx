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
#include "AliPWGFunc.h"


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

void DoExtrap(
  TString Func = "BlastWave",
  TString inputfile = "",
  Double_t * percentile = 0x0,
  Long_t binnumber = 1,
  TString outputfile = "",
  TString lCascType = "Xi",
  TString fWhichEstimator = "",
  Double_t lFixedLo = 0.0,
  Double_t lFixedHi = 100.0
  );

  void DoExtrapSystematics(
    TString lCascType = "Xi",
    Double_t lFixedLo = 40.0,
    Double_t lFixedHi = 50.0,
    TString lWhichVarEstimator = "SPDClusters",
    TString lWhichFixedEstimator = "V0M",
    Bool_t DoMB = kFALSE
  ){

   //gROOT->SetBatch(kTRUE);
    //Percentile
    Double_t * percentile;
    Long_t percbinnumb;
    Double_t pmb[] = {0,100};
    Long_t nmb = sizeof(pmb)/sizeof(Double_t) - 1;
    Double_t p0[] = {0,5,10,15,20,30,40,50,70,100};
    Long_t n0 = sizeof(p0)/sizeof(Double_t) - 1;
    Double_t p1[] = {0,5,10,20,30,40,50,100};
    Long_t n1 = sizeof(p1)/sizeof(Double_t) - 1;
    Double_t p2[] = {0,20,30,40,50,60,70,100};
    Long_t n2 = sizeof(p2)/sizeof(Double_t) - 1; 
    Double_t p4[] = {0,5,10,20,30,40,50,100};
    Long_t n4 = sizeof(p4)/sizeof(Double_t) - 1;
    Double_t p5[] = {0,10,20,30,40,50,60,70,100};
    Long_t n5 = sizeof(p5)/sizeof(Double_t) - 1;
    Double_t pOmega[] = {0,5,10,30,50,100};
    Long_t npOmega = sizeof(pOmega)/sizeof(Double_t) - 1;
    Double_t pOmega2[] = {0,40,70,100};
    Long_t npOmega2 = sizeof(pOmega2)/sizeof(Double_t) - 1;

    if (DoMB){
      if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 0. && lFixedHi == 100.){
        percentile = pmb;
        percbinnumb = nmb;
      }
    } else{
      if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 0. && lFixedHi == 100.){
        if (lCascType.Contains("Omega")){
          percentile = pOmega;
          percbinnumb = npOmega;
        } else{
          percentile = p0;
          percbinnumb = n0;
        }
      }
      if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
        if (lCascType.Contains("Omega")){
          percentile = pOmega;
          percbinnumb = npOmega;
        } else{
          percentile = p1;
          percbinnumb = n1;
        }
      }
      if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        if (lCascType.Contains("Omega")){
          percentile = pOmega2;
          percbinnumb = npOmega2;
        } else{
        percentile = p2;
        percbinnumb = n2;
        }
      }
      if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        percentile = p4;
        percbinnumb = n4;
      }
      if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
        percentile = p5;
        percbinnumb = n5;
      }
    }

    TString inputfilename = "";
    if (DoMB) {
      inputfilename = Form("../correctedresults/CorrSpectra-%s-13TeV_INELgt0.root",lCascType.Data());
    } else {
      inputfilename = Form("../correctedresults/CorrSpectra-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", 
                                lCascType.Data(), lWhichVarEstimator.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi);
    }

    TString outputfilename = "";
    if (DoMB) {
      outputfilename = Form("ExtrSyst-%s-13TeV_INELgt0.root",lCascType.Data());
    } else {
      outputfilename = Form("ExtrSyst-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", 
                                lCascType.Data(), lWhichVarEstimator.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi);
    }

    if (1){
    DoExtrap( "Levy", inputfilename, percentile, percbinnumb, outputfilename, lCascType, lWhichVarEstimator ,lFixedLo,lFixedHi  );
    DoExtrap( "BlastWave", inputfilename, percentile, percbinnumb, outputfilename, lCascType, lWhichVarEstimator ,lFixedLo,lFixedHi );
    DoExtrap( "Boltz", inputfilename, percentile, percbinnumb, outputfilename, lCascType, lWhichVarEstimator ,lFixedLo,lFixedHi );
    DoExtrap( "MTexpo", inputfilename, percentile, percbinnumb, outputfilename, lCascType, lWhichVarEstimator ,lFixedLo,lFixedHi );
    DoExtrap( "FermiDirac", inputfilename, percentile, percbinnumb, outputfilename, lCascType, lWhichVarEstimator ,lFixedLo,lFixedHi );
    DoExtrap( "BoseEinstein", inputfilename, percentile, percbinnumb, outputfilename, lCascType, lWhichVarEstimator ,lFixedLo,lFixedHi );
    } else {
      DoExtrap( "BlastWave", inputfilename, percentile, percbinnumb, outputfilename, lCascType, lWhichVarEstimator ,lFixedLo,lFixedHi );
    }
 
  }

void DoExtrap(
  TString Func,
  TString inputfile,
  Double_t * percentile,
  Long_t binnumber,
  TString outputfile,
  TString lCascType,
  TString fWhichEstimator ,
  Double_t lFixedLo = 0.0,
    Double_t lFixedHi = 100.0
  ) {

    

    TString fWhichSelEstimator = "";
    if (fWhichEstimator.Contains("V0M"))fWhichSelEstimator = "SPDClusters";
    if (fWhichEstimator.Contains("SPDClusters"))fWhichSelEstimator = "V0M";

    // Load common libraries
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    gSystem->Load("libOADB");
    gSystem->Load("libPWGLFSTRANGENESS");
    // Use AliRoot includes to compile our task
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

    //Argument: System analysed - one of "XiMinus", "XiPlus", "OmegaMinus"
    cout<<" ---> Macro to compute systematics on extrapolation... "<<endl;
    cout<<"----------------------------------------------------"<<endl;
    cout<<endl;
    cout<<"----------------------------------------------------"<<endl;
    cout<<" ---> Compiling needed class, please wait... "<<endl;
  
    //Load Class
    gSystem->Load("AliPWGFunc_cxx.so");
    //

    double minfit, maxfit;
    Double_t mass = 1.32171;

    if (lCascType.Contains("Xi")){
        minfit = 0.6;
        maxfit = 6.5;
        //mass is default
    }
    if (lCascType.Contains("Omega")){
        minfit = 0.9;
        maxfit = 5.5;
        mass = 1.67245;
    }
    if (lCascType.Contains("Lambda")){
        minfit = 0.4;
        maxfit = 8.;
        mass = 1.115683;
    }   
   
    //Systematics Files
    TFile* lResultsFile = TFile::Open(inputfile);
    TH1D* lHistPt[binnumber];
    TH1D* lSystPt[binnumber];
    TH1D* lHistPtClone[binnumber];
    TString folder = "NormCorr/AfterSignalLossCorr";
    TFile* filesyst = TFile::Open(Form("SystematicsFinalResults-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root",lCascType.Data(),0.,100.,0.,100.));
    TH1F* hsyst = (TH1F*)filesyst->Get("hSystTot");
    //Get syst on signal and event loss 
    TFile* filesgnlosssyst = TFile::Open(
      //Form("/home/fercoles/strg_analysis/Analisi2022/Systematics/SgnLossSystematics-%s-13TeV_INELgt0.root",lCascType.Data()));
      Form("/home/fercoles/strg_analysis/Analisi2022/Systematics/SgnLossSystematics-%s_%s_Fixed%sin%03.0f_%03.0f.root", lCascType.Data(), fWhichEstimator.Data(), fWhichSelEstimator.Data(), lFixedLo, lFixedHi));
    TH1F* hsyssgnloss[binnumber];
    for (int i = 0; i<binnumber; i++){
      if (lCascType.Contains("Xi") ) hsyssgnloss[i] = (TH1F*)filesgnlosssyst->Get(Form("var_%.0f-%0.f",percentile[i],percentile[i+1]));
      if (lCascType.Contains("Omega") ) hsyssgnloss[i] = (TH1F*)filesgnlosssyst->Get(Form("var_%.0f-%0.f",0.,100.));
    }
      
    for (int i = 0; i < binnumber; i++){           
      double multmin, multmax, eemin, eemax;
      //
      if (fWhichEstimator.Contains("SPD")){
          multmin = percentile[i];
          multmax = percentile[i+1];
          eemin = lFixedLo; 
          eemax = lFixedHi;
      }
      //
      if (fWhichEstimator.Contains("V0M") || fWhichEstimator.Contains("ZDC")){
          eemin = percentile[i];
          eemax = percentile[i+1];
          multmin = lFixedLo; 
          multmax = lFixedHi;
      }
      lHistPt[i] = (TH1D*)lResultsFile->Get(Form("%s/PtSpectrumCorr_%s_%.0f_%.0f_%s_%.0f_%.0f",folder.Data(),"SPDClusters",multmin,multmax,"V0M",eemin,eemax)); 
      lHistPt[i]->SetLineColor(kBlack);
      lHistPt[i]->SetMarkerColor(kBlack);
      lHistPt[i]->SetMarkerStyle(24);
      lHistPt[i]->SetMarkerSize(1.2);
      lHistPt[i]->SetName(Form("lHistPt%i",i));
      lHistPtClone[i] = (TH1D*) lHistPt[i]->Clone(Form("lHistClonePt%i",i));      
      lSystPt[i] = (TH1D*)lResultsFile->Get(Form("%s/PtSpectrumCorr_%s_%.0f_%.0f_%s_%.0f_%.0f",folder.Data(),"SPDClusters",multmin,multmax,"V0M",eemin,eemax));
      lSystPt[i]->SetName(Form("lSystPt%i",i));
      int start = hsyssgnloss[i]->FindBin(minfit)-1;
      for (int b = 1; b<(lSystPt[i]->GetNbinsX()+1); b++){
          double s1 = hsyst->GetBinContent(b);
          double s2 = hsyssgnloss[i]->GetBinContent(b+start);
          lSystPt[i]->SetBinError(b, TMath::Sqrt(s1*s1+s2*s2)*lSystPt[i]->GetBinContent(b));
      }
    }

    //Initialize Object
    AliPWGFunc *PWGFunc = new AliPWGFunc();
    //chose the variable you want to use calling AliPWGFunc::SetVarType with one of the elements of the VarType_t enum
    PWGFunc->SetVarType(AliPWGFunc::kdNdpt);

    //Functions to Fit
    TF1* FFit[binnumber];
    //
    for (int i = 0; i < binnumber; i++){ //loop over selection classes   
        if (Func.Contains("Boltz")) FFit[i] = PWGFunc->GetBoltzmann( mass, 0.3, 1.,  Form("%s_%.0f-%.0f","Boltz",percentile[i],percentile[i+1]));
        if (Func.Contains("BlastWave")) FFit[i] = PWGFunc->GetBGBW( mass, .8, 0.3, 1., 5.,  Form("%s_%.0f-%.0f","BlastWave",percentile[i],percentile[i+1]));//};
        if (Func.Contains("MTexpo")) FFit[i] = PWGFunc->GetMTExp( mass, 1., 1., Form("%s_%.0f-%.0f","MTexpo",percentile[i],percentile[i+1]));
        if (Func.Contains("FermiDirac")) FFit[i] = PWGFunc->GetFermiDirac( mass,  .1,  1.,  Form("%s_%.0f-%.0f","FermiDirac",percentile[i],percentile[i+1]));
        if (Func.Contains("BoseEinstein")) FFit[i] = PWGFunc->GetBoseEinstein( mass, .1, 1.,  Form("%s_%.0f-%.0f","BoseEinstein",percentile[i],percentile[i+1]));
        if (Func.Contains("Levy")) FFit[i] = LevyTsallis(Form("%s_%.0f-%.0f","Levy",percentile[i],percentile[i+1]), mass);

        FFit[i]->SetName(Form("%s_%.0f-%.0f",Func.Data(),percentile[i],percentile[i+1]));

    }

    //Do the magic
    TH1D* hout[binnumber];
    Double_t Yield[binnumber], YieldStat[binnumber], Chi[binnumber]; //Levy is [i][0]
    double min[binnumber], max[binnumber];
    for (int i = 0; i < binnumber; i++){ //defaults fit
      min[i] = minfit;
      max[i] = maxfit;
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.) max[2] = 4.;
    }

    if (Func.Contains("Levy")) {
      for (int i = 0; i < binnumber; i++){ 
        if (lCascType.Contains("Xi") && fWhichEstimator.Contains("SPDClusters") && lFixedLo == 40. && lFixedHi == 50.){
        max[0] = 6.5;
        max[1] = 6.5;
        max[2] = 6.5;
        max[3] = 6.5;
        max[4] = 6.5;
        max[5] = 6.5;
        max[6] = 6.5;      
      }
      }
    }

    if (Func.Contains("BlastWave")) {
      for (int i = 0; i < binnumber; i++){ 
        if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
        max[0] = 5.;
        max[1] = 2.5;
        max[2] = 2.5;
        /*max[3] = 5.;
        max[4] = 5.;
        max[5] = 5;
        max[6] = 5.;*/    
      }
        if (lCascType.Contains("Xi") && fWhichEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        max[0] = 6.5;
        max[1] = 6.5;
        max[2] = 6.5;
        max[3] = 6.5;
        max[4] = 6.5;
        max[5] = 6.5;
        max[6] = 6.5;
        max[7] = 5.;      
      }
      }
    }

    if (Func.Contains("Boltz")) {
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.) {
        max[0] = 3.;
        max[1] = 3.;
      }
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 0. && lFixedHi == 100.){
        max[0] = 5;  
        max[1] = 3;   
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 0. && lFixedHi == 100.){
        max[0] = 3;   
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 2.5;
        max[1] = 3.;
        max[2] = 3.;
        max[3] = 3.;
        max[4] = 3.;
        max[5] = 3.;
        max[6] = 3.;      
      }
       if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
        max[0] = 2.;
        max[1] = 2.;
        max[2] = 2.;
        max[3] = 2.;
        max[4] = 2.;
        max[5] = 2.;
        max[6] = 2.;      
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        max[0] = 2.5;
        max[1] = 2.;
        max[2] = 2.5;
        max[3] = 2.5;
        max[4] = 2.;
        max[5] = 2.;
        max[6] = 2.;
        max[7] = 2.;      
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 3.5;
        max[1] = 3.5;
        max[2] = 3.5;
        max[3] = 3.;
        max[4] = 3.;
        max[5] = 3.;
        max[6] = 2.5;
        max[7] = 2.5;      
      }
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 2.5;
        max[1] = 5.;
        max[2] = 5.;
        max[3] = 4.;
        max[4] = 4.;
        max[5] = 4.;
      }
    }
    
    if (Func.Contains("MTexpo")) {
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.) {
        max[0] = 3.;
        max[1] = 3.;
      }
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 0. && lFixedHi == 100.){
        max[0] = 3;
        max[1] = 3;   
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 0. && lFixedHi == 100.){
        max[0] = 3;   
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 4.;
        max[1] = 4.;
        max[2] = 4.;
        max[3] = 4.;
        max[4] = 3.;
        max[5] = 3.;
        max[6] = 4.;      
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
        max[0] = 2.;
        max[1] = 2.;
        max[2] = 2.;
        max[3] = 2.;
        max[4] = 2.;
        max[5] = 2.;
        max[6] = 2.;      
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        max[0] = 2.5;
        max[1] = 2.;
        max[2] = 2.5;
        max[3] = 2.5;
        max[4] = 2.5;
        max[5] = 2.5;
        max[6] = 2.5;
        max[7] = 2.2;      
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 3.5;
        max[1] = 3.5;
        max[2] = 3.5;
        max[3] = 3.5;
        max[4] = 3.5;
        max[5] = 3.5;
        max[6] = 2.5;
        max[7] = 3.;      
      }
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 2.5;
        max[1] = 5.;
        max[2] = 5.;
        max[3] = 4.;
        max[4] = 4.;
      }
    }

    if (Func.Contains("FermiDirac")) {
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.) {
        max[0] = 3.;
        max[1] = 3.;
      }
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 0. && lFixedHi == 100.){
        max[0] = 4; 
        max[1] = 4;   
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 0. && lFixedHi == 100.){
        max[0] = 3;   
      }
        if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 4.;
        max[1] = 4.;
        max[2] = 4.;
        max[3] = 4.;
        max[4] = 3.5;
        max[5] = 3.;
        max[6] = 4.;      
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
        max[0] = 2.5;
        max[1] = 2.5;
        max[2] = 2.5;
        max[3] = 2.5;
        max[4] = 2.;
        max[5] = 2.;
        max[6] = 2.;      
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        max[0] = 2.5;
        max[1] = 2.5;
        max[2] = 2.5;
        max[3] = 2.5;
        max[4] = 2.5;
        max[5] = 2.5;
        max[6] = 2.5;
        max[7] = 2.2;      
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 3.5;
        max[1] = 3.5;
        max[2] = 3.5;
        max[3] = 3.5;
        max[4] = 3.5;
        max[5] = 3.5;
        max[6] = 2.5;
        max[7] = 2.2;      
      }
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 2.5;
        max[2] = 5.;
        max[3] = 4.;
        max[4] = 4.;
        max[5] = 4.;
      }
    }

    if (Func.Contains("BoseEinstein")) {
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.) {
        max[0] = 3.5;
        max[1] = 3.;
      }
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 0. && lFixedHi == 100.){
        max[0] = 4;
        max[1] = 4;    
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 0. && lFixedHi == 100.){
        max[0] = 3.5;   
      }
        if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 4.;
        max[1] = 4.;
        max[2] = 4.;
        max[3] = 4.;
        max[4] = 3.5;
        max[5] = 3.;
        max[6] = 4.;      
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
        max[0] = 2.5;
        max[1] = 2.5;
        max[2] = 2.5;
        max[3] = 2.5;
        max[4] = 2.;
        max[5] = 2.;
        max[6] = 2.;      
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        max[0] = 2.5;
        max[1] = 2.5;
        max[2] = 2.5;
        max[3] = 2.5;
        max[4] = 2.5;
        max[5] = 2.5;
        max[6] = 2.5;
        max[7] = 2.5;      
      }
      if (lCascType.Contains("Xi") && fWhichEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 3.5;
        max[1] = 3.5;
        max[2] = 3.5;
        max[3] = 3.5;
        max[4] = 3.5;
        max[5] = 3.5;
        max[6] = 2.5;
        max[7] = 2.5;      
      }
      if (lCascType.Contains("Omega") && fWhichEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        max[0] = 2.5;
        max[1] = 5.;
        max[2] = 5.;
        max[3] = 4.;
        max[4] = 4.;
      }
    }
    
    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
   
      hout[i] = (TH1D*)YieldMean(lHistPt[i],lSystPt[i],FFit[i],min[i],max[i]);
      // Get yield
      Yield[i] = hout[i]->GetBinContent(kYield);
      YieldStat[i] = hout[i]->GetBinContent(kYieldStat);
      // Get Chi/ndf from fit
      Chi[i] = hout[i]->GetBinContent(kChi);
    }//end loop over classes

    for (int i = 0; i < binnumber; i++){ //loop over selection classes 
        if (Func.Contains("Boltz")) FFit[i]->SetLineColor(kOrange);
        if (Func.Contains("BlastWave")) FFit[i]->SetLineColor(kYellow+1);
        if (Func.Contains("MTexpo"))FFit[i] ->SetLineColor(kGreen+3);
        if (Func.Contains("FermiDirac"))FFit[i]->SetLineColor(kBlue);
        if (Func.Contains("BoseEinstein")) FFit[i]->SetLineColor(kRed);
        if (Func.Contains("Levy")) FFit[i]->SetLineColor(kBlack);
        FFit[i]->SetLineWidth(2);
    }

    Double_t n[binnumber];
    TF1* funcdraw[binnumber];
    for (int i = 0; i < binnumber; i++){  
        n[i] = std::pow(2,binnumber-1-i);
        double par = FFit[i]->GetParameter("norm");
        FFit[i]->SetParameter("norm",par*n[i]);
        FFit[i]->SetRange(0.,max[i]);
        lHistPtClone[i]->Scale(n[i]);
        cout << max[i] << endl;
    }
   
    TCanvas* levy = new TCanvas("levy", "", 800,1500);
    levy->SetRightMargin(0.09);
    levy->SetLeftMargin(0.15);
    levy->SetBottomMargin(0.15);
    levy->SetFillColor(kWhite);
    levy->SetLogy();
    lHistPtClone[0]->GetYaxis()->SetRangeUser(1E-7,1000);
    lHistPtClone[0]->GetXaxis()->SetRangeUser(0.,6.5);
    lHistPtClone[0]->GetYaxis()->SetTitle("1/N_{ev} d^{2}N/(dp_{T}dy) [(GeV/c)^{-1}]");    
    lHistPtClone[0]->SetStats(0);
    for (int i = 0; i < binnumber; i++){           
        lHistPtClone[i]->Draw("SAME");    
        FFit[i]->Draw("SAME");
    }
    FFit[0]->Draw("SAME");
  
    //   
    double errpercentile[binnumber], centrpercentile[binnumber];
    double err0[binnumber];
    for (Int_t n = 0; n < binnumber; n++) {
      centrpercentile[n] = (percentile[n] + percentile[ n + 1 ]) / 2;
      errpercentile[n] = (percentile[ n + 1 ] - percentile[n]) / 2;
      err0[n] = 0.;
    }
   
    TGraphErrors *ChiGraph = new TGraphErrors(binnumber,centrpercentile,Chi,errpercentile,err0); 
  
    TH1D* Yields = new TH1D("Yields", Form(";percentile %s;#GT dN/dy #LT", fWhichEstimator.Data()), binnumber, percentile);
    for (int bin = 1; bin < Yields->GetNbinsX()+1; bin++){
      Yields->SetBinContent(bin,Yield[bin-1]);
      Yields->SetBinError(bin,YieldStat[bin-1]);
    }

    //Chi graph
    TCanvas* chi = new TCanvas("chi","",1100,500);
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
    TH1D* bkg = new TH1D("bkg", " ", 102, -10., 101.);
    bkg->SetBinContent(1,0.);
    bkg->SetBinContent(50,2.3);
    bkg->SetLineColor(kWhite);
    bkg->SetStats(0);
    bkg->GetYaxis()->SetTitle("#Chi^{2} / ndf");
    bkg->GetYaxis()->SetRangeUser(0.,10.);
    bkg->GetXaxis()->SetTitle(Form("percentile %s (%)",fWhichEstimator.Data()));
    bkg->Draw();
    ChiGraph->Draw("SAME EP");
    if (Func.Contains("Boltz")) ChiGraph->SetName(Form("Chi_%s","Boltz"));
    if (Func.Contains("BlastWave"))  ChiGraph->SetName(Form("Chi_%s","BlastWave"));
    if (Func.Contains("MTexpo")) ChiGraph->SetName(Form("Chi_%s","MTexpo"));
    if (Func.Contains("FermiDirac"))  ChiGraph->SetName(Form("Chi_%s","FermiDirac"));
    if (Func.Contains("BoseEinstein"))  ChiGraph->SetName(Form("Chi_%s","BoseEinstein"));
    if (Func.Contains("Levy"))  ChiGraph->SetName(Form("Chi_%s","Levy"));

    if (Func.Contains("Boltz")) ChiGraph->SetLineColor(kOrange);
    if (Func.Contains("BlastWave")) ChiGraph->SetLineColor(kYellow+3);
    if (Func.Contains("MTexpo"))ChiGraph ->SetLineColor(kGreen+1);
    if (Func.Contains("FermiDirac"))ChiGraph->SetLineColor(kBlue);
    if (Func.Contains("BoseEinstein")) ChiGraph->SetLineColor(kBlack);
    if (Func.Contains("Levy")) ChiGraph->SetLineColor(kRed);
    if (Func.Contains("Boltz")) ChiGraph->SetMarkerColor(kOrange);
    if (Func.Contains("BlastWave")) ChiGraph->SetMarkerColor(kYellow+3);
    if (Func.Contains("MTexpo"))ChiGraph ->SetMarkerColor(kGreen+1);
    if (Func.Contains("FermiDirac"))ChiGraph->SetMarkerColor(kBlue);
    if (Func.Contains("BoseEinstein")) ChiGraph->SetMarkerColor(kBlack);
    if (Func.Contains("Levy")) ChiGraph->SetMarkerColor(kRed);


    TFile* Write = new TFile(outputfile, "UPDATE");
    TDirectoryFile *lDirCutVar = new TDirectoryFile(Form("%s",Func.Data()),Form("Function %s",Func.Data()));
    lDirCutVar->cd();
    for (int i = 0; i < binnumber; i++){  
        FFit[i]->Write();
    }
    ChiGraph->Write();
    Yields->Write();
    Write->cd();
    Write->Close();
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
    if(trials > 10) {
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
    
    double lSigmaDelta = TMath::Sqrt( TMath::Abs( TMath::Power(Aerr,2) - TMath::Power(Berr,2) ) );
    //Computation of relationship to h2 for plotting in ratio plot 
    if ( B > 1e-12 ){ 
      lSigmaDelta /= B; 
    }else{ 
      lSigmaDelta = 0; 
    }
  return lSigmaDelta;
}
