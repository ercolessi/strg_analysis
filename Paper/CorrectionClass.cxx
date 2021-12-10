/***********************************************
  Cascade Correction Module
  ----------------------
    This version: 2021
    --- Francesca Ercolessi
    francesca.ercolessi@cern.ch
 ***********************************************/

//--- For C++ ----
#include <TApplication.h>
#include <TMatrix.h>
#include <TMatrixD.h>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLine.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <THashList.h>
#include <TDirectoryFile.h>
#include <cstdlib>
using namespace std;

//--- For ROOT ---
#include "AliVTrack.h"
#include "CorrectionClass.h"

CorrectionClass::CorrectionClass()
{
    fWhichParticle = "Xi"; //Default
    fOutputDataFile = "";
    fMCFilepTshape = "";
    fMCFileNormCorr = "";
    fWhichVarEstimator = "";
    fWhichFixedEstimator = "";
    fWhichEnergyEstimator = "";
    fWhichMultEstimator = "";
    fLowFixed = -1.;
    fHighFixed = -1.;
    fDoMCShape = kFALSE;
    fDoEventLoss = kFALSE;
    fDoSignalLoss = kFALSE;
}

CorrectionClass::CorrectionClass(TString fParticleType)
{
    // Allows definition of Particle Type in analysis.
    // Possible Options are "Xi", "Omega", "Lambda".
    // If some other string is given, this constructor will
    // default to "Xi".
    fWhichParticle = fParticleType;
    fOutputDataFile  = "";
    fMCFilepTshape = "";
    fMCFileNormCorr = "";
    fWhichVarEstimator = "";
    fWhichFixedEstimator = "";
    fWhichEnergyEstimator = "";
    fWhichMultEstimator = "";
    fLowFixed = -1.;
    fHighFixed = -1.;
    fDoMCShape = kFALSE;
    fDoEventLoss = kFALSE;
    fDoSignalLoss = kFALSE;
}

/***********************************************
  --- Setters For Configuration ---
************************************************/

// Filename Setters
void CorrectionClass::SetOutputDataFile  ( TString DataFilename   ){
    fOutputDataFile = DataFilename;
}
void CorrectionClass::SetMCFilepTshape  ( TString MCFilename   ){
    fMCFilepTshape = MCFilename;
}
void CorrectionClass::SetMCFileNormCorr  ( TString MCNormCorr   ){
    fMCFileNormCorr = MCNormCorr;
}

// Bin Limit Setters
void CorrectionClass::SetPtBinLimits(Long_t got_fptbinnumb,const Double_t *got_fptbinlimits){
    //Function to set pt binning. First argument is the number of pt bins, second is
    //an array with bin limits.
    fptbinnumb = got_fptbinnumb;
    for(int ix = 0;ix<fptbinnumb+1;ix++){
        fptbinlimits[ix] = got_fptbinlimits[ix];
    }
    for(int ix = 0;ix<fptbinnumb;ix++){
        fptX[ix] = (fptbinlimits[ix+1] + fptbinlimits[ix])/2.;
    }
    cout<<"[CorrectionClass] Received "<<fptbinnumb<<" pt bins, set accordingly."<<endl;
}
//
void CorrectionClass::SetPercBinLimits(Long_t got_fpercbinnumb, const Double_t *got_fpercbinlimits){
    //Function to set pt binning. First argument is the number of pt bins, second is
    //an array with bin limits.
    fpercbinnumb = got_fpercbinnumb;
    for(int ix = 0;ix<fpercbinnumb+1;ix++){
        fpercbinlimits[ix] = got_fpercbinlimits[ix];
    }
    for(int ix = 0;ix<fpercbinnumb;ix++){
        fpercX[ix] = (fpercbinlimits[ix+1] + fpercbinlimits[ix])/2.;
    }
    cout<<"[CorrectionClass] Received "<<fpercbinnumb<<" perc bins, set accordingly."<<endl;
}

//Set Fixed Bin Limits
void CorrectionClass::SetFixedPercLimits(Double_t got_Low, Double_t got_High){
    //Sets the fixed selection variables limits
    fLowFixed = got_Low;
    fHighFixed = got_High;
}

//Set Estimators
void CorrectionClass::SetVarEstimator(TString got_var){
    fWhichVarEstimator = got_var;
}
void CorrectionClass::SetFixedEstimator(TString got_fixed){
    fWhichFixedEstimator = got_fixed;
}
void CorrectionClass::SetEnergyEstimator(TString got_energyvar){
    fWhichEnergyEstimator = got_energyvar;
}
void CorrectionClass::SetMultEstimator(TString got_multvar){
    fWhichMultEstimator = got_multvar;
}

//Do corrections
void CorrectionClass::DoMCShape(Bool_t got_DoMCShape){
    fDoMCShape = got_DoMCShape;
}
//
void CorrectionClass::DoEventLoss(Bool_t got_DoEventLoss){
    fDoEventLoss = got_DoEventLoss;
}
//
void CorrectionClass::DoSignalLoss(Bool_t got_DoSignalLoss){
    fDoSignalLoss = got_DoSignalLoss;
}

//Levy-Tsallis function
Double_t CorrectionClass::LevyTsallis_Func(const Double_t *x, const Double_t *p)
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

TF1 * CorrectionClass::LevyTsallis(const Char_t *name, Double_t mass, Double_t par2min , Double_t par2max ,Double_t par3min, Double_t par3max, Double_t n, Double_t C, Double_t norm)
{ 
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 25., 4.);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e+3);
  fLevyTsallis->SetParLimits(2, par2min, par2max);
  fLevyTsallis->SetParLimits(3, par3min, par3max);
  return fLevyTsallis;
}

Double_t CorrectionClass::ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ){
    //Error in a Ratio
    if(B!=0){
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( errorfromtop + errorfrombottom );
    }
    return 1;
}

//----------------------------------------------------------------------------------------------------
void CorrectionClass::ErrorInRatioCorr ( TH1F* h1, TH1F *h2 ){ 
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
    if ( h2->GetBinContent(i) != 0 ){ 
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


//Printer
void CorrectionClass::PrintConfiguration() {
    //Print current analysis configuration
    Double_t lParticleMass = 1.32171; //Default
    if(fWhichParticle == "Omega") lParticleMass = 1.67245;
    if(fWhichParticle == "Lambda") lParticleMass = 1.115683;

    //
    cout << "\n\n";
    cout<<"--------------- Configuration --------------------------"<<endl;
    cout<<" Analysed Particle.............: "<<fWhichParticle<<endl;
    cout<<" This Particle Mass............: "<<lParticleMass<<endl;
    cout<<" Variable estimator............: "<<fWhichVarEstimator<<endl;
    cout<<" Fixed estimator...............: "<<fWhichFixedEstimator<<" in ["<< fLowFixed << " - " << fHighFixed << "] " << endl;
    cout << "\n";
    cout<<"------------- Corrections to apply----------------------"<<endl;
    cout<<" MC pT shape ...................: "<<(fDoMCShape ? "YES" : "NO")<<endl;
    cout<<" Event Loss ....................: "<<(fDoEventLoss ? "YES" : "NO")<<endl;
    cout<<" Signal Loss ...................: "<<(fDoSignalLoss ? "YES" : "NO")<<endl;
    cout << "\n";

}

void CorrectionClass::DoAnalysis(){

    Double_t lParticleMass = 1.32171; //Default
    if(fWhichParticle == "Omega") lParticleMass = 1.67245;
    if(fWhichParticle == "Lambda") lParticleMass = 1.115683;


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
      fAntiParticle = "AntiLambda";
    }  
    //
    TString folder = "../FullStatistics/RawSpectra/18i";

    TFile* filepart[fpercbinnumb], *fileantipart[fpercbinnumb];
    TH1F* hpart[fpercbinnumb], *hantipart[fpercbinnumb], *hspectra[fpercbinnumb];

    for(int nmult = 0; nmult < fpercbinnumb; nmult++) { //loop over perc selections

        double multmin, multmax, eemin, eemax;
        //
        if (fWhichVarEstimator.Contains("SPD")){
            multmin = fpercbinlimits[nmult];
            multmax = fpercbinlimits[nmult+1];
            eemin = fLowFixed; 
            eemax = fHighFixed;
        }
        //
        if (fWhichVarEstimator.Contains("V0M") || fWhichVarEstimator.Contains("ZDC")){
            eemin = fpercbinlimits[nmult];
            eemax = fpercbinlimits[nmult+1];
            multmin = fLowFixed; 
            multmax = fHighFixed;
        }

        //Get spectra from files
        //particle
        filepart[nmult] = new TFile(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root", folder.Data(), fParticle.Data(), fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        hpart[nmult] = (TH1F *) filepart[nmult]->Get(Form("fHistPt%s", fParticle.Data()));
        //anti-particle
        fileantipart[nmult] = new TFile(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root", folder.Data(), fAntiParticle.Data(), fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        hantipart[nmult] = (TH1F *) fileantipart[nmult]->Get(Form("fHistPt%s", fAntiParticle.Data()));
     
        //Sum spectrum
        hspectra[nmult] = (TH1F*)hpart[nmult]->Clone(Form("PtSpectrumPercBin%i",nmult));
        hspectra[nmult]->SetTitle(Form("%s_%03.0f_%03.0f - %s_%03.0f_%03.0f", fWhichVarEstimator.Data(), fpercbinlimits[nmult], fpercbinlimits[nmult+1], fWhichFixedEstimator.Data(), fLowFixed, fHighFixed));
        hspectra[nmult]->Reset();
        for (int bin = 1 ; bin <=  hpart[nmult]->GetNbinsX(); bin ++ ){ // loop over bins
            hspectra[nmult]->SetBinContent(bin, hpart[nmult]->GetBinContent(bin) + hantipart[nmult]->GetBinContent(bin));
            hspectra[nmult]->SetBinError(bin, TMath::Sqrt(hpart[nmult]->GetBinError(bin)*hpart[nmult]->GetBinError(bin)+hantipart[nmult]->GetBinError(bin)*hantipart[nmult]->GetBinError(bin)));
        } // end loop over bins
    } // end loop over perc selections

    //===========================================================================================
    // ============================= pT shape MC correction =====================================
    //===========================================================================================

    TH1D* HistRecoPart = 0x0;
    TH1D* HistRecoAntiPart = 0x0;
    TH1D* HistReco = 0x0;
    TH1F* HistRecoMultIT1[fpercbinnumb],* HistRecoMultIT2[fpercbinnumb],* HistRecoMultIT3[fpercbinnumb];
    TH1F* HistGenMultIT1[fpercbinnumb],* HistGenMultIT2[fpercbinnumb],* HistGenMultIT3[fpercbinnumb];
    Double_t par_n[fpercbinnumb];
    Double_t par_C[fpercbinnumb];
    Double_t par_norm[fpercbinnumb];
    TF1* Levy[fpercbinnumb], *LevyIT1[fpercbinnumb], *LevyIT2[fpercbinnumb], *LevyIT3[fpercbinnumb];
    TF1* LevyMCIT1[fpercbinnumb], *LevyMCIT2[fpercbinnumb], *LevyMCIT3[fpercbinnumb];
    TF1* LevyMC = LevyTsallis(Form("LevyMCIT0FitMB"), lParticleMass);
    double minfit = 0.6, maxfit = 6.5;
    TH1F *hPtSpectrumGen[fpercbinnumb];
    TH1F *hPtSpectrumGenMCMB = new TH1F("hPtSpectrumGenMCIT0MB", "Cascade generated from Levy-Tsallis Fit;#it{p}_{T} (GeV/#it{c});Counts", 250, 0.,25);
    hPtSpectrumGenMCMB->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), 0., 100., fWhichEnergyEstimator.Data(), 0., 100.));
    TH1F *hPtSpectrumGenMCIT1[fpercbinnumb], *hPtSpectrumGenMCIT2[fpercbinnumb], *hPtSpectrumGenMCIT3[fpercbinnumb];
    TFile* fileMC;
    TList* clistMC;
    TH3D* h3DGeneratedPart, *h3DGeneratedAntiPart;
    Int_t rapBin_min;
    Int_t rapBin_max;
    Int_t multBin_min;
    Int_t multBin_max;    
    TH1D* HistGenPart = 0x0, *HistGenAntiPart = 0x0, *HistGen = 0x0, *HistGenClone = 0x0;
    TH1F* hRatioPtShape[fpercbinnumb], * hRatioPtShapeIT1[fpercbinnumb], * hRatioPtShapeIT2[fpercbinnumb], * hRatioPtShapeIT3[fpercbinnumb];
    TH1F* HistRecoMult[fpercbinnumb], *HistGenMult[fpercbinnumb], *HistRecoMult_rebin[fpercbinnumb], *HistGenMult_rebin[fpercbinnumb], *CorrEfficiency[fpercbinnumb], *CorrEfficiencyIT1[fpercbinnumb], *CorrEfficiencyIT2[fpercbinnumb], *CorrEfficiencyIT3[fpercbinnumb];
    TH1F* RatioEff[fpercbinnumb],* RatioEffIT1[fpercbinnumb],* RatioEffIT2[fpercbinnumb],* RatioEffIT3[fpercbinnumb]; 
    TH1F* HistReco_rebin;
    TH1F* HistGen_rebin;  
    TH1F* EfficiencyMB;  
    TH1F* hspectraIT1[fpercbinnumb],* hspectraIT2[fpercbinnumb],* hspectraIT3[fpercbinnumb];
    TH1F *hPtSpectrumGenIT1[fpercbinnumb], *hPtSpectrumGenIT2[fpercbinnumb], *hPtSpectrumGenIT3[fpercbinnumb];
    TH1F *hPtSpectrumGenMCMBIT1, *hPtSpectrumGenMCMBIT2, *hPtSpectrumGenMCMBIT3;
        
    if (fDoMCShape){

        cout << "\n\n";
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"-------------------- Doing MC pT shape correction --------------------"<<endl;
        cout<<"----------------------------------------------------------------------"<<endl;
        cout << "\n";
        cout<<"--------------------------- Iteration 0 ------------------------------"<<endl;
        

        //Get Generated MC Histo (input MC pT-shape)
        fileMC = new TFile(fMCFilepTshape,"READ");
        clistMC = (TList*)fileMC->Get("PWGLF_StrVsMult_MC/cList");
        h3DGeneratedPart = (TH3D*)clistMC->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", fParticle.Data()));
        h3DGeneratedAntiPart = (TH3D*)clistMC->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", fAntiParticle.Data()));
        // find projection bins in rapidity
        rapBin_min = h3DGeneratedPart->GetYaxis()->FindBin( -0.5+1.e-6 );
        rapBin_max = h3DGeneratedPart->GetYaxis()->FindBin( +0.5-1.e-6 );
        // 
        // find projection bins in multiplicity
        multBin_min = h3DGeneratedPart->GetZaxis()->FindBin( 0.+1.e-6 );
        multBin_max = h3DGeneratedPart->GetZaxis()->FindBin( 100.-1.e-6 );
        
        HistGenPart = (TH1D*)h3DGeneratedPart->ProjectionX(Form("HistGen%s", fParticle.Data()),rapBin_min,rapBin_max, multBin_min, multBin_max);
        HistGenAntiPart = (TH1D*)h3DGeneratedAntiPart->ProjectionX(Form("HistGen%s", fAntiParticle.Data()),rapBin_min,rapBin_max, multBin_min, multBin_max);
        HistGen = (TH1D*)HistGenPart->Clone("HistGen");
        HistGen->Reset();
        for (int bin = 1; bin <= HistGen->GetNbinsX(); bin++){   
        HistGen->SetBinContent(bin, HistGenPart->GetBinContent(bin) + HistGenAntiPart->GetBinContent(bin));
        HistGen->SetBinError(bin, TMath::Sqrt(HistGenPart->GetBinError(bin)*HistGenPart->GetBinError(bin)+
            HistGenAntiPart ->GetBinError(bin)*HistGenAntiPart->GetBinError(bin)));
        }    
        HistGenClone = (TH1D*)HistGen->Clone("hMCPtshape0");
        HistGenClone->Reset();
        
        //Rebin
        for (int ibin = 1; ibin <= HistGen->GetNbinsX(); ibin++){
            HistGenClone->SetBinContent(ibin, HistGen->GetBinContent(ibin)/HistGen->GetBinWidth(ibin));
            HistGenClone->SetBinError(ibin, HistGen->GetBinError(ibin)/HistGen->GetBinWidth(ibin));
        }
        double integral = HistGenClone->Integral(0,-1);
        HistGenClone->Scale(1./integral);
        // Fit this wight away
        HistGenClone->Fit(LevyMC,"0q","",0.,10.);

        //Get Reco MC histo
        HistRecoPart = (TH1D*)filepart[0]->Get("fHistReco");
        HistRecoAntiPart = (TH1D*)fileantipart[0]->Get("fHistReco");
        HistReco = (TH1D*) HistRecoPart->Clone("HistReco");
        HistReco->Reset();
        for (int bin = 1 ; bin <=  HistRecoPart->GetNbinsX(); bin ++ ){
            HistReco->SetBinContent(bin, HistRecoPart->GetBinContent(bin) + HistRecoAntiPart->GetBinContent(bin));
            HistReco->SetBinError(bin, TMath::Sqrt(HistRecoPart->GetBinError(bin)*HistRecoPart->GetBinError(bin)+HistRecoAntiPart->GetBinError(bin)*HistRecoAntiPart->GetBinError(bin)));
        }

        //Rebin Gen e Reco per fare efficienze default
        HistReco_rebin   = new TH1F("HistReco_rebin",";p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);
        HistGen_rebin  = new TH1F("HistGen_rebin",";p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);  
        //Let's define a total efficiency, it will be usefull aftewards     
        EfficiencyMB  = new TH1F("EfficiencyMB","Eff #Xi+#bar{#Xi};p_{T} (GeV/c);Eff", fptbinnumb, fptbinlimits);  

        Double_t tempptreco;
        for(long i = 1; i<HistReco->GetNbinsX()+1;i++){
            tempptreco = HistReco->GetXaxis()->GetBinCenter(i);
            for(long filling = 0; filling<HistReco->GetBinContent(i); filling++){
                HistReco_rebin->Fill(tempptreco);
            }
        }
        Double_t tempptgen;
        for(long i = 1; i<HistGen->GetNbinsX()+1;i++){
            tempptgen = HistGen->GetXaxis()->GetBinCenter(i);
            for(long filling = 0; filling<HistGen->GetBinContent(i); filling++){
                HistGen_rebin->Fill(tempptgen);
            }
        }   
    
        for (int bin = 1; bin <= EfficiencyMB->GetNbinsX(); bin ++){
            if (HistGen_rebin->GetBinContent(bin) != 0){
                EfficiencyMB->SetBinContent(bin , HistReco_rebin->GetBinContent(bin)/HistGen_rebin->GetBinContent(bin));
                EfficiencyMB->SetBinError(bin, ErrorInRatio(HistReco_rebin->GetBinContent(bin),HistReco_rebin->GetBinError(bin),
                    HistGen_rebin->GetBinContent(bin),HistGen_rebin->GetBinError(bin)));
            }
        }

        for(int nmult = 0; nmult < fpercbinnumb; nmult++){ //loop over perc selections

            double multmin, multmax, eemin, eemax;
            //
            if (fWhichVarEstimator.Contains("SPD")){
                multmin = fpercbinlimits[nmult];
                multmax = fpercbinlimits[nmult+1];
                eemin = fLowFixed; 
                eemax = fHighFixed;
            }
            //
            if (fWhichVarEstimator.Contains("V0M") || fWhichVarEstimator.Contains("ZDC")){
                eemin = fpercbinlimits[nmult];
                eemax = fpercbinlimits[nmult+1];
                multmin = fLowFixed; 
                multmax = fHighFixed;
            }

            //Fit function 
            Levy[nmult] = LevyTsallis(Form("LevyIT0FitPercBin%i",nmult),lParticleMass);
            Levy[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

            Int_t fitres;
            Int_t trials = 0;
            trials = 0;
            do {
                fitres = hspectra[nmult]->Fit(Levy[nmult], "0q","",minfit,maxfit);
                Printf("Fit trial: %d", trials++);
                if(trials > 10) {
                    Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
                    break;
                }
            }
            while (fitres != 0);
            par_n[nmult] = Levy[nmult]->GetParameter(1);
            par_C[nmult] = Levy[nmult]->GetParameter(2);
            par_norm[nmult] = Levy[nmult]->GetParameter(3);

            //Define generated histos
            hPtSpectrumGen[nmult] = new TH1F(Form("hPtSpectrumGenIT0PercBin%i",nmult), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", 250, 0, 25);
            hPtSpectrumGen[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        
            //Define ratio shape
            hRatioPtShape[nmult] = (TH1F*)hPtSpectrumGen[nmult]->Clone(Form("hRatioPtShapeIT0PercBin%i",nmult));
            hRatioPtShape[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            hRatioPtShape[nmult]->Reset();

        }// end loop over perc selections

        //++++++++++++++++++++++++++++++++++++++++
        // 0) Generate histos from Fit Levy-Tsallis
        //++++++++++++++++++++++++++++++++++++++++
        //The pT-spectra for both data and Monte Carlo were fitted by a Levy-Tsallis function. 
        //Here pT-spectra are “generated” according to the fitted shapes using an high granularity. 

        const int n = 50000000; 
        
        //MB spectra
        for (Int_t i = 0; i < n; i++) {
            double xMC = (double)LevyMC->GetRandom(0., 25.);
            hPtSpectrumGenMCMB->Fill(xMC);
        }
        // /dpT
        for (Int_t i = 0; i <= hPtSpectrumGenMCMB->GetNbinsX(); i++) {
            double entryMC = hPtSpectrumGenMCMB->GetBinContent(i + 1) /
                            hPtSpectrumGenMCMB->GetBinWidth(i + 1);
            hPtSpectrumGenMCMB->SetBinContent(i + 1, entryMC);
        }

        // multiplicity spectra RD
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
            for (Int_t i = 0; i < n; i++) {
                double x = (double)Levy[nmult]->GetRandom(0., 25.);
                hPtSpectrumGen[nmult]->Fill(x);
                }
                // /dpT
                for (Int_t i = 0; i <hPtSpectrumGen[nmult]->GetNbinsX() ; i++) {
                double entry = hPtSpectrumGen[nmult]->GetBinContent(i + 1) /
                                hPtSpectrumGen[nmult]->GetBinWidth(i + 1);
                hPtSpectrumGen[nmult]->SetBinContent(i + 1, entry);
                hPtSpectrumGen[nmult]->SetMarkerStyle(8);
                hPtSpectrumGen[nmult]->SetMarkerSize(0.9);
            }
        }

        //++++++++++++++++++++++++++++++++++++++++
        // 1) Iterative process for correction
        //++++++++++++++++++++++++++++++++++++++++ 
        //The ratios of the “fitted” shapes obtained from the data over the Monte Carlo input pT-shape
        //are computed with the high granularity

        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
           
            for (int bin = 1; bin <= hRatioPtShape[0]->GetNbinsX(); bin ++){
                if (hPtSpectrumGenMCMB->GetBinContent(bin) != 0 ){      
                    hRatioPtShape[nmult]->SetBinContent(bin, hPtSpectrumGen[nmult]->GetBinContent(bin)/hPtSpectrumGenMCMB->GetBinContent(bin));
                    hRatioPtShape[nmult]->SetBinError(bin, ErrorInRatio(hPtSpectrumGen[nmult]->GetBinContent(bin),hPtSpectrumGen[nmult]->GetBinError(bin),
                                    hPtSpectrumGenMCMB->GetBinContent(bin),hPtSpectrumGenMCMB->GetBinError(bin)));
                }
                else hRatioPtShape[nmult]->SetBinContent(bin,1); //x*1=x
            }
        }

        //++++++++++++++++++++++++++++++++++++++++
        // 2) Re-weight reco and gen spectra 
        //++++++++++++++++++++++++++++++++++++++++ 
        //The ratios at the previous point were used to “re-weight” the reconstructed and generated pT spectra 
        //used to compute efficiencies (same high granularity)

        for (int nmult = 0; nmult < fpercbinnumb; nmult++) { //loop pver perc selections
            
            double multmin, multmax, eemin, eemax;
            //
            if (fWhichVarEstimator.Contains("SPD")){
                multmin = fpercbinlimits[nmult];
                multmax = fpercbinlimits[nmult+1];
                eemin = fLowFixed; 
                eemax = fHighFixed;
            }
            //
            if (fWhichVarEstimator.Contains("V0M") || fWhichVarEstimator.Contains("ZDC")){
                eemin = fpercbinlimits[nmult];
                eemax = fpercbinlimits[nmult+1];
                multmin = fLowFixed; 
                multmax = fHighFixed;
            }

            HistRecoMult_rebin[nmult]   = new TH1F(Form("HistRecoMult_rebin%i",nmult),"Cascade MC count;p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);
            HistGenMult_rebin[nmult]   = new TH1F(Form("HistGenMult_rebin%i",nmult),"Cascade MC count;p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);

            HistRecoMult[nmult] = (TH1F*) HistReco->Clone(Form("HistRecoPercBin%i",nmult));
            HistRecoMult[nmult] ->Reset();
            HistGenMult[nmult] = (TH1F*) HistGen->Clone(Form("HistGenPercBin%i",nmult));
            HistGenMult[nmult] ->Reset();
            for (int bin = 1; bin <= HistRecoMult[nmult]->GetNbinsX(); bin ++){
                //
                HistRecoMult[nmult]-> SetBinContent(bin,HistReco->GetBinContent(bin)*hRatioPtShape[nmult]->GetBinContent(bin));
                HistRecoMult[nmult]-> SetBinError(bin,ErrorInRatio(HistReco->GetBinContent(bin), HistReco->GetBinError(bin),
                    hRatioPtShape[nmult]->GetBinContent(bin), hRatioPtShape[nmult]->GetBinError(bin))
                    );
                //
                HistGenMult[nmult]-> SetBinContent(bin, HistGen->GetBinContent(bin)*hRatioPtShape[nmult]->GetBinContent(bin));
                HistGenMult[nmult]-> SetBinError(bin, ErrorInRatio(HistGen->GetBinContent(bin), HistGen->GetBinError(bin),
                    hRatioPtShape[nmult]->GetBinContent(bin), hRatioPtShape[nmult]->GetBinError(bin)
                    ));
            }

            Double_t tempptreco;
            for(long i = 1; i<HistRecoMult[nmult]->GetNbinsX()+1;i++){
                tempptreco = HistRecoMult[nmult]->GetXaxis()->GetBinCenter(i);
                for(long filling = 0; filling<HistRecoMult[nmult]->GetBinContent(i); filling++){
                    HistRecoMult_rebin[nmult]->Fill(tempptreco);
                }
            }
            Double_t tempptgen;
            for(long i = 1; i<HistGenMult[nmult]->GetNbinsX()+1;i++){
                tempptgen = HistGenMult[nmult]->GetXaxis()->GetBinCenter(i);
                for(long filling = 0; filling<HistGenMult[nmult]->GetBinContent(i); filling++){
                    HistGenMult_rebin[nmult]->Fill(tempptgen);
                }
            }
        
            //++++++++++++++++++++++++++++++++++++++++
            // 3) Correct Efficiency
            //++++++++++++++++++++++++++++++++++++++++  
            //New efficiencies are recomputed by “re-binning” the pT-spectra at the point 2) in the same pT
            //bins used for the analysis. In this way the “corrected” efficiencies are obtained.

        
            CorrEfficiency[nmult] = (TH1F*)HistRecoMult_rebin[nmult]->Clone(Form("hCorrEfficiencyIT0PercBin%i",nmult));
            CorrEfficiency[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            CorrEfficiency[nmult]->Reset();
            for (int bin = 1; bin <= HistRecoMult_rebin[nmult]->GetNbinsX(); bin ++){
                if (HistGenMult_rebin[nmult]->GetBinContent(bin) != 0){
                    CorrEfficiency[nmult]->SetBinContent(bin , HistRecoMult_rebin[nmult]->GetBinContent(bin)/HistGenMult_rebin[nmult]->GetBinContent(bin));
                    CorrEfficiency[nmult]->SetBinError(bin, ErrorInRatio(HistRecoMult_rebin[nmult]->GetBinContent(bin),HistRecoMult_rebin[nmult]->GetBinError(bin),
                    HistGenMult_rebin[nmult]->GetBinContent(bin),HistGenMult_rebin[nmult]->GetBinError(bin)));
                }
            }

            RatioEff[nmult] = (TH1F*)CorrEfficiency[nmult]->Clone(Form("CorrFactorsIT0PercBin%i",nmult));
            RatioEff[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        
        } //end loop over perc selections

        //++++++++++++++++++++++++++++++++++++++++
        // 4) Apply Correction Factors
        //++++++++++++++++++++++++++++++++++++++++  
        //The correction factors, which are the ratios of the new efficiencies w.r.t. the older ones, were594
        //finally applied to the measured spectra.
        //Iterate 3 times

        //First iteration works with default efficiency  
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

            ErrorInRatioCorr(RatioEff[nmult],EfficiencyMB);
            hspectraIT1[nmult] = (TH1F*)hspectra[nmult]->Clone();
            hspectraIT1[nmult]->SetName(Form("PtSpectrumIT1PercBin%i",nmult));
            //
            for (int bin = 1 ; bin <=  hspectraIT1[nmult]->GetNbinsX(); bin ++ ){
                //
                double content = hspectraIT1[nmult]->GetBinContent(bin)/RatioEff[nmult]->GetBinContent(bin);
                double error = ErrorInRatio(hspectraIT1[nmult]->GetBinContent(bin), hspectraIT1[nmult]->GetBinError(bin), RatioEff[nmult]->GetBinContent(bin), RatioEff[nmult]->GetBinError(bin));
                    hspectraIT1[nmult]->SetBinContent(bin, content);
                    hspectraIT1[nmult]->SetBinError(bin,error);
            }

            //Reset parameters to avoid problems
            par_n[nmult] = 0.;
            par_C[nmult] = 0.;
            par_norm[nmult] = 0.;
        }     
        
        cout<<"------------------------------ Done! ----------------------------------"<<endl;
        cout << "\n";
        cout<<"--------------------------- Iteration 1 ------------------------------"<<endl;

        //Fit spectra obtained in previous iteration with a Levy Tsallis

        for(int nmult = 0; nmult < fpercbinnumb; nmult++){ //loop over perc selections

            double multmin, multmax, eemin, eemax;
            //
            if (fWhichVarEstimator.Contains("SPD")){
                multmin = fpercbinlimits[nmult];
                multmax = fpercbinlimits[nmult+1];
                eemin = fLowFixed; 
                eemax = fHighFixed;
            }
            //
            if (fWhichVarEstimator.Contains("V0M") || fWhichVarEstimator.Contains("ZDC")){
                eemin = fpercbinlimits[nmult];
                eemax = fpercbinlimits[nmult+1];
                multmin = fLowFixed; 
                multmax = fHighFixed;
            }

            //Fit function 
            LevyIT1[nmult] = LevyTsallis(Form("LevyIT1FitPercBin%i",nmult),lParticleMass);
            LevyIT1[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            LevyMCIT1[nmult] = LevyTsallis(Form("LevyMCIT1FitPercBin%i",nmult),lParticleMass);
            LevyMCIT1[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

            //Data----
            Int_t fitres;
            Int_t trials = 0;
            trials = 0;
            do {
                fitres = hspectraIT1[nmult]->Fit(LevyIT1[nmult], "0q","",minfit,maxfit);
                Printf("Fit trial: %d", trials++);
                if(trials > 10) {
                    Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
                    break;
                }
            }
            while (fitres != 0);
            par_n[nmult] = LevyIT1[nmult]->GetParameter(1);
            par_C[nmult] = LevyIT1[nmult]->GetParameter(2);
            par_norm[nmult] = LevyIT1[nmult]->GetParameter(3);

            //MC shape ---
            //
            Int_t fitres2;
            Int_t trials2 = 0;
            trials2 = 0;
            do {
                fitres =  HistGenMult[nmult]->Fit(LevyMCIT1[nmult], "0q","",minfit,maxfit);
                Printf("Fit trial: %d", trials++);
                if(trials > 10) {
                    Printf("FIT DOES NOT CONVERGE MC IN LINE %d",__LINE__);
                    break;
                }
            }
            while (fitres != 0);

            //Define generated histos
            hPtSpectrumGenIT1[nmult] = new TH1F(Form("hPtSpectrumGenIT1PercBin%i",nmult), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", 250, 0, 25);
            hPtSpectrumGenIT1[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            hPtSpectrumGenMCIT1[nmult] = new TH1F(Form("hPtSpectrumGenMCIT1PercBin%i",nmult), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", 250, 0, 25);
            hPtSpectrumGenMCIT1[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            
            //Define ratio shape
            hRatioPtShapeIT1[nmult] = (TH1F*)hPtSpectrumGenIT1[nmult]->Clone(Form("hRatioPtShapeIT1PercBin%i",nmult));
            hRatioPtShapeIT1[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            hRatioPtShapeIT1[nmult]->Reset();

        }// end loop over perc selections

        //++++++++++++++++++++++++++++++++++++++++++++++
        // 0_IT1) Generate histos from Fit Levy-Tsallis
        //++++++++++++++++++++++++++++++++++++++++++++++

        // multiplicity spectra 
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

            //MC spectra
            for (Int_t i = 0; i < n; i++) {
                double xMC = (double)LevyMCIT1[nmult]->GetRandom(0., 25.);
                hPtSpectrumGenMCIT1[nmult]->Fill(xMC);
            }
            // /dpT
            for (Int_t i = 0; i <= hPtSpectrumGenMCIT1[nmult]->GetNbinsX(); i++) {
                double entryMC = hPtSpectrumGenMCIT1[nmult]->GetBinContent(i + 1) /
                                hPtSpectrumGenMCIT1[nmult]->GetBinWidth(i + 1);
                hPtSpectrumGenMCIT1[nmult]->SetBinContent(i + 1, entryMC);
            }

            //Data
            for (Int_t i = 0; i < n; i++) {
                double x = (double)LevyIT1[nmult]->GetRandom(0., 25.);
                hPtSpectrumGenIT1[nmult]->Fill(x);
                }
                // /dpT
                for (Int_t i = 0; i <hPtSpectrumGenIT1[nmult]->GetNbinsX() ; i++) {
                double entry = hPtSpectrumGenIT1[nmult]->GetBinContent(i + 1) /
                                hPtSpectrumGenIT1[nmult]->GetBinWidth(i + 1);
                hPtSpectrumGenIT1[nmult]->SetBinContent(i + 1, entry);
                hPtSpectrumGenIT1[nmult]->SetMarkerStyle(8);
                hPtSpectrumGenIT1[nmult]->SetMarkerSize(0.9);
            }
        }

        //++++++++++++++++++++++++++++++++++++++++
        // 1_IT1) Iterative process for correction
        //++++++++++++++++++++++++++++++++++++++++ 
        
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
           
            for (int bin = 1; bin <= hRatioPtShapeIT1[0]->GetNbinsX(); bin ++){
                if (hPtSpectrumGenMCIT1[nmult]->GetBinContent(bin) != 0 ){      
                    hRatioPtShapeIT1[nmult]->SetBinContent(bin, hPtSpectrumGenIT1[nmult]->GetBinContent(bin)/hPtSpectrumGenMCIT1[nmult]->GetBinContent(bin));
                    hRatioPtShapeIT1[nmult]->SetBinError(bin, ErrorInRatio(hPtSpectrumGenIT1[nmult]->GetBinContent(bin),hPtSpectrumGenIT1[nmult]->GetBinError(bin),
                                    hPtSpectrumGenMCIT1[nmult]->GetBinContent(bin),hPtSpectrumGenMCIT1[nmult]->GetBinError(bin)));
                }
                else hRatioPtShapeIT1[nmult]->SetBinContent(bin,1); //x*1=x
            }
        }
        
        //++++++++++++++++++++++++++++++++++++++++
        // 2_IT1) Re-weight reco and gen spectra 
        //++++++++++++++++++++++++++++++++++++++++ 
       
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) { //loop pver perc selections
            
            double multmin, multmax, eemin, eemax;
            //
            if (fWhichVarEstimator.Contains("SPD")){
                multmin = fpercbinlimits[nmult];
                multmax = fpercbinlimits[nmult+1];
                eemin = fLowFixed; 
                eemax = fHighFixed;
            }
            //
            if (fWhichVarEstimator.Contains("V0M") || fWhichVarEstimator.Contains("ZDC")){
                eemin = fpercbinlimits[nmult];
                eemax = fpercbinlimits[nmult+1];
                multmin = fLowFixed; 
                multmax = fHighFixed;
            }

            HistRecoMult_rebin[nmult]->Reset();
            HistGenMult_rebin[nmult]->Reset();

            HistRecoMultIT1[nmult] = (TH1F*) HistRecoMult[nmult]->Clone(Form("HistRecoIT1PercBin%i",nmult));
            HistRecoMultIT1[nmult] ->Reset();
            HistGenMultIT1[nmult] = (TH1F*) HistGenMult[nmult]->Clone(Form("HistGenIT1PercBin%i",nmult));
            HistGenMultIT1[nmult] ->Reset();
            for (int bin = 1; bin <= HistRecoMultIT1[nmult]->GetNbinsX(); bin ++){
                //
                HistRecoMultIT1[nmult]-> SetBinContent(bin,HistRecoMult[nmult]->GetBinContent(bin)*hRatioPtShapeIT1[nmult]->GetBinContent(bin));
                HistRecoMultIT1[nmult]-> SetBinError(bin,ErrorInRatio(HistRecoMult[nmult]->GetBinContent(bin), HistRecoMult[nmult]->GetBinError(bin),
                    hRatioPtShapeIT1[nmult]->GetBinContent(bin), hRatioPtShapeIT1[nmult]->GetBinError(bin))
                    );
                //
                HistGenMultIT1[nmult]-> SetBinContent(bin, HistGenMult[nmult]->GetBinContent(bin)*hRatioPtShapeIT1[nmult]->GetBinContent(bin));
                HistGenMultIT1[nmult]-> SetBinError(bin, ErrorInRatio(HistGenMult[nmult]->GetBinContent(bin), HistGenMult[nmult]->GetBinError(bin),
                    hRatioPtShapeIT1[nmult]->GetBinContent(bin), hRatioPtShapeIT1[nmult]->GetBinError(bin)
                    ));
            }

            Double_t tempptreco;
            for(long i = 1; i<HistRecoMultIT1[nmult]->GetNbinsX()+1;i++){
                tempptreco = HistRecoMultIT1[nmult]->GetXaxis()->GetBinCenter(i);
                for(long filling = 0; filling<HistRecoMultIT1[nmult]->GetBinContent(i); filling++){
                    HistRecoMult_rebin[nmult]->Fill(tempptreco);
                }
            }
            Double_t tempptgen;
            for(long i = 1; i<HistGenMultIT1[nmult]->GetNbinsX()+1;i++){
                tempptgen = HistGenMultIT1[nmult]->GetXaxis()->GetBinCenter(i);
                for(long filling = 0; filling<HistGenMultIT1[nmult]->GetBinContent(i); filling++){
                    HistGenMult_rebin[nmult]->Fill(tempptgen);
                }
            }

            //++++++++++++++++++++++++++++++++++++++++
            // 3) Correct Efficiency
            //++++++++++++++++++++++++++++++++++++++++  

            //
            CorrEfficiencyIT1[nmult] = (TH1F*)HistRecoMult_rebin[nmult]->Clone(Form("hCorrEfficiencyIT1PercBin%i",nmult));
            CorrEfficiencyIT1[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            CorrEfficiencyIT1[nmult]->Reset();
            for (int bin = 1; bin <= HistRecoMult_rebin[nmult]->GetNbinsX(); bin ++){
                if (HistGenMult_rebin[nmult]->GetBinContent(bin) != 0){
                    CorrEfficiencyIT1[nmult]->SetBinContent(bin , HistRecoMult_rebin[nmult]->GetBinContent(bin)/HistGenMult_rebin[nmult]->GetBinContent(bin));
                    CorrEfficiencyIT1[nmult]->SetBinError(bin, ErrorInRatio(HistRecoMult_rebin[nmult]->GetBinContent(bin),HistRecoMult_rebin[nmult]->GetBinError(bin),
                    HistGenMult_rebin[nmult]->GetBinContent(bin),HistGenMult_rebin[nmult]->GetBinError(bin)));
                }
            }

            RatioEffIT1[nmult] = (TH1F*)CorrEfficiencyIT1[nmult]->Clone(Form("CorrFactorsIT1PercBin%i",nmult));
            RatioEffIT1[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        
        }

        //++++++++++++++++++++++++++++++++++++++++
        // 4_IT1) Apply Correction Factors
        //++++++++++++++++++++++++++++++++++++++++  
 
        //First iteration works with default efficiency  
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

            ErrorInRatioCorr(RatioEffIT1[nmult],CorrEfficiency[nmult]);
            hspectraIT2[nmult] = (TH1F*)hspectraIT1[nmult]->Clone();
            hspectraIT2[nmult]->SetName(Form("PtSpectrumIT1PercBin%i",nmult));
            //
            for (int bin = 1 ; bin <=  hspectraIT2[nmult]->GetNbinsX(); bin ++ ){
                //
                double content = hspectraIT2[nmult]->GetBinContent(bin)/RatioEffIT1[nmult]->GetBinContent(bin);
                double error = ErrorInRatio(hspectraIT2[nmult]->GetBinContent(bin), hspectraIT2[nmult]->GetBinError(bin), RatioEffIT1[nmult]->GetBinContent(bin), RatioEffIT1[nmult]->GetBinError(bin));
                    hspectraIT2[nmult]->SetBinContent(bin, content);
                    hspectraIT2[nmult]->SetBinError(bin,error);
            }

            //Reset parameters to avoid problems
            par_n[nmult] = 0.;
            par_C[nmult] = 0.;
            par_norm[nmult] = 0.;
        }     
        
        cout<<"------------------------------ Done! ----------------------------------"<<endl;
        cout << "\n\n";
        cout<<"--------------------------- Iteration 2 ------------------------------"<<endl;

        //Fit spectra obtained in previous iteration with a Levy Tsallis

        for(int nmult = 0; nmult < fpercbinnumb; nmult++){ //loop over perc selections

            double multmin, multmax, eemin, eemax;
            //
            if (fWhichVarEstimator.Contains("SPD")){
                multmin = fpercbinlimits[nmult];
                multmax = fpercbinlimits[nmult+1];
                eemin = fLowFixed; 
                eemax = fHighFixed;
            }
            //
            if (fWhichVarEstimator.Contains("V0M") || fWhichVarEstimator.Contains("ZDC")){
                eemin = fpercbinlimits[nmult];
                eemax = fpercbinlimits[nmult+1];
                multmin = fLowFixed; 
                multmax = fHighFixed;
            }

            //Fit function 
            LevyIT2[nmult] = LevyTsallis(Form("LevyIT2FitPercBin%i",nmult),lParticleMass);
            LevyIT2[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            LevyMCIT2[nmult] = LevyTsallis(Form("LevyMCIT2FitPercBin%i",nmult),lParticleMass);
            LevyMCIT2[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

            //Data----
            Int_t fitres;
            Int_t trials = 0;
            trials = 0;
            do {
                fitres = hspectraIT2[nmult]->Fit(LevyIT2[nmult], "0q","",minfit,maxfit);
                Printf("Fit trial: %d", trials++);
                if(trials > 10) {
                    Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
                    break;
                }
            }
            while (fitres != 0);
            par_n[nmult] = LevyIT2[nmult]->GetParameter(1);
            par_C[nmult] = LevyIT2[nmult]->GetParameter(2);
            par_norm[nmult] = LevyIT2[nmult]->GetParameter(3);

            //MC shape ---
            //
            Int_t fitres2;
            Int_t trials2 = 0;
            trials2 = 0;
            do {
                fitres2 =  HistGenMultIT1[nmult]->Fit(LevyMCIT2[nmult], "0q","",minfit,maxfit);
                Printf("Fit trial: %d", trials2++);
                if(trials2 > 10) {
                    Printf("FIT DOES NOT CONVERGE MC IN LINE %d",__LINE__);
                    break;
                }
            }
            while (fitres2 != 0);

            //Define generated histos
            hPtSpectrumGenIT2[nmult] = new TH1F(Form("hPtSpectrumGenIT2PercBin%i",nmult), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", 250, 0, 25);
            hPtSpectrumGenIT2[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            hPtSpectrumGenMCIT2[nmult] = new TH1F(Form("hPtSpectrumGenMCIT2PercBin%i",nmult), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", 250, 0, 25);
            hPtSpectrumGenMCIT2[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            
            //Define ratio shape
            hRatioPtShapeIT2[nmult] = (TH1F*)hPtSpectrumGenIT2[nmult]->Clone(Form("hRatioPtShapeIT2PercBin%i",nmult));
            hRatioPtShapeIT2[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            hRatioPtShapeIT2[nmult]->Reset();

        }// end loop over perc selections

        //++++++++++++++++++++++++++++++++++++++++++++++
        // 0_IT2) Generate histos from Fit Levy-Tsallis
        //++++++++++++++++++++++++++++++++++++++++++++++

        // multiplicity spectra 
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

            //MC spectra
            for (Int_t i = 0; i < n; i++) {
                double xMC = (double)LevyMCIT2[nmult]->GetRandom(0., 25.);
                hPtSpectrumGenMCIT2[nmult]->Fill(xMC);
            }
            // /dpT
            for (Int_t i = 0; i <= hPtSpectrumGenMCIT2[nmult]->GetNbinsX(); i++) {
                double entryMC = hPtSpectrumGenMCIT2[nmult]->GetBinContent(i + 1) /
                                hPtSpectrumGenMCIT2[nmult]->GetBinWidth(i + 1);
                hPtSpectrumGenMCIT2[nmult]->SetBinContent(i + 1, entryMC);
            }

            //Data
            for (Int_t i = 0; i < n; i++) {
                double x = (double)LevyIT2[nmult]->GetRandom(0., 25.);
                hPtSpectrumGenIT2[nmult]->Fill(x);
                }
                // /dpT
                for (Int_t i = 0; i <hPtSpectrumGenIT2[nmult]->GetNbinsX() ; i++) {
                double entry = hPtSpectrumGenIT2[nmult]->GetBinContent(i + 1) /
                                hPtSpectrumGenIT2[nmult]->GetBinWidth(i + 1);
                hPtSpectrumGenIT2[nmult]->SetBinContent(i + 1, entry);
                hPtSpectrumGenIT2[nmult]->SetMarkerStyle(8);
                hPtSpectrumGenIT2[nmult]->SetMarkerSize(0.9);
            }
        }

        //++++++++++++++++++++++++++++++++++++++++
        // 1_IT2) Iterative process for correction
        //++++++++++++++++++++++++++++++++++++++++ 
        
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
           
            for (int bin = 1; bin <= hRatioPtShapeIT2[0]->GetNbinsX(); bin ++){
                if (hPtSpectrumGenMCIT2[nmult]->GetBinContent(bin) != 0 ){      
                    hRatioPtShapeIT2[nmult]->SetBinContent(bin, hPtSpectrumGenIT2[nmult]->GetBinContent(bin)/hPtSpectrumGenMCIT2[nmult]->GetBinContent(bin));
                    hRatioPtShapeIT2[nmult]->SetBinError(bin, ErrorInRatio(hPtSpectrumGenIT2[nmult]->GetBinContent(bin),hPtSpectrumGenIT2[nmult]->GetBinError(bin),
                                    hPtSpectrumGenMCIT2[nmult]->GetBinContent(bin),hPtSpectrumGenMCIT2[nmult]->GetBinError(bin)));
                }
                else hRatioPtShapeIT2[nmult]->SetBinContent(bin,1); //x*1=x
            }
        }
    
        //++++++++++++++++++++++++++++++++++++++++
        // 2_IT2) Re-weight reco and gen spectra 
        //++++++++++++++++++++++++++++++++++++++++ 
       
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) { //loop pver perc selections
            
            double multmin, multmax, eemin, eemax;
            //
            if (fWhichVarEstimator.Contains("SPD")){
                multmin = fpercbinlimits[nmult];
                multmax = fpercbinlimits[nmult+1];
                eemin = fLowFixed; 
                eemax = fHighFixed;
            }
            //
            if (fWhichVarEstimator.Contains("V0M") || fWhichVarEstimator.Contains("ZDC")){
                eemin = fpercbinlimits[nmult];
                eemax = fpercbinlimits[nmult+1];
                multmin = fLowFixed; 
                multmax = fHighFixed;
            }

            HistRecoMult_rebin[nmult]->Reset();
            HistGenMult_rebin[nmult]->Reset();

            HistRecoMultIT2[nmult] = (TH1F*) HistRecoMultIT1[nmult]->Clone(Form("HistRecoIT2PercBin%i",nmult));
            HistRecoMultIT2[nmult] ->Reset();
            HistGenMultIT2[nmult] = (TH1F*) HistGenMultIT1[nmult]->Clone(Form("HistGenIT2PercBin%i",nmult));
            HistGenMultIT2[nmult] ->Reset();
            for (int bin = 1; bin <= HistRecoMultIT2[nmult]->GetNbinsX(); bin ++){
                //
                HistRecoMultIT2[nmult]-> SetBinContent(bin,HistRecoMultIT1[nmult]->GetBinContent(bin)*hRatioPtShapeIT2[nmult]->GetBinContent(bin));
                HistRecoMultIT2[nmult]-> SetBinError(bin,ErrorInRatio(HistRecoMultIT1[nmult]->GetBinContent(bin), HistRecoMultIT1[nmult]->GetBinError(bin),
                    hRatioPtShapeIT2[nmult]->GetBinContent(bin), hRatioPtShapeIT2[nmult]->GetBinError(bin))
                    );
                //
                HistGenMultIT2[nmult]-> SetBinContent(bin, HistGenMultIT1[nmult]->GetBinContent(bin)*hRatioPtShapeIT2[nmult]->GetBinContent(bin));
                HistGenMultIT2[nmult]-> SetBinError(bin, ErrorInRatio(HistGenMultIT1[nmult]->GetBinContent(bin), HistGenMultIT1[nmult]->GetBinError(bin),
                    hRatioPtShapeIT2[nmult]->GetBinContent(bin), hRatioPtShapeIT2[nmult]->GetBinError(bin)
                    ));
            }

            Double_t tempptreco;
            for(long i = 1; i<HistRecoMultIT2[nmult]->GetNbinsX()+1;i++){
                tempptreco = HistRecoMultIT2[nmult]->GetXaxis()->GetBinCenter(i);
                for(long filling = 0; filling<HistRecoMultIT2[nmult]->GetBinContent(i); filling++){
                    HistRecoMult_rebin[nmult]->Fill(tempptreco);
                }
            }
            Double_t tempptgen;
            for(long i = 1; i<HistGenMultIT2[nmult]->GetNbinsX()+1;i++){
                tempptgen = HistGenMultIT2[nmult]->GetXaxis()->GetBinCenter(i);
                for(long filling = 0; filling<HistGenMultIT2[nmult]->GetBinContent(i); filling++){
                    HistGenMult_rebin[nmult]->Fill(tempptgen);
                }
            }

            //++++++++++++++++++++++++++++++++++++++++
            // 3) Correct Efficiency
            //++++++++++++++++++++++++++++++++++++++++  

            //
            CorrEfficiencyIT2[nmult] = (TH1F*)HistRecoMult_rebin[nmult]->Clone(Form("hCorrEfficiencyIT2PercBin%i",nmult));
            CorrEfficiencyIT2[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            CorrEfficiencyIT2[nmult]->Reset();
            for (int bin = 1; bin <= HistRecoMult_rebin[nmult]->GetNbinsX(); bin ++){
                if (HistGenMult_rebin[nmult]->GetBinContent(bin) != 0){
                    CorrEfficiencyIT2[nmult]->SetBinContent(bin , HistRecoMult_rebin[nmult]->GetBinContent(bin)/HistGenMult_rebin[nmult]->GetBinContent(bin));
                    CorrEfficiencyIT2[nmult]->SetBinError(bin, ErrorInRatio(HistRecoMult_rebin[nmult]->GetBinContent(bin),HistRecoMult_rebin[nmult]->GetBinError(bin),
                    HistGenMult_rebin[nmult]->GetBinContent(bin),HistGenMult_rebin[nmult]->GetBinError(bin)));
                }
            }

            RatioEffIT2[nmult] = (TH1F*)CorrEfficiencyIT2[nmult]->Clone(Form("CorrFactorsIT2PercBin%i",nmult));
            RatioEffIT2[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        
        }

       
        //++++++++++++++++++++++++++++++++++++++++
        // 4_IT2) Apply Correction Factors
        //++++++++++++++++++++++++++++++++++++++++  
 
        //First iteration works with default efficiency  
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

            ErrorInRatioCorr(RatioEffIT2[nmult],CorrEfficiencyIT1[nmult]);
            hspectraIT3[nmult] = (TH1F*)hspectraIT2[nmult]->Clone();
            hspectraIT3[nmult]->SetName(Form("PtSpectrumIT3PercBin%i",nmult));
            //
            for (int bin = 1 ; bin <=  hspectraIT3[nmult]->GetNbinsX(); bin ++ ){
                //
                double content = hspectraIT3[nmult]->GetBinContent(bin)/RatioEffIT2[nmult]->GetBinContent(bin);
                double error = ErrorInRatio(hspectraIT3[nmult]->GetBinContent(bin), hspectraIT3[nmult]->GetBinError(bin), RatioEffIT2[nmult]->GetBinContent(bin), RatioEffIT2[nmult]->GetBinError(bin));
                hspectraIT3[nmult]->SetBinContent(bin, content);
                hspectraIT3[nmult]->SetBinError(bin,error);
            }
            
            //Reset parameters to avoid problems
            par_n[nmult] = 0.;
            par_C[nmult] = 0.;
            par_norm[nmult] = 0.;
        }     
        
        cout<<"------------------------------ Done! ----------------------------------"<<endl;
        cout << "\n";
        cout<<"------------ MC pT shape correction applied succesfully! --------------"<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl;
    }

    //===========================================================================================
    //============================== END pT shape MC correction =================================
    //===========================================================================================

    //===========================================================================================
    //=============================== Event Loss correction =====================================
    //===========================================================================================

    TFile* filenorm;
    TH1F* heventloss;
    TH1F* hspectraEvtLoss[fpercbinnumb];
	
    //se fai la MC pt fagli prendere gli spettri suoi, altrimenti fallo partire da spettro 1 
    if (fDoEventLoss){

        cout << "\n\n";
        cout<<"-----------------------------------------------------------------------"<<endl;
        cout<<"-------------------- Doing Event Loss correction ---------------------"<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl;
        cout << "\n";
        cout<<" Event Loss Corr File ..........: "<< fMCFileNormCorr.Data() <<endl;
        cout << "\n";
        cout<<" !!!!!!!!!!! Please make sure event loss percentile bins match your analisis!!! No check here !!!!!!!!!!!!!!!"<<endl;
        cout << "\n";
        
        filenorm = TFile::Open(fMCFileNormCorr);
        //
        heventloss = (TH1F*)filenorm->Get("EventLoss/hevtloss");

        if (fDoMCShape){
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
                hspectraEvtLoss[nmult] = (TH1F*) hspectraIT3[nmult]->Clone(Form("PtSpectrumEvtLossCorrPercBin%i",nmult));
            }
        } else {
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
                hspectraEvtLoss[nmult] = (TH1F*) hspectra[nmult]->Clone(Form("PtSpectrumEvtLossCorrPercBin%i",nmult));
            }
        }

        //Apply correction
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
            for (int bin = 1 ; bin <=  hspectraEvtLoss[nmult]->GetNbinsX(); bin ++ ){
                //
                double content = hspectraEvtLoss[nmult]->GetBinContent(bin)*heventloss->GetBinContent(nmult+1);
                double error = ErrorInRatio(hspectraEvtLoss[nmult]->GetBinContent(bin), hspectraEvtLoss[nmult]->GetBinError(bin), heventloss->GetBinContent(nmult+1), heventloss->GetBinError(nmult+1));
                hspectraEvtLoss[nmult]->SetBinContent(bin, content);
                hspectraEvtLoss[nmult]->SetBinError(bin,error);
            }    
        }    
    }
    cout<<"------------------------------ Done! ----------------------------------"<<endl;
    cout << "\n";
    cout<<"------------- Event Loss correction applied succesfully! --------------"<<endl;
    cout<<"-----------------------------------------------------------------------"<<endl;

    //===========================================================================================
    // ============================== END Event Loss correction =================================
    //===========================================================================================

    
    //===========================================================================================
    //=============================== Signal Loss correction ====================================
    //===========================================================================================

    TH1D* hsgnloss[fpercbinnumb]; //these are TH2D!! As from definition in the AliPhyiscs task
    TH1F* hspectraSgnLoss[fpercbinnumb];

    //se fai la MC pt fagli prendere gli spettri suoi, altrimenti fallo partire da spettro 1 
    if (fDoSignalLoss){

        cout << "\n\n";
        cout<<"-----------------------------------------------------------------------"<<endl;
        cout<<"-------------------- Doing Signal Loss correction ---------------------"<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl;
        cout << "\n";
        cout<<" Signal Loss Corr File ..........: "<< fMCFileNormCorr.Data() <<endl;
        cout << "\n\n";
        
        if (!fDoEventLoss) filenorm = TFile::Open(fMCFileNormCorr);

        if (fDoEventLoss){
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
                hspectraSgnLoss[nmult] = (TH1F*) hspectraEvtLoss[nmult]->Clone(Form("PtSpectrumEvtLossCorrPercBin%i",nmult));
            }
        } else if (fDoMCShape){
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
                hspectraSgnLoss[nmult] = (TH1F*) hspectraIT3[nmult]->Clone(Form("PtSpectrumEvtLossCorrPercBin%i",nmult));
            }
        } else {
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
                hspectraSgnLoss[nmult] = (TH1F*) hspectra[nmult]->Clone(Form("PtSpectrumEvtLossCorrPercBin%i",nmult));
            }
        }

        //Apply correction
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

            hsgnloss[nmult] = (TH1D*)filenorm->Get(Form("SignalLoss/hsgnloss_%i_%i", (int)fpercbinlimits[nmult], (int)fpercbinlimits[nmult+1] ));
            //
            for (int bin = 1 ; bin <=  hspectraSgnLoss[nmult]->GetNbinsX(); bin ++ ){
                //
                double content = hspectraSgnLoss[nmult]->GetBinContent(bin)/hsgnloss[nmult]->GetBinContent(bin);
                double error = ErrorInRatio(hspectraSgnLoss[nmult]->GetBinContent(bin), hspectraSgnLoss[nmult]->GetBinError(bin), hsgnloss[nmult]->GetBinContent(bin), hsgnloss[nmult]->GetBinError(bin));
                hspectraSgnLoss[nmult]->SetBinContent(bin, content);
                hspectraSgnLoss[nmult]->SetBinError(bin,error);
            }    
        }          
    }

    cout<<"------------------------------ Done! ----------------------------------"<<endl;
    cout << "\n";
    cout<<"------------ Signal Loss correction applied succesfully! --------------"<<endl;
    cout<<"-----------------------------------------------------------------------"<<endl;

    //===========================================================================================
    //============================== END Signal Loss correction =================================
    //===========================================================================================
    
    //===========================================================================================
    //================================== Write on file ==========================================
    //===========================================================================================

    cout << "\n\n";
    cout<<"--------------- Result Output --------------------------"<<endl;
    cout<<" ---> Writing information to "<<fOutputDataFile<<endl;
    // Open an output file
    TFile* lResultsFile = TFile::Open(fOutputDataFile, "RECREATE");
    if (!lResultsFile || !lResultsFile->IsOpen()){
        cout<<"Error! Couldn't open file!"<<endl;
        return;
    }

    TDirectoryFile *UncorrSpectra = new TDirectoryFile("NoCorrSpectra","Spectra before final corrections");
    UncorrSpectra->cd();
    for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
        hspectra[nmult]->Write();
        if (fDoMCShape) Levy[nmult]->Write();
    }
    lResultsFile->cd();

    TDirectoryFile *MCpTshape;
    TDirectoryFile *Iteration[3];
    TDirectoryFile *MCshape[3];
    TDirectoryFile *AfterIt[3];
    TDirectoryFile *Gen[3];
    if (fDoMCShape){
        MCpTshape = new TDirectoryFile("MCpTshape","MC pT shape correction");
        MCpTshape->cd();
        Iteration[0] = new TDirectoryFile("Iteration0","MC pT shape correction");
        Iteration[0]->cd();
        MCshape[0] = new TDirectoryFile("MCshape","MC pT shape correction");
        MCshape[0]->cd();
        HistGenClone->Write();
        LevyMC->Write();
        Iteration[0]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hRatioPtShape[nmult]->Write();
            RatioEff[nmult]->Write();
        }    
        Gen[0] = new TDirectoryFile("GeneratedSpectra","MC pT shape correction");
        Gen[0]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hPtSpectrumGen[nmult]->Write();
            hPtSpectrumGenMCMB->Write();
        }        
        AfterIt[0] = new TDirectoryFile("SpectraAfterIt0","MC pT shape correction");
        AfterIt[0]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hspectraIT1[nmult]->Write();
            LevyIT1[nmult]->Write(); //to do mettere i fit 
        }    
        Iteration[0]->cd();
        MCpTshape->cd();
        //
        Iteration[1] = new TDirectoryFile("Iteration1","MC pT shape correction");
        Iteration[1]->cd();
        MCshape[1] = new TDirectoryFile("MCshape","MC pT shape correction");
        MCshape[1]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            LevyMCIT1[nmult]->Write();
            HistGenMult[nmult]->Write();
        }
        Iteration[1]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hRatioPtShapeIT1[nmult]->Write();
            RatioEffIT1[nmult]->Write();
        } 
        Gen[1] = new TDirectoryFile("GeneratedSpectra","MC pT shape correction");
        Gen[1]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hPtSpectrumGenIT1[nmult]->Write();
            hPtSpectrumGenMCIT1[nmult]->Write();
        } 
        AfterIt[1] = new TDirectoryFile("SpectraAfterIt1","MC pT shape correction");
        AfterIt[1]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hspectraIT2[nmult]->Write();
            LevyIT2[nmult]->Write(); //to do mettere i fit 
        }    
        Iteration[1]->cd();
        MCpTshape->cd();
        //
        Iteration[2] = new TDirectoryFile("Iteration2","MC pT shape correction");
        Iteration[2]->cd();
        MCshape[2] = new TDirectoryFile("MCshape","MC pT shape correction");
        MCshape[2]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            LevyMCIT2[nmult]->Write();
            HistGenMultIT1[nmult]->Write();
        }
        Iteration[2]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hRatioPtShapeIT2[nmult]->Write();
            RatioEffIT2[nmult]->Write();
        }  
        Gen[2] = new TDirectoryFile("GeneratedSpectra","MC pT shape correction");
        Gen[2]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hPtSpectrumGenIT2[nmult]->Write();
            hPtSpectrumGenMCIT2[nmult]->Write();
        } 
        AfterIt[2] = new TDirectoryFile("SpectraAfterIt2","MC pT shape correction");
        AfterIt[2]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hspectraIT3[nmult]->Write();
        }    
        Iteration[2]->cd();
        MCpTshape->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hspectraIT3[nmult]->Write();
        }    
    }
    lResultsFile->cd();
    
    TDirectoryFile *NormCorr, *EventLoss, *SignalLoss;
    NormCorr = new TDirectoryFile("NormCorr","EventLoss + Signal Loss");
    NormCorr->cd();
    if (fDoEventLoss){    
        EventLoss = new TDirectoryFile("AfterEventLossCorr","EventLoss");
        EventLoss->cd();
        heventloss->Write();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hspectraEvtLoss[nmult]->Write(); 
        }
        NormCorr->cd();
    }
    if (fDoSignalLoss){  
        SignalLoss = new TDirectoryFile("AfterSignalLossCorr","SignalLoss");
        SignalLoss->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) { 
            hspectraSgnLoss[nmult]->Write();
            hsgnloss[nmult]->Write();
        }      
    }
    lResultsFile->cd();

    //Exit Batch Mode
    gROOT->SetBatch(kFALSE);

    cout<<"--------------------------------------------------------"<<endl;
    cout<<" There, done! "<<endl;
    cout<<endl;
    cout << "\n\n";
}