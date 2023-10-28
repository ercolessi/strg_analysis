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
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <cstdlib>
using namespace std;

//--- For ROOT ---
#include "CorrectionClass.h"

CorrectionClass::CorrectionClass()
{
    fWhichParticle = "Xi"; //Default
    fOutputDataFile = "";
    fMCFilepTshape0 = "";
    fMCFilepTshape2 = "";
    fMCFilepTshape1 = "";
    fMCFileSignalLoss = "";
    fMCFileEventLoss = "";
    fWhichVarEstimator = "";
    fWhichFixedEstimator = "";
    fWhichEnergyEstimator = "";
    fWhichMultEstimator = "";
    fLowFixed = -1.;
    fHighFixed = -1.;
    fDoMCShape = kFALSE;
    fDoEventLoss = kFALSE;
    fDoSignalLoss = kFALSE;
    fDoSystematics = kFALSE;
}

CorrectionClass::CorrectionClass(TString fParticleType)
{
    // Allows definition of Particle Type in analysis.
    // Possible Options are "Xi", "Omega", "Lambda".
    // If some other string is given, this constructor will
    // default to "Xi".
    fWhichParticle = fParticleType;
    fOutputDataFile  = "";
    fMCFilepTshape0 = "";
    fMCFilepTshape2 = "";
    fMCFilepTshape1 = "";
    fMCFileSignalLoss = "";
    fMCFileEventLoss = "";
    fWhichVarEstimator = "";
    fWhichFixedEstimator = "";
    fWhichEnergyEstimator = "";
    fWhichMultEstimator = "";
    fLowFixed = -1.;
    fHighFixed = -1.;
    fDoMCShape = kFALSE;
    fDoEventLoss = kFALSE;
    fDoSignalLoss = kFALSE;
    fDoSystematics = kFALSE;
}

/***********************************************
  --- Setters For Configuration ---
************************************************/

// Filename Setters
void CorrectionClass::SetOutputDataFile  ( TString DataFilename   ){
    fOutputDataFile = DataFilename;
}
void CorrectionClass::SetMCFilepTshape  ( TString MCFilename0, TString MCFilename1, TString MCFilename2  ){
    fMCFilepTshape0 = MCFilename0;
    fMCFilepTshape1 = MCFilename1;
    fMCFilepTshape2 = MCFilename2;
}
void CorrectionClass::SetMCFileEventLoss  ( TString MCFilename  ){
    fMCFileEventLoss = MCFilename;
}
void CorrectionClass::SetMCFileSignalLoss  ( TString MCFilename  ){
    fMCFileSignalLoss = MCFilename;
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
void CorrectionClass::SetPercBinLimits(Long_t got_percbinnumb, std::vector<Double_t> got_percbinlimitsMultlow, std::vector<Double_t> got_percbinlimitsMulthigh, std::vector<Double_t> got_percbinlimitsEnergylow, std::vector<Double_t> got_percbinlimitsEnergyhigh)
{
    //Function to set pt binning. First argument is the number of pt bins, second is
    //an array with bin limits.
    fpercbinnumb = got_percbinnumb;
    for(int ix = 0;ix<fpercbinnumb;ix++){
        fpercentileMultlow[ix] = got_percbinlimitsMultlow[ix];
        fpercentileMulthigh[ix] = got_percbinlimitsMulthigh[ix];
        fpercentileEnergylow[ix] = got_percbinlimitsEnergylow[ix];
        fpercentileEnergyhigh[ix] = got_percbinlimitsEnergyhigh[ix];
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
//
void CorrectionClass::DoSystematics(Bool_t got_DoSystematics){
    fDoSystematics = got_DoSystematics;
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

Double_t CorrectionClass::ErrorInRatio( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ){
    //Error in a Ratio
    if(B!=0){
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( errorfromtop + errorfrombottom );
    }
    return 1;
}

void CorrectionClass::ErrorInRatioUncorr( TH1F* h1, TH1F *h2 ) {
    Double_t lh1NBins = h1->GetNbinsX();
    Double_t lh2NBins = h2->GetNbinsX();

    if( lh1NBins != lh2NBins ){
        cout<<"Problem! Number of bins doesn't match! "<<endl;
        return;
    }

    //Regular Division
    h1 -> Divide( h2 );
    //Replace Erorrs
    for( Int_t i=1; i<h1->GetNbinsX()+1; i++){
        h1->SetBinError(i, ErrorInRatio(h1->GetBinContent(i),h1->GetBinError(i),h2->GetBinContent(i),h2->GetBinError(i)));
    }
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

Float_t  CorrectionClass::DoWeightedMean(Bool_t err , Float_t h1, Float_t h2, Float_t h3, Float_t e1, Float_t e2, Float_t e3, Float_t w1, Float_t w2, Float_t w3) {

    Float_t wmean,  werror;

    Float_t num = h1*w1 + h2*w2 + h3*w3;
    Float_t den = w1 + w2 + w3;

    wmean = num/den;
    werror =  1./den * (TMath::Sqrt(w1*w1*TMath::Power(e1,2) + w2*w2*TMath::Power(e2,2) + w3*w3*TMath::Power(e3,2)));

    if (err) {
        return werror;
    } else {
        return wmean;
    }

}

//Printer
void CorrectionClass::PrintConfiguration() {
    //Print current analysis configuration
    Double_t lParticleMass = 1.32171; //Default
    if(fWhichParticle.Contains("Omega")) lParticleMass = 1.67245;
    if(fWhichParticle.Contains("Lambda")) lParticleMass = 1.115683;
    if(fWhichParticle.Contains("K0Short")) lParticleMass = 0.497611;


    //
    cout << "\n\n";
    cout<<"--------------- Configuration --------------------------"<<endl;
    cout<<" Analysed Particle.............: "<<fWhichParticle<<endl;
    cout<<" This Particle Mass............: "<<lParticleMass<<endl;
    cout<<" Mult estimator............: "<<fWhichMultEstimator;
    for (int i = 0; i < fpercbinnumb; i++){
        cout << " [" << fpercentileMultlow[i] << "-" << fpercentileMulthigh[i] << "]";
    }
    cout << endl;
    cout<<" Energy estimator...............: "<<fWhichEnergyEstimator;
    for (int i = 0; i < fpercbinnumb; i++){
        cout << " [" << fpercentileEnergylow[i] << "-" << fpercentileEnergyhigh[i] << "]";
    }
    cout << endl;
    cout << "\n";
    cout<<"------------- Corrections to apply----------------------"<<endl;
    cout<<" MC pT shape ...................: "<<(fDoMCShape ? "YES" : "NO")<<endl;
    cout<<" Event Loss ....................: "<<(fDoEventLoss ? "YES" : "NO")<<endl;
    cout<<" Signal Loss ...................: "<<(fDoSignalLoss ? "YES" : "NO")<<endl;
    cout<<" Sistematics ...................: "<<(fDoSystematics ? "YES" : "NO")<<endl;
    cout << "\n";

}

//----------------------------------------------------------------------------------------------------
Double_t CorrectionClass::PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin){
  Double_t dev = h->GetBinContent(bin);
  Double_t RBsigma = h->GetBinError(bin);

  if (dev>(nsigmas*RBsigma)) {return dev;}
  else {return 0.;}
}

void CorrectionClass::DoAnalysis(){


    Double_t lParticleMass = 1.32171; //Default
    if(fWhichParticle.Contains("Omega")) lParticleMass = 1.67245;
    if(fWhichParticle.Contains("Lambda")) lParticleMass = 1.115683;
    if(fWhichParticle.Contains("K0Short")) lParticleMass = 0.497611;

    TString fParticle = "";
    TString fAntiParticle = "";
    double minfit = 0., maxfit = 0.;
    if (fWhichParticle.Contains("Xi")){
      fParticle = "XiMinus";
      fAntiParticle = "XiPlus";
      minfit = 0.6;
      maxfit = 6.5;
    }
    if (fWhichParticle.Contains("Omega")){
      fParticle = "OmegaMinus";
      fAntiParticle = "OmegaPlus";
      minfit = .9;
      maxfit = 5.5;
    }
    if (fWhichParticle.Contains("Lambda")){
      fParticle = "Lambda";
      fAntiParticle = "AntiLambda";
      minfit = 0.4;
      maxfit = 8.0;
    }
    if (fWhichParticle.Contains("K0Short")){
      fParticle = "K0Short";
      fAntiParticle = "";
      minfit = 0.;
      maxfit = 10.0;
    }
    const int nbins = (maxfit-minfit)*10;
    const int plus = (minfit)*10;
    //
    TString folder = "V0Analysis/results";
    if(fParticle.Contains("Minus")) folder = "CascadeAnalysis/results";
    if(fParticle.Contains("Lambda")) folder = "V0Analysis/resultsFD";

    TFile* filepart[fpercbinnumb], *fileantipart[fpercbinnumb], *efficiencypart, *efficiencyantipart;
    TH1F* hpart[fpercbinnumb], *hantipart[fpercbinnumb], *hpartraw[fpercbinnumb], *hantipartraw[fpercbinnumb], *hspectra[fpercbinnumb], *heffpart, *heffantipart;

    TH1F *hFinalSpectraStat[fpercbinnumb], *hFinalSpectraSyst[fpercbinnumb];

    TF1 *LevyFit[fpercbinnumb];

    // Needed for pT shape correction
    TString sfilepart, sfileantipart;
    if (fDoMCShape){
        sfilepart = Form("%s/Results-%s-13TeV-SPDClusters_000_100_V0M_000_100.root", folder.Data(), fParticle.Data());
        sfileantipart = Form("%s/Results-%s-13TeV-SPDClusters_000_100_V0M_000_100.root", folder.Data(), fAntiParticle.Data());
        efficiencypart = TFile::Open(Form("%s", sfilepart.Data()));
        efficiencyantipart = TFile::Open(Form("%s", sfileantipart.Data()));
    }

    for(int nmult = 0; nmult < fpercbinnumb; nmult++) { //loop over perc selections

        double multmin, multmax, eemin, eemax;
        multmin = fpercentileMultlow[nmult];
        multmax = fpercentileMulthigh[nmult];
        eemin = fpercentileEnergylow[nmult];
        eemax = fpercentileEnergyhigh[nmult];

        //Get spectra from files
        //particle
        if (fWhichParticle.Contains("K0Short")) filepart[nmult] = new TFile(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f_FDUseMCRatio.root", folder.Data(), fParticle.Data(), fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        if (fWhichParticle.Contains("Lambda")) filepart[nmult] = new TFile(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f_FDUseMCRatio.root", folder.Data(), fParticle.Data(), fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        if (fWhichParticle.Contains("Xi") || fWhichParticle.Contains("Omega")) filepart[nmult] = new TFile(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root", folder.Data(), fParticle.Data(), fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        cout << "Opening file: " << Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root", folder.Data(), fParticle.Data(), fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax) << endl;

        hpartraw[nmult] = (TH1F *)filepart[nmult]->Get(
        Form("fHistPt%s", fParticle.Data()));
        //"lInvMassReal/lInvMassRealRawData/fHistPtRaw");
        hpart[nmult] = (TH1F*)hpartraw[nmult]->Clone(Form("hpart%i",nmult));
        //hpart[nmult]->Reset();
        //heffpart = (TH1F*) efficiencypart->Get("fHistEfficiency");
        /*for(Int_t ibin=0;ibin<hpart[nmult]->GetNbinsX();ibin++){
            hpart[nmult]->SetBinContent(ibin+1, hpartraw[nmult]->GetBinContent(ibin+1)/heffpart->GetBinContent(ibin+1));
            hpart[nmult]->SetBinError(ibin+1, ErrorInRatio(
                hpartraw[nmult]->GetBinContent(ibin+1),hpartraw[nmult]->GetBinError(ibin+1),
                heffpart->GetBinContent(ibin+1),heffpart->GetBinError(ibin+1)));
        }*/
        //anti-particle
        if (!fWhichParticle.Contains("K0Short")){
            if (fWhichParticle.Contains("Lambda")) fileantipart[nmult] = new TFile(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f_FDUseMCRatio.root", folder.Data(), fAntiParticle.Data(), fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            if (fWhichParticle.Contains("Xi") || fWhichParticle.Contains("Omega")) fileantipart[nmult] = new TFile(Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root", folder.Data(), fAntiParticle.Data(), fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            cout << "Opening file: " << Form("%s/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root", folder.Data(), fAntiParticle.Data(), fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax) << endl;

            hantipartraw[nmult] = (TH1F *) fileantipart[nmult]->Get(
            Form("fHistPt%s",fAntiParticle.Data()));
            //"lInvMassReal/lInvMassRealRawData/fHistPtRaw");
            hantipart[nmult] = (TH1F*)hantipartraw[nmult]->Clone(Form("hantipart%i",nmult));
            //hantipart[nmult]->Reset();
            //heffantipart = (TH1F*) efficiencyantipart->Get("fHistEfficiency");
            /*for(Int_t ibin=0;ibin<hantipart[nmult]->GetNbinsX();ibin++){
                hantipart[nmult]->SetBinContent(ibin+1, hantipartraw[nmult]->GetBinContent(ibin+1)/heffantipart->GetBinContent(ibin+1));
                hantipart[nmult]->SetBinError(ibin+1, ErrorInRatio(
                    hantipartraw[nmult]->GetBinContent(ibin+1),hantipartraw[nmult]->GetBinError(ibin+1),
                    heffantipart->GetBinContent(ibin+1),heffantipart->GetBinError(ibin+1)));
            }*/
        } else { // if K0Short no antiparticle
            fileantipart[nmult] = 0x0;
            hantipart[nmult] = (TH1F *)hpartraw[nmult]->Clone(Form("hantipart%i", nmult));
            for (int bin = 1 ; bin <=  hantipart[nmult]->GetNbinsX(); bin ++ ){ // loop over bins
                hantipart[nmult]->SetBinContent(bin, 0.);
                hantipart[nmult]->SetBinError(bin, 0.);
            }
        }


        // Sum spectrum
        hspectra[nmult] = (TH1F*)hpart[nmult]->Clone(Form("PtSpectrumPercBin%i",nmult));
        hspectra[nmult]->SetTitle(Form("%s_%03.0f_%03.0f - %s_%03.0f_%03.0f", fWhichMultEstimator.Data(), fpercentileMultlow[nmult], fpercentileMulthigh[nmult], fWhichEnergyEstimator.Data(), fpercentileEnergylow[nmult], fpercentileEnergyhigh[nmult]));
        hspectra[nmult]->Reset();
        for (int bin = 1 ; bin <=  hpart[nmult]->GetNbinsX(); bin ++ ){ // loop over bins
            hspectra[nmult]->SetBinContent(bin, hpart[nmult]->GetBinContent(bin) + hantipart[nmult]->GetBinContent(bin));
            hspectra[nmult]->SetBinError(bin, TMath::Sqrt(hpart[nmult]->GetBinError(bin)*hpart[nmult]->GetBinError(bin)+hantipart[nmult]->GetBinError(bin)*hantipart[nmult]->GetBinError(bin)));
        } // end loop over bins


        // Final spectra
        hFinalSpectraStat[nmult] = (TH1F *)hspectra[nmult]->Clone();
        hFinalSpectraSyst[nmult] = (TH1F *)hspectra[nmult]->Clone();
        hFinalSpectraSyst[nmult]->Reset();
        hFinalSpectraStat[nmult]->SetName(Form("FinalPtSpectrumStat_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        hFinalSpectraSyst[nmult]->SetName(Form("FinalPtSpectrumSyst_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

        // Fit function
        if (!fDoSystematics){
            LevyFit[nmult] = LevyTsallis(Form("LevyFit%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax), lParticleMass);
            LevyFit[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

            // Fit spectra with syst + stat errors
            Int_t fitres;
            Int_t trials = 0;
            trials = 0;
            do
            {
                fitres = hFinalSpectraStat[nmult]->Fit(LevyFit[nmult], "0q", "", minfit, maxfit);
                Printf("Fit trial: %d", trials++);
                if (trials > 10)
                {
                    Printf("FIT DOES NOT CONVERGE IN LINE %d", __LINE__);
                    break;
                }
            } while (fitres != 0);
        }
    } // end loop over perc selections

    //===========================================================================================
    // ============================= pT shape MC correction =====================================
    //===========================================================================================

    TH1D* HistGenPart[3], *HistGenAntiPart[3], *HistGen[3], *HistGenClone[3];
    TH1D* HistRecoPart[3];
    TH1D* HistRecoAntiPart[3];
    TH1D* HistReco[3];
    TH1F* HistRecoMultIT1[fpercbinnumb][3],* HistRecoMultIT2[fpercbinnumb][3],* HistRecoMultIT3[fpercbinnumb][3];
    TH1F* HistGenMultIT1[fpercbinnumb][3],* HistGenMultIT2[fpercbinnumb][3],* HistGenMultIT3[fpercbinnumb][3], * HistGenMultCloneIT1[fpercbinnumb][3],* HistGenMultCloneIT2[fpercbinnumb][3],* HistGenMultCloneIT3[fpercbinnumb][3];
    Double_t par_n[fpercbinnumb];
    Double_t par_C[fpercbinnumb];
    Double_t par_norm[fpercbinnumb];
    TF1* Levy[fpercbinnumb], *LevyIT1[fpercbinnumb], *LevyIT2[fpercbinnumb], *LevyIT3[fpercbinnumb];
    TF1* LevyMCIT1[fpercbinnumb][3], *LevyMCIT2[fpercbinnumb][3], *LevyMCIT3[fpercbinnumb][3];
    TF1* LevyMC[3];
    TH1F *hPtSpectrumGen[fpercbinnumb];
    TH1F *hPtSpectrumGenMCMB[3];
    TH1F *hPtSpectrumGenMCIT1[fpercbinnumb][3], *hPtSpectrumGenMCIT2[fpercbinnumb][3], *hPtSpectrumGenMCIT3[fpercbinnumb][3];
    TFile* fileMC0, *fileMC1, *fileMC2;
    TList* clistMC0, * clistMC1, * clistMC2;
    TH3D* h3DGeneratedPart[3], *h3DGeneratedAntiPart[3];
    Int_t rapBin_min;
    Int_t rapBin_max;
    Int_t multBin_min;
    Int_t multBin_max;
    TH1F* hRatioPtShape[fpercbinnumb][3], * hRatioPtShapeIT1[fpercbinnumb][3], * hRatioPtShapeIT2[fpercbinnumb][3], * hRatioPtShapeIT3[fpercbinnumb][3];
    TH1F* HistRecoMult[fpercbinnumb][3], *HistGenMult[fpercbinnumb][3], *HistGenMultClone[fpercbinnumb][3], *HistRecoMult_rebin[fpercbinnumb][3], *HistGenMult_rebin[fpercbinnumb][3], *CorrEfficiencyparz[fpercbinnumb][3], *CorrEfficiencyIT1parz[fpercbinnumb][3], *CorrEfficiencyIT2parz[fpercbinnumb][3], *CorrEfficiencyIT3parz[fpercbinnumb][3], *CorrEfficiency[fpercbinnumb], *CorrEfficiencyIT1[fpercbinnumb], *CorrEfficiencyIT2[fpercbinnumb], *CorrEfficiencyIT3[fpercbinnumb];
    TH1F* RatioEff[fpercbinnumb],* RatioEffIT1[fpercbinnumb],* RatioEffIT2[fpercbinnumb],* RatioEffIT3[fpercbinnumb];
    TH1F* RatioEffparz[fpercbinnumb][3], * RatioEffparzIT1[fpercbinnumb][3],* RatioEffparzIT2[fpercbinnumb][3],* RatioEffparzIT3[fpercbinnumb][3];
    TH1F* HistReco_rebin[3];
    TH1F* HistGen_rebin[3];
    TH1F* EfficiencyMB[3];
    TH1F* EfficiencyIT0;
    TH1F* hspectraIT1[fpercbinnumb],* hspectraIT2[fpercbinnumb],* hspectraIT3[fpercbinnumb], * hspectraIT4[fpercbinnumb];
    TH1F *hPtSpectrumGenIT1[fpercbinnumb], *hPtSpectrumGenIT2[fpercbinnumb], *hPtSpectrumGenIT3[fpercbinnumb];
    TH1F *hPtSpectrumGenMCMBIT1, *hPtSpectrumGenMCMBIT2, *hPtSpectrumGenMCMBIT3;
    Int_t colors[] = {kRed+1, kRed-4, kOrange+7, kOrange-3, kYellow+1, kSpring-7, kGreen+2, kAzure+8, kBlue-4, kBlue+3};
    Double_t w0 = 3.591657E+7;
    Double_t w1 = 2.984665E+7;
    Double_t w2 = 5.014868E+7;

    if (fDoMCShape){

        cout << "\n\n";
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"-------------------- Doing MC pT shape correction --------------------"<<endl;
        cout<<"----------------------------------------------------------------------"<<endl;
        cout << "\n";
        cout<<"--------------------------- Iteration 0 ------------------------------"<<endl;

        //Get Generated MC Histo (input MC pT-shape)
        fileMC0 = new TFile(fMCFilepTshape0,"READ");
        fileMC1 = new TFile(fMCFilepTshape1,"READ");
        fileMC2 = new TFile(fMCFilepTshape2,"READ");
        //
        clistMC0 = (TList*)fileMC0->Get("PWGLF_StrVsMult_MC/cList");
        clistMC1 = (TList*)fileMC1->Get("PWGLF_StrVsMult/cList");
        clistMC2 = (TList*)fileMC2->Get("PWGLF_StrVsMult/cList");
        //
        h3DGeneratedPart[0] = (TH3D*)clistMC0->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", fParticle.Data()));
        h3DGeneratedAntiPart[0] = (TH3D*)clistMC0->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", fAntiParticle.Data()));
        h3DGeneratedPart[1] = (TH3D*)clistMC1->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", fParticle.Data()));
        h3DGeneratedAntiPart[1] = (TH3D*)clistMC1->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", fAntiParticle.Data()));
        h3DGeneratedPart[2] = (TH3D*)clistMC2->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", fParticle.Data()));
        h3DGeneratedAntiPart[2] = (TH3D*)clistMC2->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", fAntiParticle.Data()));
        // find projection bins in rapidity
        rapBin_min = h3DGeneratedPart[0]->GetYaxis()->FindBin( -0.5+1.e-6 );
        rapBin_max = h3DGeneratedPart[0]->GetYaxis()->FindBin( +0.5-1.e-6 );
        //
        // find projection bins in multiplicity
        multBin_min = h3DGeneratedPart[0]->GetZaxis()->FindBin( 0.+1.e-6 );
        multBin_max = h3DGeneratedPart[0]->GetZaxis()->FindBin( 100.-1.e-6 );



        for(int imc = 0; imc < 3; imc ++){ //loop over mc files

            //===================================GENERATED====================================
            HistGenPart[imc] = (TH1D*)h3DGeneratedPart[imc]->ProjectionX(Form("HistGen%s_MC%i", fParticle.Data(),imc),rapBin_min,rapBin_max, multBin_min, multBin_max);
            HistGenAntiPart[imc] = (TH1D*)h3DGeneratedAntiPart[imc]->ProjectionX(Form("HistGen%s_MC%i", fAntiParticle.Data(),imc),rapBin_min,rapBin_max, multBin_min, multBin_max);
            HistGen[imc] = (TH1D*)HistGenPart[imc]->Clone(Form("HistGen%i",imc));
            HistGen[imc]->Reset();


            for (int bin = 1; bin <= HistGen[imc]->GetNbinsX(); bin++){
            HistGen[imc]->SetBinContent(bin, HistGenPart[imc]->GetBinContent(bin) + HistGenAntiPart[imc]->GetBinContent(bin));
            HistGen[imc]->SetBinError(bin, TMath::Sqrt(HistGenPart[imc]->GetBinError(bin)*HistGenPart[imc]->GetBinError(bin)+
                HistGenAntiPart[imc]->GetBinError(bin)*HistGenAntiPart[imc]->GetBinError(bin)));
            }

            // Rebin Gen
            HistGen_rebin[imc] = new TH1F(Form("HistGen%i_rebin", imc), ";p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);
            Double_t tempptgen;
            for(long i = 1; i<HistGen[imc]->GetNbinsX()+1;i++){
                tempptgen = HistGen[imc]->GetXaxis()->GetBinCenter(i);
                for(long filling = 0; filling<HistGen[imc]->GetBinContent(i); filling++){
                    HistGen_rebin[imc]->Fill(tempptgen);
                }
            }

            //HERE
            //Fit this already (clone HistGen for high granularity or HistGen_rebin for low granularity)
            //My hypothesis is that Fiorella used low granularity for the fit
            //But maybe the correct way would be use directly the input shape --> try to be consistent with Fiorella
            HistGenClone[imc] = (TH1D *)HistGen_rebin[imc]->Clone(Form("HistGen%i", imc));
            HistGenClone[imc]->Reset();
            //Rebin
            for (int ibin = 1; ibin <= HistGen_rebin[imc]->GetNbinsX(); ibin++){
                HistGenClone[imc]->SetBinContent(ibin, HistGen_rebin[imc]->GetBinContent(ibin) / HistGen_rebin[imc]->GetBinWidth(ibin));
                HistGenClone[imc]->SetBinError(ibin, HistGen_rebin[imc]->GetBinError(ibin) / HistGen_rebin[imc]->GetBinWidth(ibin));
            }
            double integral = HistGenClone[imc]->Integral(0,-1);
            HistGenClone[imc]->Scale(1./integral);

            //
            LevyMC[imc] = LevyTsallis(Form("LevyMCIT0FitMB%i",imc), lParticleMass,1E-3,1E+3);
            Bool_t okfit = HistGenClone[imc]->Fit(LevyMC[imc],"0q","", minfit, maxfit);
            if (okfit) cout << "FIT DOES NOT CONVERGE IN LINE " << __LINE__ << endl;

            //=====================================RECO=========================================
            //Get Reco MC histo
            HistRecoPart[imc] = (TH1D*)efficiencypart->Get(Form("lInvMassMC/fHistReco%i",imc));
            HistRecoAntiPart[imc] = (TH1D*)efficiencyantipart->Get(Form("lInvMassMC/fHistReco%i",imc));
            HistReco[imc] = (TH1D*) HistRecoPart[imc]->Clone(Form("fHistReco%i",imc));
            HistReco[imc]->Reset();

            for (int bin = 1 ; bin <=  HistRecoPart[imc]->GetNbinsX(); bin ++ ){
                HistReco[imc]->SetBinContent(bin, HistRecoPart[imc]->GetBinContent(bin) + HistRecoAntiPart[imc]->GetBinContent(bin));
                HistReco[imc]->SetBinError(bin, TMath::Sqrt(HistRecoPart[imc]->GetBinError(bin)*HistRecoPart[imc]->GetBinError(bin)+HistRecoAntiPart[imc]->GetBinError(bin)*HistRecoAntiPart[imc]->GetBinError(bin)));
            }


            //====================================================================================================
            //Rebin Reco per fare efficienze default
            HistReco_rebin[imc] = new TH1F(Form("HistReco%i_rebin", imc), ";p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);
            //Let's define a total efficiency, it will be usefull aftewards
            EfficiencyMB[imc]  = new TH1F(Form("EfficiencyMB%i",imc),"Eff #Xi+#bar{#Xi};p_{T} (GeV/c);Eff", fptbinnumb, fptbinlimits);

            Double_t tempptreco;
            for(long i = 1; i<HistReco[imc]->GetNbinsX()+1;i++){
                tempptreco = HistReco[imc]->GetXaxis()->GetBinCenter(i);
                for(long filling = 0; filling<HistReco[imc]->GetBinContent(i); filling++){
                    HistReco_rebin[imc]->Fill(tempptreco);
                }
            }

            for (int bin = 1; bin <= EfficiencyMB[imc]->GetNbinsX(); bin ++){
                if (HistGen_rebin[imc]->GetBinContent(bin) != 0){
                    EfficiencyMB[imc]->SetBinContent(bin , HistReco_rebin[imc]->GetBinContent(bin)/HistGen_rebin[imc]->GetBinContent(bin));
                    EfficiencyMB[imc]->SetBinError(bin, ErrorInRatio(HistReco_rebin[imc]->GetBinContent(bin),HistReco_rebin[imc]->GetBinError(bin),
                        HistGen_rebin[imc]->GetBinContent(bin),HistGen_rebin[imc]->GetBinError(bin)));
                }
            }


        }//end loop over mc files

        EfficiencyIT0 = (TH1F*)EfficiencyMB[0]->Clone("EfficiencyIT0");
        EfficiencyIT0->Reset();
        for (int bin = 1; bin <= EfficiencyIT0->GetNbinsX(); bin ++){
            EfficiencyIT0->SetBinContent(bin,DoWeightedMean(0,
                    EfficiencyMB[0]->GetBinContent(bin), EfficiencyMB[1]->GetBinContent(bin), EfficiencyMB[2]->GetBinContent(bin),
                    EfficiencyMB[0]->GetBinError(bin), EfficiencyMB[1]->GetBinError(bin), EfficiencyMB[2]->GetBinError(bin),
                    w0,w1,w2));
            EfficiencyIT0->SetBinError(bin,DoWeightedMean(1,
                    EfficiencyMB[0]->GetBinContent(bin), EfficiencyMB[1]->GetBinContent(bin), EfficiencyMB[2]->GetBinContent(bin),
                    EfficiencyMB[0]->GetBinError(bin), EfficiencyMB[1]->GetBinError(bin), EfficiencyMB[2]->GetBinError(bin),
                    w0,w1,w2));
        }

        for(int nmult = 0; nmult < fpercbinnumb; nmult++){ //loop over perc selections

            double multmin, multmax, eemin, eemax;
            //
            multmin = fpercentileMultlow[nmult];
            multmax = fpercentileMulthigh[nmult];
            eemin = fpercentileEnergylow[nmult];
            eemax = fpercentileEnergyhigh[nmult];

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
            hPtSpectrumGen[nmult] = new TH1F(Form("hPtSpectrumGenIT0PercBin%i",nmult), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", nbins, minfit, maxfit);
            hPtSpectrumGen[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

            //Define ratio shape
            for(int imc = 0; imc < 3; imc ++){ //loop over mc files
                hRatioPtShape[nmult][imc] = (TH1F*)hPtSpectrumGen[nmult]->Clone(Form("hRatioPtShapeIT0PercBin%i_MC%i",nmult,imc));
                hRatioPtShape[nmult][imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
                hRatioPtShape[nmult][imc]->Reset();
            }
        }// end loop over perc selections

        //++++++++++++++++++++++++++++++++++++++++
        // 0) Generate histos from Fit Levy-Tsallis
        //++++++++++++++++++++++++++++++++++++++++
        //The pT-spectra for both data and Monte Carlo were fitted by a Levy-Tsallis function.
        //Here pT-spectra are “generated” according to the fitted shapes using an high granularity.

        const int n = 100000000;

        // multiplicity spectra RD
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
            for (Int_t i = 0; i < n; i++) {
                double x = (double)Levy[nmult]->GetRandom(minfit, maxfit);
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

        for(int imc = 0; imc < 3; imc ++){ //loop over mc files

            hPtSpectrumGenMCMB[imc] = new TH1F(Form("hPtSpectrumGenMCIT0MB%i",imc), "Cascade generated from Levy-Tsallis Fit;#it{p}_{T} (GeV/#it{c});Counts", nbins, minfit,maxfit);
            hPtSpectrumGenMCMB[imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), 0., 100., fWhichEnergyEstimator.Data(), 0., 100.));

            //MB spectra
            for (Int_t i = 0; i < n; i++) { //HERE
                //double xMC = (double)HistGenClone[imc]->GetRandom();
                double xMC = (double)LevyMC[imc]->GetRandom(minfit, maxfit);
                while (xMC<minfit || xMC>maxfit)
                {
                    xMC = (double)HistGenClone[imc]->GetRandom();
                }
                hPtSpectrumGenMCMB[imc]->Fill(xMC);
            }
            // /dpT
            for (Int_t i = 0; i <= hPtSpectrumGenMCMB[imc]->GetNbinsX(); i++) {
                double entryMC = hPtSpectrumGenMCMB[imc]->GetBinContent(i + 1) /
                                hPtSpectrumGenMCMB[imc]->GetBinWidth(i + 1);
                hPtSpectrumGenMCMB[imc]->SetBinContent(i + 1, entryMC);
            }

            //++++++++++++++++++++++++++++++++++++++++
            // 1) Iterative process for correction
            //++++++++++++++++++++++++++++++++++++++++
            //The ratios of the “fitted” shapes obtained from the data over the Monte Carlo input pT-shape
            //are computed with the high granularity

            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

                for (int bin = 1; bin <= hRatioPtShape[0][imc]->GetNbinsX(); bin ++){
                    if (hPtSpectrumGenMCMB[imc]->GetBinContent(bin) != 0 ){
                        hRatioPtShape[nmult][imc]->SetBinContent(bin, hPtSpectrumGen[nmult]->GetBinContent(bin)/hPtSpectrumGenMCMB[imc]->GetBinContent(bin));
                        hRatioPtShape[nmult][imc]->SetBinError(bin, ErrorInRatio(hPtSpectrumGen[nmult]->GetBinContent(bin),hPtSpectrumGen[nmult]->GetBinError(bin),
                                        hPtSpectrumGenMCMB[imc]->GetBinContent(bin),hPtSpectrumGenMCMB[imc]->GetBinError(bin)));
                    } else hRatioPtShape[nmult][imc]->SetBinContent(bin,1); //x*1=x
                }
            }

            //++++++++++++++++++++++++++++++++++++++++
            // 2) Re-weight reco and gen spectra
            //++++++++++++++++++++++++++++++++++++++++
            //The ratios at the previous point were used to “re-weight” the reconstructed and generated pT spectra
            //used to compute efficiencies (same high granularity)

            for (int nmult = 0; nmult < fpercbinnumb; nmult++) { //loop over perc selections

                double multmin, multmax, eemin, eemax;
                //
                multmin = fpercentileMultlow[nmult];
                multmax = fpercentileMulthigh[nmult];
                eemin = fpercentileEnergylow[nmult];
                eemax = fpercentileEnergyhigh[nmult];

                HistRecoMult_rebin[nmult][imc]   = new TH1F(Form("HistRecoMult_rebin%i_MC%i",nmult,imc),"Cascade MC count;p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);
                HistGenMult_rebin[nmult][imc]   = new TH1F(Form("HistGenMult_rebin%i_MC%i",nmult,imc),"Cascade MC count;p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);
                HistRecoMult[nmult][imc] = new TH1F(Form("HistRecoPercBin%i_MC%i",nmult,imc), ";#it{p}_{T} (GeV/#it{c});Counts", nbins, minfit,maxfit);
                HistGenMult[nmult][imc] = new TH1F(Form("HistGenPercBin%i_MC%i",nmult,imc), ";#it{p}_{T} (GeV/#it{c});Counts", nbins, minfit,maxfit);

                for (int bin = 1; bin <= HistRecoMult[nmult][imc]->GetNbinsX(); bin ++){
                    //
                    HistRecoMult[nmult][imc]-> SetBinContent(bin,HistReco[imc]->GetBinContent(bin+plus)*hRatioPtShape[nmult][imc]->GetBinContent(bin));
                    HistRecoMult[nmult][imc]-> SetBinError(bin,ErrorInRatio(HistReco[imc]->GetBinContent(bin+plus), HistReco[imc]->GetBinError(bin+plus),
                        hRatioPtShape[nmult][imc]->GetBinContent(bin), hRatioPtShape[nmult][imc]->GetBinError(bin))
                        );

                    //
                    HistGenMult[nmult][imc]-> SetBinContent(bin, HistGen[imc]->GetBinContent(bin+plus)*hRatioPtShape[nmult][imc]->GetBinContent(bin));
                    HistGenMult[nmult][imc]-> SetBinError(bin, ErrorInRatio(HistGen[imc]->GetBinContent(bin+plus), HistGen[imc]->GetBinError(bin+plus),
                        hRatioPtShape[nmult][imc]->GetBinContent(bin), hRatioPtShape[nmult][imc]->GetBinError(bin)
                        ));

                }

                Double_t tempptreco;
                for(long i = 1; i<HistRecoMult[nmult][imc]->GetNbinsX()+1;i++){
                    tempptreco = HistRecoMult[nmult][imc]->GetXaxis()->GetBinCenter(i);
                    for(long filling = 0; filling<HistRecoMult[nmult][imc]->GetBinContent(i); filling++){
                        HistRecoMult_rebin[nmult][imc]->Fill(tempptreco);
                    }
                }
                Double_t tempptgen;
                for(long i = 1; i<HistGenMult[nmult][imc]->GetNbinsX()+1;i++){
                    tempptgen = HistGenMult[nmult][imc]->GetXaxis()->GetBinCenter(i);
                    for(long filling = 0; filling<HistGenMult[nmult][imc]->GetBinContent(i); filling++){
                        HistGenMult_rebin[nmult][imc]->Fill(tempptgen);
                    }
                }

                //HERE --> vedi sopra se usare rebin o high granularity
                //Clone for Fit
                HistGenMultClone[nmult][imc] = (TH1F *)HistGenMult_rebin[nmult][imc]->Clone(Form("HistGenMultClone%i_MC%i", nmult, imc));
                HistGenMultClone[nmult][imc]->Reset();
                //
                for (int ibin = 1; ibin <= HistGenMult_rebin[nmult][imc]->GetNbinsX(); ibin++){
                    HistGenMultClone[nmult][imc]->SetBinContent(ibin, HistGenMult_rebin[nmult][imc]->GetBinContent(ibin) / HistGenMult_rebin[nmult][imc]->GetBinWidth(ibin));
                    HistGenMultClone[nmult][imc]->SetBinError(ibin, HistGenMult_rebin[nmult][imc]->GetBinError(ibin) / HistGenMult_rebin[nmult][imc]->GetBinWidth(ibin));
                }
                double integral = HistGenMultClone[nmult][imc]->Integral(0,-1);
                HistGenMultClone[nmult][imc]->Scale(1./integral);
                //

                //++++++++++++++++++++++++++++++++++++++++
                // 3) Correct Efficiency
                //++++++++++++++++++++++++++++++++++++++++
                //New efficiencies are recomputed by “re-binning” the pT-spectra at the point 2) in the same pT
                //bins used for the analysis. In this way the “corrected” efficiencies are obtained.


                CorrEfficiencyparz[nmult][imc] = (TH1F*)HistRecoMult_rebin[nmult][imc]->Clone(Form("hCorrEfficiencyparzIT0PercBin%i_MC%i",nmult,imc));
                CorrEfficiencyparz[nmult][imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
                CorrEfficiencyparz[nmult][imc]->Reset();

                for (int bin = 1; bin <= HistRecoMult_rebin[nmult][imc]->GetNbinsX(); bin ++){
                    if (HistGenMult_rebin[nmult][imc]->GetBinContent(bin) != 0){
                        CorrEfficiencyparz[nmult][imc]->SetBinContent(bin , HistRecoMult_rebin[nmult][imc]->GetBinContent(bin)/HistGenMult_rebin[nmult][imc]->GetBinContent(bin));
                        CorrEfficiencyparz[nmult][imc]->SetBinError(bin, ErrorInRatio(HistRecoMult_rebin[nmult][imc]->GetBinContent(bin),HistRecoMult_rebin[nmult][imc]->GetBinError(bin),
                        HistGenMult_rebin[nmult][imc]->GetBinContent(bin),HistGenMult_rebin[nmult][imc]->GetBinError(bin)));
                    }
                }
                //
                RatioEffparz[nmult][imc] = (TH1F*)CorrEfficiencyparz[nmult][imc]->Clone(Form("CorrFactorsParzIT0PercBin%i_MC%i",nmult,imc));
                RatioEffparz[nmult][imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

            } //end loop over perc selections


        }//end loop over mc files



        for (int nmult = 0; nmult < fpercbinnumb; nmult++) { //loop pver perc selections

            double multmin, multmax, eemin, eemax;
            //
            multmin = fpercentileMultlow[nmult];
            multmax = fpercentileMulthigh[nmult];
            eemin = fpercentileEnergylow[nmult];
            eemax = fpercentileEnergyhigh[nmult];

            //Weighted efficiency after IT0
            CorrEfficiency[nmult] = (TH1F*)CorrEfficiencyparz[nmult][0]->Clone(Form("CorrEfficiencyIT0PercBin%i",nmult));
            CorrEfficiency[nmult]->Reset();
            for (int bin = 1; bin <= CorrEfficiency[nmult]->GetNbinsX(); bin ++){
                CorrEfficiency[nmult]->SetBinContent(bin,DoWeightedMean(0,
                        CorrEfficiencyparz[nmult][0]->GetBinContent(bin), CorrEfficiencyparz[nmult][1]->GetBinContent(bin), CorrEfficiencyparz[nmult][2]->GetBinContent(bin),
                        CorrEfficiencyparz[nmult][0]->GetBinError(bin), CorrEfficiencyparz[nmult][1]->GetBinError(bin), CorrEfficiencyparz[nmult][2]->GetBinError(bin),
                        w0,w1,w2));
                CorrEfficiency[nmult]->SetBinError(bin,DoWeightedMean(1,
                        CorrEfficiencyparz[nmult][0]->GetBinContent(bin), CorrEfficiencyparz[nmult][1]->GetBinContent(bin), CorrEfficiencyparz[nmult][2]->GetBinContent(bin),
                        CorrEfficiencyparz[nmult][0]->GetBinError(bin), CorrEfficiencyparz[nmult][1]->GetBinError(bin), CorrEfficiencyparz[nmult][2]->GetBinError(bin),
                        w0,w1,w2));
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
        for(int imc = 0; imc < 3; imc ++){ //loop over mc files
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
            ErrorInRatioCorr(RatioEffparz[nmult][imc],EfficiencyMB[imc]);
            }
        }

        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

            ErrorInRatioCorr(RatioEff[nmult],EfficiencyIT0);

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

        for(int imc = 0; imc < 3; imc ++){ //loop over mc files

            for(int nmult = 0; nmult < fpercbinnumb; nmult++){ //loop over perc selections

                double multmin, multmax, eemin, eemax;
                //
                multmin = fpercentileMultlow[nmult];
                multmax = fpercentileMulthigh[nmult];
                eemin = fpercentileEnergylow[nmult];
                eemax = fpercentileEnergyhigh[nmult];

                LevyMCIT1[nmult][imc] = LevyTsallis(Form("LevyMCIT1FitPercBin%i_MC%i",nmult,imc),lParticleMass);
                //,1e-3 , 1e+6, 1e-8, 1e+8,  5.,  0.1, 0.1);
                LevyMCIT1[nmult][imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

                //MC shape ---
                //
                Int_t fitres;
                Int_t trials = 0;
                trials = 0;
                do {
                    fitres =  HistGenMultClone[nmult][imc]->Fit(LevyMCIT1[nmult][imc], "0q","",minfit,maxfit);
                    Printf("Fit trial: %d", trials++);
                    if(trials > 10) {
                        Printf("FIT DOES NOT CONVERGE MC IN LINE %d",__LINE__);
                        break;
                    }
                }
                while (fitres != 0);

                if(imc==0){
                    //Define generated histos
                    hPtSpectrumGenIT1[nmult] = new TH1F(Form("hPtSpectrumGenIT1PercBin%i",nmult), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", nbins, minfit,maxfit);
                    hPtSpectrumGenIT1[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
                }

                //Define ratio shape
                hRatioPtShapeIT1[nmult][imc] = (TH1F*)hPtSpectrumGenIT1[nmult]->Clone(Form("hRatioPtShapeIT1PercBin%i_MC%i",nmult,imc));
                hRatioPtShapeIT1[nmult][imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
                hRatioPtShapeIT1[nmult][imc]->Reset();

                hPtSpectrumGenMCIT1[nmult][imc] = new TH1F(Form("hPtSpectrumGenMCIT1PercBin%i_MC%i",nmult,imc), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts",nbins, minfit,maxfit);
                hPtSpectrumGenMCIT1[nmult][imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            }// end loop over perc selections
        }


        for(int nmult = 0; nmult < fpercbinnumb; nmult++){ //loop over perc selections

            double multmin, multmax, eemin, eemax;
            //
            multmin = fpercentileMultlow[nmult];
            multmax = fpercentileMulthigh[nmult];
            eemin = fpercentileEnergylow[nmult];
            eemax = fpercentileEnergyhigh[nmult];

            //Fit function
            LevyIT1[nmult] = LevyTsallis(Form("LevyIT1FitPercBin%i",nmult),lParticleMass);
            LevyIT1[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

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

        }// end loop over perc selections

        //++++++++++++++++++++++++++++++++++++++++++++++
        // 0_IT1) Generate histos from Fit Levy-Tsallis
        //++++++++++++++++++++++++++++++++++++++++++++++

        for(int imc = 0; imc < 3; imc ++){ //loop over mc files

            // multiplicity spectra
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

                //MC spectra
                for (Int_t i = 0; i < n; i++) { //HERE
                    double xMC = (double)LevyMCIT1[nmult][imc]->GetRandom(minfit, maxfit);
                    //double xMC = (double)HistGenMultClone[nmult][imc]->GetRandom();
                    hPtSpectrumGenMCIT1[nmult][imc]->Fill(xMC);
                }
                // /dpT
                for (Int_t i = 0; i <= hPtSpectrumGenMCIT1[nmult][imc]->GetNbinsX(); i++) {
                    double entryMC = hPtSpectrumGenMCIT1[nmult][imc]->GetBinContent(i + 1) /
                                    hPtSpectrumGenMCIT1[nmult][imc]->GetBinWidth(i + 1);
                    hPtSpectrumGenMCIT1[nmult][imc]->SetBinContent(i + 1, entryMC);
                }
            }
        } // end loop over mc files

        // multiplicity spectra
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

            //Data
            for (Int_t i = 0; i < n; i++) {
                double x = (double)LevyIT1[nmult]->GetRandom(minfit,maxfit);
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
        for(int imc = 0; imc < 3; imc ++){ //loop over mc files

            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

                for (int bin = 1; bin <= hRatioPtShapeIT1[0][imc]->GetNbinsX(); bin ++){
                    if (hPtSpectrumGenMCIT1[nmult][imc]->GetBinContent(bin) != 0 ){
                        hRatioPtShapeIT1[nmult][imc]->SetBinContent(bin, hPtSpectrumGenIT1[nmult]->GetBinContent(bin)/hPtSpectrumGenMCIT1[nmult][imc]->GetBinContent(bin));
                        hRatioPtShapeIT1[nmult][imc]->SetBinError(bin, ErrorInRatio(hPtSpectrumGenIT1[nmult]->GetBinContent(bin),hPtSpectrumGenIT1[nmult]->GetBinError(bin),
                                        hPtSpectrumGenMCIT1[nmult][imc]->GetBinContent(bin),hPtSpectrumGenMCIT1[nmult][imc]->GetBinError(bin)));
                    }
                    else hRatioPtShapeIT1[nmult][imc]->SetBinContent(bin,1); //x*1=x
                }
            }

            //++++++++++++++++++++++++++++++++++++++++
            // 2_IT1) Re-weight reco and gen spectra
            //++++++++++++++++++++++++++++++++++++++++

            for (int nmult = 0; nmult < fpercbinnumb; nmult++) { //loop pver perc selections

                double multmin, multmax, eemin, eemax;
                //
                multmin = fpercentileMultlow[nmult];
                multmax = fpercentileMulthigh[nmult];
                eemin = fpercentileEnergylow[nmult];
                eemax = fpercentileEnergyhigh[nmult];

                HistRecoMult_rebin[nmult][imc]->Reset();
                HistGenMult_rebin[nmult][imc]->Reset();

                HistRecoMultIT1[nmult][imc] = (TH1F*) HistRecoMult[nmult][imc]->Clone(Form("HistRecoIT1PercBin%i_MC%i",nmult,imc));
                HistRecoMultIT1[nmult][imc] ->Reset();
                HistGenMultIT1[nmult][imc] = (TH1F*) HistGenMult[nmult][imc]->Clone(Form("HistGenIT1PercBin%i_MC%i",nmult,imc));
                HistGenMultIT1[nmult][imc] ->Reset();
                for (int bin = 1; bin <= HistRecoMultIT1[nmult][imc]->GetNbinsX(); bin ++){
                    //
                    HistRecoMultIT1[nmult][imc]-> SetBinContent(bin,HistRecoMult[nmult][imc]->GetBinContent(bin)*hRatioPtShapeIT1[nmult][imc]->GetBinContent(bin));
                    HistRecoMultIT1[nmult][imc]-> SetBinError(bin,ErrorInRatio(HistRecoMult[nmult][imc]->GetBinContent(bin), HistRecoMult[nmult][imc]->GetBinError(bin),
                        hRatioPtShapeIT1[nmult][imc]->GetBinContent(bin), hRatioPtShapeIT1[nmult][imc]->GetBinError(bin))
                        );
                    //
                    HistGenMultIT1[nmult][imc]-> SetBinContent(bin, HistGenMult[nmult][imc]->GetBinContent(bin)*hRatioPtShapeIT1[nmult][imc]->GetBinContent(bin));
                    HistGenMultIT1[nmult][imc]-> SetBinError(bin, ErrorInRatio(HistGenMult[nmult][imc]->GetBinContent(bin), HistGenMult[nmult][imc]->GetBinError(bin),
                        hRatioPtShapeIT1[nmult][imc]->GetBinContent(bin), hRatioPtShapeIT1[nmult][imc]->GetBinError(bin)
                        ));
                }

                Double_t tempptreco;
                for(long i = 1; i<HistRecoMultIT1[nmult][imc]->GetNbinsX()+1;i++){
                    tempptreco = HistRecoMultIT1[nmult][imc]->GetXaxis()->GetBinCenter(i);
                    for(long filling = 0; filling<HistRecoMultIT1[nmult][imc]->GetBinContent(i); filling++){
                        HistRecoMult_rebin[nmult][imc]->Fill(tempptreco);
                    }
                }
                Double_t tempptgen;
                for(long i = 1; i<HistGenMultIT1[nmult][imc]->GetNbinsX()+1;i++){
                    tempptgen = HistGenMultIT1[nmult][imc]->GetXaxis()->GetBinCenter(i);
                    for(long filling = 0; filling<HistGenMultIT1[nmult][imc]->GetBinContent(i); filling++){
                        HistGenMult_rebin[nmult][imc]->Fill(tempptgen);
                    }
                }

                //HERE --> vedi sopra se usare rebinned o high granularity
                //Clone for Fit
                HistGenMultCloneIT1[nmult][imc] = (TH1F*) HistGenMult_rebin[nmult][imc]->Clone(Form("HistGenMultCloneIT1%i_MC%i", nmult,imc));
                HistGenMultCloneIT1[nmult][imc]->Reset();
                //
                for (int ibin = 1; ibin <= HistGenMult[nmult][imc]->GetNbinsX(); ibin++){
                    HistGenMultCloneIT1[nmult][imc]->SetBinContent(ibin, HistGenMult_rebin[nmult][imc]->GetBinContent(ibin)/HistGenMult_rebin[nmult][imc]->GetBinWidth(ibin));
                    HistGenMultCloneIT1[nmult][imc]->SetBinError(ibin, HistGenMult_rebin[nmult][imc]->GetBinError(ibin)/HistGenMult_rebin[nmult][imc]->GetBinWidth(ibin));
                }
                double integral = HistGenMultCloneIT1[nmult][imc]->Integral(0,-1);
                HistGenMultCloneIT1[nmult][imc]->Scale(1./integral);
                //

                //++++++++++++++++++++++++++++++++++++++++
                // 3) Correct Efficiency
                //++++++++++++++++++++++++++++++++++++++++

                //
                CorrEfficiencyIT1parz[nmult][imc] = (TH1F*)HistRecoMult_rebin[nmult][imc]->Clone(Form("hCorrEfficiencyIT1parzPercBin%i_MC%i",nmult,imc));
                CorrEfficiencyIT1parz[nmult][imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
                CorrEfficiencyIT1parz[nmult][imc]->Reset();

                for (int bin = 1; bin <= HistRecoMult_rebin[nmult][imc]->GetNbinsX(); bin ++){
                    if (HistGenMult_rebin[nmult][imc]->GetBinContent(bin) != 0){
                        CorrEfficiencyIT1parz[nmult][imc]->SetBinContent(bin , HistRecoMult_rebin[nmult][imc]->GetBinContent(bin)/HistGenMult_rebin[nmult][imc]->GetBinContent(bin));
                        CorrEfficiencyIT1parz[nmult][imc]->SetBinError(bin, ErrorInRatio(HistRecoMult_rebin[nmult][imc]->GetBinContent(bin),HistRecoMult_rebin[nmult][imc]->GetBinError(bin),
                        HistGenMult_rebin[nmult][imc]->GetBinContent(bin),HistGenMult_rebin[nmult][imc]->GetBinError(bin)));
                    }
                }
            }//end loop over perc
        } // end loop over mc files


        for (int nmult = 0; nmult < fpercbinnumb; nmult++) { //loop pver perc selections

            double multmin, multmax, eemin, eemax;
            //
            multmin = fpercentileMultlow[nmult];
            multmax = fpercentileMulthigh[nmult];
            eemin = fpercentileEnergylow[nmult];
            eemax = fpercentileEnergyhigh[nmult];

            //Weighted efficiency after IT0
            CorrEfficiencyIT1[nmult] = (TH1F*)CorrEfficiencyIT1parz[nmult][0]->Clone(Form("CorrEfficiencyIT1PercBin%i",nmult));
            CorrEfficiencyIT1[nmult]->Reset();
            for (int bin = 1; bin <= CorrEfficiency[nmult]->GetNbinsX(); bin ++){
                CorrEfficiencyIT1[nmult]->SetBinContent(bin,DoWeightedMean(0,
                        CorrEfficiencyIT1parz[nmult][0]->GetBinContent(bin), CorrEfficiencyIT1parz[nmult][1]->GetBinContent(bin), CorrEfficiencyIT1parz[nmult][2]->GetBinContent(bin),
                        CorrEfficiencyIT1parz[nmult][0]->GetBinError(bin), CorrEfficiencyIT1parz[nmult][1]->GetBinError(bin), CorrEfficiencyIT1parz[nmult][2]->GetBinError(bin),
                        w0,w1,w2));
                CorrEfficiencyIT1[nmult]->SetBinError(bin,DoWeightedMean(1,
                    CorrEfficiencyIT1parz[nmult][0]->GetBinContent(bin), CorrEfficiencyIT1parz[nmult][1]->GetBinContent(bin), CorrEfficiencyIT1parz[nmult][2]->GetBinContent(bin),
                    CorrEfficiencyIT1parz[nmult][0]->GetBinError(bin), CorrEfficiencyIT1parz[nmult][1]->GetBinError(bin), CorrEfficiencyIT1parz[nmult][2]->GetBinError(bin),
                    w0,w1,w2));
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
        cout << "\n";
        cout<<"--------------------------- Iteration 2 ------------------------------"<<endl;

        for(int imc = 0; imc < 3; imc ++){ //loop over mc files

            for(int nmult = 0; nmult < fpercbinnumb; nmult++){ //loop over perc selections

                double multmin, multmax, eemin, eemax;
                //
                multmin = fpercentileMultlow[nmult];
                multmax = fpercentileMulthigh[nmult];
                eemin = fpercentileEnergylow[nmult];
                eemax = fpercentileEnergyhigh[nmult];

                LevyMCIT2[nmult][imc] = LevyTsallis(Form("LevyMCIT2FitPercBin%i_MC%i",nmult,imc),lParticleMass);
                //,1e-3 , 1e+6, 1e-8, 1e+8,  5.,  0.1, 0.1);
                LevyMCIT2[nmult][imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

                //MC shape ---
                //
                Int_t fitres;
                Int_t trials = 0;
                trials = 0;
                do {
                    fitres =  HistGenMultCloneIT1[nmult][imc]->Fit(LevyMCIT2[nmult][imc], "0q","",minfit,maxfit);
                    Printf("Fit trial: %d", trials++);
                    if(trials > 10) {
                        Printf("FIT DOES NOT CONVERGE MC IN LINE %d",__LINE__);
                        break;
                    }
                }
                while (fitres != 0);

                if(imc==0){
                    //Define generated histos
                    hPtSpectrumGenIT2[nmult] = new TH1F(Form("hPtSpectrumGenIT2PercBin%i",nmult), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", nbins, minfit,maxfit);
                    hPtSpectrumGenIT2[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
                }

                //Define ratio shape
                hRatioPtShapeIT2[nmult][imc] = (TH1F*)hPtSpectrumGenIT2[nmult]->Clone(Form("hRatioPtShapeIT2PercBin%i_MC%i",nmult,imc));
                hRatioPtShapeIT2[nmult][imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
                hRatioPtShapeIT2[nmult][imc]->Reset();

                hPtSpectrumGenMCIT2[nmult][imc] = new TH1F(Form("hPtSpectrumGenMCIT2PercBin%i_MC%i",nmult,imc), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts",nbins, minfit,maxfit);
                hPtSpectrumGenMCIT2[nmult][imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            }// end loop over perc selections
        }


        for(int nmult = 0; nmult < fpercbinnumb; nmult++){ //loop over perc selections

            double multmin, multmax, eemin, eemax;
            //
            multmin = fpercentileMultlow[nmult];
            multmax = fpercentileMulthigh[nmult];
            eemin = fpercentileEnergylow[nmult];
            eemax = fpercentileEnergyhigh[nmult];

            //Fit function
            LevyIT2[nmult] = LevyTsallis(Form("LevyIT2FitPercBin%i",nmult),lParticleMass);
            LevyIT2[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

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

        }// end loop over perc selections

        //++++++++++++++++++++++++++++++++++++++++++++++
        // 0_IT2) Generate histos from Fit Levy-Tsallis
        //++++++++++++++++++++++++++++++++++++++++++++++

        for(int imc = 0; imc < 3; imc ++){ //loop over mc files

            // multiplicity spectra
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

                //MC spectra
                for (Int_t i = 0; i < n; i++) { //HERE
                    double xMC = (double)LevyMCIT2[nmult][imc]->GetRandom(minfit, maxfit);
                    //double xMC = (double)HistGenMultCloneIT1[nmult][imc]->GetRandom();
                    hPtSpectrumGenMCIT2[nmult][imc]->Fill(xMC);
                }
                // /dpT
                for (Int_t i = 0; i <= hPtSpectrumGenMCIT2[nmult][imc]->GetNbinsX(); i++) {
                    double entryMC = hPtSpectrumGenMCIT2[nmult][imc]->GetBinContent(i + 1) /
                                    hPtSpectrumGenMCIT2[nmult][imc]->GetBinWidth(i + 1);
                    hPtSpectrumGenMCIT2[nmult][imc]->SetBinContent(i + 1, entryMC);
                }
            }
        } // end loop over mc files

        // multiplicity spectra
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

            //Data
            for (Int_t i = 0; i < n; i++) {
                double x = (double)LevyIT2[nmult]->GetRandom(minfit,maxfit);
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
        for(int imc = 0; imc < 3; imc ++){ //loop over mc files

            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

                for (int bin = 1; bin <= hRatioPtShapeIT2[0][imc]->GetNbinsX(); bin ++){
                    if (hPtSpectrumGenMCIT2[nmult][imc]->GetBinContent(bin) != 0 ){
                        hRatioPtShapeIT2[nmult][imc]->SetBinContent(bin, hPtSpectrumGenIT2[nmult]->GetBinContent(bin)/hPtSpectrumGenMCIT2[nmult][imc]->GetBinContent(bin));
                        hRatioPtShapeIT2[nmult][imc]->SetBinError(bin, ErrorInRatio(hPtSpectrumGenIT2[nmult]->GetBinContent(bin),hPtSpectrumGenIT2[nmult]->GetBinError(bin),
                                        hPtSpectrumGenMCIT2[nmult][imc]->GetBinContent(bin),hPtSpectrumGenMCIT2[nmult][imc]->GetBinError(bin)));
                    }
                    else hRatioPtShapeIT2[nmult][imc]->SetBinContent(bin,1); //x*1=x
                }
            }

            //++++++++++++++++++++++++++++++++++++++++
            // 2_IT2) Re-weight reco and gen spectra
            //++++++++++++++++++++++++++++++++++++++++

            for (int nmult = 0; nmult < fpercbinnumb; nmult++) { //loop pver perc selections

                double multmin, multmax, eemin, eemax;
                //
                multmin = fpercentileMultlow[nmult];
                multmax = fpercentileMulthigh[nmult];
                eemin = fpercentileEnergylow[nmult];
                eemax = fpercentileEnergyhigh[nmult];

                HistRecoMult_rebin[nmult][imc]->Reset();
                HistGenMult_rebin[nmult][imc]->Reset();

                HistRecoMultIT2[nmult][imc] = (TH1F*) HistRecoMultIT1[nmult][imc]->Clone(Form("HistRecoIT2PercBin%i_MC%i",nmult,imc));
                HistRecoMultIT2[nmult][imc] ->Reset();
                HistGenMultIT2[nmult][imc] = (TH1F*) HistGenMultIT1[nmult][imc]->Clone(Form("HistGenIT2PercBin%i_MC%i",nmult,imc));
                HistGenMultIT2[nmult][imc] ->Reset();
                for (int bin = 1; bin <= HistRecoMultIT2[nmult][imc]->GetNbinsX(); bin ++){
                    //
                    HistRecoMultIT2[nmult][imc]-> SetBinContent(bin,HistRecoMultIT1[nmult][imc]->GetBinContent(bin)*hRatioPtShapeIT2[nmult][imc]->GetBinContent(bin));
                    HistRecoMultIT2[nmult][imc]-> SetBinError(bin,ErrorInRatio(HistRecoMultIT1[nmult][imc]->GetBinContent(bin), HistRecoMultIT1[nmult][imc]->GetBinError(bin),
                        hRatioPtShapeIT2[nmult][imc]->GetBinContent(bin), hRatioPtShapeIT2[nmult][imc]->GetBinError(bin))
                        );
                    //
                    HistGenMultIT2[nmult][imc]-> SetBinContent(bin, HistGenMultIT1[nmult][imc]->GetBinContent(bin)*hRatioPtShapeIT2[nmult][imc]->GetBinContent(bin));
                    HistGenMultIT2[nmult][imc]-> SetBinError(bin, ErrorInRatio(HistGenMultIT1[nmult][imc]->GetBinContent(bin), HistGenMultIT1[nmult][imc]->GetBinError(bin),
                        hRatioPtShapeIT2[nmult][imc]->GetBinContent(bin), hRatioPtShapeIT2[nmult][imc]->GetBinError(bin)
                        ));
                }

                Double_t tempptreco;
                for(long i = 1; i<HistRecoMultIT2[nmult][imc]->GetNbinsX()+1;i++){
                    tempptreco = HistRecoMultIT2[nmult][imc]->GetXaxis()->GetBinCenter(i);
                    for(long filling = 0; filling<HistRecoMultIT2[nmult][imc]->GetBinContent(i); filling++){
                        HistRecoMult_rebin[nmult][imc]->Fill(tempptreco);
                    }
                }
                Double_t tempptgen;
                for(long i = 1; i<HistGenMultIT2[nmult][imc]->GetNbinsX()+1;i++){
                    tempptgen = HistGenMultIT2[nmult][imc]->GetXaxis()->GetBinCenter(i);
                    for(long filling = 0; filling<HistGenMultIT2[nmult][imc]->GetBinContent(i); filling++){
                        HistGenMult_rebin[nmult][imc]->Fill(tempptgen);
                    }
                }

                //Clone for Fit
                HistGenMultCloneIT2[nmult][imc] = (TH1F*) HistGenMultIT2[nmult][imc]->Clone(Form("HistGenMultCloneIT2%i_MC%i", nmult,imc));
                HistGenMultCloneIT2[nmult][imc]->Reset();
                //
                for (int ibin = 1; ibin <= HistGenMultIT2[nmult][imc]->GetNbinsX(); ibin++){
                    HistGenMultCloneIT2[nmult][imc]->SetBinContent(ibin, HistGenMultIT2[nmult][imc]->GetBinContent(ibin)/HistGenMultIT2[nmult][imc]->GetBinWidth(ibin));
                    HistGenMultCloneIT2[nmult][imc]->SetBinError(ibin, HistGenMultIT2[nmult][imc]->GetBinError(ibin)/HistGenMultIT2[nmult][imc]->GetBinWidth(ibin));
                }
                double integral = HistGenMultCloneIT2[nmult][imc]->Integral(0,-1);
                HistGenMultCloneIT2[nmult][imc]->Scale(1./integral);
                //

                //++++++++++++++++++++++++++++++++++++++++
                // 3) Correct Efficiency
                //++++++++++++++++++++++++++++++++++++++++

                //
                CorrEfficiencyIT2parz[nmult][imc] = (TH1F*)HistRecoMult_rebin[nmult][imc]->Clone(Form("hCorrEfficiencyIT2parzPercBin%i_MC%i",nmult,imc));
                CorrEfficiencyIT2parz[nmult][imc]->SetTitle(Form("%s_%.0f_%.0f_%s_%03.0f_%03.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
                CorrEfficiencyIT2parz[nmult][imc]->Reset();

                for (int bin = 1; bin <= HistRecoMult_rebin[nmult][imc]->GetNbinsX(); bin ++){
                    if (HistGenMult_rebin[nmult][imc]->GetBinContent(bin) != 0){
                        CorrEfficiencyIT2parz[nmult][imc]->SetBinContent(bin , HistRecoMult_rebin[nmult][imc]->GetBinContent(bin)/HistGenMult_rebin[nmult][imc]->GetBinContent(bin));
                        CorrEfficiencyIT2parz[nmult][imc]->SetBinError(bin, ErrorInRatio(HistRecoMult_rebin[nmult][imc]->GetBinContent(bin),HistRecoMult_rebin[nmult][imc]->GetBinError(bin),
                        HistGenMult_rebin[nmult][imc]->GetBinContent(bin),HistGenMult_rebin[nmult][imc]->GetBinError(bin)));
                    }
                }
            }//end loop over perc
        } // end loop over mc files


        for (int nmult = 0; nmult < fpercbinnumb; nmult++) { //loop pver perc selections

            double multmin, multmax, eemin, eemax;
            //
            multmin = fpercentileMultlow[nmult];
            multmax = fpercentileMulthigh[nmult];
            eemin = fpercentileEnergylow[nmult];
            eemax = fpercentileEnergyhigh[nmult];

            //Weighted efficiency after IT0
            CorrEfficiencyIT2[nmult] = (TH1F*)CorrEfficiencyIT2parz[nmult][0]->Clone(Form("CorrEfficiencyIT2PercBin%i",nmult));
            CorrEfficiencyIT2[nmult]->Reset();
            for (int bin = 1; bin <= CorrEfficiencyIT2[nmult]->GetNbinsX(); bin ++){
                CorrEfficiencyIT2[nmult]->SetBinContent(bin,DoWeightedMean(0,
                        CorrEfficiencyIT2parz[nmult][0]->GetBinContent(bin), CorrEfficiencyIT2parz[nmult][1]->GetBinContent(bin), CorrEfficiencyIT2parz[nmult][2]->GetBinContent(bin),
                        CorrEfficiencyIT2parz[nmult][0]->GetBinError(bin), CorrEfficiencyIT2parz[nmult][1]->GetBinError(bin), CorrEfficiencyIT2parz[nmult][2]->GetBinError(bin),
                        w0,w1,w2));
                CorrEfficiencyIT2[nmult]->SetBinError(bin,DoWeightedMean(1,
                    CorrEfficiencyIT2parz[nmult][0]->GetBinContent(bin), CorrEfficiencyIT2parz[nmult][1]->GetBinContent(bin), CorrEfficiencyIT2parz[nmult][2]->GetBinContent(bin),
                    CorrEfficiencyIT2parz[nmult][0]->GetBinError(bin), CorrEfficiencyIT2parz[nmult][1]->GetBinError(bin), CorrEfficiencyIT2parz[nmult][2]->GetBinError(bin),
                    w0,w1,w2));
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
            hspectraIT3[nmult]->SetName(Form("PtSpectrumIT2PercBin%i",nmult));
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

            //Final spectra
            double multmin, multmax, eemin, eemax;
            //
            multmin = fpercentileMultlow[nmult];
            multmax = fpercentileMulthigh[nmult];
            eemin = fpercentileEnergylow[nmult];
            eemax = fpercentileEnergyhigh[nmult];

            //
            hFinalSpectraStat[nmult] = (TH1F*)hspectraIT3[nmult]->Clone();
            hFinalSpectraSyst[nmult] = (TH1F*)hspectraIT3[nmult]->Clone();
            hFinalSpectraSyst[nmult]->Reset();
            hFinalSpectraStat[nmult]->SetName(Form("FinalPtSpectrumStat_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            hFinalSpectraSyst[nmult]->SetName(Form("FinalPtSpectrumSyst_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        }

        cout<<"------------------------------ Done! ----------------------------------"<<endl;
        cout << "\n";

        cout << "\n\n";

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
        cout<<" Event Loss Corr File ..........: "<< fMCFileEventLoss.Data() <<endl;
        cout << "\n";
        cout<<" !!!!!!!!!!! Please make sure event loss percentile bins match your analisis!!! No check here !!!!!!!!!!!!!!!"<<endl;
        cout << "\n";

        filenorm = TFile::Open(fMCFileEventLoss);
        //
        heventloss = (TH1F*)filenorm->Get("EventLoss/hevtloss");

        if (fDoMCShape){
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
                hspectraEvtLoss[nmult] = (TH1F *)hspectraIT3[nmult]->Clone(Form("PtSpectrumEvLossCorr_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), fpercentileMultlow[nmult], fpercentileMulthigh[nmult], fWhichEnergyEstimator.Data(), fpercentileEnergylow[nmult], fpercentileEnergyhigh[nmult]));
            }
        } else {
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
                hspectraEvtLoss[nmult] = (TH1F *)hspectra[nmult]->Clone(Form("PtSpectrumEvLossCorr_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), fpercentileMultlow[nmult], fpercentileMulthigh[nmult], fWhichEnergyEstimator.Data(), fpercentileEnergylow[nmult], fpercentileEnergyhigh[nmult]));
            }
        }

        //Apply correction
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

            int binevtloss = -1;
            binevtloss = heventloss->GetXaxis()->FindBin(Form("%s_%.0f-%.0f_%s_%.0f-%.0f", fWhichMultEstimator.Data(), fpercentileMultlow[nmult], fpercentileMulthigh[nmult], fWhichEnergyEstimator.Data(), fpercentileEnergylow[nmult], fpercentileEnergyhigh[nmult]));
            if (binevtloss == -1) {
                cout << "ERROR: Event loss bin not found! Check your percentile bins!" << endl;
                return;
            }

            for (int bin = 1 ; bin <=  hspectraEvtLoss[nmult]->GetNbinsX(); bin ++ ){
                //
                double content = hspectraEvtLoss[nmult]->GetBinContent(bin) * heventloss->GetBinContent(binevtloss);
                double error = ErrorInRatio(hspectraEvtLoss[nmult]->GetBinContent(bin), hspectraEvtLoss[nmult]->GetBinError(bin), heventloss->GetBinContent(binevtloss), heventloss->GetBinError(binevtloss));
                hspectraEvtLoss[nmult]->SetBinContent(bin, content);
                hspectraEvtLoss[nmult]->SetBinError(bin,error);
            }

            //Final spectra
            double multmin, multmax, eemin, eemax;
            //
            multmin = fpercentileMultlow[nmult];
            multmax = fpercentileMulthigh[nmult];
            eemin = fpercentileEnergylow[nmult];
            eemax = fpercentileEnergyhigh[nmult];

            hFinalSpectraStat[nmult] = 0x0;
            hFinalSpectraSyst[nmult] = 0x0;
            hFinalSpectraStat[nmult] = (TH1F*)hspectraEvtLoss[nmult]->Clone();
            hFinalSpectraSyst[nmult] = (TH1F*)hspectraEvtLoss[nmult]->Clone();
            hFinalSpectraSyst[nmult]->Reset();
            hFinalSpectraStat[nmult]->SetName(Form("FinalPtSpectrumStat_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            hFinalSpectraSyst[nmult]->SetName(Form("FinalPtSpectrumSyst_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        }

        cout<<"------------------------------ Done! ----------------------------------"<<endl;
        cout << "\n";
        cout<<"------------- Event Loss correction applied succesfully! --------------"<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl;
    }

    //===========================================================================================
    // ============================== END Event Loss correction =================================
    //===========================================================================================


    //===========================================================================================
    //=============================== Signal Loss correction ====================================
    //===========================================================================================

    TH1D* hsgnloss[fpercbinnumb]; //these are TH2D!! As from definition in the AliPhyiscs task
    TH1F* hspectraSgnLoss[fpercbinnumb];
    TFile* filenorm2;

    //se fai la MC pt fagli prendere gli spettri suoi, altrimenti fallo partire da spettro 1
    if (fDoSignalLoss){

        cout << "\n\n";
        cout<<"-----------------------------------------------------------------------"<<endl;
        cout<<"-------------------- Doing Signal Loss correction ---------------------"<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl;
        cout << "\n";
        cout<<" Signal Loss Corr File ..........: "<< fMCFileSignalLoss.Data() <<endl;
        cout << "\n\n";

        filenorm2 = TFile::Open(fMCFileSignalLoss);

        if (fDoEventLoss){
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
                double multmin, multmax, eemin, eemax;
                //
                multmin = fpercentileMultlow[nmult];
                multmax = fpercentileMulthigh[nmult];
                eemin = fpercentileEnergylow[nmult];
                eemax = fpercentileEnergyhigh[nmult];
                hspectraSgnLoss[nmult] = (TH1F*) hspectraEvtLoss[nmult]->Clone(Form("PtSpectrumSgnLossCorr_%s_%.0f_%.0f_%s_%.0f_%.0f",fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            }
        } else if (fDoMCShape){
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
                //
                double multmin, multmax, eemin, eemax;
                multmin = fpercentileMultlow[nmult];
                multmax = fpercentileMulthigh[nmult];
                eemin = fpercentileEnergylow[nmult];
                eemax = fpercentileEnergyhigh[nmult];
                hspectraSgnLoss[nmult] = (TH1F*) hspectraIT3[nmult]->Clone(Form("PtSpectrumSgnLossCorr_%s_%.0f_%.0f_%s_%.0f_%.0f",fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            }
        } else {
            for (int nmult = 0; nmult < fpercbinnumb; nmult++) {
                //
                double multmin, multmax, eemin, eemax;
                multmin = fpercentileMultlow[nmult];
                multmax = fpercentileMulthigh[nmult];
                eemin = fpercentileEnergylow[nmult];
                eemax = fpercentileEnergyhigh[nmult];
                hspectraSgnLoss[nmult] = (TH1F*) hspectra[nmult]->Clone(Form("PtSpectrumSgnLossCorr_%s_%.0f_%.0f_%s_%.0f_%.0f",fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            }
        }

        //Apply correction
        cout << " List of applied correction factors: " << endl;
        for (int nmult = 0; nmult < fpercbinnumb; nmult++) {

            double multmin, multmax, eemin, eemax;
            //
            multmin = fpercentileMultlow[nmult];
            multmax = fpercentileMulthigh[nmult];
            eemin = fpercentileEnergylow[nmult];
            eemax = fpercentileEnergyhigh[nmult];

            hsgnloss[nmult] = (TH1D *)filenorm2->Get(Form("hsgnloss_%s_%.0f-%.0f_%s_%.0f-%.0f", fWhichMultEstimator.Data(), fpercentileMultlow[nmult], fpercentileMulthigh[nmult], fWhichEnergyEstimator.Data(), fpercentileEnergylow[nmult], fpercentileEnergyhigh[nmult]));
            int start = hsgnloss[nmult]->FindBin(minfit) -1;
            //
            for (int bin = 1 ; bin <=  hspectraSgnLoss[nmult]->GetNbinsX(); bin ++ ){
                //
                if(hsgnloss[nmult]->GetBinContent(start + bin)!=0){
                    double content = hspectraSgnLoss[nmult]->GetBinContent(bin)/hsgnloss[nmult]->GetBinContent(start + bin);
                    double error = ErrorInRatio(hspectraSgnLoss[nmult]->GetBinContent(bin), hspectraSgnLoss[nmult]->GetBinError(bin), hsgnloss[nmult]->GetBinContent(start + bin), hsgnloss[nmult]->GetBinError(start + bin));
                    hspectraSgnLoss[nmult]->SetBinContent(bin, content);
                    hspectraSgnLoss[nmult]->SetBinError(bin,error);
                } else {
                    double content = hspectraSgnLoss[nmult]->GetBinContent(bin);
                    double error = hspectraSgnLoss[nmult]->GetBinError(bin);
                    hspectraSgnLoss[nmult]->SetBinContent(bin, content);
                    hspectraSgnLoss[nmult]->SetBinError(bin,error);
                }
            }

            //Final spectra
            hFinalSpectraStat[nmult] = 0x0;
            hFinalSpectraSyst[nmult] = 0x0;
            hFinalSpectraStat[nmult] = (TH1F*)hspectraSgnLoss[nmult]->Clone();
            hFinalSpectraSyst[nmult] = (TH1F*)hspectraSgnLoss[nmult]->Clone();
            hFinalSpectraSyst[nmult]->Reset();
            hFinalSpectraStat[nmult]->SetName(Form("FinalPtSpectrumStat_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            hFinalSpectraSyst[nmult]->SetName(Form("FinalPtSpectrumSyst_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
        }

        cout<<"------------------------------ Done! ----------------------------------"<<endl;
        cout << "\n";
        cout<<"------------ Signal Loss correction applied succesfully! --------------"<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl;
    }

    //===========================================================================================
    //============================== END Signal Loss correction =================================
    //===========================================================================================

    //===========================================================================================
    //=============================== Systematics Uncertainties =================================
    //===========================================================================================

    TH1F *hspectraStat[fpercbinnumb];
    TH1F *hspectraSystTot[fpercbinnumb];
    TH1F *hspectraTotErrors[fpercbinnumb];
    TH1F *hSysSgnLoss[fpercbinnumb];
    TFile *fsystCutVar, *fsystOOB, *fsystInBunch, *fsystMatBudget, *fsystSgnLoss;
    TH1F *hSystTop, *hSystSel, *hSystSgnLoss, *hTotalSyst, *hMatBudget, *hInBunchPileUp, *hpTindip, *hMatBudgetK0s;
    TH1D *hSystOOB, *hInBunch040, *hInBunch40100;
    TF1 *funMatBudget;
    TCanvas *cSyst;
    TLegend *leg;
    TLatex* text;

    TString folder_ = "V0Analysis";
    if (fParticle.Contains("Minus")) folder_ = "CascadeAnalysis";

    if (fDoSystematics){

        cout << "\n\n";
        cout<<"-----------------------------------------------------------------------"<<endl;
        cout<<"--------------------- Applying Systematic Errors ----------------------"<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl;
        cout << "\n\n";

        fsystCutVar = TFile::Open(Form("%s/Systematics/Systematics-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root", folder_.Data(), fWhichParticle.Data(), 0., 100., 0., 100.));
        fsystOOB = TFile::Open(Form("%s/Systematics/OOBSyst%s.root", folder_.Data(), fWhichParticle.Data()));
        fsystInBunch = TFile::Open(Form("%s/Systematics/SystematicsPileUp-%s.root", folder_.Data(), fWhichParticle.Data()));
        fsystMatBudget = TFile::Open(Form("%s/Systematics/systMatBudget-%s.root", folder_.Data(), fWhichParticle.Data()));
        /*if (fpercbinnumb > 1) {
            fsystSgnLoss = TFile::Open(Form("%s/Systematics/SgnLoss/SgnLossSystematics-%s_%s_Fixed%sin%03.0f_%03.0f.root", folder_.Data(), fWhichParticle.Data(), fWhichVarEstimator.Data(), fWhichFixedEstimator.Data(), fLowFixed, fHighFixed));
        } else {
            fsystSgnLoss = TFile::Open(Form("%s/Systematics/SgnLoss/SgnLossSystematics-%s-13TeV_INELgt0.root", folder_.Data(), fWhichParticle.Data()));
        }*/
        hSystTop = (TH1F *)fsystCutVar->Get("hSystTopological");
        hSystSel = (TH1F *)fsystCutVar->Get("hSystOthers");
        hSystOOB = (TH1D *)fsystOOB->Get("hMaxDev");
        if (!fWhichParticle.Contains("Lambda")) hSystOOB->Scale(1. / 2);
        if (fWhichParticle.Contains("K0Short")){
            hInBunch040 = (TH1D *)fsystInBunch->Get("hsystPileUp");
        } else {
            hInBunch040 = (TH1D *)fsystInBunch->Get(Form("fHistPileUp%s_0_40", fWhichParticle.Data()));
            hInBunch40100 = (TH1D *)fsystInBunch->Get(Form("fHistPileUp%s_40_100", fWhichParticle.Data()));
        }
        if (fWhichParticle.Contains("K0Short")){
            hMatBudgetK0s = (TH1F *)fsystMatBudget->Get("hMatBudget");
        } else {
            funMatBudget = (TF1 *)fsystMatBudget->Get("funcFit");
        }

        hTotalSyst = (TH1F *)hFinalSpectraStat[0]->Clone("hTotalSyst");
        hTotalSyst->Reset();
        hInBunchPileUp = (TH1F *)hFinalSpectraStat[0]->Clone("hInBunchPileUp");
        hInBunchPileUp->Reset();
        hMatBudget = (TH1F *)hFinalSpectraStat[0]->Clone("hMatBudget");
        hMatBudget->Reset();
        hpTindip = (TH1F *)hFinalSpectraStat[0]->Clone("hpTindip");
        hpTindip->Reset();
        double top, sel, oob, ibp, matbud;

        for (int bin = 1 ; bin <=  hFinalSpectraStat[0]->GetNbinsX(); bin ++ ){
            top = hSystTop->GetBinContent(bin);
            sel = hSystSel->GetBinContent(bin);
            oob = hSystOOB->GetBinContent(bin);
            if (fWhichParticle.Contains("K0Short")){
                ibp = hInBunch040->GetBinContent(bin);
            } else {
                ibp = hInBunch040->GetBinError(bin) * 40. / 100 + hInBunch40100->GetBinError(bin) * 60. / 100;
            }
            hInBunchPileUp->SetBinContent(bin, ibp);
            if (fWhichParticle.Contains("K0Short")){
                matbud = w0 / (w0 + w1 + w2) * hMatBudgetK0s->GetBinContent(bin);
                hMatBudget->SetBinContent(bin, matbud);
            } else {
                matbud = w0 / (w0 + w1 + w2) * TMath::Abs((1 - funMatBudget->Eval(hFinalSpectraStat[0]->GetBinCenter(bin)))) / 2;
                hMatBudget->SetBinContent(bin, matbud);
            }
            hpTindip->SetBinContent(bin, 0.02);

            hTotalSyst->SetBinContent(bin, TMath::Sqrt( top*top  +
                                                        sel*sel +
                                                        0.02*0.02 +
                                                        oob*oob +
                                                        ibp*ibp +
                                                        matbud*matbud
                                                        ));
        }

        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            double multmin, multmax, eemin, eemax;
            //
            multmin = fpercentileMultlow[nmult];
            multmax = fpercentileMulthigh[nmult];
            eemin = fpercentileEnergylow[nmult];
            eemax = fpercentileEnergyhigh[nmult];

            // sgn loss syst
            //hSysSgnLoss[nmult] = (TH1F *)fsystSgnLoss->Get(Form("var_%s_%.0f-%0.f_%s_%.0f-%0.f", fWhichMultEstimator.Data(), fpercentileMultlow[nmult], fpercentileMulthigh[nmult], fWhichEnergyEstimator.Data(), fpercentileEnergylow[nmult], fpercentileEnergyhigh[nmult]));

            hspectraStat[nmult] = (TH1F*)hFinalSpectraStat[nmult]->Clone(Form("PtSpectrumCorrStat_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            hspectraSystTot[nmult] = (TH1F *)hFinalSpectraStat[nmult]->Clone(Form("PtSpectrumCorrSystTot_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

            for (int bin = 1 ; bin <=  hspectraStat[nmult]->GetNbinsX(); bin ++ ){
                double sgnloss = 0.;
                /*if (hSysSgnLoss[nmult]->FindBin(hspectraStat[nmult]->GetBinCenter(bin)) == 0){
                    sgnloss = hSysSgnLoss[nmult]->GetBinContent(1);
                } else {
                    sgnloss = hSysSgnLoss[nmult]->GetBinContent(hSysSgnLoss[nmult]->FindBin(hspectraStat[nmult]->GetBinCenter(bin)));
                }*/
                top = hSystTop->GetBinContent(bin);
                sel = hSystSel->GetBinContent(bin);
                oob = hSystOOB->GetBinContent(bin);
                if (fWhichParticle.Contains("K0Short")){
                    ibp = hInBunch040->GetBinContent(bin);
                } else {
                    ibp = hInBunch040->GetBinError(bin) * 40. / 100 + hInBunch40100->GetBinError(bin) * 60. / 100;
                }
                if (fWhichParticle.Contains("K0Short"))
                {
                    matbud = w0 / (w0 + w1 + w2) * hMatBudgetK0s->GetBinContent(bin);
                } else {
                    matbud = w0 / (w0 + w1 + w2) * TMath::Abs((1 - funMatBudget->Eval(hspectraSgnLoss[0]->GetBinCenter(bin)))) / 2;
                }
                hspectraSystTot[nmult]->SetBinError(bin, hspectraSystTot[nmult]->GetBinContent(bin) *
                                                                TMath::Sqrt(
                                                                    top * top +
                                                                    sel * sel +
                                                                    0.02 * 0.02 +
                                                                    oob * oob +
                                                                    ibp * ibp +
                                                                    matbud * matbud +
                                                                    sgnloss * sgnloss));
            }

            // Final spectra
            hFinalSpectraSyst[nmult] = 0x0;
            hFinalSpectraSyst[nmult] = (TH1F *)hspectraSystTot[nmult]->Clone();

            // Stat + syst
            hspectraTotErrors[nmult] = (TH1F *)hspectraSystTot[nmult]->Clone(Form("PtSpectrumCorrTotErrors_%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));
            for (int bin = 1; bin <= hspectraTotErrors[nmult]->GetNbinsX(); bin++)
            {
                hspectraTotErrors[nmult]->SetBinError(bin, TMath::Sqrt(hspectraStat[nmult]->GetBinError(bin) * hspectraStat[nmult]->GetBinError(bin) +
                                                                       hspectraSystTot[nmult]->GetBinError(bin) * hspectraSystTot[nmult]->GetBinError(bin)));
            }

            // Fit function
            LevyFit[nmult] = LevyTsallis(Form("LevyFit%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax), lParticleMass);
            LevyFit[nmult]->SetTitle(Form("%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), multmin, multmax, fWhichEnergyEstimator.Data(), eemin, eemax));

            // Fit spectra with syst + stat errors
            Int_t fitres;
            Int_t trials = 0;
            trials = 0;
            do
            {
                fitres = hspectraTotErrors[nmult]->Fit(LevyFit[nmult], "0q", "", minfit, maxfit);
                Printf("Fit trial: %d", trials++);
                if (trials > 10)
                {
                    Printf("FIT DOES NOT CONVERGE IN LINE %d", __LINE__);
                    break;
                }
            } while (fitres != 0);
        }

        cSyst = new TCanvas("cSyst", "cSyst", 800, 600);
        cSyst->cd();
        cSyst->SetBottomMargin(0.15);
        cSyst->SetLeftMargin(0.15);
        cSyst->SetRightMargin(0.05);
        cSyst->SetTopMargin(0.05);
        cSyst->SetTicks(1,1);
        hTotalSyst->SetLineColor(kBlack);
        hTotalSyst->SetLineWidth(2);
        hSystTop->SetLineColor(kRed);
        hSystTop->SetLineWidth(2);
        hSystSel->SetLineColor(kBlue);
        hSystSel->SetLineWidth(2);
        hSystOOB->SetLineColor(kGreen);
        hSystOOB->SetLineWidth(2);
        hInBunchPileUp->SetLineColor(kMagenta);
        hInBunchPileUp->SetLineWidth(2);
        hMatBudget->SetLineColor(kAzure+7);
        hMatBudget->SetLineWidth(2);
        hpTindip->SetLineColor(kOrange-3);
        hpTindip->SetLineWidth(2);
        hTotalSyst->GetYaxis()->SetRangeUser(0, hTotalSyst->GetMaximum()*1.5);
        hTotalSyst->GetYaxis()->SetTitle("Relative uncertainty");
        if (fWhichParticle.Contains("K0Short")) hTotalSyst->GetXaxis()->SetRangeUser(0.,10.);
        if (fWhichParticle.Contains("Lambda")) hTotalSyst->GetXaxis()->SetRangeUser(0.4,8.);
        if (fWhichParticle.Contains("Lambda")) hTotalSyst->GetYaxis()->SetRangeUser(0.,0.15);
        hTotalSyst->SetStats(0);
        hTotalSyst->SetTitle("");
        hTotalSyst->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hTotalSyst->GetXaxis()->SetTitleSize(0.05);
        hTotalSyst->GetXaxis()->SetTitleOffset(0.9);
        hTotalSyst->GetXaxis()->SetLabelSize(0.04);
        hTotalSyst->GetYaxis()->SetTitleSize(0.05);
        hTotalSyst->GetYaxis()->SetTitleOffset(0.9);
        hTotalSyst->GetYaxis()->SetLabelSize(0.04);
        hTotalSyst->Draw("hist");
        hSystTop->Draw("hist same");
        hSystSel->Draw("hist same");
        hSystOOB->Draw("hist same");
        hInBunchPileUp->Draw("hist same");
        hMatBudget->Draw("hist same");
        hpTindip->Draw("hist same");

        leg = new TLegend(0.5, 0.6, 0.88, 0.88);
        leg->AddEntry(hTotalSyst, "Total", "l");
        leg->AddEntry(hSystTop, "Topological", "l");
        leg->AddEntry(hSystSel, "Other selections", "l");
        leg->AddEntry(hSystOOB, "Out-of-bunch pile-up", "l");
        leg->AddEntry(hInBunchPileUp, "In-bunch pile-up", "l");
        leg->AddEntry(hMatBudget, "Material budget", "l");
        leg->AddEntry(hpTindip, "p_{T} independent", "l");
        leg->SetBorderSize(0);
        leg->Draw();

        text = new TLatex();
        text->SetNDC();
        text->SetTextSize(0.08);
        if (fWhichParticle.Contains("Xi")) text->DrawLatex(0.25, 0.82, Form("#Xi"));
        if (fWhichParticle.Contains("Lambda")) text->DrawLatex(0.25, 0.82, Form("#Lambda"));
        if (fWhichParticle.Contains("Omega")) text->DrawLatex(0.25, 0.82, Form("#Omega"));
        if (fWhichParticle.Contains("K0S")) text->DrawLatex(0.25, 0.82, Form("K^{0}_{S}"));

        cout<<"------------------------------ Done! ----------------------------------"<<endl;
        cout << "\n";
        cout<<"-------------------- End systematics computations  --------------------"<<endl;
        cout<<"-----------------------------------------------------------------------"<<endl;
    }

    //===========================================================================================
    //=========================== END Systematics Uncertainties =================================
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

    TDirectoryFile *UncorrSpectra = new TDirectoryFile("RecoEffCorrSpectra","Spectra after reco efficiency correction");
    UncorrSpectra->cd();
    for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
        hspectra[nmult]->Write();
        if (fDoMCShape) Levy[nmult]->Write();
    }
    lResultsFile->cd();

    TDirectoryFile *MCpTshape;
    TDirectoryFile *Iteration[4];
    TDirectoryFile *MCshape[4];
    TDirectoryFile *Ratios[4];
    TDirectoryFile *AfterIt[4];
    TDirectoryFile *Gen[4];
    if (fDoMCShape){
        MCpTshape = new TDirectoryFile("MCpTshape","MC pT shape correction");
        MCpTshape->cd();
        Iteration[0] = new TDirectoryFile("Iteration0","MC pT shape correction");
        Iteration[0]->cd();
        MCshape[0] = new TDirectoryFile("MCshape","MC pT shape correction");
        MCshape[0]->cd();
        for(int imc = 0; imc<3; imc++){
            HistGenClone[imc]->Write();
            LevyMC[imc]->Write();
        }
        for(int imc = 0; imc<3; imc++){
            EfficiencyMB[imc]->Write();
        }
        EfficiencyIT0->Write();
        Iteration[0]->cd();
        Ratios[0] = new TDirectoryFile("RatioPtShapeAndCorrFactors","Ratio Data/MC pT shape");
        Ratios[0]->cd();
        for(int imc = 0; imc<3; imc++){
            for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
                hRatioPtShape[nmult][imc]->Write();
            }
        }
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            RatioEff[nmult]->Write();
            for(int imc = 0; imc<3; imc++){
                RatioEffparz[nmult][imc]->Write();
            }
        }
        Iteration[0]->cd();
        Gen[0] = new TDirectoryFile("GeneratedSpectra","MC pT shape correction");
        Gen[0]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            hPtSpectrumGen[nmult]->Write();
        }
        for(int imc = 0; imc<3; imc++){
            hPtSpectrumGenMCMB[imc]->Write();
        }
        Iteration[0]->cd();
        AfterIt[0] = new TDirectoryFile("SpectraAfterIt0","MC pT shape correction");
        AfterIt[0]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            hspectraIT1[nmult]->Write();
          //  LevyIT1[nmult]->Write(); //to do mettere i fit
        }
        Iteration[0]->cd();
        MCpTshape->cd();
        //
        Iteration[1] = new TDirectoryFile("Iteration1","MC pT shape correction");
        Iteration[1]->cd();
        MCshape[1] = new TDirectoryFile("MCshape","MC pT shape correction");
        MCshape[1]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            for(int imc = 0; imc<3; imc++){
                LevyMCIT1[nmult][imc]->Write();
                HistGenMultClone[nmult][imc]->Write();
            }
        }
        Iteration[1]->cd();
        Ratios[1] = new TDirectoryFile("RatioPtShapeAndCorrFactors","Ratio Data/MC pT shape");
        Ratios[1]->cd();
        for(int imc = 0; imc<3; imc++){
            for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
                hRatioPtShapeIT1[nmult][imc]->Write();
            }
        }
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            RatioEffIT1[nmult]->Write();
        }
        Iteration[1]->cd();
        Gen[1] = new TDirectoryFile("GeneratedSpectra","MC pT shape correction");
        Gen[1]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            hPtSpectrumGenIT1[nmult]->Write();
            for(int imc = 0; imc<3; imc++){
                hPtSpectrumGenMCIT1[nmult][imc]->Write();
            }
        }
        Iteration[1]->cd();
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
            for(int imc = 0; imc<3; imc++){
                LevyMCIT2[nmult][imc]->Write();
                HistGenMultIT1[nmult][imc]->Write();
            }
        }
        Iteration[2]->cd();
        Ratios[2] = new TDirectoryFile("RatioPtShapeAndCorrFactors","Ratio Data/MC pT shape");
        Ratios[2]->cd();
        for(int imc = 0; imc<3; imc++){
            for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
                hRatioPtShapeIT2[nmult][imc]->Write();
            }
        }
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            RatioEffIT2[nmult]->Write();
        }
        Iteration[2]->cd();
        Gen[2] = new TDirectoryFile("GeneratedSpectra","MC pT shape correction");
        Gen[2]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            hPtSpectrumGenIT2[nmult]->Write();
            for(int imc = 0; imc<3; imc++){
                hPtSpectrumGenMCIT2[nmult][imc]->Write();
            }
        }
        Iteration[2]->cd();
        AfterIt[2] = new TDirectoryFile("SpectraAfterIt2","MC pT shape correction");
        AfterIt[2]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            hspectraIT3[nmult]->Write();
        }
        Iteration[2]->cd();
        MCpTshape->cd();
        //
        Iteration[3] = new TDirectoryFile("Iteration3","MC pT shape correction");
        Iteration[3]->cd();
        MCshape[3] = new TDirectoryFile("MCshape","MC pT shape correction");
        MCshape[3]->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            for(int imc = 0; imc<3; imc++){
                HistGenMultIT2[nmult][imc]->Write();
            }
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
            hsgnloss[nmult]->Write();
        }
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            hspectraSgnLoss[nmult]->Write();
        }
    }
    NormCorr->cd();
    lResultsFile->cd();
    if (fDoSystematics){
        TDirectoryFile *Systematics;
        Systematics = new TDirectoryFile("SpectraWtSystematics","Systematics");
        Systematics->cd();
        for(int nmult = 0; nmult < fpercbinnumb; nmult++) {
            hspectraStat[nmult]->Write();
            hspectraSystTot[nmult]->Write();
            hspectraTotErrors[nmult]->Write();
            LevyFit[nmult]->Write();
        }
        Systematics->cd();
        hTotalSyst->Write();
        cSyst->Write();
        cSyst->SaveAs(Form("%s/Systematics/images/spectra/TotSystematics%s.pdf", folder_.Data(), fWhichParticle.Data()));
        cSyst->SaveAs(Form("%s/Systematics/images/spectra/TotSystematics%s.png", folder_.Data(), fWhichParticle.Data()));
    }
    lResultsFile->cd();
    TDirectoryFile *Final;
    Final = new TDirectoryFile("FinalSpectra", "Final");
    Final->cd();
    for (int nmult = 0; nmult < fpercbinnumb; nmult++)
    {
        hFinalSpectraStat[nmult]->Write();
        hFinalSpectraSyst[nmult]->Write();
        LevyFit[nmult]->Write();
        if (fDoSystematics) {
            hspectraTotErrors[nmult]->Write();
        }
    }
    Final->cd();

    //Exit Batch Mode
    gROOT->SetBatch(kFALSE);

    cout<<"--------------------------------------------------------"<<endl;
    cout<<" There, done! "<<endl;
    cout<<endl;
    cout << "\n\n";
}
