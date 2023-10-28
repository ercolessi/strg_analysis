/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/***********************************************

  Lambda Analysis Module
  ----------------------

This Class enables the analysis of output files
created with AliAnalysisTaskExtractV0 and
AliAnalysisTaskExtractPerformanceV0 grid tasks.

It constructs corrected Lambda, AntiLambda and
K0Short spectra.

This version: 27th April 2012

To be compiled as a library
(.L AliV0Module.cxx++ or something to that effect)

--- David Dobrigkeit Chinellato
    daviddc@ifi.unicamp.br

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
#include "TMinuit.h"
#include "TGraph.h"
using namespace std;

//--- For ROOT ---
#include "AliV0Module.h"
//#include "AliVTrack.h"

AliV0Module::AliV0Module()
{
    //Dummy Constructor that sets some default values
    //so that nothing too horrible will happen (hopefully)
    fWhichParticle           = "Lambda"; //Default
    fCINT1BoverINELratio     = 0.852;
    fRapidityBoundary        = 0.5;
    fPosRap = kFALSE;
    fNegRap = kFALSE;
    fRealDataFile            = "";
    fMCDataFile0             = "";
    fMCDataFile1              = "";
    fMCDataFile2              = "";
    fFeedDownDataFile        = "";
    fOutputDataFile          = "";

    fTOFpercFilename = "";
    fpercFilenameN0815 = "";

    kDoEfficiency = kFALSE;
    fEventsWeights[0] = 0.;
    fEventsWeights[1] = 0.;
    fEventsWeights[2] = 0.;

    fWhichMultEstimator          = "V0M";
    fWhichEffEnergyEstimator     = "ZDC";


    //Selections
    fCutV0Radius                     = -1;
    fCutDCANegToPV                   = -1;
    fCutDCAPosToPV                   = -1;
    fCutDCAV0Daughters               = 1000;
    fCutV0CosPA                      = -2;
    fCutProperLifetime               = 1e+6;
    fCutTPCPIDNSigmas                = 1e+6;
    fCutNSigmasForSignalExtraction   = 5;
    fCutLeastNumberOfCrossedRows     = 70;
    fCutDaughterEta                  = 0.8;
    fCutCompetingV0Rejection         = -1;

    //Set Feeddown Treatment
    // - - fFDSwitch = "UseMCRatio";
    fFDSwitch = "NoFD";
    fListName = "WoV0refit";

    //Default: Use bin counting
    fFitBackgroundSwitch = kFALSE;
    // switch-on / off G3/Fluka correction
    fG3FlukaCorrOn = kFALSE;

    //Default: Min-Bias
    fPerformMultiplicityStudy = kFALSE;
    fUseIntegratedEfficiencies = kTRUE;
    fLoMultBound       = -1;
    fHiMultBound       = 10000;

    //Default: no cut in Armenteros
    fSpecialArmenterosCutK0s = kFALSE;

    //Pt Bins: undefined
    fptbinnumb = -1;

    for( Int_t ipt = 0; ipt<200; ipt++) fptbinlimits[ipt] = 0;
    for( Int_t i=0; i<3; i++){ for(int j=0; j<3; j++) fMatrix[i][j]=0.;}
    //Primary selection criteria
    // --- IsPhysicalPrimary....: Usual method
    // --- DistToPV.............: Iouri, Luke
    fPrimarySelection = "IsPhysicalPrimary";
    fDistToPVCut = 0.001;

    ///// TF1
    fLevyFitXiFeedDown = 0x0;
    for(int k=0; k<2; k++) fLevyFitXiFeedDownErr[k] = 0x0;
    fPtLegProtCut  = 0.3;
    fRequestITSTOF = -1;
    fPtMinITSTOF = 1.;
}

AliV0Module::AliV0Module(TString fParticleType)
{
    // Allows definition of Particle Type in analysis.
    // Possible Options are "Lambda", "AntiLambda" and "K0Short".
    // If some other string is given, this constructor will
    // default to "Lambda".
    fWhichParticle = fParticleType;
    if(fWhichParticle!="Lambda"&&fWhichParticle!="AntiLambda"&&fWhichParticle!="K0Short") {
        cout<<"Particle Type "<<fParticleType<<" unknown. Set to lambda."<<endl;
        fWhichParticle = "Lambda";
    }
    fCINT1BoverINELratio     = 0.852;
    fRapidityBoundary        = 0.5;
    fPosRap = kFALSE;
    fNegRap = kFALSE;
    fRealDataFile            = "";
    fMCDataFile0             = "";
    fMCDataFile1              = "";
    fMCDataFile2              = "";
    fFeedDownDataFile        = "";
    fOutputDataFile          = "";

    fTOFpercFilename = "";
    fpercFilenameN0815 = "";

    kDoEfficiency = kFALSE;
    fEventsWeights[0] = 0.;
    fEventsWeights[1] = 0.;
    fEventsWeights[2] = 0.;

    fWhichMultEstimator          = "V0M";
    fWhichEffEnergyEstimator     = "ZDC";

    //Selections
    fCutV0Radius                     = -1;
    fCutDCANegToPV                   = -1;
    fCutDCAPosToPV                   = -1;
    fCutDCAV0Daughters               = 1000;
    fCutV0CosPA                      = -2;
    fCutProperLifetime               = 1e+6;
    fCutTPCPIDNSigmas                = 1e+6;
    fCutNSigmasForSignalExtraction   = 5;
    fCutLeastNumberOfCrossedRows     = 70;
    fCutDaughterEta                  = 0.8;
    fCutCompetingV0Rejection         = -1;

    //Set Feeddown Treatment
    // - - fFDSwitch = "UseMCRatio";
    fFDSwitch = "NoFD";
    fListName = "WoV0refit";

    //Default: Use bin counting
    fFitBackgroundSwitch = kFALSE;

    // switch-on / off G3/Fluka correction
    fG3FlukaCorrOn = kFALSE;

    // fill QA for candidates
    fFillQA = kFALSE;

    //Default: Min-Bias
    fPerformMultiplicityStudy = kFALSE;
    fUseIntegratedEfficiencies = kTRUE;
    fLoMultBound       = -1;
    fHiMultBound       = 10000;

    //Default: no cut in Armenteros
    fSpecialArmenterosCutK0s = kFALSE;

    //Pt Bins: undefined
    fptbinnumb = -1;

    for( Int_t ipt = 0; ipt<200; ipt++) fptbinlimits[ipt] = 0;
    for( Int_t i=0; i<3; i++){ for(int j=0; j<3; j++) fMatrix[i][j]=0.;}

    //Primary selection criteria
    // --- IsPhysicalPrimary....: Usual method
    // --- DistToPV.............: Iouri, Luke
    fPrimarySelection = "IsPhysicalPrimary";
    fDistToPVCut = 0.001;
    ///// TF1
    fLevyFitXiFeedDown = 0x0;
    for(int k=0; k<2; k++) fLevyFitXiFeedDownErr[k] = 0x0;
    fPtLegProtCut  = 0.3;
    fRequestITSTOF = -1;
    fPtMinITSTOF = 1.;
}

/***********************************************
 --- Setters For Configuration ---
***********************************************/

// Filename Setters
void AliV0Module::SetRealDataFile    ( TString RealDataFilename     ) {
    //Set root file containing real data candidates.
    fRealDataFile = RealDataFilename;
}
void AliV0Module::SetMCDataFile      (  TString MCDataFilename0, TString MCDataFilename1, TString MCDataFilename2  ) {
    //Set root file containing Monte Carlo data (for efficiency computation).
    fMCDataFile0   = MCDataFilename0;
    fMCDataFile1   = MCDataFilename1;
    fMCDataFile2   = MCDataFilename2;
}
void AliV0Module::SetFeedDownDataFile( TString FeedDownDataFilename ) {
    //Set root file containing Monte Carlo data (for feeddown computation).
    fFeedDownDataFile = FeedDownDataFilename;
}
void AliV0Module::SetOutputFile      ( TString OutputFilename       ) {
    //Set root filename for the analysis output.
    fOutputDataFile = OutputFilename;
    cout<<"[AliV0Module] Set output file to \'"<<OutputFilename<<"\'."<<endl;
}

// Bin Limit Setter
void AliV0Module::SetPtBinLimits(Long_t got_fptbinnumb,const Double_t *got_fptbinlimits) {
    //Function to set pt binning. First argument is the number of pt bins, second is
    //an array with bin limits.
    fptbinnumb = got_fptbinnumb;
    for(int ix = 0; ix<fptbinnumb+1; ix++) {
        fptbinlimits[ix] = got_fptbinlimits[ix];
    }
    for(int ix = 0; ix<fptbinnumb; ix++) {
        fptX[ix] = (fptbinlimits[ix+1] + fptbinlimits[ix])/2.;
    }
    cout<<"[AliV0Module] Received "<<fptbinnumb<<" pt bins, set accordingly."<<endl;
}

// Rapidity Window Setter
void AliV0Module::SetRapidityWindow(Double_t got_fRapidityBoundary) {
    //Set Rapidity Boundary used in analysis.
    //Value provided will be used as upper limit for the modulus of y (|y|< given value)
    fRapidityBoundary = got_fRapidityBoundary;
    cout<<"[AliV0Module] Received "<<got_fRapidityBoundary<<" as rapidity limits, set accordingly."<<endl;
}

void AliV0Module::SetEventsWeights( Double_t w0, Double_t w1, Double_t w2){
    //Set the number of events per dataset to get weighted final efficiency
    fEventsWeights[0] = w0;
    fEventsWeights[1] = w1;
    fEventsWeights[2] = w2;

}

// CINT1B/INEL Setter for normalization to yields
void AliV0Module::SetCINT1BoverINEL(Double_t got_fCINT1BoverINEL) {
    //Set CINT1B/INEL ratio (determined by collaboration).
    fCINT1BoverINELratio = got_fCINT1BoverINEL;
    cout<<"[AliV0Module] Received CINT1B/INEL = "<<got_fCINT1BoverINEL<<" to normalize with, set accordingly."<<endl;
}

// Topological Selection Setters
void AliV0Module::SetCutV0Radius(Double_t cut) {
    //Set minimum decay radius for the V0 in centimeters.
    //Note: this is (R_{2D}, cylindrical, centered around ALICE detector)
    fCutV0Radius = cut;
    cout<<"[AliV0Module] Received V0 Radius (min value) = "<<cut<<endl;
}
void AliV0Module::SetCutDCANegToPV(Double_t cut) {
    //Set minimum distance of closest approach between V0 negative daughter
    //track and primary vertex (in centimeters).
    fCutDCANegToPV = cut;
    cout<<"[AliV0Module] Received DCA Negative track to PV (min value) = "<<cut<<endl;
}
void AliV0Module::SetCutDCAPosToPV(Double_t cut) {
    //Set minimum distance of closest approach between V0 positive daughter
    //track and primary vertex (in centimeters).
    fCutDCAPosToPV = cut;
    cout<<"[AliV0Module] Received DCA Positive track to PV (min value) = "<<cut<<endl;
}
void AliV0Module::SetCutDCAV0Daughters(Double_t cut) {
    //Set minimum distance of closest approach between V0 daughter
    //tracks, in sigmas. This is not in centimeters because if the
    //tracks have been determined with ITS refit the resolution is
    //greatly improved; thus, the cut may be tighter in that case.
    //Using sigmas will therefore take the tracking method in
    //consideration.
    fCutDCAV0Daughters = cut;
    cout<<"[AliV0Module] Received DCA V0 Daughters (max value) = "<<cut<<endl;
}
void AliV0Module::SetCutV0CosPA(Double_t cut) {
    //Set minimum value for the cosine of pointing angle of the V0.
    fCutV0CosPA = cut;
    cout<<"[AliV0Module] Received V0 Cosine of Pointing Angle (min value) = "<<cut<<endl;
}

// Other Selection Setters
void AliV0Module::SetCutProperLifetime(Double_t cut) {
    //Set maximum value for m*L/p variable for the V0.
    //This is the "proper lifetime" selection and is usually called a
    //"c*tau cut". Should be set to a value larger than the c*tau for
    //the V0 considered.
    fCutProperLifetime = cut;
    cout<<"[AliV0Module] Received proper lifetime cut (max value) = "<<cut<<endl;
}
void AliV0Module::SetCutTPCPIDNSigmas(Double_t cut) {
    //Set maximum deviation from the expected energy loss in
    //the TPC, in multiples of sigmas as computed from the AliPIDResponse
    //object. Selection is only used in real data and should thus
    //be very loose to ensure negligible signal loss.
    fCutTPCPIDNSigmas = cut;
    cout<<"[AliV0Module] Received TPC N-sigmas selection (max dist from BB curve) = "<<cut<<endl;
}
void AliV0Module::SetCutSigmaForSignalExtraction(Double_t cut) {
    //Set number of sigmas for the signal extraction method. The value
    //provided is the number of sigmas from the peak position used for
    //the peak sampling, i.e. the peak will be from [<m>-cut*sigma,<m>+cut*sigma]
    //while the background sampling regions will be from
    //[<m>-2*cut*sigma,<m>+2*cut*sigma.
    fCutNSigmasForSignalExtraction = cut;
    cout<<"[AliV0Module] Received N-sigmas for sig. ext.: peak is (-"<<cut<<",+"<<cut<<") in sigmas"<<endl;
}
void AliV0Module::SetCutLeastNumberOfCrossedRows(Double_t cut) {
    //Set smallest allowed number of TPC clusters for the V0 daughter tracks.
    fCutLeastNumberOfCrossedRows = cut;
    cout<<"[AliV0Module] Received Least Nbr of crossed rows (min value) = "<<cut<<endl;
}
void AliV0Module::SetCutLeastNumberOfCrossedRowsOverFindable(Double_t cut) {
    //Set smallest allowed number of TPC clusters for the V0 daughter tracks.
    fCutLeastNumberOfCrossedRowsOverFindable = cut;
    cout<<"[AliV0Module] Received Least Nbr of crossed rows over findable clusters (min value) = "<<cut<<endl;
}

void AliV0Module::SetCutDaughterEta(Double_t cut) {
    //Set maximum eta value allowed for the V0 daughter tracks (|eta|<cut).
    fCutDaughterEta = cut;
    cout<<"[AliV0Module] Received Daughter |eta| cut (max value) = "<<cut<<endl;
}

void AliV0Module::SetCutCompetingV0Rejection(Double_t cut) {
    //Set rejection window around invariant mass of competing V0 species.
    //If negative, no rejection will occur. If analysis is for Lambdas, K0s
    //will be rejected and vice-versa. Typical values revolve around
    //0.003 - 0.010 GeV/c^2.
    fCutCompetingV0Rejection = cut;
    cout<<"[AliV0Module] Received Competing V0 Rejection window of +/- (in GeV/c^2) = "<<cut<<endl;
}

void AliV0Module::SetFeeddownTreatment ( TString fFDMethod ) {
    //Set method used to compute charged and neutral Xi baryon feeddown to Lambda.
    //Methods allowed are:
    //"NoFD" - No Feeddown subtraction performed.
    //"DoubleChargedXi" - Feeddown is computed for the charged Xi, and is then
    //subtracted twice from real data. Assumes charged and neutral Xi enter
    //the analysis in the same way and are produced at a ratio 1:1.
    //"UseMCRatio" - The Feeddown matrix F_{ij} is filled not only with Lambdas
    //coming from charged Xis but also from neutral ones. The scaling is performed
    //according to the measured charged Xis, and since the matrix is more populated
    //a larger subtraction will occurr. This method assumes that MC correctly
    //reproduces the ratio between charged and neutral Xi baryons. Warning: this
    //should not be used with charged Xi triggered Monte Carlo data.
    fFDSwitch = fFDMethod;
    cout<<"[AliV0Module] Received Feeddown treatment method: "<<fFDMethod<<endl;
    if( fFDMethod == "NoFD" )
        cout<<"[AliV0Module] ---> No Feeddown correction will be performed."<<endl;
    if( fFDMethod == "DoubleChargedXi" )
        cout<<"[AliV0Module] ---> Feeddown performed by doubling charged Xi correction."<<endl;
    if( fFDMethod == "UseMCRatio" )
        cout<<"[AliV0Module] ---> Feeddown performed by using MC neutral/charged Xi."<<endl;
}

void AliV0Module::SetPrimarySelection ( TString lPrimSelMethod ) {
    cout<<"[AliV0ModulePbPb] Received Primary Definition Criteria = "<<lPrimSelMethod<<endl;
    fPrimarySelection = lPrimSelMethod;
}

void AliV0Module::SetDistToPVPrimary ( Double_t fDistToPVSelection ) {
    cout<<"[AliV0ModulePbPb] Received max dist to PV for primary = "<<fDistToPVSelection<<endl;
    fDistToPVCut = fDistToPVSelection;
}

void AliV0Module::SetFitBackground ( Bool_t fitBgSwitch ) {
    //Turns on background fitting for signal extraction instead of pure
    //bin counting. Useful, among other things, for systematics.
    fFitBackgroundSwitch = fitBgSwitch;
}

void AliV0Module::SetGeant3FlukaCorr ( Bool_t switchOnG3Fk ) {
    //Turns on background fitting for signal extraction instead of pure
    //bin counting. Useful, among other things, for systematics.
    fG3FlukaCorrOn = switchOnG3Fk;
}

void AliV0Module::SetPosRap(Bool_t lPosRap){
    fPosRap = lPosRap;
}

void AliV0Module::SetNegRap(Bool_t lNegRap){
    fNegRap = lNegRap;
}

Bool_t AliV0Module::GetPosRap(){
    return fPosRap;
}

Bool_t AliV0Module::GetNegRap(){
    return fNegRap;
}


void AliV0Module::SetPerformMultiplicityStudy ( Bool_t lPerformMultStudy ) {
    //Turns on selection according to multiplicity of the event.
    //Boundaries are set in charged track multiplicity (pp) or in
    //centrality percentiles (PbPb). This is a requirement for studying
    //PbPb.
    fPerformMultiplicityStudy = lPerformMultStudy;
}

void AliV0Module::SetLowMultValue ( Double_t lLoMultBound       ) {
    //Lower boundary (inclusive) in integer number for mult selection.
    //Note: If in PbPb and you want, say, 10-20% centrality, set this
    //to 10.
    fLoMultBound = lLoMultBound;
}

void AliV0Module::SetHighMultValue ( Double_t lHiMultBound       ) {
    //Lower boundary (inclusive) in integer number for mult selection.
    //Note: If in PbPb and you want, say, 10-20% centrality, set this
    //to 10.
    fHiMultBound = lHiMultBound;
}

void AliV0Module::SetSpecialArmenterosCutK0s ( Bool_t lSpecialArmenterosCutK0s ) {
    //Special armenteros cut: |alpha|<5*pt_{arm}
    fSpecialArmenterosCutK0s = lSpecialArmenterosCutK0s;
}

void AliV0Module::SetDefaultCuts() {
    //Sets Default cuts for analysis. (adjusted for adequate pp analysis)
    cout<<"[AliV0Module] Setting default cuts for particle species: "<<fWhichParticle<<endl;
    //Set Cuts - topological
    SetCutV0Radius       (0.500);
    SetCutDCANegToPV     (0.060);
    SetCutDCAPosToPV     (0.060);
    SetCutDCAV0Daughters (1.000);
    if ( fWhichParticle == "K0Short")
        SetCutV0CosPA      (0.970);
    if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
        SetCutV0CosPA      (0.995);
        //SetCutV0CosPA      (0.970);


    //Set Cuts - other
    if ( fWhichParticle == "K0Short")
        SetCutProperLifetime                          ( 20);
    if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
        SetCutProperLifetime                          ( 30);

    SetCutTPCPIDNSigmas                           (  5);
    SetCutSigmaForSignalExtraction                (  6);
    SetCutLeastNumberOfCrossedRows                ( 70);
    SetCutLeastNumberOfCrossedRowsOverFindable    (0.8);
    SetCutDaughterEta                             (0.8);
    if ( fWhichParticle == "K0Short")
        SetCutCompetingV0Rejection      (0.005);
    if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
        SetCutCompetingV0Rejection      (0.010);

    //Primary selection criteria
    // --- IsPhysicalPrimary....: Usual method
    // --- DistToPV.............: Iouri, Luke
    fPrimarySelection = "IsPhysicalPrimary";
    fDistToPVCut = 0.001;
    if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ) SetPtProCut(0.3);
    if( fWhichParticle == "K0Short") SetPtProCut(0.0);
    fRequestITSTOF = -1;
    if(fWhichParticle == "K0Short") fPtMinITSTOF = 1.;
    if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )  fPtMinITSTOF = 1.;
}


void AliV0Module::SetMinimumCuts() {
    //Sets Default cuts for analysis. (adjusted for adequate pp analysis)
    cout<<"[AliV0Module] Setting default cuts for particle species: "<<fWhichParticle<<endl;
    //Set Cuts - topological
    SetCutV0Radius       (0.00);
    SetCutDCANegToPV     (0.00);
    SetCutDCAPosToPV     (0.00);
    SetCutDCAV0Daughters (5000.000);
    if ( fWhichParticle == "K0Short")
        SetCutV0CosPA      (0.0);
    if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
        SetCutV0CosPA      (0.0);
        //SetCutV0CosPA      (0.970);


    //Set Cuts - other
    if ( fWhichParticle == "K0Short")
        SetCutProperLifetime                          ( 2000);
    if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
        SetCutProperLifetime                          ( 3000);

    SetCutTPCPIDNSigmas                           (  50);
    SetCutSigmaForSignalExtraction                (  6);
    SetCutLeastNumberOfCrossedRows                ( 70);
    SetCutLeastNumberOfCrossedRowsOverFindable    (0.8);
    SetCutDaughterEta                             (0.8);
    if ( fWhichParticle == "K0Short")
        SetCutCompetingV0Rejection      (0.000);
    if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
        SetCutCompetingV0Rejection      (0.000);

    //Primary selection criteria
    // --- IsPhysicalPrimary....: Usual method
    // --- DistToPV.............: Iouri, Luke
    fPrimarySelection = "IsPhysicalPrimary";
    fDistToPVCut = 0.001;
    fRequestITSTOF = -1;
    fPtMinITSTOF = 1.;
}


TString AliV0Module::IntToString(int input) {
    //Integer to TString Converter
    char dummyChar[50];
    sprintf(dummyChar, "%d", input);
    TString outputstring = dummyChar;
    return outputstring;
}

TString AliV0Module::DoubleToString(double input) {
    //Double to TString Converter
    char dummyChar[50];
    sprintf(dummyChar, "%.3f", input);
    TString outputstring = dummyChar;
    return outputstring;
}

Double_t AliV0Module::ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( errorfromtop + errorfrombottom );
    }
    return 1;
}

Double_t AliV0Module::MyGeant3FlukaCorrectionForProtons(const Double_t *x, const Double_t *par) {
    //Parametrization Used for Geant3/Fluka Correction for protons
    //Credit: Antonin Maire
    return 1 - par[0]*TMath::Exp(par[1]*x[0]) + par[2];
}


Double_t AliV0Module::MyGeant3FlukaCorrectionForAntiProtons(const Double_t *x, const Double_t *par) {
    // Parametrization Used for Geant3/Fluka Correction for antiprotons
    // Credit: Antonin Maire
    //
    // La fonction A*TMath::Erf( B.x ) ne marche pas à bas pt.
    // Différentes fonctions ont été testées.
    // La fonction suivante semble donner de bons résultats
    // On peut jouer sur la puissance n du terme "1/x^n" pour repousser le pt auquel la fonction prend la valeur 1.0
    // Ici, pour n = 0.2, la fonction prend la valeur 0.9990 en pt = 10 GeV/c

    return 1 - par[0]*TMath::Exp(par[1]*x[0]) + par[2] + par[3]*1/TMath::Power(x[0], 0.2)*TMath::Log(x[0]);
}

Double_t AliV0Module::MyLevyPtXi(const Double_t *pt, const Double_t *par)
{
    //Levy Fit Function
    Double_t lMass  = 1.32171; //Xi Mass
    Double_t ldNdy  = par[0];
    Double_t lTemp = par[1];
    Double_t lPower = par[2];

    Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
    Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);

    return ldNdy * pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
    //return ldNdy * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
}

Double_t AliV0Module::MyBgPol1(const Double_t *x, const Double_t *par)
{
    //Function for background fitting, rejects peak region
    if ( x[0] > par[2] && x[0] < par[3]) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0];
}

Double_t AliV0Module::MyBgPolToEval1(const Double_t *x, const Double_t *par)
{
    //Just a plain linear function.
    return par[0] + par[1]*x[0];
}

Double_t AliV0Module::RoundToThousandth( const Double_t lToRound ) {
    //Round any number to a hundredth...
    //Well, within machine precision...
    return TMath::Nint( 1000 * lToRound ) * 0.001;
}

Bool_t AliV0Module::CheckITSTOF(ULong64_t lPosTrackStatus, ULong64_t lNegTrackStatus, Int_t lPosTOFBCID, Int_t lNegTOFBCID){

    Bool_t ITSrefitAllPtOneLeg = kFALSE;
    if ((lPosTrackStatus & 4) || (lNegTrackStatus & 4) )
        ITSrefitAllPtOneLeg = kTRUE;
    //
    Bool_t TOFmatchAllPtOneLeg = kFALSE;
    if ((TMath::Abs(lPosTOFBCID+100.)>1e-6) || (TMath::Abs(lNegTOFBCID+100.)>1e-6))
        TOFmatchAllPtOneLeg = kTRUE;
    //
    if( ITSrefitAllPtOneLeg || TOFmatchAllPtOneLeg ) {
        return kTRUE;
    }
    else {
        return kFALSE;
    }
}

Bool_t AliV0Module::CheckITSTOFOne(Bool_t DoITS, Double_t pt, Double_t thr, ULong64_t lPosTrackStatus, ULong64_t lNegTrackStatus,
                                                        Int_t lPosTOFBCID, Int_t lNegTOFBCID) {
    //Check whether ITS or TOF is available for one of the daughters

    Bool_t boolvar = kFALSE;

    if (DoITS && (pt < thr)){
        if ((lPosTrackStatus & 4) || (lNegTrackStatus & 4) ) boolvar = kTRUE;
    }
    else if (DoITS && (pt > thr)){
        boolvar =  kTRUE;
    }
    else if (!DoITS && (pt > thr)){
        if ((TMath::Abs(lPosTOFBCID+100.)>1e-6) || (TMath::Abs(lNegTOFBCID+100.)>1e-6) ) boolvar = kTRUE;
    }
    else if (!DoITS && (pt < thr)){
        boolvar =  kTRUE;
    }

    return boolvar;
}

Double_t AliV0Module::DoWeightedMean(Bool_t err , Double_t h1, Double_t h2, Double_t h3, Double_t w1, Double_t w2, Double_t w3) {
    //#eventi nel period i = 1./wi^2
    //wi = 1./sqrt(eventi)

    Double_t wmean,  werror;

    Double_t num = h1/(w1*w1) + h2/(w2*w2) + h3/(w3*w3);
    Double_t den = (1./(w1*w1)) + (1./(w2*w2)) + (1./(w3*w3));

    wmean = num/den;
    werror = TMath::Sqrt(1./((1./(w1*w1)) + (1./(w2*w2)) + (1./(w3*w3))));

    if (err) {
        return werror;
    } else {
        return wmean;
    }
}

Float_t  AliV0Module::DoWeightedMeanXBin(Bool_t err , Float_t h1, Float_t h2, Float_t h3, Float_t e1, Float_t e2, Float_t e3, Float_t w1, Float_t w2, Float_t w3) {

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


void AliV0Module::SetWhichEstimators ( TString got_multestimator = "SPDClusters", TString got_effenergyestimator = "V0M" ){
    //Set estimators (one of V0M, ZDC, SPDClusters)
    fWhichMultEstimator = got_multestimator;
    fWhichEffEnergyEstimator = got_effenergyestimator;
    cout<<"[AliV0Module] Received Multiplicity Estimator: "<<got_multestimator<<endl;
    cout<<"[AliV0Module] Received Effective Energy Estimator: "<<got_effenergyestimator<<endl;
}

void AliV0Module::SetLowEEValue ( Double_t lLoEEBound       ){
    //Lower boundary (inclusive) in integer number for effective energy selection.
    fLoEEBound = lLoEEBound;
}

void AliV0Module::SetHighEEValue ( Double_t lHiEEBound       ){
    //Lower boundary (inclusive) in integer number for effective energy selection.
    fHiEEBound = lHiEEBound;
}

void AliV0Module::SetTOFpercFileName ( TString TOFpercFilename ) {
	//Set root file containing TOF calibrations for multiplicity
    fTOFpercFilename = TOFpercFilename;
}

void AliV0Module::SetN0815percFileName ( TString percFilenameN0815 ) {
	fpercFilenameN0815 = percFilenameN0815;
}

void AliV0Module::SetDoEfficiency(Bool_t DoIt){
    //Do you want to compute weighted efficiency?
    kDoEfficiency = DoIt;
}

Double_t AliV0Module::GetTOFpercentile ( TFile* lfilename , Int_t lNTOFtrgPads, Int_t lRun ){ //set in the run func
	//Converts TOF Pads to TOF percentile.

    TH1F * fhCumulative = (TH1F *)lfilename->Get(Form("%i/hCumulative_%i",lRun, lRun));
	Double_t fTOFpercentile;

 	fTOFpercentile = 100*(1-fhCumulative->Interpolate(lNTOFtrgPads));

    return fTOFpercentile;
}

Double_t AliV0Module::GetPercentilefromValue ( TFile* lfilename, Int_t lRun, Int_t lValue, TString lEstimator ){ //set in the run func

    //Getting percentile
    TString name = Form("hcum%s_%i", lEstimator.Data(), lRun);
    TH1D *hcum = (TH1D *)lfilename->Get(name);

    Double_t percentile = 100 * (hcum->GetBinContent(hcum->FindBin(lValue)));

    return percentile;
}


void AliV0Module::DoAnalysis() {
    //----------------------
    // V0 Analysis Function
    //----------------------
    //
    //Consists of the following steps:
    //
    // (1) Loop over Real data candidates, acquire peak position and widths.
    // (2) Loop over Real data candidates, extract signal with variable extraction
    //     areas. Two loops: totally avoids binning data, allowing for sub-MeV/c^2
    //     granularity in signal extraction areas.
    // (3) Loop over MC data reconstructed candidates, find associated-to-MC primary
    //     candidates for efficiency numerator.
    // (4) (if Lambda or AntiLambda) Perform Geant3/Fluka correction.
    // (5) Get generated primary V0 histograms from MC file for efficiency denominator.
    // (6) (if Lambda or AntiLambda + FD correction enabled) Open MC file for feeddown
    //     subtraction. Loop over candidates, find Lambda associated with primary Xi
    //     to fill feeddown matrix. Scale the feeddown contribution according to real-life
    //     measured Xi- production under specified feeddown subtraction scheme (see
    //     function SetFeeddownTreatment). Perform subtraction of raw Lambda or AntiLambda
    //     count estimated to be coming from charged or neutral Xi.
    // (7) Perform detection efficiency correction and compute final corrected spectra.
    //
    //
    // Normalization:
    //  --- Number of Inelastic Events estimated to be:
    //      <Number of CINT1B triggers> / <CINT1B/INEL ratio>
    //
    //  --- Efficiency denominator: filled at pre-physics selection stage. All signal
    //      losses are therefore estimated from Monte Carlo.
    //
    //  --- Pileup: not included in MC. Fraction of events removed by IsPileupFromSPD
    //      is artificially removed from Number of Inelastic Events as well.
    //
    // Output: written to specified file (with SetOutputFile(...)).
    //   --- includes spectra and efficiencies in base directory.
    //   --- includes a number of different subdirectories storing plots such as
    //       such as invariant mass histos, signal extraction, pt resolution, G3/F
    //       related distributions, feeddown related distributions, and so on.

    //Set Batch Mode: Ignore all canvases in this section
    gROOT->SetBatch (kTRUE);

    cout<<"======================================"<<endl;
    cout<<" -------- Spectra Extraction -------- "<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;
    if(fptbinnumb == -1) {
        cout<<"[AliV0Module] It's unclear what kind of pt binning you want."<<endl;
        cout<<"[AliV0Module] Most likely you forgot to set it using SetPtBinLimits..."<<endl;
        cout<<"[AliV0Module] Analysis will NOT be done. Returning."<<endl;
        return;
    }
    cout<<"--------------- Configuration --------------------------"<<endl;
    cout<<" Analysed Particle.............: "<<fWhichParticle<<endl;
    cout<<" Rapidity Window...............: "<<fRapidityBoundary<<endl;
    cout<<" CINT1B/INEL Ratio used........: "<<fCINT1BoverINELratio<<endl;
    cout<<" V0 Decay Radius...............: "<<fCutV0Radius<<endl;
    cout<<" DCA Negative track to PV......: "<<fCutDCANegToPV<<endl;
    cout<<" DCA Positive track to PV......: "<<fCutDCAPosToPV<<endl;
    cout<<" DCA V0 Daughters..............: "<<fCutDCAV0Daughters<<endl;
    cout<<" Cosine of Pointing Angle V0...: "<<fCutV0CosPA<<endl;
    //cout<<" Cosine of Pointing Angle V0...: parametric"<<endl;
    cout<<" Proper Lifetime cut (cm)......: "<<fCutProperLifetime<<endl;
    cout<<" TPC dE/dx sigmas cut (Real)...: "<<fCutTPCPIDNSigmas<<endl;
    cout<<" CutNSigmasForSignalExtraction.: "<<fCutNSigmasForSignalExtraction<<endl;
    cout<<" Least # of Crossed Rows.......: "<<fCutLeastNumberOfCrossedRows<<endl;
    cout<<" Least # of C. R. over findable: "<<fCutLeastNumberOfCrossedRowsOverFindable<<endl;
    cout<<" Daughter Track |eta| < .......: "<<fCutDaughterEta<<endl;
    cout<<" Competing V0 Reject. (GeV/c^2): "<<fCutCompetingV0Rejection<<endl;
    cout<<"--------------- File Names -----------------------------"<<endl;
    cout<<" Real Data File................: "<<fRealDataFile<<endl;
    cout<<" MC File 1.......................: "<<fMCDataFile0<<endl;
    cout<<" MC File 2.......................: "<<fMCDataFile1<<endl;
    cout<<" MC File 3.......................: "<<fMCDataFile2<<endl;
    cout<<" File for feeddown..........: "<<fFeedDownDataFile<<endl;
    cout<<" Analysis output file..........: "<<fOutputDataFile<<endl;
    cout<<" Setting multiplicity estimator........:" << fWhichMultEstimator.Data() << endl;
    cout<<" Setting eff energy estimator........:" << fWhichEffEnergyEstimator.Data() << endl;
    cout << " Ntrk0815 percentile file......: " << fpercFilenameN0815 << endl;
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;

    //=== Real data loop 1: Acquire peak positions, widths ===
    //--- Preparing... ---
    //defining helping histogram - only used to FindBin Index!========
    TH1F* fHistPt 		= new TH1F("fHistPt","Dummy;p_{T} (GeV/c);Counts",fptbinnumb,fptbinlimits);
    TH1F* fHistPtRaw 		= new TH1F("fHistPtRaw","Raw Counts;#it{p}_{T} (GeV/it{c});Counts",fptbinnumb,fptbinlimits);
    //Peak Position Histograms
    TH1F* fHistPeakPosition	= new TH1F("fHistPeakPosition","Peak Position (Real);p_{T} (GeV/c);Peak Position",fptbinnumb,fptbinlimits);
    TH1F* fHistPeakPositionMC	= new TH1F("fHistPeakPositionMC","Peak Position (MC);p_{T} (GeV/c);Peak Position",fptbinnumb,fptbinlimits);
    //Signal to Noise Histograms
    TH1F* fHistSigToNoise	  = new TH1F("fHistSigToNoise","Signal to Noise Ratio (real);p_{T} (GeV/c);Sig / Bg",fptbinnumb,fptbinlimits);
    TH1F* fHistSigToNoiseMC	= new TH1F("fHistSigToNoiseMC","Signal to Noise Ratio (MC);p_{T} (GeV/c);Sig / Bg",fptbinnumb,fptbinlimits);

    //Signal Extraction Range Histogram
    TH1F* fHistSignalExtractionRange	= new TH1F("fHistSignalExtractionRange","Sig. Ext. Range;p_{T} (GeV/c);Range",fptbinnumb,fptbinlimits);

    //Resolution Histogram (filled with MC)
    TH2F* f2dHistPtResolution 		= new TH2F("f2dHistPtResolution","p_{t} Resolution;p_{t} (reco);p_{t} (mc)",fptbinnumb,fptbinlimits, fptbinnumb, fptbinlimits);

    TH2F* f2dHistPtBlur 		= new TH2F("f2dHistPtBlur","p_{t} blurring matrix;p_{t} (reco);p_{t} (mc)",fptbinnumb,fptbinlimits, fptbinnumb, fptbinlimits);
    TH2F* f2dHistPtSharpen 		= new TH2F("f2dHistPtSharpen","p_{t} sharpening matrix;p_{t} (reco);p_{t} (mc)",fptbinnumb,fptbinlimits, fptbinnumb, fptbinlimits);

    //Efficiency Denominator + Numerator
    TH1F* fHistEffNumerator 		= new TH1F("fHistEffNumerator","Efficiency Numerator;p_{T} (GeV/c);Reconst.+Assoc.",fptbinnumb,fptbinlimits);
    TH1F* fHistEffDenominator 		= new TH1F("fHistEffDenominator","Efficiency Denominator;p_{T} (GeV/c);Generated",fptbinnumb,fptbinlimits);

    TH1F* fHistReco[3];
    for (int i = 0; i<3; i++){
        fHistReco[i] = new TH1F(Form("fHistReco%i",i), "HistReco ;p_{T} (GeV/c);Counts", 250, 0., 25.);
    }

    //Average Multiplicity
    TH1F* fHistAverageMult          = new TH1F("fHistAverageMult","Average Multiplicity;Multiplicity Bin;Average Multiplicity",1,0,1);

     ///
    TH2F *etaDaughters = new TH2F("fEtaDaughs","yMoth; etaD1; etaD2",200,-2.,2.,200,-2.,2.);

    TH1F *lProtonMomentumData[100];
    TH1F *lProtonMomentumMC[100];
     char lNameOne[100];
    for(Int_t ibin=0; ibin<100; ibin++) {
        sprintf(lNameOne,"lProtonMomentumBin%i",ibin);
        lProtonMomentumData[ibin] = new TH1F(Form("%s_data",lNameOne),"",800,0,20);
        lProtonMomentumMC[ibin] = new TH1F(Form("%s_mc",lNameOne),"",800,0,20);
        }
        //
    /********************************************************

      ---> Let's Remember the limits of the data we're analyzing!
      ---> Important so that we don't try signal extraction outside
      ---> the bundaries of available data in the Tree object.

      From AliAnalysisTaskExtractV0.cxx

      //Second Selection: rough 20-sigma band, parametric.
      //K0Short: Enough to parametrize peak broadening with linear function.
      Double_t UpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*tree_lPt;
      Double_t LowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*tree_lPt;

      //Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
      //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
      Double_t UpperLimitLambda = (1.13688e+00) + (5.27838e-03)*tree_lPt + (8.42220e-02)*TMath::Exp(-(3.80595e+00)*tree_lPt);
      Double_t LowerLimitLambda = (1.09501e+00) - (5.23272e-03)*tree_lPt - (7.52690e-02)*TMath::Exp(-(3.46339e+00)*tree_lPt);

    ********************************************************/

    TF1 *fKDataUpper = new TF1("fKDataUpper","[0]+[1]*x",0,20);
    TF1 *fKDataLower = new TF1("fKDataLower","[0]+[1]*x",0,20);
    TF1 *fLDataUpper = new TF1("fLDataUpper","[0]+[1]*x+[2]*TMath::Exp([3]*x)",0,20);
    TF1 *fLDataLower = new TF1("fLDataLower","[0]+[1]*x+[2]*TMath::Exp([3]*x)",0,20);

    fKDataUpper->SetParameter(0, 5.63707e-01);
    fKDataUpper->SetParameter(1, 1.14979e-02);
    fKDataLower->SetParameter(0, 4.30006e-01);
    fKDataLower->SetParameter(1,-1.10029e-02);

    fLDataUpper->SetParameter(0, 1.13688e+00);
    fLDataUpper->SetParameter(1, 5.27838e-03);
    fLDataUpper->SetParameter(2, 8.42220e-02);
    fLDataUpper->SetParameter(3,-3.80595e+00);

    fLDataLower->SetParameter(0, 1.09501e+00);
    fLDataLower->SetParameter(1,-5.23272e-03);
    fLDataLower->SetParameter(2,-7.52690e-02);
    fLDataLower->SetParameter(3,-3.46339e+00);

    //7 Double_t weightMass = 0.;
    Int_t lWeAreAtBin = 0;
    Int_t lWeAreAtBin_pTCorr = 0;
    Double_t lParticleMass = 1.115683;          //needed for proper lifetime selection
    Double_t lCompetingParticleMass = -1; //needed for Competing V0 Rejection
    Double_t lHistoLowerBoundary = -1;
    Double_t lHistoUpperBoundary = -1;
    Long_t lHistoNBins = -1;
    if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ) {
        lParticleMass          = 1.115683;
        lCompetingParticleMass = 0.4976;
        lHistoLowerBoundary = 1.116 - 0.1;
        lHistoUpperBoundary = 1.116 + 0.1;
        lHistoNBins = 200; //1MeV/c^2 binning
    }
    if ( fWhichParticle == "K0Short" ) {
        lParticleMass          = 0.4976;
        lCompetingParticleMass = 1.115683;
        lHistoLowerBoundary = 0.498 - 0.15;
        lHistoUpperBoundary = 0.498 + 0.15;
        lHistoNBins = 300; //1MeV/c^2 binning
    }

    //Setting Up QA plots to use===============================
    TH1F *lHistoV0Radius;
    TH1F *lHistDCAnegToPrim;
    TH1F *lHistDCAposToPrim;
    TH1F *lHistDCAV0Daughters;
    TH1F *lProperLifeTime;
    TH1F *lCosPointingAngle;
    TH2F *lHistoNsigmaTPCpos;
    TH2F *lHistoNsigmaTPCneg;

    TH1F *lHistoV0RadiusMC;
    TH1F *lHistDCAnegToPrimMC;
    TH1F *lHistDCAposToPrimMC;
    TH1F *lHistDCAV0DaughtersMC;
    TH1F *lProperLifeTimeMC;
    TH1F *lCosPointingAngleMC;
    TH2F *lHistoNsigmaTPCposMC;
    TH2F *lHistoNsigmaTPCnegMC;

    if(fFillQA){
    lHistoV0Radius = new TH1F("lHistoV0Radius","V0 radius; V0radius[cm]; Entries", 300, 0., 30.);
    lHistDCAnegToPrim = new TH1F("lHistDCAnegToPrim","DCAnegToPrim; DCAnegToPrim[cm]; Entries", 500, 0., 5.);
    lHistDCAposToPrim = new TH1F("lHistDCAposToPrim","DCAposToPrim; DCAposToPrim[cm]; Entries", 500, 0., 5.);
    lHistDCAV0Daughters = new TH1F("lHistDCAV0Daughters","DCAV0Daughters; DCAV0Daughters[cm]; Entries", 200, 0., 2.);
    if(fWhichParticle == "K0Short") lProperLifeTime = new TH1F("lProperLifeTime","ProperLifeTime; ProperLifeTime[cm]; Entries", 5000, 0., 50.);
    else  lProperLifeTime = new TH1F("lProperLifeTime","ProperLifeTime; ProperLifeTime[cm]; Entries", 1000, 0., 50.);
    lCosPointingAngle = new TH1F("lCosPointingAngle","CosPointingAngle; CosPointingAngle; Entries", 1500, 0.95, 1.);
    lHistoNsigmaTPCpos = new TH2F("NsigmaTPCpos","NsigmaTPCpos; p_{T}^{V0}(GeV/c); NsigmaTPCpos", 100,0.,10., 400, -20., 20.);
    lHistoNsigmaTPCneg = new TH2F("NsigmaTPCneg","NsigmaTPCneg; p_{T}^{V0}(GeV/c); NsigmaTPCneg", 100,0.,10., 400, -20., 20.);
    ///
    lHistoV0RadiusMC = new TH1F("lHistoV0RadiusMC","V0 radius; V0radius[cm]; Entries", 300, 0., 30.);
    lHistDCAnegToPrimMC = new TH1F("lHistDCAnegToPrimMC","DCAnegToPrim; DCAnegToPrim[cm]; Entries", 500, 0., 5.);
    lHistDCAposToPrimMC = new TH1F("lHistDCAposToPrimMC","DCAposToPrim; DCAposToPrim[cm]; Entries", 500, 0., 5.);
    lHistDCAV0DaughtersMC = new TH1F("lHistDCAV0DaughtersMC","DCAV0Daughters; DCAV0Daughters[cm]; Entries", 200, 0., 2.);
    if(fWhichParticle == "K0Short") lProperLifeTimeMC = new TH1F("lProperLifeTimeMC","ProperLifeTime; ProperLifeTime[cm]; Entries", 5000, 0., 50.);
    else lProperLifeTimeMC = new TH1F("lProperLifeTimeMC","ProperLifeTime; ProperLifeTime[cm]; Entries", 1000, 0., 50.);
    lCosPointingAngleMC = new TH1F("lCosPointingAngleMC","CosPointingAngle; CosPointingAngle; Entries", 1500, 0.95, 1.);
    lHistoNsigmaTPCposMC = new TH2F("NsigmaTPCposMC","NsigmaTPCpos; p_{T}^{V0}(GeV/c); NsigmaTPCpos", 100,0.,10., 400, -20., 20.);
    lHistoNsigmaTPCnegMC = new TH2F("NsigmaTPCnegMC","NsigmaTPCneg; p_{T}^{V0}(GeV/c); NsigmaTPCneg", 100,0.,10., 400, -20., 20.);
    }

    //Setting Up Base Histograms to use===============================
    TH1F* lHistoFullV0[100];
    TH1F* lHistoSelectedV0[100];
    TH1F* lHistoFullV0MC[100];
    TCanvas* lCanvasHistoFullV0[100];
    TCanvas* lCanvasHistoFullV0MC[100];
    TH1F *lHistResolution[100];
    TH1F *lHistPtDistr[100];
    TH1F *lHistPtDistrMCtruth[100];
    TF1 *lHistResolutionGaussian[100];
    char lHistResolutionName[100];
    TH1F *fHistMCtruth    = new TH1F("fHistMCSpectra","V0 MC count;p_{T} (GeV/c);Counts",fptbinnumb,fptbinlimits);
    TH1F* fHistResolutionVsPt 		= new TH1F("fHistResolutionVsPt","Resolution vs p_{t};p_{t} (GeV/c);#Delta(p_{t}) = p_{t}-p_{t}^{true} (GeV/c)",fptbinnumb,fptbinlimits);
    TH1F* fHistResolutionVsPtDivByBinWidth 		= new TH1F("fHistResolutionVsPtDivByBinWidth","Resolution vs p_{t};p_{t} (GeV/c);#Delta(p_{t})/#Delta^{bin}(p_{t}) = (p_{t}-p_{t}^{true})/#Delta^{bin}(p_{t}) (GeV/c)",fptbinnumb,fptbinlimits);
    TH1F* fHistResolutionVsPtWithGaussians 		= new TH1F("fHistResolutionVsPtWithGaussians","Resolution vs p_{t};p_{t} (GeV/c);#Delta(p_{t}) = p_{t}-p_{t}^{true} (GeV/c)",fptbinnumb,fptbinlimits);
    TH1F* fHistResolutionVsPtDivByBinWidthWithGaussians 		= new TH1F("fHistResolutionVsPtDivByBinWidthWithGaussians","Resolution vs p_{t};p_{t} (GeV/c);#Delta(p_{t})/#Delta^{bin}(p_{t}) = (p_{t}-p_{t}^{true})/#Delta^{bin}(p_{t}) (GeV/c)",fptbinnumb,fptbinlimits);
    for(Int_t ibin=0; ibin<fptbinnumb; ibin++) {
    }
    char histname[80];
    char histnameSelected[80];
    TString bindescription = "";
    for(Int_t ihist=0; ihist<100; ihist++) {
        //Histo For Real Data
        sprintf(histname,"lHistoFullV0%i",ihist);
        sprintf(histnameSelected,"lHistoSelectedV0%i",ihist);
        if(fWhichParticle == "Lambda")     bindescription="#Lambda, bin #";
        if(fWhichParticle == "AntiLambda") bindescription="#bar{#Lambda}, bin #";
        if(fWhichParticle == "K0Short")    bindescription="K^{0}_{S}, bin #";
        bindescription.Append(IntToString( ihist ));
        if ( ihist < fptbinnumb ) {
            bindescription.Append(", ");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist ]));
            bindescription.Append("-");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist+1 ]));
            bindescription.Append("GeV/c");
        }
        lHistoFullV0[ihist]   =	new TH1F(histname,"Candidates;M(appropriate) (GeV/c^{2});Counts", lHistoNBins,lHistoLowerBoundary,lHistoUpperBoundary);
        lHistoFullV0[ihist]->SetTitle(bindescription);
        lHistoSelectedV0[ihist]   =	new TH1F(histnameSelected,"Candidates;M(appropriate) (GeV/c^{2});Counts", lHistoNBins,lHistoLowerBoundary,lHistoUpperBoundary);
        lHistoSelectedV0[ihist]->SetTitle(bindescription);
        sprintf(histname,"lCanvasHistoFullV0%i",ihist);
        lCanvasHistoFullV0[ihist] = new TCanvas(histname, bindescription, 800, 600);
        lCanvasHistoFullV0[ihist] -> SetFillColor(kWhite);
        lCanvasHistoFullV0[ihist] -> SetLeftMargin( 0.175 );
        lCanvasHistoFullV0[ihist] -> SetBottomMargin( 0.175 );

        //Histo for MC
        sprintf(histname,"lHistoFullV0MC%i",ihist);
        if(fWhichParticle == "Lambda")     bindescription="#Lambda, MC, bin #";
        if(fWhichParticle == "AntiLambda") bindescription="#bar{#Lambda}, MC, bin #";
        if(fWhichParticle == "K0Short")    bindescription="K^{0}_{S}, MC, bin #";
        bindescription.Append(IntToString( ihist ));
        if ( ihist < fptbinnumb ) {
            bindescription.Append(", ");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist ]));
            bindescription.Append("-");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist+1 ]));
            bindescription.Append("GeV/c");
        }
        lHistoFullV0MC[ihist]   =	new TH1F(histname,"Candidates;M(appropriate) (GeV/c^{2});Counts", lHistoNBins,lHistoLowerBoundary,lHistoUpperBoundary);
        lHistoFullV0MC[ihist]->SetTitle(bindescription);
        sprintf(histname,"lCanvasHistoFullV0MC%i",ihist);
        lCanvasHistoFullV0MC[ihist] = new TCanvas(histname, bindescription, 800, 600);
        lCanvasHistoFullV0MC[ihist] -> SetFillColor(kWhite);
        lCanvasHistoFullV0MC[ihist] -> SetLeftMargin( 0.175 );
        lCanvasHistoFullV0MC[ihist] -> SetBottomMargin( 0.175 );

        //Histo for resolution tests
        sprintf(lHistResolutionName,"lHistResolution%i",ihist);
        Long_t lNumberOfBinsResolution = 5000;
        if ( (fptbinlimits[ihist+1]+fptbinlimits[ihist])*0.5 > 5 )
            lNumberOfBinsResolution = 500;
        if ( (fptbinlimits[ihist+1]+fptbinlimits[ihist])*0.5 > 10 )
            lNumberOfBinsResolution = 50;

        lHistResolution[ihist] = new TH1F ( lHistResolutionName,bindescription,lNumberOfBinsResolution, -5, 5); //histo with 5MeV/c precision!
        lHistPtDistr[ihist] = new TH1F ( Form("%s_ptDistr",lHistResolutionName),bindescription,1000, -1., 9.); //histo with 5MeV/c precision!
        lHistPtDistrMCtruth[ihist] = new TH1F ( Form("%s_ptDistrMCtruth",lHistResolutionName),bindescription,1000, -1., 9.); //histo with 5MeV/c precision!
        sprintf(lHistResolutionName,"lHistResolutionGaussian%i",ihist);
        lHistResolutionGaussian[ihist] = new TF1(lHistResolutionName, "[0]*TMath::Gaus(x,[1],[2])",-5,5);

    }
    //================================================================
    cout<<endl;

    cout<<"--------------- Open Real Data File --------------------"<<endl;
    TFile* file = TFile::Open(fRealDataFile, "READ");
    file->cd("PWGLF_StrVsMult");
    TList* v0list  = (TList*)file->FindObjectAny("cList");
    TTree* lTree;
    TTree* lTreeEvent;
    lTree = (TTree*)file->FindObjectAny("fTreeV0");
    lTreeEvent = (TTree*)file->Get("PWGLF_StrVsMult/fTreeEvent");
    TH1F *HistEvents = (TH1F*)v0list->FindObject("fHistEventCounter");
    Float_t nEventsMC  = HistEvents->GetBinContent(1);

    Float_t fMultCentrality = 0.;
    Float_t fEnergyCentrality = 0.;
    Bool_t  fMVPileupFlag = 0;
    Int_t   fRun = 0;
    Int_t   fTOFPads = 0;
    Float_t ZDCFired = -1.;
    Int_t fNTracksGlobal0815 = 0.;
    Int_t fSPDtracklets0815 = 0;

    if (fWhichMultEstimator.Contains("SPDCl") || fWhichMultEstimator.Contains("V0M")) {
        lTreeEvent->SetBranchAddress(Form("fCentrality_%s",fWhichMultEstimator.Data()), &fMultCentrality);
    }
    lTreeEvent->SetBranchAddress(Form("fCentrality_%s",fWhichEffEnergyEstimator.Data()), &fEnergyCentrality);
    lTreeEvent->SetBranchAddress("fNTOFtrgPads",&fTOFPads);
    lTreeEvent->SetBranchAddress("fRun", &fRun);
    lTreeEvent->SetBranchAddress("fMVPileupFlag", &fMVPileupFlag);
    lTreeEvent->SetBranchAddress(Form("fCentrality_%s","ZDCFired"), &ZDCFired);
    lTreeEvent->SetBranchAddress("fNTracksGlobal0815", &fNTracksGlobal0815);
    lTreeEvent->SetBranchAddress("fSPDtracklets0815", &fSPDtracklets0815);

    TFile* Read =  0x0;
    if (fWhichMultEstimator.Contains("TOF")) Read = new TFile (fTOFpercFilename);
    TFile* Read0815 =  0x0;
    if (fWhichMultEstimator.Contains("0815")) Read0815 = new TFile (fpercFilenameN0815);

    Long_t lNInelasticEvents = 0;

    cout<<"--------------------------------------------------------"<<endl;
    cout<<" Will now loop over events, please wait..."<<endl;

    lNInelasticEvents = 0;
    for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries(); iEv++) {

            lTreeEvent->GetEntry(iEv);
            // printf("cent = %f \n",fCentrality_V0M); getchar();
            if ( iEv % ( lTreeEvent->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;

            //If you want to use TOF percentile
            if (fWhichMultEstimator.Contains("TOF")) {
                cout << "\n........ Setting TOF clusters as multiplicity estimator ........\n" << endl;
                fMultCentrality = GetTOFpercentile(Read, fTOFPads, fRun);
            }

            //If you want to use 0815 percentile
            if (fWhichMultEstimator.Contains("0815")) {
                if (fWhichMultEstimator.Contains("SPDtrk"))
                    fMultCentrality = GetPercentilefromValue(Read0815, fRun, fSPDtracklets0815, fWhichMultEstimator.Data());
                if (fWhichMultEstimator.Contains("NtrkGlobal"))
                    fMultCentrality = GetPercentilefromValue(Read0815, fRun, fNTracksGlobal0815, fWhichMultEstimator.Data());
                if (iEv==0)
                    cout << "INFO: Multiplicity percentile = " << fWhichMultEstimator.Data() << endl;
            }

            if (ZDCFired>100) continue;

            //Count events
            if( fPerformMultiplicityStudy == kTRUE &&  //inside mult bin
                fMultCentrality>=fLoMultBound &&
                fMultCentrality<fHiMultBound &&
                fEnergyCentrality>=fLoEEBound &&
                fEnergyCentrality<fHiEEBound
            ) lNInelasticEvents++;
        }
    cout<<"--------------------------------------------------------"<<endl;

    //Variable Definition=============================================
    //Kinematic
    Float_t lPt, lRap, lPtMC, lNegEta, lPosEta;
    Float_t lPosPx, lPosPy, lPosPz;
    Float_t lNegPx, lNegPy, lNegPz;
    //Invariant Masses
    Float_t lInvariantMass, lInvariantMassCompetingOne, lInvariantMassCompetingTwo; //wildcard
    //DCA Variables
    Float_t lDcaV0Daughters;
    Float_t lDcaPosToPrimVertex,  lDcaNegToPrimVertex;
    //Cosine of Pointing Angle variable
    Float_t lV0CosinePointingAngle;
    Float_t lDCAV0ToPrimVertex;
    //Decay Radius and distance over total momentum
    Float_t lV0Radius, lDistOverTotMom;
    //Least Number of TPC Clusters
    Int_t lLeastNbrCrossedRows;
    Float_t lLeastNbrCrossedRowsOverFindable;
    //TPC dE/dx acquired with AliPIDResponse class
    Float_t lNSigmasPosProton,lNSigmasNegProton,lNSigmasPosPion,lNSigmasNegPion;
    Float_t lArmPt,lArmAlpha;
    //Multiplicity Variable
    Float_t lMultiplicity = 0; //for ease of handling later
    Float_t lEnergyPercentile = -1.;
    ///
    Float_t fTreeVariablePtLegPos;
    Float_t fTreeVariablePtLegNeg;
    Float_t lPtXv0Pos;
    Float_t lPtYv0Pos;
    Float_t lPtZv0Pos;
    Float_t lPtXv0Neg;
    Float_t lPtYv0Neg;
    Float_t lPtZv0Neg;
    Int_t pidFromTrackNeg;
    Int_t pidFromTrackPos;
    Bool_t lITSrefitPos = kFALSE;
    Bool_t lITSrefitNeg = kFALSE;
    Int_t fTreeVariableNegTOFBCid;
    Int_t fTreeVariablePosTOFBCid;
    Bool_t ITSrefitAllPtOneLeg;
    Bool_t TOFmatchAllPtOneLeg;
    Bool_t ITSrefitLowPtBothLegs;
    Bool_t TOFmatchHighPtBothLegs;
    Bool_t TOFmatchHighPtOneLeg;
    //ITS
    ULong64_t lNegTrackStatus, lPosTrackStatus;
    //================================================================

    //Linking to Tree=================================================
    //--- Base Variables ----------------------------------------------
    lTree->SetBranchAddress("fTreeVariablePosEta",&lPosEta);
    lTree->SetBranchAddress("fTreeVariableNegEta",&lNegEta);
    lTree->SetBranchAddress("fTreeVariableNegTOFBCid",&fTreeVariableNegTOFBCid);
    lTree->SetBranchAddress("fTreeVariablePosTOFBCid",&fTreeVariablePosTOFBCid);
    lTree->SetBranchAddress("fTreeVariablePt",&lPt);
    //--- ITS flag -----------------------------------------------------
    lTree->SetBranchAddress("fTreeVariablePosTrackStatus", &lPosTrackStatus);
    lTree->SetBranchAddress("fTreeVariableNegTrackStatus", &lNegTrackStatus);
    if ( fWhichParticle == "Lambda"      )  lTree->SetBranchAddress("fTreeVariableInvMassLambda",&lInvariantMass);
    if ( fWhichParticle == "AntiLambda"  ) 	lTree->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMass);
    if ( fWhichParticle == "K0Short"     ) 	lTree->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMass);
    if ( fWhichParticle == "Lambda"      ) { //For symmetry of computation...
        lTree->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
        lInvariantMassCompetingTwo = -1; //will always be far enough
    }
    if ( fWhichParticle == "AntiLambda"  ) { //For symmetry of computation...
        lTree->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
        lInvariantMassCompetingTwo = -1; //will always be far enough
    }
    if ( fWhichParticle == "K0Short"     ) {
        lTree->SetBranchAddress("fTreeVariableInvMassLambda"    ,&lInvariantMassCompetingOne);
        lTree->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMassCompetingTwo);
    }
    if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" )
        lTree->SetBranchAddress("fTreeVariableRapLambda",&lRap);
    if ( fWhichParticle == "K0Short" )
        lTree->SetBranchAddress("fTreeVariableRapK0Short",&lRap);
    lTree->SetBranchAddress("fTreeVariableDistOverTotMom",&lDistOverTotMom);
    // lTree->SetBranchAddress("fTreeVariableIsCowboy",&fTreeVariableIsCowboy);
    lTree->SetBranchAddress("fTreeVariableLeastNbrCrossedRows",&lLeastNbrCrossedRows);
    lTree->SetBranchAddress("fTreeVariableLeastRatioCrossedRowsOverFindable",&lLeastNbrCrossedRowsOverFindable);
    //--- TPC dEdx Variables ------------------------------------------
    lTree->SetBranchAddress("fTreeVariableNSigmasPosProton",&lNSigmasPosProton);
    lTree->SetBranchAddress("fTreeVariableNSigmasNegProton",&lNSigmasNegProton);
    lTree->SetBranchAddress("fTreeVariableNSigmasPosPion",&lNSigmasPosPion);
    lTree->SetBranchAddress("fTreeVariableNSigmasNegPion",&lNSigmasNegPion);
    //--- Topological selection variables -----------------------------
    lTree->SetBranchAddress("fTreeVariableV0Radius",&lV0Radius);
    lTree->SetBranchAddress("fTreeVariableDcaNegToPrimVertex",&lDcaNegToPrimVertex);
    lTree->SetBranchAddress("fTreeVariableDcaPosToPrimVertex",&lDcaPosToPrimVertex);
    lTree->SetBranchAddress("fTreeVariableDcaV0Daughters",&lDcaV0Daughters);
    lTree->SetBranchAddress("fTreeVariableDcaV0ToPrimVertex",&lDCAV0ToPrimVertex);
    lTree->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle",&lV0CosinePointingAngle);
    //--- Multiplicity Variable ----------------------------------------
    if (fWhichMultEstimator.Contains("SPDCl") || fWhichMultEstimator.Contains("V0M")){
        lTree->SetBranchAddress(Form("fTreeVariableCentrality_%s",fWhichMultEstimator.Data()), &fMultCentrality);
    }
    lTree->SetBranchAddress(Form("fTreeVariableCentrality_%s",fWhichEffEnergyEstimator.Data()), &fEnergyCentrality);
    lTree->SetBranchAddress("fTreeVariablePosInnerP",&fTreeVariablePtLegPos);
    lTree->SetBranchAddress("fTreeVariableNegInnerP",&fTreeVariablePtLegNeg);
    /*lTree->SetBranchAddress("fTreeVariablePIDTrackingPos",&pidFromTrackPos);
    lTree->SetBranchAddress("fTreeVariablePIDTrackingNeg",&pidFromTrackNeg);*/
    lTree->SetBranchAddress("fTreeVariableRun", &fRun);
    lTree->SetBranchAddress("fTreeVariablePosPx"     , &lPosPx);
    lTree->SetBranchAddress("fTreeVariablePosPy"     , &lPosPy);
    lTree->SetBranchAddress("fTreeVariablePosPz"     , &lPosPz);
    lTree->SetBranchAddress("fTreeVariableNegPx"     , &lNegPx);
    lTree->SetBranchAddress("fTreeVariableNegPy"     , &lNegPy);
    lTree->SetBranchAddress("fTreeVariableNegPz"     , &lNegPz);
    lTree->SetBranchAddress(Form("fTreeVariableCentrality_%s","ZDCFired"), &ZDCFired);
    lTree->SetBranchAddress("fTreeVariableSPDtracklets0815", &fSPDtracklets0815);

    //================================================================
    Long_t lNCandidates = lTree->GetEntries();
    cout<<"--------------- Real Data File Loop 1 ------------------"<<endl;
    Long_t lOneTenthOfNCandidates = ((double)(lNCandidates) / 10. );
    for(Long_t icand = 0; icand<lNCandidates; icand++) {
        lTree->GetEntry(icand);
        if( icand % lOneTenthOfNCandidates == 0 )
            cout<<" Currently at candidate........: "<<icand<<" / "<<lNCandidates<<" ( "<<(long)(((double)(icand)/(double)(lNCandidates))*(100.+1e-3))<<"% )"<<endl;

        //Multiplicity Switch -- use integrated sample for peak finding
        lMultiplicity = (Double_t)fMultCentrality;

        if (ZDCFired>100) continue;

        if (fWhichMultEstimator.Contains("TOF")) {
            cout << "\n........ Setting TOF clusters as multiplicity estimator ........\n" << endl;
            lMultiplicity = GetTOFpercentile(Read, fTOFPads, fRun);
        }

        //If you want to use 0815 percentile
        if (fWhichMultEstimator.Contains("0815")) {
            if (fWhichMultEstimator.Contains("SPDtrk"))
                lMultiplicity = GetPercentilefromValue(Read0815, fRun, fSPDtracklets0815, fWhichMultEstimator.Data());
            if (fWhichMultEstimator.Contains("NtrkGlobal"))
                lMultiplicity = GetPercentilefromValue(Read0815, fRun, fNTracksGlobal0815, fWhichMultEstimator.Data());
            if (icand == 0)
                cout << "INFO: Multiplicity percentile = %s" << fWhichMultEstimator.Data() << endl;
        }

        if(TMath::Abs(lRap)<fRapidityBoundary &&
                //((!fPosRap && !fNegRap) || (fPosRap && lRap > 0. && lNegEta > 0. && lPosEta >0) || (fNegRap && lRap<0. && lNegEta < 0. && lPosEta <0)) &&
                TMath::Abs(lNegEta)       <= fCutDaughterEta               &&
                TMath::Abs(lPosEta)       <= fCutDaughterEta               &&
                lV0Radius                 >= fCutV0Radius                  &&
                lDcaNegToPrimVertex       >= fCutDCANegToPV                &&
                lDcaPosToPrimVertex       >= fCutDCAPosToPV                &&
                lDcaV0Daughters           <= fCutDCAV0Daughters            &&
                lV0CosinePointingAngle    >= fCutV0CosPA                   &&
                lParticleMass*lDistOverTotMom    <= fCutProperLifetime     &&
                lLeastNbrCrossedRows             >= fCutLeastNumberOfCrossedRows             &&
                lLeastNbrCrossedRowsOverFindable >= fCutLeastNumberOfCrossedRowsOverFindable &&
                TMath::Abs(lInvariantMassCompetingOne - lCompetingParticleMass) > fCutCompetingV0Rejection &&
                TMath::Abs(lInvariantMassCompetingTwo - lCompetingParticleMass) > fCutCompetingV0Rejection &&
                //( fSpecialArmenterosCutK0s == kFALSE || ( fSpecialArmenterosCutK0s == kTRUE && lArmPt*5>TMath::Abs(lArmAlpha) ) ) &&
                ( //official response code
                    ( fWhichParticle == "Lambda"
                      && TMath::Abs(lNSigmasNegPion)   <= fCutTPCPIDNSigmas
                      && TMath::Abs(lNSigmasPosProton) <= fCutTPCPIDNSigmas) ||
                    ( fWhichParticle == "AntiLambda"
                      && TMath::Abs(lNSigmasPosPion)   <= fCutTPCPIDNSigmas
                      && TMath::Abs(lNSigmasNegProton) <= fCutTPCPIDNSigmas) ||
                    ( fWhichParticle == "K0Short"
                      && TMath::Abs(lNSigmasPosPion)   <= fCutTPCPIDNSigmas
                      && TMath::Abs(lNSigmasNegPion)   <= fCutTPCPIDNSigmas) )
                && //OOB condition
                 (
                 (fWhichParticle == "Lambda" && fTreeVariablePtLegPos > fPtLegProtCut ) ||
                 (fWhichParticle == "AntiLambda" && fTreeVariablePtLegNeg > fPtLegProtCut ) ||
                 (fWhichParticle == "K0Short")
                )  &&
                CheckITSTOF(lPosTrackStatus,  lNegTrackStatus,  fTreeVariablePosTOFBCid,  fTreeVariableNegTOFBCid)
                //CheckITSTOFOne( kFALSE, lPt, 2., lPosTrackStatus,  lNegTrackStatus,  fTreeVariablePosTOFBCid,  fTreeVariableNegTOFBCid)

          ) {

            lWeAreAtBin = fHistPt->FindBin(lPt)-1;
            if(lWeAreAtBin == -1) lWeAreAtBin = 99; //UnderFlow, special treatment

            lHistoFullV0[lWeAreAtBin]->Fill(lInvariantMass); //fill with specific inv mass

            if(lMultiplicity>=fLoMultBound &&  lMultiplicity<fHiMultBound) lHistoSelectedV0[lWeAreAtBin]->Fill(lInvariantMass); //fill with specific inv mass

            if(fFillQA){
                lHistoV0Radius->Fill(lV0Radius);
                lHistDCAnegToPrim->Fill(lDcaNegToPrimVertex);
                lHistDCAposToPrim->Fill(lDcaPosToPrimVertex);
                lHistDCAV0Daughters->Fill(lDcaV0Daughters);
                lProperLifeTime->Fill(lParticleMass*lDistOverTotMom);
                lCosPointingAngle->Fill(lV0CosinePointingAngle);
                if(fWhichParticle == "Lambda") {
                    lHistoNsigmaTPCpos->Fill(lPt, lNSigmasPosProton);
                    lHistoNsigmaTPCneg->Fill(lPt, lNSigmasNegPion);
                }
                if(fWhichParticle == "AntiLambda") {
                    lHistoNsigmaTPCpos->Fill(lPt, lNSigmasPosPion);
                    lHistoNsigmaTPCneg->Fill(lPt, lNSigmasNegProton);
                }
                   if(fWhichParticle == "K0Short") {
                    lHistoNsigmaTPCpos->Fill(lPt, lNSigmasPosPion);
                    lHistoNsigmaTPCneg->Fill(lPt, lNSigmasNegPion);
                }
            }
          }
    }
    cout<<"--------------- Loop Completed -------------------------"<<endl;
    cout<<endl;

    cout<<"--------------- Peak Finding (gauss+linear) ------------"<<endl;
    //== Variables for holding peak position, width ==============
    Double_t lPeakPosition[100];
    Double_t lPeakWidth[100];
    Double_t lLeftBgLeftLimit[100];
    Double_t lLeftBgRightLimit[100];
    Double_t lPeakLeftLimit[100];
    Double_t lPeakRightLimit[100];
    Double_t lRightBgLeftLimit[100];
    Double_t lRightBgRightLimit[100];
    //May be needed if bg areas are different in the future

    TLine *lLineLeftMost[100];
    TLine *lLineLeft[100];
    TLine *lLineRight[100];
    TLine *lLineRightMost[100];

    TLine *lLineLeftMostMC[100];
    TLine *lLineLeftMC[100];
    TLine *lLineRightMC[100];
    TLine *lLineRightMostMC[100];

    Double_t lMiddle = 0;
    Double_t lUpperFit = 0;
    Double_t lLowerFit = 0;

    char fgausname[100];
    TF1 *fgausPt[100];
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        cout<<"---> Peak Finding, bin #"<<ibin<<"..."<<endl;
        sprintf(fgausname,"fGausPt%i",ibin);
        if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ) {
            lMiddle = 0.5 * ( fLDataLower->Eval( fHistPt->GetBinCenter(ibin+1) ) + fLDataUpper->Eval( fHistPt->GetBinCenter(ibin+1) ) );
            lUpperFit = lMiddle + 0.7 * ( fLDataUpper->Eval( fHistPt->GetBinCenter(ibin+1) ) - lMiddle );
            lLowerFit = lMiddle - 0.7 * ( lMiddle - fLDataLower->Eval( fHistPt->GetBinCenter(ibin+1) ) );
            fgausPt[ibin]= new TF1(fgausname,"[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]", lLowerFit, lUpperFit );
            fgausPt[ibin]->SetParameter(1,1.116);
            fgausPt[ibin]->SetParameter(2,0.0025);
            fgausPt[ibin]->SetParLimits(1,1,1.2);
            fgausPt[ibin]->SetParLimits(2,0.001,0.01);
        }
        if ( fWhichParticle == "K0Short") {
            fgausPt[ibin]= new TF1(fgausname,"[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]",0.498-0.06,0.498+0.060);
            fgausPt[ibin]->SetParameter(1,0.498);
            fgausPt[ibin]->SetParameter(2,0.004);
        }
        fgausPt[ibin]->SetParameter(0,lHistoFullV0[ibin]->GetMaximum() * 0.9);


        fgausPt[ibin]->SetParameter(3,0);
        fgausPt[ibin]->SetParameter(4,lHistoFullV0[ibin]->GetMaximum() * 0.1);
        //if(ibin == 0) lHistoFullV0[ibin]->Fit(fgausname,"QREM0","",1.105,1.118);
        //else
        lHistoFullV0[ibin]->Fit(fgausname,"QREM0");
        lPeakPosition[ibin] = fgausPt[ibin]->GetParameter(1);
        lPeakWidth[ibin] = TMath::Abs( fgausPt[ibin]->GetParameter(2) );
        cout<<"---> ["<<fptbinlimits[ibin]<<" - "<<fptbinlimits[ibin+1]<<" GeV/c]\tPeak at: "<<lPeakPosition[ibin]<<", sigma = "<<lPeakWidth[ibin]<<endl;
        //Find Corresponding Limits In this bin
        lLeftBgLeftLimit[ibin]  = lPeakPosition[ibin] - 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
        lLeftBgRightLimit[ibin] = lPeakPosition[ibin] - 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
        lPeakLeftLimit[ibin]    = lPeakPosition[ibin] - 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
        lPeakRightLimit[ibin]   = lPeakPosition[ibin] + 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
        lRightBgLeftLimit[ibin] = lPeakPosition[ibin] + 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
        lRightBgRightLimit[ibin]= lPeakPosition[ibin] + 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
        fHistPeakPosition->SetBinContent(ibin+1, lPeakPosition[ibin]);
        fHistPeakPosition->SetBinError(ibin+1, lPeakWidth[ibin]);
        //Create Signal Extraction Range Histogram
        fHistSignalExtractionRange->SetBinContent(ibin+1, lPeakPosition[ibin]);
        fHistSignalExtractionRange->SetBinError(ibin+1, 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin] );
        //Create appropriate TLine Objects for Canvases
        lLineLeftMost[ibin]  = new TLine( lLeftBgLeftLimit[ibin],   0, lLeftBgLeftLimit[ibin],   lHistoFullV0[ibin]->GetMaximum() * 0.95 );
        lLineLeft[ibin]      = new TLine( lLeftBgRightLimit[ibin],  0, lLeftBgRightLimit[ibin],  lHistoFullV0[ibin]->GetMaximum() * 0.95 );
        lLineRight[ibin]     = new TLine( lRightBgLeftLimit[ibin],  0, lRightBgLeftLimit[ibin],  lHistoFullV0[ibin]->GetMaximum() * 0.95 );
        lLineRightMost[ibin] = new TLine( lRightBgRightLimit[ibin], 0, lRightBgRightLimit[ibin], lHistoFullV0[ibin]->GetMaximum() * 0.95 );

        //Preparing Canvas for storing...
        lCanvasHistoFullV0[ibin]->cd();
        lHistoFullV0[ibin]->Draw();
        lLineLeftMost[ibin]->Draw();
        lLineLeft[ibin]->Draw();
        lLineRight[ibin]->Draw();
        lLineRightMost[ibin]->Draw();
        fgausPt[ibin]->Draw("same");

    }
    cout<<"--------------- Peak Finding Finished ------------------"<<endl;

    //Defining Signal holding variables===============================
    Double_t lSigRealV0[100];
    Double_t lSigErrRealV0[100];
    Double_t lSigMCV0[100];
    Double_t lSigErrMCV0[100];
    Double_t lSigMCV0_pTCorr[100];
    Double_t lSigErrMCV0_pTCorr[100];
    //================================================================
    Double_t lLeftPlusRightBgV0[100];
    Double_t lSigPlusCenterBgV0[100];
    Double_t lLeftPlusRightBgV0MC[100];
    Double_t lSigPlusCenterBgV0MC[100];
    for(Long_t n1 = 0; n1<100; n1++) {
        lLeftPlusRightBgV0[n1]=0;
        lSigPlusCenterBgV0[n1]=0;
        lLeftPlusRightBgV0MC[n1]=0;
        lSigPlusCenterBgV0MC[n1]=0;
    }
    //MC pT-shape correction
    Double_t lLeftPlusRightBgV0MC_pTCorr[250];
    Double_t lSigPlusCenterBgV0MC_pTCorr[250];
    for(Long_t n1 = 0; n1<250; n1++){
        lLeftPlusRightBgV0MC_pTCorr[n1]=0;
        lSigPlusCenterBgV0MC_pTCorr[n1]=0;
    }
    //================================================================
    cout<<endl;
    cout<<"--------------- Real Data File Loop 2 ------------------"<<endl;
    for(Long_t icand = 0; icand<lNCandidates; icand++) {
        lTree->GetEntry(icand);

        if( icand % lOneTenthOfNCandidates == 0 )
            cout<<" Currently at candidate........: "<<icand<<" / "<<lNCandidates<<" ( "<<(long)(((double)(icand)/(double)(lNCandidates))*(100.+1e-3))<<"% )"<<endl;

        //Multiplicity Switch -- use integrated sample for peak finding
        lMultiplicity = (Double_t)fMultCentrality;
        //Eff energy Switch
        lEnergyPercentile = (Double_t)fEnergyCentrality;

        if (ZDCFired>100) continue;

         if (fWhichMultEstimator.Contains("TOF")) {
            cout << "\n........ Setting TOF clusters as multiplicity estimator ........\n" << endl;
            lMultiplicity = GetTOFpercentile(Read, fTOFPads, fRun);
        }

        //If you want to use 0815 percentile
        if (fWhichMultEstimator.Contains("0815")) {
            if (fWhichMultEstimator.Contains("SPDtrk"))
                lMultiplicity = GetPercentilefromValue(Read0815, fRun, fSPDtracklets0815, fWhichMultEstimator.Data());
            if (fWhichMultEstimator.Contains("NtrkGlobal"))
                lMultiplicity = GetPercentilefromValue(Read0815, fRun, fNTracksGlobal0815, fWhichMultEstimator.Data());
            if (icand == 0)
                cout << "INFO: Multiplicity percentile = %s" << fWhichMultEstimator.Data() << endl;
        }

        if( fPerformMultiplicityStudy && (lMultiplicity<fLoMultBound || lMultiplicity>=fHiMultBound)) continue;
        if  (fPerformMultiplicityStudy && (lEnergyPercentile<fLoEEBound || lEnergyPercentile>=fHiEEBound) ) continue;

        if(TMath::Abs(lRap)<fRapidityBoundary &&
              //  ((!fPosRap && !fNegRap) || (fPosRap && lRap > 0. && lNegEta > 0. && lPosEta >0 ) || (fNegRap && lRap<0. && lNegEta < 0. && lPosEta < 0)) &&
                TMath::Abs(lNegEta)       <= fCutDaughterEta               &&
                TMath::Abs(lPosEta)       <= fCutDaughterEta               &&
                lV0Radius                 >= fCutV0Radius                  &&
                lDcaNegToPrimVertex       >= fCutDCANegToPV                &&
                lDcaPosToPrimVertex       >= fCutDCAPosToPV                &&
                lDcaV0Daughters           <= fCutDCAV0Daughters            &&
                lV0CosinePointingAngle    >= fCutV0CosPA                   &&
                lParticleMass*lDistOverTotMom    <= fCutProperLifetime     &&
                lLeastNbrCrossedRows             >= fCutLeastNumberOfCrossedRows             &&
                lLeastNbrCrossedRowsOverFindable >= fCutLeastNumberOfCrossedRowsOverFindable &&
                TMath::Abs(lInvariantMassCompetingOne - lCompetingParticleMass) > fCutCompetingV0Rejection &&
                TMath::Abs(lInvariantMassCompetingTwo - lCompetingParticleMass) > fCutCompetingV0Rejection &&
               // ( fSpecialArmenterosCutK0s == kFALSE || ( fSpecialArmenterosCutK0s == kTRUE && lArmPt*5>TMath::Abs(lArmAlpha) ) ) &&
                ( //official response code
                    ( fWhichParticle == "Lambda"
                      && TMath::Abs(lNSigmasNegPion)   <= fCutTPCPIDNSigmas
                      && TMath::Abs(lNSigmasPosProton) <= fCutTPCPIDNSigmas) ||
                    ( fWhichParticle == "AntiLambda"
                      && TMath::Abs(lNSigmasPosPion)   <= fCutTPCPIDNSigmas
                      && TMath::Abs(lNSigmasNegProton) <= fCutTPCPIDNSigmas) ||
                    ( fWhichParticle == "K0Short"
                      && TMath::Abs(lNSigmasPosPion)   <= fCutTPCPIDNSigmas
                      && TMath::Abs(lNSigmasNegPion)   <= fCutTPCPIDNSigmas)
                )  &&
                (
                 (fWhichParticle == "Lambda" && fTreeVariablePtLegPos > fPtLegProtCut ) ||
                 (fWhichParticle == "AntiLambda" && fTreeVariablePtLegNeg > fPtLegProtCut ) ||
                 (fWhichParticle == "K0Short")
                )  &&
                //OOB condition
                CheckITSTOF(lPosTrackStatus,  lNegTrackStatus,  fTreeVariablePosTOFBCid,  fTreeVariableNegTOFBCid)
                //CheckITSTOFOne( kFALSE, lPt, 2., lPosTrackStatus,  lNegTrackStatus,  fTreeVariablePosTOFBCid,  fTreeVariableNegTOFBCid)

          ) { // Start Entry Loop

            lWeAreAtBin = fHistPt->FindBin(lPt)-1;
            if(lWeAreAtBin == -1) lWeAreAtBin = 99; //UnderFlow, special treatment

            //Extract left and right background areas and peak
            //--- Left Background Sample
            if(lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]    && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]   ) {
                lLeftPlusRightBgV0[lWeAreAtBin]++;
            }

            if(lInvariantMass>lRightBgLeftLimit[lWeAreAtBin]   && lInvariantMass<lRightBgRightLimit[lWeAreAtBin]  ) {
                lLeftPlusRightBgV0[lWeAreAtBin]++;
            }

            //--- Peak Region
            if(lInvariantMass>lPeakLeftLimit[lWeAreAtBin]      && lInvariantMass<lPeakRightLimit[lWeAreAtBin]     ) {
                lSigPlusCenterBgV0[lWeAreAtBin]++;
            }

            if ( fWhichParticle == "Lambda"  /* && pidFromTrackPos != 4 */ ) lProtonMomentumData[lWeAreAtBin]->Fill( fTreeVariablePtLegPos);
            if ( fWhichParticle == "AntiLambda" /* && pidFromTrackNeg != 4 */ ) lProtonMomentumData[lWeAreAtBin]->Fill( fTreeVariablePtLegNeg );
        } // End Entry Loop
    }
    cout<<"--------------- Loop Completed -------------------------"<<endl;
    cout<<endl;
    cout<<"--------------- Memory Cleanup -------------------------"<<endl;
    if ( v0list ) {
        v0list->Delete();
        delete v0list;
    }

    file->Close("R");
    file->Delete();
    delete file;
    cout<<endl;

    TF1 *lfitNoise[50];
    TF1 *lSampleNoise[50];
    char lFitNameOne[50];

    if( fFitBackgroundSwitch ) {
        cout<<"--------------- Backgrounds: fitting with linear -------"<<endl;
        for(long i=0; i<fptbinnumb; i++) {
            //Define Function to Fit Background
            sprintf(lFitNameOne,"lfitNoise%i",((int)(i)));
            //cout<<"creating fitnoise, named "<<lFitNameOne<<endl;
            lfitNoise[i] = new TF1(lFitNameOne, this, &AliV0Module::MyBgPol1,
                                   RoundToThousandth ( lLeftBgLeftLimit[i]   ),
                                   RoundToThousandth ( lRightBgRightLimit[i] ), 4 , "AliV0Module", "MyBgPol1");
            lfitNoise[i] -> FixParameter(2,RoundToThousandth ( lLeftBgRightLimit[i] ) );
            lfitNoise[i] -> FixParameter(3,RoundToThousandth ( lRightBgLeftLimit[i] ) );
            lfitNoise[i] -> SetParameter(0,lLeftPlusRightBgV0[i] * lHistoFullV0[i]->GetBinWidth(5) / (lRightBgLeftLimit[i]-lLeftBgRightLimit[i] + 1e-6 ) );
            lfitNoise[i] -> SetParameter(1,0);
            cout<<"Guessed Parameter 0 fot "<<lFitNameOne<<" to be "<<lLeftPlusRightBgV0[i] * lHistoFullV0[i]->GetBinWidth(5) / (lRightBgLeftLimit[i]-lLeftBgRightLimit[i] + 1e-6 )<<endl;
            sprintf(lFitNameOne,"lSampleNoise%i",((int)(i)));

            //Define Function to Sample Background
            //cout<<"creating sample "<<i<<endl;
            lSampleNoise[i] = new TF1(lFitNameOne, this, &AliV0Module::MyBgPolToEval1,
                                      RoundToThousandth ( lLeftBgLeftLimit[i]   ),
                                      RoundToThousandth ( lRightBgRightLimit[i] ), 2 , "AliV0Module", "MyBgPolToEval1");
        }
        for(long i=0; i<fptbinnumb; i++) {
            //cout<<"Fitting function for bin "<<i<<", get name = "<<lfitNoise[i]->GetName()<<endl;
            sprintf(lFitNameOne,"lfitNoise%i",((int)(i)));
            lHistoFullV0[i] -> Fit(lFitNameOne,"LLrie+0");
            lSampleNoise[i]->SetParameter(0, lfitNoise[i]->GetParameter(0) );
            lSampleNoise[i]->SetParameter(1, lfitNoise[i]->GetParameter(1) );
        }
        for(long i=0; i<fptbinnumb; i++) {
            cout<<"Overriding Background info: Was "<<lLeftPlusRightBgV0[i]<<", is now "<<lSampleNoise[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoFullV0[i]->GetBinWidth(5)<<endl;
            lLeftPlusRightBgV0[i] = lSampleNoise[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoFullV0[i]->GetBinWidth(5);
        }
        cout<<"--------------- Fitting Finished! ----------------------"<<endl;
    }

    //=============================================================
    // Compute Signal + Sig to noise
    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
        //Signal computation
        lSigRealV0[ipoint] = lSigPlusCenterBgV0[ipoint] - lLeftPlusRightBgV0[ipoint];
        lSigErrRealV0[ipoint] = TMath::Sqrt(lSigPlusCenterBgV0[ipoint]+lLeftPlusRightBgV0[ipoint]);
        //Error Formula: Equivalent to Sqrt(S+B+B) = Sqrt(S+2B)

        //Attribute Raw Counts
        fHistPtRaw->SetBinContent(ipoint+1, lSigRealV0[ipoint]/((fptbinlimits[ipoint+1]-fptbinlimits[ipoint])*lNInelasticEvents) );
        fHistPtRaw->SetBinError(ipoint+1, lSigErrRealV0[ipoint]/((fptbinlimits[ipoint+1]-fptbinlimits[ipoint])*lNInelasticEvents) );

        //Signal-to-noise computation
        if( lLeftPlusRightBgV0[ipoint] != 0 ) {
            fHistSigToNoise->SetBinContent(ipoint+1, (lSigPlusCenterBgV0[ipoint] - lLeftPlusRightBgV0[ipoint])/lLeftPlusRightBgV0[ipoint] );
        } else {
            fHistSigToNoise->SetBinContent(ipoint+1, -1) ; //-1 means: no background
        }
    }
    //=============================================================

    cout<<"--------------- Extracted Signal -----------------------"<<endl;
    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
        cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tSignal: "<<lSigRealV0[ipoint]<<" +/- "<<lSigErrRealV0[ipoint]<<endl;
    }
    cout<<"--------------------------------------------------------"<<endl;

    //=============================================================
    // Preparations for MC loop
    //Extra info available only in MC
    Int_t lPID = 0;
    Float_t lNegTransvMomentumMC, lPosTransvMomentumMC;
    Int_t lPIDPositive = 0;
    Int_t lPIDNegative = 0;
    Int_t lPIDMother = -1;
    Float_t lPtMother = -1;
    Int_t lPrimaryStatus = 0;
    Int_t lPrimaryStatusMother = 0;
    //=============================================================
    cout<<endl;
    //Prepare containers for Geant3/fluka
    //Will only be used and saved if we're dealing with Lambdas
    TH1F *lProtonMomentum[100][3];
    // char lNameOne[100];
    for(Int_t ibin=0; ibin<100; ibin++) {
        for(int ifile = 0; ifile < 3; ifile++){
            sprintf(lNameOne,"lProtonMomentumBin%imc%i",ibin,ifile);
            lProtonMomentum[ibin][ifile] = new TH1F(lNameOne,"",800,0,20);
            if(fWhichParticle == "Lambda")     bindescription="#Lambda, bin #";
            if(fWhichParticle == "AntiLambda") bindescription="#bar{#Lambda}, bin #";
            if(fWhichParticle == "K0Short")    bindescription="K^{0}_{S}, bin #";
            bindescription.Append(IntToString( ibin ));
            if ( ibin < fptbinnumb ) {
                bindescription.Append(", ");
                bindescription.Append(DoubleToString(fptbinlimits[ ibin ]));
                bindescription.Append("-");
                bindescription.Append(DoubleToString(fptbinlimits[ ibin+1 ]));
                bindescription.Append("GeV/c");
            }
            lProtonMomentum[ibin][ifile]->SetTitle(bindescription);
        }
    }

    //MC File Acquisition=============================================
    Double_t lEfficiencyStep[100][3];
    Double_t lEfficiencyErrorStep[100][3];
    Double_t lPureEfficiencyStep[100][3];
    Double_t lPureEfficiencyErrorStep[100][3];
    Double_t lEfficiency[100];
    Double_t lEfficiencyError[100];
    Double_t lPureEfficiency[100];
    Double_t lPureEfficiencyError[100];
    TString filename;
    Long_t lNEventsMC;
    //== Variables for holding peak position, width ==============
    Double_t lPeakPositionMC[100];
    Double_t lPeakWidthMC[100];

    char fgausnameMC[100];
    TF1 *fgausPtMC[100];


    TH1F* fHistEfficiencyParz[3];
    for (int ifile = 0; ifile<3; ifile++){
        fHistEfficiencyParz[ifile] = new TH1F(Form("fHistEfficiencyParz%i",ifile)    , "G3/G4 corrected efficiency;p_{T} (GeV/c);Efficiency"          , fptbinnumb, fptbinlimits);
    }
    TH1F* fHistPureEfficiencyParz[3];
    for (int ifile = 0; ifile<3; ifile++){
        fHistPureEfficiencyParz[ifile] = new TH1F(Form("fHistPureEfficiencyParz%i",ifile)    , "Pure efficiency;p_{T} (GeV/c);Efficiency"          , fptbinnumb, fptbinlimits);
    }

    TH3D* h3DGenerated = 0x0;
    TH1F *fHistMCCountbyptV0 = 0x0;
    TH1D* fHistGen = 0x0;
    //If it's Lambda or AntiLambda, do Geant3/fluka correction
    //Keep fit function outside 'if' scope
    TF1 *fitGeant3FlukaCorr = 0x0;
    TH1D* h1dPanosCorrections =0x0;
    TFile *fileCorrGeanT3FlukaPanos = new TFile("antiPrCorrFunc.root","READ");
    fitGeant3FlukaCorr = (TF1*)fileCorrGeanT3FlukaPanos->Get("funcCorrAntiProtonFLUKA");

    Bool_t fEvSel_AllSelections;

    if (kDoEfficiency) {

        for (int ifile = 0; ifile < 3; ifile ++){ //0=15f, 1=17j, 2=18i

            //put all to 0 to be clean from iteration before
            for (int i = 0; i < 100; i++){
                lEfficiencyStep[i][ifile] = 0;
                lEfficiencyErrorStep[i][ifile]=0;
                lPureEfficiencyStep[i][ifile] = 0;
                lPureEfficiencyErrorStep[i][ifile]=0;
                lNEventsMC = 0;
                //
                //lHistoFullV0MC[i]->Reset();
                //f2dHistPtResolution->Reset();
                lLeftPlusRightBgV0MC[i] = 0.;
                lLeftPlusRightBgV0MC[i] = 0.;
                lSigPlusCenterBgV0MC[i] = 0.;
                //if (ifile>0) lProtonMomentum[i][ifile]->Reset();
                //lHistResolution[i]->Reset();
                //lHistPtDistr[i]->Reset();
                //lHistPtDistrMCtruth[i]->Reset();
                //
                lPeakPositionMC[i]= 0;
                lPeakWidthMC[i]= 0;
                //
                lSigMCV0[i] = 0;
                lSigErrMCV0[i] = 0;
                lSigMCV0_pTCorr[i] = 0;
                lSigErrMCV0_pTCorr[i] = 0;
            }

            for(Long_t n1 = 0; n1<250; n1++){
                lLeftPlusRightBgV0MC_pTCorr[n1]=0;
                lSigPlusCenterBgV0MC_pTCorr[n1]=0;
            }

            if(ifile == 0) filename =  fMCDataFile0;
            if(ifile == 1) filename =  fMCDataFile1;
            if(ifile == 2) filename =  fMCDataFile2;

            cout<<"--------------- Opening MC file " << ifile+1 << " : " << filename.Data() << " ------------------------"<<endl;
            TFile* fileMC;
            TList* v0listMC;
            TTree *lTreeMC;
            TTree* lTreeEventMC;
            fileMC = TFile::Open(filename, "READ");
            if (ifile==0) {
                v0listMC      = (TList*)fileMC->Get("PWGLF_StrVsMult_MC/cList");
                lTreeMC      = (TTree*)fileMC->Get("PWGLF_StrVsMult_MC/fTreeV0");
                lTreeEventMC = (TTree*)fileMC->Get("PWGLF_StrVsMult_MC/fTreeEvent");
            } else {
                v0listMC      = (TList*)fileMC->Get("PWGLF_StrVsMult/cList");
                lTreeMC      = (TTree*)fileMC->Get("PWGLF_StrVsMult/fTreeV0");
                lTreeEventMC = (TTree*)fileMC->Get("PWGLF_StrVsMult/fTreeEvent");
            }

            Long_t lNCandidatesMC = lTreeMC->GetEntries();
            cout<<"--------------------------------------------------------"<<endl;
            //================================================================

            //Linking to Tree=================================================
            //--- Base Variables ----------------------------------------------

            lTreeMC->SetBranchAddress("fTreeVariablePosTrackStatus", &lPosTrackStatus);
            lTreeMC->SetBranchAddress("fTreeVariableNegTrackStatus", &lNegTrackStatus);
            lTreeMC->SetBranchAddress("fTreeVariableNegTOFBCid",&fTreeVariableNegTOFBCid);
            lTreeMC->SetBranchAddress("fTreeVariablePosTOFBCid",&fTreeVariablePosTOFBCid);
            lTreeMC->SetBranchAddress("fTreeVariablePosEta",&lPosEta);
            lTreeMC->SetBranchAddress("fTreeVariableNegEta",&lNegEta);
            lTreeMC->SetBranchAddress("fTreeVariablePrimaryStatus",&lPrimaryStatus);
            lTreeMC->SetBranchAddress("fTreeVariablePrimaryStatusMother",&lPrimaryStatusMother);
            lTreeMC->SetBranchAddress("fTreeVariablePt",&lPt);
            lTreeMC->SetBranchAddress("fTreeVariablePtMother",&lPtMother);
            lTreeMC->SetBranchAddress("fTreeVariablePtMC",&lPtMC);
            if ( fWhichParticle == "Lambda"      )  lTreeMC->SetBranchAddress("fTreeVariableInvMassLambda",&lInvariantMass);
            if ( fWhichParticle == "AntiLambda"  ) 	lTreeMC->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMass);
            if ( fWhichParticle == "K0Short"     ) 	lTreeMC->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMass);
            if ( fWhichParticle == "Lambda"      ) { //For symmetry of computation...
                lTreeMC->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
                lInvariantMassCompetingTwo = -1; //will always be far enough
            }
            if ( fWhichParticle == "AntiLambda"  ) { //For symmetry of computation...
                lTreeMC->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
                lInvariantMassCompetingTwo = -1; //will always be far enough
            }
            if ( fWhichParticle == "K0Short"     ) {
                lTreeMC->SetBranchAddress("fTreeVariableInvMassLambda"    ,&lInvariantMassCompetingOne);
                lTreeMC->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMassCompetingTwo);
            }
            lTreeMC->SetBranchAddress("fTreeVariableRapMC",&lRap);
            lTreeMC->SetBranchAddress("fTreeVariableLeastNbrCrossedRows",&lLeastNbrCrossedRows);
            lTreeMC->SetBranchAddress("fTreeVariableLeastRatioCrossedRowsOverFindable",&lLeastNbrCrossedRowsOverFindable);
            //--- Topological selection variables -----------------------------
            lTreeMC->SetBranchAddress("fTreeVariableV0Radius",&lV0Radius);
            lTreeMC->SetBranchAddress("fTreeVariableDcaNegToPrimVertex",&lDcaNegToPrimVertex);
            lTreeMC->SetBranchAddress("fTreeVariableDcaPosToPrimVertex",&lDcaPosToPrimVertex);
            lTreeMC->SetBranchAddress("fTreeVariableDcaV0Daughters",&lDcaV0Daughters);
            lTreeMC->SetBranchAddress("fTreeVariableDcaV0ToPrimVertex",&lDCAV0ToPrimVertex);
            lTreeMC->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle",&lV0CosinePointingAngle);
            lTreeMC->SetBranchAddress("fTreeVariableDistOverTotMom",&lDistOverTotMom);
            //---- Extra  Information Available in MC only --------------------
            lTreeMC->SetBranchAddress("fTreeVariablePID",&lPID);
            lTreeMC->SetBranchAddress("fTreeVariablePIDMother",&lPIDMother);
            lTreeMC->SetBranchAddress("fTreeVariablePIDPositive",&lPIDPositive);
            lTreeMC->SetBranchAddress("fTreeVariablePIDNegative",&lPIDNegative);
            //--- Multiplicity ------------------------------------------------
            lTreeMC->SetBranchAddress(Form("fTreeVariableCentrality_%s","V0M"), &fMultCentrality);
            lTreeMC->SetBranchAddress(Form("fTreeVariableCentrality_%s",fWhichEffEnergyEstimator.Data()), &fEnergyCentrality);
            lTreeMC->SetBranchAddress("fTreeVariablePosInnerP",&fTreeVariablePtLegPos);
            lTreeMC->SetBranchAddress("fTreeVariableNegInnerP",&fTreeVariablePtLegNeg);

            //--- Armenteros-Podolansky ----------------------------------------
            // lTreeMC->SetBranchAddress("fTreeVariableAlphaV0",&lArmAlpha);
            // lTreeMC->SetBranchAddress("fTreeVariablePtArmV0",&lArmPt);
            //================================================================
            // PID
            //--- TPC dEdx Variables ------------------------------------------
            lTreeMC->SetBranchAddress("fTreeVariableNSigmasPosProton",&lNSigmasPosProton);
            lTreeMC->SetBranchAddress("fTreeVariableNSigmasNegProton",&lNSigmasNegProton);
            lTreeMC->SetBranchAddress("fTreeVariableNSigmasPosPion",&lNSigmasPosPion);
            lTreeMC->SetBranchAddress("fTreeVariableNSigmasNegPion",&lNSigmasNegPion);
            // lTreeMC->SetBranchAddress("fTreeVariablePIDTrackingPos",&pidFromTrackPos);
            // lTreeMC->SetBranchAddress("fTreeVariablePIDTrackingNeg",&pidFromTrackNeg);
            lTreeMC->SetBranchAddress("fTreeVariablePosPx"           , &lPosPx);
            lTreeMC->SetBranchAddress("fTreeVariablePosPy"           , &lPosPy);
            lTreeMC->SetBranchAddress("fTreeVariablePosPz"           , &lPosPz);
            lTreeMC->SetBranchAddress("fTreeVariableNegPx"           , &lNegPx);
            lTreeMC->SetBranchAddress("fTreeVariableNegPy"           , &lNegPy);
            lTreeMC->SetBranchAddress("fTreeVariableNegPz"           , &lNegPz);
            lTreeMC->SetBranchAddress("fTreeVariableEvSel_AllSelections"  , &fEvSel_AllSelections  );
            lTreeMC->SetBranchAddress("fTreeVariableRun"  , &fRun  );

            //================================================================
            cout<<endl;
            cout<<"--------------- MC Data File Loop  ---------------------"<<endl;
            Long_t lOneTenthOfNCandidatesMC = ((double)(lNCandidatesMC) / 10. );
            for(Long_t icand = 0; icand<lNCandidatesMC; icand++) {
                lTreeMC->GetEntry(icand);

                if(fEvSel_AllSelections==0) continue;

                // Multiplicity Switch
                lMultiplicity = (Double_t)fMultCentrality;

                if( icand % lOneTenthOfNCandidatesMC == 0 )
                    cout<<" Currently at candidate........: "<<icand<<" / "<<lNCandidatesMC<<" ( "<<(long)(((double)(icand)/(double)(lNCandidatesMC))*(100.+1e-3))<<"% )"<<endl;
                if(TMath::Abs(lRap)<fRapidityBoundary &&
                        // ((!fPosRap && !fNegRap) || (fPosRap && lRap > 0. && lNegEta > 0. && lPosEta >0) || (fNegRap && lRap<0. && lNegEta < 0. && lPosEta <0)) &&
                        TMath::Abs(lNegEta)       <= fCutDaughterEta               &&
                        TMath::Abs(lPosEta)       <= fCutDaughterEta               &&
                        lV0Radius                 >= fCutV0Radius                  &&
                        lDcaNegToPrimVertex       >= fCutDCANegToPV                &&
                        lDcaPosToPrimVertex       >= fCutDCAPosToPV                &&
                        lDcaV0Daughters           <= fCutDCAV0Daughters            &&
                        lV0CosinePointingAngle    >= fCutV0CosPA                   &&
                        lParticleMass*lDistOverTotMom    <= fCutProperLifetime     &&
                        lLeastNbrCrossedRows             >= fCutLeastNumberOfCrossedRows             &&
                        lLeastNbrCrossedRowsOverFindable >= fCutLeastNumberOfCrossedRowsOverFindable &&
                        TMath::Abs(lInvariantMassCompetingOne - lCompetingParticleMass) > fCutCompetingV0Rejection &&
                        TMath::Abs(lInvariantMassCompetingTwo - lCompetingParticleMass) > fCutCompetingV0Rejection &&
                        // ( fSpecialArmenterosCutK0s == kFALSE || ( fSpecialArmenterosCutK0s == kTRUE && lArmPt*5>TMath::Abs(lArmAlpha) ) ) &&
                        ( //perfect PID association, IsPhysicalPrimary association
                            ( fWhichParticle == "Lambda"
                            && lPID           == 3122 //V0 is a Lambda
                            && lPIDPositive   == 2212 //Pos Daughter is p
                            && lPIDNegative   == -211 //Neg Daughter is pi-
                            &&lPrimaryStatus==1 //Physical Primary Criterion
                            ) ||
                            ( fWhichParticle == "AntiLambda"
                            && lPID           == -3122  //V0 is an AntiLambda
                            && lPIDPositive   == 211    //Pos Daughter is pi+
                            && lPIDNegative   == -2212  //Neg Daughter is antiproton
                            && lPrimaryStatus==1 //Physical Primary Criterion
                            ) ||
                            ( fWhichParticle == "K0Short"
                            && lPID           == 310   //V0 is an AntiLambda
                            && lPIDPositive   == 211   //Pos Daughter is pi+
                            && lPIDNegative   == -211  //Neg Daughter is pi-
                            && lPrimaryStatus==1 //Physical Primary Criterion
                            ) )
                        && //OOB condition
                        CheckITSTOF(lPosTrackStatus,  lNegTrackStatus,  fTreeVariablePosTOFBCid,  fTreeVariableNegTOFBCid)
                        //CheckITSTOFOne( kFALSE, lPt, 2., lPosTrackStatus,  lNegTrackStatus,  fTreeVariablePosTOFBCid,  fTreeVariableNegTOFBCid)
                        &&
                        (
                        (fWhichParticle == "Lambda" && fTreeVariablePtLegPos > fPtLegProtCut  ) ||
                        (fWhichParticle == "AntiLambda" && fTreeVariablePtLegNeg > fPtLegProtCut ) ||
                        (fWhichParticle == "K0Short")
                        )
                    ) { // Start Entry Loop

                    lWeAreAtBin = fHistPt->FindBin(lPt)-1;
                    if(lWeAreAtBin == -1) lWeAreAtBin = 99; //UnderFlow, special treatment

                    lWeAreAtBin_pTCorr = fHistReco[ifile]->FindBin( lPt ) - 1;
                    if(lWeAreAtBin_pTCorr == -1) lWeAreAtBin_pTCorr = 250; //UnderFlow, special treatment

                    //1 - IsPhysicalPrimary Stuff
                    if ( lPrimaryStatus == 1) {
                        lHistoFullV0MC[lWeAreAtBin]->Fill(lInvariantMass); //fill with specific inv mass
                        f2dHistPtResolution->Fill(lPt,lPtMC);
                        //Extract left and right background areas and peak
                        //--- Left Background Sample
                        if(lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]    && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]   ) {
                            // lLeftPlusRightBgV0MC[lWeAreAtBin]++;
                            lLeftPlusRightBgV0MC[lWeAreAtBin] += 1.;//*weightMass;
                        }
                        if(lInvariantMass>lRightBgLeftLimit[lWeAreAtBin]   && lInvariantMass<lRightBgRightLimit[lWeAreAtBin]  ) {
                            // lLeftPlusRightBgV0MC[lWeAreAtBin]++;
                            lLeftPlusRightBgV0MC[lWeAreAtBin] += 1.;//*weightMass;
                            lLeftPlusRightBgV0MC_pTCorr[lWeAreAtBin_pTCorr]++;//*weightMass;
                        }
                        //--- Peak Region
                        if(lInvariantMass>lPeakLeftLimit[lWeAreAtBin]      && lInvariantMass<lPeakRightLimit[lWeAreAtBin]     ) {
                            // lSigPlusCenterBgV0MC[lWeAreAtBin]++;
                            lSigPlusCenterBgV0MC[lWeAreAtBin] += 1.; // *weightMass;
                            lSigPlusCenterBgV0MC_pTCorr[lWeAreAtBin_pTCorr]++;
                        }

                        //--- Info may be needed for geant3/fluka
                        if ( fWhichParticle == "Lambda"     ) lProtonMomentum[lWeAreAtBin][ifile]->Fill( TMath::Sqrt( lPosPx*lPosPx + lPosPy*lPosPy ) );
                        if ( fWhichParticle == "AntiLambda" ) lProtonMomentum[lWeAreAtBin][ifile]->Fill( TMath::Sqrt( lNegPx*lNegPx + lNegPy*lNegPy ) );

                        //--- Resolution tests
                        lHistResolution[lWeAreAtBin]->Fill(lPt - lPtMC);
                        lHistPtDistr[lWeAreAtBin]->Fill(lPt);
                        lHistPtDistrMCtruth[lWeAreAtBin]->Fill(lPtMC);
                        //
                        if(fFillQA){
                        lHistoV0RadiusMC->Fill(lV0Radius);
                        lHistDCAnegToPrimMC->Fill(lDcaNegToPrimVertex);
                        lHistDCAposToPrimMC->Fill(lDcaPosToPrimVertex);
                        lHistDCAV0DaughtersMC->Fill(lDcaV0Daughters);
                        lProperLifeTimeMC->Fill(lParticleMass*lDistOverTotMom);
                        lCosPointingAngleMC->Fill(lV0CosinePointingAngle);
                        if(fWhichParticle == "Lambda") {
                                lHistoNsigmaTPCposMC->Fill(lPt, lNSigmasPosProton); lHistoNsigmaTPCnegMC->Fill(lPt, lNSigmasNegPion);
                                }
                        if(fWhichParticle == "AntiLambda") {
                                lHistoNsigmaTPCposMC->Fill(lPt, lNSigmasPosPion); lHistoNsigmaTPCnegMC->Fill(lPt, lNSigmasNegProton);
                                }
                        if(fWhichParticle == "K0Short") {
                                lHistoNsigmaTPCposMC->Fill(lPt, lNSigmasPosPion); lHistoNsigmaTPCnegMC->Fill(lPt, lNSigmasNegPion);
                                }
                        }

                }

                } // End Entry Loop
            }
            cout<<"--------------- Loop Completed -------------------------"<<endl;
            cout<<endl;

            cout<<"--------------- X-Check Peak Finding (gauss+linear) ----"<<endl;
            for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
                cout<<"---> Peak Finding, bin #"<<ibin<<" (perfect count = "<<lHistoFullV0MC[ibin]->GetEntries()<<")..."<<endl;
                sprintf(fgausnameMC,"fGausPtMC%i",ibin);
                if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ) {
                    // fgausPtMC[ibin]= new TF1(fgausnameMC,"[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]",1.116-0.040,1.116+0.040);
                    fgausPtMC[ibin]= new TF1(fgausnameMC,"[0]*TMath::Gaus(x,[1],[2])",1.116-0.040,1.116+0.040);
                    fgausPtMC[ibin]->SetParameter(1,1.116);
                    fgausPtMC[ibin]->SetParameter(2,0.0025);
                    fgausPtMC[ibin]->SetParLimits(1,1.1,1.2);
                    fgausPtMC[ibin]->SetParLimits(2,0.001,0.02);
                }
                if ( fWhichParticle == "K0Short") {
                    fgausPtMC[ibin]= new TF1(fgausnameMC,"[0]*TMath::Gaus(x,[1],[2])",0.498-0.06,0.498+0.060);
                    fgausPtMC[ibin]->SetParameter(1,0.498);
                    fgausPtMC[ibin]->SetParameter(2,0.004);
                }
                fgausPtMC[ibin]->SetParameter(0,lHistoFullV0MC[ibin]->GetMaximum() * 0.9);
                if(fUseIntegratedEfficiencies) lHistoFullV0MC[ibin]->Fit(fgausnameMC,"QREM0");
                else lHistoFullV0MC[ibin]->Fit(fgausnameMC);
                lPeakPositionMC[ibin] = TMath::Abs( fgausPtMC[ibin]->GetParameter(1) );
                lPeakWidthMC[ibin] = fgausPtMC[ibin]->GetParameter(2);
                cout<<"---> ["<<fptbinlimits[ibin]<<" - "<<fptbinlimits[ibin+1]<<" GeV/c]\tPeak at: "<<lPeakPositionMC[ibin]<<", sigma = "<<lPeakWidthMC[ibin]<<endl;
                fHistPeakPositionMC->SetBinContent(ibin+1, lPeakPositionMC[ibin]);
                fHistPeakPositionMC->SetBinError(ibin+1, lPeakWidthMC[ibin]);

                //Create appropriate TLine Objects for Canvases
                lLineLeftMostMC[ibin]  = new TLine( lLeftBgLeftLimit[ibin],   0, lLeftBgLeftLimit[ibin],   lHistoFullV0MC[ibin]->GetMaximum() * 0.95 );
                lLineLeftMC[ibin]      = new TLine( lLeftBgRightLimit[ibin],  0, lLeftBgRightLimit[ibin],  lHistoFullV0MC[ibin]->GetMaximum() * 0.95 );
                lLineRightMC[ibin]     = new TLine( lRightBgLeftLimit[ibin],  0, lRightBgLeftLimit[ibin],  lHistoFullV0MC[ibin]->GetMaximum() * 0.95 );
                lLineRightMostMC[ibin] = new TLine( lRightBgRightLimit[ibin], 0, lRightBgRightLimit[ibin], lHistoFullV0MC[ibin]->GetMaximum() * 0.95 );

                //Preparing Canvas for storing...
                lCanvasHistoFullV0MC[ibin]->cd();
                lHistoFullV0MC[ibin]->Draw();
                lLineLeftMostMC[ibin]->Draw();
                lLineLeftMC[ibin]->Draw();
                lLineRightMC[ibin]->Draw();
                lLineRightMostMC[ibin]->Draw();
                fgausPtMC[ibin]->Draw("same");

            }
            cout<<"--------------- Peak Finding Finished (in MC) ----------"<<endl;
            cout<<endl;

            TF1 *lfitNoiseMC[50];
            TF1 *lSampleNoiseMC[50];
            char lFitNameOneMC[50];

            if( fFitBackgroundSwitch ) {
                cout<<"--------------- Backgrounds: fitting with linear -------"<<endl;
                for(long i=0; i<fptbinnumb; i++) {
                    //Define Function to Fit Background
                    sprintf(lFitNameOneMC,"lfitNoiseMC%i",(int)i);
                    cout<<"creating fitnoise, named "<<lFitNameOneMC<<endl;
                    lfitNoiseMC[i] = new TF1(lFitNameOneMC, this, &AliV0Module::MyBgPol1,
                                            RoundToThousandth ( lLeftBgLeftLimit[i]   ),
                                            RoundToThousandth ( lRightBgRightLimit[i] ), 4 , "AliV0Module", "MyBgPol1");
                    lfitNoiseMC[i] -> FixParameter(2,RoundToThousandth ( lLeftBgRightLimit[i] ) );
                    lfitNoiseMC[i] -> FixParameter(3,RoundToThousandth ( lRightBgLeftLimit[i] ) );
                    lfitNoiseMC[i] -> SetParameter(0,lLeftPlusRightBgV0MC[i]*lHistoFullV0MC[i]->GetMaximum() / (lSigPlusCenterBgV0MC[i]+1e-6) );
                    lfitNoiseMC[i] -> SetParameter(1,0);
                    sprintf(lFitNameOneMC,"lSampleNoiseMC%i",(int)i);

                    //Define Function to Sample Background
                    cout<<"creating sample "<<i<<endl;
                    lSampleNoiseMC[i] = new TF1(lFitNameOneMC, this, &AliV0Module::MyBgPolToEval1,
                                                RoundToThousandth ( lLeftBgLeftLimit[i]   ),
                                                RoundToThousandth ( lRightBgRightLimit[i] ), 2, "AliV0Module", "MyBgPolToEval1" );
                }
                for(long i=0; i<fptbinnumb; i++) {
                    cout<<"Fitting function for bin "<<i<<", get name = "<<lfitNoiseMC[i]->GetName()<<endl;
                    sprintf(lFitNameOneMC,"lfitNoiseMC%i",(int)i);
                    lHistoFullV0MC[i] -> Fit(lFitNameOneMC,"LLrie+0");
                    lSampleNoiseMC[i]->SetParameter(0, lfitNoiseMC[i]->GetParameter(0) );
                    lSampleNoiseMC[i]->SetParameter(1, lfitNoiseMC[i]->GetParameter(1) );
                }
                for(long i=0; i<fptbinnumb; i++) {
                    cout<<"Overriding Background info: Was "<<lLeftPlusRightBgV0MC[i]<<", is now "<<lSampleNoiseMC[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoFullV0MC[i]->GetBinWidth(5)<<endl;
                    lLeftPlusRightBgV0MC[i] = lSampleNoiseMC[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoFullV0MC[i]->GetBinWidth(5);
                }
                cout<<"--------------- Fitting Finished! ----------------------"<<endl;
            }

            //=============================================================
            // Compute Signal
            //MC pT shape Correction
            for(Int_t ipoint = 0; ipoint<250; ipoint++){
            //MC pT Corr
                lSigMCV0_pTCorr[ipoint] = lSigPlusCenterBgV0MC_pTCorr[ipoint] - lLeftPlusRightBgV0MC_pTCorr[ipoint];
                lSigErrMCV0_pTCorr[ipoint] = TMath::Sqrt(lSigPlusCenterBgV0MC_pTCorr[ipoint]+lLeftPlusRightBgV0MC_pTCorr[ipoint]);

                fHistReco[ifile]->SetBinContent(ipoint+1, lSigMCV0_pTCorr[ipoint]);
            }

            //=============================================================
            // Compute Signal
            for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
                lSigMCV0[ipoint] = lSigPlusCenterBgV0MC[ipoint] - lLeftPlusRightBgV0MC[ipoint];
                lSigErrMCV0[ipoint] = TMath::Sqrt(lSigPlusCenterBgV0MC[ipoint]+lLeftPlusRightBgV0MC[ipoint]);
                //Error Formula: Equivalent to Sqrt(S+B+B) = Sqrt(S+2B)

                //Signal-to-noise computation
                if( lLeftPlusRightBgV0MC[ipoint] != 0 ) {
                    fHistSigToNoiseMC->SetBinContent(ipoint+1, (lSigPlusCenterBgV0MC[ipoint] - lLeftPlusRightBgV0MC[ipoint])/lLeftPlusRightBgV0MC[ipoint] );
                } else {
                    fHistSigToNoiseMC->SetBinContent(ipoint+1, -1) ; //-1 means: no background
                }
            }
            //=============================================================
            cout<<"--------------- Extracted Signal (MC) ------------------"<<endl;
            for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
                cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tSignal: "<<lSigMCV0[ipoint]<<" +/- "<<lSigErrMCV0[ipoint]<<endl;
            }
            cout<<"--------------------------------------------------------"<<endl;
            cout<<endl;


            cout<<"--------------- Process Generated V0s ------------------"<<endl;
            //Filling histogram with original MC particles====================
            TH1F* fHistGen = 0x0;
            TH3F* fHistGen2D = 0x0;
            cout<<"Getting Histogram..."<<endl;

            fHistGen2D  = (TH3F*)v0listMC->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s",fWhichParticle.Data()));

            Int_t rapBin_min = fHistGen2D->GetYaxis()->FindBin( -0.5+1.e-6 );
            Int_t rapBin_max = fHistGen2D->GetYaxis()->FindBin( 0.5-1.e-6 );
            //
            Int_t multBin_min = fHistGen2D->GetZaxis()->FindBin( 0.+1.e-6 );
            Int_t multBin_max = fHistGen2D->GetZaxis()->FindBin( 100.-1.e-6 );
            //fHistGen
            fHistGen = (TH1F*)fHistGen2D->ProjectionX("fHistPt_Gen", rapBin_min, rapBin_max, multBin_min, multBin_max);

            TH1F *fHistMCCountbyptV0	= new TH1F("fHistMCCountbyptV0","V0 MC count;p_{T} (GeV/c);Counts",fptbinnumb,fptbinlimits);

            Double_t temppt; //  weight; // Int_t binMCw = -1;
            if(!fHistGen) { printf("non esiste1!!!! \n"); getchar();}
            for(long i = 1; i<fHistGen->GetNbinsX()+1; i++) {

                temppt = fHistGen->GetXaxis()->GetBinCenter(i);
                Double_t maxFillBin = fHistGen->GetBinContent(i);
                if(fPosRap || fNegRap) maxFillBin = fHistGen->GetBinContent(i)/2.;
                for(long filling = 0; filling<maxFillBin; filling++) {
                    fHistMCCountbyptV0->Fill(temppt);
                    fHistMCtruth->Fill(temppt);
                }
            }

            fHistMCtruth->Scale(1./nEventsMC);

            cout<<"--------------- Generated V0 Dump ----------------------"<<endl;
            for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
                cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tSignal: "<<fHistMCCountbyptV0->GetBinContent(ipoint+1)<<" +/- "<<TMath::Sqrt(fHistMCCountbyptV0->GetBinContent(ipoint+1))<<endl;
            }
            cout<<"--------------------------------------------------------"<<endl;

            //=============================================================
            // Compute Efficiency
            for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
                Double_t lEffNumerator = lSigMCV0[ipoint];
                Double_t lEffDenominator = fHistMCCountbyptV0->GetBinContent(ipoint+1);
                lEfficiencyStep[ipoint][ifile] = lEffNumerator / lEffDenominator;
                fHistEffNumerator->SetBinContent(ipoint+1,lEffNumerator);
                fHistEffDenominator->SetBinContent(ipoint+1,lEffDenominator);
                lEfficiencyErrorStep[ipoint][ifile] = ErrorInRatio(
                                            lSigMCV0[ipoint],
                                            lSigErrMCV0[ipoint],
                                            fHistMCCountbyptV0->GetBinContent(ipoint+1),
                                            TMath::Sqrt(fHistMCCountbyptV0->GetBinContent(ipoint+1))
                                        );

                lPureEfficiencyStep[ipoint][ifile] = lEfficiencyStep[ipoint][ifile];
                lPureEfficiencyErrorStep[ipoint][ifile] = lEfficiencyErrorStep[ipoint][ifile];
            }

            for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
                cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tEff.: "<<lEfficiencyStep[ipoint][ifile] <<" +/- "<<lEfficiencyErrorStep[ipoint][ifile] <<endl;
                fHistPureEfficiencyParz[ifile]->SetBinContent(ipoint+1, lEfficiencyStep[ipoint][ifile]);
                fHistPureEfficiencyParz[ifile]->SetBinError(ipoint+1, lEfficiencyErrorStep[ipoint][ifile]);
            }

            //=============================================================
            cout<<endl;
            cout<<"--------------- Memory Cleanup -------------------------"<<endl;
            //fHistGen->Delete();
            fHistMCCountbyptV0->Delete();
            v0listMC->Delete();
            fileMC->Close("R");
            fileMC->Delete();
            delete fileMC;
            cout<<endl;

            if( (fWhichParticle == "AntiLambda") && fG3FlukaCorrOn && ifile == 0) {

                cout<<"--------------- Geant3/Fluka Correction ----------------"<<endl;

                Printf("Test fit function...");
                Printf(" - fitGeant3FlukaCorr for pt =  .25 GeV/c : %f", fitGeant3FlukaCorr->Eval( .25) );
                Printf(" - fitGeant3FlukaCorr for pt =  .5  GeV/c : %f", fitGeant3FlukaCorr->Eval( .5) );
                Printf(" - fitGeant3FlukaCorr for pt =  1   GeV/c : %f", fitGeant3FlukaCorr->Eval( 1.0) );
                Printf(" - fitGeant3FlukaCorr for pt =  2   GeV/c : %f", fitGeant3FlukaCorr->Eval( 2.0) );
                Printf(" - fitGeant3FlukaCorr for pt =  5   GeV/c : %f", fitGeant3FlukaCorr->Eval( 5.0) );
                Printf(" - fitGeant3FlukaCorr for pt =  7   GeV/c : %f", fitGeant3FlukaCorr->Eval( 7.0) );
                Printf(" - fitGeant3FlukaCorr for pt =  8   GeV/c : %f", fitGeant3FlukaCorr->Eval( 8.0) );
                Printf(" - fitGeant3FlukaCorr for pt = 10   GeV/c : %f", fitGeant3FlukaCorr->Eval(10.0) );

                cout<<"--------------- Embedding G3/F in Efficiencies ---------"<<endl;
                for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
                    lEfficiencyStep[ipoint][ifile]       *= ( fitGeant3FlukaCorr->Eval( lProtonMomentum[ipoint][ifile]->GetMean() ) ) ;
                    lEfficiencyErrorStep[ipoint][ifile]  *= ( fitGeant3FlukaCorr->Eval( lProtonMomentum[ipoint][ifile]->GetMean() ) ) ;
                }
                cout<<"--------------- G3/F Corrected Efficiencies ------------"<<endl;
                for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
                    cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tEff.: "<<lEfficiencyStep[ipoint][ifile] <<" +/- "<<lEfficiencyErrorStep[ipoint][ifile] <<endl;
                    fHistEfficiencyParz[ifile]->SetBinContent(ipoint+1, lEfficiencyStep[ipoint][ifile]);
                    fHistEfficiencyParz[ifile]->SetBinError(ipoint+1, lEfficiencyErrorStep[ipoint][ifile]);
                }
                cout<<"--------------------------------------------------------"<<endl;
                cout<<endl;

            } //end Geant3/fluka correction if
            else {
                for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
                    cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tEff.: "<<lEfficiencyStep[ipoint][ifile] <<" +/- "<<lEfficiencyErrorStep[ipoint][ifile] <<endl;
                    fHistEfficiencyParz[ifile]->SetBinContent(ipoint+1, lEfficiencyStep[ipoint][ifile]);
                    fHistEfficiencyParz[ifile]->SetBinError(ipoint+1, lEfficiencyErrorStep[ipoint][ifile]);
                }
                cout<<"--------------------------------------------------------"<<endl;
                cout<<endl;
            }
        } // end loop on files

        for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
            //Do weighted mean of all efficiencies
            lEfficiency[ipoint] = DoWeightedMeanXBin(0, lEfficiencyStep[ipoint][0], lEfficiencyStep[ipoint][1], lEfficiencyStep[ipoint][2], lEfficiencyErrorStep[ipoint][0], lEfficiencyErrorStep[ipoint][1], lEfficiencyErrorStep[ipoint][2], fEventsWeights[0], fEventsWeights[1], fEventsWeights[2]);
            lEfficiencyError[ipoint] = DoWeightedMeanXBin(1, lEfficiencyStep[ipoint][0], lEfficiencyStep[ipoint][1], lEfficiencyStep[ipoint][2], lEfficiencyErrorStep[ipoint][0], lEfficiencyErrorStep[ipoint][1], lEfficiencyErrorStep[ipoint][2], fEventsWeights[0], fEventsWeights[1], fEventsWeights[2]);

            lPureEfficiency[ipoint] = DoWeightedMeanXBin(0, lPureEfficiencyStep[ipoint][0], lPureEfficiencyStep[ipoint][1], lPureEfficiencyStep[ipoint][2], lPureEfficiencyErrorStep[ipoint][0], lPureEfficiencyErrorStep[ipoint][1], lPureEfficiencyErrorStep[ipoint][2], fEventsWeights[0], fEventsWeights[1], fEventsWeights[2]);
            lPureEfficiencyError[ipoint] = DoWeightedMeanXBin(1, lPureEfficiencyStep[ipoint][0], lPureEfficiencyStep[ipoint][1], lPureEfficiencyStep[ipoint][2], lPureEfficiencyErrorStep[ipoint][0], lPureEfficiencyErrorStep[ipoint][1], lPureEfficiencyErrorStep[ipoint][2], fEventsWeights[0], fEventsWeights[1], fEventsWeights[2]);
        }
    } else {
        for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
            lEfficiency[ipoint] = 1.;
            lEfficiencyError[ipoint] = 0.00001;
            //
            lPureEfficiency[ipoint] = 1.;
            lPureEfficiencyError[ipoint] = 0.00001;
        }
    } // end if-else for reco efficiency

    fileCorrGeanT3FlukaPanos->Close();
    fileCorrGeanT3FlukaPanos->Delete();
    delete fileCorrGeanT3FlukaPanos;

    //=============================================================
    cout<<endl;
    cout<<"--------------- Memory Cleanup -------------------------"<<endl;

    cout<<"--------------- Pure Efficiency Numbers ----------------"<<endl;
    TH1F* fHistPureEfficiency 		= new TH1F("fHistPureEfficiency","Pure Efficiency;p_{T} (GeV/c);Pure Efficiency",fptbinnumb,fptbinlimits);
    TH1F* fHistEfficiency 		    = new TH1F("fHistEfficiency","Efficiency;p_{T} (GeV/c);Efficiency",fptbinnumb,fptbinlimits);
    if( fWhichParticle == "K0Short"      ) fHistPureEfficiency->SetTitle("K^{0}_{S} Efficiency");
    if( fWhichParticle == "Lambda"       ) fHistPureEfficiency->SetTitle("#Lambda Efficiency (no G3/F corr.)");
    if( fWhichParticle == "AntiLambda"   ) fHistPureEfficiency->SetTitle("#bar{#Lambda} Efficiency (no G3/F corr.)");

    if( fWhichParticle == "Lambda"       ) fHistEfficiency->SetTitle("#Lambda Efficiency (with G3/F corr.)");
    if( fWhichParticle == "AntiLambda"   ) fHistEfficiency->SetTitle("#bar{#Lambda} Efficiency (with G3/F corr.)");

    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
        cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tEff.: "<<lEfficiency[ipoint]<<" +/- "<<lEfficiencyError[ipoint]<<endl;
        fHistEfficiency->SetBinContent(ipoint+1, lEfficiency[ipoint]);
        fHistEfficiency->SetBinError(ipoint+1, lEfficiencyError[ipoint]);
        fHistPureEfficiency->SetBinContent(ipoint+1, lPureEfficiency[ipoint]);
        fHistPureEfficiency->SetBinError(ipoint+1, lPureEfficiencyError[ipoint]);
    }
    cout<<"--------------------------------------------------------"<<endl;

    cout<<endl;

    //Declaring outside Feeddown scope to keep in memory
    TH2F *f2dFeedDownMatrix[3];
    TH2F *f2dFeedDownEfficiency[3];
    TH2F *f2dFeedDownEfficiencyGFCorrected[3];
    for(int jj=0;jj<3;jj++){
        f2dFeedDownMatrix[jj] = 0x0;
        f2dFeedDownEfficiency[jj] = 0x0;
        f2dFeedDownEfficiencyGFCorrected[jj] = 0x0;
    }
    TH1F *fHistFeeddownSubtraction = 0x0;
    TH1F *XiSpectra = 0x0;
    TH2F* f2dFeedDownEfficiencyTOT = 0x0;
    TH2F* f2dFeedDownEfficiencyGFCorrectedTOT = 0x0;

   // === for systematics due to Xi spectra
    TH1F *feedDownSyst[2];
    TH1F *feedDownSystSub[2];

    for(int jj=0;jj<2;jj++){
        feedDownSyst[jj]  = new TH1F(Form("fHistFeeddownSubtraction_Syst%d",jj),"FD Subtraction;p_{T} (GeV/c);Ratio removed",fptbinnumb,fptbinlimits);
        feedDownSyst[jj]->SetDirectory(0);
        feedDownSystSub[jj]  = new TH1F(Form("fHistFeeddownSubtraction_SystSub%d",jj),"FD Subtraction;p_{T} (GeV/c);Ratio removed",fptbinnumb,fptbinlimits);
        feedDownSystSub[jj]->SetDirectory(0);
    }

    Double_t lFeedDownEfficiency[100][100][3];
    Double_t lFeedDownEfficiencyGFCorrected[100][100][3];
    Double_t lFeedDownEfficiencyError[100][100][3];
    Double_t lFeedDownEfficiencyGFCorrectedError[100][100][3];
    Double_t lFeedDownEfficiencyTOT[100][100];
    Double_t lFeedDownEfficiencyTOTError[100][100];
    Double_t lFeedDownMatrix[100][100];

    TFile *fileXiSpectra;
    TF1 *fLevyFit;
    TF1 *fLevyFitXiFeedDown = new TF1("LevyFitXiFeedDown",this ,&AliV0Module::MyLevyPtXi, 0.0,15,3, "AliV0Module","MyLevyPXi");
    fileXiSpectra = TFile::Open(fFeedDownDataFile,"READ");
    Double_t lFDParameter_dNdy;
    Double_t lFDParameter_Temp;
    Double_t lFDParameter_n;

    //=========================================================================
    // --- Feeddown correction section
    if( fFDSwitch != "NoFD" && fWhichParticle != "K0Short") {

        fLevyFit = (TF1 *)(fileXiSpectra->Get(Form("FinalSpectra/LevyFit%s_%.0f_%.0f_%s_%.0f_%.0f", fWhichMultEstimator.Data(), fLoMultBound, fHiMultBound, fWhichEffEnergyEstimator.Data(), fLoEEBound, fHiEEBound)));

        lFDParameter_dNdy = fLevyFit->GetParameter("norm");
        lFDParameter_Temp = fLevyFit->GetParameter("C");
        lFDParameter_n = fLevyFit->GetParameter("n");

        //////
        fLevyFitXiFeedDown->SetParameter(0, lFDParameter_dNdy);
        fLevyFitXiFeedDown->SetParameter(1, lFDParameter_Temp);
        fLevyFitXiFeedDown->SetParameter(2, lFDParameter_n);
        fLevyFitXiFeedDown->SetNpx(5000);

        cout<<"--------------- Feeddown Correction --------------------"<<endl;

        //Define Xi Binning
        //Result should be ~invariant with respect to this choice...
        Double_t xibinlimits[] = {
            0.00,  0.20,  0.40,  0.60,  0.80,  0.90,
            1.00,  1.10,  1.20,  1.30,  1.40,  1.50,
            1.70,  1.90,  2.20,  2.60,  3.10,  3.90,
            4.90,  6.00,  7.20,  8.50 ,10.00, 12.00
        };
        Long_t xibinnumb = sizeof(xibinlimits)/sizeof(Double_t) - 1;

        cout<<"--------------------------------------------------------"<<endl;
        cout<<" Feeddown Parameters used for Tsallis "<<endl;
        cout<<"--------------------------------------------------------"<<endl;
        cout<<" dN/dy ...............: "<<fLevyFitXiFeedDown->GetParameter(0)<<endl;
        cout<<" Temp ................: "<<fLevyFitXiFeedDown->GetParameter(1)<<endl;
        cout<<" n ...................: "<<fLevyFitXiFeedDown->GetParameter(2)<<endl;

        Double_t lProducedXi[100];

        //If you want me to double charged Xi feeddown, I'll do it here
        for(Int_t ixb = 0; ixb<xibinnumb; ixb++) {
        if(fFDSwitch == "DoubleChargedXi") lProducedXi[ixb] = lNInelasticEvents * (fLevyFitXiFeedDown->Integral(xibinlimits[ixb],xibinlimits[ixb+1]));
        else if(fFDSwitch == "UseMCRatio") lProducedXi[ixb] = lNInelasticEvents * (fLevyFitXiFeedDown->Integral(xibinlimits[ixb],xibinlimits[ixb+1])/2.);
        }
        if(fWhichParticle == "Lambda")
            cout<<"--------------- Generated Xi- Dump (real-corrected) ----"<<endl;
        if(fWhichParticle == "AntiLambda")
            cout<<"--------------- Generated Xi+ Dump (real-corrected) ----"<<endl;
        for(Int_t ixb = 0; ixb<xibinnumb; ixb++) {
            cout<<"Xi bin "<<ixb<<"\t"<<lProducedXi[ixb]<<endl;
        }
        cout<<"--------------------------------------------------------"<<endl;
        cout<<endl;

        for (int ifile = 0; ifile < 3; ifile ++){ //0=15f, 1=17j, 2=18i

            cout<<"--------------- Feeddown Correction --------------------"<<endl;
            //MC File Acquisition=============================================

            if(ifile == 0) filename =  fMCDataFile0;
            if(ifile == 1) filename =  fMCDataFile1;
            if(ifile == 2) filename =  fMCDataFile2;

            cout<<"--------------- Opening MC file " << ifile+1 << " : " << filename.Data() << " ------------------------"<<endl;

            TFile* fileMCFD = TFile::Open(filename, "READ");
            TList* v0listMCFD     = 0x0;
            TTree* lTreeMCFD      = 0x0;
            TTree* lTreeEventMCFD = 0x0;
            if (ifile == 0)
            {
                v0listMCFD = (TList *)fileMCFD->Get("PWGLF_StrVsMult_MC/cList");
                lTreeMCFD = (TTree *)fileMCFD->Get("PWGLF_StrVsMult_MC/fTreeV0");
                lTreeEventMCFD = (TTree *)fileMCFD->Get("PWGLF_StrVsMult_MC/fTreeEvent");
            }
            else
            {
                v0listMCFD = (TList *)fileMCFD->Get("PWGLF_StrVsMult/cList");
                lTreeMCFD = (TTree *)fileMCFD->Get("PWGLF_StrVsMult/fTreeV0");
                lTreeEventMCFD = (TTree *)fileMCFD->Get("PWGLF_StrVsMult/fTreeEvent");
            }

            Long_t lNCandidatesMCFD = lTreeMCFD->GetEntries();
            Long_t lNEventsMCFD = lTreeEventMCFD->GetEntries();
            cout<<" MC Events (before trig sel)...: "<<lNEventsMCFD<<endl;
            cout<<" Cascade MC Candidates.........: "<<lNCandidatesMCFD<<endl;
            cout<<"--------------------------------------------------------"<<endl;
            //================================================================

            //Linking to Tree=================================================
            //--- Base Variables ----------------------------------------------
            lTreeMCFD->SetBranchAddress("fTreeVariablePosTrackStatus", &lPosTrackStatus);
            lTreeMCFD->SetBranchAddress("fTreeVariableNegTrackStatus", &lNegTrackStatus);
            lTreeMCFD->SetBranchAddress("fTreeVariableNegTOFBCid",&fTreeVariableNegTOFBCid);
            lTreeMCFD->SetBranchAddress("fTreeVariablePosTOFBCid",&fTreeVariablePosTOFBCid);
            lTreeMCFD->SetBranchAddress("fTreeVariablePosEta",&lPosEta);
            lTreeMCFD->SetBranchAddress("fTreeVariableNegEta",&lNegEta);
            lTreeMCFD->SetBranchAddress("fTreeVariablePrimaryStatus",&lPrimaryStatus);
            lTreeMCFD->SetBranchAddress("fTreeVariablePrimaryStatusMother",&lPrimaryStatusMother);
            lTreeMCFD->SetBranchAddress("fTreeVariablePt",&lPt);
            lTreeMCFD->SetBranchAddress("fTreeVariablePtMC",&lPtMC);
            lTreeMCFD->SetBranchAddress("fTreeVariablePtMother",&lPtMother);
            if ( fWhichParticle == "Lambda"      )  lTreeMCFD->SetBranchAddress("fTreeVariableInvMassLambda",&lInvariantMass);
            if ( fWhichParticle == "AntiLambda"  ) 	lTreeMCFD->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMass);
            if ( fWhichParticle == "K0Short"     ) 	lTreeMCFD->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMass);
            if ( fWhichParticle == "Lambda"      ) { //For symmetry of computation...
                lTreeMCFD->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
                lInvariantMassCompetingTwo = -1.;
            }
            if ( fWhichParticle == "AntiLambda"  ) { //For symmetry of computation...
                lTreeMCFD->SetBranchAddress("fTreeVariableInvMassK0s",&lInvariantMassCompetingOne);
                lInvariantMassCompetingTwo = -1.;
            }
            if ( fWhichParticle == "K0Short"     ) {
                lTreeMCFD->SetBranchAddress("fTreeVariableInvMassLambda"    ,&lInvariantMassCompetingOne);
                lTreeMCFD->SetBranchAddress("fTreeVariableInvMassAntiLambda",&lInvariantMassCompetingTwo);
            }
            lTreeMCFD->SetBranchAddress("fTreeVariableRapMC",&lRap);
            lTreeMCFD->SetBranchAddress("fTreeVariableLeastNbrCrossedRows",&lLeastNbrCrossedRows);
            lTreeMCFD->SetBranchAddress("fTreeVariableLeastRatioCrossedRowsOverFindable",&lLeastNbrCrossedRowsOverFindable);
            //--- Topological selection variables -----------------------------
            lTreeMCFD->SetBranchAddress("fTreeVariableV0Radius",&lV0Radius);
            lTreeMCFD->SetBranchAddress("fTreeVariableDcaNegToPrimVertex",&lDcaNegToPrimVertex);
            lTreeMCFD->SetBranchAddress("fTreeVariableDcaPosToPrimVertex",&lDcaPosToPrimVertex);
            lTreeMCFD->SetBranchAddress("fTreeVariableDcaV0Daughters",&lDcaV0Daughters);
            lTreeMCFD->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle",&lV0CosinePointingAngle);
            lTreeMCFD->SetBranchAddress("fTreeVariableDistOverTotMom",&lDistOverTotMom);
            //---- Extra  Information Available in MC only --------------------
            lTreeMCFD->SetBranchAddress("fTreeVariableNegTransvMomentumMC",&lNegTransvMomentumMC);
            lTreeMCFD->SetBranchAddress("fTreeVariablePosTransvMomentumMC",&lPosTransvMomentumMC);
            lTreeMCFD->SetBranchAddress("fTreeVariablePID",&lPID);
            lTreeMCFD->SetBranchAddress("fTreeVariablePIDMother",&lPIDMother);
            lTreeMCFD->SetBranchAddress("fTreeVariablePIDPositive",&lPIDPositive);
            lTreeMCFD->SetBranchAddress("fTreeVariablePIDNegative",&lPIDNegative);
            //---- Multiplicity info -------------------------------------------
            lTreeMCFD->SetBranchAddress("fTreeVariableCentV0M",&lMultiplicity);
            //--- Armenteros-Podolansky ----------------------------------------
            lTreeMCFD->SetBranchAddress("fTreeVariableAlphaV0",&lArmAlpha);
            lTreeMCFD->SetBranchAddress("fTreeVariablePtArmV0",&lArmPt);
            lTreeMCFD->SetBranchAddress("fTreeVariablePosInnerP",&fTreeVariablePtLegPos);
            lTreeMCFD->SetBranchAddress("fTreeVariableNegInnerP",&fTreeVariablePtLegNeg);
            //--- TPC dEdx Variables ------------------------------------------
            lTreeMCFD->SetBranchAddress("fTreeVariableNSigmasPosProton",&lNSigmasPosProton);
            lTreeMCFD->SetBranchAddress("fTreeVariableNSigmasNegProton",&lNSigmasNegProton);
            lTreeMCFD->SetBranchAddress("fTreeVariableNSigmasPosPion",&lNSigmasPosPion);
            lTreeMCFD->SetBranchAddress("fTreeVariableNSigmasNegPion",&lNSigmasNegPion);
            lTreeMCFD->SetBranchAddress("fTreeVariableRun", &fRun);
            lTreeMCFD->SetBranchAddress("fTreeVariableEvSel_AllSelections", &fEvSel_AllSelections);

            //================================================================

            //================================================================
            // Define Feeddown matrix
            //       FeedDownMatrix [Lambda Bin][Xi Bin];
            for(Int_t ilb = 0; ilb<100; ilb++) {
                for(Int_t ixb = 0; ixb<100; ixb++) {
                    lFeedDownMatrix[ilb][ixb]=0;
                }
            }
            //Define Xi Binning
            //Result should be ~invariant with respect to this choice...
            Double_t xibinlimits[] = {
                0.00,  0.20,  0.40,  0.60,  0.80,  0.90,
                1.00,  1.10,  1.20,  1.30,  1.40,  1.50,
                1.70,  1.90,  2.20,  2.60,  3.10,  3.90,
                4.90,  6.00,  7.20,  8.50 ,10.00, 12.00
            };
            Long_t xibinnumb = sizeof(xibinlimits)/sizeof(Double_t) - 1;

            TH1F* fHistPtXiReference =
                new TH1F("fHistPtXiReference","#Xi candidates uncorrected p_{T};p_{T} (GeV/c);Counts",xibinnumb,xibinlimits);

            //Feeddown: will use a double-coordinate system:
            // ( lambda bin, xi bin ) all the time!
            lWeAreAtBin = 0;         // Lambda bin
            Int_t lWeAreAtXiBin = 0; // Xi Bin

            cout<<"--------------- MC Data File Loop, Feeddown  -----------"<<endl;
            Long_t lOneTenthOfNCandidatesMCFD = ((double)(lNCandidatesMCFD) / 10. );
            for(Long_t icand = 0; icand<lNCandidatesMCFD; icand++) {
                // Int_t lBinProfile = cutProfCosPA->FindBin(lPt);
                lTreeMCFD->GetEntry(icand);

                if( icand % lOneTenthOfNCandidatesMCFD == 0 )
                    cout<<" Currently at candidate........: "<<icand<<" / "<<lNCandidatesMCFD<<" ( "<<(long)(((double)(icand)/(double)(lNCandidatesMCFD))*(100.+1e-3))<<"% )"<<endl;

                if (fEvSel_AllSelections == 0) continue;

                if(TMath::Abs(lRap)<fRapidityBoundary &&
                    // ((!fPosRap && !fNegRap) || (fPosRap && lRap > 0. && lNegEta > 0. && lPosEta >0) || (fNegRap && lRap<0. && lNegEta < 0. && lPosEta <0)) &&
                        TMath::Abs(lNegEta)       <= fCutDaughterEta               &&
                        TMath::Abs(lPosEta)       <= fCutDaughterEta               &&
                        lV0Radius                 >= fCutV0Radius                  &&
                        lDcaNegToPrimVertex       >= fCutDCANegToPV                &&
                        lDcaPosToPrimVertex       >= fCutDCAPosToPV                &&
                        lDcaV0Daughters           <= fCutDCAV0Daughters            &&
                        lV0CosinePointingAngle    >= fCutV0CosPA                   &&
                        lParticleMass*lDistOverTotMom    <= fCutProperLifetime     &&
                        lLeastNbrCrossedRows             >= fCutLeastNumberOfCrossedRows             &&
                        lLeastNbrCrossedRowsOverFindable >= fCutLeastNumberOfCrossedRowsOverFindable   &&
                        TMath::Abs(lInvariantMassCompetingOne - lCompetingParticleMass) > fCutCompetingV0Rejection &&
                        TMath::Abs(lInvariantMassCompetingTwo - lCompetingParticleMass) > fCutCompetingV0Rejection &&
                        ( fSpecialArmenterosCutK0s == kFALSE || ( fSpecialArmenterosCutK0s == kTRUE && lArmPt*5>TMath::Abs(lArmAlpha) ) )  &&
                        ( //perfect PID association, IsPhysicalPrimary association
                            ( fWhichParticle == "Lambda"
                            && lPID           == 3122 //V0 is a Lambda
                            && lPIDPositive   == 2212 //Pos Daughter is p
                            && lPIDNegative   == -211 //Neg Daughter is pi-
                            && (lPIDMother     == 3312 || (lPIDMother == 3322 && fFDSwitch=="UseMCRatio") )
                            && lPrimaryStatusMother == 1 //Xi is actually a primary (should matter little)
                            ) ||
                            ( fWhichParticle == "AntiLambda"
                            && lPID           == -3122  //V0 is an AntiLambda
                            && lPIDPositive   == 211    //Pos Daughter is pi+
                            && lPIDNegative   == -2212  //Neg Daughter is antiproton
                            && (lPIDMother     ==-3312 || (lPIDMother ==-3322 && fFDSwitch=="UseMCRatio") )
                            && lPrimaryStatusMother == 1 //Xi is actually a primary (should matter little)
                            )
                        )  &&
                     (
                        (fWhichParticle == "Lambda" && fTreeVariablePtLegPos > fPtLegProtCut ) ||
                        (fWhichParticle == "AntiLambda" && fTreeVariablePtLegNeg > fPtLegProtCut) ||
                        (fWhichParticle == "K0Short")) &&

                        //OOB condition
                        CheckITSTOF(lPosTrackStatus,  lNegTrackStatus,  fTreeVariablePosTOFBCid,  fTreeVariableNegTOFBCid)
                      //  CheckITSTOFOne( kFALSE, lPt, 2., lPosTrackStatus,  lNegTrackStatus,  fTreeVariablePosTOFBCid,  fTreeVariableNegTOFBCid)

                ) { // Start Entry Loop
                    lWeAreAtBin   = fHistPt->FindBin(lPtMC)-1;
                    if(lWeAreAtBin == -1) lWeAreAtBin   = 99; //UnderFlow, special treatment
                    lWeAreAtXiBin = fHistPtXiReference->FindBin(lPtMother)-1;
                    if(lWeAreAtXiBin == -1) lWeAreAtXiBin = 99; //UnderFlow, special treatment
                    //Populate Feeddown Matrix (only for peak region counts)
                    // cout<<" Populate at coordinates "<<lWeAreAtBin<<" vs "<<lWeAreAtXiBin<<" ....."<<endl;
                    if(lInvariantMass>lPeakLeftLimit[lWeAreAtBin]      && lInvariantMass<lPeakRightLimit[lWeAreAtBin]     )
                        { lFeedDownMatrix[lWeAreAtBin][lWeAreAtXiBin]++; /*printf("qui arrivo!!!!!! \n");*/}
                } // End Entry Loop
            }
            cout<<"--------------- Loop Completed -------------------------"<<endl;
            cout<<endl;

            //Filling histogram with original MC particles====================
            TH1F* fHistDummyV0FeedDown = 0x0;
            TH3F* f3DHistV0FeedDown = 0x0;


            if (fWhichParticle == "Lambda" ){
                f3DHistV0FeedDown = (TH3F*)v0listMCFD->FindObject("fHistGeneratedPtVsYVsCentralityXiMinus");
                fHistDummyV0FeedDown = (TH1F*)f3DHistV0FeedDown->ProjectionX("fHistPt_GenXiMinus",
                        f3DHistV0FeedDown->GetYaxis()->FindBin(-0.5+1e-2),
                        f3DHistV0FeedDown->GetYaxis()->FindBin(+0.5-1e-2),
                        f3DHistV0FeedDown->GetZaxis()->FindBin(0.),
                        f3DHistV0FeedDown->GetZaxis()->FindBin(+100-1e-2)
                        );
                //fHistDummyV0FeedDown = (TH1F*)v0listMCFD->FindObject("fHistPt_GenXiMinus");
            }
            if (fWhichParticle == "AntiLambda" ){
                f3DHistV0FeedDown = (TH3F*)v0listMCFD->FindObject("fHistGeneratedPtVsYVsCentralityXiPlus");
                fHistDummyV0FeedDown = (TH1F*)f3DHistV0FeedDown->ProjectionX("fHistPt_GenXiPlus",
                        f3DHistV0FeedDown->GetYaxis()->FindBin(-0.5+1e-2),
                        f3DHistV0FeedDown->GetYaxis()->FindBin(+0.5-1e-2),
                        f3DHistV0FeedDown->GetZaxis()->FindBin(0.),
                        f3DHistV0FeedDown->GetZaxis()->FindBin(+100-1e-2)
                        );
                //fHistDummyV0FeedDown = (TH1F*)v0listMCFD->FindObject("fHistPt_GenXiPlus");
            }

            TH1F *fHistMCCountbyptXiFeedDown	= new TH1F("fHistMCCountbyptXiFeedDown","#Xi MC count;p_{T} (GeV/c);Counts",xibinnumb,xibinlimits);
            Double_t temppt;//, y;// --- declared before
            for(long i = 1; i<fHistDummyV0FeedDown->GetNbinsX()+1; i++) {
                temppt = fHistDummyV0FeedDown->GetXaxis()->GetBinCenter(i);
                for(long filling = 0; filling<fHistDummyV0FeedDown->GetBinContent(i); filling++) {
                    fHistMCCountbyptXiFeedDown->Fill(temppt);
                }
            }
            if(fWhichParticle == "Lambda")
                cout<<"--------------- Generated Xi- Dump ---------------------"<<endl;
            if(fWhichParticle == "AntiLambda")
                cout<<"--------------- Generated Xi+ Dump ---------------------"<<endl;
            for(Int_t ipoint = 0; ipoint<xibinnumb; ipoint++) {
                cout<<"---> ["<<xibinlimits[ipoint]<<" - "<<xibinlimits[ipoint+1]<<" GeV/c]\tSignal: "<<fHistMCCountbyptXiFeedDown->GetBinContent(ipoint+1)<<" +/- "<<TMath::Sqrt(fHistMCCountbyptXiFeedDown->GetBinContent(ipoint+1))<<endl;
            }
            cout<<"--------------------------------------------------------"<<endl;
            cout<<endl;
            cout<<"Computing actual feeddown matrix"<<endl;

            for(Int_t ilb = 0; ilb<fptbinnumb; ilb++) {
                for(Int_t ixb = 0; ixb<xibinnumb; ixb++) {
                    if( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1)!= 0 ) { //avoid div by zero
                        lFeedDownEfficiency[ilb][ixb][ifile]=((double)(lFeedDownMatrix[ilb][ixb])) / ((double)( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1) )) ;
                        lFeedDownEfficiencyError[ilb][ixb][ifile]=ErrorInRatio(
                                                            ((double)(lFeedDownMatrix[ilb][ixb])),
                                                            TMath::Sqrt((double)(lFeedDownMatrix[ilb][ixb])),
                                                            ((double)( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1) )),
                                                            TMath::Sqrt((double)( fHistMCCountbyptXiFeedDown->GetBinContent(ixb+1) ))
                                                        );
                    } else {
                        lFeedDownEfficiency[ilb][ixb][ifile] = 0;
                        lFeedDownEfficiencyError[ilb][ixb][ifile] = 0;
                    }
                }
            }

            //Beware: the feeddown correction efficiency is computed with MC,
            //Which has the geant3/fluka problem. Thus, we actually geant3/fluka
            //correct the Feeddown efficiencies, too...

            for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
                for(Int_t ixb = 0; ixb<xibinnumb; ixb++) {
                    if ((fWhichParticle == "AntiLambda") && fG3FlukaCorrOn && ifile == 0) {
                        lFeedDownEfficiencyGFCorrected[ipoint][ixb][ifile] = lFeedDownEfficiency[ipoint][ixb][ifile]*(fitGeant3FlukaCorr->Eval( lProtonMomentum[ipoint][ifile]->GetMean() ) ) ;
                        lFeedDownEfficiencyGFCorrectedError[ipoint][ixb][ifile] = lFeedDownEfficiencyError[ipoint][ixb][ifile]*( fitGeant3FlukaCorr->Eval( lProtonMomentum[ipoint][ifile]->GetMean() ) ) ;
                    }
                    else
                    {
                        lFeedDownEfficiencyGFCorrected[ipoint][ixb][ifile] = lFeedDownEfficiency[ipoint][ixb][ifile];
                        lFeedDownEfficiencyGFCorrectedError[ipoint][ixb][ifile] = lFeedDownEfficiencyError[ipoint][ixb][ifile];
                    }
                }
            }

            //Create FD TH2Fs for storing
            f2dFeedDownMatrix[ifile] = new TH2F(Form("f2dFeedDownMatrix_mc%i",ifile),"",fptbinnumb,fptbinlimits,xibinnumb,xibinlimits);
            f2dFeedDownEfficiency[ifile]  = new TH2F(Form("f2dFeedDownEfficiency_mc%i",ifile),"",fptbinnumb,fptbinlimits,xibinnumb,xibinlimits);
            f2dFeedDownEfficiencyGFCorrected[ifile] = new TH2F(Form("f2dFeedDownEfficiencyGFCorrected_mc%i",ifile),"",fptbinnumb,fptbinlimits,xibinnumb,xibinlimits);
            f2dFeedDownMatrix[ifile]->SetDirectory(0);
            f2dFeedDownEfficiency[ifile]->SetDirectory(0);
            f2dFeedDownEfficiencyGFCorrected[ifile]->SetDirectory(0);
            for(Int_t ilb = 0; ilb<fptbinnumb; ilb++) {
                for(Int_t ixb = 0; ixb<xibinnumb; ixb++) {
                    f2dFeedDownMatrix[ifile]->SetBinContent(ilb+1,ixb+1,lFeedDownMatrix[ilb][ixb]);
                    f2dFeedDownEfficiency[ifile]->SetBinContent(ilb+1,ixb+1,lFeedDownEfficiency[ilb][ixb][ifile]);
                    f2dFeedDownEfficiency[ifile]->SetBinError(ilb+1,ixb+1,lFeedDownEfficiencyError[ilb][ixb][ifile]);
                    f2dFeedDownEfficiencyGFCorrected[ifile]->SetBinContent(ilb+1,ixb+1,lFeedDownEfficiencyGFCorrected[ilb][ixb][ifile]);
                    f2dFeedDownEfficiencyGFCorrected[ifile]->SetBinError(ilb+1,ixb+1,lFeedDownEfficiencyGFCorrectedError[ilb][ixb][ifile]);
                }
            }


            cout<<"--------------- Memory Cleanup -------------------------"<<endl;
            //fLevyFitXiFeedDown->Delete();
            v0listMCFD->Delete();
            delete v0listMCFD;
            fileMCFD->Close("R");
            fileMCFD->Delete();
            delete fileMCFD;
            cout<<endl;
           // fHistDummyV0FeedDown->Delete();
           // fHistMCCountbyptXiFeedDown->Delete();

        }//end loop files


        for(Int_t ilb = 0; ilb<fptbinnumb; ilb++) {
            for(Int_t ixb = 0; ixb<xibinnumb; ixb++) {
                lFeedDownEfficiencyTOT[ilb][ixb] = DoWeightedMeanXBin(0, lFeedDownEfficiencyGFCorrected[ilb][ixb][0], lFeedDownEfficiencyGFCorrected[ilb][ixb][1], lFeedDownEfficiencyGFCorrected[ilb][ixb][2],
                                                                     lFeedDownEfficiencyGFCorrectedError[ilb][ixb][0], lFeedDownEfficiencyGFCorrectedError[ilb][ixb][1], lFeedDownEfficiencyGFCorrectedError[ilb][ixb][2],
                                                                     fEventsWeights[0], fEventsWeights[1], fEventsWeights[2]);
                lFeedDownEfficiencyTOTError[ilb][ixb] = DoWeightedMeanXBin(1, lFeedDownEfficiencyGFCorrected[ilb][ixb][0], lFeedDownEfficiencyGFCorrected[ilb][ixb][1], lFeedDownEfficiencyGFCorrected[ilb][ixb][2],
                                                                          lFeedDownEfficiencyGFCorrectedError[ilb][ixb][0], lFeedDownEfficiencyGFCorrectedError[ilb][ixb][1], lFeedDownEfficiencyGFCorrectedError[ilb][ixb][2],
                                                                          fEventsWeights[0], fEventsWeights[1], fEventsWeights[2]);
            }
        }

        fHistFeeddownSubtraction = new TH1F("fHistFeeddownSubtraction","FD Subtraction;p_{T} (GeV/c);Ratio removed",fptbinnumb,fptbinlimits);
        fHistFeeddownSubtraction->SetDirectory(0);
        f2dFeedDownEfficiencyTOT  = new TH2F("f2dFeedDownEfficiencyTOT","",fptbinnumb,fptbinlimits,xibinnumb,xibinlimits);
        f2dFeedDownEfficiencyGFCorrectedTOT = new TH2F("f2dFeedDownEfficiencyGFCorrectedTOT","",fptbinnumb,fptbinlimits,xibinnumb,xibinlimits);

        for(Int_t ilb = 0; ilb<fptbinnumb; ilb++) {
            for(Int_t ixb = 0; ixb<xibinnumb; ixb++) {
                f2dFeedDownEfficiencyTOT->SetBinContent(ilb+1,ixb+1,lFeedDownEfficiencyTOT[ilb][ixb]);
                f2dFeedDownEfficiencyTOT->SetBinError(ilb+1,ixb+1,lFeedDownEfficiencyTOTError[ilb][ixb]);
                f2dFeedDownEfficiencyGFCorrectedTOT->SetBinContent(ilb+1,ixb+1,lFeedDownEfficiencyTOT[ilb][ixb]);
                f2dFeedDownEfficiencyGFCorrectedTOT->SetBinError(ilb+1,ixb+1,lFeedDownEfficiencyTOTError[ilb][ixb]);
            }
        }

        /// close Xi spectra file
       /* if(kFALSE){
        //
    Double_t lFeedDownToSubtractSyst[100];
    Double_t lProducedXiSyst[100];

        gMinuit->mnemat(&fMatrix[0][0],3);

        //Int_t attemptN = 0;
    for(int i=0; i<2; i++){

        cout<<"--------------------------------------------------------"<<endl;
            cout<<" Feeddown Parameters used for Tsallis (systematics) "<<endl;
            cout<<"--------------------------------------------------------"<<endl;

                for(Int_t ixb = 0; ixb<xibinnumb; ixb++) {
                if(fFDSwitch == "DoubleChargedXi") lProducedXiSyst[ixb] = lNInelasticEvents * fLevyFitXiFeedDownErr[i]->GetBinContent(ixb+1);
                else if(fFDSwitch == "UseMCRatio") lProducedXiSyst[ixb] = lNInelasticEvents * (fLevyFitXiFeedDownErr[i]->GetBinContent(ixb+1)/2.);
                cout<<"Yield Xi bin["<<ixb<<"] = "<< lProducedXiSyst[ixb] <<endl;
                }
            for(Int_t ilb = 0; ilb < fptbinnumb; ilb++) {
                lFeedDownToSubtractSyst[ilb] = 0;
                for(Int_t ixb = 0; ixb < xibinnumb; ixb++) {
                lFeedDownToSubtractSyst[ilb] += lProducedXiSyst[ixb]*lFeedDownEfficiencyGFCorrected[ilb][ixb];
                }
                if(lSigRealV0[ilb]>0){
                cout<<"FeeDown Lambda bin["<<ilb<<"] = "<< ((double)(lFeedDownToSubtractSyst[ilb])) / ((double)(lSigRealV0[ilb])) <<endl;
                cout<<"signal " << (double)(lSigRealV0[ilb]);
                feedDownSyst[i]->SetBinContent(ilb+1, ((double)(lFeedDownToSubtractSyst[ilb])) / ((double)(lSigRealV0[ilb])) );
                feedDownSyst[i]->SetBinError(ilb+1,0.);
                feedDownSystSub[i]->SetBinContent(ilb+1, 1.- (((double)(lFeedDownToSubtractSyst[ilb])) / ((double)(lSigRealV0[ilb]))) );
                feedDownSystSub[i]->SetBinError(ilb+1,0.);
                }else{
                cout<<"signal is zero!!!!" << endl;
                }
                }
                feedDownSyst[i]->SetTitle(Form("fHistFeeddownSubtraction_Syst_%d",i));
                feedDownSystSub[i]->SetTitle(Form("fHistFeeddownSubtraction_SystSub_%d",i));
    }

        } */
        fileXiSpectra->Close();


        //=== Computing actual feeddown numbers... ===
        Double_t lFeedDownToSubtract[100];
        Double_t lFeedDownToSubtractError[100];
        // === for systematics due to Xi spectra
        for(Int_t ilb = 0; ilb < fptbinnumb; ilb++) {
            lFeedDownToSubtract[ilb] = 0;
            lFeedDownToSubtractError[ilb] = 0;
            for(Int_t ixb = 0; ixb < xibinnumb; ixb++) {
                lFeedDownToSubtract[ilb] += lProducedXi[ixb]*lFeedDownEfficiencyTOT[ilb][ixb];
                lFeedDownToSubtractError[ilb] += TMath::Power( (lProducedXi[ixb]*lFeedDownEfficiencyTOTError[ilb][ixb]), 2 );
            }
            lFeedDownToSubtractError[ilb] = TMath::Sqrt(lFeedDownToSubtractError[ilb]);
        }

        for(Int_t ilb = 0; ilb<fptbinnumb; ilb++) {
            if(lSigRealV0[ilb] > 0.){
            fHistFeeddownSubtraction->SetBinContent( ilb+1 , ((double)(lFeedDownToSubtract[ilb])) / ((double)(lSigRealV0[ilb]) ) );
            fHistFeeddownSubtraction->SetBinError  ( ilb+1 ,
                    ErrorInRatio(
                        ((double)(lFeedDownToSubtract[ilb])),
                        ((double)(lFeedDownToSubtractError[ilb])),
                        ((double)(lSigRealV0[ilb]) ),
                        ((double)(lSigErrRealV0[ilb]) ) )
                                                );
                }else{
                fHistFeeddownSubtraction->SetBinContent( ilb+1 , 0. );
                fHistFeeddownSubtraction->SetBinError( ilb+1 , 0. );
            }

        }


        cout<<"--------------- FD Subtraction Fraction ----------------"<<endl;
        for(Int_t ilb = 0; ilb<fptbinnumb; ilb++) {
            cout<<"---> ["<<fptbinlimits[ilb]<<" - "<<fptbinlimits[ilb+1]<<" GeV/c]\t"<<fHistFeeddownSubtraction->GetBinContent(ilb+1)<<"\t+/-\t"<<fHistFeeddownSubtraction->GetBinError(ilb+1)<<endl;
        }
        cout<<"--------------------------------------------------------"<<endl;
        cout<<endl;
        cout<<"--------------------------------------------------------"<<endl;
        cout<<" Performing Actual Correction... "<<endl;
        for(Int_t ilb = 0; ilb<fptbinnumb; ilb++) {
            lSigRealV0   [ilb] = lSigRealV0[ilb] - lFeedDownToSubtract[ilb];
            lSigErrRealV0[ilb] = TMath::Sqrt( TMath::Power(lSigErrRealV0[ilb],2) + TMath::Power(lFeedDownToSubtractError[ilb],2) );
        }
        cout<<"--------------------------------------------------------"<<endl;
        cout<<endl;

    }//--- End Feeddown Correction section
    //=========================================================================


    //========================================================================
    //---> At this stage, everthing's just ready for the actual spectrum
    //---> computation to occur!
    Double_t lSpectrum[100];
    Double_t lSpectrumError[100];

    for(Int_t ibin=0; ibin<fptbinnumb; ibin++) {
        if(lEfficiency[ibin] > 0){
        lSpectrum[ibin] = lSigRealV0[ibin] / lEfficiency [ibin];
        lSpectrumError[ibin] = ErrorInRatio(
                                   lSigRealV0    [ibin],
                                   lSigErrRealV0 [ibin],
                                   lEfficiency      [ibin],
                                   lEfficiencyError [ibin]);
	} else{
         lSpectrum[ibin] = 999.;
         lSpectrumError[ibin] = 0.;
        }
    }
    //Divide by: Bin Width, Rapidity Window, N_{Inel}
    for(Int_t ibin=0; ibin<fptbinnumb; ibin++) {
        lSpectrum[ibin] /= (fptbinlimits[ibin+1]-fptbinlimits[ibin]);
        lSpectrum[ibin] /= lNInelasticEvents;
        if(fPosRap || fNegRap) lSpectrum[ibin] /= fRapidityBoundary;
        else lSpectrum[ibin] /= 2*fRapidityBoundary;
        lSpectrumError[ibin] /= (fptbinlimits[ibin+1]-fptbinlimits[ibin]);
        lSpectrumError[ibin] /= lNInelasticEvents;
        if(fPosRap || fNegRap) lSpectrumError[ibin] /= fRapidityBoundary;
        else lSpectrumError[ibin] /= 2*fRapidityBoundary;
    }

    TH1F* fHistPtLambda = new TH1F("fHistPtLambda","#Lambda Corrected Spectrum;p_{T} (GeV/c);1/N_{INEL} #frac{d^{2}N}{dydp_{t}}",fptbinnumb,fptbinlimits);
    TH1F* fHistPtAntiLambda = new TH1F("fHistPtAntiLambda","#bar{#Lambda} Corrected Spectrum;p_{T} (GeV/c);1/N_{INEL} #frac{d^{2}N}{dydp_{t}}",fptbinnumb,fptbinlimits);
    TH1F* fHistPtK0Short = new TH1F("fHistPtK0Short","K^{0}_{S} Corrected Spectrum;p_{T} (GeV/c);1/N_{INEL} #frac{d^{2}N}{dydp_{t}}",fptbinnumb,fptbinlimits);

    //Copy to Histogram
    for(Int_t ibin=0; ibin<fptbinnumb; ibin++) {
        if(fWhichParticle == "Lambda" ) {
            fHistPtLambda->SetBinContent( ibin+1, lSpectrum[ibin]      );
            fHistPtLambda->SetBinError  ( ibin+1, lSpectrumError[ibin] );
        }
        if(fWhichParticle == "AntiLambda" ) {
            fHistPtAntiLambda->SetBinContent( ibin+1, lSpectrum[ibin]      );
            fHistPtAntiLambda->SetBinError  ( ibin+1, lSpectrumError[ibin] );
        }
        if(fWhichParticle == "K0Short" ) {
            fHistPtK0Short->SetBinContent( ibin+1, lSpectrum[ibin]      );
            fHistPtK0Short->SetBinError  ( ibin+1, lSpectrumError[ibin] );
        }
    }


    //=========================================================================

    cout<<"--------------- Result Output --------------------------"<<endl;
    cout<<" ---> Writing information to "<<fOutputDataFile<<endl;
    // Open an output file
    TFile* lResultsFile = TFile::Open(fOutputDataFile, "RECREATE");
    if (!lResultsFile || !lResultsFile->IsOpen()) {
        cout<<"Error! Couldn't open file!"<<endl;
        return;
    }

    //Preparing Signal Extraction Range Canvas
    TCanvas *cSigExtRange = new TCanvas("cSigExtRange","Extraction Range",900,600);
    cSigExtRange->SetFillColor(kWhite);
    cSigExtRange->SetLeftMargin(0.17);
    cSigExtRange->SetRightMargin(0.17);
    cSigExtRange->cd();
    fHistSignalExtractionRange->SetFillColor(18);
    fHistSignalExtractionRange->SetMarkerStyle(20);
    if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ) {
        fHistSignalExtractionRange->GetYaxis()->SetRangeUser(1.115-0.08,1.115+0.08);
    }
    if( fWhichParticle == "K0Short" ) {
        fHistSignalExtractionRange->GetYaxis()->SetRangeUser(0.498-0.15,0.498+0.15);
    }
    fHistSignalExtractionRange->Draw("E2");
    if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ) {
        fLDataUpper->Draw("same");
        fLDataLower->Draw("same");
    }
    if( fWhichParticle == "K0Short" ) {
        fKDataUpper->Draw("same");
        fKDataLower->Draw("same");
    }
    //Saving Invariant Mass Plots (real data)
    lResultsFile->cd();
    TDirectoryFile *lInvMassReal = new TDirectoryFile("lInvMassReal","Invariant Mass Plots (Real Data)");
    lInvMassReal->cd();
    cSigExtRange->Write();
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        lCanvasHistoFullV0[ibin] -> Write();
    }

    TDirectoryFile *lInvMassRealRawData = new TDirectoryFile("lInvMassRealRawData","Objects for Inv Mass Plots (Real Data)");
    lInvMassRealRawData->cd();
    fHistSignalExtractionRange->Write();
    fHistPeakPosition->Write();
    fHistSigToNoise->Write();
    fHistPtRaw->Write();
    if(fFillQA){
    lHistoV0Radius->Write();
    lHistDCAnegToPrim->Write();
    lHistDCAposToPrim->Write();
    lHistDCAV0Daughters->Write();
    lProperLifeTime->Write();
    lCosPointingAngle->Write();
    lHistoNsigmaTPCpos->Write();
    lHistoNsigmaTPCneg->Write();
    //
    lHistoV0RadiusMC->Write();
    lHistDCAnegToPrimMC->Write();
    lHistDCAposToPrimMC->Write();
    lHistDCAV0DaughtersMC->Write();
    lProperLifeTimeMC->Write();
    lCosPointingAngleMC->Write();
    lHistoNsigmaTPCposMC->Write();
    lHistoNsigmaTPCnegMC->Write();
    }
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        lHistoFullV0[ibin] -> Write();
        lHistoSelectedV0[ibin] -> Write();
        fgausPt[ibin]         -> Write();
        if( fFitBackgroundSwitch ) lfitNoise[ibin] -> Write();
    }
    if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ) {
        fLDataUpper->Write();
        fLDataLower->Write();
    }
    if( fWhichParticle == "K0Short" ) {
        fKDataUpper->Write();
        fKDataLower->Write();
    }


    //Saving Invariant Mass Plots (MC)
    lResultsFile->cd();
    TDirectoryFile *lInvMassMC = new TDirectoryFile("lInvMassMC","Invariant Mass Plots (Monte Carlo)");
    lInvMassMC->cd();

    if(kDoEfficiency){
        fHistEfficiencyParz[0]->Write();
        fHistEfficiencyParz[1]->Write();
        fHistEfficiencyParz[2]->Write();
        fHistPureEfficiencyParz[0]->Write();
        fHistPureEfficiencyParz[1]->Write();
        fHistPureEfficiencyParz[2]->Write();
        fHistReco[0]->Write();
        fHistReco[1]->Write();
        fHistReco[2]->Write();
    }
/*
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        lCanvasHistoFullV0MC[ibin] -> Write();
    }

    TDirectoryFile *lInvMassMCRawData = new TDirectoryFile("lInvMassMCRawData","Objects for Inv Mass Plots (MC)");
    lInvMassMCRawData->cd();
    fHistPeakPositionMC->Write();
    fHistSigToNoiseMC->Write();
    f2dHistPtResolution->Write();
    f2dHistPtBlur->Write();
    f2dHistPtSharpen->Write();
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        lHistoFullV0MC[ibin] ->Write();
    //    lProtonMomentumData[ibin]  ->Write();
        lProtonMomentumMC[ibin]  ->Write();
        fgausPtMC[ibin]         ->Write();
        if( fFitBackgroundSwitch ) lfitNoiseMC[ibin] -> Write();
    }


    //Saving Geant3/Fluka Correction Data (MC)
    if( (fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda")  && fG3FlukaCorrOn) {
        lResultsFile->cd();
        TDirectoryFile *lGFCorrection = new TDirectoryFile("lGFCorrection","Geant3/Fluka Correction Histograms");
        lGFCorrection->cd();
        // h1dPanosCorrections->Write();
        fitGeant3FlukaCorr->Write();
        for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
            lProtonMomentum[ibin]  ->Write();
        }
    }

    fHistEffNumerator->Write();
    fHistEffDenominator->Write();*/


    //Saving Feeddown Correction information, if needed
    if( (fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda") && fFDSwitch!="NoFD" ) {
        lResultsFile->cd();
        TDirectoryFile *lFeeddown = new TDirectoryFile("lFeeddown","Feeddown Subtraction Information");
        lFeeddown->cd();
        for(int ifile = 0; ifile < 3; ifile++){
            f2dFeedDownMatrix[ifile]->Write();
            f2dFeedDownEfficiency[ifile]->Write();
            f2dFeedDownEfficiencyGFCorrected[ifile]->Write();
        }
        f2dFeedDownEfficiencyTOT->Write();
        f2dFeedDownEfficiencyGFCorrectedTOT->Write();
        fHistFeeddownSubtraction->Write();
         for(Int_t ibin = 0; ibin<2; ibin++) {
         /*feedDownSyst[ibin]->Write();
         feedDownSystSub[ibin]->Write();*/
         }
    }

    //Saving Average Multiplicity
    fHistAverageMult->Write();

    //Saving Resolution Information

    //preparing...
    for(Int_t ibin=0; ibin<fptbinnumb; ibin++) {
        fHistResolutionVsPt->SetBinContent(ibin+1,lHistResolution[ibin]->GetMean());
        fHistResolutionVsPt->SetBinError(ibin+1,lHistResolution[ibin]->GetRMS());
        fHistResolutionVsPtDivByBinWidth->SetBinContent(ibin+1,lHistResolution[ibin]->GetMean() / (fptbinlimits[ibin+1]-fptbinlimits[ibin]) );
        fHistResolutionVsPtDivByBinWidth->SetBinError(ibin+1,lHistResolution[ibin]->GetRMS() / (fptbinlimits[ibin+1]-fptbinlimits[ibin]) );

        fHistResolutionVsPtWithGaussians->SetBinContent(ibin+1,lHistResolutionGaussian[ibin]->GetParameter(1));
        fHistResolutionVsPtWithGaussians->SetBinError(ibin+1,lHistResolutionGaussian[ibin]->GetParameter(2));
        fHistResolutionVsPtDivByBinWidthWithGaussians->SetBinContent(ibin+1,lHistResolutionGaussian[ibin]->GetParameter(1) / (fptbinlimits[ibin+1]-fptbinlimits[ibin]) );
        fHistResolutionVsPtDivByBinWidthWithGaussians->SetBinError(ibin+1,lHistResolutionGaussian[ibin]->GetParameter(2) / (fptbinlimits[ibin+1]-fptbinlimits[ibin]) );
    }


    lResultsFile->cd();
    TDirectoryFile *lDirResolution = new TDirectoryFile("lInfoResolution","Resolution Information");
    lDirResolution->cd();
    fHistResolutionVsPt->Write();
    fHistResolutionVsPtDivByBinWidth->Write();
    fHistResolutionVsPtWithGaussians->Write();
    fHistResolutionVsPtDivByBinWidthWithGaussians->Write();
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        lHistResolution[ibin] ->Write();
        lHistPtDistr[ibin] ->Write();
        lHistPtDistrMCtruth[ibin] ->Write();
        lHistResolutionGaussian[ibin]->Write();
     }


    lResultsFile->cd();
    fHistMCtruth->Write();
    // HistEvents->Write();

    if( fWhichParticle == "K0Short") fHistPureEfficiency->Write();
    if( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ) {
        fHistPureEfficiency->Write();
        fHistEfficiency->Write();
    }
    if(fWhichParticle == "Lambda" ) {
        fHistPtLambda     ->  Write();
    }
    if(fWhichParticle == "AntiLambda" ) {
        fHistPtAntiLambda -> Write();
    }
    if(fWhichParticle == "K0Short" ) {
        fHistPtK0Short    -> Write();
    }
    // fHistMCCountbyptV0->Write();
    lResultsFile->Write();
    lResultsFile->Close();
    delete lResultsFile;

    //================================================
    //Manual Cleanup of all created Histograms
    //================================================
    // needed if you want to re-run the whole thing without
    // memory leaks (systematics, etc)

    //switch on if you want large amounts of printout
    Bool_t lDebugCleaningProcess = kTRUE;

    if (lDebugCleaningProcess) cout<<"fHistPt->Delete()"<<endl;
    fHistPt->Delete();
    fHistPtRaw->Delete();
    if (lDebugCleaningProcess) cout<<"fHistEffNumerator->Delete()"<<endl;
    fHistEffNumerator->Delete();
    if (lDebugCleaningProcess) cout<<"fHistEffDenominator->Delete()"<<endl;
    fHistEffDenominator->Delete();
    if (lDebugCleaningProcess) cout<<"fHistPeakPosition->Delete()"<<endl;
    fHistPeakPosition->Delete();
    if (lDebugCleaningProcess) cout<<"fHistPeakPositionMC->Delete()"<<endl;
    fHistPeakPositionMC->Delete();
    if (lDebugCleaningProcess) cout<<"fHistSigToNoise->Delete()"<<endl;
    fHistSigToNoise->Delete();
    if (lDebugCleaningProcess) cout<<"fHistSigToNoiseMC->Delete()"<<endl;
    fHistSigToNoiseMC->Delete();
    if (lDebugCleaningProcess) cout<<"fHistSignalExtractionRange->Delete()"<<endl;
    fHistSignalExtractionRange->Delete();
    if (lDebugCleaningProcess) cout<<"fKData*->Delete()"<<endl;
    fKDataUpper->Delete();
    fKDataLower->Delete();
    fLDataUpper->Delete();
    fLDataLower->Delete();
    if(lDebugCleaningProcess) cout<<"f2dHistPtResolution->Delete()"<<endl;
    f2dHistPtResolution->Delete();

    if(lDebugCleaningProcess) cout<<"lfitNoise*[*]->Delete(); lSampleNoise*[*]->Delete()"<<endl;
    if( fFitBackgroundSwitch ) {
        /*for(long i=0; i<fptbinnumb; i++) {
            lfitNoise[i]    -> Delete();
            lSampleNoise[i] -> Delete();
            lfitNoiseMC[i]    -> Delete();
            lSampleNoiseMC[i] -> Delete();
        }*/
    }

    //average multiplicity histo
    fHistAverageMult->Delete();


    //pt-by-pt histos
    fHistResolutionVsPt->Delete();
    fHistResolutionVsPtDivByBinWidth->Delete();
    fHistResolutionVsPtWithGaussians->Delete();
    fHistResolutionVsPtDivByBinWidthWithGaussians->Delete();
    f2dHistPtBlur->Delete();
    f2dHistPtSharpen->Delete();
    if (lDebugCleaningProcess) cout<<"lHistoFullV0*[*]->Delete()"<<endl;
    for(Int_t ihist=0; ihist<100; ihist++) {
        lHistoFullV0[ihist]->Delete();
        lHistoSelectedV0[ihist]->Delete();
        lHistoFullV0MC[ihist]->Delete();
    //    lProtonMomentumData[ihist]  ->Delete();
        lProtonMomentumMC[ihist]  ->Delete();
        lHistResolution[ihist]->Delete();
        lHistPtDistr[ihist]->Delete();
        lHistPtDistrMCtruth[ihist]->Delete();
        lHistResolutionGaussian[ihist]->Delete();
    }

    if (lDebugCleaningProcess) cout<<"lLine*[*]->Delete()"<<endl;
    //gaussian fit functions, drawing lines
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        /*lLineLeftMost[ibin] ->Delete();
        lLineLeft[ibin]     ->Delete();
        lLineRight[ibin]    ->Delete();
        lLineRightMost[ibin]->Delete();
        lLineLeftMostMC[ibin] ->Delete();
        lLineLeftMC[ibin]     ->Delete();
        lLineRightMC[ibin]    ->Delete();
        lLineRightMostMC[ibin]->Delete();

        if ( fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda" ) {
            fgausPt[ibin]->Delete();
            fgausPtMC[ibin]->Delete();
        }
        if ( fWhichParticle == "K0Short") {
            fgausPt[ibin]->Delete();
            fgausPtMC[ibin]->Delete();
        }*/
    }
    if (lDebugCleaningProcess) cout<<"lCanvasHistoFullV0*[*]->Delete()"<<endl;
    /*for(Int_t ihist=0; ihist<100; ihist++) {
        lCanvasHistoFullV0[ihist]->Close();
        lCanvasHistoFullV0MC[ihist]->Close();
        delete lCanvasHistoFullV0[ihist];
        delete lCanvasHistoFullV0MC[ihist];
    }*/

    if (lDebugCleaningProcess) cout<<"cSigExtRange->Delete()"<<endl;
    cSigExtRange->Close();
    delete cSigExtRange;

    if (lDebugCleaningProcess) cout<<"lProtonMomentum[*]->Delete()"<<endl;
    //histograms for G3/F correction
    for(Int_t ibin=0; ibin<100; ibin++) {
        for(Int_t ifile=0; ifile<3; ifile++) {
            lProtonMomentum[ibin][ifile]->Delete();
        }
    }

    if (lDebugCleaningProcess) cout<<"fHist*Efficiency->Delete()"<<endl;
    //MC info, efficiencies
    fHistPureEfficiency->Delete();
    fHistEfficiency->Delete();
    //data for G3/F, continued

    if( (fWhichParticle == "Lambda" || fWhichParticle == "AntiLambda") && fG3FlukaCorrOn ) {
        fitGeant3FlukaCorr->Delete();
    }

    //histograms, feeddown
    if (lDebugCleaningProcess) cout<<"f2dFeedDown*->Delete()"<<endl;
    if( fFDSwitch != "NoFD" && fWhichParticle != "K0Short") {
        for (int i = 0; i<3; i++){
            f2dFeedDownMatrix[i]->Delete();
            f2dFeedDownEfficiency[i]->Delete();
            f2dFeedDownEfficiencyGFCorrected[i]->Delete();
        }
        fHistFeeddownSubtraction->Delete();
        f2dFeedDownEfficiencyTOT->Delete();
        f2dFeedDownEfficiencyGFCorrectedTOT->Delete();

        for(Int_t ibin = 0; ibin<2; ibin++) {
           // feedDownSyst[ibin]->Delete();
           // feedDownSystSub[ibin]->Delete();
        }
    }
    if (lDebugCleaningProcess) cout<<"fHistPt*->Delete()"<<endl;
    //Corrected Spectra Histograms
    fHistPtLambda->Delete();
    fHistPtAntiLambda->Delete();
    fHistPtK0Short->Delete();

    for (int i = 0; i<3; i++){
        fHistEfficiencyParz[i]->Delete();
        fHistPureEfficiencyParz[i]->Delete();
        fHistReco[i]->Delete();
    }

    //Exit Batch Mode
    gROOT->SetBatch (kFALSE);

    cout<<"--------------------------------------------------------"<<endl;
    cout<<" There, done! "<<endl;
    cout<<endl;




}

Double_t AliV0Module::FuncPlusErr(Double_t *xx, Double_t *p) {

  Double_t value = fLevyFitXiFeedDown->Eval(xx[0]);
  Double_t error = ErrorFunction(xx, p);
// printf("error %f \n",error);
  if(p[0] > 0) {
    return value + error;
  }
  return value -error;
}


Double_t AliV0Module::ErrorFunction(Double_t *xx, Double_t *p ) {

  //TF1 * fLevyFitXiFeedDown = (TF1*)TVirtualFitter::GetFCN();

  Double_t func[3];
  //  In general, one can compute the derivative numerically using root.
  for(Int_t ipar = 0; ipar < 3; ipar++){
         func[ipar] =  fLevyFitXiFeedDown->GradientPar(ipar, xx);
      // printf("der = %f \n", func[ipar]);
  }



  Double_t variance = 0;

  for(Int_t ipar = 0; ipar < 3; ipar++){
    for(Int_t jpar = 0; jpar < 3; jpar++){
      variance = variance + func[ipar]*func[jpar] * fMatrix[ipar][jpar];
    }

  }
  Double_t error = TMath::Sqrt(variance);
  //return error;
  Double_t ndf = fLevyFitXiFeedDown->GetNDF();
  Double_t chidf =  TMath::Sqrt(fLevyFitXiFeedDown->GetChisquare()/fLevyFitXiFeedDown->GetNDF());
  Double_t t = TMath::StudentQuantile(0.5 + 0.68/2, fLevyFitXiFeedDown->GetNDF());
  return error*t*chidf;

}


