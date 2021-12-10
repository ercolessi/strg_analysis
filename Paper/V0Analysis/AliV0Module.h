#ifndef ALIV0MODULE_H
#define ALIV0MODULE_H
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

  Lambda Analysis Module - Header
  -------------------------------

This version: 27th April 2012 

--- David Dobrigkeit Chinellato
    daviddc@ifi.unicamp.br

***********************************************/

class AliV0Module{
public: 
  //Constructor
  AliV0Module(); 
  AliV0Module(TString ParticleType);

  //Set Files to Use
  void SetRealDataFile      ( TString RealDataFilename     );
  void SetMCDataFile        ( TString MCDataFilename       );
  void SetFeedDownDataFile  ( TString FeedDownDataFilename );
  void SetOutputFile        ( TString OutputFilename       );
  void SetListName          ( TString ListName            ) {fListName = ListName;};
  void SetTOFpercFileName   ( TString TOFpercFilename  );

  //Set Pt Bin Limits
  void SetPtBinLimits(Long_t got_ptbinnumb, const Double_t *got_ptbinlimits);

  //Set Rapidity Window
  void SetRapidityWindow(Double_t got_RapidityBoundary);

  //Set CINT1B/INEL to normalize to yield
  void SetCINT1BoverINEL(Double_t got_CINT1BoverINEL);

  //Set Cuts - topological
  void SetCutV0Radius       (Double_t cut);
  void SetCutDCANegToPV     (Double_t cut);
  void SetCutDCAPosToPV     (Double_t cut);
  void SetCutDCAV0Daughters (Double_t cut);
  void SetCutV0CosPA        (Double_t cut);

  //Set Cuts - other
  void SetCutProperLifetime                         (Double_t cut);
  void SetCutTPCPIDNSigmas                          (Double_t cut);
  void SetCutSigmaForSignalExtraction               (Double_t cut);
  void SetCutLeastNumberOfCrossedRows               (Double_t cut);
  void SetCutLeastNumberOfCrossedRowsOverFindable   (Double_t cut);
  void SetCutDaughterEta                            (Double_t cut);
  void SetCutCompetingV0Rejection                   (Double_t cut);
  void SetITSTOFrequest (Int_t itstofreq) {fRequestITSTOF = itstofreq;};
  void SetPtLimitForTOF (Double_t ptmin) {fPtMinITSTOF = ptmin;};
  //Set Feeddown treatment
  void SetFeeddownTreatment ( TString FDMethod );

  //Set Fit Background or not 
  void SetFitBackground ( Bool_t fitBgSwitch );

  //Switvh on/off G3/Fluka correction
  void SetGeant3FlukaCorr ( Bool_t switchOnG3Fk );
  void SetPosRap( Bool_t lPosRap );
  void SetNegRap( Bool_t lNegRap );
  Bool_t GetPosRap();
  Bool_t GetNegRap();

  //Multiplicity Study Setters
  void SetPerformMultiplicityStudy ( Bool_t lPerformMultStudy );
  void SetLowMultValue             ( Double_t lLoMultBound    );
  void SetHighMultValue            ( Double_t lHiMultBound    );
  void SetLowEEValue             ( Double_t lLoEEBound    );
  void SetHighEEValue            ( Double_t lHiEEBound    );
  
  //Set Estimators
  void SetWhichEstimators ( TString got_multestimator , TString got_effenergyestimator  );

  //TOFpercentile
  Double_t GetTOFpercentile ( TFile* lfilename , Int_t lNTOFtrgPads, Int_t lRun );

  //Set Fit Background or not 
  void SetSpecialArmenterosCutK0s ( Bool_t lSpecialArmenterosCutK0s );

  void SetUseIntegratedEff(Bool_t useIntegrEff) {fUseIntegratedEfficiencies = useIntegrEff;}
  //Do Analysis
  void DoAnalysis();

  //Set Default Cuts
  void SetDefaultCuts();
  void SetMinimumCuts();
  void SetPtProCut(Double_t cut){fPtLegProtCut = cut;}

  //Auxiliary Functions 
  TString IntToString(int input);
  TString DoubleToString(double input);
  Double_t ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
  Double_t MyGeant3FlukaCorrectionForProtons(const Double_t *x, const Double_t *par);
  Double_t MyGeant3FlukaCorrectionForAntiProtons(const Double_t *x, const Double_t *par);
  Double_t MyLevyPtXi(const Double_t *pt, const Double_t *par);
  Double_t MyBgPol1(const Double_t *x, const Double_t *par);
  Double_t MyBgPolToEval1(const Double_t *x, const Double_t *par);
  Double_t RoundToThousandth( const Double_t lToRound );
  Double_t FuncPlusErr(Double_t *xx, Double_t *p);
  Double_t ErrorFunction (Double_t *xx, Double_t *p ); 
  Bool_t CheckITSTOF ( ULong64_t lPosTrackStatus, ULong64_t lNegTrackStatus, Int_t lPosTOFBCID, Int_t lNegTOFBCID );
 

  //Set primary definitions (MC)
  void SetPrimarySelection ( TString lPrimSelMethod );
  void SetDistToPVPrimary ( Double_t fDistToPVSelection );

private:
  //root file names
  TString fRealDataFile;
  TString fMCDataFile;
  TString fFeedDownDataFile;
  TString fOutputDataFile;
  TString fTOFpercFilename;

  //Store Pt Bin Limits
  //Max Number of Pt Bins Set here (100)
  Double_t fptbinlimits[200];
  Double_t fptX[200];
  Long_t fptbinnumb;
  
  /// Xi spectra for feeddown
  TF1 *fLevyFitXiFeedDown; 
  TH1D *fLevyFitXiFeedDownErr[2];
  Double_t fMatrix[3][3];

  //Rapidity Range Window
  Double_t fRapidityBoundary;

  //CINT1B/INEL for normalization
  Double_t fCINT1BoverINELratio; 
  
  //Main Analysis Parameters
  //--- 5 Topological Selections
  Double_t fCutV0Radius;
  Double_t fCutDCANegToPV;
  Double_t fCutDCAPosToPV;
  Double_t fCutDCAV0Daughters;
  Double_t fCutV0CosPA;
  //--- Proper Lifetime
  Double_t fCutProperLifetime;
  //--- TPC dE/dx N-sigmas
  Double_t fCutTPCPIDNSigmas;
  //--- Sigmas for signal extraction
  Double_t fCutNSigmasForSignalExtraction;
  //--- Smallest Number of Crossed Rows in TPC accepted for tracks
  Double_t fCutLeastNumberOfCrossedRows;
  Double_t fCutLeastNumberOfCrossedRowsOverFindable;
  // pt prot
  Double_t fPtLegProtCut;
  //--- Daughter Track eta cut 
  Double_t fCutDaughterEta;
  //--- Competing V0 Species Rejection
  Double_t fCutCompetingV0Rejection;

  //WhichParticle switch: "Lambda", "AntiLambda" or "K0Short"
  TString fWhichParticle;

  //Mult Estimator 
  TString fWhichMultEstimator;
  TString fWhichEffEnergyEstimator;
  
  //Do Fitting instead of bin counting switch
  Bool_t fFitBackgroundSwitch;

  // switch on/off G3/Fluka correction
  Bool_t fG3FlukaCorrOn;

  // switch on QA
  Bool_t fFillQA;
  Bool_t fPosRap;
  Bool_t fNegRap;

  //Special Armenteros Selection For K0s
  Bool_t fSpecialArmenterosCutK0s;

  Bool_t fUseIntegratedEfficiencies;
  //Perform multiplicity selection 
  // --- (mult estimator in pp, centrality in PbPb)
  Bool_t fPerformMultiplicityStudy; 
  Double_t fLoMultBound;
  Double_t fHiMultBound;
  Double_t fLoEEBound;
  Double_t fHiEEBound;


  //FeedDownTreatment switch:  
  // --- NoFD.............: Doesn't feeddown subtract at all 
  // --- DoubleChargedXi..: Multiply Charged Xi subtraction by 2
  // --- UseMCRatio.......: Fill FD Matrix with Xi- and Xi0 
  //  (effectively use Xi0/Xi- from MC, should not be used with 
  //   Xi- and Xi+ enhanced datasets!!)
  TString fFDSwitch;
  TString fListName;

  //Primary Selection Switch: 
  // --- IsPhysicalPrimary....: Usual method
  // --- DistToPV.............: Iouri, Luke
  TString fPrimarySelection;
  Double_t fDistToPVCut;
  Int_t fRequestITSTOF;
  Double_t  fPtMinITSTOF;

};
#endif
