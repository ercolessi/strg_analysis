#ifndef CorrectionClass_H
#define CorrectionClass_H

/***********************************************
  Cascade Correction Module
  ----------------------
    This version: 2021
    --- Francesca Ercolessi
    francesca.ercolessi@cern.ch
 ***********************************************/

class CorrectionClass{
public: 
  //Constructor
  CorrectionClass(); 
  CorrectionClass(TString ParticleType);

  //Set Files to Use
  void SetOutputDataFile  ( TString DataFilename );
  void SetMCFilepTshape  ( TString MCFilename );
  void SetMCFileNormCorr  ( TString MCNormCorr );

  //Do Analysis
  void DoAnalysis();

  //Set Pt Bin Limits
  void SetPtBinLimits(Long_t got_ptbinnumb, const Double_t *got_ptbinlimits);

  //Set Variable Bin Limits
  void SetPercBinLimits(Long_t got_percbinnumb, const Double_t *got_percbinlimits);

  //Set Fixed Bin Limits
  void SetFixedPercLimits(Double_t got_Low, Double_t got_High);

  //Set Estimators
  void SetVarEstimator(TString got_var);
  void SetFixedEstimator(TString got_fixed);
  void SetEnergyEstimator(TString got_energyvar);
  void SetMultEstimator(TString got_multvar);

  //Print Configuration
  void PrintConfiguration();

  //Do corrections
  void DoMCShape(Bool_t got_DoMCShape);
  void DoEventLoss(Bool_t got_DoEventLoss);
  void DoSignalLoss(Bool_t got_DoEventLoss);

  //Compute errors
  Double_t ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
  void ErrorInRatioCorr ( TH1F* h1, TH1F *h2 );

  //Levy-Tsallis function
  static Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p);
  TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t par2min = 1e-3, Double_t par2max = 1e+3,Double_t par3min = 1e-6, Double_t par3max = 1e+6, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);

private:
  //root file names
  TString fOutputDataFile;
  TString fMCFilepTshape;
  TString fMCFileNormCorr;
  //WhichParticle switch
  TString fWhichParticle;
  //
  Double_t fLowFixed;
  Double_t fHighFixed;
  //
  TString fWhichVarEstimator;
  TString fWhichFixedEstimator;
  TString fWhichMultEstimator;
  TString fWhichEnergyEstimator;

  //Setters for pT and percentiles
  Double_t fptbinlimits[200];
  Double_t fptX[200];
  Long_t fptbinnumb;
  Double_t fpercbinlimits[200];
  Double_t fpercX[200];
  Long_t fpercbinnumb;

  //Booleans for corrections
  Bool_t fDoMCShape;
  Bool_t fDoEventLoss;
  Bool_t fDoSignalLoss;

};

#endif