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
  void SetMCFilepTshape  ( TString MCFilename0, TString MCFilename1, TString MCFilename2 );
  void SetMCFileEventLoss  ( TString MCFilename );
  void SetMCFileSignalLoss  ( TString MCFilename );

  //Do Analysis
  void DoAnalysis();

  //Set Pt Bin Limits
  void SetPtBinLimits(Long_t got_ptbinnumb, const Double_t *got_ptbinlimits);

  //Set Variable Bin Limits
  void SetPercBinLimits(Long_t got_percbinnumb, std::vector<Double_t> got_percbinlimitsMultlow, std::vector<Double_t> got_percbinlimitsMulthigh, std::vector<Double_t> got_percbinlimitsEnergylow, std::vector<Double_t> got_percbinlimitsEnergyhigh);

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
  void DoSystematics(Bool_t got_DoSystematics);

  //Compute errors
  Double_t ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
  void ErrorInRatioCorr ( TH1F* h1, TH1F *h2 );
  void ErrorInRatioUncorr( TH1F* h1, TH1F *h2 ) ;
  Double_t PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin);

  //Levy-Tsallis function
  static Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p);
  TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t par2min = 1e-3, Double_t par2max = 1e+3,Double_t par3min = 1e-5, Double_t par3max = 1e+5, Double_t n = 5., Double_t C = 0.1, Double_t norm = 0.1);

  //Weighted mean
  Float_t DoWeightedMean(Bool_t err , Float_t h1, Float_t h2, Float_t h3, Float_t e1, Float_t e2, Float_t e3, Float_t w1, Float_t w2, Float_t w3);

private:
  //root file names
  TString fOutputDataFile;
  TString fMCFilepTshape0;
  TString fMCFilepTshape1;
  TString fMCFilepTshape2;
  TString fMCFileSignalLoss;
  TString fMCFileEventLoss;
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
  Double_t fpercentileMultlow[200];
  Double_t fpercentileMulthigh[200];
  Double_t fpercentileEnergylow[200];
  Double_t fpercentileEnergyhigh[200];
  Long_t fpercbinnumb;

  //Booleans for corrections
  Bool_t fDoMCShape;
  Bool_t fDoEventLoss;
  Bool_t fDoSignalLoss;
  Bool_t fDoSystematics;

};

#endif