#include "AliV0Module.h"

void runV0Analysis(
  Double_t lMultBoundLo = 0.0,
  Double_t lMultBoundHi = 100.0,
  Double_t lEEBoundLo = 0.0,
  Double_t lEEBoundHi = 100.0,
  TString lWhichMultEstimator = "SPDClusters",
  TString lWhichEffEnergyEstimator = "V0M",
  TString lV0Type = "Lambda", 
  TString fd = "NoFD",
  Bool_t lDoSystematics = kFALSE
  ){
  
  cout<<"Macro to test V0 analysis module"<<endl;

  cout<<"----------------------------------------------------"<<endl;
  cout<<"               V0 Analysis Macro "<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<" ---> Compiling needed class, please wait... "<<endl;
  //Compile Macro
  //Int_t workedornot = gSystem->CompileMacro("AliV0Module.cxx","-kfo");
 /* cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
   if( workedornot == 0 ){ 
      cout<<"*************************************"<<endl;
      cout<<" AliV0Module.cxx compilation failed! "<<endl;
      cout<<"*************************************"<<endl;
      return;
   }*/

  //Load Class
  gSystem->Load("AliV0Module_cxx");

  //Initialize Analysis Object
  AliV0Module *v0 = new AliV0Module(lV0Type);

  //Set data files
  TString lDataFilename = "../../AODAnalysis/grid/esdFile.root";
  //"../data/LHC17j.root";
  //"/home/fercoles/strg_analysis/FullStatistics/RootFiles/LHC18i.root" ; //0x0
  TString lMCFilename  = "../../AODAnalysis/MCESDmerged.root";
  //"../data/LHC20i2c_17j.root" ;//0x0;
  TString lTOFpercFilename  = "0x0";  //"ExtractedTOFPercentile_18i_TOFinfo.root";
  TString lFDFilename  = "0x0";  
  
  //PT
  Double_t ptbinlimits[] = {0.4,8.};
  // {0.4,  0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.5, 2.9, 3.4, 4, 5, 6.5, 8, 10};
  Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;;
  //Set Pt Bins Used in Analysis
  v0->SetPtBinLimits( ptbinnumb, ptbinlimits );

  //Estimators
  v0->SetPerformMultiplicityStudy(kTRUE);
  if (lWhichMultEstimator.Contains("TOF")) v0->SetTOFpercFileName ( lTOFpercFilename.Data()  );
  v0->SetLowMultValue(lMultBoundLo);
  v0->SetHighMultValue(lMultBoundHi);
  v0->SetLowEEValue(lEEBoundLo);
  v0->SetHighEEValue(lEEBoundHi);
  v0->SetWhichEstimators( lWhichMultEstimator, lWhichEffEnergyEstimator );
  
  //Set Default Cuts - note: Particle dependent
  v0->SetDefaultCuts(); 
  
  //v0->SetPtProCut(ptCutProt); ??????????????????????????????
  //if(lV0Type.Contains("Lambda")) v0->SetPtLimitForTOF(2.0); ???????????????????????????
  //
 
  //Feeddown treatment switch (applies only to Lambda, AntiLambda) 
  v0->SetFeeddownTreatment (fd);
  
  v0->SetGeant3FlukaCorr(kFALSE);
  if(lV0Type == "AntiLambda") v0->SetGeant3FlukaCorr(kTRUE);

  v0->SetUseIntegratedEff(kTRUE);
  v0->SetRealDataFile( lDataFilename.Data() );
  v0->SetMCDataFile ( lMCFilename.Data()   );
  v0->SetFeedDownDataFile  (lFDFilename.Data() );
  
  v0->SetOutputFile( Form("results/Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root", lV0Type.Data(), lWhichMultEstimator.Data(), 
    (lMultBoundLo), (lMultBoundHi), lWhichEffEnergyEstimator.Data(), (lEEBoundLo), (lEEBoundHi)) );

  //Run the analysis 
  if (lDoSystematics == kFALSE) v0->DoAnalysis();
  if (lDoSystematics == kTRUE){
  }
}
