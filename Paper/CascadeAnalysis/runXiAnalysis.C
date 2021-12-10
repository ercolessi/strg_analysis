#include <iostream>
#include <fstream>
#include "AliCascadeModule.h"

void runXiAnalysis(
  Double_t lMultBoundLo = 0.0,
  Double_t lMultBoundHi = 100.0,
  Double_t lEEBoundLo = 0.0,
  Double_t lEEBoundHi = 100.0,
  TString lCascType = "XiMinus",
  TString lWhichMultEstimator = "SPDClusters",
  TString lWhichEffEnergyEstimator = "V0M",
  Bool_t lDoSystematics = kFALSE,
  Bool_t  lSweepFlag = kFALSE){

  TString lWhichSystVar = "";
  TString lWhichMult = Form("%05.0lfto%05.0lf", 100.*lMultBoundLo, 100.*lMultBoundHi); 

  //Set data files
  TString lDataFilename = "../data/LHC17j.root";
  TString lMCFilename  = "../data/LHC20i2c_17j.root";
  TString lTOFpercFilename  = "0x0";  //"ExtractedTOFPercentile_18i_TOFinfo.root";
  
  // Load common libraries
  /*gSystem->Load("libCore.so");
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
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");*/

  //Argument: System analysed - one of "XiMinus", "XiPlus", "OmegaMinus"
  cout<<"Macro to test Cascade analysis module"<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<"               Cascade Analysis Macro "<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<" ---> Compiling needed class, please wait... "<<endl;

  //Compile Macro
  //gSystem->CompileMacro("AliCascadeModule.cxx","-kfo");

  //Load Class
  gSystem->Load("AliCascadeModule_cxx.so");

  //Initialize Analysis Object
  AliCascadeModule *casc = new AliCascadeModule(lCascType);

  //Set analysis files
  casc->SetRealDataFile( lDataFilename.Data() );
  casc->SetMCDataFile ( lMCFilename.Data()   );
  
  //Omega
  //Double_t ptbinlimits[]   = {0.90, 1.60, 2.20, 2.60, 3.00, 3.80, 5.50, 8.00, 12.0 }; 
  //Xi
  Double_t ptbinlimits[] = {0.6, 1.0, 1.5, 2.0, 2.9, 3.4, 6.5};
  //{0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
  Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;
  
  //Set Pt Bins 
  casc->SetPtBinLimits( ptbinnumb, ptbinlimits );
  casc->SetRapidityWindow(-0.5, 0.5);
 
  // Default cuts configuration (Fiorella's)
  casc->SetDefaultCuts();
 
  //Set estimators configurations
  if (lWhichMultEstimator.Contains("TOF")) casc->SetTOFpercFileName ( lTOFpercFilename.Data()  );
  casc->SetLowMultValue(lMultBoundLo);
  casc->SetHighMultValue(lMultBoundHi);
  casc->SetLowEEValue(lEEBoundLo);
  casc->SetHighEEValue(lEEBoundHi);
  casc->SetPerformMultiplicityStudy( kTRUE );
  casc->SetPerformEEStudy( kTRUE );
  casc->SetWhichEstimators( lWhichMultEstimator, lWhichEffEnergyEstimator );
  casc->SetUseIntegratedEfficiency( kTRUE );

  // Geant-Fluka correction for anti-protons 
  if( lCascType.Contains("Plus") ) {
      TFile* fGeantFluka = new TFile("antiPrCorrFunc.root", "READ");
      TF1* funcGeantFlukaCorr = (TF1*) fGeantFluka->Get("funcCorrAntiProtonGEANT4");
      //TF1* funcGeantFlukaCorr = (TF1*)fGeantFluka->Get("funcCorrAntiProtonFLUKA"); // used for systematics
      casc->SetGeantFlukaCorrection ( funcGeantFlukaCorr );
  }

  //Set output file name
  cout << "[runCascadeAnalysis] Performing multiplicity selection: [" << lMultBoundLo << "," << lMultBoundHi << "[ (" << lWhichMultEstimator << ")" << endl;
  cout << "[runCascadeAnalysis] Performing effective energy selection: [" << lEEBoundLo << "," << lEEBoundHi << "[ (" << lWhichEffEnergyEstimator << ")" << endl;
  
  casc->SetOutputFile( Form("Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f.root", lCascType.Data(), lWhichMultEstimator.Data(), 
    (lMultBoundLo), (lMultBoundHi), lWhichEffEnergyEstimator.Data(), (lEEBoundLo), (lEEBoundHi)) );

    //Run the analysis 
    if (lDoSystematics == kFALSE) casc->DoAnalysis();
    if (lDoSystematics == kTRUE){

        TString lSystFile = "Results-Systematics";
        lSystFile.Append( Form("-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f-",lCascType.Data(),lWhichMultEstimator.Data(),
           (lMultBoundLo),(lMultBoundHi),lWhichEffEnergyEstimator.Data(),(lEEBoundLo),(lEEBoundHi)) );

        //---------------------------------------------------------------------------------------------
        // Topological Selection Variables Systematics

        ////////////////////////////////////////////////////////////////////////////
        // V0RADIUS ////////////////////////////////////////////////////////////////
        lWhichSystVar = "V0Radius";
      	cout<<endl<<"---> Performing Systematics Studies: V0 Radius"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutV0Radius(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Radius(1.10); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Radius(2.50); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Radius(5.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutV0Radius(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Radius(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Radius(2.50); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Radius(6.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
          
        
      	////////////////////////////////////////////////////////////////////////////
      	// CASCRADIUS //////////////////////////////////////////////////////////////
      	lWhichSystVar = "CascRadius";
        casc->SetDefaultCuts(); //Do this before starting again
      	cout<<endl<<"---> Performing Systematics Studies: Cascade Radius"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutCascRadius(0.40); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutCascRadius(0.50); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutCascRadius(0.80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutCascRadius(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutCascRadius(0.30); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutCascRadius(0.40); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutCascRadius(0.60); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutCascRadius(0.70); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
              
        ////////////////////////////////////////////////////////////////////////////
        // DCABACHTOPV /////////////////////////////////////////////////////////////
        lWhichSystVar = "DCABachToPV";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: DCA Bachelor to Primary Vertex"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutDCABachToPV(0.03); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCABachToPV(0.03); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCABachToPV(0.10); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCABachToPV(0.17); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutDCABachToPV(0.03); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCABachToPV(0.03); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCABachToPV(0.07); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCABachToPV(0.10); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
              
        ////////////////////////////////////////////////////////////////////////////
        // DCAV0TOPV ///////////////////////////////////////////////////////////////
        lWhichSystVar = "DCAV0ToPV";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: DCA V0 to Primary Vertex"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutDCAV0ToPV(0.05); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0ToPV(0.05); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0ToPV(0.10); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0ToPV(0.15); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutDCAV0ToPV(0.05); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0ToPV(0.05); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0ToPV(0.08); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0ToPV(0.10); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
              
        ////////////////////////////////////////////////////////////////////////////
        // DCANEGTOPV //////////////////////////////////////////////////////////////
        lWhichSystVar = "DCANegToPV";
        casc->SetDefaultCuts(); //Do this before starting again
      	cout<<endl<<"---> Performing Systematics Studies: DCA Negative to Primary Vertex"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.03 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.15 : 0.09 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.30 : 0.11 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.03 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.10 : 0.05 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.30 : 0.10 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
              
        ////////////////////////////////////////////////////////////////////////////
        // DCAPOSTOPV //////////////////////////////////////////////////////////////
        lWhichSystVar = "DCAPosToPV";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: DCA Positive to Primary Vertex"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.03 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.09 : 0.15 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.11 : 0.30 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.03 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.05 : 0.10 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.10 : 0.30 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
              
        ////////////////////////////////////////////////////////////////////////////
        // DCAV0DAUGHTERS //////////////////////////////////////////////////////////
        lWhichSystVar = "DCAV0Daughters";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: DCA V0 Daughters"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutDCAV0Daughters(2.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0Daughters(1.80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0Daughters(1.20); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0Daughters(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutDCAV0Daughters(2.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0Daughters(1.80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0Daughters(1.30); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCAV0Daughters(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
              
        ////////////////////////////////////////////////////////////////////////////
        // DCACASCDAUGHTERS ////////////////////////////////////////////////////////
        lWhichSystVar = "DCACascDaughters";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: DCA Cascade Daughters"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutDCACascDaughters(2.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCACascDaughters(1.80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCACascDaughters(1.20); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCACascDaughters(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutDCACascDaughters(2.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutDCACascDaughters(1.80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutDCACascDaughters(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutDCACascDaughters(0.60); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
          
        ////////////////////////////////////////////////////////////////////////////
        // V0COSPA /////////////////////////////////////////////////////////////////
        lWhichSystVar = "V0CosPA";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: Cosine of Pointing Angle of the V0"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutV0CosPA(0.95); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutV0CosPA(0.96); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutV0CosPA(0.98); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutV0CosPA(0.99); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutV0CosPA(0.95); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutV0CosPA(0.96); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutV0CosPA(0.98); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutV0CosPA(0.99); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
          
        ////////////////////////////////////////////////////////////////////////////
        // CASCCOSPA //////////////////////////////////////////////////////////////
        lWhichSystVar = "CascCosPA";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: Cosine of Pointing Angle of the Cascade"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutCascCosPA(0.95); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutCascCosPA(0.96); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutCascCosPA(0.98); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutCascCosPA(0.99); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutCascCosPA(0.95); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutCascCosPA(0.96); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutCascCosPA(0.98); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutCascCosPA(0.99); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
              
        ////////////////////////////////////////////////////////////////////////////
        // V0MASS //////////////////////////////////////////////////////////////////
        lWhichSystVar = "V0Mass";
        casc->SetDefaultCuts(); //Do this before starting again
      	cout<<endl<<"---> Performing Systematics Studies: V0 Mass Window"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutV0Mass(0.010); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Mass(0.009); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Mass(0.007); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Mass(1.006); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutV0Mass(0.010); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Mass(0.009); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Mass(0.007); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutV0Mass(0.006); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	}
              
              
        ////////////////////////////////////////////////////////////////////////////
        //---------------------------------------------------------------------------------------------
        // Other Selection Variables Systematics
        ////////////////////////////////////////////////////////////////////////////
        // TPCPIDNSIGMAS ///////////////////////////////////////////////////////////
        lWhichSystVar = "TPCPIDNSigmas";
        casc->SetDefaultCuts(); //Do this before starting again
      	cout<<endl<<"---> Performing Systematics Studies: TPC PID Nsigmas"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutTPCPIDNSigmas(7.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutTPCPIDNSigmas(6.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutTPCPIDNSigmas(4.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  //casc->SetCutTPCPIDNSigmas(4.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutTPCPIDNSigmas(7.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutTPCPIDNSigmas(6.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutTPCPIDNSigmas(4.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  //casc->SetCutTPCPIDNSigmas(4.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
      	}
          
        ////////////////////////////////////////////////////////////////////////////
        // TPCNCLUSTERS ///////////////////////////////////////////////////////////
        lWhichSystVar = "TPCNClusters";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: TPC Least Number of Clusters"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutLeastNumberOfClusters(70); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  //casc->SetCutLeastNumberOfClusters(70); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutLeastNumberOfClusters(75); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutLeastNumberOfClusters(80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutLeastNumberOfClusters(70); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  //casc->SetCutLeastNumberOfClusters(70); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutLeastNumberOfClusters(73); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  casc->SetCutLeastNumberOfClusters(76); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
      	}
              
        ////////////////////////////////////////////////////////////////////////////
        // PROPERLIFETIME //////////////////////////////////////////////////////////
        lWhichSystVar = "ProperLifetime";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: Proper Lifetime"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutProperLifetime(5.0*4.91); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutProperLifetime(4.0*4.91); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutProperLifetime(2.5*4.91); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  //casc->SetCutProperLifetime(2.5*4.91); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutProperLifetime(5.0*2.461); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  casc->SetCutProperLifetime(4.0*2.461); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutProperLifetime(3.0*2.461); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  //casc->SetCutProperLifetime(3.0*2.461); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
      	}
              
        ////////////////////////////////////////////////////////////////////////////
        // COMPETINGSPECIES ////////////////////////////////////////////////////////
        lWhichSystVar = "CompetingSpecies";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: Competing Species"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutCompetingSpecies(-1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  //casc->SetCutCompetingSpecies(-1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  //casc->SetCutCompetingSpecies(-1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  //casc->SetCutCompetingSpecies(-1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
      	  gSystem->Exec(Form("cp %s-1.root %s-3.root", systname.Data(), systname.Data()));
      	  gSystem->Exec(Form("cp %s-1.root %s-4.root", systname.Data(), systname.Data()));
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutCompetingSpecies( -1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  //casc->SetCutCompetingSpecies( -1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutCompetingSpecies(0.008); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  //casc->SetCutCompetingSpecies(0.008); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
      	  gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
      	}
          
        ////////////////////////////////////////////////////////////////////////////
        // SIGEXTBINCOUNT //////////////////////////////////////////////////////////
        lWhichSystVar = "SigExtBinCount";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: Sigma Cut for Signal Extraction (bin count)"<<endl;
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutSigmaForSignalExtraction(5.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  //casc->SetCutSigmaForSignalExtraction(5.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutSigmaForSignalExtraction(3.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  //casc->SetCutSigmaForSignalExtraction(3.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
      	  gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutSigmaForSignalExtraction(4.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  //casc->SetCutSigmaForSignalExtraction(4.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutSigmaForSignalExtraction(3.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  //casc->SetCutSigmaForSignalExtraction(3.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
      	  gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
      	}
              
        ////////////////////////////////////////////////////////////////////////////
        // SIGEXTFITTING ///////////////////////////////////////////////////////////
        lWhichSystVar = "SigExtFitting";
      	casc->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: Sigma Cut for Signal Extraction (fitting)"<<endl;
      	casc->SetFitBackground( kTRUE );
      	if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
      	  casc->SetCutSigmaForSignalExtraction(5.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  //casc->SetCutSigmaForSignalExtraction(5.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutSigmaForSignalExtraction(3.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  //casc->SetCutSigmaForSignalExtraction(3.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
      	  gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
      	}
      	if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
      	  casc->SetCutSigmaForSignalExtraction(4.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
      	  //casc->SetCutSigmaForSignalExtraction(4.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
      	  casc->SetCutSigmaForSignalExtraction(3.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
      	  //casc->SetCutSigmaForSignalExtraction(3.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
      	  //
      	  TString systname = lSystFile + lWhichSystVar;
      	  gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
      	  gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
      	}
    
        ////////////////////////////////////////////////////////////////////////////
        // GEANTFLUKACORRECTION ////////////////////////////////////////////////////
        if( lWhichSystVar == "GeantFlukaCorrection" ) {
            cout<<endl<<"---> Performing Systematics Studies: Geant-Fluka correction"<<endl;
            if( lCascType == "XiPlus" || lCascType == "OmegaPlus" ) {
                TFile* fGeantFluka = new TFile("antiPrCorrFunc.root", "READ");
                TF1* funcGeantFlukaCorr = (TF1*)fGeantFluka->Get("funcCorrAntiProtonFLUKA"); // used for systematics
                casc->SetGeantFlukaCorrection( funcGeantFlukaCorr ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                //casc->SetGeantFlukaCorrection( funcGeantFlukaCorr ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                //casc->SetGeantFlukaCorrection( funcGeantFlukaCorr ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetGeantFlukaCorrection( funcGeantFlukaCorr ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
                gSystem->Exec(Form("cp %s-1.root %s-3.root", systname.Data(), systname.Data()));
                gSystem->Exec(Form("cp %s-1.root %s-4.root", systname.Data(), systname.Data()));
            }
        }
        ////////////////////////////////////////////////////////////////////////////
    }
    //---------------------------------------------------------------------------------------------
  }
