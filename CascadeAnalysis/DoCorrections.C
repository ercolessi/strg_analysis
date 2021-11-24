//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Funziona per SPDClusters vs ZDC o vs V0M, SPD come mult estimator e ZDC o V0M come energy estimators.
// Basta mettere i due estimatori var o fixed (uno dei due DEVE essere SPDClusters) e fa tutto da solo.
// Ricordarsi di controllare i bin usati nell'event e signal loss correction, che matchino quelli dell'analisi
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

void DoCorrections(
	TString lCascType = "Xi",
	Double_t lFixedLo = 0.0,
	Double_t lFixedHi = 100.0,
	TString lWhichVarEstimator = "SPDClusters",
	TString lWhichFixedEstimator = "V0M")
	{
	
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
	cout<<"Macro to test Cascade analysis module"<<endl;
	cout<<"----------------------------------------------------"<<endl;
	cout<<"               Cascade Correction Macro "<<endl;
	cout<<"----------------------------------------------------"<<endl;
	cout<<endl;
	cout<<"----------------------------------------------------"<<endl;
	cout<<" ---> Compiling needed class, please wait... "<<endl;
	
	//Load Class
	gSystem->Load("CorrectionClass_cxx.so");

	//Initialize Analysis Object
	CorrectionClass *casc = new CorrectionClass(lCascType);

	casc->SetOutputDataFile("testcorrclass.root");
	casc->SetMCFilepTshape("~/strg_analysis/CascadesResults/ANALISI/MC15g3b1_ZDCRuns.root");
	casc->SetMCFileNormCorr( Form("Norm-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", 
                                lCascType.Data(), lWhichVarEstimator.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi)
							);

	//=================
	//===== pT ========
	//=================
	
	//Omega
	//Double_t ptbinlimits[]   = {0.90, 1.60, 2.20, 2.60, 3.00, 3.80, 5.50, 8.00, 12.0 }; 
	//Xi
	Double_t ptbinlimits[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
	Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;

	//Set Pt Bins 
	casc->SetPtBinLimits( ptbinnumb, ptbinlimits );

	//=================
	//== Estimators ===
	//=================

	//Set variables
	casc->SetVarEstimator(lWhichVarEstimator);
	casc->SetFixedEstimator(lWhichFixedEstimator);
	casc->SetMultEstimator("SPDClusters");
	if (lWhichVarEstimator.Contains("V0M") || lWhichFixedEstimator.Contains("V0M")) casc->SetEnergyEstimator("V0M");
	if (lWhichVarEstimator.Contains("ZDC") || lWhichFixedEstimator.Contains("ZDC")) casc->SetEnergyEstimator("ZDC");
	
	//Set percentile classes
	Double_t percentile[] = {0,100};
	const Long_t percentilebinnumb = sizeof(percentile)/sizeof(Double_t) - 1;
	//
	casc->SetPercBinLimits(percentilebinnumb, percentile);

	//Set fixed percentile class
	casc->SetFixedPercLimits(lFixedLo,lFixedHi);

	//Which corrections?
	casc->DoMCShape(kTRUE);
	casc->DoEventLoss(kTRUE);
	casc->DoSignalLoss(kTRUE);
	
	casc->PrintConfiguration();
	//
	casc->DoAnalysis();
}
