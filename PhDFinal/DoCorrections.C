//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Funziona per SPDtrk0815trk0815 vs ZDC o vs V0M, SPDtrk0815 come mult estimator e ZDC o V0M come energy estimators.
// Basta mettere i due estimatori var o fixed (uno dei due DEVE essere SPDtrk0815trk0815) e fa tutto da solo.
// Ricordarsi di controllare i bin usati nell'event e signal loss correction, che matchino quelli dell'analisi
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include "CorrectionClass.h"

void Init(Int_t lClassCode,
		  TString lWhichEnergyEstimator,
		  TString lWhichMultEstimator,
		  Bool_t lDoMB,
		  TString fWhichParticle,
		  vector<Double_t> &percentileMult_low,
		  vector<Double_t> &percentileMult_high,
		  vector<Double_t> &percentileEnergy_low,
		  vector<Double_t> &percentileEnergy_high);

// classes
enum classname
{
	kStandalone = 0,
	kHighMult,
	kLowMult,
	kHighZN,
	kLowZN,
	kVeryLowZN
};

void DoCorrections(
	TString lCascType = "Omega",
	Int_t lClassCode = 0,
	TString lWhichEnergyEstimator = "V0M",
	TString lWhichMultEstimator = "SPDClusters",
	Bool_t DoMB = kTRUE)
{

	//Argument: System analysed - one of "XiMinus", "XiPlus", "OmegaMinus"
	cout<<"Macro to test Cascade analysis module"<<endl;
	cout<<"----------------------------------------------------"<<endl;
	cout<<"               Cascade Correction Macro "<<endl;
	cout<<"----------------------------------------------------"<<endl;
	cout<<endl;
	cout<<"----------------------------------------------------"<<endl;

	//Load Class
	gSystem->Load("CorrectionClass_cxx.so");

	//Initialize Analysis Object
	CorrectionClass *casc = new CorrectionClass(lCascType);

	if (!DoMB) {
		casc->SetOutputDataFile(Form("correctedspectra/CorrSpectra-%s-13TeV_%s_%s_class%i.root",
									 lCascType.Data(), lWhichMultEstimator.Data(), lWhichEnergyEstimator.Data(), lClassCode));
		casc->SetMCFileSignalLoss(Form("NormCorrections/SignalLoss-%s-13TeV_class%i.root", lCascType.Data(), lClassCode));
		casc->SetMCFileEventLoss(Form("NormCorrections/EventLoss-13TeV_class%i.root", lClassCode));
	} else{
		casc->SetOutputDataFile(Form("correctedspectra/CorrSpectra-%s-13TeV_INELgt0.root", lCascType.Data()));
		casc->SetMCFileSignalLoss(Form("NormCorrections/SignalLoss-%s-13TeV_INELgt0.root", lCascType.Data()));
		casc->SetMCFileEventLoss("NormCorrections/EventLoss-13TeV_INELgt0.root");
	}
	if (lCascType.Contains("Xi")) casc->SetMCFilepTshape("/home/fercoles/strg_analysis/PhDWork/Data/LHC15g3b1.root","/home/fercoles/strg_analysis/PhDWork/Data/LHC22e1_17jRunList.root", "/home/fercoles/strg_analysis/PhDWork/Data/LHC22e1_18iRunList.root");
	if (lCascType.Contains("Omega")) casc->SetMCFilepTshape("/home/fercoles/strg_analysis/PhDWork/Data/LHC15g3b2.root","/home/fercoles/strg_analysis/PhDWork/Data/LHC22e1_17jRunList.root", "/home/fercoles/strg_analysis/PhDWork/Data/LHC22e1_18iRunList.root");

	//=================
	//===== pT ========
	//=================
	Double_t* ptbinlimits;
	Long_t ptbinnumb;

	Double_t ptLambda[] = {0.4,  0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.5, 2.9, 3.4, 4, 5, 6.5, 8.};
    Long_t nLambda = sizeof(ptLambda)/sizeof(Double_t) - 1;
	Double_t ptXi[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
    Long_t nXi = sizeof(ptXi)/sizeof(Double_t) - 1;
	Double_t ptOmega[] = {0.90, 1.60, 2.20, 2.60, 3.00, 3.80, 5.50 };
    Long_t nOmega = sizeof(ptOmega)/sizeof(Double_t) - 1;
	Double_t ptK0S[] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.2,2.4,2.6,2.8,3.0,
						3.3,3.6,3.9,4.2,4.6,5,5.4,5.9, 6.5,7,7.5,8,8.5,9.2,10,11,12};
	Long_t nK0S = sizeof(ptK0S)/sizeof(Double_t) - 1;

	if (lCascType.Contains("Lambda")){
		ptbinlimits = ptLambda;
		ptbinnumb = nLambda;
	}
	if (lCascType.Contains("Xi")){
		ptbinlimits = ptXi;
		ptbinnumb = nXi;
	}
	if (lCascType.Contains("Omega")){
		ptbinlimits = ptOmega;
		ptbinnumb = nOmega;
	}
	if (lCascType.Contains("K0Short")){
		ptbinlimits = ptK0S;
		ptbinnumb = nK0S;
	}

	//Set Pt Bins
	casc->SetPtBinLimits( ptbinnumb, ptbinlimits );

	//=================
	//== Estimators ===
	//=================

	casc->SetMultEstimator(lWhichMultEstimator);
	casc->SetEnergyEstimator(lWhichEnergyEstimator);

	//Percentile
	vector<Double_t> percentileMult_low;
	vector<Double_t> percentileMult_high;
	vector<Double_t> percentileEnergy_low;
	vector<Double_t> percentileEnergy_high;

	Init(lClassCode, lWhichEnergyEstimator, lWhichMultEstimator, DoMB, lCascType, percentileMult_low, percentileMult_high, percentileEnergy_low, percentileEnergy_high);
	const Long_t percbinnumb = percentileMult_low.size();
	//
	casc->SetPercBinLimits(percbinnumb, percentileMult_low, percentileMult_high, percentileEnergy_low, percentileEnergy_high);

	//Which corrections?
	if (lCascType.Contains("Xi")){// || lCascType.Contains("Omega")){
		casc->DoMCShape(kTRUE);
	} else {
		casc->DoMCShape(kFALSE);
	}
	casc->DoEventLoss(kTRUE);
	casc->DoSignalLoss(kTRUE);
	casc->DoSystematics(kFALSE);

	casc->PrintConfiguration();
	//
	casc->DoAnalysis();
}

void Init(Int_t lClassCode,
		  TString lWhichEnergyEstimator,
		  TString lWhichMultEstimator,
		  Bool_t lDoMB,
		  TString fWhichParticle,
		  vector<Double_t> &percentileMult_low,
		  vector<Double_t> &percentileMult_high,
		  vector<Double_t> &percentileEnergy_low,
		  vector<Double_t> &percentileEnergy_high)
{
	// class 0 --> standalone
	Double_t percentileEnergy_low_0[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
	Double_t percentileEnergy_high_0[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
	Double_t percentileMult_low_0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t percentileMult_high_0[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
	Long_t n0 = sizeof(percentileEnergy_low_0) / sizeof(Double_t);

	// class 10 --> kStandalone for Omegas
	Double_t percentileEnergy_low_10[] = {0, 5, 15, 30, 50};
	Double_t percentileEnergy_high_10[] = {5, 15, 30, 50, 100};
	Double_t percentileMult_low_10[] = {0, 0, 0, 0, 0};
	Double_t percentileMult_high_10[] = {100, 100, 100, 100, 100};
	Long_t n10 = sizeof(percentileEnergy_low_10) / sizeof(Double_t);

	// class 1 --> kHighMult
	Double_t percentileMult_low_1[] = {10, 10, 10, 10, 10, 10, 10};
	Double_t percentileMult_high_1[] = {20, 20, 20, 20, 20, 20, 20};
	Double_t percentileEnergy_low_1[] = {0, 5, 10, 20, 30, 40, 50};
	Double_t percentileEnergy_high_1[] = {5, 10, 20, 30, 40, 50, 100};
	Long_t n1 = sizeof(percentileMult_low_1) / sizeof(Double_t);

	// class 2 --> kLowMult
	Double_t percentileMult_low_2[] = {40, 40, 40, 40, 40, 40, 40};
	Double_t percentileMult_high_2[] = {50, 50, 50, 50, 50, 50, 50};
	Double_t percentileEnergy_low_2[] = {0, 20, 30, 40, 50, 60, 70};
	Double_t percentileEnergy_high_2[] = {20, 30, 40, 50, 60, 70, 100};
	Long_t n2 = sizeof(percentileMult_low_2) / sizeof(Double_t);

	// class 3 --> fixed high ZN
	Double_t percentileMult_low_3[] = {10, 40, 60, 70, 80};
	Double_t percentileMult_high_3[] = {40, 60, 70, 80, 100};
	Double_t percentileEnergy_low_3[] = {70, 60, 40, 40, 40};
	Double_t percentileEnergy_high_3[] = {100, 100, 100, 80, 70};
	Long_t n3 = sizeof(percentileMult_low_3) / sizeof(Double_t);

	// class 4 --> fixed low ZN
	Double_t percentileMult_low_4[] = {0, 10, 20, 30, 50};
	Double_t percentileMult_high_4[] = {20, 30, 40, 50, 100};
	Double_t percentileEnergy_low_4[] = {40, 30, 30, 20, 0};
	Double_t percentileEnergy_high_4[] = {60, 70, 50, 50, 30};
	Long_t n4 = sizeof(percentileMult_low_4) / sizeof(Double_t);

	// class 5 --> fixed very low ZN
	Double_t percentileMult_low_5[] = {0, 10, 20, 30};
	Double_t percentileMult_high_5[] = {10, 20, 30, 50};
	Double_t percentileEnergy_low_5[] = {20, 10, 0, 0};
	Double_t percentileEnergy_high_5[] = {30, 30, 20, 10};
	Long_t n5 = sizeof(percentileMult_low_5) / sizeof(Double_t);

	if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 0) {
			if (!lDoMB) {
				for (Int_t i = 0; i < n0; i++) {
					percentileMult_low.push_back(percentileMult_low_0[i]);
					percentileMult_high.push_back(percentileMult_high_0[i]);
					percentileEnergy_low.push_back(percentileEnergy_low_0[i]);
					percentileEnergy_high.push_back(percentileEnergy_high_0[i]);
				}
			} else {
			for (Int_t i = 0; i < 2; i++){
					percentileMult_low.push_back(0.);
					percentileMult_high.push_back(100.);
					percentileEnergy_low.push_back(0.);
					percentileEnergy_high.push_back(100.);
			}
		}
	}

	if (lClassCode == 10)
	{
		for (Int_t i = 0; i < n10; i++)
		{
			percentileMult_low.push_back(percentileMult_low_10[i]);
			percentileMult_high.push_back(percentileMult_high_10[i]);
			percentileEnergy_low.push_back(percentileEnergy_low_10[i]);
			percentileEnergy_high.push_back(percentileEnergy_high_10[i]);
		}
		cout << "------------------------------------------" << endl;
		cout << " Initializing kStandalone class for Omegas..." << endl;
	}

	if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 1)
	{
		for (Int_t i = 0; i < n1; i++)
		{
			percentileMult_low.push_back(percentileMult_low_1[i]);
			percentileMult_high.push_back(percentileMult_high_1[i]);
			percentileEnergy_low.push_back(percentileEnergy_low_1[i]);
			percentileEnergy_high.push_back(percentileEnergy_high_1[i]);
		}
	}

	if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 2)
	{
		for (Int_t i = 0; i < n2; i++)
		{
			percentileMult_low.push_back(percentileMult_low_2[i]);
			percentileMult_high.push_back(percentileMult_high_2[i]);
			percentileEnergy_low.push_back(percentileEnergy_low_2[i]);
			percentileEnergy_high.push_back(percentileEnergy_high_2[i]);
		}
	}

	if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 3) {
		for (Int_t i = 0; i < n3; i++) {
			percentileMult_low.push_back(percentileMult_low_3[i]);
			percentileMult_high.push_back(percentileMult_high_3[i]);
			percentileEnergy_low.push_back(percentileEnergy_low_3[i]);
			percentileEnergy_high.push_back(percentileEnergy_high_3[i]);
		}
	}

	if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 4) {
		for (Int_t i = 0; i < n4; i++) {
			percentileMult_low.push_back(percentileMult_low_4[i]);
			percentileMult_high.push_back(percentileMult_high_4[i]);
			percentileEnergy_low.push_back(percentileEnergy_low_4[i]);
			percentileEnergy_high.push_back(percentileEnergy_high_4[i]);
		}
	}

	if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 5) {
		for (Int_t i = 0; i < n5; i++) {
			percentileMult_low.push_back(percentileMult_low_5[i]);
			percentileMult_high.push_back(percentileMult_high_5[i]);
			percentileEnergy_low.push_back(percentileEnergy_low_5[i]);
			percentileEnergy_high.push_back(percentileEnergy_high_5[i]);
		}
	}
}