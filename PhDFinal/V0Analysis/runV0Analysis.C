#include "AliV0Module.h"

void runV0Analysis(
    Double_t lMultBoundLo = 0.0,
    Double_t lMultBoundHi = 100.0,
    Double_t lEEBoundLo = 0.0,
    Double_t lEEBoundHi = 100.0,
    TString lV0Type = "K0Short",
    TString lWhichMultEstimator = "SPDClusters",
    TString lWhichEffEnergyEstimator = "V0M",
    Int_t lClassCode = 0,
    TString fd = "UseMCRatio", //"UseMCRatio",//"DoubleChargedXi", //NoFD
    Bool_t lDoSystematics = kFALSE)
{

  cout<<"Macro to test V0 analysis module"<<endl;

  cout<<"----------------------------------------------------"<<endl;
  cout<<"               V0 Analysis Macro "<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<" ---> Compiling needed class, please wait... "<<endl;

  //Load Class
  gSystem->Load("AliV0Module_cxx");

  //Initialize Analysis Object
  AliV0Module *v0 = new AliV0Module(lV0Type);

  //Set data files
  TString lDataFilename = "~/strg_analysis/AnalisiFinale/data/real/LHCFull_pass2.root";
  TString lMCFilename0  =
      "/home/fercoles/strg_analysis/PhDWork/Data/LHC15g3a3.root" ;//0x0;
  TString lMCFilename1 =
      "/home/fercoles/strg_analysis/PhDWork/Data/LHC22k1_17jRunList.root"; // 0x0;
  TString lMCFilename2 =
      "/home/fercoles/strg_analysis/PhDWork/Data/LHC22k1_18iRunList.root"; // 0x0;
  TString lTOFpercFilename  = "0x0";  //"ExtractedTOFPercentile_18i_TOFinfo.root";
  TString lpercnameN0815 = "/home/fercoles/strg_analysis/PhDWork/Selections/PercentileCalibration.root";

  TString lFDFilename = Form("/home/fercoles/strg_analysis/PhDFinal/correctedspectra/CorrSpectra-%s-13TeV_%s_%s_class%i.root",
                             "Xi", lWhichMultEstimator.Data(), lWhichEffEnergyEstimator.Data(), lClassCode);
  if (lMultBoundLo < 1E-10 && lMultBoundHi > 99 && lEEBoundLo < 1E-10 && lEEBoundHi > 99)
    lFDFilename = "/home/fercoles/strg_analysis/PhDFinal/correctedspectra/CorrSpectra-Xi-13TeV_INELgt0.root";

  //PT
  Double_t* ptbinlimits;
	Long_t ptbinnumb;

	Double_t ptLambda[] = {0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.5, 2.9, 3.4, 4, 5, 6.5, 8, 10};
  Long_t nLambda = sizeof(ptLambda)/sizeof(Double_t) - 1;

  Double_t ptK0S[] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6., 8., 10.};
  //{0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.2,2.4,2.6,2.8,3.0,
   // 3.3,3.6,3.9,4.2,4.6,5,5.4,5.9, 6.5,7,7.5,8,8.5,9.2,10,11,12,13.5,15};
  Long_t nK0S = sizeof(ptK0S)/sizeof(Double_t) - 1;

	if (lV0Type.Contains("Lambda")){
		ptbinlimits = ptLambda;
		ptbinnumb = nLambda;
	}

	if (lV0Type.Contains("K0Short")){
		ptbinlimits = ptK0S;
		ptbinnumb = nK0S;
	}

  //Set Pt Bins Used in Analysis
  v0->SetPtBinLimits( ptbinnumb, ptbinlimits );

  //Estimators
  v0->SetPerformMultiplicityStudy(kTRUE);
  if (lWhichMultEstimator.Contains("TOF")) v0->SetTOFpercFileName ( lTOFpercFilename.Data()  );
  if (lWhichMultEstimator.Contains("0815")) v0->SetN0815percFileName ( lpercnameN0815.Data()  );
  v0->SetLowMultValue(lMultBoundLo);
  v0->SetHighMultValue(lMultBoundHi);
  v0->SetLowEEValue(lEEBoundLo);
  v0->SetHighEEValue(lEEBoundHi);
  v0->SetWhichEstimators( lWhichMultEstimator, lWhichEffEnergyEstimator );

  //Set Default Cuts - note: Particle dependent
  v0->SetDefaultCuts();
  v0->SetDoEfficiency(kTRUE);

  //Feeddown treatment switch (applies only to Lambda, AntiLambda)
  v0->SetFeeddownTreatment (fd);

  v0->SetGeant3FlukaCorr(kFALSE);
  if(lV0Type == "AntiLambda") v0->SetGeant3FlukaCorr(kTRUE);

  v0->SetUseIntegratedEff(kTRUE);
  v0->SetRealDataFile( lDataFilename.Data() );
  v0->SetMCDataFile ( lMCFilename0.Data(), lMCFilename1.Data(), lMCFilename2.Data()   );
  v0->SetFeedDownDataFile  (lFDFilename.Data() );

  //weights for efficiency
  Double_t w0 = 3.591657E+7;
  Double_t w1 = 2.984665E+7;
  Double_t w2 = 5.014868E+7;
  v0->SetEventsWeights(w0,w1,w2);

  TString lWhichSystVar = "";

 // v0->SetOutputFile( Form("Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f_FD%s_TOF2GeV.root", lV0Type.Data(), lWhichMultEstimator.Data(),
  //  (lMultBoundLo), (lMultBoundHi), lWhichEffEnergyEstimator.Data(), (lEEBoundLo), (lEEBoundHi), fd.Data()) );

  v0->SetOutputFile( Form("Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f_FD%s.root", lV0Type.Data(), lWhichMultEstimator.Data(),
    (lMultBoundLo), (lMultBoundHi), lWhichEffEnergyEstimator.Data(), (lEEBoundLo), (lEEBoundHi), fd.Data()) );

  //Run the analysis
  if (lDoSystematics == kFALSE) v0->DoAnalysis();
  if (lDoSystematics == kTRUE){

    v0->DoAnalysis();

    TString lSystFile = "Results-Systematics";
        lSystFile.Append( Form("-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f-",lV0Type.Data(),lWhichMultEstimator.Data(),
           (lMultBoundLo),(lMultBoundHi),lWhichEffEnergyEstimator.Data(),(lEEBoundLo),(lEEBoundHi)) );

        //---------------------------------------------------------------------------------------------
        // Topological Selection Variables Systematics

        ////////////////////////////////////////////////////////////////////////////
        // V0RADIUS ////////////////////////////////////////////////////////////////
        lWhichSystVar = "V0Radius";
      	cout<<endl<<"---> Performing Systematics Studies: V0 Radius"<<endl;
      	v0->SetCutV0Radius(0.30); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
      	v0->SetCutV0Radius(0.40); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
      	v0->SetCutV0Radius(0.60); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
      	v0->SetCutV0Radius(0.70); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();

        ////////////////////////////////////////////////////////////////////////////
        // V0COSPA /////////////////////////////////////////////////////////////////
        lWhichSystVar = "V0CosPA";
      	v0->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: Cosine of Pointing Angle of the V0"<<endl;
    	  if (lV0Type.Contains("Lambda")){
          v0->SetCutV0CosPA(0.993); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
          v0->SetCutV0CosPA(0.994); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
          v0->SetCutV0CosPA(0.996); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
          v0->SetCutV0CosPA(0.997); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();
        }
        if (lV0Type.Contains("K0Short")){
          v0->SetCutV0CosPA(0.95); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
          v0->SetCutV0CosPA(0.96); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
          v0->SetCutV0CosPA(0.98); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
          v0->SetCutV0CosPA(0.99); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();
        }

        lWhichSystVar = "DCAPosToPV";
        v0->SetDefaultCuts(); //Do this before starting again
      	cout<<endl<<"---> Performing Systematics Studies: DCA Negative to Primary Vertex"<<endl;
    	  v0->SetCutDCAPosToPV( 0.05 ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
    	  v0->SetCutDCAPosToPV( 0.055 ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
    	  v0->SetCutDCAPosToPV( 0.07 ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
    	  v0->SetCutDCAPosToPV( 0.08 ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();

      	lWhichSystVar = "DCANegToPV";
        v0->SetDefaultCuts(); //Do this before starting again
      	cout<<endl<<"---> Performing Systematics Studies: DCA Negative to Primary Vertex"<<endl;
    	  v0->SetCutDCANegToPV( 0.05 ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
    	  v0->SetCutDCANegToPV( 0.055 ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
    	  v0->SetCutDCANegToPV( 0.07 ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
    	  v0->SetCutDCANegToPV( 0.08 ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();

        lWhichSystVar = "DCAV0Daughters";
      	v0->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: DCA V0 Daughters"<<endl;
    	  v0->SetCutDCAV0Daughters(1.50); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
    	  v0->SetCutDCAV0Daughters(1.25); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
    	  v0->SetCutDCAV0Daughters(0.75); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
    	  v0->SetCutDCAV0Daughters(0.50); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();

        ////////////////////////////////////////////////////////////////////////////
        // PROPERLIFETIME //////////////////////////////////////////////////////////
        lWhichSystVar = "ProperLifetime";
      	v0->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: Proper Lifetime"<<endl;
      	if (lV0Type.Contains("Lambda")){
          v0->SetCutProperLifetime(40); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
          v0->SetCutProperLifetime(20); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
        }
        if (lV0Type.Contains("K0Short")){
          v0->SetCutProperLifetime(40); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
          v0->SetCutProperLifetime(30); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
          v0->SetCutProperLifetime(12); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
        }

        //v0->SetCutProperLifetime(2.5*4.91); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
    	  //v0->SetCutProperLifetime(2.5*4.91); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();
    	  //
    	  TString systname = lSystFile + lWhichSystVar;
    	  if (lV0Type.Contains("Lambda")){
          gSystem->Exec(Form("cp %s-2.root %s-1.root", systname.Data(), systname.Data()));
          gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
        }
        if (lV0Type.Contains("K0Short")){
          gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
        }

        ////////////////////////////////////////////////////////////////////////////
        // TPCNCLUSTERS ///////////////////////////////////////////////////////////
        lWhichSystVar = "TPCNClusters";
      	v0->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: TPC Least Number of Clusters"<<endl;
     	  //v0->SetCutLeastNumberOfClusters(70); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
     	  //v0->SetCutLeastNumberOfClusters(70); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
     	  v0->SetCutLeastNumberOfCrossedRows(75); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
     	  v0->SetCutLeastNumberOfCrossedRows(80); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();
     	  //
     	  systname = lSystFile + lWhichSystVar;
     	  gSystem->Exec(Form("cp %s-3.root %s-1.root", systname.Data(), systname.Data()));
        gSystem->Exec(Form("cp %s-3.root %s-2.root", systname.Data(), systname.Data()));

         // TPCPIDNSIGMAS ///////////////////////////////////////////////////////////
        lWhichSystVar = "TPCPIDNSigmas";
        v0->SetDefaultCuts(); //Do this before starting again
      	cout<<endl<<"---> Performing Systematics Studies: TPC PID Nsigmas"<<endl;
    	  v0->SetCutTPCPIDNSigmas(7.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
    	  v0->SetCutTPCPIDNSigmas(6.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
    	  v0->SetCutTPCPIDNSigmas(4.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
    	  //v0->SetCutTPCPIDNSigmas(4.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();
    	  //
    	  systname = lSystFile + lWhichSystVar;
    	  gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));


        ////////////////////////////////////////////////////////////////////////////
        // SIGEXTBINCOUNT //////////////////////////////////////////////////////////
        lWhichSystVar = "SigExtBinCount";
      	v0->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: Sigma Cut for Signal Extraction (bin count)"<<endl;
     	  //v0->SetCutSigmaForSignalExtraction(5.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
     	  v0->SetCutSigmaForSignalExtraction(7.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
     	  v0->SetCutSigmaForSignalExtraction(5.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
     	  v0->SetCutSigmaForSignalExtraction(4.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();
     	  //
     	  systname = lSystFile + lWhichSystVar;
     	  gSystem->Exec(Form("cp %s-2.root %s-1.root", systname.Data(), systname.Data()));


        lWhichSystVar = "CompetingRejection";
      	v0->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: Competing V0 rejection "<<endl;
     	  if (lV0Type.Contains("Lambda")){
          v0->SetCutCompetingV0Rejection(0.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
     	  }
        if (lV0Type.Contains("K0Short")){
          v0->SetCutCompetingV0Rejection(3.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
          v0->SetCutCompetingV0Rejection(5.5); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
     	  }
        //v0->SetCutCompetingV0Rejection(5.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
     	  //v0->SetCutCompetingV0Rejection(3.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
     	  //v0->SetCutCompetingV0Rejection(3.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();
     	  //
     	  systname = lSystFile + lWhichSystVar;
     	  if (lV0Type.Contains("Lambda")){
          gSystem->Exec(Form("cp %s-2.root %s-1.root", systname.Data(), systname.Data()));
          gSystem->Exec(Form("cp %s-2.root %s-3.root", systname.Data(), systname.Data()));
          gSystem->Exec(Form("cp %s-2.root %s-4.root", systname.Data(), systname.Data()));
        }
        if (lV0Type.Contains("K0Short")){
          gSystem->Exec(Form("cp %s-2.root %s-1.root", systname.Data(), systname.Data()));
          gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
        }

        lWhichSystVar = "NCrossedOverFindable";
      	v0->SetDefaultCuts(); //Do this before starting again
        cout<<endl<<"---> Performing Systematics Studies: "<<endl;
     	  //v0->SetCutLeastNumberOfCrossedRowsOverFindable(.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
     	  //v0->SetCutLeastNumberOfCrossedRowsOverFindable(5.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
     	  v0->SetCutLeastNumberOfCrossedRowsOverFindable(0.95); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
     	  //v0->SetCutLeastNumberOfCrossedRowsOverFindable(3.0); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();
     	  //
     	  systname = lSystFile + lWhichSystVar;
     	  gSystem->Exec(Form("cp %s-3.root %s-1.root", systname.Data(), systname.Data()));
     	  gSystem->Exec(Form("cp %s-3.root %s-2.root", systname.Data(), systname.Data()));
        gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));

         ////////////////////////////////////////////////////////////////////////////
        // GEANTFLUKACORRECTION ////////////////////////////////////////////////////
        lWhichSystVar = "GeantFlukaCorrection";
        cout<<endl<<"---> Performing Systematics Studies: Geant-Fluka correction"<<endl;
        if( lV0Type == "AntiLambda" ) {
            TFile* fGeantFluka = new TFile("antiPrCorrFunc.root", "READ");
            TF1* funcGeantFlukaCorr = (TF1*)fGeantFluka->Get("funcCorrAntiProtonFLUKA"); // used for systematics
            v0->SetGeant3FlukaCorr( funcGeantFlukaCorr ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); v0->DoAnalysis();
            //v0->SetGeantFlukaCorrection( funcGeantFlukaCorr ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); v0->DoAnalysis();
            //v0->SetGeantFlukaCorrection( funcGeantFlukaCorr ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); v0->DoAnalysis();
            //v0->SetGeantFlukaCorrection( funcGeantFlukaCorr ); v0->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); v0->DoAnalysis();
            //
            TString systname = lSystFile + lWhichSystVar;
            gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
            gSystem->Exec(Form("cp %s-1.root %s-3.root", systname.Data(), systname.Data()));
            gSystem->Exec(Form("cp %s-1.root %s-4.root", systname.Data(), systname.Data()));
        }
        if( lV0Type == "Lambda" ) {
        gSystem->Exec(Form("cp %s %s-1.root", Form("Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f_FD%s.root", lV0Type.Data(), lWhichMultEstimator.Data(),
    (lMultBoundLo), (lMultBoundHi), lWhichEffEnergyEstimator.Data(), (lEEBoundLo), (lEEBoundHi), fd.Data()), systname.Data()));
        gSystem->Exec(Form("cp %s %s-2.root", Form("Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f_FD%s.root", lV0Type.Data(), lWhichMultEstimator.Data(),
      (lMultBoundLo), (lMultBoundHi), lWhichEffEnergyEstimator.Data(), (lEEBoundLo), (lEEBoundHi), fd.Data()), systname.Data()));
        gSystem->Exec(Form("cp %s %s-2.root", Form("Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f_FD%s.root", lV0Type.Data(), lWhichMultEstimator.Data(),
    (lMultBoundLo), (lMultBoundHi), lWhichEffEnergyEstimator.Data(), (lEEBoundLo), (lEEBoundHi), fd.Data()), systname.Data()));
        gSystem->Exec(Form("cp %s %s-2.root", Form("Results-%s-13TeV-%s_%03.0f_%03.0f_%s_%03.0f_%03.0f_FD%s.root", lV0Type.Data(), lWhichMultEstimator.Data(),
      (lMultBoundLo), (lMultBoundHi), lWhichEffEnergyEstimator.Data(), (lEEBoundLo), (lEEBoundHi), fd.Data()), systname.Data()));
        }

  }
}
