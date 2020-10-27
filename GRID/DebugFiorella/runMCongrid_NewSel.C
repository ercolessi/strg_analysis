
////////////////////////////////////////////////////////////////////////////////////
//                                                                                //
// Macro to run strangeness MC analysis vs multiplicity & effective energy on grid//
// author: Francesca Ercolessi                                                    //
//         francesca.ercolessi@cern.ch                                            //
//                                                                                //
//                                                                                //
// To Run:                                                                        //
// testmode) gridtest=KTRUE                                                       //
// gridmode) gridtest=KFALSE                                                      //
// 1)run with runmode=full and setMergeViaJDL(kTRUE) to run on grid               //
// 2)run with runmode=terminate and SetMergeViaJDL(kTRUE) to merge                //
// 3)run with runmode=terminate and SetMergeViaJDL(kFALSE) to download in local   //
//    the merged root file                                                        //
//                                                                                //
// in aliroot (ROOT5)                                                             //
// >> .L runMCongrid.C                                                            //
// >> runMCongrid(Period)                                                         //
//                         example for LHC12b  Period=122                         //
//                               (year (12) + b second letter of the alphabet)    // 
//                                                                                //
//                                                                                //              
////////////////////////////////////////////////////////////////////////////////////

#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskLeading.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AddTaskPhysicsSelection.h"
#include "AliAnalysisTask.h"
#include "AliPhysicsSelectionTask.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"
#include "AliAnalysisTaskStrangenessVsMultiplicityEEMCRun222.h"
#endif

void runMCongrid_NewSel(Int_t Period = 156)
  
{
  // if you need to access data remotely from the GRID
  Bool_t grid = 1;  
  // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
  Bool_t gridTest = kFALSE;
  // if the data are MC
  Bool_t isMCdata=kFALSE;
  Bool_t isMC=kTRUE; // if you have the kinematics information
  if(isMC) isMCdata=kTRUE;
 
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
   
  // include the path you need to compile and the library you need
  gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");   
  //Tell root where to look for headers
  #if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
  #else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
  #endif

  // Create the analysis manager and ESD handler
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  // add interface to MC if requested
  if(isMC){
    AliMCEventHandler* mcH = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcH);
  }  

  //PhysicsSelection Configuration
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* ps = AddTaskPhysicsSelection(kTRUE);
  
  //MultSelection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask* ms = AddTaskMultSelection();
  ms->SetAddInfo(kTRUE);
  
  //PID Response
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse* pid = AddTaskPIDResponse(kTRUE);
  
  //Weak Decay Vertexer
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexer.C");
  AliAnalysisTaskWeakDecayVertexer * wv = AddTaskWeakDecayVertexer();
  //Fiorella's configuration
  wv->SetPreselectDedx(kFALSE);
  wv->SetPreselectDedxLambda(kFALSE);
  wv->SetExtraCleanup(kFALSE);
  wv->SetUseExtraEvSels(kFALSE);
  wv->SetRunV0Vertexer(kTRUE);
  wv->SetRunCascadeVertexer(kTRUE);
  wv->SetDoV0Refit(kTRUE);
  wv->SetDoCascadeRefit(kTRUE);
  wv->SetDoImprovedDCAV0DauPropagation(kTRUE);
  wv->SetXYCase1Preoptimization(kTRUE);
  wv->SetXYCase2Preoptimization(kTRUE);
  wv->SetDoImprovedDCACascDauPropagation(kTRUE);
  wv->SetDoPureGeometricMinimization(kTRUE);
  wv->SetUseImprovedFinding();
  //V0
  wv->SetV0VertexerDCAFirstToPV(0.02);
  wv->SetV0VertexerDCASecondtoPV(0.02);
  wv->SetV0VertexerDCAV0Daughters(2.0) ;
  wv->SetV0VertexerCosinePA(0.0) ;
  wv->SetV0VertexerMinRadius(0.2) ;
  wv->SetV0VertexerMaxRadius(200);
  //Cascade
  wv->SetCascVertexerMinV0ImpactParameter(0.05) ; 
  wv->SetCascVertexerV0MassWindow(0.01) ;
  wv->SetCascVertexerDCABachToPV(0.03) ;
  wv->SetCascVertexerDCACascadeDaughters(2.0) ;
  wv->SetCascVertexerCascadeMinRadius(.3) ;
  wv->SetCascVertexerCascadeCosinePA(.95) ;

  //Strangness Task
  gROOT->LoadMacro("AliAnalysisTaskStrangenessVsMultiplicityEEMCRun22.cxx++g");
  gROOT->LoadMacro("AddTaskStrangenessVsMultiplicityEEMCRun22.C");
  AliAnalysisTaskStrangenessVsMultiplicityEEMCRun22* st = AddTaskStrangenessVsMultiplicityEEMCRun22(kTRUE,kTRUE,kTRUE,"ABC");
  st->SetDownScaleV0(kFALSE, 1.0000);
  st->SetDownScaleCascade(kFALSE, 1.0000);
  st->SetRunVertexers(kFALSE);  
  st->SetSelectedTriggerClass(AliVEvent::kINT7);
  st->SetApplySPDClsVsTrackletsCut(kTRUE);  
  st->SetUseExtraEvSels(kFALSE);
   

  if(!mgr->InitAnalysis()) return;
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();

  if (grid)
    {
      AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
      alienHandler->SetOverwriteMode();

      alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        
      // Set versions of used packages
      alienHandler->SetAPIVersion("V1.1x");
      //alienHandler->SetROOTVersion("v5-26-00b-6");
      // alienHandler->SetAliROOTVersion("v4-19-21-AN");
      //Please keep this version updated
      alienHandler->SetAliPhysicsVersion("vAN-20200902-1");   
      alienHandler->SetAnalysisMacro("AnalysisStrangenessMC.C"); 

      // number of files per subjob
      //alienHandler->SetSplitMaxInputFileNumber(10);
      //set Executable
      alienHandler->SetExecutable("MCStrangenessESD.sh");
      //specify how many seconds your job may take
      alienHandler->SetTTL(36000);
      //set jdl name
      alienHandler->SetJDLName("MCStrangenessESD.jdl");
      
      alienHandler->SetOutputToRunNo(kTRUE);
      alienHandler->SetKeepLogs(kTRUE);
   
      alienHandler->SetTerminateFiles("event_stat.root");  //to have the output file of the Physics Selection class
      alienHandler->SetInputFormat("xml-single");
      alienHandler->SetPrice(1);      
      // Optionally modify split mode (default 'se')    
      alienHandler->SetSplitMode("se");
       
      // make sure your source files get copied to grid
      alienHandler->SetAdditionalLibs("AliAnalysisTaskStrangenessVsMultiplicityEEMCRun22.cxx AliAnalysisTaskStrangenessVsMultiplicityEEMCRun22.h");
      alienHandler->SetAnalysisSource("AliAnalysisTaskStrangenessVsMultiplicityEEMCRun22.cxx");
      
      //Declare input data to be processed.
      //Method 1: Create automatically XML collections using alien 'find' command.
          
      //LHC15f---------------------------------------------
      if (Period==156){
	
      	alienHandler->SetGridDataDir("/alice/sim/2015/LHC15g3b1/");
      	//alienHandler->SetRunPrefix("000");
      	alienHandler->SetDataPattern("/*/AliESDs.root");
      	// runnumber
      	Int_t runList[51] = {226500, 226495, 226483, 226476, 226472, 226468, 226466, 226452, 226445, 226444, 226225, 226220, 226170, 226062, 225768, 225766, 225763, 225762, 225757, 225753, 225719, 225717, 225716, 225710, 225709, 225708, 225707, 225705, 225587, 225586, 225579, 225578, 225576, 225322, 225315, 225314, 225313, 225310, 225309, 225307, 225305, 225106, 225052, 225051, 225050, 225043, 225041, 225037, 225035, 225031, 225026};
      	for (Int_t i = 0;i <5; i++) alienHandler->AddRunNumber(runList[i]);
      	
      	alienHandler->SetGridWorkingDir("NewEvSel/MC/LHC15f");
      }
              
    	//LHC17j
    	if (Period==17){

        alienHandler->SetGridDataDir("/alice/sim/2020/LHC20i2c/");
        //alienHandler->SetRunPrefix("000");
        alienHandler->SetDataPattern("/*/AliESDs.root");
        // runnumber
        Int_t runList[10] = {274593, 274594, 274595, 274596, 274601, 274653, 274657, 274667, 274669, 274671};
        for (Int_t i = 0;i <10; i++) alienHandler->AddRunNumber(runList[i]);

        alienHandler->SetGridWorkingDir("test/LHC17j");
      }


      //LHC18i---------------------------------------------
      if (Period==18){
      	// select the input data 
      	alienHandler->SetGridDataDir("/alice/sim/2020/LHC20i2b/");
      	//alienHandler->SetRunPrefix("000");
      	alienHandler->SetDataPattern("/*/AliESDs.root");
      	// runnumber
      	Int_t runList[10] = {288861, 288862, 288863, 288864, 288868, 288897, 288902, 288903, 288908, 288909};
      	for (Int_t i= 1;i <10; i++) alienHandler->AddRunNumber(runList[i]);
      	 
      	alienHandler->SetGridWorkingDir("MC/LHC18i");
      }                 
          
      // define the output folder
      alienHandler->SetGridOutputDir("OutputDir");
      alienHandler->SetDefaultOutputs();
      
      //number of times to merge (if a lot of data need a higher number)
      alienHandler->SetMaxMergeStages(5);
      alienHandler->SetSplitMaxInputFileNumber(100);
      // we can specify that we want to, later on, use Grid to also merge
      // our output. to enable this, we will set 'SetMergeViaJDL' to kTRUE
      alienHandler->SetMergeViaJDL(kTRUE);
      //When all your merging jobs have finished running,
      //there is one step to be taken still, which is downloading the output
      //of the merging jobs. This is done by changing SetMergeViaJDL(kFALSE)
      //and running one last time
       
      // Connect plug-in to the analysis manager
      mgr->SetGridHandler(alienHandler);
      
      if(gridTest) {
      	// speficy on how many files you want to run
      	alienHandler->SetNtestFiles(2);
      	// and launch the analysis
      	// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
      	alienHandler->SetRunMode("test");
      	mgr->StartAnalysis("grid");
    	
      } else
    	{
    	  //full grid      
    	  alienHandler->SetRunMode("terminate");
    	  mgr->StartAnalysis("grid");
    	}
      
    }

}

