#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AliAnalysisTask.h"
#include "AliPhysicsSelectionTask.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"
#include "AliAnalysisTaskStrangenessVsMultVsEffEnergyAODRun2bis.h"
#include "AliPPVsMultUtils.h"
#include "AliAODTracklets.h"

#endif

void runongrid_default(Int_t Period = 15)

{
  // if you need to access data remotely from the GRID
  Bool_t grid = 1;
  // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
  Bool_t gridTest = kFALSE;

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

  // Create the analysis manager and AOD handler
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  AliAODInputHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);

  //PhysicsSelection Configuration
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* ps =  reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(kFALSE)"));

  //MultSelection
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask* ms =  reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
  ms->SetAddInfo(kTRUE);

  //PID Response
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse* pid = reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(kFALSE)"));

  //Strangness Task
  gROOT->LoadMacro("AliAnalysisTaskStrangenessVsMultVsEffEnergyAODRun2bis.cxx++g");
  AliAnalysisTaskStrangenessVsMultVsEffEnergyAODRun2bis *st = reinterpret_cast<AliAnalysisTaskStrangenessVsMultVsEffEnergyAODRun2bis *>(gInterpreter->ExecuteMacro("AddTaskStrangenessVsMultVsEffEnergyAODRun2bis.C(1,1,1,\"ABC\")"));
  st->SetSaveV0s(kTRUE);
  st->SetSaveCascades(kTRUE);
  st->SetDownScaleEvent(kFALSE, 1.0);
  st->SetDownScaleV0(kFALSE, 1.0);
  st->SetDownScaleCascade(kFALSE, 1.0);
  st->SetRunVertexers(kFALSE); //can remove
  st->SetSelectedTriggerClass(AliVEvent::kINT7);
  st->SetApplySPDClsVsTrackletsCut(kTRUE);
  st->SetUseExtraEvSels(kFALSE);
  st->SetUseLightVertexers(kFALSE); //can also remove

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
      //alienHandler->SetAliROOTVersion("v4-19-21-AN");
      //Please keep this version updated
      alienHandler->SetAliPhysicsVersion("vAN-20230307_O2-1");
      alienHandler->SetAnalysisMacro("AnalysisLeading.C");

      // number of files per subjob
      //alienHandler->SetSplitMaxInputFileNumber(5);
      //set Executable
      alienHandler->SetExecutable("runLeadingAOD.sh");
      //specify how many seconds your job may take
      alienHandler->SetTTL(36000);
      //set jdl name
      alienHandler->SetJDLName("runLeadingAOD.jdl");

      alienHandler->SetOutputToRunNo(kTRUE);
      alienHandler->SetKeepLogs(kTRUE);

      alienHandler->SetTerminateFiles("event_stat.root");  //to have the output file of the Physics Selection class
      alienHandler->SetInputFormat("xml-single");
      alienHandler->SetPrice(1);
      // Optionally modify split mode (default 'se')
      alienHandler->SetSplitMode("se");

      // make sure your source files get copied to grid
      alienHandler->SetAdditionalLibs("AliAnalysisTaskStrangenessVsMultVsEffEnergyAODRun2bis.cxx AliAnalysisTaskStrangenessVsMultVsEffEnergyAODRun2bis.h");
      alienHandler->SetAnalysisSource("AliAnalysisTaskStrangenessVsMultVsEffEnergyAODRun2bis.cxx ");

      //Declare input data to be processed.
      //Method 1: Create automatically XML collections using alien 'find' command.

      //LHC17j---------------------------------------------
     	// select the input data
      if (Period==17){
        alienHandler->SetGridDataDir("/alice/data/2017/LHC17j/");
        alienHandler->SetRunPrefix("000");
        alienHandler->SetDataPattern("/pass1/AOD234/*/AliAOD.root");
        // runnumber
        Int_t runList[9] = {274593,  274595, 274596, 274601, 274653, 274657, 274667, 274669, 274671};
        for (Int_t i= 0;i<9; i++) alienHandler->AddRunNumber(runList[i]);
        alienHandler->SetGridWorkingDir("Calibration/LHC17j");
      }
      //LHC15f---------------------------------------------
      if (Period==15){
    	// select the input data
    	alienHandler->SetGridDataDir("/alice/data/2015/LHC15f/");
    	alienHandler->SetRunPrefix("000");
    	alienHandler->SetDataPattern("/pass2/AOD234/*/AliAOD.root");
    	// runnumber
    	Int_t runList[43] = {
        226500, 226495, 226483, 226472, 226468, 226466, 226452,
        226445, 226444, 226225, 226220, 226062, 225719, 225717,
        225716, 225710, 225709, 225708, 225707, 225705, 225587,
        225586, 225579, 225578, 225576, 225322, 225315, 225314,
        225313, 225310, 225309, 225307, 225305, 225052, 225051,
        225050, 225043, 225041, 225037, 225035, 225031, 225026,
        225106};
        //226500, 226495, 226483, 226476, 226472, 226468, 226466, 226452, 226445, 226444, 226225, 226220, 226170, 226062, 225768, 225766, 225763, 225762, 225757, 225753, 225719, 225717, 225716, 225710, 225709, 225708, 225707, 225705, 225587, 225586, 225579, 225578, 225576, 225322, 225315, 225314, 225313, 225310, 225309, 225307, 225305, 225106, 225052, 225051, 225050, 225043, 225041, 225037, 225035, 225031, 225026};
    	for (Int_t i= 0;i <1; i++) alienHandler->AddRunNumber(runList[i]);
    	  alienHandler->SetGridWorkingDir("PlotsTesi");
      }


      //LHC18i---------------------------------------------
      if (Period==18){
        // select the input data
        alienHandler->SetGridDataDir("/alice/data/2018/LHC18i/");
        alienHandler->SetRunPrefix("000");
        alienHandler->SetDataPattern("pass1/AOD264/AOD/*/AliAOD.root");
        // runnumber
        Int_t runList[9] = {288861, 288862, 288863, 288868, 288902, 288903, 288908, 288864, 288909};
        for (Int_t i= 0;i <9; i++) alienHandler->AddRunNumber(runList[i]);
          alienHandler->SetGridWorkingDir("Calibration/LHC18i");
      }


      // define the output folder
      alienHandler->SetGridOutputDir("OutputDir");
      //alienHandler->SetDefaultOutputs();

      //number of times to merge (if a lot of data need a higher number)
      alienHandler->SetMaxMergeStages(10);
      alienHandler->SetSplitMaxInputFileNumber(3);
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
      	alienHandler->SetNtestFiles(1);
      	// and launch the analysis
      	// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
      	alienHandler->SetRunMode("test");
      	mgr->StartAnalysis("grid");

      } else
    	{
    	  //full grid
    	  alienHandler->SetRunMode("full");
    	  mgr->StartAnalysis("grid");
    	}

    }

}

