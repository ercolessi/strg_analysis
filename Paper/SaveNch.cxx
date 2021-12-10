using namespace std;
#include <iostream>
#include <fstream>

#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLine.h>
#include <TFile.h>
#include <TTree.h>

void SaveNch(TString period = "17j"){

	TFile* file = new TFile(Form("data/LHC%sevents.root",period.Data()),"READ");

	cout<<"--------------- Open Real Data File --------------------"<<endl;
    TList* clist      = (TList*)file->Get("PWGLF_StrVsMult/cList");
    TTree* lTreeEvent = (TTree*)file->Get("PWGLF_StrVsMult/fTreeEvent");

    Float_t fCentrality_V0M = 0.;
    Float_t fCentrality_SPDClusters = 0.;
    Float_t fCentrality_ZDC = 0.;
    Int_t fRun = 0;
    Int_t fNTracksGlobal2015 = 0;
    Int_t fSPDTracklets = 0;
    // 
    lTreeEvent->SetBranchAddress("fRun",&fRun);
    lTreeEvent->SetBranchAddress("fCentrality_V0M", &fCentrality_V0M);
    lTreeEvent->SetBranchAddress("fCentrality_ZDC", &fCentrality_ZDC);
    lTreeEvent->SetBranchAddress("fCentrality_SPDClusters", &fCentrality_SPDClusters);
    lTreeEvent->SetBranchAddress("fNTracksGlobal", &fNTracksGlobal2015);
    lTreeEvent->SetBranchAddress("fSPDtracklets", &fSPDTracklets);
    
    TH3F* hspd_spdv0m = new TH3F("hspd_spdv0m","SPD-V0M vs SPD Tracklets; SPDTracklets; SPD [%]; V0M [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    TH3F* hntracks_spdv0m = new TH3F("hntracks_spdv0m","SPD-V0M vs N tracks global ; NTracksGlobal2015; SPD [%]; V0M [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    //
    TH3F* hspd_spdzdc = new TH3F("hspd_spdzdc","SPD-ZDC vs SPD Tracklets; SPDTracklets; SPD [%]; ZDC [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    TH3F* hntracks_spdzdc = new TH3F("hntracks_spdzdc","SPD-ZDC vs N tracks global ; NTracksGlobal2015; SPD [%]; ZDC [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    //
    TH3F* hspd_v0mzdc = new TH3F("hspd_v0mzdc","V0M-ZDC vs SPD Tracklets; SPDTracklets; V0M [%]; ZDC [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    TH3F* hntracks_v0mzdc = new TH3F("hntracks_v0mzdc","V0M-ZDC vs N tracks global ; NTracksGlobal2015; V0M [%]; ZDC [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    
 
    cout<<"\n\n--------------------------------------------------------"<<endl;
	cout<<" Will now loop over events, please wait..."<<endl;
	Long_t lNEvents = 0;
	for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries(); iEv++) {
	    lTreeEvent->GetEntry(iEv);
	    if( iEv % ( lTreeEvent->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;

        hspd_spdv0m->Fill(fSPDTracklets, fCentrality_SPDClusters, fCentrality_V0M);
        hntracks_spdv0m->Fill(fNTracksGlobal2015, fCentrality_SPDClusters, fCentrality_V0M);
        //
        hspd_spdzdc->Fill(fSPDTracklets, fCentrality_SPDClusters, fCentrality_ZDC);
        hntracks_spdzdc->Fill(fNTracksGlobal2015, fCentrality_SPDClusters, fCentrality_ZDC);
        //
        hspd_v0mzdc->Fill(fSPDTracklets, fCentrality_V0M, fCentrality_ZDC);
        hntracks_v0mzdc->Fill(fNTracksGlobal2015, fCentrality_V0M, fCentrality_ZDC);      
	}
	cout<<"--------------------------------------------------------\n\n"<<endl;

    TFile* Write = new TFile(Form("NchRawContainer_%s.root",period.Data()), "recreate");
    hspd_spdv0m->Write();
    hntracks_spdv0m->Write();
    hspd_spdzdc->Write();
    hntracks_spdzdc->Write();
    hspd_v0mzdc->Write();
    hntracks_v0mzdc->Write();
}

