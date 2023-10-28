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

// ZDC calibrations per period
float refZDCN = 140; 
float refZDCP = 35; 

float cutZDCN = -100;
float cutZDCP = -100;

void SaveNch(TString period = "18i"){
    TFile* file = new TFile(Form("LHC%s.root",period.Data()),"READ");

    TFile *fcal = TFile::Open("calZDC.root");
    TProfile *calna, *calnc, *calpa, *calpc;
    if(fcal){
      calna = (TProfile *) fcal->Get("zdcna");
      calna->SetName("calna");
      calnc = (TProfile *) fcal->Get("zdcnc");
      calnc->SetName("calnc");
      calpa = (TProfile *) fcal->Get("zdcpa");
      calpa->SetName("calpa");
      calpc = (TProfile *) fcal->Get("zdcpc");
      calpc->SetName("calpc");
    }

    TFile *fped = TFile::Open("pedZDC.root");
    TProfile *pedna, *pednc, *pedpa, *pedpc;
    if(fped){
      pedna = (TProfile *) fped->Get("zdcnaPed");
      pedna->SetName("pedna");
      pednc = (TProfile *) fped->Get("zdcncPed");
      pednc->SetName("pednc");
      pedpa = (TProfile *) fped->Get("zdcpaPed");
      pedpa->SetName("pedpa");
      pedpc = (TProfile *) fped->Get("zdcpcPed");
      pedpc->SetName("pedpc");
    }

    cout<<"--------------- Open Real Data File --------------------"<<endl;
    TList* clist      = (TList*)file->Get("PWGLF_StrVsMult/cList");
    TTree* lTreeEvent = (TTree*)file->Get("PWGLF_StrVsMult/fTreeEvent");

    Float_t fCentrality_V0M = 0.;
    Float_t fCentrality_SPDClusters = 0.;
    Float_t fCentrality_ZDC = 0.;
    Int_t fRun = 0;
    Int_t fNTracksGlobal2015 = 0;
    Int_t fSPDTracklets = 0;
    Float_t fZPApp = 0.;
    Float_t fZPCpp = 0.;
    Float_t fZNApp = 0.;
    Float_t fZNCpp = 0.;
    Float_t fCentrality_ZDCFired = 0.;

    // 
    lTreeEvent->SetBranchAddress("fRun",&fRun);
    lTreeEvent->SetBranchAddress("fCentrality_V0M", &fCentrality_V0M);
    lTreeEvent->SetBranchAddress("fCentrality_ZDC", &fCentrality_ZDC);
    lTreeEvent->SetBranchAddress("fCentrality_SPDClusters", &fCentrality_SPDClusters);
    lTreeEvent->SetBranchAddress("fNTracksGlobal", &fNTracksGlobal2015);
    lTreeEvent->SetBranchAddress("fSPDtracklets", &fSPDTracklets);
    lTreeEvent->SetBranchAddress("fZPApp", &fZPApp);
    lTreeEvent->SetBranchAddress("fZPCpp", &fZPCpp);
    lTreeEvent->SetBranchAddress("fZNApp", &fZNApp);
    lTreeEvent->SetBranchAddress("fZNCpp", &fZNCpp);
    lTreeEvent->SetBranchAddress("fCentrality_ZDCFired", &fCentrality_ZDCFired);
    
    TH3F* hspd_spdv0m = new TH3F("hspd_spdv0m","SPD-V0M vs SPD Tracklets; SPDTracklets; SPD [%]; V0M [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    TH3F* hzdcsum_spdv0m = new TH3F("hzdcsum_spdv0m","SPD-V0M vs ZDC Sum (a.u.); ZDC Sum (a.u.); SPD [%]; V0M [%]",4100,-100.,4000.,100,0.,100.,100,0.,100.);
    TH3F* hznsum_spdv0m = new TH3F("hznsum_spdv0m","SPD-V0M vs ZN (a.u.); ZN (a.u.); SPD [%]; V0M [%]",4100,-100.,4000.,100,0.,100.,100,0.,100.);
    TH3F* hzpsum_spdv0m = new TH3F("hzpsum_spdv0m","SPD-V0M vs ZP Sum (a.u.); ZP (a.u.); SPD [%]; V0M [%]",4100,-100.,4000.,100,0.,100.,100,0.,100.);
    
    TH3F* hntracks_spdv0m = new TH3F("hntracks_spdv0m","SPD-V0M vs N tracks global ; NTracksGlobal2015; SPD [%]; V0M [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    //
    TH3F* hspd_spdzdc = new TH3F("hspd_spdzdc","SPD-ZDC vs SPD Tracklets; SPDTracklets; SPD [%]; ZDC [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    TH3F* hntracks_spdzdc = new TH3F("hntracks_spdzdc","SPD-ZDC vs N tracks global ; NTracksGlobal2015; SPD [%]; ZDC [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    //
    TH3F* hspd_v0mzdc = new TH3F("hspd_v0mzdc","V0M-ZDC vs SPD Tracklets; SPDTracklets; V0M [%]; ZDC [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    TH3F* hntracks_v0mzdc = new TH3F("hntracks_v0mzdc","V0M-ZDC vs N tracks global ; NTracksGlobal2015; V0M [%]; ZDC [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    //
    TH3F* hspdclspdclv0m = new TH3F("hspdclspdclv0m","SPDcl vs SPDcl vs V0M; SPD [%]; SPD [%]; V0M [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    TH3F* hspdclspdclzdc = new TH3F("hspdclspdclzdc","SPDcl vs SPDcl vs ZDC; SPD [%]; SPD [%]; ZDC [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    TH3F* hv0mv0mspdcl = new TH3F("hv0mv0mspdcl","V0M vs V0M vs SPDcl; V0M [%]; V0M [%]; SPD [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    TH3F* hzdczdcspdcl = new TH3F("hzdczdcspdcl","ZDC vs ZDC vs SPDcl; ZDC [%]; ZDC [%]; SPD [%]",100,0.,100.,100,0.,100.,100,0.,100.);
    
 
    TProfile *zdcnc = new TProfile("zdcnc",";run;<ZDCNC> (a.u.)",100000,200000,300000);
    TProfile *zdcna = new TProfile("zdcna",";run;<ZDCNC> (a.u.)",100000,200000,300000);
    TProfile *zdcpc = new TProfile("zdcpc",";run;<ZDCNC> (a.u.)",100000,200000,300000);
    TProfile *zdcpa = new TProfile("zdcpa",";run;<ZDCNC> (a.u.)",100000,200000,300000);

    TProfile *zdcncPed = new TProfile("zdcncPed",";run;<ZDCNC> (a.u.)",100000,200000,300000);
    TProfile *zdcnaPed = new TProfile("zdcnaPed",";run;<ZDCNC> (a.u.)",100000,200000,300000);
    TProfile *zdcpcPed = new TProfile("zdcpcPed",";run;<ZDCNC> (a.u.)",100000,200000,300000);
    TProfile *zdcpaPed = new TProfile("zdcpaPed",";run;<ZDCNC> (a.u.)",100000,200000,300000);

    cout<<"\n\n--------------------------------------------------------"<<endl;
	cout<<" Will now loop over events, please wait..."<<endl;
	Long_t lNEvents = 0;
	for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries(); iEv++) {
	    lTreeEvent->GetEntry(iEv);
	    if( iEv % ( lTreeEvent->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;


         if(fRun == 225052 || fRun == 225051 || fRun == 225050 || fRun == 225043 || fRun == 225041 || fRun ==  225037 || fRun == 225035 || fRun == 225031 || fRun == 225026){
//printf("run %d\n",fRun);
           continue;
         } else {
//         continue;
         }

//printf("good\n");
        // event selection
        if(fCentrality_ZDCFired > 100) continue;
//printf("goodx2\n");

//        fZNApp=0;
//        fZPApp=0;

        if(fped){
          fZNApp -= pedna->GetBinContent(fRun - 199999);
          fZNCpp -= pednc->GetBinContent(fRun - 199999);
          fZPApp -= pedpa->GetBinContent(fRun - 199999);
          fZPCpp -= pedpc->GetBinContent(fRun - 199999);
        }

        if(TMath::Abs(fZNCpp) < 20) zdcncPed->Fill(fRun, fZNCpp);
        if(TMath::Abs(fZNApp) < 20) zdcnaPed->Fill(fRun, fZNApp);
        if(TMath::Abs(fZPCpp) < 20) zdcpcPed->Fill(fRun, fZPCpp);
        if(TMath::Abs(fZPApp) < 20) zdcpaPed->Fill(fRun, fZPApp);


        // apply calibrations
        if(fcal){
          fZNApp /= calna->GetBinContent(fRun - 199999);
          fZNCpp /= calnc->GetBinContent(fRun - 199999);
          fZPApp /= calpa->GetBinContent(fRun - 199999);
          fZPCpp /= calpc->GetBinContent(fRun - 199999);
        }

        if(fZNCpp > cutZDCN) zdcnc->Fill(fRun, fZNCpp/refZDCN);
        if(fZNApp > cutZDCN) zdcna->Fill(fRun, fZNApp/refZDCN);
        if(fZPCpp > cutZDCP) zdcpc->Fill(fRun, fZPCpp/refZDCP);
        if(fZPApp > cutZDCP) zdcpa->Fill(fRun, fZPApp/refZDCP);

        Float_t lZDCSum = fZNApp + fZNCpp + fZPApp + fZPCpp;

        hspd_spdv0m->Fill(fSPDTracklets, fCentrality_SPDClusters, fCentrality_V0M);
        hntracks_spdv0m->Fill(fNTracksGlobal2015, fCentrality_SPDClusters, fCentrality_V0M);
        hzdcsum_spdv0m->Fill(lZDCSum, fCentrality_SPDClusters, fCentrality_V0M);
        hznsum_spdv0m->Fill(fZNApp + fZNCpp, fCentrality_SPDClusters, fCentrality_V0M);
        hzpsum_spdv0m->Fill(fZPApp + fZPCpp, fCentrality_SPDClusters, fCentrality_V0M);
        //
        hspd_spdzdc->Fill(fSPDTracklets, fCentrality_SPDClusters, fCentrality_ZDC);
        hntracks_spdzdc->Fill(fNTracksGlobal2015, fCentrality_SPDClusters, fCentrality_ZDC);
        //
        hspd_v0mzdc->Fill(fSPDTracklets, fCentrality_V0M, fCentrality_ZDC);
        hntracks_v0mzdc->Fill(fNTracksGlobal2015, fCentrality_V0M, fCentrality_ZDC);      
        //
        hspdclspdclv0m->Fill(fCentrality_SPDClusters, fCentrality_SPDClusters, fCentrality_V0M);
        hspdclspdclzdc->Fill(fCentrality_SPDClusters, fCentrality_SPDClusters, fCentrality_ZDC);
        //
        hv0mv0mspdcl->Fill(fCentrality_V0M, fCentrality_V0M, fCentrality_SPDClusters);
        hzdczdcspdcl->Fill(fCentrality_ZDC, fCentrality_ZDC, fCentrality_SPDClusters);
	}
	cout<<"--------------------------------------------------------\n\n"<<endl;

    TFile* Write = new TFile(Form("NchRawContainer_%s.root",period.Data()), "recreate");
    hspd_spdv0m->Write();
    hntracks_spdv0m->Write();
    hzdcsum_spdv0m->Write();
    hspd_spdzdc->Write();
    hntracks_spdzdc->Write();
    hspd_v0mzdc->Write();
    hntracks_v0mzdc->Write();
    hspdclspdclv0m->Write();
    hspdclspdclzdc->Write();
    hv0mv0mspdcl->Write();
    hzdczdcspdcl->Write();
    hznsum_spdv0m->Write();
    hzpsum_spdv0m->Write();
    zdcnc->Write();
    zdcna->Write();
    zdcpc->Write();
    zdcpa->Write();
    zdcncPed->Write();
    zdcnaPed->Write();
    zdcpcPed->Write();
    zdcpaPed->Write();
}
