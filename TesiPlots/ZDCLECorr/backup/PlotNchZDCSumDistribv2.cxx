void dostudy(TString f = "", TString period = "18i", TString sel = "SPD", float sel1 = 0., float sel2 = 100.,  float var1 = 0., float var2 = 100., Int_t colors=0);

void PlotNchZDCSumDistribv2(){

    Int_t color[] = {kBlack, kRed+1, kGreen+2, kBlue+1};

    TString f = "ZDCDistributionCheck.root";

    dostudy(f,"18i","SPD",0.,100.,0,100,color[0]);
    dostudy(f,"18i","SPD",10.,20.,0,5,color[1]);
    dostudy(f,"18i","SPD",10.,20.,20,30,color[2]);
    dostudy(f,"18i","SPD",10.,20.,50,100,color[3]);
    dostudy(f,"18i","SPD",40.,50.,0,20,color[1]);
    dostudy(f,"18i","SPD",40.,50.,40,50,color[2]);
    dostudy(f,"18i","SPD",40.,50.,70,100,color[3]);

}

void dostudy(TString f = "", TString period = "18i", TString sel = "SPD", float sel1 = 0., float sel2 = 100., float var1 = 0., float var2 = 100., Int_t colors){

    TFile* file = new TFile(Form("LHC%s.root",period.Data()),"READ");

	cout<<"--------------- Open Real Data File --------------------"<<endl;
    TTree* lTreeEvent = (TTree*)file->Get("PWGLF_StrVsMult/fTreeEvent");

    Float_t fCentrality_V0M = 0.;
    Float_t fCentrality_SPDClusters = 0.;
    Float_t fCentrality_ZDCFired = 0.;
    Int_t fRun = 0;
    Int_t fSPDTracklets = 0;
    Float_t fZPApp = 0.;
    Float_t fZPCpp = 0.;
    Float_t fZNApp = 0.;
    Float_t fZNCpp = 0.;

    // 
    lTreeEvent->SetBranchAddress("fRun",&fRun);
    lTreeEvent->SetBranchAddress("fCentrality_V0M", &fCentrality_V0M);
    lTreeEvent->SetBranchAddress("fCentrality_SPDClusters", &fCentrality_SPDClusters);
    lTreeEvent->SetBranchAddress("fCentrality_ZDCFired", &fCentrality_ZDCFired);
    lTreeEvent->SetBranchAddress("fSPDtracklets", &fSPDTracklets);
    lTreeEvent->SetBranchAddress("fZPApp", &fZPApp);
    lTreeEvent->SetBranchAddress("fZPCpp", &fZPCpp);
    lTreeEvent->SetBranchAddress("fZNApp", &fZNApp);
    lTreeEvent->SetBranchAddress("fZNCpp", &fZNCpp);

    TH1F* hnch  = new TH1F(Form("hnch"), "; SPD tracklets; Counts (normalised)", 35,0.,35.); 
    TH1F* hnchmb = new TH1F(Form("hnch"), "; SPD tracklets; Counts (normalised)", 35,0.,35.);
    TH1F* hzdcsum = new TH1F(Form("hzdcsum"), "; ZDC Energy Sum (a.u.); Counts (normalised)", 200,0.,2000.);
    TH1F *hzdcsummb = new TH1F(Form("hzdcsum"), "; ZDC Energy Sum (a.u.); Counts (normalised)", 200,0.,2000.);

    cout<<"\n\n--------------------------------------------------------"<<endl;
    cout<<" Will now loop over events, please wait..."<<endl;
    Long_t lNEvents = 0;
    for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries(); iEv++) {

        lTreeEvent->GetEntry(iEv);
        if( iEv % ( lTreeEvent->GetEntries() / 100 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;

        Float_t lZDCSum = fZNApp + fZNCpp + fZPApp + fZPCpp;

        if(fCentrality_ZDCFired>100) continue;

        if(sel.Contains("SPD")) {
            if(fCentrality_SPDClusters >= sel1 && fCentrality_SPDClusters <sel2 && fCentrality_V0M >= var1 && fCentrality_V0M < var2){
                hnch->Fill(fSPDTracklets);
                hzdcsum->Fill(lZDCSum);
            }
        }
        if(sel.Contains("V0M")) {
                if(fCentrality_V0M >= sel1 && fCentrality_V0M <sel2 && fCentrality_SPDClusters >= var1 && fCentrality_SPDClusters < var2){
                hnch->Fill(fSPDTracklets);
                hzdcsum->Fill(lZDCSum);
            }             
        }
         
        hnchmb->Fill(fSPDTracklets);        
        hzdcsummb->Fill(lZDCSum);
    }

    file->Close();

    hnch->SetLineColor(colors);
    hnch->SetMarkerColor(colors);
    hnch->SetMarkerStyle(20);
    hnch->SetLineWidth(2);
    hzdcsum->SetLineColor(colors);
    hzdcsum->SetMarkerColor(colors);
    hzdcsum->SetMarkerStyle(20);     
    hzdcsum->SetLineWidth(2);     

    TFile* f1 = TFile::Open(f,"update");
    if(sel.Contains("SPD"))hnch->SetName(Form("hnch_%s_%.0f_%.0f_%s_%.0f_%.0f",sel.Data(),sel1,sel2,"V0M",var1,var2));
    if(sel.Contains("V0M"))hnch->SetName(Form("hnch_%s_%.0f_%.0f_%s_%.0f_%.0f",sel.Data(),sel1,sel2,"SPDClusters",var1,var2));
    hnch->Write();
    if(sel.Contains("SPD"))hzdcsum->SetName(Form("hzdcsum_%s_%.0f_%.0f_%s_%.0f_%.0f",sel.Data(),sel1,sel2,"V0M",var1,var2));
    if(sel.Contains("V0M"))hzdcsum->SetName(Form("hzdcsum_%s_%.0f_%.0f_%s_%.0f_%.0f",sel.Data(),sel1,sel2,"SPDClusters",var1,var2));
    hzdcsum->Write();
    f1->Close();

}
