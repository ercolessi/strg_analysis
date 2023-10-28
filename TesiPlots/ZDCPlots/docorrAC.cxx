void histobeauty(TH2D* h);

void docorrAC(){

    TFile* f = TFile::Open("../../PhDWork/Data/FullStat_pass2.root");
    cout<<"--------------- Open Data File --------------------\n"<<endl;
    TTree* lTreeEvent = (TTree *)f->Get("PWGLF_StrVsMult/fTreeEvent");

    Float_t fCentrality_V0M;
    Float_t fCentrality_SPDClusters;
    Float_t fZNApp;
    Float_t fZNCpp;
    Float_t fZPApp;
    Float_t fZPCpp;

    //Get from Tree
    lTreeEvent->SetBranchAddress("fCentrality_V0M", &fCentrality_V0M);
    lTreeEvent->SetBranchAddress("fCentrality_SPDClusters", &fCentrality_SPDClusters);
    lTreeEvent->SetBranchAddress("fZNApp", &fZNApp);
    lTreeEvent->SetBranchAddress("fZNCpp", &fZNCpp);
    lTreeEvent->SetBranchAddress("fZPApp", &fZPApp);
    lTreeEvent->SetBranchAddress("fZPCpp", &fZPCpp);

    TH2D* hZNAC = new TH2D("hZNAC", ";ZNA (a.u.);ZNC (a.u.)", 200, 0, 2000, 200, 0, 2000);
    TH2D *hZPAC = new TH2D("hZPAC", ";ZPA (a.u.);ZPC (a.u.)", 150, 0, 1500, 150, 0, 1500);
    histobeauty(hZNAC);
    histobeauty(hZPAC);
    cout<<" \nWill now loop over events, please wait...\n"<<endl;
    for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries()/10; iEv++) {

        lTreeEvent->GetEntry(iEv);
        if( iEv % ( lTreeEvent->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;

        hZNAC->Fill(fZNApp, fZNCpp);
        hZPAC->Fill(fZPApp, fZPCpp);
    }

    TLatex *tex = new TLatex(0.45, 0.89, "ALICE, pp #sqrt{s} = 13 TeV");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    TLatex *tex2 = new TLatex(0.45, 0.84, "This work");
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(2);

    TCanvas *c1 = new TCanvas("c1", "c1", 1500, 700);
    c1->Divide(2,1);
    c1->cd(1);
    c1->cd(1)->SetLogz();
    c1->cd(1)->SetBottomMargin(0.15);
    c1->cd(1)->SetRightMargin(0.15);
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetTopMargin(0.05);
    hZNAC->Draw("colz");
    tex->Draw();
    tex2->Draw();
    c1->cd(2)->SetLogz();
    c1->cd(2)->SetBottomMargin(0.15);
    c1->cd(2)->SetRightMargin(0.15);
    c1->cd(2)->SetLeftMargin(0.15);
    c1->cd(2)->SetTopMargin(0.05);
    c1->cd(2);
    hZPAC->Draw("colz");
    tex->Draw();
    tex2->Draw();

    c1->SaveAs("ZNACvsZPAC.pdf");
}

void histobeauty(TH2D* h){
    h->SetStats(0);
    h->SetTitle("");
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.035);
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleOffset(1.6);
    h->GetYaxis()->SetLabelSize(0.035);
    h->SetMaximum(100000);
}