void histobeauty(TH2D* h);

void plot()
{

    TFile *f = TFile::Open("../../PhDWork/Data/FullStat_pass2.root");
    cout << "--------------- Open Data File --------------------\n"
         << endl;
    TTree *lTreeEvent = (TTree *)f->Get("PWGLF_StrVsMult/fTreeEvent");

    Float_t fCentrality_V0M;
    Float_t fCentrality_SPDClusters;
    Float_t fZNApp;
    Float_t fZNCpp;
    Float_t fZPApp;
    Float_t fZPCpp;
    Float_t fAmplitudeV0A;
    Float_t fAmplitudeV0C;

    // Get from Tree
    lTreeEvent->SetBranchAddress("fCentrality_V0M", &fCentrality_V0M);
    lTreeEvent->SetBranchAddress("fCentrality_SPDClusters", &fCentrality_SPDClusters);
    lTreeEvent->SetBranchAddress("fZNApp", &fZNApp);
    lTreeEvent->SetBranchAddress("fZNCpp", &fZNCpp);
    lTreeEvent->SetBranchAddress("fZPApp", &fZPApp);
    lTreeEvent->SetBranchAddress("fZPCpp", &fZPCpp);
    lTreeEvent->SetBranchAddress("fAmplitudeV0A", &fAmplitudeV0A);
    lTreeEvent->SetBranchAddress("fAmplitudeV0C", &fAmplitudeV0C);

    TH2D *hV0MSPD = new TH2D("hV0MSPD", ";VZEROM percentile; SPDClusters percentile", 100, 0, 100, 100, 0, 100);
    histobeauty(hV0MSPD);
    cout << " \nWill now loop over events, please wait...\n"
         << endl;
    for (Long_t iEv = 0; iEv < lTreeEvent->GetEntries()/10 ; iEv++)
    {

        lTreeEvent->GetEntry(iEv);
        if (iEv % (lTreeEvent->GetEntries() / 10) == 0)
            cout << " At Event " << iEv << " out of " << lTreeEvent->GetEntries() << endl;

        hV0MSPD->Fill(fCentrality_V0M, fCentrality_SPDClusters);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 700);
    c1->SetLogz();
    c1->SetBottomMargin(0.15);
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.15);
    c1->SetTopMargin(0.05);
    hV0MSPD->Draw("colz");
    c1->SaveAs("V0MSPD.pdf");
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
    h->SetLineColor(kBlack);
    h->SetMaximum(1E+4);
    //h->SetLineWidth(2);
}