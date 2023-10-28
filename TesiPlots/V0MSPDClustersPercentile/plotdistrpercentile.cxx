void histobeauty(TH1D* h);

void plotdistrpercentile()
{

    TFile* f = TFile::Open("AnalysisResults.root");
    cout<<"--------------- Open Data File --------------------\n"<<endl;
    TTree* lTreeEvent = (TTree *)f->Get("PWGLF_StrVsMult/fTreeEvent");

    Float_t fCentrality_V0M;
    Float_t fCentrality_SPDClusters;
    Float_t fAmplitudeV0A;
    Float_t fAmplitudeV0C;
    Int_t fNSPDClusters;

    //Get from Tree
    lTreeEvent->SetBranchAddress("fCentrality_V0M", &fCentrality_V0M);
    lTreeEvent->SetBranchAddress("fCentrality_SPDClusters", &fCentrality_SPDClusters);
    lTreeEvent->SetBranchAddress("fAmplitudeV0A", &fAmplitudeV0A);
    lTreeEvent->SetBranchAddress("fAmplitudeV0C", &fAmplitudeV0C);
    lTreeEvent->SetBranchAddress("fNSPDClusters", &fNSPDClusters);

    Float_t percentile[10] = {0, 1, 5, 10, 20, 30, 40, 50, 70, 100};
    Int_t colors[10] = {kRed+2,kRed, kOrange, kYellow, kGreen, kGreen+1, kAzure+1, kBlue, kBlue+2};

    TH1D* hV0Mdist = new TH1D("hV0Mdist", ";VZEROM amplitude;counts", 800, 0, 800);
    histobeauty(hV0Mdist);
    TH1D *hV0Mdist_perc[9];
    for (int i = 0; i < 9; i++)
    {
        hV0Mdist_perc[i] = new TH1D(Form("hV0Mdist_perc_%d", i), ";VZEROM amplitude;counts", 800, 0, 800);
        histobeauty(hV0Mdist_perc[i]);
    }

    TH1D *hSPDCldist = new TH1D("hSPDCldist", ";SPD clusters;counts", 450, 0, 450);
    histobeauty(hSPDCldist);
    TH1D *hSPDCldist_perc[9];
    for (int i = 0; i < 9; i++)
    {
        hSPDCldist_perc[i] = new TH1D(Form("hSPDCldist_perc_%d", i), ";SPD clusters;counts", 450, 0, 450);
        histobeauty(hSPDCldist_perc[i]);
    }

    Int_t Nev = 0;
    for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries(); iEv++) {

        lTreeEvent->GetEntry(iEv);
        if( iEv % ( lTreeEvent->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;

        hV0Mdist->Fill(fAmplitudeV0A + fAmplitudeV0C);
        hSPDCldist->Fill(fNSPDClusters);
        for (int i = 0; i < 9; i++)
        {
            if (fCentrality_V0M < percentile[i + 1] && fCentrality_V0M >= percentile[i])
            {
                hV0Mdist_perc[i]->Fill(fAmplitudeV0A + fAmplitudeV0C);
            }
            if (fCentrality_SPDClusters < percentile[i+1] && fCentrality_SPDClusters >= percentile[i])
            {
                hSPDCldist_perc[i]->Fill(fNSPDClusters);
            }
        }
    }
    for (int i = 0; i < 9; i++)
    {
        hV0Mdist_perc[i]->SetLineColor(colors[i]);
        hSPDCldist_perc[i]->SetLineColor(colors[i]);
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
    c1->cd(1)->SetLogy();
    c1->cd(1)->SetBottomMargin(0.15);
    c1->cd(1)->SetRightMargin(0.05);
    c1->cd(1)->SetLeftMargin(0.15);
    c1->cd(1)->SetTopMargin(0.05);
    c1->cd(1)->SetTicks(1, 1);

    TLegend *leg = new TLegend(0.7, 0.5, 0.9, 0.8);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetHeader("VZEROM percentile");

    hV0Mdist->Draw("HIST");
    for (int i = 0; i < 9; i++)
    {
        hV0Mdist_perc[i]->Draw("HIST SAME");
        leg->AddEntry(hV0Mdist_perc[i], Form("%.0f-%.0f %s", percentile[i], percentile[i+1], "%"), "l");
    }
    leg->Draw();
    tex->Draw();
    tex2->Draw();
    c1->cd(2)->SetLogy();
    c1->cd(2)->SetBottomMargin(0.15);
    c1->cd(2)->SetRightMargin(0.05);
    c1->cd(2)->SetLeftMargin(0.15);
    c1->cd(2)->SetTopMargin(0.05);
    c1->cd(2)->SetTicks(1,1);
    c1->cd(2);
    TLegend *leg2 = new TLegend(0.65, 0.5, 0.9, 0.8);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.03);
    leg2->SetHeader("SPDClusters percentile");

    hSPDCldist->Draw("HIST");
    for (int i = 0; i < 9; i++)
    {
        hSPDCldist_perc[i]->Draw("HIST SAME");
        leg2->AddEntry(hSPDCldist_perc[i], Form("%.0f-%.0f %s", percentile[i], percentile[i+1], "%"), "l");
    }
    leg2->Draw();
    tex->Draw();
    tex2->Draw();

    c1->SaveAs("V0MSPDCldist.pdf");
}

void histobeauty(TH1D* h){
    h->SetStats(0);
    h->SetTitle("");
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.035);
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleOffset(1.6);
    h->GetYaxis()->SetLabelSize(0.035);
    h->SetLineColor(kBlack);
    //h->SetLineWidth(2);
}