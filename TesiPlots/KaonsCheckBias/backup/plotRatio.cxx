void histobeauty(TGraphErrors *h);

void plotRatio()
{

    TFile *f = TFile::Open("LHC20i2a-Pythia8_Monash2013.root");
    cout << "--------------- Open Data File --------------------\n" << endl;
    TTree *lTreeEvent = (TTree *)f->Get("fTree");

    Float_t fCentrality_V0M;
    Float_t fCentrality_SPDClusters;
    Float_t fCentrality_SPDTracklets;
    Int_t nchEta;
    Int_t nKchEta;
    Int_t nK0Eta;

    // Get from Tree
    lTreeEvent->SetBranchAddress("v0mPerc", &fCentrality_V0M);
    lTreeEvent->SetBranchAddress("multSPDcl", &fCentrality_SPDClusters);
    lTreeEvent->SetBranchAddress("multSPDtr", &fCentrality_SPDTracklets);
    lTreeEvent->SetBranchAddress("nchEta", &nchEta);
    lTreeEvent->SetBranchAddress("nKchEta", &nKchEta);
    lTreeEvent->SetBranchAddress("nK0Eta", &nK0Eta);

    Float_t percentile[] = {0, 10, 20, 30, 40, 50, 70, 100};
    const int nbins = sizeof(percentile) / sizeof(Float_t) - 1;

    Float_t Kch[3][nbins];
    Float_t K0[3][nbins];
    Float_t ratio[3][nbins];
    Float_t errratio[3][nbins];
    Float_t KchErr[3][nbins];
    Float_t K0Err[3][nbins];
    Float_t nch[3][nbins];
    Float_t errnch[3][nbins];
    Int_t nev[3][nbins];

    for (int j = 0; j < 3; j++){
        for (int i = 0; i < nbins; i++)
        {
            Kch[j][i] = 0;
            K0[j][i] = 0;
            KchErr[j][i] = 0;
            K0Err[j][i] = 0;
            nch[j][i] = 0;
            errnch[j][i] = 0;
            ratio[j][i] = 0;
            errratio[j][i] = 0;
            nev[j][i] = 0;
        }
    }


    for (Long_t iEv = 0; iEv < lTreeEvent->GetEntries()/10; iEv++)
    {
        lTreeEvent->GetEntry(iEv);
        if (iEv % (lTreeEvent->GetEntries() / 10) == 0)
            cout << " At Event " << iEv << " out of " << lTreeEvent->GetEntries() << endl;

        for (int i = 0; i < nbins; i++)
        {
            if (fCentrality_V0M >= percentile[i] && fCentrality_V0M < percentile[i + 1])
            {
                Kch[0][i] += nKchEta;
                K0[0][i] += nK0Eta;
                nch[0][i] += nchEta;
                nev[0][i]++;
            }
            if (fCentrality_SPDClusters >= percentile[i] && fCentrality_SPDClusters < percentile[i + 1])
            {
                Kch[1][i] += nKchEta;
                K0[1][i] += nK0Eta;
                nch[1][i] += nchEta;
                nev[1][i]++;
            }
            if (fCentrality_SPDTracklets >= percentile[i] && fCentrality_SPDTracklets < percentile[i + 1])
            {
                Kch[2][i] += nKchEta;
                K0[2][i] += nK0Eta;
                nch[2][i] += nchEta;
                nev[2][i]++;
            }
        }
    }

    for (int j = 0; j < 3; j++){
    for (int i = 0; i < nbins; i++)
        {
            errnch[j][i] = sqrt(sqrt(nch[j][i])/nch[j][i] * sqrt(nch[j][i])/nch[j][i] + sqrt(nev[j][i])/nev[j][i] * sqrt(nev[j][i])/nev[j][i]);
            KchErr[j][i] = sqrt(sqrt(Kch[j][i])/Kch[j][i] * sqrt(Kch[j][i])/Kch[j][i] + sqrt(nev[j][i])/nev[j][i] * sqrt(nev[j][i])/nev[j][i]);
            K0Err[j][i] = sqrt(sqrt(K0[j][i])/K0[j][i] * sqrt(K0[j][i])/K0[j][i] + sqrt(nev[j][i])/nev[j][i] * sqrt(nev[j][i])/nev[j][i]);

            Kch[j][i] /= nev[j][i];
            K0[j][i] /= nev[j][i];
            nch[j][i] /= nev[j][i];
            KchErr[j][i] *= Kch[j][i];
            K0Err[j][i] *= K0[j][i];
            errnch[j][i] *= nch[j][i];

            ratio[j][i] = Kch[j][i] / K0[j][i];
            errratio[j][i] = sqrt(ratio[j][i] * ratio[j][i] * (KchErr[j][i] * KchErr[j][i] / (Kch[j][i] * Kch[j][i]) + K0Err[j][i] * K0Err[j][i] / (K0[j][i] * K0[j][i])));
        }
    }

    TGraphErrors *gRatio[3];
    Int_t colors[] = {kRed, kBlue, kGreen+1};
    for (int j = 0; j < 3; j++){
        gRatio[j] = new TGraphErrors(nbins, nch[j], ratio[j], errnch[j], errratio[j]);
        gRatio[j]->SetMarkerStyle(kFullCircle);
        gRatio[j]->SetMarkerColor(j+1);
        gRatio[j]->SetLineColor(j+1);
        gRatio[j]->SetLineWidth(2);
        gRatio[j]->SetMarkerSize(2.0);
        gRatio[j]->GetYaxis()->SetRangeUser(0.0, 2);
        gRatio[j]->SetTitle(";n_{ch} (|#eta|<0.5);K^{#pm}/K^{0}");
        histobeauty(gRatio[j]);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
    c1->SetLogz();
    c1->SetBottomMargin(0.15);
    c1->SetRightMargin(0.05);
    c1->SetLeftMargin(0.15);
    c1->SetTopMargin(0.05);
    c1->SetTicks(1, 1);
    gRatio[0]->Draw("AP");
    gRatio[1]->Draw("EP SAME");
    gRatio[2]->Draw("EP SAME");

    TLegend *leg = new TLegend(0.2, 0.75, 0.4, 0.9);
    leg->SetBorderSize(0);
    leg->AddEntry(gRatio[0], "VZEROM", "lpe");
    leg->AddEntry(gRatio[1], "SPDClusters", "lpe");
    leg->AddEntry(gRatio[2], "SPDTracklets", "lpe");
    leg->Draw("SAME");

    TLatex *tex = new TLatex(0.45, 0.86, "ALICE, pp #sqrt{s} = 13 TeV");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    TLatex *tex2 = new TLatex(0.45, 0.82, "This work");
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(2);

    tex->Draw();
    tex2->Draw();

    c1->SaveAs("Kch_K0_SPDcl.pdf");
}

void histobeauty(TGraphErrors* h){
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.035);
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleOffset(1.6);
    h->GetYaxis()->SetLabelSize(0.035);
}