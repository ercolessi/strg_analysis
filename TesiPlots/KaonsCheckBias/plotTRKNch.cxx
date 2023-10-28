void histobeauty(TGraphErrors *h);

void plotTRKNch()
{

    TFile *f = TFile::Open("LHC20i2a-Pythia8_Monash2013.root");
    cout << "--------------- Open Data File --------------------\n"
         << endl;
    TTree *lTreeEvent = (TTree *)f->Get("fTree");

    Float_t fCentrality_V0M;
    Float_t fCentrality_SPDClusters;
    Float_t fCentrality_RefMultEta5;
    Int_t SPDtracklets;
    Int_t nchEta;
    Int_t nKchEta;
    Int_t nK0Eta;

    // Get from Tree
    lTreeEvent->SetBranchAddress("v0mPerc", &fCentrality_V0M);
    lTreeEvent->SetBranchAddress("multSPDcl", &fCentrality_SPDClusters);
    lTreeEvent->SetBranchAddress("multRef5", &fCentrality_RefMultEta5);
    lTreeEvent->SetBranchAddress("nchEta", &nchEta);
    lTreeEvent->SetBranchAddress("SPDtracklets", &SPDtracklets);
    lTreeEvent->SetBranchAddress("nKchEta", &nKchEta);
    lTreeEvent->SetBranchAddress("nK0Eta", &nK0Eta);

    Float_t percentile[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 100};
    const int nbins = sizeof(percentile) / sizeof(Float_t) - 1;

    Float_t Kch[3][nbins];
    Float_t K0[3][nbins];
    Float_t ratio[3][nbins];
    Float_t errratio[3][nbins];
    Float_t KchErr[3][nbins];
    Float_t K0Err[3][nbins];
    Float_t nch[3][nbins];
    Float_t errnch[3][nbins];
    Float_t spdtrk[3][nbins];
    Float_t errspdtrk[3][nbins];
    Int_t nev[3][nbins];

    for (int j = 0; j < 3; j++)
    {
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

    for (Long_t iEv = 0; iEv < lTreeEvent->GetEntries() / 10; iEv++)
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
            if (fCentrality_RefMultEta5 >= percentile[i] && fCentrality_RefMultEta5 < percentile[i + 1])
            {
                Kch[2][i] += nKchEta;
                K0[2][i] += nK0Eta;
                nch[2][i] += nchEta;
                nev[2][i]++;
            }
        }
    }

    for (int i = 0; i < nbins; i++)
    {
        errnch[i] = sqrt(sqrt(nch[i])/nch[i] * sqrt(nch[i])/nch[i] + sqrt(nev[i])/nev[i] * sqrt(nev[i])/nev[i]);
        KchErr[i] = sqrt(sqrt(Kch[i])/Kch[i] * sqrt(Kch[i])/Kch[i] + sqrt(nev[i])/nev[i] * sqrt(nev[i])/nev[i]);
        K0Err[i] = sqrt(sqrt(K0[i])/K0[i] * sqrt(K0[i])/K0[i] + sqrt(nev[i])/nev[i] * sqrt(nev[i])/nev[i]);

        Kch[i] /= nev[i];
        K0[i] /= nev[i];
        nch[i] /= nev[i];
        KchErr[i] *= Kch[i];
        K0Err[i] *= K0[i];
        errnch[i] *= nch[i];

        ratio[i] = Kch[i] / K0[i];
        errratio[i] = sqrt(ratio[i] * ratio[i] * (KchErr[i] * KchErr[i] / (Kch[i] * Kch[i]) + K0Err[i] * K0Err[i] / (K0[i] * K0[i])));

        cout << "nch: " << nch[i] << " +/- " << errnch[i] << endl;
    }

    KchErrMB = sqrt(sqrt(KchMB) / KchMB * sqrt(KchMB) / KchMB + sqrt(nevMB) / nevMB * sqrt(nevMB) / nevMB);
    K0ErrMB = sqrt(sqrt(K0MB) / K0MB * sqrt(K0MB) / K0MB + sqrt(nevMB) / nevMB * sqrt(nevMB) / nevMB);
    errnchMB = sqrt(sqrt(nchMB) / nchMB * sqrt(nchMB) / nchMB + sqrt(nevMB) / nevMB * sqrt(nevMB) / nevMB);

    KchMB /= nevMB;
    K0MB /= nevMB;
    nchMB /= nevMB;
    KchErrMB *= KchMB;
    K0ErrMB *= K0MB;
    errnchMB *= nchMB;

    for (int i = 0; i < nbins; i++)
    {
        Kch[i] /= KchMB;
        K0[i] /= K0MB;
        KchErr[i] /= KchMB;
        K0Err[i] /= K0MB;
    }

    TGraphErrors *gKch = new TGraphErrors(nbins, nch, Kch, errnch, KchErr);
    gKch->SetMarkerStyle(kFullCircle);
    gKch->SetMarkerColor(kBlack);
    gKch->SetLineColor(kBlack);
    gKch->SetLineWidth(2);
    gKch->SetMarkerSize(1.2);
    gKch->SetTitle(";n_{ch} (|#eta|<0.5);");
    histobeauty(gKch);
    gKch->GetYaxis()->SetRangeUser(0, 6);

    TGraphErrors *gK0 = new TGraphErrors(nbins, nch, K0, errnch, K0Err);
    gK0->SetMarkerStyle(24);
    gK0->SetMarkerColor(kBlack);
    gK0->SetLineColor(kBlack);
    gK0->SetLineWidth(2);
    gK0->SetMarkerSize(1.2);
    gK0->SetTitle("K^{0}_{S};n_{ch};d#textit{N}d#textit{y}");
    histobeauty(gK0);

    TGraphErrors *gRatio = new TGraphErrors(nbins, nch, ratio, errnch, errratio);
    gRatio->SetMarkerStyle(kFullCircle);
    gRatio->SetMarkerColor(kBlack);
    gRatio->SetLineColor(kBlack);
    gRatio->SetLineWidth(2);
    gRatio->SetMarkerSize(0.8);
    gRatio->SetTitle(";n_{ch} (|#eta|<0.5);");
    histobeauty(gRatio);
    gRatio->GetYaxis()->SetRangeUser(0, 2);
    gRatio->GetXaxis()->SetRangeUser(0, 26);

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
    c1->SetLogz();
    c1->SetBottomMargin(0.15);
    c1->SetRightMargin(0.05);
    c1->SetLeftMargin(0.15);
    c1->SetTopMargin(0.05);
    c1->SetTicks(1, 1);
    gKch->Draw("AP");
    gK0->Draw("EP SAME");

    TLegend *leg = new TLegend(0.2, 0.8, 0.5, 0.9);
    leg->AddEntry(gKch, "K^{#pm}/#LTK^{#pm}#GT_{MB}", "lep");
    leg->AddEntry(gK0, "K^{0}/#LTK^{0}#GT_{MB}", "lep");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->Draw("SAME");

    TLatex *ltx = new TLatex();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.035);
    ltx->SetTextAlign(12);
    if (est.Contains("V0M")) ltx->DrawLatexNDC(0.22, 0.75, "VZEROM classes ");
    if (est.Contains("Cl")) ltx->DrawLatexNDC(0.22, 0.75, "SPDClusters classes ");
    if (est.Contains("Tr")) ltx->DrawLatexNDC(0.22, 0.75, "RefMultEta5 classes ");
    if (est.Contains("dd")) ltx->DrawLatexNDC(0.22, 0.75, "VZEROM + SPDClusters classes ");


    TLatex *tex = new TLatex(0.5, 0.88, "ALICE, pp #sqrt{s} = 13 TeV");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.035);
    tex->SetLineWidth(2);
    TLatex *tex2 = new TLatex(0.5, 0.82, "This work");
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.035);
    tex2->SetLineWidth(2);

    tex->Draw();
    tex2->Draw();

    c1->SaveAs(Form("KchK0_%s.pdf",est.Data()));

}

void histobeauty(TGraphErrors* h){
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.035);
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleOffset(1.6);
    h->GetYaxis()->SetLabelSize(0.035);
}