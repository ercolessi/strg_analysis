void histobeauty(TGraphErrors *h);

void plot(TString est = "dd")
{

    TFile *f = TFile::Open("LHC20i2a-Pythia8_Monash2013.root");
    cout << "--------------- Open Data File --------------------\n" << endl;
    TTree *lTreeEvent = (TTree *)f->Get("fTree");

    Float_t fCent;
    Float_t fCentrality_V0M;
    Float_t fCentrality_SPDCl;
    Int_t nchEta;
    Int_t nKchEta;
    Int_t nK0Eta;

    // Get from Tree
    if (est.Contains("V0M")) lTreeEvent->SetBranchAddress("v0mPerc", &fCent);
    else if (est.Contains("Cl")) lTreeEvent->SetBranchAddress("multSPDcl", &fCent);
    //else if (est.Contains("Tr")) lTreeEvent->SetBranchAddress("multSPDtr", &fCent);
    else if (est.Contains("Tr")) lTreeEvent->SetBranchAddress("multRef5", &fCent);
    else if (est.Contains("dd")){
        lTreeEvent->SetBranchAddress("v0mPerc", &fCentrality_V0M);
        lTreeEvent->SetBranchAddress("multSPDcl", &fCentrality_SPDCl);
    }
    // lTreeEvent->SetBranchAddress("SPDtracklets", &nchEta);
    lTreeEvent->SetBranchAddress("nchEta", &nchEta);
    lTreeEvent->SetBranchAddress("nKchEta", &nKchEta);
    lTreeEvent->SetBranchAddress("nK0Eta", &nK0Eta);

    Float_t percentile[] = {0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100};

    Float_t percentile1_low[] = { 1,30, 30, 20, 0,
    20, 10, 0,
     10, 20, 30, 40,
    50,  30, 40,
    50, 60, 70,70};
    Float_t percentile1_high[] ={ 5,70, 50, 50, 30,
    30, 30, 20,
     20, 30, 40, 50,
    100,  40, 50,
    60, 70, 100, 100};

    Float_t percentile2_low[] = { 0,10, 20, 30, 50,
    0, 10, 20,
     10, 10, 10, 10,
    10,  40, 40,
    40, 40, 40,70};
    Float_t percentile2_high[] ={ 5,30, 40, 50, 100,
    10, 20, 30,
     20, 20, 20, 20,
    20,  50, 50,
    50, 50, 50,100};

    int nb = sizeof(percentile) / sizeof(Float_t) - 1;
    int nc;
    if (est.Contains("dd")) {
        nb = sizeof(percentile1_low) / sizeof(Float_t);
        nc = nb;
     //   nb *= nb;
    }
    const int nbins = nb;

    Float_t Kch[nbins];
    Float_t K0[nbins];
    Float_t ratio[nbins];
    Float_t errratio[nbins];
    Float_t KchErr[nbins];
    Float_t K0Err[nbins];
    Float_t nch[nbins];
    Float_t errnch[nbins];
    Float_t KchMB = 0;
    Float_t K0MB = 0;
    Float_t KchErrMB = 0;
    Float_t K0ErrMB = 0;
    Float_t nchMB = 0;
    Float_t errnchMB = 0;
    Int_t nev[nbins];

    for (int i = 0; i < nbins; i++)
    {
        Kch[i] = 0;
        K0[i] = 0;
        KchErr[i] = 0;
        K0Err[i] = 0;
        nch[i] = 0;
        errnch[i] = 0;
        ratio[i] = 0;
        errratio[i] = 0;
        nev[i] = 0;
    }


    Int_t nevMB = 0;
    for (Long_t iEv = 0; iEv < lTreeEvent->GetEntries(); iEv++)
    {
        lTreeEvent->GetEntry(iEv);
        if (iEv % (lTreeEvent->GetEntries() / 10) == 0)
            cout << " At Event " << iEv << " out of " << lTreeEvent->GetEntries() << endl;

        if (nchEta < 1) continue; //INEL>0

        if (est.Contains("dd")){
            int count = 0;
            for (int i = 0; i < nbins; i++)
            {
                if (fCentrality_V0M >= percentile1_low[i] && fCentrality_V0M < percentile1_high[i])
                {
                    if (fCentrality_SPDCl >= percentile2_low[i] && fCentrality_SPDCl < percentile2_high[i])
                    {
                        Kch[i] += nKchEta;
                        K0[i] += nK0Eta;
                        nch[i] += nchEta;
                        nev[i]++;
                    }
                }
            }
        } else {
            for (int i = 0; i < nbins; i++)
            {
                if (fCent >= percentile[i] && fCent <= percentile[i + 1])
                {
                    Kch[i] += nKchEta;
                    K0[i] += nK0Eta;
                    nch[i] += nchEta;
                    nev[i]++;
                }
            }
        }
        KchMB += nKchEta;
        K0MB += nK0Eta;
        nchMB += nchEta;
        nevMB++;
    }

    for (int i = 0; i < nbins; i++)
    {
        if (nev[i] == 0)
            nev[i] = -1;

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

        cout << percentile1_low[i] << " - " << percentile1_high[i] << " " << percentile2_low[i] << " - " << percentile2_high[i] <<
        " nch: " << nch[i] << " +/- " << errnch[i] << endl;
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