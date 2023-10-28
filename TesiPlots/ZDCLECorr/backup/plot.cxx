void histobeauty(TGraphErrors *h);

TFile *f = TFile::Open("LHC20i2a-Pythia8_Monash2013.root");
TTree *lTreeEvent = (TTree *)f->Get("fTree");

float GetExpression(const char *formula);

void plot(TString est = "dd")
{

    TFile *f = TFile::Open("LHC20i2a-Pythia8_Monash2013.root");
    cout << "--------------- Open Data File --------------------\n" << endl;
    TTree *lTreeEvent = (TTree *)f->Get("fTree");

    Float_t fCentrality_V0M;
    Float_t fCentrality_SPDClusters;
    Float_t adcZDCN1;
    Float_t adcZDCN2;
    Float_t effEnergy;
    Float_t E_n_leadC;
    Float_t E_n_leadA;
    Float_t E_p_leadC;
    Float_t E_p_leadA;
    Int_t SPDtracklets;

    // Get from Tree
    lTreeEvent->SetBranchAddress("v0mPerc", &fCentrality_V0M);
    lTreeEvent->SetBranchAddress("multSPDcl", &fCentrality_SPDClusters);
    lTreeEvent->SetBranchAddress("adcZDCN1", &adcZDCN1);
    lTreeEvent->SetBranchAddress("adcZDCN2", &adcZDCN2);
    lTreeEvent->SetBranchAddress("effEnergy", &effEnergy);
    lTreeEvent->SetBranchAddress("E_n_leadC", &E_n_leadC);
    lTreeEvent->SetBranchAddress("E_n_leadA", &E_n_leadA);
    lTreeEvent->SetBranchAddress("E_p_leadC", &E_p_leadC);
    lTreeEvent->SetBranchAddress("E_p_leadA", &E_p_leadA);
    lTreeEvent->SetBranchAddress("SPDtracklets", &SPDtracklets);

    Float_t percentile[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    int nbins = sizeof(percentile) / sizeof(Float_t) - 1;

    Float_t ZN[nbins];
    Float_t ErrZN[nbins];
    Float_t LE[nbins];
    Float_t ErrLE[nbins];
    Float_t ZN2[nbins];
    Float_t ErrZN2[nbins];
    Float_t LE2[nbins];
    Float_t ErrLE2[nbins];
    Int_t nev[nbins];
    TH1F* hEE[nbins];
    TH1F *hZN[nbins];
    TH1F *hEE2[nbins];
    TH1F *hZN2[nbins];

    for (int i = 0; i < nbins; i++)
    {
        ZN[i] = 0;
        ErrZN[i] = 0;
        LE[i] = 0;
        ErrLE[i] = 0;
        nev[i] = 0;

        hEE[i] = new TH1F(Form("hEE%d", i), ";Effective Energy [GeV]", 200, 0, 13000);
        hZN[i] = new TH1F(Form("hZN%d", i), ";Leading Energy (E_{|#eta|>8}) [GeV]", 200, 0, 13000);
        hEE2[i] = new TH1F(Form("hEE2%d", i), ";Effective Energy (GeV)", 200, 0, 13000);
        hZN2[i] = new TH1F(Form("hZN2%d", i), ";Leading Energy (GeV)", 200, 0, 13000);
    }

    float x, y;

    for (Long_t iEv = 0; iEv < lTreeEvent->GetEntries()/10; iEv++)
    {
        lTreeEvent->GetEntry(iEv);
        if (iEv % (lTreeEvent->GetEntries() / 10) == 0)
            cout << " At Event " << iEv << " out of " << lTreeEvent->GetEntries() << endl;

        if (SPDtracklets < 1)
            continue; // INEL>0

        if (fCentrality_V0M < 0 || fCentrality_V0M >= 100)
            continue;

        for (int i = 0; i < nbins; i++)
        {
            if (fCentrality_V0M >= percentile[i] && fCentrality_V0M <= percentile[i + 1])
            {
                LE[i]+=TMath::Abs(lTreeEvent->GetLeaf("effEnergy")->GetValue());
                ZN[i]+=adcZDCN1 + adcZDCN2;
                nev[i]++;
                hEE[i]->Fill(TMath::Abs(lTreeEvent->GetLeaf("effEnergy")->GetValue()));
                hZN[i]->Fill(adcZDCN1 + adcZDCN2);
            }
            if (fCentrality_SPDClusters >= percentile[i] && fCentrality_SPDClusters <= percentile[i + 1])
            {
                hEE2[i]->Fill(TMath::Abs(lTreeEvent->GetLeaf("effEnergy")->GetValue()));
                hZN2[i]->Fill(adcZDCN1 + adcZDCN2);
            }
        }
    }

    for (int i = 0; i < nbins; i++)
    {
        ZN[i] = hZN[i]->GetMean();
        LE[i] = hEE[i]->GetMean();
        ErrZN[i] = hZN[i]->GetMeanError();
        ErrLE[i] = hEE[i]->GetMeanError();
        ZN2[i] = hZN2[i]->GetMean();
        LE2[i] = hEE2[i]->GetMean();
        ErrZN2[i] = hZN2[i]->GetMeanError();
        ErrLE2[i] = hEE2[i]->GetMeanError();
    }

    TGraphErrors *gZN = new TGraphErrors(nbins, LE, ZN, ErrLE, ErrZN);
    gZN->SetMarkerStyle(kFullCircle);
    gZN->SetMarkerColor(kRed);
    gZN->SetLineColor(kRed);
    gZN->SetLineWidth(2);
    gZN->SetMarkerSize(2.);
    gZN->SetTitle(";Leading energy (GeV); ZN (a.u.) ");
    histobeauty(gZN);

    TGraphErrors *gZN2 = new TGraphErrors(nbins, LE2, ZN2, ErrLE2, ErrZN2);
    gZN2->SetMarkerStyle(kFullCircle);
    gZN2->SetMarkerColor(kBlue);
    gZN2->SetLineColor(kBlue);
    gZN2->SetLineWidth(2);
    gZN2->SetMarkerSize(2.);
    gZN2->SetTitle(";Leading energy (GeV); ZN (a.u.) ");
    histobeauty(gZN2);

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
    c1->SetBottomMargin(0.15);
    c1->SetRightMargin(0.05);
    c1->SetLeftMargin(0.15);
    c1->SetTopMargin(0.05);
    c1->SetTicks(1, 1);
    TH1D* hframe = new TH1D("hframe", ";Leading energy (GeV); ZN (a.u.) ", 100, 0, 13000);
    hframe->GetYaxis()->SetRangeUser(0, 700);
    hframe->SetStats(0);
    hframe->Draw();
    gZN->Draw("EP SAME");
    gZN2->Draw("EP SAME");
    TLatex *ltx = new TLatex();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.035);
    ltx->SetTextAlign(12);

    TLine *line = new TLine(0, 0, 12000, 700);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->SetLineColor(kBlack);
    line->Draw("SAME");


    TLatex *tex = new TLatex(0.2, 0.88, "ALICE, pp #sqrt{s} = 13 TeV");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.035);
    tex->SetLineWidth(2);
    TLatex *tex2 = new TLatex(0.2, 0.82, "This work");
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.035);
    tex2->SetLineWidth(2);

    tex->Draw();
    tex2->Draw();

    c1->SaveAs(Form("ZNvsLE.pdf"));

}

void histobeauty(TGraphErrors* h){
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.035);
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetTitleOffset(1.6);
    h->GetYaxis()->SetLabelSize(0.035);
}

float GetExpression(const char *formula)
{
    istringstream iss(formula);

    float result = 0, coef;

    do
    {
        string subs;

        // Get the word from the istringstream
        iss >> subs;

        if (subs[0] == '-')
        {
            sscanf(subs.c_str(), "%f", &coef);
            iss >> subs;
            if (subs.find("effEnergy") < 10000)
                coef *= lTreeEvent->GetLeaf(subs.c_str())->GetValue() + 13000;
            else
                coef *= lTreeEvent->GetLeaf(subs.c_str())->GetValue();
            result += coef;
        }
        else if (subs[0] == '+')
        {
            sscanf(subs.c_str(), "%f", &coef);
            iss >> subs;
            if (subs.find("effEnergy") < 10000)
                coef *= lTreeEvent->GetLeaf(subs.c_str())->GetValue() + 13000;
            else
                coef *= lTreeEvent->GetLeaf(subs.c_str())->GetValue();
            result += coef;
        }
        else if (subs[0] != '\0')
        {
            if (subs.find("effEnergy") < 10000)
                coef = lTreeEvent->GetLeaf(subs.c_str())->GetValue() + 13000;
            else
                coef = lTreeEvent->GetLeaf(subs.c_str())->GetValue();
            result += coef;
        }
    } while (iss);
    return result;
}
