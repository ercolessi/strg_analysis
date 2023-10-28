void histobeauty(TGraphErrors *h);

TFile *f = TFile::Open("LHC20i2a-Pythia8_Monash2013.root");
TTree *lTreeEvent = (TTree *)f->Get("fTree");

float GetExpression(const char *formula);

void plotAll()
{

    TFile *f1 = TFile::Open("LHC20i2a-Pythia8_Monash2013.root");
    TFile *f2 = TFile::Open("LHC16d3-EPOS.root");
    TFile *f3 = TFile::Open("LHC17h7a-Pythia6_Perugia2011.root");
    TFile *f4 = TFile::Open("LHC17h7b-Phojet.root");
    cout << "--------------- Open Data File --------------------\n" << endl;
    TTree *lt1 = (TTree *)f1->Get("fTree");
    TTree *lt2 = (TTree *)f2->Get("fTree");
    TTree *lt3 = (TTree *)f3->Get("fTree");
    TTree *lt4 = (TTree *)f4->Get("fTree");

    Float_t percentile[] = {0, 10, 20, 30, 40, 50, 60, 70,  100};
    int nbins = sizeof(percentile) / sizeof(Float_t) - 1;

    TH1D *hframe = new TH1D("hframe", ";Leading Energy E_{|#eta|>8} (GeV); ZN (a.u.)", 100, 0, 13000);
    hframe->GetYaxis()->SetRangeUser(0, 700);
    hframe->GetXaxis()->SetTitleSize(0.04);
    hframe->GetXaxis()->SetTitleOffset(1.1);
    hframe->GetYaxis()->SetTitleSize(0.04);
    hframe->SetStats(0);

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

    Float_t fCentrality_V0M;
    Float_t fCentrality_SPDClusters;
    Float_t adcZDCN1;
    Float_t adcZDCN2;
    Float_t effEnergy;
    Int_t SPDtracklets;

    // Get from Tree
    lt1->SetBranchAddress("v0mPerc", &fCentrality_V0M);
    lt1->SetBranchAddress("multSPDcl", &fCentrality_SPDClusters);
    lt1->SetBranchAddress("adcZDCN1", &adcZDCN1);
    lt1->SetBranchAddress("adcZDCN2", &adcZDCN2);
    lt1->SetBranchAddress("effEnergy", &effEnergy);
    lt1->SetBranchAddress("SPDtracklets", &SPDtracklets);

    Float_t ZNv0m_1[nbins];
    Float_t ErrZNv0m_1[nbins];
    Float_t LEv0m_1[nbins];
    Float_t ErrLEv0m_1[nbins];
    Float_t ZNspd_1[nbins];
    Float_t ErrZNspd_1[nbins];
    Float_t LEspd_1[nbins];
    Float_t ErrLEspd_1[nbins];
    TH1F* hLEv0m_1[nbins];
    TH1F *hZNv0m_1[nbins];
    TH1F *hLEspd_1[nbins];
    TH1F *hZNspd_1[nbins];

    for (int i = 0; i < nbins; i++)
    {
        ZNv0m_1[i]=0;
        ErrZNv0m_1[i]=0;
        LEv0m_1[i]=0;
        ErrLEv0m_1[i]=0;
        ZNspd_1[i]=0;
        ErrZNspd_1[i]=0;
        LEspd_1[i]=0;
        ErrLEspd_1[i]=0;

        hLEv0m_1[i] = new TH1F(Form("hLEv0m_1%d", i), ";Leading Energy (GeV)", 200, 0, 13000);
        hZNv0m_1[i] = new TH1F(Form("hZNv0m_1%d", i), ";ZN Energy (a.u.)", 200, 0, 13000);
        hLEspd_1[i] = new TH1F(Form("hLEspd_1%d", i), ";Leading Energy (GeV)", 200, 0, 13000);
        hZNspd_1[i] = new TH1F(Form("hZNspd_1%d", i), ";ZN Energy (a.u.)", 200, 0, 13000);
    }

    for (Long_t iEv = 0; iEv < lt1->GetEntries()/100; iEv++)
    {
        lt1->GetEntry(iEv);
        if (iEv % (lt1->GetEntries() / 10) == 0)
            cout << " At Event " << iEv << " out of " << lt1->GetEntries() << endl;

        if (SPDtracklets < 1)
            continue; // INEL>0

        if (fCentrality_V0M < 0 || fCentrality_V0M >= 100)
            continue;

        for (int i = 0; i < nbins; i++)
        {
            if (fCentrality_V0M >= percentile[i] && fCentrality_V0M <= percentile[i + 1])
            {
                hLEv0m_1[i]->Fill(TMath::Abs(lt1->GetLeaf("effEnergy")->GetValue()));
                hZNv0m_1[i]->Fill(adcZDCN1 + adcZDCN2);
            }
            if (fCentrality_SPDClusters >= percentile[i] && fCentrality_SPDClusters <= percentile[i + 1])
            {
                hLEspd_1[i]->Fill(TMath::Abs(lt1->GetLeaf("effEnergy")->GetValue()));
                hZNspd_1[i]->Fill(adcZDCN1 + adcZDCN2);
            }
        }
    }

    for (int i = 0; i < nbins; i++)
    {
        ZNv0m_1[i] = hZNv0m_1[i]->GetMean();
        LEv0m_1[i] = hLEv0m_1[i]->GetMean();
        ErrZNv0m_1[i] = hZNv0m_1[i]->GetMeanError();
        ErrLEv0m_1[i] = hLEv0m_1[i]->GetMeanError();
        ZNspd_1[i] = hZNspd_1[i]->GetMean();
        LEspd_1[i] = hLEspd_1[i]->GetMean();
        ErrZNspd_1[i] = hZNspd_1[i]->GetMeanError();
        ErrLEspd_1[i] = hLEspd_1[i]->GetMeanError();
    }

    TGraphErrors *gZNv0m1 = new TGraphErrors(nbins, LEv0m_1, ZNv0m_1, ErrLEv0m_1, ErrZNv0m_1);
    gZNv0m1->SetMarkerStyle(kFullCircle);
    gZNv0m1->SetMarkerColor(kRed);
    gZNv0m1->SetLineColor(kRed);
    gZNv0m1->SetLineWidth(2);
    gZNv0m1->SetMarkerSize(2.);
    gZNv0m1->SetTitle(";Leading energy (GeV); ZN (a.u.) ");
    histobeauty(gZNv0m1);

    TGraphErrors *gZNspd1 = new TGraphErrors(nbins, LEspd_1, ZNspd_1, ErrLEspd_1, ErrZNspd_1);
    gZNspd1->SetMarkerStyle(kOpenCircle);
    gZNspd1->SetMarkerColor(kRed);
    gZNspd1->SetLineColor(kRed);
    gZNspd1->SetLineWidth(2);
    gZNspd1->SetMarkerSize(2.);
    gZNspd1->SetTitle(";Leading energy (GeV); ZN (a.u.) ");
    histobeauty(gZNspd1);

    /////////////////////////////////////////
    // Get from Tree
    lt2->SetBranchAddress("v0mPerc", &fCentrality_V0M);
    lt2->SetBranchAddress("multSPDcl", &fCentrality_SPDClusters);
    lt2->SetBranchAddress("adcZDCN1", &adcZDCN1);
    lt2->SetBranchAddress("adcZDCN2", &adcZDCN2);
    lt2->SetBranchAddress("effEnergy", &effEnergy);
    lt2->SetBranchAddress("SPDtracklets", &SPDtracklets);

    Float_t ZNv0m_2[nbins];
    Float_t ErrZNv0m_2[nbins];
    Float_t LEv0m_2[nbins];
    Float_t ErrLEv0m_2[nbins];
    Float_t ZNspd_2[nbins];
    Float_t ErrZNspd_2[nbins];
    Float_t LEspd_2[nbins];
    Float_t ErrLEspd_2[nbins];
    TH1F *hLEv0m_2[nbins];
    TH1F *hZNv0m_2[nbins];
    TH1F *hLEspd_2[nbins];
    TH1F *hZNspd_2[nbins];

    for (int i = 0; i < nbins; i++)
    {
        ZNv0m_2[i] = 0;
        ErrZNv0m_2[i] = 0;
        LEv0m_2[i] = 0;
        ErrLEv0m_2[i] = 0;
        ZNspd_2[i] = 0;
        ErrZNspd_2[i] = 0;
        LEspd_2[i] = 0;
        ErrLEspd_2[i] = 0;

        hLEv0m_2[i] = new TH1F(Form("hLEv0m_2%d", i), ";Leading Energy (GeV)", 200, 0, 13000);
        hZNv0m_2[i] = new TH1F(Form("hZNv0m_2%d", i), ";ZN Energy (a.u.)", 200, 0, 13000);
        hLEspd_2[i] = new TH1F(Form("hLEspd_2%d", i), ";Leading Energy (GeV)", 200, 0, 13000);
        hZNspd_2[i] = new TH1F(Form("hZNspd_2%d", i), ";ZN Energy (a.u.)", 200, 0, 13000);
    }

    for (Long_t iEv = 0; iEv < lt2->GetEntries() / 100; iEv++)
    {
        lt2->GetEntry(iEv);
        if (iEv % (lt2->GetEntries() / 10) == 0)
            cout << " At Event " << iEv << " out of " << lt2->GetEntries() << endl;

        if (SPDtracklets < 1)
            continue; // INEL>0

        if (fCentrality_V0M < 0 || fCentrality_V0M >= 100)
            continue;

        for (int i = 0; i < nbins; i++)
        {
            if (fCentrality_V0M >= percentile[i] && fCentrality_V0M <= percentile[i + 1])
            {
                hLEv0m_2[i]->Fill(TMath::Abs(lt2->GetLeaf("effEnergy")->GetValue()));
                hZNv0m_2[i]->Fill(adcZDCN1 + adcZDCN2);
            }
            if (fCentrality_SPDClusters >= percentile[i] && fCentrality_SPDClusters <= percentile[i + 1])
            {
                hLEspd_2[i]->Fill(TMath::Abs(lt2->GetLeaf("effEnergy")->GetValue()));
                hZNspd_2[i]->Fill(adcZDCN1 + adcZDCN2);
            }
        }
    }

    for (int i = 0; i < nbins; i++)
    {
        ZNv0m_2[i] = hZNv0m_2[i]->GetMean();
        LEv0m_2[i] = hLEv0m_2[i]->GetMean();
        ErrZNv0m_2[i] = hZNv0m_2[i]->GetMeanError();
        ErrLEv0m_2[i] = hLEv0m_2[i]->GetMeanError();
        ZNspd_2[i] = hZNspd_2[i]->GetMean();
        LEspd_2[i] = hLEspd_2[i]->GetMean();
        ErrZNspd_2[i] = hZNspd_2[i]->GetMeanError();
        ErrLEspd_2[i] = hLEspd_2[i]->GetMeanError();
    }

    TGraphErrors *gZNv0m2 = new TGraphErrors(nbins, LEv0m_2, ZNv0m_2, ErrLEv0m_2, ErrZNv0m_2);
    gZNv0m2->SetMarkerStyle(kFullCircle);
    gZNv0m2->SetMarkerColor(kBlue);
    gZNv0m2->SetLineColor(kBlue);
    gZNv0m2->SetLineWidth(2);
    gZNv0m2->SetMarkerSize(2.);
    gZNv0m2->SetTitle(";Leading energy (GeV); ZN (a.u.) ");
    histobeauty(gZNv0m2);

    TGraphErrors *gZNspd2 = new TGraphErrors(nbins, LEspd_2, ZNspd_2, ErrLEspd_2, ErrZNspd_2);
    gZNspd2->SetMarkerStyle(kOpenCircle);
    gZNspd2->SetMarkerColor(kBlue);
    gZNspd2->SetLineColor(kBlue);
    gZNspd2->SetLineWidth(2);
    gZNspd2->SetMarkerSize(2.);
    gZNspd2->SetTitle(";Leading energy (GeV); ZN (a.u.) ");
    histobeauty(gZNspd2);

    /////////////////////////////////////////
    // Get from Tree
    lt3->SetBranchAddress("v0mPerc", &fCentrality_V0M);
    lt3->SetBranchAddress("multSPDcl", &fCentrality_SPDClusters);
    lt3->SetBranchAddress("adcZDCN1", &adcZDCN1);
    lt3->SetBranchAddress("adcZDCN2", &adcZDCN2);
    lt3->SetBranchAddress("effEnergy", &effEnergy);
    lt3->SetBranchAddress("SPDtracklets", &SPDtracklets);

    Float_t ZNv0m_3[nbins];
    Float_t ErrZNv0m_3[nbins];
    Float_t LEv0m_3[nbins];
    Float_t ErrLEv0m_3[nbins];
    Float_t ZNspd_3[nbins];
    Float_t ErrZNspd_3[nbins];
    Float_t LEspd_3[nbins];
    Float_t ErrLEspd_3[nbins];
    TH1F *hLEv0m_3[nbins];
    TH1F *hZNv0m_3[nbins];
    TH1F *hLEspd_3[nbins];
    TH1F *hZNspd_3[nbins];

    for (int i = 0; i < nbins; i++)
    {
        ZNv0m_3[i] = 0;
        ErrZNv0m_3[i] = 0;
        LEv0m_3[i] = 0;
        ErrLEv0m_3[i] = 0;
        ZNspd_3[i] = 0;
        ErrZNspd_3[i] = 0;
        LEspd_3[i] = 0;
        ErrLEspd_3[i] = 0;

        hLEv0m_3[i] = new TH1F(Form("hLEv0m_3%d", i), ";Leading Energy (GeV)", 200, 0, 13000);
        hZNv0m_3[i] = new TH1F(Form("hZNv0m_3%d", i), ";ZN Energy (a.u.)", 200, 0, 13000);
        hLEspd_3[i] = new TH1F(Form("hLEspd_3%d", i), ";Leading Energy (GeV)", 200, 0, 13000);
        hZNspd_3[i] = new TH1F(Form("hZNspd_3%d", i), ";ZN Energy (a.u.)", 200, 0, 13000);
    }

    for (Long_t iEv = 0; iEv < lt3->GetEntries() / 100; iEv++)
    {
        lt3->GetEntry(iEv);
        if (iEv % (lt3->GetEntries() / 10) == 0)
            cout << " At Event " << iEv << " out of " << lt3->GetEntries() << endl;

        if (SPDtracklets < 1)
            continue; // INEL>0

        if (fCentrality_V0M < 0 || fCentrality_V0M >= 100)
            continue;

        for (int i = 0; i < nbins; i++)
        {
            if (fCentrality_V0M >= percentile[i] && fCentrality_V0M <= percentile[i + 1])
            {
                hLEv0m_3[i]->Fill(TMath::Abs(lt3->GetLeaf("effEnergy")->GetValue()));
                hZNv0m_3[i]->Fill(adcZDCN1 + adcZDCN2);
            }
            if (fCentrality_SPDClusters >= percentile[i] && fCentrality_SPDClusters <= percentile[i + 1])
            {
                hLEspd_3[i]->Fill(TMath::Abs(lt3->GetLeaf("effEnergy")->GetValue()));
                hZNspd_3[i]->Fill(adcZDCN1 + adcZDCN2);
            }
        }
    }

    for (int i = 0; i < nbins; i++)
    {
        ZNv0m_3[i] = hZNv0m_3[i]->GetMean();
        LEv0m_3[i] = hLEv0m_3[i]->GetMean();
        ErrZNv0m_3[i] = hZNv0m_3[i]->GetMeanError();
        ErrLEv0m_3[i] = hLEv0m_3[i]->GetMeanError();
        ZNspd_3[i] = hZNspd_3[i]->GetMean();
        LEspd_3[i] = hLEspd_3[i]->GetMean();
        ErrZNspd_3[i] = hZNspd_3[i]->GetMeanError();
        ErrLEspd_3[i] = hLEspd_3[i]->GetMeanError();
    }

    TGraphErrors *gZNv0m3 = new TGraphErrors(nbins, LEv0m_3, ZNv0m_3, ErrLEv0m_3, ErrZNv0m_3);
    gZNv0m3->SetMarkerStyle(kFullCircle);
    gZNv0m3->SetMarkerColor(kGreen+1);
    gZNv0m3->SetLineColor(kGreen+1);
    gZNv0m3->SetLineWidth(2);
    gZNv0m3->SetMarkerSize(2.);
    gZNv0m3->SetTitle(";Leading energy (GeV); ZN (a.u.) ");
    histobeauty(gZNv0m3);

    TGraphErrors *gZNspd3 = new TGraphErrors(nbins, LEspd_3, ZNspd_3, ErrLEspd_3, ErrZNspd_3);
    gZNspd3->SetMarkerStyle(kOpenCircle);
    gZNspd3->SetMarkerColor(kGreen+1);
    gZNspd3->SetLineColor(kGreen+1);
    gZNspd3->SetLineWidth(2);
    gZNspd3->SetMarkerSize(2.);
    gZNspd3->SetTitle(";Leading energy (GeV); ZN (a.u.) ");
    histobeauty(gZNspd3);

    /////////////////////////////////////////
    // Get from Tree
    lt4->SetBranchAddress("v0mPerc", &fCentrality_V0M);
    lt4->SetBranchAddress("multSPDcl", &fCentrality_SPDClusters);
    lt4->SetBranchAddress("adcZDCN1", &adcZDCN1);
    lt4->SetBranchAddress("adcZDCN2", &adcZDCN2);
    lt4->SetBranchAddress("effEnergy", &effEnergy);
    lt4->SetBranchAddress("SPDtracklets", &SPDtracklets);

    Float_t ZNv0m_4[nbins];
    Float_t ErrZNv0m_4[nbins];
    Float_t LEv0m_4[nbins];
    Float_t ErrLEv0m_4[nbins];
    Float_t ZNspd_4[nbins];
    Float_t ErrZNspd_4[nbins];
    Float_t LEspd_4[nbins];
    Float_t ErrLEspd_4[nbins];
    TH1F *hLEv0m_4[nbins];
    TH1F *hZNv0m_4[nbins];
    TH1F *hLEspd_4[nbins];
    TH1F *hZNspd_4[nbins];

    for (int i = 0; i < nbins; i++)
    {
        ZNv0m_4[i] = 0;
        ErrZNv0m_4[i] = 0;
        LEv0m_4[i] = 0;
        ErrLEv0m_4[i] = 0;
        ZNspd_4[i] = 0;
        ErrZNspd_4[i] = 0;
        LEspd_4[i] = 0;
        ErrLEspd_4[i] = 0;

        hLEv0m_4[i] = new TH1F(Form("hLEv0m_4%d", i), ";Leading Energy (GeV)", 200, 0, 13000);
        hZNv0m_4[i] = new TH1F(Form("hZNv0m_4%d", i), ";ZN Energy (a.u.)", 200, 0, 13000);
        hLEspd_4[i] = new TH1F(Form("hLEspd_4%d", i), ";Leading Energy (GeV)", 200, 0, 13000);
        hZNspd_4[i] = new TH1F(Form("hZNspd_4%d", i), ";ZN Energy (a.u.)", 200, 0, 13000);
    }

    for (Long_t iEv = 0; iEv < lt4->GetEntries() / 100; iEv++)
    {
        lt4->GetEntry(iEv);
        if (iEv % (lt4->GetEntries() / 10) == 0)
            cout << " At Event " << iEv << " out of " << lt4->GetEntries() << endl;

        if (SPDtracklets < 1)
            continue; // INEL>0

        if (fCentrality_V0M < 0 || fCentrality_V0M >= 100)
            continue;

        for (int i = 0; i < nbins; i++)
        {
            if (fCentrality_V0M >= percentile[i] && fCentrality_V0M <= percentile[i + 1])
            {
                hLEv0m_4[i]->Fill(TMath::Abs(lt4->GetLeaf("effEnergy")->GetValue()));
                hZNv0m_4[i]->Fill(adcZDCN1 + adcZDCN2);
            }
            if (fCentrality_SPDClusters >= percentile[i] && fCentrality_SPDClusters <= percentile[i + 1])
            {
                hLEspd_4[i]->Fill(TMath::Abs(lt4->GetLeaf("effEnergy")->GetValue()));
                hZNspd_4[i]->Fill(adcZDCN1 + adcZDCN2);
            }
        }
    }

    for (int i = 0; i < nbins; i++)
    {
        ZNv0m_4[i] = hZNv0m_4[i]->GetMean();
        LEv0m_4[i] = hLEv0m_4[i]->GetMean();
        ErrZNv0m_4[i] = hZNv0m_4[i]->GetMeanError();
        ErrLEv0m_4[i] = hLEv0m_4[i]->GetMeanError();
        ZNspd_4[i] = hZNspd_4[i]->GetMean();
        LEspd_4[i] = hLEspd_4[i]->GetMean();
        ErrZNspd_4[i] = hZNspd_4[i]->GetMeanError();
        ErrLEspd_4[i] = hLEspd_4[i]->GetMeanError();
    }

    TGraphErrors *gZNv0m4 = new TGraphErrors(nbins, LEv0m_4, ZNv0m_4, ErrLEv0m_4, ErrZNv0m_4);
    gZNv0m4->SetMarkerStyle(kFullCircle);
    gZNv0m4->SetMarkerColor(kBlack);
    gZNv0m4->SetLineColor(kBlack);
    gZNv0m4->SetLineWidth(2);
    gZNv0m4->SetMarkerSize(2.);
    gZNv0m4->SetTitle(";Leading energy (GeV); ZN (a.u.) ");
    histobeauty(gZNv0m4);

    TGraphErrors *gZNspd4 = new TGraphErrors(nbins, LEspd_4, ZNspd_4, ErrLEspd_4, ErrZNspd_4);
    gZNspd4->SetMarkerStyle(kOpenCircle);
    gZNspd4->SetMarkerColor(kBlack);
    gZNspd4->SetLineColor(kBlack);
    gZNspd4->SetLineWidth(2);
    gZNspd4->SetMarkerSize(2.);
    gZNspd4->SetTitle(";Leading energy (GeV); ZN (a.u.) ");
    histobeauty(gZNspd4);

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
    c1->SetBottomMargin(0.15);
    c1->SetRightMargin(0.05);
    c1->SetLeftMargin(0.15);
    c1->SetTopMargin(0.05);
    c1->SetTicks(1, 1);
    hframe->Draw();
    gZNv0m1->Draw("EP SAME");
    gZNspd1->Draw("EP SAME");
    gZNv0m2->Draw("EP SAME");
    gZNspd2->Draw("EP SAME");
    gZNv0m3->Draw("EP SAME");
    gZNspd3->Draw("EP SAME");
    gZNv0m4->Draw("EP SAME");
    gZNspd4->Draw("EP SAME");
    TLine *line = new TLine(0, 0, 12000, 700);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->SetLineColor(kRed);
    //line->Draw("SAME");
    TLine *line2 = new TLine(0, 0, 13000, 650);
    line2->SetLineStyle(2);
    line2->SetLineWidth(1);
    line2->SetLineColor(kBlue);
    //line2->Draw("SAME");
    tex->Draw();
    tex2->Draw();

    TLegend *leg = new TLegend(0.5, 0.2, 0.85, 0.4);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(gZNv0m1, "Pythia8 Monash 2013", "L");
    leg->AddEntry(gZNv0m2, "EPOS", "L");
    leg->AddEntry(gZNv0m3, "Pythia6 Perugia 2011", "L");
    leg->AddEntry(gZNv0m4, "Phojet", "L");
    leg->Draw("SAME");

    TLegend *leg2 = new TLegend(0.2, 0.65, 0.4, 0.78);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.03);
    leg2->AddEntry(gZNv0m4, "VZEROM classes", "P");
    leg2->AddEntry(gZNspd4, "SPDClusters classes", "P");
    leg2->Draw("SAME");

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
