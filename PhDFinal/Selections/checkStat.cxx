Double_t getval(TFile *filecalib, Int_t run, Int_t value);

void checkStat(){

    TFile *file = TFile::Open("~/strg_analysis/AnalisiFinale/data/real/LHCFull_pass2.root");

    TFile *filecontainer = TFile::Open("NchRawContainer_reference.root");
    if (filecontainer == 0x0) {
        cout << "File not found" << endl;
        return;
    }

    TH3D *hnch = (TH3D *)filecontainer->Get("hspd_spdv0m");   // y -> spd, z -> v0m
    TH3D *hZN = (TH3D *)filecontainer->Get("hznsum_spdv0m");  // y -> spd, z -> v0m

    const int nbins = 10;
    TH1D *pnch[nbins][nbins], *pZN[nbins][nbins], *pznMB, *pnchMB;
    double miny, maxy, minz, maxz;
    double meanZN[nbins][nbins], meanErrZN[nbins][nbins], meanNch[nbins][nbins], meanErrNch[nbins][nbins];
    double meannormZN[nbins][nbins], negmeannormZN[nbins][nbins], meannormErrZN[nbins][nbins], meannormNch[nbins][nbins], meannormErrNch[nbins][nbins];
    double meanZNMB, meanErrZNMB, meanNchMB, meanErrNchMB;

    pnchMB = hnch->ProjectionX("pnchMB", hnch->GetYaxis()->FindBin(0. + 1e-10), hnch->GetYaxis()->FindBin(100.), hnch->GetZaxis()->FindBin(0. + 1e-10), hnch->GetZaxis()->FindBin(100.));
    pznMB = hZN->ProjectionX("pznMB", hZN->GetYaxis()->FindBin(0. + 1e-10), hZN->GetYaxis()->FindBin(100.), hZN->GetZaxis()->FindBin(0. + 1e-10), hZN->GetZaxis()->FindBin(100.));
    meanZNMB = pznMB->GetMean();
    meanNchMB = pnchMB->GetMean();

    TH2D *hMult = new TH2D("hMult", "#LT n_{ch}  #GT w.r.t. MB (| #eta | < 0.5); SPDClusters(%); V0M(%) ", nbins, 0, 100, nbins, 0, 100);
    TH2D *hZDC = new TH2D("hZDC", "#LT ZN #GT w.r.t. MB (|#eta|>8);SPDClusters (%);V0M (%)", nbins, 0, 100, nbins, 0, 100);
    TH2D *hStat = new TH2D("hStat", "Candidates;SPDClusters (%);V0M (%)", nbins, 0, 100, nbins, 0, 100);

    // double-differential
    double perc1[nbins + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    double perc2[nbins + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

    for (int c1 = 0; c1 < nbins; c1++)
    { // percentile 1
        for (int c2 = 0; c2 < nbins; c2++)
        { // percentile 2
            miny = hnch->GetYaxis()->FindBin(perc1[c1]);
            maxy = hnch->GetYaxis()->FindBin(perc1[c1 + 1]);
            minz = hnch->GetZaxis()->FindBin(perc2[c2]);
            maxz = hnch->GetZaxis()->FindBin(perc2[c2 + 1]);
            pnch[c1][c2] = hnch->ProjectionX(Form("pnch_%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1[c1], perc1[c1 + 1], "V0M", perc2[c2], perc2[c2 + 1]), miny, maxy, minz, maxz);
            pZN[c1][c2] = hZN->ProjectionX(Form("pZN_%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1[c1], perc1[c1 + 1], "V0M", perc2[c2], perc2[c2 + 1]), miny, maxy, minz, maxz);

            // nch
            meanNch[c1][c2] = pnch[c1][c2]->GetMean();
            meannormNch[c1][c2] = meanNch[c1][c2] / meanNchMB;

            // ZN
            meanZN[c1][c2] = pZN[c1][c2]->GetMean();
            meannormZN[c1][c2] = meanZN[c1][c2] / meanZNMB;
            negmeannormZN[c1][c2] = -meannormZN[c1][c2];

            hMult->SetBinContent(c1 + 1, c2 + 1, meannormNch[c1][c2]);
            hZDC->SetBinContent(c1 + 1, c2 + 1, meannormZN[c1][c2]);
        } // end percentile 2
    }     // end percentile 1

    cout << "--------------- Open Real Data File --------------------" << endl;
    TList *clist = (TList *)file->Get("PWGLF_StrVsMult/cList");
    TTree *lTree = (TTree *)file->Get("PWGLF_StrVsMult/fTreeCascade");

    Float_t fCentrality_V0M = 0.;
    Float_t fCentrality_SPDClusters = 0.;
    Float_t fCentrality_ZDC = 0.;
    Int_t fRun = 0;
    Int_t fNTracksGlobal0815 = 0;
    Int_t fSPDTracklets0815 = 0;
    Int_t fSPDTracklets = 0;
    Float_t fZPApp = 0.;
    Float_t fZPCpp = 0.;
    Float_t fZNApp = 0.;
    Float_t fZNCpp = 0.;
    Float_t lInvariantMass = 0.;
    Float_t lRap = 0.;

    //Double_t perc1_low[] = {0,10,20,30,50};
    //Double_t perc1_high[] = {20,30,40,50,100};
    //Double_t perc2_low[] = {40,30,30,20,0};
    //Double_t perc2_high[] = {60,70,50,50,30};

    //Double_t perc1_low[] = {10,40,60,70,80};
    //Double_t perc1_high[] = {40,60,70,80,100};
    //Double_t perc2_low[] = {70,60,40,40,40};
    //Double_t perc2_high[] = {100,100,100,80,70};

    Double_t perc1_low[] = {0,10,20,30};//{0,30};
    Double_t perc1_high[] = {10,20,30,50};//{10,60};
    Double_t perc2_low[] = {20,10,0,0};//{30,0};
    Double_t perc2_high[] = {30,30,20,10};//{50,20};

    const int size_ = sizeof(perc1_low) / sizeof(Double_t);

    Int_t ncand[size_] = {};

    lTree->SetBranchAddress("fTreeCascVarCentrality_V0M", &fCentrality_V0M);
    lTree->SetBranchAddress("fTreeCascVarCentrality_SPDClusters", &fCentrality_SPDClusters);
    lTree->SetBranchAddress("fTreeCascVarMassAsXi", &lInvariantMass);
    lTree->SetBranchAddress("fTreeCascVarRapXi", &lRap);

    cout << "\n\n--------------------------------------------------------" << endl;
    cout << " Will now loop over events, please wait..." << endl;
    cout << "\n\n--------------------------------------------------------" << endl;

    Float_t thrs = 3000; //minimum candidates to perform the analysis
    for (Long_t iEv = 0; iEv < lTree->GetEntries(); iEv++)
    {
        lTree->GetEntry(iEv);
        if (iEv % (lTree->GetEntries() / 10) == 0)
            cout << " At Event " << iEv << " out of " << lTree->GetEntries() << endl;

        if (TMath::Abs(lRap) < 0.5 && lInvariantMass > 1.305 && lInvariantMass < 1.335)
        {
            hStat->Fill(fCentrality_SPDClusters, fCentrality_V0M);
            for (int c1 = 0; c1 < size_; c1++){
                if (fCentrality_SPDClusters > perc1_low[c1] && fCentrality_SPDClusters <= perc1_high[c1] && fCentrality_V0M > perc2_low[c1] && fCentrality_V0M <= perc2_high[c1]){
                    ncand[c1]++;
                }
            }
        }
    }

    TCanvas *c2 = new TCanvas("c2", "c2", 2100, 700);
    c2->Divide(3, 1);
    c2->cd(1);
    hMult->SetStats(0);
    hMult->Draw("colz, text");
    c2->cd(2);
    hZDC->SetStats(0);
    hZDC->Draw("colz, text");
    c2->cd(3);
    hStat->SetStats(0);
    hStat->Draw("colz, text");
    c2->SaveAs("images/statcheck.png");

    TGraph *gtot[10][10];
    Int_t colors[] = {kRed + 3, kRed + 1, kOrange + 10, kOrange + 1, kYellow + 1, kGreen + 2, kSpring - 1, kAzure + 7, kBlue - 4, kBlue + 2};
    Int_t marker[] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullStar, kFullCrossX, kFullCross, 45, 39};
    Int_t ntot = 0;
    float xtot[1000], ytot[1000], xrev[10][10], yrev[10][10];

    for (int c1 = 0; c1 < nbins; c1++)
    {
        for (int c2 = 0; c2 < nbins; c2++)
        {
            gtot[c1][c2] = new TGraph(1, &meannormNch[c1][c2], &meannormZN[c1][c2]);
            gtot[c1][c2]->SetMarkerColor(colors[c1]);
            gtot[c1][c2]->SetMarkerStyle(marker[c2]);
            gtot[c1][c2]->SetMarkerColor(colors[c1]);
            gtot[c1][c2]->SetMarkerSize(2.);

            xtot[ntot] = meannormNch[c1][c2];
            ytot[ntot] = meannormZN[c1][c2];

            xrev[c2][c1] = meannormNch[c1][c2];
            yrev[c2][c1] = meannormZN[c1][c2];

            ntot++;
        }
    }

    TCanvas *c3 = new TCanvas("c3", "", 1500, 900);
    TGraph *g = new TGraph(ntot, xtot, ytot);
    g->GetXaxis()->SetTitle("#LT n_{ch} #GT / #LT n_{ch} #GT_{MB}");
    g->GetYaxis()->SetTitle("#LT ZN #GT / #LT ZN #GT_{MB}");
    g->GetXaxis()->SetTitleSize(0.04);
    g->GetYaxis()->SetTitleSize(0.04);
    g->SetTitle("");

    TH1D *hh = new TH1D("hh", ";n_{ch} w.r.t. to MB (|#eta|<0.5);ZN w.r.t. MB (|#eta>8|)", 10, 0.2, 3.5);
    hh->SetStats(0);
    hh->GetYaxis()->SetRangeUser(0.2, 1.5);

    hh->Draw();
    g->Draw("P SAME");
    TGraph *g1[10];
    TGraph *g2[10];
    for (int i = 0; i < 10; i++)
    {
        g1[i] = new TGraph(10, meannormNch[i], meannormZN[i]);
        g2[i] = new TGraph(10, xrev[i], yrev[i]);
        g1[i]->SetLineColor(colors[i]); // linea verticale
        g1[i]->SetMarkerStyle(0.);
        g2[i]->SetLineColor(kBlack); // linea orizzontale
        g2[i]->SetMarkerStyle(0);
    }

    for (int i = 0; i < 10; i++)
    {
        g2[i]->Draw("LEP SAME");
        g1[i]->Draw("LEP SAME");
    }

    for (int c1 = 0; c1 < nbins; c1++)
    {
        for (int c2 = 0; c2 < nbins; c2++)
        {
            gtot[c1][c2]->Draw("P SAME");
        }
    }

    // Standalone classe
    // spd
    Double_t percstd1_low[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t percstd1_high[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    // v0m
    Double_t percstd2_low[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
    Double_t percstd2_high[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    const int size_std = sizeof(percstd1_low) / sizeof(Double_t);

    TH1D *pnch_std[size_std], *pZN_std[size_std], *pspd0815_std[size_std];
    double meanZN_std[nbins], meanNch_std[nbins];

    for (int c1 = 0; c1 < size_std; c1++)
    {
        miny = hnch->GetYaxis()->FindBin(percstd1_low[c1]);
        maxy = hnch->GetYaxis()->FindBin(percstd1_high[c1] - 1E-10);
        minz = hnch->GetZaxis()->FindBin(percstd2_low[c1]);
        maxz = hnch->GetZaxis()->FindBin(percstd2_high[c1] - 1E-10);
        pnch_std[c1] = hnch->ProjectionX(Form("pnch_std%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percstd1_low[c1], percstd1_high[c1], "V0M", percstd2_low[c1], percstd2_high[c1]), miny, maxy, minz, maxz);
        pZN_std[c1] = hZN->ProjectionX(Form("pZN_std%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percstd1_low[c1], percstd1_high[c1], "V0M", percstd2_low[c1], percstd2_high[c1]), miny, maxy, minz, maxz);

        // nch
        meanNch_std[c1] = pnch_std[c1]->GetMean() / meanNchMB;

        // ZN
        meanZN_std[c1] = pZN_std[c1]->GetMean() / meanZNMB;

        cout << "v0m: " << percstd2_low[c1] << "-" << percstd2_high[c1] << "  " << meanNch_std[c1] << endl;
    }
    TGraph *g10[size_std];
    for (int i = 0; i < size_std; i++)
    {
        g10[i] = new TGraph(1, &meanNch_std[i], &meanZN_std[i]);
        g10[i]->SetMarkerStyle(kFullCircle);
        g10[i]->SetMarkerSize(2.3);
        g10[i]->SetMarkerColor(kBlack);
        g10[i]->Draw("P SAME");
    }

    TLegend *leg = new TLegend(0.7, 0.5, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.03);
    TMarker *m0[nbins];
    for (int i = 0; i < nbins; i++)
    {
        m0[i] = new TMarker(0, 0, marker[i]);
        m0[i]->SetMarkerSize(2.);
        m0[i]->SetMarkerColor(kBlack);
        leg->AddEntry(m0[i], Form("%d-%d%% V0M", i * 10, (i + 1) * 10), "p");
    }
    leg->Draw();

    TLegend *leg2 = new TLegend(0.5, 0.5, 0.65, 0.85);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->SetTextSize(0.03);
    TMarker *m[nbins];
    for (int i = 0; i < nbins; i++)
    {
        m[i] = new TMarker(0, 0, kFullCircle);
        m[i]->SetMarkerSize(2.);
        m[i]->SetMarkerColor(colors[i]);
        leg2->AddEntry(m[i], Form("%d-%d%% SPDClusters", i * 10, (i + 1) * 10), "p");
    }
    leg2->Draw();

    c3->SaveAs("images/ddsel2D.png");

    TH1D *pnch_[size_], *pZN_[size_];
    double meanZN_[nbins], meanNch_[nbins];

    for (int c1 = 0; c1 < size_; c1++)
    {
        cout << "ncand[" << c1 << "] = " << ncand[c1] << endl;
        miny = hnch->GetYaxis()->FindBin(perc1_low[c1]);
        maxy = hnch->GetYaxis()->FindBin(perc1_high[c1]);
        minz = hnch->GetZaxis()->FindBin(perc2_low[c1]);
        maxz = hnch->GetZaxis()->FindBin(perc2_high[c1]);
        pnch_[c1] = hnch->ProjectionX(Form("pnch_%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low[c1], perc1_high[c1], "V0M", perc2_low[c1], perc2_high[c1]), miny, maxy, minz, maxz);
        pZN_[c1] = hZN->ProjectionX(Form("pZN_%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low[c1], perc1_high[c1], "V0M", perc2_low[c1], perc2_high[c1]), miny, maxy, minz, maxz);

        // nch
        meanNch_[c1] = pnch_[c1]->GetMean() / meanNchMB;
        // ZN
        meanZN_[c1] = pZN_[c1]->GetMean() / meanZNMB;

    } // end percentile 1

    Int_t colors2[] = {kRed + 3, kOrange + 1, kSpring - 1, kAzure + 7, kBlue + 2};

    new TCanvas;
    for (int c1 = 0; c1 < size_; c1++){
        pnch_[c1]->Scale(1. / pnch_[c1]->Integral());
        pnch_[c1]->SetLineColor(colors2[c1]);
        pnch_[c1]->SetLineWidth(2);
        //pnch_[c1]->Divide(pnch_[0]);
        pnch_[c1]->Draw("SAME HIST");
    }

    new TCanvas;
    pZN_[0]->Scale(1. / pZN_[0]->Integral());
    pZN_[0]->Rebin(10);

    for (int c1 = 1; c1 < size_; c1++)
    {
        pZN_[c1]->Scale(1. / pZN_[c1]->Integral());
        pZN_[c1]->Rebin(10);
        pZN_[c1]->SetLineColor(colors2[c1]);
        pZN_[c1]->SetLineWidth(2);
        pZN_[c1]->SetLineStyle(1+c1);
        pZN_[c1]->Divide(pZN_[0]);
        pZN_[c1]->Draw("SAME HIST");
    }

    TGraph *g4[size_];

    TCanvas *c10 = new TCanvas();
    TH1D *h = new TH1D("h", "", 10, 0, 4);
    h->GetYaxis()->SetRangeUser(0., 1.6);
    h->Draw();
    for (int i = 0; i < size_; i++)
    {
        g4[i] = new TGraph(1, &meanNch_[i], &meanZN_[i]);
        g4[i]->SetMarkerStyle(kFullCircle);
        g4[i]->SetMarkerSize(2.3);
        g4[i]->SetMarkerColor(colors2[i]);
        g4[i]->Draw("P SAME");
    }

    TLegend *leg3 = new TLegend(0.5, 0.7, 0.65, 0.85);
    leg3->SetBorderSize(0);
    leg3->SetFillColor(0);
    leg3->SetTextSize(0.025);
    TMarker *m2[size_];
    for (int i = 0; i < size_; i++)
    {
        m2[i] = new TMarker(0, 0, kFullCircle);
        m2[i]->SetMarkerSize(2.);
        m2[i]->SetMarkerColor(colors2[i]);
        leg3->AddEntry(m2[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low[i], (int)perc1_high[i], (int)perc2_low[i], (int)perc2_high[i]), "p");
    }
    leg3->Draw();
}

Double_t getval(TFile *filecalib, Int_t run, Int_t value)
{
    // Return the percentile given the value of estimator
    TString name = Form("hcumSPDClusters_%i", run);

    TH1D *hcum = (TH1D *)filecalib->Get(name);
    Double_t percentile = 100 * (hcum->GetBinContent(hcum->FindBin(value)));

    return percentile;
}
