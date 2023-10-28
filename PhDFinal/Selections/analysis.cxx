Double_t getval(TFile *filecalib, Int_t run, Int_t value);

void analysis(){

    TFile *filecontainer = TFile::Open("NchRawContainer_reference.root");
    if (filecontainer == 0x0)
    {
        cout << "File not found" << endl;
        return;
    }

    TH3D *hnch = (TH3D *)filecontainer->Get("hspd_spdv0m");  // y -> spd, z -> v0m
    TH3D *hZN = (TH3D *)filecontainer->Get("hznsum_spdv0m"); // y -> spd, z -> v0m

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

    TH2D* hMult = new TH2D("hMult", "#LT n_{ch}  #GT w.r.t. MB (| #eta | < 0.5); SPDClusters(%); V0M(%) ", nbins, 0, 100, nbins, 0, 100);
    TH2D *hZDC = new TH2D("hZDC", "#LT ZN #GT w.r.t. MB (|#eta|>8);SPDClusters (%);V0M (%)", nbins, 0, 100, nbins, 0, 100);
    TH2D *hStat = new TH2D("hStat", "Candidates;SPDClusters (%);V0M (%)", nbins, 0, 100, nbins, 0, 100);

    // double-differential
    double perc1[nbins + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    double perc2[nbins + 1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

    for (int c1 = 0; c1 < nbins; c1++)
    { // percentile 1
        for (int c2 = 0; c2 < nbins; c2++)
        { // percentile 2
            miny = hnch->GetYaxis()->FindBin(perc1[c1] + 1e-10);
            maxy = hnch->GetYaxis()->FindBin(perc1[c1 + 1]);
            minz = hnch->GetZaxis()->FindBin(perc2[c2] + 1e-10);
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
    } // end percentile 1

    Double_t perc1_low_1[] = {10, 10, 10, 10, 10, 10, 10};
    Double_t perc1_high_1[] = {20, 20, 20, 20, 20, 20, 20};
    Double_t perc2_low_1[] = {0, 5, 10, 20, 30, 40, 50};
    Double_t perc2_high_1[] = {5, 10, 20, 30, 40, 50, 100};
    const int size1 = sizeof(perc1_low_1) / sizeof(Double_t);

    Double_t perc1_low_2[] = {40, 40, 40, 40, 40, 40, 40};
    Double_t perc1_high_2[] = {50, 50, 50, 50, 50, 50, 50};
    Double_t perc2_low_2[] = {0, 20, 30, 40, 50, 60, 70};
    Double_t perc2_high_2[] = {20, 30, 40, 50, 60, 70, 100};
    const int size2 = sizeof(perc1_low_2) / sizeof(Double_t);

    Double_t perc1_low_3[] = {0, 10, 20, 30, 50};
    Double_t perc1_high_3[] = {20, 30, 40, 50, 100};
    Double_t perc2_low_3[] = {40, 30, 30, 20, 0};
    Double_t perc2_high_3[] = {60, 70, 50, 50, 30};
    const int size3 = sizeof(perc1_low_3) / sizeof(Double_t);

    Double_t perc1_low_4[] = {0, 10, 20, 30};
    Double_t perc1_high_4[] = {10, 20, 30, 50};
    Double_t perc2_low_4[] = {20, 10, 0, 0};
    Double_t perc2_high_4[] = {30, 30, 20, 10};
    const int size4 = sizeof(perc1_low_4) / sizeof(Double_t);

    Double_t percstd1_low[] = {0,0,0,0,0,0,0,0,0,0};
    Double_t percstd1_high[] = {100,100,100,100,100,100,100,100,100,100};
    Double_t percstd2_low[] = {0,1,5,10,15,20,30,40,50,70};
    Double_t percstd2_high[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    const int size_std = sizeof(percstd1_low) / sizeof(Double_t);

    TGraph *gtot[10][10];
    Int_t colors[] = {kRed + 3, kRed + 1, kOrange + 10, kOrange + 1, kYellow + 1, kGreen + 2, kSpring - 1, kAzure + 7, kBlue - 4, kBlue + 2};
    Int_t marker[] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullStar, kFullCrossX, kFullCross, 45, 39};
    Int_t ntot = 0;
    float xtot[1020], ytot[1020],xrev[10][10],yrev[10][10];
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

    TH1D *hh = new TH1D("hh", ";#LT n_{ch} #GT / #LT n_{ch} #GT_{MB};#LT ZN #GT / #LT ZN #GT_{MB}", 10, 0.2, 3.5);
    hh->SetStats(0);
    hh->GetYaxis()->SetRangeUser(.1, 1.5);
    hh->GetYaxis()->SetTitleSize(0.035);
    hh->GetXaxis()->SetTitleSize(0.035);

    TH1D *pnch_std[size_std], *pZN_std[size_std], *pspd0815_std[size_std];
    double meanZN_std[nbins], negmeanZN_std[nbins], meanNch_std[nbins];

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
        negmeanZN_std[c1] = -meanZN_std[c1];
    }

    TGraph *g0[size_std];
    for (int i = 0; i < size_std; i++)
    {
        g0[i] = new TGraph(1, &meanNch_std[i], &meanZN_std[i]);
        g0[i]->SetMarkerStyle(kFullDiamond);
        g0[i]->SetMarkerSize(3.3);
        g0[i]->SetMarkerColor(kBlack);
        g0[i]->Draw("P SAME");
    }

    TCanvas *c = new TCanvas("c", "", 1000, 900);
    c->SetBottomMargin(0.12);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.02);
    c->SetTopMargin(0.02);
    c->SetTicky();
    c->SetTickx();

    hh->Draw();

    TH1D *pnch_1[size1], *pZN_1[size1], *pspd0815_1[size1];
    double meanZN_1[nbins], meanNch_1[nbins];

    TH1D *pnch_2[size2], *pZN_2[size2], *pspd0815_2[size2];
    double meanZN_2[nbins], meanNch_2[nbins];

    TH1D *pnch_3[size3], *pZN_3[size3], *pspd0815_3[size3];
    double meanZN_3[nbins], meanNch_3[nbins];

    TH1D *pnch_4[size4], *pZN_4[size4], *pspd0815_4[size4];
    double meanZN_4[nbins], meanNch_4[nbins];

    for (int c1 = 0; c1 < size1; c1++)
    {
        miny = hnch->GetYaxis()->FindBin(perc1_low_1[c1] );
        maxy = hnch->GetYaxis()->FindBin(perc1_high_1[c1]-1E-10);
        minz = hnch->GetZaxis()->FindBin(perc2_low_1[c1] );
        maxz = hnch->GetZaxis()->FindBin(perc2_high_1[c1]-1E-10);
        pnch_1[c1] = hnch->ProjectionX(Form("pnch_1%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_1[c1], perc1_high_1[c1], "V0M", perc2_low_1[c1], perc2_high_1[c1]), miny, maxy, minz, maxz);
        pZN_1[c1] = hZN->ProjectionX(Form("pZN_1%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_1[c1], perc1_high_1[c1], "V0M", perc2_low_1[c1], perc2_high_1[c1]), miny, maxy, minz, maxz);

        // nch
        meanNch_1[c1] = pnch_1[c1]->GetMean() / meanNchMB;

        // ZN
        meanZN_1[c1] = pZN_1[c1]->GetMean() / meanZNMB;

    }
    for (int c1 = 0; c1 < size2; c1++)
    {
        miny = hnch->GetYaxis()->FindBin(perc1_low_2[c1]);
        maxy = hnch->GetYaxis()->FindBin(perc1_high_2[c1] - 1E-10);
        minz = hnch->GetZaxis()->FindBin(perc2_low_2[c1]);
        maxz = hnch->GetZaxis()->FindBin(perc2_high_2[c1] - 1E-10);
        pnch_2[c1] = hnch->ProjectionX(Form("pnch_2%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_2[c1], perc1_high_2[c1], "V0M", perc2_low_2[c1], perc2_high_2[c1]), miny, maxy, minz, maxz);
        pZN_2[c1] = hZN->ProjectionX(Form("pZN_2%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_2[c1], perc1_high_2[c1], "V0M", perc2_low_2[c1], perc2_high_2[c1]), miny, maxy, minz, maxz);

        // nch
        meanNch_2[c1] = pnch_2[c1]->GetMean() / meanNchMB;

        // ZN
        meanZN_2[c1] = pZN_2[c1]->GetMean() / meanZNMB;
    }
    for (int c1 = 0; c1 < size3; c1++)
    {
        miny = hnch->GetYaxis()->FindBin(perc1_low_3[c1] );
        maxy = hnch->GetYaxis()->FindBin(perc1_high_3[c1]-1E-10);
        minz = hnch->GetZaxis()->FindBin(perc2_low_3[c1] );
        maxz = hnch->GetZaxis()->FindBin(perc2_high_3[c1]-1E-10);
        pnch_3[c1] = hnch->ProjectionX(Form("pnch_3%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_3[c1], perc1_high_3[c1], "V0M", perc2_low_3[c1], perc2_high_3[c1]), miny, maxy, minz, maxz);
        pZN_3[c1] = hZN->ProjectionX(Form("pZN_3%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_3[c1], perc1_high_3[c1], "V0M", perc2_low_3[c1], perc2_high_3[c1]), miny, maxy, minz, maxz);

        // nch
        meanNch_3[c1] = pnch_3[c1]->GetMean() / meanNchMB;

        // ZN
        meanZN_3[c1] = pZN_3[c1]->GetMean() / meanZNMB;

    }
    for (int c1 = 0; c1 < size4; c1++)
    {
        miny = hnch->GetYaxis()->FindBin(perc1_low_4[c1]);
        maxy = hnch->GetYaxis()->FindBin(perc1_high_4[c1] - 1E-10);
        minz = hnch->GetZaxis()->FindBin(perc2_low_4[c1]);
        maxz = hnch->GetZaxis()->FindBin(perc2_high_4[c1] - 1E-10);
        pnch_4[c1] = hnch->ProjectionX(Form("pnch_4%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_4[c1], perc1_high_4[c1], "V0M", perc2_low_4[c1], perc2_high_4[c1]), miny, maxy, minz, maxz);
        pZN_4[c1] = hZN->ProjectionX(Form("pZN_4%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_4[c1], perc1_high_4[c1], "V0M", perc2_low_4[c1], perc2_high_4[c1]), miny, maxy, minz, maxz);

        // nch
        meanNch_4[c1] = pnch_4[c1]->GetMean() / meanNchMB;

        // ZN
        meanZN_4[c1] = pZN_4[c1]->GetMean() / meanZNMB;
    } // end percentile

    TGraph *g1[size1];
    TGraph *g2[size2];
    TGraph *g3[size3];
    TGraph *g4[size4];

    Int_t colors1[] = {kRed + 2, kOrange + 1, kSpring - 1, kAzure + 7, kBlue + 2};
    Int_t colors2[] = {kRed + 2, kOrange + 1, kYellow+1, kSpring - 1, kGreen+2, kAzure + 7, kBlue + 2};

    for (int i = 0; i < size1; i++)
    {
        g1[i] = new TGraph(1, &meanNch_1[i], &meanZN_1[i]);
        g1[i]->SetMarkerStyle(kFullCircle);
        g1[i]->SetMarkerSize(2.8);
        g1[i]->SetMarkerColor(colors2[i]);
        g1[i]->Draw("P SAME");
    }
    for (int i = 0; i < size2; i++)
    {
        g2[i] = new TGraph(1, &meanNch_2[i], &meanZN_2[i]);
        g2[i]->SetMarkerStyle(kOpenCircle);
        g2[i]->SetMarkerSize(2.8);
        g2[i]->SetMarkerColor(colors2[i]);
        g2[i]->Draw("P SAME");
    }
    for (int i = 0; i < size3; i++)
    {
        g3[i] = new TGraph(1, &meanNch_3[i], &meanZN_3[i]);
        g3[i]->SetMarkerStyle(kOpenSquare);
        g3[i]->SetMarkerSize(2.8);
        g3[i]->SetMarkerColor(colors1[i]);
        g3[i]->Draw("P SAME");
    }
    for (int i = 0; i < size4; i++)
    {
        g4[i] = new TGraph(1, &meanNch_4[i], &meanZN_4[i]);
        g4[i]->SetMarkerStyle(kFullSquare);
        g4[i]->SetMarkerSize(2.8);
        g4[i]->SetMarkerColor(colors1[i]);
        g4[i]->Draw("P SAME");
    }
    for (int i = 0; i < size_std; i++)
    {
        g0[i]->Draw("P SAME");
    }

    TLegend *leg1 = new TLegend(0.55, 0.8, 0.8, 0.95);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetTextSize(0.025);
    TLegend *leg2 = new TLegend(0.15, 0.38, 0.35, 0.52);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->SetTextSize(0.025);
    TLegend *leg3 = new TLegend(0.15, 0.15, 0.35, 0.35);
    leg3->SetBorderSize(0);
    leg3->SetFillColor(0);
    leg3->SetTextSize(0.025);
    TLegend *leg4 = new TLegend(0.55, 0.65, 0.8, 0.78);
    leg4->SetBorderSize(0);
    leg4->SetFillColor(0);
    leg4->SetTextSize(0.025);
    TMarker *m1[size1], *m2[size2], *m3[size3], *m4[size4];
    for (int i = 0; i < size1; i++){
        m1[i] = new TMarker(0, 0, kFullCircle);
        m1[i]->SetMarkerSize(2.2);
        m1[i]->SetMarkerColor(colors2[i]);
        leg1->AddEntry(m1[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_1[i], (int)perc1_high_1[i], (int)perc2_low_1[i], (int)perc2_high_1[i]), "p");
    }
    for (int i = 0; i < size2; i++)
    {
        m2[i] = new TMarker(0, 0, kOpenCircle);
        m2[i]->SetMarkerSize(2.2);
        m2[i]->SetMarkerColor(colors2[i]);
        leg2->AddEntry(m2[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_2[i], (int)perc1_high_2[i], (int)perc2_low_2[i], (int)perc2_high_2[i]), "p");
    }
    for (int i = 0; i < size3; i++){
        m3[i] = new TMarker(0, 0, kOpenSquare);
        m3[i]->SetMarkerSize(2.2);
        m3[i]->SetMarkerColor(colors1[i]);
        leg3->AddEntry(m3[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_3[i], (int)perc1_high_3[i], (int)perc2_low_3[i], (int)perc2_high_3[i]), "p");
    }
    for (int i = 0; i < size4; i++){
        m4[i] = new TMarker(0, 0, kFullSquare);
        m4[i]->SetMarkerSize(2.2);
        m4[i]->SetMarkerColor(colors1[i]);
        leg4->AddEntry(m4[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_4[i], (int)perc1_high_4[i], (int)perc2_low_4[i], (int)perc2_high_4[i]), "p");
    }
    leg2->Draw();
    leg1->Draw();
    leg3->Draw();
    leg4->Draw();

    c->SaveAs(Form("images/tutteleclassi.png"));

    TCanvas *c_0 = new TCanvas("c_0", "", 1000, 900);
    c_0->SetBottomMargin(0.12);
    c_0->SetLeftMargin(0.12);
    c_0->SetRightMargin(0.02);
    c_0->SetTopMargin(0.02);
    c_0->SetTicky();
    c_0->SetTickx();
    hh->Draw();
    for (int i = 0; i < size_std; i++)
    {
        g0[i]->Draw("P SAME");
    }

    TCanvas *c_1 = new TCanvas("c_1", "", 1000, 900);
    c_1->SetBottomMargin(0.12);
    c_1->SetLeftMargin(0.12);
    c_1->SetRightMargin(0.02);
    c_1->SetTopMargin(0.02);
    c_1->SetTicky();
    c_1->SetTickx();
    hh->Draw();
    for (int i = 0; i < size4; i++)
    {
        g4[i]->Draw("P SAME");
    }
    for (int i = 0; i < size_std; i++)
    {
        g0[i]->Draw("P SAME");
    }

    TCanvas *c_2 = new TCanvas("c_2", "", 1000, 900);
    c_2->SetBottomMargin(0.12);
    c_2->SetLeftMargin(0.12);
    c_2->SetRightMargin(0.02);
    c_2->SetTopMargin(0.02);
    c_2->SetTicky();
    c_2->SetTickx();
    hh->Draw();
    for (int i = 0; i < size_std; i++)
    {
        g0[i]->Draw("P SAME");
    }
    for (int i = 0; i < size3; i++)
    {
        g3[i]->Draw("P SAME");
    }
    for (int i = 0; i < size4; i++)
    {
        g4[i]->Draw("P SAME");
    }
    c_0->SaveAs(Form("images/tutteleclassi_0.png"));
    c_1->SaveAs(Form("images/tutteleclassi_04.png"));
    c_2->SaveAs(Form("images/tutteleclassi_043.png"));

    /*for (int i = 0; i < size1; i++)
    {
        g1[i]->Draw("P SAME");
    }
    for (int i = 0; i < size2; i++)
    {
        g2[i]->Draw("P SAME");
    }
    for (int i = 0; i < size3; i++)
    {
        g3[i]->Draw("P SAME");
    }
    for (int i = 0; i < size4; i++)
    {
        g4[i]->Draw("P SAME");
    }
    for (int i = 0; i < size_std; i++)
    {
        g0[i]->Draw("P SAME");
    }*/

    //leg2->Draw();
    //leg1->Draw();
    //leg3->Draw();
    //leg4->Draw();



    TCanvas *c5 = new TCanvas("c5", "class1", 1300, 800);
    c5->Divide(2, 1);

    c5->cd(1);
    // c5->cd(1)->SetLogy();
    TPad *pad1 = new TPad("pad1", "", 0, 0.3, 1, 1);
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    pad1->SetBottomMargin(0.001);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.02);
    pad1->SetTopMargin(0.02);
    pad1->Draw();
    pad2->SetBottomMargin(0.3);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.02);
    pad2->SetTopMargin(0.001);
    pad2->Draw();

    pad1->cd();
    TH1D *rationch[size1];
    for (int i = 0; i < size1; i++)
    {
        pnch_1[i]->SetStats(0);
        pnch_1[i]->SetTitle("");
        pnch_1[i]->Scale(1. / pnch_1[i]->Integral(0, -1));
        pnch_1[i]->GetYaxis()->SetRangeUser(0, pnch_1[i]->GetMaximum() * 1.2);
        pnch_1[i]->GetXaxis()->SetRangeUser(0, 30);
        pnch_1[i]->SetLineColor(colors2[i]);
        pnch_1[i]->SetLineWidth(2);
        pnch_1[i]->SetMarkerColor(colors2[i]);
        pnch_1[i]->SetMarkerStyle(kFullCircle);
        pnch_1[i]->Draw("LEP SAME");
        pnch_1[i]->Draw("HISTs SAME");
    }
    TLegend *leg5 = new TLegend(0.5, 0.7, 0.65, 0.9);
    leg5->SetBorderSize(0);
    leg5->SetFillColor(0);
    leg5->SetTextSize(0.03);
    for (int i = 0; i < size1; i++)
    {
        leg5->AddEntry(pnch_1[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_1[i], (int)perc1_high_1[i], (int)perc2_low_1[i], (int)perc2_high_1[i]), "LEP");
        //leg5->AddEntry(pnch_1[i], Form("#LT n_{ch} #GT = %.2f", meanNch_1[i] * meanNchMB), "LEP");
    }
    leg5->Draw();

    pad2->cd();
    //pad2->SetLogy();
    for (int i = 0; i < size1; i++)
    {
        rationch[i] = (TH1D *)pnch_1[i]->Clone(Form("rationch%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_1[i], perc1_high_1[i], "V0M", perc2_low_1[i], perc2_high_1[i]));
        rationch[i]->Divide(pnch_1[2]);
        rationch[i]->SetLineColor(colors2[i]);
        rationch[i]->SetLineWidth(2);
        rationch[i]->GetYaxis()->SetRangeUser(0.55, 1.45);
        rationch[i]->GetXaxis()->SetRangeUser(0, 30);
        rationch[i]->GetYaxis()->SetTitle("Ratio to III");
        rationch[i]->GetXaxis()->SetTitleSize(0.1);
        rationch[i]->GetYaxis()->SetTitleSize(0.08);
        rationch[i]->GetXaxis()->SetLabelSize(0.09);
        rationch[i]->GetYaxis()->SetLabelSize(0.09);
        rationch[i]->GetYaxis()->SetTitleOffset(0.7);

        rationch[i]->Draw("LEP SAME");
        rationch[i]->Draw("HIST SAME");
    }

    c5->cd(2);
    //c5->cd(2)->SetLogy();

    TPad *pad1_ = new TPad("pad1_", "", 0, 0.3, 1, 1);
    TPad *pad2_ = new TPad("pad2_", "", 0, 0, 1, 0.3);
    pad1_->SetBottomMargin(0.001);
    pad1_->SetLeftMargin(0.12);
    pad1_->SetRightMargin(0.02);
    pad1_->SetTopMargin(0.02);
    pad1_->Draw();
    pad2_->SetBottomMargin(0.3);
    pad2_->SetLeftMargin(0.12);
    pad2_->SetRightMargin(0.02);
    pad2_->SetTopMargin(0.001);
    pad2_->Draw();
    pad1_->Draw();
    pad1_->cd();
    pad1_->SetLogy();
    TH1D *ratiozn[size1];
    for (int i = 0; i < size1; i++)
    {
        pZN_1[i]->SetStats(0);
        pZN_1[i]->SetTitle("");
        pZN_1[i]->Rebin(10);
        pZN_1[i]->Scale(1. / pZN_1[i]->Integral(0, -1));
        pZN_1[i]->GetYaxis()->SetRangeUser(2E-4, 0.2);
        pZN_1[i]->GetXaxis()->SetRangeUser(0.0, 1020);
        pZN_1[i]->SetLineColor(colors2[i]);
        pZN_1[i]->SetMarkerColor(colors2[i]);
        pZN_1[i]->SetMarkerStyle(kFullCircle);
        pZN_1[i]->SetLineWidth(2);
        pZN_1[i]->Draw("LEP SAME");
        pZN_1[i]->Draw("HIST SAME");
    }
    TLegend *leg4bis = new TLegend(0.4, 0.7, 0.65, 0.9);
    leg4bis->SetBorderSize(0);
    leg4bis->SetFillColor(0);
    leg4bis->SetTextSize(0.03);
    for (int i = 0; i < size1; i++)
    {
        leg4bis->AddEntry(pnch_1[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_1[i], (int)perc1_high_1[i], (int)perc2_low_1[i], (int)perc2_high_1[i]), "LEP");
        //leg4bis->AddEntry(pZN_1[i], Form("#LT ZN #GT (a.u.) = %.2f", TMath::Abs(meanZN_1[i]) * meanZNMB), "LEP");
    }
    leg4bis->Draw();

    pad2_->cd();
    pad2_->SetLogy();
    for (int i = 0; i < size1; i++)
    {
        ratiozn[i] = (TH1D *)pZN_1[i]->Clone(Form("ratiozn%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_1[i], perc1_high_1[i], "V0M", perc2_low_1[i], perc2_high_1[i]));
        ratiozn[i]->Divide(pZN_1[2]);
        ratiozn[i]->SetLineColor(colors2[i]);
        ratiozn[i]->SetLineWidth(2);
        ratiozn[i]->GetYaxis()->SetRangeUser(0.2, 4.);
        ratiozn[i]->GetXaxis()->SetRangeUser(0, 1020);
        ratiozn[i]->GetYaxis()->SetTitle("Ratio to III");
        ratiozn[i]->GetXaxis()->SetTitleSize(0.1);
        ratiozn[i]->GetYaxis()->SetTitleSize(0.08);
        ratiozn[i]->GetXaxis()->SetLabelSize(0.09);
        ratiozn[i]->GetYaxis()->SetLabelSize(0.09);
        ratiozn[i]->GetYaxis()->SetTitleOffset(0.7);

        ratiozn[i]->Draw("LEP SAME");
        ratiozn[i]->Draw("HIST SAME");
    }
    c5->SaveAs(Form("images/distributions%i.png", 1));



    TCanvas *c6 = new TCanvas("c6", "class2", 1300, 800);
    c6->Divide(2, 1);
    c6->cd(1);

    TPad *pad1_2 = new TPad("pad1_2", "", 0, 0.3, 1, 1);
    TPad *pad2_2 = new TPad("pad2_2", "", 0, 0, 1, 0.3);
    pad1_2->SetBottomMargin(0.001);
    pad1_2->SetLeftMargin(0.12);
    pad1_2->SetRightMargin(0.02);
    pad1_2->SetTopMargin(0.02);
    pad1_2->Draw();
    pad2_2->SetBottomMargin(0.3);
    pad2_2->SetLeftMargin(0.12);
    pad2_2->SetRightMargin(0.02);
    pad2_2->SetTopMargin(0.001);
    pad2_2->Draw();

    pad1_2->cd();
    TH1D *rationch_2[size2];
    // c6->cd(1)->SetLogy();
    for (int i = 0; i < size2; i++)
    {
        pnch_2[i]->SetStats(0);
        pnch_2[i]->SetTitle("");
        pnch_2[i]->Scale(1. / pnch_2[i]->Integral(0, -1));
        pnch_2[i]->GetYaxis()->SetRangeUser(0., 0.2);
        pnch_2[i]->GetXaxis()->SetRangeUser(0, 20);
        pnch_2[i]->SetLineColor(colors2[i]);
        pnch_2[i]->SetLineWidth(2);
        pnch_2[i]->SetMarkerColor(colors2[i]);
        pnch_2[i]->SetMarkerStyle(kFullCircle);
        pnch_2[i]->Draw("LEP SAME");
        pnch_2[i]->Draw("HIST SAME");
    }
    TLegend *leg6 = new TLegend(0.4, 0.7, 0.65, 0.9);
    leg6->SetBorderSize(0);
    leg6->SetFillColor(0);
    leg6->SetTextSize(0.03);
    for (int i = 0; i < size2; i++)
    {
        leg6->AddEntry(pnch_2[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_2[i], (int)perc1_high_2[i], (int)perc2_low_2[i], (int)perc2_high_2[i]), "LEP");
        //leg6->AddEntry(pnch_2[i], Form("#LT n_{ch} #GT = %.2f", meanNch_2[i] * meanNchMB), "LEP");
    }

    leg6->Draw();

    pad2_2->cd();
    // pad2_2->SetLogy();
    for (int i = 0; i < size2; i++)
    {
        rationch_2[i] = (TH1D *)pnch_2[i]->Clone(Form("rationch_2%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_2[i], perc1_high_2[i], "V0M", perc2_low_2[i], perc2_high_2[i]));
        rationch_2[i]->Divide(pnch_2[2]);
        rationch_2[i]->SetLineColor(colors2[i]);
        rationch_2[i]->SetLineWidth(2);
        rationch_2[i]->GetYaxis()->SetRangeUser(0.55, 1.45);
        rationch_2[i]->GetXaxis()->SetRangeUser(0, 20);
        rationch_2[i]->GetYaxis()->SetTitle("Ratio to III");
        rationch_2[i]->GetXaxis()->SetTitleSize(0.1);
        rationch_2[i]->GetYaxis()->SetTitleSize(0.08);
        rationch_2[i]->GetXaxis()->SetLabelSize(0.09);
        rationch_2[i]->GetYaxis()->SetLabelSize(0.09);
        rationch_2[i]->GetYaxis()->SetTitleOffset(0.7);

        rationch_2[i]->Draw("LEP SAME");
        rationch_2[i]->Draw("HIST SAME");
    }

    c6->cd(2);
    //c6->cd(2)->SetLogy();

    TPad *pad1_22 = new TPad("pad1_22", "", 0, 0.3, 1, 1);
    TPad *pad2_22 = new TPad("pad2_22", "", 0, 0, 1, 0.3);
    pad1_22->SetBottomMargin(0.001);
    pad1_22->SetLeftMargin(0.12);
    pad1_22->SetRightMargin(0.02);
    pad1_22->SetTopMargin(0.02);
    pad1_22->Draw();
    pad2_22->SetBottomMargin(0.3);
    pad2_22->SetLeftMargin(0.12);
    pad2_22->SetRightMargin(0.02);
    pad2_22->SetTopMargin(0.001);
    pad2_22->Draw();
    pad1_22->Draw();
    pad1_22->cd();
    pad1_22->SetLogy();
    TH1D *ratiozn_2[size2];
    for (int i = 0; i < size2; i++)
    {
        pZN_2[i]->SetStats(0);
        pZN_2[i]->SetTitle("");
        pZN_2[i]->Rebin(10);
        pZN_2[i]->Scale(1. / pZN_2[i]->Integral(0, -1));
        pZN_2[i]->GetYaxis()->SetRangeUser(8E-4, 0.1);
        pZN_2[i]->GetXaxis()->SetRangeUser(0.0, 1020);
        pZN_2[i]->SetLineColor(colors2[i]);
        pZN_2[i]->SetLineWidth(2);
        pZN_2[i]->SetMarkerColor(colors2[i]);
        pZN_2[i]->SetMarkerStyle(kFullCircle);
        pZN_2[i]->Draw("LEP SAME");
        pZN_2[i]->Draw("HIST SAME");
    }

    TLegend *leg6bis = new TLegend(0.4, 0.7, 0.65, 0.9);
    leg6bis->SetBorderSize(0);
    leg6bis->SetFillColor(0);
    leg6bis->SetTextSize(0.03);
    for (int i = 0; i < size2; i++)
    {
        leg6bis->AddEntry(pnch_2[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_2[i], (int)perc1_high_2[i], (int)perc2_low_2[i], (int)perc2_high_2[i]), "LEP");
        //leg6bis->AddEntry(pZN_2[i], Form("#LT ZN #GT (a.u.) = %.2f", TMath::Abs(meanZN_2[i]) * meanZNMB), "LEP");
    }
    leg6bis->Draw();

    pad2_22->cd();
    pad2_22->SetLogy();
    for (int i = 0; i < size2; i++)
    {
        ratiozn_2[i] = (TH1D *)pZN_2[i]->Clone(Form("ratiozn_2%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_2[i], perc1_high_2[i], "V0M", perc2_low_2[i], perc2_high_2[i]));
        ratiozn_2[i]->Divide(pZN_2[2]);
        ratiozn_2[i]->SetLineColor(colors2[i]);
        ratiozn_2[i]->SetLineWidth(2);
        ratiozn_2[i]->GetYaxis()->SetRangeUser(0.45, 2.1);
        ratiozn_2[i]->GetXaxis()->SetRangeUser(0, 1020);
        ratiozn_2[i]->GetYaxis()->SetTitle("Ratio to III");
        ratiozn_2[i]->GetXaxis()->SetTitleSize(0.1);
        ratiozn_2[i]->GetYaxis()->SetTitleSize(0.08);
        ratiozn_2[i]->GetXaxis()->SetLabelSize(0.09);
        ratiozn_2[i]->GetYaxis()->SetLabelSize(0.09);
        ratiozn_2[i]->GetYaxis()->SetTitleOffset(0.7);

        ratiozn_2[i]->Draw("LEP SAME");
        ratiozn_2[i]->Draw("HIST SAME");
    }

    c6->SaveAs(Form("images/distributions%i.png", 2));

    TCanvas *c7 = new TCanvas("c7", "class3", 1300, 800);
    c7->Divide(2, 1);
    c7->cd(1);

    TPad *pad1_3 = new TPad("pad1_3", "", 0, 0.3, 1, 1);
    TPad *pad2_3 = new TPad("pad2_3", "", 0, 0, 1, 0.3);
    pad1_3->SetBottomMargin(0.001);
    pad1_3->SetLeftMargin(0.12);
    pad1_3->SetRightMargin(0.02);
    pad1_3->SetTopMargin(0.02);
    pad1_3->Draw();
    pad2_3->SetBottomMargin(0.3);
    pad2_3->SetLeftMargin(0.12);
    pad2_3->SetRightMargin(0.02);
    pad2_3->SetTopMargin(0.001);
    pad2_3->Draw();

    pad1_3->cd();
    TH1D *rationch_3[size3];
    // c7->cd(1)->SetLogy();
    for (int i = 0; i < size3; i++)
    {
        pnch_3[i]->SetStats(0);
        pnch_3[i]->SetTitle("");
        pnch_3[i]->Scale(1. / pnch_3[i]->Integral(0, -1));
        pnch_3[i]->GetYaxis()->SetRangeUser(0., 0.25);
        pnch_3[i]->GetXaxis()->SetRangeUser(0, 30);
        pnch_3[i]->SetLineColor(colors1[i]);
        pnch_3[i]->SetLineWidth(2);
        pnch_3[i]->SetMarkerColor(colors1[i]);
        pnch_3[i]->SetMarkerStyle(kFullCircle);
        pnch_3[i]->Draw("LEP SAME");
        pnch_3[i]->Draw("HIST SAME");
    }
    TLegend *leg8 = new TLegend(0.4, 0.7, 0.65, 0.9);
    leg8->SetBorderSize(0);
    leg8->SetFillColor(0);
    leg8->SetTextSize(0.03);
    for (int i = 0; i < size3; i++)
    {
        leg8->AddEntry(pnch_3[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_3[i], (int)perc1_high_3[i], (int)perc2_low_3[i], (int)perc2_high_3[i]), "LEP");
        //leg8->AddEntry(pnch_3[i], Form("#LT n_{ch} #GT = %.2f", meanNch_3[i] * meanNchMB), "LEP");
    }

    leg8->Draw();

    pad2_3->cd();
    pad2_3->SetLogy();
    for (int i = 0; i < size3; i++)
    {
        rationch_3[i] = (TH1D *)pnch_3[i]->Clone(Form("rationch_3%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_3[i], perc1_high_3[i], "V0M", perc2_low_3[i], perc2_high_3[i]));
        rationch_3[i]->Divide(pnch_3[2]);
        rationch_3[i]->SetLineColor(colors1[i]);
        rationch_3[i]->SetLineWidth(2);
        rationch_3[i]->GetYaxis()->SetRangeUser(0.002, 1000);
        rationch_3[i]->GetXaxis()->SetRangeUser(0, 30);
        rationch_3[i]->GetYaxis()->SetTitle("Ratio to III");
        rationch_3[i]->GetXaxis()->SetTitleSize(0.1);
        rationch_3[i]->GetYaxis()->SetTitleSize(0.08);
        rationch_3[i]->GetXaxis()->SetLabelSize(0.09);
        rationch_3[i]->GetYaxis()->SetLabelSize(0.09);
        rationch_3[i]->GetYaxis()->SetTitleOffset(0.7);

        rationch_3[i]->Draw("LEP SAME");
        rationch_3[i]->Draw("HIST SAME");
    }

    c7->cd(2);
//    c7->cd(2)->SetLogy();

    TPad *pad1_32 = new TPad("pad1_32", "", 0, 0.3, 1, 1);
    TPad *pad2_32 = new TPad("pad2_32", "", 0, 0, 1, 0.3);
    pad1_32->SetBottomMargin(0.001);
    pad1_32->SetLeftMargin(0.12);
    pad1_32->SetRightMargin(0.02);
    pad1_32->SetTopMargin(0.02);
    pad1_32->Draw();
    pad2_32->SetBottomMargin(0.3);
    pad2_32->SetLeftMargin(0.12);
    pad2_32->SetRightMargin(0.02);
    pad2_32->SetTopMargin(0.001);
    pad2_32->Draw();
    pad1_32->Draw();
    pad1_32->cd();
    pad1_32->SetLogy();
    TH1D *ratiozn_3[size3];
    for (int i = 0; i < size3; i++)
    {
        pZN_3[i]->SetStats(0);
        pZN_3[i]->SetTitle("");
        pZN_3[i]->Rebin(10);
        pZN_3[i]->Scale(1. / pZN_3[i]->Integral(0, -1));
        pZN_3[i]->GetYaxis()->SetRangeUser(8E-4, 0.2);
        pZN_3[i]->GetXaxis()->SetRangeUser(0.0, 1020);
        pZN_3[i]->SetLineColor(colors1[i]);
        pZN_3[i]->SetLineWidth(2);
        pZN_3[i]->SetMarkerColor(colors1[i]);
        pZN_3[i]->SetMarkerStyle(kFullCircle);
        pZN_3[i]->Draw("LEP SAME");
        pZN_3[i]->Draw("HIST SAME");
    }

    TLegend *leg8bis = new TLegend(0.4, 0.7, 0.65, 0.9);
    leg8bis->SetBorderSize(0);
    leg8bis->SetFillColor(0);
    leg8bis->SetTextSize(0.03);
    for (int i = 0; i < size3; i++)
    {
        leg8bis->AddEntry(pnch_3[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_3[i], (int)perc1_high_3[i], (int)perc2_low_3[i], (int)perc2_high_3[i]), "LEP");
        //leg8bis->AddEntry(pZN_3[i], Form("#LT ZN #GT (a.u.) = %.2f", TMath::Abs(meanZN_3[i]) * meanZNMB), "LEP");
    }
    leg8bis->Draw();

    pad2_32->cd();
    //pad2_32->SetLogy();
    for (int i = 0; i < size3; i++)
    {
        ratiozn_3[i] = (TH1D *)pZN_3[i]->Clone(Form("ratiozn_3%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_3[i], perc1_high_3[i], "V0M", perc2_low_3[i], perc2_high_3[i]));
        ratiozn_3[i]->Divide(pZN_3[2]);
        ratiozn_3[i]->SetLineColor(colors1[i]);
        ratiozn_3[i]->SetLineWidth(2);
        ratiozn_3[i]->GetYaxis()->SetRangeUser(0.45, 1.55);
        ratiozn_3[i]->GetXaxis()->SetRangeUser(0, 1020);
        ratiozn_3[i]->GetYaxis()->SetTitle("Ratio to III");
        ratiozn_3[i]->GetXaxis()->SetTitleSize(0.1);
        ratiozn_3[i]->GetYaxis()->SetTitleSize(0.08);
        ratiozn_3[i]->GetXaxis()->SetLabelSize(0.09);
        ratiozn_3[i]->GetYaxis()->SetLabelSize(0.09);
        ratiozn_3[i]->GetYaxis()->SetTitleOffset(0.7);

        ratiozn_3[i]->Draw("LEP SAME");
        ratiozn_3[i]->Draw("HIST SAME");

    }

    c7->SaveAs(Form("images/distributions%i.png", 3));

    TCanvas *c8 = new TCanvas("c8", "class4", 1300, 800);
    c8->Divide(2, 1);
    c8->cd(1);

    TPad *pad1_4 = new TPad("pad1_4", "", 0, 0.3, 1, 1);
    TPad *pad2_4 = new TPad("pad2_4", "", 0, 0, 1, 0.3);
    pad1_4->SetBottomMargin(0.001);
    pad1_4->SetLeftMargin(0.12);
    pad1_4->SetRightMargin(0.02);
    pad1_4->SetTopMargin(0.02);
    pad1_4->Draw();
    pad2_4->SetBottomMargin(0.3);
    pad2_4->SetLeftMargin(0.12);
    pad2_4->SetRightMargin(0.02);
    pad2_4->SetTopMargin(0.001);
    pad2_4->Draw();

    pad1_4->cd();
    TH1D *rationch_4[size4];
    // c8->cd(1)->SetLogy();
    for (int i = 0; i < size4; i++)
    {
        pnch_4[i]->SetStats(0);
        pnch_4[i]->SetTitle("");
        pnch_4[i]->Scale(1. / pnch_4[i]->Integral(0, -1));
        pnch_4[i]->GetYaxis()->SetRangeUser(0., 0.20);
        pnch_4[i]->GetXaxis()->SetRangeUser(0, 20);
        pnch_4[i]->SetLineColor(colors1[i]);
        pnch_4[i]->SetLineWidth(2);
        pnch_4[i]->SetMarkerColor(colors1[i]);
        pnch_4[i]->SetMarkerStyle(kFullCircle);
        pnch_4[i]->Draw("LEP SAME");
        pnch_4[i]->Draw("HIST SAME");
    }
    TLegend *leg11 = new TLegend(0.4, 0.7, 0.65, 0.9);
    leg11->SetBorderSize(0);
    leg11->SetFillColor(0);
    leg11->SetTextSize(0.03);
    for (int i = 0; i < size4; i++)
    {
        leg11->AddEntry(pnch_4[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_4[i], (int)perc1_high_4[i], (int)perc2_low_4[i], (int)perc2_high_4[i]), "LEP");
        //leg11->AddEntry(pnch_4[i], Form("#LT n_{ch} #GT = %.2f", meanNch_4[i] * meanNchMB), "LEP");
    }

    leg11->Draw();

    pad2_4->cd();
    pad2_4->SetLogy();
    for (int i = 0; i < size4; i++)
    {
        rationch_4[i] = (TH1D *)pnch_4[i]->Clone(Form("rationch_4%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_4[i], perc1_high_4[i], "V0M", perc2_low_4[i], perc2_high_4[i]));
        rationch_4[i]->Divide(pnch_4[2]);
        rationch_4[i]->SetLineColor(colors1[i]);
        rationch_4[i]->SetLineWidth(2);
        rationch_4[i]->GetYaxis()->SetRangeUser(0.01, 100);
        rationch_4[i]->GetXaxis()->SetRangeUser(0, 20);
        rationch_4[i]->GetYaxis()->SetTitle("Ratio to III");
        rationch_4[i]->GetXaxis()->SetTitleSize(0.1);
        rationch_4[i]->GetYaxis()->SetTitleSize(0.08);
        rationch_4[i]->GetXaxis()->SetLabelSize(0.09);
        rationch_4[i]->GetYaxis()->SetLabelSize(0.09);
        rationch_4[i]->GetYaxis()->SetTitleOffset(0.7);

        rationch_4[i]->Draw("LEP SAME");
        rationch_4[i]->Draw("HIST SAME");
    }

    c8->cd(2);
    c8->cd(2)->SetLogy();

    TPad *pad1_42 = new TPad("pad1_42", "", 0, 0.3, 1, 1);
    TPad *pad2_42 = new TPad("pad2_42", "", 0, 0, 1, 0.3);
    pad1_42->SetBottomMargin(0.001);
    pad1_42->SetLeftMargin(0.12);
    pad1_42->SetRightMargin(0.02);
    pad1_42->SetTopMargin(0.02);
    pad1_42->Draw();
    pad2_42->SetBottomMargin(0.3);
    pad2_42->SetLeftMargin(0.12);
    pad2_42->SetRightMargin(0.02);
    pad2_42->SetTopMargin(0.001);
    pad2_42->Draw();
    pad1_42->Draw();
    pad1_42->cd();
    pad1_42->SetLogy();
    TH1D *ratiozn_4[size4];
    for (int i = 0; i < size4; i++)
    {
        pZN_4[i]->SetStats(0);
        pZN_4[i]->SetTitle("");
        pZN_4[i]->Rebin(10);
        pZN_4[i]->Scale(1. / pZN_4[i]->Integral(0, -1));
        pZN_4[i]->GetYaxis()->SetRangeUser(1E-4, 0.1);
        pZN_4[i]->GetXaxis()->SetRangeUser(0.0, 1020);
        pZN_4[i]->SetLineColor(colors1[i]);
        pZN_4[i]->SetLineWidth(2);
        pZN_4[i]->SetMarkerColor(colors1[i]);
        pZN_4[i]->SetMarkerStyle(kFullCircle);
        pZN_4[i]->Draw("LEP SAME");
        pZN_4[i]->Draw("HIST SAME");
    }

    TLegend *leg11bis = new TLegend(0.4, 0.7, 0.65, 0.9);
    leg11bis->SetBorderSize(0);
    leg11bis->SetFillColor(0);
    leg11bis->SetTextSize(0.03);
    for (int i = 0; i < size4; i++)
    {
        leg11bis->AddEntry(pnch_4[i], Form("%d-%d%% SPDClusters, %d-%d%% V0M", (int)perc1_low_4[i], (int)perc1_high_4[i], (int)perc2_low_4[i], (int)perc2_high_4[i]), "LEP");
        //leg11bis->AddEntry(pZN_4[i], Form("#LT ZN #GT (a.u.) = %.2f", TMath::Abs(meanZN_4[i]) * meanZNMB), "LEP");
    }
    leg11bis->Draw();

    pad2_42->cd();
    // pad2_42->SetLogy();
    for (int i = 0; i < size4; i++)
    {
        ratiozn_4[i] = (TH1D *)pZN_4[i]->Clone(Form("ratiozn_4%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_4[i], perc1_high_4[i], "V0M", perc2_low_4[i], perc2_high_4[i]));
        ratiozn_4[i]->Divide(pZN_4[2]);
        ratiozn_4[i]->SetLineColor(colors1[i]);
        ratiozn_4[i]->SetLineWidth(2);
        ratiozn_4[i]->GetYaxis()->SetRangeUser(0.55, 1.45);
        ratiozn_4[i]->GetXaxis()->SetRangeUser(0, 1020);
        ratiozn_4[i]->GetYaxis()->SetTitle("Ratio to III");
        ratiozn_4[i]->GetXaxis()->SetTitleSize(0.1);
        ratiozn_4[i]->GetYaxis()->SetTitleSize(0.08);
        ratiozn_4[i]->GetXaxis()->SetLabelSize(0.09);
        ratiozn_4[i]->GetYaxis()->SetLabelSize(0.09);
        ratiozn_4[i]->GetYaxis()->SetTitleOffset(0.7);

        ratiozn_4[i]->Draw("LEP SAME");
        ratiozn_4[i]->Draw("HIST SAME");
    }

    c8->SaveAs(Form("images/distributions4.png"));
}

Double_t getval(TFile *filecalib, Int_t run, Int_t value)
{
    // Return the percentile given the value of estimator
    TString name = Form("hcumSPDClusters_%i", run);

    TH1D *hcum = (TH1D *)filecalib->Get(name);
    Double_t percentile = 100 * (hcum->GetBinContent(hcum->FindBin(value)));

    return percentile;
}
