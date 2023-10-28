Double_t getval(TFile *filecalib, Int_t run, Int_t value);

void prepcanvas(TCanvas* c);
void preppad1(TPad *pad);
void preppad2(TPad* pad);
void ratiohisto(TH1D* h);
void histo(TH1D *h);

void drawdistributions(){

    TFile *filecontainer = TFile::Open("../Selections/NchRawContainer_reference.root");
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

    // 1-->HighMult
    Double_t perc1_low_1[] = {10, 10, 10, 10, 10, 10, 10};
    Double_t perc1_high_1[] = {20, 20, 20, 20, 20, 20, 20};
    Double_t perc2_low_1[] = {0, 5, 10, 20, 30, 40, 50};
    Double_t perc2_high_1[] = {5, 10, 20, 30, 40, 50, 100};
    const int size1 = sizeof(perc1_low_1) / sizeof(Double_t);

    // 2-->LowMult
    Double_t perc1_low_2[] = {40, 40, 40, 40, 40, 40, 40};
    Double_t perc1_high_2[] = {50, 50, 50, 50, 50, 50, 50};
    Double_t perc2_low_2[] = {0, 20, 30, 40, 50, 60, 70};
    Double_t perc2_high_2[] = {20, 30, 40, 50, 60, 70, 100};
    const int size2 = sizeof(perc1_low_2) / sizeof(Double_t);

    // 3-->HighZN
    Double_t perc1_low_3[] = {0, 10, 20, 30, 50};
    Double_t perc1_high_3[] = {20, 30, 40, 50, 100};
    Double_t perc2_low_3[] = {40, 30, 30, 20, 0};
    Double_t perc2_high_3[] = {60, 70, 50, 50, 30};
    const int size3 = sizeof(perc1_low_3) / sizeof(Double_t);
    // 4-->LowZN
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

    Int_t colors1[] = {kRed + 2, kOrange + 1, kSpring - 1, kAzure + 7, kBlue + 2};
    Int_t colors2[] = {kRed + 2, kOrange + 1, kYellow + 1, kSpring - 1, kGreen + 2, kAzure + 7, kBlue + 2};

    TString classe[] = {"I", "II", "III", "IV", "V","VI","VII","VIII","IX","X"};

    TCanvas *c5 = new TCanvas("c5", "class1", 1300, 800);
    prepcanvas(c5);
    c5->cd(1);

    TPad *pad1 = new TPad("pad1", "", 0, 0.3, 1, 1);
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    preppad1(pad1);
    preppad2(pad2);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();

    TLatex* l = new TLatex();
    l->SetTextSize(0.04);
    l->SetTextFont(42);
    l->SetNDC();
    l->SetTextAlign(12);



    TH1D *rationch[size1];
    for (int i = 0; i < size1; i++)
    {
        histo(pnch_1[i]);
        pnch_1[i]->SetLineColor(colors2[i]);
        pnch_1[i]->SetMarkerColor(colors2[i]);
        pnch_1[i]->GetXaxis()->SetRangeUser(0, 30);
        pnch_1[i]->GetYaxis()->SetRangeUser(-0.003, pnch_1[i]->GetMaximum() * 1.2);
        pnch_1[i]->Draw("LEP SAME");
        pnch_1[i]->Draw("HIST SAME");
    }

    TLegend *leg5 = new TLegend(0.62, 0.6, 0.8, 0.91);
    leg5->SetBorderSize(0);
    leg5->SetFillColor(0);
    leg5->SetTextSize(0.04);
    leg5->SetHeader("High Multiplicity Class");
    for (int i = 0; i < size1; i++)
    {
        leg5->AddEntry(pnch_1[i], Form("%s",classe[i].Data()), "LEP");
    }
    leg5->Draw();

    l->DrawLatex(0.18, 0.9, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
    //l->DrawLatex(0.6, 0.84, "High Multiplicity Class");

    pad2->cd();
    for (int i = 0; i < size1; i++) {
        rationch[i] = (TH1D *)pnch_1[i]->Clone(Form("rationch%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_1[i], perc1_high_1[i], "V0M", perc2_low_1[i], perc2_high_1[i]));
        rationch[i]->Divide(pnch_1[2]);
        rationch[i]->SetLineColor(colors2[i]);
        ratiohisto(rationch[i]);
        rationch[i]->GetXaxis()->SetRangeUser(0, 30);
        rationch[i]->GetYaxis()->SetRangeUser(0, 1.99);
        rationch[i]->GetXaxis()->SetTitle("tracklets");
        if (i != 2)
            rationch[i]->Draw("LEP SAME");
    }

    c5->cd(2);
    //c5->cd(2)->SetLogy();

    TPad *pad1_ = new TPad("pad1_", "", 0, 0.3, 1, 1);
    TPad *pad2_ = new TPad("pad2_", "", 0, 0, 1, 0.3);
    preppad1(pad1_);
    preppad2(pad2_);
    pad2_->Draw();
    pad1_->Draw();

    pad1_->cd();
    pad1_->SetLogy();
    TH1D *ratiozn[size1];
    for (int i = 0; i < size1; i++)
    {
        pZN_1[i]->Rebin(15);
        histo(pZN_1[i]);
        pZN_1[i]->GetXaxis()->SetRangeUser(0, 1020);
        pZN_1[i]->SetLineColor(colors2[i]);
        pZN_1[i]->SetMarkerColor(colors2[i]);
        pZN_1[i]->Draw("LEP SAME");
        pZN_1[i]->Draw("HIST SAME");
    }

    l->DrawLatex(0.18, 0.9, "ALICE, pp #sqrt{#it{s}} = 13 TeV");

    TLegend *leg4bis = new TLegend(0.75, 0.55, 0.88, 0.88);
    leg4bis->SetBorderSize(0);
    leg4bis->SetFillColor(0);
    leg4bis->SetTextSize(0.04);
    for (int i = 0; i < size1; i++)
    {
        leg4bis->AddEntry(pnch_1[i], Form("%s", classe[i].Data()), "LEP");
    }
    leg5->Draw();

    pad2_->cd();
    pad2_->SetLogy();
    for (int i = 0; i < size1; i++)
    {
        ratiozn[i] = (TH1D *)pZN_1[i]->Clone(Form("ratiozn%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_1[i], perc1_high_1[i], "V0M", perc2_low_1[i], perc2_high_1[i]));
        ratiozn[i]->Divide(pZN_1[2]);
        ratiozn[i]->SetLineColor(colors2[i]);
        ratiozn[i]->GetXaxis()->SetRangeUser(0, 1020);
        ratiohisto(ratiozn[i]);
        ratiozn[i]->GetYaxis()->SetRangeUser(0.3, 6);
        ratiozn[i]->GetXaxis()->SetTitle("ZN (a.u.)");
        if (i != 2)
            ratiozn[i]->Draw("LEP SAME");
    }
    c5->SaveAs(Form("images/distributions%i.png", 1));
    c5->SaveAs(Form("images/distributions_highmultclass.eps"));

    TCanvas *c6 = new TCanvas("c6", "class2", 1300, 800);
    c6->Divide(2, 1);
    c6->cd(1);

    TPad *pad1_2 = new TPad("pad1_2", "", 0, 0.3, 1, 1);
    TPad *pad2_2 = new TPad("pad2_2", "", 0, 0, 1, 0.3);
    preppad1(pad1_2);
    preppad2(pad2_2);
    pad1_2->Draw();
    pad2_2->Draw();

    pad1_2->cd();
    TH1D *rationch_2[size2];
    // c6->cd(1)->SetLogy();
    for (int i = 0; i < size2; i++)
    {
        pnch_2[i]->SetStats(0);
        pnch_2[i]->SetTitle("");
        pnch_2[i]->Scale(1. / pnch_2[i]->Integral(0, -1));
        pnch_2[i]->GetYaxis()->SetRangeUser(-0.003, 0.23);
        pnch_2[i]->GetXaxis()->SetRangeUser(0, 15);
        pnch_2[i]->SetLineColor(colors2[i]);
        //pnch_2[i]->SetLineWidth(2);
        pnch_2[i]->SetMarkerColor(colors2[i]);
        pnch_2[i]->SetMarkerStyle(kFullCircle);
        pnch_2[i]->GetYaxis()->SetTitle("Normalized counts");
        pnch_2[i]->GetYaxis()->SetTitleSize(0.045);
        pnch_2[i]->GetYaxis()->SetTitleOffset(1.2);
        pnch_2[i]->Draw("LEP SAME");
        pnch_2[i]->Draw("HIST SAME");
    }

    l->DrawLatex(0.18, 0.9, "ALICE, pp #sqrt{#it{s}} = 13 TeV");

    TLegend *leg6 = new TLegend(0.62, 0.6, 0.8, 0.91);
    leg6->SetBorderSize(0);
    leg6->SetFillColor(0);
    leg6->SetTextSize(0.04);
    leg6->SetHeader("Low Multiplicity Class");
    for (int i = 0; i < size2; i++)
    {
        leg6->AddEntry(pnch_2[i], Form("%s", classe[i].Data()), "LEP");
    }

    leg6->Draw();

    pad2_2->cd();
    // pad2_2->SetLogy();
    for (int i = 0; i < size2; i++)
    {
        rationch_2[i] = (TH1D *)pnch_2[i]->Clone(Form("rationch_2%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_2[i], perc1_high_2[i], "V0M", perc2_low_2[i], perc2_high_2[i]));
        rationch_2[i]->Divide(pnch_2[2]);
        rationch_2[i]->SetLineColor(colors2[i]);
        ratiohisto(rationch_2[i]);
        //rationch_2[i]->SetLineWidth(2);
        rationch_2[i]->GetXaxis()->SetRangeUser(0, 15);
        rationch_2[i]->GetYaxis()->SetRangeUser(0, 1.99);
        rationch_2[i]->GetXaxis()->SetTitle("tracklets");
        if (i!=2) rationch_2[i]->Draw("LEP SAME");
        //rationch_2[i]->Draw("HIST SAME");
    }

    c6->cd(2);
    //c6->cd(2)->SetLogy();

    TPad *pad1_22 = new TPad("pad1_22", "", 0, 0.3, 1, 1);
    TPad *pad2_22 = new TPad("pad2_22", "", 0, 0, 1, 0.3);
    preppad1(pad1_22);
    preppad2(pad2_22);
    pad2_22->Draw();
    pad1_22->Draw();

    pad1_22->cd();
    pad1_22->SetLogy();
    TH1D *ratiozn_2[size2];
    for (int i = 0; i < size2; i++)
    {
        pZN_2[i]->Rebin(15);
        histo(pZN_2[i]);
        pZN_2[i]->GetXaxis()->SetRangeUser(0, 1020);
        pZN_2[i]->SetLineColor(colors2[i]);
        pZN_2[i]->SetMarkerColor(colors2[i]);
        pZN_2[i]->Draw("LEP SAME");
        pZN_2[i]->Draw("HIST SAME");
    }

    leg6->Draw();
    l->DrawLatex(0.18, 0.9, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
    pad2_22->cd();
    pad2_22->SetLogy();
    for (int i = 0; i < size2; i++)
    {
        ratiozn_2[i] = (TH1D *)pZN_2[i]->Clone(Form("ratiozn_2%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_2[i], perc1_high_2[i], "V0M", perc2_low_2[i], perc2_high_2[i]));
        ratiozn_2[i]->Divide(pZN_2[2]);
        ratiozn_2[i]->SetLineColor(colors2[i]);
        ratiozn_2[i]->GetXaxis()->SetRangeUser(0, 1020);
        ratiohisto(ratiozn_2[i]);
        ratiozn_2[i]->GetYaxis()->SetRangeUser(0.4, 3);
        ratiozn_2[i]->GetXaxis()->SetTitle("ZN (a.u.)");
        if (i != 2)
            ratiozn_2[i]->Draw("LEP SAME");
        //ratiozn_2[i]->Draw("HIST SAME");
    }

    c6->SaveAs(Form("images/distributions%i.png", 2));
    c6->SaveAs(Form("images/distributions_lowmultclass.eps"));

    TCanvas *c7 = new TCanvas("c7", "class3", 1300, 800);
    c7->Divide(2, 1);
    c7->cd(1);

    TPad *pad1_3 = new TPad("pad1_3", "", 0, 0.3, 1, 1);
    TPad *pad2_3 = new TPad("pad2_3", "", 0, 0, 1, 0.3);
    preppad1(pad1_3);
    preppad2(pad2_3);
    pad1_3->Draw();
    pad2_3->Draw();

    pad1_3->cd();
    TH1D *rationch_3[size3];
    // c7->cd(1)->SetLogy();
    for (int i = 0; i < size3; i++)
    {
        pnch_3[i]->SetStats(0);
        pnch_3[i]->SetTitle("");
        pnch_3[i]->Scale(1. / pnch_3[i]->Integral(0, -1));
        pnch_3[i]->GetYaxis()->SetRangeUser(-0.003, 0.25);
        pnch_3[i]->GetXaxis()->SetRangeUser(0, 30);
        pnch_3[i]->SetLineColor(colors1[i]);
        //pnch_3[i]->SetLineWidth(2);
        pnch_3[i]->SetMarkerColor(colors1[i]);
        pnch_3[i]->SetMarkerStyle(kFullCircle);
        pnch_3[i]->GetYaxis()->SetTitle("Normalized counts");
        pnch_3[i]->GetYaxis()->SetTitleSize(0.045);
        pnch_3[i]->GetYaxis()->SetTitleOffset(1.2);
        pnch_3[i]->Draw("LEP SAME");
        pnch_3[i]->Draw("HIST SAME");
    }
    TLegend *leg8 = new TLegend(0.62, 0.6, 0.8, 0.91);
    leg8->SetBorderSize(0);
    leg8->SetFillColor(0);
    leg8->SetTextSize(0.04);
    leg8->SetHeader("High ZN Class");
    for (int i = 0; i < size3; i++)
    {
        leg8->AddEntry(pnch_3[i], Form("%s", classe[i].Data()), "LEP");
    }

    l->DrawLatex(0.18, 0.9, "ALICE, pp #sqrt{#it{s}} = 13 TeV");

    leg8->Draw();

    pad2_3->cd();
    pad2_3->SetLogy();
    for (int i = 0; i < size3; i++)
    {
        rationch_3[i] = (TH1D *)pnch_3[i]->Clone(Form("rationch_3%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_3[i], perc1_high_3[i], "V0M", perc2_low_3[i], perc2_high_3[i]));
        rationch_3[i]->Divide(pnch_3[2]);
        rationch_3[i]->SetLineColor(colors1[i]);
        //rationch_3[i]->SetLineWidth(2);
        rationch_3[i]->GetYaxis()->SetRangeUser(0.0001, 999);
        rationch_3[i]->GetXaxis()->SetRangeUser(0, 30);
        rationch_3[i]->GetYaxis()->SetTitle("Ratio to III");
        rationch_3[i]->GetXaxis()->SetTitle("tracklets");
        rationch_3[i]->GetXaxis()->SetTitleSize(0.1);
        rationch_3[i]->GetYaxis()->SetTitleSize(0.08);
        rationch_3[i]->GetXaxis()->SetLabelSize(0.09);
        rationch_3[i]->GetYaxis()->SetLabelSize(0.09);
        rationch_3[i]->GetYaxis()->SetTitleOffset(0.7);

        if (i != 2)
            rationch_3[i]->Draw("LEP SAME");
        //rationch_3[i]->Draw("HIST SAME");
    }

    c7->cd(2);
//    c7->cd(2)->SetLogy();

    TPad *pad1_32 = new TPad("pad1_32", "", 0, 0.3, 1, 1);
    TPad *pad2_32 = new TPad("pad2_32", "", 0, 0, 1, 0.3);
    preppad1(pad1_32);
    preppad2(pad2_32);

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
        pZN_3[i]->GetYaxis()->SetRangeUser(3E-4, 0.2);
        pZN_3[i]->GetXaxis()->SetRangeUser(0.0, 1020);
        pZN_3[i]->SetLineColor(colors1[i]);
        //pZN_3[i]->SetLineWidth(2);
        pZN_3[i]->SetMarkerColor(colors1[i]);
        pZN_3[i]->SetMarkerStyle(kFullCircle);
        pZN_3[i]->GetYaxis()->SetTitle("Normalized counts");
        pZN_3[i]->GetYaxis()->SetTitleSize(0.045);
        pZN_3[i]->GetYaxis()->SetTitleOffset(1.2);
        pZN_3[i]->Draw("LEP SAME");
        pZN_3[i]->Draw("HIST SAME");
    }

    l->DrawLatex(0.18, 0.9, "ALICE, pp #sqrt{#it{s}} = 13 TeV");

    leg8->Draw();


    pad2_32->cd();
    //pad2_32->SetLogy();
    for (int i = 0; i < size3; i++)
    {
        ratiozn_3[i] = (TH1D *)pZN_3[i]->Clone(Form("ratiozn_3%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_3[i], perc1_high_3[i], "V0M", perc2_low_3[i], perc2_high_3[i]));
        ratiozn_3[i]->Divide(pZN_3[2]);
        ratiozn_3[i]->SetLineColor(colors1[i]);
        //ratiozn_3[i]->SetLineWidth(2);
        ratiozn_3[i]->GetYaxis()->SetRangeUser(0.45, 1.55);
        ratiozn_3[i]->GetXaxis()->SetRangeUser(0, 1020);
        ratiozn_3[i]->GetYaxis()->SetTitle("Ratio to III");
        ratiozn_3[i]->GetXaxis()->SetTitleSize(0.1);
        ratiozn_3[i]->GetYaxis()->SetTitleSize(0.08);
        ratiozn_3[i]->GetXaxis()->SetLabelSize(0.08);
        ratiozn_3[i]->GetYaxis()->SetLabelSize(0.08);
        ratiozn_3[i]->GetYaxis()->SetTitleOffset(0.7);

        if (i != 2)
            ratiozn_3[i]->Draw("LEP SAME");
        //ratiozn_3[i]->Draw("HIST SAME");

    }

    c7->SaveAs(Form("images/distributions%i.png", 3));
    c7->SaveAs(Form("images/distributions_higznclass.eps"));

    TCanvas *c8 = new TCanvas("c8", "class3", 1300, 800);
    c8->Divide(2, 1);
    c8->cd(1);

    TPad *pad1_4 = new TPad("pad1_4", "", 0, 0.3, 1, 1);
    TPad *pad2_4 = new TPad("pad2_4", "", 0, 0, 1, 0.3);
    preppad1(pad1_4);
    preppad2(pad2_4);
    pad1_4->Draw();
    pad2_4->Draw();

    pad1_4->cd();
    TH1D *rationch_4[size4];
    // c8->cd(1)->SetLogy();
    for (int i = 0; i < size4; i++)
    {
        pnch_4[i]->SetStats(0);
        pnch_4[i]->SetTitle("");
        pnch_4[i]->Scale(1. / pnch_4[i]->Integral(0, -1));
        pnch_4[i]->GetYaxis()->SetRangeUser(-0.003, 0.20);
        pnch_4[i]->GetXaxis()->SetRangeUser(0, 30);
        pnch_4[i]->SetLineColor(colors1[i]);
        //pnch_4[i]->SetLineWidth(2);
        pnch_4[i]->SetMarkerColor(colors1[i]);
        pnch_4[i]->SetMarkerStyle(kFullCircle);
        pnch_4[i]->GetYaxis()->SetTitle("Normalized counts");
        pnch_4[i]->GetYaxis()->SetTitleSize(0.045);
        pnch_4[i]->GetYaxis()->SetTitleOffset(1.2);
        pnch_4[i]->Draw("LEP SAME");
        pnch_4[i]->Draw("HIST SAME");
    }
    TLegend *leg11 = new TLegend(0.62, 0.6, 0.8, 0.91);
    leg11->SetBorderSize(0);
    leg11->SetFillColor(0);
    leg11->SetTextSize(0.04);
    leg11->SetHeader("Low ZN Class");
    for (int i = 0; i < size4; i++)
    {
        leg11->AddEntry(pnch_4[i], Form("%s", classe[i].Data()), "LEP");
    }
    leg11->Draw();

    l->DrawLatex(0.18, 0.9, "ALICE, pp #sqrt{#it{s}} = 13 TeV");

    pad2_4->cd();
    pad2_4->SetLogy();
    for (int i = 0; i < size4; i++)
    {
        rationch_4[i] = (TH1D *)pnch_4[i]->Clone(Form("rationch_4%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_4[i], perc1_high_4[i], "V0M", perc2_low_4[i], perc2_high_4[i]));
        rationch_4[i]->Divide(pnch_4[2]);
        rationch_4[i]->SetLineColor(colors1[i]);
        ratiohisto(rationch_4[i]);
        //rationch_4[i]->SetLineWidth(2);
        rationch_4[i]->GetXaxis()->SetRangeUser(0, 30);
        rationch_4[i]->GetYaxis()->SetRangeUser(0.009, 3E+4);
        rationch_4[i]->GetXaxis()->SetTitle("tracklets");
        if (i!=2) rationch_4[i]->Draw("LEP SAME");
    }

    c8->cd(2);
    c8->cd(2)->SetLogy();

    TPad *pad1_42 = new TPad("pad1_42", "", 0, 0.3, 1, 1);
    TPad *pad2_42 = new TPad("pad2_42", "", 0, 0, 1, 0.3);
    preppad1(pad1_42);
    preppad2(pad2_42);

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
        pZN_4[i]->GetYaxis()->SetRangeUser(9E-5, 0.1);
        pZN_4[i]->GetXaxis()->SetRangeUser(0.0, 1020);
        pZN_4[i]->SetLineColor(colors1[i]);
        //pZN_4[i]->SetLineWidth(2);
        pZN_4[i]->SetMarkerColor(colors1[i]);
        pZN_4[i]->SetMarkerStyle(kFullCircle);
        pZN_4[i]->Draw("LEP SAME");
        pZN_4[i]->Draw("HIST SAME");
    }

    l->DrawLatex(0.18, 0.9, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
    leg11->Draw();

    pad2_42->cd();
    // pad2_42->SetLogy();
    for (int i = 0; i < size4; i++)
    {
        ratiozn_4[i] = (TH1D *)pZN_4[i]->Clone(Form("ratiozn_4%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", perc1_low_4[i], perc1_high_4[i], "V0M", perc2_low_4[i], perc2_high_4[i]));
        ratiozn_4[i]->Divide(pZN_4[2]);
        ratiozn_4[i]->SetLineColor(colors1[i]);
        //ratiozn_4[i]->SetLineWidth(2);
        ratiozn_4[i]->GetYaxis()->SetRangeUser(0., 1.99);
        ratiozn_4[i]->GetXaxis()->SetRangeUser(0, 1020);
        ratiozn_4[i]->GetYaxis()->SetTitle("Ratio to III");
        ratiozn_4[i]->GetXaxis()->SetTitleSize(0.1);
        ratiozn_4[i]->GetYaxis()->SetTitleSize(0.08);
        ratiozn_4[i]->GetXaxis()->SetLabelSize(0.09);
        ratiozn_4[i]->GetYaxis()->SetLabelSize(0.09);
        ratiozn_4[i]->GetYaxis()->SetTitleOffset(0.7);

        if (i != 2)
            ratiozn_4[i]->Draw("LEP SAME");
        // ratiozn_4[i]->Draw("HIST SAME");
    }

    c8->SaveAs(Form("images/distributions4.png"));
    c8->SaveAs(Form("images/distributions_lowznclass.eps"));
}

Double_t getval(TFile *filecalib, Int_t run, Int_t value)
{
    // Return the percentile given the value of estimator
    TString name = Form("hcumSPDClusters_%i", run);

    TH1D *hcum = (TH1D *)filecalib->Get(name);
    Double_t percentile = 100 * (hcum->GetBinContent(hcum->FindBin(value)));

    return percentile;
}

void prepcanvas(TCanvas *c){
    c->Divide(2, 1);
}

void preppad1(TPad* pad){
    pad->SetBottomMargin(0.001);
    pad->SetLeftMargin(0.12);
    pad->SetRightMargin(0.02);
    pad->SetTopMargin(0.02);
    pad->SetTicks();
}

void preppad2(TPad* pad){
    pad->SetBottomMargin(0.3);
    pad->SetLeftMargin(0.12);
    pad->SetRightMargin(0.02);
    pad->SetTopMargin(0.001);
    pad->SetTicks();
}

void ratiohisto(TH1D* h){
    // rationch[i]->SetLineWidth(2);
    h->GetYaxis()->SetRangeUser(0., 1.9);
    //h->GetXaxis()->SetRangeUser(0., 30);
    h->GetYaxis()->SetTitle("Ratio to III");
    h->GetXaxis()->SetTitleSize(0.1);
    h->GetYaxis()->SetTitleSize(0.08);
    h->GetXaxis()->SetLabelSize(0.08);
    h->GetYaxis()->SetLabelSize(0.08);
    h->GetYaxis()->SetTitleOffset(0.7);
}

void histo(TH1D* h){
    h->SetStats(0);
    h->SetTitle("");
    h->GetYaxis()->SetTitle("Normalized counts");
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->Scale(1. / h->Integral(0, -1));
    //h->GetYaxis()->SetRangeUser(-0.003, h->GetMaximum() * 1.2);
    h->SetMarkerStyle(kFullCircle);
    //h->SetLineWidth(2);
}
