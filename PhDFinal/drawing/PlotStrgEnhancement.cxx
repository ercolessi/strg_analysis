void preppad(TPad *pad1, TPad *pad2);
void prepcanvas(TCanvas* c);
void prepgraph(TGraphErrors* g, int marker, float size, int color);

void PlotStrgEnhancement(TString num = "Xi", TString den = "K0Short"){

    TString sel[] = {"V0M_SelectedWithSPDClusters0100", "V0M_SelectedWithSPDClusters1020", "V0M_SelectedWithSPDClusters4050", "SPDClusters_SelectedWithV0M1020", "SPDClusters_SelectedWithV0M4050"};
    const int nsel = sizeof(sel)/sizeof(TString);
    TFile *fNum[nsel], *fDen[nsel];

    for (int i = 0; i < nsel; i++){
        fNum[i] = TFile::Open(Form("/home/fercoles/strg_analysis/AnalisiFinale/yields/%sYields%s.root", num.Data(), sel[i].Data()));
        fDen[i] = TFile::Open(Form("/home/fercoles/strg_analysis/AnalisiFinale/yields/%sYields%s.root", den.Data(), sel[i].Data()));
    }

    TGraphErrors *gNumNchStat[nsel], *gNumZDCStat[nsel], *gDenNchStat[nsel];
    TGraphErrors *gNumNchSyst[nsel], *gNumZDCSyst[nsel], *gDenNchSyst[nsel];

    for (int i = 0; i < nsel; i++){
        gNumNchStat[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormNchStat");
        gNumZDCStat[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormZDCSumStat");
        gDenNchStat[i] = (TGraphErrors *)fDen[i]->Get("NormYieldsNormNchStat");
        gNumNchSyst[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormNchSyst");
        gNumZDCSyst[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormZDCSumSyst");
        gDenNchSyst[i] = (TGraphErrors *)fDen[i]->Get("NormYieldsNormNchSyst");
    }

    std::vector<Double_t> Nch[nsel], NchErr[nsel], NchErrSyst[nsel], ZDC[nsel], ZDCErr[nsel], ZDCErrSyst[nsel], YNum[nsel], YNumErr[nsel], YNumErrSyst[nsel], YDen[nsel], YDenErr[nsel], YDenErrSyst[nsel];
    std::vector<Double_t> RatioStat[nsel], RatioStatErr[nsel], RatioSyst[nsel], RatioSystErr[nsel];
    TGraphErrors *gNchRatioStat[nsel], *gZDCRatioStat[nsel], *gNchRatioSyst[nsel], *gZDCRatioSyst[nsel];

    for (int i = 0; i < nsel; i++){
        if (gNumNchStat[i]->GetN() != gDenNchStat[i]->GetN()) {
            cout << "ERROR" << endl;
            return;
        }

        for (int j = 0; j < gNumNchStat[i]->GetN(); j++){
            Nch[i].push_back(gNumNchStat[i]->GetX()[j]);
            NchErr[i].push_back(gNumNchStat[i]->GetEX()[j]);
            NchErrSyst[i].push_back(gNumNchSyst[i]->GetEXhigh()[j]);
            ZDC[i].push_back(gNumZDCStat[i]->GetX()[j]);
            ZDCErr[i].push_back(gNumZDCStat[i]->GetEX()[j]);
            ZDCErrSyst[i].push_back(gNumZDCSyst[i]->GetEXhigh()[j]);
            YNum[i].push_back(gNumNchStat[i]->GetY()[j]);
            YNumErr[i].push_back(gNumNchStat[i]->GetEY()[j]);
            YNumErrSyst[i].push_back(gNumNchSyst[i]->GetEYhigh()[j]);
            YDen[i].push_back(gDenNchStat[i]->GetY()[j]);
            YDenErr[i].push_back(gDenNchStat[i]->GetEY()[j]);
            YDenErrSyst[i].push_back(gDenNchSyst[i]->GetEYhigh()[j]);
        }

        for (int j = 0; j < gNumNchStat[i]->GetN(); j++){
            RatioStat[i].push_back(YNum[i][j]/YDen[i][j]);
            RatioStatErr[i].push_back(RatioStat[i][j] * TMath::Sqrt(TMath::Power(YNumErr[i][j]/YNum[i][j], 2) + TMath::Power(YDenErr[i][j]/YDen[i][j], 2)));
            RatioSyst[i].push_back(YNum[i][j]/YDen[i][j]);
            RatioSystErr[i].push_back(RatioSyst[i][j] * TMath::Sqrt(YNumErrSyst[i][j]/YNum[i][j]*YNumErrSyst[i][j]/YNum[i][j] + YDenErrSyst[i][j]/YDen[i][j]*YDenErrSyst[i][j]/YDen[i][j]));
        }

        gNchRatioStat[i] = new TGraphErrors(gNumNchStat[i]->GetN(), &Nch[i][0], &RatioStat[i][0], &NchErr[i][0], &RatioStatErr[i][0]);
        gZDCRatioStat[i] = new TGraphErrors(gNumNchStat[i]->GetN(), &ZDC[i][0], &RatioStat[i][0], &ZDCErr[i][0], &RatioStatErr[i][0]);
        gNchRatioSyst[i] = new TGraphErrors(gNumNchSyst[i]->GetN(), &Nch[i][0], &RatioSyst[i][0], &NchErrSyst[i][0], &RatioSystErr[i][0]);
        gZDCRatioSyst[i] = new TGraphErrors(gNumNchSyst[i]->GetN(), &ZDC[i][0], &RatioSyst[i][0], &ZDCErrSyst[i][0], &RatioSystErr[i][0]);
    }
    prepgraph(gNchRatioStat[0], kFullDiamond, 3., kBlack);
    prepgraph(gNchRatioStat[1], kFullCircle, 2.3, kRed);
    prepgraph(gNchRatioStat[2], kOpenCircle, 2.3, kRed);
    prepgraph(gNchRatioStat[3], kFullSquare, 2.3, kBlue);
    prepgraph(gNchRatioStat[4], kOpenSquare, 2.3, kBlue);
    prepgraph(gZDCRatioStat[0], kFullDiamond, 3., kBlack);
    prepgraph(gZDCRatioStat[1], kFullCircle, 2.3, kRed);
    prepgraph(gZDCRatioStat[2], kOpenCircle, 2.3, kRed);
    prepgraph(gZDCRatioStat[3], kFullSquare, 2.3, kBlue);
    prepgraph(gZDCRatioStat[4], kOpenSquare, 2.3, kBlue);
    prepgraph(gNchRatioSyst[0], kFullDiamond,0, kBlack);
    prepgraph(gNchRatioSyst[1], kFullCircle, 0, kRed);
    prepgraph(gNchRatioSyst[2], kOpenCircle, 0, kRed);
    prepgraph(gNchRatioSyst[3], kFullSquare, 0, kBlue);
    prepgraph(gNchRatioSyst[4], kOpenSquare, 0, kBlue);
    prepgraph(gZDCRatioSyst[0], kFullDiamond,0, kBlack);
    prepgraph(gZDCRatioSyst[1], kFullCircle, 0, kRed);
    prepgraph(gZDCRatioSyst[2], kOpenCircle, 0, kRed);
    prepgraph(gZDCRatioSyst[3], kFullSquare, 0, kBlue);
    prepgraph(gZDCRatioSyst[4], kOpenSquare, 0, kBlue);

    Int_t colors[] = {kRed + 1, kOrange + 10, kOrange + 1, kYellow + 1, kSpring - 1, kAzure + 7, kBlue - 4, kBlue + 3};

    const int nsel1 = gNchRatioStat[1]->GetN();
    const int nsel2 = gNchRatioStat[2]->GetN();
    const int nsel3 = gNchRatioStat[3]->GetN();
    const int nsel4 = gNchRatioStat[4]->GetN();

    TGraphErrors *gNchRatioStatClone_sel1[nsel1], *gZDCRatioStatClone_sel1[nsel1], *gNchRatioSystClone_sel1[nsel1], *gZDCRatioSystClone_sel1[nsel1];
    TGraphErrors *gNchRatioStatClone_sel2[nsel2], *gZDCRatioStatClone_sel2[nsel2], *gNchRatioSystClone_sel2[nsel2], *gZDCRatioSystClone_sel2[nsel2];
    TGraphErrors *gNchRatioStatClone_sel3[nsel3], *gZDCRatioStatClone_sel3[nsel3], *gNchRatioSystClone_sel3[nsel3], *gZDCRatioSystClone_sel3[nsel3];
    TGraphErrors *gNchRatioStatClone_sel4[nsel4], *gZDCRatioStatClone_sel4[nsel4], *gNchRatioSystClone_sel4[nsel4], *gZDCRatioSystClone_sel4[nsel4];

    // sel1
    for (int i = 0; i < nsel1; i++)
    {
        gNchRatioStatClone_sel1[i] = (TGraphErrors *)gNchRatioStat[1]->Clone(Form("gNchRatioStat_sel1%i", i));
        gZDCRatioStatClone_sel1[i] = (TGraphErrors *)gZDCRatioStat[1]->Clone(Form("gZDCRatioStat_sel1%i", i));
        gNchRatioSystClone_sel1[i] = (TGraphErrors *)gNchRatioSyst[1]->Clone(Form("gNchRatioSyst_sel1%i", i));
        gZDCRatioSystClone_sel1[i] = (TGraphErrors *)gZDCRatioSyst[1]->Clone(Form("gZDCRatioSyst_sel1%i", i));

        gNchRatioStatClone_sel1[i]->SetLineColor(colors[i]);
        gNchRatioStatClone_sel1[i]->SetMarkerColor(colors[i]);
        gNchRatioStatClone_sel1[i]->SetMarkerStyle(kFullCircle);
        gNchRatioStatClone_sel1[i]->SetMarkerSize(2.8);
        gNchRatioSystClone_sel1[i]->SetLineColor(colors[i]);
        gNchRatioSystClone_sel1[i]->SetMarkerColor(colors[i]);
        gNchRatioSystClone_sel1[i]->SetFillStyle(0);
        gZDCRatioStatClone_sel1[i]->SetLineColor(colors[i]);
        gZDCRatioStatClone_sel1[i]->SetMarkerColor(colors[i]);
        gZDCRatioStatClone_sel1[i]->SetMarkerStyle(kFullCircle);
        gZDCRatioStatClone_sel1[i]->SetMarkerSize(2.8);
        gZDCRatioSystClone_sel1[i]->SetLineColor(colors[i]);
        gZDCRatioSystClone_sel1[i]->SetMarkerColor(colors[i]);
        gZDCRatioSystClone_sel1[i]->SetFillStyle(0);
        for (int j = 0; j < nsel1; j++)
        {
            if (j != i)
            {
                gNchRatioStatClone_sel1[i]->SetPoint(j, 0., 0.);
                gNchRatioSystClone_sel1[i]->SetPoint(j, 0., 0.);
                gZDCRatioStatClone_sel1[i]->SetPoint(j, 0., 0.);
                gZDCRatioSystClone_sel1[i]->SetPoint(j, 0., 0.);
            }
        }
    }
    // sel2
    for (int i = 0; i < nsel2; i++)
    {
        gNchRatioStatClone_sel2[i] = (TGraphErrors *)gNchRatioStat[2]->Clone(Form("gNchRatioStat_sel2%i", i));
        gZDCRatioStatClone_sel2[i] = (TGraphErrors *)gZDCRatioStat[2]->Clone(Form("gZDCRatioStat_sel2%i", i));
        gNchRatioSystClone_sel2[i] = (TGraphErrors *)gNchRatioSyst[2]->Clone(Form("gNchRatioSyst_sel2%i", i));
        gZDCRatioSystClone_sel2[i] = (TGraphErrors *)gZDCRatioSyst[2]->Clone(Form("gZDCRatioSyst_sel2%i", i));

        gNchRatioStatClone_sel2[i]->SetLineColor(colors[i]);
        gNchRatioStatClone_sel2[i]->SetMarkerColor(colors[i]);
        gNchRatioStatClone_sel2[i]->SetMarkerStyle(kOpenCircle);
        gNchRatioStatClone_sel2[i]->SetMarkerSize(2.8);
        gNchRatioSystClone_sel2[i]->SetLineColor(colors[i]);
        gNchRatioSystClone_sel2[i]->SetMarkerColor(colors[i]);
        gNchRatioSystClone_sel2[i]->SetFillStyle(0);
        gZDCRatioStatClone_sel2[i]->SetLineColor(colors[i]);
        gZDCRatioStatClone_sel2[i]->SetMarkerColor(colors[i]);
        gZDCRatioStatClone_sel2[i]->SetMarkerStyle(kOpenCircle);
        gZDCRatioStatClone_sel2[i]->SetMarkerSize(2.8);
        gZDCRatioSystClone_sel2[i]->SetLineColor(colors[i]);
        gZDCRatioSystClone_sel2[i]->SetMarkerColor(colors[i]);
        gZDCRatioSystClone_sel2[i]->SetFillStyle(0);
        for (int j = 0; j < nsel2; j++)
        {
            if (j != i)
            {
                gNchRatioStatClone_sel2[i]->SetPoint(j, 0., 0.);
                gNchRatioSystClone_sel2[i]->SetPoint(j, 0., 0.);
                gZDCRatioStatClone_sel2[i]->SetPoint(j, 0., 0.);
                gZDCRatioSystClone_sel2[i]->SetPoint(j, 0., 0.);
            }
        }
    }
    // sel3
    for (int i = 0; i < nsel3; i++)
    {
        gNchRatioStatClone_sel3[i] = (TGraphErrors *)gNchRatioStat[3]->Clone(Form("gNchRatioStat_sel3%i", i));
        gZDCRatioStatClone_sel3[i] = (TGraphErrors *)gZDCRatioStat[3]->Clone(Form("gZDCRatioStat_sel3%i", i));
        gNchRatioSystClone_sel3[i] = (TGraphErrors *)gNchRatioSyst[3]->Clone(Form("gNchRatioSyst_sel3%i", i));
        gZDCRatioSystClone_sel3[i] = (TGraphErrors *)gZDCRatioSyst[3]->Clone(Form("gZDCRatioSyst_sel3%i", i));

        gNchRatioStatClone_sel3[i]->SetLineColor(colors[i]);
        gNchRatioStatClone_sel3[i]->SetMarkerColor(colors[i]);
        gNchRatioStatClone_sel3[i]->SetMarkerStyle(kFullSquare);
        gNchRatioStatClone_sel3[i]->SetMarkerSize(2.8);
        gNchRatioSystClone_sel3[i]->SetLineColor(colors[i]);
        gNchRatioSystClone_sel3[i]->SetMarkerColor(colors[i]);
        gNchRatioSystClone_sel3[i]->SetFillStyle(0);
        gZDCRatioStatClone_sel3[i]->SetLineColor(colors[i]);
        gZDCRatioStatClone_sel3[i]->SetMarkerColor(colors[i]);
        gZDCRatioStatClone_sel3[i]->SetMarkerStyle(kFullSquare);
        gZDCRatioStatClone_sel3[i]->SetMarkerSize(2.8);
        gZDCRatioSystClone_sel3[i]->SetLineColor(colors[i]);
        gZDCRatioSystClone_sel3[i]->SetMarkerColor(colors[i]);
        gZDCRatioSystClone_sel3[i]->SetFillStyle(0);
        for (int j = 0; j < nsel3; j++)
        {
            if (j != i)
            {
                gNchRatioStatClone_sel3[i]->SetPoint(j, 0., 0.);
                gNchRatioSystClone_sel3[i]->SetPoint(j, 0., 0.);
                gZDCRatioStatClone_sel3[i]->SetPoint(j, 0., 0.);
                gZDCRatioSystClone_sel3[i]->SetPoint(j, 0., 0.);
            }
        }
    }
    // sel4
    for (int i = 0; i < nsel4; i++)
    {
        gNchRatioStatClone_sel4[i] = (TGraphErrors *)gNchRatioStat[4]->Clone(Form("gNchRatioStat_sel4%i", i));
        gZDCRatioStatClone_sel4[i] = (TGraphErrors *)gZDCRatioStat[4]->Clone(Form("gZDCRatioStat_sel4%i", i));
        gNchRatioSystClone_sel4[i] = (TGraphErrors *)gNchRatioSyst[4]->Clone(Form("gNchRatioSyst_sel4%i", i));
        gZDCRatioSystClone_sel4[i] = (TGraphErrors *)gZDCRatioSyst[4]->Clone(Form("gZDCRatioSyst_sel4%i", i));

        gNchRatioStatClone_sel4[i]->SetLineColor(colors[i]);
        gNchRatioStatClone_sel4[i]->SetMarkerColor(colors[i]);
        gNchRatioStatClone_sel4[i]->SetMarkerStyle(kOpenSquare);
        gNchRatioStatClone_sel4[i]->SetMarkerSize(2.8);
        gNchRatioSystClone_sel4[i]->SetLineColor(colors[i]);
        gNchRatioSystClone_sel4[i]->SetMarkerColor(colors[i]);
        gNchRatioSystClone_sel4[i]->SetFillStyle(0);
        gZDCRatioStatClone_sel4[i]->SetLineColor(colors[i]);
        gZDCRatioStatClone_sel4[i]->SetMarkerColor(colors[i]);
        gZDCRatioStatClone_sel4[i]->SetMarkerStyle(kOpenSquare);
        gZDCRatioStatClone_sel4[i]->SetMarkerSize(2.8);
        gZDCRatioSystClone_sel4[i]->SetLineColor(colors[i]);
        gZDCRatioSystClone_sel4[i]->SetMarkerColor(colors[i]);
        gZDCRatioSystClone_sel4[i]->SetFillStyle(0);
        for (int j = 0; j < nsel4; j++)
        {
            if (j != i)
            {
                gNchRatioStatClone_sel4[i]->SetPoint(j, 0., 0.);
                gNchRatioSystClone_sel4[i]->SetPoint(j, 0., 0.);
                gZDCRatioStatClone_sel4[i]->SetPoint(j, 0., 0.);
                gZDCRatioSystClone_sel4[i]->SetPoint(j, 0., 0.);
            }
        }
    }

    TH1D *hnch = new TH1D("hnch", "", 10, 0, 4.3);
    hnch->SetStats(0);
    hnch->GetXaxis()->SetTitle("n_{ch} / #LT n_{ch} #GT_{MB} ");
    hnch->GetXaxis()->SetTitleSize(0.06);
    hnch->GetYaxis()->SetTitleSize(0.06);
    hnch->GetYaxis()->SetLabelSize(0.04);
    hnch->GetXaxis()->SetLabelSize(0.04);
    hnch->GetYaxis()->SetTitleOffset(1.35);
    if (num.Contains("Xi") && den.Contains("K0Short")) hnch->GetYaxis()->SetRangeUser(0.55, 1.25);
    if (num.Contains("Xi") && den.Contains("Lambda")) hnch->GetYaxis()->SetRangeUser(0.65, 1.25);
    if (num.Contains("Lambda")) hnch->GetYaxis()->SetRangeUser(0.7, 1.2);

    hnch->SetTitle("");
    if (num.Contains("Xi") && den.Contains("K0Short"))
    {
        // hnch->GetYaxis()->SetTitle("#frac{#Xi}{#LT #Xi #GT_{MB}} / #frac{K^{0}_{S}}{#LT K^{0}_{S} #GT_{MB}}");
        hnch->GetYaxis()->SetTitle("#frac{#Xi / #LT #Xi #GT_{MB}}{K^{0}_{S} / #LT K^{0}_{S} #GT_{MB}}");
    }
    if (num.Contains("Xi") && den.Contains("Lambda"))
    {
        // hnch->GetYaxis()->SetTitle("#frac{#Xi}{#LT #Xi #GT_{MB}} / #frac{#Lambda}{#LT #Lambda #GT_{MB}}");
        hnch->GetYaxis()->SetTitle("#frac{#Xi / #LT #Xi #GT_{MB}}{#Lambda / #LT #Lambda #GT_{MB}}");
    }
    if (num.Contains("Lambda") && den.Contains("K0Short"))
    {
        // hnch->GetYaxis()->SetTitle("#frac{#Xi}{#LT #Xi #GT_{MB}} / #frac{K^{0}_{S}}{#LT K^{0}_{S} #GT_{MB}}");
        hnch->GetYaxis()->SetTitle("#frac{#Lambda / #LT #Lambda #GT_{MB}}{K^{0}_{S} / #LT K^{0}_{S} #GT_{MB}}");
    }

    TH1D *hzn = new TH1D("hzn", "", 10, 0.1, 1.5);
    hzn->SetStats(0);
    hzn->GetXaxis()->SetTitle("ZN / #LT ZN #GT_{MB}");
    hzn->GetYaxis()->SetTitle("");
    hzn->GetXaxis()->SetTitleSize(0.06);
    hzn->GetYaxis()->SetTitleSize(0.06);
    hzn->GetYaxis()->SetLabelSize(0.04);
    hzn->GetXaxis()->SetLabelSize(0.04);
    hzn->GetYaxis()->SetTitleOffset(1.35);
    if (num.Contains("Xi") && den.Contains("K0Short"))
        hzn->GetYaxis()->SetRangeUser(0.55, 1.25);
    if (num.Contains("Xi") && den.Contains("Lambda"))
        hzn->GetYaxis()->SetRangeUser(0.65, 1.25);
    if (num.Contains("Lambda")) hzn->GetYaxis()->SetRangeUser(0.7, 1.2);
    hzn->SetTitle("");
    hzn->GetYaxis()->SetLabelColor(kWhite);

    TCanvas *c = new TCanvas("c", "", 1800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 0.54, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0.54, 0, 1., 1);
    prepcanvas(c);
    preppad(pad1, pad2);
    pad1->Draw();
    pad2->Draw();

    //
    pad1->cd();
    hnch->Draw();
    gNchRatioSyst[0]->Draw("SAME E2");
    gNchRatioStat[0]->Draw("SAME EP");
    for (int j = 0; j < nsel1; j++)
    {
        gNchRatioStatClone_sel1[j]->Draw("SAME EP");
        gNchRatioSystClone_sel1[j]->Draw("SAME E2");
    }
    for (int j = 0; j < nsel2; j++)
    {
        gNchRatioStatClone_sel2[j]->Draw("SAME EP");
        gNchRatioSystClone_sel2[j]->Draw("SAME E2");
    }

    TMarker *fullcircle = new TMarker(12., 0.0074509863, 2);
    fullcircle->SetMarkerStyle(kFullCircle);
    fullcircle->SetMarkerSize(2.5);
    fullcircle->Draw("SAME");
    //
    TMarker *opencircle = new TMarker(12., 0.0074509863, 2);
    opencircle->SetMarkerStyle(kOpenCircle);
    opencircle->SetMarkerSize(2.5);
    opencircle->Draw("SAME");

    TLegend *lXi1 = new TLegend(0.65, 0.18, 0.87, 0.47);
    lXi1->SetBorderSize(0);
    lXi1->AddEntry(gNchRatioStat[0], "V0M standalone", "P");
    lXi1->AddEntry(fullcircle, "SPDcl fixed [10#font[122]{-}20]%", "P");
    lXi1->AddEntry(opencircle, "SPDcl fixed [40#font[122]{-}50]%", "P");
    lXi1->AddEntry(gNchRatioStatClone_sel1[0], "high V0M activity", "LEP");
    lXi1->AddEntry(gNchRatioStatClone_sel1[nsel1-1], "low V0M activity", "LEP");
    lXi1->SetTextSize(0.03);
    lXi1->SetTextFont(42);
    lXi1->Draw("SAME");

    TLatex *rap = new TLatex();
    rap->SetTextFont(42);
    rap->SetNDC();
    rap->SetTextColor(1);
    rap->SetTextSize(0.04);
    rap->SetTextAlign(22);
    rap->SetTextAngle(0);
    rap->DrawLatex(0.75, 0.5, Form("%s", "|#it{y}|<0.5"));

    //
    pad2->cd();
    hzn->Draw();
    gZDCRatioSyst[0]->Draw("SAME E2");
    gZDCRatioStat[0]->Draw("SAME EP");

    for (int j = 0; j < nsel1; j++)
    {
        gZDCRatioStatClone_sel1[j]->Draw("SAME EP");
        gZDCRatioSystClone_sel1[j]->Draw("SAME E2");
    }
    for (int j = 0; j < nsel2; j++)
    {
        gZDCRatioStatClone_sel2[j]->Draw("SAME EP");
        gZDCRatioSystClone_sel2[j]->Draw("SAME E2");
    }

    TCanvas *c1 = new TCanvas("c1", "", 1800, 800);
    TPad *pad1_ = new TPad("pad1_", "pad1_", 0., 0., 0.54, 1);
    TPad *pad2_ = new TPad("pad2_", "pad2_", 0.54, 0, 1., 1);
    prepcanvas(c1);
    preppad(pad1_, pad2_);
    pad1_->Draw();
    pad2_->Draw();

    //
    pad1_->cd();
    hnch->Draw();
    gNchRatioSyst[0]->Draw("SAME E2");
    gNchRatioStat[0]->Draw("SAME EP");

    for (int j = 0; j < nsel3; j++)
    {
        gNchRatioStatClone_sel3[j]->Draw("SAME EP");
        gNchRatioSystClone_sel3[j]->Draw("SAME E2");
    }
    for (int j = 0; j < nsel4; j++)
    {
        gNchRatioStatClone_sel4[j]->Draw("SAME EP");
        gNchRatioSystClone_sel4[j]->Draw("SAME E2");
    }
    TLegend *lXi2 = new TLegend(0.65, 0.18, 0.87, 0.47);
    lXi2->SetBorderSize(0);
    lXi2->AddEntry(gNchRatioStat[0], "V0M standalone", "P");
    lXi2->AddEntry(fullcircle, "V0M fixed [10#font[122]{-}20]%", "P");
    lXi2->AddEntry(opencircle, "V0M fixed [40#font[122]{-}50]%", "P");
    lXi2->AddEntry(gNchRatioStatClone_sel3[0], "high SPDcl activity", "LEP");
    lXi2->AddEntry(gNchRatioStatClone_sel3[nsel3 - 1], "low SPDcl activity", "LEP");
    lXi2->SetTextSize(0.03);
    lXi2->SetTextFont(42);
    lXi2->Draw("SAME");
    rap->DrawLatex(0.75, 0.56, Form("%s", "|#it{y}|<0.5"));
    //
    pad2_->cd();
    hzn->Draw();
    gZDCRatioSyst[0]->Draw("SAME E2");
    gZDCRatioStat[0]->Draw("SAME EP");

    for (int j = 0; j < nsel3; j++)
    {
        gZDCRatioStatClone_sel3[j]->Draw("SAME EP");
        gZDCRatioSystClone_sel3[j]->Draw("SAME E2");
    }
    for (int j = 0; j < nsel4; j++)
    {
        gZDCRatioStatClone_sel4[j]->Draw("SAME EP");
        gZDCRatioSystClone_sel4[j]->Draw("SAME E2");
    }


    c->SaveAs(Form("Ratio%sover%s_SPDfixed.pdf", num.Data(), den.Data()));
    c1->SaveAs(Form("Ratio%sover%s_V0Mfixed.pdf", num.Data(), den.Data()));
    c->SaveAs(Form("Ratio%sover%s_SPDfixed.png", num.Data(), den.Data()));
    c1->SaveAs(Form("Ratio%sover%s_V0Mfixed.png", num.Data(), den.Data()));
}

void prepcanvas(TCanvas* c){
    c->SetLeftMargin(0.2);
    c->SetBottomMargin(0.18);
    c->SetRightMargin(0.1);
    c->SetTopMargin(0.1);
    c->SetTicky();
    c->SetTickx();
}

void preppad(TPad *pad1, TPad *pad2){
    pad1->SetBorderMode(0);
    pad1->SetTopMargin(0.05);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.15);
    pad1->SetRightMargin(0.0);
    pad2->SetRightMargin(0.05);
    pad1->SetLeftMargin(0.2);
    pad2->SetLeftMargin(0.0);
    pad1->SetBottomMargin(0.15);
    pad2->SetBorderMode(0);
    pad1->SetFillColor(kWhite);
    pad1->SetTicky();
    pad1->SetTickx();
    pad2->SetFillColor(kWhite);
    pad2->SetTicky();
    pad2->SetTickx();
}

void prepgraph(TGraphErrors* g, int marker, float size, int color){
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(size);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetFillStyle(0);
}