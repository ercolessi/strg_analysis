void Init(Int_t lClassCode,
          vector<Double_t> &percentileSPDCl_low,
          vector<Double_t> &percentileSPDCl_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          vector<Int_t> &colors);

void prepcanvas(TCanvas *c);
void preppad1(TPad *pad);
void preppad2(TPad *pad);
void ratiohisto(TH1D *h);
void histonch(TH1D *h);
void histozn(TH1D* h);

void plotLog10(Int_t lClassCode = 1)
{

    vector<Double_t> percentileSPDCl_low;
    vector<Double_t> percentileSPDCl_high;
    vector<Double_t> percentileV0M_low;
    vector<Double_t> percentileV0M_high;
    vector<Int_t> colors;

    Init(lClassCode, percentileSPDCl_low, percentileSPDCl_high, percentileV0M_low, percentileV0M_high, colors);

    TFile *filecontainer = TFile::Open("../../PhDFinal/Selections/NchRawContainer_reference.root");
    if (filecontainer == 0x0) {
        cout << "File not found" << endl;
        return;
    }

    TH3D *hnch = (TH3D *)filecontainer->Get("hspd_spdv0m"); // x -> nch, y -> spd, z -> v0m
    TH3D *hZN = (TH3D *)filecontainer->Get("hznsum_spdv0m"); // x -> zn, y -> spd, z -> v0m

    //INEL>0
    TH1D* pnchMB = hnch->ProjectionX("pnchMB", hnch->GetYaxis()->FindBin(0. + 1e-10), hnch->GetYaxis()->FindBin(100.), hnch->GetZaxis()->FindBin(0. + 1e-10), hnch->GetZaxis()->FindBin(100.));
    TH1D* pznMB = hZN->ProjectionX("pznMB", hZN->GetYaxis()->FindBin(0. + 1e-10), hZN->GetYaxis()->FindBin(100.), hZN->GetZaxis()->FindBin(0. + 1e-10), hZN->GetZaxis()->FindBin(100.));
    Double_t meanZNMB = pznMB->GetMean();
    Double_t meanNchMB = pnchMB->GetMean();

    //Other classes
    const int sizePerc = percentileSPDCl_low.size();
    cout << "N percentiles for class " << lClassCode << " = " << sizePerc << endl;
    for (int i = 0; i<sizePerc; i++){
        cout << "SPDCl " << percentileSPDCl_low[i] << "-" << percentileSPDCl_high[i] << " + "
        << "V0M " << percentileV0M_low[i] << "-" << percentileV0M_high[i] << endl;
    }

    TH1D *pnch[sizePerc], *pZN[sizePerc], *pspd0815[sizePerc];
    TH1D* pZNLog10[sizePerc];
    double meanZN[sizePerc], meanNch[sizePerc];
    double miny, maxy, minz, maxz;

    for (int icent = 0; icent < sizePerc; icent++)
    {
        miny = hnch->GetYaxis()->FindBin(percentileSPDCl_low[icent] + 1E-10);
        maxy = hnch->GetYaxis()->FindBin(percentileSPDCl_high[icent] - 1E-10);
        minz = hnch->GetZaxis()->FindBin(percentileV0M_low[icent] + 1E-10);
        maxz = hnch->GetZaxis()->FindBin(percentileV0M_high[icent] - 1E-10);
        pnch[icent] = hnch->ProjectionX(Form("pnch%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDCl_low[icent], percentileSPDCl_high[icent], "V0M", percentileV0M_low[icent], percentileV0M_high[icent]), miny, maxy, minz, maxz);
        pZN[icent] = hZN->ProjectionX(Form("pZN%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDCl_low[icent], percentileSPDCl_high[icent], "V0M", percentileV0M_low[icent], percentileV0M_high[icent]), miny, maxy, minz, maxz);
        pZNLog10[icent] = new TH1D(Form("pZNLog10%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDCl_low[icent], percentileSPDCl_high[icent], "V0M", percentileV0M_low[icent], percentileV0M_high[icent]), "",
                                    30,1,4);

        for (int b = 1; b <= pZN[icent]->GetNbinsX(); b++)
        {
            double temp = pZN[icent]->GetBinContent(b);
            for (int i = 0; i < temp; i++)
            {
                pZNLog10[icent]->Fill(TMath::Log10(pZN[icent]->GetBinCenter(b)+0.1));
            }
        }

        // nch
        meanNch[icent] = pnch[icent]->GetMean() / meanNchMB;

        // ZN
        meanZN[icent] = pZN[icent]->GetMean() / meanZNMB;
    }

    TString classe[] = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"};
    //
    TLatex *tex = new TLatex(0.2, 0.89, "ALICE, pp #sqrt{s} = 13 TeV");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    //
    TLatex *tex2 = new TLatex(0.2, 0.84, "This work");
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(2);
    TCanvas *c5 = new TCanvas("c5", "class1", 1300, 800);
    prepcanvas(c5);
    c5->cd(1);
    //
    TPad *pad1 = new TPad("pad1", "", 0, 0.3, 1, 1);
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    preppad1(pad1);
    preppad2(pad2);
    pad1->Draw();
    pad2->Draw();
    //
    pad1->cd();

    pnch[sizePerc-1]->Draw("LEP");
    pnch[sizePerc-1]->Draw("HIST SAME");
    for (int i = 0; i < sizePerc; i++)
    {
        histonch(pnch[i]);
        pnch[i]->SetLineColor(colors[i]);
        pnch[i]->SetMarkerColor(colors[i]);
        if (lClassCode == 1){
            pnch[i]->GetXaxis()->SetRangeUser(0, 25);
        }
        if (lClassCode == 2){
            pnch[i]->GetXaxis()->SetRangeUser(0, 15);
        }
        if (lClassCode == 3) {
            pnch[i]->GetXaxis()->SetRangeUser(0, 25);
        }
        if (lClassCode == 4)
        {
            pnch[i]->GetXaxis()->SetRangeUser(0.0, 28);
        }
        pnch[i]->Draw("LEP SAME");
        pnch[i]->Draw("HIST SAME");
    }

    TLegend *leg5 = new TLegend(0.62, 0.5, 0.8, 0.91);
    leg5->SetBorderSize(0);
    leg5->SetFillColor(0);
    leg5->SetTextSize(0.04);
    leg5->SetHeader("Standalone class:");
    if (lClassCode == 1) leg5->SetHeader("High multiplicity class:");
    if (lClassCode == 2) leg5->SetHeader("Low multiplicity class:");
    if (lClassCode == 3) leg5->SetHeader("High ZN class:");
    if (lClassCode == 4) leg5->SetHeader("Low ZN class:");
    for (int i = 0; i < sizePerc; i++) {
        leg5->AddEntry(pnch[i], Form("%s", classe[i].Data()), "LEP");
    }
    leg5->Draw();
    tex->Draw();
    tex2->Draw();
    //
    pad2->cd();
    if (lClassCode == 0 || lClassCode == 3 || lClassCode == 4)
        pad2->cd()->SetLogy();
    TH1D* rationch[sizePerc];
    int pos = sizePerc/2;
    for (int i = 0; i < sizePerc; i++)
    {
        rationch[i] = (TH1D *)pnch[i]->Clone(Form("rationch%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDCl_low[i], percentileSPDCl_high[i], "V0M", percentileV0M_low[i], percentileV0M_high[i]));
        rationch[i]->Divide(pnch[pos]);
        rationch[i]->SetLineColor(colors[i]);
        ratiohisto(rationch[i]);
        rationch[i]->GetYaxis()->SetTitle(Form("Ratio to %s",classe[pos].Data()));
        rationch[i]->GetXaxis()->SetTitle("tracklets");
        rationch[i]->GetXaxis()->SetRangeUser(0, 30);
        rationch[i]->GetYaxis()->SetRangeUser(0.0001, 3000);
        if (lClassCode == 1){
            rationch[i]->GetXaxis()->SetRangeUser(0, 25);
            rationch[i]->GetYaxis()->SetRangeUser(0.0, 2.65);
        }
        if (lClassCode == 2){
            rationch[i]->GetXaxis()->SetRangeUser(0, 15);
            rationch[i]->GetYaxis()->SetRangeUser(0.0, 2.65);
        }
        if (lClassCode == 3) {
            rationch[i]->GetXaxis()->SetRangeUser(0, 25);
        }
        if (lClassCode == 4){
            rationch[i]->GetXaxis()->SetRangeUser(0.0, 28);
            rationch[i]->GetYaxis()->SetRangeUser(0.001, 50000);
        }

        if (i != pos)
            rationch[i]->Draw("LEP SAME");
    }

    c5->cd(2);

    TPad *pad1_ = new TPad("pad1_", "", 0, 0.3, 1, 1);
    TPad *pad2_ = new TPad("pad2_", "", 0, 0, 1, 0.3);
    preppad1(pad1_);
    preppad2(pad2_);
    pad2_->Draw();
    pad1_->Draw();

    pad1_->cd();
   // pad1_->SetLogy();
    for (int i = 0; i < sizePerc; i++)
    {
        //pZNLog10[i]->Rebin(10);
        histozn(pZNLog10[i]);
        pZNLog10[i]->GetYaxis()->SetRangeUser(0.00004, pZNLog10[i]->GetMaximum() * 18);
        pZNLog10[i]->GetXaxis()->SetRangeUser(1,4);
        if (lClassCode == 2)
        {
            pZNLog10[i]->GetYaxis()->SetRangeUser(0.00008, pZNLog10[i]->GetMaximum() * 1.5);
        }
        if (lClassCode == 3)
        {
            pZNLog10[i]->GetYaxis()->SetRangeUser(0.0002, 0.2);
        }
        if (lClassCode == 4)
        {
            pZNLog10[i]->GetYaxis()->SetRangeUser(0.00008, 0.2);
        }
        pZNLog10[i]->Scale(1. / pZNLog10[i]->Integral(0, -1));
        pZNLog10[i]->SetLineColor(colors[i]);
        pZNLog10[i]->SetMarkerColor(colors[i]);
        pZNLog10[i]->Draw("LEP SAME");
        pZNLog10[i]->Draw("HIST SAME");
    }
    leg5->Draw();
    tex->Draw();
    tex2->Draw();

    pad2_->cd();
    if (lClassCode == 0 || lClassCode == 1 || lClassCode == 2){
        pad2_->SetLogy();
    }
    TH1D *ratiozn[sizePerc];
    for (int i = 0; i < sizePerc; i++)
    {
        ratiozn[i] = (TH1D *)pZN[i]->Clone(Form("ratiozn%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDCl_low[i], percentileSPDCl_high[i], "V0M", percentileV0M_low[i], percentileV0M_high[i]));
        ratiozn[i]->Divide(pZN[pos]);
        ratiozn[i]->SetLineColor(colors[i]);
        ratiozn[i]->GetXaxis()->SetRangeUser(0, 1100);
        ratiohisto(ratiozn[i]);
        ratiozn[i]->GetYaxis()->SetTitle(Form("Ratio to %s", classe[pos].Data()));
        ratiozn[i]->GetYaxis()->SetRangeUser(0.1, 9);
        if (lClassCode == 2)
        {
            ratiozn[i]->GetYaxis()->SetRangeUser(0.2,3);
        }
        if (lClassCode == 4)
        {
            ratiozn[i]->GetYaxis()->SetRangeUser(0., 2.1);
        }
        if (lClassCode == 3)
        {
            ratiozn[i]->GetYaxis()->SetRangeUser(0., 2.1);
        }
        ratiozn[i]->GetXaxis()->SetTitle("ZN (a.u.)");
        if (i != pos)
            ratiozn[i]->Draw("LEP SAME");
    }

   c5->SaveAs(Form("distrLog10nchZN_class%i.pdf",lClassCode));

}

void Init(Int_t lClassCode,
          vector<Double_t> &percentileSPDCl_low,
          vector<Double_t> &percentileSPDCl_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          vector<Int_t> &colors)
{
    // class 0 --> standalone
    Double_t percentileV0M_low_0[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
    Double_t percentileV0M_high_0[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    Double_t percentileSPDCl_low_0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t percentileSPDCl_high_0[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Int_t colors_0[] = {kRed + 2, kRed, kOrange + 7, kYellow + 1, kSpring - 1, kGreen + 2, kTeal, kAzure + 7, kBlue, kBlue + 2};

    Long_t n0 = sizeof(percentileV0M_low_0) / sizeof(Double_t);

    // class 1 --> kHighMult
    Double_t percentileSPDCl_low_1[] = {10, 10, 10, 10, 10, 10, 10};
    Double_t percentileSPDCl_high_1[] = {20, 20, 20, 20, 20, 20, 20};
    Double_t percentileV0M_low_1[] = {0, 5, 10, 20, 30, 40, 50};
    Double_t percentileV0M_high_1[] = {5, 10, 20, 30, 40, 50, 100};
    Int_t colors_1[] = {kRed + 2, kOrange + 7, kYellow + 1, kSpring - 1, kTeal, kAzure + 7, kBlue + 1};
    Long_t n1 = sizeof(percentileSPDCl_low_1) / sizeof(Double_t);

    // class 2 --> kLowMult
    Double_t percentileSPDCl_low_2[] = {40, 40, 40, 40, 40, 40, 40};
    Double_t percentileSPDCl_high_2[] = {50, 50, 50, 50, 50, 50, 50};
    Double_t percentileV0M_low_2[] = {0, 20, 30, 40, 50, 60, 70};
    Double_t percentileV0M_high_2[] = {20, 30, 40, 50, 60, 70, 100};
    Int_t colors_2[] = {kRed + 1, kOrange + 7, kYellow + 1, kSpring - 1, kTeal, kAzure + 7, kBlue + 1};
    Long_t n2 = sizeof(percentileSPDCl_low_2) / sizeof(Double_t);

    // class 3 --> fixed high ZN
    Double_t percentileSPDCl_low_3[] = {0, 10, 20, 30, 50};
    Double_t percentileSPDCl_high_3[] = {20, 30, 40, 50, 100};
    Double_t percentileV0M_low_3[] = {40, 30, 30, 20, 0};   // 40
    Double_t percentileV0M_high_3[] = {60, 70, 50, 50, 30}; // 60
    Int_t colors_3[] = {kRed + 1, kOrange + 7, kSpring - 1, kAzure + 7, kBlue + 1};
    Long_t n3 = sizeof(percentileSPDCl_low_3) / sizeof(Double_t);

    // class 4 --> fixed low ZN
    Double_t percentileSPDCl_low_4[] = {0, 10, 20, 30};
    Double_t percentileSPDCl_high_4[] = {10, 20, 30, 50};
    Double_t percentileV0M_low_4[] = {20, 10, 0, 0};
    Double_t percentileV0M_high_4[] = {30, 30, 20, 10};
    Int_t colors_4[] = {kRed + 1, kOrange + 7, kAzure + 7, kBlue + 1};
    Long_t n4 = sizeof(percentileSPDCl_low_4) / sizeof(Double_t);

    if (lClassCode == 0)
    {
        for (Int_t i = 0; i < n0; i++)
        {
            percentileSPDCl_low.push_back(percentileSPDCl_low_0[i]);
            percentileSPDCl_high.push_back(percentileSPDCl_high_0[i]);
            percentileV0M_low.push_back(percentileV0M_low_0[i]);
            percentileV0M_high.push_back(percentileV0M_high_0[i]);
            colors.push_back(colors_0[i]);
        }
    }

    if (lClassCode == 1)
    {
        for (Int_t i = 0; i < n1; i++)
        {
            percentileSPDCl_low.push_back(percentileSPDCl_low_1[i]);
            percentileSPDCl_high.push_back(percentileSPDCl_high_1[i]);
            percentileV0M_low.push_back(percentileV0M_low_1[i]);
            percentileV0M_high.push_back(percentileV0M_high_1[i]);
            colors.push_back(colors_1[i]);
        }
    }

    if (lClassCode == 2)
    {
        for (Int_t i = 0; i < n2; i++)
        {
            percentileSPDCl_low.push_back(percentileSPDCl_low_2[i]);
            percentileSPDCl_high.push_back(percentileSPDCl_high_2[i]);
            percentileV0M_low.push_back(percentileV0M_low_2[i]);
            percentileV0M_high.push_back(percentileV0M_high_2[i]);
            colors.push_back(colors_2[i]);
        }
    }

    if (lClassCode == 3)
    {
        for (Int_t i = 0; i < n3; i++)
        {
            percentileSPDCl_low.push_back(percentileSPDCl_low_3[i]);
            percentileSPDCl_high.push_back(percentileSPDCl_high_3[i]);
            percentileV0M_low.push_back(percentileV0M_low_3[i]);
            percentileV0M_high.push_back(percentileV0M_high_3[i]);
            colors.push_back(colors_3[i]);
        }
    }

    if (lClassCode == 4)
    {
        for (Int_t i = 0; i < n4; i++)
        {
            percentileSPDCl_low.push_back(percentileSPDCl_low_4[i]);
            percentileSPDCl_high.push_back(percentileSPDCl_high_4[i]);
            percentileV0M_low.push_back(percentileV0M_low_4[i]);
            percentileV0M_high.push_back(percentileV0M_high_4[i]);
            colors.push_back(colors_4[i]);
        }
    }
}

void prepcanvas(TCanvas *c)
{
    c->Divide(2, 1);
}

void preppad1(TPad *pad)
{
    pad->SetBottomMargin(0.001);
    pad->SetLeftMargin(0.12);
    pad->SetRightMargin(0.02);
    pad->SetTopMargin(0.02);
    pad->SetTicks();
}

void preppad2(TPad *pad)
{
    pad->SetBottomMargin(0.3);
    pad->SetLeftMargin(0.12);
    pad->SetRightMargin(0.02);
    pad->SetTopMargin(0.001);
    pad->SetTicks();
}

void ratiohisto(TH1D *h)
{
    // rationch[i]->SetLineWidth(2);
    h->GetYaxis()->SetRangeUser(0., 1.9);
    // h->GetXaxis()->SetRangeUser(0., 30);
    h->GetXaxis()->SetTitleSize(0.1);
    h->GetYaxis()->SetTitleSize(0.08);
    h->GetXaxis()->SetLabelSize(0.08);
    h->GetYaxis()->SetLabelSize(0.08);
    h->GetYaxis()->SetTitleOffset(0.7);
}

void histonch(TH1D *h)
{
    h->SetStats(0);
    h->SetTitle("");
    h->GetYaxis()->SetTitle("Normalized counts");
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->Scale(1. / h->Integral(0, -1));
    // h->GetYaxis()->SetRangeUser(-0.003, h->GetMaximum() * 1.2);
    h->SetMarkerStyle(kFullCircle);
    // h->SetLineWidth(2);
    h->GetXaxis()->SetRangeUser(0, 30);
    h->GetYaxis()->SetRangeUser(-0.005, h->GetMaximum() * 1.3);
}

void histozn(TH1D *h)
{
    h->SetStats(0);
    h->SetTitle("");
    h->GetYaxis()->SetTitle("Normalized counts");
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.2);
    //h->Scale(1. / h->Integral(0, -1));
    // h->GetYaxis()->SetRangeUser(-0.003, h->GetMaximum() * 1.2);
    h->SetMarkerStyle(kFullCircle);
    // h->SetLineWidth(2);
    h->GetXaxis()->SetRangeUser(0, 1100);
}
