void histobeauty(TH1D* h);

void drawSgnBkg(TString part = "K0Short"){

    Float_t percentileV0M_low[] = {0,0,1,5,10,15,20,30,40,50,70};
    Float_t percentileV0M_high[] = {100,1,5,10,15,20,30,40,50,70,100};
    Float_t percentileSPDCl_low[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    Float_t percentileSPDCl_high[] = {100,100,100,100,100,100,100,100,100,100,100,100};

    /*Float_t percentileV0M_low[] = { 30, 30, 20, 0,
                                 20, 10, 0,
                                 10, 20, 30, 40,
                                 50, 30, 40,
                                 50, 60, 70};
    Float_t percentileV0M_high[] = { 70, 50, 50, 30,
                                  30, 30, 20,
                                  20, 30, 40, 50,
                                  100, 40, 50,
                                  60, 70, 100};

    Float_t percentileSPDCl_low[] = { 10, 20, 30, 50,
                                 0, 10, 20,
                                 10, 10, 10, 10,
                                 10, 40, 40,
                                 40, 40, 40};
    Float_t percentileSPDCl_high[] = { 30, 40, 50, 100,
                                  10, 20, 30,
                                  20, 20, 20, 20,
                                  20, 50, 50,
                                  50, 50, 50};*/

    Int_t colors[] = {kBlack, kRed + 2, kRed, kOrange + 7, kYellow + 1, kSpring - 1, kGreen + 2, kTeal, kAzure + 7, kBlue, kBlue + 2,
                      kBlack, kRed + 2, kRed, kOrange + 7, kYellow + 1, kSpring - 1, kGreen + 2, kTeal, kAzure + 7, kBlue, kBlue + 2};

    const int nbins = sizeof(percentileV0M_low) / sizeof(Float_t) ;

    TFile *f[nbins];
    for (int i = 0; i < nbins; i++) {
        if (part.Contains("Lambda"))
        {
            if (part.Contains("Anti"))
                f[i] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/resultsFD/Results-AntiLambda-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_FDUseMCRatio.root", percentileSPDCl_low[i], percentileSPDCl_high[i], percentileV0M_low[i], percentileV0M_high[i]));
            else f[i] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/resultsFD/Results-Lambda-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_FDUseMCRatio.root", percentileSPDCl_low[i], percentileSPDCl_high[i], percentileV0M_low[i], percentileV0M_high[i]));

        }
        else if (part.Contains("K0Short"))
        {
            f[i] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/results/Results-K0Short-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_FDUseMCRatio.root", percentileSPDCl_low[i], percentileSPDCl_high[i], percentileV0M_low[i], percentileV0M_high[i]));
        }
        else if (part.Contains("XiMinus"))
        {
            f[i] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/CascadeAnalysis/results/Results-XiMinus-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root",percentileSPDCl_low[i], percentileSPDCl_high[i], percentileV0M_low[i], percentileV0M_high[i]));
        }
        else if (part.Contains("XiPlus"))
        {
            f[i] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/CascadeAnalysis/results/Results-XiPlus-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root",percentileSPDCl_low[i], percentileSPDCl_high[i], percentileV0M_low[i], percentileV0M_high[i]));
        }
        else return;
    }

    TH1D *h[nbins];
    TH1D *hclone[nbins];
    TH1D* hraw[nbins];

    for (int i = 0; i < nbins; i++)
    {
        h[i] = (TH1D *)f[i]->Get(Form("lInvMassReal/lInvMassRealRawData/fHistSigToNoise"));
        hraw[i] = (TH1D *)f[i]->Get(Form("lInvMassReal/lInvMassRealRawData/fHistPtRaw"));
        hclone[i] = (TH1D*)h[i]->Clone(Form("hclone_%d",i));
        hclone[i]->Reset();
        for (int b = 1; b <= h[i]->GetNbinsX(); b++)
        {
            if (h[i]->GetBinContent(b) > 0)
                hclone[i]->SetBinContent(b, 1. / ((h[i]->GetBinContent(b) + 1) * 1. / h[i]->GetBinContent(b)));
            if (hclone[i]->GetBinContent(b) == 0) {
                cout << "bin " << b << " is zero" << endl;
                cout << "perc " << percentileSPDCl_low[i] << " " << percentileSPDCl_high[i] << " + " << percentileV0M_low[i] << " " << percentileV0M_high[i] << endl;
            }
            hclone[i]->SetBinError(b,
            //hclone[i]->GetBinContent(b)*(hraw[i]->GetBinError(b)/hraw[i]->GetBinContent(b)));
             0);
            histobeauty(hclone[i]);
            if (part.Contains("K0Short")) {
                hclone[i]->GetXaxis()->SetRangeUser(0.,10.0);
                hclone[i]->GetYaxis()->SetRangeUser(0.7,1.1);
            }
            if (part.Contains("Lambda")) {
                hclone[i]->GetXaxis()->SetRangeUser(0.6,8.0);
                hclone[i]->GetYaxis()->SetRangeUser(0.7,1.1);
            }
            if (part.Contains("AntiLambda")) {
                hclone[i]->GetXaxis()->SetRangeUser(0.6,8.0);
                hclone[i]->GetYaxis()->SetRangeUser(0.7, 1.1);
            }
            if (part.Contains("Xi")) hclone[i]->GetYaxis()->SetRangeUser(0.5,1.1);
            hclone[i]->SetLineColor(colors[i]);
            hclone[i]->SetMarkerColor(colors[i]);
        }
    }

    TCanvas* c = new TCanvas("c", "c", 900, 800);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.15);
    c->SetTopMargin(0.05);
    c->SetTicks();

    TLegend *leg = new TLegend(0.3, 0.18, 0.75, 0.5);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetNColumns(2);
    leg->SetHeader("VZERO classes");
    for (int i = 0; i < nbins; i++)
    {
        hclone[i]->Draw("SAME");
        leg->AddEntry(hclone[i], Form("%d-%d %%", (int)percentileV0M_low[i], (int)percentileV0M_high[i]), "lep");
    }
    leg->Draw("SAME");

    TLatex *ltx = new TLatex();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.035);
    ltx->SetTextAlign(12);
    TLatex *ltxbig = new TLatex();
    ltxbig->SetTextFont(42);
    ltxbig->SetTextSize(0.055);
    ltxbig->SetTextAlign(12);
    if (part.Contains("K0Short"))
    {
        ltxbig->DrawLatexNDC(0.65, 0.88, "K^{0}_{S} #rightarrow #pi^{+} + #pi^{-}");
    }
    if (part.Contains("AntiLambda"))
    {
        ltxbig->DrawLatexNDC(0.65, 0.88, "#bar{#Lambda} #rightarrow #bar{p} + #pi^{+}");
    }
    else if (part.Contains("Lambda"))
    {
        ltxbig->DrawLatexNDC(0.65, 0.88, "#Lambda #rightarrow p + #pi^{-}");
    }
    if (part.Contains("XiMinus"))
    {
        ltxbig->DrawLatexNDC(0.65, 0.88, "#Xi^{-} #rightarrow #Lambda + #pi^{-}");
    }
    if (part.Contains("XiPlus"))
    {
        ltxbig->DrawLatexNDC(0.65, 0.88, "#bar{#Xi}^{+} #rightarrow #bar{#Lambda} + #pi^{+}");
    }

    TLatex *tex = new TLatex(0.2, 0.88, "ALICE, pp #sqrt{s} = 13 TeV");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    TLatex *tex2 = new TLatex(0.2, 0.82, "This work");
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(2);
    tex->Draw();
    tex2->Draw();



    c->SaveAs(Form("SgnToNoise_%s.pdf",part.Data()));
}

void histobeauty(TH1D* h){
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("S/(S+B)");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.3);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetRangeUser(0, 1.2 * h->GetMaximum());
    h->SetLineWidth(2);
    h->SetTitle("");
    h->SetStats(0);
}