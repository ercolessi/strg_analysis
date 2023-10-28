void histobeauty(TH1D *h);
void histobeautyratio(TH1D *h);

void draweff(TString part = "XiMinus"){

    Double_t percentileV0M_low[4] = {0, 0, 30, 70};
    Double_t percentileV0M_high[4] = {100, 10, 40, 100};

    TFile *f[4];
    if (part.Contains("K0Short"))
    {
        f[0] = TFile::Open("Results-K0Short-13TeV-SPDClusters_000_100_V0M_000_100_FDNoFD.root");
        f[1] = TFile::Open("Results-K0Short-13TeV-SPDClusters_000_010_V0M_000_100_FDNoFD.root");
        f[2] = TFile::Open("Results-K0Short-13TeV-SPDClusters_030_040_V0M_000_100_FDNoFD.root");
        f[3] = TFile::Open("Results-K0Short-13TeV-SPDClusters_070_100_V0M_000_100_FDNoFD.root");
    }
    else if (part.Contains("Lambda"))
    {
        f[0] = TFile::Open("Results-Lambda-13TeV-SPDClusters_000_100_V0M_000_100_FDNoFD.root");
        f[1] = TFile::Open("Results-Lambda-13TeV-SPDClusters_000_010_V0M_000_100_FDNoFD.root");
        f[2] = TFile::Open("Results-Lambda-13TeV-SPDClusters_030_040_V0M_000_100_FDNoFD.root");
        f[3] = TFile::Open("Results-Lambda-13TeV-SPDClusters_070_100_V0M_000_100_FDNoFD.root");
    }
    else if (part.Contains("XiMinus"))
    {
        f[0] = TFile::Open("Results-XiMinus-13TeV-SPDClusters_000_100_V0M_000_100.root");
        f[1] = TFile::Open("Results-XiMinus-13TeV-SPDClusters_000_010_V0M_000_100.root");
        f[2] = TFile::Open("Results-XiMinus-13TeV-SPDClusters_030_040_V0M_000_100.root");
        f[3] = TFile::Open("Results-XiMinus-13TeV-SPDClusters_070_100_V0M_000_100.root");
    }
    else return;

    Int_t colors[] = {kBlack, kRed, kGreen + 1, kBlue};

    TH1D *h[4];
    for (int i = 0; i < 4; i++){
        if (part.Contains("K0Short")) h[i] = (TH1D *)f[i]->Get("fHistPureEfficiency");
        else if (part.Contains("Lambda")) h[i] = (TH1D *)f[i]->Get("fHistEfficiency");
        else if (part.Contains("XiMinus")) h[i] = (TH1D *)f[i]->Get("fHistEfficiency");
        else return;
        histobeauty(h[i]);
        h[i]->SetLineColor(colors[i]);
        h[i]->SetMarkerColor(colors[i]);
    }

    TCanvas *c = new TCanvas("c", "c", 1000, 1300);
    c->SetLeftMargin(0.05);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.05);
    c->SetTopMargin(0.05);
    c->SetTicks();

    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
    pad1->SetBottomMargin(0.01);
    pad1->SetTopMargin(0.05);
    pad1->SetTicks();
    pad1->Draw();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1, 0.3);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->SetBottomMargin(0.25);
    pad2->SetTopMargin(0.01);
    pad2->SetTicks();
    pad2->Draw();

    pad1->cd();

    TH1D* hratio[4];
    for (int i = 0; i < 4; i++) {
        h[i]->Draw("SAME");
        hratio[i] = (TH1D*)h[i]->Clone(Form("hratio_%d", i));
        hratio[i]->Divide(h[0]);
    }

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

    TLegend *leg = new TLegend(0.6, 0.12, 0.88, 0.4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetHeader("VZERO classes");
    for (int i = 0; i < 4; i++)
    {
        leg->AddEntry(h[i], Form("%d-%d %%", (int)percentileV0M_low[i], (int)percentileV0M_high[i]), "lep");
    }
    leg->Draw("SAME");

    pad2->cd();
    for (int i = 0; i < 4; i++) {
        histobeautyratio(hratio[i]);
        hratio[i]->SetLineColor(colors[i]);
        hratio[i]->SetMarkerColor(colors[i]);
        hratio[i]->Draw("SAME");
    }

    c->SaveAs(Form("EffMult_%s.pdf", part.Data()));
}

void histobeauty(TH1D * h)
{
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("#varepsilon");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.07);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetRangeUser(0, 1.5 * h->GetMaximum());
    h->SetTitle("");
    h->SetStats(0);
}

void histobeautyratio(TH1D *h)
{
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("Ratio to 0-100%");
    h->GetXaxis()->SetTitleSize(0.09);
    h->GetYaxis()->SetTitleSize(0.1);
   // h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(0.6);
    h->GetXaxis()->SetLabelSize(0.08);
    h->GetYaxis()->SetLabelSize(0.08);
    h->GetYaxis()->SetRangeUser(0.5, 1.52);
    h->SetTitle("");
    h->SetStats(0);
}