void histobeauty(TH1D* h);

void drawpeakposition(TString part = "AntiLambda"){
    TFile *f;
    if (part.Contains("AntiLambda"))
    {
        f = TFile::Open("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/resultsFD/Results-AntiLambda-13TeV-SPDClusters_000_100_V0M_000_100_FDUseMCRatio.root");
    }
    else if (part.Contains("K0Short")) {
        f = TFile::Open("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/results/Results-K0Short-13TeV-SPDClusters_000_100_V0M_000_100_FDNoFD.root");
    }
    else if (part.Contains("Lambda")) {
        f = TFile::Open("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/resultsFD/Results-Lambda-13TeV-SPDClusters_000_100_V0M_000_100_FDUseMCRatio.root");
    }
    else if (part.Contains("XiMinus")) {
        f = TFile::Open("/home/fercoles/strg_analysis/22Sett2023/CascadeAnalysis/CascadeAnalysis/results/Results-XiMinus-13TeV-SPDClusters_000_100_V0M_000_100.root");
    }
    else if (part.Contains("XiPlus")) {
        f = TFile::Open("/home/fercoles/strg_analysis/22Sett2023/CascadeAnalysis/CascadeAnalysis/results/Results-XiPlus-13TeV-SPDClusters_000_100_V0M_000_100.root");
    }
    else return;

    TH1D *h = (TH1D *)f->Get(Form("lInvMassReal/lInvMassRealRawData/fHistPeakPosition"));

    TCanvas* c = new TCanvas("c", "c", 900, 800);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.15);
    c->SetTopMargin(0.05);
    c->SetTicks();
    histobeauty(h);
    if (part.Contains("K0Short")) {
        h->GetYaxis()->SetRangeUser(0.48,0.515);
        h->GetXaxis()->SetRangeUser(0.,10.0);
    }
    if (part.Contains("Lambda")) {
        h->GetYaxis()->SetRangeUser(1.11,1.1225);
        h->GetXaxis()->SetRangeUser(0.6,8.0);
    }
    if (part.Contains("AntiLambda")) {
        h->GetYaxis()->SetRangeUser(1.11,1.1225);
        h->GetXaxis()->SetRangeUser(0.6,8.0);
    }
    if (part.Contains("Xi")) h->GetYaxis()->SetRangeUser(1.315,1.33);
    h->Draw();

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

    c->SaveAs(Form("PeakPos_%s.pdf",part.Data()));
}

void histobeauty(TH1D* h){
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetYaxis()->SetTitle("Peak position (GeV/#it{c}^{2})");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.6);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetRangeUser(0, 1.2 * h->GetMaximum());
    h->SetTitle("");
    h->SetStats(0);
}