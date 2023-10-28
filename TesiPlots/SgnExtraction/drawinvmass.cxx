void histobeauty(TH1D* h);

void drawinvmass(TString part = "AntiLambda"){
    TFile *f;
    int nbin = 3;
    int nsigma = 6;
    if (part.Contains("AntiLambda"))
    {
        f = TFile::Open("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/resultsFD/Results-AntiLambda-13TeV-SPDClusters_000_100_V0M_000_100_FDUseMCRatio.root");
        nsigma = 6;
        nbin = 3;
    }
    else if (part.Contains("K0Short")) {
        f = TFile::Open("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/results/Results-K0Short-13TeV-SPDClusters_000_100_V0M_000_100_FDNoFD.root");
        nbin = 3;
        nsigma = 6;
    }
    else if (part.Contains("Lambda")) {
        f = TFile::Open("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/resultsFD/Results-Lambda-13TeV-SPDClusters_000_100_V0M_000_100_FDUseMCRatio.root");
        nsigma = 6;
        nbin = 3;
    }
    else if (part.Contains("XiMinus")) {
        f = TFile::Open("/home/fercoles/strg_analysis/22Sett2023/CascadeAnalysis/CascadeAnalysis/results/Results-XiMinus-13TeV-SPDClusters_000_100_V0M_000_100.root");
        nbin = 5;
        nsigma = 4;
    }
    else if (part.Contains("XiPlus")) {
        f = TFile::Open("/home/fercoles/strg_analysis/22Sett2023/CascadeAnalysis/results/Results-XiPlus-13TeV-SPDClusters_000_100_V0M_000_100.root");
        nbin = 5;
        nsigma = 4;
    }
    else return;

    TH1D *h;
    if (part.Contains("K0") || part.Contains("Lambda")) h = (TH1D*)f->Get(Form("lInvMassReal/lInvMassRealRawData/lHistoSelectedV0%i", nbin));
    else if (part.Contains("Xi")) h = (TH1D*)f->Get(Form("lInvMassReal/lInvMassRealRawData/lHistoSelectedCasc%i", nbin));
    else return;
    TF1 *f1 = (TF1 *)f->Get(Form("lInvMassReal/lInvMassRealRawData/fGausPt%i", nbin));
    TCanvas* c = new TCanvas("c", "c", 900, 800);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.15);
    c->SetTopMargin(0.05);
    c->SetTicks();
    histobeauty(h);
    if (part.Contains("K0Short")) h->GetXaxis()->SetRangeUser(0.43,0.56);
    if (part.Contains("Lambda")) h->GetXaxis()->SetRangeUser(1.09,1.145);
    if (part.Contains("AntiLambda")) h->GetXaxis()->SetRangeUser(1.09,1.145);
    if (part.Contains("Xi")) h->GetXaxis()->SetRangeUser(1.295,1.35);
    h->Fit(f1, "NWWR");
    f1->SetLineStyle(2);
    f1->SetLineWidth(2);
    f1->SetLineColor(kRed);
    TH1D *hclone = (TH1D *)h->Clone();
    hclone->Reset();
    TH1D *hclonebkg1 = (TH1D *)h->Clone();
    TH1D *hclonesgn1 = (TH1D *)h->Clone();
    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
        if (h->GetBinLowEdge(i+1) > (f1->GetParameter(1) - 2 * nsigma * f1->GetParameter(2)) && h->GetBinLowEdge(i+1) < (f1->GetParameter(1) - nsigma * f1->GetParameter(2)))
            hclonebkg1->SetBinContent(i, h->GetBinContent(i));
        else if (h->GetBinLowEdge(i+1) > (f1->GetParameter(1) + nsigma * f1->GetParameter(2)) && h->GetBinLowEdge(i+1) < (f1->GetParameter(1) + 2 * nsigma * f1->GetParameter(2)))
            hclonebkg1->SetBinContent(i, h->GetBinContent(i));
        else
            hclonebkg1->SetBinContent(i, 0);

        if (h->GetBinLowEdge(i+1) > (f1->GetParameter(1) - nsigma * f1->GetParameter(2)) && h->GetBinLowEdge(i+1) < (f1->GetParameter(1) + nsigma * f1->GetParameter(2)))
            hclonesgn1->SetBinContent(i, h->GetBinContent(i));
        else
            hclonesgn1->SetBinContent(i, 0);
    }
    hclonebkg1->SetFillColor(kGray);
    hclonesgn1->SetFillColor(kRed - 10);
    hclone->Draw();
    hclonebkg1->Draw("SAME hist");
    hclonesgn1->Draw("SAME hist");
    h->Draw("EPsame");
    f1->Draw("same");

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
        ltx->DrawLatexNDC(0.2, 0.76, "0.3 < #it{p}_{T} < 0.4 GeV/#it{c}");
        ltxbig->DrawLatexNDC(0.65, 0.88, "K^{0}_{S} #rightarrow #pi^{+} + #pi^{-}");
    }
    if (part.Contains("AntiLambda"))
    {
        ltx->DrawLatexNDC(0.2, 0.76, "1.0 < #it{p}_{T} < 1.2 GeV/#it{c}");
        ltxbig->DrawLatexNDC(0.65, 0.88, "#bar{#Lambda} #rightarrow #bar{p} + #pi^{+}");
    }
    else if (part.Contains("Lambda"))
    {
        ltx->DrawLatexNDC(0.2, 0.76, "1.0 < #it{p}_{T} < 1.2 GeV/#it{c}");
        ltxbig->DrawLatexNDC(0.65, 0.88, "#Lambda #rightarrow p + #pi^{-}");
    }
    if (part.Contains("XiMinus"))
    {
        ltx->DrawLatexNDC(0.2, 0.76, "1.8 < #it{p}_{T} < 2.0 GeV/#it{c}");
        ltxbig->DrawLatexNDC(0.65, 0.88, "#Xi^{-} #rightarrow #Lambda + #pi^{-}");
    }
    if (part.Contains("XiPlus"))
    {
        ltx->DrawLatexNDC(0.2, 0.76, "1.8 < #it{p}_{T} < 2.0 GeV/#it{c}");
        ltxbig->DrawLatexNDC(0.65, 0.88, "#bar{#Xi}^{+} #rightarrow #bar{#Lambda} + #pi^{+}");
    }


    TLine *line = new TLine(f1->GetParameter(1) - nsigma*f1->GetParameter(2), 0, f1->GetParameter(1) - nsigma*f1->GetParameter(2), h->GetMaximum()*2/3);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw("same");
    TLine *line2 = new TLine(f1->GetParameter(1) + nsigma * f1->GetParameter(2), 0, f1->GetParameter(1) + nsigma * f1->GetParameter(2), h->GetMaximum()*2/3);
    line2->SetLineColor(kBlack);
    line2->SetLineStyle(2);
    line2->Draw("same");
    TLine *line3 = new TLine(f1->GetParameter(1) - 2*nsigma* f1->GetParameter(2), 0, f1->GetParameter(1) - 2*nsigma* f1->GetParameter(2), h->GetMaximum() * 2 / 3);
    line3->SetLineColor(kBlack);
    line3->SetLineStyle(2);
    line3->Draw("same");
    TLine *line4 = new TLine(f1->GetParameter(1) + 2*nsigma* f1->GetParameter(2), 0, f1->GetParameter(1) + 2*nsigma* f1->GetParameter(2), h->GetMaximum() * 2 / 3);
    line4->SetLineColor(kBlack);
    line4->SetLineStyle(2);
    line4->Draw("same");

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

    TLegend *leg = new TLegend(0.68, 0.6, 0.88, 0.82);
    leg->AddEntry(h, "Data", "lep");
    leg->AddEntry(hclonebkg1, "Background", "f");
    leg->AddEntry(hclonesgn1, "Signal", "f");
    leg->AddEntry(f1, "Gaus+pol0 fit", "l");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.032);
    leg->Draw("SAME");
    gPad->RedrawAxis();

    c->SaveAs(Form("InvMass_%s.pdf",part.Data()));
}

void histobeauty(TH1D* h){
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    h->GetXaxis()->SetTitle("m_{inv} (GeV/#it{c}^{2})");
    h->GetYaxis()->SetTitle(Form("Counts / %.0f MeV/#it{c}^{2}", h->GetBinWidth(1) * 1000));
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