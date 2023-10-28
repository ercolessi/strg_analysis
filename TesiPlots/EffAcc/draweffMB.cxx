void histobeauty(TH1D *h);
void histobeautyratio(TH1D *h);

void draweffMB(TString part = "XiMinus"){

    TLatex *ltx = new TLatex();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.035);
    ltx->SetTextAlign(12);
    TLatex *ltxbig = new TLatex();
    ltxbig->SetTextFont(42);
    ltxbig->SetTextSize(0.055);
    ltxbig->SetTextAlign(12);
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

    TFile *f[5];
    f[2] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/resultsFD/Results-AntiLambda-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_FDUseMCRatio.root", 0., 100., 0., 100.));
    f[1] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/resultsFD/Results-Lambda-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_FDUseMCRatio.root", 0., 100., 0., 100.));
    f[0] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/results/Results-K0Short-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_FDUseMCRatio.root", 0., 100., 0., 100.));
    f[3] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/CascadeAnalysis/results/Results-XiMinus-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root", 0., 100., 0., 100.));
    f[4] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/CascadeAnalysis/results/Results-XiPlus-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root", 0., 100., 0., 100.));

    TH1D *h[5];
    for (int i = 0; i < 5; i++){
        if (i==0) h[i] = (TH1D *)f[i]->Get("fHistPureEfficiency");
        else h[i] = (TH1D *)f[i]->Get("fHistEfficiency");
        histobeauty(h[i]);
        h[i]->SetLineColor(kBlack);
        h[i]->SetMarkerColor(kBlack);
    }

    TCanvas *c = new TCanvas("c", "c", 900, 800);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.15);
    c->SetTopMargin(0.05);
    c->SetTicks();

    h[0]->Draw();
    h[0]->GetXaxis()->SetRangeUser(0., 10.);

    TLegend* leg = new TLegend(0.75, 0.8, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.06);
    leg->AddEntry(h[0], "K^{0}_{S}", "lep");
    leg->Draw();
    tex->Draw();
    tex2->Draw();

    c->SaveAs(Form("EffMB_%s.pdf", "K0Short"));

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 800);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.05);
    c1->SetBottomMargin(0.15);
    c1->SetTopMargin(0.05);
    c1->SetTicks();

    h[1]->Draw();
    h[2]->SetMarkerStyle(kOpenCircle);
    h[2]->Draw("SAME");
    h[1]->GetXaxis()->SetRangeUser(0., 8.);

    TLegend *leg2 = new TLegend(0.75, 0.75, 0.89, 0.89);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.06);
    leg2->AddEntry(h[1], "#Lambda", "lep");
    leg2->AddEntry(h[2], "#bar{#Lambda}", "lep");
    leg2->Draw();
    tex->Draw();
    tex2->Draw();

    c1->SaveAs(Form("EffMB_%s.pdf", "LambdaAntiLambda"));

    TCanvas *c2 = new TCanvas("c2", "c2", 900, 800);
    c2->SetLeftMargin(0.15);
    c2->SetRightMargin(0.05);
    c2->SetBottomMargin(0.15);
    c2->SetTopMargin(0.05);
    c2->SetTicks();

    h[3]->Draw();
    h[4]->SetMarkerStyle(kOpenCircle);
    h[4]->Draw("SAME");
    h[3]->GetXaxis()->SetRangeUser(0., 6.5);

    TLegend *leg3 = new TLegend(0.75, 0.75, 0.89, 0.89);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.06);
    leg3->AddEntry(h[1], "#Xi^{-}", "lep");
    leg3->AddEntry(h[2], "#bar{#Xi}^{+}", "lep");
    leg3->Draw();
    tex->Draw();
    tex2->Draw();

    c2->SaveAs(Form("EffMB_%s.pdf", "XiMinusXiPlus"));
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