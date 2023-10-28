void preppad(TPad *pad1, TPad *pad2);
void prepcanvas(TCanvas *c);
void prepgraph(TGraphErrors* g = 0x0, int marker = 20, float size = 2., int color = kBlack, int linewidth = 1, int linestyle = 1, int fillstyle = 0);
void prephisto(TH1D* h = 0x0, bool donch = 1);
void prephisto2(TH1D *h = 0x0, bool donch = 1);

void PlotFinalStandaloneWModels(TString pythia = "Ropes")
{
    //1 con 5, 2 con 4
    TString sel[] = {"_SPDClustersV0M_class0", "_SPDClustersV0M_class1", "_SPDClustersV0M_class5", "_SPDClustersV0M_class2", "_SPDClustersV0M_class4"};
    const int nsel = sizeof(sel) / sizeof(TString);
    Int_t imarker[] = {kFullDiamond, kFullCircle, kFullSquare, kFullCircle, kFullSquare};
    Double_t imarkersize[] = {2.8, 2., 1.8, 2., 1.8};
    Int_t icolor[] = {kBlack, kRed, kGreen+1, kBlue, kViolet};


    TFile *fXi[nsel], *fLambda[nsel], *fK0s[nsel];
    for (int i = 0; i < nsel; i++){
        fXi[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYields%s_June23.root", "Xi", sel[i].Data()));
        fLambda[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYields%s_June23.root", "Lambda", sel[i].Data()));
        fK0s[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYields%s_June23.root", "K0Short", sel[i].Data()));
    }

    TGraphErrors *gNumNchStatXi[nsel], *gNumZDCStatXi[nsel];
    TGraphErrors *gNumNchSystXi[nsel], *gNumZDCSystXi[nsel];
    TGraphErrors *gNumNchStatLambda[nsel], *gNumZDCStatLambda[nsel];
    TGraphErrors *gNumNchSystLambda[nsel], *gNumZDCSystLambda[nsel];
    TGraphErrors *gNumNchStatK0s[nsel], *gNumZDCStatK0s[nsel];
    TGraphErrors *gNumNchSystK0s[nsel], *gNumZDCSystK0s[nsel];

    for (int i = 0; i < nsel; i++){
        gNumNchStatXi[i] = (TGraphErrors *)fXi[i]->Get("NormYieldsNormNchStat");
        gNumZDCStatXi[i] = (TGraphErrors *)fXi[i]->Get("NormYieldsNormZDCSumStat");
        gNumNchSystXi[i] = (TGraphErrors *)fXi[i]->Get("NormYieldsNormNchSyst");
        gNumZDCSystXi[i] = (TGraphErrors *)fXi[i]->Get("NormYieldsNormZDCSumSyst");
        gNumNchStatLambda[i] = (TGraphErrors *)fLambda[i]->Get("NormYieldsNormNchStat");
        gNumZDCStatLambda[i] = (TGraphErrors *)fLambda[i]->Get("NormYieldsNormZDCSumStat");
        gNumNchSystLambda[i] = (TGraphErrors *)fLambda[i]->Get("NormYieldsNormNchSyst");
        gNumZDCSystLambda[i] = (TGraphErrors *)fLambda[i]->Get("NormYieldsNormZDCSumSyst");
        gNumNchStatK0s[i] = (TGraphErrors *)fK0s[i]->Get("NormYieldsNormNchStat");
        gNumZDCStatK0s[i] = (TGraphErrors *)fK0s[i]->Get("NormYieldsNormZDCSumStat");
        gNumNchSystK0s[i] = (TGraphErrors *)fK0s[i]->Get("NormYieldsNormNchSyst");
        gNumZDCSystK0s[i] = (TGraphErrors *)fK0s[i]->Get("NormYieldsNormZDCSumSyst");
    }

    std::vector<Double_t> Nch[nsel], NchErr[nsel], NchErrSyst[nsel], ZDC[nsel], ZDCErr[nsel], ZDCErrSyst[nsel];
    std::vector<Double_t> RatioXi[nsel], RatioXiStatErr[nsel], RatioXiSystErr[nsel], RatioLambda[nsel], RatioLambdaStatErr[nsel], RatioLambdaSystErr[nsel], RatioK0s[nsel], RatioK0sStatErr[nsel], RatioK0sSystErr[nsel];
    TGraphErrors *gNchRatioStatXi[nsel], *gZDCRatioStatXi[nsel], *gNchRatioSystXi[nsel], *gZDCRatioSystXi[nsel];
    TGraphErrors *gNchRatioStatLambda[nsel], *gZDCRatioStatLambda[nsel], *gNchRatioSystLambda[nsel], *gZDCRatioSystLambda[nsel];
    TGraphErrors *gNchRatioStatK0s[nsel], *gZDCRatioStatK0s[nsel], *gNchRatioSystK0s[nsel], *gZDCRatioSystK0s[nsel];

    for (int i = 0; i < nsel; i++){
        for (int j = 0; j < gNumNchStatXi[i]->GetN(); j++){
            Nch[i].push_back(gNumNchStatXi[i]->GetX()[j]);
            NchErr[i].push_back(gNumNchStatXi[i]->GetEX()[j]);
            NchErrSyst[i].push_back(gNumNchSystXi[i]->GetEXhigh()[j]);
            ZDC[i].push_back(gNumZDCStatXi[i]->GetX()[j]);
            ZDCErr[i].push_back(gNumZDCStatXi[i]->GetEX()[j]);
            ZDCErrSyst[i].push_back(gNumZDCSystXi[i]->GetEXhigh()[j]);
            RatioXi[i].push_back(gNumNchStatXi[i]->GetY()[j]);
            RatioXiStatErr[i].push_back(gNumNchStatXi[i]->GetEY()[j]);
            RatioXiSystErr[i].push_back(gNumNchSystXi[i]->GetEYhigh()[j]);
            RatioLambda[i].push_back(gNumNchStatLambda[i]->GetY()[j]);
            RatioLambdaStatErr[i].push_back(gNumNchStatLambda[i]->GetEY()[j]);
            RatioLambdaSystErr[i].push_back(gNumNchSystLambda[i]->GetEYhigh()[j]);
            RatioK0s[i].push_back(gNumNchStatK0s[i]->GetY()[j]);
            RatioK0sStatErr[i].push_back(gNumNchStatK0s[i]->GetEY()[j]);
            RatioK0sSystErr[i].push_back(gNumNchSystK0s[i]->GetEYhigh()[j]);
        }

        gNchRatioStatXi[i] = new TGraphErrors(gNumNchStatXi[i]->GetN(), &Nch[i][0], &RatioXi[i][0], &NchErr[i][0], &RatioXiStatErr[i][0]);
        gZDCRatioStatXi[i] = new TGraphErrors(gNumNchStatXi[i]->GetN(), &ZDC[i][0], &RatioXi[i][0], &ZDCErr[i][0], &RatioXiStatErr[i][0]);
        gNchRatioSystXi[i] = new TGraphErrors(gNumNchSystXi[i]->GetN(), &Nch[i][0], &RatioXi[i][0], &NchErrSyst[i][0], &RatioXiSystErr[i][0]);
        gZDCRatioSystXi[i] = new TGraphErrors(gNumNchSystXi[i]->GetN(), &ZDC[i][0], &RatioXi[i][0], &ZDCErrSyst[i][0], &RatioXiSystErr[i][0]);
        gNchRatioStatLambda[i] = new TGraphErrors(gNumNchStatLambda[i]->GetN(), &Nch[i][0], &RatioLambda[i][0], &NchErr[i][0], &RatioLambdaStatErr[i][0]);
        gZDCRatioStatLambda[i] = new TGraphErrors(gNumNchStatLambda[i]->GetN(), &ZDC[i][0], &RatioLambda[i][0], &ZDCErr[i][0], &RatioLambdaStatErr[i][0]);
        gNchRatioSystLambda[i] = new TGraphErrors(gNumNchSystLambda[i]->GetN(), &Nch[i][0], &RatioLambda[i][0], &NchErrSyst[i][0], &RatioLambdaSystErr[i][0]);
        gZDCRatioSystLambda[i] = new TGraphErrors(gNumNchSystLambda[i]->GetN(), &ZDC[i][0], &RatioLambda[i][0], &ZDCErrSyst[i][0], &RatioLambdaSystErr[i][0]);
        gNchRatioStatK0s[i] = new TGraphErrors(gNumNchStatK0s[i]->GetN(), &Nch[i][0], &RatioK0s[i][0], &NchErr[i][0], &RatioK0sStatErr[i][0]);
        gZDCRatioStatK0s[i] = new TGraphErrors(gNumNchStatK0s[i]->GetN(), &ZDC[i][0], &RatioK0s[i][0], &ZDCErr[i][0], &RatioK0sStatErr[i][0]);
        gNchRatioSystK0s[i] = new TGraphErrors(gNumNchSystK0s[i]->GetN(), &Nch[i][0], &RatioK0s[i][0], &NchErrSyst[i][0], &RatioK0sSystErr[i][0]);
        gZDCRatioSystK0s[i] = new TGraphErrors(gNumNchSystK0s[i]->GetN(), &ZDC[i][0], &RatioK0s[i][0], &ZDCErrSyst[i][0], &RatioK0sSystErr[i][0]);

        prepgraph(gNchRatioStatXi[i], imarker[i], imarkersize[i], kGreen+1);
        prepgraph(gZDCRatioStatXi[i], imarker[i], imarkersize[i], kGreen+1);
        prepgraph(gNchRatioSystXi[i], imarker[i], 0, kGreen+1);
        prepgraph(gZDCRatioSystXi[i], imarker[i], 0, kGreen+1);
        prepgraph(gNchRatioStatLambda[i], kFullCircle, 2., kBlue);
        prepgraph(gZDCRatioStatLambda[i], kFullCircle, 2., kBlue);
        prepgraph(gNchRatioSystLambda[i], kFullCircle, 0, kBlue);
        prepgraph(gZDCRatioSystLambda[i], kFullCircle, 0, kBlue);
        prepgraph(gNchRatioStatK0s[i], kFullSquare,1.8, kRed);
        prepgraph(gZDCRatioStatK0s[i], kFullSquare,1.8, kRed);
        prepgraph(gNchRatioSystK0s[i], kFullSquare, 0, kRed);
        prepgraph(gZDCRatioSystK0s[i], kFullSquare, 0, kRed);
    }

    TLatex* l = new TLatex();
    l->SetTextFont(42);
    l->SetNDC();
    l->SetTextColor(1);
    l->SetTextSize(0.04);
    l->SetTextAlign(22);
    l->SetTextAngle(0);

    TLatex *lbig = new TLatex();
    lbig->SetTextFont(42);
    lbig->SetNDC();
    lbig->SetTextColor(1);
    lbig->SetTextSize(0.07);
    lbig->SetTextAlign(22);
    lbig->SetTextAngle(0);

    TLatex *l2 = new TLatex();
    l2->SetTextFont(42);
    l2->SetNDC();
    l2->SetTextColor(1);
    l2->SetTextSize(0.06);
    l2->SetTextAlign(22);
    l2->SetTextAngle(0);

    TLatex *lbig2 = new TLatex();
    lbig2->SetTextFont(42);
    lbig2->SetNDC();
    lbig2->SetTextColor(1);
    lbig2->SetTextSize(0.095);
    lbig2->SetTextAlign(22);
    lbig2->SetTextAngle(0);


    TH1D *hnch = new TH1D("hnch", "", 10, 0.1, 4.3);
    prephisto(hnch, 1);
    TH1D *hzdc = new TH1D("hzdc", "", 10, 0.1, 1.6);
    prephisto2(hzdc, 0);


    //models
    TString ddname[] = {"Standalone", "kHighMult", "kLowMult", "kHighZN", "kLowZN"};

    TFile *fileXim = TFile::Open(Form("../models/ParticleRatios/Ratio%s%s_%s.root", "XiMinusXiPlus", "nch", "Monash"));
    TFile *fileXir = TFile::Open(Form("../models/ParticleRatios/Ratio%s%s_%s.root", "XiMinusXiPlus", "nch", "Ropes"));
    TFile *fileLambdam = TFile::Open(Form("../models/ParticleRatios/Ratio%s%s_%s.root", "LambdaAntiLambda", "nch", "Monash"));
    TFile *fileLambdar = TFile::Open(Form("../models/ParticleRatios/Ratio%s%s_%s.root", "LambdaAntiLambda", "nch", "Ropes"));
    TFile *fileK0sm = TFile::Open(Form("../models/ParticleRatios/Ratio%s%s_%s.root", "K0Short", "nch", "Monash"));
    TFile *fileK0sr = TFile::Open(Form("../models/ParticleRatios/Ratio%s%s_%s.root", "K0Short", "nch", "Ropes"));

    TGraphErrors *gXiNchm, *gXiNchr, *gXiLEm, *gXiLEr;
    TGraphErrors *gLambdaNchm, *gLambdaNchr, *gLambdaLEm, *gLambdaLEr;
    TGraphErrors *gK0sNchm, *gK0sNchr, *gK0sLEm, *gK0sLEr;

    TGraphErrors *hblack = (TGraphErrors *)fileXim->Get(Form("gNchRatio_%s", "Standalone"));
    gXiNchm = (TGraphErrors *)fileXim->Get(Form("gNchRatio_%s", "Standalone"));
    gXiLEm = (TGraphErrors *)fileXim->Get(Form("gZDCRatio_%s", "Standalone"));
    gLambdaNchm = (TGraphErrors *)fileLambdam->Get(Form("gNchRatio_%s", "Standalone"));
    gLambdaLEm = (TGraphErrors *)fileLambdam->Get(Form("gZDCRatio_%s", "Standalone"));
    gK0sNchm = (TGraphErrors *)fileK0sm->Get(Form("gNchRatio_%s", "Standalone"));
    gK0sLEm = (TGraphErrors *)fileK0sm->Get(Form("gZDCRatio_%s", "Standalone"));
    gXiNchr = (TGraphErrors *)fileXir->Get(Form("gNchRatio_%s", "Standalone"));
    gXiLEr = (TGraphErrors *)fileXir->Get(Form("gZDCRatio_%s", "Standalone"));
    gLambdaNchr = (TGraphErrors *)fileLambdar->Get(Form("gNchRatio_%s", "Standalone"));
    gLambdaLEr = (TGraphErrors *)fileLambdar->Get(Form("gZDCRatio_%s", "Standalone"));
    gK0sNchr = (TGraphErrors *)fileK0sr->Get(Form("gNchRatio_%s", "Standalone"));
    gK0sLEr = (TGraphErrors *)fileK0sr->Get(Form("gZDCRatio_%s", "Standalone"));

    TCanvas *c = new TCanvas("c", "", 1400, 800);
    prepcanvas(c);
    TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 0.54, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0.54, 0, 1., 1);
    preppad(pad1, pad2);

    pad1->Draw();
    pad2->Draw();
    //
    pad1->cd();
    hnch->Draw();
    gXiNchm->SetLineWidth(2);
    gXiNchm->SetLineColor(kGreen + 2);
    gLambdaNchm->SetLineWidth(2);
    gLambdaNchm->SetLineColor(kBlue);
    gK0sNchm->SetLineWidth(2);
    gK0sNchm->SetLineColor(kRed);

    if (pythia.Contains("Monash")){
        gXiNchm->Draw("SAME");
        gLambdaNchm->Draw("SAME");
        gK0sNchm->Draw("SAME");
    } else if (pythia.Contains("Ropes")){
        gXiNchr->SetLineWidth(2);
        gXiNchr->SetLineColor(kGreen + 2);
        gLambdaNchr->SetLineWidth(2);
        gLambdaNchr->SetLineColor(kBlue);
        gK0sNchr->SetLineWidth(2);
        gK0sNchr->SetLineColor(kRed);

        gXiNchr->Draw("SAME");
        gLambdaNchr->Draw("SAME");
        gK0sNchr->Draw("SAME");
    }
    gNchRatioSystK0s[0]->Draw("SAME E2");
    gNchRatioStatK0s[0]->Draw("SAME EP");
    gNchRatioSystLambda[0]->Draw("SAME E2");
    gNchRatioStatLambda[0]->Draw("SAME EP");
    gNchRatioSystXi[0]->Draw("SAME E2");
    gNchRatioStatXi[0]->Draw("SAME EP");

    l->DrawLatex(0.43, 0.89, Form("%s", "ALICE, pp #sqrt{s} = 13 TeV"));

    TLegend *legm = new TLegend(0.5, 0.3, 0.85, 0.5);
    legm->SetBorderSize(0);
    legm->SetFillStyle(0);
    legm->SetTextSize(0.04);
    legm->AddEntry(hblack, Form("Pythia %s", pythia.Data()), "L");
    legm->Draw("same");
    //
    pad2->cd();
    hzdc->Draw();
    if (pythia.Contains("Monash")){
        gXiLEm->SetLineWidth(2);
        gXiLEm->SetLineColor(kGreen + 2);
        gLambdaLEm->SetLineWidth(2);
        gLambdaLEm->SetLineColor(kBlue);
        gK0sLEm->SetLineWidth(2);
        gK0sLEm->SetLineColor(kRed);

        gXiLEm->Draw("SAME");
        gLambdaLEm->Draw("SAME");
        gK0sLEm->Draw("SAME");
    } else if (pythia.Contains("Ropes")){
        gXiLEr->SetLineWidth(2);
        gXiLEr->SetLineColor(kGreen + 2);
        gLambdaLEr->SetLineWidth(2);
        gLambdaLEr->SetLineColor(kBlue);
        gK0sLEr->SetLineWidth(2);
        gK0sLEr->SetLineColor(kRed);

        gXiLEr->Draw("SAME");
        gLambdaLEr->Draw("SAME");
        gK0sLEr->Draw("SAME");
    }
    gZDCRatioSystK0s[0]->Draw("SAME E2");
    gZDCRatioStatK0s[0]->Draw("SAME EP");
    gZDCRatioSystLambda[0]->Draw("SAME E2");
    gZDCRatioStatLambda[0]->Draw("SAME EP");
    gZDCRatioSystXi[0]->Draw("SAME E2");
    gZDCRatioStatXi[0]->Draw("SAME EP");

    l->DrawLatex(0.5, 0.8, Form("%s", "|#it{y}|<0.5"));

    TLegend *leg = new TLegend(0.65, 0.7, 0.85, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.055);
    leg->AddEntry(gZDCRatioStatK0s[0], "K^{0}_{S}", "P");
    leg->AddEntry(gZDCRatioStatLambda[0], "#Lambda", "P");
    leg->AddEntry(gZDCRatioStatXi[0], "#Xi", "P");
    leg->Draw("same");

    hblack->SetLineColor(kBlack);

    c->SaveAs(Form("images/FinalPlotStandalone_wmodels_%s.pdf",pythia.Data()));

 }

void prepcanvas(TCanvas* c){
    c->SetLeftMargin(0.2);
    c->SetBottomMargin(0.25);
    c->SetRightMargin(0.1);
    c->SetTopMargin(0.1);
    c->SetTicky();
    c->SetTickx();
}

void preppad(TPad *pad1, TPad *pad2)
{
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
void prepgraph(TGraphErrors* g = 0x0, int marker = 20, float size = 2., int color = kBlack, int linewidth = 1, int linestyle = 1, int fillstyle = 0){
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(size);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetFillColor(color);
    g->SetLineWidth(linewidth);
    g->SetLineStyle(linestyle);
    g->SetFillStyle(fillstyle);
}

void prephisto(TH1D* h = 0x0, bool donch = 1){
    h->SetStats(0);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.035);
    h->GetXaxis()->SetLabelSize(0.035);
    h->GetYaxis()->SetTitleOffset(1.5);
    h->GetYaxis()->SetRangeUser(0.45, 1.4);
    h->SetTitle("");
    if (donch) h->GetXaxis()->SetTitle("n_{ch} / #LT n_{ch} #GT_{MB} ");
    else h->GetXaxis()->SetTitle("ZN / #LT ZN #GT_{MB}");
    h->GetYaxis()->SetTitle("#frac{#it{h}/#LT #it{h} #GT_{INEL>0}}{#it{n}_{ch}/#LT #it{n}_{ch} #GT_{INEL>0}}");
}

void prephisto2(TH1D* h = 0x0, bool donch = 1){
    h->SetStats(0);
    h->GetXaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTitleOffset(0.96);
    h->GetYaxis()->SetRangeUser(0.45, 1.4);
    h->SetTitle("");
    if (donch) h->GetXaxis()->SetTitle("n_{ch} / #LT n_{ch} #GT_{MB} ");
    else h->GetXaxis()->SetTitle("ZN / #LT ZN #GT_{MB}");
    h->GetYaxis()->SetTitle("#frac{#it{h}/#LT #it{h} #GT_{INEL>0}}{#it{n}_{ch}/#LT #it{n}_{ch} #GT_{INEL>0}}");
}