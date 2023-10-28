void preppadright(TPad *pad1, TPad *pad2, TPad *pad3);
void preppadleft(TPad *pad1, TPad *pad2, TPad *pad3);
void prepcanvas(TCanvas* c);
void prepgraph(TGraphErrors* g = 0x0, int marker = 20, float size = 2., int color = kBlack, int linewidth = 1, int linestyle = 1, int fillstyle = 0);
void prephisto(TH1D* h = 0x0, bool donch = 1);
void prephisto2(TH1D *h = 0x0, bool donch = 1);

void PlotFinalVertical(bool donch = 0, bool domodels = 1, TString pythia = "Monash")
{
    //1 con 5, 2 con 4
    TString sel[] = {"_SPDClustersV0M_class0", "_SPDClustersV0M_class1", "_SPDClustersV0M_class5", "_SPDClustersV0M_class2", "_SPDClustersV0M_class4"};
    const int nsel = sizeof(sel) / sizeof(TString);
    Int_t imarker[] = {kFullDiamond, kFullCircle, kFullSquare, kFullCircle, kFullSquare};
    Double_t imarkersize[] = {2.2, 1.7, 1.5, 1.7, 1.5};
    Int_t icolor[] = {kBlack, kRed, kGreen+1, kBlue, kViolet};
    Int_t icolormodel[] = {kGray + 3, kRed + 1, kGreen - 2, kBlue - 3, kMagenta + 1};

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

        prepgraph(gNchRatioStatXi[i], imarker[i], imarkersize[i], icolor[i]);
        prepgraph(gZDCRatioStatXi[i], imarker[i], imarkersize[i], icolor[i]);
        prepgraph(gNchRatioSystXi[i], imarker[i], 0, icolor[i]);
        prepgraph(gZDCRatioSystXi[i], imarker[i], 0, icolor[i]);
        prepgraph(gNchRatioStatLambda[i], imarker[i], imarkersize[i], icolor[i]);
        prepgraph(gZDCRatioStatLambda[i], imarker[i], imarkersize[i], icolor[i]);
        prepgraph(gNchRatioSystLambda[i], imarker[i], 0, icolor[i]);
        prepgraph(gZDCRatioSystLambda[i], imarker[i], 0, icolor[i]);
        prepgraph(gNchRatioStatK0s[i], imarker[i], imarkersize[i], icolor[i]);
        prepgraph(gZDCRatioStatK0s[i], imarker[i], imarkersize[i], icolor[i]);
        prepgraph(gNchRatioSystK0s[i], imarker[i], 0, icolor[i]);
        prepgraph(gZDCRatioSystK0s[i], imarker[i], 0, icolor[i]);
    }

    TLatex* l = new TLatex();
    l->SetTextFont(42);
    l->SetNDC();
    l->SetTextColor(1);
    l->SetTextSize(0.045);
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

    double min = 0.1;
    double max = 4.3;
    if (!donch) {
        min = 0.1;
        max = 1.6;
    }
    TH1D *h = new TH1D("h", "", 10, min, max);
    prephisto(h, donch);
    TH1D *h2 = new TH1D("h2", "", 10, min, max);
    prephisto2(h2, donch);

    // models
    TString ddname[] = {"Standalone", "kHighMult", "kLowZN", "kLowMult", "kHighZN"};

    TFile *fileModelXi = TFile::Open(Form("../models/ParticleRatios/Ratio%s%s_%s.root", "XiMinusXiPlus", "nch", pythia.Data()));
    TFile *fileModelLambda = TFile::Open(Form("../models/ParticleRatios/Ratio%s%s_%s.root", "LambdaAntiLambda", "nch", pythia.Data()));
    TFile *fileModelK0s = TFile::Open(Form("../models/ParticleRatios/Ratio%s%s_%s.root", "K0Short", "nch", pythia.Data()));

    TGraphErrors *gXiNch[nsel], *gXiLE[nsel];
    TGraphErrors *gLambdaNch[nsel], *gLambdaLE[nsel];
    TGraphErrors *gK0sNch[nsel], *gK0sLE[nsel];

    for (int i = 0; i<nsel; i++){
        gXiNch[i] = (TGraphErrors *)fileModelXi->Get(Form("gNchRatio_%s", ddname[i].Data()));
        gXiLE[i] = (TGraphErrors *)fileModelXi->Get(Form("gZDCRatio_%s", ddname[i].Data()));
        gLambdaNch[i] = (TGraphErrors *)fileModelLambda->Get(Form("gNchRatio_%s", ddname[i].Data()));
        gLambdaLE[i] = (TGraphErrors *)fileModelLambda->Get(Form("gZDCRatio_%s", ddname[i].Data()));
        gK0sNch[i] = (TGraphErrors *)fileModelK0s->Get(Form("gNchRatio_%s", ddname[i].Data()));
        gK0sLE[i] = (TGraphErrors *)fileModelK0s->Get(Form("gZDCRatio_%s", ddname[i].Data()));

        gXiNch[i]->SetLineWidth(2);
        gXiNch[i]->SetLineColor(icolormodel[i]);
        gXiLE[i]->SetLineWidth(2);
        gXiLE[i]->SetLineColor(icolormodel[i]);
        gLambdaNch[i]->SetLineWidth(2);
        gLambdaNch[i]->SetLineColor(icolormodel[i]);
        gLambdaLE[i]->SetLineWidth(2);
        gLambdaLE[i]->SetLineColor(icolormodel[i]);
        gK0sNch[i]->SetLineWidth(2);
        gK0sNch[i]->SetLineColor(icolormodel[i]);
        gK0sLE[i]->SetLineWidth(2);
        gK0sLE[i]->SetLineColor(icolormodel[i]);
    }

    TCanvas *c = new TCanvas("c", "", 900, 1000);
    prepcanvas(c);
    TPad *pad1left = new TPad("padleft1", "padleft1", 0., 0.0, 0.56, 0.4);
    TPad *pad2left = new TPad("padleft2", "padleft2", 0., 0.4, 0.56, 0.7);
    TPad *pad3left = new TPad("padleft3", "padleft3", 0., 0.7, 0.56, 1.);
    preppadleft(pad1left, pad2left, pad3left);
    TPad *pad1right = new TPad("padright1", "padright1", 0.56, 0.0, 1., 0.4);
    TPad *pad2right = new TPad("padright2", "padright2", 0.56, 0.4, 1., 0.7);
    TPad *pad3right = new TPad("padright3", "padright3", 0.56, 0.7, 1., 1.);
    preppadright(pad1right, pad2right, pad3right);

    pad1left->Draw();
    pad2left->Draw();
    pad3left->Draw();
    pad1right->Draw();
    pad2right->Draw();
    pad3right->Draw();

    TLegend *legm = new TLegend(0.6, 0.25, 0.85, 0.4);
    legm->SetBorderSize(0);
    legm->SetFillStyle(0);
    legm->SetTextSize(0.05);
    legm->AddEntry(gK0sNch[0], Form("Pythia %s",pythia.Data()), "L");
    TLegend *legm2 = new TLegend(0.5, 0.25, 0.75, 0.4);
    legm2->SetBorderSize(0);
    legm2->SetFillStyle(0);
    legm2->SetTextSize(0.05);
    legm2->AddEntry(gK0sNch[0], Form("Pythia %s", pythia.Data()), "L");

    pad1left->cd();
    h->Draw();
    if (donch)
    {
        if (domodels){
            gK0sNch[0]->Draw("C SAME");
            gK0sNch[1]->Draw("C SAME");
            gK0sNch[2]->Draw("C SAME");
        }

        gNchRatioSystK0s[0]->Draw("SAME E2");
        gNchRatioStatK0s[0]->Draw("SAME EP");
        gNchRatioSystK0s[1]->Draw("SAME E2");
        gNchRatioStatK0s[1]->Draw("SAME EP");
        gNchRatioSystK0s[2]->Draw("SAME E2");
        gNchRatioStatK0s[2]->Draw("SAME EP");
    }
    else
    {
        if (domodels){
            gK0sLE[0]->Draw("C SAME");
            gK0sLE[1]->Draw("C SAME");
            gK0sLE[2]->Draw("C SAME");
        }

        gZDCRatioSystK0s[0]->Draw("SAME E2");
        gZDCRatioStatK0s[0]->Draw("SAME EP");
        gZDCRatioSystK0s[1]->Draw("SAME E2");
        gZDCRatioStatK0s[1]->Draw("SAME EP");
        gZDCRatioSystK0s[2]->Draw("SAME E2");
        gZDCRatioStatK0s[2]->Draw("SAME EP");
    }
    l->DrawLatex(0.89, 0.82, Form("%s", "|#it{y}|<0.5"));
    lbig->DrawLatex(0.89, 0.91, Form("%s", "K^{0}_{S}"));
    legm->Draw("same");

    TLegend *leg = new TLegend(0.3, 0.75, 0.5, 0.95);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.05);
    leg->AddEntry(gZDCRatioStatK0s[0], "V0M Standalone ", "p");
    leg->AddEntry(gZDCRatioStatK0s[2], "Low Multiplicity ", "p");
    leg->AddEntry(gZDCRatioStatK0s[1], "High Multiplicity ", "p");
    leg->Draw("same");

    pad2left->cd();
    h2->Draw();
    if (donch)
    {
        if (domodels)
        {
            gLambdaNch[0]->Draw("C SAME");
            gLambdaNch[1]->Draw("C SAME");
            gLambdaNch[2]->Draw("C SAME");
        }
        gNchRatioSystLambda[0]->Draw("SAME E2");
        gNchRatioStatLambda[0]->Draw("SAME EP");
        gNchRatioSystLambda[1]->Draw("SAME E2");
        gNchRatioStatLambda[1]->Draw("SAME EP");
        gNchRatioSystLambda[2]->Draw("SAME E2");
        gNchRatioStatLambda[2]->Draw("SAME EP");
    }
    else
    {
        if (domodels)
        {
            gLambdaLE[0]->Draw("C SAME");
            gLambdaLE[1]->Draw("C SAME");
            gLambdaLE[2]->Draw("C SAME");
        }
        gZDCRatioSystLambda[0]->Draw("SAME E2");
        gZDCRatioStatLambda[0]->Draw("SAME EP");
        gZDCRatioSystLambda[1]->Draw("SAME E2");
        gZDCRatioStatLambda[1]->Draw("SAME EP");
        gZDCRatioSystLambda[2]->Draw("SAME E2");
        gZDCRatioStatLambda[2]->Draw("SAME EP");
    }
    l2->DrawLatex(0.89, 0.82, Form("%s", "|#it{y}|<0.5"));
    lbig2->DrawLatex(0.89, 0.91, Form("%s", "#Lambda"));

    pad3left->cd();
    h2->Draw();
    if (donch)
    {
        if (domodels)
        {
            gXiNch[0]->Draw("C SAME");
            gXiNch[1]->Draw("C SAME");
            gXiNch[2]->Draw("C SAME");
        }
        gNchRatioSystXi[0]->Draw("SAME E2");
        gNchRatioStatXi[0]->Draw("SAME EP");
        gNchRatioSystXi[1]->Draw("SAME E2");
        gNchRatioStatXi[1]->Draw("SAME EP");
        gNchRatioSystXi[2]->Draw("SAME E2");
        gNchRatioStatXi[2]->Draw("SAME EP");
    }
    else
    {
        if (domodels)
        {
            gXiLE[0]->Draw("C SAME");
            gXiLE[1]->Draw("C SAME");
            gXiLE[2]->Draw("C SAME");
        }
        gZDCRatioSystXi[0]->Draw("SAME E2");
        gZDCRatioStatXi[0]->Draw("SAME EP");
        gZDCRatioSystXi[1]->Draw("SAME E2");
        gZDCRatioStatXi[1]->Draw("SAME EP");
        gZDCRatioSystXi[2]->Draw("SAME E2");
        gZDCRatioStatXi[2]->Draw("SAME EP");
    }
    l2->DrawLatex(0.89, 0.77, Form("%s", "|#it{y}|<0.5"));
    lbig2->DrawLatex(0.89, 0.86, Form("%s", "#Xi"));

    l2->DrawLatex(0.47, 0.87, Form("%s", "ALICE, pp #sqrt{s} = 13 TeV"));

    pad1right->cd();
    h->Draw();

    if (donch)
    {
        if (domodels)
        {
            gK0sNch[0]->Draw("C SAME");
            gK0sNch[3]->Draw("C SAME");
            gK0sNch[4]->Draw("C SAME");
        }
        gNchRatioSystK0s[0]->Draw("SAME E2");
        gNchRatioStatK0s[0]->Draw("SAME EP");
        gNchRatioSystK0s[3]->Draw("SAME E2");
        gNchRatioStatK0s[3]->Draw("SAME EP");
        gNchRatioSystK0s[4]->Draw("SAME E2");
        gNchRatioStatK0s[4]->Draw("SAME EP");
    }
    else
    {
        if (domodels)
        {
            gK0sLE[0]->Draw("C SAME");
            gK0sLE[3]->Draw("C SAME");
            gK0sLE[4]->Draw("C SAME");
        }
        gZDCRatioSystK0s[0]->Draw("SAME E2");
        gZDCRatioStatK0s[0]->Draw("SAME EP");
        gZDCRatioSystK0s[3]->Draw("SAME E2");
        gZDCRatioStatK0s[3]->Draw("SAME EP");
        gZDCRatioSystK0s[4]->Draw("SAME E2");
        gZDCRatioStatK0s[4]->Draw("SAME EP");
    }

    l->DrawLatex(0.85, 0.82, Form("%s", "|#it{y}|<0.5"));
    lbig->DrawLatex(0.85, 0.91, Form("%s", "K^{0}_{S}"));

    legm2->Draw("same");

    TLegend *leg2 = new TLegend(0.12, 0.75, 0.3, 0.95);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.05);
    leg2->AddEntry(gZDCRatioStatK0s[0], "V0M Standalone ", "p");
    leg2->AddEntry(gZDCRatioStatK0s[4], "Low ZN ", "p");
    leg2->AddEntry(gZDCRatioStatK0s[3], "High ZN ", "p");
    leg2->Draw("same");

    pad2right->cd();
    h2->Draw();

    if (donch)
    {
        if (domodels)
        {
            gLambdaNch[0]->Draw("C SAME");
            gLambdaNch[3]->Draw("C SAME");
            gLambdaNch[4]->Draw("C SAME");
        }
        gNchRatioSystLambda[0]->Draw("SAME E2");
        gNchRatioStatLambda[0]->Draw("SAME EP");
        gNchRatioSystLambda[3]->Draw("SAME E2");
        gNchRatioStatLambda[3]->Draw("SAME EP");
        gNchRatioSystLambda[4]->Draw("SAME E2");
        gNchRatioStatLambda[4]->Draw("SAME EP");
    }
    else
    {
        if (domodels)
        {
            gLambdaLE[0]->Draw("C SAME");
            gLambdaLE[3]->Draw("C SAME");
            gLambdaLE[4]->Draw("C SAME");
        }
        gZDCRatioSystLambda[0]->Draw("SAME E2");
        gZDCRatioStatLambda[0]->Draw("SAME EP");
        gZDCRatioSystLambda[3]->Draw("SAME E2");
        gZDCRatioStatLambda[3]->Draw("SAME EP");
        gZDCRatioSystLambda[4]->Draw("SAME E2");
        gZDCRatioStatLambda[4]->Draw("SAME EP");
    }

    l2->DrawLatex(0.85, 0.82, Form("%s", "|#it{y}|<0.5"));
    lbig2->DrawLatex(0.85, 0.91, Form("%s", "#Lambda"));

    pad3right->cd();
    h2->Draw();

    if (donch)
    {
        if (domodels)
        {
            gXiNch[0]->Draw("C SAME");
            gXiNch[3]->Draw("C SAME");
            gXiNch[4]->Draw("C SAME");
        }
        gNchRatioSystXi[0]->Draw("SAME E2");
        gNchRatioStatXi[0]->Draw("SAME EP");
        gNchRatioSystXi[3]->Draw("SAME E2");
        gNchRatioStatXi[3]->Draw("SAME EP");
        gNchRatioSystXi[4]->Draw("SAME E2");
        gNchRatioStatXi[4]->Draw("SAME EP");
    }
    else
    {
        if (domodels)
        {
            gXiLE[0]->Draw("C SAME");
            gXiLE[3]->Draw("C SAME");
            gXiLE[4]->Draw("C SAME");
        }
        gZDCRatioSystXi[0]->Draw("SAME E2");
        gZDCRatioStatXi[0]->Draw("SAME EP");
        gZDCRatioSystXi[3]->Draw("SAME E2");
        gZDCRatioStatXi[3]->Draw("SAME EP");
        gZDCRatioSystXi[4]->Draw("SAME E2");
        gZDCRatioStatXi[4]->Draw("SAME EP");
    }

    l2->DrawLatex(0.85, 0.77, Form("%s", "|#it{y}|<0.5"));
    lbig2->DrawLatex(0.85, 0.86, Form("%s", "#Xi"));


    c->SaveAs(Form("images/FinalPlotVertical%s_%s.pdf", donch ? "Nch" : "ZDCSum",pythia.Data()));

 }

void prepcanvas(TCanvas* c){
    c->SetLeftMargin(0.2);
    c->SetBottomMargin(0.25);
    c->SetRightMargin(0.1);
    c->SetTopMargin(0.1);
    c->SetTicky();
    c->SetTickx();
}

void preppadleft(TPad *pad1, TPad *pad2, TPad *pad3){
    pad1->SetBorderMode(0);
    pad1->SetTopMargin(0.0);
    pad1->SetRightMargin(0.0);
    pad1->SetFillColor(kWhite);
    pad1->SetTicky();
    pad1->SetTickx();
    pad1->SetLeftMargin(0.25);
    pad1->SetBottomMargin(0.2);

    pad2->SetBorderMode(0);
    pad2->SetTopMargin(0.0);
    pad2->SetRightMargin(0.0);
    pad2->SetFillColor(kWhite);
    pad2->SetTicky();
    pad2->SetTickx();
    pad2->SetLeftMargin(0.25);
    pad2->SetBottomMargin(0.0);

    pad3->SetBorderMode(0);
    pad3->SetTopMargin(0.05);
    pad3->SetRightMargin(0.0);
    pad3->SetFillColor(kWhite);
    pad3->SetTicky();
    pad3->SetTickx();
    pad3->SetLeftMargin(0.25);
    pad3->SetBottomMargin(0.0);
}

void preppadright(TPad *pad1, TPad *pad2, TPad *pad3)
{
    pad1->SetBorderMode(0);
    pad1->SetTopMargin(0.0);
    pad1->SetRightMargin(0.05);
    pad1->SetFillColor(kWhite);
    pad1->SetTicky();
    pad1->SetTickx();
    pad1->SetLeftMargin(0.0);
    pad1->SetBottomMargin(0.2);

    pad2->SetBorderMode(0);
    pad2->SetTopMargin(0.0);
    pad2->SetRightMargin(0.05);
    pad2->SetFillColor(kWhite);
    pad2->SetTicky();
    pad2->SetTickx();
    pad2->SetLeftMargin(0.0);
    pad2->SetBottomMargin(0.0);

    pad3->SetBorderMode(0);
    pad3->SetTopMargin(0.05);
    pad3->SetRightMargin(0.05);
    pad3->SetFillColor(kWhite);
    pad3->SetTicky();
    pad3->SetTickx();
    pad3->SetLeftMargin(0.0);
    pad3->SetBottomMargin(0.0);
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
    h->GetXaxis()->SetTitleSize(0.075);
    h->GetYaxis()->SetTitleSize(0.075);
    h->GetYaxis()->SetLabelSize(0.055);
    h->GetXaxis()->SetLabelSize(0.055);
    h->GetYaxis()->SetTitleOffset(1.3);
    h->GetYaxis()->SetRangeUser(0.45, 1.68);
    h->SetTitle("");
    if (donch) h->GetXaxis()->SetTitle("n_{ch} / #LT n_{ch} #GT_{MB} ");
    else h->GetXaxis()->SetTitle("ZN / #LT ZN #GT_{MB}");
    h->GetYaxis()->SetTitle("#frac{#it{h}/#LT #it{h} #GT_{INEL>0}}{#it{n}_{ch}/#LT #it{n}_{ch} #GT_{INEL>0}}");
}

void prephisto2(TH1D* h = 0x0, bool donch = 1){
    h->SetStats(0);
    h->GetXaxis()->SetTitleSize(0.1);
    h->GetYaxis()->SetTitleSize(0.1);
    h->GetYaxis()->SetLabelSize(0.075);
    h->GetXaxis()->SetLabelSize(0.075);
    h->GetYaxis()->SetTitleOffset(0.96);
    h->GetYaxis()->SetRangeUser(0.45, 1.68);
    h->SetTitle("");
    if (donch) h->GetXaxis()->SetTitle("n_{ch} / #LT n_{ch} #GT_{MB} ");
    else h->GetXaxis()->SetTitle("ZN / #LT ZN #GT_{MB}");
    h->GetYaxis()->SetTitle("#frac{#it{h}/#LT #it{h} #GT_{INEL>0}}{#it{n}_{ch}/#LT #it{n}_{ch} #GT_{INEL>0}}");
}