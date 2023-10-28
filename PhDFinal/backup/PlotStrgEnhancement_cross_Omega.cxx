void preppad(TPad *pad1, TPad *pad2);
void prepcanvas(TCanvas* c);
void prepgraph(TGraphErrors* g, int marker, float size, int color);

void PlotStrgEnhancement_cross_Omega(int ishighmult = 0,TString num = "Omega", TString den = "nch")
{
    int nclass1 = 11, nclass2 = 5;
    if (!ishighmult) {
        nclass1 = 12;
        nclass2 = 4;
    }
    Int_t imarker1 = kOpenCircle;
    Int_t imarker2 = kFullCircle;
    if (!ishighmult)
    {
        imarker1 = kOpenSquare;
        imarker2 = kFullSquare;
    }

    TString sel[] = {"_SPDClustersV0M_class10", Form("_SPDClustersV0M_class%i", nclass1), Form("_SPDClustersV0M_class%i", nclass2)};
    const int nsel = sizeof(sel)/sizeof(TString);

    TFile *fNum[nsel], *fDen[nsel];

    for (int i = 0; i < nsel; i++){
        fNum[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYieldsRAW%s.root", num.Data(), sel[i].Data()));
        if (den.Contains("nch")){
            fDen[i] = fNum[i];
        } else {
            fDen[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYieldsRAW%s.root", den.Data(), sel[i].Data()));
        }
    }

    TGraphErrors *gNumNchStat[nsel], *gNumZDCStat[nsel], *gDenNchStat[nsel];
    TGraphErrors *gNumNchSyst[nsel], *gNumZDCSyst[nsel], *gDenNchSyst[nsel];

    for (int i = 0; i < nsel; i++){
        if (den.Contains("nch")){
            gNumNchStat[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormNchStat");
            gNumZDCStat[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormZDCSumStat");
            gDenNchStat[i] = (TGraphErrors *)fDen[i]->Get("NormYieldsNormNchStat");
        }
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
            ZDC[i].push_back(gNumZDCStat[i]->GetX()[j]);
            ZDCErr[i].push_back(gNumZDCStat[i]->GetEX()[j]);
            YNum[i].push_back(gNumNchStat[i]->GetY()[j]);
            YNumErr[i].push_back(gNumNchStat[i]->GetEY()[j]);
            YDen[i].push_back(gDenNchStat[i]->GetY()[j]);
            YDenErr[i].push_back(gDenNchStat[i]->GetEY()[j]);
        }

        for (int j = 0; j < gNumNchStat[i]->GetN(); j++){
            if (den.Contains("nch")){
                RatioStat[i].push_back(YNum[i][j]);
                RatioStatErr[i].push_back(YNumErr[i][j]);
            } else {
                RatioStat[i].push_back(YNum[i][j] / YDen[i][j]);
                RatioStatErr[i].push_back(RatioStat[i][j] * TMath::Sqrt(TMath::Power(YNumErr[i][j] / YNum[i][j], 2) + TMath::Power(YDenErr[i][j] / YDen[i][j], 2)));
            }
        }

        gNchRatioStat[i] = new TGraphErrors(gNumNchStat[i]->GetN(), &Nch[i][0], &RatioStat[i][0], &NchErr[i][0], &RatioStatErr[i][0]);
        gZDCRatioStat[i] = new TGraphErrors(gNumNchStat[i]->GetN(), &ZDC[i][0], &RatioStat[i][0], &ZDCErr[i][0], &RatioStatErr[i][0]);
    }
    prepgraph(gNchRatioStat[0], kFullDiamond, 3.2, kBlack);
    prepgraph(gZDCRatioStat[0], kFullDiamond, 3.2, kBlack);
    prepgraph(gNchRatioStat[1], imarker1, 2.5, kBlack);
    prepgraph(gZDCRatioStat[1], imarker1, 2.5, kBlack);
    prepgraph(gNchRatioStat[2], imarker2, 2.5, kBlack);
    prepgraph(gZDCRatioStat[2], imarker2, 2.5, kBlack);

    TH1D *hnch = new TH1D("hnch", "", 10, 0, 4.3);
    hnch->SetStats(0);
    hnch->GetXaxis()->SetTitle("n_{ch} / #LT n_{ch} #GT_{MB} ");
    hnch->GetXaxis()->SetTitleSize(0.06);
    hnch->GetYaxis()->SetTitleSize(0.06);
    hnch->GetYaxis()->SetLabelSize(0.04);
    hnch->GetXaxis()->SetLabelSize(0.04);
    hnch->GetYaxis()->SetTitleOffset(1.35);
    hnch->GetYaxis()->SetRangeUser(gNchRatioStat[0]->GetY()[gNchRatioStat[0]->GetN() - 1] * 0.6, gNchRatioStat[0]->GetY()[0] * 1.2);

    hnch->SetTitle("");
    if (num.Contains("Xi") && den.Contains("K0Short"))
        hnch->GetYaxis()->SetTitle("#frac{#Xi}{K^{0}_{S}}");
    if (num.Contains("Xi") && den.Contains("Lambda"))
        hnch->GetYaxis()->SetTitle("#frac{#Xi}{#Lambda}");
    if (num.Contains("Lambda") && den.Contains("K0Short"))
        hnch->GetYaxis()->SetTitle("#frac{#Lambda}{K^{0}_{S}}");
    if (num.Contains("Xi") && den.Contains("nch"))
        hnch->GetYaxis()->SetTitle("#frac{#Xi / #LT#Xi#GT_{MB}}{n_{ch} / #LTn_{ch}#GT_{MB}}");
    if (num.Contains("K0Short") && den.Contains("nch"))
        hnch->GetYaxis()->SetTitle("#frac{K^{0}_{S} / #LTK^{0}_{S}#GT_{MB}}{n_{ch} / #LTn_{ch}#GT_{MB}}");
    if (num.Contains("Lambda") && den.Contains("nch"))
        hnch->GetYaxis()->SetTitle("#frac{#Lambda / #LT#Lambda#GT_{MB}}{n_{ch} / #LTn_{ch}#GT_{MB}}");
    if (num.Contains("Omega") && den.Contains("nch"))
        hnch->GetYaxis()->SetTitle("#frac{#Omega / #LT#Omega#GT_{MB}}{n_{ch} / #LTn_{ch}#GT_{MB}}");

    TH1D *hzn = new TH1D("hzn", "", 10, 0.1, 1.5);
    hzn->SetStats(0);
    hzn->GetXaxis()->SetTitle("ZN / #LT ZN #GT_{MB}");
    hzn->GetYaxis()->SetTitle("");
    hzn->GetXaxis()->SetTitleSize(0.06);
    hzn->GetYaxis()->SetTitleSize(0.06);
    hzn->GetYaxis()->SetLabelSize(0.04);
    hzn->GetXaxis()->SetLabelSize(0.04);
    hzn->GetYaxis()->SetTitleOffset(1.35);
    hzn->SetTitle("");
    hzn->GetYaxis()->SetLabelColor(kWhite);
    hzn->GetYaxis()->SetRangeUser(gNchRatioStat[0]->GetY()[gNchRatioStat[0]->GetN() - 1]*0.6, gNchRatioStat[0]->GetY()[0]*1.2);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    gNchRatioStat[0]->Draw("SAME EP");

    const int numb1 = gNumNchStat[1]->GetN();
    const int numb2 = gNumNchStat[2]->GetN();

    TGraphErrors *gclone[numb1];
    TGraphErrors *gclone2[numb2];
    TGraphErrors *gclonesyst[numb1];
    TGraphErrors *gclone2syst[numb2];

    Int_t colors[] = {kRed + 1, kOrange+1, kYellow+1, kSpring - 1, kGreen+2, kAzure + 7, kBlue + 2};
    Int_t colors2[] = {kRed+1, kOrange+1, kSpring-1, kAzure+7, kBlue+2};

    TLegend *lXi2 = new TLegend(0.55, 0.18, 0.65, 0.52);
    lXi2->SetBorderSize(0);
    lXi2->SetTextSize(0.04);
    lXi2->SetTextFont(42);
    TString multlabel[] = {"I", "II", "III", "IV", "V", "VI", "VII"};

    for (int i = 0; i < numb1; i++)
    {
        gclone[i] = new TGraphErrors(1, &Nch[1][i], &RatioStat[1][i], &NchErr[1][i], &RatioStatErr[1][i]);
        gclone[i]->SetMarkerStyle(imarker1);
        gclone[i]->SetMarkerSize(2.8);
        gclone[i]->SetMarkerColor(colors[i]);
        gclone[i]->SetLineColor(colors[i]);
        gclone[i]->Draw("P SAME");
        lXi2->AddEntry(gclone[i], Form("%s", multlabel[i].Data()), "P");
    }
    for (int i = 0; i < numb2; i++)
    {
        gclone2[i] = new TGraphErrors(1, &Nch[2][i], &RatioStat[2][i], &NchErr[2][i], &RatioStatErr[2][i]);
        gclone2[i]->SetMarkerStyle(imarker2);
        gclone2[i]->SetMarkerSize(2.8);
        gclone2[i]->SetMarkerColor(colors2[i]);
        gclone2[i]->SetLineColor(colors2[i]);
        gclone2[i]->Draw("P SAME");
        //      lXi2->AddEntry(gclone2[i], Form("%s", multlabel[i].Data()), "P");
    }
    //lXi2->Draw("SAME");

    TMarker *fullcircle = new TMarker(12., 0.0074509863, 2);
    fullcircle->SetMarkerStyle(imarker1);
    fullcircle->SetMarkerSize(2.5);
    fullcircle->Draw("SAME");

    TLegend *lXi1 = new TLegend(0.65, 0.18, 0.87, 0.37);
    lXi1->SetBorderSize(0);
    lXi1->AddEntry(gNchRatioStat[0], "V0M standalone", "P");
    //lXi1->AddEntry(gNchRatioStat[1], Form("class %i",nclass1), "P");
    //lXi1->AddEntry(gNchRatioStat[2], Form("class %i", nclass2), "P");
    lXi1->SetTextSize(0.04);
    lXi1->SetTextFont(42);
    lXi1->Draw("SAME");

    TLatex *rap = new TLatex();
    rap->SetTextFont(42);
    rap->SetNDC();
    rap->SetTextColor(1);
    rap->SetTextSize(0.045);
    rap->SetTextAlign(22);
    rap->SetTextAngle(0);
    rap->DrawLatex(0.75, 0.45, Form("%s", "|#it{y}|<0.5"));

    //
    pad2->cd();
    hzn->Draw();
    gZDCRatioStat[0]->Draw("SAME EP");

    TGraphErrors *gclone_zdc[numb1];
    TGraphErrors *gclone_zdc2[numb2];
    TGraphErrors *gclone_zdcsyst[numb1];
    TGraphErrors *gclone_zdc2syst[numb2];

    for (int i = 0; i < numb1; i++)
    {
        gclone_zdc[i] = new TGraphErrors(1, &ZDC[1][i], &RatioStat[1][i], &ZDCErr[1][i], &RatioStatErr[1][i]);
        gclone_zdc[i]->SetMarkerStyle(imarker1);
        gclone_zdc[i]->SetMarkerSize(2.8);
        gclone_zdc[i]->SetMarkerColor(colors[i]);
        gclone_zdc[i]->SetLineColor(colors[i]);
        gclone_zdc[i]->Draw("P SAME");
    }
    for (int i = 0; i < numb2; i++)
    {
        gclone_zdc2[i] = new TGraphErrors(1, &ZDC[2][i], &RatioStat[2][i], &ZDCErr[2][i], &RatioStatErr[2][i]);
        gclone_zdc2[i]->SetMarkerStyle(imarker2);
        gclone_zdc2[i]->SetMarkerSize(2.8);
        gclone_zdc2[i]->SetMarkerColor(colors2[i]);
        gclone_zdc2[i]->SetLineColor(colors2[i]);
        gclone_zdc2[i]->Draw("P SAME");
    }

    c->SaveAs(Form("plots/Ratio%s%s_%i_cross.pdf", num.Data(), den.Data(), ishighmult));
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
    g->SetFillColor(color);
    g->SetLineWidth(1);
    g->SetFillStyle(0);
}