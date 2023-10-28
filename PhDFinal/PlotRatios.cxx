void preppad(TPad *pad1, TPad *pad2);
void prepcanvas(TCanvas* c);
void prepgraph(TGraphErrors* g = 0x0, int marker = 20, float size = 2., int color = kBlack, int linewidth = 1, int linestyle = 1, int fillstyle = 0);

void PlotRatios(int ishighmult = 1, TString num = "Xi", TString den = "nch", bool domodel = 1)
{
    int nclass1 = 1, nclass2 = 5;
    if (!ishighmult) {
        nclass1 = 2;
        nclass2 = 4;
    }
    Int_t imarker1 = kFullCircle;
    Int_t imarker2 = kFullSquare;
    Int_t icolor1 = kRed;
    Int_t icolor2 = kGreen+1;
    if (!ishighmult)
    {
        imarker1 = kFullCircle;
        imarker2 = kFullSquare;
        icolor1 = kBlue;
        icolor2 = kViolet;
    }

    int step = 3; //1 only data, 2 monash, 3 ropes

    TString sel[] = {"_SPDClustersV0M_class0", Form("_SPDClustersV0M_class%i", nclass1), Form("_SPDClustersV0M_class%i", nclass2)};
    const int nsel = sizeof(sel)/sizeof(TString);

    TFile *fNum[nsel], *fDen[nsel];

    for (int i = 0; i < nsel; i++){
        fNum[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYields%s_June23.root", num.Data(), sel[i].Data()));
        if (den.Contains("nch")){
            fDen[i] = fNum[i];
        } else {
            fDen[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYields%s_June23.root", den.Data(), sel[i].Data()));
        }
    }

    TGraphErrors *gNumNchStat[nsel], *gNumZDCStat[nsel], *gDenNchStat[nsel];
    TGraphErrors *gNumNchSyst[nsel], *gNumZDCSyst[nsel], *gDenNchSyst[nsel];

    for (int i = 0; i < nsel; i++){
        if (den.Contains("nch")){
            gNumNchStat[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormNchStat");
            gNumZDCStat[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormZDCSumStat");
            gDenNchStat[i] = (TGraphErrors *)fDen[i]->Get("NormYieldsNormNchStat");
            gNumNchSyst[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormNchSyst");
            gNumZDCSyst[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormZDCSumSyst");
            gDenNchSyst[i] = (TGraphErrors *)fDen[i]->Get("NormYieldsNormNchSyst");
        } else {
            gNumNchStat[i] = (TGraphErrors *)fNum[i]->Get("YieldsToMBNormNchStat");
            gNumZDCStat[i] = (TGraphErrors *)fNum[i]->Get("YieldsToMBNormZDCSumStat");
            gDenNchStat[i] = (TGraphErrors *)fDen[i]->Get("YieldsToMBNormNchStat");
            gNumNchSyst[i] = (TGraphErrors *)fNum[i]->Get("YieldsToMBNormNchSyst");
            gNumZDCSyst[i] = (TGraphErrors *)fNum[i]->Get("YieldsToMBNormZDCSumSyst");
            gDenNchSyst[i] = (TGraphErrors *)fDen[i]->Get("YieldsToMBNormNchSyst");
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
            if (den.Contains("nch")){
                RatioStat[i].push_back(YNum[i][j]);
                RatioStatErr[i].push_back(YNumErr[i][j]);
                RatioSyst[i].push_back(YNum[i][j]);
                RatioSystErr[i].push_back(YNumErrSyst[i][j]);
            } else {
                RatioStat[i].push_back(YNum[i][j] / YDen[i][j]);
                RatioStatErr[i].push_back(RatioStat[i][j] * TMath::Sqrt(TMath::Power(YNumErr[i][j] / YNum[i][j], 2) + TMath::Power(YDenErr[i][j] / YDen[i][j], 2)));
                RatioSyst[i].push_back(YNum[i][j] / YDen[i][j]);
                RatioSystErr[i].push_back(RatioSyst[i][j] * TMath::Sqrt(TMath::Power(YNumErrSyst[i][j] / YNum[i][j], 2) + TMath::Power(YDenErrSyst[i][j] / YDen[i][j], 2)));
            }
        }

        gNchRatioStat[i] = new TGraphErrors(gNumNchStat[i]->GetN(), &Nch[i][0], &RatioStat[i][0], &NchErr[i][0], &RatioStatErr[i][0]);
        gZDCRatioStat[i] = new TGraphErrors(gNumNchStat[i]->GetN(), &ZDC[i][0], &RatioStat[i][0], &ZDCErr[i][0], &RatioStatErr[i][0]);
        gNchRatioSyst[i] = new TGraphErrors(gNumNchSyst[i]->GetN(), &Nch[i][0], &RatioSyst[i][0], &NchErrSyst[i][0], &RatioSystErr[i][0]);
        gZDCRatioSyst[i] = new TGraphErrors(gNumNchSyst[i]->GetN(), &ZDC[i][0], &RatioSyst[i][0], &ZDCErrSyst[i][0], &RatioSystErr[i][0]);
    }
    prepgraph(gNchRatioStat[0], kFullDiamond, 3.2, kBlack);
    prepgraph(gZDCRatioStat[0], kFullDiamond, 3.2, kBlack);
    prepgraph(gNchRatioSyst[0], kFullDiamond, 0, kBlack);
    prepgraph(gZDCRatioSyst[0], kFullDiamond, 0, kBlack);
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
    //hnch->GetYaxis()->SetRangeUser(gNchRatioStat[0]->GetY()[gNchRatioStat[0]->GetN() - 1] * 0.9, gNchRatioStat[0]->GetY()[0] * 1.1);
    hnch->GetYaxis()->SetRangeUser(0.4,1.5);

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
    //hzn->GetYaxis()->SetRangeUser(gNchRatioStat[0]->GetY()[gNchRatioStat[0]->GetN() - 1]*0.9, gNchRatioStat[0]->GetY()[0]*1.1);
    hzn->GetYaxis()->SetRangeUser(0.4,1.5);

    TFile *file1, *file2;
    TString num_ = num;
    TString den_ = den;
    if (num.Contains("Xi")) num_ = "XiMinusXiPlus";
    if (den.Contains("Xi")) den_ = "XiMinusXiPlus";
    if (num.Contains("Lambda")) num_ = "LambdaAntiLambda";
    if (den.Contains("Lambda")) den_ = "LambdaAntiLambda";
    if (num.Contains("K0Short")) num_ = "KaPlusKaMinus";
    if (den.Contains("K0Short")) den_ = "KaPlusKaMinus";

    TString classmc1 = "kHighMult";
    TString classmc2 = "kLowZN";
    if (!ishighmult) {
        classmc1 = "kLowMult";
        classmc2 = "kHighZN";
    }
    TGraphErrors *gnchclass0;
    TGraphErrors *gnchclass1;
    TGraphErrors *gnchclass2;
    TGraphErrors *gzdcclass0;
    TGraphErrors *gzdcclass1;
    TGraphErrors *gzdcclass2;

    TGraphErrors *bisgnchclass0;
    TGraphErrors *bisgnchclass1;
    TGraphErrors *bisgnchclass2;
    TGraphErrors *bisgzdcclass0;
    TGraphErrors *bisgzdcclass1;
    TGraphErrors *bisgzdcclass2;

    TGraphErrors *legendgnchclass0;
    TGraphErrors *legendgnchclass1;
    TGraphErrors *legendgnchclass2;

    if (domodel){

        file1 = TFile::Open(Form("models/ParticleRatios/Ratio%s%s_%s.root", num_.Data(), den_.Data(), "Ropes"));
        gnchclass0 = (TGraphErrors *)file1->Get(Form("gNchRatio_%s", "Standalone"));
        gnchclass1 = (TGraphErrors *)file1->Get(Form("gNchRatio_%s", classmc1.Data()));
        gnchclass2 = (TGraphErrors *)file1->Get(Form("gNchRatio_%s", classmc2.Data()));
        gzdcclass0 = (TGraphErrors *)file1->Get(Form("gZDCRatio_%s", "Standalone"));
        gzdcclass1 = (TGraphErrors *)file1->Get(Form("gZDCRatio_%s", classmc1.Data()));
        gzdcclass2 = (TGraphErrors *)file1->Get(Form("gZDCRatio_%s", classmc2.Data()));
        file2 = TFile::Open(Form("models/ParticleRatios/Ratio%s%s_%s.root", num_.Data(), den_.Data(), "Monash"));
        bisgnchclass0 = (TGraphErrors *)file2->Get(Form("gNchRatio_%s", "Standalone"));
        bisgnchclass1 = (TGraphErrors *)file2->Get(Form("gNchRatio_%s", classmc1.Data()));
        bisgnchclass2 = (TGraphErrors *)file2->Get(Form("gNchRatio_%s", classmc2.Data()));
        bisgzdcclass0 = (TGraphErrors *)file2->Get(Form("gZDCRatio_%s", "Standalone"));
        bisgzdcclass1 = (TGraphErrors *)file2->Get(Form("gZDCRatio_%s", classmc1.Data()));
        bisgzdcclass2 = (TGraphErrors *)file2->Get(Form("gZDCRatio_%s", classmc2.Data()));

        legendgnchclass0 = (TGraphErrors *)file1->Get(Form("gNchRatio_%s", "Standalone"));
        legendgnchclass1 = (TGraphErrors *)file1->Get(Form("gNchRatio_%s", classmc1.Data()));
        legendgnchclass2 = (TGraphErrors *)file1->Get(Form("gNchRatio_%s", classmc2.Data()));

        prepgraph(gnchclass0, 0, 0, kBlack, 2, 1, 0);
        prepgraph(gnchclass1, 0, 0, icolor1, 2, 1, 0);
        prepgraph(gnchclass2, 0, 0, icolor2, 2, 1, 0);
        prepgraph(gzdcclass0, 0, 0, kBlack, 2, 1, 0);
        prepgraph(gzdcclass1, 0, 0, icolor1, 2, 1, 0);
        prepgraph(gzdcclass2, 0, 0, icolor2, 2, 1, 0);

        prepgraph(bisgnchclass0, 0, 0, kBlack, 2, 3, 0);
        prepgraph(bisgnchclass1, 0, 0, icolor1, 2, 3, 0);
        prepgraph(bisgnchclass2, 0, 0, icolor2, 2, 3, 0);
        prepgraph(bisgzdcclass0, 0, 0, kBlack, 2, 3, 0);
        prepgraph(bisgzdcclass1, 0, 0, icolor1, 2, 3, 0);
        prepgraph(bisgzdcclass2, 0, 0, icolor2, 2, 3, 0);

        prepgraph(legendgnchclass0, 0, 0, kBlack, 2, 1, 0);
        prepgraph(legendgnchclass1, 0, 0, kBlack, 2, 7, 0);
        prepgraph(legendgnchclass2, 0, 0, kBlack, 2, 3, 0);
    }
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

    if (domodel)
    {
        //gnchclass0->Draw("SAME C");
        //gnchclass1->Draw("SAME C");
        //gnchclass2->Draw("SAME C");

        //bisgnchclass0->Draw("SAME C");
        //bisgnchclass1->Draw("SAME C");
        //bisgnchclass2->Draw("SAME C");
    }

    gNchRatioSyst[0]->Draw("SAME E2");
    gNchRatioStat[0]->Draw("SAME EP");

    const int numb1 = gNumNchStat[1]->GetN();
    const int numb2 = gNumNchStat[2]->GetN();

    TGraphErrors *gclone[numb1];
    TGraphErrors *gclone2[numb2];
    TGraphErrors *gclonesyst[numb1];
    TGraphErrors *gclone2syst[numb2];

    //Int_t colors[] = {kRed + 1, kOrange+1, kYellow+1, kSpring - 1, kGreen+2, kAzure + 7, kBlue + 2};
    //Int_t colors2[] = {kRed+1, kOrange+1, kYellow-3, kAzure+7, kBlue+2};
    Int_t colors[] = {icolor1, icolor1, icolor1, icolor1, icolor1, icolor1, icolor1, icolor1};
    Int_t colors2[] = {icolor2, icolor2, icolor2, icolor2, icolor2, icolor2, icolor2, icolor2};

    TLegend *lXi2 = new TLegend(0.55, 0.18, 0.65, 0.52);
    lXi2->SetBorderSize(0);
    lXi2->SetTextSize(0.04);
    lXi2->SetTextFont(42);
    TString multlabel[] = {"I", "II", "III", "IV", "V", "VI", "VII"};

    for (int i = 0; i < numb1; i++)
    {
        gclone[i] = new TGraphErrors(1, &Nch[1][i], &RatioStat[1][i], &NchErr[1][i], &RatioStatErr[1][i]);
        gclone[i]->SetMarkerStyle(imarker1);
        gclone[i]->SetMarkerSize(2.5);
        gclone[i]->SetMarkerColor(colors[i]);
        gclone[i]->SetLineColor(colors[i]);
        gclonesyst[i] = new TGraphErrors(1, &Nch[1][i], &RatioSyst[1][i], &NchErrSyst[1][i], &RatioSystErr[1][i]);
        gclonesyst[i]->SetMarkerStyle(imarker1);
        gclonesyst[i]->SetMarkerSize(0);
        gclonesyst[i]->SetMarkerColor(colors[i]);
        gclonesyst[i]->SetLineColor(colors[i]);
        gclonesyst[i]->SetFillStyle(0);
        gclonesyst[i]->Draw("E2 SAME");
        gclone[i]->Draw("P SAME");
        lXi2->AddEntry(gclone[i], Form("%s", multlabel[i].Data()), "P");
    }
    for (int i = 0; i < numb2; i++)
    {
        gclone2[i] = new TGraphErrors(1, &Nch[2][i], &RatioStat[2][i], &NchErr[2][i], &RatioStatErr[2][i]);
        gclone2[i]->SetMarkerStyle(imarker2);
        gclone2[i]->SetMarkerSize(2.5);
        gclone2[i]->SetMarkerColor(colors2[i]);
        gclone2[i]->SetLineColor(colors2[i]);
        gclone2syst[i] = new TGraphErrors(1, &Nch[2][i], &RatioSyst[2][i], &NchErrSyst[2][i], &RatioSystErr[2][i]);
        gclone2syst[i]->SetMarkerStyle(imarker2);
        gclone2syst[i]->SetMarkerSize(0);
        gclone2syst[i]->SetMarkerColor(colors2[i]);
        gclone2syst[i]->SetLineColor(colors2[i]);
        gclone2syst[i]->SetFillStyle(0);
        gclone2syst[i]->Draw("E2 SAME");
        gclone2[i]->Draw("P SAME");
        //      lXi2->AddEntry(gclone2[i], Form("%s", multlabel[i].Data()), "P");
    }
    //lXi2->Draw("SAME");

    if (domodel)
    {
        if (step > 1) {
            bisgnchclass0->Draw("SAME C");
            bisgnchclass1->Draw("SAME C");
            bisgnchclass2->Draw("SAME C");
            if (step>2){
            gnchclass0->Draw("SAME C");
            gnchclass1->Draw("SAME C");
            gnchclass2->Draw("SAME C");
            }
        }
    }

    TMarker *marker1 = new TMarker(12., 0.0074509863, 2);
    marker1->SetMarkerStyle(imarker1);
    marker1->SetMarkerColor(icolor1);
    marker1->SetMarkerSize(2.5);
    marker1->Draw("SAME");

    TMarker *marker2 = new TMarker(12., 0.0074509863, 2);
    marker2->SetMarkerStyle(imarker2);
    marker2->SetMarkerColor(icolor2);
    marker2->SetMarkerSize(2.5);
    marker2->Draw("SAME");

    TLegend *lXi1 = new TLegend(0.55, 0.18, 0.87, 0.4);
    lXi1->SetBorderSize(0);
    lXi1->SetHeader("Data:");
    lXi1->AddEntry(gNchRatioStat[0], "V0M standalone", "P");
    if (ishighmult){
    lXi1->AddEntry(marker1, "SPDcl+V0M, high #LT n_{ch} #GT", "P");
    lXi1->AddEntry(marker2, "SPDcl+V0M, low #LT ZN #GT", "P");
    } else {
    lXi1->AddEntry(marker1, "SPDcl+V0M, low #LT n_{ch} #GT", "P");
    lXi1->AddEntry(marker2, "SPDcl+V0M, high #LT ZN #GT", "P");
    }
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
    rap->DrawLatex(0.41, 0.89, Form("%s", "ALICE, pp #sqrt{s} = 13 TeV"));

    //
    pad2->cd();
    hzn->Draw();

    gZDCRatioSyst[0]->Draw("SAME E2");
    gZDCRatioStat[0]->Draw("SAME EP");

    TGraphErrors *gclone_zdc[numb1];
    TGraphErrors *gclone_zdc2[numb2];
    TGraphErrors *gclone_zdcsyst[numb1];
    TGraphErrors *gclone_zdc2syst[numb2];

    for (int i = 0; i < numb1; i++)
    {
        gclone_zdc[i] = new TGraphErrors(1, &ZDC[1][i], &RatioStat[1][i], &ZDCErr[1][i], &RatioStatErr[1][i]);
        gclone_zdc[i]->SetMarkerStyle(imarker1);
        gclone_zdc[i]->SetMarkerSize(2.5);
        gclone_zdc[i]->SetMarkerColor(colors[i]);
        gclone_zdc[i]->SetLineColor(colors[i]);
        gclone_zdcsyst[i] = new TGraphErrors(1, &ZDC[1][i], &RatioSyst[1][i], &ZDCErrSyst[1][i], &RatioSystErr[1][i]);
        gclone_zdcsyst[i]->SetMarkerStyle(imarker1);
        gclone_zdcsyst[i]->SetMarkerSize(0);
        gclone_zdcsyst[i]->SetMarkerColor(colors[i]);
        gclone_zdcsyst[i]->SetLineColor(colors[i]);
        gclone_zdcsyst[i]->SetFillStyle(0);
        gclone_zdcsyst[i]->Draw("E2 SAME");
        gclone_zdc[i]->Draw("P SAME");
    }
    for (int i = 0; i < numb2; i++)
    {
        gclone_zdc2[i] = new TGraphErrors(1, &ZDC[2][i], &RatioStat[2][i], &ZDCErr[2][i], &RatioStatErr[2][i]);
        gclone_zdc2[i]->SetMarkerStyle(imarker2);
        gclone_zdc2[i]->SetMarkerSize(2.5);
        gclone_zdc2[i]->SetMarkerColor(colors2[i]);
        gclone_zdc2[i]->SetLineColor(colors2[i]);
        gclone_zdc2syst[i] = new TGraphErrors(1, &ZDC[2][i], &RatioSyst[2][i], &ZDCErrSyst[2][i], &RatioSystErr[2][i]);
        gclone_zdc2syst[i]->SetMarkerStyle(imarker2);
        gclone_zdc2syst[i]->SetMarkerSize(0);
        gclone_zdc2syst[i]->SetMarkerColor(colors2[i]);
        gclone_zdc2syst[i]->SetLineColor(colors2[i]);
        gclone_zdc2syst[i]->SetFillStyle(0);
        gclone_zdc2syst[i]->Draw("E2 SAME");
        gclone_zdc2[i]->Draw("P SAME");
    }

    if (domodel)
    {
        if (step > 1)
        {
            bisgzdcclass0->Draw("SAME C");
            bisgzdcclass1->Draw("SAME C");
            bisgzdcclass2->Draw("SAME C");

            if (step > 2){
            gzdcclass0->Draw("SAME C");
            gzdcclass1->Draw("SAME C");
            gzdcclass2->Draw("SAME C");
            }
        }
    }

    TLegend *lmodel = new TLegend(0.03, 0.18, 0.32, 0.5);
    lmodel->SetBorderSize(0);
    lmodel->AddEntry(legendgnchclass2, "Pythia8 Monash", "L");
    if (step > 2)
        lmodel->AddEntry(legendgnchclass0, "Pythia8 + Ropes", "L");
    lmodel->AddEntry(gnchclass0, "V0M standalone", "L");
    lmodel->AddEntry(gnchclass1, "SPDcl+V0M, fixed #LT n_{ch} #GT", "L");
    lmodel->AddEntry(gnchclass2, "SPDcl+V0M, fixed #LT E_{lead} #GT", "L");
    lmodel->SetTextSize(0.04);
    lmodel->SetTextFont(42);
    if (step > 1) lmodel->Draw("SAME");

    c->SaveAs(Form("plots/Ratio%s%s_%i_cross_model_step%i.pdf", num.Data(), den.Data(), ishighmult,step));
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