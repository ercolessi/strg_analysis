void preppad(TPad *pad1, TPad *pad2);
void prepcanvas(TCanvas* c);
void prepgraph(TGraphErrors* g, int marker, float size, int color);

void CompareResults0815(TString num = "K0Short", TString den = "nch")
{

    int iclass1 = 8;
    int iclass2 = 10;
    int imarker1 = kFullCircle;
    int imarker2 = kFullCircle;
    TString sel[] = {"_SPDtrk0815V0M_class0", Form("_SPDtrk0815V0M_class%i", 10), Form("_SPDClustersV0M_class%i", 12)};
    const int nsel = sizeof(sel)/sizeof(TString);

    TFile *fNum[nsel], *fDen[nsel];

    for (int i = 0; i < nsel; i++){
        fNum[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDWork/yields/%sYieldsRAW%s.root", num.Data(), sel[i].Data()));
        if (den.Contains("nch")){
            fDen[i] = fNum[i];
        } else {
            fDen[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDWork/yields/%sYieldsRAW%s.root", den.Data(), sel[i].Data()));
        }
    }

    TGraphErrors *gNumNchStat[nsel], *gNumZDCStat[nsel], *gDenNchStat[nsel];
    TGraphErrors *gNumNchSyst[nsel], *gNumZDCSyst[nsel], *gDenNchSyst[nsel];

    for (int i = 0; i < nsel; i++){
        gNumNchStat[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormNchStat");
        gNumZDCStat[i] = (TGraphErrors *)fNum[i]->Get("NormYieldsNormZDCSumStat");
        gDenNchStat[i] = (TGraphErrors *)fDen[i]->Get("NormYieldsNormNchStat");
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

    prepgraph(gNchRatioStat[1], imarker1, 2.5, kBlue);
    prepgraph(gZDCRatioStat[1], imarker1, 2.5, kBlue);

    prepgraph(gNchRatioStat[2], imarker2, 2.5, kRed);
    prepgraph(gZDCRatioStat[2], imarker2, 2.5, kRed);

    TH1D *hnch = new TH1D("hnch", "", 10, 0, 4.3);
    hnch->SetStats(0);
    hnch->GetXaxis()->SetTitle("n_{ch} / #LT n_{ch} #GT_{MB} ");
    hnch->GetXaxis()->SetTitleSize(0.06);
    hnch->GetYaxis()->SetTitleSize(0.06);
    hnch->GetYaxis()->SetLabelSize(0.04);
    hnch->GetXaxis()->SetLabelSize(0.04);
    hnch->GetYaxis()->SetTitleOffset(1.35);
    hnch->GetYaxis()->SetRangeUser(gNchRatioStat[0]->GetY()[gNchRatioStat[0]->GetN() - 1] * 0.95, gNchRatioStat[0]->GetY()[0] * 1.05);

    hnch->SetTitle("");
    if (num.Contains("Xi") && den.Contains("K0Short"))
        hnch->GetYaxis()->SetTitle("#frac{#Xi}{K^{0}_{S}}");
    if (num.Contains("Xi") && den.Contains("Lambda"))
        hnch->GetYaxis()->SetTitle("#frac{#Xi}{#Lambda}");
    if (num.Contains("Lambda") && den.Contains("K0Short"))
        hnch->GetYaxis()->SetTitle("#frac{#Lambda}{K^{0}_{S}}");
    if (num.Contains("Xi") && den.Contains("nch"))
        hnch->GetYaxis()->SetTitle("#frac{#Xi}{n_{ch}}");
    if (num.Contains("K0Short") && den.Contains("nch"))
        hnch->GetYaxis()->SetTitle("#frac{K^{0}_{S}}{n_{ch}}");
    if (num.Contains("Lambda") && den.Contains("nch"))
        hnch->GetYaxis()->SetTitle("#frac{#Lambda}{n_{ch}}");

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
    hzn->GetYaxis()->SetRangeUser(gNchRatioStat[0]->GetY()[gNchRatioStat[0]->GetN() - 1]*0.95, gNchRatioStat[0]->GetY()[0]*1.05);


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
    for (int j = 0; j < 1; j++) {
        gNchRatioStat[j]->Draw("SAME EP");
    }

    const int numb1 = gNumNchStat[1]->GetN();
    const int numb2 = gNumNchStat[2]->GetN();
    /*const int numb3 = gNumNchStat[3]->GetN();
    const int numb4 = gNumNchStat[4]->GetN();
    const int numb5 = gNumNchStat[5]->GetN();*/

    TGraphErrors *gclone[numb1];
    TGraphErrors *gclone2[numb2];
    /*TGraphErrors *gclone3[numb3];
    TGraphErrors *gclone4[numb4];
    TGraphErrors *gclone5[numb5];*/

    Int_t colors[] = {kRed + 1, kOrange+1, kSpring - 1, kAzure + 7, kBlue + 2};
    Int_t colors2[] = {kRed + 2, kRed, kOrange + 1, kYellow+1, kSpring - 1, kAzure + 7, kBlue, kBlue + 2};

    //Int_t colors[] = {kBlue + 2, kAzure + 7, kSpring - 1, kOrange + 1, kRed + 1};

    TLegend *lXi2 = new TLegend(0.55, 0.18, 0.65, 0.62);
    lXi2->SetBorderSize(0);
    lXi2->SetTextSize(0.04);
    lXi2->SetTextFont(42);
    TString multlabel[] = {"I", "II", "III", "IV", "V"};

    for (int i = 0; i < numb1; i++)
    {
        gclone[i] = new TGraphErrors(1, &Nch[1][i], &RatioStat[1][i], &NchErr[1][i], &RatioStatErr[1][i]);
        gclone[i]->SetMarkerStyle(imarker1);
        gclone[i]->SetMarkerSize(2.8);
        gclone[i]->SetMarkerColor(kBlue);
        gclone[i]->SetLineColor(kBlue);
        gclone[i]->Draw("P SAME");
        lXi2->AddEntry(gclone[i], Form("%s", multlabel[i].Data()), "P");
    }

    for (int i = 0; i < numb2; i++){
        gclone2[i] = new TGraphErrors(1, &Nch[2][i], &RatioStat[2][i], &NchErr[2][i], &RatioStatErr[2][i]);
        gclone2[i]->SetMarkerStyle(imarker2);
        gclone2[i]->SetMarkerSize(2.8);
        gclone2[i]->SetMarkerColor(kRed);
        gclone2[i]->SetLineColor(kRed);
        gclone2[i]->Draw("P SAME");
    }

  /*  for (int i = 0; i < numb3; i++){
        gclone3[i] = new TGraphErrors(1, &Nch[3][i], &RatioStat[3][i], &NchErr[3][i], &RatioStatErr[3][i]);
        gclone3[i]->SetMarkerStyle(33);
        gclone3[i]->SetMarkerSize(2.8);
        gclone3[i]->SetMarkerColor(colors2[i]);
        gclone3[i]->SetLineColor(colors2[i]);
        gclone3[i]->Draw("P SAME");
    }

    for (int i = 0; i < numb4; i++) {
        gclone4[i] = new TGraphErrors(1, &Nch[4][i], &RatioStat[4][i], &NchErr[4][i], &RatioStatErr[4][i]);
        gclone4[i]->SetMarkerStyle(34);
        gclone4[i]->SetMarkerSize(2.8);
        gclone4[i]->SetMarkerColor(colors2[i]);
        gclone4[i]->SetLineColor(colors2[i]);
        gclone4[i]->Draw("P SAME");
    }

    for (int i = 0; i < numb5; i++) {
        gclone5[i] = new TGraphErrors(1, &Nch[5][i], &RatioStat[5][i], &NchErr[5][i], &RatioStatErr[5][i]);
        gclone5[i]->SetMarkerStyle(47);
        gclone5[i]->SetMarkerSize(2.8);
        gclone5[i]->SetMarkerColor(colors2[i]);
        gclone5[i]->SetLineColor(colors2[i]);
        gclone5[i]->Draw("P SAME");
    }*/

    //lXi2->Draw("SAME");

    TMarker *fullcircle = new TMarker(12., 0.0074509863, 2);
    fullcircle->SetMarkerStyle(kFullCircle);
    fullcircle->SetMarkerSize(2.5);
    fullcircle->Draw("SAME");
    //
    TMarker *opencircle = new TMarker(12., 0.0074509863, 2);
    opencircle->SetMarkerStyle(kOpenCircle);
    opencircle->SetMarkerSize(2.5);
    opencircle->Draw("SAME");

    TMarker *fullcross = new TMarker(12., 0.0074509863, 2);
    fullcross->SetMarkerStyle(33);
    fullcross->SetMarkerSize(3);
    fullcross->Draw("SAME");
    //
    TMarker *opencross = new TMarker(12., 0.0074509863, 2);
    opencross->SetMarkerStyle(34);
    opencross->SetMarkerSize(3);
    opencross->Draw("SAME");

    TMarker *opencross2 = new TMarker(12., 0.0074509863, 2);
    opencross2->SetMarkerStyle(47);
    opencross2->SetMarkerSize(3);
    opencross2->Draw("SAME");

    TLegend *lXi1 = new TLegend(0.65, 0.18, 0.87, 0.47);
    lXi1->SetBorderSize(0);
    lXi1->AddEntry(gNchRatioStat[0], "V0M standalone", "P");
    lXi1->AddEntry(gNchRatioStat[1], "SPDtrk0815+V0M", "P");
    lXi1->AddEntry(gNchRatioStat[2], "SPDclusters+V0M", "P");
    /*lXi1->AddEntry(gNchRatioStat[5], "SPDcl[10-20]+V0M", "P");
    lXi1->AddEntry(gNchRatioStat[3], "SPDcl[20-30]+V0M", "P");
    lXi1->AddEntry(gNchRatioStat[4], "SPDcl[50-70]+V0M", "P");*/
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
  //  rap->DrawLatex(0.75, 0.45, Form("%s", "|#it{y}|<0.5"));

    //
    pad2->cd();
    hzn->Draw();
    gZDCRatioStat[0]->Draw("SAME EP");

    TGraphErrors *gclone_zdc[numb1];
    TGraphErrors *gclone_zdc2[numb2];
    /*TGraphErrors *gclone_zdc3[numb3];
    TGraphErrors *gclone_zdc4[numb4];
    TGraphErrors *gclone_zdc5[numb5];
*/
    for (int i = 0; i < numb1; i++)
    {
        gclone_zdc[i] = new TGraphErrors(1, &ZDC[1][i], &RatioStat[1][i], &ZDCErr[1][i], &RatioStatErr[1][i]);
        gclone_zdc[i]->SetMarkerStyle(imarker1);
        gclone_zdc[i]->SetMarkerSize(2.8);
        gclone_zdc[i]->SetMarkerColor(kBlue);
        gclone_zdc[i]->SetLineColor(kBlue);
        gclone_zdc[i]->Draw("P SAME");
    }

    for (int i = 0; i < numb2; i++){
        gclone_zdc2[i] = new TGraphErrors(1, &ZDC[2][i], &RatioStat[2][i], &ZDCErr[2][i], &RatioStatErr[2][i]);
        gclone_zdc2[i]->SetMarkerStyle(imarker2);
        gclone_zdc2[i]->SetMarkerSize(2.8);
        gclone_zdc2[i]->SetMarkerColor(kRed);
        gclone_zdc2[i]->SetLineColor(kRed);
        gclone_zdc2[i]->Draw("P SAME");
    }

    /*for (int i = 0; i < numb3; i++)
    {
        gclone_zdc3[i] = new TGraphErrors(1, &ZDC[3][i], &RatioStat[3][i], &ZDCErr[3][i], &RatioStatErr[3][i]);
        gclone_zdc3[i]->SetMarkerStyle(33);
        gclone_zdc3[i]->SetMarkerSize(2.8);
        gclone_zdc3[i]->SetMarkerColor(colors2[i]);
        gclone_zdc3[i]->SetLineColor(colors2[i]);
        gclone_zdc3[i]->Draw("P SAME");
    }

    for (int i = 0; i < numb4; i++)
    {
        gclone_zdc4[i] = new TGraphErrors(1, &ZDC[4][i], &RatioStat[4][i], &ZDCErr[4][i], &RatioStatErr[4][i]);
        gclone_zdc4[i]->SetMarkerStyle(34);
        gclone_zdc4[i]->SetMarkerSize(2.8);
        gclone_zdc4[i]->SetMarkerColor(colors2[i]);
        gclone_zdc4[i]->SetLineColor(colors2[i]);
        gclone_zdc4[i]->Draw("P SAME");
    }
    for (int i = 0; i < numb5; i++)
    {
        gclone_zdc5[i] = new TGraphErrors(1, &ZDC[5][i], &RatioStat[5][i], &ZDCErr[5][i], &RatioStatErr[5][i]);
        gclone_zdc5[i]->SetMarkerStyle(47);
        gclone_zdc5[i]->SetMarkerSize(2.8);
        gclone_zdc5[i]->SetMarkerColor(colors2[i]);
        gclone_zdc5[i]->SetLineColor(colors2[i]);
        gclone_zdc5[i]->Draw("P SAME");
    }
*/
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
    g->SetLineWidth(3);
    g->SetFillStyle(3001);
}