void preppad(TPad *pad1, TPad *pad2);
void prepcanvas(TCanvas* c);
void prepgraph(TGraph2DErrors* g, int marker, float size, int color);

void DoGraph2D(int nclass = 1, TString num = "Xi", TString den = "nch")
{

    TString sel[] = {"_SPDClustersV0M_class0", Form("_SPDClustersV0M_class%i", nclass), Form("_SPDClustersV0M_class%i", 2), Form("_SPDClustersV0M_class%i", 3), Form("_SPDClustersV0M_class%i", 4), Form("_SPDClustersV0M_class%i", 5)};
    const int nsel = sizeof(sel)/sizeof(TString);

    TFile *fNum[nsel], *fDen[nsel];

    for (int i = 0; i < nsel; i++){
        fNum[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYields%s.root", num.Data(), sel[i].Data()));
        if (den.Contains("nch")){
            fDen[i] = fNum[i];
        } else {
            fDen[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYields%s.root", den.Data(), sel[i].Data()));
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
    TGraph2DErrors *gRatioStat[nsel];

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

        gRatioStat[i] = new TGraph2DErrors(gNumNchStat[i]->GetN(), &Nch[i][0], &ZDC[i][0], &RatioStat[i][0], &NchErr[i][0], &ZDCErr[i][0], &RatioStatErr[i][0]);
        gRatioStat[i]->SetName(Form("RatioStat%s", sel[i].Data()));
    }
    prepgraph(gRatioStat[0], kFullDiamond, 3.2, kBlack);
    prepgraph(gRatioStat[1], kFullCircle, 2.5, kRed);

    const int n = gRatioStat[0]->GetN() + gRatioStat[1]->GetN() + gRatioStat[2]->GetN() + gRatioStat[3]->GetN() + gRatioStat[4]->GetN() + gRatioStat[5]->GetN();
    Double_t x[n], y[n], z[n], ex[n], ey[n], ez[n];
    for (int i = 0; i< gRatioStat[0]->GetN(); i++){
        x[i] = gRatioStat[0]->GetX()[i];
        y[i] = - gRatioStat[0]->GetY()[i];
        z[i] = gRatioStat[0]->GetZ()[i];
    }
    for (int i = 0; i< gRatioStat[1]->GetN(); i++){
        x[i + gRatioStat[0]->GetN()] = gRatioStat[1]->GetX()[i];
        y[i + gRatioStat[0]->GetN()] = - gRatioStat[1]->GetY()[i];
        z[i + gRatioStat[0]->GetN()] = gRatioStat[1]->GetZ()[i];
    }
    for (int i = 0; i< gRatioStat[2]->GetN(); i++){
        x[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN()] = gRatioStat[2]->GetX()[i];
        y[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN()] = - gRatioStat[2]->GetY()[i];
        z[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN()] = gRatioStat[2]->GetZ()[i];
    }
    for (int i = 0; i< gRatioStat[3]->GetN(); i++){
        x[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN() + gRatioStat[2]->GetN()] = gRatioStat[3]->GetX()[i];
        y[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN() + gRatioStat[2]->GetN()] = - gRatioStat[3]->GetY()[i];
        z[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN() + gRatioStat[2]->GetN()] = gRatioStat[3]->GetZ()[i];
    }
    for (int i = 0; i< gRatioStat[4]->GetN(); i++){
        x[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN() + gRatioStat[2]->GetN() + gRatioStat[3]->GetN()] = gRatioStat[4]->GetX()[i];
        y[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN() + gRatioStat[2]->GetN() + gRatioStat[3]->GetN()] = - gRatioStat[4]->GetY()[i];
        z[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN() + gRatioStat[2]->GetN() + gRatioStat[3]->GetN()] = gRatioStat[4]->GetZ()[i];
    }
    for (int i = 0; i< gRatioStat[5]->GetN(); i++){
        x[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN() + gRatioStat[2]->GetN() + gRatioStat[3]->GetN() + gRatioStat[4]->GetN()] = gRatioStat[5]->GetX()[i];
        y[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN() + gRatioStat[2]->GetN() + gRatioStat[3]->GetN() + gRatioStat[4]->GetN()] = - gRatioStat[5]->GetY()[i];
        z[i + gRatioStat[0]->GetN() + gRatioStat[1]->GetN() + gRatioStat[2]->GetN() + gRatioStat[3]->GetN() + gRatioStat[4]->GetN()] = gRatioStat[5]->GetZ()[i];
    }

    TGraph2D* gRatioStatbis = new TGraph2D(n, x, y, z);

    new TCanvas;
    TH2D *hFrame = new TH2D("hFrame", "", 10, 0, 4,10,0,4);
    hFrame->SetStats(0);
    hFrame->SetMaximum(2);
    hFrame->Draw();
    gRatioStat[0]->GetXaxis()->SetTitle("n_{ch} / #LT n_{ch} #GT_{MB} ");
    gRatioStat[0]->GetYaxis()->SetTitle("ZN_{sum} / #LT ZN_{sum} #GT_{MB} ");
    gRatioStat[0]->GetZaxis()->SetTitle(Form("%s / %s", num.Data(), den.Data()));
    gRatioStat[0]->SetTitle("");
    gRatioStatbis->Draw("colz");

    TFile* write = new TFile("test.root", "RECREATE");
    gRatioStat[0]->Write();
    //gRatioStat[0]->Draw("EP SAME");
    //gRatioStat[1]->Draw("EP SAME");


/*    TH1D *hnch = new TH1D("hnch", "", 10, 0, 4.3);
    hnch->SetStats(0);
    hnch->GetXaxis()->SetTitle("n_{ch} / #LT n_{ch} #GT_{MB} ");
    hnch->GetXaxis()->SetTitleSize(0.06);
    hnch->GetYaxis()->SetTitleSize(0.06);
    hnch->GetYaxis()->SetLabelSize(0.04);
    hnch->GetXaxis()->SetLabelSize(0.04);
    hnch->GetYaxis()->SetTitleOffset(1.35);
    hnch->GetYaxis()->SetRangeUser(gNchRatioStat[0]->GetY()[gNchRatioStat[0]->GetN() - 1] * 0.6, gNchRatioStat[0]->GetY()[0] * 1.1);

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
    hzn->GetYaxis()->SetRangeUser(gNchRatioStat[0]->GetY()[gNchRatioStat[0]->GetN() - 1]*0.6, gNchRatioStat[0]->GetY()[0]*1.1);


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
    TGraphErrors *gclone[numb1];

    Int_t colors[] = {kRed + 1, kOrange+1, kYellow+1, kSpring - 1, kGreen+1, kAzure + 7, kBlue + 2};
    //Int_t colors[] = {kBlue + 2, kAzure + 7, kSpring - 1, kOrange + 1, kRed + 1};

    TLegend *lXi2 = new TLegend(0.55, 0.18, 0.65, 0.52);
    lXi2->SetBorderSize(0);
    lXi2->SetTextSize(0.04);
    lXi2->SetTextFont(42);
    TString multlabel[] = {"I", "II", "III", "IV", "V"};

    for (int i = 0; i < numb1; i++)
    {
        gclone[i] = new TGraphErrors(1, &Nch[1][i], &RatioStat[1][i], &NchErr[1][i], &RatioStatErr[1][i]);
        gclone[i]->SetMarkerStyle(kFullCircle);
        gclone[i]->SetMarkerSize(2.8);
        gclone[i]->SetMarkerColor(colors[i]);
        gclone[i]->SetLineColor(colors[i]);
        gclone[i]->Draw("P SAME");
        lXi2->AddEntry(gclone[i], Form("%s", multlabel[i].Data()), "P");
    }
    //lXi2->Draw("SAME");

    TMarker *fullcircle = new TMarker(12., 0.0074509863, 2);
    fullcircle->SetMarkerStyle(kFullCircle);
    fullcircle->SetMarkerSize(2.5);
    fullcircle->Draw("SAME");

    TLegend *lXi1 = new TLegend(0.65, 0.18, 0.87, 0.37);
    lXi1->SetBorderSize(0);
    lXi1->AddEntry(gNchRatioStat[0], "V0M standalone", "P");
    lXi1->AddEntry(gNchRatioStat[1], Form("class %i",nclass), "P");
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

    for (int i = 0; i < numb1; i++)
    {
        gclone_zdc[i] = new TGraphErrors(1, &ZDC[1][i], &RatioStat[1][i], &ZDCErr[1][i], &RatioStatErr[1][i]);
        gclone_zdc[i]->SetMarkerStyle(kFullCircle);
        gclone_zdc[i]->SetMarkerSize(2.8);
        gclone_zdc[i]->SetMarkerColor(colors[i]);
        gclone_zdc[i]->SetLineColor(colors[i]);
        gclone_zdc[i]->Draw("P SAME");
    }

    c->SaveAs(Form("plots/Ratio%s%s_%i.png", num.Data(), den.Data(), nclass));*/
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

void prepgraph(TGraph2DErrors* g, int marker, float size, int color){
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(size);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetFillColor(color);
    g->SetLineWidth(3);
    g->SetFillStyle(3001);
}