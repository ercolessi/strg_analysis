void canvas_hestetics(TCanvas* c);
void graph(TGraphErrors* g, const char* x, const char*y);
void dostudy(TString period = "", TString sel = "SPD", float sel1 = 0., float sel2 = 100., Double_t* perc=0x0 , Long_t nbins=0 );

void DrawZDCVsNch(int domodel = 0){

    TFile *file[5];
    file[0] = TFile::Open("~/strg_analysis/PhDFinal/yields/XiYields_SPDClustersV0M_class0_June23.root");
    file[1] = TFile::Open("~/strg_analysis/PhDFinal/yields/XiYields_SPDClustersV0M_class1_June23.root");
    file[2] = TFile::Open("~/strg_analysis/PhDFinal/yields/XiYields_SPDClustersV0M_class2_June23.root");
    file[3] = TFile::Open("~/strg_analysis/PhDFinal/yields/XiYields_SPDClustersV0M_class4_June23.root");
    file[4] = TFile::Open("~/strg_analysis/PhDFinal/yields/XiYields_SPDClustersV0M_class5_June23.root");

    TGraphAsymmErrors* gNch[5];
    TGraphAsymmErrors* gZDC[5];

    for (int i = 0; i < 5; i++){
        gNch[i] = (TGraphAsymmErrors *)file[i]->Get("NormYieldsNormNchSyst");
        gZDC[i] = (TGraphAsymmErrors *)file[i]->Get("NormYieldsNormZDCSumSyst");
    }

    std::vector<Double_t> Nch[5], NchErr[5], ZDC[5], ZDCErr[5];
    TGraphErrors *gZDCVsNch[5];

    for (int i = 0; i < 5; i++){

        for (int j = 0; j < gNch[i]->GetN(); j++)
        {
            Nch[i].push_back(gNch[i]->GetX()[j]);
            NchErr[i].push_back(gNch[i]->GetEXhigh()[j]);
            ZDC[i].push_back(gZDC[i]->GetX()[j]);
            ZDCErr[i].push_back(gZDC[i]->GetEXhigh()[j]);
        }

        gZDCVsNch[i] = new TGraphErrors(Nch[i].size(), &Nch[i][0], &ZDC[i][0], &NchErr[i][0], &ZDCErr[i][0]);

    }

    for (int i = 0; i < 5; i++){
        graph(gZDCVsNch[i], "#LT n_{ch} #GT / #LT n_{ch} #GT_{INEL>0}", "#LT ZN #GT / #LT ZN #GT_{INEL>0}");
    }

    gZDCVsNch[0]->SetMarkerStyle(kFullDiamond);
    gZDCVsNch[0]->SetMarkerSize(2.8);
    gZDCVsNch[0]->SetMarkerColor(kBlack);
    gZDCVsNch[0]->SetLineColor(kBlack);

    gZDCVsNch[1]->SetMarkerStyle(kFullCircle);
    gZDCVsNch[1]->SetMarkerSize(2.3);
    gZDCVsNch[1]->SetMarkerColor(kRed);
    gZDCVsNch[1]->SetLineColor(kRed);

    gZDCVsNch[2]->SetMarkerStyle(kOpenCircle);
    gZDCVsNch[2]->SetMarkerSize(2.3);
    gZDCVsNch[2]->SetMarkerColor(kOrange+1);
    gZDCVsNch[2]->SetLineColor(kOrange+1);

    gZDCVsNch[3]->SetMarkerStyle(kOpenSquare);
    gZDCVsNch[3]->SetMarkerSize(2.3);
    gZDCVsNch[3]->SetMarkerColor(kAzure+7);
    gZDCVsNch[3]->SetLineColor(kAzure+7);

    gZDCVsNch[4]->SetMarkerStyle(kFullSquare);
    gZDCVsNch[4]->SetMarkerSize(2.3);
    gZDCVsNch[4]->SetMarkerColor(kBlue+1);
    gZDCVsNch[4]->SetLineColor(kBlue+1);

    TCanvas* c = new TCanvas("c", "", 1000, 900);
    canvas_hestetics(c);
    TH1F* h = new TH1F("h", "", 10, 0., 4.2);
    h->SetStats(0);
    h->GetXaxis()->SetTitle(" n_{ch}  / #LT n_{ch} #GT_{INEL>0}");
    h->GetYaxis()->SetTitle(" ZN  / #LT ZN #GT_{INEL>0}");
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.);
    h->GetXaxis()->SetTitleOffset(1.);
    h->GetYaxis()->SetRangeUser(0., 1.8);
    h->Draw();

    for (int i = 0; i < 5; i++){
        gZDCVsNch[i]->SetLineWidth(1);
        gZDCVsNch[i]->Draw("EP same");
    }

    TLegend *leg = new TLegend(0.6, 0.65, 0.86, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(gZDCVsNch[0], "Standalone ", "ep");
    leg->AddEntry(gZDCVsNch[2], "Low Multiplicity ", "ep");
    leg->AddEntry(gZDCVsNch[1], "High Multiplicity ", "ep");
    leg->AddEntry(gZDCVsNch[4], "Low ZN ", "ep");
    leg->AddEntry(gZDCVsNch[3], "High ZN ", "ep");
    leg->Draw("same");

    TLatex* label = new TLatex();
    label->SetTextFont(42);
    label->SetNDC();
    label->SetTextSize(0.04);
    label->SetTextAlign(22);
    label->SetTextAngle(0);
    label->DrawLatex(0.4, 0.88, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
    label->DrawLatex(0.3, 0.83, "This work");

    // add a cross marker
    TMarker *marker = new TMarker(1.0, 1.0, 2);
    marker->SetMarkerSize(2.5);
    marker->SetMarkerColor(kBlack);

    TLegend *leg2 = new TLegend(0.25, 0.3, 0.4, 0.4);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.035);
    leg2->AddEntry(marker, "syst.", "p");
    //leg2->Draw("same");


    c->SaveAs("ZDCVsNch.pdf");

    TCanvas *c1 = new TCanvas("c1", "", 1000, 900);
    canvas_hestetics(c1);
    h->Draw();

    for (int i = 0; i < 1; i++)
    {
        gZDCVsNch[i]->SetLineWidth(1);
        gZDCVsNch[i]->Draw("EP same");
    }

    TLegend *leg1 = new TLegend(0.6, 0.8, 0.86, 0.88);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.035);
    leg1->AddEntry(gZDCVsNch[0], "Standalone ", "ep");
    leg1->Draw("same");

     label->DrawLatex(0.4, 0.88, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
    label->DrawLatex(0.3, 0.83, "This work");


    c1->SaveAs("ZDCVsNch_std.pdf");
}

void graph(TGraphErrors* g, const char* x, const char*y){
    g->SetMarkerStyle(kFullCircle);
    g->SetMarkerSize(2.);
    g->SetMarkerColor(kBlack);
    g->SetLineColor(kBlack);
    g->SetLineWidth(2);
    g->GetYaxis()->SetTitle(y);
    g->GetXaxis()->SetTitle(x);
    g->SetTitle("");
}

void canvas_hestetics(TCanvas* c) {
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.16);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.05);
    c->SetTicky();
    c->SetTickx();
    //c->Divide(2,2);
}