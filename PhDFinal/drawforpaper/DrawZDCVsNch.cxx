void canvas_hestetics(TCanvas* c);
void graph(TGraphErrors* g, const char* x, const char*y);
void dostudy(TString period = "", TString sel = "SPD", float sel1 = 0., float sel2 = 100., Double_t* perc=0x0 , Long_t nbins=0 );

void DrawZDCVsNch(int domodel = 1){

    TFile *file[5];
    file[0] = TFile::Open("../yields/XiYields_SPDClustersV0M_class0_June23.root");
    file[1] = TFile::Open("../yields/XiYields_SPDClustersV0M_class1_June23.root");
    file[2] = TFile::Open("../yields/XiYields_SPDClustersV0M_class2_June23.root");
    file[3] = TFile::Open("../yields/XiYields_SPDClustersV0M_class4_June23.root");
    file[4] = TFile::Open("../yields/XiYields_SPDClustersV0M_class5_June23.root");

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

    gZDCVsNch[2]->SetMarkerStyle(kFullCircle);
    gZDCVsNch[2]->SetMarkerSize(2.3);
    gZDCVsNch[2]->SetMarkerColor(kBlue);
    gZDCVsNch[2]->SetLineColor(kBlue);

    gZDCVsNch[3]->SetMarkerStyle(kFullSquare);
    gZDCVsNch[3]->SetMarkerSize(2.3);
    gZDCVsNch[3]->SetMarkerColor(kViolet);
    gZDCVsNch[3]->SetLineColor(kViolet);

    gZDCVsNch[4]->SetMarkerStyle(kFullSquare);
    gZDCVsNch[4]->SetMarkerSize(2.3);
    gZDCVsNch[4]->SetMarkerColor(kGreen+1);
    gZDCVsNch[4]->SetLineColor(kGreen+1);

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
    h->GetYaxis()->SetRangeUser(0., 1.7);
    h->Draw();

    for (int i = 0; i < 5; i++){
        gZDCVsNch[i]->SetLineWidth(1);
        gZDCVsNch[i]->Draw("EP same");
    }

    TLegend *leg = new TLegend(0.6, 0.65, 0.86, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(gZDCVsNch[0], "V0M Standalone ", "p");
    leg->AddEntry(gZDCVsNch[2], "Low Multiplicity ", "p");
    leg->AddEntry(gZDCVsNch[1], "High Multiplicity ", "p");
    leg->AddEntry(gZDCVsNch[4], "Low ZN ", "p");
    leg->AddEntry(gZDCVsNch[3], "High ZN ", "p");
    leg->Draw("same");

    TLatex* label = new TLatex();
    label->SetTextFont(42);
    label->SetNDC();
    label->SetTextSize(0.04);
    label->SetTextAlign(22);
    label->SetTextAngle(0);
    label->DrawLatex(0.4, 0.88, "ALICE, pp #sqrt{#it{s}} = 13 TeV");

    // add a cross marker
    TMarker *marker = new TMarker(1.0, 1.0, 2);
    marker->SetMarkerSize(2.5);
    marker->SetMarkerColor(kBlack);

    TLegend *leg2 = new TLegend(0.25, 0.3, 0.4, 0.4);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.035);
    leg2->AddEntry(marker, "syst.", "p");
    leg2->Draw("same");


    c->SaveAs("images/ZDCVsNch.eps");


    //models


    TString ddname[] = {"Standalone", "kHighMult", "kLowMult", "kHighZN", "kLowZN"};
    const int nsel = sizeof(ddname) / sizeof(TString);

    TFile *filem = TFile::Open("../models/Results_Monash_SPDClustersV0M.root");
    TFile *filer = TFile::Open("../models/Results_Ropes_SPDClustersV0M.root");

    TGraphErrors *gNchLEm[nsel], *gNchLEr[nsel];
    for (int i = 0; i < nsel; i++)
    {
        gNchLEm[i] = (TGraphErrors *)filem->Get(Form("%s/NormLEVsNch", ddname[i].Data()));
        gNchLEr[i] = (TGraphErrors *)filer->Get(Form("%s/NormLEVsNch", ddname[i].Data()));
    }

    gNchLEm[0]->SetMarkerStyle(0);
    gNchLEm[0]->SetMarkerSize(0);
    gNchLEm[0]->SetMarkerColor(kGray+2);
    gNchLEm[0]->SetLineColor(kGray+2);
    gNchLEm[0]->SetLineStyle(7);
    gNchLEm[0]->SetLineWidth(2);

    gNchLEm[1]->SetMarkerStyle(0);
    gNchLEm[1]->SetMarkerSize(0);
    gNchLEm[1]->SetMarkerColor(kRed+1);
    gNchLEm[1]->SetLineColor(kRed+1);
    gNchLEm[1]->SetLineStyle(7);
    gNchLEm[1]->SetLineWidth(2);

    gNchLEm[2]->SetMarkerStyle(0);
    gNchLEm[2]->SetMarkerSize(0);
    gNchLEm[2]->SetMarkerColor(kBlue-3);
    gNchLEm[2]->SetLineColor(kBlue-3);
    gNchLEm[2]->SetLineStyle(7);
    gNchLEm[2]->SetLineWidth(2);

    gNchLEm[3]->SetMarkerStyle(0);
    gNchLEm[3]->SetMarkerSize(0);
    gNchLEm[3]->SetMarkerColor(kMagenta+1);
    gNchLEm[3]->SetLineColor(kMagenta+1);
    gNchLEm[3]->SetLineStyle(7);
    gNchLEm[3]->SetLineWidth(2);

    gNchLEm[4]->SetMarkerStyle(0);
    gNchLEm[4]->SetMarkerSize(0);
    gNchLEm[4]->SetMarkerColor(kGreen -2);
    gNchLEm[4]->SetLineColor(kGreen -2);
    gNchLEm[4]->SetLineStyle(7);
    gNchLEm[4]->SetLineWidth(2);

    gNchLEr[0]->SetMarkerStyle(0);
    gNchLEr[0]->SetMarkerSize(0);
    gNchLEr[0]->SetMarkerColor(kGray+2);
    gNchLEr[0]->SetLineColor(kGray+2);
    gNchLEr[0]->SetLineStyle(1);
    gNchLEr[0]->SetLineWidth(2);

    gNchLEr[1]->SetMarkerStyle(0);
    gNchLEr[1]->SetMarkerSize(0);
    gNchLEr[1]->SetMarkerColor(kRed+1);
    gNchLEr[1]->SetLineColor(kRed+1);
    gNchLEr[1]->SetLineStyle(1);
    gNchLEr[1]->SetLineWidth(2);

    gNchLEr[2]->SetMarkerStyle(0);
    gNchLEr[2]->SetMarkerSize(0);
    gNchLEr[2]->SetMarkerColor(kBlue-3);
    gNchLEr[2]->SetLineColor(kBlue-3);
    gNchLEr[2]->SetLineStyle(1);
    gNchLEr[2]->SetLineWidth(2);

    gNchLEr[3]->SetMarkerStyle(0);
    gNchLEr[3]->SetMarkerSize(0);
    gNchLEr[3]->SetMarkerColor(kMagenta+1);
    gNchLEr[3]->SetLineColor(kMagenta+1);
    gNchLEr[3]->SetLineStyle(1);
    gNchLEr[3]->SetLineWidth(2);

    gNchLEr[4]->SetMarkerStyle(0);
    gNchLEr[4]->SetMarkerSize(0);
    gNchLEr[4]->SetMarkerColor(kGreen -2);
    gNchLEr[4]->SetLineColor(kGreen -2);
    gNchLEr[4]->SetLineStyle(1);
    gNchLEr[4]->SetLineWidth(2);

    TCanvas *c1 = new TCanvas("c1", "", 1000, 900);
    canvas_hestetics(c1);
    h->Draw();

    for (int i = 0; i < 5; i++)
    {
        gZDCVsNch[i]->Draw("EP same");
    }

    leg->Draw("same");
    label->DrawLatex(0.4, 0.88, "ALICE, pp #sqrt{#it{s}} = 13 TeV");
    leg2->Draw("same");

    for (int i = 0; i < 5; i++)
    {
        gNchLEm[i]->Draw("C same");
        gNchLEr[i]->Draw("C same");
    }

    for (int i = 0; i < 5; i++)
    {
        gZDCVsNch[i]->Draw("EP same");
    }

    TLegend* legm = new TLegend(0.2, 0.2, 0.4, 0.3);
    legm->SetBorderSize(0);
    legm->SetFillStyle(0);
    legm->SetTextSize(0.035);
    legm->AddEntry(gNchLEm[0], "Pythia 8 Monash", "l");
    legm->AddEntry(gNchLEr[0], "Pythia 8 + Ropes", "l");
    legm->Draw("same");

    c1->SaveAs("images/ZDCVsNch_models.pdf");
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