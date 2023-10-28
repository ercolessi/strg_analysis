void preppad(TPad *pad1, TPad *pad2);
void prepcanvas(TCanvas *c);
void prepgraph(TGraphErrors *g, int marker, float size, int color);

void PostProcessing(TString num = "Proton", TString den = "nch", TString pythia = "Monash")
{

    TFile *file;
    if (pythia.Contains("Ropes")) file = TFile::Open("Results_Ropes_SPDClustersV0M.root");
    else if (pythia.Contains("Monash")) file = TFile::Open("Results_Monash_SPDClustersV0M.root");

    TString outname = Form("ParticleRatios/Ratio%s%s_%s.root", num.Data(), den.Data(), pythia.Data());

    TString ddname[] = {"Standalone", "kLowMult", "kHighMult", "kHighZN", "kLowZN"};
    const int nsel = sizeof(ddname) / sizeof(TString);

    TGraphErrors *gNumNch[nsel], *gDenNch[nsel], *gNumZDC[nsel], *gDenZDC[nsel];
    for (int i = 0; i < nsel; i++)
    {
        gNumNch[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsNchNorm_%s", ddname[i].Data(), num.Data()));
        gNumZDC[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsLeadENorm_%s", ddname[i].Data(), num.Data()));
        if (!den.Contains("nch")){
            gDenNch[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsNchNorm_%s", ddname[i].Data(), den.Data()));
            gDenZDC[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsLeadENorm_%s", ddname[i].Data(), den.Data()));
        } else {
            gDenNch[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsNchNorm_%s", ddname[i].Data(), num.Data()));
            gDenZDC[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsLeadENorm_%s", ddname[i].Data(), num.Data()));
        }
    }

    std::vector<Double_t> Nch[nsel], NchErr[nsel], ZDC[nsel], ZDCErr[nsel], YNum[nsel], YNumErr[nsel], YDen[nsel], YDenErr[nsel];
    std::vector<Double_t> Ratio[nsel], RatioErr[nsel];
    TGraphErrors *gNchRatio[nsel], *gZDCRatio[nsel];

    for (int i = 0; i < nsel; i++){
        if (gNumNch[i]->GetN() != gDenNch[i]->GetN()) {
            cout << "ERROR" << endl;
            return;
        }

        for (int j = 0; j < gNumNch[i]->GetN(); j++){
            Nch[i].push_back(gNumNch[i]->GetX()[j]);
            NchErr[i].push_back(gNumNch[i]->GetEX()[j]);
            ZDC[i].push_back(gNumZDC[i]->GetX()[j]);
            ZDCErr[i].push_back(gNumZDC[i]->GetEX()[j]);
            YNum[i].push_back(gNumNch[i]->GetY()[j]);
            YNumErr[i].push_back(gNumNch[i]->GetEY()[j]);
            YDen[i].push_back(gDenNch[i]->GetY()[j]);
            YDenErr[i].push_back(gDenNch[i]->GetEY()[j]);
        }

        for (int j = 0; j < gNumNch[i]->GetN(); j++){
            if (den.Contains("nch")){
                Ratio[i].push_back(YNum[i][j] / Nch[i][j]);
                RatioErr[i].push_back(Ratio[i][j] * TMath::Sqrt(TMath::Power(YNumErr[i][j] / YNum[i][j], 2) + TMath::Power(NchErr[i][j] / Nch[i][j], 2)));
            } else {
                Ratio[i].push_back(YNum[i][j]/YDen[i][j]);
                RatioErr[i].push_back(Ratio[i][j] * TMath::Sqrt(TMath::Power(YNumErr[i][j]/YNum[i][j], 2) + TMath::Power(YDenErr[i][j]/YDen[i][j], 2)));
            }
        }

        gNchRatio[i] = new TGraphErrors(gNumNch[i]->GetN(), &Nch[i][0], &Ratio[i][0], &NchErr[i][0], &RatioErr[i][0]);
        gZDCRatio[i] = new TGraphErrors(gNumNch[i]->GetN(), &ZDC[i][0], &Ratio[i][0], &ZDCErr[i][0], &RatioErr[i][0]);

        gNchRatio[i]->SetName(Form("gNchRatio_%s", ddname[i].Data()));
        gZDCRatio[i]->SetName(Form("gZDCRatio_%s", ddname[i].Data()));
    }

    prepgraph(gNchRatio[0], kFullDiamond, 3., kBlack);
    prepgraph(gZDCRatio[0], kFullDiamond, 3., kBlack);
    for (int i = 1; i < nsel; i++){
        prepgraph(gNchRatio[i], kFullCircle, 2.3, kRed);
        prepgraph(gZDCRatio[i], kFullCircle, 2.3, kRed);
    }

    TCanvas *c = new TCanvas("c", "", 1800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 0.54, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0.54, 0, 1., 1);
    prepcanvas(c);
    preppad(pad1, pad2);
    pad1->Draw();
    pad2->Draw();
    //
    pad1->cd();
    gNchRatio[0]->Draw();
    if (num.Contains("Xi") && den.Contains("K0Short"))
        gNchRatio[0]->GetYaxis()->SetTitle("#Xi/K^{0}_{S} w.r.t. MB");
    else if (num.Contains("Xi") && den.Contains("Lambda"))
        gNchRatio[0]->GetYaxis()->SetTitle("#Xi/#Lambda w.r.t. MB");
    if (num.Contains("Lambda") && den.Contains("K0Short"))
        gNchRatio[0]->GetYaxis()->SetTitle("#Lambda/K^{0}_{S} w.r.t. MB");
    else if (num.Contains("Xi") && den.Contains("nch"))
        gNchRatio[0]->GetYaxis()->SetTitle("#Xi/n_{ch} w.r.t. MB");
    else if (num.Contains("Lambda") && den.Contains("nch"))
        gNchRatio[0]->GetYaxis()->SetTitle("#Lambda/n_{ch} w.r.t. MB");
    else if (num.Contains("K0Short") && den.Contains("nch"))
        gNchRatio[0]->GetYaxis()->SetTitle("K^{0}_{S}/n_{ch} w.r.t. MB");
    else if (num.Contains("AntiProton") && den.Contains("nch"))
        gNchRatio[0]->GetYaxis()->SetTitle("#bar{p}/n_{ch} w.r.t. MB");
    else if (num.Contains("Phi") && den.Contains("nch"))
        gNchRatio[0]->GetYaxis()->SetTitle("#phi/n_{ch} w.r.t. MB");
    else if (num.Contains("KaMinus") && den.Contains("nch"))
        gNchRatio[0]->GetYaxis()->SetTitle("K^{-}/n_{ch} w.r.t. MB");
    else if (num.Contains("Proton") && den.Contains("nch"))
        gNchRatio[0]->GetYaxis()->SetTitle("p/n_{ch} w.r.t. MB");

    gNchRatio[0]->GetXaxis()->SetTitle("n_{ch} / #LT n_{ch} #GT_{MB}");
    gNchRatio[0]->GetYaxis()->SetTitleSize(0.05);
    gNchRatio[0]->GetXaxis()->SetTitleSize(0.05);
    gNchRatio[0]->GetYaxis()->SetRangeUser(0.5,1.5);
    gNchRatio[0]->SetStats(0);
    gNchRatio[0]->SetTitle("");

    TLegend *leg = new TLegend(0.25, 0.75, 0.5, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gNchRatio[0], "V0 Standalone", "lep");

    Int_t colors[] = {kBlue, kRed, kViolet, kGreen+1};
    for (int i = 1; i < nsel; i++){
        gNchRatio[i]->SetLineColor(colors[i-1]);
        gNchRatio[i]->SetMarkerColor(colors[i - 1]);
        leg->AddEntry(gNchRatio[i], Form("%s", ddname[i].Data()), "lep");

        gNchRatio[i]->Draw("same LEP");
    }
    leg->Draw();

    //
    pad2->cd();
    gZDCRatio[0]->Draw();
    gZDCRatio[0]->SetStats(0);
    gZDCRatio[0]->SetTitle("");
    for (int i = 1; i < nsel; i++){
        gZDCRatio[i]->SetLineColor(colors[i-1]);
        gZDCRatio[i]->SetMarkerColor(colors[i - 1]);

        gZDCRatio[i]->Draw("same LEP");
    }
    gZDCRatio[0]->GetXaxis()->SetTitle("E_{lead} / #LT E_{lead} #GT_{MB}");
    gZDCRatio[0]->GetYaxis()->SetTitleSize(0.05);
    gZDCRatio[0]->GetXaxis()->SetTitleSize(0.05);
    gZDCRatio[0]->GetYaxis()->SetRangeUser(0.5, 1.5);

    TLatex *tex = new TLatex(0.2, 0.85, Form("Pythia %s", pythia.Data()));
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.05);
    tex->Draw();

    c->SaveAs(Form("DoubleRatio_%s_Pythia%s_FixedNchVsZN.png",num.Data(), pythia.Data()));

    TFile* write = new TFile(outname.Data(), "recreate");
    write->cd();
    for (int i = 0; i < nsel; i++){
        gNchRatio[i]->Write();
        gZDCRatio[i]->Write();
    }
}

void prepgraph(TGraphErrors* g, int marker, float size, int color){
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(size);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetFillStyle(0);
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