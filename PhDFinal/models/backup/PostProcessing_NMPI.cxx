void preppad(TPad *pad1, TPad *pad2, TPad *pad3);
void prepcanvas(TCanvas *c);
void prepgraph(TGraphErrors *g, int marker, float size, int color);

void PostProcessing_NMPI(TString num = "XiMinusXiPlus", TString den = "nch", TString pythia = "PythiaRopes")
{

    TFile *file = TFile::Open("Results_PythiaRopes_Train2628_SPDclV0M_3040.root");
    TString train = "Train2628";
    //if (pythia.Contains("PythiaRopes")) file = TFile::Open(Form("Test_Results_%s_%s_Nch0815V0M.root", pythia.Data(), train.Data()));
    //else if (pythia.Contains("PythiaMonash")) file = TFile::Open(Form("Test_Results_%s_%s_Nch0815V0M.root", pythia.Data(), "Train2627"));
    TString outname = Form("Ratio%s%s_%s.root", num.Data(), den.Data(), pythia.Data());

    TString ddname[] =
    {"Standalone", "SPDCl1020", "SPDCl3040", "V0M1020", "V0M4050"};
    //{"Standalone", "kLowMult", "kHighMult", "kHighZN", "kLowZN"};
    const int nsel = sizeof(ddname) / sizeof(TString);

    TGraphErrors *gNumNch[nsel], *gDenNch[nsel], *gNumZDC[nsel], *gDenZDC[nsel], *gNumNMPI[nsel];
    for (int i = 0; i < nsel; i++)
    {
        cout << i << endl;
        gNumNch[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsNchNorm_%s", ddname[i].Data(), num.Data()));
        gNumZDC[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsLeadENorm_%s", ddname[i].Data(), num.Data()));
        gNumNMPI[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsNMPI_%s", ddname[i].Data(), "XiMinus"));
        if (!den.Contains("nch")){
            gDenNch[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsNchNorm_%s", ddname[i].Data(), den.Data()));
            gDenZDC[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsLeadENorm_%s", ddname[i].Data(), den.Data()));
        } else {
            gDenNch[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsNchNorm_%s", ddname[i].Data(), num.Data()));
            gDenZDC[i] = (TGraphErrors *)file->Get(Form("%s/YieldRatioToMB/YieldRatioToMB_VsLeadENorm_%s", ddname[i].Data(), num.Data()));
        }
    }

    std::vector<Double_t> Nch[nsel], NchErr[nsel], ZDC[nsel], ZDCErr[nsel], NMPI[nsel], NMPIErr[nsel], YNum[nsel], YNumErr[nsel], YDen[nsel], YDenErr[nsel];
    std::vector<Double_t> Ratio[nsel], RatioErr[nsel];
    TGraphErrors *gNchRatio[nsel], *gZDCRatio[nsel], *gNMPIRatio[nsel];

    for (int i = 0; i < nsel; i++){
        cout << i << endl;
        if (gNumNch[i]->GetN() != gDenNch[i]->GetN()) {
            cout << "ERROR" << endl;
            return;
        }

        for (int j = 0; j < gNumNch[i]->GetN(); j++){
            Nch[i].push_back(gNumNch[i]->GetX()[j]);
            NchErr[i].push_back(gNumNch[i]->GetEX()[j]);
            ZDC[i].push_back(gNumZDC[i]->GetX()[j]);
            ZDCErr[i].push_back(gNumZDC[i]->GetEX()[j]);
            NMPI[i].push_back(gNumNMPI[i]->GetX()[j]);
            NMPIErr[i].push_back(gNumNMPI[i]->GetEX()[j]);
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
        gNMPIRatio[i] = new TGraphErrors(gNumNMPI[i]->GetN(), &NMPI[i][0], &Ratio[i][0], &NMPIErr[i][0], &RatioErr[i][0]);
    }

    prepgraph(gNchRatio[0], kFullDiamond, 2.5, kBlack);
    prepgraph(gZDCRatio[0], kFullDiamond, 2.5, kBlack);
    prepgraph(gNMPIRatio[0], kFullDiamond, 2.5, kBlack);

    const int numb1 = gNchRatio[1]->GetN();
    TGraphErrors *gclone1[numb1];
    const int numb2 = gNchRatio[2]->GetN();
    TGraphErrors *gclone2[numb2];
    const int numb3 = gNchRatio[3]->GetN();
    TGraphErrors *gclone3[numb3];
    const int numb4 = gNchRatio[4]->GetN();
    TGraphErrors *gclone4[numb4];

    TGraphErrors *zclone1[numb1];
    TGraphErrors *zclone2[numb2];
    TGraphErrors *zclone3[numb3];
    TGraphErrors *zclone4[numb4];

    TGraphErrors *mpiclone1[numb1];
    TGraphErrors *mpiclone2[numb2];
    TGraphErrors *mpiclone3[numb3];
    TGraphErrors *mpiclone4[numb4];

    Int_t colors[] = {kRed + 2, kOrange + 1, kYellow + 1, kSpring - 1, kAzure + 7, kBlue, kBlue + 2};

    /*for (int i = 1; i < nsel; i++){
        prepgraph(gNchRatio[i], kFullCircle, 2.3, kRed);
        prepgraph(gZDCRatio[i], kFullCircle, 2.3, kRed);
        prepgraph(gNMPIRatio[i], kFullCircle, 2.3, kRed);
    }*/

    TCanvas *c = new TCanvas("c", "", 2100, 600);
    TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 0.33, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0.33, 0, 0.66, 1);
    TPad *pad3 = new TPad("pad3", "pad3", 0.66, 0, 1., 1);
    prepcanvas(c);
    preppad(pad1, pad2, pad3);
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    //
    pad1->cd();
    gNchRatio[0]->Draw("AEP");
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

    gNchRatio[0]->GetXaxis()->SetTitle("n_{ch} / #LT n_{ch} #GT_{MB}");
    gNchRatio[0]->GetYaxis()->SetTitleSize(0.05);
    gNchRatio[0]->GetXaxis()->SetTitleSize(0.05);
    gNchRatio[0]->GetYaxis()->SetRangeUser(0.5,1.5);
    gNchRatio[0]->SetStats(0);
    gNchRatio[0]->SetTitle("");
    for (int i = 1; i < nsel; i++){
        //gNchRatio[i]->Draw("same LEP");
    }
    for (int i = 0; i < numb1; i++)
    {
        gclone1[i] = new TGraphErrors(1, &Nch[1][i], &Ratio[1][i], &NchErr[1][i], &RatioErr[1][i]);
        gclone1[i]->SetMarkerStyle(kFullCircle);
        gclone1[i]->SetMarkerSize(1.8);
        gclone1[i]->SetMarkerColor(colors[i]);
        gclone1[i]->SetLineColor(colors[i]);
        //gclone1[i]->Draw("EP SAME");
    }
    for (int i = 0; i < numb2; i++)
    {
        gclone2[i] = new TGraphErrors(1, &Nch[2][i], &Ratio[2][i], &NchErr[2][i], &RatioErr[2][i]);
        gclone2[i]->SetMarkerStyle(kOpenCircle);
        gclone2[i]->SetMarkerSize(1.8);
        gclone2[i]->SetMarkerColor(colors[i]);
        gclone2[i]->SetLineColor(colors[i]);
        //gclone2[i]->Draw("EP SAME");
    }
    for (int i = 0; i < numb3; i++)
    {
        gclone3[i] = new TGraphErrors(1, &Nch[3][i], &Ratio[3][i], &NchErr[3][i], &RatioErr[3][i]);
        gclone3[i]->SetMarkerStyle(kFullSquare);
        gclone3[i]->SetMarkerSize(1.8);
        gclone3[i]->SetMarkerColor(colors[i]);
        gclone3[i]->SetLineColor(colors[i]);
        gclone3[i]->Draw("EP SAME");
    }
    for (int i = 0; i < numb1; i++)
    {
        gclone4[i] = new TGraphErrors(1, &Nch[4][i], &Ratio[4][i], &NchErr[4][i], &RatioErr[4][i]);
        gclone4[i]->SetMarkerStyle(kOpenSquare);
        gclone4[i]->SetMarkerSize(1.8);
        gclone4[i]->SetMarkerColor(colors[i]);
        gclone4[i]->SetLineColor(colors[i]);
        gclone4[i]->Draw("EP SAME");
    }

    //
    pad2->cd();
    gZDCRatio[0]->Draw("AEP");
    gZDCRatio[0]->SetStats(0);
    gZDCRatio[0]->SetTitle("");
    for (int i = 1; i < nsel; i++){
        //gZDCRatio[i]->Draw("same LEP");
    }
    for (int i = 0; i < numb1; i++)
    {
        zclone1[i] = new TGraphErrors(1, &ZDC[1][i], &Ratio[1][i], &ZDCErr[1][i], &RatioErr[1][i]);
        zclone1[i]->SetMarkerStyle(kFullCircle);
        zclone1[i]->SetMarkerSize(1.8);
        zclone1[i]->SetMarkerColor(colors[i]);
        zclone1[i]->SetLineColor(colors[i]);
        //zclone1[i]->Draw("EP SAME");
    }
    for (int i = 0; i < numb2; i++)
    {
        zclone2[i] = new TGraphErrors(1, &ZDC[2][i], &Ratio[2][i], &ZDCErr[2][i], &RatioErr[2][i]);
        zclone2[i]->SetMarkerStyle(kOpenCircle);
        zclone2[i]->SetMarkerSize(1.8);
        zclone2[i]->SetMarkerColor(colors[i]);
        zclone2[i]->SetLineColor(colors[i]);
        //zclone2[i]->Draw("EP SAME");
    }
    for (int i = 0; i < numb3; i++)
    {
        zclone3[i] = new TGraphErrors(1, &ZDC[3][i], &Ratio[3][i], &ZDCErr[3][i], &RatioErr[3][i]);
        zclone3[i]->SetMarkerStyle(kFullSquare);
        zclone3[i]->SetMarkerSize(1.8);
        zclone3[i]->SetMarkerColor(colors[i]);
        zclone3[i]->SetLineColor(colors[i]);
        zclone3[i]->Draw("EP SAME");
    }
    for (int i = 0; i < numb1; i++)
    {
        zclone4[i] = new TGraphErrors(1, &ZDC[4][i], &Ratio[4][i], &ZDCErr[4][i], &RatioErr[4][i]);
        zclone4[i]->SetMarkerStyle(kOpenSquare);
        zclone4[i]->SetMarkerSize(1.8);
        zclone4[i]->SetMarkerColor(colors[i]);
        zclone4[i]->SetLineColor(colors[i]);
        zclone4[i]->Draw("EP SAME");
    }
    gZDCRatio[0]->GetXaxis()->SetTitle("E_{lead} / #LT E_{lead} #GT_{MB}");
    gZDCRatio[0]->GetYaxis()->SetTitleSize(0.05);
    gZDCRatio[0]->GetXaxis()->SetTitleSize(0.05);
    gZDCRatio[0]->GetYaxis()->SetRangeUser(0.5, 1.5);

    //
    pad3->cd();
    gNMPIRatio[0]->Draw("AEP");
    gNMPIRatio[0]->SetStats(0);
    gNMPIRatio[0]->SetTitle("");
    for (int i = 1; i < nsel; i++)
    {
      //  gNMPIRatio[i]->Draw("same LEP");
    }
    gNMPIRatio[0]->GetXaxis()->SetTitle("#LT MPI #GT");
    gNMPIRatio[0]->GetYaxis()->SetTitleSize(0.05);
    gNMPIRatio[0]->GetXaxis()->SetTitleSize(0.05);
    gNMPIRatio[0]->GetYaxis()->SetRangeUser(0.5, 1.5);
    for (int i = 0; i < numb1; i++)
    {
        mpiclone1[i] = new TGraphErrors(1, &NMPI[1][i], &Ratio[1][i], &NMPIErr[1][i], &RatioErr[1][i]);
        mpiclone1[i]->SetMarkerStyle(kFullCircle);
        mpiclone1[i]->SetMarkerSize(1.8);
        mpiclone1[i]->SetMarkerColor(colors[i]);
        mpiclone1[i]->SetLineColor(colors[i]);
        //mpiclone1[i]->Draw("EP SAME");
    }
    for (int i = 0; i < numb2; i++)
    {
        mpiclone2[i] = new TGraphErrors(1, &NMPI[2][i], &Ratio[2][i], &NMPIErr[2][i], &RatioErr[2][i]);
        mpiclone2[i]->SetMarkerStyle(kOpenCircle);
        mpiclone2[i]->SetMarkerSize(1.8);
        mpiclone2[i]->SetMarkerColor(colors[i]);
        mpiclone2[i]->SetLineColor(colors[i]);
        //mpiclone2[i]->Draw("EP SAME");
    }
    for (int i = 0; i < numb3; i++)
    {
        mpiclone3[i] = new TGraphErrors(1, &NMPI[3][i], &Ratio[3][i], &NMPIErr[3][i], &RatioErr[3][i]);
        mpiclone3[i]->SetMarkerStyle(kFullSquare);
        mpiclone3[i]->SetMarkerSize(1.8);
        mpiclone3[i]->SetMarkerColor(colors[i]);
        mpiclone3[i]->SetLineColor(colors[i]);
        mpiclone3[i]->Draw("EP SAME");
    }
    for (int i = 0; i < numb1; i++)
    {
        mpiclone4[i] = new TGraphErrors(1, &NMPI[4][i], &Ratio[4][i], &NMPIErr[4][i], &RatioErr[4][i]);
        mpiclone4[i]->SetMarkerStyle(kOpenSquare);
        mpiclone4[i]->SetMarkerSize(1.8);
        mpiclone4[i]->SetMarkerColor(colors[i]);
        mpiclone4[i]->SetLineColor(colors[i]);
        mpiclone4[i]->Draw("EP SAME");
    }

    c->SaveAs(Form("%sropes1_fixZN.png",num.Data()));
    /* TFile* write = new TFile(outname.Data(), "recreate");
     write->cd();
     gNchRatio[0]->Write("gNchRatio0");
     gNchRatio[1]->Write("gNchRatio1"); */
    // gNchRatio[2]->Write("gNchRatio2");
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

void preppad(TPad *pad1, TPad *pad2, TPad *pad3){
    pad1->SetBorderMode(0);
    pad1->SetTopMargin(0.05);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.15);
    pad1->SetRightMargin(0.0);
    pad2->SetRightMargin(0.0);
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
    pad3->SetBorderMode(0);
    pad3->SetTopMargin(0.05);
    pad3->SetBottomMargin(0.15);
    pad3->SetRightMargin(0.05);
    pad3->SetLeftMargin(0.0);
    pad3->SetFillColor(kWhite);
    pad3->SetTicky();
    pad3->SetTickx();
}