// classes
enum classname
{
    kStandalone = 0,
    kHighMult,
    kLowMult,
    kHighZN,
    kLowZN,
    kVeryLowZN
};

void drawsgnloss(
    TString lCascType = "Xi",
    Int_t lClassCode = kStandalone,
    Bool_t DoMB = kFALSE,
    TString fWhichMultEstimator = "SPDClusters",
    TString fWhichEnergyEstimator = "V0M");

void init(int lClassCode,
          vector<Double_t> &percentileSPDClusters_low,
          vector<Double_t> &percentileSPDClusters_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          int &nbins,
          Bool_t DoMB);

void DrawNormCorrections(TString lPart = "Omega"){

    //drawsgnloss(lPart.Data(), 10);
    //drawsgnloss(lPart.Data(), kStandalone);
    //drawsgnloss(lPart.Data(), kHighMult);
    //drawsgnloss(lPart.Data(), kLowMult);
    //drawsgnloss(lPart.Data(), kHighZN);
    //drawsgnloss(lPart.Data(), kLowZN);
    //drawsgnloss(lPart.Data(), kStandalone, 1);

    drawsgnloss(lPart.Data(), kVeryLowZN);
}

void drawsgnloss(
    TString lCascType = "Xi",
    Int_t lClassCode = kStandalone,
    Bool_t DoMB = kFALSE,
    TString fWhichMultEstimator = "SPDClusters",
    TString fWhichEnergyEstimator = "V0M")
{

    TFile *file;
    if (!DoMB) {
        file = TFile::Open(Form("SignalLoss-%s-13TeV_class%i.root", lCascType.Data(), lClassCode));
    } else {
        file = TFile::Open(Form("SignalLoss-%s-13TeV_INELgt0.root", lCascType.Data()));
    }

    // Percentile
    vector<Double_t> percentileSPDClusters_low;
    vector<Double_t> percentileSPDClusters_high;
    vector<Double_t> percentileV0M_low;
    vector<Double_t> percentileV0M_high;
    int nbins;
    // initialize
    init(lClassCode, percentileSPDClusters_low, percentileSPDClusters_high, percentileV0M_low, percentileV0M_high, nbins, DoMB);
    const int percbinnumb = nbins;

    Int_t colors5[] = {kRed + 2, kOrange + 1, kSpring - 1, kAzure + 7, kBlue + 2};
    Int_t colors7[] = {kRed + 2, kOrange + 1, kYellow + 1, kSpring - 1, kGreen + 2, kAzure + 7, kBlue + 2};
    Int_t colors10[] = {kRed + 3, kRed + 1, kRed - 4, kOrange + 7, kYellow + 1, kSpring - 7, kGreen + 2, kAzure + 8, kBlue - 4, kBlue + 3};

    TString labels[] = {"I", "II","III", "IV", "V", "VI", "VII", "VIII", "IX", "X"};

    Int_t* colors;
    if (lClassCode == kStandalone){
        colors = colors10;
    } else if (lClassCode == 10){
        colors = colors5;
    } else if (lClassCode == kHighMult){
        colors = colors7;
    } else if (lClassCode == kLowMult){
        colors = colors7;
    } else if (lClassCode == kHighZN){
        colors = colors5;
    } else if (lClassCode == kLowZN){
        colors = colors5;
    }  else if (lClassCode == kVeryLowZN){
        colors = colors5;
    }

    TH1D* hsgnloss[percbinnumb];
    for (int iperc = 0; iperc < percbinnumb; iperc ++){
        cout << Form("hsgnloss_%s_%.0f-%.0f_%s_%.0f-%.0f", fWhichMultEstimator.Data(), percentileSPDClusters_low[iperc], percentileSPDClusters_high[iperc], fWhichEnergyEstimator.Data(), percentileV0M_low[iperc], percentileV0M_high[iperc]) << endl;

        hsgnloss[iperc] = (TH1D *)file->Get(Form("hsgnloss_%s_%.0f-%.0f_%s_%.0f-%.0f", fWhichMultEstimator.Data(), percentileSPDClusters_low[iperc], percentileSPDClusters_high[iperc], fWhichEnergyEstimator.Data(), percentileV0M_low[iperc], percentileV0M_high[iperc]));
        hsgnloss[iperc]->SetMarkerStyle(20);
        hsgnloss[iperc]->SetMarkerColor(colors[iperc]);
        hsgnloss[iperc]->SetLineColor(colors[iperc]);
        hsgnloss[iperc]->SetMarkerSize(2.);
    }

    TCanvas *c2 = new TCanvas(Form("c2"), "", 1200, 1100);
    c2->SetLeftMargin(0.15);
    c2->SetBottomMargin(0.12);
    c2->SetRightMargin(0.1);
    c2->SetTopMargin(0.1);
    c2->SetTicky();
    c2->SetTickx();

    TH1D *h = new TH1D("h", ";#it{p}_{T} (GeV/#it{c}); #varepsilon_{part}", 100, 0., 10.);
    h->SetStats(0);
    if(lCascType.Contains("Xi")) h->GetXaxis()->SetRangeUser(0.6,7.0);
    if(lCascType.Contains("Lambda")) h->GetXaxis()->SetRangeUser(0.4,8.0);
    if(lCascType.Contains("K0Short")) h->GetXaxis()->SetRangeUser(0.0,10.0);
    if(lCascType.Contains("Omega")) h->GetXaxis()->SetRangeUser(.9,6.0);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->Draw();
    h->SetMinimum(hsgnloss[percbinnumb - 1]->GetMinimum() * 0.9);
    h->SetMaximum(1.01);

    TLegend *leg = new TLegend(0.7, 0.15, 0.88, 0.45);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    TString classnames[5] = {"kStandalone", "kHighMult", "kLowMult", "kHighZN", "kLowZN"};
    if(lClassCode!=10) leg->SetHeader(classnames[lClassCode]);
    else leg->SetHeader("kStandalone");

    for (int iperc = 0; iperc < percbinnumb; iperc ++){
        hsgnloss[iperc]->Draw("same");
        leg->AddEntry(hsgnloss[iperc], labels[iperc], "LEP");
    }

    leg->Draw();

    TLatex *labbig = new TLatex();
    labbig->SetTextFont(42);
    labbig->SetNDC();
    labbig->SetTextColor(1);
    labbig->SetTextSize(0.07);
    labbig->SetTextAlign(22);
    labbig->SetTextAngle(0);

    if (lCascType.Contains("Xi")) labbig->DrawLatex(0.75, 0.55, "#Xi");
    if (lCascType.Contains("Lambda")) labbig->DrawLatex(0.75, 0.55, "#Lambda");
    if (lCascType.Contains("K0Short")) labbig->DrawLatex(0.75, 0.55, "K^{0}_{S}");
    if (lCascType.Contains("Omega")) labbig->DrawLatex(0.75, 0.55, "#Omega");

    if (!DoMB) c2->SaveAs(Form("images/sgnloss_%s_class%i.png", lCascType.Data(), lClassCode));
    else c2->SaveAs(Form("images/sgnloss_%s_INELgt0.png", lCascType.Data()));
}

void init(int lClassCode,
          vector<Double_t> &percentileSPDClusters_low,
          vector<Double_t> &percentileSPDClusters_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          int &nbins,
          Bool_t DoMB)
{

    // class 0 --> kStandalone
    Double_t percentileV0M_low_0[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
    Double_t percentileV0M_high_0[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    Double_t percentileSPDClusters_low_0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t percentileSPDClusters_high_0[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Long_t n0 = sizeof(percentileV0M_low_0) / sizeof(Double_t);

    // class 10 --> kStandalone for Omegas
    Double_t percentileV0M_low_10[] = {0, 5, 15, 30, 50};
    Double_t percentileV0M_high_10[] = {5, 15, 30, 50, 100};
    Double_t percentileSPDClusters_low_10[] = {0, 0, 0, 0, 0};
    Double_t percentileSPDClusters_high_10[] = {100, 100, 100, 100, 100};
    Long_t n10 = sizeof(percentileV0M_low_10) / sizeof(Double_t);

    // class 1 --> kHighMult
    Double_t percentileSPDClusters_low_1[] = {10, 10, 10, 10, 10, 10, 10};
    Double_t percentileSPDClusters_high_1[] = {20, 20, 20, 20, 20, 20, 20};
    Double_t percentileV0M_low_1[] = {0, 5, 10, 20, 30, 40, 50};
    Double_t percentileV0M_high_1[] = {5, 10, 20, 30, 40, 50, 100};
    Long_t n1 = sizeof(percentileSPDClusters_low_1) / sizeof(Double_t);

    // class 2 --> kLowMult
    Double_t percentileSPDClusters_low_2[] = {40, 40, 40, 40, 40, 40, 40};
    Double_t percentileSPDClusters_high_2[] = {50, 50, 50, 50, 50, 50, 50};
    Double_t percentileV0M_low_2[] = {0, 20, 30, 40, 50, 60, 70};
    Double_t percentileV0M_high_2[] = {20, 30, 40, 50, 60, 70, 100};
    Long_t n2 = sizeof(percentileSPDClusters_low_2) / sizeof(Double_t);

    // class 3 --> kHighZN
    Double_t percentileSPDClusters_low_3[] = {10, 40, 60, 70, 80};
    Double_t percentileSPDClusters_high_3[] = {40, 60, 70, 80, 100};
    Double_t percentileV0M_low_3[] = {70, 60, 40, 40, 40};
    Double_t percentileV0M_high_3[] = {100, 100, 100, 80, 70};
    Long_t n3 = sizeof(percentileSPDClusters_low_3) / sizeof(Double_t);

    // class 4 --> kLowZN
    Double_t percentileSPDClusters_low_4[] = {0, 10, 20, 30, 50};
    Double_t percentileSPDClusters_high_4[] = {20, 30, 40, 50, 100};
    Double_t percentileV0M_low_4[] = {40, 30, 30, 20, 0};
    Double_t percentileV0M_high_4[] = {60, 70, 50, 50, 30};
    Long_t n4 = sizeof(percentileSPDClusters_low_4) / sizeof(Double_t);

    // class 5 --> kVeryLowZN
    Double_t percentileSPDClusters_low_5[] = {0, 10, 20, 30};
    Double_t percentileSPDClusters_high_5[] = {10, 20, 30, 50};
    Double_t percentileV0M_low_5[] = {20, 10, 0, 0};
    Double_t percentileV0M_high_5[] = {30, 30, 20, 10};
    Long_t n5 = sizeof(percentileSPDClusters_low_5) / sizeof(Double_t);

    Double_t percentileSPDClusters_low_mb[] = {0};
    Double_t percentileSPDClusters_high_mb[] = {100};
    Double_t percentileV0M_low_mb[] = {0};
    Double_t percentileV0M_high_mb[] = {100};
    Long_t nmb = sizeof(percentileSPDClusters_low_mb) / sizeof(Double_t);

    if (DoMB)
    {
        nbins = nmb;
        for (Int_t i = 0; i < nmb; i++)
        {
            percentileSPDClusters_low.push_back(percentileSPDClusters_low_mb[i]);
            percentileSPDClusters_high.push_back(percentileSPDClusters_high_mb[i]);
            percentileV0M_low.push_back(percentileV0M_low_mb[i]);
            percentileV0M_high.push_back(percentileV0M_high_mb[i]);
        }
        cout << "------------------------------------------" << endl;
        cout << " Initializing MB class..." << endl;
    }
    else
    {

        if (lClassCode == kStandalone)
        {
            nbins = n0;
            for (Int_t i = 0; i < n0; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_0[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_0[i]);
                percentileV0M_low.push_back(percentileV0M_low_0[i]);
                percentileV0M_high.push_back(percentileV0M_high_0[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kStandalone class..." << endl;
        }

        if (lClassCode == 10)
        {
            nbins = n10;
            for (Int_t i = 0; i < n10; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_10[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_10[i]);
                percentileV0M_low.push_back(percentileV0M_low_10[i]);
                percentileV0M_high.push_back(percentileV0M_high_10[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kStandalone class for Omegas..." << endl;
        }

        if (lClassCode == kHighMult)
        {
            nbins = n1;
            for (Int_t i = 0; i < n1; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_1[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_1[i]);
                percentileV0M_low.push_back(percentileV0M_low_1[i]);
                percentileV0M_high.push_back(percentileV0M_high_1[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kHighMult class..." << endl;
        }

        if (lClassCode == kLowMult)
        {
            nbins = n2;
            for (Int_t i = 0; i < n2; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_2[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_2[i]);
                percentileV0M_low.push_back(percentileV0M_low_2[i]);
                percentileV0M_high.push_back(percentileV0M_high_2[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kLowMult class..." << endl;
        }

        if (lClassCode == kHighZN)
        {
            nbins = n3;
            for (Int_t i = 0; i < n3; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_3[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_3[i]);
                percentileV0M_low.push_back(percentileV0M_low_3[i]);
                percentileV0M_high.push_back(percentileV0M_high_3[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kHighZN class..." << endl;
        }

        if (lClassCode == kLowZN)
        {
            nbins = n4;
            for (Int_t i = 0; i < n4; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_4[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_4[i]);
                percentileV0M_low.push_back(percentileV0M_low_4[i]);
                percentileV0M_high.push_back(percentileV0M_high_4[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kLowZN class..." << endl;
        }

        if (lClassCode == kVeryLowZN)
        {
            nbins = n5;
            for (Int_t i = 0; i < n5; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_5[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_5[i]);
                percentileV0M_low.push_back(percentileV0M_low_5[i]);
                percentileV0M_high.push_back(percentileV0M_high_5[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kVeryLowZN class..." << endl;
        }
    }
}