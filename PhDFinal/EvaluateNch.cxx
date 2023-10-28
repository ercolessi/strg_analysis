void Init(Int_t lClassCode,
          TString lWhichEnergyEstimator,
          TString lWhichMultEstimator,
          vector<Double_t> &percentileMult_low,
          vector<Double_t> &percentileMult_high,
          vector<Double_t> &percentileEnergy_low,
          vector<Double_t> &percentileEnergy_high,
          vector<Double_t> &nch,
          vector<Double_t> &nchErr);

void EvaluateNch(){


    //Get official Nch
    // Percentile
    vector<Double_t> percentileMult_low;
    vector<Double_t> percentileMult_high;
    vector<Double_t> percentileEnergy_low;
    vector<Double_t> percentileEnergy_high;
    vector<Double_t> offNch;
    vector<Double_t> offNchErr;


    Init(0, "V0M", "SPDClusters", percentileMult_low, percentileMult_high, percentileEnergy_low, percentileEnergy_high, offNch, offNchErr);
    Init(1, "V0M", "SPDClusters", percentileMult_low, percentileMult_high, percentileEnergy_low, percentileEnergy_high, offNch, offNchErr);
    Init(2, "V0M", "SPDClusters", percentileMult_low, percentileMult_high, percentileEnergy_low, percentileEnergy_high, offNch, offNchErr);
    Init(3, "V0M", "SPDClusters", percentileMult_low, percentileMult_high, percentileEnergy_low, percentileEnergy_high, offNch, offNchErr);
    Init(4, "V0M", "SPDClusters", percentileMult_low, percentileMult_high, percentileEnergy_low, percentileEnergy_high, offNch, offNchErr);
    const int binnumber = percentileMult_low.size();

    //Get My Nch
    // Nch ZDC file
    TFile *NchZDCfile = TFile::Open("Selections/NchRawContainer_reference.root");

    Double_t myNch[100];
    Double_t myNchErr[100];
    TH3D *NchHisto = (TH3D *)NchZDCfile->Get("hspd_spdv0m");

    int miny = 1, maxy = 1, minz = 1, maxz = 1;
    for (int i = 0; i < binnumber; i++)
    {
        miny = NchHisto->GetYaxis()->FindBin(percentileMult_low[i] + 1e-10);
        maxy = NchHisto->GetYaxis()->FindBin(percentileMult_high[i]);
        minz = NchHisto->GetZaxis()->FindBin(percentileEnergy_low[i] + 1e-10);
        maxz = NchHisto->GetZaxis()->FindBin(percentileEnergy_high[i]);

        TH1D *dummynch = NchHisto->ProjectionX(Form("dummynch%i", i), miny, maxy, minz, maxz);
        myNch[i] = dummynch->GetMean();
        myNchErr[i] = dummynch->GetMeanError();
    }

    //Make a correlation graph
    TGraphErrors *gNchCorr = new TGraphErrors(binnumber, &myNch[0], &offNch[0], &myNchErr[0], &offNchErr[0]);
    gNchCorr->SetTitle("");
    gNchCorr->SetMarkerStyle(20);
    gNchCorr->SetMarkerSize(1.5);
    gNchCorr->SetMarkerColor(kBlack);
    gNchCorr->SetLineColor(kBlack);
    gNchCorr->GetYaxis()->SetTitle("N_{ch}^{official}");
    gNchCorr->GetXaxis()->SetTitle("N_{ch}^{this work}");
    TCanvas* c = new TCanvas();
    TH1F* h = new TH1F("h", "", 100, 0, 27);
    h->SetStats(0);
    h->GetYaxis()->SetTitle("N_{ch}^{official}");
    h->GetXaxis()->SetTitle("N_{ch}^{this work}");
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetRangeUser(0, 27);
    h->Draw();
    gNchCorr->Draw("SAME EP");

    TF1* fitf = new TF1("fitf", "[0]+[1]*x", 0, 100);
    fitf->SetLineColor(kRed);
    fitf->SetLineWidth(2);
    gNchCorr->Fit("fitf");
    fitf->Draw("same");

    gStyle->SetOptFit(111);

    vector<Double_t> percMult_low_5;
    vector<Double_t> percMult_high_5;
    vector<Double_t> percEnergy_low_5;
    vector<Double_t> percEnergy_high_5;

    vector<Double_t> percMult_low_6;
    vector<Double_t> percMult_high_6;
    vector<Double_t> percEnergy_low_6;
    vector<Double_t> percEnergy_high_6;

    vector<Double_t> nptr5;
    vector<Double_t> nptr6;

    Init(5, "V0M", "SPDClusters", percMult_low_5, percMult_high_5, percEnergy_low_5, percEnergy_high_5, nptr5, nptr5);
    Init(6, "V0M", "SPDClusters", percMult_low_6, percMult_high_6, percEnergy_low_6, percEnergy_high_6, nptr6, nptr6);

    Double_t myNch4[100];
    Double_t myNch5[100];

    for (int i = 0; i < percMult_low_5.size(); i++)
    {
        miny = NchHisto->GetYaxis()->FindBin(percMult_low_5[i] + 1e-10);
        maxy = NchHisto->GetYaxis()->FindBin(percMult_high_5[i]);
        minz = NchHisto->GetZaxis()->FindBin(percEnergy_low_5[i] + 1e-10);
        maxz = NchHisto->GetZaxis()->FindBin(percEnergy_high_5[i]);

        TH1D *dummynch = NchHisto->ProjectionX(Form("dummynch%i", i), miny, maxy, minz, maxz);
        myNch4[i] = dummynch->GetMean();
    }
    for (int i = 0; i < percMult_low_6.size(); i++)
    {
        miny = NchHisto->GetYaxis()->FindBin(percMult_low_6[i] + 1e-10);
        maxy = NchHisto->GetYaxis()->FindBin(percMult_high_6[i]);
        minz = NchHisto->GetZaxis()->FindBin(percEnergy_low_6[i] + 1e-10);
        maxz = NchHisto->GetZaxis()->FindBin(percEnergy_high_6[i]);

        TH1D *dummynch = NchHisto->ProjectionX(Form("dummynch%i", i), miny, maxy, minz, maxz);
        myNch5[i] = dummynch->GetMean();
    }
    cout << "class4" << endl;
    cout << "Mine: {";
    for (int i = 0; i < percMult_low_5.size(); i++)
    {
        cout << myNch4[i] << ",";
    }
    cout << "}" << endl;
    cout << "Evaluated: {" ;
    for (int i = 0; i < percMult_low_5.size(); i++)
    {
        cout << fitf->Eval(myNch4[i]) << ",";
    }
    cout << "}" << endl;

    cout << "class5" << endl;
    cout << "Mine: {";
    for (int i = 0; i < percMult_low_6.size(); i++)
    {
        cout << myNch5[i] << ",";
    }
    cout << "}" << endl;
    cout << "Evaluated: {";
    for (int i = 0; i < percMult_low_6.size(); i++)
    {
        cout << fitf->Eval(myNch5[i]) << ",";
    }
    cout << "}" << endl;
}

void Init(Int_t lClassCode,
          TString lWhichEnergyEstimator,
          TString lWhichMultEstimator,
          vector<Double_t> &percentileMult_low,
          vector<Double_t> &percentileMult_high,
          vector<Double_t> &percentileEnergy_low,
          vector<Double_t> &percentileEnergy_high,
          vector<Double_t> &nch,
          vector<Double_t> &nchErr)
{
	// class 0 --> standalone
	Double_t percentileEnergy_low_0[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
	Double_t percentileEnergy_high_0[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
	Double_t percentileMult_low_0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	Double_t percentileMult_high_0[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Long_t n0 = sizeof(percentileMult_low_0) / sizeof(Double_t);
    Double_t nch0[] = {25.75, 19.83, 16.12, 13.76, 12.06, 10.11, 8.07, 6.48, 4.64, 2.52};
    Double_t nchErr0[] = {0.06, 0.06, 0.08, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09};

    // class 1 --> kHighMult
    Double_t percentileMult_low_1[] = {10, 10, 10, 10, 10, 10, 10};
    Double_t percentileMult_high_1[] = {20, 20, 20, 20, 20, 20, 20};
    Double_t percentileEnergy_low_1[] = {0, 5, 10, 20, 30, 40, 50};
    Double_t percentileEnergy_high_1[] = {5, 10, 20, 30, 40, 50, 100};
    Long_t n1 = sizeof(percentileMult_low_1) / sizeof(Double_t);
    Double_t nch1[] = {13.97, 13.79, 13.65, 13.48, 13.35, 13.24, 13.15};
    Double_t nchErr1[] = {0.16, 0.17, 0.17, 0.17, 0.17, 0.17, 0.16};

    // class 2 --> kLowMult
    Double_t percentileMult_low_2[] = {40, 40, 40, 40, 40, 40, 40};
    Double_t percentileMult_high_2[] = {50, 50, 50, 50, 50, 50, 50};
    Double_t percentileEnergy_low_2[] = {0, 20, 30, 40, 50, 60, 70};
    Double_t percentileEnergy_high_2[] = {20, 30, 40, 50, 60, 70, 100};
    Long_t n2 = sizeof(percentileMult_low_2) / sizeof(Double_t);
    Double_t nch2[] = {6.19, 6.15, 6.14, 6.13, 6.09, 6.07, 6.07};
    Double_t nchErr2[] = {0.07, 0.07, 0.07, 0.08, 0.08, 0.09, 0.09};

    Double_t percentileMult_low_3[] = {0, 5, 10, 20, 30, 40, 50};
    Double_t percentileMult_high_3[] = {5, 10, 20, 30, 40, 50, 100};
    Double_t percentileEnergy_low_3[] = {10,10,10,10,10,10,10};
    Double_t percentileEnergy_high_3[] = {20,20,20,20,20,20,20};
    Long_t n3 = sizeof(percentileMult_low_3) / sizeof(Double_t);
    Double_t nch3[] = {22.87, 17.66, 13.65, 10.29, 7.94, 6.19, 4.19};
    Double_t nchErr3[] = {0.3, 0.22, 0.17, 0.12, 0.1, 0.07, 0.05};

    Double_t percentileMult_low_4[] = {10, 20, 30, 40, 50, 60, 70};
    Double_t percentileMult_high_4[] = {20, 30, 40, 50, 60, 70, 100};
    Double_t percentileEnergy_low_4[] = {40, 40, 40, 40, 40, 40, 40};
    Double_t percentileEnergy_high_4[] = {50, 50, 50, 50, 50, 50, 50};
    Long_t n4 = sizeof(percentileMult_low_4) / sizeof(Double_t);
    Double_t nch4[] = {13.25, 10.11, 7.84, 6.13, 4.73, 3.59, 2.19};
    Double_t nchErr4[] = {0.17, 0.13, 0.10, 0.08, 0.06, 0.05, 0.03};

    // class 4 --> fixed low ZN
    Double_t percentileMult_low_5[] = {0, 10, 20, 30, 50};
    Double_t percentileMult_high_5[] = {20, 30, 40, 50, 100};
    Double_t percentileEnergy_low_5[] = {40, 30, 30, 20, 0};
    Double_t percentileEnergy_high_5[] = {60, 70, 50, 50, 30};
    Long_t n5 = sizeof(percentileMult_low_5) / sizeof(Double_t);
    Double_t nch5[] = {};
    Double_t nchErr5[] = {};

    // class 5 --> fixed very low ZN
    Double_t percentileMult_low_6[] = {0, 10, 20, 30};
    Double_t percentileMult_high_6[] = {10, 20, 30, 50};
    Double_t percentileEnergy_low_6[] = {20, 10, 0, 0};
    Double_t percentileEnergy_high_6[] = {30, 30, 20, 10};
    Long_t n6 = sizeof(percentileMult_low_6) / sizeof(Double_t);
    Double_t nch6[] = {};
    Double_t nchErr6[] = {};

    if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 0)
    {
        for (Int_t i = 0; i < n0; i++){
            percentileMult_low.push_back(percentileMult_low_0[i]);
            percentileMult_high.push_back(percentileMult_high_0[i]);
            percentileEnergy_low.push_back(percentileEnergy_low_0[i]);
            percentileEnergy_high.push_back(percentileEnergy_high_0[i]);
            nch.push_back(nch0[i]);
            nchErr.push_back(nchErr0[i]);
        }
    }
    if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 1)
    {
        for (Int_t i = 0; i < n1; i++)
        {
            percentileMult_low.push_back(percentileMult_low_1[i]);
            percentileMult_high.push_back(percentileMult_high_1[i]);
            percentileEnergy_low.push_back(percentileEnergy_low_1[i]);
            percentileEnergy_high.push_back(percentileEnergy_high_1[i]);
            nch.push_back(nch1[i]);
            nchErr.push_back(nchErr1[i]);
        }
    }
    if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 2)
    {
        for (Int_t i = 0; i < n2; i++)
        {
            percentileMult_low.push_back(percentileMult_low_2[i]);
            percentileMult_high.push_back(percentileMult_high_2[i]);
            percentileEnergy_low.push_back(percentileEnergy_low_2[i]);
            percentileEnergy_high.push_back(percentileEnergy_high_2[i]);
            nch.push_back(nch2[i]);
            nchErr.push_back(nchErr2[i]);
        }
    }
    if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 3)
    {
        for (Int_t i = 0; i < n3; i++)
        {
            percentileMult_low.push_back(percentileMult_low_3[i]);
            percentileMult_high.push_back(percentileMult_high_3[i]);
            percentileEnergy_low.push_back(percentileEnergy_low_3[i]);
            percentileEnergy_high.push_back(percentileEnergy_high_3[i]);
            nch.push_back(nch3[i]);
            nchErr.push_back(nchErr3[i]);
        }
    }
    if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 4)
    {
        for (Int_t i = 0; i < n4; i++)
        {
            percentileMult_low.push_back(percentileMult_low_4[i]);
            percentileMult_high.push_back(percentileMult_high_4[i]);
            percentileEnergy_low.push_back(percentileEnergy_low_4[i]);
            percentileEnergy_high.push_back(percentileEnergy_high_4[i]);
            nch.push_back(nch4[i]);
            nchErr.push_back(nchErr4[i]);
        }
    }
    if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 5)
    {
        for (Int_t i = 0; i < n5; i++)
        {
            percentileMult_low.push_back(percentileMult_low_5[i]);
            percentileMult_high.push_back(percentileMult_high_5[i]);
            percentileEnergy_low.push_back(percentileEnergy_low_5[i]);
            percentileEnergy_high.push_back(percentileEnergy_high_5[i]);
            nch.push_back(nch5[i]);
            nchErr.push_back(nchErr5[i]);
        }
    }
    if (lWhichEnergyEstimator.Contains("V0M") && lClassCode == 6)
    {
        for (Int_t i = 0; i < n6; i++)
        {
            percentileMult_low.push_back(percentileMult_low_6[i]);
            percentileMult_high.push_back(percentileMult_high_6[i]);
            percentileEnergy_low.push_back(percentileEnergy_low_6[i]);
            percentileEnergy_high.push_back(percentileEnergy_high_6[i]);
            nch.push_back(nch6[i]);
            nchErr.push_back(nchErr6[i]);
        }
    }
}

Double_t pSPD0100[] = {0,1,5,10,15,20,30,40,50,70,100};
    Long_t nSPD0100 = sizeof(pSPD0100)/sizeof(Double_t) - 1;
    Double_t nchStandalone[] = {25.75, 19.83, 16.12, 13.76, 12.06, 10.11, 8.07, 6.48, 4.64, 2.52};
    Double_t nchSysStandalone[] = {0.06,0.06,0.08,0.09,0.09,0.09,0.09,0.09,0.09,0.09};
