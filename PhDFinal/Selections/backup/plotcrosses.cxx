void init(int lClassCode,
          vector<Double_t> &percentileSPDClusters_low,
          vector<Double_t> &percentileSPDClusters_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          int &nbins);

void plotcrosses()
{
    // Container for ZN and raw multiplicity
    TFile *filecontainer = TFile::Open("NchRawContainer_reference.root");
    if (filecontainer == 0x0)
    {
        cout << "File not found" << endl;
        return;
    }

    TH3D *hnch = (TH3D *)filecontainer->Get("hspd_spdv0m");  // y -> spd, z -> v0m
    TH3D *hZN = (TH3D *)filecontainer->Get("hznsum_spdv0m"); // y -> spd, z -> v0m

    // classes
    enum classname {kStandalone = 0, kHighMult, kLowMult, kHighZN, kLowZN, kVeryLowZN};
    const int nclasses = 6;
    // Percentile
    vector<Double_t> percentileSPDClusters_low[nclasses];
    vector<Double_t> percentileSPDClusters_high[nclasses];
    vector<Double_t> percentileV0M_low[nclasses];
    vector<Double_t> percentileV0M_high[nclasses];
    int nbins[nclasses];
    TString label[] = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"};
    // initiatlize
    init(kStandalone, percentileSPDClusters_low[kStandalone], percentileSPDClusters_high[kStandalone], percentileV0M_low[kStandalone], percentileV0M_high[kStandalone], nbins[kStandalone]);
    cout << " Class code enum: " << kStandalone << endl;
    cout << " Found " << nbins[kStandalone] << " selections:" << endl;
    for (int i = 0; i < nbins[kStandalone]; i++)
    {
        cout << label[i] << " SPDClusters [" << percentileSPDClusters_low[kStandalone][i] << "-" << percentileSPDClusters_high[kStandalone][i] << "] + V0M [" << percentileV0M_low[kStandalone][i] << "-" << percentileV0M_high[kStandalone][i] << "]" << endl;
    }
    cout << "------------------------------------------\n" << endl;

    init(kHighMult, percentileSPDClusters_low[kHighMult], percentileSPDClusters_high[kHighMult], percentileV0M_low[kHighMult], percentileV0M_high[kHighMult], nbins[kHighMult]);
    cout << " Class code enum: " << kHighMult << endl;
    cout << " Found " << nbins[kHighMult] << " selections:" << endl;
    for (int i = 0; i < nbins[kHighMult]; i++)
    {
        cout << label[i] << " SPDClusters [" << percentileSPDClusters_low[kHighMult][i] << "-" << percentileSPDClusters_high[kHighMult][i] << "] + V0M [" << percentileV0M_low[kHighMult][i] << "-" << percentileV0M_high[kHighMult][i] << "]" << endl;
    }
    cout << "------------------------------------------\n" << endl;

    init(kLowMult, percentileSPDClusters_low[kLowMult], percentileSPDClusters_high[kLowMult], percentileV0M_low[kLowMult], percentileV0M_high[kLowMult], nbins[kLowMult]);
    cout << " Class code enum: " << kLowMult << endl;
    cout << " Found " << nbins[kLowMult] << " selections:" << endl;
    for (int i = 0; i < nbins[kLowMult]; i++)
    {
        cout << label[i] << " SPDClusters [" << percentileSPDClusters_low[kLowMult][i] << "-" << percentileSPDClusters_high[kLowMult][i] << "] + V0M [" << percentileV0M_low[kLowMult][i] << "-" << percentileV0M_high[kLowMult][i] << "]" << endl;
    }
    cout << "------------------------------------------\n" << endl;

    init(kHighZN, percentileSPDClusters_low[kHighZN], percentileSPDClusters_high[kHighZN], percentileV0M_low[kHighZN], percentileV0M_high[kHighZN], nbins[kHighZN]);
    cout << " Class code enum: " << kHighZN << endl;
    cout << " Found " << nbins[kHighZN] << " selections:" << endl;
    for (int i = 0; i < nbins[kHighZN]; i++)
    {
        cout << label[i] << " SPDClusters [" << percentileSPDClusters_low[kHighZN][i] << "-" << percentileSPDClusters_high[kHighZN][i] << "] + V0M [" << percentileV0M_low[kHighZN][i] << "-" << percentileV0M_high[kHighZN][i] << "]" << endl;
    }
    cout << "------------------------------------------\n" << endl;

    init(kLowZN, percentileSPDClusters_low[kLowZN], percentileSPDClusters_high[kLowZN], percentileV0M_low[kLowZN], percentileV0M_high[kLowZN], nbins[kLowZN]);
    cout << " Class code enum: " << kLowZN << endl;
    cout << " Found " << nbins[kLowZN] << " selections:" << endl;
    for (int i = 0; i < nbins[kLowZN]; i++)
    {
        cout << label[i] << " SPDClusters [" << percentileSPDClusters_low[kLowZN][i] << "-" << percentileSPDClusters_high[kLowZN][i] << "] + V0M [" << percentileV0M_low[kLowZN][i] << "-" << percentileV0M_high[kLowZN][i] << "]" << endl;
    }
    cout << "------------------------------------------\n" << endl;

    init(kVeryLowZN, percentileSPDClusters_low[kVeryLowZN], percentileSPDClusters_high[kVeryLowZN], percentileV0M_low[kVeryLowZN], percentileV0M_high[kVeryLowZN], nbins[kVeryLowZN]);
    cout << " Class code enum: " << kVeryLowZN << endl;
    cout << " Found " << nbins[kVeryLowZN] << " selections:" << endl;
    for (int i = 0; i < nbins[kVeryLowZN]; i++)
    {
        cout << label[i] << " SPDClusters [" << percentileSPDClusters_low[kVeryLowZN][i] << "-" << percentileSPDClusters_high[kVeryLowZN][i] << "] + V0M [" << percentileV0M_low[kVeryLowZN][i] << "-" << percentileV0M_high[kVeryLowZN][i] << "]" << endl;
    }
    cout << "------------------------------------------\n" << endl;

    TH1D *pnch0[nbins[0]], *pZN0[nbins[0]];
    TH1D *pnch1[nbins[1]], *pZN1[nbins[1]];
    TH1D *pnch2[nbins[2]], *pZN2[nbins[2]];
    TH1D *pnch3[nbins[3]], *pZN3[nbins[3]];
    TH1D *pnch4[nbins[4]], *pZN4[nbins[4]];
    TH1D *pnch5[nbins[5]], *pZN5[nbins[5]];
    double meanZN0[nbins[0]], negmeanZN0[nbins[0]], meanNch0[nbins[0]];
    double meanZN1[nbins[1]], negmeanZN1[nbins[1]], meanNch1[nbins[1]];
    double meanZN2[nbins[2]], negmeanZN2[nbins[2]], meanNch2[nbins[2]];
    double meanZN3[nbins[3]], negmeanZN3[nbins[3]], meanNch3[nbins[3]];
    double meanZN4[nbins[4]], negmeanZN4[nbins[4]], meanNch4[nbins[4]];
    double meanZN5[nbins[5]], negmeanZN5[nbins[5]], meanNch5[nbins[5]];

    //MB
    TH1D* pnchMB = hnch->ProjectionX("pnchMB", hnch->GetYaxis()->FindBin(0. + 1e-10), hnch->GetYaxis()->FindBin(100.), hnch->GetZaxis()->FindBin(0. + 1e-10), hnch->GetZaxis()->FindBin(100.));
    TH1D* pznMB = hZN->ProjectionX("pznMB", hZN->GetYaxis()->FindBin(0. + 1e-10), hZN->GetYaxis()->FindBin(100.), hZN->GetZaxis()->FindBin(0. + 1e-10), hZN->GetZaxis()->FindBin(100.));
    Double_t meanZNMB = pznMB->GetMean();
    Double_t meanNchMB = pnchMB->GetMean();

    double miny, maxy, minz, maxz;
    // kStandalone
    for (int isel = 0; isel < nbins[kStandalone]; isel++){
        miny = hnch->GetYaxis()->FindBin(percentileSPDClusters_low[kStandalone][isel]);
        maxy = hnch->GetYaxis()->FindBin(percentileSPDClusters_high[kStandalone][isel] - 1E-10);
        minz = hnch->GetZaxis()->FindBin(percentileV0M_low[kStandalone][isel]);
        maxz = hnch->GetZaxis()->FindBin(percentileV0M_high[kStandalone][isel] - 1E-10);

        pnch0[isel] = hnch->ProjectionX(Form("pnch0%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kStandalone][isel], percentileSPDClusters_high[kStandalone][isel], "V0M", percentileV0M_low[kStandalone][isel], percentileV0M_high[kStandalone][isel]), miny, maxy, minz, maxz);
        pZN0[isel] = hZN->ProjectionX(Form("pZN0%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kStandalone][isel], percentileSPDClusters_high[kStandalone][isel], "V0M", percentileV0M_low[kStandalone][isel], percentileV0M_high[kStandalone][isel]), miny, maxy, minz, maxz);
        // nch
        meanNch0[isel] = pnch0[isel]->GetMean() / meanNchMB;
        // ZN
        meanZN0[isel] = pZN0[isel]->GetMean() / meanZNMB;
        negmeanZN0[isel] = -meanZN0[isel];
    }
    // kHighMult
    for (int isel = 0; isel < nbins[kHighMult]; isel++){
        miny = hnch->GetYaxis()->FindBin(percentileSPDClusters_low[kHighMult][isel]);
        maxy = hnch->GetYaxis()->FindBin(percentileSPDClusters_high[kHighMult][isel] - 1E-10);
        minz = hnch->GetZaxis()->FindBin(percentileV0M_low[kHighMult][isel]);
        maxz = hnch->GetZaxis()->FindBin(percentileV0M_high[kHighMult][isel] - 1E-10);

        pnch1[isel] = hnch->ProjectionX(Form("pnch1%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kHighMult][isel], percentileSPDClusters_high[kHighMult][isel], "V0M", percentileV0M_low[kHighMult][isel], percentileV0M_high[kHighMult][isel]), miny, maxy, minz, maxz);
        pZN1[isel] = hZN->ProjectionX(Form("pZN1%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kHighMult][isel], percentileSPDClusters_high[kHighMult][isel], "V0M", percentileV0M_low[kHighMult][isel], percentileV0M_high[kHighMult][isel]), miny, maxy, minz, maxz);
        // nch
        meanNch1[isel] = pnch1[isel]->GetMean() / meanNchMB;
        // ZN
        meanZN1[isel] = pZN1[isel]->GetMean() / meanZNMB;
        negmeanZN1[isel] = -meanZN1[isel];
    }
    // kLowMult
    for (int isel = 0; isel < nbins[kLowMult]; isel++){
        miny = hnch->GetYaxis()->FindBin(percentileSPDClusters_low[kLowMult][isel]);
        maxy = hnch->GetYaxis()->FindBin(percentileSPDClusters_high[kLowMult][isel] - 1E-10);
        minz = hnch->GetZaxis()->FindBin(percentileV0M_low[kLowMult][isel]);
        maxz = hnch->GetZaxis()->FindBin(percentileV0M_high[kLowMult][isel] - 1E-10);

        pnch2[isel] = hnch->ProjectionX(Form("pnch2%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kLowMult][isel], percentileSPDClusters_high[kLowMult][isel], "V0M", percentileV0M_low[kLowMult][isel], percentileV0M_high[kLowMult][isel]), miny, maxy, minz, maxz);
        pZN2[isel] = hZN->ProjectionX(Form("pZN2%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kLowMult][isel], percentileSPDClusters_high[kLowMult][isel], "V0M", percentileV0M_low[kLowMult][isel], percentileV0M_high[kLowMult][isel]), miny, maxy, minz, maxz);
        // nch
        meanNch2[isel] = pnch2[isel]->GetMean() / meanNchMB;
        // ZN
        meanZN2[isel] = pZN2[isel]->GetMean() / meanZNMB;
        negmeanZN2[isel] = -meanZN2[isel];
    }
    // kHighZN
    for (int isel = 0; isel < nbins[kHighZN]; isel++){
        miny = hnch->GetYaxis()->FindBin(percentileSPDClusters_low[kHighZN][isel]);
        maxy = hnch->GetYaxis()->FindBin(percentileSPDClusters_high[kHighZN][isel] - 1E-10);
        minz = hnch->GetZaxis()->FindBin(percentileV0M_low[kHighZN][isel]);
        maxz = hnch->GetZaxis()->FindBin(percentileV0M_high[kHighZN][isel] - 1E-10);

        pnch3[isel] = hnch->ProjectionX(Form("pnch3%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kHighZN][isel], percentileSPDClusters_high[kHighZN][isel], "V0M", percentileV0M_low[kHighZN][isel], percentileV0M_high[kHighZN][isel]), miny, maxy, minz, maxz);
        pZN3[isel] = hZN->ProjectionX(Form("pZN3%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kHighZN][isel], percentileSPDClusters_high[kHighZN][isel], "V0M", percentileV0M_low[kHighZN][isel], percentileV0M_high[kHighZN][isel]), miny, maxy, minz, maxz);
        // nch
        meanNch3[isel] = pnch3[isel]->GetMean() / meanNchMB;
        // ZN
        meanZN3[isel] = pZN3[isel]->GetMean() / meanZNMB;
        negmeanZN3[isel] = -meanZN3[isel];
    }
    // kLowZN
    for (int isel = 0; isel < nbins[kLowZN]; isel++){
        miny = hnch->GetYaxis()->FindBin(percentileSPDClusters_low[kLowZN][isel]);
        maxy = hnch->GetYaxis()->FindBin(percentileSPDClusters_high[kLowZN][isel] - 1E-10);
        minz = hnch->GetZaxis()->FindBin(percentileV0M_low[kLowZN][isel]);
        maxz = hnch->GetZaxis()->FindBin(percentileV0M_high[kLowZN][isel] - 1E-10);

        pnch4[isel] = hnch->ProjectionX(Form("pnch4%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kLowZN][isel], percentileSPDClusters_high[kLowZN][isel], "V0M", percentileV0M_low[kLowZN][isel], percentileV0M_high[kLowZN][isel]), miny, maxy, minz, maxz);
        pZN4[isel] = hZN->ProjectionX(Form("pZN4%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kLowZN][isel], percentileSPDClusters_high[kLowZN][isel], "V0M", percentileV0M_low[kLowZN][isel], percentileV0M_high[kLowZN][isel]), miny, maxy, minz, maxz);
        // nch
        meanNch4[isel] = pnch4[isel]->GetMean() / meanNchMB;
        // ZN
        meanZN4[isel] = pZN4[isel]->GetMean() / meanZNMB;
        negmeanZN4[isel] = -meanZN4[isel];
    }
    // kVeryLowZN
    for (int isel = 0; isel < nbins[kVeryLowZN]; isel++)
    {
        miny = hnch->GetYaxis()->FindBin(percentileSPDClusters_low[kVeryLowZN][isel]);
        maxy = hnch->GetYaxis()->FindBin(percentileSPDClusters_high[kVeryLowZN][isel] - 1E-10);
        minz = hnch->GetZaxis()->FindBin(percentileV0M_low[kVeryLowZN][isel]);
        maxz = hnch->GetZaxis()->FindBin(percentileV0M_high[kVeryLowZN][isel] - 1E-10);

        pnch5[isel] = hnch->ProjectionX(Form("pnch5%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kVeryLowZN][isel], percentileSPDClusters_high[kVeryLowZN][isel], "V0M", percentileV0M_low[kVeryLowZN][isel], percentileV0M_high[kVeryLowZN][isel]), miny, maxy, minz, maxz);
        pZN5[isel] = hZN->ProjectionX(Form("pZN5%s%.0f_%.0f_%s%.0f_%.0f", "SPDClusters", percentileSPDClusters_low[kVeryLowZN][isel], percentileSPDClusters_high[kVeryLowZN][isel], "V0M", percentileV0M_low[kVeryLowZN][isel], percentileV0M_high[kVeryLowZN][isel]), miny, maxy, minz, maxz);
        // nch
        meanNch5[isel] = pnch5[isel]->GetMean() / meanNchMB;
        // ZN
        meanZN5[isel] = pZN5[isel]->GetMean() / meanZNMB;
        negmeanZN5[isel] = -meanZN5[isel];
    }

    //Draw stuff
    TH1D *hh = new TH1D("hh", ";#LT n_{ch} #GT w.r.t. to MB (|#eta|<0.5);- #LT ZN #GT w.r.t. MB (|#eta>8|)", 10, 0.3, 3.2);
    hh->SetStats(0);
    hh->GetYaxis()->SetRangeUser(-1.6, -0.2);
    Int_t colors[] = {kRed+2, kOrange+1, kYellow+1, kSpring-1, kGreen+2, kAzure+7, kBlue+2};
    Int_t colors2[] = {kRed + 2, kOrange + 1,kSpring - 1, kAzure + 7, kBlue + 2};

    TGraph *graph0[nbins[kStandalone]];
    for (int isel = 0; isel < nbins[kStandalone]; isel++){
        graph0[isel] = new TGraph(1, &meanNch0[isel], &negmeanZN0[isel]);
        graph0[isel]->SetMarkerStyle(kFullDiamond);
        graph0[isel]->SetMarkerSize(3.8);
        graph0[isel]->SetMarkerColor(kBlack);
        //graph0[isel]->Draw("P SAME");
    }
    TGraph *graph1[nbins[kHighMult]];
    for (int isel = 0; isel < nbins[kHighMult]; isel++){
        graph1[isel] = new TGraph(1, &meanNch1[isel], &negmeanZN1[isel]);
        graph1[isel]->SetMarkerStyle(kOpenCircle);
        graph1[isel]->SetMarkerSize(3.);
        graph1[isel]->SetMarkerColor(colors[isel]);
        //graph1[isel]->Draw("P SAME");
    }
    TGraph *graph2[nbins[kLowMult]];
    for (int isel = 0; isel < nbins[kLowMult]; isel++){
        graph2[isel] = new TGraph(1, &meanNch2[isel], &negmeanZN2[isel]);
        graph2[isel]->SetMarkerStyle(kOpenSquare);
        graph2[isel]->SetMarkerSize(3.);
        graph2[isel]->SetMarkerColor(colors[isel]);
        //graph2[isel]->Draw("P SAME");
    }

    TLegend *legend1 = new TLegend(0.75, 0.18, 0.9, 0.52);
    legend1->SetBorderSize(0);
    legend1->SetFillColor(0);
    legend1->SetTextSize(0.035);
    legend1->SetHeader("class kHighMult:");

    TMarker *m1[nbins[kHighMult]];
    for (int i = 0; i < nbins[kHighMult]; i++) {
        m1[i] = new TMarker(0, 0, kOpenCircle);
        m1[i]->SetMarkerSize(3.0);
        m1[i]->SetMarkerColor(colors[i]);
        legend1->AddEntry(m1[i], label[i], "p");
    }
    //legend1->Draw();

    TLegend *legend2 = new TLegend(0.75, 0.18, 0.9, 0.52);
    legend2->SetBorderSize(0);
    legend2->SetFillColor(0);
    legend2->SetTextSize(0.035);
    legend2->SetHeader("class kLowMult:");

    TMarker *m2[nbins[kLowMult]];
    for (int i = 0; i < nbins[kLowMult]; i++) {
        m2[i] = new TMarker(0, 0, kOpenSquare);
        m2[i]->SetMarkerSize(3.0);
        m2[i]->SetMarkerColor(colors[i]);
        legend2->AddEntry(m2[i], label[i], "p");
    }
    //legend2->Draw();

    for (int isel = 0; isel < nbins[kStandalone]; isel++){
        graph0[isel]->Draw("P SAME");
    }
    TGraph *graph3[nbins[kHighZN]];
    for (int isel = 0; isel < nbins[kHighZN]; isel++){
        graph3[isel] = new TGraph(1, &meanNch3[isel], &negmeanZN3[isel]);
        graph3[isel]->SetMarkerStyle(kOpenSquare);
        graph3[isel]->SetMarkerSize(3.);
        graph3[isel]->SetMarkerColor(colors2[isel]);
        //graph3[isel]->Draw("P SAME");
    }
    TGraph *graph4[nbins[kLowZN]];
    for (int isel = 0; isel < nbins[kLowZN]; isel++){
        graph4[isel] = new TGraph(1, &meanNch4[isel], &negmeanZN4[isel]);
        graph4[isel]->SetMarkerStyle(kFullSquare);
        graph4[isel]->SetMarkerSize(3.);
        graph4[isel]->SetMarkerColor(colors2[isel]);
        //graph4[isel]->Draw("P SAME");
    }
    TGraph *graph5[nbins[kVeryLowZN]];
    for (int isel = 0; isel < nbins[kVeryLowZN]; isel++)
    {
        graph5[isel] = new TGraph(1, &meanNch5[isel], &negmeanZN5[isel]);
        graph5[isel]->SetMarkerStyle(kFullCircle);
        graph5[isel]->SetMarkerSize(3.);
        graph5[isel]->SetMarkerColor(colors2[isel]);
        //graph5[isel]->Draw("P SAME");
    }

    TLegend *legend3 = new TLegend(0.55, 0.18, 0.7, 0.45);
    legend3->SetBorderSize(0);
    legend3->SetFillColor(0);
    legend3->SetTextSize(0.035);
    legend3->SetHeader("class kHighZN:");

    TMarker *m3[nbins[kHighZN]];
    for (int i = 0; i < nbins[kHighZN]; i++) {
        m3[i] = new TMarker(0, 0, kOpenSquare);
        m3[i]->SetMarkerSize(3.);
        m3[i]->SetMarkerColor(colors2[i]);
        legend3->AddEntry(m3[i], label[i], "p");
    }
    //legend3->Draw();

    TLegend *legend4 = new TLegend(0.55, 0.18, 0.7, 0.45);
    legend4->SetBorderSize(0);
    legend4->SetFillColor(0);
    legend4->SetTextSize(0.035);
    legend4->SetHeader("class kHighZN:");

    TMarker *m4[nbins[kLowZN]];
    for (int i = 0; i < nbins[kLowZN]; i++) {
        m4[i] = new TMarker(0, 0, kFullSquare);
        m4[i]->SetMarkerSize(3.);
        m4[i]->SetMarkerColor(colors2[i]);
        legend4->AddEntry(m4[i], label[i], "p");
    }
    //legend4->Draw();

    TLegend *legend5 = new TLegend(0.55, 0.18, 0.7, 0.45);
    legend5->SetBorderSize(0);
    legend5->SetFillColor(0);
    legend5->SetTextSize(0.035);
    legend5->SetHeader("class kLowZN:");

    TMarker *m5[nbins[kVeryLowZN]];
    for (int i = 0; i < nbins[kVeryLowZN]; i++)
    {
        m5[i] = new TMarker(0, 0, kFullCircle);
        m5[i]->SetMarkerSize(3.);
        m5[i]->SetMarkerColor(colors2[i]);
        legend5->AddEntry(m5[i], label[i], "p");
    }
    //legend5->Draw();

    for (int isel = 0; isel < nbins[kStandalone]; isel++){
        graph0[isel]->Draw("P SAME");
    }
    for (int isel = 0; isel < nbins[kHighZN]; isel++){
        graph3[isel]->Draw("P SAME");
    }
    for (int isel = 0; isel < nbins[kLowZN]; isel++){
        graph4[isel]->Draw("P SAME");
    }
    for (int isel = 0; isel < nbins[kHighMult]; isel++){
        graph1[isel]->Draw("P SAME");
    }
    for (int isel = 0; isel < nbins[kLowMult]; isel++){
        graph2[isel]->Draw("P SAME");
    }
    for (int isel = 0; isel < nbins[kVeryLowZN]; isel++)
    {
        graph5[isel]->Draw("P SAME");
    }

    legend1->Draw();
    legend2->Draw();
    legend3->Draw();
    legend4->Draw();
    legend5->Draw();

    //Cross low multiplicity
    TCanvas *canvas4 = new TCanvas("canvas4", "", 1100, 900);
    canvas4->SetBottomMargin(0.12);
    canvas4->SetLeftMargin(0.12);
    canvas4->SetRightMargin(0.02);
    canvas4->SetTopMargin(0.02);
    canvas4->SetTicky();
    canvas4->SetTickx();
    hh->Draw();

    TLegend *legend0 = new TLegend(0.15, 0.85, 0.35, 0.95);
    legend0->SetBorderSize(0);
    legend0->SetFillColor(0);
    legend0->SetTextSize(0.035);
    legend0->AddEntry(graph0[0], "V0M standalone", "p");

    for (int isel = 0; isel < nbins[kStandalone]; isel++)
    {
        graph0[isel]->Draw("P SAME");
    }
    for (int isel = 0; isel < nbins[kLowMult]; isel++)
    {
        graph2[isel]->Draw("P SAME");
    }
    for (int isel = 0; isel < nbins[kLowZN]; isel++)
    {
        graph4[isel]->Draw("P SAME");
    }
    legend2->Draw();
    legend4->Draw();
    legend0->Draw();

    // Cross low multiplicity
    TCanvas *canvas4_bis = new TCanvas("canvas4_bis", "", 1100, 900);
    canvas4_bis->SetBottomMargin(0.12);
    canvas4_bis->SetLeftMargin(0.12);
    canvas4_bis->SetRightMargin(0.02);
    canvas4_bis->SetTopMargin(0.02);
    canvas4_bis->SetTicky();
    canvas4_bis->SetTickx();
    hh->Draw();

    for (int isel = 0; isel < nbins[kStandalone]; isel++)
    {
        graph0[isel]->Draw("P SAME");
    }
    for (int isel = 0; isel < nbins[kLowMult]; isel++)
    {
        graph2[isel]->Draw("P SAME");
    }
    legend2->Draw();
    legend0->Draw();

    canvas4->SaveAs("crosslow.pdf");
    canvas4_bis->SaveAs("crosslow_1.pdf");

    // Cross low multiplicity
    TCanvas *canvas5 = new TCanvas("canvas5", "", 1100, 900);
    canvas5->SetBottomMargin(0.12);
    canvas5->SetLeftMargin(0.12);
    canvas5->SetRightMargin(0.02);
    canvas5->SetTopMargin(0.02);
    canvas5->SetTicky();
    canvas5->SetTickx();
    hh->Draw();


    for (int isel = 0; isel < nbins[kStandalone]; isel++)
    {
        graph0[isel]->Draw("P SAME");
    }
    for (int isel = 0; isel < nbins[kLowMult]; isel++)
    {
        graph1[isel]->Draw("P SAME");
    }
    for (int isel = 0; isel < nbins[kVeryLowZN]; isel++)
    {
        graph5[isel]->Draw("P SAME");
    }
    legend1->Draw();
    legend5->Draw();
    legend0->Draw();

    // Cross low multiplicity
    TCanvas *canvas5_bis = new TCanvas("canvas5_bis", "", 1100, 900);
    canvas5_bis->SetBottomMargin(0.12);
    canvas5_bis->SetLeftMargin(0.12);
    canvas5_bis->SetRightMargin(0.02);
    canvas5_bis->SetTopMargin(0.02);
    canvas5_bis->SetTicky();
    canvas5_bis->SetTickx();
    hh->Draw();

    for (int isel = 0; isel < nbins[kStandalone]; isel++)
    {
        graph0[isel]->Draw("P SAME");
    }
    for (int isel = 0; isel < nbins[kLowMult]; isel++)
    {
        graph1[isel]->Draw("P SAME");
    }
    legend1->Draw();
    legend0->Draw();

    canvas5->SaveAs("crosshigh.pdf");
    canvas5_bis->SaveAs("crosshigh_1.pdf");
}

void init(int lClassCode,
          vector<Double_t> &percentileSPDClusters_low,
          vector<Double_t> &percentileSPDClusters_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          int &nbins){

    // kStandalone
    Double_t percentileV0M_low_0[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
    Double_t percentileV0M_high_0[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    Double_t percentileSPDClusters_low_0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t percentileSPDClusters_high_0[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Long_t n0 = sizeof(percentileV0M_low_0) / sizeof(Double_t);

    // class 8 --> kHighMult
    Double_t percentileSPDClusters_low_1[] = {10, 10, 10, 10, 10, 10, 10};
    Double_t percentileSPDClusters_high_1[] = {20, 20, 20, 20, 20, 20, 20};
    Double_t percentileV0M_low_1[] = {0, 5, 10, 20, 30, 40, 50};
    Double_t percentileV0M_high_1[] = {5, 10, 20, 30, 40, 50, 100};
    Long_t n1 = sizeof(percentileSPDClusters_low_1) / sizeof(Double_t);

    // class 10 --> kLowMult
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

    if (lClassCode == 0)
    {
        nbins = n0;
        for (Int_t i = 0; i < n0; i++)
        {
            percentileSPDClusters_low.push_back(percentileSPDClusters_low_0[i]);
            percentileSPDClusters_high.push_back(percentileSPDClusters_high_0[i]);
            percentileV0M_low.push_back(percentileV0M_low_0[i]);
            percentileV0M_high.push_back(percentileV0M_high_0[i]);
        }
        cout << "\n------------------------------------------" << endl;
        cout << " Initializing kStandalone class..." << endl;
    }

    if (lClassCode == 1)
    {
        nbins = n1;
        for (Int_t i = 0; i < n1; i++)
        {
            percentileSPDClusters_low.push_back(percentileSPDClusters_low_1[i]);
            percentileSPDClusters_high.push_back(percentileSPDClusters_high_1[i]);
            percentileV0M_low.push_back(percentileV0M_low_1[i]);
            percentileV0M_high.push_back(percentileV0M_high_1[i]);
        }
        cout << "\n------------------------------------------" << endl;
        cout << " Initializing kHighMult class..." << endl;
    }

    if (lClassCode == 2)
    {
        nbins = n2;
        for (Int_t i = 0; i < n2; i++)
        {
            percentileSPDClusters_low.push_back(percentileSPDClusters_low_2[i]);
            percentileSPDClusters_high.push_back(percentileSPDClusters_high_2[i]);
            percentileV0M_low.push_back(percentileV0M_low_2[i]);
            percentileV0M_high.push_back(percentileV0M_high_2[i]);
        }
        cout << "\n------------------------------------------" << endl;
        cout << " Initializing kLowMult class..." << endl;
    }

    if (lClassCode == 3)
    {
        nbins = n3;
        for (Int_t i = 0; i < n3; i++)
        {
            percentileSPDClusters_low.push_back(percentileSPDClusters_low_3[i]);
            percentileSPDClusters_high.push_back(percentileSPDClusters_high_3[i]);
            percentileV0M_low.push_back(percentileV0M_low_3[i]);
            percentileV0M_high.push_back(percentileV0M_high_3[i]);
        }
        cout << "\n------------------------------------------" << endl;
        cout << " Initializing kHighZN class..." << endl;
    }

    if (lClassCode == 4)
    {
        nbins = n4;
        for (Int_t i = 0; i < n4; i++)
        {
            percentileSPDClusters_low.push_back(percentileSPDClusters_low_4[i]);
            percentileSPDClusters_high.push_back(percentileSPDClusters_high_4[i]);
            percentileV0M_low.push_back(percentileV0M_low_4[i]);
            percentileV0M_high.push_back(percentileV0M_high_4[i]);
        }
        cout << "\n------------------------------------------" << endl;
        cout << " Initializing kLowZN class..." << endl;
    }

    if (lClassCode == 5)
    {
        nbins = n5;
        for (Int_t i = 0; i < n5; i++)
        {
            percentileSPDClusters_low.push_back(percentileSPDClusters_low_5[i]);
            percentileSPDClusters_high.push_back(percentileSPDClusters_high_5[i]);
            percentileV0M_low.push_back(percentileV0M_low_5[i]);
            percentileV0M_high.push_back(percentileV0M_high_5[i]);
        }
        cout << "\n------------------------------------------" << endl;
        cout << " Initializing kVeryLowZN class..." << endl;
    }
}