void canvas_hestetics(TCanvas* c);
void graph(TGraphErrors* g, const char* x, const char*y);
void dostudy(TString period = "", TString sel = "SPD", float sel1 = 0., float sel2 = 100., Double_t* perc=0x0 , Long_t nbins=0 );

void SaveToFileNchZDC(){

    TFile *file[5];
    file[0] = TFile::Open("../yields/XiYields_SPDClustersV0M_class0_June23.root");
    file[1] = TFile::Open("../yields/XiYields_SPDClustersV0M_class1_June23.root");
    file[2] = TFile::Open("../yields/XiYields_SPDClustersV0M_class2_June23.root");
    file[3] = TFile::Open("../yields/XiYields_SPDClustersV0M_class4_June23.root");
    file[4] = TFile::Open("../yields/XiYields_SPDClustersV0M_class5_June23.root");

    TGraphAsymmErrors* gNch[5];
    TGraphAsymmErrors* gZDC[5];

    for (int i = 0; i < 5; i++){
        gNch[i] = (TGraphAsymmErrors *)file[i]->Get("NormYieldsNchSyst");
        gZDC[i] = (TGraphAsymmErrors *)file[i]->Get("NormYieldsZDCSumSyst");
    }

    std::vector<Double_t> Nch[5], NchErr[5], ZDC[5], ZDCErr[5];

    for (int i = 0; i < 5; i++){

        for (int j = 0; j < gNch[i]->GetN(); j++)
        {
            Nch[i].push_back(gNch[i]->GetX()[j]);
            NchErr[i].push_back(gNch[i]->GetEXhigh()[j]);
            ZDC[i].push_back(gZDC[i]->GetX()[j]);
            ZDCErr[i].push_back(gZDC[i]->GetEXhigh()[j]);
        }
    }

    // 1-->HighMult
    Double_t perc1_low_1[] = {10, 10, 10, 10, 10, 10, 10};
    Double_t perc1_high_1[] = {20, 20, 20, 20, 20, 20, 20};
    Double_t perc2_low_1[] = {0, 5, 10, 20, 30, 40, 50};
    Double_t perc2_high_1[] = {5, 10, 20, 30, 40, 50, 100};
    const int size1 = sizeof(perc1_low_1) / sizeof(Double_t);

    // 2-->LowMult
    Double_t perc1_low_2[] = {40, 40, 40, 40, 40, 40, 40};
    Double_t perc1_high_2[] = {50, 50, 50, 50, 50, 50, 50};
    Double_t perc2_low_2[] = {0, 20, 30, 40, 50, 60, 70};
    Double_t perc2_high_2[] = {20, 30, 40, 50, 60, 70, 100};
    const int size2 = sizeof(perc1_low_2) / sizeof(Double_t);

    // 3-->HighZN
    Double_t perc1_low_3[] = {0, 10, 20, 30, 50};
    Double_t perc1_high_3[] = {20, 30, 40, 50, 100};
    Double_t perc2_low_3[] = {40, 30, 30, 20, 0};
    Double_t perc2_high_3[] = {60, 70, 50, 50, 30};
    const int size3 = sizeof(perc1_low_3) / sizeof(Double_t);
    // 4-->LowZN
    Double_t perc1_low_4[] = {0, 10, 20, 30};
    Double_t perc1_high_4[] = {10, 20, 30, 50};
    Double_t perc2_low_4[] = {20, 10, 0, 0};
    Double_t perc2_high_4[] = {30, 30, 20, 10};
    const int size4 = sizeof(perc1_low_4) / sizeof(Double_t);

    Double_t percstd1_low[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t percstd1_high[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Double_t percstd2_low[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
    Double_t percstd2_high[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    const int size_std = sizeof(percstd1_low) / sizeof(Double_t);

    // Write numbers
    std::ofstream outfile;
    TString namefile = Form("NchZDCValues.txt");
    TString classi[] = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"};
    outfile.open(namefile, std::ofstream::out | std::ofstream::trunc); // std::ios_base::app
    outfile << "\nStandalone Class:\n ";
    for (int j = 0; j < gNch[0]->GetN(); j++)
    {
        outfile << "\\hline\n"
                << " & " << classi[j].Data() << " & " << percstd2_low[j] << "--" << percstd2_high[j] << "\\% & " << percstd1_low[j] << "--" << percstd1_high[j] << "\\% & "
                << Nch[0].at(j) << "$\\pm$" << NchErr[0].at(j)
                << " & " << ZDC[0].at(j) << "$\\pm$" << ZDCErr[0].at(j) << " \\\\ " << endl;
    }
    outfile << " \n\n" << endl;

    outfile << "\nHigh Multiplicity Class:\n ";
    for (int j = 0; j < gNch[1]->GetN(); j++)
    {
        outfile << "\\hline\n"
                << " & " << classi[j].Data() << " & " << perc2_low_1[j] << "--" << perc2_high_1[j] << "\\% & " << perc1_low_1[j] << "--" << perc1_high_1[j] << "\\% & "
                << Nch[1].at(j) << "$\\pm$" << NchErr[1].at(j)
                << " & " << ZDC[1].at(j) << "$\\pm$" << ZDCErr[1].at(j) << " \\\\ " << endl;
    }
    outfile << " \n\n"
            << endl;

    outfile << "\nLow Multiplicity Class:\n ";
    for (int j = 0; j < gNch[2]->GetN(); j++)
    {
        outfile << "\\hline\n"
                << " & " << classi[j].Data() << " & " << perc2_low_2[j] << "--" << perc2_high_2[j] << "\\% & " << perc1_low_2[j] << "--" << perc1_high_2[j] << "\\% & "
                << Nch[2].at(j) << "$\\pm$" << NchErr[2].at(j)
                << " & " << ZDC[2].at(j) << "$\\pm$" << ZDCErr[2].at(j) << " \\\\ " << endl;
    }
    outfile << " \n\n"
            << endl;

    outfile << "\nHigh ZN Class:\n ";
    for (int j = 0; j < gNch[3]->GetN(); j++)
    {
        outfile << "\\hline\n"
                << " & " << classi[j].Data() << " & " << perc2_low_3[j] << "--" << perc2_high_3[j] << "\\% & " << perc1_low_3[j] << "--" << perc1_high_3[j] << "\\% & "
                << Nch[3].at(j) << "$\\pm$" << NchErr[3].at(j)
                << " & " << ZDC[3].at(j) << "$\\pm$" << ZDCErr[3].at(j) << " \\\\ " << endl;
    }
    outfile << " \n\n"
            << endl;

    outfile << "\nLow ZN Class:\n ";
    for (int j = 0; j < gNch[4]->GetN(); j++)
    {
        outfile << "\\hline\n"
                << " & " << classi[j].Data() << " & " << perc2_low_4[j] << "--" << perc2_high_4[j] << "\\% & " << perc1_low_4[j] << "--" << perc1_high_4[j] << "\\% & "
                << Nch[4].at(j) << "$\\pm$" << NchErr[4].at(j)
                << " & " << ZDC[4].at(j) << "$\\pm$" << ZDCErr[4].at(j) << " \\\\ " << endl;
    }
    outfile << " \n\n"
            << endl;
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