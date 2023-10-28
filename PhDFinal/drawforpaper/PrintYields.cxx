void PrintYields(int nclass = 0){
    // 0 standalone
    // 1 highmult
    // 2 lowmult
    // 3 lowZN
    // 4 highZN
    TString classname[] = {"Standalone", "HighMultiplicity", "LowMultiplicity", "LowZN", "HighZN"};
    TString sel[] = {"_SPDClustersV0M_class0", "_SPDClustersV0M_class1", "_SPDClustersV0M_class5", "_SPDClustersV0M_class2", "_SPDClustersV0M_class4"};
    const int nsel = sizeof(sel) / sizeof(TString);
    Int_t imarker[] = {kFullDiamond, kFullCircle, kFullSquare, kFullCircle, kFullSquare};
    Double_t imarkersize[] = {2.8, 2.3, 2., 2.3, 2.};
    Int_t icolor[] = {kBlack, kRed, kGreen+1, kBlue, kViolet};

    TFile *fXi[nsel], *fLambda[nsel], *fK0s[nsel];
    for (int i = 0; i < nsel; i++){
        fXi[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYields%s_June23.root", "Xi", sel[i].Data()));
        fLambda[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYields%s_June23.root", "Lambda", sel[i].Data()));
        fK0s[i] = TFile::Open(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYields%s_June23.root", "K0Short", sel[i].Data()));
    }

    TGraphErrors *gNumNchStatXi[nsel], *gNumZDCStatXi[nsel];
    TGraphErrors *gNumNchSystXi[nsel], *gNumZDCSystXi[nsel];
    TGraphErrors *gNumNchStatLambda[nsel], *gNumZDCStatLambda[nsel];
    TGraphErrors *gNumNchSystLambda[nsel], *gNumZDCSystLambda[nsel];
    TGraphErrors *gNumNchStatK0s[nsel], *gNumZDCStatK0s[nsel];
    TGraphErrors *gNumNchSystK0s[nsel], *gNumZDCSystK0s[nsel];

    for (int i = 0; i < nsel; i++){
        gNumNchStatXi[i] = (TGraphErrors *)fXi[i]->Get("NormYieldsNormNchStat");
        gNumZDCStatXi[i] = (TGraphErrors *)fXi[i]->Get("NormYieldsNormZDCSumStat");
        gNumNchSystXi[i] = (TGraphErrors *)fXi[i]->Get("NormYieldsNormNchSyst");
        gNumZDCSystXi[i] = (TGraphErrors *)fXi[i]->Get("NormYieldsNormZDCSumSyst");
        gNumNchStatLambda[i] = (TGraphErrors *)fLambda[i]->Get("NormYieldsNormNchStat");
        gNumZDCStatLambda[i] = (TGraphErrors *)fLambda[i]->Get("NormYieldsNormZDCSumStat");
        gNumNchSystLambda[i] = (TGraphErrors *)fLambda[i]->Get("NormYieldsNormNchSyst");
        gNumZDCSystLambda[i] = (TGraphErrors *)fLambda[i]->Get("NormYieldsNormZDCSumSyst");
        gNumNchStatK0s[i] = (TGraphErrors *)fK0s[i]->Get("NormYieldsNormNchStat");
        gNumZDCStatK0s[i] = (TGraphErrors *)fK0s[i]->Get("NormYieldsNormZDCSumStat");
        gNumNchSystK0s[i] = (TGraphErrors *)fK0s[i]->Get("NormYieldsNormNchSyst");
        gNumZDCSystK0s[i] = (TGraphErrors *)fK0s[i]->Get("NormYieldsNormZDCSumSyst");
    }

    std::vector<Double_t> Nch[nsel], NchErr[nsel], NchErrSyst[nsel], ZDC[nsel], ZDCErr[nsel], ZDCErrSyst[nsel];
    std::vector<Double_t> RatioXi[nsel], RatioXiStatErr[nsel], RatioXiSystErr[nsel], RatioLambda[nsel], RatioLambdaStatErr[nsel], RatioLambdaSystErr[nsel], RatioK0s[nsel], RatioK0sStatErr[nsel], RatioK0sSystErr[nsel];
    TGraphErrors *gNchRatioStatXi[nsel], *gZDCRatioStatXi[nsel], *gNchRatioSystXi[nsel], *gZDCRatioSystXi[nsel];
    TGraphErrors *gNchRatioStatLambda[nsel], *gZDCRatioStatLambda[nsel], *gNchRatioSystLambda[nsel], *gZDCRatioSystLambda[nsel];
    TGraphErrors *gNchRatioStatK0s[nsel], *gZDCRatioStatK0s[nsel], *gNchRatioSystK0s[nsel], *gZDCRatioSystK0s[nsel];

    for (int i = 0; i < nsel; i++)
    {
        for (int j = 0; j < gNumNchStatXi[i]->GetN(); j++)
        {
            Nch[i].push_back(gNumNchStatXi[i]->GetX()[j]);
            NchErr[i].push_back(gNumNchStatXi[i]->GetEX()[j]);
            NchErrSyst[i].push_back(gNumNchSystXi[i]->GetEXhigh()[j]);
            ZDC[i].push_back(gNumZDCStatXi[i]->GetX()[j]);
            ZDCErr[i].push_back(gNumZDCStatXi[i]->GetEX()[j]);
            ZDCErrSyst[i].push_back(gNumZDCSystXi[i]->GetEXhigh()[j]);
            RatioXi[i].push_back(gNumNchStatXi[i]->GetY()[j]);
            RatioXiStatErr[i].push_back(gNumNchStatXi[i]->GetEY()[j]);
            RatioXiSystErr[i].push_back(gNumNchSystXi[i]->GetEYhigh()[j]);
            RatioLambda[i].push_back(gNumNchStatLambda[i]->GetY()[j]);
            RatioLambdaStatErr[i].push_back(gNumNchStatLambda[i]->GetEY()[j]);
            RatioLambdaSystErr[i].push_back(gNumNchSystLambda[i]->GetEYhigh()[j]);
            RatioK0s[i].push_back(gNumNchStatK0s[i]->GetY()[j]);
            RatioK0sStatErr[i].push_back(gNumNchStatK0s[i]->GetEY()[j]);
            RatioK0sSystErr[i].push_back(gNumNchSystK0s[i]->GetEYhigh()[j]);
        }
    }

    // Write numbers
    std::ofstream outfile;
    TString namefile = Form("class%i_%s.txt",nclass,classname[nclass].Data());
    TString classi[] = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"};
    outfile.open(namefile, std::ofstream::out | std::ofstream::trunc); // std::ios_base::app
    outfile << Form("\n Class %i: %s\n", nclass, classname[nclass].Data());
    for (int j = 0; j < gNumNchStatXi[nclass]->GetN(); j++)
    {
        outfile << " Selection: " << classi[j].Data() << ":..........................\n" <<
        " Yield K0Short: " << RatioK0s[nclass].at(j) << " +- " << RatioK0sStatErr[nclass].at(j) << " +- " << RatioK0sSystErr[nclass].at(j) << " \n " <<
        "Yield Lambda: " << RatioLambda[nclass].at(j) << " +- " << RatioLambdaStatErr[nclass].at(j) << " +- " << RatioLambdaSystErr[nclass].at(j) << " \n " <<
        "Yield Xi: " << RatioXi[nclass].at(j) << " +- " << RatioXiStatErr[nclass].at(j) << " +- " << RatioXiSystErr[nclass].at(j) << " \n " <<
        "Norm Nch: " << Nch[nclass].at(j) << " +- " << NchErr[nclass].at(j) << " +- " << NchErrSyst[nclass].at(j) << " \n " <<
        "Norm ZDC: " << ZDC[nclass].at(j) << " +- " << ZDCErr[nclass].at(j) << " +- " << ZDCErrSyst[nclass].at(j) << " \n " << endl;
    }
    outfile << " \n\n"  << endl;
}