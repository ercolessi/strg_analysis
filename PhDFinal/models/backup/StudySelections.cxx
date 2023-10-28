TH1D* dopercentile(TH1D* hvar);
Double_t getval(TH1D* h, Double_t perc);
Double_t getperc(TH1D *hcum, Int_t value);

void StudySelections()
{

    //Open file
    TFile *file = TFile::Open(Form("Files/%s.root", "PythiaMonash_Train2627"));
    file->cd("PWGLF_MCPredictionsStrgVsMultVsZDC");
    TList* list  = (TList*)file->FindObjectAny("cList");

    const int nPart = 22;
    //Particles
    TString lPartNames[nPart] = {
        "PiPlus", "PiMinus",
        "KaPlus", "KaMinus",
        "Proton", "AntiProton",
        "K0Short",
        "Lambda", "AntiLambda",
        "XiMinus", "XiPlus",
        "OmegaMinus", "OmegaPlus",
        "Phi",
        "D0", "AntiD0",
        "DPlus", "DMinus",
        "Lambdac", "AntiLambdac",
        "JPsi",
        "Pi0"
    };

    // Definition of histograms
    TH1D *fHistV0MMult;
    TH1D *fHistSPDClusters;
    TH1D *fHistMult08to15;
    TH1D *fHistPt[nPart];
    TH2D *f2DHistPartSPDV0M[nPart];
    TH2D *f2DHistPartRecoPercSPDV0M[nPart];
    TH2D *f2DHistAvPtSPDV0M[nPart];
    TH2D *f2DHistINELgt0SPDV0M;
    TH2D *f2DHistLeadingESPDV0M;
    TH2D *f2DHistEffEnergySPDV0M;
    TH2D *f2DHistNchSPDV0M;
    TH2D *f2DHistNMPISPDV0M;

    // Get histograms
    fHistV0MMult = (TH1D *)list->FindObject("fHistV0MMult");
    fHistSPDClusters = (TH1D *)list->FindObject("fHistSPDClusters");
    fHistMult08to15   = (TH1D*)list->FindObject("fHistMult08to15");
    f2DHistINELgt0SPDV0M = (TH2D *)list->FindObject("f2DHistINELgt0Nch0815V0M");
    f2DHistLeadingESPDV0M = (TH2D *)list->FindObject("f2DHistLeadingENch0815V0M");
    f2DHistEffEnergySPDV0M = (TH2D *)list->FindObject("f2DHistEffEnergyNch0815V0M");
    f2DHistNchSPDV0M = (TH2D *)list->FindObject("f2DHistNchNch0815V0M");
    f2DHistNMPISPDV0M = (TH2D *)list->FindObject("f2DHistNMPINch0815V0M");
    for (Int_t ih = 0; ih < nPart; ih++) {
        fHistPt[ih] = (TH1D *)list->FindObject(Form("fHistPt_%s", lPartNames[ih].Data()));
        f2DHistAvPtSPDV0M[ih] = (TH2D *)list->FindObject(Form("f2DHistAvPtSPDV0M_%s", lPartNames[ih].Data()));
        f2DHistPartSPDV0M[ih] = (TH2D *)list->FindObject(Form("f2DHistPartNch0815V0M_%s", lPartNames[ih].Data()));
    }

    // Percentile calibrations
    TH1D *hcalibV0M = dopercentile(fHistV0MMult);
    hcalibV0M->SetName("hcalibV0M");
    TH1D *hcalibSPDCl = dopercentile(fHistMult08to15);
    hcalibSPDCl->SetName("hcalibNch0815");

    Double_t percentile[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100};
    const int nperc = sizeof(percentile) / sizeof(Double_t);

    TH2D *hnch_ = (TH2D *)f2DHistNchSPDV0M->Clone("hnch_");
    hnch_->Divide(f2DHistINELgt0SPDV0M);

    TH2D* hnch = (TH2D*)f2DHistNchSPDV0M->Clone("hnch");
    hnch->Reset();
    TH2D *hpartk0 = (TH2D *)f2DHistPartSPDV0M[6]->Clone("hpartk0");
    hpartk0->Reset();
    TH2D *hpartL = (TH2D *)f2DHistPartSPDV0M[7]->Clone("hpartL");
    hpartL->Reset();
    TH2D *hzdc = (TH2D *)f2DHistLeadingESPDV0M->Clone("hzdc");
    hzdc->Reset();

    for (int i = 1; i <= hnch_->GetNbinsX(); i++){
        for (int j = 1; j <= hnch_->GetNbinsY(); j++){
            if (hnch_->GetBinContent(i, j) < 7 && hnch_->GetBinContent(i, j) > 6.5) {

                hnch->SetBinContent(i, j, f2DHistNchSPDV0M->GetBinContent(i, j) / f2DHistINELgt0SPDV0M->GetBinContent(i, j));
                hpartk0->SetBinContent(i, j, f2DHistPartSPDV0M[6]->GetBinContent(i, j) / f2DHistINELgt0SPDV0M->GetBinContent(i, j));
                hpartL->SetBinContent(i, j, f2DHistPartSPDV0M[7]->GetBinContent(i, j) / f2DHistINELgt0SPDV0M->GetBinContent(i, j));
                hzdc->SetBinContent(i, j, f2DHistLeadingESPDV0M->GetBinContent(i, j) / f2DHistINELgt0SPDV0M->GetBinContent(i, j));
                if(hzdc->GetBinContent(i, j)!=0) cout << "i = " << i << ", j = " << j << "leading e = " << hzdc->GetBinContent(i, j) << endl;
            }
        }
    }

    TCanvas* c2 = new TCanvas("c2","c2",800,600);
    hnch->Draw("colz");

    TCanvas* c3 = new TCanvas("c3","c3",800,600);
    hpartk0->Draw("colz");

    TCanvas* c4 = new TCanvas("c4","c4",800,600);
    hzdc->Draw("colz");

}

//==========================================================================================
Double_t getval(TH1D* h, Double_t perc){
    //Return the bin of the cumulative once given a percentile

    Double_t diff = 10000.;
    Double_t val = h->GetBinCenter(1);

    for(int i = 2; i <= h->GetNbinsX(); i++){
        if(diff > TMath::Abs(h->GetBinContent(i) - perc/100)) {
            val = h->GetBinCenter(i);
            diff = TMath::Abs(h->GetBinContent(i) - perc/100);
        }
    }

    return val;
}

//==========================================================================================
Double_t getperc(TH1D *hcum, Int_t value)
{
    Double_t percentile = 100 * (hcum->GetBinContent(hcum->FindBin(value)));

    return percentile;
}

//==========================================================================================
TH1D* dopercentile(TH1D* hvar){
    //Creates the cumulative for percentile

    TH1D* hcum = new TH1D(*hvar);
    hcum->SetName("hcum");
    hcum->Reset();

    //Fill Cumulative
    Double_t integral = hvar->Integral(1,hvar->GetNbinsX());
    Double_t val = 0;
    for (Int_t k=1; k<hvar->GetNbinsX(); k++){
        val+= (hvar->GetBinContent(k))/integral;
        hcum->SetBinContent(k,1-val);
    }

    return hcum;
}
