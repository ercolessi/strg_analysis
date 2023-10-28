TH1D* dopercentile(TH1D* hvar);
Double_t getval(TH1D* h, Double_t perc);
void process(TString directory, TString outputfile , TString option , TString mcname , Double_t* SPDClperc, Int_t nSPDCl, Double_t* V0Mperc, Int_t nV0M, Bool_t IsFullSim);
void doerrorAB(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr);
void doerrorABC(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr, Double_t C, Double_t Cerr);
void doerror4(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr, Double_t C, Double_t Cerr, Double_t D, Double_t Derr);

void studyselections(TString mcname = "PythiaRopes_Train2628")
{
    //Open file
    TFile* file = TFile::Open(Form("Files/%s.root",mcname.Data()));
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

    cout << "-----------------------------------------------\n" << endl;
    cout << " Processing MC: ....." << mcname.Data() << endl;
    cout << "\n-----------------------------------------------" << endl;

    // Definition of histograms
    TH1D *fHistV0MMult;
    TH1D *fHistMultSPDClusters;
    TH1D *fHistPt[nPart];
    TH2D *f2DHistPartSPDClustersV0M[nPart];
    TH2D *f2DHistPartRecoPercSPDClustersV0M[nPart];
    TH2D *f2DHistAvPtSPDClustersV0M[nPart];
    TH2D *f2DHistINELgt0SPDClustersV0M;
    TH2D *f2DHistLeadingESPDClustersV0M;
    TH2D *f2DHistEffEnergySPDClustersV0M;
    TH2D *f2DHistNchSPDClustersV0M;
    TH2D *f2DHistNMPISPDClustersV0M;

    // Get histograms
    fHistV0MMult = (TH1D *)list->FindObject("fHistV0MMult");
    fHistMultSPDClusters = (TH1D *)list->FindObject("fHistSPDClusters");
    f2DHistINELgt0SPDClustersV0M = (TH2D *)list->FindObject("f2DHistINELgt0SPDV0M");
    f2DHistLeadingESPDClustersV0M = (TH2D *)list->FindObject("f2DHistLeadingESPDV0M");
    f2DHistEffEnergySPDClustersV0M = (TH2D *)list->FindObject("f2DHistEffEnergySPDV0M");
    f2DHistNchSPDClustersV0M = (TH2D *)list->FindObject("f2DHistNchSPDV0M");
    f2DHistNMPISPDClustersV0M = (TH2D *)list->FindObject("f2DHistNMPISPDV0M");
    for (Int_t ih = 0; ih < nPart; ih++) {
        fHistPt[ih] = (TH1D *)list->FindObject(Form("fHistPt_%s", lPartNames[ih].Data()));
        f2DHistPartSPDClustersV0M[ih] = (TH2D *)list->FindObject(Form("f2DHistPartSPDV0M_%s", lPartNames[ih].Data()));
    }

    // Percentile calibrations
    TH1D *hcalibV0M = dopercentile(fHistV0MMult);
    hcalibV0M->SetName("hcalibV0M");
    TH1D *hcalibSPDClusters = dopercentile(fHistMultSPDClusters);
    hcalibSPDClusters->SetName("hcalibSPDClusters");

    // Convert to values V0M and SPDCl

    Double_t percentileV0M[] = {0,10,20,30,40,50,60,70,80,90,100};
    const int npercV0M = sizeof(percentileV0M) / sizeof(Double_t);
    Double_t percentileSPDClusters[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    const int npercSPDClusters = sizeof(percentileSPDClusters) / sizeof(Double_t);

    Double_t V0Mval[npercV0M], SPDClustersval[npercSPDClusters];
    for(int i = 0; i < npercV0M; i++){
        V0Mval[i] = getval(hcalibV0M, percentileV0M[i]);
    }
    for(int i = 0; i < npercSPDClusters; i++){
        SPDClustersval[i] = getval(hcalibSPDClusters, percentileSPDClusters[i]);
    }

    //Variable definition
    Double_t EvtMB, EvtMB_e; // errors must be double for the method TH1::IntegralAndError()
    Double_t NchMB, NchMB_e;
    Double_t EffEnergyMB, EffEnergyMB_e;
    Double_t LeadingEMB, LeadingEMB_e;
    Double_t NMPIMB, NMPIMB_e;
    Double_t YieldMB[nPart], YieldMB_e[nPart];
    Int_t lPAP[9] = {0, 2, 4, 7, 9, 11, 14, 16, 18}; // places in the array which are particle and antiparticle
    Double_t PAPYieldMB[9], PAPYieldMB_e[9];
    Double_t Evt[npercSPDClusters][npercV0M],Evt_e[npercSPDClusters][npercV0M];
    Double_t Nch[npercSPDClusters][npercV0M], Nch_e[npercSPDClusters][npercV0M];
    Double_t NchNorm[npercSPDClusters][npercV0M], NchNorm_e[npercSPDClusters][npercV0M];
    Double_t NMPI[npercSPDClusters][npercV0M], NMPI_e[npercSPDClusters][npercV0M];
    Double_t NMPINorm[npercSPDClusters][npercV0M], NMPINorm_e[npercSPDClusters][npercV0M];
    Double_t EffEnergy[npercSPDClusters][npercV0M], EffEnergy_e[npercSPDClusters][npercV0M];
    Double_t LeadingE[npercSPDClusters][npercV0M], LeadingE_e[npercSPDClusters][npercV0M];
    Double_t LeadingENorm[npercSPDClusters][npercV0M], LeadingENorm_e[npercSPDClusters][npercV0M];
    Double_t Yield[nPart][npercSPDClusters][npercV0M], Yield_e[nPart][npercSPDClusters][npercV0M];
    Double_t YieldNorm[nPart][npercSPDClusters][npercV0M], YieldNorm_e[nPart][npercSPDClusters][npercV0M];
    Double_t Npart[nPart][npercSPDClusters][npercV0M];
    Double_t PAPYieldNorm[9][npercSPDClusters][npercV0M], PAPYieldNorm_e[9][npercSPDClusters][npercV0M];
    Double_t PAPYield[9][npercSPDClusters][npercV0M], PAPYield_e[9][npercSPDClusters][npercV0M];

    TH2D *hSPDClustersV0M_Nch = new TH2D("hSPDClustersV0M_Nch", "hSPDClustersV0M_Nch;SPDClusters;V0M", 10, 0, 100, 10, 0, 100);
    TH2D *hSPDClustersV0M_LeadingE = new TH2D("hSPDClustersV0M_LeadingE", "hSPDClustersV0M_LeadingE;SPDClusters;V0M", 10, 0, 100, 10, 0, 100);

    //-------------------------------------------------------
    //--------------------- INEL>0 (MB)----------------------
    //-------------------------------------------------------

    EvtMB = f2DHistINELgt0SPDClustersV0M -> IntegralAndError(
        1, f2DHistINELgt0SPDClustersV0M->GetNbinsX(),
        1, f2DHistINELgt0SPDClustersV0M->GetNbinsY(),
        EvtMB_e
    );
    //
    NchMB = f2DHistNchSPDClustersV0M -> IntegralAndError(
        1, f2DHistNchSPDClustersV0M->GetNbinsX(),
        1, f2DHistNchSPDClustersV0M->GetNbinsY(),
        NchMB_e
    );
    doerrorAB(NchMB, NchMB_e, EvtMB, EvtMB_e);
    NchMB /= EvtMB;
    //
    EffEnergyMB = f2DHistEffEnergySPDClustersV0M -> IntegralAndError(
        1, f2DHistEffEnergySPDClustersV0M->GetNbinsX(),
        1, f2DHistEffEnergySPDClustersV0M->GetNbinsY(),
        EffEnergyMB_e
        );
    doerrorAB(EffEnergyMB, EffEnergyMB_e, EvtMB, EvtMB_e);
    EffEnergyMB /= EvtMB;
    //
    LeadingEMB = f2DHistLeadingESPDClustersV0M -> IntegralAndError(
        1, f2DHistLeadingESPDClustersV0M->GetNbinsX(),
        1, f2DHistLeadingESPDClustersV0M->GetNbinsY(),
        LeadingEMB_e
        );
    doerrorAB(LeadingEMB, LeadingEMB_e, EvtMB, EvtMB_e);
    LeadingEMB /= EvtMB;
    //
    NMPIMB = f2DHistNMPISPDClustersV0M -> IntegralAndError(
        1, f2DHistNMPISPDClustersV0M->GetNbinsX(),
        1, f2DHistNMPISPDClustersV0M->GetNbinsY(),
        NMPIMB_e
        );
    doerrorAB(NMPIMB, NMPIMB_e, EvtMB, EvtMB_e);
    NMPIMB /= EvtMB;
    //
    for(Int_t ih = 0; ih < nPart; ih++){
        YieldMB[ih]  = f2DHistPartSPDClustersV0M[ih] -> IntegralAndError(
            1, f2DHistPartSPDClustersV0M[ih]->GetNbinsX(),
            1, f2DHistPartSPDClustersV0M[ih]->GetNbinsY(),
            YieldMB_e[ih]
            );
        doerrorAB(YieldMB[ih], YieldMB_e[ih],EvtMB, EvtMB_e);
        YieldMB[ih] /=  EvtMB;
    }
    //
    for (Int_t i = 0; i < 9; i++){
        PAPYieldMB[i] = YieldMB[lPAP[i]] + YieldMB[lPAP[i]+1];
        PAPYieldMB_e[i] = TMath::Sqrt(YieldMB_e[lPAP[i]]*YieldMB_e[lPAP[i]] + YieldMB_e[lPAP[i]+1]*YieldMB_e[lPAP[i]+1]);
    }


    cout << "\nNch(MB) = ......... " << NchMB << " +- " << NchMB_e << endl;
    cout << "EffEnergy(MB) = ... " << EffEnergyMB << " +- " << EffEnergyMB_e << endl;
    cout << "LeadingE(MB) = .... " << LeadingEMB << " +- " << LeadingEMB_e << endl;
    cout << "NMPI(MB) = ........ " << NMPIMB << " +- " << NMPIMB_e << endl;
    cout << "\n-----------------------------------------------" << endl;


    //-------------------------------------------------------
    //--------------------- Selections ----------------------
    //-------------------------------------------------------

    for(int j = 0; j < npercSPDClusters-1; j++){ // SPDClusters loop
        for(int i = 0; i < npercV0M-1; i++){ // V0M loop

            Evt[j][i] = f2DHistINELgt0SPDClustersV0M -> IntegralAndError(
                f2DHistINELgt0SPDClustersV0M->GetXaxis()->FindBin(SPDClustersval[j+1]), f2DHistINELgt0SPDClustersV0M->GetXaxis()->FindBin(SPDClustersval[j]),
                f2DHistINELgt0SPDClustersV0M->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistINELgt0SPDClustersV0M->GetYaxis()->FindBin(V0Mval[i]),
                Evt_e[j][i]
            );
            //
            NchNorm[j][i] = f2DHistNchSPDClustersV0M -> IntegralAndError(
                f2DHistNchSPDClustersV0M->GetXaxis()->FindBin(SPDClustersval[j+1]), f2DHistNchSPDClustersV0M->GetXaxis()->FindBin(SPDClustersval[j]),
                f2DHistNchSPDClustersV0M->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistNchSPDClustersV0M->GetYaxis()->FindBin(V0Mval[i]),
                NchNorm_e[j][i]
            );
            Nch_e[j][i] = NchNorm_e[j][i];
            doerrorAB(NchNorm[j][i], Nch_e[j][i], Evt[j][i], Evt_e[j][i]);
            doerrorABC(NchNorm[j][i], NchNorm_e[j][i], Evt[j][i], Evt_e[j][i], NchMB, NchMB_e);
            NchNorm[j][i] /= Evt[j][i];
            Nch[j][i] = NchNorm[j][i];
            NchNorm[j][i] /= NchMB;
            //
            EffEnergy[j][i] = f2DHistEffEnergySPDClustersV0M -> IntegralAndError(
                f2DHistEffEnergySPDClustersV0M->GetXaxis()->FindBin(SPDClustersval[j+1]), f2DHistEffEnergySPDClustersV0M->GetXaxis()->FindBin(SPDClustersval[j]),
                f2DHistEffEnergySPDClustersV0M->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistEffEnergySPDClustersV0M->GetYaxis()->FindBin(V0Mval[i]),
                EffEnergy_e[j][i]
            );
            doerrorAB(EffEnergy[j][i], EffEnergy_e[j][i], Evt[j][i], Evt_e[j][i]);
            EffEnergy[j][i] /= Evt[j][i];
            //
            LeadingENorm[j][i] = f2DHistLeadingESPDClustersV0M -> IntegralAndError(
                f2DHistLeadingESPDClustersV0M->GetXaxis()->FindBin(SPDClustersval[j+1]), f2DHistLeadingESPDClustersV0M->GetXaxis()->FindBin(SPDClustersval[j]),
                f2DHistLeadingESPDClustersV0M->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistLeadingESPDClustersV0M->GetYaxis()->FindBin(V0Mval[i]),
                LeadingENorm_e[j][i]
            );
            LeadingE_e[j][i] = LeadingENorm_e[j][i];
            doerrorAB(LeadingENorm[j][i], LeadingE_e[j][i], Evt[j][i], Evt_e[j][i]);
            doerrorABC(LeadingENorm[j][i], LeadingENorm_e[j][i], Evt[j][i], Evt_e[j][i], LeadingEMB, LeadingEMB_e);
            LeadingENorm[j][i] /= Evt[j][i];
            LeadingE[j][i] = LeadingENorm[j][i];
            LeadingENorm[j][i] /= LeadingEMB;
            //
            NMPINorm[j][i] = f2DHistNMPISPDClustersV0M -> IntegralAndError(
                f2DHistNMPISPDClustersV0M->GetXaxis()->FindBin(SPDClustersval[j+1]), f2DHistNMPISPDClustersV0M->GetXaxis()->FindBin(SPDClustersval[j]),
                f2DHistNMPISPDClustersV0M->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistNMPISPDClustersV0M->GetYaxis()->FindBin(V0Mval[i]),
                NMPINorm_e[j][i]
            );
            NMPI_e[j][i] = NMPINorm_e[j][i];
            doerrorAB(NMPINorm[j][i], NMPI_e[j][i], Evt[j][i], Evt_e[j][i]);
            doerrorABC(NMPINorm[j][i], NMPINorm_e[j][i], Evt[j][i], Evt_e[j][i], NMPIMB, NMPIMB_e);
            NMPINorm[j][i] /= Evt[j][i];
            NMPI[j][i] = NMPINorm[j][i];
            NMPINorm[j][i] /= NMPIMB;
            //
            for(Int_t ih = 0; ih < nPart; ih++){
                YieldNorm[ih][j][i]  = f2DHistPartSPDClustersV0M[ih] -> IntegralAndError(
                    f2DHistPartSPDClustersV0M[ih]->GetXaxis()->FindBin(SPDClustersval[j+1]), f2DHistPartSPDClustersV0M[ih]->GetXaxis()->FindBin(SPDClustersval[j]),
                    f2DHistPartSPDClustersV0M[ih]->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistPartSPDClustersV0M[ih]->GetYaxis()->FindBin(V0Mval[i]),
                    YieldNorm_e[ih][j][i]
                );
                Yield_e[ih][j][i] = YieldNorm_e[ih][j][i];
                doerrorAB(YieldNorm[ih][j][i], Yield_e[ih][j][i], Evt[j][i], Evt_e[j][i]);
                doerrorABC(YieldNorm[ih][j][i], YieldNorm_e[ih][j][i], Evt[j][i], Evt_e[j][i], YieldMB[ih], YieldMB_e[ih]);
                Npart[ih][j][i] = YieldNorm[ih][j][i];
                YieldNorm[ih][j][i] /= Evt[j][i];
                Yield[ih][j][i] = YieldNorm[ih][j][i];
                YieldNorm[ih][j][i] /= YieldMB[ih];
            }

            for(Int_t ih = 0; ih < 9; ih++){
                PAPYield[ih][j][i] = Yield[lPAP[ih]][j][i] + Yield[lPAP[ih] + 1][j][i];
                PAPYield_e[ih][j][i] = TMath::Sqrt((Yield_e[lPAP[ih]][j][i] * Yield_e[lPAP[ih]][j][i] + Yield_e[lPAP[ih] + 1][j][i] * Yield_e[lPAP[ih] + 1][j][i]));
                PAPYieldNorm[ih][j][i] = PAPYield[ih][j][i];
                PAPYieldNorm_e[ih][j][i] = PAPYield_e[ih][j][i];
                doerrorAB(PAPYieldNorm[ih][j][i], PAPYieldNorm_e[ih][j][i], PAPYieldMB[ih], PAPYieldMB_e[ih]);
                PAPYieldNorm[ih][j][i] /= PAPYieldMB[ih];
            }

            if (LeadingENorm[j][i] < .99 && LeadingENorm[j][i] > 0.80)
            {
                hSPDClustersV0M_Nch->SetBinContent(j+1, i+1, NchNorm[j][i]);
                hSPDClustersV0M_LeadingE->SetBinContent(j + 1, i + 1, LeadingENorm[j][i]);
            }

        } // end V0M loop
    } // end SPDClusters Cluster loop

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    hSPDClustersV0M_Nch->SetStats(0);
    hSPDClustersV0M_Nch->Draw("colz, text");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    hSPDClustersV0M_LeadingE->SetStats(0);
    hSPDClustersV0M_LeadingE->Draw("colz, text");

    //Monash tuned classes
    //Double_t percentileSPDClusters_low[] = {0, 10, 20, 30, 40, 80};
    //Double_t percentileSPDClusters_high[] = {10, 20, 30, 40, 80, 90};
    //Double_t percentileV0M_low[] = {40, 30, 20, 10, 10, 0};
    //Double_t percentileV0M_high[] = {50, 40, 40, 30, 20, 10};

    Double_t percentileSPDClusters_low[] = {};
    Double_t percentileSPDClusters_high[] = {};
    Double_t percentileV0M_low[] = {};
    Double_t percentileV0M_high[] = {};
    Long_t n = sizeof(percentileV0M_low) / sizeof(Double_t);


    Double_t V0Mval_low[n], SPDClustersval_low[n];
    Double_t V0Mval_high[n], SPDClustersval_high[n];
    for(int i = 0; i < n; i++){
        V0Mval_low[i] = getval(hcalibV0M, percentileV0M_low[i]);
        V0Mval_high[i] = getval(hcalibV0M, percentileV0M_high[i]);
        SPDClustersval_low[i] = getval(hcalibSPDClusters, percentileSPDClusters_low[i]);
        SPDClustersval_high[i] = getval(hcalibSPDClusters, percentileSPDClusters_high[i]);
    }

    Double_t Evt_[n], Evt_e_[n];
    Double_t LeadingENorm_[n], LeadingENorm_e_[n], LeadingE_[n], LeadingE_e_[n];

    for (int k = 0; k < n; k++)
    { // SPD Cluster loop

        Evt_[k] = f2DHistINELgt0SPDClustersV0M->IntegralAndError(
            f2DHistINELgt0SPDClustersV0M->GetXaxis()->FindBin(SPDClustersval_high[k]), f2DHistINELgt0SPDClustersV0M->GetXaxis()->FindBin(SPDClustersval_low[k]),
            f2DHistINELgt0SPDClustersV0M->GetYaxis()->FindBin(V0Mval_high[k]), f2DHistINELgt0SPDClustersV0M->GetYaxis()->FindBin(V0Mval_low[k]),
            Evt_e_[k]);
        //
        LeadingENorm_[k] = f2DHistLeadingESPDClustersV0M->IntegralAndError(
            f2DHistLeadingESPDClustersV0M->GetXaxis()->FindBin(SPDClustersval_high[k]), f2DHistLeadingESPDClustersV0M->GetXaxis()->FindBin(SPDClustersval_low[k]),
            f2DHistLeadingESPDClustersV0M->GetYaxis()->FindBin(V0Mval_high[k]), f2DHistLeadingESPDClustersV0M->GetYaxis()->FindBin(V0Mval_low[k]),
            LeadingENorm_e_[k]);
        LeadingE_e_[k] = LeadingENorm_e_[k];
        doerrorAB(LeadingENorm_[k], LeadingE_e_[k], Evt_[k], Evt_e_[k]);
        doerrorABC(LeadingENorm_[k], LeadingENorm_e_[k], Evt_[k], Evt_e_[k], LeadingEMB, LeadingEMB_e);
        LeadingENorm_[k] /= Evt_[k];
        LeadingE_[k] = LeadingENorm_[k];
        LeadingENorm_[k] /= LeadingEMB;

        cout << "SPDcl " << percentileSPDClusters_low[k] << "-" << percentileSPDClusters_high[k] << " V0M :" << percentileV0M_low[k] << "-" << percentileV0M_high[k] << " --> LE " << LeadingENorm_[k] << endl;
    }

        // Yields
        /*TGraphErrors *gYield_VsNchNorm[nPart];
        TGraphErrors *gYield_VsNch[nPart];
        TGraphErrors *gYield_VsPerc[nPart];
        TGraphErrors *gYield_VsLeadENorm[nPart];
        TGraphErrors *gYield_VsLeadE[nPart];
        TGraphErrors *gYield_VsEE[nPart];
        TGraphErrors *gYield_VsNMPINorm[nPart];
        TGraphErrors *gYield_VsNMPI[nPart];

        // Yields norm
        TGraphErrors *gYieldRatioToMB_VsNchNorm[nPart];
        TGraphErrors *gYieldRatioToMB_VsNch[nPart];
        TGraphErrors *gYieldRatioToMB_VsPerc[nPart];
        TGraphErrors *gYieldRatioToMB_VsLeadENorm[nPart];
        TGraphErrors *gYieldRatioToMB_VsLeadE[nPart];
        TGraphErrors *gYieldRatioToMB_VsEE[nPart];
        TGraphErrors *gYieldRatioToMB_VsNMPINorm[nPart];
        TGraphErrors *gYieldRatioToMB_VsNMPI[nPart];

        // Particle + antiparticle
        TGraphErrors *gPAPYield_VsNchNorm[9];
        TGraphErrors *gPAPYield_VsNch[9];
        TGraphErrors *gPAPYield_VsPerc[9];
        TGraphErrors *gPAPYield_VsLeadENorm[9];
        TGraphErrors *gPAPYield_VsLeadE[9];
        TGraphErrors *gPAPYield_VsEE[9];
        TGraphErrors *gPAPYield_VsNMPINorm[9];
        TGraphErrors *gPAPYield_VsNMPI[9];
        TGraphErrors *gPAPYieldRatioToMB_VsNchNorm[9];
        TGraphErrors *gPAPYieldRatioToMB_VsNch[9];
        TGraphErrors *gPAPYieldRatioToMB_VsPerc[9];
        TGraphErrors *gPAPYieldRatioToMB_VsLeadENorm[9];
        TGraphErrors *gPAPYieldRatioToMB_VsLeadE[9];
        TGraphErrors *gPAPYieldRatioToMB_VsEE[9];
        TGraphErrors *gPAPYieldRatioToMB_VsNMPINorm[9];
        TGraphErrors *gPAPYieldRatioToMB_VsNMPI[9];

        // Mean pT
        TGraphErrors *gAvPt_VsNch[nPart];
        TGraphErrors *gAvPt_VsNchNorm[nPart];
        TGraphErrors *gAvPt_VsLeadE[nPart];
        TGraphErrors *gAvPt_VsLeadENorm[nPart];
        TGraphErrors *gAvPt_VsEE[nPart];
        TGraphErrors *gAvPt_VsNMPI[nPart];
        TGraphErrors *gAvPt_VsNMPINorm[nPart];

        TGraphErrors *gNchVsNMPI;
        TGraphErrors *gLEVsNMPI;
        TGraphErrors *gEEVsNMPI;
        if (isPythia) {
            gNchVsNMPI = new TGraphErrors(nDiff, Nch, NMPI, Nch_e, NMPI_e);
            gNchVsNMPI->GetYaxis()->SetTitle("NMPI");
            gNchVsNMPI->GetXaxis()->SetTitle("Nch");
            gNchVsNMPI->SetTitle("");
            gNchVsNMPI->SetMarkerStyle(kFullCircle);
            gNchVsNMPI->SetName("NchVsNMPI");
            //
            gLEVsNMPI = new TGraphErrors(nDiff, LeadingE, NMPI, LeadingE_e, NMPI_e);
            gLEVsNMPI->GetYaxis()->SetTitle("NMPI");
            gLEVsNMPI->GetXaxis()->SetTitle("Leading Energy");
            gLEVsNMPI->SetTitle("");
            gLEVsNMPI->SetMarkerStyle(kFullCircle);
            gLEVsNMPI->SetName("LEVsNMPI");
            //
            gEEVsNMPI = new TGraphErrors(nDiff, EffEnergy, NMPI, EffEnergy_e, NMPI_e);
            gEEVsNMPI->GetYaxis()->SetTitle("NMPI");
            gEEVsNMPI->GetXaxis()->SetTitle("Effective Energy");
            gEEVsNMPI->SetTitle("");
            gEEVsNMPI->SetMarkerStyle(kFullCircle);
            gEEVsNMPI->SetName("EEVsNMPI");
        }
        TGraphErrors *gLEVsNch = new TGraphErrors(nDiff, Nch, LeadingE, Nch_e, LeadingE_e);
        gLEVsNch->GetYaxis()->SetTitle("Leading Energy");
        gLEVsNch->GetXaxis()->SetTitle("Nch");
        gLEVsNch->SetTitle("");
        gLEVsNch->SetMarkerStyle(kFullCircle);
        gLEVsNch->SetName("LEVsNch");
        TGraphErrors *gNormLEVsNch = new TGraphErrors(nDiff, NchNorm, LeadingENorm, NchNorm_e, LeadingENorm_e);
        gNormLEVsNch->GetYaxis()->SetTitle("E_{leading}/E_{leading}^{MB}");
        gNormLEVsNch->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
        gNormLEVsNch->SetTitle("");
        gNormLEVsNch->SetMarkerStyle(kFullCircle);
        gNormLEVsNch->SetName("NormLEVsNch");

        for(Int_t ih = 0; ih < nPart; ih++){ //Single particle

            // Vs Nch
            gYieldRatioToMB_VsNchNorm[ih] = new TGraphErrors(nDiff,NchNorm,YieldNorm[ih],NchNorm_e,YieldNorm_e[ih]);
            gYieldRatioToMB_VsNchNorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{#LT %s #GT_{MB}}",lPartNames[ih].Data(),lPartNames[ih].Data()));
            gYieldRatioToMB_VsNchNorm[ih]->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
            gYieldRatioToMB_VsNchNorm[ih]->SetTitle("");
            gYieldRatioToMB_VsNchNorm[ih]->SetMarkerStyle(kFullCircle);
            gYieldRatioToMB_VsNchNorm[ih]->SetName(Form("YieldRatioToMB_VsNchNorm_%s",lPartNames[ih].Data()));
            //
            gYieldRatioToMB_VsNch[ih] = new TGraphErrors(nDiff,Nch,YieldNorm[ih],Nch_e,YieldNorm_e[ih]);
            gYieldRatioToMB_VsNch[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{#LT %s #GT_{MB}}",lPartNames[ih].Data(),lPartNames[ih].Data()));
            gYieldRatioToMB_VsNch[ih]->GetXaxis()->SetTitle("n_{ch}");
            gYieldRatioToMB_VsNch[ih]->SetTitle("");
            gYieldRatioToMB_VsNch[ih]->SetMarkerStyle(kFullCircle);
            gYieldRatioToMB_VsNch[ih]->SetName(Form("YieldRatioToMB_VsNch_%s",lPartNames[ih].Data()));

            // Vs Lead Energy
            gYieldRatioToMB_VsLeadENorm[ih] = new TGraphErrors(nDiff,LeadingENorm,YieldNorm[ih],LeadingENorm_e,YieldNorm_e[ih]);
            gYieldRatioToMB_VsLeadENorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{#LT %s #GT_{MB}}",lPartNames[ih].Data(),lPartNames[ih].Data()));
            gYieldRatioToMB_VsLeadENorm[ih]->GetXaxis()->SetTitle("E_{leading}/E_{leading}^{MB}");
            gYieldRatioToMB_VsLeadENorm[ih]->SetTitle("");
            gYieldRatioToMB_VsLeadENorm[ih]->SetMarkerStyle(kFullCircle);
            gYieldRatioToMB_VsLeadENorm[ih]->SetName(Form("YieldRatioToMB_VsLeadENorm_%s",lPartNames[ih].Data()));
            //
            gYieldRatioToMB_VsLeadE[ih] = new TGraphErrors(nDiff,LeadingE,YieldNorm[ih],LeadingE_e,YieldNorm_e[ih]);
            gYieldRatioToMB_VsLeadE[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{#LT %s #GT_{MB}}",lPartNames[ih].Data(),lPartNames[ih].Data()));
            gYieldRatioToMB_VsLeadE[ih]->GetXaxis()->SetTitle("E_{leading}");
            gYieldRatioToMB_VsLeadE[ih]->SetTitle("");
            gYieldRatioToMB_VsLeadE[ih]->SetMarkerStyle(kFullCircle);
            gYieldRatioToMB_VsLeadE[ih]->SetName(Form("YieldRatioToMB_VsLeadE_%s",lPartNames[ih].Data()));

            // Vs Eff Energy
            gYieldRatioToMB_VsEE[ih] = new TGraphErrors(nDiff,EffEnergy,YieldNorm[ih],EffEnergy_e,YieldNorm_e[ih]);
            gYieldRatioToMB_VsEE[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{#LT %s #GT_{MB}}",lPartNames[ih].Data(),lPartNames[ih].Data()));
            gYieldRatioToMB_VsEE[ih]->GetXaxis()->SetTitle("E_{effective}");
            gYieldRatioToMB_VsEE[ih]->SetTitle("");
            gYieldRatioToMB_VsEE[ih]->SetMarkerStyle(kFullCircle);
            gYieldRatioToMB_VsEE[ih]->SetName(Form("YieldRatioToMB_VsEE_%s",lPartNames[ih].Data()));

            // Vs Nch
            gYield_VsNchNorm[ih] = new TGraphErrors(nDiff, NchNorm, Yield[ih], NchNorm_e, Yield_e[ih]);
            gYield_VsNchNorm[ih]->GetYaxis()->SetTitle(Form("%s", lPartNames[ih].Data()));
            gYield_VsNchNorm[ih]->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
            gYield_VsNchNorm[ih]->SetTitle("");
            gYield_VsNchNorm[ih]->SetMarkerStyle(kFullCircle);
            gYield_VsNchNorm[ih]->SetName(Form("Yield_VsNchNorm_%s", lPartNames[ih].Data()));
            //
            gYield_VsNch[ih] = new TGraphErrors(nDiff, Nch, Yield[ih], Nch_e, Yield_e[ih]);
            gYield_VsNch[ih]->GetYaxis()->SetTitle(Form("%s", lPartNames[ih].Data()));
            gYield_VsNch[ih]->GetXaxis()->SetTitle("n_{ch}");
            gYield_VsNch[ih]->SetTitle("");
            gYield_VsNch[ih]->SetMarkerStyle(kFullCircle);
            gYield_VsNch[ih]->SetName(Form("Yield_VsNch_%s", lPartNames[ih].Data()));

            // Vs Lead Energy
            gYield_VsLeadENorm[ih] = new TGraphErrors(nDiff, LeadingENorm, Yield[ih], LeadingENorm_e, Yield_e[ih]);
            gYield_VsLeadENorm[ih]->GetYaxis()->SetTitle(Form("%s", lPartNames[ih].Data()));
            gYield_VsLeadENorm[ih]->GetXaxis()->SetTitle("E_{leading}/E_{leading}^{MB}");
            gYield_VsLeadENorm[ih]->SetTitle("");
            gYield_VsLeadENorm[ih]->SetMarkerStyle(kFullCircle);
            gYield_VsLeadENorm[ih]->SetName(Form("Yield_VsLeadENorm_%s", lPartNames[ih].Data()));
            //
            gYield_VsLeadE[ih] = new TGraphErrors(nDiff, LeadingE, Yield[ih], LeadingE_e, Yield_e[ih]);
            gYield_VsLeadE[ih]->GetYaxis()->SetTitle(Form("%s", lPartNames[ih].Data()));
            gYield_VsLeadE[ih]->GetXaxis()->SetTitle("E_{leading}");
            gYield_VsLeadE[ih]->SetTitle("");
            gYield_VsLeadE[ih]->SetMarkerStyle(kFullCircle);
            gYield_VsLeadE[ih]->SetName(Form("Yield_VsLeadE_%s", lPartNames[ih].Data()));

            // Vs Eff Energy
            gYield_VsEE[ih] = new TGraphErrors(nDiff, EffEnergy, Yield[ih], EffEnergy_e, Yield_e[ih]);
            gYield_VsEE[ih]->GetYaxis()->SetTitle(Form("%s", lPartNames[ih].Data()));
            gYield_VsEE[ih]->GetXaxis()->SetTitle("E_{effective}");
            gYield_VsEE[ih]->SetTitle("");
            gYield_VsEE[ih]->SetMarkerStyle(kFullCircle);
            gYield_VsEE[ih]->SetName(Form("Yield_VsEE_%s", lPartNames[ih].Data()));

            //Mean Pt
            gAvPt_VsNchNorm[ih] = new TGraphErrors(nDiff,NchNorm,AvPt[ih],NchNorm_e,AvPt_e[ih]);
            gAvPt_VsNchNorm[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
            gAvPt_VsNchNorm[ih]->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
            gAvPt_VsNchNorm[ih]->SetTitle("");
            gAvPt_VsNchNorm[ih]->SetMarkerStyle(kFullCircle);
            gAvPt_VsNchNorm[ih]->SetName(Form("AvPt_VsNchNorm_%s",lPartNames[ih].Data()));
            //
            gAvPt_VsNch[ih] = new TGraphErrors(nDiff,Nch,AvPt[ih],Nch_e,AvPt_e[ih]);
            gAvPt_VsNch[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
            gAvPt_VsNch[ih]->GetXaxis()->SetTitle("n_{ch}");
            gAvPt_VsNch[ih]->SetTitle("");
            gAvPt_VsNch[ih]->SetMarkerStyle(kFullCircle);
            gAvPt_VsNch[ih]->SetName(Form("AvPt_VsNch_%s",lPartNames[ih].Data()));
            //
            gAvPt_VsLeadENorm[ih] = new TGraphErrors(nDiff,LeadingENorm,AvPt[ih],LeadingENorm_e,AvPt_e[ih]);
            gAvPt_VsLeadENorm[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
            gAvPt_VsLeadENorm[ih]->GetXaxis()->SetTitle("E_{leading}/E_{leading}^{MB}");
            gAvPt_VsLeadENorm[ih]->SetTitle("");
            gAvPt_VsLeadENorm[ih]->SetMarkerStyle(kFullCircle);
            gAvPt_VsLeadENorm[ih]->SetName(Form("AvPt_VsLeadENorm_%s",lPartNames[ih].Data()));
            //
            gAvPt_VsLeadE[ih] = new TGraphErrors(nDiff,LeadingE,AvPt[ih],LeadingE_e,AvPt_e[ih]);
            gAvPt_VsLeadE[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
            gAvPt_VsLeadE[ih]->GetXaxis()->SetTitle("E_{leading}");
            gAvPt_VsLeadE[ih]->SetTitle("");
            gAvPt_VsLeadE[ih]->SetMarkerStyle(kFullCircle);
            gAvPt_VsLeadE[ih]->SetName(Form("AvPt_VsLeadE_%s",lPartNames[ih].Data()));
            //
            gAvPt_VsEE[ih] = new TGraphErrors(nDiff,EffEnergy,AvPt[ih],EffEnergy_e,AvPt_e[ih]);
            gAvPt_VsEE[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
            gAvPt_VsEE[ih]->GetXaxis()->SetTitle("E_{effective}");
            gAvPt_VsEE[ih]->SetTitle("");
            gAvPt_VsEE[ih]->SetMarkerStyle(kFullCircle);
            gAvPt_VsEE[ih]->SetName(Form("AvPt_VsEE_%s",lPartNames[ih].Data()));

            // Vs NMPI
            if (isPythia) {
                gYieldRatioToMB_VsNMPINorm[ih] = new TGraphErrors(nDiff,NMPINorm,YieldNorm[ih],NMPINorm_e,YieldNorm_e[ih]);
                gYieldRatioToMB_VsNMPINorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{#LT %s #GT_{MB}}",lPartNames[ih].Data(),lPartNames[ih].Data()));
                gYieldRatioToMB_VsNMPINorm[ih]->GetXaxis()->SetTitle("NMPI/NMPI^{MB}");
                gYieldRatioToMB_VsNMPINorm[ih]->SetTitle("");
                gYieldRatioToMB_VsNMPINorm[ih]->SetMarkerStyle(kFullCircle);
                gYieldRatioToMB_VsNMPINorm[ih]->SetName(Form("YieldRatioToMB_VsNMPINorm_%s",lPartNames[ih].Data()));
                //
                gYieldRatioToMB_VsNMPI[ih] = new TGraphErrors(nDiff,NMPI,YieldNorm[ih],NMPI_e,YieldNorm_e[ih]);
                gYieldRatioToMB_VsNMPI[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{#LT %s #GT_{MB}}",lPartNames[ih].Data(),lPartNames[ih].Data()));
                gYieldRatioToMB_VsNMPI[ih]->GetXaxis()->SetTitle("NMPI");
                gYieldRatioToMB_VsNMPI[ih]->SetTitle("");
                gYieldRatioToMB_VsNMPI[ih]->SetMarkerStyle(kFullCircle);
                gYieldRatioToMB_VsNMPI[ih]->SetName(Form("YieldRatioToMB_VsNMPI_%s",lPartNames[ih].Data()));
                //
                gYield_VsNMPINorm[ih] = new TGraphErrors(nDiff, NMPINorm, Yield[ih], NMPINorm_e, Yield_e[ih]);
                gYield_VsNMPINorm[ih]->GetYaxis()->SetTitle(Form("%s", lPartNames[ih].Data()));
                gYield_VsNMPINorm[ih]->GetXaxis()->SetTitle("NMPI/NMPI^{MB}");
                gYield_VsNMPINorm[ih]->SetTitle("");
                gYield_VsNMPINorm[ih]->SetMarkerStyle(kFullCircle);
                gYield_VsNMPINorm[ih]->SetName(Form("Yield_VsNMPINorm_%s", lPartNames[ih].Data()));
                //
                gYield_VsNMPI[ih] = new TGraphErrors(nDiff, NMPI, Yield[ih], NMPI_e, Yield_e[ih]);
                gYield_VsNMPI[ih]->GetYaxis()->SetTitle(Form("%s", lPartNames[ih].Data()));
                gYield_VsNMPI[ih]->GetXaxis()->SetTitle("NMPI");
                gYield_VsNMPI[ih]->SetTitle("");
                gYield_VsNMPI[ih]->SetMarkerStyle(kFullCircle);
                gYield_VsNMPI[ih]->SetName(Form("Yield_VsNMPI_%s", lPartNames[ih].Data()));
                //
                gAvPt_VsNMPINorm[ih] = new TGraphErrors(nDiff,NMPINorm,AvPt[ih],NMPINorm_e,AvPt_e[ih]);
                gAvPt_VsNMPINorm[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
                gAvPt_VsNMPINorm[ih]->GetXaxis()->SetTitle("NMPI/NMPI^{MB}");
                gAvPt_VsNMPINorm[ih]->SetTitle("");
                gAvPt_VsNMPINorm[ih]->SetMarkerStyle(kFullCircle);
                gAvPt_VsNMPINorm[ih]->SetName(Form("AvPt_VsNMPINorm_%s",lPartNames[ih].Data()));
                //
                gAvPt_VsNMPI[ih] = new TGraphErrors(nDiff,NMPI,AvPt[ih],NMPI_e,AvPt_e[ih]);
                gAvPt_VsNMPI[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
                gAvPt_VsNMPI[ih]->GetXaxis()->SetTitle("NMPI");
                gAvPt_VsNMPI[ih]->SetTitle("");
                gAvPt_VsNMPI[ih]->SetMarkerStyle(kFullCircle);
                gAvPt_VsNMPI[ih]->SetName(Form("AvPt_VsNMPI_%s",lPartNames[ih].Data()));
            }
        }

        for (int ih = 0; ih < 9; ih++){
            gPAPYieldRatioToMB_VsNchNorm[ih] = new TGraphErrors(nDiff, NchNorm, PAPYieldNorm[ih], NchNorm_e, PAPYieldNorm_e[ih]);
            gPAPYieldRatioToMB_VsNchNorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s%s}{#LT %s%s #GT_{MB}}", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data(), lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            gPAPYieldRatioToMB_VsNchNorm[ih]->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
            gPAPYieldRatioToMB_VsNchNorm[ih]->SetTitle("");
            gPAPYieldRatioToMB_VsNchNorm[ih]->SetMarkerStyle(kFullCircle);
            gPAPYieldRatioToMB_VsNchNorm[ih]->SetName(Form("YieldRatioToMB_VsNchNorm_%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            //
            gPAPYieldRatioToMB_VsNch[ih] = new TGraphErrors(nDiff, Nch, PAPYieldNorm[ih], Nch_e, PAPYieldNorm_e[ih]);
            gPAPYieldRatioToMB_VsNch[ih]->GetYaxis()->SetTitle(Form("#frac{%s%s}{#LT %s%s #GT_{MB}}", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data(), lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            gPAPYieldRatioToMB_VsNch[ih]->GetXaxis()->SetTitle("n_{ch}");
            gPAPYieldRatioToMB_VsNch[ih]->SetTitle("");
            gPAPYieldRatioToMB_VsNch[ih]->SetMarkerStyle(kFullCircle);
            gPAPYieldRatioToMB_VsNch[ih]->SetName(Form("YieldRatioToMB_VsNch_%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));

            // Vs Lead Energy
            gPAPYieldRatioToMB_VsLeadENorm[ih] = new TGraphErrors(nDiff, LeadingENorm, PAPYieldNorm[ih], LeadingENorm_e, PAPYieldNorm_e[ih]);
            gPAPYieldRatioToMB_VsLeadENorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s%s}{#LT %s%s #GT_{MB}}", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data(), lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            gPAPYieldRatioToMB_VsLeadENorm[ih]->GetXaxis()->SetTitle("E_{leading}/E_{leading}^{MB}");
            gPAPYieldRatioToMB_VsLeadENorm[ih]->SetTitle("");
            gPAPYieldRatioToMB_VsLeadENorm[ih]->SetMarkerStyle(kFullCircle);
            gPAPYieldRatioToMB_VsLeadENorm[ih]->SetName(Form("YieldRatioToMB_VsLeadENorm_%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            //
            gPAPYieldRatioToMB_VsLeadE[ih] = new TGraphErrors(nDiff, LeadingE, PAPYieldNorm[ih], LeadingE_e, PAPYieldNorm_e[ih]);
            gPAPYieldRatioToMB_VsLeadE[ih]->GetYaxis()->SetTitle(Form("#frac{%s%s}{#LT %s%s #GT_{MB}}", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data(), lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            gPAPYieldRatioToMB_VsLeadE[ih]->GetXaxis()->SetTitle("E_{leading}");
            gPAPYieldRatioToMB_VsLeadE[ih]->SetTitle("");
            gPAPYieldRatioToMB_VsLeadE[ih]->SetMarkerStyle(kFullCircle);
            gPAPYieldRatioToMB_VsLeadE[ih]->SetName(Form("YieldRatioToMB_VsLeadE_%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));

            // Vs Eff Energy
            gPAPYieldRatioToMB_VsEE[ih] = new TGraphErrors(nDiff, EffEnergy, PAPYieldNorm[ih], EffEnergy_e, PAPYieldNorm_e[ih]);
            gPAPYieldRatioToMB_VsEE[ih]->GetYaxis()->SetTitle(Form("#frac{%s%s}{#LT %s%s #GT_{MB}}", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data(), lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            gPAPYieldRatioToMB_VsEE[ih]->GetXaxis()->SetTitle("E_{effective}");
            gPAPYieldRatioToMB_VsEE[ih]->SetTitle("");
            gPAPYieldRatioToMB_VsEE[ih]->SetMarkerStyle(kFullCircle);
            gPAPYieldRatioToMB_VsEE[ih]->SetName(Form("YieldRatioToMB_VsEE_%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));

            // Vs Nch
            gPAPYield_VsNchNorm[ih] = new TGraphErrors(nDiff, NchNorm, PAPYield[ih], NchNorm_e, PAPYield_e[ih]);
            gPAPYield_VsNchNorm[ih]->GetYaxis()->SetTitle(Form("%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            gPAPYield_VsNchNorm[ih]->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
            gPAPYield_VsNchNorm[ih]->SetTitle("");
            gPAPYield_VsNchNorm[ih]->SetMarkerStyle(kFullCircle);
            gPAPYield_VsNchNorm[ih]->SetName(Form("Yield_VsNchNorm_%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            //
            gPAPYield_VsNch[ih] = new TGraphErrors(nDiff, Nch, PAPYield[ih], Nch_e, PAPYield_e[ih]);
            gPAPYield_VsNch[ih]->GetYaxis()->SetTitle(Form("%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            gPAPYield_VsNch[ih]->GetXaxis()->SetTitle("n_{ch}");
            gPAPYield_VsNch[ih]->SetTitle("");
            gPAPYield_VsNch[ih]->SetMarkerStyle(kFullCircle);
            gPAPYield_VsNch[ih]->SetName(Form("Yield_VsNch_%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));

            // Vs Lead Energy
            gPAPYield_VsLeadENorm[ih] = new TGraphErrors(nDiff, LeadingENorm, PAPYield[ih], LeadingENorm_e, PAPYield_e[ih]);
            gPAPYield_VsLeadENorm[ih]->GetYaxis()->SetTitle(Form("%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            gPAPYield_VsLeadENorm[ih]->GetXaxis()->SetTitle("E_{leading}/E_{leading}^{MB}");
            gPAPYield_VsLeadENorm[ih]->SetTitle("");
            gPAPYield_VsLeadENorm[ih]->SetMarkerStyle(kFullCircle);
            gPAPYield_VsLeadENorm[ih]->SetName(Form("Yield_VsLeadENorm_%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            //
            gPAPYield_VsLeadE[ih] = new TGraphErrors(nDiff, LeadingE, PAPYield[ih], LeadingE_e, PAPYield_e[ih]);
            gPAPYield_VsLeadE[ih]->GetYaxis()->SetTitle(Form("%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            gPAPYield_VsLeadE[ih]->GetXaxis()->SetTitle("E_{leading}");
            gPAPYield_VsLeadE[ih]->SetTitle("");
            gPAPYield_VsLeadE[ih]->SetMarkerStyle(kFullCircle);
            gPAPYield_VsLeadE[ih]->SetName(Form("Yield_VsLeadE_%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));

            // Vs Eff Energy
            gPAPYield_VsEE[ih] = new TGraphErrors(nDiff, EffEnergy, PAPYield[ih], EffEnergy_e, PAPYield_e[ih]);
            gPAPYield_VsEE[ih]->GetYaxis()->SetTitle(Form("%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
            gPAPYield_VsEE[ih]->GetXaxis()->SetTitle("E_{effective}");
            gPAPYield_VsEE[ih]->SetTitle("");
            gPAPYield_VsEE[ih]->SetMarkerStyle(kFullCircle);
            gPAPYield_VsEE[ih]->SetName(Form("Yield_VsEE_%s%s", lPartNames[lPAP[ih]].Data(), lPartNames[lPAP[ih]+1].Data()));
        }

        TFile * write = new TFile(outputfile, option);

        TDirectoryFile *dir = new TDirectoryFile(directory.Data(),directory.Data());
        dir->cd();
        hcalibV0M->Write();
        hcalibSPDCl->Write();

        if (isPythia) gNchVsNMPI->Write();
        if (isPythia) gLEVsNMPI->Write();
        gLEVsNch->Write();
        gNormLEVsNch->Write();
        if (isPythia)gEEVsNMPI->Write();

        TDirectoryFile *lYield = new TDirectoryFile("Yield", "h");
        lYield->cd();
        for (Int_t ih = 0; ih < nPart; ih++)
        {
            gYield_VsNchNorm[ih]->Write();
            gYield_VsNch[ih]->Write();
            gYield_VsLeadENorm[ih]->Write();
            gYield_VsLeadE[ih]->Write();
            gYield_VsEE[ih]->Write();
            if (isPythia)
            {
                gYield_VsNMPINorm[ih]->Write();
                gYield_VsNMPI[ih]->Write();
            }
        }
        for (Int_t ih = 0; ih < 9; ih++)
        {
            gPAPYield_VsNchNorm[ih]->Write();
            gPAPYield_VsNch[ih]->Write();
            gPAPYield_VsLeadENorm[ih]->Write();
            gPAPYield_VsLeadE[ih]->Write();
            gPAPYield_VsEE[ih]->Write();
        }
        dir->cd();

        TDirectoryFile *lYieldRatioToMB = new TDirectoryFile("YieldRatioToMB","h/(h)_{MB}");
        lYieldRatioToMB->cd();
        for(Int_t ih = 0; ih < nPart; ih++){
            gYieldRatioToMB_VsNchNorm[ih]->Write();
            gYieldRatioToMB_VsNch[ih]->Write();
            gYieldRatioToMB_VsLeadENorm[ih]->Write();
            gYieldRatioToMB_VsLeadE[ih]->Write();
            gYieldRatioToMB_VsEE[ih]->Write();
            if (isPythia) {
                gYieldRatioToMB_VsNMPINorm[ih]->Write();
                gYieldRatioToMB_VsNMPI[ih]->Write();
            }
        }
        for(Int_t ih = 0; ih < 9; ih++){
            gPAPYieldRatioToMB_VsNchNorm[ih]->Write();
            gPAPYieldRatioToMB_VsNch[ih]->Write();
            gPAPYieldRatioToMB_VsLeadENorm[ih]->Write();
            gPAPYieldRatioToMB_VsLeadE[ih]->Write();
            gPAPYieldRatioToMB_VsEE[ih]->Write();
        }
        dir->cd();

        TDirectoryFile *lMeanPt = new TDirectoryFile("MeanPt","<p_{T}>");
        lMeanPt->cd();
        for(Int_t ih = 0; ih < nPart; ih++){
            gAvPt_VsNch[ih]->Write();
            gAvPt_VsNchNorm[ih]->Write();
            gAvPt_VsLeadE[ih]->Write();
            gAvPt_VsLeadENorm[ih]->Write();
            gAvPt_VsEE[ih]->Write();
            if (isPythia) {
                gAvPt_VsNMPI[ih]->Write();
                gAvPt_VsNMPINorm[ih]->Write();}
        }
        dir->cd();
        write->cd();
        write->Close();*/
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

//==========================================================================================
void doerrorAB(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr){
    // Compute the error of a division

    Double_t out =
        (A/B) *
        TMath::Sqrt(
            (Aerr/A)*(Aerr/A) + (Berr/B)*(Berr/B)
        );
    Aerr = out;

    return;

}

//==========================================================================================
void doerrorABC(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr, Double_t C, Double_t Cerr){
    // Compute the error of a division

    Double_t out =
        (A/B/C) *
        TMath::Sqrt(
            (Aerr/A)*(Aerr/A) + (Berr/B)*(Berr/B) + (Cerr/C)*(Cerr/C)
        );
    Aerr = out;

    return;

}

//==========================================================================================
void doerror4(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr, Double_t C, Double_t Cerr, Double_t D, Double_t Derr){
    // Compute the error of a division

    Double_t out =
        (A/B/C/D) *
        TMath::Sqrt(
            (Aerr/A)*(Aerr/A) + (Berr/B)*(Berr/B) + (Cerr/C)*(Cerr/C) + (Derr/D)*(Derr/D)
        );
    Aerr = out;

    return;

}
