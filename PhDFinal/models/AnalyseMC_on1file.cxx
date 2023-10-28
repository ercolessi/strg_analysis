TH1D* dopercentile(TH1D* hvar);
Double_t getval(TH1D* h, Double_t perc);
void process(TString directory, TString outputfile, TString option, TString mcname, Double_t *SPDClustersperclow, Double_t *SPDClustersperchigh, Double_t *V0Mperclow, Double_t *V0Mperchigh, Int_t n);
void doerrorAB(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr);
void doerrorABC(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr, Double_t C, Double_t Cerr);
void doerror4(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr, Double_t C, Double_t Cerr, Double_t D, Double_t Derr);

void AnalyseMC_on1file(TString mcname = "Ropes")
{

    // class 0 --> standalone
    Double_t percentileV0M_low_0[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
    Double_t percentileV0M_high_0[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    Double_t percentileSPDClusters_low_0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t percentileSPDClusters_high_0[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Long_t n0 = sizeof(percentileV0M_low_0) / sizeof(Double_t);

    Double_t percentileSPDClusters_low_1[] = {10,10,10,10,10,10};
    Double_t percentileSPDClusters_high_1[] = {20,20,20,20,20,20};
    Double_t percentileV0M_low_1[] = {0, 5, 10, 20, 30, 50};
    Double_t percentileV0M_high_1[] = {5, 10, 20, 30, 50, 100};
    Long_t n1 = sizeof(percentileV0M_low_1) / sizeof(Double_t);

    Double_t percentileSPDClusters_low_2[] = {40,40,40,40,40,40,40};
    Double_t percentileSPDClusters_high_2[] = {50,50,50,50,50,50,50};
    Double_t percentileV0M_low_2[] = {0, 20, 30, 40, 50, 60, 70};
    Double_t percentileV0M_high_2[] = {20, 30, 40, 50, 60, 70, 100};
    Long_t n2 = sizeof(percentileV0M_low_2) / sizeof(Double_t);

    //Monash tuned classes
    /*Double_t percentileSPDClusters_low_3[] = {0, 10, 20, 30, 40, 80};
    Double_t percentileSPDClusters_high_3[] = {10, 20, 30, 40, 80, 90};
    Double_t percentileV0M_low_3[] = {40, 30, 20, 10, 10, 0};
    Double_t percentileV0M_high_3[] = {50, 40, 40, 30, 20, 10};
    Long_t n3 = sizeof(percentileV0M_low_3) / sizeof(Double_t);*/

    // fixed low ZN
    Double_t percentileSPDClusters_low_3[] = {0, 10, 20, 30, 50};
    Double_t percentileSPDClusters_high_3[] = {20, 30, 40, 50, 100};
    Double_t percentileV0M_low_3[] = {40, 30, 30, 20, 0};
    Double_t percentileV0M_high_3[] = {60, 70, 50, 50, 30};
    Long_t n3 = sizeof(percentileSPDClusters_low_3) / sizeof(Double_t);

    // fixed very low ZN
    Double_t percentileSPDClusters_low_4[] = {0, 10, 20, 30};
    Double_t percentileSPDClusters_high_4[] = {10, 20, 30, 50};
    Double_t percentileV0M_low_4[] = {20, 10, 0, 0};
    Double_t percentileV0M_high_4[] = {30, 30, 20, 10};
    Long_t n4 = sizeof(percentileV0M_low_4) / sizeof(Double_t);

    TString outputfile = Form("Results_%s_SPDClustersV0M.root",mcname.Data());

    TString mcname_ = "PythiaMonash_FN_1340";
    if (mcname.Contains("Ropes")){
        mcname_ = "PythiaRopes_FN_1388";
    }

    process("Standalone", outputfile, "RECREATE", mcname_, percentileSPDClusters_low_0, percentileSPDClusters_high_0, percentileV0M_low_0, percentileV0M_high_0, n0);
    process("kHighMult", outputfile, "UPDATE", mcname_, percentileSPDClusters_low_1, percentileSPDClusters_high_1, percentileV0M_low_1, percentileV0M_high_1, n1);
    process("kLowMult", outputfile, "UPDATE", mcname_, percentileSPDClusters_low_2, percentileSPDClusters_high_2, percentileV0M_low_2, percentileV0M_high_2, n2);
    process("kHighZN", outputfile, "UPDATE", mcname_, percentileSPDClusters_low_3, percentileSPDClusters_high_3, percentileV0M_low_3, percentileV0M_high_3, n3);
    process("kLowZN", outputfile, "UPDATE", mcname_, percentileSPDClusters_low_4, percentileSPDClusters_high_4, percentileV0M_low_4, percentileV0M_high_4, n4);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void process(TString directory = "", TString outputfile = "", TString option = "RECREATE", TString mcname = "EPOSLHC", Double_t *SPDClustersperclow = 0x0, Double_t *SPDClustersperchigh = 0x0, Double_t *V0Mperclow = 0x0, Double_t *V0Mperchigh = 0x0, Int_t n = 0)
{

    //Open file
    TFile* file = TFile::Open(Form("Files/%s.root",mcname.Data()));
    file->cd("PWGLF_MCPredictionsStrgVsMultVsZDC");
    TList* list  = (TList*)file->FindObjectAny("cList");
    Bool_t isPythia = kFALSE;
    if (mcname.Contains("Pythia")) isPythia = kTRUE;

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
    cout << " Is Pythia ? ........" << isPythia << endl;
    cout << "\n-----------------------------------------------" << endl;

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
    f2DHistINELgt0SPDV0M = (TH2D *)list->FindObject("f2DHistINELgt0SPDV0M");
    f2DHistLeadingESPDV0M = (TH2D *)list->FindObject("f2DHistLeadingESPDV0M");
    f2DHistEffEnergySPDV0M = (TH2D *)list->FindObject("f2DHistEffEnergySPDV0M");
    f2DHistNchSPDV0M = (TH2D *)list->FindObject("f2DHistNchSPDV0M");
    f2DHistNMPISPDV0M = (TH2D *)list->FindObject("f2DHistNMPISPDV0M");
    for (Int_t ih = 0; ih < nPart; ih++) {
        fHistPt[ih] = (TH1D *)list->FindObject(Form("fHistPt_%s", lPartNames[ih].Data()));
        f2DHistAvPtSPDV0M[ih] = (TH2D *)list->FindObject(Form("f2DHistAvPtSPDV0M_%s", lPartNames[ih].Data()));
        f2DHistPartSPDV0M[ih] = (TH2D *)list->FindObject(Form("f2DHistPartSPDV0M_%s", lPartNames[ih].Data()));
    }

    // Percentile calibrations
    TH1D *hcalibV0M = dopercentile(fHistV0MMult);
    hcalibV0M->SetName("hcalibV0M");
    TH1D *hcalibSPDCl = dopercentile(fHistSPDClusters);
    hcalibSPDCl->SetName("hcalibSPDCl");

    // Convert to values V0M and SPDCl
    const int nDiff = n;
    Double_t V0Mvallow[nDiff], V0Mvalhigh[nDiff];
    Double_t SPDClvallow[nDiff], SPDClvalhigh[nDiff];
    for(int i = 0; i < n; i++){
        V0Mvallow[i] = getval(hcalibV0M,V0Mperclow[i]);
        V0Mvalhigh[i] = getval(hcalibV0M,V0Mperchigh[i]);
        SPDClvallow[i] = getval(hcalibSPDCl, SPDClustersperclow[i]);
        SPDClvalhigh[i] = getval(hcalibSPDCl, SPDClustersperchigh[i]);
    }

    //Variable definition
    Double_t EvtMB, EvtMB_e; // errors must be double for the method TH1::IntegralAndError()
    Double_t NchMB, NchMB_e;
    Double_t EffEnergyMB, EffEnergyMB_e;
    Double_t LeadingEMB, LeadingEMB_e;
    Double_t NMPIMB, NMPIMB_e;
    Double_t YieldMB[nPart], YieldMB_e[nPart];
    Double_t Evt[nDiff], Evt_e[nDiff];
    Double_t Nch[nDiff], Nch_e[nDiff];
    Double_t NchNorm[nDiff], NchNorm_e[nDiff];
    Double_t NMPI[nDiff], NMPI_e[nDiff];
    Double_t NMPINorm[nDiff], NMPINorm_e[nDiff];
    Double_t EffEnergy[nDiff], EffEnergy_e[nDiff];
    Double_t LeadingE[nDiff], LeadingE_e[nDiff];
    Double_t LeadingENorm[nDiff], LeadingENorm_e[nDiff];
    Double_t Yield[nPart][nDiff], Yield_e[nPart][nDiff];
    Double_t YieldNorm[nPart][nDiff], YieldNorm_e[nPart][nDiff];
    Double_t AvPt[nPart][nDiff], AvPt_e[nPart][nDiff];
    Double_t Npart[nPart][nDiff];
    Int_t lPAP[9] = {0, 2, 4, 7, 9, 11, 14, 16, 18}; // places in the array which are particle and antiparticle
    Double_t PAPYieldNorm[9][nDiff], PAPYieldNorm_e[9][nDiff];
    Double_t PAPYield[9][nDiff], PAPYield_e[9][nDiff];
    Double_t PAPYieldMB[9], PAPYieldMB_e[9];
    //-------------------------------------------------------
    //--------------------- INEL>0 (MB)----------------------
    //-------------------------------------------------------

    EvtMB = f2DHistINELgt0SPDV0M -> IntegralAndError(
        1, f2DHistINELgt0SPDV0M->GetNbinsX(),
        1, f2DHistINELgt0SPDV0M->GetNbinsY(),
        EvtMB_e
    );
    //
    NchMB = f2DHistNchSPDV0M -> IntegralAndError(
        1, f2DHistNchSPDV0M->GetNbinsX(),
        1, f2DHistNchSPDV0M->GetNbinsY(),
        NchMB_e
    );
    doerrorAB(NchMB, NchMB_e, EvtMB, EvtMB_e);
    NchMB /= EvtMB;
    //
    EffEnergyMB = f2DHistEffEnergySPDV0M -> IntegralAndError(
        1, f2DHistEffEnergySPDV0M->GetNbinsX(),
        1, f2DHistEffEnergySPDV0M->GetNbinsY(),
        EffEnergyMB_e
        );
    doerrorAB(EffEnergyMB, EffEnergyMB_e, EvtMB, EvtMB_e);
    EffEnergyMB /= EvtMB;
    //
    LeadingEMB = f2DHistLeadingESPDV0M -> IntegralAndError(
        1, f2DHistLeadingESPDV0M->GetNbinsX(),
        1, f2DHistLeadingESPDV0M->GetNbinsY(),
        LeadingEMB_e
        );
    doerrorAB(LeadingEMB, LeadingEMB_e, EvtMB, EvtMB_e);
    LeadingEMB /= EvtMB;
    //
    if (isPythia){
    NMPIMB = f2DHistNMPISPDV0M -> IntegralAndError(
        1, f2DHistNMPISPDV0M->GetNbinsX(),
        1, f2DHistNMPISPDV0M->GetNbinsY(),
        NMPIMB_e
        );
    doerrorAB(NMPIMB, NMPIMB_e, EvtMB, EvtMB_e);
    NMPIMB /= EvtMB;
    }
    //
    for(Int_t ih = 0; ih < nPart; ih++){
        YieldMB[ih]  = f2DHistPartSPDV0M[ih] -> IntegralAndError(
            1, f2DHistPartSPDV0M[ih]->GetNbinsX(),
            1, f2DHistPartSPDV0M[ih]->GetNbinsY(),
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

    //-------------------------------------------------------
    //--------------------- Selections ----------------------
    //-------------------------------------------------------

    for(int k = 0; k < nDiff; k++){ // SPD Cluster loop

            Evt[k] = f2DHistINELgt0SPDV0M -> IntegralAndError(
                f2DHistINELgt0SPDV0M->GetXaxis()->FindBin(SPDClvalhigh[k]), f2DHistINELgt0SPDV0M->GetXaxis()->FindBin(SPDClvallow[k]),
                f2DHistINELgt0SPDV0M->GetYaxis()->FindBin(V0Mvalhigh[k]), f2DHistINELgt0SPDV0M->GetYaxis()->FindBin(V0Mvallow[k]),
                Evt_e[k]
            );
            //
            NchNorm[k] = f2DHistNchSPDV0M -> IntegralAndError(
                f2DHistNchSPDV0M->GetXaxis()->FindBin(SPDClvalhigh[k]), f2DHistNchSPDV0M->GetXaxis()->FindBin(SPDClvallow[k]),
                f2DHistNchSPDV0M->GetYaxis()->FindBin(V0Mvalhigh[k]), f2DHistNchSPDV0M->GetYaxis()->FindBin(V0Mvallow[k]),
                NchNorm_e[k]
            );
            Nch_e[k] = NchNorm_e[k];
            doerrorAB(NchNorm[k], Nch_e[k], Evt[k], Evt_e[k]);
            doerrorABC(NchNorm[k], NchNorm_e[k], Evt[k], Evt_e[k], NchMB, NchMB_e);
            NchNorm[k] /= Evt[k];
            Nch[k] = NchNorm[k];
            NchNorm[k] /= NchMB;

            cout << "Nch[" << k << "] = " << Nch[k] << " +- " << Nch_e[k] << endl;
            //
            EffEnergy[k] = f2DHistEffEnergySPDV0M -> IntegralAndError(
                f2DHistEffEnergySPDV0M->GetXaxis()->FindBin(SPDClvalhigh[k]), f2DHistEffEnergySPDV0M->GetXaxis()->FindBin(SPDClvallow[k]),
                f2DHistEffEnergySPDV0M->GetYaxis()->FindBin(V0Mvalhigh[k]), f2DHistEffEnergySPDV0M->GetYaxis()->FindBin(V0Mvallow[k]),
                EffEnergy_e[k]
            );
            doerrorAB(EffEnergy[k], EffEnergy_e[k], Evt[k], Evt_e[k]);
            EffEnergy[k] /= Evt[k];
            //
            LeadingENorm[k] = f2DHistLeadingESPDV0M -> IntegralAndError(
                f2DHistLeadingESPDV0M->GetXaxis()->FindBin(SPDClvalhigh[k]), f2DHistLeadingESPDV0M->GetXaxis()->FindBin(SPDClvallow[k]),
                f2DHistLeadingESPDV0M->GetYaxis()->FindBin(V0Mvalhigh[k]), f2DHistLeadingESPDV0M->GetYaxis()->FindBin(V0Mvallow[k]),
                LeadingENorm_e[k]
            );
            LeadingE_e[k] = LeadingENorm_e[k];
            doerrorAB(LeadingENorm[k], LeadingE_e[k], Evt[k], Evt_e[k]);
            doerrorABC(LeadingENorm[k], LeadingENorm_e[k], Evt[k], Evt_e[k], LeadingEMB, LeadingEMB_e);
            LeadingENorm[k] /= Evt[k];
            LeadingE[k] = LeadingENorm[k];
            LeadingENorm[k] /= LeadingEMB;
            //
            if (isPythia){
            NMPINorm[k] = f2DHistNMPISPDV0M -> IntegralAndError(
                f2DHistNMPISPDV0M->GetXaxis()->FindBin(SPDClvalhigh[k]), f2DHistNMPISPDV0M->GetXaxis()->FindBin(SPDClvallow[k]),
                f2DHistNMPISPDV0M->GetYaxis()->FindBin(V0Mvalhigh[k]), f2DHistNMPISPDV0M->GetYaxis()->FindBin(V0Mvallow[k]),
                NMPINorm_e[k]
            );
            NMPI_e[k] = NMPINorm_e[k];
            doerrorAB(NMPINorm[k], NMPI_e[k], Evt[k], Evt_e[k]);
            doerrorABC(NMPINorm[k], NMPINorm_e[k], Evt[k], Evt_e[k], NMPIMB, NMPIMB_e);
            NMPINorm[k] /= Evt[k];
            NMPI[k] = NMPINorm[k];
            NMPINorm[k] /= NMPIMB;
            }
            //
            for(Int_t ih = 0; ih < nPart; ih++){
                YieldNorm[ih][k]  = f2DHistPartSPDV0M[ih] -> IntegralAndError(
                    f2DHistPartSPDV0M[ih]->GetXaxis()->FindBin(SPDClvalhigh[k]), f2DHistPartSPDV0M[ih]->GetXaxis()->FindBin(SPDClvallow[k]),
                    f2DHistPartSPDV0M[ih]->GetYaxis()->FindBin(V0Mvalhigh[k]), f2DHistPartSPDV0M[ih]->GetYaxis()->FindBin(V0Mvallow[k]),
                    YieldNorm_e[ih][k]
                );
                Yield_e[ih][k] = YieldNorm_e[ih][k];
                doerrorAB(YieldNorm[ih][k], Yield_e[ih][k], Evt[k], Evt_e[k]);
                doerrorABC(YieldNorm[ih][k], YieldNorm_e[ih][k], Evt[k], Evt_e[k], YieldMB[ih], YieldMB_e[ih]);
                Npart[ih][k] = YieldNorm[ih][k];
                YieldNorm[ih][k] /= Evt[k];
                Yield[ih][k] = YieldNorm[ih][k];
                YieldNorm[ih][k] /= YieldMB[ih];

                AvPt[ih][k] = f2DHistAvPtSPDV0M[ih]->IntegralAndError(
                        f2DHistAvPtSPDV0M[ih]->GetXaxis()->FindBin(SPDClvalhigh[k]), f2DHistAvPtSPDV0M[ih]->GetXaxis()->FindBin(SPDClvallow[k]),
                        f2DHistAvPtSPDV0M[ih]->GetYaxis()->FindBin(V0Mvalhigh[k]), f2DHistAvPtSPDV0M[ih]->GetYaxis()->FindBin(V0Mvallow[k]),
                        AvPt_e[ih][k]
                );
                AvPt[ih][k] /= Npart[ih][k];
                AvPt_e[ih][k] /= Npart[ih][k];
            }

            for(Int_t ih = 0; ih < 9; ih++){
                PAPYield[ih][k] = Yield[lPAP[ih]][k] + Yield[lPAP[ih] + 1][k];
                PAPYield_e[ih][k] = TMath::Sqrt((Yield_e[lPAP[ih]][k] * Yield_e[lPAP[ih]][k] + Yield_e[lPAP[ih] + 1][k] * Yield_e[lPAP[ih] + 1][k]));
                PAPYieldNorm[ih][k] = PAPYield[ih][k];
                PAPYieldNorm_e[ih][k] = PAPYield_e[ih][k];
                doerrorAB(PAPYieldNorm[ih][k], PAPYieldNorm_e[ih][k], PAPYieldMB[ih], PAPYieldMB_e[ih]);
                PAPYieldNorm[ih][k] /= PAPYieldMB[ih];
            }
    } // end  loop

    // Yields
    TGraphErrors *gYield_VsNchNorm[nPart];
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
    write->Close();
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
