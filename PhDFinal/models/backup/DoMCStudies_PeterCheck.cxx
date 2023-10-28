TH1D* dopercentile(TH1D* hvar);
Double_t getval(TH1D* h, Double_t perc);
void process(TString directory, TString outputfile , TString option , TString mcname , Double_t* SPDClperc, Int_t nSPDCl, Double_t* V0Mperc, Int_t nV0M, Bool_t IsFullSim);
void doerror(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr);
void doerror3(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr, Double_t C, Double_t Cerr);
void doerror4(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr, Double_t C, Double_t Cerr, Double_t D, Double_t Derr);

void DoMCStudies_PeterCheck(TString mcname = "PythiaMonash_Train2627", Bool_t IsFullSim = kFALSE)
{

    TString Type = "Full";
    if (!IsFullSim) Type = "Fast";

    Double_t V0Mstandalone[]   = {0,1,5,10,15,20,30,40,50,70,100};
    const int nV0Mstandalone = sizeof(V0Mstandalone)/sizeof(Double_t);
    Double_t V0M_SPD1020[]   = {0,5,10,20,30,50,100};//{0,5,10,20,30,40,50,100};
    const int nV0M_SPD1020 = sizeof(V0M_SPD1020)/sizeof(Double_t);
    Double_t V0M_SPD4050[] = {0,20,30,40,50,60,70,100};
    const int nV0M_SPD4050 = sizeof(V0M_SPD4050)/sizeof(Double_t);
    Double_t V0M_V0M1020[] = {10,20};
    const int nV0M_V0M1020 = sizeof(V0M_V0M1020)/sizeof(Double_t);
    Double_t V0M_V0M4050[] = {40,50};
    const int nV0M_V0M4050 = sizeof(V0M_V0M4050)/sizeof(Double_t);
    //
    Double_t SPDCl_V0Mstandalone[] = {0,100};
    const int nSPDCl_V0Mstandalone = sizeof(SPDCl_V0Mstandalone)/sizeof(Double_t);
    Double_t SPDCl_SPD1020[] = {10,20};
    const int nSPDCl_SPD1020 = sizeof(SPDCl_SPD1020)/sizeof(Double_t);
    Double_t SPDCl_SPD4050[] = {40,50};
    const int nSPDCl_SPD4050 = sizeof(SPDCl_SPD4050)/sizeof(Double_t);
    Double_t SPDCl_V0M1020[] = {0,5,10,20,30,40,50,100};
    const int nSPDCl_V0M1020 = sizeof(SPDCl_V0M1020)/sizeof(Double_t);
    Double_t SPDCl_V0M4050[] = {0,10,20,30,40,50,60,70,100};
    const int nSPDCl_V0M4050 = sizeof(SPDCl_V0M4050)/sizeof(Double_t);

    TString outputfile = Form("PeterCheck_ResultsMC%sSim_%s.root",Type.Data(),mcname.Data());

    process("Standalone",outputfile,"RECREATE",mcname, SPDCl_V0Mstandalone, nSPDCl_V0Mstandalone, V0Mstandalone, nV0Mstandalone, IsFullSim);
    process("SPDCl1020",outputfile,"UPDATE",mcname, SPDCl_SPD1020, nSPDCl_SPD1020, V0M_SPD1020, nV0M_SPD1020, IsFullSim);
    process("SPDCl4050",outputfile,"UPDATE",mcname, SPDCl_SPD4050, nSPDCl_SPD4050, V0M_SPD4050, nV0M_SPD4050, IsFullSim);
    process("V0M1020",outputfile,"UPDATE",mcname, SPDCl_V0M1020, nSPDCl_V0M1020, V0M_V0M1020, nV0M_V0M1020, IsFullSim);
    process("V0M4050",outputfile,"UPDATE",mcname, SPDCl_V0M4050, nSPDCl_V0M4050, V0M_V0M4050, nV0M_V0M4050, IsFullSim);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void process(TString directory = "", TString outputfile = "", TString option = "RECREATE", TString mcname = "EPOSLHC", Double_t* SPDClperc = 0x0, Int_t nSPDCl = 0, Double_t* V0Mperc = 0x0, Int_t nV0M = 0, Bool_t IsFullSim = kTRUE){

    Bool_t FixMult = kTRUE; // flag to understand in which double differential analysis we are
    if (nV0M < nSPDCl) FixMult = kFALSE;

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

    //Particles
    Int_t lPAP[9] = { // places in the array which are particle and antiparticle
        0,
        2,
        4,
        7,
        9,
        11,
        14,
        16,
        18
    };

    //Histograms
    TH1D *fHistEventCounter;
    TH1D *fHistV0MMult;
    TH1D *fHistV0AMult;
    TH1D *fHistV0CMult;
    TH1D *fHistMult05;
    TH1D *fHistMult08;
    TH1D *fHistMult08to15;
    TH1D *fHistSPDClusters;
    TH1D *fHistNMPI;
    TH1D *fHistQ2;
    TH1D *fHistb;
    TH1D *fHistLeadingE;
    TH1D *fHistEffEnergy;
    TH1D *fHistPt[nPart];
    TH2D *f2DHistPartSPDV0M[nPart];
    TH2D *f2DHistPartRecoPercSPDV0M[nPart];
    TH2D *f2DHistAvPtSPDV0M[nPart];
    TH2D *f2DHistINELgt0SPDV0M;
    TH2D *f2DHistLeadingESPDV0M;
    TH2D *f2DHistEffEnergySPDV0M;
    TH2D *f2DHistNchSPDV0M;
    TH2D *f2DHistNMPISPDV0M;
    TH2D *f2DHistQ2SPDV0M;
    TH2D *f2DHistbSPDV0M;

    //Get histograms
    fHistEventCounter = (TH1D*)list->FindObject("fHistEventCounter");
    fHistV0MMult      = (TH1D*)list->FindObject("fHistV0MMult");
    fHistV0AMult      = (TH1D*)list->FindObject("fHistV0AMult");
    fHistV0CMult      = (TH1D*)list->FindObject("fHistV0CMult");
    fHistMult05       = (TH1D*)list->FindObject("fHistMult05");
    fHistMult08       = (TH1D*)list->FindObject("fHistMult08");
    fHistMult08to15   = (TH1D*)list->FindObject("fHistMult08to15");
    fHistSPDClusters  = (TH1D*)list->FindObject("fHistSPDClusters");
    fHistNMPI         = (TH1D*)list->FindObject("fHistNMPI");
    fHistQ2           = (TH1D*)list->FindObject("fHistQ2");
    fHistb            = (TH1D*)list->FindObject("fHistb");
    fHistLeadingE     = (TH1D*)list->FindObject("fHistLeadingE");
    fHistEffEnergy    = (TH1D*)list->FindObject("fHistEffEnergy");
    f2DHistINELgt0SPDV0M   = (TH2D*)list->FindObject("f2DHistINELgt0Nch0815V0M");
    f2DHistLeadingESPDV0M  = (TH2D*)list->FindObject("f2DHistLeadingENch0815V0M");
    f2DHistEffEnergySPDV0M = (TH2D*)list->FindObject("f2DHistEffEnergyNch0815V0M");
    f2DHistNchSPDV0M       = (TH2D*)list->FindObject("f2DHistNchNch0815V0M");
    f2DHistNMPISPDV0M      = (TH2D*)list->FindObject("f2DHistNMPINch0815V0M");
    f2DHistQ2SPDV0M        = (TH2D*)list->FindObject("f2DHistQ2Nch0815V0M");
    f2DHistbSPDV0M         = (TH2D*)list->FindObject("f2DHistbNch0815V0M");
    for(Int_t ih = 0; ih < nPart; ih++){
        fHistPt[ih]                = (TH1D*)list->FindObject(Form("fHistPt_%s",lPartNames[ih].Data()));
        if (IsFullSim){
            f2DHistPartSPDV0M[ih]      = (TH2D*)list->FindObject(Form("f2DHistPartRecoPercSPDV0M_%s",lPartNames[ih].Data()));
        } else{
            f2DHistPartSPDV0M[ih]      = (TH2D*)list->FindObject(Form("f2DHistPartNch0815V0M_%s",lPartNames[ih].Data()));
        }
        f2DHistAvPtSPDV0M[ih]      = (TH2D*)list->FindObject(Form("f2DHistAvPtSPDV0M_%s",lPartNames[ih].Data()));
    }

    //Getpercentile calibrations
    TH1D* hcalibV0M = dopercentile(fHistV0MMult);
    TH1D* hcalibSPDCl = dopercentile(fHistMult08to15);
    TH1D* hcalibMult08to15 = dopercentile(fHistMult08to15);

    new TCanvas;
    hcalibSPDCl->Draw();

    Double_t V0Mval[nV0M];
    for(int i = 0; i < nV0M; i++){
        V0Mval[i] = getval(hcalibV0M,V0Mperc[i]);
    }
    Double_t SPDClval[nSPDCl];
    for(int i = 0; i < nSPDCl; i++){
        SPDClval[i] = getval(hcalibSPDCl,SPDClperc[i]);
        cout <<      SPDClperc[i] << "  " << SPDClval[i] << endl;
    }
    if (IsFullSim){
        for(int i = 0; i < nV0M; i++){
            V0Mval[i] = V0Mperc[i];
       }
        for(int i = 0; i < nSPDCl; i++){
            SPDClval[i] = SPDClperc[i];
        }
    }
    const int nV0M_ = nV0M - 1;
    const int nSPDCl_ = nSPDCl - 1;
    int n_;
    if (FixMult) n_ = nV0M_;
    if (!FixMult) n_ = nSPDCl_;
    const int nDiff = n_;

    Double_t Evt[nDiff];
    Double_t Evt_e[nDiff]; //errors must be double for the method TH1::IntegralAndError()
    Double_t Nch[nDiff];
    Double_t Nch_e[nDiff];
    Double_t NchNorm[nDiff];
    Double_t NchNorm_e[nDiff];
    Double_t NMPI[nDiff];
    Double_t NMPI_e[nDiff];
    Double_t NMPINorm[nDiff];
    Double_t NMPINorm_e[nDiff];
    Double_t EffEnergy[nDiff];
    Double_t EffEnergy_e[nDiff];
    Double_t LeadingE[nDiff];
    Double_t LeadingE_e[nDiff];
    Double_t LeadingENorm[nDiff];
    Double_t LeadingENorm_e[nDiff];
    Double_t Yield[nPart][nDiff];
    Double_t Yield_forPi[nPart][nDiff];
    Double_t Yield_forPi_e[nPart][nDiff];
    Double_t Yield_e[nPart][nDiff];
    Double_t PAPYield[9][nDiff];
    Double_t PAPYield_e[9][nDiff];
    Double_t YieldNorm[nPart][nDiff];
    Double_t YieldNorm_e[nPart][nDiff];
    Double_t PAPYieldNorm[9][nDiff];
    Double_t PAPYieldNorm_e[9][nDiff];
    Double_t PAPYield_forPi[9][nDiff];
    Double_t PAPYield_forPi_e[9][nDiff];
    Double_t PAPYieldNormPich[9][nDiff];
    Double_t PAPYieldNormPich_e[9][nDiff];
    Double_t PAPYieldNormPi0[9][nDiff];
    Double_t PAPYieldNormPi0_e[9][nDiff];
    Double_t AvPt[nPart][nDiff];
    Double_t AvPt_e[nPart][nDiff];
    Double_t PAPAvPt[9][nDiff];
    Double_t PAPAvPt_e[9][nDiff];
    Double_t Npart[nPart][nDiff];

    // INEL>0 (MB)
    Double_t EvtMB_e, NchMB_e, EffEnergyMB_e, LeadingEMB_e, NMPIMB_e;
    Double_t EvtMB =
    f2DHistINELgt0SPDV0M -> IntegralAndError(
        1, f2DHistINELgt0SPDV0M->GetNbinsX(),
        1, f2DHistINELgt0SPDV0M->GetNbinsY(),
        EvtMB_e
        );
    //
    Double_t NchMB =
    f2DHistNchSPDV0M -> IntegralAndError(
        1, f2DHistNchSPDV0M->GetNbinsX(),
        1, f2DHistNchSPDV0M->GetNbinsY(),
        NchMB_e
        );
    doerror(NchMB, NchMB_e, EvtMB, EvtMB_e);
    NchMB /= EvtMB;
    cout << " MB multiplicity... " << NchMB << endl;
    //
    Double_t EffEnergyMB =
    f2DHistEffEnergySPDV0M -> IntegralAndError(
        1, f2DHistEffEnergySPDV0M->GetNbinsX(),
        1, f2DHistEffEnergySPDV0M->GetNbinsY(),
        EffEnergyMB_e
        );
    doerror(EffEnergyMB, EffEnergyMB_e, EvtMB, EvtMB_e);
    EffEnergyMB /= EvtMB;
    //
    Double_t LeadingEMB =
    f2DHistLeadingESPDV0M -> IntegralAndError(
        1, f2DHistLeadingESPDV0M->GetNbinsX(),
        1, f2DHistLeadingESPDV0M->GetNbinsY(),
        LeadingEMB_e
        );
    doerror(LeadingEMB, LeadingEMB_e, EvtMB, EvtMB_e);
    LeadingEMB /= EvtMB;
    //
    Double_t NMPIMB;
    if (isPythia){
    NMPIMB =
    f2DHistNMPISPDV0M -> IntegralAndError(
        1, f2DHistNMPISPDV0M->GetNbinsX(),
        1, f2DHistNMPISPDV0M->GetNbinsY(),
        NMPIMB_e
        );
    doerror(NMPIMB, NMPIMB_e, EvtMB, EvtMB_e);
    NMPIMB /= EvtMB;
    }
    //
    Double_t YieldMB[nPart];
    Double_t YieldMB_e[nPart];
    Double_t YieldPAPMB[nPart];
    for(Int_t ih = 0; ih < nPart; ih++){
        YieldMB[ih]  =
        f2DHistPartSPDV0M[ih] -> IntegralAndError(
            1, f2DHistPartSPDV0M[ih]->GetNbinsX(),
            1, f2DHistPartSPDV0M[ih]->GetNbinsY(),
            YieldMB_e[ih]
            );
        doerror3(YieldMB[ih], YieldMB_e[ih],EvtMB, EvtMB_e, NchMB, NchMB_e);
        YieldMB[ih] /=  EvtMB;
        YieldPAPMB[ih] = YieldMB[ih];
        YieldMB[ih] /=  NchMB;
    }

    for(int j = 0; j < nSPDCl_; j++){
        for(int i = 0; i < nV0M_; i++){

            int k = j;
            if(FixMult) k = i;

            Evt[k] =
            f2DHistINELgt0SPDV0M -> IntegralAndError(
                f2DHistINELgt0SPDV0M->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistINELgt0SPDV0M->GetXaxis()->FindBin(SPDClval[j]),
                f2DHistINELgt0SPDV0M->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistINELgt0SPDV0M->GetYaxis()->FindBin(V0Mval[i]),
                Evt_e[k]
                );

            // Nch
            NchNorm[k] =
            f2DHistNchSPDV0M -> IntegralAndError(
                f2DHistNchSPDV0M->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistNchSPDV0M->GetXaxis()->FindBin(SPDClval[j]),
                f2DHistNchSPDV0M->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistNchSPDV0M->GetYaxis()->FindBin(V0Mval[i]),
                NchNorm_e[k]
                );
            Nch_e[k] = NchNorm_e[k];
            doerror(NchNorm[k], Nch_e[k], Evt[k], Evt_e[k]);
            doerror3(NchNorm[k], NchNorm_e[k], Evt[k], Evt_e[k], NchMB, NchMB_e);
            NchNorm[k] /= Evt[k];
            Nch[k] = NchNorm[k];
            NchNorm[k] /= NchMB;

            // Nch
            if (isPythia){
            NMPINorm[k] =
            f2DHistNMPISPDV0M -> IntegralAndError(
                f2DHistNMPISPDV0M->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistNMPISPDV0M->GetXaxis()->FindBin(SPDClval[j]),
                f2DHistNMPISPDV0M->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistNMPISPDV0M->GetYaxis()->FindBin(V0Mval[i]),
                NMPINorm_e[k]
                );
            NMPI_e[k] = NMPINorm_e[k];
            doerror(NMPINorm[k], NMPI_e[k], Evt[k], Evt_e[k]);
            doerror3(NMPINorm[k], NMPINorm_e[k], Evt[k], Evt_e[k], NMPIMB, NMPIMB_e);
            NMPINorm[k] /= Evt[k];
            NMPI[k] = NMPINorm[k];
            NMPINorm[k] /= NMPIMB;
            }

            // Effective energy
            EffEnergy[k] =
            f2DHistEffEnergySPDV0M -> IntegralAndError(
                f2DHistEffEnergySPDV0M->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistEffEnergySPDV0M->GetXaxis()->FindBin(SPDClval[j]),
                f2DHistEffEnergySPDV0M->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistEffEnergySPDV0M->GetYaxis()->FindBin(V0Mval[i]),
                EffEnergy_e[k]
                );
            doerror(EffEnergy[k], EffEnergy_e[k], Evt[k], Evt_e[k]);
            EffEnergy[k] /= Evt[k];

            // Leading energy
            LeadingENorm[k] =
            f2DHistLeadingESPDV0M -> IntegralAndError(
                f2DHistLeadingESPDV0M->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistLeadingESPDV0M->GetXaxis()->FindBin(SPDClval[j]),
                f2DHistLeadingESPDV0M->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistLeadingESPDV0M->GetYaxis()->FindBin(V0Mval[i]),
                LeadingENorm_e[k]
                );
            LeadingE_e[k] = LeadingENorm_e[k];
            doerror(LeadingENorm[k], LeadingE_e[k], Evt[k], Evt_e[k]);
            doerror3(LeadingENorm[k], LeadingENorm_e[k], Evt[k], Evt_e[k], LeadingEMB, LeadingEMB_e);
            LeadingENorm[k] /= Evt[k];
            LeadingE[k] = LeadingENorm[k];
            LeadingENorm[k] /= LeadingEMB;

            // Particle yield
            for(Int_t ih = 0; ih < nPart; ih++){
                YieldNorm[ih][k]  =
                f2DHistPartSPDV0M[ih] -> IntegralAndError(
                    f2DHistPartSPDV0M[ih]->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistPartSPDV0M[ih]->GetXaxis()->FindBin(SPDClval[j]),
                    f2DHistPartSPDV0M[ih]->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistPartSPDV0M[ih]->GetYaxis()->FindBin(V0Mval[i]),
                    YieldNorm_e[ih][k]
                    );
                Yield_e[ih][k] = YieldNorm_e[ih][k];
                Yield_forPi_e[ih][k] = YieldNorm_e[ih][k];
                doerror(YieldNorm[ih][k], Yield_forPi_e[ih][k], Evt[k], Evt_e[k]);
                doerror3(YieldNorm[ih][k], Yield_e[ih][k], Evt[k], Evt_e[k], Nch[k], Nch_e[k]);
                doerror4(YieldNorm[ih][k], YieldNorm_e[ih][k], Evt[k], Evt_e[k], Nch[k], Nch_e[k], YieldMB[ih], YieldMB_e[ih]);
                Npart[ih][k] = YieldNorm[ih][k];
                YieldNorm[ih][k] /= Evt[k];
                Yield_forPi[ih][k] = YieldNorm[ih][k];
                YieldNorm[ih][k] /=  Nch[k];
                Yield[ih][k] = YieldNorm[ih][k];
                YieldNorm[ih][k] /= YieldMB[ih];
            }

            //Particle + AntiParticle yield
            for(Int_t ih = 0; ih < 9; ih++){
                Double_t parterr, parterr_norm, antiparterr, antiparterr_norm;
                PAPYieldNorm[ih][k]  =
                f2DHistPartSPDV0M[lPAP[ih]] -> IntegralAndError(
                    f2DHistPartSPDV0M[lPAP[ih]]->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistPartSPDV0M[lPAP[ih]]->GetXaxis()->FindBin(SPDClval[j]),
                    f2DHistPartSPDV0M[lPAP[ih]]->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistPartSPDV0M[lPAP[ih]]->GetYaxis()->FindBin(V0Mval[i]),
                    parterr_norm
                    );
                PAPYieldNorm[ih][k] +=
                f2DHistPartSPDV0M[lPAP[ih]+1] -> IntegralAndError(
                    f2DHistPartSPDV0M[lPAP[ih]+1]->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistPartSPDV0M[lPAP[ih]+1]->GetXaxis()->FindBin(SPDClval[j]),
                    f2DHistPartSPDV0M[lPAP[ih]+1]->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistPartSPDV0M[lPAP[ih]+1]->GetYaxis()->FindBin(V0Mval[i]),
                    antiparterr
                    );
                PAPYieldNorm_e[ih][k] = TMath::Sqrt(parterr_norm*parterr_norm + antiparterr*antiparterr);
                PAPYield_e[ih][k] = PAPYieldNorm_e[ih][k];
                doerror(PAPYield_forPi[ih][k], PAPYield_forPi_e[ih][k], Evt[k], Evt_e[k]);
                doerror3(PAPYieldNorm[ih][k], PAPYield_e[ih][k], Evt[k], Evt_e[k], Nch[k], Nch_e[k]);
                doerror4(PAPYieldNorm[ih][k], PAPYieldNorm_e[ih][k], Evt[k], Evt_e[k], Nch[k], Nch_e[k], (YieldMB[lPAP[ih]]+YieldMB[lPAP[ih]+1]), TMath::Sqrt(YieldMB_e[lPAP[ih]]*YieldMB_e[lPAP[ih]]+YieldMB_e[lPAP[ih]+1]*YieldMB_e[lPAP[ih]+1]));
                PAPYieldNorm[ih][k] /=  Evt[k];
                PAPYield_forPi[ih][k] = PAPYieldNorm[ih][k];
                PAPYieldNorm[ih][k] /=  Nch[k];
                PAPYield[ih][k] = PAPYieldNorm[ih][k];
                PAPYieldNorm[ih][k] /=  ((YieldPAPMB[lPAP[ih]]+YieldPAPMB[lPAP[ih]+1])/NchMB);

                //Ratio to Pions +-
                Double_t parterrPich, parterrPich_norm, antiparterrPich, antiparterrPich_norm;
                PAPYieldNormPich[ih][k]  =
                f2DHistPartSPDV0M[lPAP[ih]] -> IntegralAndError(
                    f2DHistPartSPDV0M[lPAP[ih]]->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistPartSPDV0M[lPAP[ih]]->GetXaxis()->FindBin(SPDClval[j]),
                    f2DHistPartSPDV0M[lPAP[ih]]->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistPartSPDV0M[lPAP[ih]]->GetYaxis()->FindBin(V0Mval[i]),
                    parterrPich_norm
                    );
                PAPYieldNormPich[ih][k] +=
                f2DHistPartSPDV0M[lPAP[ih]+1] -> IntegralAndError(
                    f2DHistPartSPDV0M[lPAP[ih]+1]->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistPartSPDV0M[lPAP[ih]+1]->GetXaxis()->FindBin(SPDClval[j]),
                    f2DHistPartSPDV0M[lPAP[ih]+1]->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistPartSPDV0M[lPAP[ih]+1]->GetYaxis()->FindBin(V0Mval[i]),
                    antiparterrPich
                    );
                PAPYieldNormPich_e[ih][k] = TMath::Sqrt(parterrPich_norm*parterrPich_norm + antiparterrPich*antiparterrPich);
                doerror4(PAPYieldNormPich[ih][k], PAPYieldNormPich_e[ih][k], Evt[k], Evt_e[k], PAPYield_forPi[0][k], PAPYield_forPi_e[0][k], (YieldMB[lPAP[ih]]+YieldMB[lPAP[ih]+1]), TMath::Sqrt(YieldMB_e[lPAP[ih]]*YieldMB_e[lPAP[ih]]+YieldMB_e[lPAP[ih]+1]*YieldMB_e[lPAP[ih]+1]));
                PAPYieldNormPich[ih][k] /=  Evt[k];
                PAPYieldNormPich[ih][k] /=  PAPYield_forPi[0][k];
                PAPYieldNormPich[ih][k] /=  ((YieldPAPMB[lPAP[ih]]+YieldPAPMB[lPAP[ih]+1])/(YieldPAPMB[lPAP[0]]+YieldPAPMB[lPAP[0]+1]));

                //Ratio to Pions 0
                Double_t parterrPi0, parterrPi0_norm, antiparterrPi0, antiparterrPi0_norm;
                PAPYieldNormPi0[ih][k]  =
                f2DHistPartSPDV0M[lPAP[ih]] -> IntegralAndError(
                    f2DHistPartSPDV0M[lPAP[ih]]->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistPartSPDV0M[lPAP[ih]]->GetXaxis()->FindBin(SPDClval[j]),
                    f2DHistPartSPDV0M[lPAP[ih]]->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistPartSPDV0M[lPAP[ih]]->GetYaxis()->FindBin(V0Mval[i]),
                    parterrPi0_norm
                    );
                PAPYieldNormPi0[ih][k] +=
                f2DHistPartSPDV0M[lPAP[ih]+1] -> IntegralAndError(
                    f2DHistPartSPDV0M[lPAP[ih]+1]->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistPartSPDV0M[lPAP[ih]+1]->GetXaxis()->FindBin(SPDClval[j]),
                    f2DHistPartSPDV0M[lPAP[ih]+1]->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistPartSPDV0M[lPAP[ih]+1]->GetYaxis()->FindBin(V0Mval[i]),
                    antiparterrPi0
                    );
                PAPYieldNormPi0_e[ih][k] = TMath::Sqrt(parterrPi0_norm*parterrPi0_norm + antiparterrPi0*antiparterrPi0);
                doerror4(PAPYieldNormPi0[ih][k], PAPYieldNormPi0_e[ih][k], Evt[k], Evt_e[k], 2*Yield_forPi[21][k], 2*Yield_forPi_e[21][k], (YieldMB[lPAP[ih]]+YieldMB[lPAP[ih]+1]), TMath::Sqrt(YieldMB_e[lPAP[ih]]*YieldMB_e[lPAP[ih]]+YieldMB_e[lPAP[ih]+1]*YieldMB_e[lPAP[ih]+1]));
                PAPYieldNormPi0[ih][k] /=  Evt[k];
                PAPYieldNormPi0[ih][k] /=  2*Yield_forPi[21][k];
                PAPYieldNormPi0[ih][k] /=  ((YieldPAPMB[lPAP[ih]]+YieldPAPMB[lPAP[ih]+1])/(2*YieldPAPMB[21]));
               // cout << PAPYieldNormPi0_e[ih][k] << "---"<<  Yield_forPi_e[21][k] << endl;


            }

            // Particle AvPt
            for(Int_t ih = 0; ih < nPart; ih++){
                AvPt[ih][k]  =
                f2DHistAvPtSPDV0M[ih] -> IntegralAndError(
                    f2DHistAvPtSPDV0M[ih]->GetXaxis()->FindBin(SPDClval[j+1]), f2DHistAvPtSPDV0M[ih]->GetXaxis()->FindBin(SPDClval[j]),
                    f2DHistAvPtSPDV0M[ih]->GetYaxis()->FindBin(V0Mval[i+1]), f2DHistAvPtSPDV0M[ih]->GetYaxis()->FindBin(V0Mval[i]),
                    AvPt_e[ih][k]
                    );
                AvPt[ih][k] /=  Npart[ih][k];
                AvPt_e[ih][k] /= Npart[ih][k];
            }

            // Particle + Antiparticle AvPt
            for(Int_t ih = 0; ih < 9; ih++){
                Double_t err;
                TH2D* hclone = (TH2D*)f2DHistAvPtSPDV0M[lPAP[ih]]->Clone("hclone");
                hclone->Add(f2DHistAvPtSPDV0M[lPAP[ih]+1]);
                PAPAvPt[ih][k] =
                hclone -> IntegralAndError(
                    hclone->GetXaxis()->FindBin(SPDClval[j+1]), hclone->GetXaxis()->FindBin(SPDClval[j]),
                    hclone->GetYaxis()->FindBin(V0Mval[i+1]), hclone->GetYaxis()->FindBin(V0Mval[i]),
                    err
                    );
                PAPAvPt[ih][k] /= (Npart[lPAP[ih]][k]+Npart[lPAP[ih]+1][k]);
                PAPAvPt_e[ih][k] = err/(Npart[lPAP[ih]][k]+Npart[lPAP[ih]+1][k]);
            }
        }
    }

    Double_t dNchNorm_e[nDiff];
    Double_t dNch_e[nDiff];
    Double_t dNMPINorm_e[nDiff];
    Double_t dNMPI_e[nDiff];
    Double_t dEffEnergy_e[nDiff];
    Double_t dLeadingE_e[nDiff];
    Double_t dLeadingENorm_e[nDiff];
    Double_t dYield_e[nPart][nDiff];
    Double_t dPAPYield_e[9][nDiff];
    Double_t dYieldNorm_e[nPart][nDiff];
    Double_t dPAPYieldNorm_e[9][nDiff];
    Double_t dPAPYieldNormPich_e[9][nDiff];
    Double_t dPAPYieldNormPi0_e[9][nDiff];
    Double_t dAvPt_e[nPart][nDiff];
    Double_t dPAPAvPt_e[9][nDiff];

    for(Int_t j = 0; j < nDiff; j++){

        dNchNorm_e[j] = (Double_t)NchNorm_e[j];
        dNch_e[j] = (Double_t)Nch_e[j];

        if (isPythia){
            dNMPI_e[j] = (Double_t)NMPI_e[j];
            dNMPINorm_e[j] = (Double_t)NMPINorm_e[j];
        }

        dEffEnergy_e[j] = (Double_t)EffEnergy_e[j];

        dLeadingE_e[j] = (Double_t)LeadingE_e[j];
        dLeadingENorm_e[j] = (Double_t)LeadingENorm_e[j];

        for(Int_t ih = 0; ih < nPart; ih++){
            dYield_e[ih][j] = (Double_t)Yield_e[ih][j];
            if(ih<9) dPAPYield_e[ih][j] = (Double_t)PAPYield_e[ih][j];
            dYieldNorm_e[ih][j] = (Double_t)YieldNorm_e[ih][j];
            if(ih<9) dPAPYieldNorm_e[ih][j] = (Double_t)PAPYieldNorm_e[ih][j];
            if(ih<9) dPAPYieldNormPich_e[ih][j] = (Double_t)PAPYieldNormPich_e[ih][j];
            if(ih<9) dPAPYieldNormPi0_e[ih][j] = (Double_t)PAPYieldNormPi0_e[ih][j];
            dAvPt_e[ih][j] = (Double_t)AvPt_e[ih][j];
            if(ih<9) dPAPAvPt_e[ih][j] = (Double_t)PAPAvPt_e[ih][j];
        }
    }

    // Vs Nch
    TGraphErrors* gNormRatioToNch_VsNchNorm[nPart];
    TGraphErrors* gNormPAPRatioToNch_VsNchNorm[9];
    TGraphErrors* gNormRatioToNch_VsNch[nPart];
    TGraphErrors* gNormPAPRatioToNch_VsNch[9];
    TGraphErrors* gRatioToNch_VsNch[nPart];
    TGraphErrors* gPAPRatioToNch_VsNch[9];
    TGraphErrors* gNormPAPRatioToPich_VsNchNorm[9];
    TGraphErrors* gNormPAPRatioToPi0_VsNchNorm[9];
    // percentile
    TGraphErrors* gNormRatioToNch_VsPerc[nPart];
    TGraphErrors* gNormPAPRatioToNch_VsPerc[9];
    // Vs Lead Energy
    TGraphErrors* gNormRatioToNch_VsLeadENorm[nPart];
    TGraphErrors* gNormPAPRatioToNch_VsLeadENorm[9];
    TGraphErrors* gNormRatioToNch_VsLeadE[nPart];
    TGraphErrors* gNormPAPRatioToNch_VsLeadE[9];
    TGraphErrors* gRatioToNch_VsLeadE[nPart];
    TGraphErrors* gPAPRatioToNch_VsLeadE[9];
    // Vs Eff Energy
    TGraphErrors* gNormRatioToNch_VsEE[nPart];
    TGraphErrors* gNormPAPRatioToNch_VsEE[9];
    TGraphErrors* gRatioToNch_VsEE[nPart];
    TGraphErrors* gPAPRatioToNch_VsEE[9];
    // Vs N MPIs
    TGraphErrors* gNormRatioToNch_VsNMPINorm[nPart];
    TGraphErrors* gNormPAPRatioToNch_VsNMPINorm[9];
    TGraphErrors* gNormRatioToNch_VsNMPI[nPart];
    TGraphErrors* gNormPAPRatioToNch_VsNMPI[9];
    TGraphErrors* gRatioToNch_VsNMPI[nPart];
    TGraphErrors* gPAPRatioToNch_VsNMPI[9];
    // Mean pT
    TGraphErrors* gAvPt_VsNch[nPart];
    TGraphErrors* gPAPAvPt_VsNch[9];
    TGraphErrors* gAvPt_VsNchNorm[nPart];
    TGraphErrors* gPAPAvPt_VsNchNorm[9];
    TGraphErrors* gAvPt_VsLeadE[nPart];
    TGraphErrors* gPAPAvPt_VsLeadE[9];
    TGraphErrors* gAvPt_VsLeadENorm[nPart];
    TGraphErrors* gPAPAvPt_VsLeadENorm[9];
    TGraphErrors* gAvPt_VsEE[nPart];
    TGraphErrors* gPAPAvPt_VsEE[9];
    TGraphErrors* gAvPt_VsNMPI[nPart];
    TGraphErrors* gPAPAvPt_VsNMPI[9];
    TGraphErrors* gAvPt_VsNMPINorm[nPart];
    TGraphErrors* gPAPAvPt_VsNMPINorm[9];

    TGraphErrors* gNchVsNMPI;
    if (isPythia){
        gNchVsNMPI = new TGraphErrors(nDiff,Nch,NMPI,Nch_e,NMPI_e);
        gNchVsNMPI->GetYaxis()->SetTitle("NMPI");
        gNchVsNMPI->GetXaxis()->SetTitle("Nch");
        gNchVsNMPI->SetTitle("");
        gNchVsNMPI->SetMarkerStyle(kFullCircle);
        gNchVsNMPI->SetName("NchVsNMPI");
        //
    }
    TGraphErrors* gLEVsNMPI;
    if (isPythia){
        gLEVsNMPI = new TGraphErrors(nDiff,LeadingE,NMPI,LeadingE_e,NMPI_e);
        gLEVsNMPI->GetYaxis()->SetTitle("NMPI");
        gLEVsNMPI->GetXaxis()->SetTitle("Leading Energy");
        gLEVsNMPI->SetTitle("");
        gLEVsNMPI->SetMarkerStyle(kFullCircle);
        gLEVsNMPI->SetName("LEVsNMPI");
    }
    TGraphErrors* gLEVsNch = new TGraphErrors(nDiff,Nch,LeadingE,Nch_e,LeadingE_e);
    gLEVsNch->GetYaxis()->SetTitle("Leading Energy");
    gLEVsNch->GetXaxis()->SetTitle("Nch");
    gLEVsNch->SetTitle("");
    gLEVsNch->SetMarkerStyle(kFullCircle);
    gLEVsNch->SetName("LEVsNch");
    TGraphErrors* gNormLEVsNch = new TGraphErrors(nDiff,NchNorm,LeadingENorm,NchNorm_e,LeadingENorm_e);
    gNormLEVsNch->GetYaxis()->SetTitle("E_{leading}/E_{leading}^{MB}");
    gNormLEVsNch->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
    gNormLEVsNch->SetTitle("");
    gNormLEVsNch->SetMarkerStyle(kFullCircle);
    gNormLEVsNch->SetName("NormLEVsNch");
    TGraphErrors* gEEVsNMPI;
    if (isPythia){
        gEEVsNMPI = new TGraphErrors(nDiff,EffEnergy,NMPI,EffEnergy_e,NMPI_e);
        gEEVsNMPI->GetYaxis()->SetTitle("NMPI");
        gEEVsNMPI->GetXaxis()->SetTitle("Effective Energy");
        gEEVsNMPI->SetTitle("");
        gEEVsNMPI->SetMarkerStyle(kFullCircle);
        gEEVsNMPI->SetName("EEVsNMPI");
    }

    for(Int_t ih = 0; ih < nPart; ih++){ //Single particle

        // Vs Nch
        gNormRatioToNch_VsNchNorm[ih] = new TGraphErrors(nDiff,NchNorm,YieldNorm[ih],dNchNorm_e,dYieldNorm_e[ih]);
        gNormRatioToNch_VsNchNorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{n_{ch}} / (#frac{%s}{n_{ch}})_{INEL>0}",lPartNames[ih].Data(),lPartNames[ih].Data()));
        gNormRatioToNch_VsNchNorm[ih]->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
        gNormRatioToNch_VsNchNorm[ih]->SetTitle("");
        gNormRatioToNch_VsNchNorm[ih]->SetMarkerStyle(kFullCircle);
        gNormRatioToNch_VsNchNorm[ih]->SetName(Form("NormRatioToNch_VsNchNorm_%s",lPartNames[ih].Data()));
        //
        gNormRatioToNch_VsNch[ih] = new TGraphErrors(nDiff,Nch,YieldNorm[ih],dNch_e,dYieldNorm_e[ih]);
        gNormRatioToNch_VsNch[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{n_{ch}} / (#frac{%s}{n_{ch}})_{INEL>0}",lPartNames[ih].Data(),lPartNames[ih].Data()));
        gNormRatioToNch_VsNch[ih]->GetXaxis()->SetTitle("n_{ch}");
        gNormRatioToNch_VsNch[ih]->SetTitle("");
        gNormRatioToNch_VsNch[ih]->SetMarkerStyle(kFullCircle);
        gNormRatioToNch_VsNch[ih]->SetName(Form("NormRatioToNch_VsNch_%s",lPartNames[ih].Data()));
        //
        gRatioToNch_VsNch[ih] = new TGraphErrors(nDiff,Nch,Yield[ih],dNch_e,dYield_e[ih]);
        gRatioToNch_VsNch[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{n_{ch}}",lPartNames[ih].Data()));
        gRatioToNch_VsNch[ih]->GetXaxis()->SetTitle("n_{ch}");
        gRatioToNch_VsNch[ih]->SetTitle("");
        gRatioToNch_VsNch[ih]->SetMarkerStyle(kFullCircle);
        gRatioToNch_VsNch[ih]->SetName(Form("RatioToNch_VsNch_%s",lPartNames[ih].Data()));

        // Vs Lead Energy
        gNormRatioToNch_VsLeadENorm[ih] = new TGraphErrors(nDiff,LeadingENorm,YieldNorm[ih],dLeadingENorm_e,dYieldNorm_e[ih]);
        gNormRatioToNch_VsLeadENorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{n_{ch}} / (#frac{%s}{n_{ch}})_{INEL>0}",lPartNames[ih].Data(),lPartNames[ih].Data()));
        gNormRatioToNch_VsLeadENorm[ih]->GetXaxis()->SetTitle("E_{leading}/E_{leading}^{MB}");
        gNormRatioToNch_VsLeadENorm[ih]->SetTitle("");
        gNormRatioToNch_VsLeadENorm[ih]->SetMarkerStyle(kFullCircle);
        gNormRatioToNch_VsLeadENorm[ih]->SetName(Form("NormRatioToNch_VsLeadENorm_%s",lPartNames[ih].Data()));
        //
        gNormRatioToNch_VsLeadE[ih] = new TGraphErrors(nDiff,LeadingE,YieldNorm[ih],dLeadingE_e,dYieldNorm_e[ih]);
        gNormRatioToNch_VsLeadE[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{n_{ch}} / (#frac{%s}{n_{ch}})_{INEL>0}",lPartNames[ih].Data(),lPartNames[ih].Data()));
        gNormRatioToNch_VsLeadE[ih]->GetXaxis()->SetTitle("E_{leading}");
        gNormRatioToNch_VsLeadE[ih]->SetTitle("");
        gNormRatioToNch_VsLeadE[ih]->SetMarkerStyle(kFullCircle);
        gNormRatioToNch_VsLeadE[ih]->SetName(Form("NormRatioToNch_VsLeadE_%s",lPartNames[ih].Data()));
        //
        gRatioToNch_VsLeadE[ih] = new TGraphErrors(nDiff,LeadingE,Yield[ih],dLeadingE_e,dYield_e[ih]);
        gRatioToNch_VsLeadE[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{n_{ch}}",lPartNames[ih].Data()));
        gRatioToNch_VsLeadE[ih]->GetXaxis()->SetTitle("E_{leading}");
        gRatioToNch_VsLeadE[ih]->SetTitle("");
        gRatioToNch_VsLeadE[ih]->SetMarkerStyle(kFullCircle);
        gRatioToNch_VsLeadE[ih]->SetName(Form("RatioToNch_VsLeadE_%s",lPartNames[ih].Data()));

        // Vs Eff Energy
        gNormRatioToNch_VsEE[ih] = new TGraphErrors(nDiff,EffEnergy,YieldNorm[ih],dEffEnergy_e,dYieldNorm_e[ih]);
        gNormRatioToNch_VsEE[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{n_{ch}} / (#frac{%s}{n_{ch}})_{INEL>0}",lPartNames[ih].Data(),lPartNames[ih].Data()));
        gNormRatioToNch_VsEE[ih]->GetXaxis()->SetTitle("E_{effective}");
        gNormRatioToNch_VsEE[ih]->SetTitle("");
        gNormRatioToNch_VsEE[ih]->SetMarkerStyle(kFullCircle);
        gNormRatioToNch_VsEE[ih]->SetName(Form("NormRatioToNch_VsEE_%s",lPartNames[ih].Data()));
        //
        gRatioToNch_VsEE[ih] = new TGraphErrors(nDiff,EffEnergy,Yield[ih],dEffEnergy_e,dYield_e[ih]);
        gRatioToNch_VsEE[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{n_{ch}}",lPartNames[ih].Data()));
        gRatioToNch_VsEE[ih]->GetXaxis()->SetTitle("E_{effective}");
        gRatioToNch_VsEE[ih]->SetTitle("");
        gRatioToNch_VsEE[ih]->SetMarkerStyle(kFullCircle);
        gRatioToNch_VsEE[ih]->SetName(Form("RatioToNch_VsEE_%s",lPartNames[ih].Data()));

        //Mean Pt
        gAvPt_VsNchNorm[ih] = new TGraphErrors(nDiff,NchNorm,AvPt[ih],dNchNorm_e,dAvPt_e[ih]);
        gAvPt_VsNchNorm[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
        gAvPt_VsNchNorm[ih]->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
        gAvPt_VsNchNorm[ih]->SetTitle("");
        gAvPt_VsNchNorm[ih]->SetMarkerStyle(kFullCircle);
        gAvPt_VsNchNorm[ih]->SetName(Form("AvPt_VsNchNorm_%s",lPartNames[ih].Data()));
        //
        gAvPt_VsNch[ih] = new TGraphErrors(nDiff,Nch,AvPt[ih],dNch_e,dAvPt_e[ih]);
        gAvPt_VsNch[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
        gAvPt_VsNch[ih]->GetXaxis()->SetTitle("n_{ch}");
        gAvPt_VsNch[ih]->SetTitle("");
        gAvPt_VsNch[ih]->SetMarkerStyle(kFullCircle);
        gAvPt_VsNch[ih]->SetName(Form("AvPt_VsNch_%s",lPartNames[ih].Data()));
        //
        gAvPt_VsLeadENorm[ih] = new TGraphErrors(nDiff,LeadingENorm,AvPt[ih],dLeadingENorm_e,dAvPt_e[ih]);
        gAvPt_VsLeadENorm[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
        gAvPt_VsLeadENorm[ih]->GetXaxis()->SetTitle("E_{leading}/E_{leading}^{MB}");
        gAvPt_VsLeadENorm[ih]->SetTitle("");
        gAvPt_VsLeadENorm[ih]->SetMarkerStyle(kFullCircle);
        gAvPt_VsLeadENorm[ih]->SetName(Form("AvPt_VsLeadENorm_%s",lPartNames[ih].Data()));
        //
        gAvPt_VsLeadE[ih] = new TGraphErrors(nDiff,LeadingE,AvPt[ih],dLeadingE_e,dAvPt_e[ih]);
        gAvPt_VsLeadE[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
        gAvPt_VsLeadE[ih]->GetXaxis()->SetTitle("E_{leading}");
        gAvPt_VsLeadE[ih]->SetTitle("");
        gAvPt_VsLeadE[ih]->SetMarkerStyle(kFullCircle);
        gAvPt_VsLeadE[ih]->SetName(Form("AvPt_VsLeadE_%s",lPartNames[ih].Data()));
        //
        gAvPt_VsEE[ih] = new TGraphErrors(nDiff,EffEnergy,AvPt[ih],dEffEnergy_e,dAvPt_e[ih]);
        gAvPt_VsEE[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
        gAvPt_VsEE[ih]->GetXaxis()->SetTitle("E_{effective}");
        gAvPt_VsEE[ih]->SetTitle("");
        gAvPt_VsEE[ih]->SetMarkerStyle(kFullCircle);
        gAvPt_VsEE[ih]->SetName(Form("AvPt_VsEE_%s",lPartNames[ih].Data()));

        // Vs NMPI
        if (isPythia) {
            gNormRatioToNch_VsNMPINorm[ih] = new TGraphErrors(nDiff,NMPINorm,YieldNorm[ih],dNMPINorm_e,dYieldNorm_e[ih]);
            gNormRatioToNch_VsNMPINorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{n_{ch}} / (#frac{%s}{n_{ch}})_{INEL>0}",lPartNames[ih].Data(),lPartNames[ih].Data()));
            gNormRatioToNch_VsNMPINorm[ih]->GetXaxis()->SetTitle("NMPI/NMPI^{MB}");
            gNormRatioToNch_VsNMPINorm[ih]->SetTitle("");
            gNormRatioToNch_VsNMPINorm[ih]->SetMarkerStyle(kFullCircle);
            gNormRatioToNch_VsNMPINorm[ih]->SetName(Form("NormRatioToNch_VsNMPINorm_%s",lPartNames[ih].Data()));
            //
            gNormRatioToNch_VsNMPI[ih] = new TGraphErrors(nDiff,NMPI,YieldNorm[ih],dNMPI_e,dYieldNorm_e[ih]);
            gNormRatioToNch_VsNMPI[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{n_{ch}} / (#frac{%s}{n_{ch}})_{INEL>0}",lPartNames[ih].Data(),lPartNames[ih].Data()));
            gNormRatioToNch_VsNMPI[ih]->GetXaxis()->SetTitle("NMPI");
            gNormRatioToNch_VsNMPI[ih]->SetTitle("");
            gNormRatioToNch_VsNMPI[ih]->SetMarkerStyle(kFullCircle);
            gNormRatioToNch_VsNMPI[ih]->SetName(Form("NormRatioToNch_VsNMPI_%s",lPartNames[ih].Data()));
            //
            gRatioToNch_VsNMPI[ih] = new TGraphErrors(nDiff,NMPI,Yield[ih],dNMPI_e,dYield_e[ih]);
            gRatioToNch_VsNMPI[ih]->GetYaxis()->SetTitle(Form("#frac{%s}{n_{ch}}",lPartNames[ih].Data()));
            gRatioToNch_VsNMPI[ih]->GetXaxis()->SetTitle("NMPI");
            gRatioToNch_VsNMPI[ih]->SetTitle("");
            gRatioToNch_VsNMPI[ih]->SetMarkerStyle(kFullCircle);
            gRatioToNch_VsNMPI[ih]->SetName(Form("RatioToNch_VsNMPI_%s",lPartNames[ih].Data()));
            //
            gAvPt_VsNMPINorm[ih] = new TGraphErrors(nDiff,NMPINorm,AvPt[ih],dNMPINorm_e,dAvPt_e[ih]);
            gAvPt_VsNMPINorm[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
            gAvPt_VsNMPINorm[ih]->GetXaxis()->SetTitle("NMPI/NMPI^{MB}");
            gAvPt_VsNMPINorm[ih]->SetTitle("");
            gAvPt_VsNMPINorm[ih]->SetMarkerStyle(kFullCircle);
            gAvPt_VsNMPINorm[ih]->SetName(Form("AvPt_VsNMPINorm_%s",lPartNames[ih].Data()));
            //
            gAvPt_VsNMPI[ih] = new TGraphErrors(nDiff,NMPI,AvPt[ih],dNMPI_e,dAvPt_e[ih]);
            gAvPt_VsNMPI[ih]->GetYaxis()->SetTitle(Form("%s #LT p_{T} #GT (GeV/c) }",lPartNames[ih].Data()));
            gAvPt_VsNMPI[ih]->GetXaxis()->SetTitle("NMPI");
            gAvPt_VsNMPI[ih]->SetTitle("");
            gAvPt_VsNMPI[ih]->SetMarkerStyle(kFullCircle);
            gAvPt_VsNMPI[ih]->SetName(Form("AvPt_VsNMPI_%s",lPartNames[ih].Data()));
        }
    }

    for(Int_t ih = 0; ih < 9; ih++){ //Part + AntiPart
        // Vs Nch
        gNormPAPRatioToNch_VsNchNorm[ih] = new TGraphErrors(nDiff,NchNorm,PAPYieldNorm[ih],dNchNorm_e,dPAPYieldNorm_e[ih]);
        gNormPAPRatioToNch_VsNchNorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}} / (#frac{%s+%s}{n_{ch}})_{INEL>0}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data(),lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gNormPAPRatioToNch_VsNchNorm[ih]->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
        gNormPAPRatioToNch_VsNchNorm[ih]->SetTitle("");
        gNormPAPRatioToNch_VsNchNorm[ih]->SetMarkerStyle(kFullCircle);
        gNormPAPRatioToNch_VsNchNorm[ih]->SetName(Form("NormRatioToNch_VsNchNorm_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        //
        gNormPAPRatioToNch_VsPerc[ih] = new TGraphErrors(nDiff,V0Mperc,PAPYieldNorm[ih],dNch_e,dPAPYieldNorm_e[ih]);
        gNormPAPRatioToNch_VsPerc[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}} / (#frac{%s+%s}{n_{ch}})_{INEL>0}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data(),lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gNormPAPRatioToNch_VsPerc[ih]->GetXaxis()->SetTitle("V0M Percentile");
        gNormPAPRatioToNch_VsPerc[ih]->SetTitle("");
        gNormPAPRatioToNch_VsPerc[ih]->SetMarkerStyle(kFullCircle);
        gNormPAPRatioToNch_VsPerc[ih]->SetName(Form("NormRatioToNch_VsPerc_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        //
        gNormPAPRatioToNch_VsNch[ih] = new TGraphErrors(nDiff,Nch,PAPYieldNorm[ih],dNch_e,dPAPYieldNorm_e[ih]);
        gNormPAPRatioToNch_VsNch[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}} / (#frac{%s+%s}{n_{ch}})_{INEL>0}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data(),lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gNormPAPRatioToNch_VsNch[ih]->GetXaxis()->SetTitle("n_{ch}");
        gNormPAPRatioToNch_VsNch[ih]->SetTitle("");
        gNormPAPRatioToNch_VsNch[ih]->SetMarkerStyle(kFullCircle);
        gNormPAPRatioToNch_VsNch[ih]->SetName(Form("NormRatioToNch_VsNch_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        //
        gPAPRatioToNch_VsNch[ih] = new TGraphErrors(nDiff,Nch,PAPYield[ih],dNch_e,dPAPYield_e[ih]);
        gPAPRatioToNch_VsNch[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gPAPRatioToNch_VsNch[ih]->GetXaxis()->SetTitle("n_{ch}");
        gPAPRatioToNch_VsNch[ih]->SetTitle("");
        gPAPRatioToNch_VsNch[ih]->SetMarkerStyle(kFullCircle);
        gPAPRatioToNch_VsNch[ih]->SetName(Form("RatioToNch_VsNch_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));

        // Vs Lead Energy
        gNormPAPRatioToNch_VsLeadENorm[ih] = new TGraphErrors(nDiff,LeadingENorm,PAPYieldNorm[ih],dLeadingENorm_e,dPAPYieldNorm_e[ih]);
        gNormPAPRatioToNch_VsLeadENorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}} / (#frac{%s+%s}{n_{ch}})_{INEL>0}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data(),lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gNormPAPRatioToNch_VsLeadENorm[ih]->GetXaxis()->SetTitle("E_{leading}/E_{leading}^{MB}");
        gNormPAPRatioToNch_VsLeadENorm[ih]->SetTitle("");
        gNormPAPRatioToNch_VsLeadENorm[ih]->SetMarkerStyle(kFullCircle);
        gNormPAPRatioToNch_VsLeadENorm[ih]->SetName(Form("NormRatioToNch_VsLeadENorm_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        //
        gNormPAPRatioToNch_VsLeadE[ih] = new TGraphErrors(nDiff,LeadingE,PAPYieldNorm[ih],dLeadingE_e,dPAPYieldNorm_e[ih]);
        gNormPAPRatioToNch_VsLeadE[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}} / (#frac{%s+%s}{n_{ch}})_{INEL>0}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data(),lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gNormPAPRatioToNch_VsLeadE[ih]->GetXaxis()->SetTitle("E_{leading}");
        gNormPAPRatioToNch_VsLeadE[ih]->SetTitle("");
        gNormPAPRatioToNch_VsLeadE[ih]->SetMarkerStyle(kFullCircle);
        gNormPAPRatioToNch_VsLeadE[ih]->SetName(Form("NormRatioToNch_VsLeadE_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        //
        gPAPRatioToNch_VsLeadE[ih] = new TGraphErrors(nDiff,LeadingE,PAPYield[ih],dLeadingE_e,dPAPYield_e[ih]);
        gPAPRatioToNch_VsLeadE[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gPAPRatioToNch_VsLeadE[ih]->GetXaxis()->SetTitle("E_{leading}");
        gPAPRatioToNch_VsLeadE[ih]->SetTitle("");
        gPAPRatioToNch_VsLeadE[ih]->SetMarkerStyle(kFullCircle);
        gPAPRatioToNch_VsLeadE[ih]->SetName(Form("RatioToNch_VsLeadE_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));

        // Vs Eff Energy
        gNormPAPRatioToNch_VsEE[ih] = new TGraphErrors(nDiff,EffEnergy,PAPYieldNorm[ih],dEffEnergy_e,dPAPYieldNorm_e[ih]);
        gNormPAPRatioToNch_VsEE[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}} / (#frac{%s+%s}{n_{ch}})_{INEL>0}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data(),lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gNormPAPRatioToNch_VsEE[ih]->GetXaxis()->SetTitle("E_{effective}");
        gNormPAPRatioToNch_VsEE[ih]->SetTitle("");
        gNormPAPRatioToNch_VsEE[ih]->SetMarkerStyle(kFullCircle);
        gNormPAPRatioToNch_VsEE[ih]->SetName(Form("NormRatioToNch_VsEE_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        //
        gPAPRatioToNch_VsEE[ih] = new TGraphErrors(nDiff,EffEnergy,PAPYield[ih],dEffEnergy_e,dPAPYield_e[ih]);
        gPAPRatioToNch_VsEE[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gPAPRatioToNch_VsEE[ih]->GetXaxis()->SetTitle("E_{effective}");
        gPAPRatioToNch_VsEE[ih]->SetTitle("");
        gPAPRatioToNch_VsEE[ih]->SetMarkerStyle(kFullCircle);
        gPAPRatioToNch_VsEE[ih]->SetName(Form("RatioToNch_VsEE_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));

        // Mean Pt
        gPAPAvPt_VsNchNorm[ih] = new TGraphErrors(nDiff,NchNorm,PAPAvPt[ih],dNchNorm_e,dPAPAvPt_e[ih]);
        gPAPAvPt_VsNchNorm[ih]->GetYaxis()->SetTitle(Form("%s+%s #LT p_{T} #GT (GeV/c) }",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gPAPAvPt_VsNchNorm[ih]->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
        gPAPAvPt_VsNchNorm[ih]->SetTitle("");
        gPAPAvPt_VsNchNorm[ih]->SetMarkerStyle(kFullCircle);
        gPAPAvPt_VsNchNorm[ih]->SetName(Form("AvPt_VsNchNorm_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        //
        gPAPAvPt_VsNch[ih] = new TGraphErrors(nDiff,Nch,PAPAvPt[ih],dNch_e,dPAPAvPt_e[ih]);
        gPAPAvPt_VsNch[ih]->GetYaxis()->SetTitle(Form("%s+%s #LT p_{T} #GT (GeV/c) }",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gPAPAvPt_VsNch[ih]->GetXaxis()->SetTitle("n_{ch}");
        gPAPAvPt_VsNch[ih]->SetTitle("");
        gPAPAvPt_VsNch[ih]->SetMarkerStyle(kFullCircle);
        gPAPAvPt_VsNch[ih]->SetName(Form("AvPt_VsNch_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        //
        gPAPAvPt_VsLeadENorm[ih] = new TGraphErrors(nDiff,LeadingENorm,PAPAvPt[ih],dLeadingENorm_e,dPAPAvPt_e[ih]);
        gPAPAvPt_VsLeadENorm[ih]->GetYaxis()->SetTitle(Form("%s+%s #LT p_{T} #GT (GeV/c) }",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gPAPAvPt_VsLeadENorm[ih]->GetXaxis()->SetTitle("E_{leading}/E_{leading}^{MB}");
        gPAPAvPt_VsLeadENorm[ih]->SetTitle("");
        gPAPAvPt_VsLeadENorm[ih]->SetMarkerStyle(kFullCircle);
        gPAPAvPt_VsLeadENorm[ih]->SetName(Form("AvPt_VsLeadENorm_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        //
        gPAPAvPt_VsLeadE[ih] = new TGraphErrors(nDiff,LeadingE,PAPAvPt[ih],dLeadingE_e,dPAPAvPt_e[ih]);
        gPAPAvPt_VsLeadE[ih]->GetYaxis()->SetTitle(Form("%s+%s #LT p_{T} #GT (GeV/c) }",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gPAPAvPt_VsLeadE[ih]->GetXaxis()->SetTitle("E_{leading}");
        gPAPAvPt_VsLeadE[ih]->SetTitle("");
        gPAPAvPt_VsLeadE[ih]->SetMarkerStyle(kFullCircle);
        gPAPAvPt_VsLeadE[ih]->SetName(Form("AvPt_VsLeadE_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        //
        gPAPAvPt_VsEE[ih] = new TGraphErrors(nDiff,EffEnergy,PAPAvPt[ih],dEffEnergy_e,dPAPAvPt_e[ih]);
        gPAPAvPt_VsEE[ih]->GetYaxis()->SetTitle(Form("%s+%s #LT p_{T} #GT (GeV/c) }",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gPAPAvPt_VsEE[ih]->GetXaxis()->SetTitle("E_{effective}");
        gPAPAvPt_VsEE[ih]->SetTitle("");
        gPAPAvPt_VsEE[ih]->SetMarkerStyle(kFullCircle);
        gPAPAvPt_VsEE[ih]->SetName(Form("AvPt_VsEE_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));

        // Vs NMPI
        if (isPythia) {
            gNormPAPRatioToNch_VsNMPINorm[ih] = new TGraphErrors(nDiff,NMPINorm,PAPYieldNorm[ih],dNMPINorm_e,dPAPYieldNorm_e[ih]);
            gNormPAPRatioToNch_VsNMPINorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}} / (#frac{%s+%s}{n_{ch}})_{INEL>0}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data(),lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
            gNormPAPRatioToNch_VsNMPINorm[ih]->GetXaxis()->SetTitle("NMPI/NMPI^{MB}");
            gNormPAPRatioToNch_VsNMPINorm[ih]->SetTitle("");
            gNormPAPRatioToNch_VsNMPINorm[ih]->SetMarkerStyle(kFullCircle);
            gNormPAPRatioToNch_VsNMPINorm[ih]->SetName(Form("NormRatioToNch_VsNMPINorm_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
            //
            gNormPAPRatioToNch_VsNMPI[ih] = new TGraphErrors(nDiff,NMPI,PAPYieldNorm[ih],dNMPI_e,dPAPYieldNorm_e[ih]);
            gNormPAPRatioToNch_VsNMPI[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}} / (#frac{%s+%s}{n_{ch}})_{INEL>0}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data(),lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
            gNormPAPRatioToNch_VsNMPI[ih]->GetXaxis()->SetTitle("NMPI");
            gNormPAPRatioToNch_VsNMPI[ih]->SetTitle("");
            gNormPAPRatioToNch_VsNMPI[ih]->SetMarkerStyle(kFullCircle);
            gNormPAPRatioToNch_VsNMPI[ih]->SetName(Form("NormRatioToNch_VsNMPI_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
            //
            gPAPRatioToNch_VsNMPI[ih] = new TGraphErrors(nDiff,NMPI,PAPYield[ih],dNMPI_e,dPAPYield_e[ih]);
            gPAPRatioToNch_VsNMPI[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{n_{ch}}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
            gPAPRatioToNch_VsNMPI[ih]->GetXaxis()->SetTitle("NMPI");
            gPAPRatioToNch_VsNMPI[ih]->SetTitle("");
            gPAPRatioToNch_VsNMPI[ih]->SetMarkerStyle(kFullCircle);
            gPAPRatioToNch_VsNMPI[ih]->SetName(Form("RatioToNch_VsNMPI_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
            //
            gPAPAvPt_VsNMPINorm[ih] = new TGraphErrors(nDiff,NMPINorm,PAPAvPt[ih],dNMPINorm_e,dPAPAvPt_e[ih]);
            gPAPAvPt_VsNMPINorm[ih]->GetYaxis()->SetTitle(Form("%s+%s #LT p_{T} #GT (GeV/c) }",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
            gPAPAvPt_VsNMPINorm[ih]->GetXaxis()->SetTitle("NMPI/NMPI^{MB}");
            gPAPAvPt_VsNMPINorm[ih]->SetTitle("");
            gPAPAvPt_VsNMPINorm[ih]->SetMarkerStyle(kFullCircle);
            gPAPAvPt_VsNMPINorm[ih]->SetName(Form("AvPt_VsNMPINorm_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
            //
            gPAPAvPt_VsNMPI[ih] = new TGraphErrors(nDiff,NMPI,PAPAvPt[ih],dNMPI_e,dPAPAvPt_e[ih]);
            gPAPAvPt_VsNMPI[ih]->GetYaxis()->SetTitle(Form("%s+%s #LT p_{T} #GT (GeV/c) }",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
            gPAPAvPt_VsNMPI[ih]->GetXaxis()->SetTitle("NMPI");
            gPAPAvPt_VsNMPI[ih]->SetTitle("");
            gPAPAvPt_VsNMPI[ih]->SetMarkerStyle(kFullCircle);
            gPAPAvPt_VsNMPI[ih]->SetName(Form("AvPt_VsNMPI_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        }

        //Ratio to Pions
        // Vs Nch
        gNormPAPRatioToPich_VsNchNorm[ih] = new TGraphErrors(nDiff,NchNorm,PAPYieldNormPich[ih],dNchNorm_e,dPAPYieldNormPich_e[ih]);
        gNormPAPRatioToPich_VsNchNorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{#pi^{#pm}} / (#frac{%s+%s}{#pi^{#pm}})_{INEL>0}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data(),lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gNormPAPRatioToPich_VsNchNorm[ih]->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
        gNormPAPRatioToPich_VsNchNorm[ih]->SetTitle("");
        gNormPAPRatioToPich_VsNchNorm[ih]->SetMarkerStyle(kFullCircle);
        gNormPAPRatioToPich_VsNchNorm[ih]->SetName(Form("NormRatioToPich_VsNchNorm_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        //
        gNormPAPRatioToPi0_VsNchNorm[ih] = new TGraphErrors(nDiff,NchNorm,PAPYieldNormPi0[ih],dNchNorm_e,dPAPYieldNormPi0_e[ih]);
        gNormPAPRatioToPi0_VsNchNorm[ih]->GetYaxis()->SetTitle(Form("#frac{%s+%s}{2#pi^{0}} / (#frac{%s+%s}{2#pi^{0}})_{INEL>0}",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data(),lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));
        gNormPAPRatioToPi0_VsNchNorm[ih]->GetXaxis()->SetTitle("n_{ch}/n_{ch}^{MB}");
        gNormPAPRatioToPi0_VsNchNorm[ih]->SetTitle("");
        gNormPAPRatioToPi0_VsNchNorm[ih]->SetMarkerStyle(kFullCircle);
        gNormPAPRatioToPi0_VsNchNorm[ih]->SetName(Form("NormRatioToPi0_VsNchNorm_%s%s",lPartNames[lPAP[ih]].Data(),lPartNames[lPAP[ih]+1].Data()));

    }

    TFile * write = new TFile(outputfile, option);

    TDirectoryFile *dir = new TDirectoryFile(directory.Data(),directory.Data());
    dir->cd();

    if (isPythia) gNchVsNMPI->Write();
    if (isPythia) gLEVsNMPI->Write();
    gLEVsNch->Write();
    gNormLEVsNch->Write();
    if (isPythia)gEEVsNMPI->Write();

    TDirectoryFile *lNormRatioToNch = new TDirectoryFile("NormRatioToNch","h/nch / (h/nch)_{MB}");
    lNormRatioToNch->cd();
    for(Int_t ih = 0; ih < nPart; ih++){
        gNormRatioToNch_VsNchNorm[ih]->Write();
        gNormRatioToNch_VsNch[ih]->Write();
        gNormRatioToNch_VsLeadENorm[ih]->Write();
        gNormRatioToNch_VsLeadE[ih]->Write();
        gNormRatioToNch_VsEE[ih]->Write();
        if (isPythia) {
            gNormRatioToNch_VsNMPINorm[ih]->Write();
            gNormRatioToNch_VsNMPI[ih]->Write();
        }
    }
    for(Int_t ih = 0; ih < 9; ih++){
        gNormPAPRatioToNch_VsNchNorm[ih]->Write();
        gNormPAPRatioToNch_VsNch[ih]->Write();
        gNormPAPRatioToNch_VsPerc[ih]->Write();
        gNormPAPRatioToNch_VsLeadENorm[ih]->Write();
        gNormPAPRatioToNch_VsLeadE[ih]->Write();
        gNormPAPRatioToNch_VsEE[ih]->Write();
        if (isPythia) {
            gNormPAPRatioToNch_VsNMPINorm[ih]->Write();
            gNormPAPRatioToNch_VsNMPI[ih]->Write();
        }
    }
    dir->cd();

    TDirectoryFile *lNormRatioToPions = new TDirectoryFile("NormRatioToPions","h/pi / (h/pi)_{MB}");
    lNormRatioToPions->cd();
    for(Int_t ih = 0; ih < 9; ih++){
        gNormPAPRatioToPich_VsNchNorm[ih]->Write();
        gNormPAPRatioToPi0_VsNchNorm[ih]->Write();
    }
    dir->cd();

    TDirectoryFile *lRatioToNch = new TDirectoryFile("RatioToNch","h/nch");
    lRatioToNch->cd();
    for(Int_t ih = 0; ih < nPart; ih++){
        gRatioToNch_VsNch[ih]->Write();
        gRatioToNch_VsLeadE[ih]->Write();
        gRatioToNch_VsEE[ih]->Write();
        if (isPythia) {
            gRatioToNch_VsNMPI[ih]->Write();
        }
    }
    for(Int_t ih = 0; ih < 9; ih++){
        gPAPRatioToNch_VsNch[ih]->Write();
        gPAPRatioToNch_VsLeadE[ih]->Write();
        gPAPRatioToNch_VsEE[ih]->Write();
        if (isPythia) {
            gPAPRatioToNch_VsNMPI[ih]->Write();
        }
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
    for(Int_t ih = 0; ih < 9; ih++){
        gPAPAvPt_VsNch[ih]->Write();
        gPAPAvPt_VsNchNorm[ih]->Write();
        gPAPAvPt_VsLeadE[ih]->Write();
        gPAPAvPt_VsLeadENorm[ih]->Write();
        gPAPAvPt_VsEE[ih]->Write();
        if (isPythia) {
            gPAPAvPt_VsNMPI[ih]->Write();
            gPAPAvPt_VsNMPINorm[ih]->Write();};
    }
    dir->cd();
    write->cd();
    write->Close();

    //Save in pdf important stuff
    TCanvas* cpdf = new TCanvas("cpdf",Form("cpdf%i",0),1200,800);
    if (directory.Contains("Standalone")){
        fHistEventCounter->Draw();
        cpdf->Print(Form("QA/QA_%s.pdf(",mcname.Data()),"pdf");
        fHistV0MMult->Draw();
        cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
        fHistMult05->Draw();
        cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
        fHistSPDClusters->Draw();
        cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
        f2DHistINELgt0SPDV0M->SetStats(0);
        f2DHistINELgt0SPDV0M->Draw("colz");
        cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
        f2DHistNchSPDV0M->SetStats(0);
        f2DHistNchSPDV0M->Draw("colz");
        cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
        if (isPythia) {
            fHistNMPI->Draw();
            cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
            f2DHistNMPISPDV0M->SetStats(0);
            f2DHistNMPISPDV0M->Draw("colz");
            cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
        }
        fHistLeadingE->Draw();
        cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
        f2DHistLeadingESPDV0M->SetStats(0);
        f2DHistLeadingESPDV0M->Draw("colz");
        cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
        fHistEffEnergy->Draw();
        cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
        f2DHistEffEnergySPDV0M->SetStats(0);
        f2DHistEffEnergySPDV0M->Draw("colz");
        cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
        hcalibV0M->Draw();
        cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
        hcalibSPDCl->Draw();
        cpdf->Print(Form("QA/QA_%s.pdf",mcname.Data()),"pdf");
    }
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
void doerror(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr){
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
void doerror3(Double_t A, Double_t &Aerr, Double_t B, Double_t Berr, Double_t C, Double_t Cerr){
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