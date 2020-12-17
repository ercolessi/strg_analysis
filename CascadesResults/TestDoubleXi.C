    
void TestDoubleXi(){

    cout<<"--------------- Open Real Data File --------------------"<<endl;
    TFile* file = new TFile("LHC15f_pass2.root", "READ");
    TTree* lTree      = (TTree*)file->Get("PWGLF_StrVsMult/fTreeCascade");

    TString fWhichParticle = "XiMinus";

    //Variable Definition=============================================
    Float_t fCentrality = 0.;
    Float_t fZDCCentrality = 0.;
    Bool_t  fMVPileupFlag = 0;
    Int_t   fClosestNonEmptyBC = 0;
    Float_t fZPApp = 0.;
    Float_t fZPCpp = 0.;
    Float_t fZNApp = 0.;
    Float_t fZNCpp = 0.;
    Int_t   fRun = 0;
    //Kinematic
    Float_t lPt, lRap, lPtMC, lNegEta, lPosEta, lBachEta;
    Float_t lPosPx, lPosPy, lPosPz;
    Float_t lNegPx, lNegPy, lNegPz;
    Float_t lBachPx, lBachPy, lBachPz;
    //Invariant Masses
    Float_t lInvariantMass; //wildcard
    Float_t lCompetingParticleMass = -1; //Competing Species rejection
    //Topological variable definitions
    Float_t lDcaV0Daughters; //1
    Float_t lDcaPosToPrimVertex,  lDcaNegToPrimVertex, lDcaBachToPrimVertex; //2, 3, 4
    Float_t lDCAzPosToPrimVertex, lDCAzNegToPrimVertex, lDCAzBachToPrimVertex; // 12, 13, 14;
    Float_t lDCABachToBaryon;//15
    Float_t lDCAxyCascToPV, lDCAzCascToPV, lDCACascToPV;
    //Cosine of Pointing Angle variable
    Float_t lV0CosinePointingAngle, lCascCosinePointingAngle; //5, 6
    Float_t lV0CosinePointingAngleSpecial;
    Float_t lBBCosPA;
    //Decay Radius and distance over total momentum
    Float_t lV0Radius, lCascRadius, lDistOverTotMom; //7, 8
    Float_t lV0Mass, lDcaCascDaughters; //9, 10
    Float_t lDcaV0ToPV; //11
    //Least Number of TPC Clusters
    Int_t lLeastNbrClusters;
    Float_t lMinTrackLength;
    Float_t lMaxChi2PerCluster;
    //TPC dE/dx acquired with AliPIDResponse class
    Float_t lNSigmasPosProton,lNSigmasNegProton,lNSigmasPosPion,lNSigmasNegPion,
            lNSigmasBachPion, lNSigmasBachKaon;
    Float_t lPosdEdx, lNegdEdx, lBachdEdx;
    Float_t lPosInnerP, lNegInnerP, lBachInnerP;
    //ITS
    ULong64_t lNegTrackStatus, lPosTrackStatus, lBachTrackStatus;
    //TOF
    Float_t lNegTOFExpTDiff, lPosTOFExpTDiff, lBachTOFExpTDiff;
    Float_t lNegTOFSignal, lPosTOFSignal, lBachTOFSignal;
    //Charge
    Int_t lCharge = 0;
    ULong64_t fEventNumber = 0;
    //Multiplicity Variable
    Float_t lMultiplicity = -1.;
    Float_t lCascZDCPercentile = -1.;
    Bool_t ITSrefitLowPtBothLegs;
    Bool_t TOFmatchHighPtBothLegs;
    Bool_t TOFmatchHighPtOneLeg;
    const int kITSrefit = 4;
    //================================================================

    //Linking to Tree=================================================
    //--- Base Variables ----------------------------------------------
    lTree->SetBranchAddress("fTreeCascVarCharge"  ,&lCharge );
    lTree->SetBranchAddress("fTreeCascVarPosEta"  ,&lPosEta );
    lTree->SetBranchAddress("fTreeCascVarNegEta"  ,&lNegEta );
    lTree->SetBranchAddress("fTreeCascVarBachEta" ,&lBachEta);
    lTree->SetBranchAddress("fTreeCascVarPt",&lPt);
    if ( fWhichParticle == "XiMinus"      )  lTree->SetBranchAddress("fTreeCascVarMassAsXi"   ,&lInvariantMass);
    if ( fWhichParticle == "XiPlus"       )  lTree->SetBranchAddress("fTreeCascVarMassAsXi"   ,&lInvariantMass);
    if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" ){
        lTree->SetBranchAddress("fTreeCascVarMassAsOmega",&lInvariantMass);
        lTree->SetBranchAddress("fTreeCascVarMassAsXi"   ,&lCompetingParticleMass); //For Competing Rejection
    }
    if ( fWhichParticle == "XiMinus"    || fWhichParticle == "XiPlus"    )
        lTree->SetBranchAddress("fTreeCascVarRapXi"   ,&lRap);
    if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" )
        lTree->SetBranchAddress("fTreeCascVarRapOmega",&lRap);
    lTree->SetBranchAddress("fTreeCascVarLeastNbrClusters",&lLeastNbrClusters);
    lTree->SetBranchAddress("fTreeCascVarPosPx"     , &lPosPx);
    lTree->SetBranchAddress("fTreeCascVarPosPy"     , &lPosPy);
    lTree->SetBranchAddress("fTreeCascVarPosPz"     , &lPosPz);
    lTree->SetBranchAddress("fTreeCascVarNegPx"     , &lNegPx);
    lTree->SetBranchAddress("fTreeCascVarNegPy"     , &lNegPy);
    lTree->SetBranchAddress("fTreeCascVarNegPz"     , &lNegPz);
    lTree->SetBranchAddress("fTreeCascVarBachPx"    , &lBachPx);
    lTree->SetBranchAddress("fTreeCascVarBachPy"    , &lBachPy);
    lTree->SetBranchAddress("fTreeCascVarBachPz"    , &lBachPz);
    lTree->SetBranchAddress("fTreeCascVarPosInnerP",  &lPosInnerP);
    lTree->SetBranchAddress("fTreeCascVarNegInnerP",  &lNegInnerP);
    lTree->SetBranchAddress("fTreeCascVarBachInnerP", &lBachInnerP);
    //--- TPC Variables -----------------------------------------------
    lTree->SetBranchAddress("fTreeCascVarPosNSigmaProton",&lNSigmasPosProton);
    lTree->SetBranchAddress("fTreeCascVarNegNSigmaProton",&lNSigmasNegProton);
    lTree->SetBranchAddress("fTreeCascVarPosNSigmaPion",&lNSigmasPosPion);
    lTree->SetBranchAddress("fTreeCascVarNegNSigmaPion",&lNSigmasNegPion);
    lTree->SetBranchAddress("fTreeCascVarBachNSigmaPion",&lNSigmasBachPion);
    lTree->SetBranchAddress("fTreeCascVarBachNSigmaKaon",&lNSigmasBachKaon);
    lTree->SetBranchAddress("fTreeCascVarPosdEdx",  &lPosdEdx);
    lTree->SetBranchAddress("fTreeCascVarNegdEdx",  &lNegdEdx);
    lTree->SetBranchAddress("fTreeCascVarBachdEdx", &lBachdEdx);
    lTree->SetBranchAddress("fTreeCascVarMinTrackLength", &lMinTrackLength);
    lTree->SetBranchAddress("fTreeCascVarMaxChi2PerCluster", &lMaxChi2PerCluster);
    //--- Topological selection variables -----------------------------
    lTree->SetBranchAddress("fTreeCascVarV0Radius",&lV0Radius); //1
    lTree->SetBranchAddress("fTreeCascVarCascRadius",&lCascRadius); //2
    lTree->SetBranchAddress("fTreeCascVarV0Mass",&lV0Mass); //3
    lTree->SetBranchAddress("fTreeCascVarV0CosPointingAngle",&lV0CosinePointingAngle); //4
    lTree->SetBranchAddress("fTreeCascVarV0CosPointingAngleSpecial",&lV0CosinePointingAngleSpecial);
    lTree->SetBranchAddress("fTreeCascVarCascCosPointingAngle",&lCascCosinePointingAngle); //5
    lTree->SetBranchAddress("fTreeCascVarDCANegToPrimVtx",&lDcaNegToPrimVertex); //6
    lTree->SetBranchAddress("fTreeCascVarDCAPosToPrimVtx",&lDcaPosToPrimVertex); //7
    lTree->SetBranchAddress("fTreeCascVarDCABachToPrimVtx",&lDcaBachToPrimVertex); //8
    lTree->SetBranchAddress("fTreeCascVarDCAV0ToPrimVtx",&lDcaV0ToPV); //9
    lTree->SetBranchAddress("fTreeCascVarCascDCAtoPVxy", &lDCAxyCascToPV);
    lTree->SetBranchAddress("fTreeCascVarCascDCAtoPVz", &lDCAzCascToPV);
    lTree->SetBranchAddress("fTreeCascVarDCAV0Daughters",&lDcaV0Daughters); //10
    lTree->SetBranchAddress("fTreeCascVarDCACascDaughters",&lDcaCascDaughters); //11
    lTree->SetBranchAddress("fTreeCascVarDistOverTotMom",&lDistOverTotMom); //11
    lTree->SetBranchAddress("fTreeCascVarNegDCAz",  &lDCAzNegToPrimVertex);
    lTree->SetBranchAddress("fTreeCascVarPosDCAz",  &lDCAzPosToPrimVertex);
    lTree->SetBranchAddress("fTreeCascVarBachDCAz", &lDCAzBachToPrimVertex);
    lTree->SetBranchAddress("fTreeCascVarDCABachToBaryon", &lDCABachToBaryon);
    lTree->SetBranchAddress("fTreeCascVarWrongCosPA", &lBBCosPA);
    //--- ITS flag -----------------------------------------------------
    lTree->SetBranchAddress("fTreeCascVarPosTrackStatus", &lPosTrackStatus);
    lTree->SetBranchAddress("fTreeCascVarNegTrackStatus", &lNegTrackStatus);
    lTree->SetBranchAddress("fTreeCascVarBachTrackStatus", &lBachTrackStatus);
    //--- TOF info -----------------------------------------------------
    lTree->SetBranchAddress("fTreeCascVarNegTOFExpTDiff",  &lNegTOFExpTDiff);
    lTree->SetBranchAddress("fTreeCascVarPosTOFExpTDiff",  &lPosTOFExpTDiff);
    lTree->SetBranchAddress("fTreeCascVarBachTOFExpTDiff", &lBachTOFExpTDiff);
    lTree->SetBranchAddress("fTreeCascVarNegTOFSignal",  &lNegTOFSignal);
    lTree->SetBranchAddress("fTreeCascVarPosTOFSignal",  &lPosTOFSignal);
    lTree->SetBranchAddress("fTreeCascVarBachTOFSignal", &lBachTOFSignal);
    //--- Multiplicity Variable ----------------------------------------
    lTree->SetBranchAddress("fTreeCascVarCentrality",&fCentrality);
    lTree->SetBranchAddress("fTreeCascVarZPApp", &fZPApp);
    lTree->SetBranchAddress("fTreeCascVarZPCpp", &fZPCpp);
    lTree->SetBranchAddress("fTreeCascVarZNApp", &fZNApp);
    lTree->SetBranchAddress("fTreeCascVarZNCpp", &fZNCpp);
    lTree->SetBranchAddress("fTreeCascVarRun", &fRun);
    lTree->SetBranchAddress("fTreeCascVarEventNumber", &fEventNumber);
    //--- MV pileup flag -----------------------------------------------
    lTree->SetBranchAddress("fTreeCascVarMVPileupFlag", &fMVPileupFlag);
    lTree->SetBranchAddress("fTreeCascVarClosestNonEmptyBC", &fClosestNonEmptyBC);
    //================================================================

    Long_t lNCandidates = lTree->GetEntries();
    Long_t lOneTenthOfNCandidates = ((double)(lNCandidates) / 10. );
    ULong64_t evt = 0;
    double mass = 0.;
    int count = 1;
    TH2D* hInvMass = new TH2D("hInvMass", " ", 200, 1.24,1.4, 200, 1.24,1.4);
    
    cout<<"--------------- Real Data File ------------------"<<endl;
    for(Long_t icand = 0;icand<lNCandidates;icand++){
        lTree->GetEntry(icand);

        //Now check validity
        if( lRap<0.5 && lRap>(-0.5) &&
            (
            TMath::Abs(lNegEta)       < 0.8       &&
            TMath::Abs(lPosEta)       < 0.8       &&
            TMath::Abs(lBachEta)      < 0.8                    
		    ))
            { // Start Entry Loop
                //hEvtNumber->AddBinContent(fEventNumber,1);
                if (fEventNumber!= 0 && fEventNumber == evt) {
                    //hEvtNumber->AddBinContent(fEventNumber,1);
                    count++;
                    hInvMass->Fill(lInvariantMass, mass); //fill with specific inv mass                        
                }  
                else count = 1;

                evt = fEventNumber;
                mass = lInvariantMass;                         
        }
    }
    cout<<"--------------- Loop Completed -------------------------"<<endl;
    cout<<endl;    

    //TH2
    TCanvas* c = new TCanvas();
    hInvMass->Draw("LEGO2Z");
	hInvMass->GetXaxis()->SetRangeUser(1.26,1.4);
	hInvMass->GetYaxis()->SetRangeUser(1.26,1.4);
    hInvMass->GetXaxis()->SetTitleOffset(1.5);
    hInvMass->GetYaxis()->SetTitleOffset(1.5);
    hInvMass->GetXaxis()->SetTitle("M_{inv} (GeV/(c^{2})");
    hInvMass->GetYaxis()->SetTitle("M_{inv} (GeV/(c^{2})");
    //
    TCanvas* b = new TCanvas();
    hInvMass->Draw("colz");
    hInvMass->GetXaxis()->SetRangeUser(1.26,1.4);
    hInvMass->GetYaxis()->SetRangeUser(1.26,1.4);
    hInvMass->GetXaxis()->SetTitleOffset(1.5);
    hInvMass->GetYaxis()->SetTitleOffset(1.5);
    hInvMass->GetXaxis()->SetTitle("M_{inv} (GeV/(c^{2})");
    hInvMass->GetYaxis()->SetTitle("M_{inv} (GeV/(c^{2})");

    //Projections step-1
    TH1F* projx = (TH1F*)hInvMass->ProjectionX("projx",0,-1);
    new TCanvas;
    projx->Draw();
    TH1F* projy = (TH1F*)hInvMass->ProjectionY("projy",0,-1);
    new TCanvas;
    projy->Draw();

    //Fit projections
    //
    //X
    TF1* fgausPt_px = new TF1("fgausPt_px","[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]", 1.322-0.03, 1.322+0.03 );
    fgausPt_px->SetParameter(1,1.322);
    fgausPt_px->SetParameter(2,0.0025);
    fgausPt_px->SetParLimits(2,0.001,0.01);
    //
    projx->Fit("fgausPt_px","R");
    Double_t PeakPosition_px = fgausPt_px->GetParameter(1);
    Double_t PeakWidth_px = TMath::Abs(fgausPt_px->GetParameter(2));
    Int_t uplim_px = hInvMass->GetXaxis()->FindBin(PeakPosition_px+3*PeakWidth_px);
    Int_t downlim_px = hInvMass->GetXaxis()->FindBin(PeakPosition_px-3*PeakWidth_px); 
    //
    //Y
    TF1* fgausPt_py = new TF1("fgausPt_py","[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]", 1.322-0.03, 1.322+0.03 );
    fgausPt_py->SetParameter(1,1.322);
    fgausPt_py->SetParameter(2,0.0025);
    fgausPt_py->SetParLimits(2,0.001,0.01);
    //
    projy->Fit("fgausPt_py","R");
    Double_t PeakPosition_py = fgausPt_py->GetParameter(1);
    Double_t PeakWidth_py = TMath::Abs(fgausPt_py->GetParameter(2));
    Int_t uplim_py = hInvMass->GetYaxis()->FindBin(PeakPosition_py+3*PeakWidth_py);
    Int_t downlim_py = hInvMass->GetYaxis()->FindBin(PeakPosition_py-3*PeakWidth_py);
    //
    //New projections less bkg
    TH1F* ProjX_cut = (TH1F*)hInvMass->ProjectionX("ProjX_cut", downlim_py, uplim_py);
    TH1F* ProjY_cut = (TH1F*)hInvMass->ProjectionY("ProjY_cut", downlim_px, uplim_px);
    //
    new TCanvas;
    ProjX_cut->Draw();
    new TCanvas;
    ProjY_cut->Draw();

    //Stima veloce bkg 
    //
    TH1F* ProjX_cut_bkgl = (TH1F*)hInvMass->ProjectionX("ProjX_cut_bkgl", 0, downlim_py);
    TH1F* ProjX_cut_bkgr = (TH1F*)hInvMass->ProjectionX("ProjX_cut_bkgr", uplim_py, -1);
    TH1F* ProjY_cut_bkgl = (TH1F*)hInvMass->ProjectionY("ProjY_cut_bkgl",0., downlim_px);
    TH1F* ProjY_cut_bkgr = (TH1F*)hInvMass->ProjectionY("ProjY_cut_bkgr",uplim_py, -1);
    //X
    ProjX_cut_bkgl->Add(ProjX_cut_bkgr);
    ProjX_cut_bkgl->Scale(1./(ProjX_cut_bkgl->Integral(0,-1)));
    ProjX_cut->Scale(1./(ProjX_cut->Integral(0,-1)));
    //
    new TCanvas;
    ProjX_cut_bkgl->SetLineColor(kBlue);
    ProjX_cut_bkgr->SetLineColor(kBlue);
    ProjX_cut_bkgl->SetLineStyle(2);
    ProjX_cut_bkgr->SetLineStyle(2);
    ProjX_cut->SetLineColor(kRed);
    ProjX_cut->Draw();
    ProjX_cut_bkgl->Draw("SAME");
    //ProjX_cut_bkgr->Draw("SAME");
    //Y
    ProjY_cut_bkgl->Add(ProjY_cut_bkgr);
    ProjY_cut_bkgl->Scale(1./(ProjY_cut_bkgl->Integral(0,-1)));
    ProjY_cut->Scale(1./(ProjY_cut->Integral(0,-1)));
    //
    new TCanvas;
    ProjY_cut_bkgl->SetLineColor(kBlue);
    ProjY_cut_bkgr->SetLineColor(kBlue);
    ProjY_cut_bkgl->SetLineStyle(2);
    ProjY_cut_bkgr->SetLineStyle(2);
    ProjY_cut->SetLineColor(kRed);
    ProjY_cut->Draw();
    ProjY_cut_bkgl->Draw("SAME");
    //ProjY_cut_bkgr->Draw("SAME");

	new TCanvas;
    TH1F* hSBX = (TH1F*)ProjX_cut->Clone("hSBX");
    for(int i = 1; i<=ProjX_cut->GetNbinsX();i++){
        double entry = ProjX_cut->GetBinContent(i)-ProjX_cut_bkgl->GetBinContent(i);
        if (entry > 0) hSBX->SetBinContent(i,entry);
        else hSBX->SetBinContent(i,0.);
    }
    hSBX->Draw();
    new TCanvas;
    TH1F* hSBY = (TH1F*)ProjY_cut->Clone("hSBY");
    for(int i = 1; i<=ProjY_cut->GetNbinsX();i++){
        double entry = ProjY_cut->GetBinContent(i)-ProjY_cut_bkgl->GetBinContent(i);
        if (entry > 0) hSBY->SetBinContent(i,entry);
        else hSBY->SetBinContent(i,0.);
    }
    hSBY->Draw();




}
