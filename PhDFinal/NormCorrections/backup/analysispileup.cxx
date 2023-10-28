void analysispileup(){

    TFile *fileMC = new TFile("AnalysisResults_Bis.root", "READ");

    TTree* lTreeEventMC;
    Float_t fCentrality_SPDClusters = 0.;
    Float_t fCentrality_V0M = 0.;
    Bool_t fEvSel_AllSelections = kFALSE;
    Bool_t fEvSel_INELgtZEROtrue = kFALSE;
    Bool_t fEvSel_zVtxZMC = kFALSE;
    Bool_t fEvSel_NotPileupSPDInMultBins = kFALSE;
    Bool_t fEvSel_Pileuptrue = kFALSE; // FIXME
    Bool_t fEvSel_NoInconsSPDandTrackVrtx = kFALSE;
    Bool_t fEvSel_AcceptedVertexPosition = kFALSE;
    Int_t fSPDtracklets = 0.;
    Int_t fNtrackstrue = 0.;
    Float_t fAmplitudeV0A = 0.;
    Float_t fAmplitudeV0C = 0.;
    Float_t ZVtxReco = 0.;
    Float_t ZVtxMC = 0.;
    Int_t fRun = 0.;
 
    lTreeEventMC = (TTree*)fileMC->Get("PWGLF_StrVsMult_MC/fTreeEvent");
    
    //Get from Tree
    lTreeEventMC->SetBranchAddress("fRun",&fRun);
    lTreeEventMC->SetBranchAddress("fEvSel_zVtxZMC", &fEvSel_zVtxZMC);
    lTreeEventMC->SetBranchAddress("fEvSel_AllSelections", &fEvSel_AllSelections);
    lTreeEventMC->SetBranchAddress("fEvSel_INELgtZEROtrue", &fEvSel_INELgtZEROtrue);
    lTreeEventMC->SetBranchAddress("fEvSel_AcceptedVertexPosition", &fEvSel_AcceptedVertexPosition);
    lTreeEventMC->SetBranchAddress("fEvSel_NotPileupSPDInMultBins", &fEvSel_NotPileupSPDInMultBins);
    lTreeEventMC->SetBranchAddress("fEvSel_Pileuptrue", &fEvSel_Pileuptrue); // FIXME
    lTreeEventMC->SetBranchAddress("fEvSel_NoInconsSPDandTrackVrtx", &fEvSel_NoInconsSPDandTrackVrtx); // FIXME
    lTreeEventMC->SetBranchAddress("fCentrality_SPDClusters", &fCentrality_SPDClusters);
    lTreeEventMC->SetBranchAddress("fCentrality_V0M", &fCentrality_V0M);
    lTreeEventMC->SetBranchAddress("fSPDtracklets", &fSPDtracklets);
    lTreeEventMC->SetBranchAddress("fNtrackstrue", &fNtrackstrue);
    lTreeEventMC->SetBranchAddress("fAmplitudeV0A", &fAmplitudeV0A);
    lTreeEventMC->SetBranchAddress("fAmplitudeV0C", &fAmplitudeV0C);
    //lTreeEventMC->SetBranchAddress("ZVtxReco", &ZVtxReco);
    //lTreeEventMC->SetBranchAddress("ZVtxMC", &ZVtxMC);

    TH1F * hAmplitudeV0M[5]; 
    TH1F *hNtrackstrue[5];
    TH1F *hSPDtracklets[5];
    TH1F* hZVtxReco[5];
    TH1F* hZVtxMC[5];

    for (Int_t i = 0; i < 5; i++) {
        hAmplitudeV0M[i] = new TH1F(Form("hAmplitudeV0M_%d", i), "AmplitudeV0M", 500, 0, 1500);
        hNtrackstrue[i] = new TH1F(Form("hNtrackstrue_%d", i), "Ntrackstrue", 100, 0, 100);
        hSPDtracklets[i] = new TH1F(Form("hSPDtracklets_%d", i), "SPDtracklets", 100, 0, 100);
        hZVtxReco[i] = new TH1F(Form("hZVtxReco_%d", i), "ZVtxReco", 600, -30, 30);
        hZVtxMC[i] = new TH1F(Form("hZVtxMC_%d", i), "ZVtxMC", 600, -30, 30);
        hAmplitudeV0M[i]->SetStats(0);
        hNtrackstrue[i]->SetStats(0);
        hSPDtracklets[i]->SetStats(0);
       // hZVtxReco[i]->SetStats(0);
       // hZVtxMC[i]->SetStats(0);
    }

    TH1I * hIsPileup = new TH1I("hIsPileup", "IsPileup", 2, 0, 2);
    
   
    cout<<" \nWill now loop over events, please wait...\n"<<endl;
    for(Long_t iEv = 0; iEv<lTreeEventMC->GetEntries(); iEv++) {

        lTreeEventMC->GetEntry(iEv);
        if( iEv % ( lTreeEventMC->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEventMC->GetEntries()<<endl;

        if (fCentrality_V0M < 0 || fCentrality_V0M > 1) continue;
        //if (fNtrackstrue > 5) continue;

        if (!fEvSel_Pileuptrue && fEvSel_INELgtZEROtrue && fEvSel_zVtxZMC)
        {
            hAmplitudeV0M[0]->Fill(fAmplitudeV0A + fAmplitudeV0C);
            hNtrackstrue[0]->Fill(fNtrackstrue);
            hSPDtracklets[0]->Fill(fSPDtracklets);
           // hZVtxReco[0]->Fill(ZVtxReco);
           // hZVtxMC[0]->Fill(ZVtxMC);
        }
        if (!fEvSel_Pileuptrue && fEvSel_AllSelections)
        {
            hAmplitudeV0M[1]->Fill(fAmplitudeV0A + fAmplitudeV0C);
            hNtrackstrue[1]->Fill(fNtrackstrue);
            hSPDtracklets[1]->Fill(fSPDtracklets);
           // hZVtxReco[1]->Fill(ZVtxReco);
           // hZVtxMC[1]->Fill(ZVtxMC);
        }
        if (!fEvSel_Pileuptrue && fEvSel_INELgtZEROtrue && fEvSel_AcceptedVertexPosition)
        {
            hAmplitudeV0M[2]->Fill(fAmplitudeV0A + fAmplitudeV0C);
            hNtrackstrue[2]->Fill(fNtrackstrue);
            hSPDtracklets[2]->Fill(fSPDtracklets);
           // hZVtxReco[2]->Fill(ZVtxReco);
           // hZVtxMC[2]->Fill(ZVtxMC);
        }
        if (fEvSel_Pileuptrue) hIsPileup->Fill(1);
        if (!fEvSel_Pileuptrue) hIsPileup->Fill(0);

     //   hZVtxReco[3]->Fill(ZVtxReco);
     //   hZVtxMC[3]->Fill(ZVtxMC);
    }

   /* TCanvas *c = new TCanvas("c", "V0M_0-1", 1300, 900);
    hZVtxReco[3]->SetLineColor(kMagenta);
    hZVtxMC[3]->SetLineColor(kBlack);
    hZVtxReco[3]->Draw();
    hZVtxMC[3]->Draw("same");

    TLegend *leg1 = new TLegend(0.4, 0.7, 0.9, 0.9);
    leg1->SetHeader("All events V0M[0-1]");
    leg1->AddEntry(hZVtxReco[3], "ZVtxReco");
    leg1->AddEntry(hZVtxMC[3], "ZVtxMC");
    leg1->Draw();*/

    TCanvas *c1 = new TCanvas("c1", "V0M_0-1", 1800, 900);
    c1->Divide(3, 2);

    c1->cd(1);
    hIsPileup->Draw();

    c1->cd(2);
    //c1->cd(2)->SetLogy();
    hAmplitudeV0M[0]->SetLineColor(kRed);
    hAmplitudeV0M[1]->SetLineColor(kBlue);
    hAmplitudeV0M[1]->SetLineStyle(7);
    hAmplitudeV0M[2]->SetLineColor(kGreen);
    hAmplitudeV0M[2]->SetLineStyle(9);
    hAmplitudeV0M[0]->Draw();
    hAmplitudeV0M[1]->Draw("same");
    hAmplitudeV0M[2]->Draw("same");

    TLegend *leg = new TLegend(0.4, 0.7, 0.9, 0.9);
    leg->AddEntry(hAmplitudeV0M[0], "No Pileup(true) + INEL>0(true) + ZVtxtrue<10cm", "l");
    leg->AddEntry(hAmplitudeV0M[1], "No Pileup(true) + INEL>0(data) + all data sel", "l");
    leg->AddEntry(hAmplitudeV0M[2], "No Pileup(true) + INEL>0(true) + ZVtx(data)<10cm", "l");
    leg->Draw();

    c1->cd(4);
    //c1->cd(3)->SetLogy();
    hNtrackstrue[0]->SetLineColor(kRed);
    hNtrackstrue[1]->SetLineColor(kBlue);
    hNtrackstrue[1]->SetLineStyle(7);
    hNtrackstrue[2]->SetLineColor(kGreen);
    hNtrackstrue[2]->SetLineStyle(9);
    hNtrackstrue[0]->Draw();
    hNtrackstrue[1]->Draw("same");
    hNtrackstrue[2]->Draw("same");
    leg->Draw();

    c1->cd(5);
    //c1->cd(4)->SetLogy();
    hSPDtracklets[0]->SetLineColor(kRed);
    hSPDtracklets[1]->SetLineColor(kBlue);
    hSPDtracklets[1]->SetLineStyle(7);
    hSPDtracklets[2]->SetLineColor(kGreen);
    hSPDtracklets[2]->SetLineStyle(9);
    hSPDtracklets[0]->Draw();
    hSPDtracklets[1]->Draw("same");
    hSPDtracklets[2]->Draw("same");
    leg->Draw();

   /*/ c1->cd(3);
    hZVtxReco[0]->SetLineColor(kRed);
    hZVtxReco[1]->SetLineColor(kBlue);
    hZVtxReco[1]->SetLineStyle(7);
    hZVtxReco[2]->SetLineColor(kGreen);
    hZVtxReco[2]->SetLineStyle(9);
    hZVtxReco[0]->Draw();
    hZVtxReco[1]->Draw("same");
    hZVtxReco[2]->Draw("same");
    leg->Draw();

    c1->cd(6);
    hZVtxMC[0]->SetLineColor(kRed);
    hZVtxMC[1]->SetLineColor(kBlue);
    hZVtxMC[1]->SetLineStyle(7);
    hZVtxMC[2]->SetLineColor(kGreen);
    hZVtxMC[2]->SetLineStyle(9);
    hZVtxMC[0]->Draw();
    hZVtxMC[1]->Draw("same");
    hZVtxMC[2]->Draw("same");
    leg->Draw();*/

}
