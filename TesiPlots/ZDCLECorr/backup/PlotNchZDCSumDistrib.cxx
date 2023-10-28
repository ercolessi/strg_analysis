void dostudy(TString period = "18i", TString sel = "SPD", float sel1 = 0., float sel2 = 100., Int_t* perc=0x0 , Long_t nbins=0 , Int_t* colors =0x0 );

void PlotNchZDCSumDistrib(){

    Int_t *colors;
    Int_t c0[] = {kRed+1, kRed-4, kOrange+7, kOrange-4, kYellow+1, kSpring-7, kGreen+2, kAzure+8, kBlue-4, kBlue+3};
    Int_t c7[] = {kRed+1,  kOrange+7, kYellow+1, kSpring-7, kAzure+8, kBlue-4, kBlue+3};
    Int_t c8[] = {kRed+1, kRed-4, kOrange+7, kYellow+1, kSpring-7, kAzure+8, kBlue-4, kBlue+3};

    /*Int_t p1[] = {0,5,10,20,30,40,50,80,100};
    Long_t n1 = sizeof(p1)/sizeof(Int_t) - 1;
    //
    dostudy("18i","SPD",10.,20., p1, n1, c0);

    Int_t p2[] = {0,20,30,40,50,60,70,100};
    Long_t n2 = sizeof(p2)/sizeof(Int_t) - 1;
    //
    dostudy("18i","SPD",40.,50., p2, n2, c7);

    Int_t p4[] = {0,5,10,20,30,40,50,100};
    Long_t n4 = sizeof(p4)/sizeof(Int_t) - 1;
    //
    dostudy("18i","V0M",10.,20.,p4,n4, c7);*/

    Int_t p5[] = {0,10,20,30,40,50,60,70,100};
    Long_t n5 = sizeof(p5)/sizeof(Int_t) - 1;
    //
    dostudy("18i","V0M",40.,50.,p5,n5, c8);

}

void dostudy(TString period = "18i", TString sel = "SPD", float sel1 = 0., float sel2 = 100., Int_t* perc =0x0 , Long_t nbins = 0, Int_t* colors =0x0){

    TFile* file = new TFile(Form("LHC%s.root",period.Data()),"READ");

	cout<<"--------------- Open Real Data File --------------------"<<endl;
    TTree* lTreeEvent = (TTree*)file->Get("PWGLF_StrVsMult/fTreeEvent");

    Float_t fCentrality_V0M = 0.;
    Float_t fCentrality_SPDClusters = 0.;
    Float_t fCentrality_ZDCFired = 0.;
    Int_t fRun = 0;
    Int_t fSPDTracklets = 0;
    Float_t fZPApp = 0.;
    Float_t fZPCpp = 0.;
    Float_t fZNApp = 0.;
    Float_t fZNCpp = 0.;

    // 
    lTreeEvent->SetBranchAddress("fRun",&fRun);
    lTreeEvent->SetBranchAddress("fCentrality_V0M", &fCentrality_V0M);
    lTreeEvent->SetBranchAddress("fCentrality_SPDClusters", &fCentrality_SPDClusters);
    lTreeEvent->SetBranchAddress("fCentrality_ZDCFired", &fCentrality_ZDCFired);
    lTreeEvent->SetBranchAddress("fSPDtracklets", &fSPDTracklets);
    lTreeEvent->SetBranchAddress("fZPApp", &fZPApp);
    lTreeEvent->SetBranchAddress("fZPCpp", &fZPCpp);
    lTreeEvent->SetBranchAddress("fZNApp", &fZNApp);
    lTreeEvent->SetBranchAddress("fZNCpp", &fZNCpp);

    TH1F* hnch[nbins], *hzdcsum[nbins];
    TH1F* hnchmb = new TH1F(Form("hnch"), "; SPD tracklets; Ratio to MB (normalized)", 35,0.,35.);
    TH1F *hzdcsummb = new TH1F(Form("hzdcsum"), "; ZDC Energy Sum (a.u.); Ratio to MB (normalized)", 200,0.,2000.);

    for (int n = 0; n < nbins; n++){

        hnch[n] = new TH1F(Form("hnch%i",n), "; SPD tracklets; Counts", 35,0.,35.);
        hzdcsum[n] = new TH1F(Form("hzdcsum%i",n), "; ZDC Energy Sum (a.u.); Counts", 200,0.,2000.);


        cout<<"\n\n--------------------------------------------------------"<<endl;
        cout<<" Will now loop over events, please wait..."<<endl;
        Long_t lNEvents = 0;
        for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries()/10; iEv++) {

            lTreeEvent->GetEntry(iEv);
            if( iEv % ( lTreeEvent->GetEntries() / 100 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;

            Float_t lZDCSum = fZNApp + fZNCpp + fZPApp + fZPCpp;

            if(fCentrality_ZDCFired>100) continue;

            if(sel.Contains("SPD")) {
                if(fCentrality_SPDClusters >= sel1 && fCentrality_SPDClusters <sel2 && fCentrality_V0M >= perc[n] && fCentrality_V0M < perc[n+1]){
                    hnch[n]->Fill(fSPDTracklets);
                    hzdcsum[n]->Fill(lZDCSum);
                }
            }
            if(sel.Contains("V0M")) {
                    if(fCentrality_V0M >= sel1 && fCentrality_V0M <sel2 && fCentrality_SPDClusters >= perc[n] && fCentrality_SPDClusters < perc[n+1]){
                    hnch[n]->Fill(fSPDTracklets);
                    hzdcsum[n]->Fill(lZDCSum);
                }             
            }

            if(n==0){
                hnchmb->Fill(fSPDTracklets);        
                hzdcsummb->Fill(lZDCSum);
            }

        }
            hnch[n]->SetLineColor(colors[n]);
            hnch[n]->SetMarkerColor(colors[n]);
            hnch[n]->SetMarkerStyle(20);
            hnch[n]->SetLineWidth(2);
            hzdcsum[n]->SetLineColor(colors[n]);
            hzdcsum[n]->SetMarkerColor(colors[n]);
            hzdcsum[n]->SetMarkerStyle(20);     
            hzdcsum[n]->SetLineWidth(2);     
    }

    TLegend* leg = new TLegend(0.25,0.1,0.5,0.4);
    leg->SetTextSize(0.02);
    leg->SetBorderSize(0); 
    leg->SetHeader(Form("%s fixed [%.0f-%.0f]",sel.Data(),sel1,sel2));
    for (int i = 0; i< nbins; i++){ 
        if (sel.Contains("V0M")) leg->AddEntry(hnch[i],Form("%s [%i-%i]", "SPDcl" ,perc[i], perc[i+1]),"LEP");
        if (sel.Contains("SPD")) leg->AddEntry(hnch[i],Form("%s [%i-%i]", "V0M" ,perc[i], perc[i+1]),"LEP");
    }



    TCanvas* c = new TCanvas(Form("c_%s%.0f",sel.Data(),sel1), "", 1000, 800);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.12);
    c->SetTopMargin(0.1);
    c->SetTicky();
    c->SetTickx();
    hnch[0]->GetYaxis()->SetRangeUser(0.,10000);
    for (int n = 0; n < nbins; n++){
        hnch[n]->Divide(hnchmb);
        hnch[n]->DrawNormalized("SAME");
    }
    leg->Draw("SAME");

    TCanvas* c1 = new TCanvas(Form("c1_%s%.0f",sel.Data(),sel1), "", 1000, 800);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    c1->SetRightMargin(0.12);
    c1->SetTopMargin(0.1);
    c1->SetTicky();
    c1->SetTickx();
    c1->SetLogy();
    hzdcsum[0]->GetYaxis()->SetRangeUser(0.01,6000);
    for (int n = 0; n < nbins; n++){
   //     hzdcsum[n]->Divide(hzdcsummb);
        hzdcsum[n]->DrawNormalized("SAME");
    }
    leg->Draw("SAME");
}
