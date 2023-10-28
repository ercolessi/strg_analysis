void drawsgnloss(
    TString lCascType = "Xi",
    Double_t lFixedLo = 40.0,
	Double_t lFixedHi = 50.0,
	TString lWhichVarEstimator = "V0M",
	TString lWhichFixedEstimator = "SPDClusters"
);

void drawevtloss(
    TString lCascType = "Xi",
    Double_t lFixedLo = 40.0,
	Double_t lFixedHi = 50.0,
	TString lWhichVarEstimator = "V0M",
	TString lWhichFixedEstimator = "SPDClusters"
);

void drawsgnlossxperiod(
    TString lCascType = "Xi",
    Double_t lFixedLo = 40.0,
	Double_t lFixedHi = 50.0,
	TString lWhichVarEstimator = "V0M",
	TString lWhichFixedEstimator = "SPDClusters"
);

void drawevtlossxperiod(
    TString lCascType = "Xi",
    Double_t lFixedLo = 40.0,
	Double_t lFixedHi = 50.0,
	TString lWhichVarEstimator = "V0M",
	TString lWhichFixedEstimator = "SPDClusters"
);


void DrawNormCorrections(TString lPart = "K0Short"){

    /*
    drawevtloss(Form("%s", lPart.Data()), 0, 100, "V0M", "SPDClusters");
    drawevtloss(Form("%s", lPart.Data()), 10, 20, "V0M", "SPDClusters");
    drawevtloss(Form("%s", lPart.Data()), 40, 50, "V0M", "SPDClusters");
    drawevtloss(Form("%s", lPart.Data()), 10, 20, "SPDClusters", "V0M");
    drawevtloss(Form("%s", lPart.Data()), 40, 50, "SPDClusters", "V0M");
    */

    drawsgnloss(Form("%s", lPart.Data()), 0, 100, "V0M", "SPDClusters");
    /*drawsgnloss(Form("%s", lPart.Data()), 10, 20, "V0M", "SPDClusters");
    drawsgnloss(Form("%s", lPart.Data()), 40, 50, "V0M", "SPDClusters");
    drawsgnloss(Form("%s", lPart.Data()), 10, 20, "SPDClusters", "V0M");
    drawsgnloss(Form("%s", lPart.Data()), 40, 50, "SPDClusters", "V0M");*/
}





void drawsgnloss(
    TString lCascType = "Xi",
    Double_t lFixedLo = 40.0,
	Double_t lFixedHi = 50.0,
	TString lWhichVarEstimator = "V0M",
	TString lWhichFixedEstimator = "SPDClusters"
){

    TFile* file = TFile::Open(
        Form("SignalLoss-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root",
            lCascType.Data(), lWhichVarEstimator.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi)
    );

    //Percentile
    Float_t * percentile;
    Long_t percbinnumb;
    Float_t p0[] = {0,1,5,10,15,20,30,40,50,70,100};
    Long_t n0 = sizeof(p0)/sizeof(Float_t) - 1;
    Float_t p1[] = {0,5,10,20,30,40,50,100}; //7
    Long_t n1 = sizeof(p1)/sizeof(Float_t) - 1;
    Float_t p2[] = {0,20,30,40,50,60,70,100}; //7
    Long_t n2 = sizeof(p2)/sizeof(Float_t) - 1;
    Float_t p4[] = {0,5,10,20,30,40,50,100}; //7
    Long_t n4 = sizeof(p4)/sizeof(Float_t) - 1;
    Float_t p5[] = {0,10,20,30,40,50,60,70,100}; //8
    Long_t n5 = sizeof(p5)/sizeof(Float_t) - 1;
     Float_t pOmega[] = {0,5,10,30,50,100};
    Long_t nOmega = sizeof(pOmega)/sizeof(Float_t) - 1;
    Float_t pOmega2[] = {0,40,70,100};
    Long_t nOmega2 = sizeof(pOmega2)/sizeof(Float_t) - 1;
    Float_t pOmega3[] = {0,10,30,50,100};
    Long_t nOmega3 = sizeof(pOmega3)/sizeof(Float_t) - 1;

    Int_t *colors;
    Int_t c0[] = {kRed + 3, kRed + 1, kRed - 4, kOrange + 7, kYellow + 1, kSpring - 7, kGreen + 2, kAzure + 8, kBlue - 4, kBlue + 3};
    Int_t c7[] = {kRed+1,  kOrange+7, kYellow+1, kSpring-7, kAzure+8, kBlue-4, kBlue+3};
    Int_t c8[] = {kRed+1, kRed-4, kOrange+7, kYellow+1, kSpring-7, kAzure+8, kBlue-4, kBlue+3};
    Int_t c4[] = {kRed+1, kYellow+1, kAzure+8, kBlue+3};


   if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 0. && lFixedHi == 100.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = nOmega;
            colors = c8;
        } else{
            percentile = p0;
            percbinnumb = n0;
            colors = c0;
        }
    }
    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = nOmega;
            colors = c8;
        } else{
            percentile = p1;
            percbinnumb = n1;
            colors = c7;
        }
    }
    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega2;
            percbinnumb = nOmega2;
            colors = c4;
        } else{
            percentile = p2;
            percbinnumb = n2;
            colors = c7;
        }
    }
    if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega3;
            percbinnumb = nOmega3;
            colors = c4;
        } else{
        percentile = p4;
        percbinnumb = n4;
        colors = c7;
    }
    }
    if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega2;
            percbinnumb = nOmega2;
            colors = c4;
        } else{
        percentile = p5;
        percbinnumb = n5;
        colors = c8;
        }
    }

    TH1D* hsgnloss[percbinnumb];
    for (int iperc = 0; iperc < percbinnumb; iperc ++){
        hsgnloss[iperc] = (TH1D*)file->Get(Form("wsgnloss_%i-%i",(int)percentile[iperc],(int)percentile[iperc+1]));
        hsgnloss[iperc]->SetMarkerStyle(20);
        hsgnloss[iperc]->SetMarkerColor(colors[iperc]);
        hsgnloss[iperc]->SetLineColor(colors[iperc]);
        hsgnloss[iperc]->SetMarkerSize(2.);
    }

    TLatex *lab0 = new TLatex();
    lab0->SetTextFont(42);
    lab0->SetNDC();
    lab0->SetTextColor(1);
    lab0->SetTextSize(0.03);
    lab0->SetTextAlign(22);
    lab0->SetTextAngle(0);
    TLatex *labbig = new TLatex();
    labbig->SetTextFont(42);
    labbig->SetNDC();
    labbig->SetTextColor(1);
    labbig->SetTextSize(0.07);
    labbig->SetTextAlign(22);
    labbig->SetTextAngle(0);

    TCanvas* c2 = new TCanvas(Form("c2"),"",1200,1100);
    c2->SetLeftMargin(0.15);
    c2->SetBottomMargin(0.12);
    c2->SetRightMargin(0.1);
    c2->SetTopMargin(0.1);
    c2->SetTicky();
    c2->SetTickx();

    TLegend* leg = new TLegend(0.5,0.15,0.8,0.45);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);

    TH1D* h = new TH1D("h",";#it{p}_{T} (GeV/#it{c}); #varepsilon_{part}",10,0.,10.);
    h->SetStats(0);
    //h->GetYaxis()->SetRangeUser(0.7,1.05);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->Draw();
    h->SetMinimum(hsgnloss[percbinnumb-1]->GetMinimum()*0.95);
    h->SetMaximum(1.01);
    for (int iperc = 0; iperc < percbinnumb; iperc ++){
        hsgnloss[iperc]->Draw("SAME EP");
        leg->AddEntry(hsgnloss[iperc],Form("%s [%.0f-%.0f]",lWhichVarEstimator.Data(),percentile[iperc],percentile[iperc+1]),"LEP");
    }
    leg->Draw("SAME");
    labbig->DrawLatex(0.48, 0.45, Form("#%s",lCascType.Data()));
    if (lFixedLo != 0 && lFixedHi != 100) lab0->DrawLatex(0.36, 0.38, Form("%s fixed in [%.0f-%.0f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));


    c2->SaveAs(Form("images/SignalLoss_%sFixed%sin%03.0f_%03.0f.png", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
    c2->SaveAs(Form("images/SignalLoss_%sFixed%sin%03.0f_%03.0f.pdf", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
}


void drawsgnlossxperiod(
    TString lCascType = "Xi",
    Double_t lFixedLo = 40.0,
	Double_t lFixedHi = 50.0,
	TString lWhichVarEstimator = "V0M",
	TString lWhichFixedEstimator = "SPDClusters"
){

    TFile* file = TFile::Open(
        Form("SignalLoss-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root",
            lCascType.Data(), lWhichVarEstimator.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi)
    );

    //Percentile
    Float_t * percentile;
    Long_t percbinnumb;
    Float_t p0[] = {0,5,10,15,20,30,40,50,70,100};
    Long_t n0 = sizeof(p0)/sizeof(Float_t) - 1;
    Float_t p1[] = {0,5,10,20,30,40,50,100}; //7
    Long_t n1 = sizeof(p1)/sizeof(Float_t) - 1;
    Float_t p2[] = {0,20,30,40,50,60,70,100}; //7
    Long_t n2 = sizeof(p2)/sizeof(Float_t) - 1;
    Float_t p4[] = {0,5,10,20,30,40,50,100}; //7
    Long_t n4 = sizeof(p4)/sizeof(Float_t) - 1;
    Float_t p5[] = {0,10,20,30,40,50,60,70,100}; //8
    Long_t n5 = sizeof(p5)/sizeof(Float_t) - 1;
     Float_t pOmega[] = {0,5,10,30,50,100};
    Long_t nOmega = sizeof(pOmega)/sizeof(Float_t) - 1;
    Float_t pOmega2[] = {0,40,70,100};
    Long_t nOmega2 = sizeof(pOmega2)/sizeof(Float_t) - 1;

    Int_t *colors;
    Int_t c0[] = {kRed+1, kRed-4, kOrange+7, kYellow+1, kSpring-7, kGreen+2, kAzure+8, kBlue-4, kBlue+3};
    Int_t c7[] = {kRed+1,  kOrange+7, kYellow+1, kSpring-7, kAzure+8, kBlue-4, kBlue+3};
    Int_t c8[] = {kRed+1, kRed-4, kOrange+7, kYellow+1, kSpring-7, kAzure+8, kBlue-4, kBlue+3};
    Int_t c4[] = {kRed+1, kYellow+1, kAzure+8, kBlue+3};


   if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 0. && lFixedHi == 100.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = nOmega;
            colors = c8;
        } else{
            percentile = p0;
            percbinnumb = n0;
            colors = c0;
        }
    }
    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = nOmega;
            colors = c8;
        } else{
            percentile = p1;
            percbinnumb = n1;
            colors = c7;
        }
    }
    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega2;
            percbinnumb = nOmega2;
            colors = c4;
        } else{
            percentile = p2;
            percbinnumb = n2;
            colors = c7;
        }
    }
    if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        percentile = p4;
        percbinnumb = n4;
        colors = c7;
    }
    if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
        percentile = p5;
        percbinnumb = n5;
        colors = c8;
    }

    TH1D* hsgnloss0[percbinnumb],* hsgnloss1[percbinnumb],* hsgnloss2[percbinnumb],* hsgnloss[percbinnumb];
    for (int iperc = 0; iperc < percbinnumb; iperc ++){
        hsgnloss0[iperc] = (TH1D*)file->Get(Form("hsgnloss_%i-%i_MC0",(int)percentile[iperc],(int)percentile[iperc+1]));
        hsgnloss0[iperc]->SetMarkerStyle(34);
        hsgnloss0[iperc]->SetMarkerColor(colors[iperc]);
        hsgnloss0[iperc]->SetLineColor(colors[iperc]);
        hsgnloss0[iperc]->SetMarkerSize(1.5);
        hsgnloss1[iperc] = (TH1D*)file->Get(Form("hsgnloss_%i-%i_MC1",(int)percentile[iperc],(int)percentile[iperc+1]));
        hsgnloss1[iperc]->SetMarkerStyle(3);
        hsgnloss1[iperc]->SetMarkerColor(colors[iperc]);
        hsgnloss1[iperc]->SetLineColor(colors[iperc]);
        hsgnloss1[iperc]->SetMarkerSize(1.5);
        hsgnloss2[iperc] = (TH1D*)file->Get(Form("hsgnloss_%i-%i_MC2",(int)percentile[iperc],(int)percentile[iperc+1]));
        hsgnloss2[iperc]->SetMarkerStyle(23);
        hsgnloss2[iperc]->SetMarkerColor(colors[iperc]);
        hsgnloss2[iperc]->SetLineColor(colors[iperc]);
        hsgnloss2[iperc]->SetMarkerSize(1.5);
        hsgnloss[iperc] = (TH1D*)file->Get(Form("wsgnloss_%i-%i",(int)percentile[iperc],(int)percentile[iperc+1]));
        hsgnloss[iperc]->SetMarkerStyle(1);
        hsgnloss[iperc]->SetLineWidth(3);
        hsgnloss[iperc]->SetMarkerColor(kBlack);
        hsgnloss[iperc]->SetLineColor(kBlack);
        hsgnloss[iperc]->SetMarkerSize(.1);

    }

    TLatex *lab0 = new TLatex();
    lab0->SetTextFont(42);
    lab0->SetNDC();
    lab0->SetTextColor(1);
    lab0->SetTextSize(0.07);
    lab0->SetTextAlign(22);
    lab0->SetTextAngle(0);
    TLatex *labbig = new TLatex();
    labbig->SetTextFont(42);
    labbig->SetNDC();
    labbig->SetTextColor(1);
    labbig->SetTextSize(0.1);
    labbig->SetTextAlign(22);
    labbig->SetTextAngle(0);

    TCanvas* c2 = new TCanvas(Form("c2"),"",1300,1100);
    c2->SetLeftMargin(0.15);
    c2->SetBottomMargin(0.12);
    c2->SetRightMargin(0.1);
    c2->SetTopMargin(0.1);
    c2->SetTicky();
    c2->SetTickx();
    c2->Divide(3,3);
        //percbinnumb/2,percbinnumb/2);

    TLegend* leg = new TLegend(0.2,0.15,0.5,0.5);
    leg->SetTextSize(0.07);
    leg->SetBorderSize(0);

    TH1D* h = new TH1D("h",";#it{p}_{T} (GeV/#it{c}); #varepsilon_{part}",10,0.6,6.5);
    h->SetStats(0);
   // h->GetYaxis()->SetRangeUser(0.55,1.005);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleSize(0.04);
    h->SetMinimum(hsgnloss[percbinnumb-1]->GetMinimum()*0.95);
    h->SetMaximum(1.01);



    for (int iperc = 0; iperc < percbinnumb; iperc ++){
        c2->cd(iperc+1);
        h->Draw();
        hsgnloss0[iperc]->Draw("SAME EP");
        hsgnloss1[iperc]->Draw("SAME EP");
        hsgnloss2[iperc]->Draw("SAME EP");
        hsgnloss[iperc]->Draw("SAME EP");
        lab0->DrawLatex(0.66, 0.28, Form("%s [%.0f-%.0f] %s",lWhichVarEstimator.Data(),percentile[iperc],percentile[iperc+1],"%"));

    }
    leg->AddEntry(hsgnloss0[0],"LHC15f","P");
    leg->AddEntry(hsgnloss1[0],"LHC17j","P");
    leg->AddEntry(hsgnloss2[0],"LHC18i","P");
    leg->AddEntry(hsgnloss[0],"weighted","L");
    c2->cd(percbinnumb+1);
    leg->Draw();
    labbig->DrawLatex(0.28, 0.75, Form("#%s",lCascType.Data()));
    if (lFixedLo != 0 && lFixedHi != 100) lab0->DrawLatex(0.36, 0.58, Form("%s fixed in [%.0f-%.0f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));

    c2->SaveAs(Form("images/SignalLossXPeriod_%sFixed%sin%03.0f_%03.0f.png", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
    c2->SaveAs(Form("images/SignalLossXPeriod_%sFixed%sin%03.0f_%03.0f.pdf", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
}



void drawevtloss(
    TString lCascType = "Xi",
    Double_t lFixedLo = 40.0,
	Double_t lFixedHi = 50.0,
	TString lWhichVarEstimator = "V0M",
	TString lWhichFixedEstimator = "SPDClusters"
){

    TFile* file = TFile::Open(
        Form("EventLoss-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root",
            lCascType.Data(), lWhichVarEstimator.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi)
    );

    //Percentile
    Float_t * percentile;
    Long_t percbinnumb;
    Float_t p0[] = {0,1,5,10,15,20,30,40,50,70,100};
    Long_t n0 = sizeof(p0)/sizeof(Float_t) - 1;
    Float_t p1[] = {0,5,10,20,30,40,50,100}; //7
    Long_t n1 = sizeof(p1)/sizeof(Float_t) - 1;
    Float_t p2[] = {0,20,30,40,50,60,70,100}; //7
    Long_t n2 = sizeof(p2)/sizeof(Float_t) - 1;
    Float_t p4[] = {0,5,10,20,30,40,50,100}; //7
    Long_t n4 = sizeof(p4)/sizeof(Float_t) - 1;
    Float_t p5[] = {0,10,20,30,40,50,60,70,100}; //8
    Long_t n5 = sizeof(p5)/sizeof(Float_t) - 1;
     Float_t pOmega[] = {0,5,10,30,50,100};
    Long_t nOmega = sizeof(pOmega)/sizeof(Float_t) - 1;
    Float_t pOmega2[] = {0,40,70,100};
    Long_t nOmega2 = sizeof(pOmega2)/sizeof(Float_t) - 1;

    Int_t *colors;
    Int_t c0[] = {kRed+1, kRed-4, kOrange+7, kYellow+1, kSpring-7, kGreen+2, kAzure+8, kBlue-4, kBlue+3};
    Int_t c7[] = {kRed+1,  kOrange+7, kYellow+1, kSpring-7, kAzure+8, kBlue-4, kBlue+3};
    Int_t c8[] = {kRed+1, kRed-4, kOrange+7, kYellow+1, kSpring-7, kAzure+8, kBlue-4, kBlue+3};
    Int_t c4[] = {kRed+1, kYellow+1, kAzure+8, kBlue+3};


   if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 0. && lFixedHi == 100.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = nOmega;
            colors = c4;
        } else{
            percentile = p0;
            percbinnumb = n0;
            colors = c0;
        }
    }
    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = nOmega;
            colors = c4;
        } else{
            percentile = p1;
            percbinnumb = n1;
            colors = c7;
        }
    }
    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega2;
            percbinnumb = nOmega2;
            colors = c4;
        } else{
            percentile = p2;
            percbinnumb = n2;
            colors = c7;
        }
    }
    if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        percentile = p4;
        percbinnumb = n4;
        colors = c7;
    }
    if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
        percentile = p5;
        percbinnumb = n5;
        colors = c8;
    }

    TH1D* hevtloss = (TH1D*)file->Get("EventLoss/hevtloss");

    TLatex *lab0 = new TLatex();
    lab0->SetTextFont(42);
    lab0->SetNDC();
    lab0->SetTextColor(1);
    lab0->SetTextSize(0.03);
    lab0->SetTextAlign(22);
    lab0->SetTextAngle(0);
    TLatex *labbig = new TLatex();
    labbig->SetTextFont(42);
    labbig->SetNDC();
    labbig->SetTextColor(1);
    labbig->SetTextSize(0.07);
    labbig->SetTextAlign(22);
    labbig->SetTextAngle(0);

    TCanvas* c1 = new TCanvas(Form("c1"),"",1200,1100);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.12);
    c1->SetRightMargin(0.1);
    c1->SetTopMargin(0.1);
    c1->SetTicky();
    c1->SetTickx();
    hevtloss->SetLineColor(kBlue);
    hevtloss->SetMarkerColor(kBlue);
    hevtloss->SetMarkerStyle(20);
    hevtloss->SetMarkerSize(1.8);
    hevtloss->GetYaxis()->SetTitle("#varepsilon_{event}");
    hevtloss->GetYaxis()->SetTitleSize(0.04);
    hevtloss->SetStats(0);
    //hevtloss->GetYaxis()->SetRangeUser(0.75,1.05);
    hevtloss->GetXaxis()->SetTitleOffset(1.2);
    hevtloss->Draw();
    //labbig->DrawLatex(0.28, 0.45, Form("#%s",lCascType.Data()));
    if (lFixedLo != 0 && lFixedHi != 100) lab0->DrawLatex(0.36, 0.38, Form("%s fixed in [%.0f-%.0f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));

    c1->SaveAs(Form("images/EventLoss_%sFixed%sin%03.0f_%03.0f.png", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
    c1->SaveAs(Form("images/EventLoss_%sFixed%sin%03.0f_%03.0f.pdf", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
}

