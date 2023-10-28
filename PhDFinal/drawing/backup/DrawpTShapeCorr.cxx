void DrawpTShapeCorr(
    TString lCascType = "Xi",
	Double_t lFixedLo = 40.0,
	Double_t lFixedHi = 50.0,
	TString lWhichVarEstimator = "V0M",
	TString lWhichFixedEstimator = "SPDClusters",
    Bool_t DoMB = kFALSE){

         //Percentile
    Float_t * percentile;
    Long_t percbinnumb;
    Float_t pmb[] = {0,100};
    Long_t nmb = sizeof(pmb)/sizeof(Float_t) - 1;
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
    Float_t pOmega[] = {0,10,30,50,100}; //4
    Long_t nOmega = sizeof(pOmega)/sizeof(Float_t) - 1;

    Int_t *colors;
    Int_t c0[] = {kRed+1, kRed-4, kOrange+1, kYellow+1, kSpring-1, kGreen+2, kCyan+1,  kAzure+7, kBlue-4, kBlue+2};
    Int_t c7[] = {kRed+1,  kOrange+1, kYellow+1, kSpring-1, kAzure+7, kBlue-4, kBlue+2};
    Int_t c8[] = {kRed+1, kRed-4, kOrange+1, kYellow+1, kSpring-1, kAzure+7, kBlue-4, kBlue+2};
    Int_t c4[] = {kRed+1, kYellow+1, kAzure+7, kBlue+3};

    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 0. && lFixedHi == 100.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = nOmega;
            colors = c4;
        } else {
            percentile = p0;
            percbinnumb = n0;
            colors = c0;
        }
        if (DoMB){
            percentile = pmb;
            percbinnumb = nmb;
            colors = c0;
        }
    }
    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = nOmega;
            colors = c4;
        } else {
            percentile = p1;
            percbinnumb = n1;
            colors = c7;
        }
    }
    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        percentile = p2;
        percbinnumb = n2;
        colors = c7;
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
    }//

    TFile* file;
    if (!DoMB){
        file = TFile::Open(Form("../correctedspectra/CorrSpectra-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root",
                                lCascType.Data(), lWhichVarEstimator.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
    } else {
        file = TFile::Open(Form("../correctedspectra/CorrSpectra-%s-13TeV_INELgt0.root",
                                lCascType.Data()));
    }

        TH1F* hRatio[percbinnumb][3], *hRatioIT1[percbinnumb][3], *hRatioIT2[percbinnumb][3];
        TH1F* hGenSpectraIT0[percbinnumb];
        TH1F* hGenMCIT0[3], * hGenMCIT0Clone[3];
        TH1F* CorrFactors[percbinnumb],* CorrFactorsIT1[percbinnumb], * CorrFactorsIT2[percbinnumb];
        for(int imc = 0; imc < 3 ; imc++){
            for (int i = 0; i< percbinnumb; i++){
                hRatio[i][imc] = (TH1F*)file->Get(Form("MCpTshape/Iteration0/RatioPtShapeAndCorrFactors/hRatioPtShapeIT0PercBin%i_MC%i",i,imc));
                hRatio[i][imc]->SetMarkerColor(colors[i]);
                hRatio[i][imc]->SetLineColor(colors[i]);
                hRatio[i][imc]->SetMarkerStyle(20);
                hRatio[i][imc]->SetMarkerSize(0.8);
                //
                hRatioIT1[i][imc] = (TH1F*)file->Get(Form("MCpTshape/Iteration1/RatioPtShapeAndCorrFactors/hRatioPtShapeIT1PercBin%i_MC%i",i,imc));
                hRatioIT1[i][imc]->SetMarkerColor(colors[i]);
                hRatioIT1[i][imc]->SetLineColor(colors[i]);
                hRatioIT1[i][imc]->SetMarkerStyle(20);
                hRatioIT1[i][imc]->SetMarkerSize(0.8);
                //
                hRatioIT2[i][imc] = (TH1F*)file->Get(Form("MCpTshape/Iteration2/RatioPtShapeAndCorrFactors/hRatioPtShapeIT2PercBin%i_MC%i",i,imc));
                hRatioIT2[i][imc]->SetMarkerColor(colors[i]);
                hRatioIT2[i][imc]->SetLineColor(colors[i]);
                hRatioIT2[i][imc]->SetMarkerStyle(20);
                hRatioIT2[i][imc]->SetMarkerSize(0.8);
                //
                if (imc==0){
                    hGenSpectraIT0[i] = (TH1F*)file->Get(Form("MCpTshape/Iteration0/GeneratedSpectra/hPtSpectrumGenIT0PercBin%i",i));
                    hGenSpectraIT0[i]->SetMarkerColor(colors[i]);
                    hGenSpectraIT0[i]->SetLineColor(colors[i]);
                    hGenSpectraIT0[i]->SetMarkerStyle(20);
                    hGenSpectraIT0[i]->SetMarkerSize(1);
                    //
                    CorrFactors[i] = (TH1F*)file->Get(Form("MCpTshape/Iteration0/RatioPtShapeAndCorrFactors/CorrFactorsIT0PercBin%i",i));
                    CorrFactors[i]->SetMarkerColor(colors[i]);
                    CorrFactors[i]->SetLineColor(colors[i]);
                    CorrFactors[i]->SetMarkerStyle(20);
                    //
                    CorrFactorsIT1[i] = (TH1F*)file->Get(Form("MCpTshape/Iteration1/RatioPtShapeAndCorrFactors/CorrFactorsIT1PercBin%i",i));
                    CorrFactorsIT1[i]->SetMarkerColor(colors[i]);
                    CorrFactorsIT1[i]->SetLineColor(colors[i]);
                    CorrFactorsIT1[i]->SetMarkerStyle(20);
                    //
                    CorrFactorsIT2[i] = (TH1F*)file->Get(Form("MCpTshape/Iteration2/RatioPtShapeAndCorrFactors/CorrFactorsIT2PercBin%i",i));
                    CorrFactorsIT2[i]->SetMarkerColor(colors[i]);
                    CorrFactorsIT2[i]->SetLineColor(colors[i]);
                    CorrFactorsIT2[i]->SetMarkerStyle(20);
                }
            }
            hGenMCIT0[imc] = (TH1F*)file->Get(Form("MCpTshape/Iteration0/GeneratedSpectra/hPtSpectrumGenMCIT0MB%i",imc));
            hGenMCIT0[imc]->SetLineColor(kBlack);
            hGenMCIT0[imc]->SetLineWidth(2);
            hGenMCIT0Clone[imc] = (TH1F*)hGenMCIT0[imc]->Clone(Form("hGenMCIT0%i",imc));
        }
        hGenMCIT0Clone[1]->SetLineStyle(9);
        hGenMCIT0Clone[2]->SetLineStyle(4);

        TLatex *lab0 = new TLatex();
        lab0->SetTextFont(42);
        lab0->SetNDC();
        lab0->SetTextColor(1);
        lab0->SetTextSize(0.04);
        lab0->SetTextAlign(22);
        lab0->SetTextAngle(0);
        TLatex *labbig = new TLatex();
        labbig->SetTextFont(42);
        labbig->SetNDC();
        labbig->SetTextColor(1);
        labbig->SetTextSize(0.06);
        labbig->SetTextAlign(22);
        labbig->SetTextAngle(0);


        TCanvas* c1 = new TCanvas(Form("c1"),"",800,1100);
        c1->SetLeftMargin(0.15);
        c1->SetBottomMargin(0.08);
        c1->SetRightMargin(0.05);
        c1->SetTopMargin(0.05);
        c1->SetTicky();
        c1->SetTickx();
        c1->SetLogy();



        TH1D* h1 = new TH1D("h1", "",7, 0.6,6.5);
        h1->GetYaxis()->SetRangeUser(3E+3,5E+8);
        h1->GetYaxis()->SetTitle("dN/dp_{T} (a.u.)");
        h1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h1->GetYaxis()->SetTitleOffset(1.5);
        h1->GetXaxis()->SetRangeUser(0.6, 6.5);
        h1->SetStats(0);
        h1->SetTitle("");
        h1->Draw();
        for (int i = 0; i< percbinnumb; i++){
            hGenSpectraIT0[i]->Draw("SAME EP");
        }
        hGenMCIT0[0]->Draw("SAME L");
        hGenMCIT0[1]->Draw("SAME L");
        hGenMCIT0[1]->SetLineStyle(9);
        hGenMCIT0[2]->Draw("SAME L");
        hGenMCIT0[2]->SetLineStyle(2);
        labbig->DrawLatex(0.78, 0.85, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));

        //
        TLegend* leg = new TLegend(0.25,0.1,0.5,0.4);
        leg->SetTextSize(0.02);
        leg->SetBorderSize(0);
        leg->AddEntry(hGenMCIT0[0],"MC Input LHC15g3b1","L");
        leg->AddEntry(hGenMCIT0[1],"MC Input LHC22e1_17j","L");
        leg->AddEntry(hGenMCIT0[2],"MC Input LHC22e1_18i","L");
        for (int i = 0; i< percbinnumb; i++){
            leg->AddEntry(hGenSpectraIT0[i],Form("%s [%.0f-%.0f] corr. spectrum", lWhichVarEstimator.Data(),percentile[i], percentile[i+1]),"LEP");
        }
        leg->Draw("SAME");

        TLegend* leg2 = new TLegend(0.2,0.15,0.5,0.25);
        leg2->SetTextSize(0.02);
        leg2->SetBorderSize(0);
        leg2->AddEntry(hGenSpectraIT0[0],Form("%s [%.0f-%.0f] corr. spectrum", lWhichVarEstimator.Data(),percentile[0], percentile[1]),"LEP");
        leg2->AddEntry(hGenSpectraIT0[percbinnumb-1],Form("%s [%.0f-%.0f] corr. spectrum", lWhichVarEstimator.Data(),percentile[percbinnumb-1], percentile[percbinnumb]),"LEP");


        if (!DoMB){
        c1->SaveAs(Form("mcptshapeimages/GenSpectraIT0_%sFixed%sin%03.0f_%03.0f.png", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        c1->SaveAs(Form("mcptshapeimages/GenSpectraIT0_%sFixed%sin%03.0f_%03.0f.root", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        c1->SaveAs(Form("mcptshapeimages/GenSpectraIT0_%sFixed%sin%03.0f_%03.0f.pdf", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        } else {
        c1->SaveAs(Form("mcptshapeimages/GenSpectraIT0_%s_INELgt0.png", lCascType.Data()));
        c1->SaveAs(Form("mcptshapeimages/GenSpectraIT0_%s_INELgt0.root", lCascType.Data()));
        c1->SaveAs(Form("mcptshapeimages/GenSpectraIT0_%s_INELgt0.pdf", lCascType.Data()));
        }

        TCanvas* c = new TCanvas(Form("c"),"",1000,900);
        c->Divide(2,2);
        c->cd(1);
        c->cd(1)->SetLeftMargin(0.15);
        c->cd(1)->SetBottomMargin(0.1);
        c->cd(1)->SetRightMargin(0.05);
        c->cd(1)->SetTopMargin(0.05);
        c->cd(1)->SetTicky();
        c->cd(1)->SetTickx();
        //
        TH1D* h = new TH1D("h", "",7, 0.6,6.5);
        h->GetYaxis()->SetRangeUser(0., 2.);
        h->GetYaxis()->SetTitle("Corr.Spectrum(DATA) / Input(MC)");
        h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h->GetYaxis()->SetTitleOffset(1.5);
        h->GetXaxis()->SetRangeUser(0.6, 6.5);
        h->SetStats(0);
        h->SetTitle("");
        h->Draw();
        for (int i = 0; i< percbinnumb; i++){
            hRatio[i][0]->Draw("SAME");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 0");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.8, 0.86, "LHC15g3b1");
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
        leg2->Draw("SAME");


        c->cd(2);
        c->cd(2)->SetLeftMargin(0.15);
        c->cd(2)->SetBottomMargin(0.1);
        c->cd(2)->SetRightMargin(0.05);
        c->cd(2)->SetTopMargin(0.05);
        c->cd(2)->SetTicky();
        c->cd(2)->SetTickx();
        //
        TH1D* hh = (TH1D*)h->Clone("hh");
        hh->GetYaxis()->SetRangeUser(0., 2.);
        hh->Draw();
        for (int i = 0; i< percbinnumb; i++){
            hRatio[i][1]->Draw("SAME");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 0");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.8, 0.86, "LHC22e1_17j");
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
        leg2->Draw("SAME");



        c->cd(3);
        c->cd(3)->SetLeftMargin(0.15);
        c->cd(3)->SetBottomMargin(0.1);
        c->cd(3)->SetRightMargin(0.05);
        c->cd(3)->SetTopMargin(0.05);
        c->cd(3)->SetTicky();
        c->cd(3)->SetTickx();
        //
        hh->Draw();
        for (int i = 0; i< percbinnumb; i++){
            hRatio[i][2]->Draw("SAME");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 0");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.8, 0.86, "LHC22e1_18i");
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
        leg2->Draw("SAME");



        c->cd(4);
        c->cd(4)->SetLeftMargin(0.15);
        c->cd(4)->SetBottomMargin(0.1);
        c->cd(4)->SetRightMargin(0.05);
        c->cd(4)->SetTopMargin(0.05);
        c->cd(4)->SetTicky();
        c->cd(4)->SetTickx();
        //
        TH1D* h2 = new TH1D("h2", "",7, 0.6,6.5);
        h2->GetYaxis()->SetRangeUser(0.85, 1.15);
        h2->GetYaxis()->SetTitle("Corr. factor (#varepsilon_{it}/#varepsilon_{it-1})");
        h2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        h2->GetYaxis()->SetTitleOffset(1.5);
        h2->GetXaxis()->SetRangeUser(0.6, 6.5);
        h2->SetStats(0);
        h2->SetTitle("");
        h2->Draw();
        for (int i = 0; i< percbinnumb; i++){
            CorrFactors[i]->Draw("SAME LEP");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 0");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
        leg2->Draw("SAME");

        if (!DoMB){
        c->SaveAs(Form("mcptshapeimages/RatioIT0_%sFixed%sin%03.0f_%03.0f.png", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        c->SaveAs(Form("mcptshapeimages/RatioIT0_%sFixed%sin%03.0f_%03.0f.root", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        c->SaveAs(Form("mcptshapeimages/RatioIT0_%sFixed%sin%03.0f_%03.0f.pdf", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        } else {
        c->SaveAs(Form("mcptshapeimages/RatioIT0_%s_INELgt0.png", lCascType.Data()));
        c->SaveAs(Form("mcptshapeimages/RatioIT0_%s_INELgt0.root", lCascType.Data()));
        c->SaveAs(Form("mcptshapeimages/RatioIT0_%s_INELgt0.pdf", lCascType.Data()));
        }

        //IT1
        TCanvas* cc = new TCanvas(Form("cc"),"",1000,900);
        cc->Divide(2,2);
        cc->cd(1);
        cc->cd(1)->SetLeftMargin(0.15);
        cc->cd(1)->SetBottomMargin(0.1);
        cc->cd(1)->SetRightMargin(0.05);
        cc->cd(1)->SetTopMargin(0.05);
        cc->cd(1)->SetTicky();
        cc->cd(1)->SetTickx();
        //
        h->Draw();
        for (int i = 0; i< percbinnumb; i++){
            hRatioIT1[i][0]->Draw("SAME");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 1");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.8, 0.86, "LHC15g3b1");
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
        leg2->Draw("SAME");


        cc->cd(2);
        cc->cd(2)->SetLeftMargin(0.15);
        cc->cd(2)->SetBottomMargin(0.1);
        cc->cd(2)->SetRightMargin(0.05);
        cc->cd(2)->SetTopMargin(0.05);
        cc->cd(2)->SetTicky();
        cc->cd(2)->SetTickx();
        //
        h->Draw();
        for (int i = 0; i< percbinnumb; i++){
            hRatioIT1[i][1]->Draw("SAME");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 1");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.8, 0.86, "LHC22e1_17j");
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
        leg2->Draw("SAME");



        cc->cd(3);
        cc->cd(3)->SetLeftMargin(0.15);
        cc->cd(3)->SetBottomMargin(0.1);
        cc->cd(3)->SetRightMargin(0.05);
        cc->cd(3)->SetTopMargin(0.05);
        cc->cd(3)->SetTicky();
        cc->cd(3)->SetTickx();
        //
        h->Draw();
        for (int i = 0; i< percbinnumb; i++){
            hRatioIT1[i][2]->Draw("SAME");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 1");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.8, 0.86, "LHC22e1_18i");
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
        leg2->Draw("SAME");



        cc->cd(4);
        cc->cd(4)->SetLeftMargin(0.15);
        cc->cd(4)->SetBottomMargin(0.1);
        cc->cd(4)->SetRightMargin(0.05);
        cc->cd(4)->SetTopMargin(0.05);
        cc->cd(4)->SetTicky();
        cc->cd(4)->SetTickx();
        //
        h2->Draw();
        for (int i = 0; i< percbinnumb; i++){
            CorrFactorsIT1[i]->Draw("SAME LEP");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 1");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
        leg2->Draw("SAME");

        if (!DoMB){
        cc->SaveAs(Form("mcptshapeimages/RatioIT1_%sFixed%sin%03.0f_%03.0f.png", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        cc->SaveAs(Form("mcptshapeimages/RatioIT1_%sFixed%sin%03.0f_%03.0f.root", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        cc->SaveAs(Form("mcptshapeimages/RatioIT1_%sFixed%sin%03.0f_%03.0f.pdf", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        }
        else{
        cc->SaveAs(Form("mcptshapeimages/RatioIT1_%s_INELgt0.png", lCascType.Data()));
        cc->SaveAs(Form("mcptshapeimages/RatioIT1_%s_INELgt0.root", lCascType.Data()));
        cc->SaveAs(Form("mcptshapeimages/RatioIT1_%s_INELgt0.pdf", lCascType.Data()));
        }

        //IT2
        TCanvas* ccc = new TCanvas(Form("ccc"),"",1000,900);
        ccc->Divide(2,2);
        ccc->cd(1);
        ccc->cd(1)->SetLeftMargin(0.15);
        ccc->cd(1)->SetBottomMargin(0.1);
        ccc->cd(1)->SetRightMargin(0.05);
        ccc->cd(1)->SetTopMargin(0.05);
        ccc->cd(1)->SetTicky();
        ccc->cd(1)->SetTickx();
        //
        h->Draw();
        for (int i = 0; i< percbinnumb; i++){
            hRatioIT2[i][0]->Draw("SAME");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 2");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.8, 0.86, "LHC15g3b1");
        leg2->Draw("SAME");
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));


        ccc->cd(2);
        ccc->cd(2)->SetLeftMargin(0.15);
        ccc->cd(2)->SetBottomMargin(0.1);
        ccc->cd(2)->SetRightMargin(0.05);
        ccc->cd(2)->SetTopMargin(0.05);
        ccc->cd(2)->SetTicky();
        ccc->cd(2)->SetTickx();
        //
        h->Draw();
        for (int i = 0; i< percbinnumb; i++){
            hRatioIT2[i][1]->Draw("SAME");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 2");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.8, 0.86, "LHC22e1_17j");
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
        leg2->Draw("SAME");



        ccc->cd(3);
        ccc->cd(3)->SetLeftMargin(0.15);
        ccc->cd(3)->SetBottomMargin(0.1);
        ccc->cd(3)->SetRightMargin(0.05);
        ccc->cd(3)->SetTopMargin(0.05);
        ccc->cd(3)->SetTicky();
        ccc->cd(3)->SetTickx();
        //
        h->Draw();
        for (int i = 0; i< percbinnumb; i++){
            hRatioIT2[i][2]->Draw("SAME");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 2");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.8, 0.86, "LHC22e1_18i");
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
        leg2->Draw("SAME");



        ccc->cd(4);
        ccc->cd(4)->SetLeftMargin(0.15);
        ccc->cd(4)->SetBottomMargin(0.1);
        ccc->cd(4)->SetRightMargin(0.05);
        ccc->cd(4)->SetTopMargin(0.05);
        ccc->cd(4)->SetTicky();
        ccc->cd(4)->SetTickx();
        //
        h2->Draw();
        for (int i = 0; i< percbinnumb; i++){
            CorrFactorsIT2[i]->Draw("SAME LEP");
        }
        lab0->DrawLatex(0.3, 0.86, "Iteration 2");
        labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
        lab0->DrawLatex(0.7, 0.76, Form("%s fixed [%.0f-%0.f]",lWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
        leg2->Draw("SAME");

        if(!DoMB){
        ccc->SaveAs(Form("mcptshapeimages/RatioIT2_%sFixed%sin%03.0f_%03.0f.png", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        ccc->SaveAs(Form("mcptshapeimages/RatioIT2_%sFixed%sin%03.0f_%03.0f.root", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        ccc->SaveAs(Form("mcptshapeimages/RatioIT2_%sFixed%sin%03.0f_%03.0f.pdf", lCascType.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi));
        }
        else{
        ccc->SaveAs(Form("mcptshapeimages/RatioIT2_%s_INELgt0.png", lCascType.Data()));
        ccc->SaveAs(Form("mcptshapeimages/RatioIT2_%s_INELgt0.root", lCascType.Data()));
        ccc->SaveAs(Form("mcptshapeimages/RatioIT2_%s_INELgt0.pdf", lCascType.Data()));
        }
}


