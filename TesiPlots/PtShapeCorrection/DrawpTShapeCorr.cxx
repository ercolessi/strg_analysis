// classes
enum classname
{
    kStandalone = 0,
    kHighMult,
    kLowMult,
    kHighZN,
    kLowZN
};

void DrawpTShapeCorr(
    TString lCascType = "Xi",
    Int_t lClassCode = kStandalone,
    Bool_t DoMB = kFALSE)
{

    int convertclass = 0;
    if (lClassCode == 1)
        convertclass = kHighMult;
    if (lClassCode == 2)
        convertclass = kLowMult;
    if (lClassCode == 4)
        convertclass = kLowZN;
    if (lClassCode == 5)
        convertclass = kHighZN;

    TFile* file;
    if (!DoMB){
        file = TFile::Open(Form("~/strg_analysis/PhDFinal/correctedspectra/CorrSpectra-%s-13TeV_SPDClusters_V0M_class%i.root",
                                lCascType.Data(), convertclass));
    } else {
        file = TFile::Open(Form("~/strg_analysis/PhDFinal/correctedspectra/CorrSpectra-%s-13TeV_INELgt0.root",
                                lCascType.Data()));
    }

    Int_t nbins;
    Int_t *colors;
    TString label[] = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"};
    TString lclassname[] = {"Standalone", "High Multiplicity", "Low Multiplicity", "High ZN", "Low ZN"};
    Int_t colors5[] = {kRed + 2, kOrange + 1, kSpring - 1, kAzure + 7, kBlue + 2};
    Int_t colors10[] = {kRed + 2, kRed, kOrange + 7, kYellow + 1, kSpring - 1, kGreen + 2, kTeal, kAzure + 7, kBlue, kBlue + 2};
    Int_t colors7[] = {kRed + 2,  kOrange + 7, kYellow + 1, kSpring - 1, kAzure + 7, kBlue, kBlue + 2};

    if (lClassCode == kStandalone)
    {
        colors = colors10;
        nbins = 10;
    }
    else if (lClassCode == kHighMult || lClassCode == kLowMult)
    {
        colors = colors7;
        nbins = 7;
    }
    else if (lClassCode == kLowZN)
    {
        colors = colors5;
        nbins = 4;
    }
    else if (lClassCode == kHighZN)
    {
        colors = colors5;
        nbins = 5;
    }
    const int percbinnumb = nbins;

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
    lab0->SetTextSize(0.05);
    lab0->SetTextAlign(22);
    lab0->SetTextAngle(0);
    TLatex *labbig = new TLatex();
    labbig->SetTextFont(42);
    labbig->SetNDC();
    labbig->SetTextColor(1);
    labbig->SetTextSize(0.08);
    labbig->SetTextAlign(22);
    labbig->SetTextAngle(0);
    TLatex *tex = new TLatex(0.2, 0.9, "ALICE, pp #sqrt{s} = 13 TeV");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    TLatex *tex2 = new TLatex(0.2, 0.86, "This work");
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(2);

    TCanvas* c1 = new TCanvas(Form("c1"),"",800,1100);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.1);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.05);
    c1->SetTicky();
    c1->SetTickx();
    c1->SetLogy();


    TH1D* h1 = new TH1D("h1", "",7, 0.6,6.5);
    h1->GetYaxis()->SetRangeUser(3E+3,5E+8);
    h1->GetYaxis()->SetTitle("dN/dp_{T} (a.u.)");
    h1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleOffset(1.1);
    h1->GetXaxis()->SetTitleOffset(0.9);
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
    labbig->DrawLatex(0.83, 0.85, Form("#%s",lCascType.Data()));
    lab0->DrawLatex(0.78, 0.76, Form("%s", lclassname[lClassCode].Data()));

    //
    TLegend* leg = new TLegend(0.18,0.3,0.45,0.4);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->AddEntry(hGenMCIT0[0],"MC input shape anch. to 2015","L");
    leg->AddEntry(hGenMCIT0[1],"MC input shape anch. to 2017","L");
    leg->AddEntry(hGenMCIT0[2],"MC input shape anch. to 2018","L");
    leg->Draw("SAME");

    TLegend *lede = new TLegend(0.21, 0.12, 0.5, 0.28);
    lede->SetTextSize(0.03);
    lede->SetBorderSize(0);
    lede->SetNColumns(2);
    for (int i = 0; i < percbinnumb; i++)
    {
        lede->AddEntry(hGenSpectraIT0[i], Form("%s", label[i].Data()), "LEP");
    }
    lede->Draw("SAME");

    tex->Draw();
    tex2->Draw();

    TLegend* leg2 = new TLegend(0.2,0.22,0.4,0.32);
    leg2->SetTextSize(0.05);
    leg2->SetBorderSize(0);
    leg2->AddEntry(hGenSpectraIT0[0], Form("%s", label[0].Data()), "LEP");
    leg2->AddEntry(hGenSpectraIT0[percbinnumb-1],Form("%s", label[percbinnumb-1].Data()),"LEP");

    if (!DoMB){
    c1->SaveAs(Form("GenSpectraIT0_%s_%i.png", lCascType.Data(), lClassCode));
    c1->SaveAs(Form("GenSpectraIT0_%s_%i.root", lCascType.Data(), lClassCode));
    c1->SaveAs(Form("GenSpectraIT0_%s_%i.pdf", lCascType.Data(), lClassCode));
    } else {
    c1->SaveAs(Form("GenSpectraIT0_%s_INELgt0.png", lCascType.Data()));
    c1->SaveAs(Form("GenSpectraIT0_%s_INELgt0.root", lCascType.Data()));
    c1->SaveAs(Form("GenSpectraIT0_%s_INELgt0.pdf", lCascType.Data()));
    }

    TCanvas* c = new TCanvas(Form("c"),"",1000,900);
    c->Divide(2,2);
    c->cd(1);
    c->cd(1)->SetLeftMargin(0.15);
    c->cd(1)->SetBottomMargin(0.15);
    c->cd(1)->SetRightMargin(0.05);
    c->cd(1)->SetTopMargin(0.08);
    c->cd(1)->SetTicky();
    c->cd(1)->SetTickx();
    //
    TH1D* h = new TH1D("h", "",7, 0.6,6.5);
    h->GetYaxis()->SetRangeUser(0., 2.5);
    h->GetYaxis()->SetTitle("Corr.Spectrum(DATA) / Input(MC)");
    h->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetRangeUser(0.6, 6.5);
    h->SetStats(0);
    h->SetTitle("");
    h->Draw();
    for (int i = 0; i< percbinnumb; i++){
        hRatio[i][0]->Draw("SAME");
    }
    lab0->DrawLatex(0.3, 0.86, "Iteration 0");
    labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
    lab0->DrawLatex(0.72, 0.86, "MC anch. to 2015");
    lab0->DrawLatex(0.7, 0.78, Form("%s", lclassname[lClassCode].Data()));
    leg2->Draw("SAME");


    c->cd(2);
    c->cd(2)->SetLeftMargin(0.15);
    c->cd(2)->SetBottomMargin(0.15);
    c->cd(2)->SetRightMargin(0.05);
    c->cd(2)->SetTopMargin(0.08);
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
    lab0->DrawLatex(0.72, 0.86, "MC anch. to 2017");
    lab0->DrawLatex(0.7, 0.78, Form("%s",lclassname[lClassCode].Data()));
    leg2->Draw("SAME");


    c->cd(3);
    c->cd(3)->SetLeftMargin(0.15);
    c->cd(3)->SetBottomMargin(0.15);
    c->cd(3)->SetRightMargin(0.05);
    c->cd(3)->SetTopMargin(0.08);
    c->cd(3)->SetTicky();
    c->cd(3)->SetTickx();
    //
    hh->Draw();
    for (int i = 0; i< percbinnumb; i++){
        hRatio[i][2]->Draw("SAME");
    }
    lab0->DrawLatex(0.3, 0.86, "Iteration 0");
    labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
    lab0->DrawLatex(0.72, 0.86, "MC anch. to 2018");
    lab0->DrawLatex(0.7, 0.78, Form("%s",lclassname[lClassCode].Data()));
    leg2->Draw("SAME");


    c->cd(4);
    c->cd(4)->SetLeftMargin(0.15);
    c->cd(4)->SetBottomMargin(0.15);
    c->cd(4)->SetRightMargin(0.05);
    c->cd(4)->SetTopMargin(0.08);
    c->cd(4)->SetTicky();
    c->cd(4)->SetTickx();
    //
    TH1D* h2 = new TH1D("h2", "",7, 0.6,6.5);
    h2->GetYaxis()->SetRangeUser(0.85, 1.15);
    h2->GetYaxis()->SetTitle("Corr. factor (#varepsilon_{it}/#varepsilon_{it-1})");
    h2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h2->GetYaxis()->SetTitleOffset(1.1);
    h2->GetXaxis()->SetTitleOffset(1.1);
    h2->GetYaxis()->SetTitleSize(0.055);
    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetXaxis()->SetRangeUser(0.6, 6.5);
    h2->SetStats(0);
    gStyle->SetTitleSize(0.06,"t");
    h2->SetTitle("Total correction in this iteration");
    h2->Draw();
    for (int i = 0; i< percbinnumb; i++){
        CorrFactors[i]->Draw("SAME LEP");
    }
    lab0->DrawLatex(0.3, 0.86, "Iteration 0");
    labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
    //lab0->DrawLatex(0.7, 0.78, Form("%s",lclassname[lClassCode].Data()));
    leg2->Draw("SAME");

    TLatex *texx = new TLatex(0.5, 0.83, "ALICE, pp #sqrt{s} = 13 TeV");
    texx->SetNDC();
    texx->SetTextFont(42);
    texx->SetTextSize(0.05);
    texx->SetLineWidth(2);
    TLatex *tex2x = new TLatex(0.6, 0.75, "This work");
    tex2x->SetNDC();
    tex2x->SetTextFont(42);
    tex2x->SetTextSize(0.05);
    tex2x->SetLineWidth(2);

    texx->Draw();
    tex2x->Draw();

    if (!DoMB){
    c->SaveAs(Form("RatioIT0_%s_class%i.png", lCascType.Data(), lClassCode));
    c->SaveAs(Form("RatioIT0_%s_class%i.root", lCascType.Data(), lClassCode));
    c->SaveAs(Form("RatioIT0_%s_class%i.pdf", lCascType.Data(), lClassCode));
    } else {
    c->SaveAs(Form("RatioIT0_%s_INELgt0.png", lCascType.Data()));
    c->SaveAs(Form("RatioIT0_%s_INELgt0.root", lCascType.Data()));
    c->SaveAs(Form("RatioIT0_%s_INELgt0.pdf", lCascType.Data()));
    }

    //IT1
    TCanvas* cc = new TCanvas(Form("cc"),"",1000,900);
    cc->Divide(2,2);
    cc->cd(1);
    cc->cd(1)->SetLeftMargin(0.15);
    cc->cd(1)->SetBottomMargin(0.15);
    cc->cd(1)->SetRightMargin(0.05);
    cc->cd(1)->SetTopMargin(0.08);
    cc->cd(1)->SetTicky();
    cc->cd(1)->SetTickx();
    //
    h->Draw();
    for (int i = 0; i< percbinnumb; i++){
        hRatioIT1[i][0]->Draw("SAME");
    }
    lab0->DrawLatex(0.3, 0.86, "Iteration 1");
    labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
    lab0->DrawLatex(0.72, 0.86, "MC anch. to 2015");
    lab0->DrawLatex(0.7, 0.78, Form("%s", lclassname[lClassCode].Data()));
    leg2->Draw("SAME");


    cc->cd(2);
    cc->cd(2)->SetLeftMargin(0.15);
    cc->cd(2)->SetBottomMargin(0.15);
    cc->cd(2)->SetRightMargin(0.05);
    cc->cd(2)->SetTopMargin(0.08);
    cc->cd(2)->SetTicky();
    cc->cd(2)->SetTickx();
    //
    h->Draw();
    for (int i = 0; i< percbinnumb; i++){
        hRatioIT1[i][1]->Draw("SAME");
    }
    lab0->DrawLatex(0.3, 0.86, "Iteration 1");
    labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
    lab0->DrawLatex(0.72, 0.86, "MC anch. to 2017");
    lab0->DrawLatex(0.7, 0.78, Form("%s", lclassname[lClassCode].Data()));
    leg2->Draw("SAME");



    cc->cd(3);
    cc->cd(3)->SetLeftMargin(0.15);
    cc->cd(3)->SetBottomMargin(0.15);
    cc->cd(3)->SetRightMargin(0.05);
    cc->cd(3)->SetTopMargin(0.08);
    cc->cd(3)->SetTicky();
    cc->cd(3)->SetTickx();
    //
    h->Draw();
    for (int i = 0; i< percbinnumb; i++){
        hRatioIT1[i][2]->Draw("SAME");
    }
    lab0->DrawLatex(0.3, 0.86, "Iteration 1");
    labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
    lab0->DrawLatex(0.72, 0.86, "MC anch. to 2018");
    lab0->DrawLatex(0.7, 0.78, Form("%s", lclassname[lClassCode].Data()));
    leg2->Draw("SAME");



    cc->cd(4);
    cc->cd(4)->SetLeftMargin(0.15);
    cc->cd(4)->SetBottomMargin(0.15);
    cc->cd(4)->SetRightMargin(0.05);
    cc->cd(4)->SetTopMargin(0.08);
    cc->cd(4)->SetTicky();
    cc->cd(4)->SetTickx();
    //
    h2->Draw();
    for (int i = 0; i< percbinnumb; i++){
        CorrFactorsIT1[i]->Draw("SAME LEP");
    }
    lab0->DrawLatex(0.3, 0.86, "Iteration 1");
    labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
    //lab0->DrawLatex(0.7, 0.73, Form("%s", lclassname[lClassCode].Data()));
    leg2->Draw("SAME");

    texx->Draw();
    tex2x->Draw();

    if (!DoMB){
    cc->SaveAs(Form("RatioIT1_%s_%i.png", lCascType.Data(), lClassCode));
    cc->SaveAs(Form("RatioIT1_%s_%i.root", lCascType.Data(), lClassCode));
    cc->SaveAs(Form("RatioIT1_%s_%i.pdf", lCascType.Data(), lClassCode));
    }
    else{
    cc->SaveAs(Form("RatioIT1_%s_INELgt0.png", lCascType.Data()));
    cc->SaveAs(Form("RatioIT1_%s_INELgt0.root", lCascType.Data()));
    cc->SaveAs(Form("RatioIT1_%s_INELgt0.pdf", lCascType.Data()));
    }

    //IT2
    TCanvas* ccc = new TCanvas(Form("ccc"),"",1000,900);
    ccc->Divide(2,2);
    ccc->cd(1);
    ccc->cd(1)->SetLeftMargin(0.15);
    ccc->cd(1)->SetBottomMargin(0.15);
    ccc->cd(1)->SetRightMargin(0.05);
    ccc->cd(1)->SetTopMargin(0.08);
    ccc->cd(1)->SetTicky();
    ccc->cd(1)->SetTickx();
    //
    h->Draw();
    for (int i = 0; i< percbinnumb; i++){
        hRatioIT2[i][0]->Draw("SAME");
    }
    lab0->DrawLatex(0.3, 0.86, "Iteration 2");
    labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
    lab0->DrawLatex(0.72, 0.86, "MC anch. to 2015");
    leg2->Draw("SAME");
    lab0->DrawLatex(0.7, 0.78, Form("%s", lclassname[lClassCode].Data()));

    ccc->cd(2);
    ccc->cd(2)->SetLeftMargin(0.15);
    ccc->cd(2)->SetBottomMargin(0.15);
    ccc->cd(2)->SetRightMargin(0.05);
    ccc->cd(2)->SetTopMargin(0.08);
    ccc->cd(2)->SetTicky();
    ccc->cd(2)->SetTickx();
    //
    h->Draw();
    for (int i = 0; i< percbinnumb; i++){
        hRatioIT2[i][1]->Draw("SAME");
    }
    lab0->DrawLatex(0.3, 0.86, "Iteration 2");
    labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
    lab0->DrawLatex(0.72, 0.86, "MC anch. to 2017");
    lab0->DrawLatex(0.7, 0.78, Form("%s", lclassname[lClassCode].Data()));
    leg2->Draw("SAME");



    ccc->cd(3);
    ccc->cd(3)->SetLeftMargin(0.15);
    ccc->cd(3)->SetBottomMargin(0.15);
    ccc->cd(3)->SetRightMargin(0.05);
    ccc->cd(3)->SetTopMargin(0.08);
    ccc->cd(3)->SetTicky();
    ccc->cd(3)->SetTickx();
    //
    h->Draw();
    for (int i = 0; i< percbinnumb; i++){
        hRatioIT2[i][2]->Draw("SAME");
    }
    lab0->DrawLatex(0.3, 0.86, "Iteration 2");
    labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
    lab0->DrawLatex(0.72, 0.86, "MC anch. to 2018");
    lab0->DrawLatex(0.7, 0.78, Form("%s", lclassname[lClassCode].Data()));
    leg2->Draw("SAME");



    ccc->cd(4);
    ccc->cd(4)->SetLeftMargin(0.15);
    ccc->cd(4)->SetBottomMargin(0.15);
    ccc->cd(4)->SetRightMargin(0.05);
    ccc->cd(4)->SetTopMargin(0.08);
    ccc->cd(4)->SetTicky();
    ccc->cd(4)->SetTickx();
    //
    h2->Draw();
    for (int i = 0; i< percbinnumb; i++){
        CorrFactorsIT2[i]->Draw("SAME LEP");
    }
    lab0->DrawLatex(0.3, 0.86, "Iteration 2");
    labbig->DrawLatex(0.3, 0.79, Form("#%s",lCascType.Data()));
    //lab0->DrawLatex(0.7, 0.78, Form("%s", lclassname[lClassCode].Data()));
    leg2->Draw("SAME");

    texx->Draw();
    tex2x->Draw();

    if(!DoMB){
    ccc->SaveAs(Form("RatioIT2_%s_%i.png", lCascType.Data(), lClassCode));
    ccc->SaveAs(Form("RatioIT2_%s_%i.root", lCascType.Data(), lClassCode));
    ccc->SaveAs(Form("RatioIT2_%s_%i.pdf", lCascType.Data(), lClassCode));
    }
    else{
    ccc->SaveAs(Form("RatioIT2_%s_INELgt0.png", lCascType.Data()));
    ccc->SaveAs(Form("RatioIT2_%s_INELgt0.root", lCascType.Data()));
    ccc->SaveAs(Form("RatioIT2_%s_INELgt0.pdf", lCascType.Data()));
    }
}
