void DoComparisonSpectra(TString part = "Lambda"){

    TFile *myfile = TFile::Open("../correctedspectra/CorrSpectra-Lambda-13TeV_INELgt0.root");
    TH1D *myyields = (TH1D *)myfile->Get(Form("FinalSpectra/PtSpectrumCorrStat_%s_%.0f_%.0f_%s_%.0f_%.0f", "SPDClusters", 0., 100., "V0M", 0., 100.));

    TFile *fiorellasfile = TFile::Open("HEPData-ins1748157-v1-root.root");
    TGraphAsymmErrors *fiorellasyields;
    TH1F *fiorellahist, *fiorellahist2;
    if (part.Contains("K0Short")){
        fiorellasyields = (TGraphAsymmErrors *) fiorellasfile->Get("Table 4/Graph1D_y11");
        fiorellahist = (TH1F *)fiorellasfile->Get("Table 4/Hist1D_y11_e1");
        fiorellahist2 = (TH1F *)fiorellasfile->Get("Table 4/Hist1D_y11_e2");
    }
    else if (part.Contains("Lambda")){
        fiorellasyields = (TGraphAsymmErrors *) fiorellasfile->Get("Table 2/Graph1D_y11");
        fiorellahist = (TH1F *)fiorellasfile->Get("Table 2/Hist1D_y11_e1");
        fiorellahist2 = (TH1F *)fiorellasfile->Get("Table 2/Hist1D_y11_e2");
    }
    else if (part.Contains("Xi")){
        fiorellasyields = (TGraphAsymmErrors *) fiorellasfile->Get("Table 3/Graph1D_y11");
        fiorellahist = (TH1F *)fiorellasfile->Get("Table 3/Hist1D_y11_e1");
        fiorellahist2 = (TH1F *)fiorellasfile->Get("Table 3/Hist1D_y11_e2");
    }

    Double_t* x = fiorellasyields->GetX();
    Double_t* xerr = fiorellasyields->GetEXlow();
    Double_t* yfiorella = fiorellasyields->GetY();
    Double_t* yerrfiorella2 = fiorellasyields->GetEYlow();
    const int n = fiorellasyields->GetN();
    Double_t myy[n];
    Double_t myyerr[n];

    Double_t ratio[n], ratioerr[n];

    for (int i = 0; i < n; i++){
        double yerrfiorella = fiorellahist->GetBinContent(i+1);
        myy[i] = myyields->GetBinContent(i+1);
        myyerr[i] = myyields->GetBinError(i+1);
        ratio[i] = yfiorella[i] / myy[i];
        ratioerr[i] = ratio[i] * TMath::Sqrt(TMath::Power(yerrfiorella / yfiorella[i], 2) + TMath::Power(myyerr[i] / myy[i], 2));
    }

    new TCanvas;
    fiorellasyields->Draw();
    myyields->Draw("same");

    TGraphErrors *gr = new TGraphErrors(n, x, ratio, xerr, ratioerr);
    gr->GetYaxis()->SetRangeUser(0.8,1.2);
    gr->GetXaxis()->SetRangeUser(0,35);
    gr->GetYaxis()->SetTitle("dN/dy_{Published}/dN/dy_{ThisAnalysis}");
    gr->GetXaxis()->SetTitle("pT");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor(kRed);
    gr->SetMarkerSize(1.5);
    gr->SetLineWidth(2);
    gr->SetTitle("");
    gr->GetYaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleOffset(0.9);
    gr->GetXaxis()->SetTitleSize(0.05);

    TCanvas *c = new TCanvas("c","c",800,600);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.15);
    c->SetTopMargin(0.05);
    c->SetGridy();
    c->SetTickx();
    c->SetTicky();
    gr->Draw("AP");
    c->SaveAs(Form("Comparison%s.png",part.Data()));

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.07);
    tex->SetTextFont(42);
    if (part.Contains("K0Short"))
        tex->DrawLatex(0.2, 0.85, "K^{0}_{S}");
    else if (part.Contains("Lambda"))
        tex->DrawLatex(0.2, 0.85, "#Lambda + #bar{#Lambda}");
    else if (part.Contains("Xi"))
        tex->DrawLatex(0.2, 0.85, "#Xi^{-} + #bar{#Xi}^{+}");
}