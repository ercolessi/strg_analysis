void DoComparison(TString part = "Lambda", Bool_t isMB = kTRUE){

    TFile* myfile;
    if (!isMB) {
        myfile = TFile::Open(Form(Form("/home/fercoles/strg_analysis/PhDFinal/yields/%sYields_SPDClustersV0M_class0.root", part.Data())));
    } else {
        myfile = TFile::Open(Form("../yields/%sYields_INELgt0.root", part.Data()));
    }
    TGraphErrors *myyields = (TGraphErrors *)myfile->Get("YieldsNchStat");

    TFile *fiorellasfile = TFile::Open("HEPData-ins1748157-v1-root.root");
    TGraphAsymmErrors *fiorellasyields;
    TH1F *fiorellahist, *fiorellahist2;
    if (!isMB){
        if (part.Contains("K0Short")){
            fiorellasyields = (TGraphAsymmErrors *) fiorellasfile->Get("Table 8c/Graph1D_y1");
            fiorellahist = (TH1F *)fiorellasfile->Get("Table 8c/Hist1D_y1_e1");
            fiorellahist2 = (TH1F *)fiorellasfile->Get("Table 8c/Hist1D_y1_e2");
        }
        else if (part.Contains("Lambda")){
            fiorellasyields = (TGraphAsymmErrors *) fiorellasfile->Get("Table 8c/Graph1D_y2");
            fiorellahist = (TH1F *)fiorellasfile->Get("Table 8c/Hist1D_y2_e1");
            fiorellahist2 = (TH1F *)fiorellasfile->Get("Table 8c/Hist1D_y2_e2");
        }
        else if (part.Contains("Xi")){
            fiorellasyields = (TGraphAsymmErrors *) fiorellasfile->Get("Table 8c/Graph1D_y3");
            fiorellahist = (TH1F *)fiorellasfile->Get("Table 8c/Hist1D_y3_e1");
            fiorellahist2 = (TH1F *)fiorellasfile->Get("Table 8c/Hist1D_y3_e2");
        }
    } else {
        if (part.Contains("Xi")){
            fiorellasyields = (TGraphAsymmErrors *) fiorellasfile->Get("Table 11a/Graph1D_y3");
            fiorellahist = (TH1F *)fiorellasfile->Get("Table 11a/Hist1D_y3_e1");
            fiorellahist2 = (TH1F *)fiorellasfile->Get("Table 11a/Hist1D_y3_e2");
        }
        if (part.Contains("K0Short")){
            fiorellasyields = (TGraphAsymmErrors *)fiorellasfile->Get("Table 11a/Graph1D_y1");
            fiorellahist = (TH1F *)fiorellasfile->Get("Table 11a/Hist1D_y1_e1");
            fiorellahist2 = (TH1F *)fiorellasfile->Get("Table 11a/Hist1D_y1_e2");
        }
         if (part.Contains("Lambda")){
            fiorellasyields = (TGraphAsymmErrors *)fiorellasfile->Get("Table 11a/Graph1D_y2");
            fiorellahist = (TH1F *)fiorellasfile->Get("Table 11a/Hist1D_y2_e1");
            fiorellahist2 = (TH1F *)fiorellasfile->Get("Table 11a/Hist1D_y2_e2");
        }
    }

    Double_t* x = fiorellasyields->GetX();
    Double_t* xerr = fiorellasyields->GetEXlow();
    Double_t* yfiorella = fiorellasyields->GetY();
    Double_t* yerrfiorella2 = fiorellasyields->GetEYlow();
    Double_t *myy = myyields->GetY();
    Double_t *myyerr = myyields->GetEY();
    const int n = fiorellasyields->GetN();

    Double_t ratio[n], ratioerr[n];

    for (int i = 0; i < n; i++){
        double yerrfiorella = fiorellahist->GetBinContent(i+1);
        ratio[i] = yfiorella[i] / myy[i];
        ratioerr[i] = ratio[i] * TMath::Sqrt(TMath::Power(yerrfiorella2[i]/ yfiorella[i], 2) + TMath::Power(myyerr[i] / myy[i], 2));
        cout << yfiorella[i] - yerrfiorella << endl;
    }



    TGraphErrors *gr = new TGraphErrors(n, x, ratio, xerr, ratioerr);
    gr->GetYaxis()->SetRangeUser(0.8,1.2);
    gr->GetXaxis()->SetRangeUser(0,35);
    gr->GetYaxis()->SetTitle("dN/dy_{Published}/dN/dy_{ThisAnalysis}");
    gr->GetXaxis()->SetTitle("n_{ch}");
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

    TCanvas* s = new TCanvas("s", "s", 800, 600);
    fiorellasyields->SetMarkerStyle(20);
    fiorellasyields->SetMarkerColor(kBlue);
    fiorellasyields->SetLineColor(kBlue);
    myyields->SetMarkerStyle(kOpenCircle);
    myyields->SetMarkerColor(kRed);
    myyields->SetLineColor(kRed);
    fiorellasyields->Draw("AEP");
    myyields->Draw("SAME EP");
    fiorellasyields->GetYaxis()->SetTitle("dN/dy K^{0}_{S}");
    fiorellasyields->GetXaxis()->SetTitle("n_{ch}");

        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(myyields, "This Analysis", "LEP");
    leg->AddEntry(fiorellasyields, "Published stat + syst", "LEP");
    leg->SetBorderSize(0);
    leg->Draw();
}