void canvasprep(TCanvas* c);
void histoprep(TH1* h);

void drawcal_ZNZPdist(){

    TFile *f = new TFile("perfSel.root");
    // perfSel.root");
    TH2D *hRunZDCNC = (TH2D *)f->Get("hRunZDCdist");
    //"hRunZDCdist");
    //"hRunZDCNC");

    const int nruns = hRunZDCNC->GetNbinsX()-1;
    TH1D* projYZDCNC[nruns];

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    canvasprep(c2);

    for (int i = 1; i<=nruns; i++){
        projYZDCNC[i] = hRunZDCNC->ProjectionY(Form("projYZDCNC_%s", hRunZDCNC->GetXaxis()->GetBinLabel(i + 1)), i + 1, i + 1);
        projYZDCNC[i]->SetTitle(Form("Run %s",hRunZDCNC->GetXaxis()->GetBinLabel(i+1)));
        projYZDCNC[i]->SetStats(0);
        projYZDCNC[i]->DrawNormalized("SAME");
    }

    TFile *f1 = new TFile("perfSel.root");
    TH1D *h0PeriodZDC = (TH1D *)f1->Get("h0PeriodZDC");
    TH1D *h1PeriodZDC = (TH1D *)f1->Get("h1PeriodZDC");
    TH1D *h2PeriodZDC = (TH1D *)f1->Get("h2PeriodZDC");

    h0PeriodZDC->SetLineColor(kRed);
    h1PeriodZDC->SetLineColor(kBlue);
    h2PeriodZDC->SetLineColor(kGreen);

    TCanvas* c = new TCanvas();
    c->SetLogy();
    h0PeriodZDC->Scale(1./h0PeriodZDC->Integral());
    h1PeriodZDC->Scale(1./h1PeriodZDC->Integral());
    h2PeriodZDC->Scale(1./h2PeriodZDC->Integral());
    h0PeriodZDC->Draw();
    h1PeriodZDC->Draw("same");
    h2PeriodZDC->Draw("same");

    cout << "h0PeriodZDC->GetMean() = " << h0PeriodZDC->GetMean() << endl;
    cout << "h1PeriodZDC->GetMean() = " << h1PeriodZDC->GetMean() << endl;
    cout << "h2PeriodZDC->GetMean() = " << h2PeriodZDC->GetMean() << endl;

    TH1D* h1ratioto0 = (TH1D*)h1PeriodZDC->Clone("hratioto0");
    h1ratioto0->Divide(h0PeriodZDC);
    h1ratioto0->SetLineColor(kBlue);
    TH1D* h2ratioto0 = (TH1D*)h2PeriodZDC->Clone("hratioto0");
    h2ratioto0->Divide(h0PeriodZDC);
    h2ratioto0->SetLineColor(kGreen);

    new TCanvas;
    h1ratioto0->Draw();
    h2ratioto0->Draw("same");
}

void canvasprep(TCanvas* c){
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.1);
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.15);
    c->SetTicks();
}

void histoprep(TH1* h){
    h->SetStats(0);
    h->SetMarkerStyle(kFullSquare);
    h->SetMarkerSize(1.5);
    h->SetMarkerColor(kRed);
    h->SetLineColor(kRed);
    h->SetTitle("");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetLabelSize(0.05);
}