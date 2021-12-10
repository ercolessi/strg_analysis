double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) ;
void canvas_hestetics(TCanvas* c);


void PlotOfficialYieldTrend(){
    
    //TFile* filemb = TFile::Open("YieldsSPDClusters_SelectedWithZDC0100_17j.root");
    //TFile* filesel = TFile::Open("YieldsZDC_SelectedWithSPDClusters1015_17j.root");

    TFile* filemb = TFile::Open("../FullStatistics/YieldsSPDClusters_SelectedWithZDC0100_18i.root");
    TFile* filesel = TFile::Open("../FullStatistics/YieldsV0M_SelectedWithSPDClusters1030_18i.root");


    TGraphErrors* gmb = (TGraphErrors*)filemb->Get("YieldsNchStat");
    TGraphErrors* gselperc = (TGraphErrors*)filesel->Get("AvYieldsPercStat");
    TGraphErrors* gsel = (TGraphErrors*)filesel->Get("AvYieldsNchStat");

    double* xmb = gmb->GetX();
    double* ymb = gmb->GetY();
    double* xsel = gsel->GetX();
    double* ysel = gsel->GetY();

    double yvalmb[2], xvalsel[2];
    yvalmb[0] = gmb->Eval(xsel[0])/xsel[0];
    yvalmb[1] = gmb->Eval(xsel[7])/xsel[7];
    xvalsel[0] = 0.;
    xvalsel[1] =100.;

    TCanvas* c = new TCanvas("c","",1800,1500);
    canvas_hestetics(c);
    
    TGraph* g = new TGraph(2,xvalsel, yvalmb);
    gselperc->SetLineColor(kBlue);
    gselperc->SetMarkerColor(kBlue);
    gselperc->SetMarkerStyle(21);
    gselperc->SetMarkerSize(2.2);
    gselperc->Draw("AEP");
    gselperc->SetLineWidth(2);
    gselperc->GetXaxis()->SetRangeUser(0.,100.);
    //gselperc->GetYaxis()->SetRangeUser(0.048,.067);
    gselperc->GetYaxis()->SetTitleOffset(1.9);
    g->SetLineWidth(2);
    g->Draw("SAME L");

    new TCanvas;
    gmb->Draw();
    new TCanvas;
    gsel->Draw();
    
    TFile* Write = new TFile("OffYieldsNch.root", "recreate");

}

//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop + errorfrombottom) );
    }
    return 1.;
}

void canvas_hestetics(TCanvas* c) {
    c->SetLeftMargin(0.22);
    c->SetBottomMargin(0.16);
    c->SetRightMargin(0.12);
    c->SetTopMargin(0.12);
    c->SetTicky();
    c->SetTickx();
}

