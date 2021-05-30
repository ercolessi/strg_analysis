double domean(TH1D* h1, int bin);

void DoMobileMean(){
    TFile* f = TFile::Open("ZDCSystMobileMeanPos.root");
    TH1D* h = (TH1D*) f->Get("h");
    TH1D* hc = (TH1D*) h->Clone("hc");
    hc->Reset();
    for (int i = 1; i <= h->GetNbinsX(); i++){
        hc->SetBinContent(i, domean(h,i));
    }

    new TCanvas;
    hc->SetLineColor(kRed);
    h->Draw("SAME");
    hc->Draw("SAME");

    TFile* w = new TFile("SystMobileMeanforZDCStandalonePos.root","RECREATE");
    hc->Write();

}

double domean(TH1D* h1, int bin){
   
    if (bin == 1) {
        return (h1->GetBinContent(bin)*h1->GetBinWidth(bin) + h1->GetBinContent(bin+1)*h1->GetBinWidth(bin+1))
                /
               (h1->GetBinWidth(bin) + h1->GetBinWidth(bin+1));
            }
    else if (bin == h1->GetNbinsX()) {
        return (h1->GetBinContent(bin)*h1->GetBinWidth(bin)+h1->GetBinContent(bin-1)*h1->GetBinWidth(bin-1))
        /
        (h1->GetBinWidth(bin)+h1->GetBinWidth(bin-1));}
    else {return 
    (h1->GetBinContent(bin)*h1->GetBinWidth(bin)+h1->GetBinContent(bin+1)*h1->GetBinWidth(bin+1)+h1->GetBinContent(bin-1)*h1->GetBinWidth(bin-1))
    /
    (h1->GetBinWidth(bin-1)+h1->GetBinWidth(bin)+h1->GetBinWidth(bin+1));}
}