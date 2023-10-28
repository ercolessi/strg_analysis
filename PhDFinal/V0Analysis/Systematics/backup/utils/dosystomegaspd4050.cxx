void dosystomegaspd4050(){
    TFile* f = TFile::Open(Form("FinalSystematicsYield-%s-fixedSPDClusters_%03.0f_%03.0f.root","Lambda",0., 100.));
    TH1F* h = (TH1F*)f->Get("hSystTot");
    h->Draw();

    Double_t perc[] ={0,5,10,15,20,30,40,50,70,100};
    const int n = sizeof(perc)/sizeof(Double_t) - 1;
    TH1F* hnew = new TH1F("hnew","",n,perc);

    for(int i = 1; i<= h->GetNbinsX(); i++){
        hnew->SetBinContent(i,h->GetBinContent(i));
    }
    

    Double_t perc2[] =  {0,20,30,40,50,60,70,100};
    const int n2 = sizeof(perc2)/sizeof(Double_t) - 1;
    TH1F* hnew2 = new TH1F("hnew2","",n2,perc2);

    for(int i = 1; i<= hnew2->GetNbinsX(); i++){
        cout << (perc2[i]+perc2[i-1])/2 << endl;
        hnew2->SetBinContent(i,hnew->Interpolate((perc2[i]+perc2[i-1])/2));
    }
    hnew->Draw();
    hnew2->SetLineColor(kRed);
    hnew2->Draw("SAME");
    hnew2->SetName("hSystTot");

    TFile* f1 = new TFile(Form("FinalSystematicsYield-%s-fixedSPDClusters_%03.0f_%03.0f.root","Lambda",40., 50.),"RECREATE");
    hnew2->Write();

}