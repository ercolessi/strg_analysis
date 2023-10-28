void macrostupida(){
    Double_t ptbinlimitsK0s[] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3.0,
                                 3.3, 3.6, 3.9, 4.2, 4.6, 5, 5.4, 5.9, 6.5, 7, 7.5, 8, 8.5, 9.2, 10, 11, 12, 13.5, 15};
    Long_t ptbinnumbK0s = sizeof(ptbinlimitsK0s) / sizeof(Double_t) - 1;

    TFile *write = new TFile("SystematicsPileUp-K0Short.root", "RECREATE");
    TH1F *hsystPileUp = new TH1F("hsystPileUp", "hsystPileUp", ptbinnumbK0s, ptbinlimitsK0s);
    for (int bin = 1; bin < ptbinnumbK0s; bin++){
        if (ptbinlimitsK0s[bin-1]<1.5){
            hsystPileUp->SetBinContent(bin, 0.016);
        } else if (ptbinlimitsK0s[bin-1]>=1.5) {
            hsystPileUp->SetBinContent(bin, 0.025);
        }
    }
    hsystPileUp->Write();

    /*Double_t ptbinlimitsK0s[] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3.0,
                                 3.3, 3.6, 3.9, 4.2, 4.6, 5, 5.4, 5.9, 6.5, 7, 7.5, 8, 8.5, 9.2, 10, 11, 12, 13.5, 15};
    Long_t ptbinnumbK0s = sizeof(ptbinlimitsK0s) / sizeof(Double_t) - 1;

    TFile *write = new TFile("systMatBudget-K0Short.root", "RECREATE");
    TH1F *hMatBudget = new TH1F("hMatBudget", "hMatBudget", ptbinnumbK0s, ptbinlimitsK0s);
    for (int bin = 1; bin < ptbinnumbK0s; bin++){
        if (ptbinlimitsK0s[bin-1]<1.5){
            hMatBudget->SetBinContent(bin, 0.011);
        } else if (ptbinlimitsK0s[bin-1]>=1.5) {
            hMatBudget->SetBinContent(bin, 0.005);
        }
    }
    hMatBudget->Write();*/
 }