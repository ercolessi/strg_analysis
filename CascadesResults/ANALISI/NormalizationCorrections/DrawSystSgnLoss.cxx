void DrawSystSgnLoss(){

    TFile* file = new TFile ("SgnLossSyst.root");

    Double_t percentileV0[] = {0.,1.,5,10,15,20,30,40,50,70,100};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    //
    Double_t percentileZDC[] = {0,20,30,40,50,60,70,80,90,100};
    const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;

    Double_t percentileZDCMerged[nbinZDC] = {0,20,40,50,60,70,80,90,100};

     Double_t percentileV0Merged[nbinV0] = {0.,5,10,15,20,30,40,50,70,100};
    

    TCanvas* c = new TCanvas("c","",1300,1100);
    c->SetLeftMargin(0.19);
    c->SetRightMargin(0.12);
    c->SetBottomMargin(0.12);
    TH1D* h[nbinZDC];
    TLegend* l = new TLegend(0.57,0.7,0.87,0.89);
    l->SetTextSize(0.021);   
    l->SetBorderSize(0);
    for (int i =0; i<nbinZDC;i++){
        h[i] = (TH1D*)file->Get(Form("%s-%s/hRatioClone%i","EEsel","ZDC",i));
        h[i]->GetYaxis()->SetTitle("#frac{  #varepsilon_{part}^{Pythia62011-LHC} }{#varepsilon_{part}^{EPOS-LHC}}");
        h[i]->GetYaxis()->SetTitleOffset(2.1);
        h[i]->GetYaxis()->SetRangeUser(0.9,1.1);
        h[i]->GetXaxis()->SetRangeUser(0.6,6.5);
        h[i]->SetStats(0);
        h[i]->SetTitle("");
        h[i]->SetLineWidth(2);
        h[i]->Draw("HIST SAME");
        l->AddEntry(h[i],Form("ZDC %.0f-%.0f ",percentileZDC[i],percentileZDC[i+1]),"L");
    }
   
    l->Draw("SAME");  

    c->SaveAs("ZDCSgnLossSystematic.png");

    TText *xlabel = new TText();
    xlabel-> SetNDC();
    xlabel-> SetTextFont(42);
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.025);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);

    TCanvas* c1 = new TCanvas("c1","",1300,1100);
    c1->SetLeftMargin(0.19);
    c1->SetRightMargin(0.12);
    c1->SetBottomMargin(0.12);
    TH1D* hlowmult[nbinZDC];
    TLegend* lp = new TLegend(0.57,0.7,0.87,0.89);
    lp->SetTextSize(0.021);   
    lp->SetBorderSize(0);
    for (int i =0; i<nbinZDC-1;i++){
        hlowmult[i] = (TH1D*)file->Get(Form("%s-%s/hRatioClone%i","EEsel_fixedlowmult","ZDCFixLowmult",i));
        hlowmult[i]->GetYaxis()->SetTitle("#frac{ #varepsilon_{part}^{Pythia62011} }{#varepsilon_{part}^{EPOS-LHC}}");
        hlowmult[i]->GetYaxis()->SetTitleOffset(2.1);
        hlowmult[i]->GetYaxis()->SetRangeUser(0.8,1.2);
        hlowmult[i]->GetXaxis()->SetRangeUser(0.6,6.5);
        hlowmult[i]->SetStats(0);
        hlowmult[i]->SetTitle("");
        hlowmult[i]->SetLineWidth(2);
        hlowmult[i]->Draw("HIST SAME");
        lp->AddEntry(hlowmult[i],Form("ZDC %.0f-%.0f ",percentileZDCMerged[i],percentileZDCMerged[i+1]),"L");
    } 
    lp->Draw("SAME"); 
    xlabel-> DrawText(0.35, 0.85, "fixed V0M [70,100]%");
    c1->SaveAs("ZDCSgnLossSystematicLowMult.png");

    TCanvas* c2 = new TCanvas("c2","",1300,1100);
    c2->SetLeftMargin(0.19);
    c2->SetRightMargin(0.12);
    c2->SetBottomMargin(0.12);
    TH1D* hhighmult[nbinZDC];
    TLegend* ld = new TLegend(0.57,0.7,0.87,0.89);
    ld->SetTextSize(0.021);   
    ld->SetBorderSize(0);
    for (int i =0; i<nbinZDC;i++){
        hhighmult[i] = (TH1D*)file->Get(Form("%s-%s/hRatioClone%i","EEsel_fixedhighmult","ZDCFixHighmult",i));
        hhighmult[i]->GetYaxis()->SetTitle("#frac{  #varepsilon_{part}^{Pythia62011} }{#varepsilon_{part}^{EPOS-LHC}}");
        hhighmult[i]->GetYaxis()->SetTitleOffset(2.1);
        hhighmult[i]->GetYaxis()->SetRangeUser(0.95,1.05);
        hhighmult[i]->GetXaxis()->SetRangeUser(0.6,6.5);
        hhighmult[i]->SetStats(0);
        hhighmult[i]->SetTitle("");
        hhighmult[i]->SetLineWidth(2);
        hhighmult[i]->Draw("HIST SAME");
        ld->AddEntry(hhighmult[i],Form("ZDC %.0f-%.0f ",percentileZDC[i],percentileZDC[i+1]),"L");
   
    } 
    ld->Draw("SAME"); 
    xlabel-> DrawText(0.35, 0.85, "fixed V0M [0,30]%");
    c2->SaveAs("ZDCSgnLossSystematicHighMult.png");


    TCanvas* c11 = new TCanvas("c11","",1300,1100);
    c11->SetLeftMargin(0.19);
    c11->SetRightMargin(0.12);
    c11->SetBottomMargin(0.12);
    TH1D* hlowee[nbinZDC];
    TLegend* ldd = new TLegend(0.57,0.7,0.87,0.89);
    ldd->SetTextSize(0.021);   
    ldd->SetBorderSize(0);
    for (int i =0; i<nbinV0-1;i++){
        hlowee[i] = (TH1D*)file->Get(Form("%s-%s/hRatioClone%i","multsel_fixedlowEE","V0FixLowEE",i));
        hlowee[i]->GetYaxis()->SetTitle("#frac{ #varepsilon_{part}^{Pythia62011} }{#varepsilon_{part}^{EPOS-LHC}}");
        hlowee[i]->GetYaxis()->SetTitleOffset(2.1);
        hlowee[i]->GetYaxis()->SetRangeUser(0.9,1.1);
        hlowee[i]->GetXaxis()->SetRangeUser(0.6,6.5);
        hlowee[i]->SetStats(0);
        hlowee[i]->SetTitle("");
        hlowee[i]->SetLineWidth(2);
        hlowee[i]->Draw("HIST SAME");
        ldd->AddEntry(hlowee[i],Form("V0M %.0f-%.0f ",percentileV0Merged[i],percentileV0Merged[i+1]),"L");
   
    } 
    ldd->Draw("SAME"); 
    xlabel-> DrawText(0.35, 0.85, "fixed ZDC [70,100]%");
    c11->SaveAs("V0MSgnLossSystematicLowEE.png");

    TCanvas* c22 = new TCanvas("c22","",1300,1100);
    c22->SetLeftMargin(0.19);
    c22->SetRightMargin(0.12);
    c22->SetBottomMargin(0.12);
    TH1D* hhighee[nbinZDC];
    TLegend* dd = new TLegend(0.57,0.7,0.87,0.89);
    dd->SetTextSize(0.021);   
    dd->SetBorderSize(0);
    for (int i =0; i<nbinV0-1;i++){
        hhighee[i] = (TH1D*)file->Get(Form("%s-%s/hRatioClone%i","multsel_fixedhighEE","V0FixHighEE",i));
        hhighee[i]->GetYaxis()->SetTitle("#frac{  #varepsilon_{part}^{Pythia62011} }{#varepsilon_{part}^{EPOS-LHC}}");
        hhighee[i]->GetYaxis()->SetTitleOffset(2.1);
        hhighee[i]->GetYaxis()->SetRangeUser(0.9,1.1);
        hhighee[i]->GetXaxis()->SetRangeUser(0.6,6.5);
        hhighee[i]->SetStats(0);
        hhighee[i]->SetTitle("");
        hhighee[i]->SetLineWidth(2);
        hhighee[i]->Draw("HIST SAME");
        dd->AddEntry(hlowee[i],Form("V0M %.0f-%.0f ",percentileV0Merged[i],percentileV0Merged[i+1]),"L");
   
    } 
    dd->Draw("SAME"); 
    xlabel-> DrawText(0.35, 0.85, "fixed ZDC [0,30]%");


    c22->SaveAs("V0MSgnLossSystematicHighEE.png");
    /*c2->SaveAs();
    c11->SaveAs();
    c22->SaveAs();
*/
}
