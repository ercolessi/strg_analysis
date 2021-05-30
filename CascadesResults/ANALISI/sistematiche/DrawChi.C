void DrawChi( 
  TString fWhichEstimator = "ZDC", 
  Double_t lLoMult = 0., 
  Double_t lHiMult = 30., 
  Double_t lLoEE = 0., 
  Double_t lHiEE = 100.){

    TString inputfilename = Form("ExtSyst-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE);
    TFile* inputfile = TFile::Open(inputfilename);
    //
    TGraph* ChiF[6]; //0 is Levy
    TString func[6] = {"Levy","BlastWave","Boltz","MTexpo","FermiDirac","BoseEinstein"};
    for (int nf = 0; nf < 6 ; nf++){
        ChiF[nf] = (TGraph*)inputfile->Get(Form("%s/Chi_%s",func[nf].Data(),func[nf].Data()));   //primo bin + 0.100  
        ChiF[nf]->SetLineWidth(2);
    }

    TCanvas* c = new TCanvas("c","c",1500,1000);
    ChiF[0]->GetYaxis()->SetRangeUser(0.,2.99);
    ChiF[0]->GetXaxis()->SetRangeUser(0.,100);
    ChiF[0]->GetXaxis()->SetTitle(Form("%s percentile[%]",fWhichEstimator.Data()));
    ChiF[0]->GetYaxis()->SetTitle("#chi^{2} / NDF ");
    ChiF[0]->SetTitle("");
    ChiF[0]->Draw("AEP");
    TLegend* l = new TLegend(0.6,0.65,0.89,0.86);
    l->SetTextSize(0.034);
    l->SetBorderSize(0);
    for (int nf = 0; nf < 6 ; nf++){
        ChiF[nf]->Draw("SAME EP");
        l->AddEntry(ChiF[nf],func[nf].Data(),"LP");
    }
    l->Draw("SAME");
    TLatex *xlabel2 = new TLatex();
	xlabel2->SetTextFont(42);
    xlabel2-> SetNDC();
    xlabel2-> SetTextColor(1);
    xlabel2-> SetTextSize(0.03);
    xlabel2-> SetTextAlign(22);
    xlabel2-> SetTextAngle(0);
    if (fWhichEstimator.Contains("V0M")) {xlabel2-> DrawLatex(0.4, 0.81, Form("ZDC fixed [%.0f-%.0f] ",lLoEE,lHiEE));}
    else {xlabel2-> DrawLatex(0.4, 0.81, Form("V0M fixed [%.0f-%.0f]",lLoMult,lHiMult));}


    c->SaveAs(Form("ChiExtrapolation-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.png",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE));
   
    
   
    //
    Double_t percentileV0[] = {0.,5,10,15,20,30,40,50,70,100};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    //
    Double_t percentileZDC[] = {0,20,30,40,50,60,70,80,90,100};
    const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;
    //
    Double_t * mult, *systpercentile;
    Int_t *systcounter;
    //
    // Choose scenario and initialize correct variables
    int tempbinnumber = 0;
    if (fWhichEstimator.Contains("V0M")) {
        tempbinnumber = nbinV0;
        mult = percentileV0;
    } 
    else if (fWhichEstimator.Contains("ZDC")) {
        tempbinnumber = nbinZDC;
        mult = percentileZDC;
    } 
    else {cout << "No valid name for estimator... its V0M or ZDC" << endl; return;}
    const int multbinnumb = tempbinnumber;
    const int binn = multbinnumb;
   
    TH1D* h[multbinnumb];
    TF1* f[multbinnumb][6];
   
    for (int n = 0; n < multbinnumb ; n++){
        double zdcmin = lLoEE, zdcmax = lHiEE, v0min = lLoMult, v0max = lHiMult;
        if (fWhichEstimator.Contains("V0M")) {
            v0min = mult[n];
            v0max = mult[n+1];
        }
        if (fWhichEstimator.Contains("ZDC")) {
            zdcmin = mult[n];
            zdcmax = mult[n+1];
        }
        h[n] = (TH1D*)inputfile->Get(Form("XiSpectra_Stats_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",v0min,v0max,zdcmin,zdcmax));   //primo bin + 0.100  
        for (int i = 0; i < 6; i++){ 
            f[n][i] = (TF1*)inputfile->Get(Form("%s/%s_%.0f-%.0f",func[i].Data(),func[i].Data(),mult[n],mult[n+1]));

        }
    }

    double n[multbinnumb];
     for (int i = 0; i < multbinnumb; i++){     
         if (fWhichEstimator.Contains("ZDC") && lLoMult == 70. && lHiMult == 100.){
          h[i]->Scale(TMath::Power(5,binn-1-i));
          h[i]->SetLineWidth(2); 
          n[i] = TMath::Power(5,binn-1-i);
        }    
        else{  
        h[i]->Scale(TMath::Power(2,binn-1-i));    
        n[i] = TMath::Power(2,binn-1-i);
        h[i]->SetLineWidth(2);  
        }
    }
    TCanvas* c1 = new TCanvas("c1","c1",1600,1500);
    c1->SetLogy();
    c1->SetRightMargin(0.09);
    c1->SetLeftMargin(0.25);
    c1->SetBottomMargin(0.15);

    h[0]->GetYaxis()->SetRangeUser(1E-8,1E4);
    h[0]->GetYaxis()->SetTitle("1/N_{ev} d^{2}N/(dp_{T}dy) [(GeV/c)^{-1}]");
    h[0]->GetXaxis()->SetRangeUser(0.,6.5);
    h[0]->SetTitle("");
    h[0]->Draw();
    for (int i = 0; i < multbinnumb; i++){ 
        h[i]->Draw("SAME");
        h[i]->SetStats(0);
        f[i][0]->SetLineStyle(7);
        f[i][0]->SetLineWidth(1);
        for (int fi = 0; fi < 6; fi++){ 
            f[i][fi]->SetLineWidth(2);
            f[i][fi]->Draw("SAME");
        }
    }
    
    if (fWhichEstimator.Contains("V0M")) {xlabel2-> DrawLatex(0.4, 0.81, Form("ZDC fixed [%.0f-%.0f] ",lLoEE,lHiEE));}
    else {xlabel2-> DrawLatex(0.4, 0.81, Form("V0M fixed [%.0f-%.0f]",lLoMult,lHiMult));}

    TLegend* lf = new TLegend(0.6,0.7,0.89,0.85);
    lf->SetTextSize(0.028);
    lf->SetBorderSize(0);
    ChiF[0]->SetLineStyle(7);
    for (int nf = 0; nf < 6 ; nf++){
        lf->AddEntry(ChiF[nf],func[nf].Data(),"L");
    }
    lf->Draw("SAME");

    c1->SaveAs(Form("FitExtrapolation-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.png",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE));
    






}
