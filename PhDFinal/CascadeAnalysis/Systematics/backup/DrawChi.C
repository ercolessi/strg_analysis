void DrawChi( 
  TString fWhichParticle = "Xi", 
  TString sel = "V0M", 
  Double_t lFixedLo = 40., 
  Double_t lFixedHi = 50.){

    TString var = "SPDClusters";
    if (sel.Contains("SPD")) var ="V0M";

    TString inputfilename = Form("ExtrSyst-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", 
                                fWhichParticle.Data(), var.Data(), sel.Data(), lFixedLo, lFixedHi);
    TFile* inputfile = TFile::Open(inputfilename);
    //
    TGraph* ChiF[6]; //0 is Levy
    TString func[6] = {"Levy","BlastWave","Boltz","MTexpo","FermiDirac","BoseEinstein"};
    for (int nf = 0; nf < 6 ; nf++){
        ChiF[nf] = (TGraph*)inputfile->Get(Form("%s/Chi_%s",func[nf].Data(),func[nf].Data()));   //primo bin + 0.100  
        ChiF[nf]->SetLineWidth(2);
    }

    TCanvas* c = new TCanvas("c","c",1500,1000);
    ChiF[0]->GetYaxis()->SetRangeUser(0.,5);
    ChiF[0]->GetXaxis()->SetRangeUser(0.,100);
    ChiF[0]->GetXaxis()->SetTitle(Form("%s percentile",var.Data()));
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
    TLatex *xlabel2 = new TLatex();
	//xlabel2->SetTextFont(42);
    xlabel2-> SetNDC();
    xlabel2-> SetTextColor(1);
    xlabel2-> SetTextSize(0.03);
    xlabel2-> SetTextAlign(22);
    xlabel2-> SetTextAngle(0);
    
    TLine* el = new TLine(0.,2.5,100,2.5);
    el->SetLineWidth(2);
    el->SetLineStyle(9);
    el->SetLineColor(kRed);
    //el->Draw("SAME");
    xlabel2-> DrawLatex(0.4, 0.81, Form("%s fixed [%.0f-%.0f] ",sel.Data(),lFixedLo, lFixedHi));
    l->Draw("SAME");
   
    c->SaveAs(Form("images/Extrap/ChiExtrapolation-%s-%s_%03.0f_%03.0f.png",fWhichParticle.Data(),sel.Data(),lFixedLo, lFixedHi));
      
    //Percentile
    Double_t * percentile;
    Long_t percbinnumb;
    Double_t pmb[] = {0,100};
    Long_t nmb = sizeof(pmb)/sizeof(Double_t) - 1;
    Double_t p0[] = {0,5,10,15,20,30,40,50,70,100};
    Long_t n0 = sizeof(p0)/sizeof(Double_t) - 1;
    Double_t p1[] = {0,5,10,20,30,40,50,100};
    Long_t n1 = sizeof(p1)/sizeof(Double_t) - 1;
    Double_t p2[] = {0,20,30,40,50,60,70,100};
    Long_t n2 = sizeof(p2)/sizeof(Double_t) - 1; 
    Double_t p4[] = {0,5,10,20,30,40,50,100};
    Long_t n4 = sizeof(p4)/sizeof(Double_t) - 1;
    Double_t p5[] = {0,10,20,30,40,50,60,70,100};
    Long_t n5 = sizeof(p5)/sizeof(Double_t) - 1;
    Double_t pOmega[] = {0,5,10,30,50,100};
    Long_t npOmega = sizeof(pOmega)/sizeof(Double_t) - 1;
    Double_t pOmega2[] = {0,40,70,100};
    Long_t npOmega2 = sizeof(pOmega2)/sizeof(Double_t) - 1;

      if (sel.Contains("SPD") && lFixedLo == 0. && lFixedHi == 100.){
        if (fWhichParticle.Contains("Omega")){
          percentile = pOmega;
          percbinnumb = npOmega;
        } else{
          percentile = p0;
          percbinnumb = n0;
        }
      }
      if (sel.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
        if (fWhichParticle.Contains("Omega")){
          percentile = pOmega;
          percbinnumb = npOmega;
        } else{
          percentile = p1;
          percbinnumb = n1;
        }
      }
      if (sel.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        if (fWhichParticle.Contains("Omega")){
          percentile = pOmega2;
          percbinnumb = npOmega2;
        } else{
        percentile = p2;
        percbinnumb = n2;
        }
      }
      if (sel.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        percentile = p4;
        percbinnumb = n4;
      }
      if (sel.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
        percentile = p5;
        percbinnumb = n5;
      }
    


    TFile* resfile = TFile::Open(Form("../correctedresults/CorrSpectra-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", 
                                fWhichParticle.Data(), var.Data(), sel.Data(), lFixedLo, lFixedHi));
                                
    TH1D* h[percbinnumb];
    TF1* f[percbinnumb][6];

    TString folder = "NormCorr/AfterSignalLossCorr";
   
    for (int n = 0; n < percbinnumb ; n++){
       
        double multmin, multmax, eemin, eemax;
        //
        if (var.Contains("SPD")){
            multmin = percentile[n];
            multmax = percentile[n+1];
            eemin = lFixedLo; 
            eemax = lFixedHi;
        }
        //
        if (var.Contains("V0M") || var.Contains("ZDC")){
            eemin = percentile[n];
            eemax = percentile[n+1];
            multmin = lFixedLo; 
            multmax = lFixedHi;
        }
        h[n] = (TH1D*)resfile->Get(Form("%s/PtSpectrumCorr_%s_%.0f_%.0f_%s_%.0f_%.0f",folder.Data(),"SPDClusters",multmin,multmax,"V0M",eemin,eemax));    
        for (int i = 0; i < 6; i++){ 
          cout << Form("%s/%s_%.0f-%.0f",func[i].Data(),func[i].Data(),percentile[n],percentile[n+1]) << endl;
            f[n][i] = (TF1*)inputfile->Get(Form("%s/%s_%.0f-%.0f",func[i].Data(),func[i].Data(),percentile[n],percentile[n+1]));

        }
    }

    double n[percbinnumb];
    for (int i = 0; i < percbinnumb; i++){     
        h[i]->Scale(std::pow(2,percbinnumb-1-i));    
        n[i] = std::pow(2,percbinnumb-1-i);
        h[i]->SetLineWidth(2);  
    }
    TCanvas* c1 = new TCanvas("c1","c1",1000,1000);
    c1->SetLogy();
    c1->SetRightMargin(0.09);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.10);
    c1->SetTickx();
    c1->SetTicky();

    h[0]->GetYaxis()->SetRangeUser(1E-8,1E3);
    h[0]->GetYaxis()->SetTitle("1/N_{ev} d^{2}N/(dp_{T}dy) [(GeV/c)^{-1}]");
    h[0]->GetXaxis()->SetRangeUser(0.,6.5);
    h[0]->SetTitle("");
    h[0]->Draw();
    TLegend* lf = new TLegend(0.6,0.7,0.89,0.85);
    lf->SetTextSize(0.028);
    lf->SetBorderSize(0);
    
    for (int i = 0; i < percbinnumb; i++){ 
        //h[i]->SetMarkerStyle(kFullCircle);
        h[i]->SetMarkerSize(1.5);
        h[i]->SetLineWidth(3);
        h[i]->Draw("SAME");
        h[i]->SetStats(0);
        f[i][0]->SetLineStyle(7);
        f[i][0]->SetLineWidth(1);
        for (int fi = 0; fi < 6; fi++){ 
            f[i][fi]->SetLineWidth(2);
            f[i][fi]->Draw("SAME");
            if (i == 0) lf->AddEntry(f[0][fi],func[fi].Data(),"L");
        }
        
        h[i]->Draw("SAME");
        
    }
    xlabel2-> DrawLatex(0.4, 0.81, Form("%s fixed [%.0f-%.0f] ",sel.Data(),lFixedLo, lFixedHi));
    lf->Draw("SAME");
    TLatex *xlabel3 = new TLatex();
	//xlabel3->SetTextFont(42);
    xlabel3-> SetNDC();
    xlabel3-> SetTextColor(1);
    xlabel3-> SetTextSize(0.05);
    xlabel3-> SetTextAlign(22);
    xlabel3-> SetTextAngle(0);
    xlabel3-> DrawLatex(0.4, 0.76, Form("#%s",fWhichParticle.Data()));
    
    c1->SaveAs(Form("images/Extrap/FitExtrapolation-%s-%s_%03.0f_%03.0f.png",fWhichParticle.Data(),sel.Data(),lFixedLo, lFixedHi));
    


    //c1->SaveAs(Form("FitExtrapolation-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.png",var.Data(),lLoMult,lHiMult,lLoEE,lHiEE));
    
}
