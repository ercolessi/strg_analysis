Double_t doerror ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
TH1D* doratio(TH1D* h1, TH1D* h2);
void canvas_hestetics(TCanvas* c);

void compareESDAOD(Double_t lMultBoundLo = 0.0,
  Double_t lMultBoundHi = 100.0,
  Double_t lEEBoundLo = 0.0,
  Double_t lEEBoundHi = 100.0){

    
    //--------- spectra files ------------
    TFile* resaod = TFile::Open("~/Scaricati/resultsSpectra/Results-Lambda-13TeV-V0M_0.0_100.0_UseMCRatioFD_WithV0refitAndImprovedDCA.root");
        //Form("Results-Lambda-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_AOD.root",lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi), "READ");
    TFile* resesd = TFile::Open(Form("results/Results-Lambda-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_ESD.root",lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi), "READ");
    
   
    TH1D* aodsp = (TH1D*)resaod->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
    TH1D* esdsp = (TH1D*)resesd->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");

    TH1D* aodeff = (TH1D*)resaod->Get("fHistEfficiency");
    TH1D* esdeff = (TH1D*)resesd->Get("fHistEfficiency");

    TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
    xlabel-> SetNDC();
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.1);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);

    TLegend* l = new TLegend(0.2,0.7,0.4,0.85);
    l->SetTextSize(0.022);
    l->SetBorderSize(0); 
      
     //#Lambda

    TH1D* hratio = doratio(aodsp,esdsp);
    hratio->SetMarkerStyle(20);
    hratio->SetMarkerColor(kBlack);
    hratio->SetLineColor(kBlack);
    
    TCanvas* c = new TCanvas("c","",1000,800);
    canvas_hestetics(c);
    hratio->GetYaxis()->SetRangeUser(0.8,1.2);
    hratio->GetYaxis()->SetTitle("Raw Spectra Ratio: AOD/ESD");
    hratio->GetYaxis()->SetTitleSize(0.04);
    hratio->GetXaxis()->SetTitleSize(0.04);
    hratio->SetStats(0);
    hratio->SetTitle(Form("SPD [%.0f-%.0f], V0M [%.0f-%.0f]",lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi));
    hratio->SetMarkerSize(1.7);
    hratio->Draw();
  
    l->AddEntry(hratio,"#it{p}_{T} differential","EP");
    l->SetTextSize(0.04);
    l->Draw("SAME");
    
    TLine* l5p = new TLine(0.6,1.05,10.,1.05);
    l5p->SetLineColor(kRed);
    l5p->SetLineStyle(9);
    
    TLine* l5m = new TLine(0.6,.95,10.,.95);
    l5m->SetLineColor(kRed);
    l5m->SetLineStyle(9);

    l5p->Draw("SAME");
    l5m->Draw("SAME");
    xlabel-> DrawLatex(0.75, 0.78, Form("#Lambda"));


    TH1D* hratioeff = doratio(aodeff,esdeff);
    hratioeff->SetMarkerStyle(20);
    hratioeff->SetMarkerColor(kBlack);
    hratioeff->SetLineColor(kBlack);
    
    TCanvas* ce = new TCanvas("ce","",1000,800);
    canvas_hestetics(ce);
    hratioeff->GetYaxis()->SetRangeUser(0.8,1.2);
    hratioeff->GetYaxis()->SetTitle("Efficiency Ratio: AOD/ESD");
    hratioeff->GetYaxis()->SetTitleSize(0.04);
    hratioeff->SetStats(0);
    hratioeff->SetTitle(Form("SPD [%.0f-%.0f], V0M [%.0f-%.0f]",lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi));
    hratioeff->GetXaxis()->SetTitleSize(0.04);
    hratioeff->Draw();
    hratioeff->SetMarkerSize(1.7);
    xlabel-> DrawLatex(0.75, 0.78, Form("#Lambda"));

    TLine* l5pe = new TLine(0.6,1.05,10.,1.05);
    l5pe->SetLineColor(kRed);
    l5pe->SetLineStyle(9);
    
    TLine* l5me = new TLine(0.6,.95,10.,.95);
    l5me->SetLineColor(kRed);
    l5me->SetLineStyle(9);

    l5pe->Draw("SAME");
    l5me->Draw("SAME");


}

//=================================
Double_t doerror ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ){
    //Error in a Ratio
    if(B!=0){
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop - errorfrombottom) );
    }
    return 1;
}

//=================================
TH1D* doratio(TH1D* h1, TH1D* h2){

    Double_t lh1NBins = h1->GetNbinsX(); 
    Double_t lh2NBins = h2->GetNbinsX(); 

    if( lh1NBins != lh2NBins ){ 
        cout<<"Problem! Number of bins doesn't match! "<<endl;
    }

    Double_t lSigmaDelta[100]; 
    for( Int_t i=1; i<h1->GetNbinsX()+1; i++){ 
        //Computation of roger barlow sigma_{delta} 
        lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(h1->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
        //Computation of relationship to h2 for plotting in ratio plot 
        if ( h2->GetBinContent(i) != 0 ){ 
        lSigmaDelta[i] /= h2->GetBinContent(i); 
        }else{ 
        lSigmaDelta[i] = 0; 
        }
    }
    //Regular Division 
    h1 -> Divide( h2 ); 
    //Replace Erorrs 
    for( Int_t i=1; i<h1->GetNbinsX()+1; i++){ 
        h1->SetBinError(i, lSigmaDelta[i]);
    }  

    return h1;
}

void canvas_hestetics(TCanvas* c) {
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.12);
    c->SetTopMargin(0.1);
    c->SetTicky();
    c->SetTickx();
}