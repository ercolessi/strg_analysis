double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 );
Bool_t PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin);

void DrawCorrFactor(int it = 0, TString fWhichEstimator = "ZDC", 
Double_t lLoMult = 70., Double_t lHiMult = 100., 
Double_t lLoEE = 0., Double_t lHiEE = 100.){

 TString fWhichOtherEstimator = "ZDC";
   TString fWhichParticle = "XiMinus";
    TString fWhichAntiParticle = "XiPlus";
    
  
	// Def multiplicity
   // Def multiplicity / eff energy
    Double_t percentileV0[] = {0.,5,10,15,20,30,40,50,70,100};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    //
    Double_t percentileZDC[] = {0,20,40,50,60,70,80,90,100};
    const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;
    //
    Double_t * mult;
    // Choose scenario and initialize correct variables
    int tempbinnumber = 0;
    if (fWhichEstimator.Contains("V0M")) {
        tempbinnumber = nbinV0;
        mult = percentileV0;
    } 
    else if (fWhichEstimator.Contains("ZDC")) {
        tempbinnumber = nbinZDC;
        mult = percentileZDC;
        fWhichOtherEstimator = "V0M";
    } 
    else {cout << "No valid name for estimator... its V0M or ZDC" << endl; return;}
    const int multbinnumb = tempbinnumber;
 
  Double_t ptbinlimits[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
  Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;

 TFile* fileMC = new TFile("../MC15g3b1_ZDCRuns.root","READ");
    TList* clistMC = (TList*)fileMC->Get("PWGLF_StrVsMult_MC/cList");
    TH3D* h3DGenerated = (TH3D*)clistMC->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", fWhichParticle.Data()));
    // find projection bins in rapidity
    Int_t rapBin_min = h3DGenerated->GetYaxis()->FindBin( -0.5+1.e-6 );
    Int_t rapBin_max = h3DGenerated->GetYaxis()->FindBin( +0.5-1.e-6 );
    // 
    // find projection bins in multiplicity
    Int_t multBin_min = h3DGenerated->GetZaxis()->FindBin( 0.+1.e-6 );
    Int_t multBin_max = h3DGenerated->GetZaxis()->FindBin( 100.-1.e-6 );
    // get the projection in pt
    TH1D* HistGenPart = (TH1D*)h3DGenerated->ProjectionX(Form("HistGen%s", fWhichParticle.Data()),rapBin_min,rapBin_max, multBin_min, multBin_max);
    TH1D* HistGenAntiPart = (TH1D*)h3DGenerated->ProjectionX(Form("HistGen%s", fWhichAntiParticle.Data()),rapBin_min,rapBin_max, multBin_min, multBin_max);
    TH1D* HistGen = (TH1D*)HistGenPart->Clone("HistGen");
    HistGen->Reset();
    for (int bin = 1; bin <= HistGen->GetNbinsX(); bin++){   
      HistGen->SetBinContent(bin, HistGenPart->GetBinContent(bin) + HistGenAntiPart->GetBinContent(bin));
      HistGen->SetBinError(bin, TMath::Sqrt(HistGenPart->GetBinError(bin)*HistGenPart->GetBinError(bin)+
           HistGenAntiPart ->GetBinError(bin)*HistGenAntiPart->GetBinError(bin)));
    }    

 TFile* fileP = new TFile(Form("Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichParticle.Data(), 0.,100.,0.,100.));
    TFile* fileAntiP = new TFile(Form("Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichAntiParticle.Data(), 0.,100.,0.,100.));
    TH1D* HistRecoPart = (TH1D*)fileP->Get("fHistReco");
    TH1D* HistRecoAntiPart = (TH1D*)fileAntiP->Get("fHistReco");
    TH1F* HistReco = (TH1F*) HistRecoPart->Clone("HistReco");
    HistReco->Reset();
    for (int bin = 1 ; bin <=  HistRecoPart->GetNbinsX(); bin ++ ){
       HistReco->SetBinContent(bin, HistRecoPart->GetBinContent(bin) + HistRecoAntiPart->GetBinContent(bin));
       HistReco->SetBinError(bin, TMath::Sqrt(HistRecoPart->GetBinError(bin)*HistRecoPart->GetBinError(bin)+HistRecoAntiPart->GetBinError(bin)*HistRecoAntiPart->GetBinError(bin)));
    }
    
    //Rebin Gen e Reco per fare efficienze default
    TH1F* HistReco_rebin   = new TH1F("HistReco_rebin","Cascade MC count;p_{T} (GeV/c);Counts", ptbinnumb, ptbinlimits);
    TH1F* HistGen_rebin  = new TH1F("HistGen_rebin","Cascade MC count;p_{T} (GeV/c);Counts", ptbinnumb, ptbinlimits);  
    TH1F* EfficiencyMB  = new TH1F("EfficiencyMB","Cascade MC count;p_{T} (GeV/c);Counts", ptbinnumb, ptbinlimits);  

    Double_t tempptreco;
    for(long i = 1; i<HistReco->GetNbinsX()+1;i++){
      tempptreco = HistReco->GetXaxis()->GetBinCenter(i);
      for(long filling = 0; filling<HistReco->GetBinContent(i); filling++){
        HistReco_rebin->Fill(tempptreco);
      }
    }
    Double_t tempptgen;
    for(long i = 1; i<HistGen->GetNbinsX()+1;i++){
      tempptgen = HistGen->GetXaxis()->GetBinCenter(i);
      for(long filling = 0; filling<HistGen->GetBinContent(i); filling++){
        HistGen_rebin->Fill(tempptgen);
      }
    }   
  
    for (int bin = 1; bin <= EfficiencyMB->GetNbinsX(); bin ++){
      if (HistGen_rebin->GetBinContent(bin) != 0){
        EfficiencyMB->SetBinContent(bin , HistReco_rebin->GetBinContent(bin)/HistGen_rebin->GetBinContent(bin));
        EfficiencyMB->SetBinError(bin, ErrorInRatio(HistReco_rebin->GetBinContent(bin),HistReco_rebin->GetBinError(bin),
            HistGen_rebin->GetBinContent(bin),HistGen_rebin->GetBinError(bin)));
      }
    }




//++++++++++++++++++++++++++++++++++++++++
// Fit spectra w/ Levy-Tsallis 
//++++++++++++++++++++++++++++++++++++++++ 

  //Get from files 
 
  TH1F* Eff[multbinnumb];
  TH1F* EffCorr[multbinnumb];

  TH1F* RatioEff[multbinnumb];
  TString text = "";

  TFile* effcorr = new TFile(Form("%sMCptshape%sIT%i_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi",fWhichEstimator.Data(),it,lLoMult, lHiMult,lLoEE,lHiEE));
  TFile* eff;
  if (it>0) eff = new TFile(Form("%sMCptshape%sIT%i_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi",fWhichEstimator.Data(),it-1,lLoMult, lHiMult,lLoEE,lHiEE));
  
 
  for(int nmult = 0; nmult < multbinnumb; nmult++)
  {
   
    if (it>0){
      EffCorr[nmult] = (TH1F *) effcorr->Get(Form("Efficiencies/CorrEfficiency%i",nmult));
      Eff[nmult] = (TH1F *)  eff->Get(Form("Efficiencies/CorrEfficiency%i",nmult));
    }

    if (it==0){
      EffCorr[nmult] = (TH1F *) effcorr->Get(Form("Efficiencies/CorrEfficiency%i",nmult));
      Eff[nmult] = EfficiencyMB;
    }
     RatioEff[nmult] = (TH1F*)EffCorr[nmult]->Clone(Form("RatioEff%i",nmult));
      //RatioEff[nmult]->Reset();
      DivideAndComputeRogerBarlow(RatioEff[nmult],Eff[nmult]);

    /*for (int bin = 1; bin <= RatioEff[nmult]->GetNbinsX(); bin ++){
    	RatioEff[nmult]->SetBinContent(bin,EffCorr[nmult]->GetBinContent(bin)/Eff[nmult]->GetBinContent(bin));
    	RatioEff[nmult]->SetBinError(bin,1E-12);
      //Eff[nmult]->GetBinError(bin)/EffCorr[nmult]->GetBinContent(bin));
      //ErrorInRatio(EffCorr[nmult]->GetBinContent(bin),EffCorr[nmult]->GetBinError(bin), Eff[nmult]->GetBinContent(bin), Eff[nmult]->GetBinError(bin)));
	  }*/
  }

  //Set colors
  RatioEff[0]->SetLineColor(kRed + 1);
  RatioEff[1]->SetLineColor(kRed - 4);
  RatioEff[2]->SetLineColor(kOrange + 7);
  RatioEff[3]->SetLineColor(kOrange - 3);
  RatioEff[4]->SetLineColor(kYellow + 1);
  RatioEff[5]->SetLineColor(kSpring - 7);
  RatioEff[6]->SetLineColor(kGreen + 2);
  RatioEff[7]->SetLineColor(kAzure + 8);
  RatioEff[8]->SetLineColor(kBlue - 4);
 // if (fWhichEstimator.Contains("V0M"))RatioEff[9]->SetLineColor(kBlue + 3);

  RatioEff[0]->SetMarkerColor(kRed + 1);
  RatioEff[1]->SetMarkerColor(kRed - 4);
  RatioEff[2]->SetMarkerColor(kOrange + 7);
  RatioEff[3]->SetMarkerColor(kOrange - 3);
  RatioEff[4]->SetMarkerColor(kYellow + 1);
  RatioEff[5]->SetMarkerColor(kSpring - 7);
  RatioEff[6]->SetMarkerColor(kGreen + 2);
  RatioEff[7]->SetMarkerColor(kAzure + 8);
  RatioEff[8]->SetMarkerColor(kBlue - 4);
  //if (fWhichEstimator.Contains("V0M"))RatioEff[9]->SetMarkerColor(kBlue + 3);

  TCanvas* t = new TCanvas("t","",800,800);
  t->SetGridy();
  t->SetRightMargin(0.08);
  t->SetLeftMargin(0.15);
  t->SetBottomMargin(0.15);
  TLegend* l = new TLegend(0.6,0.75,0.9,0.89);
  l->SetTextSize(0.022);
  l->SetBorderSize(0);  
  if (fWhichEstimator.Contains("V0M")) {
    l->AddEntry(RatioEff[0],"V0M 0-1% corr. spectrum","LEP");
    l->AddEntry(RatioEff[8],"V0M 70-100% corr. spectrum","LEP");
  } else if (fWhichEstimator.Contains("ZDC")) {
    l->AddEntry(RatioEff[0],"ZDC 0-20% corr. spectrum","LEP");
    l->AddEntry(RatioEff[7],"ZDC 90-100% corr. spectrum","LEP");
  }
 
  double a = 0.9, b = 1.1;
  if (it>1) {a=0.98;b=1.02;}
  RatioEff[0]->GetXaxis()->SetRangeUser(0.6,5);
  RatioEff[0]->GetYaxis()->SetRangeUser(a,b);
  RatioEff[0]->GetYaxis()->SetTitle(Form("Corr. factor  (#varepsilon_{it%i} / #varepsilon_{it%i})",it,it-1));
  RatioEff[0]->GetYaxis()->SetTitleOffset(1.6);
  RatioEff[0]->GetXaxis()->SetTitleOffset(1.4);
  RatioEff[0]->SetStats(0);
  RatioEff[0]->SetTitle("");
  RatioEff[0]->SetMarkerStyle(8);  
  RatioEff[0]->SetMarkerSize(1.4);  
  RatioEff[0]->Draw();
  for(int nmult = 1; nmult < multbinnumb; nmult++){
  	RatioEff[nmult]->SetMarkerStyle(8);  
    RatioEff[nmult]->SetMarkerSize(1.4); 
  	RatioEff[nmult]->Draw("SAME");
  }
  l->Draw("SAME");
  TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
  xlabel-> SetNDC();
  xlabel-> SetTextColor(1);
  xlabel-> SetTextSize(0.06);
  xlabel-> SetTextAlign(22);
  xlabel-> SetTextAngle(0);
  xlabel-> DrawLatex(0.6, 0.82, text.Data());
  TLatex *xlabel2 = new TLatex();
	xlabel2->SetTextFont(42);
  xlabel2-> SetNDC();
  xlabel2-> SetTextColor(1);
  xlabel2-> SetTextSize(0.03);
  xlabel2-> SetTextAlign(22);
  xlabel2-> SetTextAngle(0);
  if (fWhichEstimator.Contains("V0M")) {xlabel2-> DrawLatex(0.78, 0.71, Form("ZDC fixed [%.0f-%.0f]",lLoEE,lHiEE));}
  else {xlabel2-> DrawLatex(0.78, 0.71, Form("V0M fixed [%.0f-%.0f]",lLoMult,lHiMult));}
  TLatex *xlabe = new TLatex();
	xlabe->SetTextFont(42);
  xlabe-> SetNDC();
  xlabe-> SetTextColor(1);
  xlabe-> SetTextSize(0.04);
  xlabe-> SetTextAlign(22);
  xlabe-> SetTextAngle(0);
  xlabe-> DrawLatex(0.32, 0.79, Form("Iteration %i",it));
  t->SaveAs(Form("images/%s/CorrFactor%s_it%i-%i_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.png",fWhichEstimator.Data(),"Xi",it,it-1,lLoMult,lHiMult,lLoEE,lHiEE));

}

//---------------------------------------------------------------------------------------------------
double ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr) {
  // Error in a Ratio
  if (B != 0) {
    Double_t errorfromtop = Aerr * Aerr / (B * B);
    Double_t errorfrombottom = ((A * A) / (B * B * B * B)) * Berr * Berr;
    return TMath::Sqrt(TMath::Abs(errorfromtop - errorfrombottom));
  }
  return 1.;
}

//----------------------------------------------------------------------------------------------------
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 ){ 
  //Use Roger Barlow "sigma_{delta}" as errors for ratios
  Double_t lh1NBins = h1->GetNbinsX(); 
  Double_t lh2NBins = h2->GetNbinsX(); 

  if( lh1NBins != lh2NBins ){ 
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return;
  }

  Double_t lSigmaDelta[100]; 
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){ 
    //Computation of roger barlow sigma_{delta} 
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(h1->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
    //Computation of relationship to h2 for plotting in ratio plot 
    if ( h2->GetBinContent(i) > 1e-12 ){ 
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
}

//----------------------------------------------------------------------------------------------------
Bool_t PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin){
  double dev = h->GetBinContent(bin);
  double RBsigma = h->GetBinError(bin);

  if( TMath::Abs(dev/RBsigma) > nsigmas ) {return kTRUE;}
  else {return kFALSE;}
}
