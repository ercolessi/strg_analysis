double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 );

void DrawCorrFactor(TString fWhichParticle = "XiPlus"){

	// Def multiplicity
  Float_t mult[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
  Long_t multbinnumb = sizeof(mult) / sizeof(Float_t) - 1;
  Float_t multOmega[] = {0, 5, 15, 30, 50, 100};
  Long_t multbinnumbOmega = sizeof(mult) / sizeof(Float_t) - 1;
  Double_t ptbinlimits[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
  Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;

//++++++++++++++++++++++++++++++++++++++++
// Fit spectra w/ Levy-Tsallis 
//++++++++++++++++++++++++++++++++++++++++ 

  //Get from files 
  TFile* filemultcorr[multbinnumb];
  TFile* filemult[multbinnumb];
  
  TH1F* Eff[multbinnumb];
  TH1F* EffCorr[multbinnumb];

  TH1F* RatioEff[multbinnumb];
  TString effname = "fHistEfficiency";
 
  for(int nmult = 0; nmult < multbinnumb; nmult++)
  {
    filemultcorr[nmult] = new TFile(Form("IT3/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichParticle.Data(), mult[nmult], mult[nmult+1], 0., 100.));
    filemult[nmult] = new TFile(Form("IT2/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichParticle.Data(), mult[nmult], mult[nmult+1], 0., 100.));
    if (fWhichParticle.Data() == "XiMinus" || fWhichParticle.Data() == "OmegaMinus"){
    	effname = "fHistPureEfficiency";
    }
    if (fWhichParticle.Data() == "XiPlus" || fWhichParticle.Data() == "OmegaPlus"){
    	effname = "fHistEfficiency";
    }
    EffCorr[nmult] = (TH1F *) filemultcorr[nmult]->Get(effname.Data());
    Eff[nmult] = (TH1F *) filemult[nmult]->Get(effname.Data());

    RatioEff[nmult] = (TH1F*)EffCorr[nmult]->Clone(Form("RatioEff%i",nmult));
    //RatioEff[nmult]->Reset();
    DivideAndComputeRogerBarlow(RatioEff[nmult],Eff[nmult]);
    /*for (int bin = 1; bin <= RatioEff[nmult]->GetNbinsX(); bin ++){
    	RatioEff[nmult]->SetBinContent(bin,EffCorr[nmult]->GetBinContent(bin)/Eff[nmult]->GetBinContent(bin));
    	RatioEff[nmult]->SetBinError(bin,ErrorInRatio(EffCorr[nmult]->GetBinContent(bin),EffCorr[nmult]->GetBinError(bin), Eff[nmult]->GetBinContent(bin), Eff[nmult]->GetBinError(bin)));
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
  RatioEff[9]->SetLineColor(kBlue + 3);

  RatioEff[0]->SetMarkerColor(kRed + 1);
  RatioEff[1]->SetMarkerColor(kRed - 4);
  RatioEff[2]->SetMarkerColor(kOrange + 7);
  RatioEff[3]->SetMarkerColor(kOrange - 3);
  RatioEff[4]->SetMarkerColor(kYellow + 1);
  RatioEff[5]->SetMarkerColor(kSpring - 7);
  RatioEff[6]->SetMarkerColor(kGreen + 2);
  RatioEff[7]->SetMarkerColor(kAzure + 8);
  RatioEff[8]->SetMarkerColor(kBlue - 4);
  RatioEff[9]->SetMarkerColor(kBlue + 3);

  new TCanvas;
  RatioEff[0]->GetYaxis()->SetRangeUser(0.9,1.1);
  RatioEff[0]->GetYaxis()->SetTitle("Corr. factor  (#varepsilon_{it3} / #varepsilon_{it2})");
  RatioEff[0]->GetYaxis()->SetTitleOffset(1.2);
  RatioEff[0]->SetStats(0);
  for(int nmult = 0; nmult < multbinnumb; nmult++){
  	RatioEff[nmult]->SetMarkerStyle(8);    
  	RatioEff[nmult]->Draw("SAME");
  }

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
