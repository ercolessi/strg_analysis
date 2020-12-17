#include <fstream>
#include <iostream>

double ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr);

void DivideAndComputeRogerBarlow(TH1F *h1, TH1F *h2);

TF1 *LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5.,
                 Double_t C = 0.1, Double_t norm = 1.);

//--------------------------------------------------------------
//------------------- MAIN FUNCTION ----------------------------
//--------------------------------------------------------------
void DoCorrectionMCpTshape(TString fWhichParticle = "XiPlus") {

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
  TFile* filemult[multbinnumb];
  TFile* fileCascMB = new TFile(Form("IT0/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichParticle.Data(), 0.,100., 0., 100.));

  TF1* Levy[multbinnumb];
  TF1* LevyMCMB = LevyTsallis(Form("LevyMC%s_%.0f-%.0f",fWhichParticle.Data(),0.,100.), 1.32131 );
  TF1* LevyRDMB = LevyTsallis(Form("LevyRD%s_%.0f-%.0f",fWhichParticle.Data(),0.,100.), 1.32131 );

  TH1F* Spectra[multbinnumb];
  TH1F* SpectraRaw[multbinnumb];
  TH1F* SpectraMB = (TH1F *) fileCascMB->Get(Form("fHistPt%s", fWhichParticle.Data()));

  for(int nmult = 0; nmult < multbinnumb; nmult++)
  {
    filemult[nmult] = new TFile(Form("IT0/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichParticle.Data(), mult[nmult], mult[nmult+1], 0., 100.));
    Levy[nmult] = LevyTsallis(Form("Levy%s_%.0f-%.0f",fWhichParticle.Data(),mult[nmult],mult[nmult+1]), 1.32131 );
    Spectra[nmult] = (TH1F *) filemult[nmult]->Get(Form("fHistPt%s", fWhichParticle.Data()));
    SpectraRaw[nmult] = (TH1F *) filemult[nmult]->Get("lInvMassRealRawData/fHistPtRaw");
  }

  //Get input MC pT-shape
  TFile* fileMC = new TFile("MC15g3b1.root","READ");
  TList* clistMC = (TList*)fileMC->Get("PWGLF_StrVsMult_MC/cList");
  TH3D* h3DGenerated = (TH3D*)clistMC->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", fWhichParticle.Data()));
  //
  // find projection bins in rapidity
  Int_t rapBin_min = h3DGenerated->GetYaxis()->FindBin( -0.5+1.e-6 );
  Int_t rapBin_max = h3DGenerated->GetYaxis()->FindBin( +0.5-1.e-6 );
  // 
  // find projection bins in multiplicity
  Int_t multBin_min = h3DGenerated->GetZaxis()->FindBin( 0.+1.e-6 );
  Int_t multBin_max = h3DGenerated->GetZaxis()->FindBin( 100.-1.e-6 );
  // get the projection in pt
  TH1D* HistGen = (TH1D*)h3DGenerated->ProjectionX(Form("HistGen%s", fWhichParticle.Data()),rapBin_min,rapBin_max, multBin_min, multBin_max);
  TH1D* HistGenClone = (TH1D*)HistGen->Clone("HistGenClone");
  HistGenClone->Reset();
  //Rebin in pT
  TH1F *fHistMCCountbyptCasc  = new TH1F("fHistMCCountbyptCasc","Cascade MC count;p_{T} (GeV/c);Counts", ptbinnumb, ptbinlimits);
  Double_t temppt;
  for(long i = 1; i<HistGen->GetNbinsX()+1;i++){
    temppt = HistGen->GetXaxis()->GetBinCenter(i);
    for(long filling = 0; filling<HistGen->GetBinContent(i); filling++){
      fHistMCCountbyptCasc->Fill(temppt);
    }
  }

  for (int ibin = 1; ibin <= HistGen->GetNbinsX(); ibin++){
     HistGenClone->SetBinContent(ibin, HistGen->GetBinContent(ibin)/HistGen->GetBinWidth(ibin));
     HistGenClone->SetBinError(ibin, HistGen->GetBinError(ibin)/HistGen->GetBinWidth(ibin));
  }

  double integral = HistGenClone->Integral(0,-1);
  HistGenClone->Scale(1./integral);

  //Fit with Tsallis 
  Double_t par_n[multbinnumb];
  Double_t par_C[multbinnumb];
  Double_t par_norm[multbinnumb];
  HistGenClone->Fit(LevyMCMB,"","",0.,8.);
  SpectraMB->Fit(LevyRDMB,"","",0.,8.);
  for(int nmult = 0; nmult < multbinnumb; nmult++){ 
    Spectra[nmult]->Fit(Levy[nmult],"","",0.8,6.5);
    par_n[nmult] = Levy[nmult]->GetParameter(1);
    par_C[nmult] = Levy[nmult]->GetParameter(2);
    par_norm[nmult] = Levy[nmult]->GetParameter(3);
    Levy[nmult]->SetParameter(0,1.32131); 
    Levy[nmult]->SetParameter(1,par_n[nmult]);
    Levy[nmult]->SetParameter(2,par_C[nmult]);
    Levy[nmult]->SetParameter(3,par_norm[nmult]);
  }

//++++++++++++++++++++++++++++++++++++++++
// Generate histos from Fit Levy-Tsallis
//++++++++++++++++++++++++++++++++++++++++ 

  TH1F *hLevyGenerated[multbinnumb];
  TH1F *hLevyGeneratedMCMB = new TH1F(Form("hLevyGeneratedMC%s_%.0f-%.0f", fWhichParticle.Data(), 0., 100.), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", 100, 0., 10.);
  TH1F *hLevyGeneratedRDMB = new TH1F(Form("hLevyGeneratedRD%s_%.0f-%.0f", fWhichParticle.Data(), 0., 100.), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", 100, 0., 10.);

  for (int nmult = 0; nmult < multbinnumb; nmult++) {
    hLevyGenerated[nmult] = new TH1F(Form("hLevyGenerated%s_%.0f-%.0f", fWhichParticle.Data(), mult[nmult], mult[nmult + 1]), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", 100, 0., 10.);
  }

  const int n = 100000000; 
  //MB spectra
  for (Int_t i = 0; i < n; i++) {
    double xMC = (double)LevyMCMB->GetRandom(0., 10.);
    double xRD = (double)LevyRDMB->GetRandom(0., 10.);
    hLevyGeneratedRDMB->Fill(xRD);
    hLevyGeneratedMCMB->Fill(xMC);
  }
  // /dpT
  for (Int_t i = 0; i < 100; i++) {
    double entryRD = hLevyGeneratedRDMB->GetBinContent(i + 1) /
                     hLevyGeneratedRDMB->GetBinWidth(i + 1);
    hLevyGeneratedRDMB->SetBinContent(i + 1, entryRD);
    double entryMC = hLevyGeneratedMCMB->GetBinContent(i + 1) /
                     hLevyGeneratedMCMB->GetBinWidth(i + 1);
    hLevyGeneratedMCMB->SetBinContent(i + 1, entryMC);
  }

  // multiplicity spectra RD
  for (int nmult = 0; nmult < multbinnumb; nmult++) {
    for (Int_t i = 0; i < n; i++) {
      double x = (double)Levy[nmult]->GetRandom(0., 10.);
      hLevyGenerated[nmult]->Fill(x);
    }
    // /dpT
    for (Int_t i = 0; i < 100; i++) {
      double entry = hLevyGenerated[nmult]->GetBinContent(i + 1) /
                     hLevyGenerated[nmult]->GetBinWidth(i + 1);
      hLevyGenerated[nmult]->SetBinContent(i + 1, entry);
      hLevyGenerated[nmult]->SetMarkerStyle(8);
      hLevyGenerated[nmult]->SetMarkerSize(0.5);
    }
  }

  //Set colors
  hLevyGenerated[0]->SetLineColor(kRed + 1);
  hLevyGenerated[1]->SetLineColor(kRed - 4);
  hLevyGenerated[2]->SetLineColor(kOrange + 7);
  hLevyGenerated[3]->SetLineColor(kOrange - 3);
  hLevyGenerated[4]->SetLineColor(kYellow + 1);
  hLevyGenerated[5]->SetLineColor(kSpring - 7);
  hLevyGenerated[6]->SetLineColor(kGreen + 2);
  hLevyGenerated[7]->SetLineColor(kAzure + 8);
  hLevyGenerated[8]->SetLineColor(kBlue - 4);
  hLevyGenerated[9]->SetLineColor(kBlue + 3);

  hLevyGenerated[0]->SetMarkerColor(kRed + 1);
  hLevyGenerated[1]->SetMarkerColor(kRed - 4);
  hLevyGenerated[2]->SetMarkerColor(kOrange + 7);
  hLevyGenerated[3]->SetMarkerColor(kOrange - 3);
  hLevyGenerated[4]->SetMarkerColor(kYellow + 1);
  hLevyGenerated[5]->SetMarkerColor(kSpring - 7);
  hLevyGenerated[6]->SetMarkerColor(kGreen + 2);
  hLevyGenerated[7]->SetMarkerColor(kAzure + 8);
  hLevyGenerated[8]->SetMarkerColor(kBlue - 4);
  hLevyGenerated[9]->SetMarkerColor(kBlue + 3);

  //Draw
  TCanvas *c = new TCanvas("c", "", 800, 1000);
  c->SetLogy();
  hLevyGeneratedRDMB->SetMarkerStyle(8);
  hLevyGeneratedRDMB->SetMarkerSize(0.5);
  hLevyGeneratedRDMB->SetMarkerColor(kBlack);
  hLevyGeneratedRDMB->SetLineColor(kBlack);
  hLevyGeneratedMCMB->SetLineColor(kBlack);
  hLevyGeneratedMCMB->SetLineWidth(2);
  hLevyGeneratedRDMB->GetYaxis()->SetRangeUser(10000, 5E+8);
  hLevyGeneratedRDMB->GetYaxis()->SetTitle("dN/dp_{T}");
  hLevyGeneratedRDMB->GetXaxis()->SetRangeUser(0., 6.5);
  hLevyGeneratedRDMB->Draw("LEP");
  hLevyGeneratedMCMB->Draw("SAME L");
  for (int nmult = 0; nmult < multbinnumb; nmult++) {
    hLevyGenerated[nmult]->Draw("SAME LEP");
    // Levy[nmult]->Draw("SAME");
  }
  hLevyGeneratedRDMB->Draw("SAME LEP");
  hLevyGeneratedMCMB->Draw("SAME L");

//++++++++++++++++++++++++++++++++++++++++
// Iterative process for correction
//++++++++++++++++++++++++++++++++++++++++ 

  TH1F* hSpectraRatio[multbinnumb];

  for (int nmult = 0; nmult < multbinnumb; nmult++) {
    hSpectraRatio[nmult] = (TH1F*)hLevyGenerated[nmult]->Clone(Form("hLevyGenerated%s_%.0f-%.0f", fWhichParticle.Data(),
                      mult[nmult], mult[nmult + 1]));
    hSpectraRatio[nmult]->Reset();

    for (int bin = 1; bin < hSpectraRatio[0]->GetNbinsX(); bin ++){
      hSpectraRatio[nmult]->SetBinContent(bin, hLevyGenerated[nmult]->GetBinContent(bin)/hLevyGeneratedMCMB->GetBinContent(bin));
      hSpectraRatio[nmult]->SetBinError(bin, ErrorInRatio(hLevyGenerated[nmult]->GetBinContent(bin),hLevyGenerated[nmult]->GetBinError(bin),
                          hLevyGeneratedMCMB->GetBinContent(bin),hLevyGeneratedMCMB->GetBinError(bin)));
    }
  }

  TCanvas *it1 = new TCanvas("it1", "", 800, 1000);
  hSpectraRatio[0]->GetYaxis()->SetRangeUser(0.,2.);
  hSpectraRatio[0]->GetYaxis()->SetTitle("Corr.Spectrum(DATA) / Input(MC)");
  hSpectraRatio[0]->GetXaxis()->SetRangeUser(0., 6.5);
  hSpectraRatio[0]->Draw();
  for (int nmult = 1; nmult < multbinnumb; nmult++) {
    hSpectraRatio[nmult]->Draw("SAME");
  }


  TFile* lResultsFile = TFile::Open(Form("%sMCptshapeIT0.root",fWhichParticle.Data()), "RECREATE");
  fHistMCCountbyptCasc->Write();  
  HistGen->Write();
  HistGenClone->Write();
  lResultsFile->cd();
  TDirectoryFile *lSpectraGen = new TDirectoryFile("GeneratedSpectra","");
  lSpectraGen->cd();
  hLevyGeneratedMCMB->Write();
  hLevyGeneratedRDMB->Write();
  for(int nmult = 0; nmult < multbinnumb; nmult++){ 
  	hLevyGenerated[nmult]->Write();
  }
  lResultsFile->cd();
  TDirectoryFile *lSpectraRatio = new TDirectoryFile("RatioGenMCRD","");
  lSpectraRatio->cd();
  for(int nmult = 0; nmult < multbinnumb; nmult++){ 
    hSpectraRatio[nmult]->Write();
  }


}

//--------------------------------------------------------------
//------------------- OTHER FUNCTIONS --------------------------
//--------------------------------------------------------------

Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{ 
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4.);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e+3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e+3);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e+6);
  return fLevyTsallis;
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


