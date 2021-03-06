#include <fstream>
#include <iostream>

double ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr);
double ErrorInRatioCorr(Double_t A, Double_t Aerr, Double_t B, Double_t Berr);
void DivideAndComputeRogerBarlow(TH1F *h1, TH1F *h2);
TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t par2min = 1e-3, Double_t par2max = 1e+3,Double_t par3min = 1e-6, Double_t par3max = 1e+6, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);

//--------------------------------------------------------------
//------------------- MAIN FUNCTION ----------------------------
//--------------------------------------------------------------
void FitGenSpectra(int iteration = 0, TString fWhichEstimator = "ZDC", Double_t lLoMult = 0., 
Double_t lHiMult = 30., Double_t lLoEE = 0., Double_t lHiEE = 100.) {

    TString fWhichParticle = "XiMinus";
    TString fWhichAntiParticle = "XiPlus";
    TString fWhichOtherEstimator = "ZDC";

    //Get input MC pT-shape
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
    TH1D* HistGenClone = (TH1D*)HistGen->Clone("HistGenClone");
    HistGenClone->Reset();
    //Rebin
    for (int ibin = 1; ibin <= HistGen->GetNbinsX(); ibin++){
        HistGenClone->SetBinContent(ibin, HistGen->GetBinContent(ibin)/HistGen->GetBinWidth(ibin));
        HistGenClone->SetBinError(ibin, HistGen->GetBinError(ibin)/HistGen->GetBinWidth(ibin));
    }
    double integral = HistGenClone->Integral(0,-1);
    HistGenClone->Scale(1./integral);

    // Def multiplicity / eff energy
    Double_t percentileV0[] = {0.,1.,5,10,15,20,30,40,50,70,100};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    //
    Double_t percentileZDC[] = {0,20,30,40,50,60,70,80,90,100};
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
    //
    Double_t ptbinlimits[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
    Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;


//++++++++++++++++++++++++++++++++++++++++
// Fit spectra w/ Levy-Tsallis 
//++++++++++++++++++++++++++++++++++++++++ 

    //Get from files 
    TFile* filePart[multbinnumb], * fileAntiPart[multbinnumb];
    TFile* fileP = new TFile(Form("Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichParticle.Data(), 0.,100.,0.,100.));
    TFile* fileAntiP = new TFile(Form("Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichAntiParticle.Data(), 0.,100.,0.,100.));

    TFile* fileMB = new TFile(Form("../ResultsFiles/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", 
    fWhichParticle.Data(), lLoMult,lHiMult,lLoEE,lHiEE));
  
    TF1* Levy[multbinnumb];
    TF1* LevyMCMB = LevyTsallis(Form("LevyMC%s_%.0f_%.0f",fWhichParticle.Data(),0.,100.), 1.32131,1e-3,1e+3,1e-6,1e+6 );
    TF1* LevyRDMB = LevyTsallis(Form("LevyRD%s_%.0f_%.0f",fWhichParticle.Data(),0.,100.), 1.32131 );
  
    TH1F* SpectraPart[multbinnumb], *SpectraAntiPart[multbinnumb], *Spectra[multbinnumb], *HistReco;
    TH1F* SpectraMB = (TH1F *) fileMB->Get(Form("fHistPt%s", fWhichParticle.Data()));
  
    Levy[0] = LevyTsallis(Form("Levy%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",fWhichParticle.Data(),mult[0],mult[0+1],lLoEE,lHiEE), 1.32131 ,1e-3,1e+3,1e-6,1e+6);
    Levy[1] = LevyTsallis(Form("Levy%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",fWhichParticle.Data(),mult[1],mult[1+1],lLoEE,lHiEE), 1.32131 ,1e-3,1e+3,1e-6,1e+6);
    Levy[2] = LevyTsallis(Form("Levy%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",fWhichParticle.Data(),mult[2],mult[2+1],lLoEE,lHiEE), 1.32131 ,1e-3,1e+3,1e-6,1e+6);
    Levy[3] = LevyTsallis(Form("Levy%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",fWhichParticle.Data(),mult[3],mult[3+1],lLoEE,lHiEE), 1.32131 ,1e-3,1e+3,1e-6,1e+6);
    Levy[4] = LevyTsallis(Form("Levy%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",fWhichParticle.Data(),mult[4],mult[4+1],lLoEE,lHiEE), 1.32131 ,1e-3,1e+3,1e-6,1e+6);
    Levy[5] = LevyTsallis(Form("Levy%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",fWhichParticle.Data(),mult[5],mult[5+1],lLoEE,lHiEE), 1.32131 ,1e-3,1e+3,1e-6,1e+6);
    Levy[6] = LevyTsallis(Form("Levy%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",fWhichParticle.Data(),mult[6],mult[6+1],lLoEE,lHiEE), 1.32131 ,1e-3,1e+3,1e-6,1e+6);
    Levy[7] = LevyTsallis(Form("Levy%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",fWhichParticle.Data(),mult[7],mult[7+1],lLoEE,lHiEE), 1.32131 ,1e-3,1e+3,1e-6,1e+6);
    Levy[8] = LevyTsallis(Form("Levy%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",fWhichParticle.Data(),mult[8],mult[8+1],lLoEE,lHiEE), 1.32131 ,1e-3,1e+3,1e-6,1e+6);
    if (fWhichEstimator.Contains("V0M")) Levy[9] = LevyTsallis(Form("Levy%s_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",fWhichParticle.Data(),mult[9],mult[9+1],lLoEE,lHiEE), 1.32131 ,1e-2,1e+3,1e-6,1e+6);
    
    for(int nmult = 0; nmult < multbinnumb; nmult++)
    {
        double zdcmin = lLoEE, zdcmax = lHiEE, v0min = lLoMult, v0max = lHiMult;
        if (fWhichEstimator.Contains("V0M")) {
            v0min = mult[nmult];
            v0max = mult[nmult+1];
        }
        if (fWhichEstimator.Contains("ZDC")) {
            zdcmin = mult[nmult];
            zdcmax = mult[nmult+1];
        }

        filePart[nmult] = new TFile(Form("../ResultsFiles/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichParticle.Data(), v0min, v0max, zdcmin, zdcmax));
        SpectraPart[nmult] = (TH1F *) filePart[nmult]->Get(Form("fHistPt%s", fWhichParticle.Data()));
        fileAntiPart[nmult] = new TFile(Form("../ResultsFiles/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichAntiParticle.Data(), v0min, v0max, zdcmin, zdcmax));
        SpectraAntiPart[nmult] = (TH1F *) fileAntiPart[nmult]->Get(Form("fHistPt%s", fWhichAntiParticle.Data()));
        
        //Sum
        Spectra[nmult] = (TH1F*)SpectraPart[nmult]->Clone(Form("fHistPtXi%i",nmult));
        Spectra[nmult]->Reset();
        for (int bin = 1 ; bin <=  SpectraPart[nmult]->GetNbinsX(); bin ++ ){
            Spectra[nmult]->SetBinContent(bin, SpectraPart[nmult]->GetBinContent(bin) + SpectraAntiPart[nmult]->GetBinContent(bin));
            Spectra[nmult]->SetBinError(bin, TMath::Sqrt(SpectraPart[nmult]->GetBinError(bin)*SpectraPart[nmult]->GetBinError(bin)+SpectraAntiPart[nmult]->GetBinError(bin)*SpectraAntiPart[nmult]->GetBinError(bin)));
        }
    }

    TH1D* HistRecoPart = (TH1D*)fileP->Get("fHistReco");
    TH1D* HistRecoAntiPart = (TH1D*)fileAntiP->Get("fHistReco");
    HistReco = (TH1F*) HistRecoPart->Clone("HistReco");
    HistReco->Reset();
    for (int bin = 1 ; bin <=  HistRecoPart->GetNbinsX(); bin ++ ){
       HistReco->SetBinContent(bin, HistRecoPart->GetBinContent(bin) + HistRecoAntiPart->GetBinContent(bin));
       HistReco->SetBinError(bin, TMath::Sqrt(HistRecoPart->GetBinError(bin)*HistRecoPart->GetBinError(bin)+HistRecoAntiPart->GetBinError(bin)*HistRecoAntiPart->GetBinError(bin)));
    }
           

    //Fit Spectra Xi
    Double_t par_n[multbinnumb];
    Double_t par_C[multbinnumb];
    Double_t par_norm[multbinnumb];
    HistGenClone->Fit(LevyMCMB,"","",0.,10.);
    SpectraMB->Fit(LevyRDMB,"","",0.6,6.5);
    for(int nmult = 0; nmult < multbinnumb; nmult++){ 
        Int_t fitres;
        Int_t trials = 0;
        trials = 0;
        do {
        fitres = Spectra[nmult]->Fit(Levy[nmult], "","",0.6,6.5);
        Printf("Trial: %d", trials++);
        if(trials > 10) {
            Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
            break;
        }
        }
        while (fitres != 0);
        //Spectra[nmult]->Fit(Levy[nmult],"","",0.,8.);
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
  TH1F *hLevyGeneratedMCMB = new TH1F(Form("hLevyGeneratedMC%s_%.0f_%.0f", fWhichParticle.Data(), 0., 100.), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", 250, 0.,25);
  TH1F *hLevyGeneratedRDMB = new TH1F(Form("hLevyGeneratedRD%s_%.0f_%.0f", fWhichParticle.Data(), 0., 100.), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", 250, 0.,25);

  for (int nmult = 0; nmult < multbinnumb; nmult++) {
    hLevyGenerated[nmult] = new TH1F(Form("hLevyGenerated%s_V0M_%.0f_%.0f_ZDC_%03.0f_%03.0f", fWhichParticle.Data(), mult[nmult], mult[nmult + 1],lLoEE,lHiEE), "Cascade generated from Levy-Tsallis Fit;p_{T} (GeV/c);Counts", 250, 0, 25);
  }

  const int n = 10000000; 
  //MB spectra
  for (Int_t i = 0; i < n; i++) {
    double xMC = (double)LevyMCMB->GetRandom(0, 25);
    double xRD = (double)LevyRDMB->GetRandom(0, 25);
    hLevyGeneratedRDMB->Fill(xRD);
    hLevyGeneratedMCMB->Fill(xMC);
  }
  // /dpT
  for (Int_t i = 0; i < hLevyGeneratedRDMB->GetNbinsX(); i++) {
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
      double x = (double)Levy[nmult]->GetRandom(0, 25);
      hLevyGenerated[nmult]->Fill(x);
    }
    // /dpT
    for (Int_t i = 0; i <hLevyGenerated[nmult]->GetNbinsX() ; i++) {
      double entry = hLevyGenerated[nmult]->GetBinContent(i + 1) /
                     hLevyGenerated[nmult]->GetBinWidth(i + 1);
      hLevyGenerated[nmult]->SetBinContent(i + 1, entry);
      hLevyGenerated[nmult]->SetMarkerStyle(8);
      hLevyGenerated[nmult]->SetMarkerSize(0.9);
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
   if (fWhichEstimator.Contains("V0M")) hLevyGenerated[9]->SetLineColor(kBlue + 3);

  hLevyGenerated[0]->SetMarkerColor(kRed + 1);
  hLevyGenerated[1]->SetMarkerColor(kRed - 4);
  hLevyGenerated[2]->SetMarkerColor(kOrange + 7);
  hLevyGenerated[3]->SetMarkerColor(kOrange - 3);
  hLevyGenerated[4]->SetMarkerColor(kYellow + 1);
  hLevyGenerated[5]->SetMarkerColor(kSpring - 7);
  hLevyGenerated[6]->SetMarkerColor(kGreen + 2);
  hLevyGenerated[7]->SetMarkerColor(kAzure + 8);
  hLevyGenerated[8]->SetMarkerColor(kBlue - 4);
   if (fWhichEstimator.Contains("V0M")) hLevyGenerated[9]->SetMarkerColor(kBlue + 3);

  //Draw
  TCanvas *c = new TCanvas("c", "", 800, 800);
  c->SetRightMargin(0.08);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);
  c->SetLogy();
  hLevyGeneratedRDMB->SetMarkerStyle(8);
  hLevyGeneratedRDMB->SetMarkerSize(1.1);
  hLevyGeneratedRDMB->SetMarkerColor(kBlack);
  hLevyGeneratedRDMB->SetLineColor(kBlack);
  hLevyGeneratedMCMB->SetLineColor(kBlack);
  hLevyGeneratedMCMB->SetLineWidth(3);
  hLevyGeneratedRDMB->GetYaxis()->SetRangeUser(1000, 5E+7);
  hLevyGeneratedRDMB->GetYaxis()->SetTitle("dN/dp_{T} (a.u.)");
  hLevyGeneratedRDMB->GetYaxis()->SetTitleOffset(1.2);
  hLevyGeneratedRDMB->GetXaxis()->SetRangeUser(0.6, 6.5);
  hLevyGeneratedRDMB->SetStats(0);
  hLevyGeneratedRDMB->SetTitle("");
  hLevyGeneratedRDMB->Draw("LEP");
  hLevyGeneratedMCMB->Draw("SAME L");
  for (int nmult = 0; nmult < multbinnumb; nmult++) {
    hLevyGenerated[nmult]->Draw("SAME LEP");
    // Levy[nmult]->Draw("SAME");
  }
  hLevyGeneratedRDMB->Draw("SAME LEP");
  hLevyGeneratedMCMB->Draw("SAME L");
  TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
  xlabel-> SetNDC();
  xlabel-> SetTextColor(1);
  xlabel-> SetTextSize(0.08);
  xlabel-> SetTextAlign(22);
  xlabel-> SetTextAngle(0);
  TString text = "#Xi";
  xlabel-> DrawLatex(0.85, 0.85, text.Data());
  TLegend* l = new TLegend(0.5,0.75,0.7,0.89);
  l->SetTextSize(0.02);
  l->SetBorderSize(0);  
  //l->SetHeader("V0M selection (percentiles):");
  l->AddEntry(hLevyGeneratedMCMB,"MC input","L");
  if (fWhichEstimator.Contains("V0M")) {
    l->AddEntry(hLevyGeneratedRDMB,"V0M 0-100% corr. spectrum","LEP");
    l->AddEntry(hLevyGenerated[0],"V0M 0-1% corr. spectrum","LEP");
    l->AddEntry(hLevyGenerated[9],"V0M 70-100% corr. spectrum","LEP");
  } else if (fWhichEstimator.Contains("ZDC")) {
    l->AddEntry(hLevyGeneratedRDMB,"ZDC 0-100% corr. spectrum","LEP");
    l->AddEntry(hLevyGenerated[0],"ZDC 0-20% corr. spectrum","LEP");
    l->AddEntry(hLevyGenerated[8],"ZDC 90-100% corr. spectrum","LEP");
  }
  l->Draw("SAME");
  TLatex *xlabe = new TLatex();
	xlabe->SetTextFont(42);
  xlabe-> SetNDC();
  xlabe-> SetTextColor(1);
  xlabe-> SetTextSize(0.03);
  xlabe-> SetTextAlign(22);
  xlabe-> SetTextAngle(0);
 xlabe-> DrawLatex(0.78, 0.71, Form("%s fixed [%.0f-%.0f]",fWhichOtherEstimator.Data(),lLoMult,lHiMult));
  c->SaveAs(Form("images/%s/SpectraMCpTCorrection%sIT%i_ZDC_%03.0f_%03.0f.png",fWhichEstimator.Data(),"Xi",iteration,lLoMult,lHiMult));

//++++++++++++++++++++++++++++++++++++++++
// Iterative process for correction
//++++++++++++++++++++++++++++++++++++++++ 

  TH1F* hSpectraRatio[multbinnumb];

  for (int nmult = 0; nmult < multbinnumb; nmult++) {
    hSpectraRatio[nmult] = (TH1F*)hLevyGenerated[nmult]->Clone(Form("hLevyGenerated%s_V0M_%.0f_%.0f_ZDC_%03.0f_%03.0f", fWhichParticle.Data(),
                      mult[nmult], mult[nmult + 1],lLoEE,lHiEE));
    hSpectraRatio[nmult]->Reset();

    for (int bin = 1; bin < hSpectraRatio[0]->GetNbinsX(); bin ++){
      //DivideAndComputeRogerBarlow(hSpectraRatio[nmult],hLevyGeneratedMCMB);
      if (hLevyGeneratedMCMB->GetBinContent(bin) != 0 )
        {
      hSpectraRatio[nmult]->SetBinContent(bin, hLevyGenerated[nmult]->GetBinContent(bin)/hLevyGeneratedMCMB->GetBinContent(bin));
      hSpectraRatio[nmult]->SetBinError(bin, ErrorInRatio(hLevyGenerated[nmult]->GetBinContent(bin),hLevyGenerated[nmult]->GetBinError(bin),
                          hLevyGeneratedMCMB->GetBinContent(bin),hLevyGeneratedMCMB->GetBinError(bin)));
        }else hSpectraRatio[nmult]->SetBinContent(bin,1);
    }
  }

//++++++++++++++++++++++++++++++++++++++++
// Correct Efficiency
//++++++++++++++++++++++++++++++++++++++++ 

  TH1F* HistRecoMult[multbinnumb], *HistGenMult[multbinnumb], *HistRecoMult_rebin[multbinnumb], *HistGenMult_rebin[multbinnumb], *CorrEfficiency[multbinnumb];
 
  for (int nmult = 0; nmult < multbinnumb; nmult++) {

    HistRecoMult_rebin[nmult]   = new TH1F(Form("HistRecoMult_rebin%i",nmult),"Cascade MC count;p_{T} (GeV/c);Counts", ptbinnumb, ptbinlimits);
    HistGenMult_rebin[nmult]   = new TH1F(Form("HistGenMult_rebin%i",nmult),"Cascade MC count;p_{T} (GeV/c);Counts", ptbinnumb, ptbinlimits);

    HistRecoMult[nmult] = (TH1F*) HistReco->Clone(Form("HistRecoMult%i",nmult));
    HistRecoMult[nmult] ->Reset();
    HistGenMult[nmult] = (TH1F*) HistGen->Clone(Form("HistGenMult%i",nmult));
    HistGenMult[nmult] ->Reset();
    for (int bin = 1; bin <= HistRecoMult[nmult]->GetNbinsX(); bin ++){
      //
      HistRecoMult[nmult]-> SetBinContent(bin,HistReco->GetBinContent(bin)*hSpectraRatio[nmult]->GetBinContent(bin));
      HistRecoMult[nmult]-> SetBinError(bin,ErrorInRatio(HistReco->GetBinContent(bin), HistReco->GetBinError(bin),
        hSpectraRatio[nmult]->GetBinContent(bin), hSpectraRatio[nmult]->GetBinError(bin))
        );
      //
      HistGenMult[nmult]-> SetBinContent(bin, HistGen->GetBinContent(bin)*hSpectraRatio[nmult]->GetBinContent(bin));
      HistGenMult[nmult]-> SetBinError(bin, ErrorInRatio(HistGen->GetBinContent(bin), HistGen->GetBinError(bin),
        hSpectraRatio[nmult]->GetBinContent(bin), hSpectraRatio[nmult]->GetBinError(bin)
        ));
    }

    Double_t tempptreco;
    for(long i = 1; i<HistRecoMult[nmult]->GetNbinsX()+1;i++){
      tempptreco = HistRecoMult[nmult]->GetXaxis()->GetBinCenter(i);
      for(long filling = 0; filling<HistRecoMult[nmult]->GetBinContent(i); filling++){
        HistRecoMult_rebin[nmult]->Fill(tempptreco);
      }
    }
    Double_t tempptgen;
    for(long i = 1; i<HistGenMult[nmult]->GetNbinsX()+1;i++){
      tempptgen = HistGenMult[nmult]->GetXaxis()->GetBinCenter(i);
      for(long filling = 0; filling<HistGenMult[nmult]->GetBinContent(i); filling++){
        HistGenMult_rebin[nmult]->Fill(tempptgen);
      }
    }
   
  
    CorrEfficiency[nmult] = (TH1F*)HistRecoMult_rebin[nmult]->Clone(Form("CorrEfficiency%i",nmult));
    CorrEfficiency[nmult]->Reset();
    for (int bin = 1; bin <= HistRecoMult_rebin[nmult]->GetNbinsX(); bin ++){
      if (HistGenMult[nmult]->GetBinContent(bin) != 0){
        CorrEfficiency[nmult]->SetBinContent(bin , HistRecoMult_rebin[nmult]->GetBinContent(bin)/HistGenMult_rebin[nmult]->GetBinContent(bin));
        CorrEfficiency[nmult]->SetBinError(bin, ErrorInRatioCorr(HistRecoMult_rebin[nmult]->GetBinContent(bin),HistRecoMult_rebin[nmult]->GetBinError(bin),
        HistGenMult_rebin[nmult]->GetBinContent(bin),HistGenMult_rebin[nmult]->GetBinError(bin)));
      }
    }
   
  }

  TFile* lResultsFile = TFile::Open(Form("%sMCptshape%sIT%i_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi",fWhichEstimator.Data(),iteration,lLoMult, lHiMult,lLoEE,lHiEE), "RECREATE");
  
  lResultsFile->cd();
  TDirectoryFile *lEff = new TDirectoryFile("Efficiencies","");
  lEff->cd();
  for(int nmult = 0; nmult < multbinnumb; nmult++){ 
    CorrEfficiency[nmult]->Write();
  }
  lResultsFile->cd();
  TDirectoryFile *lSpectraCorr = new TDirectoryFile("SpectraCorr","");
  lSpectraCorr->cd();
  for(int nmult = 0; nmult < multbinnumb; nmult++){ 
    Spectra[nmult]->Write();
  }
  lResultsFile->cd();
  TDirectoryFile *lGen = new TDirectoryFile("HistGen","");
  lGen->cd();
  for(int nmult = 0; nmult < multbinnumb; nmult++){ 
    HistGenMult[nmult]->Write();
  }
  lResultsFile->cd();
  TDirectoryFile *lReco = new TDirectoryFile("HistReco","");
  lReco->cd();
  for(int nmult = 0; nmult < multbinnumb; nmult++){ 
    HistRecoMult[nmult]->Write();
  }
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
  

  TCanvas *it1 = new TCanvas("it1", "", 800, 800);
  it1->SetRightMargin(0.08);
  it1->SetLeftMargin(0.15);
  it1->SetBottomMargin(0.15);
  it1->SetGridy();
  hSpectraRatio[0]->SetStats(0);
  hSpectraRatio[0]->SetTitle("");
  hSpectraRatio[0]->GetYaxis()->SetRangeUser(0.,3.);
  hSpectraRatio[0]->GetYaxis()->SetTitle("Corr.Spectrum(DATA) / Input(MC)");
  hSpectraRatio[0]->GetXaxis()->SetRangeUser(0.6, 6.5);
  hSpectraRatio[0]->Draw();
  for (int nmult = 1; nmult < multbinnumb; nmult++) {
    hSpectraRatio[nmult]->Draw("SAME");
  }
  xlabel-> DrawLatex(0.85, 0.85, text.Data());
  //
  TLatex *xlabel2 = new TLatex();
	xlabel2->SetTextFont(42);
  xlabel2-> SetNDC();
  xlabel2-> SetTextColor(1);
  xlabel2-> SetTextSize(0.04);
  xlabel2-> SetTextAlign(22);
  xlabel2-> SetTextAngle(0);
  xlabel2-> DrawLatex(0.32, 0.81, Form("Iteration %i",0));
  TLegend* l1 = new TLegend(0.5,0.75,0.7,0.89);
  l1->SetTextSize(0.02);
  l1->SetBorderSize(0);  
  //l1->SetHeader("V0M selection (percentiles):");
  if (fWhichEstimator.Contains("V0M")) {
    l1->AddEntry(hSpectraRatio[0],"V0M 0-1% corr. spectrum","LEP");
    l1->AddEntry(hSpectraRatio[9],"V0M 70-100% corr. spectrum","LEP");
  } else if (fWhichEstimator.Contains("ZDC")) {
    l1->AddEntry(hSpectraRatio[0],"ZDC 0-20% corr. spectrum","LEP");
    l1->AddEntry(hSpectraRatio[8],"ZDC 90-100% corr. spectrum","LEP");
  }
  
  l1->Draw("SAME");
  if (lLoEE !=0 || lHiEE !=100)xlabe-> DrawLatex(0.32, 0.71, Form("%s fixed [%.0f-%.0f]",fWhichOtherEstimator.Data(),lLoEE,lHiEE));
  it1->SaveAs(Form("images/%s/SpectraRatio%sIT%i_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.png",fWhichEstimator.Data(),"Xi",iteration,lLoMult, lHiMult, lLoEE,lHiEE));

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

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t par2min = 1e-3, Double_t par2max = 1e+3,Double_t par3min = 1e-6, Double_t par3max = 1e+6, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{ 
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 25., 4.);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e+3);
  fLevyTsallis->SetParLimits(2, par2min, par2max);
  fLevyTsallis->SetParLimits(3, par3min, par3max);
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

//---------------------------------------------------------------------------------------------------
double ErrorInRatioCorr ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop - errorfrombottom) );
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