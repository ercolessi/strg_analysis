#include <iostream>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TPad.h"
#include "TPaveText.h"
#include <TLatex.h>
//---------------------------------------------------------------------------------------------------
double ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr);
void DivideAndComputeRogerBarlow(TH1F *h1, TH1F *h2);
TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t par2min = 1e-3, Double_t par2max = 1e+3,Double_t par3min = 1e-6, Double_t par3max = 1e+6, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);

void DrawSpectraWithSyst(TString fWhichEstimator = "V0M", Double_t lLo = 0., Double_t lHi = 30.){

    TString fWhichOtherEstimator = "";
    if (fWhichEstimator.Contains("V0M")) fWhichOtherEstimator = "ZDC";
    else fWhichOtherEstimator = "V0M";
    TString outputname = Form("CorrectedSpectraV0MXiTOT_ZDC_%03.0f_%03.0f.root",lLo, lHi);

    Double_t percentileV0[] = {0.,5,10,15,20,30,40,50,70,100};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    Double_t percentileZDC[] = {0,20,30,40,50,60,70,80,90,100};
    const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;
    Double_t * mult;
    TString namesystsgnloss;
    //
    int tempbinnumber = 0;
    if (fWhichEstimator.Contains("V0M")) {
        tempbinnumber = nbinV0;
        mult = percentileV0;
        if (lLo > 0) {
           namesystsgnloss = Form("%s-%s","multsel_fixedlowEE","V0FixLowEE");  
        }
        if (lHi < 100) {
            namesystsgnloss = Form("%s-%s","multsel_fixedhighEE","V0FixHighEE");  
        }
    } 
    else if (fWhichEstimator.Contains("ZDC")) {
        tempbinnumber = nbinZDC;
        mult = percentileZDC;
        namesystsgnloss = Form("%s-%s","EEsel","ZDC");
         if (lLo > 0) {
            namesystsgnloss = Form("%s-%s","EEsel_fixedlowmult","ZDCFixLowmult");
        }
        if (lHi < 100) {
            namesystsgnloss = Form("%s-%s","EEsel_fixedhighmult","ZDCFixHighmult");
        }
    } 
    else {cout << "No valid name for estimator... its V0M or ZDC" << endl; return;}
    const int multbinnumb = tempbinnumber;
	//
	Double_t ptbinlimits[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
	Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;
    //
	TH1D* SpectraStat[multbinnumb],*hSystPart[3],*SpectraSyst[multbinnumb];
	TFile* fileStat = new TFile(Form("StatSpectra%s-Xi-%s_%03.0f_%03.0f.root",fWhichEstimator.Data(),fWhichOtherEstimator.Data(),lLo,lHi));
    TFile* fileSyst[3];
    TF1* levy[multbinnumb];
    TFile* fileSgnLossSyst = TFile::Open("NormalizationCorrections/SgnLossSyst.root");
    TH1D* hsgnloss[multbinnumb];
	
    for(int nmult = 0; nmult < multbinnumb; nmult++){
        double zdcmin = lLo, zdcmax = lHi, v0min = lLo, v0max = lHi;
		if (fWhichEstimator.Contains("V0M")) {
			v0min = mult[nmult];
			v0max = mult[nmult+1];
		}
		if (fWhichEstimator.Contains("ZDC")) {
			zdcmin = mult[nmult];
			zdcmax = mult[nmult+1];
		}		
        // levy
        levy[nmult] = LevyTsallis(Form("Levy%s_%.0f-%.0f","Xi",mult[nmult],mult[nmult+1]), 1.32131 ,1e-3,1e+2,1e-6,1e+6);
        // ---- stats -----
        SpectraStat[nmult] = (TH1D *) fileStat->Get(Form("XiSpectra_Stats_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",v0min,v0max,zdcmin,zdcmax));
        SpectraStat[nmult]->SetStats(0);
    
    }
    //

    //Syst loop
    double systmult[4] = //{0,10,30,100};
{0.,30.,60.,100.};
    for(int nmult = 0; nmult < 3; nmult++){
        double zdcmin = lLo, zdcmax = lHi, v0min = lLo, v0max = lHi;
		if (fWhichEstimator.Contains("V0M")) {
			v0min = systmult[nmult];
			v0max = systmult[nmult+1];
		}
		if (fWhichEstimator.Contains("ZDC")) {
			zdcmin = systmult[nmult];
			zdcmax = systmult[nmult+1];
		}	
         // files
	    fileSyst[nmult] = TFile::Open(Form("sistematiche/SystematicsFinalResults-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",
            "Xi",0.,100.,0.,100.));
            //v0min,v0max,zdcmin,zdcmax));	
        // histo x systematics
		hSystPart[nmult] = (TH1D*)fileSyst[nmult]->Get("hSystTot");
    }

    for(int nmult = 0; nmult < multbinnumb; nmult++){
        double zdcmin = lLo, zdcmax = lHi, v0min = lLo, v0max = lHi;
		if (fWhichEstimator.Contains("V0M")) {
			v0min = mult[nmult];
			v0max = mult[nmult+1];
		}
		if (fWhichEstimator.Contains("ZDC")) {
			zdcmin = mult[nmult];
			zdcmax = mult[nmult+1];
		}	
        // ---- syst ------
        // do spectra syst
        hsgnloss[nmult] = (TH1D*)fileSgnLossSyst->Get(Form("%s/hRatioClone%i",namesystsgnloss.Data(),nmult));
        SpectraSyst[nmult] = (TH1D*)SpectraStat[nmult]->Clone(Form("XiSpectra_Syst_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            v0min,v0max,zdcmin,zdcmax));
        //ZDC 3-6-10 V0M 4-7-11
   		for (int jbin = 1; jbin <= SpectraSyst[nmult]->GetNbinsX(); jbin ++){
            SpectraSyst[nmult]->SetBinError(jbin,(SpectraSyst[nmult]->GetBinContent(jbin)*
                                                    TMath::Sqrt( hSystPart[0]->GetBinContent(jbin)*hSystPart[0]->GetBinContent(jbin) +
                                                      TMath::Abs(hsgnloss[nmult]->GetBinContent(jbin+1)-1)*TMath::Abs(hsgnloss[nmult]->GetBinContent(jbin+1)-1)
                                                    )
                ));
           cout << TMath::Sqrt( hSystPart[0]->GetBinContent(jbin)*hSystPart[0]->GetBinContent(jbin) +
                                                      TMath::Abs(hsgnloss[nmult]->GetBinContent(jbin+1)-1)*TMath::Abs(hsgnloss[nmult]->GetBinContent(jbin+1)-1)
                                                    ) << endl;

   		/*	if (fWhichEstimator.Contains("V0M")) {
                if (nmult < 4) SpectraSyst[nmult]->SetBinError(jbin,SpectraSyst[nmult]->GetBinContent(jbin)*(hSystPart[0]->GetBinContent(jbin)));
                if (nmult > 4 && nmult < 7) SpectraSyst[nmult]->SetBinError(jbin,SpectraSyst[nmult]->GetBinContent(jbin)*(hSystPart[1]->GetBinContent(jbin)));
                if (nmult > 7 && nmult < 11) SpectraSyst[nmult]->SetBinError(jbin,SpectraSyst[nmult]->GetBinContent(jbin)*(hSystPart[2]->GetBinContent(jbin)));
            }
            if (fWhichEstimator.Contains("ZDC")) {
                if (nmult < 3) SpectraSyst[nmult]->SetBinError(jbin,SpectraSyst[nmult]->GetBinContent(jbin)*(hSystPart[0]->GetBinContent(jbin)));
                if (nmult > 3 && nmult < 6) SpectraSyst[nmult]->SetBinError(jbin,SpectraSyst[nmult]->GetBinContent(jbin)*(hSystPart[1]->GetBinContent(jbin)));
                if (nmult > 6 && nmult < 10) SpectraSyst[nmult]->SetBinError(jbin,SpectraSyst[nmult]->GetBinContent(jbin)*(hSystPart[2]->GetBinContent(jbin)));
          }*/
        }
	}

   	TH1D* LevyFitSpectra[multbinnumb];
   	for(int nmult = 0; nmult < multbinnumb; nmult++){
        double zdcmin = lLo, zdcmax = lHi, v0min = lLo, v0max = lHi;
		if (fWhichEstimator.Contains("V0M")) {
			v0min = mult[nmult];
			v0max = mult[nmult+1];
		}
		if (fWhichEstimator.Contains("ZDC")) {
			zdcmin = mult[nmult];
			zdcmax = mult[nmult+1];
		}	
		LevyFitSpectra[nmult] = (TH1D*)SpectraStat[nmult]->Clone(Form("LevyFitSpectra-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            v0min,v0max,zdcmin,zdcmax ));
        //
		Int_t fitres;
		Int_t trials = 0;
		trials = 0;
		do {
		if (fWhichEstimator.Contains("ZDC") && lLo==70){
            if (nmult == 1) {fitres = LevyFitSpectra[nmult]->Fit(levy[nmult],"","",0.6,5);}
            else fitres = LevyFitSpectra[nmult]->Fit(levy[nmult],"","",0.6,5);
        } else fitres = LevyFitSpectra[nmult]->Fit(levy[nmult],"","",0.6,6.5);
        cout << levy[nmult]->GetChisquare() << " " << levy[nmult]->GetNDF() << endl;
		Printf("Trial: %d", trials++);
		if(trials > 10) {
			Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
			break;
		}
		}
		while (fitres != 0);	
		levy[nmult]->SetLineStyle(9);
        levy[nmult]->SetLineColor(kBlack);
        levy[nmult]->SetLineWidth(1);
        SpectraStat[nmult]->SetMarkerSize(1.8);
        SpectraSyst[nmult]->SetMarkerSize(1.8);

   	}

    /*SpectraStat[0]->Scale(64);
	SpectraStat[1]->Scale(32);
	SpectraStat[2]->Scale(16);
	SpectraStat[3]->Scale(8);
	SpectraStat[4]->Scale(4);
	SpectraStat[5]->Scale(2);
	SpectraStat[6]->Scale(1);
	SpectraStat[7]->Scale(0.5);
	SpectraStat[8]->Scale(0.25);
	if (fWhichEstimator.Contains("V0M")) SpectraStat[9]->Scale(0.10);*/

    SpectraStat[0]->SetLineColor(kRed + 1);
	SpectraStat[1]->SetLineColor(kRed - 4);
	SpectraStat[2]->SetLineColor(kOrange + 7);
	SpectraStat[3]->SetLineColor(kOrange - 3);
	SpectraStat[4]->SetLineColor(kYellow + 1);
	SpectraStat[5]->SetLineColor(kSpring - 7);
	SpectraStat[6]->SetLineColor(kGreen + 2);
	SpectraStat[7]->SetLineColor(kAzure + 8);
	if ((fWhichEstimator.Contains("ZDC") && lLo == 0)||fWhichEstimator.Contains("V0M"))SpectraStat[8]->SetLineColor(kBlue - 4);
	
	SpectraStat[0]->SetMarkerColor(kRed + 1);
	SpectraStat[1]->SetMarkerColor(kRed - 4);
	SpectraStat[2]->SetMarkerColor(kOrange + 7);
	SpectraStat[3]->SetMarkerColor(kOrange - 3);
	SpectraStat[4]->SetMarkerColor(kYellow + 1);
	SpectraStat[5]->SetMarkerColor(kSpring - 7);
	SpectraStat[6]->SetMarkerColor(kGreen + 2);
	SpectraStat[7]->SetMarkerColor(kAzure + 8);
	if ((fWhichEstimator.Contains("ZDC") && lLo == 0)||fWhichEstimator.Contains("V0M"))SpectraStat[8]->SetMarkerColor(kBlue - 4);
	
    /*SpectraSyst[0]->Scale(64);
	SpectraSyst[1]->Scale(32);
	SpectraSyst[2]->Scale(16);
	SpectraSyst[3]->Scale(8);
	SpectraSyst[4]->Scale(4);
	SpectraSyst[5]->Scale(2);
	SpectraSyst[6]->Scale(1);
	SpectraSyst[7]->Scale(0.5);
	SpectraSyst[8]->Scale(0.25);
	if (fWhichEstimator.Contains("V0M")) SpectraSyst[9]->Scale(0.10);*/

    SpectraSyst[0]->SetLineColor(kRed + 1);
	SpectraSyst[1]->SetLineColor(kRed - 4);
	SpectraSyst[2]->SetLineColor(kOrange + 7);
	SpectraSyst[3]->SetLineColor(kOrange - 3);
	SpectraSyst[4]->SetLineColor(kYellow + 1);
	SpectraSyst[5]->SetLineColor(kSpring - 7);
	SpectraSyst[6]->SetLineColor(kGreen + 2);
	SpectraSyst[7]->SetLineColor(kAzure + 8);
	if ((fWhichEstimator.Contains("ZDC") && lLo == 0)||fWhichEstimator.Contains("V0M"))SpectraSyst[8]->SetLineColor(kBlue - 4);
	
	SpectraSyst[0]->SetMarkerColor(kRed + 1);
	SpectraSyst[1]->SetMarkerColor(kRed - 4);
	SpectraSyst[2]->SetMarkerColor(kOrange + 7);
	SpectraSyst[3]->SetMarkerColor(kOrange - 3);
	SpectraSyst[4]->SetMarkerColor(kYellow + 1);
	SpectraSyst[5]->SetMarkerColor(kSpring - 7);
	SpectraSyst[6]->SetMarkerColor(kGreen + 2);
	SpectraSyst[7]->SetMarkerColor(kAzure + 8);
	if ((fWhichEstimator.Contains("ZDC") && lLo == 0)||fWhichEstimator.Contains("V0M"))SpectraSyst[8]->SetMarkerColor(kBlue - 4);
	
    

    SpectraStat[0]->SetMarkerStyle(20+0);
    SpectraSyst[0]->SetMarkerStyle(20+0);
    SpectraStat[1]->SetMarkerStyle(20+1);
    SpectraSyst[1]->SetMarkerStyle(20+1);
     SpectraStat[2]->SetMarkerStyle(20+2);
    SpectraSyst[2]->SetMarkerStyle(20+2);
     SpectraStat[3]->SetMarkerStyle(20+3);
    SpectraSyst[3]->SetMarkerStyle(20+3);
     SpectraStat[4]->SetMarkerStyle(20+9);
    SpectraSyst[4]->SetMarkerStyle(20+9);
     SpectraStat[5]->SetMarkerStyle(20+13);
    SpectraSyst[5]->SetMarkerStyle(20+13);
     SpectraStat[6]->SetMarkerStyle(20+14);
    SpectraSyst[6]->SetMarkerStyle(20+14);
     SpectraStat[7]->SetMarkerStyle(47);
    SpectraSyst[7]->SetMarkerStyle(47);
    
   if ((fWhichEstimator.Contains("ZDC") && lLo == 0)||fWhichEstimator.Contains("V0M")){ SpectraSyst[8]->SetMarkerStyle(39);
     SpectraStat[8]->SetMarkerStyle(39);
     SpectraStat[8]->SetMarkerSize(2.2);
        SpectraSyst[8]->SetMarkerSize(2.2);
        SpectraStat[8]->SetMarkerSize(2.2);}
        SpectraSyst[4]->SetMarkerSize(2.2);
        SpectraStat[4]->SetMarkerSize(2.2);
        SpectraSyst[5]->SetMarkerSize(2.2);
        SpectraStat[5]->SetMarkerSize(2.2);
   

    /*if (fWhichEstimator.Contains("ZDC") && lLo==70){ 
        SpectraStat[1]->SetBinContent(1,0.);
        SpectraSyst[1]->SetBinContent(1,0.);
        }*/


  /*  levy[0]->SetLineColor(kRed+1);
    levy[1]->SetLineColor(kRed-4);
    levy[2]->SetLineColor(kOrange+7);
    levy[3]->SetLineColor(kOrange-3);
    levy[4]->SetLineColor(kYellow+1);
    levy[5]->SetLineColor(kSpring-7);
    levy[6]->SetLineColor(kGreen+2);
    levy[7]->SetLineColor(kAzure+8);
    levy[8]->SetLineColor(kBlue-4);*/
    //if (fWhichEstimator.Contains("V0M")) levy[9]->SetLineColor(kBlue+3);
 
 	Double_t n[multbinnumb];
     n[0]=64.;
     n[1]=32.;
     n[2]=16.;
     n[3]=8.;
     n[4]=4.;
     n[5]=2.;
     n[6]=1.;
     n[7]=0.5;
     if ((fWhichEstimator.Contains("ZDC") && lLo == 0)||fWhichEstimator.Contains("V0M")) n[8]=0.25;
   
    for(int nmult = 0; nmult < multbinnumb; nmult++){
        n[nmult] = TMath::Power(5,multbinnumb-1-nmult);
        SpectraSyst[nmult]->Scale(n[nmult]);
        SpectraStat[nmult]->Scale(n[nmult]);

        double par = levy[nmult]->GetParameter("norm");
        levy[nmult]->SetParameter("norm",par*(n[nmult]));
    
    }

	  TLegend* l = new TLegend(0.2,0.18,0.4,0.45);
      l->SetTextSize(0.022);
      l->SetHeader(Form("%s classes (percentiles):",fWhichEstimator.Data()));
      l->SetBorderSize(0);     

	  TCanvas* c = new TCanvas("c","",1100,1600);
	  TStyle* mcStyle = new TStyle("mcStyle","Francesca's Root Styles");  
      mcStyle->SetPadTickX(1); 
      mcStyle->SetPadTickY(1); 
      mcStyle->SetPalette(1,0); 
      mcStyle->cd();
	  c->SetRightMargin(0.09);
      c->SetLeftMargin(0.15);
      c->SetBottomMargin(0.15);
	  c->SetFillColor(kWhite);
	  c->SetLogy();
      SpectraStat[0]->GetYaxis()->SetRangeUser(1E-12,10000000);
      if (fWhichEstimator.Contains("ZDC") && lLo==70){ SpectraStat[0]->GetXaxis()->SetRangeUser(0.,5);}
      else SpectraStat[0]->GetXaxis()->SetRangeUser(0.,6.5);
      SpectraStat[0]->GetYaxis()->SetTitle("1/N_{ev} d^{2}N/(dp_{T}dy) [(GeV/c)^{-1}]");
	  SpectraStat[0]->Draw();
      levy[0]->Draw("SAME");
      SpectraSyst[0]->SetFillStyle(0.);
	  SpectraSyst[0]->Draw("SAME E2");
      l->AddEntry(SpectraStat[0],Form("0-5 (x 5^{%i})",0),"LEP");	
	  for(int nmult = 1; nmult < multbinnumb; nmult++){
       
		SpectraSyst[nmult]->SetFillStyle(0.);
		SpectraSyst[nmult]->Draw("SAME E2");
		SpectraStat[nmult]->SetTitle(Form("%.0f-%.0f (x 5^{%i})",mult[nmult],mult[nmult+1],nmult));
	  	SpectraStat[nmult]->Draw("SAME");
        levy[nmult]->Draw("SAME");

        cout << nmult << endl;
	 
		l->AddEntry(SpectraStat[nmult],"","LEP");		
	  }
	  l->AddEntry(levy[multbinnumb-1],"Levy-Tsallis","L");
	  l->Draw("SAME");  
      // draw text
      TLatex *xlabel = new TLatex();
	  xlabel->SetTextFont(42);
      xlabel-> SetNDC();
      xlabel-> SetTextColor(1);
      xlabel-> SetTextSize(0.07);
      xlabel-> SetTextAlign(22);
      xlabel-> SetTextAngle(0);
      xlabel-> DrawLatex(0.75, 0.82, Form("#Xi^{-} + #bar{#Xi}^{+}"));
      TLatex *xlabe = new TLatex();
	  xlabe->SetTextFont(42);
      xlabe-> SetNDC();
      xlabe-> SetTextColor(1);
      xlabe-> SetTextSize(0.03);
      xlabe-> SetTextAlign(22);
      xlabe-> SetTextAngle(0);
      xlabe-> DrawLatex(0.25, 0.48, "|y|<0.5");
      TLatex *xlaber = new TLatex();
      xlaber->SetTextFont(42);
      xlaber-> SetNDC();
      xlaber-> SetTextColor(1);
      xlaber-> SetTextSize(0.03);
      xlaber-> SetTextAlign(22);
      xlaber-> SetTextAngle(0);
      if (lLo !=0 || lHi !=100) xlaber-> DrawLatex(0.6, 0.25, Form("%s fixed [%.0f-%.0f]",fWhichOtherEstimator.Data(),lLo,lHi));
        
        c->SaveAs(Form("images/Spectra%sSel_%sfix_%.0f-%.0f.png",fWhichEstimator.Data(),fWhichOtherEstimator.Data(),lLo,lHi));

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
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4.);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e+3);
  fLevyTsallis->SetParLimits(2, par2min, par2max);
  fLevyTsallis->SetParLimits(3, par3min, par3max);
  return fLevyTsallis;
}