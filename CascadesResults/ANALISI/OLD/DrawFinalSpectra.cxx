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

void DrawFinalSpectra(){

    Double_t percentileV0[] = {0.,1.,5,10,15,20,30,40,50,70,100};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    Double_t percentileZDC[] = {0,20,30,40,50,60,70,80,90,100};
    const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;
    //
	Double_t ptbinlimits[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
	Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;
    //
    // ZDC
    TString namesystsgnloss = Form("%s-%s","EEsel","ZDC");
    TFile* fileSgnLossSyst = TFile::Open("NormalizationCorrections/SgnLossSyst.root");
    TH1D* hsgnloss[nbinZDC];
    
	TH1D* SpectraStatZDC[nbinZDC],*hSystPartZDC,*SpectraSystZDC[nbinZDC];
	TFile* fileStatZDC = new TFile(Form("StatSpectra%s-Xi-%s_%03.0f_%03.0f.root","ZDC","V0M",0.,100.));
    if (!fileStatZDC) cout << "NO FILE" << endl;
    TFile* fileSystZDC = TFile::Open(Form("sistematiche/SystematicsFinalResults-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",
            "Xi",0.,100.,0.,100.));;
    TF1* levy[nbinZDC];
    //
    // 
    for(int nmult = 0; nmult < nbinZDC; nmult++){
        // levy
        levy[nmult] = LevyTsallis(Form("Levy%s_%.0f-%.0f","Xi",percentileZDC[nmult],percentileZDC[nmult+1]), 1.32131 ,1e-3,1e+3,1e-6,1e+6);
        // ---- stats -----
        SpectraStatZDC[nmult] = (TH1D *) fileStatZDC->Get(Form("XiSpectra_Stats_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",0.,100.,percentileZDC[nmult],percentileZDC[nmult+1]));
        SpectraStatZDC[nmult]->SetStats(0);
        if (!SpectraStatZDC[nmult]) cout << "NO Spectra" << endl;
        // ---- syst ------
        // do spectra syst
        hSystPartZDC = (TH1D*)fileSystZDC->Get("hSystTot");
        hsgnloss[nmult] = (TH1D*)fileSgnLossSyst->Get(Form("%s/hRatioClone%i",namesystsgnloss.Data(),nmult));


        SpectraSystZDC[nmult] = (TH1D*)SpectraStatZDC[nmult]->Clone(Form("XiSpectra_Syst_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            0.,100.,percentileZDC[nmult],percentileZDC[nmult+1]));
        for (int jbin = 1; jbin <= SpectraSystZDC[nmult]->GetNbinsX(); jbin ++){
            //
            SpectraSystZDC[nmult]->SetBinError(jbin,(SpectraSystZDC[nmult]->GetBinContent(jbin)*
                                                    TMath::Sqrt( hSystPartZDC->GetBinContent(jbin)*hSystPartZDC->GetBinContent(jbin) +
                                                      TMath::Abs(hsgnloss[nmult]->GetBinContent(jbin+1)-1)*TMath::Abs(hsgnloss[nmult]->GetBinContent(jbin+1)-1)
                                                    )
                ));
        }
        SpectraStatZDC[nmult]->SetMarkerStyle(20);
        SpectraSystZDC[nmult]->SetMarkerStyle(20);
        SpectraStatZDC[nmult]->SetMarkerSize(2.);
        SpectraSystZDC[nmult]->SetMarkerSize(2.);
    }
    //Fit ZDC spectra
   	TH1D* LevyFitSpectra[nbinZDC];
   	for(int nmult = 0; nmult < nbinZDC; nmult++){
        LevyFitSpectra[nmult] = (TH1D*)SpectraStatZDC[nmult]->Clone(Form("LevyFitSpectra-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            0.,100.,percentileZDC[nmult],percentileZDC[nmult+1] ));
        //
		Int_t fitres;
		Int_t trials = 0;
		trials = 0;
		do {
        fitres = LevyFitSpectra[nmult]->Fit(levy[nmult],"","",0.6,6.5);
        Printf("Trial: %d", trials++);
		if(trials > 10) {
			Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
			break;
		}}
		while (fitres != 0);	
		levy[nmult]->SetLineStyle(7);
        levy[nmult]->SetLineColor(kBlack);
        levy[nmult]->SetLineWidth(1);
   	}
      
    //
    SpectraStatZDC[0]->SetLineColor(kRed + 1);
	SpectraStatZDC[1]->SetLineColor(kRed - 4);
	SpectraStatZDC[2]->SetLineColor(kOrange + 7);
	SpectraStatZDC[3]->SetLineColor(kOrange - 3);
	SpectraStatZDC[4]->SetLineColor(kYellow + 1);
	SpectraStatZDC[5]->SetLineColor(kSpring - 7);
	SpectraStatZDC[6]->SetLineColor(kGreen + 2);
	SpectraStatZDC[7]->SetLineColor(kAzure + 8);
	SpectraStatZDC[8]->SetLineColor(kBlue - 4);
	//
	SpectraStatZDC[0]->SetMarkerColor(kRed + 1);
	SpectraStatZDC[1]->SetMarkerColor(kRed - 4);
	SpectraStatZDC[2]->SetMarkerColor(kOrange + 7);
	SpectraStatZDC[3]->SetMarkerColor(kOrange - 3);
	SpectraStatZDC[4]->SetMarkerColor(kYellow + 1);
	SpectraStatZDC[5]->SetMarkerColor(kSpring - 7);
	SpectraStatZDC[6]->SetMarkerColor(kGreen + 2);
	SpectraStatZDC[7]->SetMarkerColor(kAzure + 8);
	SpectraStatZDC[8]->SetMarkerColor(kBlue - 4);
	//
    SpectraSystZDC[0]->SetLineColor(kRed + 1);
	SpectraSystZDC[1]->SetLineColor(kRed - 4);
	SpectraSystZDC[2]->SetLineColor(kOrange + 7);
	SpectraSystZDC[3]->SetLineColor(kOrange - 3);
	SpectraSystZDC[4]->SetLineColor(kYellow + 1);
	SpectraSystZDC[5]->SetLineColor(kSpring - 7);
	SpectraSystZDC[6]->SetLineColor(kGreen + 2);
	SpectraSystZDC[7]->SetLineColor(kAzure + 8);
	SpectraSystZDC[8]->SetLineColor(kBlue - 4);
    //
	SpectraSystZDC[0]->SetMarkerColor(kRed + 1);
	SpectraSystZDC[1]->SetMarkerColor(kRed - 4);
	SpectraSystZDC[2]->SetMarkerColor(kOrange + 7);
	SpectraSystZDC[3]->SetMarkerColor(kOrange - 3);
	SpectraSystZDC[4]->SetMarkerColor(kYellow + 1);
	SpectraSystZDC[5]->SetMarkerColor(kSpring - 7);
	SpectraSystZDC[6]->SetMarkerColor(kGreen + 2);
	SpectraSystZDC[7]->SetMarkerColor(kAzure + 8);
	SpectraSystZDC[8]->SetMarkerColor(kBlue - 4);   

   /* SpectraStatZDC[0]->SetMarkerStyle(20+0);
    SpectraSystZDC[0]->SetMarkerStyle(20+0);
    SpectraStatZDC[1]->SetMarkerStyle(20+1);
    SpectraSystZDC[1]->SetMarkerStyle(20+1);
    SpectraStatZDC[2]->SetMarkerStyle(20+2);
    SpectraSystZDC[2]->SetMarkerStyle(20+2);
    SpectraStatZDC[3]->SetMarkerStyle(20+3);
    SpectraSystZDC[3]->SetMarkerStyle(20+3);
    SpectraStatZDC[4]->SetMarkerStyle(20+9);
    SpectraSystZDC[4]->SetMarkerStyle(20+9);
    SpectraStatZDC[5]->SetMarkerStyle(20+13);
    SpectraSystZDC[5]->SetMarkerStyle(20+13);
    SpectraStatZDC[6]->SetMarkerStyle(20+14);
    SpectraSystZDC[6]->SetMarkerStyle(20+14);
    SpectraStatZDC[7]->SetMarkerStyle(47);
    SpectraSystZDC[7]->SetMarkerStyle(47);
    SpectraStatZDC[8]->SetMarkerStyle(39);
    SpectraSystZDC[8]->SetMarkerStyle(39);
    SpectraStatZDC[8]->SetMarkerSize(1.5);
    SpectraSystZDC[8]->SetMarkerSize(1.5);
    SpectraStatZDC[8]->SetMarkerSize(1.5);
    SpectraSystZDC[4]->SetMarkerSize(1.5);
    SpectraStatZDC[4]->SetMarkerSize(1.5);
    SpectraSystZDC[5]->SetMarkerSize(1.5);
    SpectraStatZDC[5]->SetMarkerSize(1.5);*/

    // Ratio to MB
    TH1D* RatioSpectraZDC[nbinZDC], * RatioSpectraZDCSyst[nbinZDC];
    TH1D* SpectraZDCMB = (TH1D*)SpectraStatZDC[0]->Clone("SpectraZDCMB");
    SpectraZDCMB->Reset();
    cout << "done" << endl;
    // Fill histos
    for (int bin = 1; bin <= SpectraZDCMB->GetNbinsX(); bin ++){ //loop over histogram bins
        // volatile variables
        double contsum = 0, staterrsum = 0;
        cout << "done" << endl;cout << "done" << endl;

        for (int i = 0; i < nbinZDC; i++){ //loop over classes
            contsum += SpectraStatZDC[i]->GetBinContent(bin) * (percentileZDC[i+1]-percentileZDC[i]) / 100.;
            staterrsum += (SpectraStatZDC[i]->GetBinError(bin)*(percentileZDC[i+1]-percentileZDC[i])/100.
                              *SpectraStatZDC[i]->GetBinError(bin)*(percentileZDC[i+1]-percentileZDC[i])/100.);
        }    

        SpectraZDCMB->SetBinContent(bin,contsum);
        SpectraZDCMB->SetBinError(bin,TMath::Sqrt(staterrsum));
        cout << "done" << endl;    
    } 


    // V0M 
    TH1D* SpectraStatV0M[nbinV0],*SpectraSystV0M[nbinV0];	
    TF1* levyV0M[nbinV0],*levyFior[nbinV0];
	TFile* FiorfileMB = new TFile("../CLEAN/SpectraVsMultiplicityXi.root","READ");
    TFile* FiorfileLT = new TFile("~/Scaricati/LTspectraFitsVsMultiplicityXi.root","READ");
	
    for(int nmult = 0; nmult < nbinV0; nmult++){
        // levy
        
        levyFior[nmult] = (TF1*)FiorfileLT->Get(Form("LTfitSpectraVsV0M%03.0f00to%03.0f00",percentileV0[nmult],percentileV0[nmult+1]));
        
        levyV0M[nmult] = LevyTsallis(Form("Levy%s_%.0f-%.0f","Xi",percentileZDC[nmult],percentileZDC[nmult+1]), 1.32131 ,1e-3,1e+3,1e-6,1e+6);
        levyV0M[nmult]->SetParameter("norm",levyFior[nmult]->GetParameter(0));
        levyV0M[nmult]->SetParameter("C",levyFior[nmult]->GetParameter(1));
        levyV0M[nmult]->SetParameter("n",levyFior[nmult]->GetParameter(2));
       
        // ---- stats -----
        SpectraStatV0M[nmult] = (TH1D*)FiorfileMB->Get(Form("hPtXiStatOnly_V0M_%03.0f00to%03.0f00-epsPart-epsEv-Corrected",percentileV0[nmult],percentileV0[nmult+1]));
        SpectraSystV0M[nmult] = (TH1D*)FiorfileMB->Get(Form("hPtXiSystUnco_V0M_%03.0f00to%03.0f00-epsPart-epsEv-Corrected",percentileV0[nmult],percentileV0[nmult+1]));
    }

    //
   
    SpectraStatV0M[0]->SetLineColor(kRed + 1);
	SpectraStatV0M[1]->SetLineColor(kRed - 4);
	SpectraStatV0M[2]->SetLineColor(kOrange + 7);
	SpectraStatV0M[3]->SetLineColor(kOrange - 3);
	SpectraStatV0M[4]->SetLineColor(kYellow + 1);
	SpectraStatV0M[5]->SetLineColor(kSpring - 7);
	SpectraStatV0M[6]->SetLineColor(kGreen + 2);
	SpectraStatV0M[7]->SetLineColor(kAzure + 8);
	SpectraStatV0M[8]->SetLineColor(kBlue - 4);
	SpectraStatV0M[9]->SetLineColor(kBlue + 3);

	SpectraStatV0M[0]->SetMarkerColor(kRed + 1);
	SpectraStatV0M[1]->SetMarkerColor(kRed - 4);
	SpectraStatV0M[2]->SetMarkerColor(kOrange + 7);
	SpectraStatV0M[3]->SetMarkerColor(kOrange - 3);
	SpectraStatV0M[4]->SetMarkerColor(kYellow + 1);
	SpectraStatV0M[5]->SetMarkerColor(kSpring - 7);
	SpectraStatV0M[6]->SetMarkerColor(kGreen + 2);
	SpectraStatV0M[7]->SetMarkerColor(kAzure + 8);
	SpectraStatV0M[8]->SetMarkerColor(kBlue - 4);
	SpectraStatV0M[9]->SetMarkerColor(kBlue + 3);



    // Syst

    SpectraSystV0M[0]->SetLineColor(kRed + 1);
	SpectraSystV0M[1]->SetLineColor(kRed - 4);
	SpectraSystV0M[2]->SetLineColor(kOrange + 7);
	SpectraSystV0M[3]->SetLineColor(kOrange - 3);
	SpectraSystV0M[4]->SetLineColor(kYellow + 1);
	SpectraSystV0M[5]->SetLineColor(kSpring - 7);
	SpectraSystV0M[6]->SetLineColor(kGreen + 2);
	SpectraSystV0M[7]->SetLineColor(kAzure + 8);
	SpectraSystV0M[8]->SetLineColor(kBlue - 4);
	SpectraSystV0M[9]->SetLineColor(kBlue + 3);

	SpectraSystV0M[0]->SetMarkerColor(kRed + 1);
	SpectraSystV0M[1]->SetMarkerColor(kRed - 4);
	SpectraSystV0M[2]->SetMarkerColor(kOrange + 7);
	SpectraSystV0M[3]->SetMarkerColor(kOrange - 3);
	SpectraSystV0M[4]->SetMarkerColor(kYellow + 1);
	SpectraSystV0M[5]->SetMarkerColor(kSpring - 7);
	SpectraSystV0M[6]->SetMarkerColor(kGreen + 2);
	SpectraSystV0M[7]->SetMarkerColor(kAzure + 8);
	SpectraSystV0M[8]->SetMarkerColor(kBlue - 4);
	SpectraSystV0M[9]->SetMarkerColor(kBlue + 3);

   	TH1D* LevyFitSpectraV0M[nbinV0];
   	for(int nmult = 0; nmult < nbinV0; nmult++){
       	LevyFitSpectraV0M[nmult] = (TH1D*)SpectraStatV0M[nmult]->Clone(Form("LevyFitSpectra-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            percentileV0[nmult],percentileV0[nmult+1],0.,100.  ));
         LevyFitSpectraV0M[nmult]->SetStats(0);

        //	
		levyV0M[nmult]->SetLineStyle(7);
        levyV0M[nmult]->SetLineColor(kBlack);
        levyV0M[nmult]->SetLineWidth(1);
        SpectraStatV0M[nmult]->SetMarkerStyle(20);
        SpectraSystV0M[nmult]->SetMarkerStyle(20);
        SpectraStatV0M[nmult]->SetMarkerSize(2);
        SpectraSystV0M[nmult]->SetMarkerSize(2);
       

   	}

   /* SpectraStatV0M[0]->SetMarkerStyle(20+0);
    SpectraSystV0M[0]->SetMarkerStyle(20+0);
    SpectraStatV0M[1]->SetMarkerStyle(20+1);
    SpectraSystV0M[1]->SetMarkerStyle(20+1);
     SpectraStatV0M[2]->SetMarkerStyle(20+2);
    SpectraSystV0M[2]->SetMarkerStyle(20+2);
     SpectraStatV0M[3]->SetMarkerStyle(20+3);
    SpectraSystV0M[3]->SetMarkerStyle(20+3);
     SpectraStatV0M[4]->SetMarkerStyle(20+9);
    SpectraSystV0M[4]->SetMarkerStyle(20+9);
     SpectraStatV0M[5]->SetMarkerStyle(20+13);
    SpectraSystV0M[5]->SetMarkerStyle(20+13);
     SpectraStatV0M[6]->SetMarkerStyle(20+14);
    SpectraSystV0M[6]->SetMarkerStyle(20+14);
     SpectraStatV0M[7]->SetMarkerStyle(47);
    SpectraSystV0M[7]->SetMarkerStyle(47);
     SpectraStatV0M[8]->SetMarkerStyle(39);
    SpectraSystV0M[8]->SetMarkerStyle(39);
     SpectraStatV0M[8]->SetMarkerSize(2.6);
        SpectraSystV0M[8]->SetMarkerSize(2.6);
        SpectraStatV0M[8]->SetMarkerSize(2.6);
        SpectraSystV0M[4]->SetMarkerSize(2.6);
        SpectraStatV0M[4]->SetMarkerSize(2.6);
        SpectraSystV0M[5]->SetMarkerSize(2.6);
        SpectraStatV0M[5]->SetMarkerSize(2.6);
    SpectraStatV0M[9]->SetMarkerStyle(41);
    SpectraSystV0M[9]->SetMarkerStyle(41);*/

     // Ratio to MB
    TH1D* RatioSpectraV0M[nbinV0],* RatioSpectraV0MSyst[nbinV0];
    TH1D* SpectraV0MMB = (TH1D*)FiorfileMB->Get(Form("hPtXiStatOnly_V0M_%03.0f00to%03.0f00-epsPart-epsEv-Corrected",0.,100.));
     TH1D* SpectraV0MMBSyst = (TH1D*)FiorfileMB->Get(Form("hPtXiSystOnly_V0M_%03.0f00to%03.0f00-epsPart-epsEv-Corrected",0.,100.));
   
   	for(int nmult = 0; nmult < nbinV0; nmult++){
        RatioSpectraV0M[nmult] = (TH1D*)SpectraStatV0M[nmult]->Clone(Form("RatioSpectraV0M-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            0.,100.,percentileV0[nmult],percentileV0[nmult+1] ));
        RatioSpectraV0M[nmult]->Reset();
        RatioSpectraV0MSyst[nmult] = (TH1D*)SpectraSystV0M[nmult]->Clone(Form("RatioSpectraSystV0M-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            0.,100.,percentileV0[nmult],percentileV0[nmult+1] ));
        RatioSpectraV0MSyst[nmult]->Reset();
        cout << "done" << endl;
        //
        for (int bin = 1; bin <= RatioSpectraV0M[0]->GetNbinsX(); bin ++){ //loop over histogram bins
            RatioSpectraV0M[nmult]->SetBinContent(bin,SpectraStatV0M[nmult]->GetBinContent(bin)/SpectraV0MMB->GetBinContent(bin));
            RatioSpectraV0M[nmult]->SetBinError(bin,ErrorInRatio(SpectraStatV0M[nmult]->GetBinContent(bin),SpectraStatV0M[nmult]->GetBinError(bin),SpectraV0MMB->GetBinContent(bin),SpectraV0MMB->GetBinError(bin)));
            RatioSpectraV0MSyst[nmult]->SetBinContent(bin,SpectraStatV0M[nmult]->GetBinContent(bin)/SpectraV0MMB->GetBinContent(bin));
            RatioSpectraV0MSyst[nmult]->SetBinError(bin,SpectraSystV0M[nmult]->GetBinError(bin)/SpectraV0MMB->GetBinContent(bin));
            
        }
    }

       	for(int nmult = 0; nmult < nbinZDC; nmult++){
        RatioSpectraZDC[nmult] = (TH1D*)SpectraStatZDC[nmult]->Clone(Form("RatioSpectra-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            0.,100.,percentileZDC[nmult],percentileZDC[nmult+1] ));
        RatioSpectraZDC[nmult]->Reset();
        RatioSpectraZDCSyst[nmult] = (TH1D*)SpectraSystZDC[nmult]->Clone(Form("RatioSpectraSyst-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            0.,100.,percentileZDC[nmult],percentileZDC[nmult+1] ));
        RatioSpectraZDCSyst[nmult]->Reset();
        cout << "done" << endl;
        //
        for (int bin = 1; bin <= RatioSpectraZDC[0]->GetNbinsX(); bin ++){ //loop over histogram bins
            RatioSpectraZDC[nmult]->SetBinContent(bin,SpectraStatZDC[nmult]->GetBinContent(bin)/SpectraV0MMB->GetBinContent(bin));
            RatioSpectraZDC[nmult]->SetBinError(bin,ErrorInRatio(SpectraStatZDC[nmult]->GetBinContent(bin),SpectraStatZDC[nmult]->GetBinError(bin),SpectraV0MMB->GetBinContent(bin),SpectraV0MMB->GetBinError(bin)));
            RatioSpectraZDCSyst[nmult]->SetBinContent(bin,SpectraStatZDC[nmult]->GetBinContent(bin)/SpectraV0MMB->GetBinContent(bin));
            RatioSpectraZDCSyst[nmult]->SetBinError(bin,
            //ErrorInRatio(SpectraSystZDC[nmult]->GetBinContent(bin),SpectraSystZDC[nmult]->GetBinError(bin),SpectraV0MMBSyst->GetBinContent(bin),SpectraV0MMBSyst->GetBinError(bin)));
            SpectraSystZDC[nmult]->GetBinError(bin)/SpectraV0MMBSyst->GetBinContent(bin));
            cout<< "new " << SpectraSystZDC[nmult]->GetBinError(bin)/SpectraV0MMBSyst->GetBinContent(bin) << endl;
            cout << "conservativo  " << ErrorInRatio(SpectraSystZDC[nmult]->GetBinContent(bin),SpectraSystZDC[nmult]->GetBinError(bin),SpectraV0MMBSyst->GetBinContent(bin),SpectraV0MMBSyst->GetBinError(bin)) << endl;
            // ErrorInRatio(SpectraSystZDC[nmult]->GetBinContent(bin),SpectraSystZDC[nmult]->GetBinError(bin),SpectraV0MMB->GetBinContent(bin),SpectraV0MMB->GetBinError(bin)));
        }
    }
    
    //
  /*  SpectraStatZDC[0]->Scale(64);
	SpectraStatZDC[1]->Scale(32);
	SpectraStatZDC[2]->Scale(16);
	SpectraStatZDC[3]->Scale(8);
	SpectraStatZDC[4]->Scale(4);
	SpectraStatZDC[5]->Scale(2);
	SpectraStatZDC[6]->Scale(1);
	SpectraStatZDC[7]->Scale(0.5);
	SpectraStatZDC[8]->Scale(0.25);
	//

    //
    SpectraSystZDC[0]->Scale(64);
	SpectraSystZDC[1]->Scale(32);
	SpectraSystZDC[2]->Scale(16);
	SpectraSystZDC[3]->Scale(8);
	SpectraSystZDC[4]->Scale(4);
	SpectraSystZDC[5]->Scale(2);
	SpectraSystZDC[6]->Scale(1);
	SpectraSystZDC[7]->Scale(0.5);
	SpectraSystZDC[8]->Scale(0.25);*/

    /*SpectraSystV0M[0]->Scale(64);
	SpectraSystV0M[1]->Scale(32);
	SpectraSystV0M[2]->Scale(16);
	SpectraSystV0M[3]->Scale(8);
	SpectraSystV0M[4]->Scale(4);
	SpectraSystV0M[5]->Scale(2);
	SpectraSystV0M[6]->Scale(1);
	SpectraSystV0M[7]->Scale(0.5);
	SpectraSystV0M[8]->Scale(0.25);
	SpectraSystV0M[9]->Scale(0.10);
    //
    SpectraStatV0M[0]->Scale(64);
	SpectraStatV0M[1]->Scale(32);
	SpectraStatV0M[2]->Scale(16);
	SpectraStatV0M[3]->Scale(8);
	SpectraStatV0M[4]->Scale(4);
	SpectraStatV0M[5]->Scale(2);
	SpectraStatV0M[6]->Scale(1);
	SpectraStatV0M[7]->Scale(0.5);
	SpectraStatV0M[8]->Scale(0.25);
	SpectraStatV0M[9]->Scale(0.10);*/

    
   

    //
 	Double_t n[nbinV0];
    /*n[0]=64.;
    n[1]=32.;
    n[2]=16.;
    n[3]=8.;
    n[4]=4.;
    n[5]=2.;
    n[6]=1.;
    n[7]=0.5;
    n[8]=0.25;   */
    for(int nmult = 0; nmult < nbinZDC; nmult++){
        n[nmult] = TMath::Power(2,8-nmult);
        double par = levy[nmult]->GetParameter("norm");
        levy[nmult]->SetParameter("norm",par*(n[nmult]));
        SpectraSystZDC[nmult]->Scale(n[nmult]);
        SpectraStatZDC[nmult]->Scale(n[nmult]);
    }
    Double_t nv0[nbinV0];
    for(int nmult = 0; nmult < nbinV0; nmult++){
        double parV0 = levyV0M[nmult]->GetParameter("norm");
       // cout << parV0 << endl;
        nv0[nmult] = TMath::Power(2,9-nmult);
        levyV0M[nmult]->SetParameter("norm",parV0*(nv0[nmult]));
        SpectraSystV0M[nmult]->Scale(nv0[nmult]);
        SpectraStatV0M[nmult]->Scale(nv0[nmult]);
        levyV0M[nmult]->SetRange(0.,7.);
       // cout << levyV0M[nmult]->GetParameter(0) << endl;
    }
  
    //
  
	TCanvas* c = new TCanvas("c","",2300,2000);

    c->Divide(2);
	c->SetRightMargin(0.09);
    c->SetLeftMargin(0.25);
    c->SetBottomMargin(0.15);
    c->SetFillColor(0);
    c->SetTickx(); 
    c->SetTicky();
	c->cd(1);
    TPad *pad1 = new TPad("pad1","pad1",0,0.37,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.37);
    pad1->SetBottomMargin(0.0001);
    pad1->SetBorderMode(0);
    pad1->SetTopMargin(0.05);
    pad2->SetTopMargin(0.0001);
    pad2->SetBottomMargin(0.3);
    pad1->SetRightMargin(0.01);
    pad2->SetRightMargin(0.01);
    pad1->SetLeftMargin(0.15);
    pad2->SetLeftMargin(0.15);
    pad1->SetBottomMargin(0.01);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();
    //
    pad1->cd();
    pad1->SetLogy();    
    pad1->SetFillColor(kWhite);
  
    
    TLegend* l = new TLegend(0.18,0.03,0.35,0.14);
    l->SetTextSize(0.03);
    l->SetTextFont(42);
   // l->SetHeader(Form("%s classes (percentiles):","(#sqrt{#it{s}} - ZDC)"));
    l->SetBorderSize(0);     
    TLegend* lw = new TLegend(0.42,0.03,0.63,0.14);
    lw->SetTextSize(0.03);
    lw->SetTextFont(42);
   // l->SetHeader(Form("%s classes (percentiles):","(#sqrt{#it{s}} - ZDC)"));
    lw->SetBorderSize(0);
    TLegend* lw1 = new TLegend(0.68,0.03,0.85,0.14);
    lw1->SetTextSize(0.03);
    lw1->SetTextFont(42);
   // l->SetHeader(Form("%s classes (percentiles):","(#sqrt{#it{s}} - ZDC)"));
    lw1->SetBorderSize(0);

    TLegend* lTS = new TLegend(0.20,0.2,0.5,0.310);
    lTS->SetTextSize(0.032);
  //  lTS->SetTextFont(42);
   // l->SetHeader(Form("%s classes (percentiles):","(#sqrt{#it{s}} - ZDC)"));
    lTS->SetBorderSize(0);   

    TH1D* h1 = new TH1D("h1","",100,0.,8.01);
    h1->GetYaxis()->SetRangeUser(1.3E-7,90);
    h1->GetYaxis()->SetLabelSize(0.04);
    h1->GetXaxis()->SetRangeUser(0.,8.);
    h1->GetYaxis()->SetTitle("1/#it{N}_{ev}  d^{2}#it{N}/(d#it{p}_{T}d#it{y}) (GeV/#it{c})^{-1}");
    h1->GetXaxis()->SetTitle("#it{p}_{T}");
    h1->GetYaxis()->SetTitleOffset(1.44);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->SetStats(0);
    h1->Draw();
    SpectraStatZDC[0]->SetTitle("");
	SpectraStatZDC[0]->Draw("SAME");
    levy[0]->SetRange(0.,7.);
    levy[0]->Draw("SAME");
    SpectraSystZDC[0]->SetFillStyle(0.);
	SpectraSystZDC[0]->Draw("SAME E2");
    l->AddEntry(SpectraStatZDC[0],"0-20 % ( #times 2^{8} )","LEP");	
	for(int nmult = 1; nmult < nbinZDC; nmult++){ 
		SpectraSystZDC[nmult]->SetFillStyle(0.);
		SpectraSystZDC[nmult]->Draw("SAME E2");
		SpectraStatZDC[nmult]->SetTitle(Form("%.0f-%.0f % ( #times 2^{%i} )",percentileZDC[nmult],percentileZDC[nmult+1],8-nmult));
	  	SpectraStatZDC[nmult]->Draw("SAME");
        levy[nmult]->SetRange(0.,7.);
        levy[nmult]->Draw("SAME");		
	}
    for(int nmult = 1; nmult < 3; nmult++){
        l->AddEntry(SpectraStatZDC[nmult],"","LEP");
    }
    for(int nmult = 3; nmult < 6; nmult++){
        lw->AddEntry(SpectraStatZDC[nmult],"","LEP");
    }
    for(int nmult = 6; nmult < nbinZDC; nmult++){
        lw1->AddEntry(SpectraStatZDC[nmult],"","LEP");
    }
	lTS->AddEntry(levy[nbinZDC-1],"Levy-Tsallis fit","L");
	l->Draw("SAME");  
    lw->Draw("SAME");
    lw1->Draw("SAME");
    lTS->Draw("SAME");
    // draw text
    TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
    xlabel-> SetNDC();
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.08);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    xlabel-> DrawLatex(0.82, 0.77, Form("#Xi^{-} + #bar{#Xi}^{+}"));
    TLatex *xlabe = new TLatex(1.,1.E-4,"|#it{y}| < 0.5");
	//xlabe->SetTextFont(42);
    xlabe-> SetTextSize(0.035);
    xlabe->Draw("SAME");
    TLatex *   texbr = new TLatex(4.7,22,"ALICE Preliminary");
    TLatex *   texbrpp = new TLatex(5.2,9,"pp #sqrt{#it{s}} = 13 TeV");
    TLatex *   texclasses = new TLatex(0.7,.65E-5,"(#sqrt{#it{s}} - ZDC) effective energy classes:");

    //, pp #sqrt{#it{s}} = 13 TeV");
    //texbr->SetTextFont(42);
    texbr->SetTextSize(0.03208798);
    texbrpp->SetTextFont(42);
    texbrpp->SetTextSize(0.03208798);
   // texb->SetLineWidth(2);
    texbr->Draw("SAME");
    texbrpp->Draw("SAME");
    texclasses->SetTextSize(0.032);
   // texb->SetLineWidth(2);
    texclasses->Draw("SAME");

      

    pad2->cd();
    pad2->SetLogy();    
    pad2->SetFillColor(kWhite);
    pad2->SetLogy();
    TH1D* h2 = new TH1D("h2","",100,0.,8.01);
    h2->Draw();

    h2->GetYaxis()->SetRangeUser(0.035,15);
    h2->GetXaxis()->SetRangeUser(0.00,8.);
    h2->GetYaxis()->SetTitle("Ratio w.r.t. INEL>0");
    h2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h2->GetYaxis()->SetTitleOffset(.7);
    h2->GetYaxis()->SetTitleSize(0.09);
    h2->GetYaxis()->SetLabelSize(0.09);
    h2->GetXaxis()->SetTitleOffset(.93);
    h2->GetXaxis()->SetTitleSize(0.1);
    h2->GetXaxis()->SetLabelSize(0.08);
    h2->SetTitle("");
    h2->SetStats(0);
    RatioSpectraZDC[0]->SetLineWidth(2);
    RatioSpectraZDCSyst[0]->SetFillStyle(0.);
	RatioSpectraZDC[0]->Draw("SAME");
    RatioSpectraZDCSyst[0]->Draw("SAME E2");
    RatioSpectraZDC[0]->SetFillStyle(0.);
	for(int nmult = 1; nmult < nbinZDC; nmult++){ 
        RatioSpectraZDCSyst[nmult]->SetFillStyle(0.);
        RatioSpectraZDCSyst[nmult]->Draw("SAME E2");
        RatioSpectraZDC[nmult]->SetStats(0);
		RatioSpectraZDC[nmult]->Draw("SAME");
    }


    c->cd(2);
    TPad *pad1v0 = new TPad("pad1v0","pad1",0,0.37,1,1);
    TPad *pad2v0 = new TPad("pad2v0","pad2",0,0,1,0.37);
    pad1v0->SetBottomMargin(0.0001);
    pad1v0->SetBorderMode(0);
    pad2v0->SetTopMargin(0.00001);
    pad1v0->SetTopMargin(0.05);
    pad2v0->SetBottomMargin(0.3);
    pad1v0->SetRightMargin(0.095);
    pad2v0->SetRightMargin(0.095);
    pad1v0->SetLeftMargin(0.008);
    pad2v0->SetLeftMargin(0.008);
    pad1v0->SetBottomMargin(0.01);
    pad2v0->SetBorderMode(0);
    pad1v0->Draw();
    pad2v0->Draw();
    //
    pad1v0->cd();
    pad1v0->SetLogy();    
    pad1v0->SetFillColor(kWhite);
    
    TLegend* lv0 = new TLegend(0.06,0.03,0.2,0.18);
    lv0->SetTextSize(0.027);
    lv0->SetTextFont(42);
   // l->SetHeader(Form("%s classes (percentiles):","(#sqrt{#it{s}} - ZDC)"));
    lv0->SetBorderSize(0);     
    TLegend* lwv0 = new TLegend(0.25,0.03,0.4,0.18);
    lwv0->SetTextSize(0.027);
    lwv0->SetTextFont(42);

    TLegend* lTSv0 = new TLegend(0.05,0.2,0.35,0.310);
    lTSv0->SetTextSize(0.027);
  //  lTS->SetTextFont(42);
   // l->SetHeader(Form("%s classes (percentiles):","(#sqrt{#it{s}} - ZDC)"));
    lTSv0->SetBorderSize(0); 
    h1->Draw();
 
    SpectraStatV0M[0]->GetYaxis()->SetRangeUser(1.1E-7,90);
    SpectraStatV0M[0]->GetXaxis()->SetRangeUser(-0.001,8.);
    SpectraStatV0M[0]->GetYaxis()->SetTitle("1/N_{ev} d^{2}N/(dp_{T}dy) [(GeV/c)^{-1}]");
    SpectraStatV0M[0]->GetYaxis()->SetTitleOffset(1.75);
    SpectraStatV0M[0]->SetTitle("");
    SpectraStatV0M[0]->SetStats(0);
	SpectraStatV0M[0]->Draw("SAME");
    levyV0M[0]->Draw("SAME");
    SpectraSystV0M[0]->SetFillStyle(0.);
	SpectraSystV0M[0]->Draw("SAME E2");
    lv0->AddEntry(SpectraStatV0M[0],"0-1 ( #times 2^{9} )","LEP");	
	for(int nmult = 1; nmult < nbinV0; nmult++){ 
		SpectraSystV0M[nmult]->SetFillStyle(0.);
		SpectraSystV0M[nmult]->Draw("SAME E2");
		SpectraStatV0M[nmult]->SetTitle(Form("%.0f-%.0f ( #times 2^{%i} )",percentileV0[nmult],percentileV0[nmult+1],9-nmult));
	  	SpectraStatV0M[nmult]->Draw("SAME");
        SpectraStatV0M[nmult]->SetStats(0);
        levyV0M[nmult]->Draw("SAME");		
	}
    for(int nmult = 1; nmult < 5; nmult++){
	//lv0->AddEntry(levyV0M[nbinV0-1],"Levy-Tsallis","L");
        lv0->AddEntry(SpectraStatV0M[nmult],"","LEP");
    }
    for(int nmult = 5; nmult < nbinV0; nmult++){
	//lv0->AddEntry(levyV0M[nbinV0-1],"Levy-Tsallis","L");
        lwv0->AddEntry(SpectraStatV0M[nmult],"","LEP");
    }
    lTSv0->AddEntry(levyV0M[nbinV0-1],"Levy-Tsallis fit","L");
	lv0->Draw("SAME");  
    lwv0->Draw("SAME");
    lTSv0->Draw("SAME");
    // draw text
    xlabel-> DrawLatex(0.75, 0.77, Form("#Xi^{-} + #bar{#Xi}^{+}"));
   // xlaber-> DrawLatex(0.66, 0.18, Form("EPJC80167(2020)"));
    TLatex *   texb = new TLatex(4.3,22,"ALICE EPJC80167(2020)");
   // texb->SetTextFont(42);
    texb->SetTextSize(0.0320798);
   // texb->SetLineWidth(2);
    texb->Draw("SAME");
    TLatex *   texbrppv0 = new TLatex(5.5,9,"pp #sqrt{#it{s}} = 13 TeV");
    texbrppv0->SetTextFont(42);
    texbrppv0->SetTextSize(0.03208798);
    texbrppv0->Draw("SAME");
    
    TLatex *   texclassesv0 = new TLatex(0.6,.65E-5,"V0M multiplicity classes:");
    texclassesv0->SetTextSize(0.025);
   // texb->SetLineWidth(2);
    texclassesv0->Draw("SAME");

   
    TLatex *xlabev0 = new TLatex(1,1E-4,"|#it{y}| < 0.5");
	//xlabe->SetTextFont(42);
    xlabev0-> SetTextSize(0.035);
    xlabev0->Draw("SAME");
   
  
    

    pad2v0->cd();
    pad2v0->SetLogy();    
    pad2v0->SetFillColor(kWhite);
    pad2v0->SetLogy();

    h2->Draw();
    RatioSpectraV0M[0]->GetYaxis()->SetRangeUser(0.035,15);
    RatioSpectraV0M[0]->GetXaxis()->SetRangeUser(0.00,8.);
    RatioSpectraV0M[0]->GetYaxis()->SetTitle("Ratio w.r.t. [0-100] %");
    RatioSpectraV0M[0]->GetYaxis()->SetTitleOffset(.7);
    RatioSpectraV0M[0]->GetYaxis()->SetTitleSize(0.08);
    RatioSpectraV0M[0]->GetYaxis()->SetLabelSize(0.08);
    RatioSpectraV0M[0]->GetXaxis()->SetTitleOffset(1.1);
    RatioSpectraV0M[0]->GetXaxis()->SetTitleSize(0.1);
    RatioSpectraV0M[0]->GetXaxis()->SetLabelSize(0.08);
    RatioSpectraV0M[0]->SetTitle("");
    RatioSpectraV0M[0]->SetStats(0);
	RatioSpectraV0M[0]->Draw("SAME");
    RatioSpectraV0MSyst[0]->Draw("SAME E2");
    RatioSpectraV0MSyst[0]->SetFillStyle(0.);
    RatioSpectraV0M[0]->Draw("SAME");
    RatioSpectraV0M[0]->SetFillStyle(0.);
	for(int nmult = 1; nmult < nbinV0; nmult++){ 
        RatioSpectraV0MSyst[nmult]->SetFillStyle(0.);
		RatioSpectraV0M[nmult]->Draw("SAME");
        RatioSpectraV0MSyst[nmult]->Draw("E2 SAME");
    }

    new TCanvas;
    RatioSpectraV0MSyst[0]->SetFillStyle(0.);
    RatioSpectraV0MSyst[0]->SetLineColor(kRed);
    RatioSpectraV0MSyst[0]->Draw("E2");
    RatioSpectraV0MSyst[0]->SetLineColor(kRed);    

      
    c->SaveAs("immaginifinali/XiCorrectedpTSpectra_ZDCV0M_LevyTsallis.eps");
   

}


//---------------------------------------------------------------------------------------------------
double ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr) {	
  // Error in a Ratio
  if (B != 0) {
    Double_t errorfromtop = Aerr * Aerr / (B * B);
    Double_t errorfrombottom = ((A * A) / (B * B * B * B)) * Berr * Berr;
    return TMath::Sqrt(TMath::Abs(errorfromtop + errorfrombottom));
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
