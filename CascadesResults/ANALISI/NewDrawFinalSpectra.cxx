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
double ErrorInRatioBarlow(Double_t A, Double_t Aerr, Double_t B, Double_t Berr);
void DivideAndComputeRogerBarlow(TH1F *h1, TH1F *h2);
TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t par2min = 1e-3, Double_t par2max = 1e+3,Double_t par3min = 1e-6, Double_t par3max = 1e+6, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);

void NewDrawFinalSpectra(){

    //percentile
    Double_t percentileV0[] = {0.,1.,5,10,15,20,30,40,50,70,100};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    Double_t percentileZDC[] = {0,20,30,40,50,60,70,80,90,100};
    const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;
    
    //ZDC ================================================================================================================
    TString namesystsgnloss = Form("%s-%s","EEsel","ZDC");
    TFile* fileSgnLossSyst = TFile::Open("NormalizationCorrections/SgnLossSyst.root");
    TH1D* hsgnloss[nbinZDC];
    
	TH1D* SpectraStatZDC[nbinZDC],*hSystPartZDC,*SpectraSystZDC[nbinZDC];
	TFile* fileStatZDC = new TFile(Form("StatSpectra%s-Xi-%s_%03.0f_%03.0f.root","ZDC","V0M",0.,100.));
    if (!fileStatZDC) cout << "NO FILE" << endl;
    TFile* fileSystZDC = TFile::Open(Form("sistematiche/SystematicsFinalResults-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",
            "Xi",0.,100.,0.,100.));;
    TF1* levy[nbinZDC];

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
        SpectraStatZDC[nmult]->SetMarkerSize(2.2);
        SpectraSystZDC[nmult]->SetMarkerSize(2.2);
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
	SpectraStatZDC[5]->SetLineColor(kGreen+1);
	SpectraStatZDC[6]->SetLineColor(kTeal - 5);
	SpectraStatZDC[7]->SetLineColor(kAzure + 7);
	SpectraStatZDC[8]->SetLineColor(kBlue - 4);
	//
	SpectraStatZDC[0]->SetMarkerColor(kRed + 1);
	SpectraStatZDC[1]->SetMarkerColor(kRed - 4);
	SpectraStatZDC[2]->SetMarkerColor(kOrange + 7);
	SpectraStatZDC[3]->SetMarkerColor(kOrange - 3);
	SpectraStatZDC[4]->SetMarkerColor(kYellow + 1);
	SpectraStatZDC[5]->SetMarkerColor(kGreen+1);
	SpectraStatZDC[6]->SetMarkerColor(kTeal - 5);
	SpectraStatZDC[7]->SetMarkerColor(kAzure + 7);
	SpectraStatZDC[8]->SetMarkerColor(kBlue - 4);
	//
    SpectraSystZDC[0]->SetLineColor(kRed + 1);
	SpectraSystZDC[1]->SetLineColor(kRed - 4);
	SpectraSystZDC[2]->SetLineColor(kOrange + 7);
	SpectraSystZDC[3]->SetLineColor(kOrange - 3);
	SpectraSystZDC[4]->SetLineColor(kYellow + 1);
	SpectraSystZDC[5]->SetLineColor(kGreen+1);
	SpectraSystZDC[6]->SetLineColor(kTeal - 5);
	SpectraSystZDC[7]->SetLineColor(kAzure + 7);
	SpectraSystZDC[8]->SetLineColor(kBlue - 4);

    //V0 ======================================================================================================

    TH1D* SpectraStatV0M[nbinV0],*SpectraSystTotV0M[nbinV0],*SpectraSystUncoV0M[nbinV0];	
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
        SpectraSystTotV0M[nmult] = (TH1D*)FiorfileMB->Get(Form("hPtXiSystOnly_V0M_%03.0f00to%03.0f00-epsPart-epsEv-Corrected",percentileV0[nmult],percentileV0[nmult+1]));
        SpectraSystUncoV0M[nmult] = (TH1D*)FiorfileMB->Get(Form("hPtXiSystUnco_V0M_%03.0f00to%03.0f00-epsPart-epsEv-Corrected",percentileV0[nmult],percentileV0[nmult+1]));
    
        levyV0M[nmult]->SetLineStyle(7);
        levyV0M[nmult]->SetLineColor(kBlack);
        levyV0M[nmult]->SetLineWidth(1);
        SpectraStatV0M[nmult]->SetMarkerStyle(20);
        SpectraSystTotV0M[nmult]->SetMarkerStyle(20);
        SpectraStatV0M[nmult]->SetMarkerSize(2.2);
        SpectraSystTotV0M[nmult]->SetMarkerSize(2.2);
    }

    // Stat
   
    SpectraStatV0M[0]->SetLineColor(kRed + 1);
	SpectraStatV0M[1]->SetLineColor(kRed - 4);
	SpectraStatV0M[2]->SetLineColor(kOrange + 7);
	SpectraStatV0M[3]->SetLineColor(kOrange - 3);
	SpectraStatV0M[4]->SetLineColor(kYellow + 1);
	SpectraStatV0M[5]->SetLineColor(kGreen+ 1);
	SpectraStatV0M[6]->SetLineColor(kTeal - 5);
	SpectraStatV0M[7]->SetLineColor(kAzure + 8);
	SpectraStatV0M[8]->SetLineColor(kAzure+7 - 4);
	SpectraStatV0M[9]->SetLineColor(kBlue + 3);

	SpectraStatV0M[0]->SetMarkerColor(kRed + 1);
	SpectraStatV0M[1]->SetMarkerColor(kRed - 4);
	SpectraStatV0M[2]->SetMarkerColor(kOrange + 7);
	SpectraStatV0M[3]->SetMarkerColor(kOrange - 3);
	SpectraStatV0M[4]->SetMarkerColor(kYellow + 1);
	SpectraStatV0M[5]->SetMarkerColor(kGreen +1);
	SpectraStatV0M[6]->SetMarkerColor(kTeal - 5);
	SpectraStatV0M[7]->SetMarkerColor(kAzure+7);
	SpectraStatV0M[8]->SetMarkerColor(kBlue - 4);
	SpectraStatV0M[9]->SetMarkerColor(kBlue + 3);

    // Syst

    SpectraSystTotV0M[0]->SetLineColor(kRed + 1);
	SpectraSystTotV0M[1]->SetLineColor(kRed - 4);
	SpectraSystTotV0M[2]->SetLineColor(kOrange + 7);
	SpectraSystTotV0M[3]->SetLineColor(kOrange - 3);
	SpectraSystTotV0M[4]->SetLineColor(kYellow + 1);
	SpectraSystTotV0M[5]->SetLineColor(kGreen +1);
	SpectraSystTotV0M[6]->SetLineColor(kTeal - 5);
	SpectraSystTotV0M[7]->SetLineColor(kAzure + 7);
	SpectraSystTotV0M[8]->SetLineColor(kBlue - 4);
	SpectraSystTotV0M[9]->SetLineColor(kBlue + 3);

    
    // Ratio to MB
    TH1D* RatioSpectraV0M[nbinV0],* RatioSpectraV0MSyst[nbinV0];


    TFile* fhep = TFile::Open("HEPData-ins1748157-v1-Table_3R.root");
    TFile* fhep1 = TFile::Open("HEPData-ins1748157-v1-Table_3.root");
    TH1D* erruncorr[nbinV0]; 
    TH1D* huncorrINEL = (TH1D*)fhep1->Get("Table 3/Hist1D_y11_e3");
    TH1D* SpectraV0MMB = (TH1D*)FiorfileMB->Get(Form("hPtXiStatOnly_V0M_%03.0f00to%03.0f00-epsPart-epsEv-Corrected",0.,100.));
    TH1D* SpectraV0MMBSyst = (TH1D*)FiorfileMB->Get(Form("hPtXiSystOnly_V0M_%03.0f00to%03.0f00-epsPart-epsEv-Corrected",0.,100.));
   
   	for(int nmult = 0; nmult < nbinV0; nmult++){

        erruncorr[nmult] = (TH1D*)fhep->Get(Form("Table 3R/Hist1D_y%i_e2",nmult+1));
        RatioSpectraV0M[nmult] = (TH1D*)SpectraStatV0M[nmult]->Clone(Form("RatioSpectraV0M-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            0.,100.,percentileV0[nmult],percentileV0[nmult+1] ));
        RatioSpectraV0M[nmult]->Reset();
        RatioSpectraV0MSyst[nmult] = (TH1D*)SpectraSystTotV0M[nmult]->Clone(Form("RatioSpectraSystV0M-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            0.,100.,percentileV0[nmult],percentileV0[nmult+1] ));
        RatioSpectraV0MSyst[nmult]->Reset();
        //
        for (int bin = 1; bin <= RatioSpectraV0M[0]->GetNbinsX(); bin ++){ //loop over histogram bins
            RatioSpectraV0M[nmult]->SetBinContent(bin,SpectraStatV0M[nmult]->GetBinContent(bin)/SpectraV0MMB->GetBinContent(bin));
            RatioSpectraV0M[nmult]->SetBinError(bin,ErrorInRatioBarlow(SpectraStatV0M[nmult]->GetBinContent(bin),SpectraStatV0M[nmult]->GetBinError(bin),SpectraV0MMB->GetBinContent(bin),SpectraV0MMB->GetBinError(bin)));
            RatioSpectraV0MSyst[nmult]->SetBinContent(bin,SpectraStatV0M[nmult]->GetBinContent(bin)/SpectraV0MMB->GetBinContent(bin));
            RatioSpectraV0MSyst[nmult]->SetBinError(bin,erruncorr[nmult]->GetBinContent(bin));
       }
    }
  
    TH1D* RatioSpectraZDC[nbinZDC], * RatioSpectraZDCSyst[nbinZDC];
    for(int nmult = 0; nmult < nbinZDC; nmult++){
        RatioSpectraZDC[nmult] = (TH1D*)SpectraStatZDC[nmult]->Clone(Form("RatioSpectra-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            0.,100.,percentileZDC[nmult],percentileZDC[nmult+1] ));
        RatioSpectraZDC[nmult]->Reset();
        RatioSpectraZDCSyst[nmult] = (TH1D*)SpectraSystZDC[nmult]->Clone(Form("RatioSpectraSyst-Xi-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f", 
            0.,100.,percentileZDC[nmult],percentileZDC[nmult+1] ));
        RatioSpectraZDCSyst[nmult]->Reset();
        //
        for (int bin = 1; bin <= RatioSpectraZDC[0]->GetNbinsX(); bin ++){ //loop over histogram bins
            RatioSpectraZDC[nmult]->SetBinContent(bin,SpectraStatZDC[nmult]->GetBinContent(bin)/SpectraV0MMB->GetBinContent(bin));
            RatioSpectraZDC[nmult]->SetBinError(bin,ErrorInRatioBarlow(SpectraStatZDC[nmult]->GetBinContent(bin),SpectraStatZDC[nmult]->GetBinError(bin),SpectraV0MMB->GetBinContent(bin),SpectraV0MMB->GetBinError(bin)));
            RatioSpectraZDCSyst[nmult]->SetBinContent(bin,SpectraStatZDC[nmult]->GetBinContent(bin)/SpectraV0MMB->GetBinContent(bin));
            RatioSpectraZDCSyst[nmult]->SetBinError(bin,
            ErrorInRatio(SpectraSystZDC[nmult]->GetBinContent(bin),SpectraSystZDC[nmult]->GetBinError(bin),SpectraV0MMB->GetBinContent(bin),huncorrINEL->GetBinContent(bin)));
            //SpectraSystZDC[nmult]->GetBinError(bin)/SpectraV0MMBSyst->GetBinContent(bin));

          //  if (nmult==0) cout << huncorrINEL->GetBinContent(bin) << endl;
        }
    }

    //Scaliamo gli spettri
    Double_t nzdc[nbinZDC];
    for(int nmult = 0; nmult < nbinZDC; nmult++){
        nzdc[nmult] = TMath::Power(2,8-nmult);
        double par = levy[nmult]->GetParameter("norm");
        levy[nmult]->SetParameter("norm",par*(nzdc[nmult]));
        SpectraSystZDC[nmult]->Scale(nzdc[nmult]);
        SpectraStatZDC[nmult]->Scale(nzdc[nmult]);
    }
    Double_t nv0[nbinV0];
    for(int nmult = 0; nmult < nbinV0; nmult++){
        double parV0 = levyV0M[nmult]->GetParameter("norm");
       // cout << parV0 << endl;
        nv0[nmult] = TMath::Power(2,9-nmult);
        levyV0M[nmult]->SetParameter("norm",parV0*(nv0[nmult]));
        SpectraSystTotV0M[nmult]->Scale(nv0[nmult]);
        SpectraStatV0M[nmult]->Scale(nv0[nmult]);
        levyV0M[nmult]->SetRange(0.,7.);
    }

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
    pad1->SetTickx(); 
    pad1->SetTicky();
    pad2->SetTickx(); 
    pad2->SetTicky();
    pad1->SetBottomMargin(0.001);
    pad1->SetBorderMode(0);
    pad1->SetTopMargin(0.01);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.25);
    pad1->SetRightMargin(0.01);
    pad2->SetRightMargin(0.01);
    pad1->SetLeftMargin(0.15);
    pad2->SetLeftMargin(0.15);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();

    //pad1
    pad1->cd();
    pad1->SetLogy();    
    pad1->SetFillColor(kWhite);
    TH1D* h1 = new TH1D("h1","",100,0.,8.01);
    h1->GetYaxis()->SetRangeUser(5E-8,900);
    h1->GetYaxis()->SetLabelSize(0.05);
    h1->GetXaxis()->SetRangeUser(0.,8.);
    h1->GetYaxis()->SetTitle("1/#it{N}_{ev}  d^{2}#it{N}/(d#it{p}_{T}d#it{y}) (GeV/#it{c})^{-1}");
    h1->GetYaxis()->SetTitleOffset(1.5);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->SetStats(0);
    h1->Draw();

    TLegend* l = new TLegend(0.17,0.03,0.35,0.17);
    l->SetTextSize(0.033);
    l->SetTextFont(42);
    l->SetBorderSize(0);     
    TLegend* lw = new TLegend(0.42,0.03,0.63,0.17);
    lw->SetTextSize(0.033);
    lw->SetTextFont(42);
    lw->SetBorderSize(0);
    TLegend* lw1 = new TLegend(0.69,0.03,0.86,0.17);
    lw1->SetTextSize(0.033);
    lw1->SetTextFont(42);
    lw1->SetBorderSize(0);

    TLegend* lTS = new TLegend(0.20,0.25,0.5,0.310);
    lTS->SetTextSize(0.035);
    lTS->SetBorderSize(0); 

    SpectraStatZDC[0]->SetTitle("");
    levy[0]->SetRange(0.,7.);
    SpectraSystZDC[0]->SetFillStyle(0.);
	SpectraSystZDC[0]->Draw("SAME E2");
    SpectraStatZDC[0]->Draw("SAME");
    levy[0]->Draw("SAME");
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

    TLatex *   texclasses = new TLatex(0.2,.6E-5,"(#sqrt{#it{s}} - ZDC) effective energy classes:");
    texclasses->SetTextSize(0.035);
    texclasses->Draw("SAME");
    TLatex *xlabe = new TLatex(1.,3.E-4,"|#it{y}| < 0.5");
    xlabe-> SetTextSize(0.035);
    xlabe->Draw("SAME");

    TLatex *   texbr = new TLatex(4.7,22,"ALICE Preliminary");
    TLatex *   texbrpp = new TLatex(5.2,9,"pp #sqrt{#it{s}} = 13 TeV");
    texbr->SetTextSize(0.03808798);
    texbrpp->SetTextFont(42);
    texbrpp->SetTextSize(0.03808798);
    texbr->Draw("SAME");
    texbrpp->Draw("SAME");

    TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
    xlabel-> SetNDC();
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.08);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    xlabel-> DrawLatex(0.82, 0.77, Form("#Xi^{-} + #bar{#Xi}^{+}"));

    TPavesText *pstr = new TPavesText(3.5,30,4,5);
    pstr->SetBorderSize(0);
    pstr->SetLineWidth(0);
    pstr->SetFillColor(kWhite);
    pstr->SetTextSize(0.03104787);
    pstr->SetTextFont(42);
    TText *pst_LaTexbr = pstr->AddText("stat.");
    pst_LaTexbr = pstr->AddText("syst.");
    pstr->Draw("SAME");
    TMarker *markerbr = new TMarker(3.,30.,2);
    markerbr->SetMarkerStyle(2);
    markerbr->SetMarkerSize(2.9);
    markerbr->Draw("SAME");
    TMarker* marker2br = new TMarker(3.,10.,25);
    marker2br->SetMarkerStyle(25);
    marker2br->SetMarkerSize(2.8);
    marker2br->Draw("SAME");

    //pad2
    pad2->cd();
    pad2->SetLogy();    
    pad2->SetFillColor(kWhite);
    TH1D* h2 = new TH1D("h2","",100,0.,8.01);
    h2->GetYaxis()->SetRangeUser(0.035,15);
    h2->GetXaxis()->SetRangeUser(0.00,8.);
    h2->GetYaxis()->SetTitle("Ratio w.r.t. INEL>0");
    h2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h2->GetYaxis()->SetTitleOffset(.78);
    h2->GetYaxis()->SetTitleSize(0.09);
    h2->GetYaxis()->SetLabelSize(0.09);
    h2->GetXaxis()->SetTitleOffset(.93);
    h2->GetXaxis()->SetTitleSize(0.1);
    h2->GetXaxis()->SetLabelSize(0.08);
    h2->SetTitle("");
    h2->SetStats(0);
    h2->Draw();

    TLine* line = new TLine(0.,1.,8.,1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(9);
    line->SetLineWidth(2);
    line->Draw("SAME");

    RatioSpectraZDCSyst[0]->SetFillStyle(0.);
	RatioSpectraZDCSyst[0]->Draw("SAME E2");
    RatioSpectraZDC[0]->Draw("SAME");
	for(int nmult = 1; nmult < nbinZDC; nmult++){ 
        RatioSpectraZDCSyst[nmult]->SetFillStyle(0.);
        RatioSpectraZDCSyst[nmult]->Draw("SAME E2");
        RatioSpectraZDC[nmult]->Draw("SAME");
    }


    c->cd(2);

    TPad *pad1bis = new TPad("pad1bis","pad1bis",0,0.37,1,1);
    TPad *pad2bis = new TPad("pad2bis","pad2bis",0,0,1,0.37);
    pad1bis->SetTickx(); 
    pad1bis->SetTicky();
    pad2bis->SetTickx(); 
    pad2bis->SetTicky();
    pad1bis->SetBottomMargin(0.001);
    pad1bis->SetBorderMode(0);
    pad1bis->SetTopMargin(0.01);
    pad2bis->SetTopMargin(0.02);
    pad2bis->SetBottomMargin(0.25);
    pad1bis->SetRightMargin(0.1);
    pad2bis->SetRightMargin(0.1);
    pad1bis->SetLeftMargin(0.01);
    pad2bis->SetLeftMargin(0.01);
    pad2bis->SetBorderMode(0);
    pad1bis->Draw();
    pad2bis->Draw();

    //pad1
    pad1bis->cd();
    pad1bis->SetLogy();    
    pad1bis->SetFillColor(kWhite);
    TH1D* h1b = new TH1D("h1","",100,0.,8.01);
    h1b->GetYaxis()->SetRangeUser(5E-8,900);
    h1b->GetYaxis()->SetLabelSize(0.05);
    h1b->GetYaxis()->SetLabelColor(kWhite);
    h1b->GetXaxis()->SetRangeUser(0.,8.);
    h1b->GetYaxis()->SetTitle("");
    h1b->SetStats(0);
    h1b->Draw();

    TLegend* lv = new TLegend(0.03,0.03,0.3,0.18);
    lv->SetTextSize(0.033);
    lv->SetTextFont(42);
    lv->SetBorderSize(0);     
    TLegend* lwv = new TLegend(0.31,0.03,0.57,0.18);
    lwv->SetTextSize(0.033);
    lwv->SetTextFont(42);
    lwv->SetBorderSize(0);
    TLegend* lw1v = new TLegend(0.6,0.03,0.83,0.18);
    lw1v->SetTextSize(0.033);
    lw1v->SetTextFont(42);
    lw1v->SetBorderSize(0);

    TLegend* lTSv = new TLegend(0.03,0.25,0.3,0.310);
    lTSv->SetTextSize(0.035);
    lTSv->SetBorderSize(0);
    
    SpectraSystTotV0M[0]->SetFillStyle(0.);
	SpectraSystTotV0M[0]->Draw("SAME E2");
    SpectraStatV0M[0]->Draw("SAME");
    levyV0M[0]->Draw("SAME");
    lv->AddEntry(SpectraStatV0M[0],"0-1 % ( #times 2^{9} )","LEP");	
	for(int nmult = 1; nmult < nbinV0; nmult++){ 
		SpectraSystTotV0M[nmult]->SetFillStyle(0.);
		SpectraSystTotV0M[nmult]->Draw("SAME E2");
		SpectraStatV0M[nmult]->SetTitle(Form("%.0f-%.0f % ( #times 2^{%i} )",percentileV0[nmult],percentileV0[nmult+1],9-nmult));
	  	SpectraStatV0M[nmult]->Draw("SAME");
        SpectraStatV0M[nmult]->SetStats(0);
        levyV0M[nmult]->Draw("SAME");		
	}

    for(int nmult = 1; nmult < 4; nmult++){
	    lv->AddEntry(SpectraStatV0M[nmult],"","LEP");
    }
    for(int nmult = 4; nmult < 7; nmult++){
	    lwv->AddEntry(SpectraStatV0M[nmult],"","LEP");
    }
    for(int nmult = 7; nmult < nbinV0; nmult++){
	    lw1v->AddEntry(SpectraStatV0M[nmult],"","LEP");
    }
    lTSv->AddEntry(levyV0M[nbinV0-1],"Levy-Tsallis fit","L");
	lv->Draw("SAME");  
    lwv->Draw("SAME");
    lw1v->Draw("SAME");
    lTSv->Draw("SAME");

    TLatex *   texclassesv0 = new TLatex(0.4,.75E-5,"V0M multiplicity classes:");
    texclassesv0->SetTextSize(0.035);
    texclassesv0->Draw("SAME");

    xlabe->Draw("SAME");

    TLatex *   texbrb = new TLatex(4.7,22,"ALICE EPJC80167(2020)");
    texbrb->SetTextSize(0.03808798);
    texbrpp->Draw("SAME");
    texbrb->Draw("SAME");
    
    xlabel-> DrawLatex(0.72, 0.81, Form("#Xi^{-} + #bar{#Xi}^{+}"));

    
    markerbr->Draw("SAME");
    marker2br->Draw("SAME");
    TPavesText *pstrb = new TPavesText(3.5,30,4,5);
    pstrb->SetBorderSize(0);
    pstrb->SetLineWidth(0);
    pstrb->SetFillColor(kWhite);
    pstrb->SetTextSize(0.03104787);
    pstrb->SetTextFont(42);
    TText *pst_LaTexbrb = pstrb->AddText("stat.");
    pst_LaTexbrb = pstrb->AddText("syst.");
    pstrb->Draw("SAME");

    //pad2
    pad2bis->cd();
    pad2bis->SetLogy();    
    pad2bis->SetFillColor(kWhite);
    TH1D* h2b = new TH1D("h2","",100,0.,8.01);
    h2b->GetYaxis()->SetRangeUser(0.035,15);
    h2b->GetYaxis()->SetLabelColor(kWhite);
    h2b->GetXaxis()->SetRangeUser(0.00,8.);
    h2b->SetTitle("");
    h2b->SetStats(0);
    h2b->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h2b->GetXaxis()->SetTitleOffset(.93);
    h2b->GetXaxis()->SetTitleSize(0.1);
    h2b->GetXaxis()->SetLabelSize(0.08);
    h2b->Draw();
    line->Draw("SAME");

    RatioSpectraV0MSyst[0]->SetFillStyle(0.);
    RatioSpectraV0MSyst[0]->Draw("SAME E2");
    RatioSpectraV0M[0]->Draw("SAME");
  	for(int nmult = 1; nmult < nbinV0; nmult++){ 
        RatioSpectraV0MSyst[nmult]->SetFillStyle(0.);
        RatioSpectraV0MSyst[nmult]->Draw("E2 SAME");
      	RatioSpectraV0M[nmult]->Draw("SAME");
    }


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
double ErrorInRatioBarlow(Double_t A, Double_t Aerr, Double_t B, Double_t Berr) {	
  // Error in a Ratio
  if (B != 0) {
    return TMath::Sqrt(TMath::Abs( Aerr*Aerr - Berr*Berr ))/B;
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