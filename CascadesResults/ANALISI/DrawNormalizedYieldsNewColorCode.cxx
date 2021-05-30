#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TH1.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TAttFill.h"

#endif
using namespace std;
#include <TString.h>

double MergeBins(Double_t* p);
double MergeBinsErrors(Double_t* e);
//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
double ErrorInRatioUncorr ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );

void DrawNormalizedYieldsNewColorCode(){
    
    TFile* fileLowMult  = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","ZDC",70.,100.,0.,100.),"READ");
    TFile* fileHighMult = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","ZDC",0.,30.,0.,100.),"READ");
   //
    TFile* fileLowEE  = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","V0M",0.,100.,70.,100.),"READ");
    TFile* fileHighEE = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","V0M",0.,100.,0.,30.),"READ");
   //
  //  TFile* fileV0  = TFile::Open("TestfullfitYield-Xi-V0M-V0M_000_100_ZDC_000_100.root","READ");
    TFile* fileZDC = TFile::Open("Test2SystTestfullfitYield-Xi-ZDC-V0M_000_100_ZDC_000_100.root","READ");
   //
    
    TGraphErrors* hLowMultStat = (TGraphErrors*)fileLowMult->Get("NormYieldsvspercentile_Stat");
    TGraphErrors* hHighMultStat = (TGraphErrors*)fileHighMult->Get("NormYieldsvspercentile_Stat");

    TGraphAsymmErrors* hLowMultSyst = (TGraphAsymmErrors*)fileLowMult->Get("NormYieldsvspercentile_Syst");
    TGraphAsymmErrors* hHighMultSyst = (TGraphAsymmErrors*)fileHighMult->Get("NormYieldsvspercentile_Syst");

    TGraphAsymmErrors* NormYieldsNchStat_EELowMult = (TGraphAsymmErrors*)fileLowMult->Get("NormYieldsNchStat");
    TGraphAsymmErrors* NormYieldsNchSyst_EELowMult = (TGraphAsymmErrors*)fileLowMult->Get("NormYieldsNchSyst");

    TGraphAsymmErrors* NormYieldsNchStat_EEHighMult = (TGraphAsymmErrors*)fileHighMult->Get("NormYieldsNchStat");
    TGraphAsymmErrors* NormYieldsNchSyst_EEHighMult = (TGraphAsymmErrors*)fileHighMult->Get("NormYieldsNchSyst");

    TGraphErrors* hLowEEStat = (TGraphErrors*)fileLowEE->Get("NormYieldsvspercentile_Stat");
    TGraphErrors* hHighEEStat = (TGraphErrors*)fileHighEE->Get("NormYieldsvspercentile_Stat");

    TGraphAsymmErrors* hLowEESyst = (TGraphAsymmErrors*)fileLowEE->Get("NormYieldsvspercentile_Syst");
    TGraphAsymmErrors* hHighEESyst = (TGraphAsymmErrors*)fileHighEE->Get("NormYieldsvspercentile_Syst");

    TGraphAsymmErrors* NormYieldsNchStat_LowEE = (TGraphAsymmErrors*)fileLowEE->Get("NormYieldsNchStat");
    TGraphAsymmErrors* NormYieldsNchSyst_LowEE = (TGraphAsymmErrors*)fileLowEE->Get("NormYieldsNchSyst");

    TGraphAsymmErrors* NormYieldsNchStat_HighEE = (TGraphAsymmErrors*)fileHighEE->Get("NormYieldsNchStat");
    TGraphAsymmErrors* NormYieldsNchSyst_HighEE = (TGraphAsymmErrors*)fileHighEE->Get("NormYieldsNchSyst");
 
    /*TStyle* mcStyle = new TStyle("mcStyle","Francesca's Root Styles");  
    mcStyle->SetPadTickX(1); 
    mcStyle->SetPadTickY(1); 
    mcStyle->SetPalette(1,0); 
    mcStyle->cd();*/
   

    TCanvas* EE = new TCanvas("EE", "", 1700,1400);
    EE->SetFillColor(kWhite);
    EE->SetTicky();
    EE->SetTickx();

    TH1D* g = new TH1D("g", " ", 12, 0., 100.);
    g->SetBinContent(1,-0.02);
    g->SetBinContent(12,2.2);
    g->SetLineColor(kWhite);
    g->SetStats(0);
    g->GetYaxis()->SetTitle("#left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{Sel} / #left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{INEL>0}");
    g->GetYaxis()->SetTitleSize(0.035);
    g->GetYaxis()->SetTitleOffset(1.8);
    g->GetXaxis()->SetTitle("(#sqrt{#it{s}} - ZDC) percentile");
    g->GetXaxis()->SetTitleSize(0.04);
    g->GetXaxis()->SetTitleOffset(1.2);
    g->GetYaxis()->SetRangeUser(0.4,1.6);
    g->SetTitle("");
    g->Draw();

    TLegend* l2 = new TLegend (0.29029276,0.680894,0.486995,0.818088);

    EE->SetRightMargin(0.09);
    EE->SetLeftMargin(0.25);
    EE->SetBottomMargin(0.15);

    hHighMultSyst ->SetFillColorAlpha(kMagenta-10,0.1);
    hHighMultSyst ->Draw("E2 P SAME");
    hHighMultStat ->Draw("E P SAME");
    hHighMultSyst->SetFillStyle(0);
    hHighMultSyst ->SetLineColor(kRed+1);
    hHighMultSyst ->SetLineWidth(1);
    hHighMultStat ->SetLineColor(kRed+1);
    hHighMultStat ->SetMarkerColor(kRed+1);
    hHighMultStat ->SetMarkerStyle(kFullCircle);
    hHighMultStat ->SetMarkerSize(3.8);
    hLowMultStat ->SetMarkerColor(kRed+1);
    hLowMultStat ->SetLineColor(kRed+1);
    hLowMultStat ->SetMarkerStyle(kOpenCircle);
    hLowMultSyst ->SetLineColor(kRed+1);
    hLowMultSyst ->SetFillColor(kGreen-8);
    hLowMultSyst ->SetFillStyle(0);
    hLowMultStat ->SetMarkerSize(3.8);
    hLowMultSyst ->SetMarkerSize(3.8);
    hLowMultSyst ->SetMarkerStyle(kOpenCircle);
    hLowMultSyst ->SetMarkerColor(kRed+1);
    hHighMultSyst ->SetMarkerSize(3.8);
    hHighMultSyst ->SetMarkerStyle(kFullCircle);

    hHighMultSyst ->SetMarkerColor(kRed+1);
    hLowMultSyst ->Draw("E2 P SAME");
    hLowMultStat ->Draw("E P SAME");

    TH1D* hcloneMulthigh = (TH1D*)hHighMultSyst->Clone("hclonehigh");
    TH1D* hcloneMultlow = (TH1D*)hLowMultSyst->Clone("hclonelow");
    hcloneMulthigh->SetFillStyle(3000);
    hcloneMultlow->SetFillStyle(3000);
    hcloneMulthigh->Draw("E2 SAME");
    hcloneMultlow->Draw("E2 SAME");

    TLatex* lel = new TLatex(0.72, 0.7, Form("#Xi^{-} + #bar{#Xi}^{+}"));
    lel->SetTextFont(42);
    lel-> SetNDC();
    lel-> SetTextColor(1);
    lel-> SetTextSize(0.08);
    lel->Draw("SAME");

    l2->SetBorderSize(0);
    l2->SetHeader("V0M multiplicity:","NDF");
    l2->AddEntry(hHighMultSyst ,"0-30%, high mult.","P");
    l2->AddEntry(hLowMultSyst ,"70-100%, low mult.","P");
    l2->SetTextSize(0.037);
    l2->SetTextFont(42);
    l2->Draw("SAME");

    TPavesText *pstr = new TPavesText(85.,1.53,90,1.74);
    pstr->SetBorderSize(0);
    pstr->SetLineWidth(0);
    pstr->SetFillColor(kWhite);
    pstr->SetTextSize(0.02904787);
    pstr->SetTextFont(42);
    TText *pst_LaTexbr = pstr->AddText("stat.");
    pst_LaTexbr = pstr->AddText("syst.");
    pstr->Draw("SAME");
    TMarker *markerbr = new TMarker(80.280374,1.689863,2);
    markerbr->SetMarkerStyle(2);
    markerbr->SetMarkerSize(3.9);
    markerbr->Draw("SAME");
    TMarker* marker2br = new TMarker(80.287895,1.589041,25);
    marker2br->SetMarkerStyle(25);
    marker2br->SetMarkerSize(3.8);
    marker2br->Draw("SAME");
    TLatex *   texbr = new TLatex(4.5,1.68,"ALICE Preliminary");//, pp #sqrt{#it{s}} = 13 TeV");
    TLatex *   texbbr = new TLatex(40,1.68,"pp #sqrt{#it{s}} = 13 TeV");
    texbbr->SetTextFont(42);
    texbr->SetTextSize(0.03808798);
   // texb->SetLineWidth(2);
    texbr->Draw("SAME");
    texbbr->SetTextSize(0.03808798);
   // texb->SetLineWidth(2);
    texbbr->Draw("SAME");

    EE->SaveAs("immaginifinali/D-DiffZDCSelection.pdf");
    EE->SaveAs("immaginifinali/D-DiffZDCSelection.eps");
    EE->SaveAs("immaginifinali/D-DiffZDCSelection.png");
    EE->Draw();

    TCanvas* Mult = new TCanvas("Mult", "", 1700,1400);
    //Mult->SetGridy();
    Mult->SetFillColor(kWhite);
     Mult->SetTicky();
    Mult->SetTickx();

    TH1D* ga = new TH1D("g", " ", 12, 0., 100.);
    ga->SetBinContent(1,-0.02);
    ga->SetBinContent(12,2.2);
    ga->SetLineColor(kWhite);
    ga->SetStats(0);
    ga->GetYaxis()->SetTitle("#left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{Sel} / #left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{INEL>0}");
    ga->GetYaxis()->SetTitleSize(0.035);
    ga->GetYaxis()->SetTitleOffset(1.8);
    ga->GetXaxis()->SetTitle("V0M percentile");
    ga->GetXaxis()->SetTitleSize(0.04);
    ga->GetXaxis()->SetTitleOffset(1.2);
    ga->GetYaxis()->SetRangeUser(0.4,1.6);
    ga->SetTitle("");
    ga->Draw();

    TLegend* l2a = new TLegend (0.29029276,0.230894,0.486995,0.38088);

    Mult->SetRightMargin(0.09);
    Mult->SetLeftMargin(0.25);
    Mult->SetBottomMargin(0.15);

  
    hHighEESyst->SetFillStyle(0);
    hHighEESyst ->SetLineColor(kBlue+1);
    hHighEESyst ->SetLineWidth(1);
    hHighEEStat ->SetLineColor(kBlue+1);
    hHighEEStat ->SetMarkerColor(kBlue+1);
    hHighEEStat ->SetMarkerStyle(21);
    hHighEEStat ->SetMarkerSize(3.7);
    hLowEEStat ->SetMarkerColor(kBlue+1);
    hLowEEStat ->SetLineColor(kBlue+1);
    hLowEEStat ->SetMarkerStyle(25);
    hLowEESyst ->SetLineColor(kBlue+1);
    hLowEESyst ->SetFillColor(kBlue+1);
    hLowEESyst ->SetFillStyle(0);
    hLowEEStat ->SetMarkerSize(3.7);
    hLowEESyst ->SetMarkerSize(3.7);
    hLowEESyst ->SetMarkerStyle(25);
    hLowEESyst ->SetMarkerColor(kBlue+1);
    hHighEESyst ->SetMarkerSize(3.7);
    hHighEESyst ->SetMarkerStyle(21);
    hHighEESyst ->SetMarkerColor(kBlue+1);
    hLowEESyst ->Draw("E2 P SAME");
    hLowEEStat ->Draw("E P SAME");
    hHighEESyst ->Draw("E2 P SAME");
    hHighEEStat ->Draw("E P SAME");

    /*TH1D* hcloneEEhigh = (TH1D*)hHighEESyst->Clone("hclonehighEE");
    TH1D* hcloneEElow = (TH1D*)hLowEESyst->Clone("hclonelowEE");
    hcloneEEhigh->SetFillStyle(3000);
    hcloneEElow->SetFillStyle(3000);
    hcloneEEhigh->Draw("E2 SAME");
    hcloneEElow->Draw("E2 SAME");*/

    l2a->SetBorderSize(0);
    l2a->SetHeader("(#sqrt{#it{s}} - ZDC) effective energy:","");
    l2a->AddEntry(hHighEESyst ,"0-30%, high eff. energy","P");
    l2a->AddEntry(hLowEESyst ,"70-100%, low eff. energy","P");
    l2a->SetTextSize(0.027);
    l2a->SetTextFont(42);
    
    l2a->Draw("SAME");
    lel->Draw("SAME");

    pstr->Draw("SAME");
    texbbr->Draw("SAME");
    markerbr->Draw("SAME");
    marker2br->Draw("SAME");
    texbr->Draw("SAME");

    Mult->SaveAs("immaginifinali/D-DiffV0MSelection.png");
    Mult->SaveAs("immaginifinali/D-DiffV0MSelection.pdf");
    Mult->SaveAs("immaginifinali/D-DiffV0MSelection.eps");


    //Mult Sel
    TCanvas* NchMultsel = new TCanvas("NchMultsel", "", 1700,1400);
     NchMultsel->SetTicky();
    NchMultsel->SetTickx();
    NchMultsel->SetFillColor(kWhite);
    //NchMultsel->SetGridy();
    NchMultsel->SetLogx();
    TH1D* nee = new TH1D("nee", " ", 12, 2., 21);
    nee->SetBinContent(1,0.4);
    nee->SetBinContent(12,1.6);
    nee->SetLineColor(kWhite);
    nee->SetStats(0);
    nee->GetYaxis()->SetTitle("#left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{Sel} / #left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{INEL>0}");
    nee->GetYaxis()->SetTitleSize(0.035);
    nee->GetYaxis()->SetTitleOffset(1.8);
    nee->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}");
    nee->GetXaxis()->SetTitleSize(0.04);
    nee->GetXaxis()->SetTitleOffset(1.2);
    nee->GetYaxis()->SetRangeUser(0.4,1.6);
    nee->SetTitle("");
    nee->Draw();

    TLegend* lnmultsel = new TLegend (0.30,0.65,0.49,0.8);

    NchMultsel->SetRightMargin(0.09);
    NchMultsel->SetLeftMargin(0.25);
    NchMultsel->SetBottomMargin(0.15);
    
    NormYieldsNchSyst_EEHighMult->SetFillColorAlpha(kRed-10,0.1);
    NormYieldsNchSyst_EEHighMult->SetLineColor(kRed+1);
    NormYieldsNchSyst_EEHighMult->SetLineWidth(1);
    NormYieldsNchStat_EEHighMult->SetLineColor(kRed+1);
    NormYieldsNchStat_EEHighMult->SetMarkerColor(kRed+1);
    NormYieldsNchStat_EEHighMult->SetMarkerStyle(kFullCircle);
    NormYieldsNchStat_EEHighMult->SetMarkerSize(3.8);
    NormYieldsNchStat_EELowMult->SetMarkerColor(kRed+1);
    NormYieldsNchStat_EELowMult->SetLineColor(kRed+1);
    NormYieldsNchStat_EELowMult->SetMarkerStyle(kOpenCircle);
    NormYieldsNchSyst_EELowMult->SetLineColor(kRed+1);
    NormYieldsNchSyst_EELowMult->SetFillColor(kBlack-10);
    NormYieldsNchSyst_EELowMult->SetFillStyle(0);
    NormYieldsNchSyst_EEHighMult->SetFillStyle(0);
    NormYieldsNchStat_EELowMult->SetMarkerSize(3.8);
    NormYieldsNchSyst_EELowMult->SetMarkerSize(3.8);
    NormYieldsNchSyst_EELowMult->SetMarkerStyle(kOpenCircle);
    NormYieldsNchSyst_EELowMult->SetMarkerColor(kRed+1);
    NormYieldsNchSyst_EEHighMult->SetMarkerSize(3.8);
    NormYieldsNchSyst_EEHighMult->SetMarkerStyle(kFullCircle);
    NormYieldsNchSyst_EEHighMult->SetMarkerColor(kRed+1);

    NormYieldsNchSyst_EEHighMult ->Draw("E2 P SAME");
    NormYieldsNchStat_EEHighMult ->Draw("E P SAME");

    NormYieldsNchSyst_EELowMult ->Draw("E2 P SAME");
    NormYieldsNchStat_EELowMult ->Draw("E P SAME");

    TH1D* hclone_LowMult = (TH1D*)NormYieldsNchSyst_EELowMult->Clone("hclone_LowMult");
    TH1D* hclone_HighMult = (TH1D*)NormYieldsNchSyst_EEHighMult->Clone("hclone_HighMult");
   
    hclone_LowMult->SetFillStyle(3000);
    hclone_LowMult->Draw("E2 SAME");
    hclone_HighMult->SetFillStyle(3000);
    hclone_HighMult->Draw("E2 SAME");

    lnmultsel->SetBorderSize(0);
    lnmultsel->SetHeader("V0 multiplicity:");
    lnmultsel->AddEntry(hHighMultSyst ,"V0M 0-30%, high mult.","P");
    lnmultsel->AddEntry(hLowMultSyst ,"V0M 70-100%, low mult.","P");
    lnmultsel->SetTextSize(0.037);
    lnmultsel->SetTextFont(42);
    
    lnmultsel->Draw("SAME");

    TMarker *arkerc = new TMarker(12.,1.709863,2);
    arkerc->SetMarkerStyle(2);
    arkerc->SetMarkerSize(3.9);
    arkerc->Draw("SAME");
    TMarker* arker2c = new TMarker(12.,1.62041,25);
    arker2c->SetMarkerStyle(25);
    arker2c->SetMarkerSize(3.8);
    arker2c->Draw("SAME");

    TPavesText *stc = new TPavesText(13.,1.58,17,1.76);
    stc->SetBorderSize(0);
    stc->SetLineWidth(0);
    stc->SetFillColor(kWhite);
    stc->SetTextSize(0.03904787);
    stc->SetTextFont(42);
    TText *st_LaTexc = stc->AddText("stat.");
    st_LaTexc = stc->AddText("syst.");
    stc->Draw("SAME");


    TLatex *   ttexc = new TLatex(2.3,1.68,"ALICE Preliminary");
    //, pp #sqrt{#it{s}} = 13 TeV");
    ttexc->SetTextFont(42);
    ttexc->SetTextSize(0.0408798);
   // texb->SetLineWidth(2);
   ttexc->Draw("SAME");
    TLatex *   ttexbbb = new TLatex(5,1.68,"pp #sqrt{#it{s}} = 13 TeV");
    ttexbbb->SetTextFont(42);
    ttexbbb->SetTextSize(0.0408798);
   // texb->SetLineWidth(2);
    ttexbbb->Draw("SAME");

    TLatex* xl = new TLatex(0.65, 0.25, Form("#Xi^{-} + #bar{#Xi}^{+}"));
    xl->SetTextFont(42);
    xl-> SetNDC();
    xl-> SetTextColor(1);
    xl-> SetTextSize(0.08);
    xl->Draw("SAME");


    NchMultsel->SaveAs("immaginifinali/D-DiffZDCSelV0M_Fixed_Nch.png");
    NchMultsel->SaveAs("immaginifinali/D-DiffZDCSelV0M_Fixed_Nch.pdf");
    NchMultsel->SaveAs("immaginifinali/D-DiffZDCSelV0M_Fixed_Nch.eps");

     //Mult Sel
    TCanvas* NchEEsel = new TCanvas("NchEEsel", "", 1700,1400);
    NchEEsel ->SetTicky();
    NchEEsel->SetTickx();
    NchEEsel->SetFillColor(kWhite);
    //NchEEsel->SetGridy();
    NchEEsel->SetLogx();
    TH1D* neea = new TH1D("nee", " ", 12, 2, 28);
    neea->SetBinContent(1,0.4);
    neea->SetBinContent(12,1.6);
    neea->SetLineColor(kWhite);
    neea->SetStats(0);
    neea->GetYaxis()->SetTitle("#left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{Sel} / #left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{INEL>0}");
    neea->GetYaxis()->SetTitleSize(0.035);
    neea->GetYaxis()->SetTitleOffset(1.8);
    neea->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}");
    neea->GetXaxis()->SetTitleSize(0.04);
    neea->GetXaxis()->SetTitleOffset(1.2);
    neea->GetYaxis()->SetRangeUser(0.4,1.6);
    neea->SetTitle("");
    neea->Draw();

    TLegend* lnEEsel = new TLegend (0.28,0.65,0.49,0.8);

    NchEEsel->SetRightMargin(0.05);
    NchEEsel->SetLeftMargin(0.2);
    NchEEsel->SetBottomMargin(0.15);
    
    NormYieldsNchSyst_HighEE->SetFillColorAlpha(kRed-10,0.1);
    NormYieldsNchSyst_HighEE->SetLineColor(kBlue+1);
    NormYieldsNchSyst_HighEE->SetLineWidth(1);
    NormYieldsNchStat_HighEE->SetLineColor(kBlue+1);
    NormYieldsNchStat_HighEE->SetMarkerColor(kBlue+1);
    NormYieldsNchStat_HighEE->SetMarkerStyle(21);
    NormYieldsNchStat_HighEE->SetMarkerSize(3.7);
    NormYieldsNchStat_LowEE->SetMarkerColor(kBlue+1);
    NormYieldsNchStat_LowEE->SetLineColor(kBlue+1);
    NormYieldsNchStat_LowEE->SetMarkerStyle(25);
    NormYieldsNchSyst_LowEE->SetLineColor(kBlue+1);
    NormYieldsNchSyst_LowEE->SetFillColor(1-10);
    NormYieldsNchSyst_LowEE->SetFillStyle(0);
    NormYieldsNchSyst_HighEE->SetFillStyle(0);
    NormYieldsNchStat_LowEE->SetMarkerSize(3.7);
    NormYieldsNchSyst_LowEE->SetMarkerSize(3.7);
    NormYieldsNchSyst_LowEE->SetMarkerStyle(25);
    NormYieldsNchSyst_LowEE->SetMarkerColor(kBlue+1);
    NormYieldsNchSyst_HighEE->SetMarkerSize(3.7);
    NormYieldsNchSyst_HighEE->SetMarkerStyle(21);
    NormYieldsNchSyst_HighEE->SetMarkerColor(kBlue+1);

    NormYieldsNchSyst_HighEE ->Draw("E2 P SAME");
    NormYieldsNchStat_HighEE ->Draw("E P SAME");

    NormYieldsNchSyst_LowEE ->Draw("E2 P SAME");
    NormYieldsNchStat_LowEE ->Draw("E P SAME");

    TH1D* hclone_LowEE = (TH1D*)NormYieldsNchSyst_LowEE->Clone("hclone_LowEE");
    TH1D* hclone_HighEE = (TH1D*)NormYieldsNchSyst_HighEE->Clone("hclone_HighEE");
   
    hclone_LowEE->SetFillStyle(3000);
    hclone_LowEE->Draw("E2 SAME");
    hclone_HighEE->SetFillStyle(3000);
    hclone_HighEE->Draw("E2 SAME");
    
    TMarker *akerc = new TMarker(15.,1.709863,2);
    akerc->SetMarkerStyle(2);
    akerc->SetMarkerSize(4.9);
    akerc->Draw("SAME");
    TMarker* aker2c = new TMarker(15.,1.62041,25);
    aker2c->SetMarkerStyle(25);
    aker2c->SetMarkerSize(4.8);
    aker2c->Draw("SAME");

    stc->Draw("SAME");
    ttexc->Draw("SAME");
    

    TLatex *   ttxbbb = new TLatex(6,1.68,"pp #sqrt{#it{s}} = 13 TeV");
    ttxbbb->SetTextFont(42);
    ttxbbb->SetTextSize(0.0408798);
   // texb->SetLineWidth(2);
    ttxbbb->Draw("SAME");

  

    lnEEsel->SetBorderSize(0);
    lnEEsel->SetHeader("(#sqrt{#it{s}} - ZDC) effective energy:");
    lnEEsel->AddEntry(hHighEESyst ,"0-30%, high eff. energy","P");
    lnEEsel->AddEntry(hLowEESyst ,"70-100%, low eff. energy","P");
    lnEEsel->SetTextSize(0.032);
    lnEEsel->SetTextFont(42);
 
    lnEEsel->Draw("SAME");
    xl->Draw("SAME");

  
 



    NchEEsel->SaveAs("immaginifinali/ppD-DiffV0MSelZDC_Fixed_Nch.png");
     NchEEsel->SaveAs("immaginifinali/ppD-DiffV0MSelZDC_Fixed_Nch.eps");
      NchEEsel->SaveAs("immaginifinali/ppD-DiffV0MSelZDC_Fixed_Nch.pdf");


    /////////////////////////////////////////////////////////////////////////////

    Double_t percentileV0Full[] = {0.,1.,5,10,15,20,30,40,50,70,100};
    const int nbinV0 = 9;
    //    
    double errpercentileFull[nbinV0+1], centrpercentileFull[nbinV0+1],err0Full[nbinV0+1];
    for (Int_t n = 0; n < nbinV0+1; n++) {
      centrpercentileFull[n] = (percentileV0Full[n] + percentileV0Full[n + 1]) / 2;
      errpercentileFull[n] = (percentileV0Full[n + 1] - percentileV0Full[n]) / 2;
      err0Full[n] = 0;
    }

    TGraphErrors* hEEStats = (TGraphErrors*)fileZDC->Get("NormYieldsvspercentile_Stat");
    TGraphAsymmErrors* hEESyst = (TGraphAsymmErrors*)fileZDC->Get("NormYieldsvspercentile_Syst");
    TGraphAsymmErrors* NormYieldsNchStat_EE = (TGraphAsymmErrors*)fileZDC->Get("NormYieldsNchStat");
    TGraphAsymmErrors* NormYieldsNchSyst_EE = (TGraphAsymmErrors*)fileZDC->Get("NormYieldsNchSyst");

    
    Double_t FYields[nbinV0+1] = {0.1284456594, 0.0924905535, 0.0750700817,0.0609966459, 0.0550614327, 0.0447087709, 0.0333413242, 0.0261796735, 0.0158834917, 0.0061615164};
    Double_t FStatYields[nbinV0+1] = { 0.0029034064, 0.0012247650, 0.0010058078, 0.0008704892, 0.0008879130, 0.0005695243, 0.0005306301, 0.0004594446, 0.0002620426, 0.0001610817 }; 
    Double_t FSystYields[nbinV0+1] = { 0.0028243299, 0.0018076214, 0.0021215272, 0.0016507798, 0.0022896177, 0.0022972079, 0.0019656341, 0.0016173249, 0.0010877748, 0.0006121448 };
    Double_t OffdNch[nbinV0+1] = {25.75, 19.83, 16.12, 13.76, 12.06, 10.11, 8.07, 6.48, 4.64, 2.52};
    Double_t OffSystdNch[nbinV0+1] = {0.4, 0.3, 0.24, 0.21, 0.18, 0.15, 0.12, 0.1, 0.07, 0.04};
    Double_t OffStatdNch[nbinV0+1] = {0.};
    Double_t FYieldMB = 0.0273555395;
    Double_t FStatYieldMB = 0.0001852769;
    Double_t OffdNchMB = 6.89;
    Double_t OffSystdNchMB = 0.11;

    Double_t FNormYields[nbinV0+1];
    Double_t FStatNormYields[nbinV0+1];
    Double_t FSystNormYields[nbinV0+1];

    for (int i = 0; i < nbinV0+1; i++){
        FNormYields[i] = ( FYields[i] / OffdNch[i] ) / ( FYieldMB / OffdNchMB );
        FStatNormYields[i] = ErrorInRatio(FYields[i],FStatYields[i], FYieldMB,FStatYieldMB)*OffdNchMB/OffdNch[i];
        FSystNormYields[i] = FNormYields[i]*
            TMath::Sqrt(
                (FSystYields[i] / FYields[i])*(FSystYields[i] / FYields[i]) + 
                (OffSystdNch[i] / OffdNch[i] )*(OffSystdNch[i] / OffdNch[i] ) +
                (OffSystdNchMB / OffdNchMB )* (OffSystdNchMB / OffdNchMB )
            );        
    }

    Double_t NY7[nbinV0+1] = { 1.298 , 1.233, 1.244, 1.179, 1.126, 1.091, 1.013, .9484, .8211,  0.5825};
    Double_t StatNY7[nbinV0+1] = { .02159, .01279, .0131 , .01316, .01378, .01161, .01197, .01422, .01191, .01198   };
    Double_t SystNY7[nbinV0+1] = { 0.04577,.0439, .04987, .05137, .05827, .05886, .06235,  .06873, .07593,  .07689 };
    Double_t Nch7[nbinV0+1] = { 21.29, 16.51, 13.46, 11.51, 10.08, 8.45, 6.72, 5.4, 3.9, 2.26 };
    Double_t SystNch7[nbinV0+1] = {0.64, 0.5,0.4, 0.35, 0.3,0.25,0.21,0.17,0.14,0.12};

    TGraphErrors *Stat7 = new TGraphErrors(nbinV0+1,Nch7,NY7,OffStatdNch,StatNY7);
    TGraphErrors *Syst7 = new TGraphErrors(nbinV0+1,Nch7,NY7,SystNch7,SystNY7);

    TGraphErrors *FiorStat = new TGraphErrors(nbinV0+1,OffdNch,FNormYields,OffStatdNch,FStatNormYields);
    TGraphErrors *FiorSyst = new TGraphErrors(nbinV0+1,OffdNch,FNormYields,OffSystdNch,FSystNormYields);

    TGraphErrors *FiorStatPerc = new TGraphErrors(nbinV0+1,centrpercentileFull,FNormYields,err0Full,FStatNormYields);
    TGraphErrors *FiorSystPerc = new TGraphErrors(nbinV0+1,centrpercentileFull,FNormYields,errpercentileFull,FSystNormYields);

    TCanvas* multV0 = new TCanvas("multV0", "", 1700,1400);
     multV0->SetTicky();
    multV0->SetTickx();
   
    TH1D* gp = new TH1D("g2", " ", 12, -1., 100.);
    gp->SetBinContent(1,0.41);
    gp->SetBinContent(12,1.61);
    gp->SetLineColor(kWhite);
    gp->SetStats(0);
    gp->GetYaxis()->SetTitle("#left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{Sel} / #left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{INEL>0}");
    gp->GetYaxis()->SetTitleSize(0.035);
    gp->GetYaxis()->SetTitleOffset(1.6);
    gp->GetXaxis()->SetTitle("percentile");
    gp->GetXaxis()->SetTitleSize(0.04);
    gp->GetXaxis()->SetTitleOffset(1.2);
    gp->GetYaxis()->SetRangeUser(0.4,1.6);
    gp->SetTitle("");
    gp->Draw();

    TLegend* lp = new TLegend (0.66,0.7,0.89,0.89);

    multV0->SetRightMargin(0.09);
    multV0->SetLeftMargin(0.25);
    multV0->SetBottomMargin(0.15);
    multV0->SetFillColor(kWhite);

    hEEStats->SetMarkerColor(kRed+1);
    hEEStats->SetLineColor(kRed+1);
    hEEStats->SetMarkerStyle(20);
    hEESyst->SetLineColor(kRed+1);
    hEESyst->SetFillColor(kBlack-10);
    hEESyst->SetFillStyle(0);
    hEEStats->SetMarkerSize(3.5);
    hEESyst->SetMarkerSize(3.5);
    hEESyst->SetMarkerStyle(20);
    hEESyst->SetMarkerColor(kRed+1);
    hEEStats->SetLineWidth(2);

    FiorSystPerc->SetFillColorAlpha(kRed-10,0.1);
    FiorSystPerc->SetLineColor(kBlue+1);
    FiorSystPerc->SetFillStyle(0);
    FiorSystPerc->SetLineWidth(1);
    FiorStatPerc->SetLineColor(kBlue+1);
    FiorStatPerc->SetMarkerColor(kBlue+1);
    FiorStatPerc->SetMarkerStyle(21);
    FiorStatPerc->SetMarkerSize(3.3);
    FiorStatPerc->SetLineWidth(2);
    FiorSystPerc->SetMarkerSize(3.3);
    FiorSystPerc->SetMarkerStyle(21);
    FiorSystPerc->SetMarkerColor(kBlue+1);
    //FiorSystPerc->SetLineWidth(2);
   // hEESyst->SetLineWidth(2);

    TH1D* hcloneee = (TH1D*)hEESyst->Clone("hcloneee");
    TH1D* hclonemult = (TH1D*)FiorSystPerc->Clone("hclonemult");
    hcloneee->SetFillStyle(3000);
    hclonemult->SetFillStyle(3000);

    FiorSystPerc->Draw("E2 P SAME");
    FiorStatPerc->Draw("EP SAME");
    hclonemult->Draw("E2 SAME");
    hEESyst->Draw("E2 P SAME");
    hEEStats->Draw("P SAME");
    hcloneee->Draw("E2 SAME");

    TPavesText *pstb = new TPavesText(85.,1.53,90,1.74);
    pstb->SetBorderSize(0);
    pstb->SetLineWidth(0);
    pstb->SetFillColor(kWhite);
    pstb->SetTextSize(0.02904787);
    pstb->SetTextFont(42);
    TText *pst_LaTexb = pstb->AddText("stat.");
    pst_LaTexb = pstb->AddText("syst.");
    pstb->Draw("SAME");
    TMarker *markerb = new TMarker(80.280374,1.689863,2);
    markerb->SetMarkerStyle(2);
    markerb->SetMarkerSize(2.9);
    markerb->Draw("SAME");
    TMarker* marker2b = new TMarker(80.287895,1.589041,25);
    marker2b->SetMarkerStyle(25);
    marker2b->SetMarkerSize(2.8);
    marker2b->Draw("SAME");
    TLatex *   texb = new TLatex(4.5,1.65,"ALICE Preliminary");//, pp #sqrt{#it{s}} = 13 TeV");
    TLatex *   texbb = new TLatex(42,1.65,"pp #sqrt{#it{s}} = 13 TeV");
    texbb->SetTextFont(42);
    texb->SetTextSize(0.0308798);
   // texb->SetLineWidth(2);
    texb->Draw("SAME");
    texbb->SetTextSize(0.0308798);
   // texb->SetLineWidth(2);
    texbb->Draw("SAME");

    TLegend* lnee = new TLegend (0.2829276,0.700894,0.56995,0.808088);
    lnee->SetBorderSize(0);
    lnee->AddEntry(hEEStats ,"(#sqrt{#it{s}} - ZDC), Preliminary","P ");
    lnee->AddEntry(FiorStatPerc, "V0M, EPJC80167(2020)","P ");
    lnee->SetTextSize(0.04);
    lnee->SetTextFont(42);

    TLatex* label = new TLatex(0.35, 0.25, Form("#Xi^{-} + #bar{#Xi}^{+}"));
    label->SetTextFont(42);
    label-> SetNDC();
    label-> SetTextColor(1);
    label-> SetTextSize(0.08);
    label->Draw("SAME");
    
    lnee->Draw("SAME");

    //multV0->SaveAs("immaginifinali/XiNormYieldsvsperc_ZDCV0M_pp13TeV.eps");
    multV0->SaveAs("immaginifinali/NormYields.png");
    multV0->SaveAs("immaginifinali/NormYields.pdf");

    TCanvas* NchEE = new TCanvas("NchEE", "", 1700,1400);
     NchEE->SetTicky();
    NchEE->SetTickx();
    NchEE->SetFillColor(kWhite);
    NchEE->SetLogx();

    TH1D* neeb = new TH1D("nee", " ", 12, 2., 30);
    neeb->SetBinContent(1,0.4);
    neeb->SetBinContent(12,1.6);
    neeb->SetLineColor(kWhite);
    neeb->SetStats(0);
    neeb->GetYaxis()->SetTitle("#left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{Sel} / #left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{INEL>0}");
    neeb->GetYaxis()->SetTitleSize(0.045);
    neeb->GetYaxis()->SetTitleOffset(1.8);
    neeb->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}");
    neeb->GetXaxis()->SetTitleSize(0.04);
    neeb->GetXaxis()->SetTitleOffset(1.2);
    neeb->GetYaxis()->SetRangeUser(0.4,1.6);
    neeb->SetTitle("");

    neeb->Draw();

    NchEE->SetRightMargin(0.05);
    NchEE->SetLeftMargin(0.22);
    NchEE->SetBottomMargin(0.15);

    NormYieldsNchSyst_EE->SetFillStyle(0);
    NormYieldsNchStat_EE ->SetMarkerColor(kRed+1);
    NormYieldsNchStat_EE ->SetLineColor(kRed+1);
    NormYieldsNchStat_EE ->SetMarkerStyle(20);
    NormYieldsNchSyst_EE ->SetLineColor(kRed+1);
   // NormYieldsNchSyst_EE ->SetFillColor(kBlue+2-10);
    NormYieldsNchStat_EE ->SetMarkerSize(3.5);
   // NormYieldsNchSyst_EE ->SetMarkerSize(2.2);
   // NormYieldsNchSyst_EE ->SetMarkerStyle(24);
 //   NormYieldsNchSyst_EE ->SetLineWidth(2);
    NormYieldsNchSyst_EE ->SetMarkerColor(kBlue+1);
   

  //  FiorSyst->SetFillColorAlpha(kRed-10,0.1);
    FiorSyst->SetLineColor(kBlue+1);
    FiorSyst->SetFillStyle(0);
    FiorSyst->SetLineWidth(1);
    FiorSyst->SetLineColor(kBlue+1);
    FiorStat->SetMarkerColor(kBlue+1);
    FiorStat->SetMarkerStyle(21);
    FiorStat->SetMarkerSize(3.3);
    FiorStat->SetLineColor(kBlue+1);
    FiorSyst->SetMarkerColor(kBlue+1);   
    FiorSyst->Draw("E2 P SAME");
    FiorStat->Draw("E P SAME");
    //FiorSyst->SetLineWidth(2);
    NormYieldsNchSyst_EE ->Draw("E2 P SAME");
    NormYieldsNchStat_EE ->Draw("E P SAME");
    TLatex* xlabel = new TLatex(0.65, 0.25, Form("#Xi^{-} + #bar{#Xi}^{+}"));
    xlabel->SetTextFont(42);
    xlabel-> SetNDC();
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.08);
    xlabel->Draw("SAME");

   /* Syst7->SetFillColorAlpha(kBlue+2,0.1);
    Syst7->SetLineColor(kBlack);
    Syst7->SetFillStyle(3000);
    Syst7->SetLineWidth(1);
    Syst7->SetLineColor(kBlack);
    Stat7->SetMarkerColor(kBlack);
    Stat7->SetMarkerStyle(33);
    Stat7->SetMarkerSize(2.6);
    Stat7->SetLineColor(kBlack);
    Syst7->SetMarkerColor(kBlack);   
    Syst7->Draw("E2 P SAME");
    Stat7->Draw("E P SAME");*/
    
   // NormYieldsNchStat_Mult->Draw("SAME");
   // NormYieldsNchSyst_Mult->Draw("E2 SAME");

    TPavesText *pstc = new TPavesText(24.,1.58,25,1.76);
    pstc->SetBorderSize(0);
    pstc->SetLineWidth(0);
    pstc->SetFillColor(kWhite);
    pstc->SetTextSize(0.02904787);
    pstc->SetTextFont(42);
    TText *pst_LaTexc = pstc->AddText("stat.");
    pst_LaTexc = pstc->AddText("syst.");
    pstc->Draw("SAME");

    TMarker *markerc = new TMarker(20.,1.509863,2);
    markerc->SetMarkerStyle(2);
    markerc->SetMarkerSize(3.5);
    markerc->Draw("SAME");
    TMarker* marker2c = new TMarker(20.,1.32041,25);
    marker2c->SetMarkerStyle(25);
    marker2c->SetMarkerSize(3.5);
    marker2c->Draw("SAME");


    TLatex *   texc = new TLatex(2.3,1.68,"ALICE Preliminary");
    //, pp #sqrt{#it{s}} = 13 TeV");
    texc->SetTextFont(42);
    texc->SetTextSize(0.0408798);
   // texb->SetLineWidth(2);
    texc->Draw("SAME");
    TLatex *   texbbb = new TLatex(6,1.68,"pp #sqrt{#it{s}} = 13 TeV");
    texbbb->SetTextFont(42);
    texbbb->SetTextSize(0.0408798);
   // texb->SetLineWidth(2);
    texbbb->Draw("SAME");

    TLatex *   ty = new TLatex(2.3,1.68,"|y|<0.5");
    ty->Draw("SAME");

   
    
    lnee->Draw("SAME");

    
    NchEE->Draw();
    NchEE->SaveAs("immaginifinali/XiNormYieldsvsNch_ZDCV0M_pp13TeV.eps");
    NchEE->SaveAs("immaginifinali/D-DiffV0MZDC_Nch.pdf");
    NchEE->SaveAs("immaginifinali/D-DiffV0MZDC_Nch.png");


}


Double_t MergeBinsErrors(Double_t* e){

    double e1 = e[0]/5;
    double e2 = e[1]*4/5;

    double out = TMath::Sqrt( e1*e1 + e2*e2 ) ;

    return out;
}

Double_t MergeBins(Double_t* p){

    double out = (p[0] + p[1]*4)/5;

    return out;
}

//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
       // return TMath::Sqrt(TMath::Abs( Aerr*Aerr - Berr*Berr ))/B;
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop + errorfrombottom) );
    }
    return 1.;
}

//---------------------------------------------------------------------------------------------------
double ErrorInRatioUncorr ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop + errorfrombottom) );
    }
    return 1.;
}

/*markers(int lineWidth=0) { // Display the table of markers with their numbers. 
    TMarker *marker = new TMarker(); 
    marker->SetMarkerSize(3); 
    TText *text = new TText(); 
    text->SetTextFont(62); 
    text->SetTextAlign(22); 
    text->SetTextSize(0.1); 
    char atext[] = " "; 
    Double_t x = 0; 
    Double_t dx = 1/12.0; 
    for (Int_t i=1;i<12;i++) { 
        x += dx; sprintf(atext,"%d",i); 
        marker->SetMarkerStyle(i+ 1000*lineWidth); 
        marker->DrawMarker(x,.35); 
        text->DrawText(x,.17,atext); 
        sprintf(atext,"%d",i+19); 
        marker->SetMarkerStyle(i+19 + 1000*lineWidth); 
        marker->DrawMarker(x,.8); 
        text->DrawText(x,.62,atext); 
    } 
    delete marker; 
    delete text; 
}*/