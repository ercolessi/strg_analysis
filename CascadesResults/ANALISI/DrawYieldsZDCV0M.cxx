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

#endif
using namespace std;
#include <TString.h>

void DrawYieldsZDCV0M(){
    TFile* fileZDC = TFile::Open("TestfullfitYield-Xi-ZDC-V0M_000_100_ZDC_000_100.root","READ");
    TGraphErrors* hZDCStat = (TGraphErrors*)fileZDC->Get("YieldsNchStat");
    TGraphAsymmErrors* hZDCSyst = (TGraphAsymmErrors*)fileZDC->Get("YieldsNchSyst");

    const int nbinV0 = 10;

    Double_t FYields[nbinV0+1] = {0.1284456594, 0.0924905535, 0.0750700817,0.0609966459, 0.0550614327, 0.0447087709, 0.0333413242, 0.0261796735, 0.0158834917, 0.0061615164};
    Double_t FStatYields[nbinV0+1] = { 0.0029034064, 0.0012247650, 0.0010058078, 0.0008704892, 0.0008879130, 0.0005695243, 0.0005306301, 0.0004594446, 0.0002620426, 0.0001610817 }; 
    Double_t FSystUncorrYields[nbinV0+1] = { 0.0028243299, 0.0018076214, 0.0021215272, 0.0016507798, 0.0022896177, 0.0022972079, 0.0019656341, 0.0016173249, 0.0010877748, 0.0006121448 };
    Double_t FSystYields[nbinV0+1] = { 0.0083107249, 0.0059671232, 0.0053410856, 0.0041266370, 0.0042970402, 0.0036712523, 0.0029716380, 0.0024463140, 0.0015786774, 0.0010696392};
    Double_t OffdNch[nbinV0+1] = {25.75, 19.83, 16.12, 13.76, 12.06, 10.11, 8.07, 6.48, 4.64, 2.52};
    Double_t OffSystdNch[nbinV0+1] = {0.4, 0.3, 0.24, 0.21, 0.18, 0.15, 0.12, 0.1, 0.07, 0.04};
    Double_t OffStatdNch[nbinV0+1] = {0.};

    TGraphErrors *FiorSyst = new TGraphErrors(10,OffdNch,FYields,OffSystdNch,FSystYields);
    TGraphErrors *FiorStat = new TGraphErrors(10,OffdNch,FYields,OffStatdNch,FStatYields);

    TStyle* mcStyle = new TStyle("mcStyle","Francesca's Root Styles");  
    mcStyle->SetPadTickX(1); 
    mcStyle->SetPadTickY(1); 
    mcStyle->SetPalette(1,0); 
    mcStyle->cd();

    TH1D* g = new TH1D("g", " ", 40, 0., 40.);
    g->SetBinContent(1,0.);
    g->SetBinContent(29,2.);
    g->SetLineColor(kWhite);
    TCanvas* c = new TCanvas("c","",1000,900);
    c->SetFillColor(kWhite);
    c->SetRightMargin(0.09);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetGridy();
    c->SetGridx();
    g->GetXaxis()->SetRangeUser(0.,30);
    g->GetYaxis()->SetRangeUser(0.,0.15);
    g->GetYaxis()->SetTitle("#LT d#it{N}/d#it{y} #GT");
    g->GetYaxis()->SetTitleSize(0.05);
    g->GetYaxis()->SetTitleOffset(1.1);
    g->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#eta #GT_{|#eta|<0.5}");
    g->GetXaxis()->SetTitleSize(0.05);
    g->GetXaxis()->SetTitleOffset(0.9);
    g->SetTitle("");
    g->SetStats(0);
    g->Draw();

    TPavesText *pst = new TPavesText(2.258112,0.03395189,5.67125,0.04311871,5,"br");
    pst->SetBorderSize(1);
    pst->SetLineWidth(0);
    pst->SetFillColor(kWhite);
    pst->SetTextFont(42);
    pst->SetTextSize(0.03104787);
    TText *pst_LaTex = pst->AddText("stat.");
    pst_LaTex = pst->AddText("syst.");
    pst->Draw("SAME");
    TMarker *marker = new TMarker(2.777503,0.04075613,2);
    marker->SetMarkerStyle(2);
    marker->SetMarkerSize(2.9);
    marker->Draw();
    marker = new TMarker(2.777503,0.03650348,25);
    marker->SetMarkerStyle(25);
    marker->SetMarkerSize(2.8);
    marker->Draw("SAME");
    TLatex* tex = new TLatex(1.6,.14378024,"ALICE pp, #sqrt{s} = 13 TeV");
    tex->SetTextSize(0.03622251);
    tex->SetTextFont(42);
    //tex->SetLineWidth(2);
    tex->Draw("SAME");


    //FiorSyst->SetFillColorAlpha(kRed-10,0.1);
    FiorSyst->SetLineColor(kRed+1);
    FiorSyst->SetFillStyle(3000);
    FiorSyst->SetLineWidth(1);
    FiorSyst->SetLineColor(kRed+1);
    FiorStat->SetMarkerColor(kRed+1);
    FiorStat->SetMarkerStyle(20);
    FiorStat->SetMarkerSize(2.);
    FiorStat->SetLineColor(kRed+1);
    FiorSyst->SetMarkerColor(kRed+1);   
    FiorSyst->Draw("E2 P SAME");
    FiorStat->Draw("E P SAME");

    hZDCSyst->SetLineColor(kBlue+1);
    hZDCSyst->SetFillStyle(3000);
    hZDCSyst->SetLineWidth(1);
    hZDCSyst->SetLineColor(kBlue+1);
    hZDCStat->SetMarkerColor(kBlue+1);
    hZDCStat->SetMarkerStyle(20);
    hZDCStat->SetMarkerSize(2.);
    hZDCStat->SetLineColor(kBlue+1);
    hZDCSyst->SetMarkerColor(kBlue+1);   
    hZDCSyst->Draw("E2 P SAME");
    hZDCStat->Draw("E P SAME");

    TLegend* lnee = new TLegend (0.3029276,0.6810894,0.506995,0.7818088);
    lnee->SetBorderSize(0);
    lnee->AddEntry(hZDCStat ,"(#sqrt{s} - ZDC), Preliminary","P ");
    lnee->AddEntry(FiorStat, "V0M, EPJC80167(2020)","P ");
    lnee->SetTextSize(0.026);
    lnee->SetTextFont(42);
    
    lnee->Draw("SAME");




}