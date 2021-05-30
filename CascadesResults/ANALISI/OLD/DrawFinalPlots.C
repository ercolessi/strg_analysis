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

void DrawFinalPlots(){

    TFile* fileLowMult  = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","ZDC",70.,100.,0.,100.),"READ");
    TFile* fileHighMult = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","ZDC",0.,30.,0.,100.),"READ");
    //
    TFile* fileLowEE  = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","V0M",0.,100.,70.,100.),"READ");
    TFile* fileHighEE = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","V0M",0.,100.,0.,30.),"READ");
   
    //Nch
  
    TGraphErrors* hAvLowMultStat = (TGraphErrors*)fileLowMult->Get("AvYieldsNchStat");
    TGraphErrors* hAvHighMultStat = (TGraphErrors*)fileHighMult->Get("AvYieldsNchStat");

    TGraphAsymmErrors* hAvLowMultSyst = (TGraphAsymmErrors*)fileLowMult->Get("AvYieldsNchSyst");
    TGraphAsymmErrors* hAvHighMultSyst = (TGraphAsymmErrors*)fileHighMult->Get("AvYieldsNchSyst");

    TGraphAsymmErrors* hAvLowMultSystCorr = (TGraphAsymmErrors*)fileLowMult->Get("AvYieldsNchSystCorr");
    TGraphAsymmErrors* hAvHighMultSystCorr = (TGraphAsymmErrors*)fileHighMult->Get("AvYieldsNchSystCorr");

    TGraphErrors* hAvLowEEStat = (TGraphErrors*)fileLowEE->Get("AvYieldsNchStat");
    TGraphErrors* hAvHighEEStat = (TGraphErrors*)fileHighEE->Get("AvYieldsNchStat");

    TGraphAsymmErrors* hAvLowEESyst = (TGraphAsymmErrors*)fileLowEE->Get("AvYieldsNchSyst");
    TGraphAsymmErrors* hAvHighEESyst = (TGraphAsymmErrors*)fileHighEE->Get("AvYieldsNchSyst");

    TGraphAsymmErrors* hAvLowEESystCorr = (TGraphAsymmErrors*)fileLowEE->Get("AvYieldsNchSystCorr");
    TGraphAsymmErrors* hAvHighEESystCorr = (TGraphAsymmErrors*)fileHighEE->Get("AvYieldsNchSystCorr");

    //percentile
    TGraphErrors* hAvLowMultStatpercentile = (TGraphErrors*)fileLowMult->Get("AvYieldsStatpercentile");
    TGraphErrors* hAvHighMultStatpercentile = (TGraphErrors*)fileHighMult->Get("AvYieldsStatpercentile");

    TGraphAsymmErrors* hAvLowMultSystpercentile = (TGraphAsymmErrors*)fileLowMult->Get("AvYieldsSystpercentile");
    TGraphAsymmErrors* hAvHighMultSystpercentile = (TGraphAsymmErrors*)fileHighMult->Get("AvYieldsSystpercentile");

    TGraphAsymmErrors* hAvLowMultSystCorrpercentile = (TGraphAsymmErrors*)fileLowMult->Get("AvYieldsSystCorrpercentile");
    TGraphAsymmErrors* hAvHighMultSystCorrpercentile = (TGraphAsymmErrors*)fileHighMult->Get("AvYieldsSystCorrpercentile");

    TGraphErrors* hAvLowEEStatpercentile = (TGraphErrors*)fileLowEE->Get("AvYieldsStatpercentile");
    TGraphErrors* hAvHighEEStatpercentile = (TGraphErrors*)fileHighEE->Get("AvYieldsStatpercentile");

    TGraphAsymmErrors* hAvLowEESystpercentile = (TGraphAsymmErrors*)fileLowEE->Get("AvYieldsSystpercentile");
    TGraphAsymmErrors* hAvHighEESystpercentile = (TGraphAsymmErrors*)fileHighEE->Get("AvYieldsSystpercentile");

    TGraphAsymmErrors* hAvLowEESystCorrpercentile = (TGraphAsymmErrors*)fileLowEE->Get("AvYieldsSystCorrpercentile");
    TGraphAsymmErrors* hAvHighEESystCorrpercentile = (TGraphAsymmErrors*)fileHighEE->Get("AvYieldsSystCorrpercentile");
    
    double val[1] = {0.0273555395/6.89};
    double nch[1] = {6.89};
    double errnch[1] = {0.};
    double errsystnch[1] = {0.11};
    double errsyst[1] = {val[0]*TMath::Sqrt(0.11/6.89*0.11/6.89 + 0.0019195267/0.0273555395*0.0019195267/0.0273555395)};
    double errstat[1] = {0.0001852769/6.89};
    const int n = 1; 
    
    TLine* lmb = new TLine(0.,val[0],100,val[0]);
    lmb->SetLineColor(kBlack);
    double upstat = val[0] + errstat[0];
    double downstat = val[0] - errstat[0];
    double upsyst = val[0] + errsyst[0];
    double downsyst = val[0] - errsyst[0];
    TLine* lupstat = new TLine(0,upstat,100,upstat);
    TLine* ldownstat = new TLine(0,downstat,100,downstat);
    lupstat->SetLineColor(kBlack);
    lupstat->SetLineStyle(7);
    ldownstat->SetLineColor(kBlack);
    ldownstat->SetLineStyle(7);

    double bandx[1] = {50.};
    double bandxx[1] = {50.};
    TGraphErrors* hMBband = new TGraphErrors(n,bandx,val,bandxx,errsyst);
    hMBband->SetLineColor(kBlack);
    hMBband->SetFillColor(kGray);
    hMBband->SetFillStyle(1001);

    TGraphErrors* hMBstat = new TGraphErrors(n,nch,val,errnch,errstat);
    TGraphErrors* hMBsyst = new TGraphErrors(n,nch,val,errsystnch,errsyst);

    hMBstat->SetMarkerColor(kBlack);
    hMBstat->SetMarkerStyle(34);
    hMBstat->SetLineColor(kBlack);
    hMBsyst->SetLineColor(kBlack);
    hMBsyst->SetFillStyle(0);
    hMBstat->SetMarkerSize(4);
  
    TCanvas* multlow = new TCanvas("mult", "", 1600,1300);
    multlow->SetLogx();
    multlow->SetFillColor(kWhite);
    multlow->SetTickx();
    multlow->SetTicky();

    TH1D* g3 = new TH1D("g1", " ", 30, 2., 30.);
    g3->SetBinContent(1,0.);
    g3->SetBinContent(10,0.006);
    g3->SetLineColor(kWhite);
    g3->SetStats(0);
    g3->GetYaxis()->SetTitle("#frac{ #LT d#it{N}/d#it{y} #GT }{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} ");
    g3->GetYaxis()->SetTitleSize(0.04);
    g3->GetYaxis()->SetTitleOffset(1.9);
    g3->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{#eta}|<0.5} ");
    g3->GetXaxis()->SetTitleSize(0.04);
    g3->GetYaxis()->SetRangeUser(0.,.008);
    g3->GetXaxis()->SetRangeUser(2.,21.);
    g3->GetXaxis()->SetTitleOffset(1.2);
    g3->Draw();

   
    multlow->SetRightMargin(0.09);
    multlow->SetLeftMargin(0.25);
    multlow->SetBottomMargin(0.15);

    TLegend* lnmultsel = new TLegend (0.29,0.6,0.49,0.8);
    
    hAvHighMultSyst->SetFillColor(kMagenta-9);
    hAvHighMultSyst->SetLineColor(kMagenta+2);
    hAvHighMultSystCorr->SetLineColor(kMagenta+2);
    hAvHighMultSystCorr->SetLineWidth(1);
    hAvHighMultStat->SetLineColor(kMagenta+2);
    hAvHighMultStat->SetMarkerColor(kMagenta+2);
    hAvHighMultStat->SetMarkerStyle(kOpenCircle);
    hAvHighMultStat->SetMarkerSize(2.9);
    hAvLowMultStat->SetMarkerColor(kGreen+3);
    hAvLowMultStat->SetLineColor(kGreen+3);
    hAvLowMultStat->SetMarkerStyle(kFullCircle);
    hAvLowMultSystCorr->SetLineColor(kGreen+3);
    hAvLowMultSystCorr->SetFillColor(kBlack-10);
    hAvLowMultSystCorr->SetFillStyle(0);
    hAvHighMultSystCorr->SetFillStyle(0);
    hAvLowMultSyst->SetLineColor(kGreen+3);
    hAvLowMultSyst->SetFillColor(kGreen-10);
    hAvLowMultSyst->SetFillStyle(1001);
    hAvHighMultSyst->SetFillStyle(1001);
    hAvLowMultStat->SetMarkerSize(2.9);
    hAvLowMultSystCorr->SetMarkerSize(2.9);
    hAvLowMultSystCorr->SetMarkerStyle(kFullCircle);
    hAvLowMultSystCorr->SetMarkerColor(kGreen+3);
    hAvHighMultSystCorr->SetMarkerSize(2.9);
    hAvHighMultSystCorr->SetMarkerStyle(kOpenCircle);
    hAvHighMultSyst->SetMarkerStyle(kOpenCircle);
    hAvHighMultSyst->SetMarkerSize(2.9);
    hAvHighMultSystCorr->SetMarkerColor(kMagenta+2);

    TH1D* hclone_LowMult = (TH1D*)hAvLowMultSyst->Clone("hclone_LowMult");
    TH1D* hclone_HighMult = (TH1D*)hAvHighMultSyst->Clone("hclone_HighMult");
   
    hclone_LowMult->SetFillStyle(3000);
    hclone_LowMult->Draw("E2 SAME");
    hclone_HighMult->SetFillStyle(3000);
    hclone_HighMult->Draw("E2 SAME");

    lnmultsel->SetBorderSize(0);
    lnmultsel->SetHeader("Multiplicity selection:");
    lnmultsel->AddEntry( hAvHighMultSystCorr,"V0M 0-30%, high multiplicity","P");
    lnmultsel->AddEntry(hAvLowMultSystCorr ,"V0M 70-100%, low multiplicity","P");
    lnmultsel->AddEntry(hMBstat ,"INEL>0","P");
    lnmultsel->SetTextSize(0.027);
    lnmultsel->SetTextFont(42);
    
    lnmultsel->Draw("SAME");

    TMarker *arkerc = new TMarker(12.,0.0074509863,2);
    arkerc->SetMarkerStyle(2);
    arkerc->SetMarkerSize(2.9);
    arkerc->Draw("SAME");
    TMarker* arker2c = new TMarker(12.,0.00695041,25);
    arker2c->SetMarkerStyle(25);
    arker2c->SetMarkerSize(2.9);
    arker2c->Draw("SAME");
    TMarker* arker3c = new TMarker(12.,0.00645041,21);
    arker3c->SetMarkerColor(kGray);
    arker3c->SetMarkerStyle(21);
    arker3c->SetMarkerSize(2.9);
    arker3c->Draw("SAME");

    TPavesText *stc = new TPavesText(13.,0.00615,17,.00768);
    stc->SetBorderSize(0);
    stc->SetLineWidth(0);
    stc->SetFillColor(kWhite);
    stc->SetTextSize(0.02904787);
    stc->SetTextFont(42);
    TText *st_LaTexc = stc->AddText("stat.");
    st_LaTexc = stc->AddText("syst. total");
    st_LaTexc = stc->AddText("syst. uncorr.");
    stc->SetTextAlign(12);
    stc->Draw("SAME");


    TLatex *   ttexc = new TLatex(2.3,.0072,"ALICE Preliminary");
    //, pp #sqrt{#it{s}} = 13 TeV");
  //  texc->SetTextFont(42);
    ttexc->SetTextSize(0.0308798);
   // texb->SetLineWidth(2);
   ttexc->Draw("SAME");
    TLatex *   ttexbbb = new TLatex(5.3,0.0072,"pp #sqrt{#it{s}} = 13 TeV");
    ttexbbb->SetTextFont(42);
    ttexbbb->SetTextSize(0.0308798);
   // texb->SetLineWidth(2);
    ttexbbb->Draw("SAME");

    TLatex* xl = new TLatex(0.65, 0.25, Form("#Xi^{-} + #bar{#Xi}^{+}"));
    xl->SetTextFont(42);
    xl-> SetNDC();
    xl-> SetTextColor(1);
    xl-> SetTextSize(0.07);
    xl->Draw("SAME");

    hAvHighMultSyst->Draw("E2 P SAME");
    hAvHighMultSystCorr->Draw("E2 P SAME");
    hAvHighMultStat->Draw("P SAME");
    //
    hAvLowMultSyst->Draw("E2 P SAME");
    hAvLowMultSystCorr->Draw("E2 P SAME");
    hAvLowMultStat->Draw("P SAME");

    hMBsyst->Draw("SAME E2");
    hMBstat->Draw("SAME P");


    multlow->Draw();
    multlow->SaveAs("immaginifinali/multselNch.png");
     multlow->SaveAs("immaginifinali/multselNch.pdf");
      multlow->SaveAs("immaginifinali/multselNch.eps");

    TCanvas*  eelow = new TCanvas(" ee", "", 1600,1300);
    eelow->SetLogx();
    eelow->SetFillColor(kWhite);
    eelow->SetTickx();
    eelow->SetTicky();

    TH1D* sg3 = new TH1D("g1", " ", 30, 2., 28.);
    sg3->SetBinContent(1,0.);
    sg3->SetBinContent(10,0.006);
    sg3->SetLineColor(kWhite);
    sg3->SetStats(0);
    sg3->GetYaxis()->SetTitle("#frac{ #LT d#it{N}/d#it{y} #GT }{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} ");
    sg3->GetYaxis()->SetTitleSize(0.04);
    sg3->GetYaxis()->SetTitleOffset(1.9);
    sg3->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT_{|#it{#eta}|<0.5} ");
    sg3->GetXaxis()->SetTitleSize(0.04);
    sg3->GetYaxis()->SetRangeUser(0.,.008);
    sg3->GetXaxis()->SetRangeUser(2.,28.);
    sg3->GetXaxis()->SetTitleOffset(1.2);
    sg3->Draw();

    TLegend* l = new TLegend (0.6,0.7,0.89,0.89);

    eelow->SetRightMargin(0.09);
    eelow->SetLeftMargin(0.25);
    eelow->SetBottomMargin(0.15);

    hAvHighEESyst->SetLineColor(kOrange+10);
    hAvHighEESyst->SetLineWidth(1);
    hAvHighEESyst->SetFillColor(kOrange-9);
    hAvHighEESyst->SetFillStyle(1001);
    hAvHighEESystCorr->SetLineColor(kOrange+10);
    hAvHighEESystCorr->SetLineWidth(1);
    hAvHighEESystCorr->SetFillColor(kOrange-9);
    hAvHighEESystCorr->SetFillStyle(0);
    hAvHighEEStat->SetLineColor(kOrange+10);
    hAvHighEEStat->SetLineWidth(1);
    hAvHighEEStat->SetMarkerColor(kOrange+10);
    hAvHighEEStat->SetMarkerStyle(25);
    hAvHighEEStat->SetMarkerSize(2.7);
    //
    hAvLowEESyst->SetLineColor(kBlue+1);
    hAvLowEESyst->SetLineWidth(1);
    hAvLowEESyst->SetFillColor(kBlue-9);
    hAvLowEESyst->SetFillStyle(1001);
    hAvLowEESystCorr->SetLineColor(kBlue+1);
    hAvLowEESystCorr->SetLineWidth(1);
    hAvLowEESystCorr->SetFillColor(kBlue-9);
    hAvLowEESystCorr->SetFillStyle(0);
    hAvLowEESystCorr->SetMarkerStyle(21);
    hAvLowEESystCorr->SetMarkerSize(2.7);
    hAvLowEESyst->SetMarkerStyle(21);
    hAvLowEESyst->SetMarkerSize(2.7);
    hAvLowEEStat->SetLineColor(kBlue+1);
    hAvLowEEStat->SetLineWidth(1);
    hAvLowEEStat->SetMarkerColor(kBlue+1);
    hAvLowEEStat->SetMarkerStyle(21);
    hAvLowEEStat->SetMarkerSize(2.7);
    

    hAvHighEESyst->Draw("E2  SAME");
    hAvHighEESystCorr->Draw("E2  SAME");
    hAvHighEEStat->Draw("P SAME");
    //
    hAvLowEESyst->Draw("E2  SAME");
    hAvLowEESystCorr->Draw("E2  SAME");
    hAvLowEEStat->Draw("P SAME");

    
    stc->Draw("SAME");
    ttexc->Draw("SAME");
    ttexbbb->Draw("SAME");
    xl->Draw("SAME");
    TMarker *arkerca = new TMarker(14.5,0.0074509863,2);
    arkerca->SetMarkerStyle(2);
    arkerca->SetMarkerSize(2.9);
    arkerca->Draw("SAME");
    TMarker* arker2ca = new TMarker(14.5,0.00695041,25);
    arker2ca->SetMarkerStyle(25);
    arker2ca->SetMarkerSize(2.9);
    arker2ca->Draw("SAME");
    TMarker* arker3ca = new TMarker(14.5,0.00645041,21);
    arker3ca->SetMarkerColor(kGray);
    arker3ca->SetMarkerStyle(21);
    arker3ca->SetMarkerSize(2.9);
    arker3ca->Draw("SAME");

    TLegend* l2a = new TLegend (0.29029276,0.580894,0.486995,0.7918088);

    l2a->SetBorderSize(0);
    l2a->SetHeader("Effective energy selection:");
    l2a->AddEntry(hAvHighEEStat ,"(#sqrt{#it{s}} - ZDC) 0-30%, high eff. energy","P");
    l2a->AddEntry(hAvLowEEStat ,"(#sqrt{#it{s}} - ZDC) 70-100%, low eff. energy","P");
    l2a->AddEntry(hMBstat ,"INEL>0","P");
    l2a->SetTextSize(0.027);
    l2a->SetTextFont(42);

    l2a->Draw();

    hMBsyst->Draw("SAME E2");
    hMBstat->Draw("SAME P");

    eelow->Draw();
    eelow->SaveAs("immaginifinali/eeselNch.png");
    eelow->SaveAs("immaginifinali/eeselNch.eps");
    eelow->SaveAs("immaginifinali/eeselNch.pdf");

    TCanvas* multlowperc = new TCanvas("multperc", "", 1600,1300);
    //multlowperc->SetGridy();
    multlowperc->SetTickx();
    multlowperc->SetTicky();
    //multlowperc->SetLogx();
    multlowperc->SetFillColor(kWhite);

    TH1D* g3perc = new TH1D("g1perc", " ", 30, 0, 100.);
    g3perc->SetBinContent(1,0.);
    g3perc->SetBinContent(10,0.006);
    g3perc->SetLineColor(kWhite);
    g3perc->SetStats(0);
    g3perc->GetYaxis()->SetTitle("#frac{ #LT d#it{N}/d#it{y} #GT }{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} ");
    g3perc->GetYaxis()->SetTitleSize(0.04);
    g3perc->GetYaxis()->SetTitleOffset(1.9);
    g3perc->GetXaxis()->SetTitle("(#sqrt{#it{s}} - ZDC) percentile");
    g3perc->GetXaxis()->SetTitleSize(0.04);
    g3perc->GetYaxis()->SetRangeUser(0.,.008);
    g3perc->GetXaxis()->SetRangeUser(0.,100.);
    g3perc->GetXaxis()->SetTitleOffset(1.2);
    g3perc->Draw();

    multlowperc->SetRightMargin(0.09);
    multlowperc->SetLeftMargin(0.25);
    multlowperc->SetBottomMargin(0.15);

    hAvHighMultSystpercentile->SetLineColor(kMagenta+2);
    hAvHighMultSystpercentile->SetLineWidth(1);
    hAvHighMultSystpercentile->SetFillColor(kMagenta-9);
    hAvHighMultSystpercentile->SetFillStyle(1001);
    hAvHighMultSystCorrpercentile->SetLineColor(kMagenta+2);
    hAvHighMultSystCorrpercentile->SetLineWidth(1);
    hAvHighMultSystCorrpercentile->SetFillColor(kMagenta-9);
    hAvHighMultSystCorrpercentile->SetFillStyle(0);
    hAvHighMultStatpercentile->SetLineColor(kMagenta+2);
    hAvHighMultStatpercentile->SetLineWidth(1);
    hAvHighMultStatpercentile->SetMarkerColor(kMagenta+2);
    hAvHighMultStatpercentile->SetMarkerStyle(kOpenCircle);
    hAvHighMultStatpercentile->SetMarkerSize(2.9);
    //
    hAvLowMultSystpercentile->SetLineColor(kGreen+3);
    hAvLowMultSystpercentile->SetLineWidth(1);
    hAvLowMultSystpercentile->SetFillColor(kGreen-10);
    hAvLowMultSystpercentile->SetFillStyle(1001);
    hAvLowMultSystCorrpercentile->SetLineColor(kGreen+3);
    hAvLowMultSystCorrpercentile->SetLineWidth(1);
    hAvLowMultSystCorrpercentile->SetFillColor(kGreen-10);
    hAvLowMultSystCorrpercentile->SetFillStyle(0);
    hAvLowMultStatpercentile->SetLineColor(kGreen+3);
    hAvLowMultStatpercentile->SetLineWidth(1);
    hAvLowMultStatpercentile->SetMarkerColor(kGreen+3);
    hAvLowMultStatpercentile->SetMarkerStyle(kFullCircle);
    hAvLowMultStatpercentile->SetMarkerSize(2.9);

    hMBband->Draw("SAME E2");
    lmb->Draw("SAME");
    lupstat->Draw("SAME");
    ldownstat->Draw("SAME");

    TPavesText *pstr = new TPavesText(78.5,0.0062,95,.0077);
    pstr->SetBorderSize(0);
    //pstr->SetLineWidth(0);
    pstr->SetFillColor(kWhite);
    pstr->SetTextSize(0.02904787);
    pstr->SetTextFont(42);
    TText *pst_LaTexbr = pstr->AddText("stat.");
    pst_LaTexbr = pstr->AddText("syst. total");
    pst_LaTexbr = pstr->AddText("syst. uncorr.");
    pstr->SetTextAlign(12);
  
    TMarker *markerbr = new TMarker(75.280374,0.00745,2);
    markerbr->SetMarkerStyle(2);
    markerbr->SetMarkerSize(2.9);
    markerbr->Draw("SAME");
    TMarker* marker2br = new TMarker(75.287895,0.00695,25);
    marker2br->SetMarkerStyle(25);
    marker2br->SetMarkerSize(2.9);
    marker2br->Draw("SAME");
    TMarker* marker3c = new TMarker(75.287895,0.00645041,21);
    marker3c->SetMarkerColor(kGray);
    marker3c->SetMarkerStyle(21);
    marker3c->SetMarkerSize(2.9);
    marker3c->Draw("SAME");
      pstr->Draw("SAME");
    TLatex *   texbr = new TLatex(4.5,0.0073,"ALICE Preliminary");//, pp #sqrt{#it{s}} = 13 TeV");
    TLatex *   texbbr = new TLatex(38,0.0073,"pp #sqrt{#it{s}} = 13 TeV");
    texbbr->SetTextFont(42);
    texbr->SetTextSize(0.03008798);
   // texb->SetLineWidth(2);
    texbr->Draw("SAME");
    texbbr->SetTextSize(0.03008798);
   // texb->SetLineWidth(2);
    texbbr->Draw("SAME");    

    hAvHighMultSystpercentile->Draw("E2  SAME");
    hAvHighMultSystCorrpercentile->Draw("E2  SAME");
    hAvHighMultStatpercentile->Draw("P SAME");
    //
    hAvLowMultSystpercentile->Draw("E2  SAME");
    hAvLowMultSystCorrpercentile->Draw("E2  SAME");
    hAvLowMultStatpercentile->Draw("P SAME");

    TLegend* l21 = new TLegend (0.28229276,0.630894,0.506995,0.8218088);

    l21->SetBorderSize(0);
    l21->SetHeader("Multiplicity selection:");
    l21->AddEntry( hAvHighMultStat,"V0M 0-30%, high multiplicity","P");
    l21->AddEntry( hAvLowMultStat,"V0M 70-100%, low multiplicity","P");
    l21->AddEntry(lmb,"INEL>0, stat. (dashed line) syst. (band)","L");
    l21->SetTextSize(0.027);
    l21->SetTextFont(42);
    l21->Draw("SAME");

    TLatex* lel = new TLatex(0.75, 0.67, Form("#Xi^{-} + #bar{#Xi}^{+}"));
    lel->SetTextFont(42);
    lel-> SetNDC();
    lel-> SetTextColor(1);
    lel-> SetTextSize(0.06);
    lel->Draw("SAME");

   

    multlowperc->Draw();
    multlowperc->SaveAs("immaginifinali/multselperc.png");
    multlowperc->SaveAs("immaginifinali/multselperc.pdf");
    multlowperc->SaveAs("immaginifinali/multselperc.eps");

    TCanvas*  eelowperc = new TCanvas(" eeperc", "", 1600,1300);
     //eelowperc->SetGridy();
     eelowperc->SetTicky();
     eelowperc->SetTickx();
    // eelowperc->SetLogx();
     eelowperc->SetFillColor(kWhite);

      TH1D* g1perc = new TH1D("g1perc", " ", 30, 0., 100.);
    g1perc->SetBinContent(1,0.);
    g1perc->SetBinContent(10,0.006);
    g1perc->SetLineColor(kWhite);
    g1perc->SetStats(0);
    g1perc->GetYaxis()->SetTitle("#frac{ #LT d#it{N}/d#it{y} #GT }{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} ");
    g1perc->GetYaxis()->SetTitleSize(0.04);
    g1perc->GetYaxis()->SetTitleOffset(1.9);
    g1perc->GetXaxis()->SetTitle("V0M percentile");
    g1perc->GetXaxis()->SetTitleSize(0.04);
    g1perc->GetYaxis()->SetRangeUser(0.,.008);
    g1perc->GetXaxis()->SetRangeUser(0.,100.);
    g1perc->GetXaxis()->SetTitleOffset(1.2);
    g1perc->Draw();




  
     eelowperc->SetRightMargin(0.09);
     eelowperc->SetLeftMargin(0.25);
     eelowperc->SetBottomMargin(0.15);

    hAvHighEESystpercentile->SetLineColor(kOrange+10);
    hAvHighEESystpercentile->SetLineWidth(1);
    hAvHighEESystpercentile->SetFillColor(kOrange-9);
    hAvHighEESystpercentile->SetFillStyle(1001);
    hAvHighEESystCorrpercentile->SetLineColor(kOrange+10);
    hAvHighEESystCorrpercentile->SetLineWidth(1);
    hAvHighEESystCorrpercentile->SetFillColor(kOrange-9);
    hAvHighEESystCorrpercentile->SetFillStyle(0);
    hAvHighEEStatpercentile->SetLineColor(kOrange+10);
    hAvHighEEStatpercentile->SetLineWidth(1);
    hAvHighEEStatpercentile->SetMarkerColor(kOrange+10);
    hAvHighEEStatpercentile->SetMarkerStyle(25);
    hAvHighEEStatpercentile->SetMarkerSize(2.7);
    //
    hAvLowEESystpercentile->SetLineColor(kBlue+1);
    hAvLowEESystpercentile->SetLineWidth(1);
    hAvLowEESystpercentile->SetFillColor(kBlue-9);
    hAvLowEESystpercentile->SetFillStyle(1001);
    hAvLowEESystCorrpercentile->SetLineColor(kBlue+1);
    hAvLowEESystCorrpercentile->SetLineWidth(1);
    hAvLowEESystCorrpercentile->SetFillColor(kBlue-9);
    hAvLowEESystCorrpercentile->SetFillStyle(0);
    hAvLowEEStatpercentile->SetLineColor(kBlue+1);
    hAvLowEEStatpercentile->SetLineWidth(1);
    hAvLowEEStatpercentile->SetMarkerColor(kBlue+1);
    hAvLowEEStatpercentile->SetMarkerStyle(21);
    hAvLowEEStatpercentile->SetMarkerSize(2.7);
    
     hMBband->Draw("SAME E2");
    lmb->Draw("SAME");
    lupstat->Draw("SAME");
    ldownstat->Draw("SAME");

   
    markerbr->Draw("SAME");
    marker2br->Draw("SAME");
    marker3c->Draw("SAME");
    pstr->Draw("SAME");
    texbr->Draw("SAME");
    texbbr->Draw("SAME");  

    hAvHighEESystpercentile->Draw("E2 P SAME");
    hAvHighEESystCorrpercentile->Draw("E2 P SAME");
    hAvHighEEStatpercentile->Draw("P SAME");
    //
    hAvLowEESystpercentile->Draw("E2 P SAME");
    hAvLowEESystCorrpercentile->Draw("E2 P SAME");
    hAvLowEEStatpercentile->Draw("P SAME");

        
     TLegend* l2a1 = new TLegend (0.29029276,0.230894,0.486995,0.42088);

    l2a1->SetBorderSize(0);
    l2a1->SetHeader("Effective energy selection:");
    l2a1->AddEntry(hAvHighEEStat ,"(#sqrt{#it{s}} - ZDC) 0-30%, high eff. energy","P");
    l2a1->AddEntry(hAvLowEEStat ,"(#sqrt{#it{s}} - ZDC) 70-100%, low eff. energy","P");
    l2a1->AddEntry(lmb,"INEL>0, stat. (dashed line) syst. (band)","L");
    l2a1->SetTextSize(0.027);
    l2a1->SetTextFont(42);

    lel->Draw("SAME");


    l2a1->Draw("SAME");

    eelowperc->Draw("");
    eelowperc->SaveAs("immaginifinali/eeselperc.png");
    eelowperc->SaveAs("immaginifinali/eeselperc.pdf");
    eelowperc->SaveAs("immaginifinali/eeselperc.eps");
    

}

//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop - errorfrombottom) );
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