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

void beautifyhisto(TH1F* h, Int_t marker, Int_t color);
void beautifyfunc(TF1* f, Int_t style, Int_t color, Int_t width);
void beautifycanvas(TCanvas* c, int n1, int n2);
void Init(Int_t lClassCode,
          vector<Double_t> &percentileSPDCl_low,
          vector<Double_t> &percentileSPDCl_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          vector<Int_t> &colors);

void DrawSpectra(Int_t lClassCode = 0, TString type = "XiPlus")
{

    vector<Double_t> percentileSPDCl_low;
    vector<Double_t> percentileSPDCl_high;
    vector<Double_t> percentileV0M_low;
    vector<Double_t> percentileV0M_high;
    vector<Int_t> colors;

    Init(lClassCode, percentileSPDCl_low, percentileSPDCl_high, percentileV0M_low, percentileV0M_high, colors);

    const Int_t percbinnumb = percentileSPDCl_low.size();

    TH1F *spectra[percbinnumb];
    TFile *file[percbinnumb];
    Double_t n[percbinnumb];
    TH1F *spectraMB;
    TFile *fileMB;

    TString folder = "V0Analysis/results";
    TString v0add = "_FDUseMCRatio";
    if (type.Contains("Xi")) {
        folder = "CascadeAnalysis/results";
        v0add = "";
    }
    if (type.Contains("Lambda")) {
        folder = "V0Analysis/resultsFD";
        v0add = "_FDUseMCRatio";
    }


    for (int i = 0; i < percbinnumb; i++){
        file[i] = TFile::Open(Form("~/strg_analysis/22Sett2023/%s/Results-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f%s.root", folder.Data(), type.Data(), percentileSPDCl_low[i], percentileSPDCl_high[i], percentileV0M_low[i], percentileV0M_high[i], v0add.Data()));
        spectra[i] = (TH1F *)file[i]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
        beautifyhisto(spectra[i], kFullCircle, colors[i]);
        n[i] = pow(2,(percbinnumb-(i+1)));
        spectra[i]->Scale(n[i]);
    }
    fileMB = TFile::Open(Form("~/strg_analysis/22Sett2023/%s/Results-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f%s.root", folder.Data(), type.Data(), 0., 100., 0., 100., v0add.Data()));
    spectraMB = (TH1F *)fileMB->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
    beautifyhisto(spectraMB, kFullCircle, kBlack);
    spectraMB->Sumw2();
    spectraMB->Scale(1./pow(2,percbinnumb));

    TCanvas* c = new TCanvas("c","",800,1000);
    beautifycanvas(c,0,2);
    c->SetLogy();

    TLatex *xlabel = new TLatex();
    xlabel->SetTextFont(42);
    xlabel-> SetNDC();
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.07);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    //
    TH1D* h1 = new TH1D("h1","",100,0.,10.);
    if (type.Contains("K0Short")){
        h1->GetYaxis()->SetRangeUser(1.E-12,50000);
        h1->GetXaxis()->SetRangeUser(0.,10.);
    }
    if (type.Contains("Lambda")){
        h1->GetYaxis()->SetRangeUser(1.E-13, 10000);
        h1->GetXaxis()->SetRangeUser(0.4,8.0);
    }
    if (type.Contains("Xi")){
        h1->GetYaxis()->SetRangeUser(1.E-14, 1000);
        h1->GetXaxis()->SetRangeUser(0.6,6.5);
    }
    h1->GetYaxis()->SetLabelSize(0.04);
    h1->GetYaxis()->SetTitle("1/#it{N}_{ev}  d^{2}#it{N}_{raw}/(d#it{p}_{T}d#it{y}) (GeV/#it{c})^{-1}");
    h1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    h1->GetYaxis()->SetTitleOffset(1.4);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetTitleOffset(0.9);
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->SetStats(0);
    h1->Draw();

    TString classe[] = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"};
    TString classname[] = {"Standalone", "High Multiplicity", "Low Multiplicity", "High ZN", "Low ZN"};

    TLegend* leg = new TLegend(0.2,0.15,0.6,0.35);
    leg->SetTextSize(0.03);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetNColumns(2);
    leg->SetHeader(Form("%s Class:",classname[lClassCode].Data()));

    spectraMB->Draw("SAME EP");
    leg->AddEntry(spectraMB, Form("INEL>0 #times 1/2^{%i}",percbinnumb),"LEP");
    for (int i = 0; i<percbinnumb; i++){
        spectra[i]->Draw("SAME EP");
        leg->AddEntry(spectra[i],Form("%s #times 2^{%i}",classe[i].Data(),percbinnumb-(i+1)),"LEP");
    }
    leg->Draw("SAME");

    TLatex *labbig = new TLatex();
    labbig->SetTextFont(42);
    labbig->SetNDC();
    labbig->SetTextColor(1);
    labbig->SetTextSize(0.07);
    labbig->SetTextAlign(22);
    labbig->SetTextAngle(0);
    if (type.Contains("XiMinus")) labbig->DrawLatex(0.85, 0.89, "#Xi^{-}");
    if (type.Contains("XiPlus")) labbig->DrawLatex(0.85, 0.89, "#bar{#Xi}^{+}");
    if (type.Contains("AntiLambda")) labbig->DrawLatex(0.85, 0.89, "#bar{#Lambda}");
    else if (type.Contains("Lambda")) labbig->DrawLatex(0.85, 0.89, "#Lambda");
    if (type.Contains("K0Short")) labbig->DrawLatex(0.85, 0.89, "K^{0}_{S}");

    TLatex *tex = new TLatex(0.2, 0.9, "ALICE, pp #sqrt{s} = 13 TeV");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    TLatex *tex2 = new TLatex(0.2, 0.86, "This work");
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(2);
    tex->Draw();
    tex2->Draw();

    c->SaveAs(Form("Spectra_%s_class%i.pdf",type.Data(),lClassCode));

}

void Init(Int_t lClassCode,
          vector<Double_t> &percentileSPDCl_low,
          vector<Double_t> &percentileSPDCl_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          vector<Int_t> &colors)
{
    // class 0 --> standalone
    Double_t percentileV0M_low_0[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
    Double_t percentileV0M_high_0[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    Double_t percentileSPDCl_low_0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t percentileSPDCl_high_0[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Int_t colors_0[] = {kRed + 2, kRed, kOrange + 7, kYellow + 1, kSpring - 1, kGreen + 2, kTeal, kAzure + 7, kBlue, kBlue + 2};

    Long_t n0 = sizeof(percentileV0M_low_0) / sizeof(Double_t);

    // class 1 --> kHighMult
    Double_t percentileSPDCl_low_1[] = {10, 10, 10, 10, 10, 10, 10};
    Double_t percentileSPDCl_high_1[] = {20, 20, 20, 20, 20, 20, 20};
    Double_t percentileV0M_low_1[] = {0, 5, 10, 20, 30, 40, 50};
    Double_t percentileV0M_high_1[] = {5, 10, 20, 30, 40, 50, 100};
    Int_t colors_1[] = {kRed + 2, kOrange + 7, kYellow + 1, kSpring - 1, kTeal, kAzure + 7, kBlue + 1};
    Long_t n1 = sizeof(percentileSPDCl_low_1) / sizeof(Double_t);

    // class 2 --> kLowMult
    Double_t percentileSPDCl_low_2[] = {40, 40, 40, 40, 40, 40, 40};
    Double_t percentileSPDCl_high_2[] = {50, 50, 50, 50, 50, 50, 50};
    Double_t percentileV0M_low_2[] = {0, 20, 30, 40, 50, 60, 70};
    Double_t percentileV0M_high_2[] = {20, 30, 40, 50, 60, 70, 100};
    Int_t colors_2[] = {kRed + 1, kOrange + 7, kYellow + 1, kSpring - 1, kTeal, kAzure + 7, kBlue + 1};
    Long_t n2 = sizeof(percentileSPDCl_low_2) / sizeof(Double_t);

    // class 3 --> fixed high ZN
    Double_t percentileSPDCl_low_3[] = {0, 10, 20, 30, 50};
    Double_t percentileSPDCl_high_3[] = {20, 30, 40, 50, 100};
    Double_t percentileV0M_low_3[] = {40, 30, 30, 20, 0};   // 40
    Double_t percentileV0M_high_3[] = {60, 70, 50, 50, 30}; // 60
    Int_t colors_3[] = {kRed + 1, kOrange + 7, kSpring - 1, kAzure + 7, kBlue + 1};
    Long_t n3 = sizeof(percentileSPDCl_low_3) / sizeof(Double_t);

    // class 4 --> fixed low ZN
    Double_t percentileSPDCl_low_4[] = {0, 10, 20, 30};
    Double_t percentileSPDCl_high_4[] = {10, 20, 30, 50};
    Double_t percentileV0M_low_4[] = {20, 10, 0, 0};
    Double_t percentileV0M_high_4[] = {30, 30, 20, 10};
    Int_t colors_4[] = {kRed + 1, kOrange + 7, kAzure + 7, kBlue + 1};
    Long_t n4 = sizeof(percentileSPDCl_low_4) / sizeof(Double_t);

    if (lClassCode == 0)
    {
        for (Int_t i = 0; i < n0; i++)
        {
            percentileSPDCl_low.push_back(percentileSPDCl_low_0[i]);
            percentileSPDCl_high.push_back(percentileSPDCl_high_0[i]);
            percentileV0M_low.push_back(percentileV0M_low_0[i]);
            percentileV0M_high.push_back(percentileV0M_high_0[i]);
            colors.push_back(colors_0[i]);
        }
    }

    if (lClassCode == 1)
    {
        for (Int_t i = 0; i < n1; i++)
        {
            percentileSPDCl_low.push_back(percentileSPDCl_low_1[i]);
            percentileSPDCl_high.push_back(percentileSPDCl_high_1[i]);
            percentileV0M_low.push_back(percentileV0M_low_1[i]);
            percentileV0M_high.push_back(percentileV0M_high_1[i]);
            colors.push_back(colors_1[i]);
        }
    }

    if (lClassCode == 2)
    {
        for (Int_t i = 0; i < n2; i++)
        {
            percentileSPDCl_low.push_back(percentileSPDCl_low_2[i]);
            percentileSPDCl_high.push_back(percentileSPDCl_high_2[i]);
            percentileV0M_low.push_back(percentileV0M_low_2[i]);
            percentileV0M_high.push_back(percentileV0M_high_2[i]);
            colors.push_back(colors_2[i]);
        }
    }

    if (lClassCode == 3)
    {
        for (Int_t i = 0; i < n3; i++)
        {
            percentileSPDCl_low.push_back(percentileSPDCl_low_3[i]);
            percentileSPDCl_high.push_back(percentileSPDCl_high_3[i]);
            percentileV0M_low.push_back(percentileV0M_low_3[i]);
            percentileV0M_high.push_back(percentileV0M_high_3[i]);
            colors.push_back(colors_3[i]);
        }
    }

    if (lClassCode == 4)
    {
        for (Int_t i = 0; i < n4; i++)
        {
            percentileSPDCl_low.push_back(percentileSPDCl_low_4[i]);
            percentileSPDCl_high.push_back(percentileSPDCl_high_4[i]);
            percentileV0M_low.push_back(percentileV0M_low_4[i]);
            percentileV0M_high.push_back(percentileV0M_high_4[i]);
            colors.push_back(colors_4[i]);
        }
    }
}

void beautifyhisto(TH1F* h, Int_t marker, Int_t color){
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetMarkerStyle(marker);
}

void beautifyfunc(TF1* f, Int_t style, Int_t color, Int_t width){
    f->SetLineWidth(width);
    f->SetLineColor(color);
    f->SetLineStyle(style);
}

void beautifycanvas(TCanvas* c1, int n1, int n2){
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.12);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.05);
    c1->SetTicky();
    c1->SetTickx();
    if (n1!=0 || n2!=0) c1->Divide(n1,n2);
}
