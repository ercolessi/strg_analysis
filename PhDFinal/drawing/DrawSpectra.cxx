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
Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p);
TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t par2min = 1e-3, Double_t par2max = 1e+3,Double_t par3min = 1e-5, Double_t par3max = 1e+5, Double_t n = 5., Double_t C = 0.1, Double_t norm = 0.1);

void DrawSpectra(
    TString type = "K0Short",
    TString var = "V0M",
    TString fixed = "SPDClusters",
    Double_t lLow = 0.,
    Double_t lHigh = 100.
    ){

        Double_t * percentile;
        Long_t percbinnumb;
        Double_t pmb[] = {0,100};
        Long_t nmb = sizeof(pmb)/sizeof(Double_t) - 1;
        Double_t p0[] = {0,1,5,10,15,20,30,40,50,70,100};
        Long_t n0 = sizeof(p0)/sizeof(Double_t) - 1;
        Double_t p1[] = {0,5,10,20,30,40,50,100};
        Long_t n1 = sizeof(p1)/sizeof(Double_t) - 1;
        Double_t p2[] = {0,20,30,40,50,60,70,100};
        Long_t n2 = sizeof(p2)/sizeof(Double_t) - 1;
        Double_t p4[] = {0,5,10,20,30,40,50,100};
        Long_t n4 = sizeof(p4)/sizeof(Double_t) - 1;
        Double_t p5[] = {0,10,20,30,40,50,60,70,100};
        Long_t n5 = sizeof(p5)/sizeof(Double_t) - 1;
        Double_t pOmega[] = {0,5,10,30,50,100};
        Long_t npOmega = sizeof(pOmega)/sizeof(Double_t) - 1;
        Double_t pOmega2[] = {0,40,70,100};
        Long_t npOmega2 = sizeof(pOmega2)/sizeof(Double_t) - 1;
        Double_t pOmega3[] = {0,5,10,30,100};
        Long_t npOmega3 = sizeof(pOmega3)/sizeof(Double_t) - 1;
        Double_t pOmega4[] = {0,30,50,100};
        Long_t npOmega4 = sizeof(pOmega4)/sizeof(Double_t) - 1;


        if (fixed.Contains("SPD") && lLow == 0. && lHigh == 100.){
            if (type.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = npOmega;
            } else{
            percentile = p0;
            percbinnumb = n0;
            }
        }
        if (fixed.Contains("SPD") && lLow == 10. && lHigh == 20.){
            if (type.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = npOmega;
            } else{
            percentile = p1;
            percbinnumb = n1;
            }
        }
        if (fixed.Contains("SPD") && lLow == 40. && lHigh == 50.){
            if (type.Contains("Omega")){
            percentile = pOmega2;
            percbinnumb = npOmega2;
            } else{
            percentile = p2;
            percbinnumb = n2;
            }
        }
        if (fixed.Contains("V0M") && lLow == 10. && lHigh == 20.){
            if (type.Contains("Omega")){
            percentile = pOmega3;
            percbinnumb = npOmega3;
            } else{
            percentile = p4;
            percbinnumb = n4;
            }
        }
        if (fixed.Contains("V0M") && lLow == 40. && lHigh == 50.){
            if (type.Contains("Omega")){
            percentile = pOmega4;
            percbinnumb = npOmega4;
            } else{
            percentile = p5;
            percbinnumb = n5;
            }
        }

        Int_t colors[] = {kRed + 2, kRed + 1, kRed - 4, kOrange + 1, kYellow + 1, kGreen + 1, kTeal - 6, kAzure + 8, kBlue - 4, kBlue + 2};
        TH1F *spectra[percbinnumb], *spectrasys[percbinnumb];

        //---------------------
        TFile* file = TFile::Open(Form("../correctedspectra/CorrSpectra-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", type.Data() , var.Data(), fixed.Data(), lLow, lHigh));
        //

        Double_t lParticleMass = 1.32171; //Default
        if(type.Contains("Omega")) lParticleMass = 1.67245;
        if(type.Contains("Lambda")) lParticleMass = 1.115683;
        if(type.Contains("K0Short")) lParticleMass = 0.497611;

        for (int i = 0; i<percbinnumb; i++){

            double multmin, multmax, eemin, eemax;
            //
            if (var.Contains("SPD")){
                multmin = percentile[i];
                multmax = percentile[i+1];
                eemin = lLow;
                eemax = lHigh;
            }
            //
            if (var.Contains("V0M")){
                eemin = percentile[i];
                eemax = percentile[i+1];
                multmin = lLow;
                multmax = lHigh;
            }

            spectra[i] = (TH1F *)file->Get(Form("FinalSpectra/PtSpectrumCorrStat_%s_%.0f_%.0f_%s_%.0f_%.0f", "SPDClusters", multmin, multmax, "V0M", eemin, eemax));
            spectra[i]->SetName(Form("spectra%i",i));

            spectrasys[i] = (TH1F *)file->Get(Form("FinalSpectra/PtSpectrumCorrSystTot_%s_%.0f_%.0f_%s_%.0f_%.0f", "SPDClusters", multmin, multmax, "V0M", eemin, eemax));
            spectrasys[i]->SetName(Form("spectrasys%i", i));

            beautifyhisto(spectra[i],kFullCircle,colors[i]);
            beautifyhisto(spectrasys[i],kFullCircle,colors[i]);
        }

        Double_t n[percbinnumb];
        for(int nmult = 0; nmult < percbinnumb; nmult++){
            n[nmult] = pow(2,(percbinnumb-(nmult+1)));
            spectra[nmult]->Sumw2();
            spectra[nmult]->Scale(n[nmult]);
            spectrasys[nmult]->Scale(n[nmult]);
        }

        TCanvas* c = new TCanvas("c","",800,1300);
        beautifycanvas(c,0,2);
        c->SetLogy();
        TPad *pad1 = new TPad("pad1","pad1",0,0.27,1,1);
        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.27);
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
       // pad1->Draw();
       // pad2->Draw();
        pad1->SetTicky();
        pad1->SetTickx();
        //
       // pad1->cd();
        pad1->SetLogy();
        pad1->SetFillColor(kWhite);
        //
        TLatex *xlabel = new TLatex();
        xlabel->SetTextFont(42);
        xlabel-> SetNDC();
        xlabel-> SetTextColor(1);
        xlabel-> SetTextSize(0.07);
        xlabel-> SetTextAlign(22);
        xlabel-> SetTextAngle(0);
        //
        TH1D* h1 = new TH1D("h1","",100,0.,10.);
        h1->GetYaxis()->SetRangeUser(1.E-10,10000);
        h1->GetYaxis()->SetLabelSize(0.04);
        h1->GetXaxis()->SetRangeUser(0.,10.);
        h1->GetYaxis()->SetTitle("1/#it{N}_{ev}  d^{2}#it{N}/(d#it{p}_{T}d#it{y}) (GeV/#it{c})^{-1}");
        h1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        h1->GetYaxis()->SetTitleOffset(1.4);
        h1->GetYaxis()->SetTitleSize(0.05);
        h1->SetStats(0);
        h1->Draw();

        //
        TLegend* leg = new TLegend(0.2,0.15,0.45,0.45);
        leg->SetTextSize(0.025);
        leg->SetTextFont(42);
        leg->SetBorderSize(0);
        leg->SetHeader(Form("%s fixed [%.0f-%.0f]%s", fixed.Data(),lLow,lHigh,"%"));
        for (int i = 0; i<percbinnumb; i++){
            spectra[i]->SetTitle(Form("%s [%.0f-%.0f]%s ( #times 2^{%i} )",var.Data(),percentile[i],percentile[i+1],"%",(int)(percbinnumb-(i+1))));
            spectra[i]->Draw("SAME EP");

            spectrasys[i]->SetFillStyle(3000);
            spectra[i]->Draw("SAME E P");
            spectrasys[i]->Draw("SAME E2 P");
            spectra[i]->Draw("SAME");
            leg->AddEntry(spectra[i],"","LEP");
        }
        leg->Draw("SAME");
        if (!type.Contains("K0Short")) xlabel-> DrawLatex(0.85, 0.85,Form("#%s",type.Data()));
        else xlabel-> DrawLatex(0.85, 0.85,Form("K^{0}_{S}"));

        c->SaveAs(Form("Spectra%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", type.Data() , var.Data(), fixed.Data(), lLow, lHigh));
        c->SaveAs(Form("Spectra%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.png", type.Data() , var.Data(), fixed.Data(), lLow, lHigh));
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

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t par2min , Double_t par2max ,Double_t par3min, Double_t par3max, Double_t n, Double_t C, Double_t norm)
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