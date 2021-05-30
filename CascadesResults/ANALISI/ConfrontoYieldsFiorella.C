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

void ConfrontoYieldsFiorella(){

    double Ymio[] = { 

        0.126981, 0.0930338, 0.0750684, 0.0605956, 0.0540262, 0.0442288, 0.0334032, 0.0261914, 0.0160898, 0.00586512

        
       /* 0.12695,
0.09314,
0.075128,
0.060503,
0.0540761,
0.0443123,
0.0334782,
0.026311,
0.0163416,
0.00637379*/};
        //0.1269,  0.0931319,  0.0751731,  0.0605998,  0.054119,  0.0443682, 0.0335511, 0.026411, 0.0164182, 0.00644213 };
    double Statmio[] = {0.00322639,0.00132684,
0.00118275,
0.00104529,
0.00099384,
0.000723653,
0.000592889,
0.000526717,
0.000348851,
0.000314879};
    double PosSystmio[] = {0.00894824,
0.00631219,
0.0051502,
0.00430984,
0.003736,
0.00268908,
0.00196119,
0.00175497,
0.000841647,
0.000356619};
    double NegSystmio[] = {0.00811975,
0.00640544,
0.00530662,
0.00411615,
0.0039349,
0.00271134,
0.00217207,
0.00147264,
0.00111481,
0.000419717
};
    double Yfior[] = {0.1284456594, 0.0924905535, 0.0750700817,0.0609966459, 0.0550614327, 0.0447087709, 0.0333413242, 0.0261796735, 0.0158834917, 0.0061615164};
        //FROM ANALYSIS NOTE 0.12241,0.09185,0.07389,0.06017,0.05439, 0.04365, 0.03239, 0.02598, 0.01606, 0.00723};
    double Statfior[] = {0.00252,0.00121,0.00101,0.00083,0.00078,0.00055,0.00047,0.00045,0.00028,0.00018};
    double PosSystfior[] = {0.00808,0.00608,0.00503,0.00419,0.00404,0.00338,0.00264,0.00234,0.00145,0.00084};
    double NegSystfior[] = {0.00828,0.00604,0.00518,0.00423,0.00425,0.00329,0.00268,0.00216,0.00145,0.00086};
    Double_t dNchV0[] = {26.02, 20.02, 16.17, 13.77, 12.04, 10.02, 7.95, 6.32, 4.50, 2.55};
    Double_t PosSystdNchV0[] = { +0.35, +0.27, +0.22, +0.19, +0.17, +0.14, +0.11, +0.09, +0.07, +0.04 };
    Double_t NegSystdNchV0[] = { 0.29, 0.22, 0.18, 0.16, 0.14, 0.11, 0.09, 0.07, 0.05, +0.03 };
    Double_t StatdNchV0[10] = {0.};
   

    double errmio[10] = {0};
    double errfior[10] = {0};

    TGraphAsymmErrors *systmio = new TGraphAsymmErrors(10,dNchV0,Ymio,PosSystdNchV0,NegSystdNchV0,PosSystmio,NegSystmio);
    TGraphAsymmErrors *systfior = new TGraphAsymmErrors(10,dNchV0,Yfior,PosSystdNchV0,NegSystdNchV0,PosSystfior,NegSystfior);
    TGraphErrors *statmio = new TGraphErrors(10,dNchV0,Ymio,StatdNchV0,Statmio);
    TGraphErrors *statfior = new TGraphErrors(10,dNchV0,Yfior,StatdNchV0,Statfior);

    statmio->SetLineColor(kRed);
    statmio->SetMarkerColor(kRed);
    statmio->SetMarkerStyle(20);
    statmio->SetMarkerSize(1.7);
    systmio->GetXaxis()->SetTitle("<dN_{ch}/d#eta>_{|#eta|<0.5}");
    systmio->GetYaxis()->SetTitle("<dN/dy>");
    systmio->GetYaxis()->SetTitleSize(0.05);
    systmio->GetYaxis()->SetTitleOffset(1.1);
    systmio->GetXaxis()->SetTitleSize(0.05);
    systmio->GetXaxis()->SetTitleOffset(0.9);

    statfior->SetLineColor(kBlue);
    statfior->SetMarkerColor(kBlue);
    statfior->SetMarkerStyle(22);
    statfior->SetMarkerSize(1.7);

    systmio->SetLineColor(kRed+1);
    //systmio->SetMarkerColor();
    //systmio->SetMarkerStyle(0);
    //systmio->SetMarkerSize(0);
    systmio->SetFillStyle(3001);
    systmio->SetFillColor(kRed);

    systfior->SetLineColor(kBlue);
    //systmio->SetMarkerColor();
    //systfior->SetMarkerStyle(0);
    //systfior->SetMarkerSize(0);
    systfior->SetFillStyle(3000);
    systfior->SetFillColor(kBlue);
  
    TLegend* l = new TLegend(0.2,0.6,0.5,0.89);
    l->SetTextSize(0.022);
    l->SetBorderSize(0);     
    l->AddEntry(statmio, "THIS","LP");
    l->AddEntry(statfior, "published pp@13TeV","LP");


    TCanvas* c = new TCanvas("c","",1200,1100);
    c->SetRightMargin(0.09);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetGridy();
    c->SetGridx();
    systmio->SetTitle("");
    systmio->GetYaxis()->SetRangeUser(0,0.15);
    systmio->GetYaxis()->SetTitle("#LT dN/dy #GT");
    systmio->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    systmio->Draw("A2");
    statmio->Draw("SAME EP");
    systfior->Draw("SAME 2");
    statfior->Draw("SAME EP");
    l->Draw("SAME");
    c->SaveAs("images/ConfrontoYields.png");



}