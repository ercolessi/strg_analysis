#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"

#endif
using namespace std;
#include <TString.h>

TCanvas* eventloss(TString which = "hevtlossV0", TString estim = "V0 multiplicity [%]");
TCanvas* signallossv0(TString part = "XiMinus", TString dir = "multsel", TString estim = "V0M", TString text = "Xi^{-}", TString histo = "V0", Int_t* mult =0x0 );
TCanvas* signallosszdc(TString part = "XiMinus", TString dir = "multsel", TString estim = "V0M", TString text = "Xi^{-}", TString histo = "V0", Int_t* mult = 0x0 );

void DrawEfficiencies_MergedBins(){

    Int_t multV0[11] = {0,5,10,15,20,30,40,50,70,100};
    Int_t multZDC[10] = {0,20,40,50,60,70,80,90,100};
    TText *xlabel = new TText();
    xlabel-> SetNDC();
    xlabel-> SetTextFont(42);
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.04);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    

    //Draw Event Loss-------------------------------------------
    /*TCanvas* eevtv0 = eventloss();
    eevtv0->Draw();
    eevtv0->SaveAs("hevtlossV0.png");
    //
    TCanvas* eevtzdc = eventloss("hevtlossZDC", "ZDC [%]");
    eevtzdc->Draw();
    eevtzdc->SaveAs("hevtlossZDC.png");*/
    //
    TCanvas* eevtlom = eventloss("hevtlossV0FixLowmult", "ZDC [%]");
    eevtlom->Draw();
    xlabel-> DrawText(0.415, 0.55, "fixed V0M [70,100]%");
    eevtlom->SaveAs("hevtlossV0FixedLowmult.png");
    //
   /* TCanvas* eevthim = eventloss("hevtlossV0FixHighmult", "ZDC [%]");
    eevthim->Draw();
    xlabel-> DrawText(0.415, 0.55, "fixed V0M [0,30]%");
    eevthim->SaveAs("hevtlossV0FixedHighmult.png");*/
    //
    TCanvas* eevtloe = eventloss("hevtlossV0FixLowEE");
    eevtloe->Draw();
    xlabel-> DrawText(0.415, 0.55, "fixed ZDC [70,100]%");
    eevtloe->SaveAs("hevtlossV0FixedLowEE.png");
    //
    TCanvas* eevthie = eventloss("hevtlossV0FixHighEE");
    eevthie->Draw();
    xlabel-> DrawText(0.415, 0.55, "fixed ZDC [0,30]%");
    eevthie->SaveAs("hevtlossV0FixedHighEE.png"); 
    
    //Signal Loss--------------------------------------------------
    /*TCanvas* epartv0 = signallossv0("XiMinus","multsel","V0","Xi^{-}","V0",multV0);
    epartv0->Draw();
    epartv0->SaveAs("SgnLossXiMinus_multsel.png");
    //
    TCanvas* epartv0p = signallossv0("XiPlus","multsel","V0","Xi^{+}","V0",multV0);
    epartv0p->Draw();
    epartv0p->SaveAs("SgnLossXiPlus_multsel.png");
    //
    TCanvas* epartzdc = signallosszdc("XiMinus","EEsel","ZDC","Xi^{-}","ZDC",multZDC);
    epartzdc->Draw();
    epartzdc->SaveAs("SgnLossXiMinus_EEsel.png");
    //
    TCanvas* epartzdcp = signallosszdc("XiPlus","EEsel","ZDC","Xi^{+}","ZDC",multZDC);
    epartzdcp->Draw();
    epartzdcp->SaveAs("SgnLossXiPlus_EEsel.png");*/
    //
    TCanvas* epartloee = signallossv0("XiMinus","multsel_fixedlowEE","V0","Xi^{-}","V0FixLowEE",multV0);
    epartloee->Draw();
    xlabel-> DrawText(0.715, 0.55, "fixed ZDC [70,100]%");
    epartloee->SaveAs("SgnLossXiMinus_multsel_fixedlowEE.png");
    //
    /*TCanvas* epartloeep = signallossv0("XiPlus","multsel_fixedlowEE","V0","Xi^{+}","V0FixLowEE",multV0);
    epartloeep->Draw();
    xlabel-> DrawText(0.715, 0.55, "fixed ZDC [70,100]%");
    epartloeep->SaveAs("SgnLossXiPlus_multsel_fixedlowEE.png");*/
    //
    TCanvas* eparthiee = signallossv0("XiMinus","multsel_fixedhighEE","V0","Xi^{-}","V0FixHighEE",multV0);
    eparthiee->Draw();
    xlabel-> DrawText(0.715, 0.55, "fixed ZDC [0,30]%");
    eparthiee->SaveAs("SgnLossXiMinus_multsel_fixedhighEE.png");
    //
   /* TCanvas* eparthieep = signallossv0("XiPlus","multsel_fixedhighEE","V0","Xi^{+}","V0FixHighEE",multV0);
    eparthieep->Draw();
    xlabel-> DrawText(0.715, 0.55, "fixed ZDC [0,30]%");
    eparthieep->SaveAs("SgnLossXiPlus_multsel_fixedhighEE.png");*/
    //
    TCanvas* epartlom = signallosszdc("XiMinus","EEsel_fixedlowmult","ZDC","Xi^{-}","ZDCFixLowmult",multZDC);
    epartlom->Draw();
    xlabel-> DrawText(0.715, 0.55, "fixed V0M [70,100]%");
    epartlom->SaveAs("SgnLossXiMinus_EEsel_fixedlowmult.png");
    //
 /*   TCanvas* epartlomp = signallosszdc("XiPlus","EEsel_fixedlowmult","ZDC","Xi^{+}","ZDCFixLowmult",multZDC);
    epartlomp->Draw();
    xlabel-> DrawText(0.715, 0.55, "fixed V0M [70,100]%");
    epartlomp->SaveAs("SgnLossXiPlus_EEsel_fixedlowmult.png");
    //
    TCanvas* eparthim = signallosszdc("XiMinus","EEsel_fixedhighmult","ZDC","Xi^{-}","ZDCFixHighmult",multZDC);
    eparthim->Draw();
    xlabel-> DrawText(0.715, 0.55, "fixed V0M [0,30]%");
    eparthim->SaveAs("SgnLossXiMinus_EEsel_fixedhighmult.png");
    //
   /* TCanvas* eparthimp = signallosszdc("XiPlus","EEsel_fixedhighmult","ZDC","Xi^{+}","ZDCFixHighmult",multZDC);
    eparthimp->Draw();
    xlabel-> DrawText(0.715, 0.55, "fixed V0M [0,30]%");
    eparthimp->SaveAs("SgnLossXiPlus_EEsel_fixedhighmult.png");*/
    //


}

TCanvas* eventloss(TString which = "hevtlossV0", TString estim = "V0 multiplicity [%]") {   
    
    TFile* file = new TFile("EventCountLoss_16d3_mergedbins.root","READ");
    
    TH1F* hEevt = (TH1F*)file->Get(Form("EventLoss/%s",which.Data()));
    hEevt->SetMarkerColor(kBlue);
    hEevt->SetLineColor(kBlue);
    hEevt->GetYaxis()->SetRangeUser(0.6,1.01);
    hEevt->SetStats(0);
    hEevt->GetYaxis()->SetTitle("#varepsilon_{evt} ");
    hEevt->GetXaxis()->SetTitle(Form("%s",estim.Data()));  
    hEevt->SetMarkerStyle(21);
    hEevt->SetMarkerSize(1.4);
    hEevt->SetLineWidth(1);
    hEevt->GetYaxis()->SetTitleOffset(1.);
    hEevt->GetXaxis()->SetTitleOffset(1.);
    hEevt->GetYaxis()->SetTitleSize(0.07);
    hEevt->GetXaxis()->SetTitleSize(0.06);
    
    TCanvas* sR = new TCanvas("sR","",900,800);
    sR->SetRightMargin(0.09);
    sR->SetLeftMargin(0.15);
    sR->SetBottomMargin(0.15);
    
    TText *xlabel = new TText();
    xlabel-> SetNDC();
    xlabel-> SetTextFont(1);
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.04);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    
    hEevt->Draw("HIST L");
    hEevt->Draw("HIST P SAME");
   // xlabel->DrawText(0.815, 0.85, Form("%s",period.Data()));
    return sR;
}

TCanvas* signallossv0(TString part = "XiMinus", TString dir = "multsel", TString estim = "V0M", TString text = "Xi^{-}", TString histo = "V0", Int_t* mult =0x0) {   
    
    TFile* file = new TFile("EventCountLoss_16d3_mergedbins.root","READ");
    TH1F* hsgnloss[9];

    for (int i =0; i<9;i++){
        hsgnloss[i] = (TH1F*)file->Get(Form("SgnLoss/%s/%s/fHistptSel%s_%i-%i_%s",part.Data(),dir.Data(), histo.Data(),mult[i],mult[i+1],part.Data()));
    }
    hsgnloss[0]->SetMarkerColor(kRed+1);
    hsgnloss[1]->SetMarkerColor(kRed-4);
    hsgnloss[2]->SetMarkerColor(kOrange+7);
    hsgnloss[3]->SetMarkerColor(kOrange-3);
    hsgnloss[4]->SetMarkerColor(kYellow+1);
    hsgnloss[5]->SetMarkerColor(kSpring-7);
    hsgnloss[6]->SetMarkerColor(kGreen+2);
    hsgnloss[7]->SetMarkerColor(kAzure+8);
    hsgnloss[8]->SetMarkerColor(kBlue-4);
  //  hsgnloss[9]->SetMarkerColor(kBlue+3);
    hsgnloss[0]->SetLineColor(kRed+1);
    hsgnloss[1]->SetLineColor(kRed-4);
    hsgnloss[2]->SetLineColor(kOrange+7);
    hsgnloss[3]->SetLineColor(kOrange-3);
    hsgnloss[4]->SetLineColor(kYellow+1);
    hsgnloss[5]->SetLineColor(kSpring-7);
    hsgnloss[6]->SetLineColor(kGreen+2);
    hsgnloss[7]->SetLineColor(kAzure+8);
    hsgnloss[8]->SetLineColor(kBlue-4);
   // hsgnloss[9]->SetLineColor(kBlue+3);
    
    TCanvas* sR = new TCanvas("sR","",900,800);

    TText *xlabel = new TText();
    xlabel-> SetNDC();
    xlabel-> SetTextFont(0);
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.03);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    
    TMathText plabel;
    plabel.SetNDC();
    plabel.SetTextFont(1);
    plabel.SetTextColor(1);
    plabel.SetTextSize(0.05);
    plabel.SetTextAlign(22);
    plabel.SetTextAngle(0);
    
    TLegend* l = new TLegend(0.57,0.2,0.87,0.5);
    l->SetTextSize(0.023);
    l->SetHeader(Form("%s selection (percentiles):",estim.Data()));
    l->SetBorderSize(0);

    for (int i = 0; i<9;i++){
        hsgnloss[i]->SetTitle(Form("%i-%i",mult[i],mult[i+1]));
        hsgnloss[i]->GetYaxis()->SetRangeUser(0.7,1.02);
        hsgnloss[i]->GetXaxis()->SetRangeUser(0.6,6.5);
        hsgnloss[i]->GetYaxis()->SetTitle(Form("#varepsilon_{part} (p_{T})  [#%s]", "Xi"));
        hsgnloss[i]->GetYaxis()->SetTitleOffset(1.5);
        hsgnloss[i]->SetStats(0);
        hsgnloss[i]->SetMarkerStyle(20);
        hsgnloss[i]->SetMarkerSize(1.6);
        hsgnloss[i]->Draw("EP SAME");
        l->AddEntry(hsgnloss[i],"","LEP");
        l->Draw("SAME");
    }
    
    hsgnloss[0]->SetTitle("");//Form("#%s", text.Data()));
    hsgnloss[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hsgnloss[0]->GetYaxis()->SetTitleOffset(.9);
    hsgnloss[0]->GetXaxis()->SetTitleOffset(.9);
    hsgnloss[0]->GetXaxis()->SetTitleSize(0.06);
    hsgnloss[0]->GetYaxis()->SetTitleSize(0.06);
    //xlabel-> DrawText(0.8, 0.35, Form("%s", period.Data()));
    //plabel.DrawMathText(0.355, 0.3, text.Data());
    
    sR->SetRightMargin(0.09);
    sR->SetLeftMargin(0.15);
    sR->SetBottomMargin(0.15);
    sR->SetTopMargin(0.1);

    return sR;   
}

TCanvas* signallosszdc(TString part = "XiMinus", TString dir = "multsel", TString estim = "V0M", TString text = "Xi^{-}", TString histo = "V0", Int_t* mult =0x0) {   
    
    TFile* file = new TFile("EventCountLoss_16d3_mergedbins.root","READ");
    TH1F* hsgnloss[8];

    for (int i =0; i<8;i++){
        hsgnloss[i] = (TH1F*)file->Get(Form("SgnLoss/%s/%s/fHistptSel%s_%i-%i_%s",part.Data(),dir.Data(), histo.Data(),mult[i],mult[i+1],part.Data()));
    }
    hsgnloss[0]->SetMarkerColor(kRed+1);
    hsgnloss[1]->SetMarkerColor(kRed-4);
    hsgnloss[2]->SetMarkerColor(kOrange+7);
    hsgnloss[3]->SetMarkerColor(kOrange-3);
    hsgnloss[4]->SetMarkerColor(kYellow+1);
    hsgnloss[5]->SetMarkerColor(kSpring-7);
    hsgnloss[6]->SetMarkerColor(kGreen+2);
    hsgnloss[7]->SetMarkerColor(kAzure+8);
  //  hsgnloss[8]->SetMarkerColor(kBlue-4);
    
    hsgnloss[0]->SetLineColor(kRed+1);
    hsgnloss[1]->SetLineColor(kRed-4);
    hsgnloss[2]->SetLineColor(kOrange+7);
    hsgnloss[3]->SetLineColor(kOrange-3);
    hsgnloss[4]->SetLineColor(kYellow+1);
    hsgnloss[5]->SetLineColor(kSpring-7);
    hsgnloss[6]->SetLineColor(kGreen+2);
    hsgnloss[7]->SetLineColor(kAzure+8);
   // hsgnloss[8]->SetLineColor(kBlue-4);
    
    TCanvas* sR = new TCanvas("sR","",900,800);

    TText *xlabel = new TText();
    xlabel-> SetNDC();
    xlabel-> SetTextFont(0);
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.03);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    
    TMathText plabel;
    plabel.SetNDC();
    plabel.SetTextFont(1);
    plabel.SetTextColor(1);
    plabel.SetTextSize(0.05);
    plabel.SetTextAlign(22);
    plabel.SetTextAngle(0);
    
    TLegend* l = new TLegend(0.57,0.2,0.87,0.5);
    l->SetTextSize(0.023);
    l->SetHeader(Form("%s selection (percentiles):",estim.Data()));
    l->SetBorderSize(0);
    
    for (int i = 0; i<8;i++){
        hsgnloss[i]->SetTitle(Form("%i-%i",mult[i],mult[i+1]));
        hsgnloss[i]->GetYaxis()->SetRangeUser(0.7,1.02);
        hsgnloss[i]->GetXaxis()->SetRangeUser(0.6,6.5);
        hsgnloss[i]->GetYaxis()->SetTitle(Form("#varepsilon_{part} (p_{T})  [#%s]", "Xi"));
        hsgnloss[i]->GetYaxis()->SetTitleOffset(1.5);
        hsgnloss[i]->SetStats(0);
        hsgnloss[i]->SetMarkerStyle(20);
        hsgnloss[i]->SetMarkerSize(1.6);
        hsgnloss[i]->Draw("EP SAME");
        l->AddEntry(hsgnloss[i],"","LEP");
        l->Draw("SAME");
    }
    
    hsgnloss[0]->SetTitle("");//Form("#%s", text.Data()));
    hsgnloss[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hsgnloss[0]->GetXaxis()->SetTitleOffset(.9);
    hsgnloss[0]->GetYaxis()->SetTitleOffset(.9);
    hsgnloss[0]->GetXaxis()->SetTitleSize(0.06);
    hsgnloss[0]->GetYaxis()->SetTitleSize(0.06);
    //xlabel-> DrawText(0.8, 0.35, Form("%s", period.Data()));
    //plabel.DrawMathText(0.355, 0.3, text.Data());
    
    sR->SetRightMargin(0.09);
    sR->SetLeftMargin(0.15);
    sR->SetBottomMargin(0.15);
    sR->SetTopMargin(0.1);

    return sR;   
}
