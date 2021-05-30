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

TH1F* ReturnMaxDev(TString inputfile, TString fWhichSystVar = "V0Radius", TString fWhichEstimator = "V0M",Double_t lLoMult = 0., 
    Double_t lHiMult = 100., 
    Double_t lLoEE = 0., 
    Double_t lHiEE = 100.);
double PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin);
double ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr);

void StudySystematicsonYieldsFullXi(
    TString fWhichEstimator = "V0M", 
    TString fWhichParticle = "Xi", 
    Double_t lLoMult = 0., 
    Double_t lHiMult = 100., 
    Double_t lLoEE = 0., 
    Double_t lHiEE = 30.){

    TString inputfilename = Form("Systematics-Yield-%s-%s-sel-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",fWhichParticle.Data(),fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE);
    TH1F* hV0Radius = ReturnMaxDev(inputfilename,"V0Radius", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hCascRadius = ReturnMaxDev(inputfilename,"CascRadius", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hNegToPV = ReturnMaxDev(inputfilename,"DCANegToPV", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hPosToPV = ReturnMaxDev(inputfilename,"DCAPosToPV", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hV0ToPV = ReturnMaxDev(inputfilename,"DCAV0ToPV", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hBachToPV = ReturnMaxDev(inputfilename,"DCABachToPV", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hV0Daught = ReturnMaxDev(inputfilename,"DCAV0Daughters", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hCascDaught = ReturnMaxDev(inputfilename,"DCACascDaughters", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hV0CosPA = ReturnMaxDev(inputfilename,"V0CosPA", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hCascCosPA = ReturnMaxDev(inputfilename,"CascCosPA", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hV0Mass = ReturnMaxDev(inputfilename,"V0Mass", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hCompetingSpecies = 0x0;
    TH1F* hPLT = ReturnMaxDev(inputfilename,"ProperLifetime", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hTPCNSigmas = ReturnMaxDev(inputfilename,"TPCPIDNSigmas", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hTPCNClusters = ReturnMaxDev(inputfilename,"TPCNClusters", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE);
    TH1F* hSigExtBinCount = ReturnMaxDev(inputfilename,"SigExtBinCount", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE); 
    TH1F* hGeantFlukaCorrection = 0;// ReturnMaxDev(inputfilename,"GeantFlukaCorrection", fWhichEstimator, lLoMult, lHiMult, lLoEE, lHiEE); 


    //Topological contribution systematics
    TH1F* hSystTopological = (TH1F*)hV0Radius->Clone("hSystTopological");
    hSystTopological->Reset();

    for (int i = 1; i <= hV0Radius->GetNbinsX(); i++){
        Double_t V0Radius = hV0Radius->GetBinContent(i);
        Double_t CascRadius = hCascRadius->GetBinContent(i);
        Double_t NegToPV = hNegToPV->GetBinContent(i);
        Double_t PosToPV= hPosToPV->GetBinContent(i);
        Double_t BachToPV = hBachToPV->GetBinContent(i);
        Double_t V0Daught = hV0Daught->GetBinContent(i);
        Double_t V0CosPA = hV0CosPA->GetBinContent(i);
        Double_t V0ToPV = hV0ToPV->GetBinContent(i);
        Double_t CascDaught = hCascDaught->GetBinContent(i);
        Double_t CascCosPA = hCascCosPA->GetBinContent(i);
        
        hSystTopological->SetBinContent(i, 
        TMath::Sqrt( V0Radius*V0Radius + NegToPV*NegToPV + PosToPV*PosToPV + BachToPV*BachToPV +
                    V0Daught*V0Daught + CascDaught*CascDaught + V0CosPA*V0CosPA + CascCosPA*CascCosPA +
                    V0CosPA*V0CosPA + V0ToPV*V0ToPV )    
        );
    }

    //Other selections contribution systematics
    TH1F* hSystOthers = (TH1F*)hV0Radius->Clone("hSystOthers");
    hSystOthers->Reset();

    for (int i = 1; i<= hV0Radius->GetNbinsX(); i++){
        Double_t V0Mass = hV0Mass->GetBinContent(i);
        Double_t PLT = hPLT->GetBinContent(i);
        Double_t CompetingSpecies = 0;
        Double_t TPCNClusters = hTPCNClusters->GetBinContent(i);
        Double_t TPCNSigmas = hTPCNSigmas->GetBinContent(i);
        Double_t SigExtBinCount = hSigExtBinCount->GetBinContent(i);
        Double_t GeantFluka = 0;//hGeantFlukaCorrection->GetBinContent(i);

        hSystOthers->SetBinContent(i, 
        TMath::Sqrt(V0Mass*V0Mass  + CompetingSpecies*CompetingSpecies + TPCNClusters*TPCNClusters +
                    TPCNSigmas*TPCNSigmas + SigExtBinCount*SigExtBinCount + PLT*PLT + GeantFluka*GeantFluka ) 
        );
    }

    //Total systematics
    TH1F* hSystTot = (TH1F*)hV0Radius->Clone("hSystTot");  
    hSystTot->Reset();
    for (int i = 1; i<= hV0Radius->GetNbinsX(); i++){
        hSystTot->SetBinContent(i, 
            TMath::Sqrt(
            hSystOthers->GetBinContent(i)*hSystOthers->GetBinContent(i) +
            hSystTopological->GetBinContent(i)*hSystTopological->GetBinContent(i) )
            );        
    }


    double low = 0., high = 100;
    if (fWhichEstimator.Contains("V0M")) {low = lLoEE; high = lHiEE;}
    if (fWhichEstimator.Contains("ZDC")) {low = lLoMult; high = lHiMult;}
    TString fWhichOtherEstimator = "ZDC";
    if (fWhichEstimator.Contains("ZDC")) fWhichOtherEstimator = "V0M";
    TLatex *xlabel = new TLatex();
 	xlabel->SetTextFont(42);
    xlabel-> SetNDC();
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.03);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    
    //Topological Displayed
    TCanvas * Top = new TCanvas("Top","",1000,1000);
	Top->SetRightMargin(0.09);
	Top->SetLeftMargin(0.2);
	Top->SetBottomMargin(0.15);
    TLegend* l = new TLegend (0.5,0.55,0.89,0.89);
    l->SetBorderSize(0);
    l->AddEntry(hV0Radius,"V0 Radius","L");
    l->AddEntry(hCascRadius,"Casc Radius","L");
    l->AddEntry(hNegToPV,"DCA Neg To PV","L");
    l->AddEntry(hPosToPV,"DCA Pos To PV","L");
    l->AddEntry(hBachToPV,"DCA Bach To PV","L");
    l->AddEntry(hV0Daught,"DCA V0 Daughters","L");
    l->AddEntry(hCascDaught,"DCA Casc Daughters","L");
    l->AddEntry(hV0CosPA,"V0 Cosine of PA","L");
    l->AddEntry(hCascCosPA,"Casc Cosine of PA","L");
    l->AddEntry(hV0ToPV,"V0 To PV","L");
    l->SetTextSize(0.03);
    hV0Radius->SetLineColor(kRed+1);
    hCascRadius->SetLineColor(kRed);
    hNegToPV->SetLineColor(kOrange+7);
    hPosToPV->SetLineColor(kOrange);
    hBachToPV->SetLineColor(kYellow);
    hV0Daught->SetLineColor(kSpring+10);
    hCascDaught->SetLineColor(kSpring);
    hV0CosPA->SetLineColor(kAzure+8);
    hCascCosPA->SetLineColor(kAzure);
    hV0ToPV->SetLineColor(kBlue+3);
    hV0Radius->SetMarkerStyle(1);
    hNegToPV->SetMarkerStyle(1);
    hCascRadius->SetMarkerStyle(1);
    hPosToPV->SetMarkerStyle(1);
    hBachToPV->SetMarkerStyle(1);
    hV0Daught->SetMarkerStyle(1);
    hCascDaught->SetMarkerStyle(1);
    hV0CosPA->SetMarkerStyle(1);
    hCascCosPA->SetMarkerStyle(1);
    hV0ToPV->SetMarkerStyle(1);
    hV0Radius->SetLineWidth(2);
    hNegToPV->SetLineWidth(2);
    hCascRadius->SetLineWidth(2);
    hPosToPV->SetLineWidth(2);
    hBachToPV->SetLineWidth(2);
    hV0Daught->SetLineWidth(2);
    hCascDaught->SetLineWidth(2);
    hV0CosPA->SetLineWidth(2);
    hCascCosPA->SetLineWidth(2);
    hV0ToPV->SetLineWidth(2);
    hV0Radius->SetTitle(Form("Contribution of Topological Variables %s",fWhichParticle.Data()));
    //hV0Radius->GetYaxis()->SetTitleOffset(1.2);
    hV0Radius->SetTitle("");
    hV0Radius->Draw("HIST ");
    hCascRadius->Draw("HIST SAME ");
    hNegToPV->Draw("HIST SAME ");
    hPosToPV->Draw("HIST SAME ");
    hBachToPV->Draw("HIST SAME ");
    hV0Daught->Draw("HIST SAME ");
    hCascDaught->Draw("HIST SAME ");
    hV0CosPA->Draw("HIST SAME ");
    hCascCosPA->Draw("HIST SAME ");
    hV0ToPV->Draw("HIST SAME ");
    l->Draw("SAME");
    xlabel-> DrawLatex(0.35, 0.82, Form("%s fixed [%.0f,%.0f]",fWhichOtherEstimator.Data(),low,high));

    Top->SaveAs(Form("images/yields/%s/TopologicalYieldsSystematics%s-%s_V0-%03.0f_%03.0f-ZDC-%03.0f_%03.0f.png",fWhichEstimator.Data(),fWhichParticle.Data(),fWhichEstimator.Data(),lLoMult, lHiMult, lLoEE, lHiEE));

    //Selection Displayed
    TCanvas * Sel = new TCanvas("Sel","",1000,1000);
    Sel->SetRightMargin(0.09);
    Sel->SetLeftMargin(0.2);
    Sel->SetBottomMargin(0.15);
    TLegend* l1 = new TLegend (0.45,0.55,0.85,0.89);
    l1->SetBorderSize(0);
    l1->AddEntry(hPLT,"Proper Life Time","L");
    l1->AddEntry(hV0Mass,"V0 Mass","L");
    l1->AddEntry(hTPCNClusters,"TPC N of Clusters","L");
    l1->AddEntry(hTPCNSigmas,"TPC N Sigmas","L");
    l1->AddEntry(hSigExtBinCount,"Sigmas for Sgn Extraction","L");
    //l1->AddEntry(hGeantFlukaCorrection,"Geant-Fluka Correction","L");
    l1->SetTextSize(0.03);
    hPLT->SetLineColor(kRed+1);
    hV0Mass->SetLineColor(kOrange+1);
    hSigExtBinCount->SetLineColor(kSpring-1);
    hTPCNClusters->SetLineColor(kAzure-4);
    hTPCNSigmas->SetLineColor(kBlue+2);
    hPLT->SetMarkerStyle(1);
    hPLT->SetLineWidth(2);
    hV0Mass->SetMarkerStyle(1);
    hV0Mass->SetLineWidth(2);
    hTPCNClusters->SetMarkerStyle(1);
    hTPCNClusters->SetLineWidth(2);
    hTPCNSigmas->SetMarkerStyle(1);
    hTPCNSigmas->SetLineWidth(2);
    hSigExtBinCount->SetMarkerStyle(1);
    hSigExtBinCount->SetLineWidth(2);
    //hGeantFlukaCorrection->SetLineWidth(2);
    hPLT->SetTitle("Contribution of Selection Variables");
    //hPLT->GetYaxis()->SetTitleOffset(1.2);
    hPLT->SetTitle("");
    hPLT->Draw("HIST ");
    hV0Mass->Draw("HIST SAME ");
    hTPCNClusters->Draw("HIST SAME ");
    hTPCNSigmas->Draw("HIST SAME ");
    hSigExtBinCount->Draw("HIST SAME ");
    //hGeantFlukaCorrection->Draw("HIST SAME");
    l1->Draw("SAME");
    xlabel-> DrawLatex(0.35, 0.82, Form("%s fixed [%.0f,%.0f]",fWhichOtherEstimator.Data(),low,high));

    Sel->SaveAs(Form("images/yields/%s/%s-YieldsSelectionSystematics%s_V0-%03.0f_%03.0f-ZDC-%03.0f_%03.0f.png",fWhichEstimator.Data(),fWhichParticle.Data(),fWhichEstimator.Data(),lLoMult, lHiMult, lLoEE, lHiEE));

    //Contributions Displayed
    TCanvas* cn = new TCanvas("cn","",900,1000);
    cn->SetRightMargin(0.09);
    cn->SetLeftMargin(0.2);
    cn->SetBottomMargin(0.15);
     
    TLegend* legend = new TLegend(0.25,0.55,0.65,0.89);
    legend->SetBorderSize(0);
    legend->AddEntry(hSystTopological,"Topological systematics","L");
    legend->AddEntry(hSystOthers,"Selection cuts systematics","L");
    legend->AddEntry(hSystTot,"Total systematics","L");
    legend->SetTextSize(0.027);

    hSystTot->SetLineWidth(2);
    hSystTot->SetLineColor(kBlack);
    hSystOthers->SetLineWidth(2);
    hSystOthers->SetLineColor(kRed);
    hSystTopological->SetLineWidth(2);
    hSystTopological->SetLineColor(kBlue);
  
    hSystTot->SetTitle("");
    hSystTot->GetYaxis()->SetRangeUser(-0.005,.1);
    //hSystTot->GetYaxis()->SetTitleOffset(1.);
    //hSystTot->SetTitle(Form("Systematics contributions for %s",fWhichParticle.Data()));
    //hSystTot->GetYaxis()->SetTitleOffset(1.2);
    hSystTot->Draw();
    hSystOthers->Draw("SAME");
    hSystTopological->Draw("SAME");
    legend->Draw("SAME");  
    xlabel-> DrawLatex(0.4, 0.52, Form("%s fixed [%.0f,%.0f]",fWhichOtherEstimator.Data(),low,high));


    cn->SaveAs(Form("images/yields/%s/%s-YieldsSystematics-Displayed-%s_V0-%03.0f_%03.0f-ZDC-%03.0f_%03.0f.png",fWhichEstimator.Data(),fWhichParticle.Data(),fWhichEstimator.Data(),lLoMult, lHiMult, lLoEE, lHiEE));


    //Total syst canvas shaded 
    TCanvas *cfinal = new TCanvas("cfinal","Total Systematics",900,600);
    cfinal->SetFillColor(kWhite);
    cfinal->SetLeftMargin(0.17);
    cfinal->SetRightMargin(0.17);
    cfinal->SetBottomMargin(0.17);
    cfinal->cd();
    TH1F* hSystTotNeg = (TH1F*)hSystTot->Clone("hSystTotNeg");
    hSystTotNeg->Reset();
    for (int i = 1;i<=hSystTot->GetNbinsX();i++){
        hSystTotNeg->SetBinError(i,(hSystTot->GetBinContent(i)));
        hSystTotNeg->SetBinContent(i, 0.);
    }
    hSystTotNeg->SetFillColor(18);
    hSystTotNeg->SetMarkerStyle(20);
    //hSystTotNeg->GetYaxis()->SetTitle("max. rel. deviation (%)");
    hSystTotNeg->GetYaxis()->SetRangeUser(0.-0.1,0.+0.1);
    hSystTotNeg->SetTitle("");
    hSystTotNeg->Draw("E2");
    xlabel-> DrawLatex(0.5, 0.82, Form("%s fixed [%.0f,%.0f]",fWhichOtherEstimator.Data(),low,high));

    cfinal->SaveAs(Form("images/yields/%s/FinalYieldsSystContribPlot%s%s-V0M_%03.0f_%03.0f-ZDC-%03.0f_%03.0f.png",fWhichEstimator.Data(),fWhichParticle.Data(),fWhichEstimator.Data(),lLoMult, lHiMult, lLoEE, lHiEE));
    
    TString outputfilename = Form("CorrSystematics-Yield-%s-%s-sel-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",fWhichParticle.Data(),fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE);
    TFile* Write = new TFile (outputfilename, "RECREATE");
    hSystTopological->Write();
    hSystOthers->Write();
    hSystTot->Write(); 
    hV0Radius->Write();
    hCascRadius->Write();
    hNegToPV->Write();
    hPosToPV->Write();
    hBachToPV->Write();
    hV0Daught->Write();
    hCascDaught->Write();
    hV0CosPA->Write();
    hCascCosPA->Write();
    hV0ToPV->Write();
    hPLT->Write();
    hV0Mass->Write();
    hTPCNClusters->Write();
    hTPCNSigmas->Write();
    hSigExtBinCount->Write();
    //hGeantFlukaCorrection->Write();
}

TH1F* ReturnMaxDev(
    TString inputfile, 
    TString fWhichSystVar = "V0Radius", 
    TString fWhichEstimator = "V0M",
    Double_t lLoMult = 0., 
    Double_t lHiMult = 100., 
    Double_t lLoEE = 0., 
    Double_t lHiEE = 100.
    ){
        double low = 0., high = 100;
        if (fWhichEstimator.Contains("V0M")) {low = lLoEE; high = lHiEE;}
        if (fWhichEstimator.Contains("ZDC")) {low = lLoMult; high = lHiMult;}
        TString fWhichOtherEstimator = "ZDC";
        if (fWhichEstimator.Contains("ZDC")) fWhichOtherEstimator = "V0M";
       

    Double_t percentileV0[] = {0,15,40,100};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    Double_t percentileZDC[] = {0,40,70,100};
    const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;
    Double_t percentileMB[] = {0,100};
    const Long_t nbinMB = sizeof(percentileMB)/sizeof(Double_t) - 1;
    Double_t * percentile;
    int tempbinnumber = 0;
    if (fWhichEstimator.Contains("V0M")) {
        tempbinnumber = nbinV0;
        percentile = percentileV0;
    } 
    else if (fWhichEstimator.Contains("ZDC")) {
        tempbinnumber = nbinZDC;
        percentile = percentileZDC;
    } 
    else if (fWhichEstimator.Contains("MB")) {
        tempbinnumber = nbinMB;
        percentile = percentileMB;
       // low = lLoEE; 
       // high = lHiEE;
        low = lLoMult; 
        high = lHiMult;
    } 
    else {cout << "No valid name for estimator... its V0M or ZDC" << endl; }
    const int binnumber = tempbinnumber;
    
    // Open File
    TFile* inputfilename = TFile::Open(inputfile);
    TFile* MBinputfilename = TFile::Open(Form("Systematics-Yield-%s-%s-sel-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","MB",0.,100.,0.,100.));

    TH1F* hMBCutVar[4], *hCutVar[4], *hMaxDev, *hCutVarClone[4];
    for (int nfile = 0; nfile < 4 ; nfile++){
        hMBCutVar[nfile] = (TH1F*)MBinputfilename->Get(Form("MBVariationCuts%s/hYield-%s-MB-%i",fWhichSystVar.Data(),fWhichSystVar.Data(),nfile+1));
        hCutVar[nfile] = (TH1F*)inputfilename->Get(Form("%sVariationCuts%s/hYield-%s-%s-%i",fWhichEstimator.Data(),fWhichSystVar.Data(),fWhichSystVar.Data(),fWhichEstimator.Data(),nfile+1));
        hCutVarClone[nfile] = (TH1F* )hCutVar[nfile]->Clone(Form("hYield-%s-%s-%i",fWhichSystVar.Data(),fWhichEstimator.Data(),nfile+1));
        hCutVarClone[nfile]->Reset();
        for (int bin = 1; bin <= hCutVarClone[nfile]->GetNbinsX(); bin++){
            hCutVarClone[nfile]->SetBinContent(bin,hCutVar[nfile]->GetBinContent(bin)/hMBCutVar[nfile]->GetBinContent(1));
            hCutVarClone[nfile]->SetBinError(bin,ErrorInRatio(hCutVar[nfile]->GetBinContent(bin),hCutVar[nfile]->GetBinError(bin),
                    hMBCutVar[nfile]->GetBinContent(1),hMBCutVar[nfile]->GetBinError(1)));

        }
    }

    TCanvas* c1 = new TCanvas("c","",1000,700);
    c1->SetGridy();
    c1->SetGridx();
    c1->SetRightMargin(0.09);
    c1->SetLeftMargin(0.17);
    c1->SetBottomMargin(0.15);
    hCutVarClone[0]->GetYaxis()->SetRangeUser(0.95,1.05);
    hCutVarClone[0]->GetYaxis()->SetTitle("R = #frac{[dN/dy^{syst-cut}/dN/dy^{def-cut}]_{sel}}{[dN/dy^{syst-cut}/dN/dy^{def-cut}]_{MB}}");
    hCutVarClone[0]->GetYaxis()->SetTitleOffset(2.4);
    hCutVarClone[0]->GetXaxis()->SetTitleOffset(1.3);
    hCutVarClone[0]->GetYaxis()->SetTitleSize(0.03);
    hCutVarClone[0]->GetXaxis()->SetTitleSize(0.03);
    for (int nfile = 0; nfile < 4 ; nfile++){
        //hCutVarClone[nfile]->SetLineStyle(7);
        //hCutVarClone[nfile]->SetLineWidth(2);
        hCutVarClone[nfile]->SetMarkerStyle(24);
        hCutVarClone[nfile]->SetMarkerSize(1.5);
        hCutVarClone[nfile]->Draw("SAME LEP");
    }
    TLegend* legend = new TLegend (0.7,.7,0.91,0.89);  
    legend->SetBorderSize(0);
    legend->AddEntry(hCutVarClone[0],"very loose","LEP");
    legend->AddEntry(hCutVarClone[1],"loose","LEP");
    legend->AddEntry(hCutVarClone[2],"tight","LEP");
    legend->AddEntry(hCutVarClone[3],"very tight","LEP");
    legend->SetTextSize(0.03);
    legend->Draw("SAME");
    TLatex *xlabel = new TLatex();
 	xlabel->SetTextFont(42);
    xlabel-> SetNDC();
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.03);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);
    xlabel-> DrawLatex(0.5, 0.82, Form("%s fixed [%.0f,%.0f]",fWhichOtherEstimator.Data(),low,high));


    c1->SaveAs(Form("images/yields/%s/Var%s%s_V0-%03.0f_%03.0f-ZDC-%03.0f_%03.0f.png",fWhichEstimator.Data(),fWhichSystVar.Data(),fWhichEstimator.Data(),lLoMult, lHiMult, lLoEE, lHiEE));
   

    //
    hMaxDev = (TH1F*)hCutVar[0]->Clone(Form("hMaxDev%s",fWhichSystVar.Data()));
    hMaxDev->Reset();  
   
    for (int nmult = 0; nmult < binnumber; nmult ++){
        double maxdev = 0;
        double valuexfile[4];
        for (int nfile = 0; nfile < 4 ; nfile++){
            //valuexfile[nfile] = TMath::Abs((hCutVar[nfile]->GetBinContent(nmult+1) - hMBCutVar[nfile]->GetBinContent(1)));
            valuexfile[nfile] = PassRogerBarlowCriterion(1,hCutVarClone[nfile],nmult+1);           
        }
        maxdev = valuexfile[0];
        int counter = 0;
        for (int k = 1; k < 4; k++){ 
            maxdev = max(maxdev,valuexfile[k]);
            if (maxdev == valuexfile[k]) counter = k;
        }
        hMaxDev->SetBinContent(nmult+1, maxdev);
        cout << Form("%s   ",fWhichSystVar.Data()) << hMBCutVar[0]->GetBinContent(1) << endl;
        hMaxDev->SetBinError(nmult+1, hCutVarClone[counter]->GetBinError(nmult+1));
        hMaxDev->SetLineColor(kBlack);
        hMaxDev->GetYaxis()->SetTitle("R^{MAX} = #frac{[dN/dy^{syst-cut}/dN/dy^{def-cut}]_{sel}}{[dN/dy^{syst-cut}/dN/dy^{def-cut}]_{MB}}");
        hMaxDev->GetYaxis()->SetRangeUser(0.,0.05);
    }
    hMaxDev->GetYaxis()->SetTitleOffset(2.4);
    hMaxDev->GetXaxis()->SetTitleOffset(1.3);
    hMaxDev->GetYaxis()->SetTitleSize(0.03);
    hMaxDev->GetXaxis()->SetTitleSize(0.03);
    hMaxDev->SetTitle("");
    
    return hMaxDev;

}

double PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin){
  double dev = h->GetBinContent(bin);
  double RBsigma = h->GetBinError(bin);

  if( (TMath::Abs(dev-1)) > RBsigma*nsigmas ) {return TMath::Abs(dev-1);}
  else {return 0.;}
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
