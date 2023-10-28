#include <iostream>
#include <fstream>
//#define max(a,b) (a>b ? a : b)

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 );
double PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin);

TH1F* makeSystPlotsV0s(
  TString filename,
  TString lType = "Xi",
	Double_t lMultBoundLo = 0.0,
	Double_t lMultBoundHi = 100.0,
	Double_t lEEBoundLo = 0.0,
	Double_t lEEBoundHi = 100.0,
	TString lWhichSystVar = "V0Radius" );


//--------------------------------------------------------------
//------------------- MAIN FUNCTION ----------------------------
//--------------------------------------------------------------
void DoCutVariationSystematics(
  TString Type = "Lambda",
  Double_t LowMult = 0.0,
  Double_t HighMult = 100,
  Double_t LowEE = 0.0,
  Double_t HighEE = 100.0)
{

  TString outputfilename = Form("Systematics-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root",Type.Data(),LowMult, HighMult, LowEE, HighEE);
  TFile* CreateFile = new TFile (outputfilename, "recreate");
  CreateFile->Close();

  TH1F* hV0Radius = makeSystPlotsV0s(outputfilename,Type,LowMult,HighMult,LowEE,HighEE,"V0Radius");
  TH1F* hNegToPV = makeSystPlotsV0s(outputfilename,Type,LowMult,HighMult,LowEE,HighEE,"DCANegToPV");
  TH1F* hPosToPV = makeSystPlotsV0s(outputfilename,Type,LowMult,HighMult,LowEE,HighEE,"DCAPosToPV");
  TH1F* hV0Daught = makeSystPlotsV0s(outputfilename,Type,LowMult,HighMult,LowEE,HighEE,"DCAV0Daughters");
  TH1F* hV0CosPA = makeSystPlotsV0s(outputfilename,Type,LowMult,HighMult,LowEE,HighEE,"V0CosPA");
  TH1F *hCompetingSpecies = 0x0;//makeSystPlotsV0s(outputfilename, Type, LowMult, HighMult, LowEE, HighEE, "CompetingRejection");
  TH1F* hPLT = makeSystPlotsV0s(outputfilename,Type,LowMult,HighMult,LowEE,HighEE,"ProperLifetime");
  TH1F* hTPCNSigmas = makeSystPlotsV0s(outputfilename,Type,LowMult,HighMult,LowEE,HighEE,"TPCPIDNSigmas");
  TH1F* hTPCNClusters = makeSystPlotsV0s(outputfilename,Type,LowMult,HighMult,LowEE,HighEE,"TPCNClusters");
  TH1F* hSigExtBinCount = makeSystPlotsV0s(outputfilename,Type,LowMult,HighMult,LowEE,HighEE,"SigExtBinCount");
  TH1F* hGeantFluka = 0x0;//makeSystPlotsV0s(outputfilename,Type,LowMult,HighMult,LowEE,HighEE,"GeantFlukaCorrection");
  TH1F* hNCr = makeSystPlotsV0s(outputfilename,Type,LowMult,HighMult,LowEE,HighEE,"NCrossedOverFindable");

  //Topological contribution systematics
  TH1F* hSystTopological = (TH1F*)hV0Radius->Clone("hSystTopological");
  hSystTopological->Reset();

  for (int i = 1; i <= hV0Radius->GetNbinsX(); i++){
    Double_t V0Radius = hV0Radius->GetBinContent(i);
    Double_t NegToPV = hNegToPV->GetBinContent(i);
    Double_t PosToPV= hPosToPV->GetBinContent(i);
    Double_t V0Daught = hV0Daught->GetBinContent(i);
    Double_t V0CosPA = hV0CosPA->GetBinContent(i);
    Double_t NCr = hNCr->GetBinContent(i);

    hSystTopological->SetBinContent(i,
      TMath::Sqrt( V0Radius*V0Radius + NegToPV*NegToPV + PosToPV*PosToPV +
                   V0Daught*V0Daught + V0CosPA*V0CosPA +
                   V0CosPA*V0CosPA + NCr*NCr)
    );
  }

  TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
  xlabel-> SetNDC();
  xlabel-> SetTextColor(1);
  xlabel-> SetTextSize(0.07);
  xlabel-> SetTextAlign(22);
  xlabel-> SetTextAngle(0);

  TLatex *xlabeltiny = new TLatex();
  xlabeltiny->SetTextFont(42);
  xlabeltiny->SetNDC();
  xlabeltiny->SetTextColor(1);
  xlabeltiny->SetTextSize(0.03);
  xlabeltiny->SetTextAlign(22);
  xlabeltiny->SetTextAngle(0);

  //Canvas
  //Topological Displayed
  TCanvas * Top = new TCanvas("Top","",1000,1000);
	Top->SetRightMargin(0.09);
	Top->SetLeftMargin(0.15);
	Top->SetBottomMargin(0.15);
  TLegend* l = new TLegend (0.5,0.55,0.89,0.89);
  l->SetBorderSize(0);
  l->AddEntry(hV0Radius,"V0 Radius","L");
  l->AddEntry(hPosToPV,"DCA Pos To PV","L");
  l->AddEntry(hV0Daught,"DCA V0 Daughters","L");
  l->AddEntry(hV0CosPA,"V0 Cosine of PA","L");
  l->AddEntry(hNCr, "N Crossed Rows/Findable");
  l->SetTextSize(0.03);
  TH1D* h = (TH1D*)hV0Radius->Clone("h");
  h->Reset();
  h->SetStats(0);
  h->GetYaxis()->SetRangeUser(0.,0.1);
  h->GetXaxis()->SetRangeUser(0.4,8.);
  if (Type.Contains("K0Short")){
    h->GetYaxis()->SetRangeUser(0., 0.05);
    h->GetXaxis()->SetRangeUser(0., 10.);
  }
  hV0Radius->SetLineColor(kRed+1);
  hPosToPV->SetLineColor(kOrange);
  hV0Daught->SetLineColor(kSpring);
  hV0CosPA->SetLineColor(kAzure+8);
  hV0Radius->SetMarkerStyle(1);
  hPosToPV->SetMarkerStyle(1);
  hV0Daught->SetMarkerStyle(1);
  hV0CosPA->SetMarkerStyle(1);
  hNCr->SetMarkerStyle(1);
  hV0Radius->SetLineWidth(2);
  hPosToPV->SetLineWidth(2);
  hV0Daught->SetLineWidth(2);
  hV0CosPA->SetLineWidth(2);
  hNCr->SetLineWidth(2);
  hNCr->SetLineColor(kBlue);
  hV0Radius->SetTitle(Form("Contribution of Topological Variables"));
  hV0Radius->GetYaxis()->SetRangeUser(0.,0.1);
  h->Draw();
  hV0Radius->Draw("HIST SAME");
  hPosToPV->Draw("HIST SAME");
  hV0Daught->Draw("HIST SAME");
  hV0CosPA->Draw("HIST SAME");
  hNCr->Draw("HIST SAME");
  l->Draw("SAME");
  if (Type.Contains("Lambda"))
    xlabel->DrawLatex(0.25, 0.8, "#Lambda + #bar{#Lambda}");
  else xlabel->DrawLatex(0.25, 0.8, "K^{0}_{S}");
  xlabeltiny->DrawLatex(0.35, 0.73, Form("V0M %.0f-%.0f SPDcl %.0f-%.0f", LowEE, HighEE, LowMult, HighMult));
  Top->SaveAs(Form("images/spectra/%s-Topological_SPDCl-%03.0f_%03.0f-V0M-%03.0f_%03.0f.png",Type.Data(),LowMult, HighMult, LowEE, HighEE));

  //Other selections contribution systematics
  TH1F* hSystOthers = (TH1F*)hV0Radius->Clone("hSystOthers");
  hSystOthers->Reset();

  for (int i = 1; i<= hV0Radius->GetNbinsX(); i++){
    Double_t PLT = hPLT->GetBinContent(i);
    Double_t CompetingSpecies = 0.;//hCompetingSpecies->GetBinContent(i);
    Double_t TPCNClusters = hTPCNClusters->GetBinContent(i);
    Double_t TPCNSigmas = hTPCNSigmas->GetBinContent(i);
    Double_t SigExtBinCount = hSigExtBinCount->GetBinContent(i);
    Double_t geant = 0; //hGeantFluka->GetBinContent(i);

    hSystOthers->SetBinContent(i,
      TMath::Sqrt(CompetingSpecies*CompetingSpecies + TPCNClusters*TPCNClusters +
                  TPCNSigmas*TPCNSigmas + SigExtBinCount*SigExtBinCount + PLT*PLT + geant*geant)
    );
  }

  //Selection Displayed
  TCanvas * Sel = new TCanvas("Sel","",1000,1000);
	Sel->SetRightMargin(0.09);
	Sel->SetLeftMargin(0.15);
	Sel->SetBottomMargin(0.15);
  TLegend* l1 = new TLegend (0.45,0.55,0.85,0.89);
  l1->SetBorderSize(0);
  l1->AddEntry(hPLT,"Proper Life Time","L");
  l1->AddEntry(hCompetingSpecies,"Competing Species","L");
  l1->AddEntry(hTPCNClusters, "TPC N of Clusters", "L");
  l1->AddEntry(hTPCNSigmas,"TPC N Sigmas","L");
  l1->AddEntry(hSigExtBinCount,"Sigmas for Sgn Extraction","L");
  l1->SetTextSize(0.03);
  hPLT->SetLineColor(kRed+1);
  hSigExtBinCount->SetLineColor(kSpring-1);
  hTPCNClusters->SetLineColor(kAzure-4);
  hTPCNSigmas->SetLineColor(kBlue+2);
  //hCompetingSpecies->SetLineColor(kBlack);
  hPLT->SetMarkerStyle(1);
  hPLT->SetLineWidth(2);
  //hCompetingSpecies->SetMarkerStyle(1);
  //hCompetingSpecies->SetLineWidth(2);
  hTPCNClusters->SetMarkerStyle(1);
  hTPCNClusters->SetLineWidth(2);
  hTPCNSigmas->SetMarkerStyle(1);
  hTPCNSigmas->SetLineWidth(2);
  hSigExtBinCount->SetMarkerStyle(1);
  hSigExtBinCount->SetLineWidth(2);
  hPLT->SetTitle("Contribution of Selection Variables");
  hPLT->GetYaxis()->SetTitleOffset(1.2);
  hPLT->GetYaxis()->SetRangeUser(0.,0.1);
  hPLT->GetXaxis()->SetRangeUser(0.4,8.);
  if (Type.Contains("K0Short")){
    hPLT->GetYaxis()->SetRangeUser(0., 0.05);
    hPLT->GetXaxis()->SetRangeUser(0., 10.);
  }
  hPLT->Draw("HIST ");
  //hCompetingSpecies->Draw("HIST SAME ");
  hTPCNClusters->Draw("HIST SAME ");
  hTPCNSigmas->Draw("HIST SAME ");
  hSigExtBinCount->Draw("HIST SAME ");
  l1->Draw("SAME");
  Sel->SaveAs(Form("images/spectra/%s-Selection_SPDCl-%03.0f_%03.0f-V0M-%03.0f_%03.0f.png",Type.Data(),LowMult, HighMult, LowEE, HighEE));
  if (Type.Contains("Lambda"))
    xlabel->DrawLatex(0.25, 0.8, "#Lambda + #bar{#Lambda}");
  else
    xlabel->DrawLatex(0.25, 0.8, "K^{0}_{S}");
  xlabeltiny->DrawLatex(0.35, 0.73, Form("V0M %.0f-%.0f SPDcl %.0f-%.0f", LowEE, HighEE, LowMult, HighMult));
  //Output File
  TFile* Write = new TFile (outputfilename, "UPDATE");
  hSystTopological->Write();
  hSystOthers->Write();
  hV0Radius->Write();
  hPosToPV->Write();
  hNegToPV->Write();
  hV0Daught->Write();
  hV0CosPA->Write();
  hPLT->Write();
  //hCompetingSpecies->Write();
  hTPCNClusters->Write();
  hTPCNSigmas->Write();
  hSigExtBinCount->Write();
  hNCr->Write();
}

//-------------------------------------------------------------------------------
//------------------------- DEFINE FUNCTIONS ------------------------------------
//-------------------------------------------------------------------------------

//---------------------------------------------------------------
TH1F* makeSystPlotsV0s(
          TString filename,
          TString lType = "Xi",
          Double_t lMultBoundLo = 0.0,
          Double_t lMultBoundHi = 100.0,
          Double_t lEEBoundLo = 0.0,
          Double_t lEEBoundHi = 100.0,
          TString lWhichSystVar = "V0Radius")
{
  cout<<"  "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"                  Make systematics for V0s                   "<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<endl;

  TString lCascType = "";
  TString lAntiCascType = "";
  if (lType.Contains("Xi")){
    lCascType = "XiMinus";
    lAntiCascType = "XiPlus";
  }
  if (lType.Contains("Omega")){
    lCascType = "OmegaMinus";
    lAntiCascType = "OmegaPlus";
  }
  if (lType.Contains("Lambda")){
    lCascType = "Lambda";
    lAntiCascType = "AntiLambda";
  }
  if (lType.Contains("K0Short")){
    lCascType = "K0Short";
    lAntiCascType = "";
  }

  TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
  xlabel-> SetNDC();
  xlabel-> SetTextColor(1);
  xlabel-> SetTextSize(0.07);
  xlabel-> SetTextAlign(22);
  xlabel-> SetTextAngle(0);

  TLatex *xlabeltiny = new TLatex();
  xlabeltiny->SetTextFont(42);
  xlabeltiny->SetNDC();
  xlabeltiny->SetTextColor(1);
  xlabeltiny->SetTextSize(0.03);
  xlabeltiny->SetTextAlign(22);
  xlabeltiny->SetTextAngle(0);

  //Definition of variables
  //Number of variation cuts
  const int nfiles = 5;
  TH1F* hMaxDev;
  TH1F* lHistPt[nfiles], *lAntiHistPt[nfiles], *lHistPtXi[nfiles];
  TH1F* lHistPtRaw[nfiles], *lAntiHistPtRaw[nfiles];
  TString lDataFilename[nfiles], lAntiDataFilename[nfiles];
  TFile* InputFile[nfiles], *AntiInputFile[nfiles];
  TH1F* hCut[nfiles-1];

  //Set data files
  TString lSystFile = Form("files/%s/Results-Systematics", lCascType.Data());
  TString lAntiSystFile = Form("files/%s/Results-Systematics", lAntiCascType.Data());
  lSystFile.Append( Form( "-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f-", lCascType.Data(), lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi ) );
  if (!lType.Contains("K0Short")) lAntiSystFile.Append(Form("-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f-", lAntiCascType.Data(), lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi));

  //Particle
  lDataFilename[0] = Form("files/%s/Results-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_FDNoFD.root", lCascType.Data(), lCascType.Data(), lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi);
  InputFile[0] = new TFile(lDataFilename[0].Data(),"READ");
  for (int i = 1; i <= nfiles-1; i++){
    lDataFilename[i] = lSystFile + lWhichSystVar + Form("-%i.root",i);
    InputFile[i] = new TFile(lDataFilename[i].Data(),"READ");
  }
  //Anti Particle
  if (!lType.Contains("K0Short")) {
    lAntiDataFilename[0] = Form("files/%s/Results-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_FDNoFD.root", lAntiCascType.Data(), lAntiCascType.Data(), lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi);
    AntiInputFile[0] = new TFile(lAntiDataFilename[0].Data(),"READ");
    for (int i = 1; i <= nfiles-1; i++){
      lAntiDataFilename[i] = lAntiSystFile + lWhichSystVar + Form("-%i.root",i);
      AntiInputFile[i] = new TFile(lAntiDataFilename[i].Data(),"READ");
    }
  }

  //Inizialize Histos
  for (int i = 0; i < nfiles; i++){
    lHistPt[i] = (TH1F *) InputFile[i]->Get(Form("fHistPt%s", lCascType.Data()));
    lHistPtRaw[i] = (TH1F *) InputFile[i]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
    if (!lType.Contains("K0Short")) {
      lAntiHistPt[i] = (TH1F *)AntiInputFile[i]->Get(Form("fHistPt%s", lAntiCascType.Data()));
      lAntiHistPtRaw[i] = (TH1F *) AntiInputFile[i]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
    }

    //
    lHistPtXi[i] = (TH1F *)lHistPt[i]->Clone(Form("HistPtXi%i",i));

    if (!lType.Contains("K0Short")) {
      lHistPtXi[i]->Reset();
      for (int bin = 1; bin <= lHistPt[i]->GetNbinsX(); bin++)
      {
        lHistPtXi[i]->SetBinContent(bin, lHistPt[i]->GetBinContent(bin) + lAntiHistPt[i]->GetBinContent(bin));
        lHistPtXi[i]->SetBinError(bin, TMath::Sqrt(lHistPt[i]->GetBinError(bin) * lHistPt[i]->GetBinError(bin) +
                                                   lAntiHistPt[i]->GetBinError(bin) * lAntiHistPt[i]->GetBinError(bin)));
      }
    }

    //
    if (i > 0) hCut[i - 1] = (TH1F *)lHistPtXi[i]->Clone(Form("hVarCut%s-%i", lWhichSystVar.Data(), i)); // will contain var cuts
  }

  hMaxDev = (TH1F*)lHistPtXi[0]->Clone("hMaxDev");
  hMaxDev->Reset();

  //Divide following Roger-Barlow prescription
  for (int k = 1; k < nfiles; k++){
    DivideAndComputeRogerBarlow(hCut[k-1],lHistPtXi[0]);
  }

  //Compute Max Deviation
  double binvalue[nfiles-1];
  double maxvalue = 0.;
  for (int i = 1; i<= lHistPtXi[0]->GetNbinsX(); i++){
    for (int k = 1; k < nfiles; k++){
      double dev = hCut[k-1]->GetBinContent(i);
      binvalue[k-1] = TMath::Abs(dev-1); //deviation from 1 x each file
    }
    maxvalue = binvalue[0];
    int counter = 0;
    for (int k = 1; k < (nfiles-1); k++){
      maxvalue = max(maxvalue,binvalue[k]);
      if (maxvalue == binvalue[k]) counter = k;
    }
    hMaxDev->SetBinError(i,hCut[counter]->GetBinError(i));
    hMaxDev->SetBinContent(i,TMath::Abs(maxvalue));
    //if (TMath::Abs(maxvalue) == 0) hMaxDev->SetBinContent(i,0.0);
  }

  for (int i = 1; i<= lHistPt[0]->GetNbinsX(); i++){
    hMaxDev->SetBinContent(i,PassRogerBarlowCriterion(1,hMaxDev,i));
    if (!PassRogerBarlowCriterion(1, hMaxDev, i)) hMaxDev->SetBinError(i, 0.0);
  }
  hMaxDev->SetLineWidth(2);
  hMaxDev->SetTitle(Form("Syst %s",lWhichSystVar.Data()));
  hMaxDev->SetName(Form("hMaxDev%s",lWhichSystVar.Data()));

  //Prepare Canvas
  //Max Deviation
  hMaxDev->GetXaxis()->SetRangeUser(.4,8.);
  if (lCascType.Contains("Lambda"))
    hMaxDev->GetXaxis()->SetRangeUser(0.4, 8.);
  else
    hMaxDev->GetXaxis()->SetRangeUser(0., 10.);
  hMaxDev->GetYaxis()->SetRangeUser(0., hMaxDev->GetMaximum()*2);
  hMaxDev->SetYTitle("max rel. dev.");
  hMaxDev->SetTitle(Form("%s",lWhichSystVar.Data()));
  hMaxDev->GetYaxis()->SetTitleSize(0.05);
  hMaxDev->GetYaxis()->SetTitleOffset(1.4);
  hMaxDev->GetXaxis()->SetTitleSize(0.05);
  hMaxDev->GetXaxis()->SetTitleOffset(1.);
  hMaxDev->SetMarkerStyle(20);
  hMaxDev->SetMarkerSize(1.1);
  hMaxDev->SetMarkerColor(kBlack);
  hMaxDev->SetLineColor(kBlack);
  TCanvas* maxdev = new TCanvas("maxdev"," ",1000,1000);
  maxdev->SetRightMargin(0.15);
  maxdev->SetLeftMargin(0.15);
  maxdev->SetBottomMargin(0.15);
  maxdev->SetGridy();
  maxdev->SetGridx();
  hMaxDev->SetStats(kFALSE);
  hMaxDev->SetName(lWhichSystVar);
  hMaxDev->Draw();
  if (lCascType.Contains("Lambda")) xlabel-> DrawLatex(0.25, 0.8,Form("#%s + #bar{#%s}",lCascType.Data(),lAntiCascType.Data()));
  else xlabel-> DrawLatex(0.25, 0.8,"K^{0}_{S}");
  xlabeltiny->DrawLatex(0.35, 0.73, Form("V0M %.0f-%.0f SPDcl %.0f-%.0f", lEEBoundLo, lEEBoundHi, lMultBoundLo, lMultBoundHi));

  maxdev->SaveAs(Form("images/spectra/%s-MaxRelDev%s-V0-%03.0f_%03.0f-V0M-%03.0f_%03.0f.png",lType.Data(),lWhichSystVar.Data(),lMultBoundLo,lMultBoundHi,lEEBoundLo,lEEBoundHi));

  //Cut Variation
  TCanvas* cutvar = new TCanvas("cutvar"," ",1000,1000);
  cutvar->SetRightMargin(0.08);
  cutvar->SetLeftMargin(0.15);
  cutvar->SetBottomMargin(0.15);
  cutvar->SetGridy();
  hCut[0]->Draw();
  TLegend* legend = new TLegend (0.7,.7,0.91,0.89);
  legend->SetBorderSize(0);
  if (lCascType.Contains("Lambda")) xlabel-> DrawLatex(0.25, 0.8,Form("#%s + #bar{#%s}",lCascType.Data(),lCascType.Data()));
  else xlabel->DrawLatex(0.25, 0.8, "K^{0}_{S}");
  xlabeltiny->DrawLatex(0.35, 0.73, Form("V0M %.0f-%.0f SPDcl %.0f-%.0f", lEEBoundLo, lEEBoundHi, lMultBoundLo, lMultBoundHi));

  hCut[0]->SetMarkerColor(kRed);
  hCut[1]->SetMarkerColor(kBlue+1);
  hCut[2]->SetMarkerColor(kGreen+3);
  hCut[3]->SetMarkerColor(kAzure+9);
  hCut[0]->SetLineColor(kRed);
  hCut[1]->SetLineColor(kBlue+1);
  hCut[2]->SetLineColor(kGreen+3);
  hCut[3]->SetLineColor(kAzure+9);
  legend->AddEntry(hCut[0],"very loose","LEP");
  legend->AddEntry(hCut[1],"loose","LEP");
  legend->AddEntry(hCut[2],"tight","LEP");
  legend->AddEntry(hCut[3],"very tight","LEP");

  for (int k = 0; k < nfiles-1; k++){
    hCut[k]->SetLineWidth(2);
    if (lCascType.Contains("Lambda")) hCut[k]->GetXaxis()->SetRangeUser(0.4,8.);
    else hCut[k]->GetXaxis()->SetRangeUser(0.,10.);
    hCut[k]->GetYaxis()->SetRangeUser(1 - hMaxDev->GetMaximum() * 2., 1 + hMaxDev->GetMaximum() * 2.0);
    hCut[k]->SetYTitle("Yield^{syst-cut} / Yield^{def-cut}");
    hCut[k]->SetTitle(Form("%s",lWhichSystVar.Data()));
    hCut[k]->GetYaxis()->SetTitleSize(0.05);
    hCut[k]->GetYaxis()->SetTitleOffset(1.1);
    hCut[k]->GetXaxis()->SetTitleSize(0.05);
    hCut[k]->GetXaxis()->SetTitleOffset(0.8);
    hCut[k]->SetMarkerStyle(20);
    hCut[k]->SetMarkerSize(1.1);
  }
  hCut[0]->SetStats(kFALSE);
  for (int k = 1; k < nfiles-1; k++){
    hCut[k]->Draw("SAME");
  }
  legend->SetTextSize(0.035);
  legend->Draw("SAME");
  cutvar->SaveAs(Form("images/spectra/%s-CutVar-%s-V0-%03.0f_%03.0f-V0M-%03.0f_%03.0f.png",lType.Data(),lWhichSystVar.Data(),lMultBoundLo,lMultBoundHi,lEEBoundLo,lEEBoundHi));

  TFile* lOutputFile = TFile::Open(filename, "UPDATE");
  TDirectoryFile *lDirCutVar = new TDirectoryFile(Form("Var%s",lWhichSystVar.Data()),Form("Variation Cuts %s",lWhichSystVar.Data()));
  lDirCutVar->cd();
  for (int k = 0; k < nfiles-1; k++){
    hCut[k]->Write();
  }
  //lOutputFile->cd();
  lOutputFile->Close();

  //Return Max Dev Histo
  return hMaxDev;

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

//----------------------------------------------------------------------------------------------------
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 ){
  //Use Roger Barlow "sigma_{delta}" as errors for ratios
  Double_t lh1NBins = h1->GetNbinsX();
  Double_t lh2NBins = h2->GetNbinsX();

  if( lh1NBins != lh2NBins ){
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return;
  }

  Double_t lSigmaDelta[100];
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){
    //Computation of roger barlow sigma_{delta}
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(h1->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
    //Computation of relationship to h2 for plotting in ratio plot
    if ( h2->GetBinContent(i) > 1e-12 ){
      lSigmaDelta[i] /= h2->GetBinContent(i);
    }else{
      lSigmaDelta[i] = 0;
    }
  }
  //Regular Division
  h1 -> Divide( h2 );
  //Replace Erorrs
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){
    h1->SetBinError(i, lSigmaDelta[i]);
  }
}

//----------------------------------------------------------------------------------------------------
double PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin){
  double dev = h->GetBinContent(bin);
  double RBsigma = h->GetBinError(bin);

  if (dev>(nsigmas*RBsigma)) {return dev;}
  else {return 0.;}
}