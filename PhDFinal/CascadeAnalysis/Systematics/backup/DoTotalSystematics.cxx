#include <iostream>
#include <fstream>
//#define max(a,b) (a>b ? a : b)

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 );
double PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin);

TH1F* makeSystPlotsV0s(
  TString filename,
  TString lType = "Xi",
	TString lWhichEstimator = "SPDClusters", 
	Double_t lMultBoundLo = 0.0, 
	Double_t lMultBoundHi = 100.0, 
	Double_t lEEBoundLo = 0.0, 
	Double_t lEEBoundHi = 100.0, 
	TString lWhichSystVar = "V0Radius" );


//--------------------------------------------------------------
//------------------- MAIN FUNCTION ----------------------------
//--------------------------------------------------------------
void DoTotalSystematics(
	Double_t LowMult = 10.0, 
	Double_t HighMult = 20, 
	Double_t LowEE = 0.0, 
	Double_t HighEE = 100.0,
	Bool_t DoMult = kTRUE){

  Bool_t DoEE = !DoMult; 
  TString Type = "Omega";
  
  TString outputfilename = Form("SystematicsFinalResults-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root",Type.Data(),LowMult, HighMult, LowEE, HighEE);
  TFile* CreateFile = new TFile (outputfilename, "recreate");
  CreateFile->Close();

  TString fWhichEstimator = "SPDClusters";

  TH1F* hV0Radius = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"V0Radius");
  TH1F* hCascRadius = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"CascRadius");
  TH1F* hNegToPV = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"DCANegToPV");
  TH1F* hPosToPV = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"DCAPosToPV");
  TH1F* hV0ToPV = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"DCAV0ToPV");
  TH1F* hBachToPV = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"DCABachToPV");
  TH1F* hV0Daught = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"DCAV0Daughters");
  TH1F* hCascDaught = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"DCACascDaughters");
  TH1F* hV0CosPA = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"V0CosPA");
  TH1F* hCascCosPA = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"CascCosPA");
  TH1F* hV0Mass = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"V0Mass");
  TH1F* hCompetingSpecies = 0x0;//makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"CompetingSpecies");
  TH1F* hPLT = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"ProperLifetime");
  TH1F* hTPCNSigmas = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"TPCPIDNSigmas");
  TH1F* hTPCNClusters = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"TPCNClusters");
  TH1F* hSigExtBinCount = 0x0;//makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"SigExtBinCount");
  //TH1F* hGeantFluka = makeSystPlotsV0s(outputfilename,Type,fWhichEstimator.Data(),LowMult,HighMult,LowEE,HighEE,"GeantFlukaCorrection");

  //=========================================================
  //========= Topological contribution systematics ==========
  //=========================================================

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

  //=========================================================
  //======== Other select. contribution systematics =========
  //=========================================================

  TH1F* hSystOthers = (TH1F*)hV0Radius->Clone("hSystOthers");
  hSystOthers->Reset();

  for (int i = 1; i<= hV0Radius->GetNbinsX(); i++){
    Double_t V0Mass = hV0Mass->GetBinContent(i);
    Double_t PLT = hPLT->GetBinContent(i);
    Double_t CompetingSpecies = 0;
    Double_t TPCNClusters = hTPCNClusters->GetBinContent(i);
    Double_t TPCNSigmas = hTPCNSigmas->GetBinContent(i);
    Double_t SigExtBinCount = 0;//hSigExtBinCount->GetBinContent(i);
    Double_t geant = 0;//hGeantFluka->GetBinContent(i);

    hSystOthers->SetBinContent(i, 
      TMath::Sqrt(V0Mass*V0Mass  + CompetingSpecies*CompetingSpecies + TPCNClusters*TPCNClusters +
                  TPCNSigmas*TPCNSigmas + SigExtBinCount*SigExtBinCount + PLT*PLT + geant*geant) 
    );
  }

  //=========================================================
  //============ Other contribution systematics =============
  //=========================================================
  
  //G4
  TFile* fileG4 = TFile::Open(Form("G4Syst.root"));
  TH1D* hG4 = (TH1D*)fileG4->Get(Form("h%sTot",Type.Data())); 

  //OOB
  TFile* fileOOBSyst = TFile::Open(Form("OOBSyst%s.root",Type.Data()));
  TH1D* hoob = (TH1D*)fileOOBSyst->Get("hMaxDev");    
  
  //INBUNCH
  TFile* inbunch = TFile::Open(Form("SystematicsPileUp-%s.root",Type.Data()));
  TH1D* inbunch040 = (TH1D*)inbunch->Get(Form("fHistPileUp%s_0_40",Type.Data()));
  TH1D* inbunch40100 = (TH1D*)inbunch->Get(Form("fHistPileUp%s_40_100",Type.Data()));

  //MATBUDG
  TFile* matbudg  = TFile::Open(Form("systMatBudget-%s.root",Type.Data()));
  TF1* fmatbudg = (TF1*)matbudg->Get("funcFit");
  TH1F* hmatbudg = (TH1F*)hV0Radius->Clone("hmatbudg");
  TH1F* hibp = (TH1F*)hV0Radius->Clone("hibp");
  TH1F* hconst = (TH1F*)hV0Radius->Clone("hconst");
  
  TH1F* hSystExtra = (TH1F*)hV0Radius->Clone("hSystExtra");
  hSystExtra->Reset();

  for (int i = 1; i<= hV0Radius->GetNbinsX(); i++){
    
    Double_t MultIndipEffic = 0.02;
    hconst->SetBinContent(i,MultIndipEffic);
    Double_t MaterialBudget = TMath::Abs((1-fmatbudg->Eval(hV0Radius->GetBinCenter(i))))/2;
    hmatbudg->SetBinContent(i,MaterialBudget);
    Double_t ibp = inbunch040->GetBinError(i)*40./100 + inbunch40100->GetBinError(i)*60./100; 
    hibp->SetBinContent(i,ibp);
    Double_t oob = hoob->GetBinContent(i)/2; 
    hoob->SetBinContent(i,hoob->GetBinContent(i)/2); 
    Double_t g4 = hG4->GetBinContent(i);

    hSystExtra->SetBinContent(i, 
      TMath::Sqrt( MultIndipEffic*MultIndipEffic + MaterialBudget*MaterialBudget + ibp*ibp + oob*oob + g4*g4) 
    );
  }

  //Total systematics
  TH1F* hSystTot = (TH1F*)hV0Radius->Clone("hSystTot");  
  hSystTot->Reset();
  for (int i = 1; i<= hV0Radius->GetNbinsX(); i++){
      
      hSystTot->SetBinContent(i, 
        TMath::Sqrt(
          hSystOthers->GetBinContent(i)*hSystOthers->GetBinContent(i) +
          hSystTopological->GetBinContent(i)*hSystTopological->GetBinContent(i) + 
          hSystExtra->GetBinContent(i)*hSystExtra->GetBinContent(i)
          )
        );        
  }

  //Canvases---------------------------------------------
  TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
  xlabel-> SetNDC();
  xlabel-> SetTextColor(1);
  xlabel-> SetTextSize(0.07);
  xlabel-> SetTextAlign(22);
  xlabel-> SetTextAngle(0);

  //Topological Displayed
  TCanvas * Top = new TCanvas("Top","",1000,1000);
	Top->SetRightMargin(0.09);
	Top->SetLeftMargin(0.15);
	Top->SetBottomMargin(0.15);
  TLegend* l = new TLegend (0.5,0.65,0.89,0.89);
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
  hV0Radius->SetTitle(Form("Contribution of Topological Variables %s",Type.Data()));
  hV0Radius->GetYaxis()->SetTitleOffset(1.2);
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
  xlabel-> DrawLatex(0.25, 0.8,Form("#%s",Type.Data()));
  Top->SaveAs(Form("images/TopologicalSystematics%s_V0-%03.0f_%03.0f-V0M-%03.0f_%03.0f.png",Type.Data(),LowMult, HighMult, LowEE, HighEE));

  //Selection Displayed
  TCanvas * Sel = new TCanvas("Sel","",1000,1000);
	Sel->SetRightMargin(0.09);
	Sel->SetLeftMargin(0.15);
	Sel->SetBottomMargin(0.15);
  TLegend* l1 = new TLegend (0.45,0.65,0.85,0.89);
  l1->SetBorderSize(0);
  l1->AddEntry(hPLT,"Proper Life Time","L");
  l1->AddEntry(hV0Mass,"V0 Mass","L");
  if (Type == "OmegaMinus" || Type == "OmegaPlus") l1->AddEntry(hCompetingSpecies,"Competing Species","L");
  l1->AddEntry(hTPCNClusters,"TPC N of Clusters","L");
  l1->AddEntry(hTPCNSigmas,"TPC N Sigmas","L");
 // l1->AddEntry(hSigExtBinCount,"Sigmas for Sgn Extraction","L");
  l1->SetTextSize(0.03);
  hPLT->SetLineColor(kRed+1);
  hV0Mass->SetLineColor(kOrange+1);
 // hSigExtBinCount->SetLineColor(kSpring-1);
  hTPCNClusters->SetLineColor(kAzure-4);
  hTPCNSigmas->SetLineColor(kBlue+2);
  if (Type == "OmegaMinus" || Type == "OmegaPlus") hCompetingSpecies->SetLineColor(kBlack);
  hPLT->SetMarkerStyle(1);
  hPLT->SetLineWidth(2);
  hV0Mass->SetMarkerStyle(1);
  hV0Mass->SetLineWidth(2);
  if (Type == "OmegaMinus" || Type == "OmegaPlus") {
    hCompetingSpecies->SetMarkerStyle(1);
    hCompetingSpecies->SetLineWidth(2);
  }
  hTPCNClusters->SetMarkerStyle(1);
  hTPCNClusters->SetLineWidth(2);
  hTPCNSigmas->SetMarkerStyle(1);
  hTPCNSigmas->SetLineWidth(2);
//  hSigExtBinCount->SetMarkerStyle(1);
//  hSigExtBinCount->SetLineWidth(2);
  hPLT->SetTitle("Contribution of Selection Variables");
  hPLT->GetYaxis()->SetTitleOffset(1.2);
  hPLT->Draw("HIST ");
  hV0Mass->Draw("HIST SAME ");
  if (Type == "OmegaMinus" || Type == "OmegaPlus") hCompetingSpecies->Draw("HIST SAME ");
  hTPCNClusters->Draw("HIST SAME ");
  hTPCNSigmas->Draw("HIST SAME ");
//  hSigExtBinCount->Draw("HIST SAME ");
  l1->Draw("SAME");
  xlabel-> DrawLatex(0.25, 0.8,Form("#%s",Type.Data()));
  Sel->SaveAs(Form("images/%s-SelectionSystematics_V0-%03.0f_%03.0f-V0M-%03.0f_%03.0f.png",Type.Data(),LowMult, HighMult, LowEE, HighEE));

  //Extra syst
  TCanvas * Extra = new TCanvas("Extra","",1000,1000);
	Extra->SetRightMargin(0.09);
	Extra->SetLeftMargin(0.15);
	Extra->SetBottomMargin(0.15);
  TLegend* lextra = new TLegend (0.3,0.65,0.89,0.89);
  lextra->SetBorderSize(0);
  lextra->AddEntry(hmatbudg,"Material Budget (from pp@13TeV)","L");
  lextra->AddEntry(hibp,"In-bunch pileup (from pp@13TeV)","L");
  lextra->AddEntry(hoob,"Out-of-bunch pileup (OOB)","L");
  lextra->AddEntry(hconst,"p_{T} indipendent (mult integr. eff.)","L");
  lextra->AddEntry(hG4,"G4 anti-proton correction","L");
  lextra->SetTextSize(0.03);
  hmatbudg->SetLineColor(kRed+1);
  hibp->SetLineColor(kSpring);
  hoob->SetLineColor(kBlue);
  hconst->SetLineColor(kBlack);
  hG4->SetLineColor(kMagenta);
  hmatbudg->SetLineWidth(2);
  hibp->SetLineWidth(2);
  hoob->SetLineWidth(2);
  hconst->SetLineWidth(2);
  hconst->SetLineStyle(7);
  hG4->SetLineWidth(2);
  hmatbudg->SetTitle(Form("Other Contributions %s",Type.Data()));
  hmatbudg->GetYaxis()->SetTitleOffset(1.2);
  hmatbudg->Draw("HIST ");
  hibp->Draw("HIST SAME ");
  hconst->Draw("HIST SAME ");
  hoob->Draw("HIST SAME");
  hG4->Draw("HIST SAME");
  lextra->Draw("SAME");
  xlabel-> DrawLatex(0.25, 0.8,Form("#%s",Type.Data()));
  Extra->SaveAs(Form("images/ExtraSystematics%s_V0-%03.0f_%03.0f-V0M-%03.0f_%03.0f.png",Type.Data(),LowMult, HighMult, LowEE, HighEE));
  

  //Contributions Displayed
  TCanvas* cn = new TCanvas("cn","",1000,1000);
  cn->SetRightMargin(0.09);
	cn->SetLeftMargin(0.15);
	cn->SetBottomMargin(0.15);
  
  TLegend* legend = new TLegend(0.45,0.65,0.85,0.89);
  legend->SetBorderSize(0);
  legend->AddEntry(hSystTopological,"Topological systematics","L");
  legend->AddEntry(hSystOthers,"Selection cuts systematics","L");
  legend->AddEntry(hSystExtra,"Other contributions systematics","L");
  legend->AddEntry(hSystTot,"Total systematics","L");
  legend->SetTextSize(0.027);

  hSystTot->SetLineWidth(2);
  hSystTot->SetLineColor(kBlack);
  hSystOthers->SetLineWidth(2);
  hSystOthers->SetLineColor(kRed);
  hSystTopological->SetLineWidth(2);
  hSystTopological->SetLineColor(kBlue);
  hSystExtra->SetLineColor(kGreen);
  hSystExtra->SetLineWidth(2);
  TLatex *xlabel2 = new TLatex();
	xlabel2->SetTextFont(42);
  xlabel2-> SetNDC();
  xlabel2-> SetTextColor(1);
  xlabel2-> SetTextSize(0.04);
  xlabel2-> SetTextAlign(22);
  xlabel2-> SetTextAngle(0);
  xlabel2-> DrawLatex(0.32, 0.81, "V0M [0,100]% V0M [0,100]%");

  hSystTot->SetYTitle("max. rel. deviation (%)");
  hSystTot->GetYaxis()->SetRangeUser(-0.005,.15);
  hSystTot->GetYaxis()->SetTitleOffset(1.);
  hSystTot->SetTitle(Form("Systematics contributions for %s",Type.Data()));
  hSystTot->GetYaxis()->SetTitleOffset(1.2);
  hSystTot->Draw();
  hSystOthers->Draw("SAME");
  hSystTopological->Draw("SAME");
  hSystExtra->Draw("SAME");
  legend->Draw("SAME");  
  xlabel-> DrawLatex(0.25, 0.8,Form("#%s",Type.Data()));
  cn->SaveAs(Form("images/%s-Systematics-Displayed_V0-%03.0f_%03.0f-V0M-%03.0f_%03.0f.png",Type.Data(),LowMult, HighMult, LowEE, HighEE));

  //Total syst canvas shaded 
  TCanvas *cfinal = new TCanvas("cfinal","Total Systematics",1000,1000);
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
  hSystTotNeg->GetYaxis()->SetTitle("max. rel. deviation (%)");
  hSystTotNeg->GetYaxis()->SetRangeUser(0.-0.2,0.+0.2);
  hSystTotNeg->Draw("E2");
  xlabel-> DrawLatex(0.25, 0.8,Form("#%s",Type.Data()));
  cfinal->SaveAs(Form("images/FinalSystContribPlot%s-%03.0f_%03.0f-V0M-%03.0f_%03.0f.png",Type.Data(),LowMult, HighMult, LowEE, HighEE));
 
  //Output File
  TFile* Write = new TFile (outputfilename, "UPDATE");
  hSystTopological->Write();
  hSystOthers->Write();
  hSystExtra->Write();
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
  if (Type == "OmegaMinus" || Type == "OmegaPlus") hCompetingSpecies->Write();
  hTPCNClusters->Write();
  hTPCNSigmas->Write();
 // hSigExtBinCount->Write();

}

//-------------------------------------------------------------------------------
//------------------------- DEFINE FUNCTIONS ------------------------------------
//-------------------------------------------------------------------------------

//---------------------------------------------------------------
TH1F* makeSystPlotsV0s(
          TString filename,
          TString lType = "Xi",
          TString lWhichEstimator = "V0M",
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

  //Definition of variables
  //Number of variation cuts
  const int nfiles = 5;  
  TH1F* hMaxDev;
  TH1F* lHistPt[nfiles], *lAntiHistPt[nfiles], *lHistPtXi[nfiles];
  TH1F* lHistPtRaw[nfiles], *lAntiHistPtRaw[nfiles];
  TString lDataFilename[nfiles], lAntiDataFilename[nfiles];
  TFile* InputFile[nfiles], *AntiInputFile[nfiles];   
  TH1F* hCut[nfiles-1];
 
  TString folder = ".";
  if (lCascType.Contains("Xi")) folder = "XiToBeUpdated";

  //Set data files
  TString lSystFile = Form("sistematiche/%s/Results-Systematics",folder.Data());
  TString lAntiSystFile = Form("sistematiche/%s/Results-Systematics",folder.Data());
  lSystFile.Append( Form( "-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f-", lCascType.Data(), lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi ) );
  lAntiSystFile.Append( Form( "-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f-", lAntiCascType.Data(), lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi ) );
  
  //Particle
  lDataFilename[0] = Form( "sistematiche/%s/Results-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root", folder.Data(), lCascType.Data(), lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi );
  InputFile[0] = new TFile(lDataFilename[0].Data(),"READ");
  for (int i = 1; i <= nfiles-1; i++){
    lDataFilename[i] = lSystFile + lWhichSystVar + Form("-%i.root",i);
    InputFile[i] = new TFile(lDataFilename[i].Data(),"READ");
  }  
  //Anti Particle
  lAntiDataFilename[0] = Form( "sistematiche/%s/Results-%s-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f.root",  folder.Data(), lAntiCascType.Data(), lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi );
  AntiInputFile[0] = new TFile(lAntiDataFilename[0].Data(),"READ");
  for (int i = 1; i <= nfiles-1; i++){
    lAntiDataFilename[i] = lAntiSystFile + lWhichSystVar + Form("-%i.root",i);
    AntiInputFile[i] = new TFile(lAntiDataFilename[i].Data(),"READ");
  }

  //Inizialize Histos
  for (int i = 0; i < nfiles; i++){
    lHistPt[i] = (TH1F *) InputFile[i]->Get(Form("fHistPt%s", lCascType.Data()));
    lHistPtRaw[i] = (TH1F *) InputFile[i]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");      
    lAntiHistPt[i] = (TH1F *) AntiInputFile[i]->Get(Form("fHistPt%s", lAntiCascType.Data()));
    lAntiHistPtRaw[i] = (TH1F *) AntiInputFile[i]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");      
    //
    lHistPtXi[i] = (TH1F *)lHistPt[i]->Clone(Form("HistPtXi%i",i));
    lHistPtXi[i]->Reset();
    for (int bin = 1; bin <=  lHistPt[i]->GetNbinsX(); bin ++){
      lHistPtXi[i]->SetBinContent(bin, lHistPt[i]->GetBinContent(bin)+lAntiHistPt[i]->GetBinContent(bin));
      lHistPtXi[i]->SetBinError(bin, TMath::Sqrt(lHistPt[i]->GetBinError(bin)*lHistPt[i]->GetBinError(bin) +
                                               lAntiHistPt[i]->GetBinError(bin)*lAntiHistPt[i]->GetBinError(bin)));       
    }
    //
    if(i>0) hCut[i-1] = (TH1F*)lHistPtXi[i]->Clone(Form("hVarCut%s-%i",lWhichSystVar.Data(),i));//will contain var cuts
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
      binvalue[k-1] = TMath::Abs(dev-1);
    }
    maxvalue = binvalue[0];
    int counter = 0;
    for (int k = 1; k < nfiles; k++){ 
      maxvalue = max(maxvalue,binvalue[k]);
      if (maxvalue == binvalue[k]) counter = k;
    }   
    hMaxDev->SetBinError(i,hCut[counter]->GetBinError(i));
    hMaxDev->SetBinContent(i,TMath::Abs(maxvalue));
    //if (TMath::Abs(maxvalue) == 0) hMaxDev->SetBinContent(i,0.0);
  }
  for (int i = 1; i<= lHistPt[0]->GetNbinsX(); i++){
    hMaxDev->SetBinContent(i,PassRogerBarlowCriterion(1,hMaxDev,i));
  }
  hMaxDev->SetTitle(Form("Syst %s",lWhichSystVar.Data()));
  hMaxDev->SetName(Form("hMaxDev%s",lWhichSystVar.Data()));

  TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
  xlabel-> SetNDC();
  xlabel-> SetTextColor(1);
  xlabel-> SetTextSize(0.07);
  xlabel-> SetTextAlign(22);
  xlabel-> SetTextAngle(0);

  //Prepare Canvas
  //Max Deviation
  if(lType.Contains("Xi")) hMaxDev->GetXaxis()->SetRangeUser(0.6,6.5);
  else hMaxDev->GetXaxis()->SetRangeUser(0.9,5.5);
  hMaxDev->GetYaxis()->SetRangeUser(-0.0005,.15);
  hMaxDev->SetYTitle("max rel. dev.");
  hMaxDev->SetTitle(Form("%s",lWhichSystVar.Data()));
  hMaxDev->GetYaxis()->SetTitleSize(0.05);
  hMaxDev->GetYaxis()->SetTitleOffset(1.0);
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
  hMaxDev->SetLineWidth(2);  
  hMaxDev->SetName(lWhichSystVar);
  hMaxDev->Draw();
  xlabel-> DrawLatex(0.4, 0.8,Form("#%s",lType.Data()));
  maxdev->SaveAs(Form("images/%s-MaxRelDev%s-V0-%03.0f_%03.0f-V0M-%03.0f_%03.0f.png",lType.Data(),lWhichSystVar.Data(),lMultBoundLo,lMultBoundHi,lEEBoundLo,lEEBoundHi));

  //Cut Variation
  TCanvas* cutvar = new TCanvas("cutvar"," ",1000,1000);
  cutvar->SetRightMargin(0.08);
  cutvar->SetLeftMargin(0.15);
  cutvar->SetBottomMargin(0.15);
  cutvar->SetGridy();
  hCut[0]->Draw();
  TLegend* legend = new TLegend (0.7,.75,0.91,0.89);  
  legend->SetBorderSize(0);

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
    if(lType.Contains("Xi")) hCut[k]->GetXaxis()->SetRangeUser(0.6,6.5);
    else hCut[k]->GetXaxis()->SetRangeUser(0.9,5.5);
    hCut[k]->GetYaxis()->SetRangeUser(0.9,1.1);
    hCut[k]->SetYTitle("Yield^{syst-cut} / Yield^{def-cut}");
    hCut[k]->SetTitle(Form("%s",lWhichSystVar.Data()));
    hCut[k]->GetYaxis()->SetTitleSize(0.05);
    hCut[k]->GetYaxis()->SetTitleOffset(1.1);
    hCut[k]->GetXaxis()->SetTitleSize(0.05);
    hCut[k]->GetXaxis()->SetTitleOffset(0.8);
    hCut[k]->SetMarkerStyle(20);
    hCut[k]->SetMarkerSize(1.1);
    hCut[k]->SetLineWidth(2);
  }
  hCut[0]->SetStats(kFALSE);
  for (int k = 1; k < nfiles-1; k++){
    hCut[k]->Draw("SAME");
  }
  legend->SetTextSize(0.035);
  xlabel-> DrawLatex(0.4, 0.8,Form("#%s",lType.Data()));
  legend->Draw("SAME");    
  cutvar->SaveAs(Form("images/%s-CutVar-%s-V0-%03.0f_%03.0f-V0M-%03.0f_%03.0f.png",lType.Data(),lWhichSystVar.Data(),lMultBoundLo,lMultBoundHi,lEEBoundLo,lEEBoundHi));

  TFile* lOutputFile = TFile::Open(filename, "UPDATE");
  TDirectoryFile *lDirCutVar = new TDirectoryFile(Form("VariationCuts%s",lWhichSystVar.Data()),Form("Variation Cuts %s",lWhichSystVar.Data()));
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
