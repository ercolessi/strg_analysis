#include <iostream>
#include <fstream>
//#define max(a,b) (a>b ? a : b)

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );

TH1F* makeSystPlotsV0s(TString lV0Type = "XiMinus", 
	TString lWhichEstimator = "V0M", 
	Double_t lMultBoundLo = 0.0, 
	Double_t lMultBoundHi = 100.0, 
	Double_t lEEBoundLo = 0.0, 
	Double_t lEEBoundHi = 100.0, 
	TString lWhichSystVar = "V0Radius" );

void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 );

//--------------------------------------------------------------
//------------------- MAIN FUNCTION ----------------------------
//--------------------------------------------------------------
void DoSystematics(
	TString Type = "XiMinus", 
	Double_t LowMult = 0.0, 
	Double_t HighMult = 100.0, 
	Double_t LowEE = 0.0, 
	Double_t HighEE = 100.0,
	Bool_t DoMult = kTRUE){

  Bool_t DoEE = !DoMult; 

  TFile* InputFile = new TFile(Form( "Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", Type.Data(), LowMult, HighMult, LowEE, HighEE ));

  TH1F* hV0Radius = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"V0Radius");
  TH1F* hCascRadius = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"CascRadius");
  TH1F* hNegToPV = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"DCANegToPV");
  TH1F* hPosToPV = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"DCAPosToPV");
  TH1F* hV0ToPV = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"DCAV0ToPV");
  TH1F* hBachToPV = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"DCABachToPV");
  TH1F* hV0Daught = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"DCAV0Daughters");
  TH1F* hCascDaught = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"DCACascDaughters");
  TH1F* hV0CosPA = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"V0CosPA");
  TH1F* hCascCosPA = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"CascCosPA");
  TH1F* hV0Mass = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"V0Mass");
  TH1F* hCompetingSpecies = 0x0;
  if (Type == "OmegaMinus" || Type == "OmegaPlus")
    hCompetingSpecies = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"CompetingSpecies");
  TH1F* hPLT = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"ProperLifetime");
  TH1F* hTPCNSigmas = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"TPCPIDNSigmas");
  TH1F* hTPCNClusters = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"TPCNClusters");
  TH1F* hSigExtBinCount = makeSystPlotsV0s(Type,"V0M",LowMult,HighMult,LowEE,HighEE,"SigExtBinCount");

  //Topological contribution systematics
  TH1F* hSystTopological = (TH1F*)hV0Radius->Clone("hSystTopological");
  hSystTopological->Reset();

  Double_t V0Radius   = 0;
  Double_t CascRadius = 0;
  Double_t NegToPV    = 0;
  Double_t PosToPV    = 0;
  Double_t BachToPV   = 0;
  Double_t V0ToPV     = 0;
  Double_t CascDaught = 0;
  Double_t CascCosPA  = 0;
  Double_t V0Daught   = 0;
  Double_t V0CosPA    = 0;

  for (int i = 1; i <= hV0Radius->GetNbinsX(); i++){

  	//If deviation > 2*sigma not included
    if ( TMath::Abs( hV0Radius->GetBinContent(i)/hV0Radius->GetBinError(i) ) > 2 ) 
    	V0Radius = hV0Radius->GetBinContent(i); 
    else V0Radius = 0;
    if ( TMath::Abs(hCascRadius->GetBinContent(i)/hCascRadius->GetBinError(i)) > 2 ) 
    	CascRadius = hCascRadius->GetBinContent(i); 
    else CascRadius = 0;   
    if ( TMath::Abs(hNegToPV->GetBinContent(i)/hNegToPV->GetBinError(i)) > 2 ) 
    	NegToPV= hNegToPV->GetBinContent(i); 
    else  NegToPV = 0;
    if ( TMath::Abs(hPosToPV->GetBinContent(i)/hPosToPV->GetBinError(i)) > 2 ) 
    	PosToPV = hPosToPV->GetBinContent(i); 
    else PosToPV = 0;
    if ( TMath::Abs(hBachToPV->GetBinContent(i)/hBachToPV->GetBinError(i)) > 2 ) 
    	BachToPV = hBachToPV->GetBinContent(i); 
    else BachToPV = 0;
    if ( TMath::Abs(hV0Daught->GetBinContent(i)/hV0Daught->GetBinError(i)) > 2 ) 
    	V0Daught = hV0Daught->GetBinContent(i); 
    else V0Daught = 0;
    if ( TMath::Abs(hV0CosPA->GetBinContent(i)/hV0CosPA->GetBinError(i)) > 2 ) 
    	V0CosPA = hV0CosPA->GetBinContent(i); 
    else V0CosPA = 0;
    if ( TMath::Abs(hV0ToPV->GetBinContent(i)/hV0ToPV->GetBinError(i)) > 2 ) 
    	V0ToPV = hV0ToPV->GetBinContent(i); 
    else V0ToPV = 0;
    if ( TMath::Abs(hCascDaught->GetBinContent(i)/hCascDaught->GetBinError(i)) > 2 ) 
    	CascDaught = hCascDaught->GetBinContent(i); 
    else CascDaught = 0;
    if ( TMath::Abs(hCascCosPA->GetBinContent(i)/hCascCosPA->GetBinError(i)) > 2 ) 
    	CascCosPA = hCascCosPA->GetBinContent(i); 
    else CascCosPA = 0;
    
    hSystTopological->SetBinContent(i, 
      TMath::Sqrt( V0Radius*V0Radius + NegToPV*NegToPV + PosToPV*PosToPV + BachToPV*BachToPV +
                   V0Daught*V0Daught + CascDaught*CascDaught + V0CosPA*V0CosPA + CascCosPA*CascCosPA +
                   V0CosPA*V0CosPA + V0ToPV*V0ToPV )    
    );
  }

  //Other selections contribution systematics
  TH1F* hSystOthers = (TH1F*)hV0Radius->Clone("hSystOthers");
  hSystOthers->Reset();

  Double_t V0Mass = 0;
  Double_t PLT = 0;
  Double_t CompetingSpecies = 0;
  Double_t TPCNClusters = 0;
  Double_t TPCNSigmas = 0;
  Double_t SigExtBinCount = 0;

  for (int i = 1; i<= hV0Radius->GetNbinsX(); i++){

  	//If deviation > 2*sigma not included
      if ( TMath::Abs( hV0Mass->GetBinContent(i)/hV0Mass->GetBinError(i)) > 2 ) 
      	V0Mass = hV0Mass->GetBinContent(i); 
      else V0Mass = 0;
      if ( TMath::Abs( hPLT->GetBinContent(i)/hPLT->GetBinError(i)) > 2 ) 
      	PLT = hPLT->GetBinContent(i); 
      else  PLT= 0;
      if (Type == "OmegaMinus" || Type == "OmegaPlus") {
        if ( TMath::Abs(hCompetingSpecies->GetBinContent(i)/hCompetingSpecies->GetBinError(i)) > 2 ) 
        	CompetingSpecies = hCompetingSpecies->GetBinContent(i); 
        else CompetingSpecies = 0;
        } 
      if ( TMath::Abs( hTPCNClusters->GetBinContent(i)/hTPCNClusters->GetBinError(i)) > 2 ) 
      	TPCNClusters = hTPCNClusters->GetBinContent(i); 
      else TPCNClusters = 0;
      if ( TMath::Abs( hTPCNSigmas->GetBinContent(i)/hTPCNSigmas->GetBinError(i)) > 2 ) 
      	TPCNSigmas = hTPCNSigmas->GetBinContent(i); 
      else TPCNSigmas = 0;
      if ( TMath::Abs( hSigExtBinCount->GetBinContent(i)/hSigExtBinCount->GetBinError(i)) > 2 ) 
      	SigExtBinCount = hSigExtBinCount->GetBinContent(i); 
      else SigExtBinCount = 0;
      
      hSystOthers->SetBinContent(i, 
      TMath::Sqrt(V0Mass*V0Mass  + CompetingSpecies*CompetingSpecies + TPCNClusters*TPCNClusters +
       TPCNSigmas*TPCNSigmas + SigExtBinCount*SigExtBinCount + PLT*PLT) 
    );
  }

  //Other contributions:
  Double_t MaterialBudget = 0.04;
  Double_t MultIndipEffic = 0.02;
  //Double_t MultIndipFeedDown = 0.02;

  //Total systematics
  TH1F* hSystTot = (TH1F*)hV0Radius->Clone("hSystTot");  
  hSystTot->Reset();
  for (int i = 1; i<= hV0Radius->GetNbinsX(); i++){
      hSystTot->SetBinContent(i, 
        TMath::Sqrt(
          hSystOthers->GetBinContent(i)*hSystOthers->GetBinContent(i) +
          hSystTopological->GetBinContent(i)*hSystTopological->GetBinContent(i) +
          MaterialBudget*MaterialBudget +
          MultIndipEffic*MultIndipEffic )
        );        
  }



  //Canvases---------------------------------------------

  //Topological Displayed
  TCanvas * Top = new TCanvas("Top","",900,900);
  TLegend* l = new TLegend (0.2,0.55,0.66,0.9);
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
  hV0ToPV->SetLineColor(kBlue+1);
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
  Top->SaveAs(Form("images/TopologicalSystematics_V0-%03.0f_%03.0f-ZDC-%03.0f_%03.0f.png",Type.Data(),LowMult, HighMult, LowEE, HighEE));

  //Selection Displayed
  TCanvas * Sel = new TCanvas("Sel","",900,900);
  TLegend* l1 = new TLegend (0.14,0.65,0.65,0.9);
  l1->AddEntry(hPLT,"Proper Life Time","L");
  l1->AddEntry(hV0Mass,"V0 Mass","L");
  if (Type == "OmegaMinus" || Type == "OmegaPlus") l1->AddEntry(hCompetingSpecies,"Competing Species","L");
  l1->AddEntry(hTPCNClusters,"TPC N of Clusters","L");
  l1->AddEntry(hTPCNSigmas,"TPC N Sigmas","L");
  l1->AddEntry(hSigExtBinCount,"Sigmas for Sgn Extraction","L");
  l1->SetTextSize(0.03);
  hPLT->SetLineColor(kRed+1);
  hV0Mass->SetLineColor(kOrange+1);
  hSigExtBinCount->SetLineColor(kSpring-1);
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
  hSigExtBinCount->SetMarkerStyle(1);
  hSigExtBinCount->SetLineWidth(2);
  hPLT->SetTitle("Contribution of Selection Variables");
  hPLT->Draw("HIST ");
  hV0Mass->Draw("HIST SAME ");
  if (Type == "OmegaMinus" || Type == "OmegaPlus") hCompetingSpecies->Draw("HIST SAME ");
  hTPCNClusters->Draw("HIST SAME ");
  hTPCNSigmas->Draw("HIST SAME ");
  hSigExtBinCount->Draw("HIST SAME ");
  l1->Draw("SAME");
  Sel->SaveAs(Form("images/%s-SelectionSystematics_V0-%03.0f_%03.0f-ZDC-%03.0f_%03.0f.png",Type.Data(),LowMult, HighMult, LowEE, HighEE));

  //Contributions Displayed
  TCanvas* cn = new TCanvas("cn","",950,900);
  Float_t constSyst = TMath::Sqrt(MaterialBudget*MaterialBudget + MultIndipEffic*MultIndipEffic);
  TLine* lconstSyst = new TLine(1.,constSyst,6.5,constSyst);
  
  TLegend* legend = new TLegend (0.14,0.7,0.6,0.9);
  legend->AddEntry(hSystTopological,"Topological systematics","L");
  legend->AddEntry(hSystOthers,"Selection cuts systematics","L");
  legend->AddEntry(lconstSyst,"p_{T} indipendent systematics","L");
  legend->AddEntry(hSystTot,"Total systematics","L");
  legend->SetTextSize(0.03);

  hSystTot->SetLineWidth(2);
  hSystTot->SetLineColor(kBlack);
  hSystOthers->SetLineWidth(2);
  hSystOthers->SetLineColor(kRed);
  hSystTopological->SetLineWidth(2);
  hSystTopological->SetLineColor(kBlue);
  lconstSyst->SetLineColor(kGreen);
  lconstSyst->SetLineWidth(2);

  hSystTot->Draw();
  hSystTot->SetYTitle("Systematics");
  hSystTot->GetYaxis()->SetRangeUser(-0.005,.1);
  hSystTot->GetYaxis()->SetTitleOffset(1.);
  hSystTot->SetTitle(Form("Systematics contributions for %s",Type.Data()));
  hSystOthers->Draw("SAME");
  hSystTopological->Draw("SAME");
  lconstSyst->Draw("SAME");
  legend->Draw("SAME");  

  cn->SaveAs(Form("images/%s-Systematics-Displayed_V0-%03.0f_%03.0f-ZDC-%03.0f_%03.0f.png",Type.Data(),LowMult, HighMult, LowEE, HighEE));

  //Get histograms from files
  TH1F* lHistPt;
  lHistPt = (TH1F *) InputFile->Get(Form("fHistPt%s", Type.Data()));
  const int nbins = lHistPt->GetNbinsX();
  TH1F* hSystPt = (TH1F*)lHistPt->Clone("hSystPt");
  for (int i = 1; i<= nbins; i++){
    hSystPt->SetBinContent(i,lHistPt->GetBinContent(i));
    hSystPt->SetBinError( i, hSystTot->GetBinContent(i) * lHistPt->GetBinContent(i));
    hSystPt->SetMarkerStyle(1);
  }

 
  //Output File
  TFile* Write = new TFile (Form("SystematicsFinalResults-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",Type.Data(),LowMult, HighMult, LowEE, HighEE), "recreate");
  lHistPt->Write();
  hSystPt->Write();
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
  if (Type == "OmegaMinus" || Type == "OmegaPlus") hCompetingSpecies->Write();
  hTPCNClusters->Write();
  hTPCNSigmas->Write();
  hSigExtBinCount->Write();

  
  

}

//-------------------------------------------------------------------------------
//------------------------- DEFINE FUNCTIONS ------------------------------------
//-------------------------------------------------------------------------------

//---------------------------------------------------------------
TH1F* makeSystPlotsV0s(
          TString lCascType = "XiMinus",
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
  

  //Definition of variables
  //Number of variation cuts
  const int nfiles = 5;  
  TH1F* hMaxDev;
  TH1F* lHistPt[nfiles];
  TH1F* lHistPtRaw[nfiles];
  TH1F* hCut[nfiles-1];
  TString lDataFilename[nfiles];
  TFile* InputFile[nfiles];   

  //Set data files
  TString lSystFile = "FilesSyst/Results-Systematics";
  lSystFile.Append( Form( "-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f-", lCascType.Data(), lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi ) );
  
  lDataFilename[0] = Form( "Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", lCascType.Data(), lMultBoundLo, lMultBoundHi, lEEBoundLo, lEEBoundHi );
  InputFile[0] = new TFile(lDataFilename[0].Data(),"READ");
  for (int i = 1; i <= nfiles-1; i++){
    lDataFilename[i] = lSystFile + lWhichSystVar + Form("-%i.root",i);
    InputFile[i] = new TFile(lDataFilename[i].Data(),"READ");
  }  

  cout << "\n---------------------------------------------------------------" << endl;
  cout<<endl<<"\n---> Do Systematics for: "<< lWhichSystVar << endl;
  cout<<endl<<"\n\n---> Number of syst files: "<< nfiles-1 << endl;
  cout << " --- Set Minimum Bias file   :  " << lDataFilename[0].Data() << endl;
  for (int i = 1; i <= nfiles-1; i++){
    cout << " --- Set Systematics file #"<< i << " :  " << lDataFilename[i].Data() << endl;
  }  
  cout << "\n---------------------------------------------------------------" << endl;


  //Inizialize Histos
  lHistPt[0] = (TH1F *) InputFile[0]->Get(Form("fHistPt%s", lCascType.Data()));
  lHistPtRaw[0] = (TH1F *) InputFile[0]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");
  hMaxDev = (TH1F*)lHistPt[0]->Clone("hMaxDev");
  hMaxDev->Reset();

  for (int i = 1; i < nfiles; i++){
    lHistPt[i] = (TH1F *) InputFile[i]->Get(Form("fHistPt%s", lCascType.Data()));
    lHistPtRaw[i] = (TH1F *) InputFile[i]->Get("lInvMassReal/lInvMassRealRawData/fHistPtRaw");      
    hCut[i-1] = (TH1F*)lHistPt[i]->Clone(Form("hCut%i",i));//will contain var cuts
  }

  //Divide following Roger-Barlow prescription
  for (int k = 1; k < nfiles; k++){
   	DivideAndComputeRogerBarlow(hCut[k-1],lHistPt[0]);
  }
        
  //Compute Max Deviation
  double binvalue[nfiles-1];
  double maxvalue = 0.;
  for (int i = 1; i<= lHistPt[0]->GetNbinsX(); i++){
    for (int k = 1; k < nfiles; k++){ 
      binvalue[k-1] = TMath::Abs(hCut[k-1]->GetBinContent(i)-1);
    }
    maxvalue = binvalue[0];
    int counter = 0;
    for (int k = 1; k < nfiles; k++){ 
      maxvalue = max(maxvalue,binvalue[k]);
      if (maxvalue == binvalue[k]) counter = k;
    }       
    hMaxDev->SetBinError(i,hCut[counter]->GetBinError(i));
    hMaxDev->SetBinContent(i,TMath::Abs(maxvalue));
    if (TMath::Abs(maxvalue) == 0)hMaxDev->SetBinContent(i,0.00000001);
  }
    

  //Prepare Canvas
  //Max Deviation
  hMaxDev->GetXaxis()->SetRangeUser(1.,6.5);
  hMaxDev->GetYaxis()->SetRangeUser(-0.0005,.1);
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
  TCanvas* maxdev = new TCanvas("maxdev"," ",1000,800);
  maxdev->SetRightMargin(0.15);
  maxdev->SetLeftMargin(0.15);
  maxdev->SetBottomMargin(0.15);
  maxdev->SetGridy();
  maxdev->SetGridx();
  hMaxDev->SetStats(kFALSE);  
  hMaxDev->SetName(lWhichSystVar);
  hMaxDev->Draw();
  maxdev->SaveAs(Form("images/%s-MaxRelDev-V0-%03.0f_%03.0f-ZDC-%03.0f_%03.0f.png",lCascType.Data(),lWhichSystVar.Data(),lMultBoundLo,lMultBoundHi,lEEBoundLo,lEEBoundHi));

  //Cut Variation
  TCanvas* cutvar = new TCanvas("cutvar"," ",1000,800);
  cutvar->SetRightMargin(0.08);
  cutvar->SetLeftMargin(0.15);
  cutvar->SetBottomMargin(0.15);
  cutvar->SetGridy();
  hCut[0]->Draw();
  TLegend* legend = new TLegend (0.8,.83,0.98,0.98);  

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
    hCut[k]->GetXaxis()->SetRangeUser(1.,6.5);
    hCut[k]->GetYaxis()->SetRangeUser(0.9,1.1);
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
  cutvar->SaveAs(Form("images/%s-CutVar-%s-V0-%03.0f_%03.0f-ZDC-%03.0f_%03.0f.png",lCascType.Data(),lWhichSystVar.Data(),lMultBoundLo,lMultBoundHi,lEEBoundLo,lEEBoundHi));

  //Return Max Dev Histo
  return hMaxDev;
}


//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop + errorfrombottom) );
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
