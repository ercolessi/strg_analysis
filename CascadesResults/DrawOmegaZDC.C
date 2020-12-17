double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p);
TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);


void DrawOmegaZDC(){

	TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);
	const int nbins = 6;
	Double_t centrality[nbins] = {0,30,50,80,90,100};
	TFile* filename[nbins];
	TH1D* lHistPt[nbins];
	TH1D* lSystPt[nbins];
	TF1* levy[nbins];
	TH1D* hClone[nbins];
	TH1D* hPlot[nbins];
	TString Type = "OmegaMinus";
	TString period = "LHC15f";


	TFile* file = new TFile("EventCountLoss_15f.root","READ");
	TH1F* hEevt = (TH1F*)file->Get("EventLoss/hevtlossZDC");
	TH1F* hsgnloss[10];
	for (int i =0; i<10;i++){
	    hsgnloss[i] = (TH1F*)file->Get(Form("SgnLoss/%s/%s/fHistptSel%s_%.0f-%.0f_%s",Type.Data(),"EEsel", "ZDC",centrality[i],centrality[i+1],Type.Data()));
	}
	double evtloss15[5];
	evtloss15[0] = (hEevt->GetBinContent(1)*1+hEevt->GetBinContent(2)*4)/5; 
	evtloss15[1] = (hEevt->GetBinContent(3)*5+hEevt->GetBinContent(4)*5)/10; 
	evtloss15[2] = (hEevt->GetBinContent(5)*5+hEevt->GetBinContent(6)*10)/15; 
	evtloss15[3] = (hEevt->GetBinContent(7)*10+hEevt->GetBinContent(8)*10)/20; 
	evtloss15[4] = (hEevt->GetBinContent(9)*20+hEevt->GetBinContent(10)*30)/50; 

	for (int i = 1; i < nbins; i++){
	     filename[i] = new TFile(Form("15f/Results-%s-13TeV-V0M_000_100_ZDC_%03.0f_%03.0f.root",Type.Data(),centrality[i-1],centrality[i]),"READ");
	         //Form("SystematicsFinalResults-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",Type.Data(), 0., 100., centrality[i-1],centrality[i]),"READ");
	    
	     lHistPt[i] = (TH1D*)filename[i]->Get(Form("fHistPt%s",Type.Data()));
	     //lHistPtA[i] = (TH1D*)filenamea[i]->Get(Form("fHistPt%s",TypeA.Data()));
	     lSystPt[i] = lHistPt[i];//(TH1D*) filename[i]->Get("hSystPt");
	     levy[i] = LevyTsallis("LevyTsallis",1.672); 
	     hClone[i] = (TH1D*)lHistPt[i]->Clone("hClone");
	     hClone[i]->Reset();
	     //hAdd[i] = (TH1D*)lHistPt[i]->Clone("hAdd");
	     //hAdd[i]->Reset();
	     for (int b = 1; b <= hClone[i]->GetNbinsX(); b++){
	        hClone[i]->SetBinContent(b,lHistPt[i]->GetBinContent(b)*
	                                 //hEevt->GetBinContent(i)/
	                                 evtloss15[i-1]/
	                                 hsgnloss[i-1]->GetBinContent(b));
	        hClone[i]->SetBinError(b,ErrorInRatio(lHistPt[i]->GetBinContent(b),lHistPt[i]->GetBinError(b),hsgnloss[i-1]->GetBinContent(b),hsgnloss[i-1]->GetBinError(b)));
	        
	     }
	}

	hClone[1]->SetLineColor(kRed+1);
	hClone[1]->SetMarkerColor(kRed+1);
	hClone[1]->SetMarkerStyle(25);
	lSystPt[1]->SetLineColor(kRed+1);
	lSystPt[1]->SetMarkerColor(kRed+1);

	hClone[2]->SetLineColor(kOrange-3);
	hClone[2]->SetMarkerColor(kOrange-3);
	lSystPt[2]->SetLineColor(kOrange-3);
	lSystPt[2]->SetMarkerColor(kOrange-3);

	hClone[3]->SetLineColor(kYellow+1);
	hClone[3]->SetMarkerColor(kYellow+1);
	lSystPt[3]->SetLineColor(kYellow+1);
	lSystPt[3]->SetMarkerColor(kYellow+1);

	hClone[4]->SetLineColor(kGreen+2);
	hClone[4]->SetMarkerColor(kGreen+2);
	lSystPt[4]->SetLineColor(kGreen+2);
	lSystPt[4]->SetMarkerColor(kGreen+2);

	hClone[5]->SetLineColor(kBlue+3);
	hClone[5]->SetMarkerColor(kBlue+3);
	lSystPt[5]->SetLineColor(kBlue+3);
	lSystPt[5]->SetMarkerColor(kBlue+3);

	for (int i = 1; i < nbins; i++){
	    hPlot[i] = (TH1D*)hClone[i]->Clone("hPlot");
	    hClone[i]->Fit(levy[i], "", "",1,5.5);
	    levy[i]->SetLineStyle(2);    
	}
	    
	//levy[0]->SetLineColor(kBlack);
	levy[1]->SetLineColor(kOrange+7);
	levy[2]->SetLineColor(kOrange-3);
	levy[3]->SetLineColor(kYellow+1);
	levy[4]->SetLineColor(kGreen+2);
	levy[5]->SetLineColor(kBlue+3);

	hPlot[1]->Scale(16);
	//hSystPt[1]->Scale(16);
	hPlot[2]->Scale(8);
	//hSystPt[2]->Scale(8);
	hPlot[3]->Scale(4);
	//hSystPt[3]->Scale(4);
	hPlot[4]->Scale(2);
	//hSystPt[4]->Scale(2);
	hPlot[5]->Scale(1);
	//hSystPt[6]->Scale(1);

	Double_t normal[nbins];
	for (int i =1; i<nbins ; i++) {normal[i] = levy[i]->GetParameter(3);}
	normal[1] *= 16.;
	normal[2] *= 8.;
	normal[3] *= 4.;
	normal[4] *= 2.;
	normal[5] *= 1.;
	for (int i =1; i<nbins ; i++) {
	    levy[i]->SetParameter(3,normal[i]);
	    levy[i]->SetLineStyle(2);
	}



	TCanvas* s = new TCanvas("s","",800,1000);
	s->SetLogy();
	s->SetGridx();
	s->SetGridy();
	TLegend* l = new TLegend(0.6,0.6,0.9,0.89);
	l->SetTextSize(0.023);
	l->SetHeader("ZDC selection (percentiles):");


	TStyle* mcStyle = new TStyle("mcStyle","Francesca's Root Styles");  
	mcStyle->SetPadTickX(1); 
	mcStyle->SetPadTickY(1); 
	mcStyle->SetPalette(1,0); 
	mcStyle->cd();
	s->SetFillColor(kWhite);
	//Set empty histo
	hPlot[1]->GetYaxis()->SetTitle("1/N_{ev} d^{2}N/(dp_{T}dy) [(GeV/c)^{-1}]");
	hPlot[1]->GetXaxis()->SetRangeUser(1.,5.5);
	hPlot[1]->GetYaxis()->SetRangeUser(6E-7,5);
	hPlot[1]->GetYaxis()->SetTitleOffset(1.9);
	hPlot[1]->SetTitle("0-1");
	hPlot[1]->SetStats(0);
	hPlot[1]->SetMarkerStyle(9);
	//hPlot[10]->Scale(2);

	TH1D* hempty = (TH1D*)hPlot[1]->Clone("hempty");
	hempty->Reset();
	hempty->SetTitle(Form("%s Corrected Spectra",Type.Data()));
	hempty->Draw();

	Double_t n[nbins];
	n[1]=16.;
	n[2]=8.;
	n[3]=4.;
	n[4]=2.;
	n[5]=1.;

	for (int i = 1; i < nbins; i++){
	    hPlot[i]->SetMarkerStyle(21);
	    //lSystPt[i]->Draw("SAME E2");
	    hPlot[i]->Draw("SAME LEP");
	    hPlot[i]->SetTitle(Form("%.0f-%.0f (x%.0f)",centrality[i-1],centrality[i],n[i]));
	    l->AddEntry(hPlot[i],"","LEP");
	    levy[i]->Draw("SAME");
	}

	l->AddEntry(levy[1],"","L");
	l->Draw("SAME");  
	// draw text
	TText *xlabel = new TText();
	xlabel-> SetNDC();
	xlabel-> SetTextFont(1);
	xlabel-> SetTextColor(1);
	xlabel-> SetTextSize(0.03);
	xlabel-> SetTextAlign(22);
	xlabel-> SetTextAngle(0);
	xlabel-> DrawText(0.35, 0.82, period.Data());

	s->SetRightMargin(0.09);
	s->SetLeftMargin(0.15);
	s->SetBottomMargin(0.15);
	//s->SetFrameFillStyle(0);

	s->Draw();
	s->SaveAs(Form("%s%sZDCOmega.png",Type.Data(),period.Data())); 


}



double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop - errorfrombottom) );
    }
    return 1.;
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

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{
  
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e6);
  return fLevyTsallis;
}