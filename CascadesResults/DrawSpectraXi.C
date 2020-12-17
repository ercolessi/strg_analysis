double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p);
TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);



void DrawSpectraXi(){


	TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);
	const int nbins = 11;
	Double_t centrality[nbins] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.};
	//{0,1,5,10,15,20,30,40,50,70,100};
	
	
	TFile* filename[nbins];
	TH1D* lHistPt[nbins];
	TH1D* lSystPt[nbins];
	TF1* levy[nbins];
	TH1D* hClone[nbins];
	TH1D* hPlot[nbins];
	TString Type = "XiPlus";
	TString period = "LHC15f";


	TFile* file = new TFile("EventCountLoss_15f.root","READ");
	TH1F* hEevt = (TH1F*)file->Get("EventLoss/hevtlossZDC");
	TH1F* hsgnloss[10];
	for (int i =0; i<10;i++){
	    hsgnloss[i] = (TH1F*)file->Get(Form("SgnLoss/%s/%s/fHistptSel%s_%.0f-%.0f_%s",Type.Data(),"EEsel", "ZDC",centrality[i],centrality[i+1],Type.Data()));
	}

for (int i = 1; i < nbins; i++){
     filename[i] = new TFile(Form("15f/Results-%s-13TeV-V0M_000_100_ZDC_%03.0f_%03.0f.root",Type.Data(),centrality[i-1],centrality[i]),"READ");
     	//Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_000_100.root",Type.Data(),centrality[i-1],centrality[i]),"READ");
     	
     	
         //Form("SystematicsFinalResults-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",Type.Data(), 0., 100., centrality[i-1],centrality[i]),"READ");
     lHistPt[i] = (TH1D*)filename[i]->Get(Form("fHistPt%s",Type.Data()));
     lSystPt[i] = lHistPt[i];//(TH1D*) filename[i]->Get("hSystPt");
     levy[i] = LevyTsallis("LevyTsallis",1.31486); 
     hClone[i] = (TH1D*)lHistPt[i]->Clone("hClone");
     hClone[i]->Reset();
     double sgn = 0,sgne=0;
     for (int b = 1; b <= hClone[i]->GetNbinsX(); b++){
     	if (i==2){
     		sgn = hsgnloss[i-2]->GetBinContent(b);
     		sgne = hsgnloss[i-2]->GetBinError(b);
     	}
     	else {
     		sgn = hsgnloss[i-1]->GetBinContent(b);
     		sgne = hsgnloss[i-1]->GetBinError(b);
     	}
        hClone[i]->SetBinContent(b,lHistPt[i]->GetBinContent(b)*
                                 hEevt->GetBinContent(i)/
                                 sgn);
        hClone[i]->SetBinError(b,ErrorInRatio(lHistPt[i]->GetBinContent(b),lHistPt[i]->GetBinError(b),sgn,sgne));
    }
}

hClone[1]->SetLineColor(kRed+1);
hClone[1]->SetMarkerColor(kRed+1);

hClone[1]->SetMarkerStyle(25);
lSystPt[1]->SetLineColor(kRed+1);
lSystPt[1]->SetMarkerColor(kRed+1);
lSystPt[1]->Scale(64.);

hClone[2]->SetLineColor(kRed-4);
hClone[2]->SetMarkerColor(kRed-4);

lSystPt[2]->SetLineColor(kRed-4);
lSystPt[2]->SetMarkerColor(kRed-4);
lSystPt[2]->Scale(32.);

hClone[3]->SetLineColor(kOrange+7);
hClone[3]->SetMarkerColor(kOrange+7);

lSystPt[3]->SetLineColor(kOrange+7);
lSystPt[3]->SetMarkerColor(kOrange+7);
lSystPt[3]->Scale(16.);

hClone[4]->SetLineColor(kOrange-3);
hClone[4]->SetMarkerColor(kOrange-3);

lSystPt[4]->SetLineColor(kOrange-3);
lSystPt[4]->SetMarkerColor(kOrange-3);
lSystPt[4]->Scale(8.);

hClone[5]->SetLineColor(kYellow+1);
hClone[5]->SetMarkerColor(kYellow+1);

lSystPt[5]->SetLineColor(kYellow+1);
lSystPt[5]->SetMarkerColor(kYellow+1);
lSystPt[5]->Scale(4.);

hClone[6]->SetLineColor(kSpring-7);
hClone[6]->SetMarkerColor(kSpring-7);
lSystPt[6]->SetLineColor(kSpring-7);
lSystPt[6]->SetMarkerColor(kSpring-7);
lSystPt[6]->Scale(2.);

hClone[7]->SetLineColor(kGreen+2);
hClone[7]->SetMarkerColor(kGreen+2);

lSystPt[7]->SetLineColor(kGreen+2);
lSystPt[7]->SetMarkerColor(kGreen+2);
lSystPt[7]->Scale(1.);

hClone[8]->SetLineColor(kAzure+8);
hClone[8]->SetMarkerColor(kAzure+8);

lSystPt[8]->SetLineColor(kAzure+8);
lSystPt[8]->SetMarkerColor(kAzure+8);
lSystPt[8]->Scale(0.5);

hClone[9]->SetLineColor(kBlue-4);
hClone[9]->SetMarkerColor(kBlue-4);

lSystPt[9]->SetLineColor(kBlue-4);
lSystPt[9]->SetMarkerColor(kBlue-4);
lSystPt[9]->Scale(0.25);

hClone[10]->SetLineColor(kBlue+3);
hClone[10]->SetMarkerColor(kBlue+3);

lSystPt[10]->SetLineColor(kBlue+3);
lSystPt[10]->SetMarkerColor(kBlue+3);
lSystPt[10]->Scale(0.10);

for (int i = 1; i < nbins; i++){
    hPlot[i] = (TH1D*)hClone[i]->Clone("hPlot");
    hClone[i]->Fit(levy[i], "", "",0.6,6.5);
    levy[i]->SetLineStyle(2);    
}

hPlot[1]->Scale(64);
//hSystPt[1]->Scale(16);
hPlot[2]->Scale(32);
//hSystPt[2]->Scale(8);
hPlot[3]->Scale(16);
//hSystPt[3]->Scale(4);
hPlot[4]->Scale(8);
//hSystPt[4]->Scale(2);
hPlot[5]->Scale(4);
hPlot[6]->Scale(2);
hPlot[7]->Scale(1);
hPlot[8]->Scale(0.5);
hPlot[9]->Scale(0.25);
hPlot[10]->Scale(0.10);


Double_t normal[nbins];
for (int i =1; i<nbins ; i++) {normal[i] = levy[i]->GetParameter(3);}
normal[1] *= 64.;
normal[2] *= 32.;
normal[3] *= 16.;
normal[4] *= 8.;
normal[5] *= 4.;
normal[6] *= 2.;
normal[7] *= 1.;
normal[8] *= 0.5;
normal[9] *= 0.25;
normal[10] *= 0.10;
for (int i =1; i<nbins ; i++) {
    levy[i]->SetParameter(3,normal[i]);
    levy[i]->SetLineStyle(2);
}

//levy[0]->SetLineColor(kBlack);
levy[1]->SetLineColor(kRed+1);
levy[2]->SetLineColor(kRed-4);
levy[3]->SetLineColor(kOrange+7);
levy[4]->SetLineColor(kOrange-3);
levy[5]->SetLineColor(kYellow+1);
levy[6]->SetLineColor(kSpring-7);
levy[7]->SetLineColor(kGreen+2);
levy[8]->SetLineColor(kAzure+8);
levy[9]->SetLineColor(kBlue-4);
levy[10]->SetLineColor(kBlue+3);

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
hPlot[1]->GetXaxis()->SetRangeUser(-0.0001,6.5);
hPlot[1]->GetYaxis()->SetRangeUser(5E-8,10000);
hPlot[1]->GetYaxis()->SetTitleOffset(1.9);
hPlot[1]->SetTitle("0-1");
hPlot[1]->SetStats(0);
hPlot[1]->SetMarkerStyle(9);

TH1D* hempty = (TH1D*)hPlot[1]->Clone("hempty");
hempty->Reset();
hempty->SetTitle(Form("%s Corrected Spectra",Type.Data()));
hempty->Draw();

Double_t n[nbins];
n[0]=0;
n[1]=64.;
n[2]=32.;
n[3]=16.;
n[4]=8.;
n[5]=4.;
n[6]=2.;
n[7]=1.;
n[8]=0.5;
n[9]=0.25;
n[10]=0.10;

for (int i = 1; i < nbins; i++){
    hPlot[i]->SetMarkerStyle(21);
    //lSystPt[i]->Draw("SAME E2");
    hPlot[i]->Draw("SAME LEP");
    hPlot[i]->SetTitle(Form("%.0f-%.0f (x%.2f)",centrality[i-1],centrality[i],n[i]));
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

//s->Draw();
s->SaveAs(Form("%s%sZDC.png",Type.Data(),period.Data())); 

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