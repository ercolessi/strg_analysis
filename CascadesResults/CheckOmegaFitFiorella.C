#include <iostream>
#include <fstream>
using namespace std;

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );

void CheckOmegaFitFiorella(){

	TFile* fileThis = new TFile("~/risultatitest/Results-OmegaPlus-13TeV-V0M_000_100_ZDC_000_100.root","READ");
	TFile* fileFior = new TFile("~/Scaricati/resultsSpectra/Results-OmegaPlus-V0M-000to100_WithV0refitAndImprovedDCA.root","READ");
	
    TF1* fitFior[6];
    TF1* fitThis[6];
    double meanFior[6], errmFior[6];
    double meanThis[6], errmThis[6];
    double sigmaFior[6], errsFior[6];
    double sigmaThis[6], errsThis[6];
    double mean[6],errmean[6],sigma[6],errsigma[6];

    for (int i = 0; i < 6; i++)
    {
    	fitFior[i] = (TF1*)fileFior->Get(Form("lInvMassReal/lInvMassRealRawData/fGausPt%i",i));
    	meanFior[i] = fitFior[i]->GetParameter(1);
    	sigmaFior[i] = fitFior[i]->GetParameter(2);
    	errmFior[i] = fitFior[i]->GetParError(1);
    	errsFior[i] = fitFior[i]->GetParError(2);

    	fitThis[i] = (TF1*)fileThis->Get(Form("lInvMassReal/lInvMassRealRawData/fGausPt%i",i));
    	meanThis[i] = fitThis[i]->GetParameter(1);
    	sigmaThis[i] = fitThis[i]->GetParameter(2);
    	errmThis[i] = fitFior[i]->GetParError(1);
    	errsThis[i] = fitFior[i]->GetParError(2);

    	cout << sigmaFior[i] << "    " << sigmaThis[i]<< endl;

    	mean[i] = (double)meanThis[i]/meanFior[i];
    	sigma[i] = (double)sigmaThis[i]/sigmaFior[i];
    	errmean[i] = ErrorInRatio(meanThis[i],errmThis[i],meanFior[i],errmFior[i]);
    	errsigma[i] = ErrorInRatio(sigmaThis[i],errsThis[i],sigmaFior[i],errsFior[i]);
    } 

    Double_t ptbinlimits[]   = {0.90, 1.60, 2.20, 2.60, 3.00, 3.80, 5.50}; 
	Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;
	Double_t ept[ptbinnumb],pt[ptbinnumb];
	for (int i = 0; i < ptbinnumb; i++){
		pt[i] = (ptbinlimits[i+1]+ptbinlimits[i])/2;
		ept[i] = (ptbinlimits[i+1]-ptbinlimits[i])/2;
	}

    TGraphErrors* meancheckThis = new TGraphErrors(6,pt,meanThis,ept,errmThis);
    TGraphErrors* meancheckFior = new TGraphErrors(6,pt,meanFior,ept,errmFior);
    TGraphErrors* sigmacheckThis = new TGraphErrors(6,pt,sigmaThis,ept,errsThis);
    TGraphErrors* sigmacheckFior = new TGraphErrors(6,pt,sigmaFior,ept,errsFior);
    TGraphErrors* sigmag = new TGraphErrors(6,pt,sigma,ept,errsigma);
    TGraphErrors* media = new TGraphErrors(6,pt,mean,ept,errmean);

    TLine* l = new TLine(0.6,1,5.9,1.);
	l->SetLineColor(kRed);
	l->SetLineStyle(2);

    TCanvas* m = new TCanvas("m","",1000,800);
	m->SetGridy();
	m->SetGridx();
	meancheckThis->SetMarkerStyle(21);
	meancheckThis->GetYaxis()->SetRangeUser(1.6,1.7);
	meancheckThis->GetYaxis()->SetTitleOffset(1.4);
	meancheckThis->GetXaxis()->SetTitle("p_{T} GeV/c");
	meancheckThis->GetYaxis()->SetTitle("mean [ Omega^{+} ]");
    meancheckFior->SetLineColor(kRed);
    meancheckFior->SetMarkerColor(kRed);
    meancheckFior->SetMarkerStyle(25);

    meancheckThis->Draw();
    meancheckFior->Draw("SAME EP");

	l->Draw("SAME");

    TCanvas* d = new TCanvas("d","",1000,800);
    d->SetGridy();
    d->SetGridx();
    sigmacheckThis->SetMarkerStyle(21);
    sigmacheckThis->GetYaxis()->SetRangeUser(0.001,0.003);
    sigmacheckThis->GetYaxis()->SetTitleOffset(1.4);
    sigmacheckThis->GetXaxis()->SetTitle("p_{T} GeV/c");
    sigmacheckThis->GetYaxis()->SetTitle("sigma [ #Omega^{+} ]");
    sigmacheckFior->SetLineColor(kRed);
    sigmacheckFior->SetMarkerColor(kRed);
    sigmacheckFior->SetMarkerStyle(25);


    sigmacheckThis->Draw();
    sigmacheckFior->Draw("SAME EP");

    l->Draw("SAME");

	TCanvas* s = new TCanvas("s","",1000,800);
	s->SetGridy();
	s->SetGridx();
	sigmag->SetMarkerStyle(21);
	sigmag->GetYaxis()->SetRangeUser(0.9,1.1);
	sigmag->GetYaxis()->SetTitleOffset(1.4);
	sigmag->GetXaxis()->SetTitle("p_{T} GeV/c");
	sigmag->GetYaxis()->SetTitle("sigma_{THIS}/sigma_{Fior} [ #Omega^{+} ]");

	sigmag->Draw();
	l->Draw("SAME");

    TCanvas* sm = new TCanvas("sm","",1000,800);
    sm->SetGridy();
    sm->SetGridx();
    media->SetMarkerStyle(21);
    media->GetYaxis()->SetRangeUser(0.95,1.05);
    media->GetXaxis()->SetRangeUser(0.9,5.5);
    media->GetYaxis()->SetTitleOffset(1.4);
    media->GetXaxis()->SetTitle("p_{T} GeV/c");
    media->GetYaxis()->SetTitle("mean_{THIS}/mean_{Fior} [ #Omega^{-} ]");

    media->Draw();
    l->Draw("SAME");
	

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