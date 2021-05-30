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
//#define max(a,b) (a>b ? a : b)

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
void DivideAndComputeRogerBarlow( TH1D* h1, TH1D *h2 );
double PassRogerBarlowCriterion(int nsigmas, TH1D* h, int bin);

void DoUncorrSystExtrapMergedBins(
    TString fWhichEstimator = "ZDC", 
    TString fWhichParticle = "Xi", 
    Double_t lLoMult = 70., 
    Double_t lHiMult = 100., 
    Double_t lLoEE = 0., 
    Double_t lHiEE = 100.){

        TString inputfilename = Form("ExtSyst-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE);
        TFile* inputfile = TFile::Open(inputfilename);
        TString fWhichOtherEstimator = "V0M";
        if (fWhichEstimator.Contains("V0M")) fWhichOtherEstimator = "ZDC";
        //
        Double_t ptbinlimits[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
        Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;
         //
        Double_t percentileV0[] = {0.,5,10,15,20,30,40,50,70,100};
        const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
        //
        Double_t percentileZDC[] = {0,20,40,50,60,70,80,90,100};
        const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;
        //
        Double_t * mult, *systpercentile;
        Int_t *systcounter;
        //
        // Choose scenario and initialize correct variables
        int tempbinnumber = 0;
        if (fWhichEstimator.Contains("V0M")) {
            tempbinnumber = nbinV0;
            mult = percentileV0;
        } 
        else if (fWhichEstimator.Contains("ZDC")) {
            tempbinnumber = nbinZDC;
            mult = percentileZDC;
        } 
        else {cout << "No valid name for estimator... its V0M or ZDC" << endl; return;}
        const int multbinnumb = tempbinnumber;

        //0 is Levy
        TH1D* YieldHistoF[6], *YieldRatioFLT[5], *YieldCorrSyst[5], *YieldMB[6],*YieldRatioMB[5];

        TString func[6] = {"Levy","BlastWave","Boltz","MTexpo","FermiDirac","BoseEinstein"};
        for (int nf = 0; nf < 6 ; nf++){
            YieldHistoF[nf] = (TH1D*)inputfile->Get(Form("%s/Yields",func[nf].Data()));       
            //
            YieldMB[nf] = (TH1D*)inputfile->Get(Form("YieldMB%s",func[nf].Data()));    
        }
        //Do Correlated Contribution fist

        for (int nf = 0; nf < 5 ; nf++){
            
            //Correlated contribution
            YieldCorrSyst[nf] = (TH1D*)YieldHistoF[nf+1]->Clone(Form("YieldCorrSyst%i",nf));
            YieldCorrSyst[nf]->Reset();
            YieldCorrSyst[nf]->GetYaxis()->SetTitle("| dN/dy^{varied-func} - dN/dy^{LT} |  / dN/dy^{LT}");
            YieldCorrSyst[nf]->GetYaxis()->SetTitle(Form("%s percentile %",fWhichEstimator.Data()));
            for (int bin = 1; bin <=YieldHistoF[nf+1]->GetNbinsX(); bin ++ ){
              YieldCorrSyst[nf]->SetBinContent(bin, TMath::Abs(YieldHistoF[nf+1]->GetBinContent(bin) - YieldHistoF[0]->GetBinContent(bin))/ YieldHistoF[0]->GetBinContent(bin));          
              YieldCorrSyst[nf]->SetBinError(bin, 0.00000001);
            }
            YieldCorrSyst[nf]->SetLineWidth(2);
            YieldCorrSyst[nf]->SetMarkerStyle(20);
            YieldCorrSyst[nf]->SetMarkerSize(1.5);

            //Uncorrelated Numerator
            YieldRatioFLT[nf] = (TH1D*)YieldHistoF[nf+1]->Clone(Form("YieldRatioFLT%i",nf));
            YieldRatioFLT[nf]->Reset();
            //YieldRatioFLT[nf]->GetYaxis()->SetTitle("| dN/dy^{varied-func} - dN/dy^{LT} |  / dN/dy^{LT}");
            YieldRatioFLT[nf]->GetYaxis()->SetTitle("ZDC percentile %");
            for (int bin = 1; bin <=YieldHistoF[nf+1]->GetNbinsX(); bin ++ ){
              YieldRatioFLT[nf]->SetBinContent(bin, YieldHistoF[nf+1]->GetBinContent(bin)/ YieldHistoF[0]->GetBinContent(bin));
              YieldRatioFLT[nf]->SetBinError(bin, 0.00000001);
            }      
            YieldRatioFLT[nf]->SetLineWidth(2);
            YieldRatioFLT[nf]->SetMarkerStyle(20);
            YieldRatioFLT[nf]->SetMarkerSize(1.5);

            //Uncorrelated Denominator
            YieldRatioMB[nf] = (TH1D*)YieldMB[nf+1]->Clone(Form("YieldRatioMBFLT%i",nf));
            YieldRatioMB[nf]->Reset();
            for (int bin = 1; bin <=YieldMB[nf+1]->GetNbinsX(); bin ++ ){
              YieldRatioMB[nf]->SetBinContent(bin,YieldMB[nf+1]->GetBinContent(bin) / YieldMB[0]->GetBinContent(bin));
              YieldRatioMB[nf]->SetBinError(bin,0.0000001);
            }
            YieldRatioMB[nf]->SetLineWidth(2);
            YieldRatioMB[nf]->SetMarkerStyle(20);
            YieldRatioMB[nf]->SetMarkerSize(1.5);
        }

        //Correlated contribution------------------------------------------------------------------
        TH1D* hMaxDevCorr = (TH1D*)YieldCorrSyst[0]->Clone("hSystCorr");
        hMaxDevCorr->Reset();     
        for (int nmult = 0; nmult < multbinnumb; nmult ++){
            double maxdev = 0;
            double valuexfile[5];
            for (int nfile = 0; nfile < 5 ; nfile++){
                valuexfile[nfile] = YieldCorrSyst[nfile]->GetBinContent(nmult+1);     
            }
            maxdev = valuexfile[0];
            int counter = 0;
            for (int k = 1; k < 5; k++){ 
                if (nmult==multbinnumb-1) cout << maxdev << "   " << valuexfile[k] << endl;
                maxdev = max(maxdev,valuexfile[k]);
                if (maxdev == valuexfile[k]) counter = k;
            }
            if (nmult==multbinnumb-1) cout << maxdev << endl;
            hMaxDevCorr->SetBinContent(nmult+1, maxdev);
            //hMaxDevCorr->SetBinError(nmult+1, YieldCorrSyst[counter]->GetBinError(nmult+1));
            hMaxDevCorr->SetLineColor(kBlack);
            
        }
        hMaxDevCorr->GetYaxis()->SetTitleOffset(2.4);
        hMaxDevCorr->GetYaxis()->SetTitle("Mav. dev. correlated");
        hMaxDevCorr->GetYaxis()->SetRangeUser(0.,0.2);
        hMaxDevCorr->GetXaxis()->SetTitleOffset(1.3);
        hMaxDevCorr->GetYaxis()->SetTitleSize(0.03);
        hMaxDevCorr->GetXaxis()->SetTitleSize(0.03);
        hMaxDevCorr->SetTitle("Max Dev Correlated");
        hMaxDevCorr->SetMarkerStyle(20);
        hMaxDevCorr->SetMarkerSize(1.5);
        hMaxDevCorr->SetLineWidth(3);
        hMaxDevCorr->SetLineStyle(7);
        hMaxDevCorr->SetName("hSystCorr");

        YieldCorrSyst[1]->SetLineColor(kOrange);
        YieldCorrSyst[0]->SetLineColor(kYellow+3);
        YieldCorrSyst[2]->SetLineColor(kGreen+3);
        YieldCorrSyst[3]->SetLineColor(kBlue);
        YieldCorrSyst[4]->SetLineColor(kRed);
        //YieldCorrSyst[0]->SetLineColor(kBlack);
        YieldCorrSyst[1]->SetMarkerColor(kOrange);
        YieldCorrSyst[0]->SetMarkerColor(kYellow+3);
        YieldCorrSyst[2]->SetMarkerColor(kGreen+3);
        YieldCorrSyst[3]->SetMarkerColor(kBlue);
        YieldCorrSyst[4]->SetMarkerColor(kRed);
        //YieldCorrSyst[0]->SetMarkerColor(kBlack);

        TCanvas* c = new TCanvas("c","",1500,1000);
        c->SetLeftMargin(0.17);
        c->SetRightMargin(0.17);
        c->SetBottomMargin(0.17);
        YieldCorrSyst[0]->GetYaxis()->SetRangeUser(0.,0.2);
        YieldCorrSyst[0]->SetStats(0);
        YieldCorrSyst[0]->GetYaxis()->SetTitle("#frac{| #LT dN/dy #GT ^{varied-func} - #LT dN/dy #GT ^{LT} |}{#LT dN/dy #GT^{LT}}");
        YieldCorrSyst[0]->GetYaxis()->SetTitleOffset(1.3);
        YieldCorrSyst[0]->GetYaxis()->SetTitleSize(0.05);
        YieldCorrSyst[0]->GetXaxis()->SetTitleSize(0.05);
        YieldCorrSyst[0]->Draw();      

        for (int nf = 0; nf < 5 ; nf++){
          
            YieldCorrSyst[nf]->SetMarkerStyle(20);
            YieldCorrSyst[nf]->Draw("SAME");
          
        } 
        hMaxDevCorr->Draw("SAME");
        //
        TLegend* l = new TLegend(0.2,0.55,0.49,0.89);
        l->SetTextSize(0.028);
        l->SetBorderSize(0);
        l->AddEntry(hMaxDevCorr,"max dev.","L");        
        for (int nf = 0; nf < 5 ; nf++){
            l->AddEntry(YieldCorrSyst[nf],func[nf+1].Data(),"LP");
        }
        l->Draw("SAME");       
          TLatex *xlabel2 = new TLatex();
       	xlabel2->SetTextFont(42);
        xlabel2-> SetNDC();
        xlabel2-> SetTextColor(1);
        xlabel2-> SetTextSize(0.03);
        xlabel2-> SetTextAlign(22);
        xlabel2-> SetTextAngle(0);
        if (fWhichEstimator.Contains("V0M")) {xlabel2-> DrawLatex(0.7, 0.81, Form("ZDC fixed [%.0f-%.0f] ",lLoEE,lHiEE));}
        else {xlabel2-> DrawLatex(0.7, 0.81, Form("V0M fixed [%.0f-%.0f]",lLoMult,lHiMult));}
  

        c->SaveAs(Form("ExtrapCorrSyst-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.png",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE));
      

        //Uncorrelated contribution
        YieldRatioFLT[1]->SetLineColor(kOrange);
        YieldRatioFLT[0]->SetLineColor(kYellow+3);
        YieldRatioFLT[2]->SetLineColor(kGreen+3);
        YieldRatioFLT[3]->SetLineColor(kBlue);
        YieldRatioFLT[4]->SetLineColor(kRed);
        //YieldRatioFLT[0]->SetLineColor(kBlack);
        YieldRatioFLT[1]->SetMarkerColor(kOrange);
        YieldRatioFLT[0]->SetMarkerColor(kYellow+3);
        YieldRatioFLT[2]->SetMarkerColor(kGreen+3);
        YieldRatioFLT[3]->SetMarkerColor(kBlue);
        YieldRatioFLT[4]->SetMarkerColor(kRed);
        //YieldRatioFLT[0]->SetMarkerColor(kBlack);

        YieldRatioMB[1]->SetLineColor(kOrange);
        YieldRatioMB[0]->SetLineColor(kYellow+3);
        YieldRatioMB[2]->SetLineColor(kGreen+3);
        YieldRatioMB[3]->SetLineColor(kBlue);
        YieldRatioMB[4]->SetLineColor(kRed);
        //YieldRatioMB[0]->SetLineColor(kBlack);
        YieldRatioMB[1]->SetMarkerColor(kOrange);
        YieldRatioMB[0]->SetMarkerColor(kYellow+3);
        YieldRatioMB[2]->SetMarkerColor(kGreen+3);
        YieldRatioMB[3]->SetMarkerColor(kBlue);
        YieldRatioMB[4]->SetMarkerColor(kRed);
  

        TH1D* SystFinal[5];
        for (int nf = 0; nf < 5 ; nf++){
            SystFinal[nf] = (TH1D*)YieldHistoF[nf+1]->Clone(Form("YieldRatioFLT%i",nf));
            SystFinal[nf]->Reset();
            for(int b =1; b <=YieldRatioFLT[0]->GetNbinsX(); b++ ){
                SystFinal[nf]->SetBinContent(b,YieldRatioFLT[nf]->GetBinContent(b)/YieldRatioMB[nf]->GetBinContent(1));
                SystFinal[nf]->SetBinError(b,0.0000000001);
            }
            SystFinal[nf]->SetMarkerStyle(20);
            SystFinal[nf]->SetMarkerSize(1.5);
        }

        SystFinal[1]->SetLineColor(kOrange);
        SystFinal[0]->SetLineColor(kYellow+3);
        SystFinal[2]->SetLineColor(kGreen+3);
        SystFinal[3]->SetLineColor(kBlue);
        SystFinal[4]->SetLineColor(kRed);
        //SystFinal[0]->SetLineColor(kBlack);
        SystFinal[1]->SetMarkerColor(kOrange);
        SystFinal[0]->SetMarkerColor(kYellow+3);
        SystFinal[2]->SetMarkerColor(kGreen+3);
        SystFinal[3]->SetMarkerColor(kBlue);
        SystFinal[4]->SetMarkerColor(kRed);

        TCanvas* r = new TCanvas("r","",1500,1000);
        SystFinal[0]->GetYaxis()->SetRangeUser(0.8,1.2);
        SystFinal[0]->SetStats(0);
        SystFinal[0]->GetYaxis()->SetTitle("#frac{[dN/dy^{varied-func}/dN/dy^{LT}]_{sel}}{[dN/dy^{varied-func}/dN/dy^{LT}]_{MB}}");
        r->SetLeftMargin(0.17);
        r->SetRightMargin(0.17);
        r->SetBottomMargin(0.17);
        SystFinal[0]->GetYaxis()->SetTitleOffset(1.3);
        SystFinal[0]->GetYaxis()->SetTitleSize(0.05);
        SystFinal[0]->GetXaxis()->SetTitleSize(0.05);
        SystFinal[0]->Draw();      

        for (int nf = 0; nf < 5 ; nf++){
            SystFinal[nf]->SetLineWidth(2);
            SystFinal[nf]->Draw("SAME");          
        }

        TLegend* lu = new TLegend(0.2,0.65,0.49,0.89);
        lu->SetTextSize(0.028);
        lu->SetBorderSize(0);    
        for (int nf = 0; nf < 5 ; nf++){
            lu->AddEntry(SystFinal[nf],func[nf+1].Data(),"LP");
        }
        lu->Draw("SAME"); 
            if (fWhichEstimator.Contains("V0M")) {xlabel2-> DrawLatex(0.7, 0.81, Form("ZDC fixed [%.0f-%.0f] ",lLoEE,lHiEE));}
        else {xlabel2-> DrawLatex(0.7, 0.81, Form("V0M fixed [%.0f-%.0f]",lLoMult,lHiMult));}
  

        r->SaveAs(Form("ContrExtrapUncorrSyst-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.png",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE));
      
      
        
        //
        TH1D* hMaxDev = (TH1D*)SystFinal[0]->Clone("hSystUncorr");
        hMaxDev->Reset();  
   
        for (int nmult = 0; nmult < multbinnumb; nmult ++){
            double maxdev = 0;
            double valuexfile[5];
            for (int nfile = 0; nfile < 5 ; nfile++){
                valuexfile[nfile] = TMath::Abs(SystFinal[nfile]->GetBinContent(nmult+1)-1);      
            }
            maxdev = valuexfile[0];
            int counter = 0;
            for (int k = 1; k < 5; k++){ 
                maxdev = max(maxdev,valuexfile[k]);
                if (maxdev == valuexfile[k]) counter = k;
            }
            hMaxDev->SetBinContent(nmult+1, maxdev);
            hMaxDev->SetBinError(nmult+1, 0.000000000001);
            hMaxDev->SetLineColor(kBlack);
            hMaxDev->SetMarkerColor(kBlack);
            hMaxDev->SetLineWidth(2);
            hMaxDev->SetMarkerStyle(24);
            hMaxDev->GetYaxis()->SetTitle("Rel. Uncorr. Systematic");
            hMaxDev->GetYaxis()->SetRangeUser(0.,0.2);
        }
        hMaxDev->GetYaxis()->SetTitleOffset(2.4);
        hMaxDev->GetXaxis()->SetTitleOffset(1.3);
        hMaxDev->GetYaxis()->SetTitleSize(0.03);
        hMaxDev->GetXaxis()->SetTitleSize(0.03);
        hMaxDev->SetTitle("");

        TCanvas* d = new TCanvas("d","",1500,1000);
        hMaxDev->GetYaxis()->SetRangeUser(0.,.2);
        hMaxDev->SetStats(0);
        hMaxDev->GetYaxis()->SetTitle("Rel. Uncorr. Systematic");
        d->SetLeftMargin(0.17);
        d->SetRightMargin(0.17);
        d->SetBottomMargin(0.17);
        hMaxDev->GetYaxis()->SetTitleOffset(1.3);
        hMaxDev->GetYaxis()->SetTitleSize(0.05);
        hMaxDev->GetXaxis()->SetTitleSize(0.05);
        hMaxDev->Draw(); 
              if (fWhichEstimator.Contains("V0M")) {xlabel2-> DrawLatex(0.4, 0.81, Form("ZDC fixed [%.0f-%.0f] ",lLoEE,lHiEE));}
        else {xlabel2-> DrawLatex(0.6, 0.81, Form("V0M fixed [%.0f-%.0f]",lLoMult,lHiMult));}
  

        d->SaveAs(Form("TotExtrapUncorrSyst-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.png",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE));
       

        TFile* Write = new TFile(Form("Final-Extrap-Syst-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root",fWhichEstimator.Data(),lLoMult,lHiMult,lLoEE,lHiEE), "RECREATE");
        for (int nmult = 0; nmult < 5; nmult ++){
        YieldRatioFLT[nmult]->Write();
        YieldRatioMB[nmult]->Write();
        }
        hMaxDevCorr->Write();
        hMaxDev->Write();

}

//----------------------------------------------------------------------------------------------------
void DivideAndComputeRogerBarlow( TH1D* h1, TH1D *h2 ){ 
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
double PassRogerBarlowCriterion(int nsigmas, TH1D* h, int bin){
  double dev = h->GetBinContent(bin);
  double RBsigma = h->GetBinError(bin);

  if (RBsigma>0 &&  (dev/RBsigma)>nsigmas ) {return dev;}
  else {return 1.;}
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