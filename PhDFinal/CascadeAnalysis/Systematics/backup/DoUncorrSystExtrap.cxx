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

void DoUncorrSystExtrap(
    TString fWhichEstimator = "V0M", 
    TString fWhichParticle = "Omega", 
    Double_t lFixedLo = 10., 
    Double_t lFixedHi = 20.){

        TString fWhichSelEstimator ="SPDClusters";
        if (fWhichEstimator.Contains("SPDClusters") )fWhichSelEstimator = "V0M";

        TString inputfilename = Form("ExtrSyst-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", 
                                fWhichParticle.Data(), fWhichEstimator.Data(), fWhichSelEstimator.Data(), lFixedLo, lFixedHi);
        
        TFile* inputfile = TFile::Open(inputfilename);
        TFile* inputfileMB = TFile::Open(Form("ExtrSyst-%s-13TeV_INELgt0.root",fWhichParticle.Data()));
        
        //
        //Percentile
        Double_t * percentile;
        Long_t percbinnumb;
        Double_t pmb[] = {0,100};
        Long_t nmb = sizeof(pmb)/sizeof(Double_t) - 1;
        Double_t p0[] = {0,5,10,15,20,30,40,50,70,100};
        Long_t n0 = sizeof(p0)/sizeof(Double_t) - 1;
        Double_t p1[] = {0,5,10,20,30,40,50,100};
        Long_t n1 = sizeof(p1)/sizeof(Double_t) - 1;
        Double_t p2[] = {0,20,30,40,50,60,70,100};
        Long_t n2 = sizeof(p2)/sizeof(Double_t) - 1; 
        Double_t p4[] = {0,5,10,20,30,40,50,100};
        Long_t n4 = sizeof(p4)/sizeof(Double_t) - 1;
        Double_t p5[] = {0,10,20,30,40,50,60,70,100};
        Long_t n5 = sizeof(p5)/sizeof(Double_t) - 1;
        Double_t pOmega[] = {0,5,10,30,50,100};
        Long_t npOmega = sizeof(pOmega)/sizeof(Double_t) - 1;
        Double_t pOmega2[] = {0,40,70,100};
        Long_t npOmega2 = sizeof(pOmega2)/sizeof(Double_t) - 1;
     
        if (fWhichSelEstimator.Contains("SPD") && lFixedLo == 0. && lFixedHi == 100.){
          if (fWhichParticle.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = npOmega;
          } else{
            percentile = p0;
            percbinnumb = n0;
          }
        }
        if (fWhichSelEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
          if (fWhichParticle.Contains("Omega")){
            percentile = pOmega;
            percbinnumb = npOmega;
          } else{
            percentile = p1;
            percbinnumb = n1;
          }
        }
        if (fWhichSelEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
           if (fWhichParticle.Contains("Omega")){
            percentile = pOmega2;
            percbinnumb = npOmega2;
          }
          percentile = p2;
          percbinnumb = n2;
        }
        if (fWhichSelEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
          percentile = p4;
          percbinnumb = n4;
        }
        if (fWhichSelEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
          percentile = p5;
          percbinnumb = n5;
        }
          
        //0 is Levy
        TH1D* YieldHistoF[6], *YieldRatioFLT[5], *YieldCorrSyst[5], *YieldMB[6],*YieldRatioMB[5];
   
        TString func[6] = {"Levy","BlastWave","Boltz","MTexpo","FermiDirac","BoseEinstein"};
        //{"Levy","BlastWave","Boltz","MTexpo","FermiDirac","BoseEinstein"};
        for (int nf = 0; nf < 6 ; nf++){
            YieldHistoF[nf] = (TH1D*)inputfile->Get(Form("%s/Yields",func[nf].Data()));       
            //
            YieldMB[nf] = (TH1D*)inputfileMB->Get(Form("%s/Yields",func[nf].Data()));    
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
        TH1D* hMaxDevTot = (TH1D*)YieldCorrSyst[0]->Clone("hSystTot");
        hMaxDevTot->Reset();     
        for (int nmult = 0; nmult < percbinnumb; nmult ++){
            double maxdev = 0;
            double valuexfile[5];
            for (int nfile = 0; nfile < 5 ; nfile++){
                valuexfile[nfile] = YieldCorrSyst[nfile]->GetBinContent(nmult+1);     
            }
            maxdev = valuexfile[0];
            int counter = 0;
            for (int k = 1; k < 5; k++){ 
                if (nmult==percbinnumb-1) cout << maxdev << "   " << valuexfile[k] << endl;
                maxdev = max(maxdev,valuexfile[k]);
                if (maxdev == valuexfile[k]) counter = k;
            }
            if (nmult==percbinnumb-1) cout << maxdev << endl;
            hMaxDevTot->SetBinContent(nmult+1, maxdev);
            //hMaxDevTot->SetBinError(nmult+1, YieldCorrSyst[counter]->GetBinError(nmult+1));
            hMaxDevTot->SetLineColor(kBlack);
            
        }
        hMaxDevTot->GetYaxis()->SetTitleOffset(2.4);
        hMaxDevTot->GetYaxis()->SetTitle("Mav. dev. correlated");
        hMaxDevTot->GetYaxis()->SetRangeUser(0.,0.1);
        hMaxDevTot->GetXaxis()->SetTitleOffset(1.3);
        hMaxDevTot->GetYaxis()->SetTitleSize(0.03);
        hMaxDevTot->GetXaxis()->SetTitleSize(0.03);
        hMaxDevTot->SetTitle("Max Dev Correlated");
        hMaxDevTot->SetMarkerStyle(20);
        hMaxDevTot->SetMarkerSize(1.5);
        hMaxDevTot->SetLineWidth(3);
        hMaxDevTot->SetLineStyle(7);
      
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

        TCanvas* c = new TCanvas("c","",1200,1000);
        c->SetLeftMargin(0.17);
        c->SetRightMargin(0.17);
        c->SetBottomMargin(0.17);
        YieldCorrSyst[0]->GetYaxis()->SetRangeUser(0.,0.15);
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
        hMaxDevTot->Draw("SAME");
        //
        TLegend* l = new TLegend(0.2,0.55,0.49,0.89);
        l->SetTextSize(0.028);
        l->SetBorderSize(0);
        l->AddEntry(hMaxDevTot,"max dev.","L");        
        for (int nf = 0; nf < 5 ; nf++){
            l->AddEntry(YieldCorrSyst[nf],func[nf+1].Data(),"LP");
        }
        l->Draw("SAME");     
         TLatex *xlabel2 = new TLatex();
       	//xlabel2->SetTextFont(42);
        xlabel2-> SetNDC();
        xlabel2-> SetTextColor(1);
        xlabel2-> SetTextSize(0.03);
        xlabel2-> SetTextAlign(22);
        xlabel2-> SetTextAngle(0);
        xlabel2-> DrawLatex(0.6, 0.81, Form("%s fixed [%.0f-%.0f]",fWhichSelEstimator.Data(),lFixedLo,lFixedHi));
        TLatex *xlabel3 = new TLatex();
        //xlabel3->SetTextFont(42);
        xlabel3-> SetNDC();
        xlabel3-> SetTextColor(1);
        xlabel3-> SetTextSize(0.05);
        xlabel3-> SetTextAlign(22);
        xlabel3-> SetTextAngle(0);
        xlabel3-> DrawLatex(0.4, 0.76, Form("#%s",fWhichParticle.Data()));
  
        c->SaveAs(Form("images/Extrap/YieldExtrapolation-%s-%s_%03.0f_%03.0f.png",fWhichParticle.Data(),fWhichSelEstimator.Data(),lFixedLo, lFixedHi));
    
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

        TCanvas* r = new TCanvas("r","",1200,1000);
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
        xlabel2-> DrawLatex(0.6, 0.81, Form("%s fixed [%.0f-%.0f]",fWhichSelEstimator.Data(),lFixedLo,lFixedHi));
        xlabel3-> DrawLatex(0.4, 0.76, Form("#%s",fWhichParticle.Data()));
        r->SaveAs(Form("images/Extrap/UncorrExtrapSyst-%s-%s_%03.0f_%03.0f.png",fWhichParticle.Data(),fWhichSelEstimator.Data(),lFixedLo, lFixedHi));
    
        //
        TH1D* hMaxDevUncorr = (TH1D*)SystFinal[0]->Clone("hSystUncorr");
        hMaxDevUncorr->Reset();  
   
        for (int nmult = 0; nmult < percbinnumb; nmult ++){
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
            hMaxDevUncorr->SetBinContent(nmult+1, maxdev);
            hMaxDevUncorr->SetBinError(nmult+1, 0.000000000001);
            hMaxDevUncorr->SetLineColor(kBlack);
            hMaxDevUncorr->SetMarkerColor(kBlack);
            hMaxDevUncorr->SetLineWidth(2);
            hMaxDevUncorr->SetMarkerStyle(24);
            hMaxDevUncorr->GetYaxis()->SetTitle("Rel. Uncorr. Systematic");
            hMaxDevUncorr->GetYaxis()->SetRangeUser(0.,0.1);
        }
        hMaxDevUncorr->GetYaxis()->SetTitleOffset(2.4);
        hMaxDevUncorr->GetXaxis()->SetTitleOffset(1.3);
        hMaxDevUncorr->GetYaxis()->SetTitleSize(0.03);
        hMaxDevUncorr->GetXaxis()->SetTitleSize(0.03);
        hMaxDevUncorr->SetTitle("");

        TCanvas* d = new TCanvas("d","",1200,1000);
        hMaxDevUncorr->GetYaxis()->SetRangeUser(0.,.2);
        hMaxDevUncorr->SetStats(0);
        hMaxDevUncorr->GetYaxis()->SetTitle("Rel. Uncorr. Systematic");
        d->SetLeftMargin(0.17);
        d->SetRightMargin(0.17);
        d->SetBottomMargin(0.17);
        hMaxDevUncorr->GetYaxis()->SetTitleOffset(1.3);
        hMaxDevUncorr->GetYaxis()->SetTitleSize(0.05);
        hMaxDevUncorr->GetXaxis()->SetTitleSize(0.05);
        hMaxDevUncorr->Draw(); 
        //hMaxDevTot->Draw("SAME");
        xlabel2-> DrawLatex(0.6, 0.81, Form("%s fixed [%.0f-%.0f]",fWhichSelEstimator.Data(),lFixedLo,lFixedHi));
        xlabel3-> DrawLatex(0.4, 0.76, Form("#%s",fWhichParticle.Data()));
        d->SaveAs(Form("images/Extrap/MaxDevExtrapSyst-%s-%s_%03.0f_%03.0f.png",fWhichParticle.Data(),fWhichSelEstimator.Data(),lFixedLo, lFixedHi));
    

        TCanvas* cmb = new TCanvas("cmb","",1200,1000);
        cmb->SetLeftMargin(0.17);
        cmb->SetRightMargin(0.17);
        cmb->SetBottomMargin(0.17);
        YieldRatioMB[0]->Draw();
        YieldRatioMB[0]->GetYaxis()->SetRangeUser(.8,1.1);
        
        for (int nf = 0; nf < 5 ; nf++){
          YieldRatioMB[nf]->Draw("SAME");
          YieldRatioFLT[nf]->Draw("SAME");
        }


        TFile* Write = new TFile(Form("Final-Extrap-Syst-%s-%s_%03.0f_%03.0f.root",fWhichParticle.Data(),fWhichSelEstimator.Data(),lFixedLo, lFixedHi), "RECREATE");
        hMaxDevTot->Write();
        hMaxDevUncorr->Write();
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
