using namespace std;
#include <iostream>
#include <fstream>

#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLine.h>
#include <TFile.h>
#include <TTree.h>

double ErrorInRatioCorr ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
void ErrorInRatioHisto ( TH1D* h1, TH1D *h2 );
Double_t DoWeightedMean(Bool_t err , Double_t h1, Double_t h2, Double_t h3, Double_t w1, Double_t w2, Double_t w3);

void ComputeNormCorr(
    Bool_t DoEventLoss = kTRUE,
    Bool_t DoSignalLoss = kTRUE,
    TString lCascType = "Lambda",
    Double_t lFixedLo = 10.0,
    Double_t lFixedHi = 20.0,
    TString lWhichVarEstimator = "SPDClusters",
    TString lWhichFixedEstimator = "V0M",
    Bool_t lDoMB = kFALSE)
{

    cout << "Please check the percentile binning is fine before running!!!\n" << endl;

    TString lEnergyEstimator = "CentV0M";

    TString particle = "", antiparticle = "";
    if (lCascType.Contains("Xi")){
        particle = "XiMinus";
        antiparticle = "XiPlus";
    }
    if (lCascType.Contains("Omega")){
        particle = "OmegaMinus";
        antiparticle = "OmegaPlus";
    }
    if (lCascType.Contains("Lambda")){
        particle = "Lambda";
        antiparticle = "AntiLambda";
    }
    if (lCascType.Contains("K0Short")){
        particle = "K0Short";
        antiparticle = "K0Short";
    }

    //Percentile
    Double_t * percentile;
    Long_t percbinnumb;
	Double_t pmb[] = {0,100};
    Long_t nmb = sizeof(pmb)/sizeof(Double_t) - 1;
    Double_t p0[] = {0,1,5,10,15,20,30,40,50,70,100};
    Long_t n0 = sizeof(p0)/sizeof(Double_t) - 1;
    Double_t p1[] = {0,5,10,20,30,40,50,100};
    Long_t n1 = sizeof(p1)/sizeof(Double_t) - 1;
    Double_t p2[] = {0,20,30,40,50,60,70,100};
    Long_t n2 = sizeof(p2)/sizeof(Double_t) - 1;
    Double_t p4[] = {0,5,10,20,30,40,50,100};
    Long_t n4 = sizeof(p4)/sizeof(Double_t) - 1;
    Double_t p5[] = {0,10,20,30,40,50,60,70,100};
    Long_t n5 = sizeof(p5)/sizeof(Double_t) - 1;
	Double_t pOmega[] = {0, 5, 15, 30, 50, 100};
	Long_t npOmega = sizeof(pOmega) / sizeof(Double_t) - 1;
	Double_t pOmega1[] = {0, 5, 10, 30, 50, 100};
	Long_t npOmega1 = sizeof(pOmega1) / sizeof(Double_t) - 1;
	Double_t pOmega2[] = {0, 40, 70, 100};
	Long_t npOmega2 = sizeof(pOmega2) / sizeof(Double_t) - 1;
	Double_t pOmega3[] = {0, 5, 10, 30, 100};
	Long_t npOmega3 = sizeof(pOmega3) / sizeof(Double_t) - 1;
	Double_t pOmega4[] = {0, 30, 50, 100};
	Long_t npOmega4 = sizeof(pOmega4) / sizeof(Double_t) - 1;

	if (lDoMB){
		if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 0. && lFixedHi == 100.){
			percentile = pmb;
			percbinnumb = nmb;
		}
	} else{
		if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 0. && lFixedHi == 100.){
			if (lCascType.Contains("Omega")){
				percentile = pOmega;
				percbinnumb = npOmega;
			} else{
				percentile = p0;
				percbinnumb = n0;
			}
		}
		if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
			if (lCascType.Contains("Omega")){
				percentile = pOmega;
				percbinnumb = npOmega;
			} else{
				percentile = p1;
				percbinnumb = n1;
			}
		}
		if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
			percentile = p2;
			percbinnumb = n2;
		}
		if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
			percentile = p4;
			percbinnumb = n4;
		}
		if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
			percentile = p5;
			percbinnumb = n5;
		}
	}

    //pT bins
    Double_t* ptbinlimits;
    Long_t ptbinnumb;
    Double_t ptbinlimitsXi[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
    Long_t ptbinnumbXi = sizeof(ptbinlimitsXi)/sizeof(Double_t) - 1;
    Double_t ptbinlimitsOmega[] = {0.90, 1.60, 2.20, 2.60, 3.00, 3.80, 5.50};
    Long_t ptbinnumbOmega = sizeof(ptbinlimitsOmega)/sizeof(Double_t) - 1;
    Double_t ptbinlimitsLambda[] = {0.4,  0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.5, 2.9, 3.4, 4, 5, 6.5, 8,10.};
    Long_t ptbinnumbLambda = sizeof(ptbinlimitsLambda)/sizeof(Double_t) - 1;
    Double_t ptbinlimitsK0s[] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3.0,
                                    3.3, 3.6, 3.9, 4.2, 4.6, 5, 5.4, 5.9, 6.5, 7, 7.5, 8, 8.5, 9.2, 10, 11, 12, 13.5, 15};
    Long_t ptbinnumbK0s = sizeof(ptbinlimitsK0s) / sizeof(Double_t) - 1;

    //MC
    TFile* fileMC[2];
    fileMC[0] = new TFile("/home/fercoles/strg_analysis/Paper/data/mc/LHC15g3a3_2.root","READ");
    fileMC[1] = new TFile("/home/fercoles/strg_analysis/Paper/data/mc/LHC15g3c3.root","READ");
    //
    TList* clistMC[2];
    clistMC[0] = (TList*)fileMC[0]->Get("PWGLF_StrVsMult_MC/cList");
    clistMC[1] = (TList*)fileMC[1]->Get("PWGLF_StrVsMult/cList");

    TTree* lTreeEvent;
    Float_t fCentrality_V0M = 0.;
    Float_t fCentrality_SPDClusters = 0.;
    TTree* lTreeEventMC;
    Float_t fMultiplicity = 0.;
    Float_t fEnergy = 0.;
    Bool_t fEvSel_AllSelections = kFALSE;
    Bool_t fEvSel_INELgtZEROtrue = kFALSE;
    Bool_t fEvSel_zVtxZMC = kFALSE;
    Int_t fRun;

    TH1D* fHistptSel[percbinnumb][2];
    TH1D* fHistpttrueSel[percbinnumb][2];
    TH3D* h3DGenSel[2];
    TH3D* h3DGenSelAntiPart[2];
    TH3D* h3DGentrueSel[2];
    TH3D* h3DGentrueSelAntiPart[2];
    TH1D* fHistGentrueSel[percbinnumb][2];
    TH1D* fHistGenSel[percbinnumb][2];


    if (DoSignalLoss){

        for (int iset = 0; iset < 2; iset ++){ //loop over 3 datasets

            cout<<"--------------- Open MC File " << iset << " --------------------\n"<<endl;

            Double_t temppt = 0, temppt2 = 0;
            Int_t minmultbin = 0, maxmultbin = 0, mineebin = 0, maxeebin = 0, minmultbintrue = 0, maxmultbintrue = 0, mineebintrue = 0, maxeebintrue = 0;

            // INEL > 0
            h3DGenSel[iset] = (TH3D*)clistMC[iset]->FindObject(Form("fHistGeneratedPtVsSPDclVs%s%s", lEnergyEstimator.Data(), particle.Data()));
            h3DGenSel[iset]->Sumw2();
            if (!lCascType.Contains("K0Short")){
                h3DGenSelAntiPart[iset] = (TH3D *)clistMC[iset]->FindObject(Form("fHistGeneratedPtVsSPDclVs%s%s", lEnergyEstimator.Data(), antiparticle.Data()));
                h3DGenSelAntiPart[iset]->Sumw2();
            }
            // true INEL > 0
            h3DGentrueSel[iset] = (TH3D*)clistMC[iset]->FindObject(Form("fHistPtVsSPDclVs%s_Gen%s", lEnergyEstimator.Data(), particle.Data()));
            h3DGentrueSel[iset]->Sumw2();
            if (!lCascType.Contains("K0Short")){
                h3DGentrueSelAntiPart[iset] = (TH3D*)clistMC[iset]->FindObject(Form("fHistPtVsSPDclVs%s_Gen%s", lEnergyEstimator.Data(), antiparticle.Data()));
                h3DGentrueSelAntiPart[iset]->Sumw2();
            }

            // Add particle + antiparticle
            if (!lCascType.Contains("K0Short")){
                h3DGenSel[iset]->Add(h3DGenSelAntiPart[iset]);
                h3DGentrueSel[iset]->Add(h3DGentrueSelAntiPart[iset]);
            }

            for (int k=0; k < percbinnumb; k++) { // loop on percentile classes

                //set the limits
                double multmin, multmax, eemin, eemax;
                //
                if (lWhichVarEstimator.Contains("SPD")){
                    multmin = percentile[k];
                    multmax = percentile[k+1];
                    eemin = lFixedLo;
                    eemax = lFixedHi;
                }
                //
                if (lWhichVarEstimator.Contains("V0M") || lWhichVarEstimator.Contains("ZDC")){
                    eemin = percentile[k];
                    eemax = percentile[k+1];
                    multmin = lFixedLo;
                    multmax = lFixedHi;
                }

                minmultbin = h3DGenSel[iset]->GetYaxis()->FindBin( multmin );
                maxmultbin = h3DGenSel[iset]->GetYaxis()->FindBin( multmax );
                mineebin = h3DGenSel[iset]->GetZaxis()->FindBin( eemin );
                maxeebin = h3DGenSel[iset]->GetZaxis()->FindBin( eemax );
                //
                minmultbintrue = h3DGentrueSel[iset]->GetYaxis()->FindBin( multmin );
                maxmultbintrue = h3DGentrueSel[iset]->GetYaxis()->FindBin( multmax );
                mineebintrue = h3DGentrueSel[iset]->GetZaxis()->FindBin( eemin );
                maxeebintrue = h3DGentrueSel[iset]->GetZaxis()->FindBin( eemax );

                fHistGenSel[k][iset] = (TH1D*)h3DGenSel[iset]->ProjectionX(Form("fHistGenSel%i%s",k,lCascType.Data()),
                        minmultbin, maxmultbin, mineebin, maxeebin );
                fHistGentrueSel[k][iset] = (TH1D*)h3DGentrueSel[iset]->ProjectionX(Form("fHistGenSelTrue%i%s",k,lCascType.Data()),
                        minmultbintrue, maxmultbintrue, mineebintrue, maxeebintrue );

                fHistptSel[k][iset] = (TH1D*)fHistGenSel[k][iset]->Clone(Form("hsgnloss_%i-%i_MC%i",(int)percentile[k],(int)percentile[k+1],iset));
                //new TH1D(Form("hsgnloss_%i-%i_MC%i",(int)percentile[k],(int)percentile[k+1],iset),"Signal loss corr in mult bins;#it{p}_{T} (GeV/#it{c});Counts", ptbinnumb, ptbinlimits);
                fHistpttrueSel[k][iset] = (TH1D*)fHistGentrueSel[k][iset]->Clone(Form("hsgnlosstrue_%i_%i_MC%i",(int)percentile[k],(int)percentile[k+1],iset));
                //= new TH1D(Form("hsgnlosstrue_%i_%i_MC%i",(int)percentile[k],(int)percentile[k+1],iset),"Signal loss corr in mult bins;#it{p}_{T} (GeV/#it{c});Counts", ptbinnumb, ptbinlimits);
            }

            //Rebinning
            /*temppt = 0;
            for (int k = 0; k<percbinnumb; k++){
                for(long i = 1; i<fHistGenSel[k][iset]->GetNbinsX()+1;i++){
                    temppt = fHistGenSel[k][iset]->GetXaxis()->GetBinCenter(i);
                    for(long filling = 0; filling<fHistGenSel[k][iset]->GetBinContent(i); filling++){
                        fHistptSel[k][iset]->Fill(temppt);
                    }
                }
            }

            temppt2 = 0;
            for (int k = 0; k<percbinnumb; k++){
                for(long i = 1; i<fHistGentrueSel[k][iset]->GetNbinsX()+1;i++){
                    temppt2 = fHistGentrueSel[k][iset]->GetXaxis()->GetBinCenter(i);
                    for(long filling = 0; filling<fHistGentrueSel[k][iset]->GetBinContent(i); filling++){
                        fHistpttrueSel[k][iset]->Fill(temppt2);
                    }
                }
            }*/

            for (int k = 0; k<percbinnumb; k++){
                ErrorInRatioHisto(fHistptSel[k][iset],fHistpttrueSel[k][iset]);
            }
        }
    }

    ////////////////////////////////////////////////////////
    // -------------- Write on file ------------------------
    ////////////////////////////////////////////////////////

    TFile* Write;
    if (lDoMB) {
        Write = new TFile(Form("SystNorm-%s-13TeV_INELgt0.root", lCascType.Data()), "RECREATE");
    } else {
        Write = new TFile(Form("SystNorm-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", lCascType.Data(), lWhichVarEstimator.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi), "RECREATE");
    }

    for (int iset = 0; iset < 2; iset ++){ //loop over 3 datasets
        for (int k = 0; k<percbinnumb; k++){
            fHistptSel[k][iset]->Write();
        }
    }
}

//---------------------------------------------------------------------------------------------------
double ErrorInRatioCorr ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0 ) {
        Double_t lSigmaDelta = TMath::Sqrt( TMath::Abs( TMath::Power(Aerr,2) - TMath::Power(Berr,2) ) );
        return lSigmaDelta/B;
    }
    return 1.;
}

//---------------------------------------------------------------------------------------------------
void ErrorInRatioHisto ( TH1D* h1, TH1D *h2 ){
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
    if ( h2->GetBinContent(i) != 0 ){
      lSigmaDelta[i] /= h2->GetBinContent(i);
    }else{
      lSigmaDelta[i] = 0;
    }
  }
  //Replace Erorrs
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){
      if (h2->GetBinContent(i)!=0){
        h1->SetBinContent(i, h1->GetBinContent(i)/h2->GetBinContent(i));
      } /*else if (h1->GetBinContent(i)==0){
          h1->SetBinContent(i, 1);
      }else {
        cout << "HELP" << endl;
      };*/
      h1->SetBinError(i, lSigmaDelta[i]);
  }
}

Double_t DoWeightedMean(Bool_t err , Double_t h1, Double_t h2, Double_t h3, Double_t w1, Double_t w2, Double_t w3) {
    //#eventi nel period i = 1./wi^2
    //wi = 1./sqrt(eventi)

    Double_t wmean,  werror;

    Double_t num = h1/(w1*w1) + h2/(w2*w2) + h3/(w3*w3);
    Double_t den = (1./(w1*w1)) + (1./(w2*w2)) + (1./(w3*w3));

    wmean = num/den;
    werror = TMath::Sqrt(1./((1./(w1*w1)) + (1./(w2*w2)) + (1./(w3*w3))));

    if (err) {
        return werror;
    } else {
        return wmean;
    }


}