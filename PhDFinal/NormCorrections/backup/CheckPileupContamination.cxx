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

double ErrorBinomial ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
void ErrorInRatioHisto ( TH1D* h1, TH1D *h2 );
TH1F* DoWeightedMean(TH1F* h1, TH1F* h2, TH1F* h3, Double_t w1, Double_t w2, Double_t w3);
Float_t DoWeightedMeanXBin(Bool_t err , Float_t h1, Float_t h2, Float_t h3, Float_t e1, Float_t e2, Float_t e3, Float_t w1, Float_t w2, Float_t w3);


void CheckPileupContamination(
    TString lCascType = "Xi",
    Double_t lFixedLo = 0.0,
    Double_t lFixedHi = 100.0,
    TString lWhichVarEstimator = "V0M",
    TString lWhichFixedEstimator = "SPDClusters",
    Bool_t DoMB = kFALSE)
{

    cout << "Please check the percentile binning is fine before running!!!\n\n" << endl;

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
    if (lCascType.Contains("K0Short")) {
        particle = "K0Short";
        antiparticle = "";
    }

    //Percentile
    Float_t * percentile;
    Long_t percbinnumb;

    Float_t pstd[] = {0,100};
    Long_t nstd = sizeof(pstd)/sizeof(Float_t) - 1;
    Float_t p0[] = {0,1,5,10,15,20,30,40,50,70,100};
    Long_t n0 = sizeof(p0)/sizeof(Float_t) - 1;
    Float_t p1[] = {0,5,10,20,30,40,50,100};
    Long_t n1 = sizeof(p1)/sizeof(Float_t) - 1;
    Float_t p2[] = {0,20,30,40,50,60,70,100};
    Long_t n2 = sizeof(p2)/sizeof(Float_t) - 1;
    Float_t p4[] = {0,5,10,20,30,40,50,100};
    Long_t n4 = sizeof(p4)/sizeof(Float_t) - 1;
    Float_t p5[] = {0,10,20,30,40,50,60,70,100};
    Long_t n5 = sizeof(p5)/sizeof(Float_t) - 1;
    Float_t pOmega[] = {0,5,15,30,50,100};
    Long_t nOmega = sizeof(pOmega)/sizeof(Float_t) - 1;
    Float_t pOmega1[] = {0,5,10,30,50,100};
    Long_t nOmega1 = sizeof(pOmega1)/sizeof(Float_t) - 1;
    Float_t pOmega2[] = {0,40,70,100};
    Long_t nOmega2 = sizeof(pOmega2)/sizeof(Float_t) - 1;
    Float_t pOmega3[] = {0,5,10,30,100};
    Long_t nOmega3 = sizeof(pOmega3)/sizeof(Float_t) - 1;
    Float_t pOmega4[] = {0,30,50,100};
    Long_t nOmega4 = sizeof(pOmega4)/sizeof(Float_t) - 1;


    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 0. && lFixedHi == 100.){
        if (DoMB){
            percentile = pstd;
            percbinnumb = nstd;
        }
        else{
            if (lCascType.Contains("Omega")){
                percentile = pOmega;
                percbinnumb = nOmega;
            } else{
                percentile = p0;
                percbinnumb = n0;
            }
        }
    }
    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 10. && lFixedHi == 20.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega1;
            percbinnumb = nOmega1;
        } else{
            percentile = p1;
            percbinnumb = n1;
        }
    }
    if (lWhichFixedEstimator.Contains("SPD") && lFixedLo == 40. && lFixedHi == 50.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega2;
            percbinnumb = nOmega2;
        }else {
        percentile = p2;
        percbinnumb = n2;
        }
    }
    if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 10. && lFixedHi == 20.){
        if (lCascType.Contains("Omega")){
            percentile = pOmega3;
            percbinnumb = nOmega3;
        } else{
            percentile = p4;
            percbinnumb = n4;
        }
    }
    if (lWhichFixedEstimator.Contains("V0M") && lFixedLo == 40. && lFixedHi == 50.){
         if (lCascType.Contains("Omega")){
            percentile = pOmega4;
            percbinnumb = nOmega4;
        } else{
        percentile = p5;
        percbinnumb = n5;
        }
    }


    //pT bins
    Double_t* ptbinlimits;
    Long_t ptbinnumb;
    Double_t ptbinlimitsXi[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
    Long_t ptbinnumbXi = sizeof(ptbinlimitsXi)/sizeof(Double_t) - 1;
    Double_t ptbinlimitsOmega[] = {0.90, 1.60, 2.20, 2.60, 3.00, 3.80, 5.50, 8.00, 10.0 };
    Long_t ptbinnumbOmega = sizeof(ptbinlimitsOmega)/sizeof(Double_t) - 1;
    Double_t ptbinlimitsLambda[] = {0.4,  0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.5, 2.9, 3.4, 4, 5, 6.5, 8, 10};
    Long_t ptbinnumbLambda = sizeof(ptbinlimitsLambda)/sizeof(Double_t) - 1;
    Double_t ptbinlimitsK0S[] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3.0,
                                    3.3, 3.6, 3.9, 4.2, 4.6, 5, 5.4, 5.9, 6.5, 7, 7.5, 8, 8.5, 9.2, 10, 11, 12, 13.5, 15};
    Long_t ptbinnumbK0S = sizeof(ptbinlimitsK0S) / sizeof(Double_t) - 1;

    //MC
    TFile *fileMC = new TFile("testsignalloss.root", "READ");

    TTree* lTreeEvent;
    Float_t fCentrality_V0M = 0.;
    Float_t fCentrality_SPDClusters = 0.;

    //Variables
    TTree* lTreeEventMC;
    Float_t fMultiplicity = 0.;
    Float_t fEnergy = 0.;
    Bool_t fEvSel_AllSelections = kFALSE;
    Bool_t fEvSel_INELgtZEROtrue = kFALSE;
    Bool_t fEvSel_zVtxZMC = kFALSE;
    Bool_t fEvSel_NotPileupSPDInMultBins = kFALSE;
    Bool_t fEvSel_AcceptedVertexPosition = kFALSE;
    Bool_t fEvSel_Pileuptrue = kFALSE; // FIXME
    Bool_t fEvSel_NoInconsSPDandTrackVrtx = kFALSE;
    Int_t fRun = 0.;
    Float_t fTimePileup = 0.;

    TH1F* hevtloss = new TH1F(Form("hevtloss"), Form(";%s percentile",lWhichVarEstimator.Data()), percbinnumb, percentile);

    lTreeEventMC = (TTree*)fileMC->Get("PWGLF_StrVsMult_MC/fTreeEvent");

    //Get from Tree
    lTreeEventMC->SetBranchAddress("fRun",&fRun);
    lTreeEventMC->SetBranchAddress("fEvSel_zVtxZMC", &fEvSel_zVtxZMC);
    lTreeEventMC->SetBranchAddress("fEvSel_AllSelections", &fEvSel_AllSelections);
    lTreeEventMC->SetBranchAddress("fEvSel_INELgtZEROtrue", &fEvSel_INELgtZEROtrue);
    lTreeEventMC->SetBranchAddress("fEvSel_AcceptedVertexPosition", &fEvSel_AcceptedVertexPosition);
    lTreeEventMC->SetBranchAddress("fEvSel_NotPileupSPDInMultBins", &fEvSel_NotPileupSPDInMultBins);
    lTreeEventMC->SetBranchAddress("fEvSel_Pileuptrue", &fEvSel_Pileuptrue); // FIXME
    lTreeEventMC->SetBranchAddress("fEvSel_NoInconsSPDandTrackVrtx", &fEvSel_NoInconsSPDandTrackVrtx); // FIXME
    lTreeEventMC->SetBranchAddress("fCentrality_SPDClusters", &fMultiplicity);
    lTreeEventMC->SetBranchAddress("fCentrality_V0M", &fEnergy);
    lTreeEventMC->SetBranchAddress("fTimePileup", &fTimePileup);

    Long_t lNEvtSelected[percbinnumb];
    Long_t lNEvttrueSelected[percbinnumb];

    if (lCascType.Contains("Xi")) {
        ptbinlimits = ptbinlimitsXi;
        ptbinnumb = ptbinnumbXi;
    }
    if (lCascType.Contains("Omega")) {
        ptbinlimits = ptbinlimitsOmega;
        ptbinnumb = ptbinnumbOmega;
    }
    if (lCascType.Contains("Lambda")) {
        ptbinlimits = ptbinlimitsLambda;
        ptbinnumb = ptbinnumbLambda;
    }
    if (lCascType.Contains("K0")) {
        ptbinlimits = ptbinlimitsK0S;
        ptbinnumb = ptbinnumbK0S;
    }

    for (int k=0; k < percbinnumb; k++) { // loop on percentile classes
        //initialize counters
        lNEvtSelected[k]=0;
        lNEvttrueSelected[k]=0;
    }


    cout<<" \nWill now loop over events, please wait...\n"<<endl;
    for(Long_t iEv = 0; iEv<lTreeEventMC->GetEntries(); iEv++) {

        lTreeEventMC->GetEntry(iEv);
        if( iEv % ( lTreeEventMC->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEventMC->GetEntries()<<endl;

        if (fEvSel_AcceptedVertexPosition == 0) continue;
        if (fEvSel_INELgtZEROtrue == 0) continue;
        //if (fEvSel_Pileuptrue) continue;
        //if (fabs(fTimePileup) > 3.) continue; //Reject OOB pileup
        //if (fabs(fTimePileup) < 3. && fEvSel_Pileuptrue) continue; //Reject in-bunch pileup

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

            if (fEnergy > eemin && fEnergy < eemax &&
                fMultiplicity > multmin && fMultiplicity < multmax)
                lNEvttrueSelected[k]++;

            if (fEvSel_AllSelections &&
                fEnergy > eemin && fEnergy < eemax &&
                fMultiplicity > multmin && fMultiplicity < multmax)
                lNEvtSelected[k]++;
        }
    }

    for (int k=0; k < percbinnumb; k++) { // loop on percentile classes
        cout << " --------------------------------------------------------" << endl;
        cout << " This selection....................................: " << lWhichVarEstimator.Data() << " [" << percentile[k] << " - " << percentile[k + 1] << "], " << lWhichFixedEstimator.Data() << " [" << lFixedLo << " - " << lFixedHi << "]" << endl;
        cout << "\n Events reco selection, this selection......: " << lNEvtSelected[k] << endl;
        cout << "\n Events true selection, this selection......: " << lNEvttrueSelected[k] << endl;
        cout << "\n Event Loss Correction.............................: " << (Float_t)lNEvtSelected[k] / lNEvttrueSelected[k] << endl;
        cout << " --------------------------------------------------------" << endl;
    }

    for (int i=1; i <= hevtloss->GetNbinsX(); i++){

        hevtloss->SetBinContent(i,(Float_t)lNEvtSelected[i-1]/lNEvttrueSelected[i-1]);
        hevtloss->SetBinError(i,ErrorBinomial(
            (Float_t)lNEvtSelected[i-1],(Float_t)TMath::Sqrt(lNEvtSelected[i-1]),
            (Float_t)lNEvttrueSelected[i-1],(Float_t)TMath::Sqrt(lNEvttrueSelected[i-1]))
            );
    }

    ////////////////////////////////////////////////////////
    // -------------- Write on file ------------------------
    ////////////////////////////////////////////////////////

    TFile* Write = 0x0;
    if (DoMB) {
        Write = new TFile(Form("TESTEventLoss-%s-13TeV_INELgt0.root", lCascType.Data()), "RECREATE");
    } else{
        Write = new TFile (Form("TESTEventLoss-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", lCascType.Data(), lWhichVarEstimator.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi)
                                ,"RECREATE");
    }

    hevtloss->Write();
    Write->cd();
}

//---------------------------------------------------------------------------------------------------
double ErrorBinomial ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    float eff = A/B;
    return TMath::Sqrt((eff*(1-eff))/B);
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

    h1->Divide(h1,h2,1,1,"B");
}

TH1F* DoWeightedMean(TH1F* h1, TH1F* h2, TH1F* h3, Double_t w1, Double_t w2, Double_t w3) {

    Float_t nh1 = h1->GetNbinsX();
    Float_t nh2 = h2->GetNbinsX();

    if( nh1 != nh2 ){
        cout<<"Problem! Number of bins doesn't match! "<<endl;
    }

    TH1F* h = (TH1F*)h1->Clone("Ratio");
    h->Reset();

    for (int bin = 1; bin <= nh1; bin ++){

        Float_t num = h1->GetBinContent(bin)*w1 + h2->GetBinContent(bin)*w2 + h3->GetBinContent(bin)*w3;
        Float_t den = w1 + w2 + w3;

        Float_t wmean,  werror;
        wmean = num/den;
        werror = 1./den * (TMath::Sqrt(w1*w1*TMath::Power(h1->GetBinError(bin),2) + w2*w2*TMath::Power(h2->GetBinError(bin),2) + w3*w3*TMath::Power(h3->GetBinError(bin),2)));
    }

    return h;
}


Float_t DoWeightedMeanXBin(Bool_t err , Float_t h1, Float_t h2, Float_t h3, Float_t e1, Float_t e2, Float_t e3, Float_t w1, Float_t w2, Float_t w3) {

    Float_t wmean,  werror;

    Float_t num = h1*w1 + h2*w2 + h3*w3;
    Float_t den = w1 + w2 + w3;

    wmean = num/den;
    werror =  1./den * (TMath::Sqrt(w1*w1*TMath::Power(e1,2) + w2*w2*TMath::Power(e2,2) + w3*w3*TMath::Power(e3,2)));

    if (err) {
        return werror;
    } else {
        return wmean;
    }

}