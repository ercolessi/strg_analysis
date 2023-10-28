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
void init(int lClassCode,
          vector<Double_t> &percentileSPDClusters_low,
          vector<Double_t> &percentileSPDClusters_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          int &nbins,
          Bool_t DoMB);
Double_t GetPercentilefromValue(TFile *lfilename, Int_t lRun, Int_t lValue);

// classes
enum classname
{
    kStandalone = 0,
    kHighMult,
    kLowMult,
    kHighZN,
    kLowZN,
    kVeryLowZN
};

void ComputeEventLoss(int lClassName = 0, Bool_t DoMB = kTRUE)
{

    // Percentile
    vector<Double_t> percentileSPDClusters_low;
    vector<Double_t> percentileSPDClusters_high;
    vector<Double_t> percentileV0M_low;
    vector<Double_t> percentileV0M_high;
    int nbins;
    // initialize
    init(lClassName, percentileSPDClusters_low, percentileSPDClusters_high, percentileV0M_low, percentileV0M_high, nbins, DoMB);
    const int percbinnumb = nbins;

    cout << " Class code enum: " << lClassName << endl;
    cout << " Found " << nbins << " selections:" << endl;
    for (int i = 0; i < nbins; i++)
    {
        cout << " SPDClusters [" << percentileSPDClusters_low[i] << "-" << percentileSPDClusters_high[i] << "] + V0M [" << percentileV0M_low[i] << "-" << percentileV0M_high[i] << "]" << endl;
    }
    cout << "------------------------------------------" << endl;

    //MC
    TFile *fileMC = new TFile("FullMC_15f17j18i_normcorr_2.root", "READ");
    TList* clistMC = (TList*)fileMC->Get("PWGLF_StrVsMult_MC/cList");

    Float_t fCentrality_V0M = 0.;
    Float_t fCentrality_SPDClusters = 0.;
    Int_t fRun = 0;
    Int_t fSPDtracklets0815 = 0;
    Long_t NeventsPerc[percbinnumb];

    //Variables
    TTree* lTreeEventMC;
    Float_t fMultiplicity = 0.;
    Float_t fEnergy = 0.;
    Bool_t fEvSel_AllSelections = kFALSE;
    Bool_t fEvSel_INELgtZEROtrue = kFALSE;
    Bool_t fEvSel_zVtxZMC = kFALSE;
    Bool_t fEvSel_NotPileupSPDInMultBins = kFALSE;
    Bool_t fEvSel_AcceptedVertexPosition = kFALSE;
    Bool_t fEvSel_Pileuptrue = kFALSE;

    //Histograms
    TH1F *hevtloss = new TH1F(Form("hevtloss"), Form("Event loss correction class %i; ",lClassName), percbinnumb, 0, percbinnumb);
    for (int i = 1; i <= percbinnumb; i++){
        hevtloss->GetXaxis()->SetBinLabel(i, Form("SPDClusters_%.0f-%.0f_V0M_%.0f-%.0f", percentileSPDClusters_low[i-1], percentileSPDClusters_high[i-1], percentileV0M_low[i-1], percentileV0M_high[i-1]));
    }

    cout<<"--------------- Open MC File --------------------\n"<<endl;

    lTreeEventMC = (TTree *)fileMC->Get("PWGLF_StrVsMult_MC/fTreeEvent");

    //Get from Tree
    lTreeEventMC->SetBranchAddress("fRun",&fRun);
    lTreeEventMC->SetBranchAddress("fEvSel_zVtxZMC", &fEvSel_zVtxZMC);
    lTreeEventMC->SetBranchAddress("fEvSel_AllSelections", &fEvSel_AllSelections);
    lTreeEventMC->SetBranchAddress("fEvSel_INELgtZEROtrue", &fEvSel_INELgtZEROtrue);
    lTreeEventMC->SetBranchAddress("fEvSel_NotPileupSPDInMultBins", &fEvSel_NotPileupSPDInMultBins);
    lTreeEventMC->SetBranchAddress("fCentrality_V0M", &fEnergy);
    lTreeEventMC->SetBranchAddress("fCentrality_SPDClusters", &fMultiplicity);
    lTreeEventMC->SetBranchAddress("fEvSel_AcceptedVertexPosition", &fEvSel_AcceptedVertexPosition);
    lTreeEventMC->SetBranchAddress("fEvSel_Pileuptrue", &fEvSel_Pileuptrue);

    Long_t lNEvtSelected[percbinnumb];
    Long_t lNEvttrueSelected[percbinnumb];

    for (int k=0; k < percbinnumb; k++) { // loop on percentile classes

        //initialize counters
        lNEvtSelected[k]=0;
        lNEvttrueSelected[k]=0;

        //set the limits
        double multmin, multmax, eemin, eemax;
        multmin = percentileSPDClusters_low[k];
        multmax = percentileSPDClusters_high[k];
        eemin = percentileV0M_low[k];
        eemax = percentileV0M_high[k];

        cout<<" \nWill now loop over events, please wait...\n"<<endl;
        for(Long_t iEv = 0; iEv<lTreeEventMC->GetEntries(); iEv++) {

            lTreeEventMC->GetEntry(iEv);
            if( iEv % ( lTreeEventMC->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEventMC->GetEntries()<<endl;

            if( fEvSel_AllSelections &&
                fEvSel_INELgtZEROtrue &&
                !fEvSel_Pileuptrue &&
                fEnergy>eemin && fEnergy<eemax &&
                fMultiplicity>multmin && fMultiplicity<multmax)
                lNEvtSelected[k]++;

            if (fEvSel_INELgtZEROtrue &&
                fEvSel_AcceptedVertexPosition &&
                !fEvSel_Pileuptrue &&
                fEnergy > eemin && fEnergy < eemax &&
                fMultiplicity > multmin && fMultiplicity < multmax)
                lNEvttrueSelected[k]++;
        }

        cout<<" --------------------------------------------------------"<<endl;
        cout<<" This selection......................................:" << " SPDClusters[" << percentileSPDClusters_low[k] << " - " << percentileSPDClusters_high[k] << "] + V0M[" << percentileV0M_low[k] << " - " << percentileV0M_high[k] << "] " << endl;
        cout<<"\n Number of events reco INEL>0, this selection......: "<<lNEvtSelected[k] <<endl;
        cout<<"\n Number of events true INEL>0, this selection......: "<<lNEvttrueSelected[k] <<endl;
        cout<<"\n Event Loss Correction.............................: " << (Float_t)lNEvtSelected[k]/lNEvttrueSelected[k] << endl;
        cout<<" --------------------------------------------------------"<<endl;

    } // end loop on percentile classes

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
        Write = new TFile ("EventLoss-13TeV_INELgt0.root"
                                ,"RECREATE");
    } else{
        Write = new TFile (Form("EventLoss-13TeV_class%i.root", lClassName)
                                ,"RECREATE");
    }

    TDirectoryFile* epsevt;
    epsevt = new TDirectoryFile("EventLoss","EventLoss");
    epsevt->cd();
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

void init(int lClassCode,
          vector<Double_t> &percentileSPDClusters_low,
          vector<Double_t> &percentileSPDClusters_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          int &nbins,
          Bool_t DoMB)
{

    // class 0 --> kStandalone
    Double_t percentileV0M_low_0[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
    Double_t percentileV0M_high_0[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    Double_t percentileSPDClusters_low_0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t percentileSPDClusters_high_0[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Long_t n0 = sizeof(percentileV0M_low_0) / sizeof(Double_t);

    // class 10 --> kStandalone for Omegas
    Double_t percentileV0M_low_10[] = {0, 5, 15, 30, 50} ;
    Double_t percentileV0M_high_10[] = {5, 15, 30, 50, 100};
    Double_t percentileSPDClusters_low_10[] = {0, 0, 0, 0, 0};
    Double_t percentileSPDClusters_high_10[] = {100, 100, 100, 100, 100};
    Long_t n10 = sizeof(percentileV0M_low_10) / sizeof(Double_t);

    // class 1 --> kHighMult
    Double_t percentileSPDClusters_low_1[] = {10, 10, 10, 10, 10, 10, 10};
    Double_t percentileSPDClusters_high_1[] = {20, 20, 20, 20, 20, 20, 20};
    Double_t percentileV0M_low_1[] = {0, 5, 10, 20, 30, 40, 50};
    Double_t percentileV0M_high_1[] = {5, 10, 20, 30, 40, 50, 100};
    Long_t n1 = sizeof(percentileSPDClusters_low_1) / sizeof(Double_t);

    // class 2 --> kLowMult
    Double_t percentileSPDClusters_low_2[] = {40, 40, 40, 40, 40, 40, 40};
    Double_t percentileSPDClusters_high_2[] = {50, 50, 50, 50, 50, 50, 50};
    Double_t percentileV0M_low_2[] = {0, 20, 30, 40, 50, 60, 70};
    Double_t percentileV0M_high_2[] = {20, 30, 40, 50, 60, 70, 100};
    Long_t n2 = sizeof(percentileSPDClusters_low_2) / sizeof(Double_t);

    // class 3 --> kHighZN
    Double_t percentileSPDClusters_low_3[] = {10, 40, 60, 70, 80};
    Double_t percentileSPDClusters_high_3[] = {40, 60, 70, 80, 100};
    Double_t percentileV0M_low_3[] = {70, 60, 40, 40, 40};
    Double_t percentileV0M_high_3[] = {100, 100, 100, 80, 70};
    Long_t n3 = sizeof(percentileSPDClusters_low_3) / sizeof(Double_t);

    // class 4 --> kLowZN
    Double_t percentileSPDClusters_low_4[] = {0, 10, 20, 30, 50};
    Double_t percentileSPDClusters_high_4[] = {20, 30, 40, 50, 100};
    Double_t percentileV0M_low_4[] = {40, 30, 30, 20, 0};
    Double_t percentileV0M_high_4[] = {60, 70, 50, 50, 30};
    Long_t n4 = sizeof(percentileSPDClusters_low_4) / sizeof(Double_t);

    // class 5 --> fixed very low ZN
    Double_t percentileSPDClusters_low_5[] = {0, 10, 20, 30};
    Double_t percentileSPDClusters_high_5[] = {10, 20, 30, 50};
    Double_t percentileV0M_low_5[] = {20, 10, 0, 0};
    Double_t percentileV0M_high_5[] = {30, 30, 20, 10};
    Long_t n5 = sizeof(percentileSPDClusters_low_5) / sizeof(Double_t);

    Double_t percentileSPDClusters_low_mb[] = {0};
    Double_t percentileSPDClusters_high_mb[] = {100};
    Double_t percentileV0M_low_mb[] = {0};
    Double_t percentileV0M_high_mb[] = {100};
    Long_t nmb = sizeof(percentileSPDClusters_low_mb) / sizeof(Double_t);

    if (DoMB){
        nbins = nmb;
        for (Int_t i = 0; i < nmb; i++)
        {
            percentileSPDClusters_low.push_back(percentileSPDClusters_low_mb[i]);
            percentileSPDClusters_high.push_back(percentileSPDClusters_high_mb[i]);
            percentileV0M_low.push_back(percentileV0M_low_mb[i]);
            percentileV0M_high.push_back(percentileV0M_high_mb[i]);
        }
        cout << "------------------------------------------" << endl;
        cout << " Initializing MB class..." << endl;
    } else {

        if (lClassCode == kStandalone)
        {
            nbins = n0;
            for (Int_t i = 0; i < n0; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_0[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_0[i]);
                percentileV0M_low.push_back(percentileV0M_low_0[i]);
                percentileV0M_high.push_back(percentileV0M_high_0[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kStandalone class..." << endl;
        }

        if (lClassCode == 10)
        {
            nbins = n10;
            for (Int_t i = 0; i < n10; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_10[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_10[i]);
                percentileV0M_low.push_back(percentileV0M_low_10[i]);
                percentileV0M_high.push_back(percentileV0M_high_10[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kStandalone class for Omegas..." << endl;
        }

        if (lClassCode == kHighMult)
        {
            nbins = n1;
            for (Int_t i = 0; i < n1; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_1[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_1[i]);
                percentileV0M_low.push_back(percentileV0M_low_1[i]);
                percentileV0M_high.push_back(percentileV0M_high_1[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kHighMult class..." << endl;
        }

        if (lClassCode == kLowMult)
        {
            nbins = n2;
            for (Int_t i = 0; i < n2; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_2[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_2[i]);
                percentileV0M_low.push_back(percentileV0M_low_2[i]);
                percentileV0M_high.push_back(percentileV0M_high_2[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kLowMult class..." << endl;
        }

        if (lClassCode == kHighZN)
        {
            nbins = n3;
            for (Int_t i = 0; i < n3; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_3[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_3[i]);
                percentileV0M_low.push_back(percentileV0M_low_3[i]);
                percentileV0M_high.push_back(percentileV0M_high_3[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kHighZN class..." << endl;
        }

        if (lClassCode == kLowZN)
        {
            nbins = n4;
            for (Int_t i = 0; i < n4; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_4[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_4[i]);
                percentileV0M_low.push_back(percentileV0M_low_4[i]);
                percentileV0M_high.push_back(percentileV0M_high_4[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kLowZN class..." << endl;
        }

        if (lClassCode == kVeryLowZN)
        {
            nbins = n5;
            for (Int_t i = 0; i < n5; i++)
            {
                percentileSPDClusters_low.push_back(percentileSPDClusters_low_5[i]);
                percentileSPDClusters_high.push_back(percentileSPDClusters_high_5[i]);
                percentileV0M_low.push_back(percentileV0M_low_5[i]);
                percentileV0M_high.push_back(percentileV0M_high_5[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kVeryLowZN class..." << endl;
        }
    }
}

Double_t GetPercentilefromValue(TFile *lfilename, Int_t lRun, Int_t lValue)
{ // set in the run func

    // Getting percentile
    TString name = Form("hcum%s_%i", "SPDClusters", lRun);

    if (lRun == 274594 || lRun == 288897){
        return 0.;
    }

    TH1D *hcum = (TH1D *)lfilename->Get(name);
    if (!hcum) cout << "Run " << lRun << " not available!" << endl;

    Double_t percentile = 100 * (hcum->GetBinContent(hcum->FindBin(lValue)));

    return percentile;
}