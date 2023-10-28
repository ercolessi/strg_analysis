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
          vector<Double_t> &percentileSPDtrk0815_low,
          vector<Double_t> &percentileSPDtrk0815_high,
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
    kLowZN
};

void CheckEvents(int lClassName = kHighMult, Bool_t DoMB = kFALSE)
{

    //MC
    TFile *fileMC = new TFile("FullMC_15f17j18i_normcorr.root", "READ");
    TList* clistMC = (TList*)fileMC->Get("PWGLF_StrVsMult_MC/cList");

    Float_t fCentrality_V0M = 0.;
    Float_t fCentrality_SPDtrk0815 = 0.;
    Int_t fRun = 0;
    Int_t fSPDtracklets0815 = 0;

    TFile *Read0815 = TFile::Open("/home/fercoles/strg_analysis/PhDWork/Selections/PercentileCalibration_MC.root");

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
    TH1F *hSPDtrkdisttrue = new TH1F(Form("hSPDtrkdisttrue"), "", 100,0,100);
    TH1F *hSPDtrkdistreco = new TH1F(Form("hSPDtrkdistreco"), "", 100, 0, 100);

    cout<<"--------------- Open MC File --------------------\n"<<endl;

    lTreeEventMC = (TTree *)fileMC->Get("PWGLF_StrVsMult_MC/fTreeEvent");

    //Get from Tree
    lTreeEventMC->SetBranchAddress("fRun",&fRun);
    lTreeEventMC->SetBranchAddress("fEvSel_zVtxZMC", &fEvSel_zVtxZMC);
    lTreeEventMC->SetBranchAddress("fEvSel_AllSelections", &fEvSel_AllSelections);
    lTreeEventMC->SetBranchAddress("fEvSel_INELgtZEROtrue", &fEvSel_INELgtZEROtrue);
    lTreeEventMC->SetBranchAddress("fEvSel_NotPileupSPDInMultBins", &fEvSel_NotPileupSPDInMultBins);
    lTreeEventMC->SetBranchAddress("fCentrality_V0M", &fEnergy);
    lTreeEventMC->SetBranchAddress("fEvSel_AcceptedVertexPosition", &fEvSel_AcceptedVertexPosition);
    lTreeEventMC->SetBranchAddress("fEvSel_Pileuptrue", &fEvSel_Pileuptrue);
    lTreeEventMC->SetBranchAddress("fSPDtracklets0815", &fSPDtracklets0815);

    cout<<" \nWill now loop over events, please wait...\n"<<endl;
    for(Long_t iEv = 0; iEv<lTreeEventMC->GetEntries()/10; iEv++) {

        lTreeEventMC->GetEntry(iEv);
        if( iEv % ( lTreeEventMC->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEventMC->GetEntries()<<endl;

        fMultiplicity = GetPercentilefromValue(Read0815, fRun, fSPDtracklets0815);

        if (fEvSel_AllSelections &&
            fEvSel_AcceptedVertexPosition &&
            !fEvSel_Pileuptrue &&
            fEnergy > 0 && fEnergy < 5)
            hSPDtrkdisttrue->Fill(fSPDtracklets0815);
    }

    new TCanvas;
    hSPDtrkdisttrue->Draw();
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
          vector<Double_t> &percentileSPDtrk0815_low,
          vector<Double_t> &percentileSPDtrk0815_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          int &nbins,
          Bool_t DoMB)
{

    // kStandalone
    Double_t percentileV0M_low_0[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
    Double_t percentileV0M_high_0[] = {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    Double_t percentileSPDtrk0815_low_0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t percentileSPDtrk0815_high_0[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Long_t n0 = sizeof(percentileV0M_low_0) / sizeof(Double_t);

    // class 8 --> kHighMult
    Double_t percentileSPDtrk0815_low_1[] = {0, 0, 10, 20, 50};
    Double_t percentileSPDtrk0815_high_1[] = {20, 30, 40, 70, 70};
    Double_t percentileV0M_low_1[] = {40, 40, 30, 20, 10};
    Double_t percentileV0M_high_1[] = {100, 50, 40, 30, 20};
    Long_t n1 = sizeof(percentileSPDtrk0815_low_1) / sizeof(Double_t);

    // class 10 --> kLowMult
    Double_t percentileSPDtrk0815_low_2[] = {10, 30, 40, 60};
    Double_t percentileSPDtrk0815_high_2[] = {50, 60, 70, 100};
    Double_t percentileV0M_low_2[] = {70, 60, 60, 35};
    Double_t percentileV0M_high_2[] = {100, 90, 70, 60};
    Long_t n2 = sizeof(percentileSPDtrk0815_low_2) / sizeof(Double_t);

    // class 9 --> kHighZN
    Double_t percentileSPDtrk0815_low_3[] = {70, 40, 20, 0};
    Double_t percentileSPDtrk0815_high_3[] = {100, 70, 40, 30};
    Double_t percentileV0M_low_3[] = {10, 30, 40, 40};
    Double_t percentileV0M_high_3[] = {50, 60, 60, 80};
    Long_t n3 = sizeof(percentileSPDtrk0815_low_3) / sizeof(Double_t);

    // class 4 --> kLowZN
    Double_t percentileSPDtrk0815_low_4[] = {50, 35, 20, 10, 0};
    Double_t percentileSPDtrk0815_high_4[] = {100, 50, 35, 20, 10};
    Double_t percentileV0M_low_4[] = {0, 20, 20, 25, 30};
    Double_t percentileV0M_high_4[] = {30, 30, 40, 40, 40};
    Long_t n4 = sizeof(percentileSPDtrk0815_low_4) / sizeof(Double_t);

    Double_t percentileSPDtrk0815_low_mb[] = {0};
    Double_t percentileSPDtrk0815_high_mb[] = {100};
    Double_t percentileV0M_low_mb[] = {0};
    Double_t percentileV0M_high_mb[] = {100};
    Long_t nmb = sizeof(percentileSPDtrk0815_low_mb) / sizeof(Double_t);

    if (DoMB){
        nbins = nmb;
        for (Int_t i = 0; i < nmb; i++)
        {
            percentileSPDtrk0815_low.push_back(percentileSPDtrk0815_low_mb[i]);
            percentileSPDtrk0815_high.push_back(percentileSPDtrk0815_high_mb[i]);
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
                percentileSPDtrk0815_low.push_back(percentileSPDtrk0815_low_0[i]);
                percentileSPDtrk0815_high.push_back(percentileSPDtrk0815_high_0[i]);
                percentileV0M_low.push_back(percentileV0M_low_0[i]);
                percentileV0M_high.push_back(percentileV0M_high_0[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kStandalone class..." << endl;
        }

        if (lClassCode == kHighMult)
        {
            nbins = n1;
            for (Int_t i = 0; i < n1; i++)
            {
                percentileSPDtrk0815_low.push_back(percentileSPDtrk0815_low_1[i]);
                percentileSPDtrk0815_high.push_back(percentileSPDtrk0815_high_1[i]);
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
                percentileSPDtrk0815_low.push_back(percentileSPDtrk0815_low_2[i]);
                percentileSPDtrk0815_high.push_back(percentileSPDtrk0815_high_2[i]);
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
                percentileSPDtrk0815_low.push_back(percentileSPDtrk0815_low_3[i]);
                percentileSPDtrk0815_high.push_back(percentileSPDtrk0815_high_3[i]);
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
                percentileSPDtrk0815_low.push_back(percentileSPDtrk0815_low_4[i]);
                percentileSPDtrk0815_high.push_back(percentileSPDtrk0815_high_4[i]);
                percentileV0M_low.push_back(percentileV0M_low_4[i]);
                percentileV0M_high.push_back(percentileV0M_high_4[i]);
            }
            cout << "------------------------------------------" << endl;
            cout << " Initializing kLowZN class..." << endl;
        }
    }
}

Double_t GetPercentilefromValue(TFile *lfilename, Int_t lRun, Int_t lValue)
{ // set in the run func

    // Getting percentile
    TString name = Form("hcum%s_%i", "SPDtrk0815", lRun);

    if (lRun == 274594 || lRun == 288897){
        return 0.;
    }

    TH1D *hcum = (TH1D *)lfilename->Get(name);
    if (!hcum) cout << "Run " << lRun << " not available!" << endl;

    Double_t percentile = 100 * (hcum->GetBinContent(hcum->FindBin(lValue)));

    return percentile;
}