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
Double_t GetValueFromPercentile(TFile *lfilename, Double_t perc);

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

void ComputeSignalLoss(
    TString lCascType = "Xi",
    int lClassName = 3,
    Bool_t DoMB = kFALSE)
{

    // Percentile
    vector<Double_t> percentileSPDClusters_low;
    vector<Double_t> percentileSPDClusters_high;
    vector<Double_t> percentileV0M_low;
    vector<Double_t> percentileV0M_high;
    int nbins;
    // initiatlize
    init(lClassName, percentileSPDClusters_low, percentileSPDClusters_high, percentileV0M_low, percentileV0M_high, nbins, DoMB);
    const int percbinnumb = nbins;

    cout << " Class code enum: " << lClassName << endl;
    cout << " Found " << nbins << " selections:" << endl;
    for (int i = 0; i < nbins; i++)
    {
        cout << " SPDClusters [" << percentileSPDClusters_low[i] << "-" << percentileSPDClusters_high[i] << "] + V0M [" << percentileV0M_low[i] << "-" << percentileV0M_high[i] << "]" << endl;
    }
    cout << "------------------------------------------" << endl;

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
    if (lCascType.Contains("K0Short")) {
        particle = "K0Short";
        antiparticle = "";
    }

    //MC
    TFile* fileMC = new TFile("FullMC_15f17j18i_normcorr_2.root", "READ");
    //
    TList* clistMC = (TList*)fileMC->Get("PWGLF_StrVsMult_MC/cList");
    //DATA
    TFile* fileDATA[3];
    fileDATA[0] = new TFile("/home/fercoles/strg_analysis/PhDWork/Data/LHC15f_pass2.root","READ");
    fileDATA[1] = new TFile("/home/fercoles/strg_analysis/PhDWork/Data/LHC17j_pass2.root","READ");
    fileDATA[2] = new TFile("/home/fercoles/strg_analysis/PhDWork/Data/LHC18i_pass2.root","READ");

    TH1D* fHistptSel[percbinnumb];
    TH1D* fHistpttrueSel[percbinnumb];
    TH3D* h3DGenSel;
    TH3D* h3DGenSelAntiPart;
    TH3D* h3DGentrueSel;
    TH3D* h3DGentrueSelAntiPart;
    TH1D* fHistGentrueSel[percbinnumb];
    TH1D* fHistGenSel[percbinnumb];
    TH1D *hsgnloss[percbinnumb];

    Double_t temppt = 0, temppt2 = 0;
    Int_t minmultbin = 0, maxmultbin = 0, mineebin = 0, maxeebin = 0, minmultbintrue = 0, maxmultbintrue = 0, mineebintrue = 0, maxeebintrue = 0;

    // INEL > 0
    h3DGenSel = (TH3D*)clistMC->FindObject(Form("fHistGeneratedPtVsSPDclVs%s%s", lEnergyEstimator.Data(), particle.Data()));
    h3DGenSel->Sumw2();
    if (!particle.Contains("K0Short")){
        h3DGenSelAntiPart = (TH3D *)clistMC->FindObject(Form("fHistGeneratedPtVsSPDclVs%s%s", lEnergyEstimator.Data(), antiparticle.Data()));
        h3DGenSelAntiPart->Sumw2();
    }
    // true INEL > 0
    h3DGentrueSel = (TH3D*)clistMC->FindObject(Form("fHistPtVsSPDclVs%s_Gen%s", lEnergyEstimator.Data(), particle.Data()));
    h3DGentrueSel->Sumw2();
    if (!particle.Contains("K0Short")){
        h3DGentrueSelAntiPart = (TH3D *)clistMC->FindObject(Form("fHistPtVsSPDclVs%s_Gen%s", lEnergyEstimator.Data(), antiparticle.Data()));
        h3DGentrueSelAntiPart->Sumw2();
    }

    TH1D *h1 = (TH1D *)h3DGentrueSel->ProjectionY("h1", 0, -1, 1, 10);
    TH1D *h2 = (TH1D *)h3DGentrueSel->ProjectionY("h2", 0, -1, 1, 100);

    new TCanvas;
    h2->Draw("HIST");
    h1->SetLineColor(kRed);
    h1->Draw("HIST same");

    // Add particle + antiparticle
    if (!particle.Contains("K0Short")) {
        h3DGenSel->Add(h3DGenSelAntiPart);
        h3DGentrueSel->Add(h3DGentrueSelAntiPart);
    }

    for (int k = 0; k < percbinnumb; k++) { // loop on percentile classes

        // set the limits
        double multmin, multmax, eemin, eemax;
        multmin = percentileSPDClusters_low[k];
        multmax = percentileSPDClusters_high[k];
        eemin = percentileV0M_low[k];
        eemax = percentileV0M_high[k];

        // get bin numbers for limits
        minmultbin = h3DGenSel->GetYaxis()->FindBin(multmin);
        maxmultbin = h3DGenSel->GetYaxis()->FindBin(multmax);
        mineebin = h3DGenSel->GetZaxis()->FindBin(eemin);
        maxeebin = h3DGenSel->GetZaxis()->FindBin(eemax);
        //
        minmultbintrue = h3DGentrueSel->GetYaxis()->FindBin(multmin);
        maxmultbintrue = h3DGentrueSel->GetYaxis()->FindBin(multmax);
        mineebintrue = h3DGentrueSel->GetZaxis()->FindBin(eemin);
        maxeebintrue = h3DGentrueSel->GetZaxis()->FindBin(eemax);

        fHistGenSel[k] = (TH1D *)h3DGenSel->ProjectionX(Form("fHistGenSel%i%s", k, lCascType.Data()),
                                                        minmultbin, maxmultbin, mineebin, maxeebin);
        fHistGentrueSel[k] = (TH1D *)h3DGentrueSel->ProjectionX(Form("fHistGenSelTrue%i%s", k, lCascType.Data()),
                                                                minmultbintrue, maxmultbintrue, mineebintrue, maxeebintrue);

        fHistptSel[k] = (TH1D *)fHistGenSel[k]->Clone(Form("INELgt0RECO_SPDClusters_%.0f-%.0f_V0M_%.0f-%.0f", percentileSPDClusters_low[k], percentileSPDClusters_high[k], percentileV0M_low[k], percentileV0M_high[k]));
        fHistpttrueSel[k] = (TH1D *)fHistGentrueSel[k]->Clone(Form("INELgt0TRUE_%.0f-%.0f_V0M_%.0f-%.0f", percentileSPDClusters_low[k], percentileSPDClusters_high[k], percentileV0M_low[k], percentileV0M_high[k]));

        //fHistptSel[k]->Sumw2();
        //fHistpttrueSel[k]->Sumw2();

        hsgnloss[k] = (TH1D *)fHistptSel[k]->Clone(Form("hsgnloss_SPDClusters_%.0f-%.0f_V0M_%.0f-%.0f", percentileSPDClusters_low[k], percentileSPDClusters_high[k], percentileV0M_low[k], percentileV0M_high[k]));

        ErrorInRatioHisto(hsgnloss[k], fHistGentrueSel[k]);
    }

    ////////////////////////////////////////////////////////
    // -------------- Write on file ------------------------
    ////////////////////////////////////////////////////////

    TFile* Write = 0x0;
    if (DoMB) {
        Write = new TFile (Form("SignalLoss-%s-13TeV_INELgt0.root", lCascType.Data())
                             ,"RECREATE");
    } else{
        Write = new TFile(Form("SignalLoss-%s-13TeV_class%i.root", lCascType.Data(), lClassName)
                            , "RECREATE");
    }

    h3DGentrueSel->Write();
    h3DGenSel->Write();
    for (int k = 0; k<percbinnumb; k++){
        hsgnloss[k]->Write();
        fHistGentrueSel[k]->Write();
        fHistptSel[k]->Write();
    }
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
    Double_t percentileV0M_low_10[] = {0, 5, 15, 30, 50};
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

    if (DoMB)
    {
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
    }
    else
    {

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
            cout << " Initializing kStandalone for Omega class..." << endl;
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

    if (lRun == 274594 || lRun == 288897)
    {
        return 0.;
    }

    TH1D *hcum = (TH1D *)lfilename->Get(name);
    if (!hcum)
        cout << "Run " << lRun << " not available!" << endl;

    Double_t percentile = 100 * (hcum->GetBinContent(hcum->FindBin(lValue)));

    return percentile;
}

//==========================================================================================
Double_t GetValueFromPercentile(TFile *lfilename, Double_t perc)
{
    // Return the bin of the cumulative once given a percentile

    // Getting percentile
    TString name = Form("hcum%s_%s", "SPDClusters", "all");
    TH1D *h = (TH1D *)lfilename->Get(name);

    Double_t diff = 10000.;
    Double_t val = 0;

    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
        if (diff > TMath::Abs(h->GetBinContent(i) - perc / 100))
        {
            val = h->GetBinCenter(i);
            diff = TMath::Abs(h->GetBinContent(i) - perc / 100);
        }
    }

    return val;
}