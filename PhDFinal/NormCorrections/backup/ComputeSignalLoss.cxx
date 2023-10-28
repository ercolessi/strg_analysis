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
          int &nbins);
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

void ComputeSignalLoss(
    TString lCascType = "Lambda",
    int lClassName = kHighMult,
    Bool_t DoMB = kFALSE)
{

    // Percentile
    vector<Double_t> percentileSPDtrk0815_low;
    vector<Double_t> percentileSPDtrk0815_high;
    vector<Double_t> percentileV0M_low;
    vector<Double_t> percentileV0M_high;
    int nbins;
    // initiatlize
    init(lClassName, percentileSPDtrk0815_low, percentileSPDtrk0815_high, percentileV0M_low, percentileV0M_high, nbins);
    const int percbinnumb = nbins;

    cout << " Class code enum: " << lClassName << endl;
    cout << " Found " << nbins << " selections:" << endl;
    for (int i = 0; i < nbins; i++)
    {
        cout << " SPDtrk0815 [" << percentileSPDtrk0815_low[i] << "-" << percentileSPDtrk0815_high[i] << "] + V0M [" << percentileV0M_low[i] << "-" << percentileV0M_high[i] << "]" << endl;
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

    //pT bins
    Double_t* ptbinlimits;
    Long_t ptbinnumb;
    Double_t ptbinlimitsXi[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
    Long_t ptbinnumbXi = sizeof(ptbinlimitsXi)/sizeof(Double_t) - 1;
    Double_t ptbinlimitsOmega[] = {0.90, 1.60, 2.20, 2.60, 3.00, 3.80, 5.50};
    Long_t ptbinnumbOmega = sizeof(ptbinlimitsOmega)/sizeof(Double_t) - 1;
    Double_t ptbinlimitsLambda[] = {0.4,  0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.5, 2.9, 3.4, 4, 5, 6.5, 8};
    Long_t ptbinnumbLambda = sizeof(ptbinlimitsLambda)/sizeof(Double_t) - 1;
    Double_t ptbinlimitsK0S[] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3.0,
                                 3.3, 3.6, 3.9, 4.2, 4.6, 5, 5.4, 5.9, 6.5, 7, 7.5, 8, 8.5, 9.2, 10, 11, 12, 13.5, 15};
    Long_t ptbinnumbK0S = sizeof(ptbinlimitsK0S) / sizeof(Double_t) - 1;

    //MC
    TFile* fileMC = new TFile("FullMC_15f17j18i_normcorr.root", "READ");
    //
    TList* clistMC = (TList*)fileMC->Get("PWGLF_StrVsMult_MC/cList");
    //DATA
    TFile* fileDATA[3];
    fileDATA[0] = new TFile("/home/fercoles/strg_analysis/PhDWork/Data/LHC15f_pass2.root","READ");
    fileDATA[1] = new TFile("/home/fercoles/strg_analysis/PhDWork/Data/LHC17j_pass2.root","READ");
    fileDATA[2] = new TFile("/home/fercoles/strg_analysis/PhDWork/Data/LHC18i_pass2.root","READ");

    TTree* lTreeEvent;
    Float_t fCentrality_V0M = 0.;
    Float_t fCentrality_SPDClusters = 0.;
    Long_t NeventsPerc[percbinnumb][3];
    Long_t NMCeventsPerc[percbinnumb][3];
    Float_t wdata[percbinnumb][3], wMC[percbinnumb][3],wtot[percbinnumb][3];
    for (int iset = 0; iset < 3; iset ++){
        for (int iperc = 0; iperc < percbinnumb; iperc ++){
            NeventsPerc[iperc][iset] = 0;
            NMCeventsPerc[iperc][iset] = 0;
            wdata[iperc][iset] = 0;
            wMC[iperc][iset] = 0;
            wtot[iperc][iset] = 0;
        }
    }

    for (int iset = 0; iset < 3; iset ++){ //loop over 3 datasets

        cout<<"--------------- Open Data File " << iset << " --------------------\n"<<endl;
        lTreeEvent = (TTree*)fileDATA[iset]->Get("PWGLF_StrVsMult/fTreeEvent");
        //
        lTreeEvent->SetBranchAddress("fCentrality_V0M", &fCentrality_V0M);
        lTreeEvent->SetBranchAddress("fCentrality_SPDClusters", &fCentrality_SPDClusters);

        Long_t lNEvents = 0;
        for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries(); iEv++) {
            lTreeEvent->GetEntry(iEv);
            if( iEv % ( lTreeEvent->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;

            for (int iperc = 0; iperc < percbinnumb; iperc ++){

                //set the limits
                double multmin, multmax, eemin, eemax;
                //
                if (lWhichVarEstimator.Contains("SPD")){
                    multmin = percentile[iperc];
                    multmax = percentile[iperc+1];
                    eemin = lFixedLo;
                    eemax = lFixedHi;
                }
                //
                if (lWhichVarEstimator.Contains("V0M")){
                    eemin = percentile[iperc];
                    eemax = percentile[iperc+1];
                    multmin = lFixedLo;
                    multmax = lFixedHi;
                }

                if(fCentrality_V0M >= eemin && fCentrality_V0M < eemax && fCentrality_SPDClusters >= multmin && fCentrality_SPDClusters < multmax){
                    NeventsPerc[iperc][iset]++;
                }
            }
        }

        lTreeEvent = 0x0;
        fCentrality_V0M = 0.;
        fCentrality_SPDClusters = 0.;

        cout << "---- Dataset " << iset << " ------------------------------------------ \n" << endl;
        for (int iperc = 0; iperc < percbinnumb; iperc ++){

            //set the limits
            double multmin, multmax, eemin, eemax;
            //
            if (lWhichVarEstimator.Contains("SPD")){
                multmin = percentile[iperc];
                multmax = percentile[iperc+1];
                eemin = lFixedLo;
                eemax = lFixedHi;
            }
            //
            if (lWhichVarEstimator.Contains("V0M")){
                eemin = percentile[iperc];
                eemax = percentile[iperc+1];
                multmin = lFixedLo;
                multmax = lFixedHi;
            }

            cout << "---- SPD [" << multmin << "-" << multmax <<"] && V0M [" << eemin << "-" << eemax <<"] : ....... N events = " <<  NeventsPerc[iperc][iset] << endl;
        }
        cout << "\n\n";
    }

    //Variables
    TTree* lTreeEventMC;
    Float_t fMultiplicity = 0.;
    Float_t fEnergy = 0.;
    Bool_t fEvSel_AllSelections = kFALSE;
    Bool_t fEvSel_INELgtZEROtrue = kFALSE;
    Bool_t fEvSel_zVtxZMC = kFALSE;
    Int_t fRun = 0.;

    for (int iset = 0; iset < 3; iset ++){ //loop over 3 datasets

        cout<<"--------------- Open MC File " << iset << " --------------------\n"<<endl;

        if (iset == 0) {
            lTreeEventMC = (TTree*)fileMC[iset]->Get("PWGLF_StrVsMult_MC/fTreeEvent");
        } else {
            lTreeEventMC = (TTree *)fileMC[iset]->Get("PWGLF_StrVsMult_MC/fTreeEvent");
        }

        //Get from Tree
        lTreeEventMC->SetBranchAddress("fRun",&fRun);
        lTreeEventMC->SetBranchAddress("fEvSel_zVtxZMC", &fEvSel_zVtxZMC);
        lTreeEventMC->SetBranchAddress("fEvSel_AllSelections", &fEvSel_AllSelections);
        lTreeEventMC->SetBranchAddress("fEvSel_INELgtZEROtrue", &fEvSel_INELgtZEROtrue);
        lTreeEventMC->SetBranchAddress("fCentrality_V0M", &fCentrality_V0M);
        lTreeEventMC->SetBranchAddress("fCentrality_SPDClusters", &fCentrality_SPDClusters);

        Long_t lNEvents = 0;
        for(Long_t iEv = 0; iEv<lTreeEventMC->GetEntries(); iEv++) {
            lTreeEventMC->GetEntry(iEv);
            if( iEv % ( lTreeEventMC->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEventMC->GetEntries()<<endl;

            if(!fEvSel_AllSelections) continue;

            for (int iperc = 0; iperc < percbinnumb; iperc ++){

                //set the limits
                double multmin, multmax, eemin, eemax;
                //
                if (lWhichVarEstimator.Contains("SPD")){
                    multmin = percentile[iperc];
                    multmax = percentile[iperc+1];
                    eemin = lFixedLo;
                    eemax = lFixedHi;
                }
                //
                if (lWhichVarEstimator.Contains("V0M")){
                    eemin = percentile[iperc];
                    eemax = percentile[iperc+1];
                    multmin = lFixedLo;
                    multmax = lFixedHi;
                }

                if(fCentrality_V0M >= eemin && fCentrality_V0M < eemax && fCentrality_SPDClusters >= multmin && fCentrality_SPDClusters < multmax){
                    NMCeventsPerc[iperc][iset]++;
                }
            }// end loop on percentile classes
        }

        lTreeEventMC = 0x0;
        fMultiplicity = 0.;
        fEnergy = 0.;
        fEvSel_AllSelections = kFALSE;
        fEvSel_INELgtZEROtrue = kFALSE;
        fEvSel_zVtxZMC = kFALSE;
        fRun = 0;

        cout << "\n -------- MC Dataset " << iset << " ------------------------------------------ \n" << endl;
        for (int iperc = 0; iperc < percbinnumb; iperc ++){

            //set the limits
            double multmin, multmax, eemin, eemax;
            //
            if (lWhichVarEstimator.Contains("SPD")){
                multmin = percentile[iperc];
                multmax = percentile[iperc+1];
                eemin = lFixedLo;
                eemax = lFixedHi;
            }
            //
            if (lWhichVarEstimator.Contains("V0M")){
                eemin = percentile[iperc];
                eemax = percentile[iperc+1];
                multmin = lFixedLo;
                multmax = lFixedHi;
            }

            cout << "---- SPD [" << multmin << "-" << multmax <<"] && V0M [" << eemin << "-" << eemax <<"] : ....... N events = " <<  NeventsPerc[iperc][iset] << endl;
        }
        cout << "\n\n";
    }

    for (int k=0; k < percbinnumb; k++) {
        for (int iset = 0; iset < 3; iset ++){ //loop over 3 datasets
            wdata[k][iset] = (float)NeventsPerc[k][iset]/(NeventsPerc[k][0]+NeventsPerc[k][1]+NeventsPerc[k][2]);
            wMC[k][iset] = (float)NMCeventsPerc[k][iset];///(NMCeventsPerc[k][0]+NMCeventsPerc[k][1]+NMCeventsPerc[k][2]);
            wtot[k][iset] =  wdata[k][iset];
            ///wMC[k][iset];
            //wdata[k][iset]/wMC[k][iset];
        }
        cout << " 0 " << wdata[k][0] << " mc " << wMC[k][0] << " tot " << wtot[k][0] << endl;
        cout << " 1 " << wdata[k][1] << " mc " << wMC[k][1] << " tot " << wtot[k][1] << endl;
        cout << " 2 " << wdata[k][2] << " mc " << wMC[k][2] << " tot " << wtot[k][2] << endl;
    }


    TH1D* fHistptSel[percbinnumb][3];
    TH1D* fHistpttrueSel[percbinnumb][3];
    TH3D* h3DGenSel[3];
    TH3D* h3DGenSelAntiPart[3];
    TH3D* h3DGentrueSel[3];
    TH3D* h3DGentrueSelAntiPart[3];
    TH1D* fHistGentrueSel[percbinnumb][3];
    TH1D* fHistGenSel[percbinnumb][3];

    for (int iset = 0; iset < 3; iset ++){ //loop over 3 datasets

        cout<<"--------------- Open MC File " << iset << " --------------------\n"<<endl;

        Double_t temppt = 0, temppt2 = 0;
        Int_t minmultbin = 0, maxmultbin = 0, mineebin = 0, maxeebin = 0, minmultbintrue = 0, maxmultbintrue = 0, mineebintrue = 0, maxeebintrue = 0;

        // INEL > 0
        h3DGenSel[iset] = (TH3D*)clistMC[iset]->FindObject(Form("fHistGeneratedPtVsSPDclVs%s%s", lEnergyEstimator.Data(), particle.Data()));
        h3DGenSel[iset]->Sumw2();
        if (!particle.Contains("K0Short")){
            h3DGenSelAntiPart[iset] = (TH3D*)clistMC[iset]->FindObject(Form("fHistGeneratedPtVsSPDclVs%s%s", lEnergyEstimator.Data(),antiparticle.Data()));
            h3DGenSelAntiPart[iset]->Sumw2();
        }
        // true INEL > 0
        h3DGentrueSel[iset] = (TH3D*)clistMC[iset]->FindObject(Form("fHistPtVsSPDclVs%s_Gen%s", lEnergyEstimator.Data(), particle.Data()));
        h3DGentrueSel[iset]->Sumw2();
        if (!particle.Contains("K0Short")){
            h3DGentrueSelAntiPart[iset] = (TH3D *)clistMC[iset]->FindObject(Form("fHistPtVsSPDclVs%s_Gen%s", lEnergyEstimator.Data(), antiparticle.Data()));
            h3DGentrueSelAntiPart[iset]->Sumw2();
        }

        // Add particle + antiparticle
        if (!particle.Contains("K0Short")) {
            h3DGenSel[iset]->Add(h3DGenSelAntiPart[iset]);
            h3DGentrueSel[iset]->Add(h3DGentrueSelAntiPart[iset]);
        }

        for (int k = 0; k < percbinnumb; k++)
        { // loop on percentile classes

            // set the limits
            double multmin, multmax, eemin, eemax;
            //
            if (lWhichVarEstimator.Contains("SPD"))
            {
                multmin = percentile[k];
                multmax = percentile[k + 1];
                eemin = lFixedLo;
                eemax = lFixedHi;
            }
            //
            if (lWhichVarEstimator.Contains("V0M") || lWhichVarEstimator.Contains("ZDC"))
            {
                eemin = percentile[k];
                eemax = percentile[k + 1];
                multmin = lFixedLo;
                multmax = lFixedHi;
            }

            minmultbin = h3DGenSel[iset]->GetYaxis()->FindBin(multmin);
            maxmultbin = h3DGenSel[iset]->GetYaxis()->FindBin(multmax);
            mineebin = h3DGenSel[iset]->GetZaxis()->FindBin(eemin);
            maxeebin = h3DGenSel[iset]->GetZaxis()->FindBin(eemax);
            //
            minmultbintrue = h3DGentrueSel[iset]->GetYaxis()->FindBin(multmin);
            maxmultbintrue = h3DGentrueSel[iset]->GetYaxis()->FindBin(multmax);
            mineebintrue = h3DGentrueSel[iset]->GetZaxis()->FindBin(eemin);
            maxeebintrue = h3DGentrueSel[iset]->GetZaxis()->FindBin(eemax);

            fHistGenSel[k][iset] = (TH1D *)h3DGenSel[iset]->ProjectionX(Form("fHistGenSel%i%s", k, lCascType.Data()),
                                                                        minmultbin, maxmultbin, mineebin, maxeebin);
            fHistGentrueSel[k][iset] = (TH1D *)h3DGentrueSel[iset]->ProjectionX(Form("fHistGenSelTrue%i%s", k, lCascType.Data()),
                                                                                minmultbintrue, maxmultbintrue, mineebintrue, maxeebintrue);

            fHistptSel[k][iset] = (TH1D *)fHistGenSel[k][iset]->Clone(Form("hsgnloss_%i-%i_MC%i", (int)percentile[k], (int)percentile[k + 1], iset));
            fHistpttrueSel[k][iset] = (TH1D *)fHistGentrueSel[k][iset]->Clone(Form("DENtrue_%i_%i_MC%i", (int)percentile[k], (int)percentile[k + 1], iset));

            fHistptSel[k][iset]->Sumw2();
            fHistpttrueSel[k][iset]->Sumw2();

            wtot[k][iset] *= 1. / (fHistpttrueSel[k][iset]->Integral(0, -1));

            fHistptSel[k][iset]->Scale(wtot[k][iset]);
            fHistpttrueSel[k][iset]->Scale(wtot[k][iset]);

            ErrorInRatioHisto(fHistGenSel[k][iset], fHistGentrueSel[k][iset]);
            fHistGenSel[k][iset]->SetName(Form("hsgnloss_%i-%i_MC%i", (int)percentile[k], (int)percentile[k + 1], iset));
        }
    }

    TH1D* fHistptSelTOT[percbinnumb], * fHistpttrueSelTOT[percbinnumb];
    //
    for (int k=0; k < percbinnumb; k++) {
        fHistptSelTOT[k] = (TH1D*) fHistptSel[k][0]->Clone(Form("wsgnloss_%i-%i",(int)percentile[k],(int)percentile[k+1]));
        fHistpttrueSelTOT[k] = (TH1D*) fHistpttrueSel[k][0]->Clone(Form("DENTOT%i",k));

        fHistptSelTOT[k]->Add(fHistptSel[k][1]);
        fHistptSelTOT[k]->Add(fHistptSel[k][2]);
        fHistpttrueSelTOT[k]->Add(fHistpttrueSel[k][1]);
        fHistpttrueSelTOT[k]->Add(fHistpttrueSel[k][2]);
    }

    for (int k = 0; k<percbinnumb; k++){
        ErrorInRatioHisto(fHistptSelTOT[k],fHistpttrueSelTOT[k]);
    }


    ////////////////////////////////////////////////////////
    // -------------- Write on file ------------------------
    ////////////////////////////////////////////////////////

    TFile* Write = 0x0;
    if (DoMB) {
        Write = new TFile (Form("SignalLoss-%s-13TeV_INELgt0.root", lCascType.Data())
                                ,"RECREATE");
    } else{
        Write = new TFile (Form("SignalLoss-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", lCascType.Data(), lWhichVarEstimator.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi)
                                ,"RECREATE");
    }
    for (int iset = 0; iset < 3; iset ++){ //loop over 3 datasets
        for (int k = 0; k<percbinnumb; k++){
            fHistGenSel[k][iset]->Write();
        }
    }
    for (int k = 0; k<percbinnumb; k++){
        fHistptSelTOT[k]->Write();
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
          vector<Double_t> &percentileSPDtrk0815_low,
          vector<Double_t> &percentileSPDtrk0815_high,
          vector<Double_t> &percentileV0M_low,
          vector<Double_t> &percentileV0M_high,
          int &nbins)
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

Double_t GetPercentilefromValue(TFile *lfilename, Int_t lRun, Int_t lValue)
{ // set in the run func

    // Getting percentile
    TString name = Form("hcum%s_%i", "SPDtrk0815", lRun);

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