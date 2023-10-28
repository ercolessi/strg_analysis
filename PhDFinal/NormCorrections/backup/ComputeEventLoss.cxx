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

void ComputeEventLoss(int lClassName = kStandalone, Bool_t DoMB = kFALSE)
{

    // Percentile
    vector<Double_t> percentileSPDtrk0815_low;
    vector<Double_t> percentileSPDtrk0815_high;
    vector<Double_t> percentileV0M_low;
    vector<Double_t> percentileV0M_high;
    int nbins;
    // initiatlize
    init(lClassName, percentileSPDtrk0815_low, percentileSPDtrk0815_high, percentileV0M_low, percentileV0M_high, nbins);

    cout << " Class code enum: " << lClassName << endl;
    cout << " Found " << nbins << " selections:" << endl;
    for (int i = 0; i < nbins; i++)
    {
        cout << " SPDtrk0815 [" << percentileSPDtrk0815_low[i] << "-" << percentileSPDtrk0815_high[i] << "] + V0M [" << percentileV0M_low[i] << "-" << percentileV0M_high[i] << "]" << endl;
    }
    cout << "------------------------------------------" << endl;

    //MC
    TFile* fileMC[3];
    fileMC[0] = new TFile("LHC15g3a3_15f_normcorr.root", "READ");
    //fileMC[1] = new TFile("LHC22k1_17j_normcorr.root", "READ");
    //fileMC[2] = new TFile("LHC22k1_18i_normcorr.root", "READ");
    //
    TList* clistMC[3];
    clistMC[0] = (TList*)fileMC[0]->Get("PWGLF_StrVsMult_MC/cList");
    //clistMC[1] = (TList*)fileMC[1]->Get("PWGLF_StrVsMult_MC/cList");
    //clistMC[2] = (TList *)fileMC[2]->Get("PWGLF_StrVsMult_MC/cList");
    //DATA
    TFile* fileDATA[3];
    fileDATA[0] = new TFile("/home/fercoles/strg_analysis/PhDWork/Data/LHC15f_pass2.root","READ");
    fileDATA[1] = new TFile("/home/fercoles/strg_analysis/PhDWork/Data/LHC17j_pass2.root","READ");
    fileDATA[2] = new TFile("/home/fercoles/strg_analysis/PhDWork/Data/LHC18i_pass2.root","READ");

    TTree* lTreeEvent;
    const int percbinnumb = nbins;
    Float_t fCentrality_V0M = 0.;
    Float_t fCentrality_SPDtrk0815 = 0.;
    Int_t fRun = 0;
    Int_t fSPDtracklets0815 = 0;
    Long_t NeventsPerc[percbinnumb][3];
    for (int iset = 0; iset < 3; iset ++){
        for (int iperc = 0; iperc < percbinnumb; iperc ++){
            NeventsPerc[iperc][iset] = 0;
        }
    }

    TFile *Read0815 = TFile::Open("/home/fercoles/strg_analysis/PhDWork/Selections/PercentileCalibration.root");

    for (int iset = 0; iset < 3; iset ++){ //loop over 3 datasets

        cout<<"--------------- Open Data File " << iset << " --------------------\n"<<endl;
        lTreeEvent = (TTree*)fileDATA[iset]->Get("PWGLF_StrVsMult/fTreeEvent");
        //
        lTreeEvent->SetBranchAddress("fRun", &fRun);
        lTreeEvent->SetBranchAddress("fCentrality_V0M", &fCentrality_V0M);
        lTreeEvent->SetBranchAddress("fSPDtracklets0815", &fSPDtracklets0815);

        Long_t lNEvents = 0;
        for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries()/100; iEv++) {
            lTreeEvent->GetEntry(iEv);
            if( iEv % ( lTreeEvent->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;

            fCentrality_SPDtrk0815 = GetPercentilefromValue(Read0815, fRun, fSPDtracklets0815);

            for (int iperc = 0; iperc < percbinnumb; iperc ++){

                //set the limits
                double multmin, multmax, eemin, eemax;
                //
                multmin = percentileSPDtrk0815_low[iperc];
                multmax = percentileSPDtrk0815_high[iperc];
                eemin = percentileV0M_low[iperc];
                eemax = percentileV0M_high[iperc];

                if(fCentrality_V0M >= eemin && fCentrality_V0M < eemax && fCentrality_SPDtrk0815 >= multmin && fCentrality_SPDtrk0815 < multmax){
                    NeventsPerc[iperc][iset]++;
                }
            }
        }

        lTreeEvent = 0x0;
        fCentrality_V0M = 0.;
        fCentrality_SPDtrk0815 = 0.;

        cout << "\n---- Dataset " << iset << " ------------------------------------------ " << endl;
        for (int iperc = 0; iperc < percbinnumb; iperc ++)
        {
            double multmin, multmax, eemin, eemax;
            multmin = percentileSPDtrk0815_low[iperc];
            multmax = percentileSPDtrk0815_high[iperc];
            eemin = percentileV0M_low[iperc];
            eemax = percentileV0M_high[iperc];

            cout << "---- SPDtrk0815 [" << multmin << "-" << multmax << "] && V0M [" << eemin << "-" << eemax << "] : ....... N events = " << NeventsPerc[iperc][iset] << endl;
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
    Bool_t fEvSel_NotPileupSPDInMultBins = kFALSE;
    Bool_t fEvSel_AcceptedVertexPosition = kFALSE;
    Bool_t fEvSel_Pileuptrue = kFALSE;

    //Histograms
    TH1F* hevtloss[3];
    //TH1F* hevtlossTOT = new TH1F(Form("hevtloss"), Form(";%s percentile","temp"), percbinnumb, percentile);

    Int_t LHC15fRunList[] = {226500, 226495, 226483, 226472, 226468, 226466, 226452, 226445, 226444, 226225, 226220, 226062, 225719, 225717, 225716, 225710, 225709, 225708, 225707, 225705, 225587, 225586, 225579, 225578, 225576, 225322, 225315, 225314, 225313, 225310, 225309, 225307, 225305, 225106, 225052, 225051, 225050, 225043, 225041, 225037, 225035, 225031, 225026};
    Int_t nruns15f = sizeof(LHC15fRunList)/sizeof(Int_t);
    Int_t LHC17jRunList[] = {274593, 274595, 274596, 274601, 274653, 274657, 274667, 274669, 274671};
    Int_t nruns17j = sizeof(LHC17jRunList)/sizeof(Int_t);
    Int_t LHC18iRunList[] = {288861, 288862, 288863, 288864, 288868, 288902, 288903, 288908, 288909};
    Int_t nruns18i = sizeof(LHC18iRunList)/sizeof(Int_t);

    for (int iset = 0; iset < 3; iset ++){ //loop over 3 datasets

        cout<<"--------------- Open MC File " << iset << " --------------------\n"<<endl;

        //hevtloss[iset] = new TH1F(Form("hevtloss_MC%i",iset), Form(";%s percentile","temp"), percbinnumb, percentile);
        //hevtloss[iset]->SetTitle(Form("Event loss correction (MC %i) in %s percentile bins, %s fixed in %03.0f-%03.0f",iset, "temp","temp",0, 100));

        lTreeEventMC = (TTree *)fileMC[0]->Get("PWGLF_StrVsMult_MC/fTreeEvent");

        //Get from Tree
        lTreeEventMC->SetBranchAddress("fRun",&fRun);
        lTreeEventMC->SetBranchAddress("fEvSel_zVtxZMC", &fEvSel_zVtxZMC);
        lTreeEventMC->SetBranchAddress("fEvSel_AllSelections", &fEvSel_AllSelections);
        lTreeEventMC->SetBranchAddress("fEvSel_INELgtZEROtrue", &fEvSel_INELgtZEROtrue);
        lTreeEventMC->SetBranchAddress("fEvSel_NotPileupSPDInMultBins", &fEvSel_NotPileupSPDInMultBins);
        lTreeEventMC->SetBranchAddress("fCentrality_V0M", &fEnergy);
        //if (iset>0){
            lTreeEventMC->SetBranchAddress("fEvSel_AcceptedVertexPosition", &fEvSel_AcceptedVertexPosition);
            lTreeEventMC->SetBranchAddress("fEvSel_Pileuptrue", &fEvSel_Pileuptrue);
        //}
        lTreeEventMC->SetBranchAddress("fSPDtracklets0815", &fSPDtracklets0815);

        Long_t lNEvtSelected[percbinnumb];
        Long_t lNEvttrueSelected[percbinnumb];

        for (int k=0; k < percbinnumb; k++) { // loop on percentile classes

            //initialize counters
            lNEvtSelected[k]=0;
            lNEvttrueSelected[k]=0;

            //set the limits
            double multmin, multmax, eemin, eemax;
            multmin = percentileSPDtrk0815_low[k];
            multmax = percentileSPDtrk0815_high[k];
            eemin = percentileV0M_low[k];
            eemax = percentileV0M_high[k];

            cout<<" \nWill now loop over events, please wait...\n"<<endl;
            for(Long_t iEv = 0; iEv<lTreeEventMC->GetEntries(); iEv++) {

                lTreeEventMC->GetEntry(iEv);
                if( iEv % ( lTreeEventMC->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEventMC->GetEntries()<<endl;

                fMultiplicity = GetPercentilefromValue(Read0815, fRun, fSPDtracklets0815);

                bool isrunok = false;
                if (iset==1) {
                    for (int irun = 0; irun < nruns15f; irun ++){
                        if (fRun == LHC15fRunList[irun]) isrunok = true;
                    }
                } else if (iset==0) {
                    for (int irun = 0; irun < nruns17j; irun ++){
                        if (fRun == LHC17jRunList[irun]) isrunok = true;
                    }
                } else if (iset==2) {
                    for (int irun = 0; irun < nruns18i; irun ++){
                        if (fRun == LHC18iRunList[irun]) isrunok = true;
                    }
                }
                if (!isrunok) continue;

                //if (iset>0){
                    if( fEvSel_AllSelections &&
                        fEvSel_INELgtZEROtrue &&
                        !fEvSel_Pileuptrue &&
                        fEnergy>=eemin && fEnergy<eemax &&
                        fMultiplicity>=multmin && fMultiplicity<multmax)
                        lNEvtSelected[k]++;

                    if (fEvSel_INELgtZEROtrue &&
                        fEvSel_AcceptedVertexPosition &&
                        !fEvSel_Pileuptrue &&
                        fEnergy >= eemin && fEnergy < eemax &&
                        fMultiplicity >= multmin && fMultiplicity < multmax)
                        lNEvttrueSelected[k]++;
                //}
                /*else {
                    if (fEvSel_AllSelections &&
                        fEvSel_INELgtZEROtrue &&
                        fEnergy >= eemin && fEnergy < eemax &&
                        fMultiplicity >= multmin && fMultiplicity < multmax)
                        lNEvtSelected[k]++;

                    if( fEvSel_INELgtZEROtrue &&
                        fEvSel_zVtxZMC &&
                        fEnergy>=eemin && fEnergy<eemax &&
                        fMultiplicity>=multmin && fMultiplicity<multmax
                    ) lNEvttrueSelected[k]++;
                }*/
        }

        cout<<" --------------------------------------------------------"<<endl;
       // cout<<" This selection....................................: "<< "temp" << " [" << percentile[k] << " - " <<  percentile[k+1] << "], " << "temp" << " ["  << 0 << " - " << 100 << "]" << endl;
        cout<<"\n Number of events reco INEL>0, this selection......: "<<lNEvtSelected[k] <<endl;
        cout<<"\n Number of events true INEL>0, this selection......: "<<lNEvttrueSelected[k] <<endl;
        cout<<"\n Event Loss Correction.............................: " << (Float_t)lNEvtSelected[k]/lNEvttrueSelected[k] << endl;
        cout<<" --------------------------------------------------------"<<endl;

        } // end loop on percentile classes
    }

        /*for (int i=1; i <= hevtloss[iset]->GetNbinsX(); i++){

            hevtloss[iset]->SetBinContent(i,(Float_t)lNEvtSelected[i-1]/lNEvttrueSelected[i-1]);
            hevtloss[iset]->SetBinError(i,ErrorBinomial(
                (Float_t)lNEvtSelected[i-1],(Float_t)TMath::Sqrt(lNEvtSelected[i-1]),
                (Float_t)lNEvttrueSelected[i-1],(Float_t)TMath::Sqrt(lNEvttrueSelected[i-1]))
                );
        }

        //---------------------------------------------------------------------------------------------------

        lTreeEventMC = 0x0;
        fMultiplicity = 0.;
        fEnergy = 0.;
        fEvSel_AllSelections = kFALSE;
        fEvSel_INELgtZEROtrue = kFALSE;
        fEvSel_zVtxZMC = kFALSE;
        fEvSel_NotPileupSPDInMultBins = kFALSE;
        fEvSel_AcceptedVertexPosition = kFALSE;
        fEvSel_Pileuptrue = kFALSE;
        fRun = 0;
    }

    for (int bin = 1; bin <= hevtlossTOT->GetNbinsX(); bin ++){
        hevtlossTOT->SetBinContent(bin,DoWeightedMeanXBin(0,
                hevtloss[0]->GetBinContent(bin), hevtloss[1]->GetBinContent(bin), hevtloss[2]->GetBinContent(bin),
                hevtloss[0]->GetBinError(bin), hevtloss[1]->GetBinError(bin), hevtloss[2]->GetBinError(bin),
                NeventsPerc[bin-1][0], NeventsPerc[bin-1][1], NeventsPerc[bin-1][2]));
        hevtlossTOT->SetBinError(bin,DoWeightedMeanXBin(1,
                hevtloss[0]->GetBinContent(bin), hevtloss[1]->GetBinContent(bin), hevtloss[2]->GetBinContent(bin),
                hevtloss[0]->GetBinError(bin), hevtloss[1]->GetBinError(bin), hevtloss[2]->GetBinError(bin),
                NeventsPerc[bin-1][0], NeventsPerc[bin-1][1], NeventsPerc[bin-1][2]));
    }

    ////////////////////////////////////////////////////////
    // -------------- Write on file ------------------------
    ////////////////////////////////////////////////////////

    TFile* Write = 0x0;
    if (DoMB) {
        Write = new TFile (Form("EventLoss-%s-13TeV_INELgt0.root", lCascType.Data())
                                ,"RECREATE");
    } else{
        Write = new TFile (Form("EventLoss-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", lCascType.Data(), "temp", "temp", 0, 100)
                                ,"RECREATE");
    }

    TDirectoryFile* epsevt;
    epsevt = new TDirectoryFile("EventLoss","EventLoss");
    epsevt->cd();
    for (int iset = 0; iset < 3; iset ++){ //loop over 3 datasets
        hevtloss[iset]->Write();
    }
    hevtlossTOT->Write();
    Write->cd();*/
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
    Double_t percentileV0M_low_0[] = {0};//, 1, 5, 10, 15, 20, 30, 40, 50, 70};
    Double_t percentileV0M_high_0[] = {1};//, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    Double_t percentileSPDtrk0815_low_0[] = {0};//, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Double_t percentileSPDtrk0815_high_0[] = {100};//, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Long_t n0 = sizeof(percentileV0M_low_0) / sizeof(Double_t);

    // class 8 --> kHighMult
    Double_t percentileSPDtrk0815_low_1[] = {0, 0, 10, 20, 50};
    Double_t percentileSPDtrk0815_high_1[] = {20, 30, 40, 70, 70};
    Double_t percentileV0M_low_1[] = {40, 40, 30, 20, 10};
    Double_t percentileV0M_high_1[] = {100, 50, 40, 30, 20};
    Long_t n1 = sizeof(percentileSPDtrk0815_low_1) / sizeof(Double_t);

    // class 10 --> kLowMult
    Double_t percentileSPDtrk0815_low_2[] = {50, 40, 30, 10};
    Double_t percentileSPDtrk0815_high_2[] = {100, 70, 60, 50};
    Double_t percentileV0M_low_2[] = {35, 60, 60, 70};
    Double_t percentileV0M_high_2[] = {60, 70, 90, 100};
    Long_t n2 = sizeof(percentileSPDtrk0815_low_2) / sizeof(Double_t);

    // class 9 --> kHighZN
    Double_t percentileSPDtrk0815_low_3[] = {0, 20, 40, 70};
    Double_t percentileSPDtrk0815_high_3[] = {30, 40, 70, 100};
    Double_t percentileV0M_low_3[] = {40, 40, 30, 10};
    Double_t percentileV0M_high_3[] = {80, 60, 60, 50};
    Long_t n3 = sizeof(percentileSPDtrk0815_low_3) / sizeof(Double_t);

    // class 4 --> kLowZN
    Double_t percentileSPDtrk0815_low_4[] = {0, 10, 20, 35, 50};
    Double_t percentileSPDtrk0815_high_4[] = {10, 20, 35, 50, 100};
    Double_t percentileV0M_low_4[] = {30, 25, 20, 20, 0};
    Double_t percentileV0M_high_4[] = {40, 40, 40, 30, 30};
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

    if (lRun == 274594 || lRun == 288897){
        return 0.;
    }

    TH1D *hcum = (TH1D *)lfilename->Get(name);
    if (!hcum) cout << "Run " << lRun << " not available!" << endl;

    Double_t percentile = 100 * (hcum->GetBinContent(hcum->FindBin(lValue)));

    return percentile;
}