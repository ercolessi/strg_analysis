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

void ComputeNormCorr(
    Bool_t DoEventLoss = kTRUE,
    Bool_t DoSignalLoss = kTRUE,
    TString lCascType = "Xi",	
    Double_t lFixedLo = 0.0,
	Double_t lFixedHi = 100.0,
	TString lWhichVarEstimator = "SPDClusters",
	TString lWhichFixedEstimator = "ZDC"
    ){

    cout << "Please check the percentile binning is fine before running!!!\n" << endl;

    TFile* file = new TFile("~/Scaricati/fffff.root","READ");

	const char* outputname = Form("Norm-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", 
                                lCascType.Data(), lWhichVarEstimator.Data(), lWhichFixedEstimator.Data(), lFixedLo, lFixedHi);

    TString lEnergyEstimator = "";
    if (lWhichVarEstimator.Contains("V0M") || lWhichFixedEstimator.Contains("V0M")) lEnergyEstimator = "V0M";
    if (lWhichVarEstimator.Contains("ZDC") || lWhichFixedEstimator.Contains("ZDC")) lEnergyEstimator = "ZDC";  
    
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

	cout<<"--------------- Open Data File --------------------"<<endl;
    TList* clist      = (TList*)file->Get("PWGLF_StrVsMult_MC/cList");
    TTree* lTreeEvent = (TTree*)file->Get("PWGLF_StrVsMult_MC/fTreeEvent");

    //Variables
    Float_t fMultiplicity = 0.;
    Float_t fEnergy = 0.;
    Bool_t fEvSel_AllSelections = kFALSE;
    Bool_t fEvSel_INELgtZEROtrue = kFALSE;
    Bool_t fEvSel_zVtxZMC = kFALSE;
    Int_t fRun;
   
    //Get from Tree
    lTreeEvent->SetBranchAddress("fRun",&fRun);
    lTreeEvent->SetBranchAddress("fEvSel_zVtxZMC", &fEvSel_zVtxZMC);
    lTreeEvent->SetBranchAddress("fEvSel_AllSelections", &fEvSel_AllSelections);
    lTreeEvent->SetBranchAddress("fEvSel_INELgtZEROtrue", &fEvSel_INELgtZEROtrue);
    lTreeEvent->SetBranchAddress("fCentrality_SPDClusters", &fMultiplicity);
    if (lWhichVarEstimator.Contains("V0M") || lWhichFixedEstimator.Contains("V0M")) lTreeEvent->SetBranchAddress("fCentrality_V0M", &fEnergy);
    if (lWhichVarEstimator.Contains("ZDC") || lWhichFixedEstimator.Contains("ZDC")) lTreeEvent->SetBranchAddress("fCentrality_ZDC", &fEnergy);
    
    //Percentile
    Float_t percentile[] = {0,10};
    Long_t percbinnumb = sizeof(percentile)/sizeof(Float_t) - 1;
    Long_t lNEvtSelected[percbinnumb];
    Long_t lNEvttrueSelected[percbinnumb];    

    //Histograms
    TH1F* hevtloss = new TH1F("hevtloss", Form(";%s percentile",lWhichVarEstimator.Data()), percbinnumb, percentile);
    hevtloss->SetTitle(Form("Event loss correction in %s percentile bins, %s fixed in %03.0f-%03.0f",lWhichVarEstimator.Data(),lWhichFixedEstimator.Data(),lFixedLo, lFixedHi));

    //pT bins 
    Double_t* ptbinlimits;
    Long_t ptbinnumb;
    Double_t ptbinlimitsXi[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5, 8.0, 10.0 };
    Long_t ptbinnumbXi = sizeof(ptbinlimitsXi)/sizeof(Double_t) - 1;
    Double_t ptbinlimitsOmega[] = {0.90, 1.60, 2.20, 2.60, 3.00, 3.80, 5.50, 8.00, 10.0 }; 
    Long_t ptbinnumbOmega = sizeof(ptbinlimitsOmega)/sizeof(Double_t) - 1;
    Double_t ptbinlimitsLambda[] = {0.4,  0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.5, 2.9, 3.4, 4, 5, 6.5, 8, 10};
    Long_t ptbinnumbLambda = sizeof(ptbinlimitsLambda)/sizeof(Double_t) - 1;;
  
    
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

    if (DoEventLoss){
   
        for (int k=0; k < percbinnumb; k++) { // loop on percentile classes 
            
            //initialize counters
            lNEvtSelected[k]=0;
            lNEvttrueSelected[k]=0;

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

            cout<<" \nWill now loop over events, please wait...\n"<<endl;
            for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries(); iEv++) {

                lTreeEvent->GetEntry(iEv);
                if( iEv % ( lTreeEvent->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;

                if( fEvSel_AllSelections && 
                    fEnergy>=eemin && fEnergy<eemax &&
                    fMultiplicity>=multmin && fMultiplicity<multmax
                ) lNEvtSelected[k]++;

                if( fEvSel_INELgtZEROtrue && 
                    fEvSel_zVtxZMC && 
                    fEnergy>=eemin && fEnergy<eemax &&
                    fMultiplicity>=multmin && fMultiplicity<multmax
                ) lNEvttrueSelected[k]++;
        }

        cout<<" --------------------------------------------------------"<<endl;
        cout<<" This selection....................................: "<< lWhichVarEstimator.Data() << " [" << percentile[k] << " - " <<  percentile[k+1] << "], " << lWhichFixedEstimator.Data() << " ["  << lFixedLo << " - " << lFixedHi << "]" << endl;  
        cout<<"\n Number of events reco INEL>0, this selection......: "<<lNEvtSelected[k] <<endl;
        cout<<"\n Number of events true INEL>0, this selection......: "<<lNEvttrueSelected[k] <<endl;
        cout<<"\n Event Loss Correction.............................: " << (Float_t)lNEvtSelected[k]/lNEvttrueSelected[k] << endl;
        cout<<" --------------------------------------------------------"<<endl;
        
        } // end loop on percentile classes 

        for (int i=1; i <= hevtloss->GetNbinsX(); i++){

            hevtloss->SetBinContent(i,(Float_t)lNEvtSelected[i-1]/lNEvttrueSelected[i-1]);
            hevtloss->SetBinError(i,ErrorInRatioCorr(
                (Float_t)lNEvtSelected[i-1],(Float_t)TMath::Sqrt(lNEvtSelected[i-1]),
                (Float_t)lNEvttrueSelected[i-1],(Float_t)TMath::Sqrt(lNEvttrueSelected[i-1]))
                );
        }
    }
    //---------------------------------------------------------------------------------------------------

    TH3D* h3DGenSel;
    TH3D* h3DGenSelAntiPart;
    TH3D* h3DGentrueSel;  
    TH3D* h3DGentrueSelAntiPart;   
    TH1D* fHistGentrueSel[percbinnumb];
    TH1D* fHistGenSel[percbinnumb];
    TH1D* fHistptSel[percbinnumb];
    TH1D* fHistpttrueSel[percbinnumb];
    Double_t temppt = 0, temppt2 = 0;
    Int_t minmultbin = 0, maxmultbin = 0, mineebin = 0, maxeebin = 0, minmultbintrue = 0, maxmultbintrue = 0, mineebintrue = 0, maxeebintrue = 0;
    
    if (DoSignalLoss){

        // INEL > 0
        h3DGenSel = (TH3D*)clist->FindObject(Form("fHistGeneratedPtVsSPDclVs%s%s", lEnergyEstimator.Data(), particle.Data()));
        h3DGenSel->Sumw2();
        h3DGenSelAntiPart = (TH3D*)clist->FindObject(Form("fHistGeneratedPtVsSPDclVs%s%s", lEnergyEstimator.Data(),antiparticle.Data()));
        h3DGenSelAntiPart->Sumw2();
        // true INEL > 0
        h3DGentrueSel = (TH3D*)clist->FindObject(Form("fHistPtVsSPDclVs%s_Gen%s", lEnergyEstimator.Data(), particle.Data()));
        h3DGentrueSel->Sumw2();
        h3DGentrueSelAntiPart = (TH3D*)clist->FindObject(Form("fHistPtVsSPDclVs%s_Gen%s", lEnergyEstimator.Data(), antiparticle.Data()));
        h3DGentrueSelAntiPart->Sumw2();

        // Add particle + antiparticle
        h3DGenSel->Add(h3DGenSelAntiPart);
        h3DGentrueSel->Add(h3DGentrueSelAntiPart);

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

            minmultbin = h3DGenSel->GetYaxis()->FindBin( multmin );
            maxmultbin = h3DGenSel->GetYaxis()->FindBin( multmax-1e-6 );
            mineebin = h3DGenSel->GetZaxis()->FindBin( eemin );
            maxeebin = h3DGenSel->GetZaxis()->FindBin( eemax-1e-6 );
            //
            minmultbintrue = h3DGentrueSel->GetYaxis()->FindBin( multmin );
            maxmultbintrue = h3DGentrueSel->GetYaxis()->FindBin( multmax-1e-6 );
            mineebintrue = h3DGentrueSel->GetZaxis()->FindBin( eemin );
            maxeebintrue = h3DGentrueSel->GetZaxis()->FindBin( eemax-1e-6 );
            
            fHistGenSel[k] = (TH1D*)h3DGenSel->ProjectionX(Form("fHistGenSel%i%s",k,lCascType.Data()),
                    minmultbin, maxmultbin, mineebin, maxeebin );
            fHistGentrueSel[k] = (TH1D*)h3DGentrueSel->ProjectionX(Form("fHistGenSelTrue%i%s",k,lCascType.Data()),
                    minmultbintrue, maxmultbintrue, mineebintrue, maxeebintrue );

            fHistptSel[k] = new TH1D(Form("hsgnloss_%i_%i",(int)percentile[k],(int)percentile[k+1]),"Signal loss corr in mult bins;#it{p}_{T} (GeV/#it{c});Counts", ptbinnumb, ptbinlimits);
            fHistpttrueSel[k] = new TH1D(Form("hsgnlosstrue_%i_%i",(int)percentile[k],(int)percentile[k+1]),"Signal loss corr in mult bins;#it{p}_{T} (GeV/#it{c});Counts", ptbinnumb, ptbinlimits);
        }

        //Rebinning
        temppt = 0;
        for (int k = 0; k<percbinnumb; k++){
            for(long i = 1; i<fHistGenSel[k]->GetNbinsX()+1;i++){
                temppt = fHistGenSel[k]->GetXaxis()->GetBinCenter(i);
                for(long filling = 0; filling<fHistGenSel[k]->GetBinContent(i); filling++){
                    fHistptSel[k]->Fill(temppt);
                }
            }
        }

        
        temppt2 = 0;
        for (int k = 0; k<percbinnumb; k++){
            for(long i = 1; i<fHistGentrueSel[k]->GetNbinsX()+1;i++){
                temppt2 = fHistGentrueSel[k]->GetXaxis()->GetBinCenter(i);
                for(long filling = 0; filling<fHistGentrueSel[k]->GetBinContent(i); filling++){
                    fHistpttrueSel[k]->Fill(temppt2);
                }
            }
        }

        for (int k = 0; k<percbinnumb; k++){
            ErrorInRatioHisto(fHistptSel[k],fHistpttrueSel[k]);
        }
    }

    ////////////////////////////////////////////////////////
    // -------------- Write on file ------------------------
    ////////////////////////////////////////////////////////

    TFile* Write = new TFile (outputname,"RECREATE");
    
    TDirectoryFile* epsevt;
    if (DoEventLoss){    
        epsevt = new TDirectoryFile("EventLoss","EventLoss");
        epsevt->cd();
        hevtloss->Write();
    }
    Write->cd();
    //
    TDirectoryFile* epspart;
    if (DoSignalLoss){    
        epspart = new TDirectoryFile("SignalLoss","SignalLoss");
        epspart->cd();
        for (int k = 0; k<percbinnumb; k++){
            fHistptSel[k]->Write();
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
  //Regular Division 
  h1 -> Divide( h2 ); 
  //Replace Erorrs 
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){ 
    h1->SetBinError(i, lSigmaDelta[i]);
  }  
}