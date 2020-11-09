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

//#define max(a,b) (a>b ? a : b)

double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
void DivideAndComputeRogerBarlow( TH1D* h1, TH1D *h2 );
Double_t GetEEfromZDC ( TFile* lfilename , Float_t lZPCpp, Float_t lZNCpp, Float_t lZPApp, Float_t lZNApp, Int_t lRun );

void ComputeEfficiencies(Bool_t kDoMult = kTRUE, Bool_t kDoEE = kTRUE, TString lCascType = "XiMinus"){

	TFile* file = new TFile("~/Scaricati/VtxCorrect.root","READ");
		//"MCLHC17j_GP.root"
	const char* outputname = "EventCountLoss_15fkINT7.root";

	cout<<"--------------- Open Real Data File --------------------"<<endl;
    TList* clist      = (TList*)file->Get("PWGLF_StrVsMult_MC/cList");
    TTree* lTreeEvent = (TTree*)file->Get("PWGLF_StrVsMult_MC/fTreeEvent");

    //Variables
    Double_t lLoMultBound =   0.;
    Double_t lHiMultBound =   100.;
    Double_t lLoEEBound =   0.;
    Double_t lHiEEBound =   0.;
    Float_t fCentrality = 0.;
    Bool_t fEvSel_INELgtZERO = kFALSE;
    Bool_t fEvSel_INELgtZEROtrue = kFALSE;
    Bool_t fIsINELgtZERO = kFALSE;
    Float_t ZPCpp, ZNCpp, ZPApp, ZNApp;
    Int_t fRun;
    Float_t fZDCCentrality = 0.;
    Bool_t fIsAcceptedVertexPosition = kFALSE;
    Bool_t fHasNoInconsistentSPDandTrackVertices = kFALSE;
    Bool_t fIsSelectedTrigger= kFALSE ;
    Bool_t fAllSelections = kFALSE;

    //Get from Tree
    lTreeEvent->SetBranchAddress("fZPCpp",&ZPCpp);
    lTreeEvent->SetBranchAddress("fZPApp",&ZPApp);
    lTreeEvent->SetBranchAddress("fZNCpp",&ZNCpp);
    lTreeEvent->SetBranchAddress("fZNApp",&ZNApp);
    lTreeEvent->SetBranchAddress("fRun",&fRun);
    lTreeEvent->SetBranchAddress("fCentrality", &fCentrality);
    lTreeEvent->SetBranchAddress("fIsINELgtZERO", &fIsINELgtZERO);
    lTreeEvent->SetBranchAddress("fEvSel_INELgtZERO", &fEvSel_INELgtZERO);
    lTreeEvent->SetBranchAddress("fEvSel_INELgtZEROtrue", &fEvSel_INELgtZEROtrue);
	lTreeEvent->SetBranchAddress("fIsAcceptedVertexPosition", &fIsAcceptedVertexPosition); 
    lTreeEvent->SetBranchAddress("fIsSelectedTrigger", &fIsSelectedTrigger); 
    lTreeEvent->SetBranchAddress("fHasNoInconsistentSPDandTrackVertices", &fHasNoInconsistentSPDandTrackVertices); 
	
    
    TH3D* h3DGenXiINELgt0V0Lambda = (TH3D*)clist->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", "Lambda"));
    TH2D* h2DGenXitrueINELgt0V0Lambda = (TH2D*)clist->FindObject(Form("fHistPtVsCentV0M_Gen%s", "Lambda"));  
    TH2D* h2DGenXiINELgt0ZDCLambda = (TH2D*)clist->FindObject(Form("fHistGeneratedPtVsZDC%s", "Lambda"));
    TH2D* h2DGenXitrueINELgt0ZDCLambda = (TH2D*)clist->FindObject(Form("fHistPtVsZDC_Gen%s", "Lambda"));
	
    TH3D* h3DGenXiINELgt0V0 = (TH3D*)clist->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", "XiMinus"));
    TH2D* h2DGenXitrueINELgt0V0 = (TH2D*)clist->FindObject(Form("fHistPtVsCentV0M_Gen%s", "XiMinus"));  
    TH2D* h2DGenXiINELgt0ZDC = (TH2D*)clist->FindObject(Form("fHistGeneratedPtVsZDC%s", "XiMinus"));
    TH2D* h2DGenXitrueINELgt0ZDC = (TH2D*)clist->FindObject(Form("fHistPtVsZDC_Gen%s", "XiMinus"));

    TH3D* h3DGenXiINELgt0V0Omega = (TH3D*)clist->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", "OmegaMinus"));
    TH2D* h2DGenXitrueINELgt0V0Omega = (TH2D*)clist->FindObject(Form("fHistPtVsCentV0M_Gen%s", "OmegaMinus"));  
    TH2D* h2DGenXiINELgt0ZDCOmega = (TH2D*)clist->FindObject(Form("fHistGeneratedPtVsZDC%s", "OmegaMinus"));
    TH2D* h2DGenXitrueINELgt0ZDCOmega = (TH2D*)clist->FindObject(Form("fHistPtVsZDC_Gen%s", "OmegaMinus"));

	//-----------------------------------------------------------------------------------------------------
	//-------------------------------- Event loss correction ---------------------------------------------- 
    //-----------------------------------------------------------------------------------------------------
	
	//Mult variables
    Long_t lNEventsINELgt0V0[10];
    Long_t lIsAcceptedVertexPosition[10];
    Long_t lHasNoInconsistentSPDandTrackVertices[10];
    Long_t lIsSelectedTrigger[10];
    Long_t lNEventstrueINELgt0V0[10];
    Long_t lAllSelections[10];
    Float_t mult[11] = {0,1,5,10,15,20,30,40,50,70,100};
    TH1F* hevtlossV0 = new TH1F("hevtlossV0", "", 10, mult);
    TH1F* hevtlossSPD = new TH1F("hevtlossSPD", "", 10, mult);
    TH1F* hevtlossVtx = new TH1F("hevtlossVtx", "", 10, mult);
    TH1F* hevtlossTrig = new TH1F("hevtlossTrig", "", 10, mult);
    TH1F* hevtlossAll = new TH1F("hevtlossAll", "", 10, mult);
    //ZDC variables
    Long_t lNEventsINELgt0ZDC[10];
    Long_t lNEventstrueINELgt0ZDC[10];
    Float_t ee[11] = {0,10,20,30,40,50,60,70,80,90,100};    
    TH1F* hevtlossZDC = new TH1F("hevtlossZDC", "", 10, ee);
    TFile* Read;
    TFile* ZDCpartfile;
   
    Double_t ptbinlimits[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
    Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;

    for (int k=0; k<10;k++){

    	if (kDoMult){
    		lLoMultBound=mult[k];
    		lHiMultBound=mult[k+1];
    		lNEventsINELgt0V0[k]=0;
       		lNEventstrueINELgt0V0[k]=0;
            lIsAcceptedVertexPosition[k]=0;
            lHasNoInconsistentSPDandTrackVertices[k]=0;
            lIsSelectedTrigger[k]=0;
            lAllSelections[k]=0;

       		cout << " Multiplicity is in range: [" << lLoMultBound << " , " << lHiMultBound << "]\n" << endl;
    	}
    	  	
        
      	cout << "Event counters are at " << lNEventsINELgt0V0[k] << "  ";
      	if (lNEventsINELgt0V0[k] == 0) cout << "as expected" << "\n" << endl;

      	cout<<" \nWill now loop over events, please wait...\n"<<endl;
    	for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries(); iEv++) {

        	lTreeEvent->GetEntry(iEv);
        	if( iEv % ( lTreeEvent->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;
        	// check distance to closest non empty BC
        
            //if (fRun == 226476) continue;
    

	 		//Count events effective energy
	    	if (kDoMult){
		        if( fIsINELgtZERO &&
		        	fCentrality>=lLoMultBound &&
		            fCentrality<lHiMultBound     
		          ) lNEventsINELgt0V0[k]++;
                if( fIsSelectedTrigger&&
                    fCentrality>=lLoMultBound &&
                    fCentrality<lHiMultBound     
                  ) lIsSelectedTrigger[k]++;
                if( fIsAcceptedVertexPosition &&
                    fCentrality>=lLoMultBound &&
                    fCentrality<lHiMultBound     
                  ) lIsAcceptedVertexPosition[k]++;
                if( fHasNoInconsistentSPDandTrackVertices &&
                    fCentrality>=lLoMultBound &&
                    fCentrality<lHiMultBound     
                  ) lHasNoInconsistentSPDandTrackVertices[k]++;
                if( fEvSel_INELgtZERO &&
                    fCentrality>=lLoMultBound &&
                    fCentrality<lHiMultBound     
                  ) lAllSelections[k]++;
		        if( fEvSel_INELgtZEROtrue &&
		        	fCentrality>=lLoMultBound &&
		            fCentrality<lHiMultBound  
		          ) lNEventstrueINELgt0V0[k]++;

		    }

	    }
	    if (kDoMult){
	        cout<<"\n Number of events reco INEL>0, this multiplicity selection....: "<<lNEventsINELgt0V0[k] <<endl;
	        cout<<" Number of events true INEL>0, this multiplicity selection....: "<<lNEventstrueINELgt0V0[k] <<endl;
	        cout << " Event Loss Correction Mult Sel... = " << (Float_t)lNEventsINELgt0V0[k]/lNEventstrueINELgt0V0[k] << endl;
	        cout<<" --------------------------------------------------------"<<endl;
	    }
         
    }
   
    if (kDoMult){
	    for (int i=1; i<= hevtlossV0->GetNbinsX(); i++){

	    	hevtlossV0->SetBinContent(i,(Float_t)lNEventsINELgt0V0[i-1]/lNEventstrueINELgt0V0[i-1]);
	    	hevtlossV0->SetBinError(i,ErrorInRatio((Float_t)lNEventsINELgt0V0[i-1],(Float_t)TMath::Sqrt(lNEventsINELgt0V0[i-1]),
	    							(Float_t)lNEventstrueINELgt0V0[i-1],(Float_t)TMath::Sqrt(lNEventstrueINELgt0V0[i-1])));
            hevtlossSPD->SetBinContent(i,(Float_t)lHasNoInconsistentSPDandTrackVertices[i-1]/lNEventstrueINELgt0V0[i-1]);
            hevtlossVtx->SetBinContent(i,(Float_t)lIsAcceptedVertexPosition[i-1]/lNEventstrueINELgt0V0[i-1]);
            hevtlossTrig->SetBinContent(i,(Float_t)lIsSelectedTrigger[i-1]/lNEventstrueINELgt0V0[i-1]);
            hevtlossAll->SetBinContent(i,(Float_t)lAllSelections[i-1]/lNEventstrueINELgt0V0[i-1]);
         }   
	}

	
	//-----------------------------------------------------------------------------------------------------
	//-------------------------------- Signal loss correction --------------------------------------------- 
    //-----------------------------------------------------------------------------------------------------
	
    TH1D* fHistGentrueINELgt0V0[10];
    TH1D* fHistGenINELgt0V0[10];
    TH1D* fHistptINELgt0V0[10];
    TH1D* fHistpttrueINELgt0V0[10];

    TH1D* fHistGentrueINELgt0ZDC[10];
    TH1D* fHistGenINELgt0ZDC[10];
    TH1D* fHistptINELgt0ZDC[10];
    TH1D* fHistpttrueINELgt0ZDC[10];

    Int_t minrapbin = h3DGenXiINELgt0V0->GetYaxis()->FindBin( -0.5+1e-6 );
  	Int_t maxrapbin = h3DGenXiINELgt0V0->GetYaxis()->FindBin( +0.5-1e-6 );
    Int_t minmultbinV0_h3 = 0;
    Int_t maxmultbinV0_h3 = 0;
    Int_t minmultbinV0_h2 = 0;
    Int_t maxmultbinV0_h2 = 0;
    Int_t mineebinZDC = 0;
    Int_t maxeebinZDC = 0;
  
  

    for (int i=1; i< 11; i++){

   	    //Multiplicity
        minmultbinV0_h3 = h3DGenXiINELgt0V0->GetZaxis()->FindBin( mult[i-1]+1e-6 );
    	maxmultbinV0_h3 = h3DGenXiINELgt0V0->GetZaxis()->FindBin( mult[i]-1e-6 );
    	minmultbinV0_h2 = h2DGenXitrueINELgt0V0->GetYaxis()->FindBin( mult[i-1]+1e-6 );
    	maxmultbinV0_h2 = h2DGenXitrueINELgt0V0->GetYaxis()->FindBin( mult[i]-1e-6 );

		fHistGenINELgt0V0[i-1] = (TH1D*)h3DGenXiINELgt0V0->ProjectionX(Form("fHistGenINELgt0V0%i",i), minrapbin, maxrapbin, minmultbinV0_h3, maxmultbinV0_h3);
    	fHistGentrueINELgt0V0[i-1] = (TH1D*)h2DGenXitrueINELgt0V0->ProjectionX(Form("fHistGentrueINELgt0V0%i",i), minmultbinV0_h2, maxmultbinV0_h2);

    	fHistptINELgt0V0[i-1]  = new TH1D(Form("fHistptINELgt0V0_%i-%i",(int)mult[i-1],(int)mult[i]),"Cascade MC count;p_{T} (GeV/c);Counts", ptbinnumb, ptbinlimits);
        fHistpttrueINELgt0V0[i-1] = new TH1D(Form("fHistpttrueINELgt0V0_%i-%i",(int)mult[i-1],(int)mult[i]),"Cascade MC count;p_{T} (GeV/c);Counts", ptbinnumb, ptbinlimits);

  } 

    //Do V0 selection
    Double_t temppt;
    for (int k = 0; k<10; k++){
	    for(long i = 1; i<fHistGenINELgt0V0[k]->GetNbinsX()+1;i++){
	        temppt = fHistGenINELgt0V0[k]->GetXaxis()->GetBinCenter(i);
	        for(long filling = 0; filling<fHistGenINELgt0V0[k]->GetBinContent(i); filling++){
	            fHistptINELgt0V0[k]->Fill(temppt);
	        }
	    }
	}

	
	Double_t temppt2;
	for (int k = 0; k<10; k++){
	    for(long i = 1; i<fHistGentrueINELgt0V0[k]->GetNbinsX()+1;i++){
	        temppt2 = fHistGentrueINELgt0V0[k]->GetXaxis()->GetBinCenter(i);
	        for(long filling = 0; filling<fHistGentrueINELgt0V0[k]->GetBinContent(i); filling++){
	            fHistpttrueINELgt0V0[k]->Fill(temppt2);
	        }
	    }
	}

	for (int k = 0; k<10; k++){
		DivideAndComputeRogerBarlow(fHistptINELgt0V0[k],fHistpttrueINELgt0V0[k]);
	}


    
	
   	TFile* Write = new TFile (outputname,"RECREATE");
    hevtlossV0->Write();
        hevtlossSPD->Write();
        hevtlossVtx->Write();
        hevtlossTrig->Write();
        hevtlossAll->Write();
            
    
    for (int k = 0; k<10; k++){
    	fHistptINELgt0V0[k]->Write();
    }
}


//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0 ) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop - errorfrombottom) );
    }
    return 1.;
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
Double_t GetEEfromZDC ( TFile* lfilename , Float_t lZPCpp, Float_t lZNCpp, Float_t lZPApp, Float_t lZNApp, Int_t lRun ){ //set in the run func
    //Converts ZDC information to effective energy percentiles. 
    
    TH1F * fhCumulative = (TH1F *)lfilename->Get(Form("%i/hCumulative_%i",lRun, lRun));
  
    Double_t Sum = TMath::Log10(TMath::Abs(lZPCpp+lZNCpp+lZPApp+lZNApp));
    Double_t fZDCCentrality;

    if (Sum!=0) {
        fZDCCentrality = 100*(fhCumulative->Interpolate(Sum));
    }
    else fZDCCentrality=-1; 

    return fZDCCentrality;
}
