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

void ComputeEfficienciesFullXiFINAL(){

	TFile* file = new TFile("MCLHC15g3c3.root","READ");
		//"MCLHC17j_GP.root"
	const char* ZDCFilename = "ExtractedZDCPercentile_MC-15g3c3.root"; //""
    const char* ZDCpartfilename = ZDCFilename; //""
	const char* outputname = "EventCountLoss_15g3c3_MB.root";
    TFile* Read = new TFile (ZDCFilename);
    TFile* ZDCpartfile = new TFile (ZDCpartfilename);
    cout << " Effective Energy File: ..." <<  ZDCFilename << "\n\n" << endl;

	cout<<"--------------- Open Real Data File --------------------"<<endl;
    TList* clist      = (TList*)file->Get("PWGLF_StrVsMult_MC/cList");
    TTree* lTreeEvent = (TTree*)file->Get("PWGLF_StrVsMult_MC/fTreeEvent");

    //Variables
    Double_t lLoMultBound =   0.;
    Double_t lHiMultBound =   100.;
    Double_t lLoEEBound =   0.;
    Double_t lHiEEBound =   0.;
    Float_t fCentrality = 0.;
    Bool_t fEvSel_AllSelections = kFALSE;
    Bool_t fEvSel_INELgtZEROtrue = kFALSE;
    Bool_t fEvSel_zVtxZMC = kFALSE;
    Float_t ZPCpp, ZNCpp, ZPApp, ZNApp;
    Int_t fRun;
    Float_t fZDCCentrality = 0.;

    //Get from Tree
    lTreeEvent->SetBranchAddress("fZPCpp",&ZPCpp);
    lTreeEvent->SetBranchAddress("fZPApp",&ZPApp);
    lTreeEvent->SetBranchAddress("fZNCpp",&ZNCpp);
    lTreeEvent->SetBranchAddress("fZNApp",&ZNApp);
    lTreeEvent->SetBranchAddress("fRun",&fRun);
    lTreeEvent->SetBranchAddress("fCentrality", &fCentrality);
    lTreeEvent->SetBranchAddress("fEvSel_zVtxZMC", &fEvSel_zVtxZMC);
    lTreeEvent->SetBranchAddress("fEvSel_AllSelections", &fEvSel_AllSelections);
    lTreeEvent->SetBranchAddress("fEvSel_INELgtZEROtrue", &fEvSel_INELgtZEROtrue);
    
 	//-----------------------------------------------------------------------------------------------------
	//-------------------------------- Event loss correction ---------------------------------------------- 
    //-----------------------------------------------------------------------------------------------------

	//Multiplicity Efficiency
    Float_t mult[] = {0,100};
    Long_t multbinnumb = sizeof(mult)/sizeof(Float_t) - 1;
    Long_t lNEvtSelectedV0[multbinnumb];
    Long_t lNEvttrueSelectedV0[multbinnumb];    

    //Effective energy Efficiency
    Float_t ee[] = {0,100}; 
    Long_t eebinnumb = sizeof(ee)/sizeof(Float_t) - 1;
    Long_t lNEvtSelectedZDC[eebinnumb];
    Long_t lNEvttrueSelectedZDC[eebinnumb];   

    //Double Differential study
    //Fixed Energy 
    Double_t FixedHighEE[2] = {0.,30.};
    Double_t FixedLowEE[2] = {70.,100};
    Long_t lNEvtSelectedFixLowEE[multbinnumb];
    Long_t lNEvttrueSelectedFixLowEE[multbinnumb];
    Long_t lNEvtSelectedFixHighEE[multbinnumb];
    Long_t lNEvttrueSelectedFixHighEE[multbinnumb];
    
    //Fixed Multiplicity
    Double_t FixedHighmult[2] = {0.,30.};
    Double_t FixedLowmult[2] = {70.,100};
    Long_t lNEvtSelectedFixLowmult[eebinnumb];
    Long_t lNEvttrueSelectedFixLowmult[eebinnumb];
    Long_t lNEvtSelectedFixHighmult[eebinnumb];
    Long_t lNEvttrueSelectedFixHighmult[eebinnumb];
    
    //Histograms
    TH1F* hevtlossV0 = new TH1F("hevtlossV0", "", multbinnumb, mult);
    hevtlossV0->SetTitle("Event loss correction in mult bins");
    TH1F* hevtlossZDC = new TH1F("hevtlossZDC", "", eebinnumb, ee);
    hevtlossZDC->SetTitle("Event loss correction in eff energy bins");
    TH1F* hevtlossFixLowEE = new TH1F("hevtlossV0FixLowEE", "", multbinnumb, mult);
    hevtlossFixLowEE->SetTitle(Form("Event loss correction in mult bins - ZDCperc is fixed in [%.0f,%.0f]",FixedLowEE[0],FixedLowEE[1]));
    TH1F* hevtlossFixHighEE = new TH1F("hevtlossV0FixHighEE", "", multbinnumb, mult);
    hevtlossFixHighEE->SetTitle(Form("Event loss correction in mult bins - ZDCperc is fixed in [%.0f,%.0f]",FixedHighEE[0],FixedHighEE[1]));   
    TH1F* hevtlossFixLowmult = new TH1F("hevtlossV0FixLowmult", "", eebinnumb, ee);
    hevtlossFixLowmult->SetTitle(Form("Event loss correction in eff energy bins - V0Mperc is fixed in [%.0f,%.0f]",FixedLowmult[0],FixedLowmult[1]));
    TH1F* hevtlossFixHighmult = new TH1F("hevtlossV0FixHighmult", "", eebinnumb, ee);
    hevtlossFixHighmult->SetTitle(Form("Event loss correction in eff energy bins - V0Mperc is fixed in [%.0f,%.0f]",FixedHighmult[0],FixedHighmult[1]));


    //Pt bins definition
    Double_t ptbinlimits[] = {0.0, 0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5, 8.0, 10.0 };  
    Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;
   
    for (int k=0; k<multbinnumb; k++)
    {
        //Mult and EE ranges
   		lLoMultBound=mult[k];
   		lHiMultBound=mult[k+1];
       	cout << " Multiplicity is in range: [" << lLoMultBound << " , " << lHiMultBound << "]\n" << endl;
    	
    	
        lLoEEBound=ee[k];
        lHiEEBound=ee[k+1];    	
        cout << " Effective energy is in range: [" << lLoEEBound << " , " << lHiEEBound << "]\n" << endl;
    

        //Counters must be initialized to 0
        lNEvtSelectedV0[k]=0;
        lNEvttrueSelectedV0[k]=0;
        lNEvtSelectedFixLowEE[k]=0;
        lNEvttrueSelectedFixLowEE[k]=0;
        lNEvtSelectedFixHighEE[k]=0;
        lNEvttrueSelectedFixHighEE[k]=0;
        
            lNEvtSelectedZDC[k]=0;
            lNEvttrueSelectedZDC[k]=0;   
            lNEvtSelectedFixLowmult[k]=0;
            lNEvttrueSelectedFixLowmult[k]=0;
            lNEvtSelectedFixHighmult[k]=0;
            lNEvttrueSelectedFixHighmult[k]=0;
        
      	 
      	cout<<" \nWill now loop over events, please wait...\n"<<endl;
    	for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries(); iEv++) {

        	lTreeEvent->GetEntry(iEv);
        	if( iEv % ( lTreeEvent->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;
  
          	//Get ZDCPercentile
        	fZDCCentrality = GetEEfromZDC(Read, ZPCpp, ZNCpp, ZPApp, ZNApp, fRun); 
            
            //Count events EE
            if( 
                fEvSel_AllSelections &&
	            fZDCCentrality>=lLoEEBound &&
	            fZDCCentrality<lHiEEBound     
	            ) 
                lNEvtSelectedZDC[k]++;
	        if( 
                fEvSel_INELgtZEROtrue && 
                fEvSel_zVtxZMC &&
	            fZDCCentrality>=lLoEEBound &&
	            fZDCCentrality<lHiEEBound  
	            ) 
                lNEvttrueSelectedZDC[k]++;
           
	 		//Count events mult
	    	if( 
                fEvSel_AllSelections &&
		        fCentrality>=lLoMultBound &&
		        fCentrality<lHiMultBound     
		      ) 
                lNEvtSelectedV0[k]++;
		    if( 
                fEvSel_INELgtZEROtrue && 
                fEvSel_zVtxZMC &&
		       	fCentrality>=lLoMultBound &&
		        fCentrality<lHiMultBound  
		      ) 
                lNEvttrueSelectedV0[k]++;

            //Double Differential
            //Fixed Mult Low
           if( 
                fEvSel_AllSelections &&
                fCentrality<FixedLowmult[1] &&
                fCentrality>=FixedLowmult[0] &&
                fZDCCentrality>=lLoEEBound &&
                fZDCCentrality<lHiEEBound     
                ) 
                lNEvtSelectedFixLowmult[k]++;
            if( 
                fEvSel_INELgtZEROtrue && 
                fEvSel_zVtxZMC &&
                fCentrality<FixedLowmult[1] &&
                fCentrality>=FixedLowmult[0] &&
                fZDCCentrality>=lLoEEBound &&
                fZDCCentrality<lHiEEBound  
                ) 
                lNEvttrueSelectedFixLowmult[k]++;
            
            //Fixed Multiplicity High
            if( 
                fEvSel_AllSelections &&
                fCentrality<FixedHighmult[1] &&
                fCentrality>=FixedHighmult[0] &&
                fZDCCentrality>=lLoEEBound &&
                fZDCCentrality<lHiEEBound     
                ) 
                lNEvtSelectedFixHighmult[k]++;
            if( 
                fEvSel_INELgtZEROtrue && 
                fEvSel_zVtxZMC &&
                fCentrality<FixedHighmult[1] &&
                fCentrality>=FixedHighmult[0] &&
                fZDCCentrality>=lLoEEBound &&
                fZDCCentrality<lHiEEBound  
                ) 
                lNEvttrueSelectedFixHighmult[k]++;
           

            //Fixed Low EE
            if( 
                fEvSel_AllSelections &&
                fZDCCentrality<FixedLowEE[1] &&
                fZDCCentrality>=FixedLowEE[0] &&
                fCentrality>=lLoMultBound &&
                fCentrality<lHiMultBound     
              ) 
                lNEvtSelectedFixLowEE[k]++;
            if( 
                fEvSel_INELgtZEROtrue && 
                fEvSel_zVtxZMC &&
                fZDCCentrality<FixedLowEE[1] &&
                fZDCCentrality>=FixedLowEE[0] &&
                fCentrality>=lLoMultBound &&
                fCentrality<lHiMultBound  
              ) 
                lNEvttrueSelectedFixLowEE[k]++;
		    //Fixed High EE
            if( 
                fEvSel_AllSelections &&
                fZDCCentrality<FixedHighEE[1] &&
                fZDCCentrality>=FixedHighEE[0] &&
                fCentrality>=lLoMultBound &&
                fCentrality<lHiMultBound     
              ) 
                lNEvtSelectedFixHighEE[k]++;
            if( 
                fEvSel_INELgtZEROtrue && 
                fEvSel_zVtxZMC &&
                fZDCCentrality<FixedHighEE[1] &&
                fZDCCentrality>=FixedHighEE[0] &&
                fCentrality>=lLoMultBound &&
                fCentrality<lHiMultBound  
              ) 
                lNEvttrueSelectedFixHighEE[k]++;

        }

	    cout<<"\n Number of events reco INEL>0, this multiplicity selection....: "<<lNEvtSelectedV0[k] <<endl;
	    cout<<" Number of events true INEL>0, this multiplicity selection....: "<<lNEvttrueSelectedV0[k] <<endl;
	    cout << " Event Loss Correction Mult Sel... = " << (Float_t)lNEvtSelectedV0[k]/lNEvttrueSelectedV0[k] << endl;
	    cout<<" --------------------------------------------------------"<<endl;
	   
	    
        cout<<"\n Number of events reco INEL>0, this effective energy selection....: "<<lNEvtSelectedZDC[k] <<endl;
        cout<<" Number of events true INEL>0, this effective energy selection....: "<<lNEvttrueSelectedZDC[k] <<endl;
        cout << " Event Loss Correction EE Sel ... = " << (Float_t)lNEvtSelectedZDC[k]/lNEvttrueSelectedZDC[k] << endl;
        cout<<" --------------------------------------------------------"<<endl;

            
    }
   
    for (int i=1; i<= hevtlossV0->GetNbinsX(); i++){

        cout << (Float_t)lNEvtSelectedV0[i-1]/lNEvttrueSelectedV0[i-1] << endl;
        cout << (Float_t)lNEvtSelectedV0[i-1] << "   " << lNEvttrueSelectedV0[i-1] << endl;
	    hevtlossV0->SetBinContent(i,(Float_t)lNEvtSelectedV0[i-1]/lNEvttrueSelectedV0[i-1]);
	    hevtlossV0->SetBinError(i,ErrorInRatio(
            (Float_t)lNEvtSelectedV0[i-1],(Float_t)TMath::Sqrt(lNEvtSelectedV0[i-1]),
            (Float_t)lNEvttrueSelectedV0[i-1],(Float_t)TMath::Sqrt(lNEvttrueSelectedV0[i-1]))
            );
        hevtlossFixLowEE->SetBinContent(i,(Float_t)lNEvtSelectedFixLowEE[i-1]/lNEvttrueSelectedFixLowEE[i-1]);
        hevtlossFixLowEE->SetBinError(i,ErrorInRatio(
            (Float_t)lNEvtSelectedFixLowEE[i-1],(Float_t)TMath::Sqrt(lNEvtSelectedFixLowEE[i-1]),
            (Float_t)lNEvttrueSelectedFixLowEE[i-1],(Float_t)TMath::Sqrt(lNEvttrueSelectedFixLowEE[i-1]))
            );
        hevtlossFixHighEE->SetBinContent(i,(Float_t)lNEvtSelectedFixHighEE[i-1]/lNEvttrueSelectedFixHighEE[i-1]);
        hevtlossFixHighEE->SetBinError(i,ErrorInRatio(
            (Float_t)lNEvtSelectedFixHighEE[i-1],(Float_t)TMath::Sqrt(lNEvtSelectedFixHighEE[i-1]),
            (Float_t)lNEvttrueSelectedFixHighEE[i-1],(Float_t)TMath::Sqrt(lNEvttrueSelectedFixHighEE[i-1]))
            );
	}   

    for (int i=1; i<= hevtlossZDC->GetNbinsX(); i++){
        
        cout << (Float_t)lNEvtSelectedZDC[i-1]/lNEvtSelectedZDC[i-1] << endl;
        cout << (Float_t)lNEvtSelectedZDC[i-1] << "   " << lNEvtSelectedZDC[i-1] << endl;
	
        hevtlossZDC->SetBinContent(i,(Float_t)lNEvtSelectedZDC[i-1]/lNEvttrueSelectedZDC[i-1]);
        hevtlossZDC->SetBinError(i,ErrorInRatio(
            (Float_t)lNEvtSelectedZDC[i-1],(Float_t)TMath::Sqrt(lNEvtSelectedZDC[i-1]),
            (Float_t)lNEvttrueSelectedZDC[i-1],(Float_t)TMath::Sqrt(lNEvttrueSelectedZDC[i-1]))
            );	  
        hevtlossFixLowmult->SetBinContent(i,(Float_t)lNEvtSelectedFixLowmult[i-1]/lNEvttrueSelectedFixLowmult[i-1]);
        hevtlossFixLowmult->SetBinError(i,ErrorInRatio(
            (Float_t)lNEvtSelectedFixLowmult[i-1],(Float_t)TMath::Sqrt(lNEvtSelectedFixLowmult[i-1]),
            (Float_t)lNEvttrueSelectedFixLowmult[i-1],(Float_t)TMath::Sqrt(lNEvttrueSelectedFixLowmult[i-1]))
            );
        hevtlossFixHighmult->SetBinContent(i,(Float_t)lNEvtSelectedFixHighmult[i-1]/lNEvttrueSelectedFixHighmult[i-1]);
        hevtlossFixHighmult->SetBinError(i,ErrorInRatio(
            (Float_t)lNEvtSelectedFixHighmult[i-1],(Float_t)TMath::Sqrt(lNEvtSelectedFixHighmult[i-1]),
            (Float_t)lNEvttrueSelectedFixHighmult[i-1],(Float_t)TMath::Sqrt(lNEvttrueSelectedFixHighmult[i-1]))
            );  
	}	 


	//-----------------------------------------------------------------------------------------------------
	//-------------------------------- Signal loss correction --------------------------------------------- 
    //-----------------------------------------------------------------------------------------------------

    //GetHistos for Eps part
    const int npart = 2;
    TString parttype[npart] = {"XiMinus","XiPlus"};

    TH3D* h3DGenSelV0[npart];
    TH2D* h2DGentrueSelV0[npart];   
    TH2D* h2DGenSelZDC[npart]; 
    TH2D* h2DGentrueSelZDC[npart];
    TH3D* h3DGenSelDD[npart]; 
    TH3D* h3DGentrueSelDD[npart]; 
    
    for (int type = 0; type < npart; type++){
       h3DGenSelV0[type] = (TH3D*)clist->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", parttype[type].Data()));
       h3DGenSelV0[type]->Sumw2();
       h2DGentrueSelV0[type] = (TH2D*)clist->FindObject(Form("fHistPtVsCentV0M_Gen%s", parttype[type].Data()));
       h2DGentrueSelV0[type]->Sumw2();
       h2DGenSelZDC[type] = (TH2D*)clist->FindObject(Form("fHistGeneratedPtVsZDC%s", parttype[type].Data()));
       h2DGenSelZDC[type]->Sumw2();
       h2DGentrueSelZDC[type] = (TH2D*)clist->FindObject(Form("fHistPtVsZDC_Gen%s", parttype[type].Data()));
       h2DGentrueSelZDC[type]->Sumw2();
       h3DGenSelDD[type] = (TH3D*)clist->FindObject(Form("fHistGeneratedPtVsZDCVsCentV0M%s", parttype[type].Data()));
       h3DGenSelDD[type]->Sumw2();
       h3DGentrueSelDD[type] = (TH3D*)clist->FindObject(Form("fHistPtVsZDCVsCentV0M_Gen%s", parttype[type].Data()));
        h3DGentrueSelDD[type]->Sumw2();
    }   
    
    h3DGenSelV0[0]->Add(h3DGenSelV0[1]);
    h2DGentrueSelV0[0]->Add(h2DGentrueSelV0[1]) ;
    h2DGenSelZDC[0]->Add(h2DGenSelZDC[1]) ;
    h2DGentrueSelZDC[0]->Add(h2DGentrueSelZDC[1]);
    h3DGenSelDD[0]->Add(h3DGenSelDD[1]) ;
    h3DGentrueSelDD[0]->Add(h3DGentrueSelDD[1]) ;
  
    //Now get things right for different particles 
    Double_t *partptbinlimits = 0;
    Long_t partptbinnumb = 0;
    Float_t *partmult = 0;
    Long_t partmultbinnumb = 0;
    Float_t *partee = 0;
    Long_t parteebinnumb = 0;
    Double_t temppt;

    //Multiplicity
    TH1D* fHistGentrueSelV0[multbinnumb][npart];
    TH1D* fHistGenSelV0[multbinnumb][npart];
    TH1D* fHistptSelV0[multbinnumb][npart];
    TH1D* fHistpttrueSelV0[multbinnumb][npart];

    //Eff energy
    TH1D* fHistGentrueSelZDC[eebinnumb][npart];
    TH1D* fHistGenSelZDC[eebinnumb][npart];
    TH1D* fHistptSelZDC[eebinnumb][npart];
    TH1D* fHistpttrueSelZDC[eebinnumb][npart];

    //Double Diff
    TH1D* fHistGentrueSelV0FixHighEE[multbinnumb][npart];
    TH1D* fHistGenSelV0FixHighEE[multbinnumb][npart];
    TH1D* fHistptSelV0FixHighEE[multbinnumb][npart];
    TH1D* fHistpttrueSelV0FixHighEE[multbinnumb][npart];

    TH1D* fHistGentrueSelV0FixLowEE[multbinnumb][npart];
    TH1D* fHistGenSelV0FixLowEE[multbinnumb][npart];
    TH1D* fHistptSelV0FixLowEE[multbinnumb][npart];
    TH1D* fHistpttrueSelV0FixLowEE[multbinnumb][npart];

    TH1D* fHistGentrueSelZDCFixHighmult[eebinnumb][npart];
    TH1D* fHistGenSelZDCFixHighmult[eebinnumb][npart];
    TH1D* fHistptSelZDCFixHighmult[eebinnumb][npart];
    TH1D* fHistpttrueSelZDCFixHighmult[eebinnumb][npart];

    TH1D* fHistGentrueSelZDCFixLowmult[eebinnumb][npart];
    TH1D* fHistGenSelZDCFixLowmult[eebinnumb][npart];
    TH1D* fHistptSelZDCFixLowmult[eebinnumb][npart];
    TH1D* fHistpttrueSelZDCFixLowmult[eebinnumb][npart];

    Int_t minrapbin = h3DGenSelV0[0]->GetYaxis()->FindBin( -0.5+1e-6 );
  	Int_t maxrapbin = h3DGenSelV0[0]->GetYaxis()->FindBin( +0.5-1e-6 );
    Int_t minmultbinV0_h3 = 0;
    Int_t maxmultbinV0_h3 = 0;
    Int_t minmultbinV0_h2 = 0;
    Int_t maxmultbinV0_h2 = 0;
    Int_t mineebinZDC = 0;
    Int_t maxeebinZDC = 0;
  
  	//Reverse Calibration for ZDC------------------------------------------------------------
    Double_t LowLogSum, HighLogSum;
    Int_t counter = 0;

    //Get Calibration function
    TH1F * hCumulative = (TH1F *)ZDCpartfile->Get("hCumulative_TOT");
    const int Nbins = hCumulative->GetNbinsX();
   	
    for (int ibin = 0; ibin < Nbins ; ibin++){    	
    	if (hCumulative->GetBinContent(ibin+1) == 0) counter++;    	  	
    }

    const int filledbins = Nbins - counter;
    Float_t RevCum_x[filledbins], RevCum_y[filledbins];

    for (int ibin = counter; ibin < Nbins ; ibin++){    
    	int atbin = ibin-counter;	
    	RevCum_y[atbin] = hCumulative->GetBinCenter(ibin+1);
    	RevCum_x[atbin] = hCumulative->GetBinContent(ibin+1);     	    	  	
    }
    TGraph* ReverseCumulative = new TGraph(filledbins,RevCum_x,RevCum_y);

    //Double diff values
    Double_t LogLowEE[2];
    LogLowEE[0] = ReverseCumulative->Eval(FixedLowEE[0]/100+1e-6);
    LogLowEE[1] = ReverseCumulative->Eval(FixedLowEE[1]/100+1e-6);
    Double_t LogHighEE[2]; 
    LogHighEE[0] = ReverseCumulative->Eval(FixedHighEE[0]/100-1e-6);
    LogHighEE[1] = ReverseCumulative->Eval(FixedHighEE[1]/100-1e-6);


    //Perform epsilon part    
    Double_t fixedlowmineebin = h3DGenSelDD[0]->GetYaxis()->FindBin( LogLowEE[0] );
    Double_t fixedlowmaxeebin = h3DGenSelDD[0]->GetYaxis()->FindBin( LogLowEE[1] );
    Double_t fixedhighmineebin = h3DGenSelDD[0]->GetYaxis()->FindBin( LogHighEE[0] );
    Double_t fixedhighmaxeebin = h3DGenSelDD[0]->GetYaxis()->FindBin( LogHighEE[1] );

    Double_t fixedlowminmultbin = h3DGenSelDD[0]->GetZaxis()->FindBin( FixedLowmult[0] );
    Double_t fixedlowmaxmultbin = h3DGenSelDD[0]->GetZaxis()->FindBin( FixedLowmult[1] );
    Double_t fixedhighminmultbin = h3DGenSelDD[0]->GetZaxis()->FindBin( FixedHighmult[0] );
    Double_t fixedhighmaxmultbin = h3DGenSelDD[0]->GetZaxis()->FindBin( FixedHighmult[1] );

    partptbinlimits = ptbinlimits;
    partptbinnumb = ptbinnumb;
    partmult = mult;
    partmultbinnumb = multbinnumb;
    partee = ee;
    parteebinnumb = eebinnumb;

    for (int i=1; i <= partmultbinnumb; i++){

   	    //Multiplicity
        minmultbinV0_h3 = h3DGenSelV0[0]->GetZaxis()->FindBin( partmult[i-1]+1e-6 );
    	maxmultbinV0_h3 = h3DGenSelV0[0]->GetZaxis()->FindBin( partmult[i]-1e-6 );
    	minmultbinV0_h2 = h2DGentrueSelV0[0]->GetYaxis()->FindBin( partmult[i-1]+1e-6 );
    	maxmultbinV0_h2 = h2DGentrueSelV0[0]->GetYaxis()->FindBin( partmult[i]-1e-6 );

		fHistGenSelV0[i-1][0] = (TH1D*)h3DGenSelV0[0]->ProjectionX(Form("fHistGenSelV0%i%s",i,parttype[0].Data()), minrapbin, maxrapbin, minmultbinV0_h3, maxmultbinV0_h3);
    	fHistGentrueSelV0[i-1][0] = (TH1D*)h2DGentrueSelV0[0]->ProjectionX(Form("fHistGentrueSelV0%i%s",i,parttype[0].Data()), minmultbinV0_h2, maxmultbinV0_h2);

    	fHistptSelV0[i-1][0]  = new TH1D(Form("fHistptSelV0_%i-%i_%s",(int)partmult[i-1],(int)partmult[i],parttype[0].Data()),"Signal loss corr in mult bins;p_{T} (GeV/c);Counts", partptbinnumb, partptbinlimits);
        fHistpttrueSelV0[i-1][0] = new TH1D(Form("fHistpttrueSelV0_%i-%i_%s",(int)partmult[i-1],(int)partmult[i],parttype[0].Data()),"Signal loss corr in mult bins;p_{T} (GeV/c);Counts", partptbinnumb, partptbinlimits);

        //Effective energy
        if (partee[i-1] == 0) LowLogSum = 0;
        else LowLogSum = ReverseCumulative->Eval(partee[i-1]/100+1e-6);
        HighLogSum = ReverseCumulative->Eval(partee[i]/100-1e-6);

        mineebinZDC = h2DGentrueSelZDC[0]->GetYaxis()->FindBin( LowLogSum );
        maxeebinZDC = h2DGentrueSelZDC[0]->GetYaxis()->FindBin( HighLogSum );

        fHistGenSelZDC[i-1][0] = (TH1D*)h2DGenSelZDC[0]->ProjectionX(Form("fHistGenSelZDC%s%i",parttype[0].Data() ,i), mineebinZDC ,maxeebinZDC );
        fHistGentrueSelZDC[i-1][0] = (TH1D*)h2DGentrueSelZDC[0]->ProjectionX(Form("fHistGentrueSelZDC%s%i",parttype[0].Data() ,i),mineebinZDC , maxeebinZDC);

        fHistptSelZDC[i-1][0]  = new TH1D(Form("fHistptSelZDC_%i-%i_%s",(int)partee[i-1],(int)partee[i], parttype[0].Data()),"Signal loss corr in eff energy bins;p_{T} (GeV/c);Counts", partptbinnumb, partptbinlimits);
        fHistpttrueSelZDC[i-1][0] = new TH1D(Form("fHistpttrueSelZDC_%i-%i_%s",(int)partee[i-1],(int)partee[i], parttype[0].Data()),"Signal loss corrin eff energy bins ;p_{T} (GeV/c);Counts", partptbinnumb, partptbinlimits);
    
        //Double diff
        Double_t minmultbinDD = h3DGenSelDD[0]->GetZaxis()->FindBin( partmult[i-1]+1e-6 );
        Double_t maxmultbinDD = h3DGenSelDD[0]->GetZaxis()->FindBin( partmult[i]-1e-6 );

        fHistGenSelV0FixLowEE[i-1][0] = (TH1D*)h3DGenSelDD[0]->ProjectionX(Form("fHistGenSelDDLowEE%s%i",parttype[0].Data(),i), fixedlowmineebin, fixedlowmaxeebin, minmultbinDD, maxmultbinDD);
        fHistGenSelV0FixHighEE[i-1][0] = (TH1D*)h3DGenSelDD[0]->ProjectionX(Form("fHistGenSelDDHighEE%s%i",parttype[0].Data(),i), fixedhighmineebin, fixedhighmaxeebin, minmultbinDD, maxmultbinDD);
        fHistGentrueSelV0FixLowEE[i-1][0] = (TH1D*)h3DGentrueSelDD[0]->ProjectionX(Form("fHistGentrueSelDDLowEE%s%i",parttype[0].Data(),i), fixedlowmineebin, fixedlowmaxeebin, minmultbinDD, maxmultbinDD);
        fHistGentrueSelV0FixHighEE[i-1][0] = (TH1D*)h3DGentrueSelDD[0]->ProjectionX(Form("fHistGentrueSelDDHighEE%s%i",parttype[0].Data(),i), fixedhighmineebin, fixedhighmaxeebin, minmultbinDD, maxmultbinDD);
    
        
        fHistGenSelZDCFixLowmult[i-1][0] = (TH1D*)h3DGenSelDD[0]->ProjectionX(Form("fHistGenSelDDLowmult%s%i",parttype[0].Data(),i), mineebinZDC , maxeebinZDC, fixedlowminmultbin, fixedlowmaxmultbin);
        fHistGenSelZDCFixHighmult[i-1][0] = (TH1D*)h3DGenSelDD[0]->ProjectionX(Form("fHistGenSelDDHighmult%s%i",parttype[0].Data(),i),mineebinZDC , maxeebinZDC,  fixedhighminmultbin, fixedhighmaxmultbin);
        fHistGentrueSelZDCFixLowmult[i-1][0] = (TH1D*)h3DGentrueSelDD[0]->ProjectionX(Form("fHistGentrueSelDDLowmult%s%i",parttype[0].Data(),i),mineebinZDC , maxeebinZDC,  fixedlowminmultbin, fixedlowmaxmultbin);
        fHistGentrueSelZDCFixHighmult[i-1][0] = (TH1D*)h3DGentrueSelDD[0]->ProjectionX(Form("fHistGentrueSelDDHighmult%s%i",parttype[0].Data(),i), mineebinZDC , maxeebinZDC, fixedhighminmultbin, fixedhighmaxmultbin);
    

        fHistptSelV0FixHighEE[i-1][0]  = new TH1D(Form("fHistptSelV0FixHighEE_%i-%i_%s",(int)partmult[i-1],(int)partmult[i],parttype[0].Data()),
            Form("Signal loss corr in mult bins - ZDCperc is fixed in [%.0f,%.0f];p_{T} (GeV/c);Counts",FixedHighEE[0],FixedHighEE[1]), partptbinnumb, partptbinlimits);
        fHistptSelV0FixLowEE[i-1][0] = new TH1D(Form("fHistptSelV0FixLowEE_%i-%i_%s",(int)partmult[i-1],(int)partmult[i],parttype[0].Data()),
            Form("Signal loss corr in mult bins - ZDCperc is fixed in [%.0f,%.0f];p_{T} (GeV/c);Counts",FixedLowEE[0],FixedLowEE[1]), partptbinnumb, partptbinlimits);
     
    
        fHistptSelZDCFixHighmult[i-1][0]  = new TH1D(Form("fHistptSelZDCFixHighmult_%i-%i_%s",(int)partee[i-1],(int)partee[i],parttype[0].Data()),
            Form("Signal loss corr in eff energy bins - multperc is fixed in [%.0f,%.0f];p_{T} (GeV/c);Counts",FixedHighmult[0],FixedHighmult[1]), partptbinnumb, partptbinlimits);
        fHistptSelZDCFixLowmult[i-1][0] = new TH1D(Form("fHistptSelZDCFixLowmult_%i-%i_%s",(int)partee[i-1],(int)partee[i],parttype[0].Data()),
            Form("Signal loss corr in eff energy bins - multperc is fixed in [%.0f,%.0f];p_{T} (GeV/c);Counts",FixedLowmult[0],FixedLowmult[1]), partptbinnumb, partptbinlimits);
    
        fHistpttrueSelV0FixHighEE[i-1][0]  = new TH1D(Form("fHistpttrueSelV0FixHighEE_%i-%i_%s",(int)partmult[i-1],(int)partmult[i],parttype[0].Data()),"Signal loss corr;p_{T} (GeV/c);Counts", partptbinnumb, partptbinlimits);
        fHistpttrueSelV0FixLowEE[i-1][0] = new TH1D(Form("fHistpttrueSelV0FixLowEE_%i-%i_%s",(int)partmult[i-1],(int)partmult[i],parttype[0].Data()),"Signal loss corr;p_{T} (GeV/c);Counts", partptbinnumb, partptbinlimits);
      
        fHistpttrueSelZDCFixHighmult[i-1][0]  = new TH1D(Form("fHistpttrueSelZDCFixHighmult_%i-%i_%s",(int)partee[i-1],(int)partee[i],parttype[0].Data()),"Signal loss corr;p_{T} (GeV/c);Counts", partptbinnumb, partptbinlimits);
        fHistpttrueSelZDCFixLowmult[i-1][0] = new TH1D(Form("fHistpttrueSelZDCFixLowmult_%i-%i_%s",(int)partee[i-1],(int)partee[i],parttype[0].Data()),"Signal loss corr;p_{T} (GeV/c);Counts", partptbinnumb, partptbinlimits);
    
    } 

    //Do V0 selection
    temppt = 0;
    for (int k = 0; k<multbinnumb; k++){
	    for(long i = 1; i<fHistGenSelV0[k][0]->GetNbinsX()+1;i++){
	        temppt = fHistGenSelV0[k][0]->GetXaxis()->GetBinCenter(i);
	        for(long filling = 0; filling<fHistGenSelV0[k][0]->GetBinContent(i); filling++){
	            fHistptSelV0[k][0]->Fill(temppt);
	        }
	    }
	}

	
	temppt = 0;
	for (int k = 0; k<multbinnumb; k++){
	    for(long i = 1; i<fHistGentrueSelV0[k][0]->GetNbinsX()+1;i++){
	        temppt = fHistGentrueSelV0[k][0]->GetXaxis()->GetBinCenter(i);
	        for(long filling = 0; filling<fHistGentrueSelV0[k][0]->GetBinContent(i); filling++){
	            fHistpttrueSelV0[k][0]->Fill(temppt);
	        }
	    }
	}

	for (int k = 0; k<multbinnumb; k++){
		DivideAndComputeRogerBarlow(fHistptSelV0[k][0],fHistpttrueSelV0[k][0]);
	}


	//Do ZDC selection
    temppt = 0;
    for (int k = 0; k<eebinnumb; k++){
	    for(long i = 1; i<fHistGenSelZDC[k][0]->GetNbinsX()+1;i++){
	        temppt = fHistGenSelZDC[k][0]->GetXaxis()->GetBinCenter(i);
	        for(long filling = 0; filling<fHistGenSelZDC[k][0]->GetBinContent(i); filling++){
	            fHistptSelZDC[k][0]->Fill(temppt);
	        }
	    }
	}

	
	temppt = 0;
	for (int k = 0; k<eebinnumb; k++){
	    for(long i = 1; i<fHistGentrueSelZDC[k][0]->GetNbinsX()+1;i++){
	        temppt = fHistGentrueSelZDC[k][0]->GetXaxis()->GetBinCenter(i);
	        for(long filling = 0; filling<fHistGentrueSelZDC[k][0]->GetBinContent(i); filling++){
	            fHistpttrueSelZDC[k][0]->Fill(temppt);
	        }
	    }
	}

	for (int k = 0; k<eebinnumb; k++){
		DivideAndComputeRogerBarlow(fHistptSelZDC[k][0],fHistpttrueSelZDC[k][0]);
	}

    //Do Double Diff Selection
    //Fix ee

    temppt = 0;
    for (int k = 0; k<multbinnumb; k++){
        for(long i = 1; i<fHistGenSelV0FixHighEE[k][0]->GetNbinsX()+1;i++){
            temppt = fHistGenSelV0FixHighEE[k][0]->GetXaxis()->GetBinCenter(i);
            for(long filling = 0; filling<fHistGenSelV0FixHighEE[k][0]->GetBinContent(i); filling++){
                fHistptSelV0FixHighEE[k][0]->Fill(temppt);
            }
        }
    }

    
    temppt = 0;
    for (int k = 0; k<multbinnumb; k++){
        for(long i = 1; i<fHistGentrueSelV0FixHighEE[k][0]->GetNbinsX()+1;i++){
            temppt = fHistGentrueSelV0FixHighEE[k][0]->GetXaxis()->GetBinCenter(i);
            for(long filling = 0; filling<fHistGentrueSelV0FixHighEE[k][0]->GetBinContent(i); filling++){
                fHistpttrueSelV0FixHighEE[k][0]->Fill(temppt);
            }
        }
    }

    for (int k = 0; k<multbinnumb; k++){
        DivideAndComputeRogerBarlow(fHistptSelV0FixHighEE[k][0],fHistpttrueSelV0FixHighEE[k][0]);
    }


    temppt = 0;
    for (int k = 0; k<multbinnumb; k++){
        for(long i = 1; i<fHistGenSelV0FixLowEE[k][0]->GetNbinsX()+1;i++){
            temppt = fHistGenSelV0FixLowEE[k][0]->GetXaxis()->GetBinCenter(i);
            for(long filling = 0; filling<fHistGenSelV0FixLowEE[k][0]->GetBinContent(i); filling++){
                fHistptSelV0FixLowEE[k][0]->Fill(temppt);
            }
        }
    }

    
    temppt = 0;
    for (int k = 0; k<multbinnumb; k++){
        for(long i = 1; i<fHistGentrueSelV0FixLowEE[k][0]->GetNbinsX()+1;i++){
            temppt = fHistGentrueSelV0FixLowEE[k][0]->GetXaxis()->GetBinCenter(i);
            for(long filling = 0; filling<fHistGentrueSelV0FixLowEE[k][0]->GetBinContent(i); filling++){
                fHistpttrueSelV0FixLowEE[k][0]->Fill(temppt);
            }
        }
    }

    for(long b = 1; b<=fHistGentrueSelV0FixLowEE[0][0]->GetNbinsX(); b++){
        cout << " reco " << fHistptSelV0FixLowEE[0][0]->GetBinContent(b) << endl;
        cout <<  " MC " << fHistpttrueSelV0FixLowEE[0][0]->GetBinContent(b) << endl;
    }


    for (int k = 0; k<multbinnumb; k++){
        DivideAndComputeRogerBarlow(fHistptSelV0FixLowEE[k][0],fHistpttrueSelV0FixLowEE[k][0]);
    }
   

    //Fix mult
    temppt = 0;
    for (int k = 0; k<eebinnumb; k++){
        for(long i = 1; i<fHistGenSelZDCFixHighmult[k][0]->GetNbinsX()+1;i++){
            temppt = fHistGenSelZDCFixHighmult[k][0]->GetXaxis()->GetBinCenter(i);
            for(long filling = 0; filling<fHistGenSelZDCFixHighmult[k][0]->GetBinContent(i); filling++){
                fHistptSelZDCFixHighmult[k][0]->Fill(temppt);
            }
        }
    }

    
    temppt = 0;
    for (int k = 0; k<eebinnumb; k++){
        for(long i = 1; i<fHistGentrueSelZDCFixHighmult[k][0]->GetNbinsX()+1;i++){
            temppt = fHistGentrueSelZDCFixHighmult[k][0]->GetXaxis()->GetBinCenter(i);
            for(long filling = 0; filling<fHistGentrueSelZDCFixHighmult[k][0]->GetBinContent(i); filling++){
                fHistpttrueSelZDCFixHighmult[k][0]->Fill(temppt);
            }
        }
    }

    for (int k = 0; k<eebinnumb; k++){
        DivideAndComputeRogerBarlow(fHistptSelZDCFixHighmult[k][0],fHistpttrueSelZDCFixHighmult[k][0]);
    }


    temppt = 0;
    for (int k = 0; k<eebinnumb; k++){
        for(long i = 1; i<fHistGenSelZDCFixLowmult[k][0]->GetNbinsX()+1;i++){
            temppt = fHistGenSelZDCFixLowmult[k][0]->GetXaxis()->GetBinCenter(i);
            for(long filling = 0; filling<fHistGenSelZDCFixLowmult[k][0]->GetBinContent(i); filling++){
                fHistptSelZDCFixLowmult[k][0]->Fill(temppt);
            }
        }
    }

    
    temppt = 0;
    for (int k = 0; k<eebinnumb; k++){
        for(long i = 1; i<fHistGentrueSelZDCFixLowmult[k][0]->GetNbinsX()+1;i++){
            temppt = fHistGentrueSelZDCFixLowmult[k][0]->GetXaxis()->GetBinCenter(i);
            for(long filling = 0; filling<fHistGentrueSelZDCFixLowmult[k][0]->GetBinContent(i); filling++){
                fHistpttrueSelZDCFixLowmult[k][0]->Fill(temppt);
            }
        }
    }

    for (int k = 0; k<eebinnumb; k++){
        DivideAndComputeRogerBarlow(fHistptSelZDCFixLowmult[k][0],fHistpttrueSelZDCFixLowmult[k][0]);
    }

     cout << "I did " << "Xi" << "!" << endl;
	
   	TFile* Write = new TFile (outputname,"RECREATE");
    TDirectoryFile* epsevt = new TDirectoryFile("EventLoss","EventLoss");
    epsevt->cd();
    hevtlossV0->Write();
    hevtlossZDC->Write();
    hevtlossFixHighEE->Write();
    hevtlossFixLowEE->Write();
    hevtlossFixLowmult->Write();
    hevtlossFixHighmult->Write();
    Write->cd();
    TDirectoryFile* epspart = new TDirectoryFile("SgnLoss","SgnLoss");
    epspart->cd();
    TDirectoryFile* dirtype[npart];
    TDirectoryFile* direstim[6][npart];   

    
    dirtype[0] = new TDirectoryFile(Form("%s",parttype[0].Data()),Form("%s",parttype[0].Data()));
    dirtype[0]->cd();
    direstim[0][0] = new TDirectoryFile("multsel","multsel");
    direstim[1][0] = new TDirectoryFile("EEsel","EEsel");
    direstim[2][0] = new TDirectoryFile("EEsel_fixedlowmult","EEsel_fixedlowmult");
    direstim[3][0] = new TDirectoryFile("EEsel_fixedhighmult","EEsel_fixedhighmult");
    direstim[4][0] = new TDirectoryFile("multsel_fixedlowEE","multsel_fixedlowEE");
    direstim[5][0] = new TDirectoryFile("multsel_fixedhighEE","multsel_fixedhighEE");
    
    for (int k = 0; k<multbinnumb; k++){

        direstim[0][0]->cd();
    	fHistptSelV0[k][0]->Write();
        dirtype[0]->cd();

        direstim[4][0]->cd();
        fHistptSelV0FixLowEE[k][0]->Write();
        dirtype[0]->cd();

        direstim[5][0]->cd();
        fHistptSelV0FixHighEE[k][0]->Write();
        dirtype[0]->cd();
    }

    for (int k = 0; k<eebinnumb; k++){

        //ZDC
        direstim[1][0]->cd();
    	fHistptSelZDC[k][0]->Write();
        dirtype[0]->cd();

        direstim[2][0]->cd();
        fHistptSelZDCFixLowmult[k][0]->Write();
        dirtype[0]->cd();

        direstim[3][0]->cd();
        fHistptSelZDCFixHighmult[k][0]->Write();
        dirtype[0]->cd();

    }
    epspart->cd();


    Write->cd();
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
    //
    if (h1->GetBinContent(i) == 0 && h2->GetBinContent(i) == 0){
        h1->SetBinContent(i, 1);
        h2->SetBinContent(i, 1);      
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
