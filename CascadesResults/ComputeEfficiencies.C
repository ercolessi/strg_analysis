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

	TFile* file = new TFile("Final15f.root","READ");
		//"MCLHC17j_GP.root"
	const char* ZDCFilename = "ExtractedZDCPercentile_final15f.root"; //""
    const char* ZDCpartfilename = ZDCFilename; //""
	const char* outputname = "EventCountLoss_Final15f.root";
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

	//Multiplicity Efficiency
    Float_t mult[] = {0,1,5,10,15,20,30,40,50,70,100};
    Long_t multbinnumb = sizeof(mult)/sizeof(Float_t) - 1;
    Long_t lNEvtSelectedV0[multbinnumb];
    Long_t lNEvttrueSelectedV0[multbinnumb];    

    //Effective energy Efficiency
    Float_t ee[] = {0,10,20,30,40,50,60,70,80,90,100}; 
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
    TH1F* hevtlossFixLowEE = new TH1F("hevtlossV0FixedLowEE", "", multbinnumb, mult);
    hevtlossFixLowEE->SetTitle(Form("Event loss correction in mult bins - ZDCperc is fixed in [%.0f,%.0f]",FixedLowEE[0],FixedLowEE[1]));
    TH1F* hevtlossFixHighEE = new TH1F("hevtlossV0FixedHighEE", "", multbinnumb, mult);
    hevtlossFixHighEE->SetTitle(Form("Event loss correction in mult bins - ZDCperc is fixed in [%.0f,%.0f]",FixedHighEE[0],FixedHighEE[1]));   
    TH1F* hevtlossFixLowmult = new TH1F("hevtlossV0FixedLowmult", "", eebinnumb, ee);
    hevtlossFixLowmult->SetTitle(Form("Event loss correction in eff energy bins - V0Mperc is fixed in [%.0f,%.0f]",FixedLowmult[0],FixedLowmult[1]));
    TH1F* hevtlossFixHighmult = new TH1F("hevtlossV0FixedHighmult", "", eebinnumb, ee);
    hevtlossFixHighmult->SetTitle(Form("Event loss correction in eff energy bins - V0Mperc is fixed in [%.0f,%.0f]",FixedHighmult[0],FixedHighmult[1]));


    //Pt bins definition
    Double_t ptbinlimits[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
    Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;

    /*for (int k=0; k<multbinnumb; k++)
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
        lNEvtSelectedZDC[k]=0;
        lNEvttrueSelectedZDC[k]=0;           
        lNEvtSelectedFixLowEE[k]=0;
        lNEvttrueSelectedFixLowEE[k]=0;
        lNEvtSelectedFixLowmult[k]=0;
        lNEvttrueSelectedFixLowmult[k]=0;
        lNEvtSelectedFixHighEE[k]=0;
        lNEvttrueSelectedFixHighEE[k]=0;
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
	}	 */


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

        //Effective energy
        if (ee[i-1] == 0) LowLogSum = 0;
        else LowLogSum = ReverseCumulative->Eval(ee[i-1]/100+1e-6);
        HighLogSum = ReverseCumulative->Eval(ee[i]/100-1e-6);

        cout << "low " << LowLogSum << "   high " << HighLogSum << endl;

        mineebinZDC = h2DGenXitrueINELgt0ZDC->GetYaxis()->FindBin( LowLogSum );
        maxeebinZDC = h2DGenXitrueINELgt0ZDC->GetYaxis()->FindBin( HighLogSum );

        fHistGenINELgt0ZDC[i-1] = (TH1D*)h2DGenXiINELgt0ZDC->ProjectionX(Form("fHistGenINELgt0ZDC%i",i), mineebinZDC ,maxeebinZDC );
        fHistGentrueINELgt0ZDC[i-1] = (TH1D*)h2DGenXitrueINELgt0ZDC->ProjectionX(Form("fHistGentrueINELgt0ZDC%i",i),mineebinZDC , maxeebinZDC);

        fHistptINELgt0ZDC[i-1]  = new TH1D(Form("fHistptINELgt0ZDC_%i-%i",(int)ee[i-1],(int)ee[i]),"Cascade MC count;p_{T} (GeV/c);Counts", ptbinnumb, ptbinlimits);
        fHistpttrueINELgt0ZDC[i-1] = new TH1D(Form("fHistpttrueINELgt0ZDC_%i-%i",(int)ee[i-1],(int)ee[i]),"Cascade MC count;p_{T} (GeV/c);Counts", ptbinnumb, ptbinlimits);

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

	//Do ZDC selection
    temppt = 0;
    for (int k = 0; k<10; k++){
	    for(long i = 1; i<fHistGenINELgt0ZDC[k]->GetNbinsX()+1;i++){
	        temppt = fHistGenINELgt0ZDC[k]->GetXaxis()->GetBinCenter(i);
	        for(long filling = 0; filling<fHistGenINELgt0ZDC[k]->GetBinContent(i); filling++){
	            fHistptINELgt0ZDC[k]->Fill(temppt);
	        }
	    }
	}

	
	temppt2 = 0;
	for (int k = 0; k<10; k++){
	    for(long i = 1; i<fHistGentrueINELgt0ZDC[k]->GetNbinsX()+1;i++){
	        temppt2 = fHistGentrueINELgt0ZDC[k]->GetXaxis()->GetBinCenter(i);
	        for(long filling = 0; filling<fHistGentrueINELgt0ZDC[k]->GetBinContent(i); filling++){
	            fHistpttrueINELgt0ZDC[k]->Fill(temppt2);
	        }
	    }
	}

	for (int k = 0; k<10; k++){
		DivideAndComputeRogerBarlow(fHistptINELgt0ZDC[k],fHistpttrueINELgt0ZDC[k]);
	}
    
	
   	TFile* Write = new TFile (outputname,"RECREATE");
    /*hevtlossV0->Write();
    hevtlossZDC->Write();
    hevtlossFixHighEE->Write();
    hevtlossFixLowEE->Write();
    hevtlossFixLowmult->Write();
    hevtlossFixHighmult->Write();*/
    for (int k = 0; k<10; k++){
    	fHistptINELgt0V0[k]->Write();
    	fHistptINELgt0ZDC[k]->Write();
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
