Double_t GetEEfromZDC ( TFile* lfilename , Float_t lZPCpp, Float_t lZNCpp, Float_t lZPApp, Float_t lZNApp, Int_t lRun );



void DoPercentile(){

    //Open file and get Tree 
    TFile* Read = new TFile ("MCLHC15g3c3.root");
    TTree * T = (TTree *)Read->Get("PWGLF_StrVsMult_MC/fTreeEvent");
    Float_t ZPCpp, ZNCpp, ZPApp, ZNApp;
    Int_t fRun;
  
    //Set branches
    T->SetBranchAddress("fZPCpp",&ZPCpp);
    T->SetBranchAddress("fZPApp",&ZPApp);
    T->SetBranchAddress("fZNCpp",&ZNCpp);
    T->SetBranchAddress("fZNApp",&ZNApp);
    T->SetBranchAddress("fRun",&fRun);
    TH1D* h = new TH1D("h","h",100,0.,100.);
    TFile* file = new TFile ("ExtractedZDCPercentile_MC-15g3c3.root");
    //Get runList
    for(Int_t i=0; i<T->GetEntries();i++)
    {
      T->GetEvent(i);//get Tree entries
      double zdc = GetEEfromZDC(file,ZPCpp,ZNCpp,ZPApp,ZNApp,fRun);
      h->Fill(zdc);



      
    }
    new TCanvas;
    h->Draw();
}

Double_t GetEEfromZDC ( TFile* lfilename , Float_t lZPCpp, Float_t lZNCpp, Float_t lZPApp, Float_t lZNApp, Int_t lRun ){ //set in the run func
	//Converts ZDC information to effective energy percentiles. 
    
    TH1F * fhCumulative = (TH1F *)lfilename->Get(Form("%i/hCumulative_%i",lRun, lRun));
  
	Double_t Sum = TMath::Log10(TMath::Abs(lZPCpp+lZNCpp+lZPApp+lZNApp)+1);
	Double_t fZDCCentrality;

    if (Sum!=0) {
    	fZDCCentrality = 100*(fhCumulative->Interpolate(Sum));
    }
    else fZDCCentrality = 0;	

    return fZDCCentrality;
}