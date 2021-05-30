double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 );
double PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin);

void CorrectMB(TString Type = "XiMinus", TString caso = "normal"){
//LASCIARE XIMINUS PERCHE' LE EFFCIENZE DI SGN LOSS FUNZIONANO SOLO COS^ RICORDA DI CAMBIARE
  const int nbins = 2;
	Double_t centrality[nbins] = {0.,100.};
  
  // Get corrections
  TFile* file = new TFile("NormalizationCorrections/EventCountLoss_15g3c3_MB.root","READ");
   // Get corrections
  TFile* filesgn = new TFile("NormalizationCorrections/EventCountLoss_15g3c3_MB.root","READ");
  // Which analysis -->Corrections
  //
  TString outputfilename = Form("MB15g3c3StatSpectra%s-Xi-%s_%03.0f_%03.0f.root","V0M","ZDC",0.,100.);
  //
  TString estimator = "", selection = "";
  //NON MODIFICARE QUI
  Double_t lLoEE = 0.;
  Double_t lHiEE = 100.;
  estimator = "V0";
  selection = "multsel"; 

	TH1F* hEevt = (TH1F*)file->Get(Form("EventLoss/hevtloss%s",estimator.Data()));
	TH1F* hsgnloss = (TH1F*)filesgn->Get(Form("SgnLoss/%s/%s/fHistptSel%s_%.0f-%.0f_%s",
                                      Type.Data(),selection.Data(),estimator.Data(),centrality[0],centrality[1],Type.Data()));
	
	
  TFile* filename = TFile::Open("MBSpectra.root", "READ");
	TH1D* lHistPt = (TH1D*)filename->Get("fHistPtXi");
	TH1D* hClone = (TH1D*)lHistPt->Clone(Form("XiSpectra_Stats_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",centrality[0],centrality[1], lLoEE, lHiEE));
  hClone->Reset();
         
  double sgn = 0, sgne = 0;
  for (int b = 1; b <= hClone->GetNbinsX(); b++){
    sgn = hsgnloss->GetBinContent(b+1);
  	sgne = hsgnloss->GetBinError(b+1);
     
    if (sgn!=0){
    hClone->SetBinContent(b,lHistPt->GetBinContent(b)*
                              hEevt->GetBinContent(1)/
                              sgn);
    hClone->SetBinError(b,ErrorInRatio(lHistPt->GetBinContent(b),lHistPt->GetBinError(b),
                             sgn,sgne));
    }
   
  }
    
        
    // Output File
    TFile* Write = new TFile (outputfilename, "RECREATE");
  
    hClone->Write();
    
}

//===========================================================================================
//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop - errorfrombottom) );
    }
    return 1.;
}

//----------------------------------------------------------------------------------------------------
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 ){ 
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
double PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin){
  double dev = h->GetBinContent(bin);
  double RBsigma = h->GetBinError(bin);

  if( TMath::Abs(dev/RBsigma) > nsigmas ) {return dev;}
  else {return 0.;}
}