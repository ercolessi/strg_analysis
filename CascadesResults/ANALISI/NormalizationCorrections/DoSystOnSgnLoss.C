void signalloss(TString fWhichEstimator = "V0M", TString part = "XiMinus", TString dir = "multsel", TString histo = "V0");
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 );
double PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin);

void DoSystOnSgnLoss(){

    Double_t percentileV0[] = {0.,1.,5,10,15,20,30,40,50,70,100};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    //
    Double_t percentileZDC[] = {0,20,30,40,50,60,70,80,90,100};
    const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;

    signalloss("ZDC","XiMinus","EEsel","ZDC");
   // signalloss("ZDC","XiMinus","EEsel_fixedlowmult","ZDCFixLowmult");
    signalloss("ZDC","XiMinus","EEsel_fixedhighmult","ZDCFixHighmult");
   // signalloss("V0M","XiMinus","multsel_fixedlowEE","V0FixLowEE");
   // signalloss("V0M","XiMinus","multsel_fixedhighEE","V0FixHighEE");
    signalloss("V0M","XiMinus","multsel","V0");

    }

void signalloss(TString fWhichEstimator = "V0M", TString part = "XiMinus", TString dir = "multsel", TString histo = "V0") {   

        TFile* MC15g3c3 = new TFile ("EventCountLoss_15g3c3.root");
        TFile* MC16d3 = new TFile ("EventCountLoss_16d3.root");
      
        Int_t percentileV0[] = {0,1,5,10,15,20,30,40,50,70,100};
        const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Int_t) - 1;
        //
        Int_t percentileZDC[] = {0,20,30,40,50,60,70,80,90,100};
        const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Int_t) - 1;
        //
        Int_t * percentile, *systpercentile;
        Int_t *systcounter;
        //
        // Choose scenario and initialize correct variables
        int tempbinnumber = 0;
        if (fWhichEstimator.Contains("V0M")) {
            tempbinnumber = nbinV0;
            percentile = percentileV0;
        } 
        else if (fWhichEstimator.Contains("ZDC")) {
            tempbinnumber = nbinZDC;
            percentile = percentileZDC;
        } 
        else {cout << "No valid name for estimator... its V0M or ZDC" << endl;}
        const int binnumber = tempbinnumber;

        TH1F* hsgnloss15[binnumber],* hsgnloss16[binnumber];
        TH1F* hRatio[binnumber],* hRatioClone[binnumber];
        for (int i =0; i<binnumber;i++){
            hsgnloss15[i] = (TH1F*)MC15g3c3->Get(Form("SgnLoss/%s/%s/fHistptSel%s_%i-%i_%s",part.Data(),dir.Data(), histo.Data(),percentile[i],percentile[i+1],part.Data()));
            hsgnloss16[i] = (TH1F*)MC16d3->Get(Form("SgnLoss/%s/%s/fHistptSel%s_%i-%i_%s",part.Data(),dir.Data(), histo.Data(),percentile[i],percentile[i+1],part.Data()));
       
            hRatio[i] = (TH1F*)hsgnloss15[i]->Clone(Form("hRatio%i",i));
            hRatioClone[i] = (TH1F*)hsgnloss16[i]->Clone(Form("hRatioClone%i",i));
            hRatioClone[i]->Reset();

            hRatioClone[i]->GetYaxis()->SetTitle("#frac{ | #varepsilon_{part}^{EPOS-LHC} - #varepsilon_{part}^{Pythia62011} |}{#varepsilon_{part}^{EPOS-LHC}}");
            DivideAndComputeRogerBarlow(hRatio[i],hsgnloss16[i]);
            for (int bin = 1; bin <=hRatio[i]->GetNbinsX(); bin ++ ){
              double dev = PassRogerBarlowCriterion(1,hRatio[i],bin);
              hRatioClone[i]->SetBinContent(bin, dev);          
              hRatioClone[i]->SetBinError(bin, .0000000001);
            }
            hRatioClone[i]->SetLineWidth(2);
        }  

        hRatioClone[0]->SetLineColor(kRed+1);
        hRatioClone[1]->SetLineColor(kRed-4);
        hRatioClone[2]->SetLineColor(kOrange+7);
        hRatioClone[3]->SetLineColor(kOrange-3);
        hRatioClone[4]->SetLineColor(kYellow+1);
        hRatioClone[5]->SetLineColor(kSpring-7);
        hRatioClone[6]->SetLineColor(kGreen+2);
        hRatioClone[7]->SetLineColor(kAzure+8);
        hRatioClone[8]->SetLineColor(kBlue-4);
        if (fWhichEstimator.Contains("V0M")) hRatioClone[9]->SetLineColor(kBlue+3);

        TFile* Write = new TFile("SgnLossSyst.root","UPDATE");
        TDirectoryFile *lDirCutVar = new TDirectoryFile(Form("%s-%s",dir.Data(),histo.Data()),Form("%s-%s",dir.Data(),histo.Data()));
        lDirCutVar->cd();
        for (int i =0; i<binnumber;i++){
            hRatioClone[i]->Write();
        }
        Write->cd();
        Write->Close();

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

    if( TMath::Abs(dev-1) > nsigmas*RBsigma ) {return dev;}
    else {return 1.;}
 
}
