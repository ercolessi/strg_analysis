void eventloss(TString which = "hevtlossV0");
void eventlossMerged(TString which = "hevtlossV0");
void EpsEventSyst(){

    Double_t percentileV0[] = {0.,1.,5,10,15,20,30,40,50,70,100};
    const Long_t nbinV0 = sizeof(percentileV0)/sizeof(Double_t) - 1;
    //
    Double_t percentileZDC[] = {0,20,30,40,50,60,70,80,90,100};
    const Long_t nbinZDC = sizeof(percentileZDC)/sizeof(Double_t) - 1;

    eventloss();
    eventloss("hevtlossZDC");
    eventloss("hevtlossV0FixHighmult");
    eventlossMerged("hevtlossV0FixLowmult");
    eventlossMerged("hevtlossV0FixLowEE");
    eventlossMerged("hevtlossV0FixHighEE");

    }

void eventloss(TString which = "hevtlossV0") {   

        TFile* MC15g3c3 = new TFile ("EventCountLoss_15g3c3.root");
        TFile* MC16d3 = new TFile ("EventCountLoss_16d3.root");
        TFile* MC16MB = new TFile ("EventCountLoss_16d3_MB.root");
        TFile* MC15MB = new TFile ("EventCountLoss_15g3c3_MB.root");
      
        TH1F* hEevt15 = (TH1F*)MC15g3c3->Get(Form("EventLoss/%s",which.Data()));
        TH1F* hEevt16 = (TH1F*)MC16d3->Get(Form("EventLoss/%s",which.Data()));
        TH1F* hEevt16MB = (TH1F*)MC16MB->Get(Form("EventLoss/%s","hevtlossV0"));
        TH1F* hEevt15MB = (TH1F*)MC15MB->Get(Form("EventLoss/%s","hevtlossV0"));
           
        TH1F* hRatio = (TH1F*)hEevt16->Clone(Form("hEvtLoss%s",which.Data()));
        hRatio->Reset();
        TH1F* hRatio2 = (TH1F*)hEevt16->Clone(Form("hEvtLossTot%s",which.Data()));
        hRatio2->Reset();
        for (int bin = 1; bin <=hRatio->GetNbinsX(); bin ++ ){

          double a = TMath::Abs((hEevt15->GetBinContent(bin)))/hEevt16->GetBinContent(bin);
          double b = TMath::Abs((hEevt15MB->GetBinContent(1)))/hEevt16MB->GetBinContent(1);
          double dev = TMath::Abs(1-a/b);

          cout << a << endl; 
          cout << b << endl;
         
          hRatio->SetBinContent(bin, dev);          
          hRatio->SetBinError(bin, .0000000001);
          hRatio2->SetBinContent(bin, TMath::Abs( 1-a));          
          hRatio2->SetBinError(bin, .0000000001);

        }
        hRatio->SetLineWidth(2);

        TFile* Write = new TFile("EvtLossSyst.root","UPDATE");
        hRatio->Write();
         hRatio2->Write();
        Write->cd();
        Write->Close();

}
void eventlossMerged(TString which = "hevtlossV0") {   

        TFile* MC15g3c3 = new TFile ("EventCountLoss_15g3c3_mergedbins.root");
        TFile* MC16d3 = new TFile ("EventCountLoss_16d3_mergedbins.root");
        TFile* MC16MB = new TFile ("EventCountLoss_16d3_MB.root");
        TFile* MC15MB = new TFile ("EventCountLoss_15g3c3_MB.root");

        TH1F* hEevt15 = (TH1F*)MC15g3c3->Get(Form("EventLoss/%s",which.Data()));
        TH1F* hEevt16 = (TH1F*)MC16d3->Get(Form("EventLoss/%s",which.Data()));
        TH1F* hEevt16MB = (TH1F*)MC16MB->Get(Form("EventLoss/%s","hevtlossV0"));
        TH1F* hEevt15MB = (TH1F*)MC15MB->Get(Form("EventLoss/%s","hevtlossV0"));
      
        TH1F* hRatio = (TH1F*)hEevt16->Clone(Form("hEvtLoss%s",which.Data()));
        hRatio->Reset();
        TH1F* hRatio2 = (TH1F*)hEevt16->Clone(Form("hEvtLossTot%s",which.Data()));
        hRatio2->Reset();
        for (int bin = 1; bin <=hRatio->GetNbinsX(); bin ++ ){
          double a = TMath::Abs((hEevt15->GetBinContent(bin)))/hEevt16->GetBinContent(bin);
          double b = TMath::Abs((hEevt15MB->GetBinContent(1)))/hEevt16MB->GetBinContent(1);
          double dev = TMath::Abs(1-a/b);
          hRatio->SetBinContent(bin, dev);          
          hRatio->SetBinError(bin, .0000000001);
          hRatio2->SetBinContent(bin,  TMath::Abs(a-1));          
          hRatio2->SetBinError(bin, .0000000001);
        }
        hRatio->SetLineWidth(2);

        TFile* Write = new TFile("EvtLossSyst.root","UPDATE");
        hRatio->Write();
        hRatio2->Write();
        Write->cd();
        Write->Close();

}
