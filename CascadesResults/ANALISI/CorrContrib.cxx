{
     TFile* fileLowMult  = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","ZDC",70.,100.,0.,100.),"READ");
    TFile* fileHighMult = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","ZDC",0.,30.,0.,100.),"READ");
   //
    TFile* fileLowEE  = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","V0M",0.,100.,70.,100.),"READ");
    TFile* fileHighEE = TFile::Open(Form("TestfullfitYield-%s-%s-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","Xi","V0M",0.,100.,0.,30.),"READ");
   //
  //  TFile* fileV0  = TFile::Open("TestfullfitYield-Xi-V0M-V0M_000_100_ZDC_000_100.root","READ");
    TFile* fileZDC = TFile::Open("TestfullfitYield-Xi-ZDC-V0M_000_100_ZDC_000_100.root","READ");
   //
    
    TH1D* h1 = (TH1D*)fileLowMult->Get("hTot_Uncorr");
    TH1D* h2 = (TH1D*)fileHighMult->Get("hTot_Uncorr");
    TH1D* h3 = (TH1D*)fileLowEE->Get("hTot_Uncorr");
    TH1D* h4 = (TH1D*)fileHighEE->Get("hTot_Uncorr");
    TH1D* h5 = (TH1D*)fileZDC->Get("hTot_Uncorr");

    h1->Add(h2);
     h1->Add(h3);
      h1->Add(h4);
       h1->Add(h5);
      
      new TCanvas;
      h1->GetXaxis()->SetTitle("Syst_{TOT} - Syst_{UNCORR}");
      h1->GetYaxis()->SetTitle("Counts");
      h1->Draw();
      h1->Fit("gaus");
}