{
    TFile *filePart = new TFile(Form("ResultsFiles/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", "XiMinus", 0.,100.,0.,100.));
    TH1F* SpectraPart = (TH1F *) filePart->Get(Form("fHistPt%s", "XiMinus"));
    TFile* fileAntiPart = new TFile(Form("ResultsFiles/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root","XiPlus", 0.,100.,0.,100.));
    TH1F* SpectraAntiPart = (TH1F *) fileAntiPart->Get(Form("fHistPt%s", "XiPlus"));
        
    //Sum
    TH1F* Spectra = (TH1F*)SpectraPart->Clone("fHistPtXi");
    Spectra->Reset();
    for (int bin = 1 ; bin <=  SpectraPart->GetNbinsX(); bin ++ ){
        Spectra->SetBinContent(bin, SpectraPart->GetBinContent(bin) + SpectraAntiPart->GetBinContent(bin));
        Spectra->SetBinError(bin, TMath::Sqrt(SpectraPart->GetBinError(bin)*SpectraPart->GetBinError(bin)+SpectraAntiPart->GetBinError(bin)*SpectraAntiPart->GetBinError(bin)));
    }
    
     // Output File
    TFile* Write = new TFile ("MBSpectra.root", "RECREATE");  
    Spectra->Write();
}