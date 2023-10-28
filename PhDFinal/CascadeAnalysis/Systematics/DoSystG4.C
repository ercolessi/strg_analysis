void DoSystG4(){
    TF1* f = new TF1("f","[0]*TMath::Power(x,[1])",0,10); //0.00327*TMath::Power(pt,0.19716);
    f->SetParameter(0,0.00327);
    f->SetParameter(1,0.19716);

    //f->Draw();

    Double_t ptLambda[] = {0.4,  0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.5, 2.9, 3.4, 4, 5, 6.5, 8};
    Long_t nLambda = sizeof(ptLambda)/sizeof(Double_t) - 1;
	Double_t ptXi[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
    Long_t nXi = sizeof(ptXi)/sizeof(Double_t) - 1;
	Double_t ptOmega[] = {0.90, 1.60, 2.20, 2.60, 3.00, 3.80, 5.50 }; 
    Long_t nOmega = sizeof(ptOmega)/sizeof(Double_t) - 1;

    TH1D* hXi = new TH1D("hXi", "Syst. uncertainty on #Xi^{+};p_{T};Rel. Syst. Unc.",nXi,ptXi);
    for(int b = 1; b<=hXi->GetNbinsX(); b++){
        hXi->SetBinContent(b,f->Eval((ptXi[b-1]+ptXi[b])/2));
    }

    TCanvas* c = new TCanvas();
    c->SetRightMargin(0.09);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    hXi->SetStats(0);
    hXi->Draw();
    f->Draw("SAME");

    TFile* fpXi = TFile::Open("../CascadeAnalysis/results/Results-XiMinus-13TeV-SPDClusters_000_100_V0M_000_100.root");
    TFile* fapXi = TFile::Open("../CascadeAnalysis/results/Results-XiPlus-13TeV-SPDClusters_000_100_V0M_000_100.root");
    TH1D* hpXi = (TH1D *) fpXi->Get(Form("fHistPt%s", "XiMinus"));
    TH1D* hapXi = (TH1D *) fapXi->Get(Form("fHistPt%s", "XiPlus"));
    TH1D* htotXi = (TH1D*)hpXi->Clone("HistPtXi");
    htotXi->Reset();
   
    for (int bin = 1; bin <=  htotXi->GetNbinsX(); bin ++){
      htotXi->SetBinContent(bin, hpXi->GetBinContent(bin)+hapXi->GetBinContent(bin));
      htotXi->SetBinError(bin, TMath::Sqrt(hpXi->GetBinError(bin)*hpXi->GetBinError(bin) +
                                               hapXi->GetBinError(bin)*hapXi->GetBinError(bin)));       
    }

    TH1D* hXiTot = new TH1D("hXiTot", "Syst. uncertainty on #Xi^{-}+#Xi^{+};p_{T};Rel. Syst. Unc.",nXi,ptXi);
    for(int b = 1; b<=hXiTot->GetNbinsX(); b++){
        double sigx = f->Eval((ptXi[b-1]+ptXi[b])/2)*hapXi->GetBinContent(b);
        double reltot = sigx/htotXi->GetBinContent(b);
        hXiTot->SetBinContent(b,reltot);
    }
    
    TCanvas* ct = new TCanvas();
    ct->SetRightMargin(0.09);
    ct->SetLeftMargin(0.15);
    ct->SetBottomMargin(0.15);
    hXiTot->SetStats(0);
    hXiTot->Draw();
    

    
    TH1D* hLambda = new TH1D("hLambda", "Syst. uncertainty on #bar{#Lambda};p_{T};Rel. Syst. Unc.",nLambda,ptLambda);
    for(int b = 1; b<=hLambda->GetNbinsX(); b++){
        hLambda->SetBinContent(b,f->Eval((ptLambda[b-1]+ptLambda[b])/2));
    }

    TCanvas* c1 = new TCanvas();
    c1->SetRightMargin(0.09);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.15);
    hLambda->SetStats(0);
    hLambda->Draw();
    f->Draw("SAME");

    TFile* fpLambda = TFile::Open("../V0Analysis/results/Results-Lambda-13TeV-SPDClusters_000_100_V0M_000_100_FDUseMCRatio.root");
    TFile* fapLambda = TFile::Open("../V0Analysis/results/Results-AntiLambda-13TeV-SPDClusters_000_100_V0M_000_100_FDUseMCRatio.root");
    TH1D* hpLambda = (TH1D *) fpLambda->Get(Form("fHistPt%s", "Lambda"));
    TH1D* hapLambda = (TH1D *) fapLambda->Get(Form("fHistPt%s", "AntiLambda"));
    TH1D* htotLambda = (TH1D*)hpLambda->Clone("HistPtLambda");
    htotLambda->Reset();
   
    for (int bin = 1; bin <=  htotLambda->GetNbinsX(); bin ++){
      htotLambda->SetBinContent(bin, hpLambda->GetBinContent(bin)+hapLambda->GetBinContent(bin));
      htotLambda->SetBinError(bin, TMath::Sqrt(hpLambda->GetBinError(bin)*hpLambda->GetBinError(bin) +
                                               hapLambda->GetBinError(bin)*hapLambda->GetBinError(bin)));       
    }

    TH1D* hLambdaTot = new TH1D("hLambdaTot", "Syst. uncertainty on #Lambda+#bar{#Lambda};p_{T};Rel. Syst. Unc.",nLambda,ptLambda);
    for(int b = 1; b<=hLambdaTot->GetNbinsX(); b++){
        double sigx = f->Eval((ptLambda[b-1]+ptLambda[b])/2)*hapLambda->GetBinContent(b);
        double reltot = sigx/htotLambda->GetBinContent(b);
        hLambdaTot->SetBinContent(b,reltot);
    }
    
    TCanvas* c1t = new TCanvas();
    c1t->SetRightMargin(0.09);
    c1t->SetLeftMargin(0.15);
    c1t->SetBottomMargin(0.15);
    hLambdaTot->SetStats(0);
    hLambdaTot->Draw();


    TH1D* hOmega = new TH1D("hOmega", "Syst. uncertainty on #Omega^{+};p_{T};Rel. Syst. Unc.",nOmega,ptOmega);
    for(int b = 1; b<=hOmega->GetNbinsX(); b++){
        hOmega->SetBinContent(b,f->Eval((ptOmega[b-1]+ptOmega[b])/2));
    }

    TCanvas* c2 = new TCanvas();
    c2->SetRightMargin(0.09);
    c2->SetLeftMargin(0.15);
    c2->SetBottomMargin(0.15);
    hOmega->SetStats(0);
    hOmega->Draw();
    f->Draw("SAME");

    TFile* fpOmega = TFile::Open("../CascadeAnalysis/results/Results-OmegaMinus-13TeV-SPDClusters_000_100_V0M_000_100.root");
    TFile* fapOmega = TFile::Open("../CascadeAnalysis/results/Results-OmegaPlus-13TeV-SPDClusters_000_100_V0M_000_100.root");
    TH1D* hpOmega = (TH1D *) fpOmega->Get(Form("fHistPt%s", "OmegaMinus"));
    TH1D* hapOmega = (TH1D *) fapOmega->Get(Form("fHistPt%s", "OmegaPlus"));
    TH1D* htotOmega = (TH1D*)hpOmega->Clone("HistPtOmega");
    htotOmega->Reset();
   
    for (int bin = 1; bin <=  htotOmega->GetNbinsX(); bin ++){
      htotOmega->SetBinContent(bin, hpOmega->GetBinContent(bin)+hapOmega->GetBinContent(bin));
      htotOmega->SetBinError(bin, TMath::Sqrt(hpOmega->GetBinError(bin)*hpOmega->GetBinError(bin) +
                                               hapOmega->GetBinError(bin)*hapOmega->GetBinError(bin)));       
    }

    TH1D* hOmegaTot = new TH1D("hOmegaTot", "Syst. uncertainty on #Omega^{-}+#Omega^{+};p_{T};Rel. Syst. Unc.",nOmega,ptOmega);
    for(int b = 1; b<=hOmegaTot->GetNbinsX(); b++){
        double sigx = f->Eval((ptOmega[b-1]+ptOmega[b])/2)*hapOmega->GetBinContent(b);
        double reltot = sigx/htotOmega->GetBinContent(b);
        hOmegaTot->SetBinContent(b,reltot);
    }
    
    TCanvas* c2t = new TCanvas();
    c2t->SetRightMargin(0.09);
    c2t->SetLeftMargin(0.15);
    c2t->SetBottomMargin(0.15);
    hOmegaTot->SetStats(0);
    hOmegaTot->Draw();


    c->SaveAs("images/SystG4XiPlus.png");
    c1->SaveAs("images/SystG4AntiLambda.png");
    c2->SaveAs("images/SystG4OmegaPlus.png");
    ct->SaveAs("images/SystG4Xi.png");
    c1t->SaveAs("images/SystG4Lambda.png");
    c2t->SaveAs("images/SystG4Omega.png");


    TFile* Write = new TFile("G4Syst.root", "RECREATE");
    hXi->Write();
    hOmega->Write();
    hLambda->Write();
    hXiTot->Write();
    hOmegaTot->Write();
    hLambdaTot->Write();
    

}