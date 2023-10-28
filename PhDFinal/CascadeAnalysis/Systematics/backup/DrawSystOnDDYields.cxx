void DrawSystOnDDYields(
    TString fWhichParticle = "Xi",
    TString fWhichFixedEstimator = "V0M",
    TString fWhichVarEstimator = "SPDClusters",
    Double_t lFixedLo = 40., 
    Double_t lFixedHi = 50.){

    TFile* fCutVar = TFile::Open(Form("FinalSystematicsYield-%s-fixed%s_%03.0f_%03.0f.root",fWhichParticle.Data(),fWhichFixedEstimator.Data(),lFixedLo,lFixedHi));
    TFile* fExtrap = TFile::Open(Form("Final-Extrap-Syst-%s-%s_%03.0f_%03.0f.root",fWhichParticle.Data(),fWhichFixedEstimator.Data(),lFixedLo, lFixedHi));

    TH1F* hCutVar = (TH1F*) fCutVar->Get("hSystTot");
    TH1F* hExtrap = (TH1F*) fExtrap->Get("hSystUncorr");
    hCutVar->SetLineColor(kBlue);
    TH1F* htot = (TH1F*)hCutVar->Clone("htot");
    TH1F* hExtrapClone = (TH1F*)hCutVar->Clone("hExtrapClone");
    hExtrapClone->Reset();
    htot->Reset();
    for(int i = 1; i<=hCutVar->GetNbinsX();i++){
        htot->SetBinContent(i,TMath::Sqrt(hExtrap->GetBinContent(i)*hExtrap->GetBinContent(i) + hCutVar->GetBinContent(i)*hCutVar->GetBinContent(i)));
        hExtrapClone->SetBinContent(i,hExtrap->GetBinContent(i));
    }
    hExtrapClone->SetLineColor(kRed);
    htot->SetLineColor(kBlack);
    htot->GetYaxis()->SetTitle("Relative Syst. Uncertainty");
    htot->GetXaxis()->SetTitle(Form("%s percentile",fWhichVarEstimator.Data()));

    TCanvas * c = new TCanvas("c", "", 1000,800);
    c->SetFillColor(kWhite);
    c->SetLeftMargin(0.17);
    c->SetRightMargin(0.17);
    c->SetBottomMargin(0.17);
    c->SetGridy();
    c->SetGridx();
    htot->GetXaxis()->SetTitleOffset(1.5);
    htot->GetXaxis()->SetTitleSize(.04);
    htot->GetYaxis()->SetTitleOffset(1.5);
    htot->GetYaxis()->SetTitleSize(.04);
    htot->SetLineWidth(3);
    htot->GetYaxis()->SetRangeUser(0.,0.2);
    htot->Draw();
    hCutVar->Draw("SAME");
    hExtrapClone->Draw("SAME HIST");

    TLatex *xlabel = new TLatex();
    xlabel->SetTextFont(42);
    xlabel-> SetNDC();
    xlabel-> SetTextColor(1);
    xlabel-> SetTextSize(0.08);
    xlabel-> SetTextAlign(22);
    xlabel-> SetTextAngle(0);

    xlabel-> DrawLatex(0.65, 0.75,Form("#%s",fWhichParticle.Data()));

    TLatex *xlabel2 = new TLatex();
   	//xlabel2->SetTextFont(42);
    xlabel2-> SetNDC();
    xlabel2-> SetTextColor(1);
    xlabel2-> SetTextSize(0.03);
    xlabel2-> SetTextAlign(22);
    xlabel2-> SetTextAngle(0);
    xlabel2-> DrawLatex(0.65, 0.85, Form("%s fixed [%.0f-%.0f]",fWhichFixedEstimator.Data(),lFixedLo,lFixedHi));

    TLegend* l = new TLegend(0.2,0.6,0.5,0.89);
    l->SetTextSize(0.03);
    l->SetBorderSize(0);     
    l->AddEntry(htot, "Total","L");
    l->AddEntry(hCutVar, "Cut variation","L");
    l->AddEntry(hExtrapClone, "Extrapolation","L");

    l->Draw("SAME");

    c->SaveAs(Form("images/DD%sFinalContribution_%s_%03.0f_%03.0f.png",fWhichParticle.Data(),fWhichVarEstimator.Data(),lFixedLo,lFixedHi));
    
}