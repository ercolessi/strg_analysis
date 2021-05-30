double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );

void FaiRapportoFinali(){

    TString Type = "XiMinus";

    //Spettri THIS
    TFile* MyfileMB = new TFile("StatSpectraV0M-Xi-ZDC_000_100.root","READ");
    //Spettri Fiorella
    TFile* FiorfileMB = new TFile("../CLEAN/SpectraVsMultiplicityXi.root","READ");
  
    Float_t mult[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
	Long_t multbinnumb = sizeof(mult) / sizeof(Float_t) - 1;

    TH1D* mHistPtMB[multbinnumb];
    TH1D* fHistPtMB[multbinnumb];
    TH1D* HistPtCloneMB[multbinnumb];
    for(int nmult = 0; nmult < multbinnumb; nmult++)
  	{
        mHistPtMB[nmult] = (TH1D*)MyfileMB->Get(Form("XiSpectra_Stats_V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f",mult[nmult],mult[nmult+1], 0., 100.));
        fHistPtMB[nmult] = (TH1D*)FiorfileMB->Get(Form("hPtXiStatOnly_V0M_%03.0f00to%03.0f00-epsPart-epsEv-Corrected",mult[nmult],mult[nmult+1]));
        //cout << Form("hPtXiStatOnly_V0M_%03.0f00to%03.0f00-epsPart-epsEv-Corrected",mult[nmult],mult[nmult+1]) << endl;
        HistPtCloneMB[nmult] = (TH1D*)mHistPtMB[nmult]->Clone(Form("fHistPtMB%03.0f_%03.0f",mult[nmult],mult[nmult+1]));  
        HistPtCloneMB[nmult]->Reset();
    

        //Fill them
        for (int i=1; i<=fHistPtMB[0]->GetNbinsX(); i++){    
            if (fHistPtMB[nmult]->GetBinContent(i)!=0) {
                HistPtCloneMB[nmult]->SetBinContent(i, mHistPtMB[nmult]->GetBinContent(i)/fHistPtMB[nmult]->GetBinContent(i));
                HistPtCloneMB[nmult]->SetBinError(i,ErrorInRatio(mHistPtMB[nmult]->GetBinContent(i),mHistPtMB[nmult]->GetBinError(i), 
                                                                fHistPtMB[nmult]->GetBinContent(i),fHistPtMB[nmult]->GetBinError(i)));
            }
            else {
                continue;
            }
        }
    }

    //Systematics
    //Total syst canvas shaded 
    TFile* SystFile = TFile::Open("sistematiche/SystematicsFinalResults-Xi-13TeV-V0M_000_100_ZDC_000_100.root"); 
    TH1F* hSystTot = (TH1F*) SystFile->Get("hSystTot");
    TH1F* hSystTotNeg = (TH1F*)hSystTot->Clone("hSystTotNeg");
    hSystTotNeg->Reset();
    for (int i = 1;i <= hSystTot->GetNbinsX();i++){
        hSystTotNeg->SetBinError(i,hSystTot->GetBinContent(i));
        hSystTotNeg->SetBinContent(i, 1.);
    }
    hSystTotNeg->SetFillColor( 18 );
    hSystTotNeg->SetMarkerStyle(0);
 
    
    TCanvas* c1 = new TCanvas("c1","",1500,1000);
    //Corrected spectra
    c1->SetGridy();
    c1->SetRightMargin(0.09);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.15);

    TLegend* l = new TLegend(0.6,0.55,0.89,0.89);
    l->SetTextSize(0.020);
    l->SetBorderSize(0);
    HistPtCloneMB[0]->SetStats(0);
    HistPtCloneMB[0]->GetYaxis()->SetTitleSize(0.04);
    HistPtCloneMB[0]->GetXaxis()->SetTitleSize(0.04);
    HistPtCloneMB[0]->GetYaxis()->SetTitleOffset(1.2);
    HistPtCloneMB[0]->SetTitle(" ");
    HistPtCloneMB[0]->GetYaxis()->SetTitle("CorrSpectra_{THIS} / CorrSpectra_{published_pp@13TeV}");
    HistPtCloneMB[0]->GetYaxis()->SetRangeUser(0.4,1.8);
    HistPtCloneMB[0]->GetXaxis()->SetRangeUser(.8,6.5);
    HistPtCloneMB[0]->SetMarkerStyle(8);
    HistPtCloneMB[0]->SetMarkerColor(2);
    HistPtCloneMB[0]->SetMarkerSize(1.8);
    HistPtCloneMB[0]->Draw();
    hSystTotNeg->Draw("E2 SAME");
    HistPtCloneMB[0]->Draw("SAME");
    l->SetHeader("V0 selections (percentiles):");
    l->AddEntry(HistPtCloneMB[0],"0-1 ","LEP");
    for(int nmult = 1; nmult < multbinnumb; nmult++){
        HistPtCloneMB[nmult]->SetMarkerSize(1.8);
        HistPtCloneMB[nmult]->Draw("SAME");
        HistPtCloneMB[nmult]->SetTitle(Form("%.0f-%.0f ",mult[nmult],mult[nmult+1]));
        l->AddEntry(HistPtCloneMB[nmult],"","LEP");
    }
    l->AddEntry(hSystTotNeg,"syst integrated mult. case (centered in 1)","F");
    l->Draw("SAME");

    c1->Draw();
    c1->SaveAs("images/RapportiSpectra.png");
    
}

//---------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop + errorfrombottom) );
    }
    return 1.;
}