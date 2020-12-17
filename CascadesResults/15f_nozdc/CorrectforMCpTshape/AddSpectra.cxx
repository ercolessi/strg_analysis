
//---------------------------------------------------------------------------------------------------
double ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr);

void AddSpectra(TString fWhichParticle = "XiMinus", TString fWhichAntiParticle = "XiPlus"){

	// Def multiplicity
	Float_t mult[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
	Long_t multbinnumb = sizeof(mult) / sizeof(Float_t) - 1;
	Float_t multOmega[] = {0, 5, 15, 30, 50, 100};
	Long_t multbinnumbOmega = sizeof(mult) / sizeof(Float_t) - 1;
	Double_t ptbinlimits[] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};//, 8.0, 10.0 };
	Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;


	TFile* filePart[multbinnumb];
  	TFile* fileAntiPart[multbinnumb];
	TH1F* SpectraPart[multbinnumb];
	TH1F* SpectraAntiPart[multbinnumb];

	TFile* file = new TFile("~/EventCountLoss_15f.root","READ");
	TH1F* hEevt = (TH1F*)file->Get("EventLoss/hevtlossV0");
	TH1F* hsgnlossPart[multbinnumb];
	for (int i =0; i<multbinnumb;i++){
	    hsgnlossPart[i] = (TH1F*)file->Get(Form("SgnLoss/%s/%s/fHistptSel%s_%.0f-%.0f_%s",fWhichParticle.Data(),"multsel", "V0",mult[i],mult[i+1],fWhichParticle.Data()));
	}
	TH1F* hsgnlossAntiPart[multbinnumb];
	for (int i =0; i<multbinnumb;i++){
	    hsgnlossAntiPart[i] = (TH1F*)file->Get(Form("SgnLoss/%s/%s/fHistptSel%s_%.0f-%.0f_%s",fWhichAntiParticle.Data(),"multsel", "V0",mult[i],mult[i+1],fWhichAntiParticle.Data()));
	}
  	
  	for(int nmult = 0; nmult < multbinnumb; nmult++)
  	{
    	filePart[nmult] = new TFile(Form("IT3/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichParticle.Data(), mult[nmult], mult[nmult+1], 0., 100.));
    	fileAntiPart[nmult] = new TFile(Form("IT3/Results-%s-13TeV-V0M_%03.0f_%03.0f_ZDC_%03.0f_%03.0f.root", fWhichAntiParticle.Data(), mult[nmult], mult[nmult+1], 0., 100.));
    	
	    SpectraPart[nmult] = (TH1F *) filePart[nmult]->Get(Form("fHistPt%s", fWhichParticle.Data()));
	    SpectraAntiPart[nmult] = (TH1F *) fileAntiPart[nmult]->Get(Form("fHistPt%s", fWhichAntiParticle.Data()));

	    for (int bin = 1; bin <= SpectraPart[nmult]->GetNbinsX(); bin ++){
    		SpectraPart[nmult]->SetBinContent(bin,SpectraPart[nmult]->GetBinContent(bin)/hsgnlossPart[nmult]->GetBinContent(bin)*hEevt->GetBinContent(nmult+1));
    		SpectraPart[nmult]->SetBinError(bin,hEevt->GetBinContent(nmult+1)*(ErrorInRatio(SpectraPart[nmult]->GetBinContent(bin),SpectraPart[nmult]->GetBinError(bin), 
    			hsgnlossPart[nmult]->GetBinContent(bin), hsgnlossPart[nmult]->GetBinError(bin))));
    		SpectraAntiPart[nmult]->SetBinContent(bin,SpectraAntiPart[nmult]->GetBinContent(bin)/hsgnlossAntiPart[nmult]->GetBinContent(bin)*hEevt->GetBinContent(nmult+1));
    		SpectraAntiPart[nmult]->SetBinError(bin,hEevt->GetBinContent(nmult+1)*(ErrorInRatio(SpectraAntiPart[nmult]->GetBinContent(bin),SpectraAntiPart[nmult]->GetBinError(bin), 
    			hsgnlossAntiPart[nmult]->GetBinContent(bin), hsgnlossAntiPart[nmult]->GetBinError(bin))));
		}

		for (int jbin = 1; jbin <= SpectraPart[nmult]->GetNbinsX(); jbin ++){
			SpectraPart[nmult]->SetBinContent(jbin,SpectraPart[nmult]->GetBinContent(jbin)+SpectraAntiPart[nmult]->GetBinContent(jbin));
			SpectraPart[nmult]->SetBinError(jbin,TMath::Sqrt(SpectraPart[nmult]->GetBinError(jbin)*SpectraPart[nmult]->GetBinError(jbin)+SpectraAntiPart[nmult]->GetBinError(jbin)*SpectraAntiPart[nmult]->GetBinError(jbin)));
		}

	}

		//Set colors
	  SpectraPart[0]->SetLineColor(kRed + 1);
	  SpectraPart[1]->SetLineColor(kRed - 4);
	  SpectraPart[2]->SetLineColor(kOrange + 7);
	  SpectraPart[3]->SetLineColor(kOrange - 3);
	  SpectraPart[4]->SetLineColor(kYellow + 1);
	  SpectraPart[5]->SetLineColor(kSpring - 7);
	  SpectraPart[6]->SetLineColor(kGreen + 2);
	  SpectraPart[7]->SetLineColor(kAzure + 8);
	  SpectraPart[8]->SetLineColor(kBlue - 4);
	  SpectraPart[9]->SetLineColor(kBlue + 3);

	  SpectraPart[0]->SetMarkerColor(kRed + 1);
	  SpectraPart[1]->SetMarkerColor(kRed - 4);
	  SpectraPart[2]->SetMarkerColor(kOrange + 7);
	  SpectraPart[3]->SetMarkerColor(kOrange - 3);
	  SpectraPart[4]->SetMarkerColor(kYellow + 1);
	  SpectraPart[5]->SetMarkerColor(kSpring - 7);
	  SpectraPart[6]->SetMarkerColor(kGreen + 2);
	  SpectraPart[7]->SetMarkerColor(kAzure + 8);
	  SpectraPart[8]->SetMarkerColor(kBlue - 4);
	  SpectraPart[9]->SetMarkerColor(kBlue + 3);

	  TCanvas* c = new TCanvas();
	  c->SetLogy();
	  SpectraPart[0]->GetYaxis()->SetRangeUser(1E-6,1);
	  SpectraPart[0]->GetYaxis()->SetTitle("");
	  SpectraPart[0]->GetYaxis()->SetTitleOffset(1.2);
	  SpectraPart[0]->SetStats(0);
	  for(int nmult = 0; nmult < multbinnumb; nmult++){
	  	SpectraPart[nmult]->SetMarkerStyle(8);    
	  	SpectraPart[nmult]->Draw("SAME");
	  	SpectraPart[nmult]->SetName(Form("XiSpectra_Stats_%.0f-%.0f",mult[nmult],mult[nmult+1]));
	  }

	  TFile* lResultsFile = TFile::Open("XiSpectra.root", "RECREATE");
	  TH1F* hClone[multbinnumb];
	  for(int nmult = 0; nmult < multbinnumb; nmult++){
	  	SpectraPart[nmult]->Write();
	  	hClone[nmult] = (TH1F*)SpectraPart[nmult]->Clone(Form("hClone%i",nmult));  
	  	hClone[nmult]->Reset();

	  }


	  TFile* fileFior = new TFile("~/Scaricati/SpectraVsMultiplicityXi.root", "READ");
	  TH1F* SpectraFior[multbinnumb];
	  for(int nmult = 0; nmult < multbinnumb; nmult++){
	  	SpectraFior[nmult] = (TH1F *)fileFior->Get(Form("hPtXiStatOnly_V0M_%03.0f00to%03.0f00-epsPart-epsEv-Corrected", mult[nmult],mult[nmult+1])); 	 
	  
	  	for (int bin = 1; bin <= SpectraPart[nmult]->GetNbinsX(); bin ++){
    		if(SpectraFior[nmult]->GetBinContent(bin)>0) hClone[nmult]->SetBinContent(bin,SpectraPart[nmult]->GetBinContent(bin)/SpectraFior[nmult]->GetBinContent(bin));
    		hClone[nmult]->SetBinError(bin,(ErrorInRatio(SpectraPart[nmult]->GetBinContent(bin),SpectraPart[nmult]->GetBinError(bin), 
    			SpectraFior[nmult]->GetBinContent(bin), SpectraFior[nmult]->GetBinError(bin))));
    		}
	  }

	  TCanvas* d = new TCanvas();
	  
	  hClone[0]->GetYaxis()->SetRangeUser(0.9,1.1);
	  hClone[0]->GetYaxis()->SetTitle("");
	  hClone[0]->GetYaxis()->SetTitleOffset(1.2);
	  hClone[0]->SetStats(0);
	  for(int nmult = 0; nmult < multbinnumb; nmult++){
	  	hClone[nmult]->SetMarkerStyle(8);    
	  	hClone[nmult]->Draw("SAME");
	  }

	  TCanvas* c1 = new TCanvas();
	  c1->SetLogy();
	  SpectraPart[9]->Draw();
	  SpectraFior[9]->Draw("SAME");


}


//---------------------------------------------------------------------------------------------------
double ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr) {	
  // Error in a Ratio
  if (B != 0) {
    Double_t errorfromtop = Aerr * Aerr / (B * B);
    Double_t errorfrombottom = ((A * A) / (B * B * B * B)) * Berr * Berr;
    return TMath::Sqrt(TMath::Abs(errorfromtop - errorfrombottom));
  }
  return 1.;
}