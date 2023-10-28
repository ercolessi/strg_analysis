void plotZDCperc(int imc=0){
  const char *input[4] = {"LHC20i2a-Pythia8_Monash2013.root", "LHC16d3-EPOS.root", "LHC17h7a-Pythia6_Perugia2011.root", "LHC17h7b-Phojet.root"};
  TFile *f = new TFile(input[imc]);
  TTree *t = (TTree *) f->Get("fTree");

  int nev = t->GetEntries();

  Float_t fCentrality_V0M;
  Float_t adcZDCN1;
  Float_t adcZDCN2;
  Float_t effEnergy;
  Int_t nchEta;

  // Get from Tree
  lTreeEvent->SetBranchAddress("v0mPerc", &fCentrality_V0M);
  lTreeEvent->SetBranchAddress("adcZDCN1", &adcZDCN1);
  lTreeEvent->SetBranchAddress("adcZDCN2", &adcZDCN2);
  lTreeEvent->SetBranchAddress("effEnergy", &effEnergy);
  lTreeEvent->SetBranchAddress("nchEta", &nchEta);

  TH1F *h[ncentr], *h2[ncentr];
  TH2F *h_2d, *h2_2d;
  h_2d = new TH2F(Form("h_2d"),";perc;E_{leading} (|#eta>8|) (GeV)",100,0,100,100,0,13000);
  h2_2d = new TH2F(Form("h2_2d"),";perc;dN_{ch}/d#eta",100,0,100,100,0,100);

  for(int i=0; i < ncentr; i++){
    h[i] = new TH1F(Form("h%d",i),";E_{leading} (|#eta>8|) (GeV)",100,0,13000);
    h2[i] = new TH1F(Form("h2%d",i),";dN_{ch}/d#eta",100,0,100);
  }

  float centr, centrCut;
  int icentr;

  for(int i=0; i < nev; i++){

    t->GetEvent(i);
    centrCut = t->GetLeaf(sel2)->GetValue();

    if(t->GetLeaf("SPDtracklets")->GetValue()< 1) continue;

    if(centrCut < cMin || centrCut > cMax) continue;

    centr = t->GetLeaf(sel)->GetValue();

    if(centr < 0 || centr >= 100) continue;

    h_2d->Fill(centr,TMath::Abs(t->GetLeaf("effEnergy")->GetValue()));
    h2_2d->Fill(centr,t->GetLeaf("nchEta")->GetValue());

    for (int n=0; n < ncentr; n++)
    {

      if (centr >= bin1[n] && centr < bin2[n])
      {
        h[n]->Fill(TMath::Abs(t->GetLeaf("effEnergy")->GetValue()));
        h2[n]->Fill(t->GetLeaf("nchEta")->GetValue());
      }

      }

    }


  h[0]->SetStats(0);
  h2[0]->SetStats(0);
  TCanvas *cc = new TCanvas("cDistr","");
  TLegend *leg = new TLegend(0.3,0.5,0.5,0.7);
  for(int i=0; i < ncentr; i++){
    h[i]->SetLineWidth(3);
    h[i]->SetLineColor(icolor[i]);
    h[i]->GetXaxis()->SetTitleSize(0.05);
    h[i]->GetYaxis()->SetTitle("Counts (normalised)");
    if(normalized){
      if(i==0) h[0]->DrawNormalized();
      else h[i]->DrawNormalized("SAME");
      cout << h[i]->GetMean() << endl;
    }
    else{
      if(i==0) h[0]->Draw();
      else h[i]->Draw("SAME");
    }
    leg->AddEntry(h[i],Form("%.0f-%.0f%c",bin1[i],bin2[i],'%'),"l");
  }
  leg->SetHeader(Form("%s scan (%s: %d-%d%c)",sel,sel2,cMin,cMax,'%'));

  if(cMin < 1 && cMax > 99)
    leg->SetHeader(Form("%s",sel));
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->Draw("SAME");

  TCanvas *cc2 = new TCanvas("cDistr2","");
  cc2->SetLogy();
  TLegend *leg2 = new TLegend(0.3,0.5,0.5,0.7);
  for(int i=0; i < ncentr; i++){
    h2[i]->SetLineWidth(3);
    h2[i]->SetLineColor(icolor[i]);
    h2[i]->GetXaxis()->SetTitleSize(0.05);
    h2[i]->GetYaxis()->SetTitle("Counts (normalised)");
    if(normalized){
      if(i==0) h2[0]->DrawNormalized();
      else h2[i]->DrawNormalized("SAME");
    }
    else{
      if(i==0) h2[0]->Draw();
      else h2[i]->Draw("SAME");
    }
    leg2->AddEntry(h2[i],Form("%d-%d%c",bin1[i],bin2[i+1],'%'),"l");
  }
  leg2->SetHeader(Form("%s scan (%s: %d-%d%c)",sel,sel2,cMin,cMax,'%'));

  if(cMin < 1 && cMax > 99)
    leg2->SetHeader(Form("%s",sel));

  leg2->Draw("SAME");

  TCanvas *cc3 = new TCanvas("cc3","",800,600);
  cc3->SetLogz();
  TH2D* h_2dclone = (TH2D*)h_2d->Clone("h_2dclone");
  h_2dclone->Rebin(10);
  TH2D* h2_2dclone = (TH2D*)h2_2d->Clone("h2_2dclone");
  h2_2dclone->Rebin(10);
  h_2d->SetStats(0);
  h_2d->Draw("colz");
  h_2d->GetXaxis()->SetTitle(Form("%s (%)",sel));
   TProfile* p_2d = h_2dclone->ProfileX();
  TProfile* p2_2d = h2_2dclone->ProfileX();
  p_2d->Draw("SAME");

  TFile *fout = new TFile(Form("distribution_%d_%s_%s_%d-%d.root",imc,sel,sel2,cMin,cMax),"RECREATE");
  cc->Write();
  cc2->Write();
  h_2d->Write();
  h2_2d->Write();


  p_2d->SetName("p_d");
  p2_2d->SetName("p_2d");
  p_2d->SetLineColor(kRed);
  p2_2d->SetLineColor(kRed);
  p_2d->SetLineWidth(2);
  p2_2d->SetLineWidth(2);
  p_2d->Write();
  p2_2d->Write();

  fout->Close();
}
