void drawEE(){
  TFile *f[4];
  TTree *t[4];
  const char *filename[4] = {"LHC20i2a-Pythia8_Monash2013.root", "LHC16d3-EPOS.root", "LHC17h7a-Pythia6_Perugia2011.root", "LHC16d3-EPOS.root"};

  int centrSPD1020min[3] = {0,20,50};
  int centrSPD1020max[3] = {5,30,100};

  int centrSPD4050min[3] = {0,40,70};
  int centrSPD4050max[3] = {20,50,100};

  int centrV01020min[3] = {0,20,50};
  int centrV01020max[3] = {5,30,100};

  int centrV04050min[3] = {0,30,70};
  int centrV04050max[3] = {10,40,100};

  TFile *fout = new TFile("eeMC.root","recreate");
  TH1F *hMBZ[4];
  TH1F *hMB[4];
  
  TH1F *hSPD1020Zlow[4];
  TH1F *hSPD1020Zmid[4];
  TH1F *hSPD1020Zhigh[4];
  TH1F *hSPD1020low[4];
  TH1F *hSPD1020mid[4];
  TH1F *hSPD1020high[4];

  TH1F *hSPD4050Zlow[4];
  TH1F *hSPD4050Zmid[4];
  TH1F *hSPD4050Zhigh[4];
  TH1F *hSPD4050low[4];
  TH1F *hSPD4050mid[4];
  TH1F *hSPD4050high[4];

  TH1F *hV01020Zlow[4];
  TH1F *hV01020Zmid[4];
  TH1F *hV01020Zhigh[4];
  TH1F *hV01020low[4];
  TH1F *hV01020mid[4];
  TH1F *hV01020high[4];

  TH1F *hV04050Zlow[4];
  TH1F *hV04050Zmid[4];
  TH1F *hV04050Zhigh[4];
  TH1F *hV04050low[4];
  TH1F *hV04050mid[4];
  TH1F *hV04050high[4];

  TH2F *hMatrix[4];
  
  for(int i=0; i < 4; i++){
    printf("MC %d\n",i);
    fout->cd();
    hMatrix[i] = new TH2F(Form("hMatrix%d",i),";<ZDC sum> (a.u.); Effective Energy (GeV)",200,0,2000,200,0,13000);
    hMBZ[i] = new TH1F(Form("hMBZ%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hMB[i] = new TH1F(Form("hMB%d",i),";Effective Energy (GeV)",200,0,13000);

    hSPD1020Zlow[i] = new TH1F(Form("hSPD1020Zlow%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hSPD1020Zmid[i] = new TH1F(Form("hSPD1020Zmid%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hSPD1020Zhigh[i] = new TH1F(Form("hSPD1020Zhigh%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hSPD1020low[i] = new TH1F(Form("hSPD1020low%d",i),";Effective Energy (GeV)",200,0,13000);
    hSPD1020mid[i] = new TH1F(Form("hSPD1020mid%d",i),";Effective Energy (GeV)",200,0,13000);
    hSPD1020high[i] = new TH1F(Form("hSPD1020high%d",i),";Effective Energy (GeV)",200,0,13000);

    hSPD4050Zlow[i] = new TH1F(Form("hSPD4050Zlow%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hSPD4050Zmid[i] = new TH1F(Form("hSPD4050Zmid%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hSPD4050Zhigh[i] = new TH1F(Form("hSPD4050Zhigh%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hSPD4050low[i] = new TH1F(Form("hSPD4050low%d",i),";Effective Energy (GeV)",200,0,13000);
    hSPD4050mid[i] = new TH1F(Form("hSPD4050mid%d",i),";Effective Energy (GeV)",200,0,13000);
    hSPD4050high[i] = new TH1F(Form("hSPD4050high%d",i),";Effective Energy (GeV)",200,0,13000);

    hV01020Zlow[i] = new TH1F(Form("hV01020Zlow%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hV01020Zmid[i] = new TH1F(Form("hV01020Zmid%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hV01020Zhigh[i] = new TH1F(Form("hV01020Zhigh%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hV01020low[i] = new TH1F(Form("hV01020low%d",i),";Effective Energy (GeV)",200,0,13000);
    hV01020mid[i] = new TH1F(Form("hV01020mid%d",i),";Effective Energy (GeV)",200,0,13000);
    hV01020high[i] = new TH1F(Form("hV01020high%d",i),";Effective Energy (GeV)",200,0,13000);

    hV04050Zlow[i] = new TH1F(Form("hV04050Zlow%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hV04050Zmid[i] = new TH1F(Form("hV04050Zmid%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hV04050Zhigh[i] = new TH1F(Form("hV04050Zhigh%d",i),";<ZDC sum> (a.u.)",200,0,2000);
    hV04050low[i] = new TH1F(Form("hV04050low%d",i),";Effective Energy (GeV)",200,0,13000);
    hV04050mid[i] = new TH1F(Form("hV04050mid%d",i),";Effective Energy (GeV)",200,0,13000);
    hV04050high[i] = new TH1F(Form("hV04050high%d",i),";Effective Energy (GeV)",200,0,13000);

    f[i] = new TFile(filename[i]);
    t[i] = (TTree *) f[i]->Get("fTree");

    for(int j=0; j < t[i]->GetEntries(); j++){
      t[i]->GetEvent(j);
      if(t[i]->GetLeaf("SPDtracklets")->GetValue() < 1) continue;
      
      float v0m = t[i]->GetLeaf("v0mPerc")->GetValue();
      float spd = t[i]->GetLeaf("multSPDcl")->GetValue();
      float ee = 13000 + t[i]->GetLeaf("effEnergy")->GetValue();
      float zdc = t[i]->GetLeaf("adcZDCN1")->GetValue(0) + t[i]->GetLeaf("adcZDCN2")->GetValue(0) + t[i]->GetLeaf("adcZDCP1")->GetValue(0) + t[i]->GetLeaf("adcZDCP2")->GetValue(0);

      hMatrix[i]->Fill(zdc, ee);
      hMBZ[i]->Fill(zdc);
      hMB[i]->Fill(ee);
      
      if(spd > 10 && spd < 20){
	if(v0m > centrSPD1020min[0] && v0m < centrSPD1020max[0]){
	  hSPD1020Zhigh[i]->Fill(zdc);
	  hSPD1020high[i]->Fill(ee);
	} else if(v0m > centrSPD1020min[1] && v0m < centrSPD1020max[1]) {
	  hSPD1020Zmid[i]->Fill(zdc);
	  hSPD1020mid[i]->Fill(ee);
	} else if(v0m > centrSPD1020min[2] && v0m < centrSPD1020max[2]) {
	  hSPD1020Zlow[i]->Fill(zdc);
	  hSPD1020low[i]->Fill(ee);
	}
      }

      if(spd > 40 && spd < 50){
	if(v0m > centrSPD4050min[0] && v0m < centrSPD4050max[0]){
	  hSPD4050Zhigh[i]->Fill(zdc);
	  hSPD4050high[i]->Fill(ee);
	} else if(v0m > centrSPD4050min[1] && v0m < centrSPD4050max[1]) {
	  hSPD4050Zmid[i]->Fill(zdc);
	  hSPD4050mid[i]->Fill(ee);
	} else if(v0m > centrSPD4050min[2] && v0m < centrSPD4050max[2]) {
	  hSPD4050Zlow[i]->Fill(zdc);
	  hSPD4050low[i]->Fill(ee);
	}
      }

      if(v0m > 10 && v0m < 20){
	if(spd > centrV01020min[0] && spd < centrV01020max[0]){
	  hV01020Zhigh[i]->Fill(zdc);
	  hV01020high[i]->Fill(ee);
	} else if(spd > centrV01020min[1] && spd < centrV01020max[1]) {
	  hV01020Zmid[i]->Fill(zdc);
	  hV01020mid[i]->Fill(ee);
	} else if(spd > centrV01020min[2] && spd < centrV01020max[2]) {
	  hV01020Zlow[i]->Fill(zdc);
	  hV01020low[i]->Fill(ee);
	}
      }

      if(v0m > 40 && v0m < 50){
	if(spd > centrV04050min[0] && spd < centrV04050max[0]){
	  hV04050Zhigh[i]->Fill(zdc);
	  hV04050high[i]->Fill(ee);
	} else if(spd > centrV04050min[1] && spd < centrV04050max[1]) {
	  hV04050Zmid[i]->Fill(zdc);
	  hV04050mid[i]->Fill(ee);
	} else if(spd > centrV04050min[2] && spd < centrV04050max[2]) {
	  hV04050Zlow[i]->Fill(zdc);
	  hV04050low[i]->Fill(ee);
	}
      }

    }
    f[i]->Close();
  }

  fout->cd();
  for(int i=0; i < 4; i++){
    hMatrix[i]->Write();
    hMBZ[i]->Write();
    hMB[i]->Write();
    hSPD1020Zlow[i]->Write();
    hSPD1020Zmid[i]->Write();
    hSPD1020Zhigh[i]->Write();
    hSPD1020low[i]->Write();
    hSPD1020mid[i]->Write();
    hSPD1020high[i]->Write();
    hSPD4050Zlow[i]->Write();
    hSPD4050Zmid[i]->Write();
    hSPD4050Zhigh[i]->Write();
    hSPD4050low[i]->Write();
    hSPD4050mid[i]->Write();
    hSPD4050high[i]->Write();
    hV01020Zlow[i]->Write();
    hV01020Zmid[i]->Write();
    hV01020Zhigh[i]->Write();
    hV01020low[i]->Write();
    hV01020mid[i]->Write();
    hV01020high[i]->Write();
    hV04050Zlow[i]->Write();
    hV04050Zmid[i]->Write();
    hV04050Zhigh[i]->Write();
    hV04050low[i]->Write();
    hV04050mid[i]->Write();
    hV04050high[i]->Write();
  }
}
