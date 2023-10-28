const char *estimator[4] ={"v0mPerc", "multSPDcl", "SPDtracklets", "TOFclusters"};
float maxVal[4] = {100, 100, 60, 400};
const char *title[4] = {";V0M percentile (%);kaons","SPDcl percentile (%);kaons",";SPD tracklets;kaons",";TOF clusters (%);kaons"};

void checkBias(int iEst=3){
  TFile *f = new TFile("LHC20i2a-Pythia8_Monash2013.root");
  TTree *t = (TTree *) f->Get("fTree");

  int nev = t->GetEntries();

  TProfile *hKaC = new TProfile("hKaC",title[iEst],100,0,maxVal[iEst]);
  TProfile *hKa0 = new TProfile("hKa0",title[iEst],100,0,maxVal[iEst]);

  for(int i=0; i < nev/100; i++){
    t->GetEvent(i);
    if(! t->GetLeaf("inelGT0")->GetValue()){
      continue;
    }
    float x = t->GetLeaf(estimator[iEst])->GetValue();
    hKaC->Fill(x, t->GetLeaf("nKchEta")->GetValue());
    hKa0->Fill(x, t->GetLeaf("nK0Eta")->GetValue());
  }

  hKaC->Draw();
  hKa0->Draw("SAME");
  hKaC->SetLineColor(2);
  hKa0->SetLineColor(1);
  hKaC->SetLineWidth(2);
  hKa0->SetLineWidth(2);
  hKaC->SetStats(0);
}
