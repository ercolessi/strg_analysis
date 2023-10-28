const char *estimator[4] ={"v0mPerc", "multSPDcl", "SPDtracklets", "TOFclusters"};
float maxVal[4] = {100, 100, 60, 400};
const char *title[4] = {";V0M percentile (%);kaons","SPDcl percentile (%);kaons",";SPD tracklets;kaons",";TOF clusters (%);kaons"};

bool percInverted[4] = {0, 0, 1, 1};

int nbinEst = 10000;
TH1F *hPerc,*hEstim;

void doEstim(int iEst){

  float maxValue;
  for(int i=1; i <= nbinEst; i++){
    int k = i;
    int kBefore = i-1;
    if(percInverted[iEst]){
      k = nbinEst - i;
      kBefore = nbinEst - i + 1;
    }
    if(i == 1){
      hEstim->SetBinContent(k, hPerc->GetBinContent(k));
    } else {
      hEstim->SetBinContent(k, hPerc->GetBinContent(k) + hEstim->GetBinContent(kBefore));
    }

    if(i == nbinEst){
      maxValue = hEstim->GetBinContent(k);
    }
  }
  hEstim->Scale(100./maxValue);
}

void checkBias2(int iEst=3){
  TFile *f = new TFile("LHC20i2a-Pythia8_Monash2013.root");
  TTree *t = (TTree *) f->Get("fTree");

  int nev = t->GetEntries();

  TProfile *hKaC = new TProfile("hKaC",title[iEst],100,0,maxVal[iEst]);
  TProfile *hKa0 = new TProfile("hKa0",title[iEst],100,0,maxVal[iEst]);

  TProfile *hKaC_2 = new TProfile("hKaC_2",title[iEst],100,0,100);
  TProfile *hKa0_2 = new TProfile("hKa0_2",title[iEst],100,0,100);

  hPerc = new TH1F("hPerc","",nbinEst,0,maxVal[iEst]);
  hEstim = new TH1F("hEstim","",10000,0,maxVal[iEst]);

  int nsel = 0;
  float kaC = 0;
  float ka0 = 0;

  for(int i=0; i < nev/100; i++){
    t->GetEvent(i);
    if(! t->GetLeaf("inelGT0")->GetValue()){
      continue;
    }

    kaC += t->GetLeaf("nKchEta")->GetValue();
    ka0 += t->GetLeaf("nK0Eta")->GetValue();
    nsel++;

    float x = t->GetLeaf(estimator[iEst])->GetValue();
    hPerc->Fill(x);
  }

  kaC /= nsel;
  ka0 /= nsel;

  doEstim(iEst);

  for(int i=0; i < nev; i++){
    t->GetEvent(i);
    if(! t->GetLeaf("inelGT0")->GetValue()){
      continue;
    }
    float x = t->GetLeaf(estimator[iEst])->GetValue();
    float x2 = hEstim->Interpolate(x);
    hKaC->Fill(x, t->GetLeaf("nKchEta")->GetValue() / kaC);
    hKa0->Fill(x, t->GetLeaf("nK0Eta")->GetValue() / ka0);
    hKaC_2->Fill(x2, t->GetLeaf("nKchEta")->GetValue() / kaC);
    hKa0_2->Fill(x2, t->GetLeaf("nK0Eta")->GetValue() / ka0);
  }

  hKaC->Draw();
  hKa0->Draw("SAME");
  hKaC->SetLineColor(2);
  hKa0->SetLineColor(1);
  hKaC->SetLineWidth(2);
  hKa0->SetLineWidth(2);
  hKaC->SetStats(0);

  new TCanvas;
  hKaC_2->Draw();
  hKa0_2->Draw("SAME");
  hKaC_2->SetLineColor(2);
  hKa0_2->SetLineColor(1);
  hKaC_2->SetLineWidth(2);
  hKa0_2->SetLineWidth(2);
  hKaC_2->SetStats(0);
}
