void makecalzdc(){
  TFile *f[3];
  f[0] = TFile::Open("NchRawContainer_15f.root");
  f[1] = TFile::Open("NchRawContainer_17j.root");
  f[2] = TFile::Open("NchRawContainer_18i.root");

  TProfile *hzna[3],*hznc[3],*hzpa[3],*hzpc[3];
  for(int i=0; i < 3; i++){
    hzna[i] = (TProfile *) f[i]->Get("zdcna");
    hznc[i] = (TProfile *) f[i]->Get("zdcnc");
    hzpa[i] = (TProfile *) f[i]->Get("zdcpa");
    hzpc[i] = (TProfile *) f[i]->Get("zdcpc");

    if(i==0){
      hzpa[0]->Scale(1.05);
      hzpc[0]->Scale(1.05);
    }

    if(i>0){
      hzna[0]->Add(hzna[i]);
      hznc[0]->Add(hznc[i]);
      hzpa[0]->Add(hzpa[i]);
      hzpc[0]->Add(hzpc[i]);
    }
  }

  TFile *fout = new TFile("calZDC.root","RECREATE");
  hzna[0]->Write();
  hznc[0]->Write();
  hzpa[0]->Write();
  hzpc[0]->Write();
  fout->Close();
}
