void makepedzdc(){
  TFile *f[3];
  f[0] = TFile::Open("NchRawContainer_15f.root");
  f[1] = TFile::Open("NchRawContainer_17j.root");
  f[2] = TFile::Open("NchRawContainer_18i.root");

  TProfile *hzna[3],*hznc[3],*hzpa[3],*hzpc[3];
  for(int i=0; i < 3; i++){
    hzna[i] = (TProfile *) f[i]->Get("zdcnaPed");
    hznc[i] = (TProfile *) f[i]->Get("zdcncPed");
    hzpa[i] = (TProfile *) f[i]->Get("zdcpaPed");
    hzpc[i] = (TProfile *) f[i]->Get("zdcpcPed");

    if(i>0){
      hzna[0]->Add(hzna[i]);
      hznc[0]->Add(hznc[i]);
      hzpa[0]->Add(hzpa[i]);
      hzpc[0]->Add(hzpc[i]);
    }
  }

  TFile *fout = new TFile("pedZDC.root","RECREATE");
  hzna[0]->Write();
  hznc[0]->Write();
  hzpa[0]->Write();
  hzpc[0]->Write();
  fout->Close();
}
