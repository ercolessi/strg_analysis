void lambda(const char *mc="1epos", const string centr="1030"){
  const int nV0cent = 9;
  //int V0centr[nV0cent+1] = {0,10,40,70,100};
  int V0centr[nV0cent+1] = {0,5,10,15,20,30,40,50,70,100};

  const int nZDCcent = 6;
  int ZDCcentr[nZDCcent+1] = {0,20,30,40,50,70,100};

  const int nSPDcent = 11;
  int SPDcentr[nSPDcent+1] = {0,1,5,10,15,20,30,40,50,60,70,100};

  TFile *fvzero = new TFile(Form("out%sLambda_spd%s_vzero.root",mc,centr.c_str()));
  TFile *fzdc = new TFile(Form("out%sLambda_spd%s_zdc.root",mc,centr.c_str()));
  TFile *fall = new TFile(Form("out%sLambda_spd.root",mc));

  TGraphErrors *gvzero = (TGraphErrors *) fvzero->Get("Graph");
  TGraphErrors *gzdc = (TGraphErrors *) fzdc->Get("Graph");
  TGraphErrors *gall = (TGraphErrors *) fall->Get("Graph");

  float sw = 0;
  float valX = 0;
  float valY = 0;

  if(centr.find("1030") == 0){
    valX += gall->GetX()[3] / (gall->GetEX()[3]*gall->GetEX()[3]);
    valY += gall->GetY()[3] / (gall->GetEX()[3]*gall->GetEX()[3]);
    sw += 1. / (gall->GetEX()[3]*gall->GetEX()[3]);

    valX += gall->GetX()[4] / (gall->GetEX()[4]*gall->GetEX()[4]);
    valY += gall->GetY()[4] / (gall->GetEX()[4]*gall->GetEX()[4]);
    sw += 1. / (gall->GetEX()[4]*gall->GetEX()[4]);

    valX += gall->GetX()[5] / (gall->GetEX()[5]*gall->GetEX()[5]);
    valY += gall->GetY()[5] / (gall->GetEX()[5]*gall->GetEX()[5]);
    sw += 1. / (gall->GetEX()[5]*gall->GetEX()[5]);
  } else if(centr.find("4050") == 0) {
    valX += gall->GetX()[7] / (gall->GetEX()[7]*gall->GetEX()[7]);
    valY += gall->GetY()[7] / (gall->GetEX()[7]*gall->GetEX()[7]);
    sw += 1. / (gall->GetEX()[7]*gall->GetEX()[7]);
  }

  valX /= sw;
  valY /= sw;

  for(int i=0; i < nV0cent; i++){
      gvzero->SetPoint(i, (V0centr[i] + V0centr[i+1])*0.5,  gvzero->GetY()[i] / valY);
      gvzero->SetPointError(i, (-V0centr[i] + V0centr[i+1])*0.5, gvzero->GetEY()[i] / valY);
  }
  for(int i=0; i < nZDCcent; i++){
      gzdc->SetPoint(i, (ZDCcentr[i] + ZDCcentr[i+1])*0.5,  gzdc->GetY()[i] / valY);
      gzdc->SetPointError(i, (-ZDCcentr[i] + ZDCcentr[i+1])*0.5, gzdc->GetEY()[i] / valY);
  }
  printf("X=%f - Y=%f\n",valX,valY);


  gvzero->Draw("APL");

  new TCanvas;
  gzdc->Draw("APL");

  TFile *fout = new TFile(Form("outLambda%s_%s.root",mc,centr.c_str()),"RECREATE");
  gvzero->Write();
  gzdc->Write();
  fout->Close();
}
