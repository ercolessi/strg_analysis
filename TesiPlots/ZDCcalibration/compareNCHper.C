void compareNCHper(){
  TFile *f[3];
  f[0] = TFile::Open("NchRawContainer_15f.root");
  f[1] = TFile::Open("NchRawContainer_17j.root");
  f[2] = TFile::Open("NchRawContainer_18i.root");

  int color[3] = {1, 2, 4};

  const int ncentrSPD = 1;
  int centrSPD[ncentrSPD+1] = {0,100};
  const int ncentrV0 = 8;
  int centrV0[ncentrV0+1] = {0,10,20,30,40,50,60,70,100};
  int nselections = ncentrSPD*ncentrV0;

  float mb[3] = {4.13995, 5.48181, 5.12462};

  float PrToNe[3][ncentrV0], xPrToNe[3][ncentrV0];

  TH3F *hzn[3],*hzp[3], *h3;
  TH1D *hznp[3],*hzpp[3];
  float x[3][100],y[3][100],yAv[100],xerr[100],yerr[100];
  for(int i=0; i < 3; i++){
    hzn[i] = (TH3F *) f[i]->Get("hspd_spdv0m");
    hzp[i] = (TH3F *) f[i]->Get("hspd_spdv0m");
    h3 = (TH3F *) f[i]->Get("hspd_spdv0m");
    hznp[i] = hzn[i]->ProjectionX(Form("hznp_%d",i));
    hzpp[i] = hzp[i]->ProjectionX(Form("hzpp_%d",i));
    hznp[i]->SetLineColor(color[i]);
    hzpp[i]->SetLineColor(color[i]);
    hznp[i]->SetLineWidth(2);
    hzpp[i]->SetLineWidth(2);
    hznp[i]->SetStats(0);
    hzpp[i]->SetStats(0);
    hznp[i]->SetTitle("");
    hzpp[i]->SetTitle("");
    hznp[i]->GetXaxis()->SetTitle("<ZN Sum> (a.u.)");
    hzpp[i]->GetXaxis()->SetTitle("<ZP Sum> (a.u.)");

    int n=0;
    for(int is=0;is < ncentrSPD; is++){
      for(int iv=0;iv < ncentrV0; iv++){
        x[i][n] = i - 1 + (centrV0[iv]+centrV0[iv+1])*0.5;
        y[i][n] = h3->ProjectionX("temp",1+centrSPD[is],centrSPD[is+1],centrV0[iv]+1,centrV0[iv+1])->GetMean() / mb[i];
        n++;
      }
    }
    for(int iv=0;iv < ncentrV0; iv++){
      xPrToNe[i][iv] = i - 1 + (centrV0[iv]+centrV0[iv+1])*0.5;
      PrToNe[i][iv] = hzp[i]->ProjectionX("temp",1,100,centrV0[iv]+1,centrV0[iv+1])->GetMean() / hzn[i]->ProjectionX("temp",1,100,centrV0[iv]+1,centrV0[iv+1])->GetMean();
    }
  }

  TCanvas *c = new TCanvas;
  c->Divide(2,1);
  c->cd(1)->SetLogy();
  for(int i=0; i < 3; i++){
    if(i==0) hznp[i]->DrawNormalized();
    else hznp[i]->DrawNormalized("SAME");
  }
  TLegend *leg = new TLegend(0.5,0.5,0.8,0.8);
  leg->AddEntry(hznp[0],"LHC15f","l");
  leg->AddEntry(hznp[1],"LHC17j","l");
  leg->AddEntry(hznp[2],"LHC18i","l");
  leg->SetFillStyle(0);
  leg->Draw("SAME");
  c->cd(2)->SetLogy();
  for(int i=0; i < 3; i++){
    if(i==0) hzpp[i]->DrawNormalized();
    else hzpp[i]->DrawNormalized("SAME");
  }

  TCanvas *c2 = new TCanvas;
  TGraph *g1 = new TGraph(nselections, x[0], y[0]);
  g1->SetMarkerStyle(20);
  g1->SetMarkerColor(color[0]);
  g1->SetLineColor(color[0]);
  g1->Draw("AP");
  TGraph *g2 = new TGraph(nselections, x[1], y[1]);
  g2->SetMarkerStyle(21);
  g2->SetMarkerColor(color[1]);
  g2->SetLineColor(color[1]);
  g2->Draw("P");
  TGraph *g3 = new TGraph(nselections, x[2], y[2]);
  g3->SetMarkerStyle(22);
  g3->SetMarkerColor(color[2]);
  g3->SetLineColor(color[2]);
  g3->Draw("P");
  g1->SetTitle("; V0 (%); N_{tracklets} / <N_{tracklets}>_{MB}");

  TLegend *leg2 = new TLegend(0.5,0.5,0.8,0.8);
  leg2->AddEntry(g1,"LHC15f","lp");
  leg2->AddEntry(g2,"LHC17j","lp");
  leg2->AddEntry(g3,"LHC18i","lp");
  leg2->SetFillStyle(0);
  leg2->Draw("SAME");

  TH1F *sigmaperc = new TH1F("sigmaperc",";#sigma for <ZDC Sum> (%)",100,0,10);
  
  for(int i=0; i < nselections; i++){
    float a = g1->GetY()[i];
    float b = g2->GetY()[i];
    float c = g3->GetY()[i];

    printf("%f %f\n",g1->GetX()[i],b/a);

    float mean = (a+b+c)/3;

    float maxd = TMath::Abs(a-mean);
    if(TMath::Abs(b-mean) > maxd){
      maxd = TMath::Abs(b-mean);
    }
    if(TMath::Abs(c-mean) > maxd){
      maxd = TMath::Abs(c-mean);
    }


    float sigma = (a-mean)*(a-mean);
    sigma += (b-mean)*(b-mean);
    sigma += (c-mean)*(c-mean);

    sigma = sqrt(sigma/2)/mean * 100;
    sigmaperc->Fill(maxd);//sigma);

    yAv[i] = mean;
    yerr[i] = maxd;//sigma;
    xerr[i] = 1;//g2->GetX()[i]*0.05;
  }

  TGraphErrors *gAv = new TGraphErrors(nselections,x[1],yAv,xerr,yerr);
  gAv->Draw("p,err2");
  gAv->SetName("gAv");
  gAv->SetFillStyle(0);
  gAv->SetLineColor(1);
  gAv->SetLineWidth(2);
  leg2->AddEntry(gAv,"mean + max difference as uncertainty","f");

  return;

  new TCanvas;
  sigmaperc->Draw();
  sigmaperc->SetStats(0);

  new TCanvas;
  TGraph *g1b = new TGraph(ncentrV0, xPrToNe[0], PrToNe[0]);
  g1b->SetMarkerStyle(20);
  g1b->SetMarkerColor(color[0]);
  g1b->SetLineColor(color[0]);
  g1b->Draw("AP");
  TGraph *g2b = new TGraph(ncentrV0, xPrToNe[1], PrToNe[1]);
  g2b->SetMarkerStyle(21);
  g2b->SetMarkerColor(color[1]);
  g2b->SetLineColor(color[1]);
  g2b->Draw("P");
  TGraph *g3b = new TGraph(ncentrV0, xPrToNe[2], PrToNe[2]);
  g3b->SetMarkerStyle(22);
  g3b->SetMarkerColor(color[2]);
  g3b->SetLineColor(color[2]);
  g3b->Draw("P");
  leg2->Draw("SAME");
}
