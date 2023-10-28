void compareZDCper(){
  TFile *f[3];
  f[0] = TFile::Open("NchRawContainer_15f.root");
  f[1] = TFile::Open("NchRawContainer_17j.root");
  f[2] = TFile::Open("NchRawContainer_18i.root");

  int color[3] = {1, 2, 4};

  const int ncentrSPD = 7;
  int centrSPD[ncentrSPD+1] = {0,10,20,30,40,50,70,100};
  const int ncentrV0 = 7;
  int centrV0[ncentrV0+1] = {0,10,20,30,40,50,70,100};
  int nselections = ncentrSPD*ncentrV0;

  float PrToNe[3][ncentrV0], xPrToNe[3][ncentrV0];

  TH3F *hzn[3],*hzp[3], *h3;
  TH1D *hznp[3],*hzpp[3];
  float x[3][100],y[3][100],yAv[100],xerr[100],yerr[100];
  float xv0_1[10][10], yv0_1[10][10], yv0err_1[10][10];
  float xv0_2[10][10], yv0_2[10][10], yv0err_2[10][10];
  float xv0_3[10][10], yv0_3[10][10], yv0err_3[10][10];
  for(int i=0; i < 3; i++){
    hzn[i] = (TH3F *) f[i]->Get("hznsum_spdv0m");
    hzp[i] = (TH3F *) f[i]->Get("hzpsum_spdv0m");
    h3 = (TH3F *) f[i]->Get("hznsum_spdv0m"); // change to hzdcsum_spdv0m for ZDCN+ZDCP
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
    hznp[i]->GetXaxis()->SetTitle("ZN (a.u.)");
    hzpp[i]->GetXaxis()->SetTitle("ZP (a.u.)");
    hznp[i]->GetYaxis()->SetTitle("Normalized counts");
    hzpp[i]->GetYaxis()->SetTitle("Normalized counts");
    hznp[i]->GetXaxis()->SetRangeUser(0, 3000.);
    hzpp[i]->GetXaxis()->SetRangeUser(0, 2500.);
    hznp[i]->GetXaxis()->SetTitleSize(0.045);
    hzpp[i]->GetXaxis()->SetTitleSize(0.045);

    int n=0;
    for(int is=0;is < ncentrSPD; is++){
      for(int iv=0;iv < ncentrV0; iv++){
        x[i][n] =  (centrSPD[is]+centrSPD[is+1])*0.5;
        y[i][n] = h3->ProjectionX("temp",1+centrSPD[is],centrSPD[is+1],centrV0[iv]+1,centrV0[iv+1])->GetMean();
        n++;

        if (i==1) {
          xv0_1[iv][is] = (centrSPD[is]+centrSPD[is+1])*0.5;
          yv0_1[iv][is] = h3->ProjectionY("temp",1+centrSPD[is],centrSPD[is+1],centrV0[iv]+1,centrV0[iv+1])->GetMean();
        } else if (i==2) {
          xv0_2[iv][is] = (centrSPD[is]+centrSPD[is+1])*0.5;
          yv0_2[iv][is] = h3->ProjectionY("temp",1+centrSPD[is],centrSPD[is+1],centrV0[iv]+1,centrV0[iv+1])->GetMean();
        } else if (i==3) {
          xv0_3[iv][is] = (centrSPD[is]+centrSPD[is+1])*0.5;
          yv0_3[iv][is] = h3->ProjectionY("temp",1+centrSPD[is],centrSPD[is+1],centrV0[iv]+1,centrV0[iv+1])->GetMean();
        }
      }
    }
    for(int iv=0;iv < ncentrV0; iv++){
      xPrToNe[i][iv] = i - 1 + (centrV0[iv]+centrV0[iv+1])*0.5;
      PrToNe[i][iv] = hzp[i]->ProjectionX("temp",1,100,centrV0[iv]+1,centrV0[iv+1])->GetMean() / hzn[i]->ProjectionX("temp",1,100,centrV0[iv]+1,centrV0[iv+1])->GetMean();
    }
  }

  TCanvas *c = new TCanvas("c","c",1400,800);
  c->Divide(2,1);
  c->cd(1)->SetLeftMargin(0.15);
  c->cd(1)->SetBottomMargin(0.15);
  c->cd(1)->SetRightMargin(0.05);
  c->cd(1)->SetTopMargin(0.05);
  c->cd(1)->SetTicks();
  c->cd(2)->SetLeftMargin(0.15);
  c->cd(2)->SetBottomMargin(0.15);
  c->cd(2)->SetRightMargin(0.05);
  c->cd(2)->SetTopMargin(0.05);
  c->cd(2)->SetTicks();
  c->cd(1)->SetLogy();
  for(int i=0; i < 3; i++){
    if(i==0) hznp[i]->DrawNormalized();
    else hznp[i]->DrawNormalized("SAME");
  }
  TLegend *leg = new TLegend(0.5,0.7,0.8,0.85);
  leg->AddEntry(hznp[0],"2015 period","l");
  leg->AddEntry(hznp[1],"2017 period","l");
  leg->AddEntry(hznp[2],"2018 period","l");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw("SAME");
  TLatex *tex = new TLatex(0.2, 0.89, "ALICE, pp #sqrt{s} = 13 TeV");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.035);
  tex->SetLineWidth(2);
  TLatex *tex2 = new TLatex(0.2, 0.85, "This work");
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.035);
  tex2->SetLineWidth(2);
  tex->Draw();
  tex2->Draw();
  c->cd(2)->SetLogy();
  for(int i=0; i < 3; i++){
    if(i==0) hzpp[i]->DrawNormalized();
    else hzpp[i]->DrawNormalized("SAME");
  }
  leg->Draw("SAME");
  tex->Draw();
  tex2->Draw();

  c->SaveAs("comparecalib.pdf");

  TCanvas *c2 = new TCanvas("c2","c2",1400,800);
  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.15);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.05);
  c2->SetTicks();

  TGraph *g1_[ncentrV0];
  TGraph *g2_[ncentrV0];
  TGraph *g3_[ncentrV0];
  for (int i=0; i<ncentrV0; i++) {
    g1_[i] = new TGraph(ncentrSPD, xv0_1[0], yv0_1[i]);
    g1_[i]->SetMarkerStyle(20);
    g1_[i]->SetMarkerColor(color[0]);
    g1_[i]->SetLineColor(color[0]);
    g1_[i]->SetLineWidth(2);
    g1_[i]->SetFillStyle(0);
    g1_[i]->SetTitle("; SPDClusters (%); #LT ZN #GT (a.u.)");
    g2_[i] = new TGraph(ncentrSPD, xv0_2[1], yv0_2[i]);
    g2_[i]->SetMarkerStyle(21);
    g2_[i]->SetMarkerColor(color[1]);
    g2_[i]->SetLineColor(color[1]);
    g2_[i]->SetLineWidth(2);
    g2_[i]->SetFillStyle(0);
    g3_[i] = new TGraph(ncentrSPD, xv0_3[2], yv0_3[i]);
    g3_[i]->SetMarkerStyle(22);
    g3_[i]->SetMarkerColor(color[2]);
    g3_[i]->SetLineColor(color[2]);
    g3_[i]->SetLineWidth(2);
    g3_[i]->SetFillStyle(0);
  }

  g1_[0]->Draw("AP");

  TGraph *g1 = new TGraph(nselections, x[0], y[0]);
  g1->SetMarkerStyle(20);
  g1->SetMarkerColor(color[0]);
  g1->SetLineColor(color[0]);
  g1->SetMarkerSize(1.3);
  TGraph *g2 = new TGraph(nselections, x[1], y[1]);
  g2->SetMarkerStyle(21);
  g2->SetMarkerColor(color[1]);
  g2->SetLineColor(color[1]);
  g2->SetMarkerSize(1.3);
  TGraph *g3 = new TGraph(nselections, x[2], y[2]);
  g3->SetMarkerStyle(22);
  g3->SetMarkerColor(color[2]);
  g3->SetLineColor(color[2]);
  g3->SetMarkerSize(1.3);
  g1->SetTitle("; SPDClusters (%); #LT ZN #GT (a.u.)"); // adjust title if you change to ZDCN+ZDCP
  g1->GetYaxis()->SetTitleSize(0.05);
  g1->GetXaxis()->SetTitleSize(0.05);
  g1->GetXaxis()->SetRangeUser(0, 100);
  TLegend *leg2 = new TLegend(0.75,0.2,0.89,0.4);
  leg2->AddEntry(g1,"2015 period","lp");
  leg2->AddEntry(g2,"2017 period","lp");
  leg2->AddEntry(g3,"2018 period","lp");
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.035);
  leg2->Draw("SAME");

  TH1F *sigmaperc = new TH1F("sigmaperc",";#sigma for <ZDC Sum> (%)",100,-0.1,0.1);

  for(int i=0; i < nselections; i++){
    float a = g1->GetY()[i];
    float b = g2->GetY()[i];
    float c = g3->GetY()[i];

    float mean = (a+b+c)/3;

    float maxa = (a-mean)/mean;
    float maxb = (b-mean)/mean;
    float maxc = (c-mean)/mean;
    float maxd = a-mean;
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
    sigmaperc->Fill(maxa);//sigma);
    sigmaperc->Fill(maxb);//sigma);
    sigmaperc->Fill(maxc);//sigma);

    yAv[i] = mean;
    yerr[i] = maxd;//sigma;
    xerr[i] = 1;//g2->GetX()[i]*0.05;
  }

  TGraphErrors *gAv = new TGraphErrors(nselections,x[1],yAv,xerr,yerr);
  gAv->SetName("gAv");
  gAv->SetLineWidth(1);
  gAv->SetFillColor(kGray);
  gAv->SetMarkerStyle(20);
  gAv->SetMarkerSize(0.1);
  gAv->SetLineColor(kBlack);
  gAv->SetFillStyle(4050);
  gAv->SetLineStyle(2);
  TGraphErrors *gAv2 = new TGraphErrors(nselections, x[1], yAv, xerr, yerr);
  gAv2->SetName("gAv2");
  gAv2->SetLineWidth(1);
  gAv2->SetMarkerStyle(20);
  gAv2->SetMarkerSize(0.1);
  gAv2->SetLineColor(kBlack);
  gAv2->SetFillStyle(0);
  g1->Draw("AP");
  leg2->AddEntry(gAv2, "mean + max dev.", "F");

  //gAv2->Draw("PE2 same");
  gAv2->Draw("PE2 same");
  g1->Draw("PE SAME");
  g2->Draw("PE SAME");
  g3->Draw("PE SAME");
  leg2->Draw("SAME");
  tex->Draw();
  tex2->Draw();

  c2->SaveAs("znsystematics.pdf");


  TCanvas* d = new TCanvas("d","d",800,600);
  sigmaperc->Draw();
  sigmaperc->SetStats(0);
  return;

  /*new TCanvas;
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
  leg2->Draw("SAME");*/
}
