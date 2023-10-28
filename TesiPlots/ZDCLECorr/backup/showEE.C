void showEE(int imc=0, const char *sel="4050", int type=0){//0=spd, 1=v0
  TFile *f = new TFile("eeMC.root");

  const char *prodAll[] = {"pythia8", "epos", "pythia6", "phojet"};

  const char *prod = prodAll[imc];
    
  TH2F *matrix;
  TH1F *lowZ,*midZ,*highZ,*mbZ;
  TH1F *low,*mid,*high,*mb;

  mb = (TH1F *) f->Get(Form("hMB%d",imc));
  mbZ = (TH1F *) f->Get(Form("hMBZ%d",imc));

  if(type == 0){
    matrix = (TH2F *) f->Get(Form("hMatrix%d",imc));
    lowZ = (TH1F *) f->Get(Form("hSPD%sZlow%d",sel,imc));
    midZ = (TH1F *) f->Get(Form("hSPD%sZmid%d",sel,imc));
    highZ = (TH1F *) f->Get(Form("hSPD%sZhigh%d",sel,imc));
    low = (TH1F *) f->Get(Form("hSPD%slow%d",sel,imc));
    mid = (TH1F *) f->Get(Form("hSPD%smid%d",sel,imc));
    high = (TH1F *) f->Get(Form("hSPD%shigh%d",sel,imc));
  } else {
    matrix = (TH2F *) f->Get(Form("hMatrix%d",imc));
    lowZ = (TH1F *) f->Get(Form("hV0%sZlow%d",sel,imc));
    midZ = (TH1F *) f->Get(Form("hV0%sZmid%d",sel,imc));
    highZ = (TH1F *) f->Get(Form("hV0%sZhigh%d",sel,imc));
    low = (TH1F *) f->Get(Form("hV0%slow%d",sel,imc));
    mid = (TH1F *) f->Get(Form("hV0%smid%d",sel,imc));
    high = (TH1F *) f->Get(Form("hV0%shigh%d",sel,imc));
  }

  matrix->SetTitle(prod);
  
  TH1D* h = new TH1D("h",";SPD tracklets;Counts (normalised)",10.,0.,35.);
  TH1D* h1 = new TH1D("h1",";ZDC Energy Sum (a.u.);Counts (normalised)",10.,0.,2000.);
  TH1D* h2 = new TH1D("h2",";#sqrt{s} - E_{|#eta|>8} (GeV);Counts (normalised)",10.,0.,13000.);
  TH1D* h1mb = new TH1D("h1mb",";ZDC Energy Sum (a.u.);Ratio to MB (normalised)",10.,0.,2000.);
  TH1D* h2mb = new TH1D("h2mb",";#sqrt{s} - E_{|#eta|>8} (GeV);Ratio to MB (normalised)",10.,0.,13000.);

  TLatex *xlabel = new TLatex();
	//xlabel->SetTextFont(42);
  xlabel-> SetNDC();
  xlabel-> SetTextColor(1);
  xlabel-> SetTextSize(0.05);
  xlabel-> SetTextAlign(22);
  xlabel-> SetTextAngle(0);

  TCanvas* c1 = new TCanvas("c1","",2000,800);
  c1->Divide(2);
  //
  c1->cd(1);
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetBottomMargin(0.15);
  c1->cd(1)->SetRightMargin(0.08);
  c1->cd(1)->SetTopMargin(0.08);
  c1->cd(1)->SetTicky();
  c1->cd(1)->SetTickx();
  c1->cd(1)->SetLogy();
  highZ->SetLineWidth(3);
  midZ->SetLineWidth(3);
  lowZ->SetLineWidth(3);
  highZ->SetLineStyle(7);
  midZ->SetLineStyle(7);
  lowZ->SetLineStyle(7);
  highZ->SetLineColor(kRed+1);
  midZ->SetLineColor(kGreen+2);
  lowZ->SetLineColor(kBlue+1);
  h1->SetStats(0);
  h1->GetYaxis()->SetRangeUser(0.0001,0.2);
  h1->GetXaxis()->SetTitleSize(0.045);
  h1->Draw();
  highZ->DrawNormalized("SAME");
  midZ->DrawNormalized("SAME");
  lowZ->DrawNormalized("SAME");
  xlabel-> DrawLatex(0.3, 0.83,Form("%s",prodAll[imc]));
  TLegend* l0 = new TLegend (0.52,0.7,0.83,0.85);
  l0->SetBorderSize(0);
  if (strncmp (sel,"1020",4) == 0){
    l0->SetHeader("SPD fixed [10,20]%");
    l0->AddEntry(highZ ,"V0M [0,5]","L");
    l0->AddEntry(midZ ,"V0M [20,30]","L");
    l0->AddEntry(lowZ ,"V0M [50,100]","L");
  }
  if (strncmp (sel,"4050",4) == 0){
    l0->SetHeader("SPD fixed [40,50]%");
    l0->AddEntry(highZ ,"V0M [0,20]","L");
    l0->AddEntry(midZ ,"V0M [40,50]","L");
    l0->AddEntry(lowZ ,"V0M [70,100]","L");
  }
  l0->SetTextSize(0.03);
  l0->SetTextFont(42);
  l0->Draw("SAME");
  //
  c1->cd(2);
  c1->cd(2)->SetLeftMargin(0.15);
  c1->cd(2)->SetBottomMargin(0.15);
  c1->cd(2)->SetRightMargin(0.08);
  c1->cd(2)->SetTopMargin(0.08);
  c1->cd(2)->SetTicky();
  c1->cd(2)->SetTickx();
  high->SetLineWidth(3);
  mid->SetLineWidth(3);
  low->SetLineWidth(3);
  high->SetLineStyle(7);
  mid->SetLineStyle(7);
  low->SetLineStyle(7);
  high->SetLineColor(kRed+1);
  mid->SetLineColor(kGreen+2);
  low->SetLineColor(kBlue+1);
  h2->SetStats(0);
  h2->GetYaxis()->SetRangeUser(0.00001,0.04);
  h2->Draw();
  h2->GetXaxis()->SetTitleSize(0.045);
  high->DrawNormalized("SAME");
  mid->DrawNormalized("SAME");
  low->DrawNormalized("SAME");
  TLegend* l1 = new TLegend (0.22,0.7,0.4,0.85);
  l1->SetBorderSize(0);
  if (strncmp (sel,"1020",4) == 0){
    l1->SetHeader("SPD fixed [10,20]%");
    l1->AddEntry(highZ ,"V0M [0,5]","L");
    l1->AddEntry(midZ ,"V0M [20,30]","L");
    l1->AddEntry(lowZ ,"V0M [50,100]","L");
  }
  if (strncmp (sel,"4050",4) == 0){
    l1->SetHeader("SPD fixed [40,50]%");
    l1->AddEntry(highZ ,"V0M [0,20]","L");
    l1->AddEntry(midZ ,"V0M [40,50]","L");
    l1->AddEntry(lowZ ,"V0M [70,100]","L");
  }
  l1->SetTextSize(0.03);
  l1->SetTextFont(42);
  l1->Draw("SAME");
  xlabel-> DrawLatex(0.7, 0.83,Form("%s",prodAll[imc]));
  //
  c1->SaveAs(Form("zdceffenergy_mc%i_esttype%i_%s.png",imc,type,sel));


  TCanvas *c2 = new TCanvas("c2","",2000,800);
  c2->Divide(2);
  //
  c2->cd(1);
  c2->cd(1)->SetLeftMargin(0.15);
  c2->cd(1)->SetBottomMargin(0.15);
  c2->cd(1)->SetRightMargin(0.08);
  c2->cd(1)->SetTopMargin(0.08);
  c2->cd(1)->SetTicky();
  c2->cd(1)->SetTickx();
  c2->cd(1)->SetLogy();
  TH1F* hZ = (TH1F*)highZ->Clone("highZ");
  TH1F* mZ = (TH1F*)midZ->Clone("midZ");
  TH1F* lZ = (TH1F*)lowZ->Clone("lowZ");
  hZ->Divide(mbZ);
  mZ->Divide(mbZ);
  lZ->Divide(mbZ);
  h1mb->SetStats(0);
  h1mb->GetYaxis()->SetRangeUser(0.0001,0.2);
  h1mb->GetXaxis()->SetTitleSize(0.045);
  h1mb->Draw();
  hZ->DrawNormalized("SAME");
  mZ->DrawNormalized("SAME");
  lZ->DrawNormalized("SAME");
  l0->Draw("SAME");
  xlabel-> DrawLatex(0.3, 0.83,Form("%s",prodAll[imc]));
  //
  c2->cd(2);
  c2->cd(2)->SetLeftMargin(0.15);
  c2->cd(2)->SetBottomMargin(0.15);
  c2->cd(2)->SetRightMargin(0.08);
  c2->cd(2)->SetTopMargin(0.08);
  c2->cd(2)->SetTicky();
  c2->cd(2)->SetTickx();
  TH1F* hh = (TH1F*)high->Clone("high");
  TH1F* m = (TH1F*)mid->Clone("mid");
  TH1F* l = (TH1F*)low->Clone("low");
  hh->Divide(mb);
  m->Divide(mb);
  l->Divide(mb);
  h2mb->SetStats(0);
  h2mb->GetYaxis()->SetRangeUser(0.00001,0.04);
  h2mb->GetXaxis()->SetTitleSize(0.045);
  h2mb->Draw();
  hh->DrawNormalized("SAME");
  m->DrawNormalized("SAME");
  l->DrawNormalized("SAME");
  l1->Draw("SAME");
  xlabel-> DrawLatex(0.7, 0.83,Form("%s",prodAll[imc]));
  //
  c2->SaveAs(Form("ratiotombzdceffenergy_mc%i_esttype%i_%s.png",imc,type,sel));



  TCanvas *c = new TCanvas("c","",2000,800);
  c->Divide(3,2);
  c->cd(1)->SetLogz();
  matrix->Draw("colz");
  c->cd(2)->SetLogy();
  highZ->DrawNormalized();
  midZ->DrawNormalized("SAME");
  lowZ->DrawNormalized("SAME");
  c->cd(3);
  high->DrawNormalized();
  mid->DrawNormalized("SAME");
  low->DrawNormalized("SAME");
  c->cd(5)->SetLogy();
  mbZ->Draw();
  c->cd(6);
  mb->Draw();
}
