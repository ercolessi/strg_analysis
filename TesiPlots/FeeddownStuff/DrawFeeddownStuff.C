void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 );
double PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin);

void DrawFeeddownStuff(){

  TLatex *tex = new TLatex(0.2, 0.8, "ALICE, pp #sqrt{s} = 13 TeV");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  TLatex *tex2 = new TLatex(0.2, 0.75, "This work");
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.04);
  tex2->SetLineWidth(2);
  TLatex *texh = new TLatex(0.2, 0.85, "ALICE, pp #sqrt{s} = 13 TeV");
  texh->SetNDC();
  texh->SetTextFont(42);
  texh->SetTextSize(0.04);
  texh->SetLineWidth(2);
  TLatex *tex2h = new TLatex(0.2, 0.8, "This work");
  tex2h->SetNDC();
  tex2h->SetTextFont(42);
  tex2h->SetTextSize(0.04);
  tex2h->SetLineWidth(2);

  TFile *fpRatio = TFile::Open("/home/fercoles/strg_analysis/Analisi2022/V0Analysis/results/backup/Results-Lambda-13TeV-SPDClusters_000_100_V0M_000_100_FDUseMCRatio.root");
  TFile *fapRatio = TFile::Open("/home/fercoles/strg_analysis/Analisi2022/V0Analysis/results/backup/Results-AntiLambda-13TeV-SPDClusters_000_100_V0M_000_100_FDUseMCRatio.root");
  TFile *fpDouble = TFile::Open("/home/fercoles/strg_analysis/Analisi2022/V0Analysis/results/backup/Results-Lambda-13TeV-SPDClusters_000_100_V0M_000_100_FDDoubleChargedXi.root");
  TFile *fapDouble = TFile::Open("/home/fercoles/strg_analysis/Analisi2022/V0Analysis/results/backup/Results-AntiLambda-13TeV-SPDClusters_000_100_V0M_000_100_FDDoubleChargedXi.root");

  TH1F *fracFDRatioP = (TH1F *)fpRatio->Get("lFeeddown/fHistFeeddownSubtraction");
  TH1F *fracFDRatioAP = (TH1F *)fapRatio->Get("lFeeddown/fHistFeeddownSubtraction");
  TH1F *fracFDDoubleP = (TH1F *)fpDouble->Get("lFeeddown/fHistFeeddownSubtraction");
  TH1F *fracFDDoubleAP = (TH1F *)fapDouble->Get("lFeeddown/fHistFeeddownSubtraction");

  TH1F *FDeffP = (TH1F *)fpRatio->Get("lFeeddown/f2dFeedDownEfficiencyGFCorrectedTOT");
  TH1F *FDeffAP = (TH1F *)fapRatio->Get("lFeeddown/f2dFeedDownEfficiencyGFCorrectedTOT");

  TCanvas *c = new TCanvas("c", "c", 900, 800);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  c->SetBottomMargin(0.15);
  c->SetTopMargin(0.1);
  c->SetTicks();
  FDeffP->SetStats(0);
  FDeffP->SetTitle("Feed-down efficiency #Lambda");
  FDeffP->GetXaxis()->SetTitle("#it{p}_{T} (#Lambda)");
  FDeffP->GetXaxis()->SetTitleSize(0.05);
  FDeffP->GetXaxis()->SetTitleOffset(1.1);
  FDeffP->GetYaxis()->SetTitle("#it{p}_{T} (#Xi^{-})");
  FDeffP->GetYaxis()->SetTitleSize(0.05);
  FDeffP->GetYaxis()->SetTitleOffset(1.1);
  FDeffP->Draw("colz");
  tex->Draw();
  tex2->Draw();

  TCanvas *c1 = new TCanvas("c1", "c1", 900, 800);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetTopMargin(0.1);
  c1->SetTicks();
  FDeffAP->SetStats(0);
  FDeffAP->SetTitle("Feed-down efficiency #bar{#Lambda}");
  FDeffAP->GetXaxis()->SetTitle("#it{p}_{T} (#bar{#Lambda})");
  FDeffAP->GetXaxis()->SetTitleSize(0.05);
  FDeffAP->GetXaxis()->SetTitleOffset(1.1);
  FDeffAP->GetYaxis()->SetTitle("#it{p}_{T} (#bar{#Xi}^{+})");
  FDeffAP->GetYaxis()->SetTitleSize(0.05);
  FDeffAP->GetYaxis()->SetTitleOffset(1.1);
  FDeffAP->Draw("colz");
  tex->Draw();
  tex2->Draw();

  c->SaveAs("FeedDownMatrix_Lambda.pdf");
  c1->SaveAs("FeedDownMatrix_AntiLambda.pdf");

  TH1D *hfd = new TH1D("hfd", ";;Fraction removed", 10, 0, 8);
  hfd->SetStats(0);
  hfd->GetYaxis()->SetRangeUser(0, 0.3);

  fracFDRatioP->SetMarkerStyle(kFullCircle);
  fracFDRatioP->SetMarkerSize(1.5);
  fracFDRatioP->SetMarkerColor(kBlack);
  fracFDRatioP->SetLineColor(kBlack);
  fracFDRatioAP->SetMarkerStyle(kFullCircle);
  fracFDRatioAP->SetMarkerSize(1.5);
  fracFDRatioAP->SetMarkerColor(kBlack);
  fracFDRatioAP->SetLineColor(kBlack);

  fracFDDoubleP->SetMarkerStyle(kOpenCircle);
  fracFDDoubleP->SetMarkerSize(1.5);
  fracFDDoubleP->SetMarkerColor(kBlack);
  fracFDDoubleP->SetLineColor(kBlack);
  fracFDDoubleAP->SetMarkerStyle(kOpenCircle);
  fracFDDoubleAP->SetMarkerSize(1.5);
  fracFDDoubleAP->SetMarkerColor(kBlack);
  fracFDDoubleAP->SetLineColor(kBlack);

  TLatex *xlabel = new TLatex();
  xlabel->SetTextFont(42);
  xlabel->SetNDC();
  xlabel->SetTextColor(1);
  xlabel->SetTextSize(0.07);
  xlabel->SetTextAlign(22);
  xlabel->SetTextAngle(0);

  TCanvas *c2 = new TCanvas(Form("c2"), "", 900, 800);
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.05);
  c2->SetBottomMargin(0.15);
  c2->SetTopMargin(0.05);
  c2->SetTicks();
  TLegend *leg = new TLegend(0.6, 0.25, 0.75, 0.45);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  hfd->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hfd->GetYaxis()->SetTitle("Fraction removed");
  hfd->GetXaxis()->SetTitleSize(0.05);
  hfd->GetYaxis()->SetTitleSize(0.05);
  hfd->GetXaxis()->SetTitleOffset(1.1);
  hfd->GetYaxis()->SetTitleOffset(1.1);
  hfd->Draw();
  fracFDRatioP->Draw("SAME EP");
  fracFDDoubleP->Draw("SAME EP");
  leg->AddEntry(fracFDRatioP, "MC Ratio", "LEP");
  leg->AddEntry(fracFDDoubleP, "Double Xi", "LEP");
  leg->SetTextSize(0.035);
  xlabel->DrawLatex(0.85, 0.75, Form("#Lambda"));
  leg->Draw("SAME");
  texh->Draw();
  tex2h->Draw();

  c2->SaveAs("FeedownLRatiovsDouble.pdf");

  TCanvas *c3 = new TCanvas(Form("c3"), "", 900, 800);
  c3->SetLeftMargin(0.15);
  c3->SetRightMargin(0.05);
  c3->SetBottomMargin(0.15);
  c3->SetTopMargin(0.05);
  c3->SetTicks();
  hfd->Draw();
  fracFDRatioAP->Draw("SAME EP");
  fracFDDoubleAP->Draw("SAME EP");
  leg->Draw("SAME");
  xlabel->DrawLatex(0.85, 0.75, Form("#bar{#Lambda}"));
  texh->Draw();
  tex2h->Draw();

  //
  c3->SaveAs("FeedownALRatiovsDouble.pdf");

  TFile *fp[100], *fap[100];
  TH1F *hp[100], *hap[100];
  Float_t percentileV0M_low[] = {0, 0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
  Float_t percentileV0M_high[] = {100, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
  Float_t percentileSPDCl_low[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Float_t percentileSPDCl_high[] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100};

  const int nbins = 11;
  Int_t color[] = {kBlack, kRed + 2, kRed, kOrange + 7, kYellow + 1, kSpring - 1, kGreen + 2, kTeal, kAzure + 7, kBlue, kBlue + 2};
  for (int i = 0; i < 11; i++)
  {
    fp[i] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/resultsFD/Results-Lambda-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_FDUseMCRatio.root", percentileSPDCl_low[i], percentileSPDCl_high[i], percentileV0M_low[i], percentileV0M_high[i]));
    hp[i] = (TH1F *)fp[i]->Get("lFeeddown/fHistFeeddownSubtraction");
    fap[i] = TFile::Open(Form("/home/fercoles/strg_analysis/22Sett2023/V0Analysis/resultsFD/Results-AntiLambda-13TeV-SPDClusters_%03.0f_%03.0f_V0M_%03.0f_%03.0f_FDUseMCRatio.root", percentileSPDCl_low[i], percentileSPDCl_high[i], percentileV0M_low[i], percentileV0M_high[i]));
    hap[i] = (TH1F *)fap[i]->Get("lFeeddown/fHistFeeddownSubtraction");
    hp[i]->SetLineColor(color[i]);
    hap[i]->SetLineColor(color[i]);
    hp[i]->SetLineWidth(2);
    hap[i]->SetLineWidth(2);
    hp[i]->SetMarkerStyle(kFullCircle);
    hap[i]->SetMarkerStyle(kFullCircle);
    hp[i]->SetMarkerSize(1.5);
    hap[i]->SetMarkerSize(1.5);
    hp[i]->SetMarkerColor(color[i]);
    hap[i]->SetMarkerColor(color[i]);
  }
  hp[0]->SetMarkerStyle(kFullSquare);
  hap[0]->SetMarkerStyle(kFullSquare);

  TCanvas *c4 = new TCanvas(Form("c4"), "", 1000, 900);
  c4->SetLeftMargin(0.15);
  c4->SetRightMargin(0.05);
  c4->SetBottomMargin(0.15);
  c4->SetTopMargin(0.05);
  c4->SetTicks();
  TLegend *lg = new TLegend(0.25, 0.2, 0.7, 0.4);
  lg->SetTextSize(0.027);
  lg->SetTextFont(42);
  lg->SetBorderSize(0);
  lg->SetNColumns(2);
  hfd->Draw();
  lg->AddEntry(hp[0], Form("VZERO %.0f-%.0f %s", percentileV0M_low[0], percentileV0M_high[0], "%"), "LPE");
  for (int i = 1; i < 11; i++)
  {
    lg->AddEntry(hp[i], Form("VZERO %.0f-%.0f %s", percentileV0M_low[i], percentileV0M_high[i], "%"), "LEP");
    hp[i]->Draw("SAME");
  }
  hp[0]->Draw("SAME");
  lg->Draw("SAME");
  xlabel->DrawLatex(0.85, 0.85, Form("#Lambda"));
  c4->SaveAs("AllFDL.pdf");
  texh->Draw();
  tex2h->Draw();

  TCanvas *c5 = new TCanvas(Form("c5"), "", 1000, 900);
  c5->SetLeftMargin(0.15);
  c5->SetRightMargin(0.05);
  c5->SetBottomMargin(0.15);
  c5->SetTopMargin(0.05);
  c5->SetTicks();
  hfd->Draw();
  for (int i = 1; i < 11; i++)
  {
    hap[i]->Draw("SAME");
  }
  hap[0]->Draw("SAME");
  lg->Draw("SAME");
  xlabel->DrawLatex(0.85, 0.85, Form("#bar{#Lambda}"));
  c5->SaveAs("AllFDAL.pdf");
  texh->Draw();
  tex2h->Draw();

  /*TCanvas* c4 = new TCanvas(Form("c4"),"",1200,1800);
  TH1F* hn = new TH1F("hn","Systematic on Feed-Down",10,0.4,8.);
  hn->SetStats(0);
  hn->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hn->GetYaxis()->SetTitle("Rel. Dev. on Yields (if > 1 #sigma_{RB})");
  hn->GetYaxis()->SetRangeUser(0.,.03);
  hn->Draw();

  TH1F* hRatioL = (TH1F*)fpRatio->Get("fHistPtLambda");
  TH1F* hRatioAL = (TH1F*)fapRatio->Get("fHistPtAntiLambda");
  TH1F* hDoubleL = (TH1F*)fpDouble->Get("fHistPtLambda");
  TH1F* hDoubleAL = (TH1F*)fapDouble->Get("fHistPtAntiLambda");

  for (int bin = 1; bin <=  hRatioL->GetNbinsX(); bin ++){
   hRatioL->SetBinContent(bin, hRatioL->GetBinContent(bin)+hRatioAL->GetBinContent(bin));
   hRatioL->SetBinError(bin, TMath::Sqrt(hRatioL->GetBinError(bin)*hRatioL->GetBinError(bin) +
                                             hRatioAL->GetBinError(bin)*hRatioAL->GetBinError(bin)));
   hDoubleL->SetBinContent(bin, hDoubleL->GetBinContent(bin)+hDoubleAL->GetBinContent(bin));
   hDoubleL->SetBinError(bin, TMath::Sqrt(hDoubleL->GetBinError(bin)*hDoubleL->GetBinError(bin) +
                                             hDoubleAL->GetBinError(bin)*hDoubleAL->GetBinError(bin)));
  }


  DivideAndComputeRogerBarlow(hDoubleL,hRatioL);
  TH1F* final = (TH1F*)hRatioL->Clone("hSystFD");
  final->Reset();
  for (int bin = 1; bin <=  hRatioL->GetNbinsX(); bin ++){
      final->SetBinContent(bin,PassRogerBarlowCriterion(1,hDoubleL,bin));
      final->SetBinError(bin,hDoubleL->GetBinError(bin));
  }
  final->SetLineWidth(2);
  final->Draw("SAME");
  xlabel-> DrawLatex(0.85, 0.75,Form("#Lambda + #bar{#Lambda}"));

  TFile* Write = new TFile ("Systematics/FDSyst.root", "RECREATE");
  final->Write();*/

}

//----------------------------------------------------------------------------------------------------
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 ){
  //Use Roger Barlow "sigma_{delta}" as errors for ratios
  Double_t lh1NBins = h1->GetNbinsX();
  Double_t lh2NBins = h2->GetNbinsX();

  if( lh1NBins != lh2NBins ){
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return;
  }

  Double_t lSigmaDelta[100];
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){
    //Computation of roger barlow sigma_{delta}
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(h1->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
    //Computation of relationship to h2 for plotting in ratio plot
    if ( h2->GetBinContent(i) > 1e-12 ){
      lSigmaDelta[i] /= h2->GetBinContent(i);
    }else{
      lSigmaDelta[i] = 0;
    }
  }
  //Regular Division
  h1 -> Divide( h2 );
  //Replace Erorrs
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){
    h1->SetBinError(i, lSigmaDelta[i]);
  }
}

//----------------------------------------------------------------------------------------------------
double PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin){
  double dev = h->GetBinContent(bin);
  double RBsigma = h->GetBinError(bin);

  if (dev>(nsigmas*RBsigma)) {return TMath::Abs(dev-1);}
  else {return 0.;}
}