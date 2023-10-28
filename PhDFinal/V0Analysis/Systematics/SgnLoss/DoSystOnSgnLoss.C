void signalloss(TString part = "Xi", TString sel = "SPD", TString var = "V0M",float sel1 = 0., float sel2 = 100., Bool_t lDoMB = kFALSE, Int_t* perc =0x0 , Long_t nbins = 0, Int_t* col = 0x0);
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 );
double PassRogerBarlowCriterion(int nsigmas, TH1F* h, int bin);

void DoSystOnSgnLoss(TString part = "K0Short"){

  Int_t c0[] = {kRed+1, kRed-4, kOrange+7, kYellow+1, kSpring-7, kGreen+2, kAzure+8, kBlue-4, kBlue+3};
  Int_t c7[] = {kRed+1,  kOrange+7, kYellow+1, kSpring-7, kAzure+8, kBlue-4, kBlue+3};
  Int_t c8[] = {kRed + 1, kOrange + 7, kYellow + 1, kSpring - 7, kGreen + 2, kAzure + 8, kBlue - 4, kBlue + 3};
  Int_t c10[] =  {kRed+1, kRed-4, kOrange+7, kOrange-3, kYellow+1, kSpring-7, kGreen+2, kAzure+8, kBlue-4, kBlue+3};
  Int_t c4[] = {kRed+1, kYellow+1, kAzure+8, kBlue+3};

  Int_t p1[] = {0, 100};
  Long_t n1 = sizeof(p1)/sizeof(Int_t) - 1;
  signalloss(part,"SPDClusters","V0M",0.,100.,kTRUE,p1,n1,c10);

  Int_t p2[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
  Long_t n2 = sizeof(p2) / sizeof(Int_t) - 1;
  signalloss(part, "SPDClusters", "V0M", 0., 100., kFALSE, p2, n2, c10);

  Int_t p3[] = {0,5,10,20,30,40,50,100};
  Long_t n3 = sizeof(p3)/sizeof(Int_t) - 1;
  //
  signalloss(part, "SPDClusters", "V0M", 10., 20., kFALSE, p3, n3, c7);

  Int_t p4[] = {0,20,30,40,50,60,70,100};
  Long_t n4 = sizeof(p4)/sizeof(Int_t) - 1;
  //
  signalloss(part, "SPDClusters", "V0M", 40., 50., kFALSE, p4, n4, c7);

  Int_t p5[] = {0, 5, 10, 20, 30, 40, 50, 100};
  Long_t n5 = sizeof(p5) / sizeof(Int_t) - 1;
  //
  signalloss(part, "V0M", "SPDClusters", 10., 20., kFALSE, p5, n5, c7);

  Int_t p6[] = {0, 10, 20, 30, 40, 50, 60, 70, 100};
  Long_t n6 = sizeof(p6) / sizeof(Int_t) - 1;
  //
  signalloss(part, "V0M", "SPDClusters", 40., 50., kFALSE, p6, n6, c8);
  }

  void signalloss(TString part = "Xi", TString sel = "SPDcl", TString var = "V0M", float sel1 = 0., float sel2 = 100., Bool_t lDoMB = kFALSE, Int_t *perc = 0x0, Long_t nbins = 0, Int_t *col)
  {

    TFile *file = 0x0;
    if (lDoMB){
      file = new TFile(Form("SystNorm-%s-13TeV_INELgt0.root", part.Data()), "READ");
    } else {
      file = new TFile(Form("SystNorm-%s-13TeV_%s_Fixed%sin%03.0f_%03.0f.root", part.Data(), var.Data(), sel.Data(), sel1, sel2), "READ");
    }

  TH1F *hsgnloss_def[nbins], *hsgnloss_var[nbins];
  TH1F *hRatio[nbins], *hRatioClone[nbins];
  for (int i = 0; i < nbins; i++)
  {
    hsgnloss_def[i] = (TH1F *)file->Get(Form("hsgnloss_%i-%i_MC%i", (int)perc[i], (int)perc[i + 1], 0));
    hsgnloss_var[i] = (TH1F *)file->Get(Form("hsgnloss_%i-%i_MC%i", (int)perc[i], (int)perc[i + 1], 1));

    hRatio[i] = (TH1F *)hsgnloss_var[i]->Clone(Form("hRatio%i", i));
    hRatioClone[i] = (TH1F *)hsgnloss_var[i]->Clone(Form("hRatioClone%i", i));
    hRatioClone[i]->Reset();

    hRatioClone[i]->GetYaxis()->SetTitle("#frac{  #varepsilon_{part}^{Pythia6}}{#varepsilon_{part}^{Pythia8} }");
    hRatioClone[i]->SetName(Form("var_%i-%i", (int)perc[i], (int)perc[i + 1]));
    DivideAndComputeRogerBarlow(hRatio[i], hsgnloss_def[i]);
    for (int bin = 1; bin <= hRatio[i]->GetNbinsX(); bin++)
    {
      double dev = PassRogerBarlowCriterion(1, hRatio[i], bin);
      hRatioClone[i]->SetBinContent(bin, TMath::Abs(dev - 1));
      hRatioClone[i]->SetBinError(bin, .0000000001);
    }
    hRatioClone[i]->SetLineWidth(2);
    hRatioClone[i]->SetLineColor(col[i]);
  }

  TFile *Write;
  if (lDoMB) {
    Write = new TFile(Form("SgnLossSystematics-%s-13TeV_INELgt0.root", part.Data()),"RECREATE");
  } else {
    Write = new TFile(Form("SgnLossSystematics-%s_%s_Fixed%sin%03.0f_%03.0f.root", part.Data(), var.Data(), sel.Data(), sel1, sel2), "RECREATE");
  }

  for (int i = 0; i < nbins; i++)
  {
    hRatioClone[i]->Write();
  }
  Write->cd();
  Write->Close();

  TLatex *xlabel = new TLatex();
  xlabel->SetTextFont(42);
  xlabel->SetNDC();
  xlabel->SetTextColor(1);
  xlabel->SetTextSize(0.07);
  xlabel->SetTextAlign(22);
  xlabel->SetTextAngle(0);

  TCanvas *c = new TCanvas("c", "", 1000, 800);
  c->SetLeftMargin(0.22);
  TH1D *h = new TH1D("h", "", 10, 0.4, 6.5);
  h->GetYaxis()->SetTitle("#frac{ | #varepsilon_{part}^{Pythia8} - #varepsilon_{part}^{Pythia6} |}{#varepsilon_{part}^{Pythia8}}");
  h->GetYaxis()->SetTitleOffset(2.);
  h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h->GetYaxis()->SetRangeUser(0., .1);
  if (sel.Contains("V0M") && sel1 == 40.)
    h->GetYaxis()->SetRangeUser(0., .2);
  h->Draw();
  h->SetStats(0);

  TLegend *lp = new TLegend(0.57, 0.65, 0.87, 0.87);
  lp->SetTextSize(0.021);
  lp->SetBorderSize(0);
  if (sel1 != 0 && sel2 != 100)
    lp->SetHeader(Form("%s%s selection fixed [%.0f-%.0f], with:", sel.Data(), "%", sel1, sel2));
  for (int i = 0; i < nbins; i++)
  {
    hRatioClone[i]->Draw("SAME HIST");
    lp->AddEntry(hRatioClone[i], Form("%s [%i-%i] %s", var.Data(), perc[i], perc[i + 1], "%"), "L");
  }
  lp->Draw("SAME");
  if (part.Contains("K0Short"))
    xlabel->DrawLatex(0.3, 0.8, Form("K^{0}_{S}", part.Data()));
  else
    xlabel->DrawLatex(0.3, 0.8, Form("#%s", part.Data()));

  c->SaveAs(Form("images/SgnLossSyst-%s_%s_Fixed%sin%03.0f_%03.0f.png", part.Data(), var.Data(), sel.Data(), sel1, sel2));
  c->SaveAs(Form("images/SgnLossSyst-%s_%s_Fixed%sin%03.0f_%03.0f.root", part.Data(), var.Data(), sel.Data(), sel1, sel2));
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

    if( TMath::Abs(dev-1) > nsigmas*RBsigma ) {return dev;}
    else {return 1.;}

}
