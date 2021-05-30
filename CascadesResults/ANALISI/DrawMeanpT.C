void DrawMeanpT(){

    TFile* fileZDCStandalone = TFile::Open("MeanpT-Xi-ZDC-V0M_000_100_ZDC_000_100.root");
    TFile* fileV0MStandalone = TFile::Open("MeanpT-Xi-V0M-V0M_000_100_ZDC_000_100.root ");
    TFile* fileZDClowmult    = TFile::Open("MeanpT-Xi-ZDC-V0M_070_100_ZDC_000_100.root ");
    TFile* fileZDChighmult   = TFile::Open("MeanpT-Xi-ZDC-V0M_000_030_ZDC_000_100.root ");
    TFile* fileV0Mlowee      = TFile::Open("MeanpT-Xi-V0M-V0M_000_100_ZDC_070_100.root ");
    TFile* fileV0Mhighee     = TFile::Open("MeanpT-Xi-V0M-V0M_000_100_ZDC_000_030.root ");


    TH1D* hZDCStandalone = (TH1D*) fileZDCStandalone->Get("MeanpTvsNch_Stat");
    hZDCStandalone->SetLineColor(kRed+1);
    hZDCStandalone->SetMarkerColor(kRed+1);
    hZDCStandalone->SetMarkerStyle(kFullCircle);
    hZDCStandalone->SetMarkerSize(2.3);
    hZDCStandalone->SetStats(0);

    TH1D* hV0MStandalone = (TH1D*) fileV0MStandalone->Get("MeanpTvsNch_Stat");
    hV0MStandalone->SetLineColor(kBlue+1);
    hV0MStandalone->SetMarkerColor(kBlue+1);
    hV0MStandalone->SetMarkerStyle(kFullSquare);
    hV0MStandalone->SetMarkerSize(2.3);
    hV0MStandalone->SetStats(0);


    TCanvas* stand = new TCanvas("stand","",1200,1000);

    stand->SetRightMargin(0.05);
    stand->SetTopMargin(0.05);
    stand->SetLeftMargin(0.15);
    stand->SetBottomMargin(0.15);

    TH1D* h = new TH1D("h","",10,0,30);
    h->SetStats(0);
    h->GetYaxis()->SetRangeUser(0.8,1.7);
    h->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/#it{c})");
    h->GetYaxis()->SetTitleOffset(1.);
    h->GetYaxis()->SetTitleSize(.05);
    h->GetXaxis()->SetTitle("#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}");
    h->GetXaxis()->SetTitleOffset(1.);
    h->GetXaxis()->SetTitleSize(.05);
    h->Draw();
    hV0MStandalone->Draw("EP SAME");
    hZDCStandalone->Draw("EP SAME");

    TLegend* l = new TLegend (0.20,0.75,0.49,0.89);
    l->SetBorderSize(0);
    l->AddEntry( hV0MStandalone,"V0M","P");
    l->AddEntry( hZDCStandalone,"(#sqrt{#it{s}} - ZDC)","P");
    l->SetTextSize(0.037);
    l->SetTextFont(42);
    l->Draw("SAME");


    TH1D* hZDClowmult = (TH1D*) fileZDClowmult->Get("MeanpTvsPerc_Stat");
     hZDClowmult->SetLineColor(kRed+1);
    hZDClowmult->SetMarkerColor(kRed+1);
    hZDClowmult->SetMarkerStyle(kOpenCircle);
    hZDClowmult->SetMarkerSize(2.3);
    hZDClowmult->SetStats(0);

    TH1D* hZDChighmult = (TH1D*) fileZDChighmult->Get("MeanpTvsPerc_Stat");
    hZDChighmult ->SetLineColor(kRed+1);
    hZDChighmult->SetMarkerColor(kRed+1);
    hZDChighmult->SetMarkerStyle(kFullCircle);
    hZDChighmult->SetMarkerSize(2.3);
    hZDChighmult->SetStats(0);

    TH1D* hZDClowmultnch = (TH1D*) fileZDClowmult->Get("MeanpTvsNch_Stat");
    hZDClowmultnch->SetLineColor(kRed+1);
    hZDClowmultnch->SetMarkerColor(kRed+1);
    hZDClowmultnch->SetMarkerStyle(kOpenCircle);
    hZDClowmultnch->SetMarkerSize(2.3);
    hZDClowmultnch->SetStats(0);

    TH1D* hZDChighmultnch = (TH1D*) fileZDChighmult->Get("MeanpTvsNch_Stat");
    hZDChighmultnch ->SetLineColor(kRed+1);
    hZDChighmultnch->SetMarkerColor(kRed+1);
    hZDChighmultnch->SetMarkerStyle(kFullCircle);
    hZDChighmultnch->SetMarkerSize(2.3);
    hZDChighmultnch->SetStats(0);

    TH1D* hV0Mlowee = (TH1D*) fileV0Mlowee->Get("MeanpTvsNch_Stat");
    hV0Mlowee ->SetLineColor(kBlue+1);
    hV0Mlowee->SetMarkerColor(kBlue+1);
    hV0Mlowee->SetMarkerStyle(kOpenSquare);
    hV0Mlowee->SetMarkerSize(2.3);
    hV0Mlowee->SetStats(0);

    TH1D* hV0Mhighee = (TH1D*) fileV0Mhighee->Get("MeanpTvsNch_Stat");
    hV0Mhighee ->SetLineColor(kBlue+1);
    hV0Mhighee->SetMarkerColor(kBlue+1);
    hV0Mhighee->SetMarkerStyle(kFullSquare);
    hV0Mhighee->SetMarkerSize(2.3);
    hV0Mhighee->SetStats(0);

    TCanvas* sel1 = new TCanvas("sel1","",1200,1000);
    sel1->SetRightMargin(0.05);
    sel1->SetTopMargin(0.05);
    sel1->SetLeftMargin(0.15);
    sel1->SetBottomMargin(0.15);

    TH1D* h1 = new TH1D("h","",10,0,100);
    h1->SetStats(0);
    h1->GetYaxis()->SetRangeUser(0.8,1.7);
    h1->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/#it{c})");
    h1->GetYaxis()->SetTitleOffset(1.);
    h1->GetYaxis()->SetTitleSize(.05);
    h1->GetXaxis()->SetTitle("(#sqrt{#it{s}} - ZDC) percentile (%)");
    h1->GetXaxis()->SetTitleOffset(1.);
    h1->GetXaxis()->SetTitleSize(.05);
    //
    h1->Draw();
    hZDClowmult->Draw("EP SAME");
    hZDChighmult->Draw("EP SAME");

    TLegend* ll = new TLegend (0.20,0.75,0.49,0.89);
    ll->SetBorderSize(0);
   // ll->SetHeader("(#sqrt{#it{s}} - ZDC) classification");
    ll->AddEntry( hZDChighmult,"V0M 0-30%, high mult. ","P");
    ll->AddEntry( hZDClowmult,"V0M 70-100%, low mult.","P");
    ll->SetTextSize(0.037);
    ll->SetTextFont(42);
    ll->Draw("SAME");

   


    TCanvas* sel2 = new TCanvas("sel2","",1200,1000);
    sel2->SetRightMargin(0.05);
    sel2->SetTopMargin(0.05);
    sel2->SetLeftMargin(0.15);
    sel2->SetBottomMargin(0.15);

    h->Draw();
    hV0Mlowee->Draw("EP SAME");
    hV0Mhighee->Draw("EP SAME");

    TCanvas* sel3 = new TCanvas("sel2","",1200,1000);
    sel3->SetRightMargin(0.05);
    sel3->SetTopMargin(0.05);
    sel3->SetLeftMargin(0.15);
    sel3->SetBottomMargin(0.15);
    h->Draw();
    hZDClowmultnch->Draw("EP SAME");
    hZDChighmultnch->Draw("EP SAME");
    hV0Mlowee->Draw("EP SAME");
    hV0Mhighee->Draw("EP SAME");






    
}