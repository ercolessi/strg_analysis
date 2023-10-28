void preppad(TPad *pad1, TPad *pad2);
void prepcanvas(TCanvas *c);
void prepgraph(TGraphErrors *g, int marker, float size, int color);

void NMPI(TString pythia = "Monash", TString var = "NMPI")
{ //var = "AvPt"
    TFile *file = TFile::Open(Form("Results_%s_SPDClustersV0M.root",pythia.Data()));

    TString ddname[] = {"Standalone", "kHighMult", "kLowMult", "kHighZN", "kLowZN"};
    const int nsel = sizeof(ddname) / sizeof(TString);

    TGraphErrors *gNMPILE[nsel], *gNMPINch[nsel];
    for (int i = 0; i < nsel; i++)
    {
        if (var.Contains("NMPI")){
            gNMPILE[i] = (TGraphErrors *)file->Get(Form("%s/LEVsNMPI", ddname[i].Data()));
            gNMPINch[i] = (TGraphErrors *)file->Get(Form("%s/NchVsNMPI", ddname[i].Data()));
        } else if (var.Contains("AvPt")){
            gNMPILE[i] = (TGraphErrors *)file->Get(Form("%s/MeanPt/AvPt_VsLeadE_PiMinus", ddname[i].Data()));
            gNMPINch[i] = (TGraphErrors *)file->Get(Form("%s/MeanPt/AvPt_VsNch_PiMinus", ddname[i].Data()));
        }
    }

    prepgraph(gNMPILE[0], kFullDiamond, 2.3, kBlack);
    prepgraph(gNMPILE[1], kFullCircle, 1.8, kRed);
    prepgraph(gNMPILE[2], kFullCircle, 1.8, kBlue);
    prepgraph(gNMPILE[3], kFullSquare, 1.8, kViolet);
    prepgraph(gNMPILE[4], kFullSquare, 1.8, kGreen + 1);

    prepgraph(gNMPINch[0], kFullDiamond, 2.3, kBlack);
    prepgraph(gNMPINch[1], kFullCircle, 1.8, kRed);
    prepgraph(gNMPINch[2], kFullCircle, 1.8, kBlue);
    prepgraph(gNMPINch[3], kFullSquare, 1.8, kViolet);
    prepgraph(gNMPINch[4], kFullSquare, 1.8, kGreen + 1);

    TH1F *hnch = new TH1F("h", "", 100, 0, 32);
    hnch->SetStats(0);
    hnch->GetYaxis()->SetTitle("#LT N_{MPI} #GT");
    hnch->GetXaxis()->SetTitle("#LT n_{ch} #GT (|#eta|<0.5)");
    hnch->GetYaxis()->SetTitleSize(0.05);
    hnch->GetXaxis()->SetTitleSize(0.05);
    hnch->GetYaxis()->SetTitleOffset(1.2);
    hnch->GetXaxis()->SetTitleOffset(1.2);
    if (pythia.Contains("Ropes")) hnch->GetYaxis()->SetRangeUser(0., 21.9);
    if (pythia.Contains("Monash")) hnch->GetYaxis()->SetRangeUser(0., 18.9);


    TH1F *hzdc = new TH1F("hzdc", "", 100, 10, 4000);
    hzdc->SetStats(0);
    hzdc->GetYaxis()->SetTitle("#LT N_{MPI} #GT");
    hzdc->GetXaxis()->SetTitle("#LT Neutron E_{leading} #GT (|#eta|>8)");
    hzdc->GetYaxis()->SetTitleSize(0.05);
    hzdc->GetXaxis()->SetTitleSize(0.05);
    hzdc->GetYaxis()->SetTitleOffset(1.2);
    hzdc->GetXaxis()->SetTitleOffset(1.2);
    if (pythia.Contains("Ropes")) hzdc->GetYaxis()->SetRangeUser(0., 21.9);
    if (pythia.Contains("Monash")) hzdc->GetYaxis()->SetRangeUser(0., 18.9);

    TCanvas *c = new TCanvas("c", "", 1800, 800);
    TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 0.54, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0.54, 0, 1., 1);
    prepcanvas(c);
    preppad(pad1, pad2);
    pad1->Draw();
    pad2->Draw();
    //
    pad1->cd();
    hnch->Draw();

    for (int i = 0; i < nsel; i++)
    {
        gNMPINch[i]->Draw("LEP same");
    }

    TLatex *label = new TLatex();
    label->SetTextFont(42);
    label->SetNDC();
    label->SetTextSize(0.04);
    label->SetTextAlign(22);
    label->SetTextAngle(0);
    label->DrawLatex(0.4, 0.88, "ALICE, pp #sqrt{#it{s}} = 13 TeV");

    //
    pad2->cd();
    hzdc->Draw();

    for (int i = 0; i < nsel; i++)
    {
        gNMPILE[i]->Draw("LEP same");
    }

    TLegend *leg = new TLegend(0.6, 0.64, 0.86, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->SetHeader(Form("Pythia %s", pythia.Data()));
    leg->AddEntry(gNMPINch[0], "V0M Standalone ", "p");
    leg->AddEntry(gNMPINch[2], "Low Multiplicity ", "p");
    leg->AddEntry(gNMPINch[1], "High Multiplicity ", "p");
    leg->AddEntry(gNMPINch[4], "Low ZN ", "p");
    leg->AddEntry(gNMPINch[3], "High ZN ", "p");
    leg->Draw("same");


    c->SaveAs(Form("images/%s_Pythia%s.eps", var.Data(), pythia.Data()));

}

void prepgraph(TGraphErrors* g, int marker, float size, int color){
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(size);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineStyle(1);
    g->SetFillStyle(0);
}

void prepcanvas(TCanvas* c){
    c->SetLeftMargin(0.2);
    c->SetBottomMargin(0.18);
    c->SetRightMargin(0.1);
    c->SetTopMargin(0.1);
    c->SetTicky();
    c->SetTickx();
}

void preppad(TPad *pad1, TPad *pad2)
{
    pad1->SetBorderMode(0);
    pad1->SetTopMargin(0.05);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.15);
    pad1->SetRightMargin(0.0);
    pad2->SetRightMargin(0.05);
    pad1->SetLeftMargin(0.2);
    pad2->SetLeftMargin(0.0);
    pad1->SetBottomMargin(0.15);
    pad2->SetBorderMode(0);
    pad1->SetFillColor(kWhite);
    pad1->SetTicky();
    pad1->SetTickx();
    pad2->SetFillColor(kWhite);
    pad2->SetTicky();
    pad2->SetTickx();
}