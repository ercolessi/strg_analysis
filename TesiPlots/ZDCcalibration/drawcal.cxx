void canvasprep(TCanvas* c);
void histoprep(TH1* h);

void drawcal(){
    TFile* f = new TFile("calZDC.root");
    TProfile* zdcna = (TProfile*)f->Get("zdcna");
    TProfile* zdcnc = (TProfile*)f->Get("zdcnc");
    TProfile* zdcpa = (TProfile*)f->Get("zdcpa");
    TProfile* zdcpc = (TProfile*)f->Get("zdcpc");

    TCanvas* c = new TCanvas("c","c",800,600);
    canvasprep(c);
    histoprep(zdcna);
    zdcna->GetYaxis()->SetTitle("#LT ZNA #GT (a.u.)");
    zdcna->GetXaxis()->SetTitle("Run");
    TProfile* zdcnaclone = (TProfile*) zdcna->Clone("zdcnaclone");
    zdcnaclone->SetMarkerColor(kBlack);
    zdcnaclone->SetLineColor(kBlack);
    zdcnaclone->SetMarkerStyle(kOpenSquare);
    zdcnaclone->SetMarkerSize(1.6);
    zdcna->GetYaxis()->SetRangeUser(0.65,1.55);
    zdcna->Draw();
    zdcnaclone->Draw("same");
    TH1F *hleg = new TH1F("hleg","hleg",1,0,1);
    hleg->SetFillStyle(4050);
    hleg->SetLineColor(kBlack);
    hleg->SetFillColor(kRed);

    TLegend* leg = new TLegend(0.55,0.82,0.67,0.865);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.039);
    leg->AddEntry(hleg," calibration constants ","F");
    leg->Draw("same");

    TLatex *tex = new TLatex(0.2, 0.88, "ALICE, pp #sqrt{s} = 13 TeV");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.045);
    tex->SetLineWidth(2);
    TLatex *tex2 = new TLatex(0.2, 0.82, "This work");
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.045);
    tex2->SetLineWidth(2);
    tex->Draw();
    tex2->Draw();
    c->SaveAs("zdcna_calibxrun.pdf");

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    canvasprep(c1);
    histoprep(zdcnc);
    zdcnc->GetYaxis()->SetTitle("#LT ZNC #GT (a.u.)");
    zdcnc->GetXaxis()->SetTitle("Run");
    TProfile *zdcncclone = (TProfile *)zdcnc->Clone("zdcncclone");
    zdcncclone->SetMarkerColor(kBlack);
    zdcncclone->SetLineColor(kBlack);
    zdcncclone->SetMarkerStyle(kOpenSquare);
    zdcncclone->SetMarkerSize(1.6);
    zdcnc->GetYaxis()->SetRangeUser(0.55, 1.45);
    zdcnc->Draw();
    zdcncclone->Draw("same");
    leg->Draw("same");
    tex->Draw();
    tex2->Draw();
    c1->SaveAs("zdcnc_calibxrun.pdf");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    canvasprep(c2);
    histoprep(zdcpa);
    zdcpa->GetYaxis()->SetTitle("#LT ZPA #GT (a.u.)");
    zdcpa->GetXaxis()->SetTitle("Run");
    TProfile *zdcpaclone = (TProfile *)zdcpa->Clone("zdcpaclone");
    zdcpaclone->SetMarkerColor(kBlack);
    zdcpaclone->SetLineColor(kBlack);
    zdcpaclone->SetMarkerStyle(kOpenSquare);
    zdcpaclone->SetMarkerSize(1.6);
    zdcpa->GetYaxis()->SetRangeUser(0.7, 1.4);
    zdcpa->Draw();
    zdcpaclone->Draw("same");
    leg->Draw("same");
    tex->Draw();
    tex2->Draw();
    c2->SaveAs("zdcpa_calibxrun.pdf");

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    canvasprep(c3);
    histoprep(zdcpc);
    zdcpc->GetYaxis()->SetTitle("#LT ZPC #GT (a.u.)");
    zdcpc->GetXaxis()->SetTitle("Run");
    TProfile *zdcpcclone = (TProfile *)zdcpc->Clone("zdcpcclone");
    zdcpcclone->SetMarkerColor(kBlack);
    zdcpcclone->SetLineColor(kBlack);
    zdcpcclone->SetMarkerStyle(kOpenSquare);
    zdcpcclone->SetMarkerSize(1.6);
    zdcpc->GetYaxis()->SetRangeUser(0.7, 1.8);
    zdcpc->Draw();
    zdcpcclone->Draw("same");
    leg->Draw("same");
    tex->Draw();
    tex2->Draw();
    c2->SaveAs("zdcpc_calibxrun.pdf");
}

void canvasprep(TCanvas* c){
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.1);
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.15);
    c->SetTicks();
}

void histoprep(TH1* h){
    h->SetStats(0);
    h->SetMarkerStyle(kFullSquare);
    h->SetMarkerSize(1.5);
    h->SetMarkerColor(kRed);
    h->SetLineColor(kRed);
    h->SetTitle("");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetLabelSize(0.05);
}