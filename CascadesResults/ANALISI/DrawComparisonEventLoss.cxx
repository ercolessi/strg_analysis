{
    const int nbins = 10;
	Double_t centrality[nbins] = {0.,20.,30.,40.,50.,60.,70.,80.,90.,100.};
    TH1D* h = new TH1D("h",";(#sqrt{s} - ZDC) percentile %;#varepsilon_{event}^{twiki} - #varepsilon_{event}^{This}",nbins-1,centrality);
            Double_t emio[]  = {0.91,0.92,0.93,0.92,0.90,0.87,0.84,0.79,0.69};
    Double_t etwiki[nbins-1] = {0.91,0.92,0.92,0.90,0.87,0.86,0.84,0.81,0.71};

    for (int i = 1; i < nbins ; i++){
        h->SetBinContent(i, emio[i-1]-etwiki[i-1]);
    }

    new TCanvas;
    h->Draw();
    




}