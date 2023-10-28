// classes
enum classname
{
    kStandalone = 0,
    kHighMult,
    kLowMult,
    kHighZN,
    kLowZN
};

void PrintEventLoss(int lClassName = kHighMult, Bool_t DoMB = kFALSE)
{

    TFile* Read = 0x0;
    if (DoMB) {
        Read = TFile::Open("EventLoss-13TeV_INELgt0.root");
    } else{
        Read = TFile::Open(Form("EventLoss-13TeV_class%i.root", lClassName));
    }

    TH1F *heventloss = (TH1F *)Read->Get("EventLoss/hevtloss");

    for (int nmult = 1; nmult <= heventloss->GetNbinsX(); nmult++)
    {
        int binevtloss = -1;
        cout << "Sel: " << heventloss->GetXaxis()->GetBinLabel(nmult) << endl;
        cout << "correction: " << Form("%.3f",heventloss->GetBinContent(nmult)) << endl;
    }
}