#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <algorithm>
using namespace std;


void ConvertZDCtoEE( const char* filename = "MCLHC15g3c3.root", const char* period = "MC-15g3c3")
{  
    //Variables
    Float_t ZPCpp, ZNCpp, ZPApp, ZNApp;
    Int_t fRun;
    const Int_t binnumb = 1000;
    const Double_t range = 4.; 

    //Open file and get Tree 
    TFile* Read = new TFile (filename);
    TTree * T = (TTree *)Read->Get("PWGLF_StrVsMult_MC/fTreeEvent");
  
    //Set branches
    T->SetBranchAddress("fZPCpp",&ZPCpp);
    T->SetBranchAddress("fZPApp",&ZPApp);
    T->SetBranchAddress("fZNCpp",&ZNCpp);
    T->SetBranchAddress("fZNApp",&ZNApp);
    T->SetBranchAddress("fRun",&fRun);
    
    //Get runList
    std::vector<int> runList;
    for(Int_t i=0; i<T->GetEntries();i++)
    {
      T->GetEvent(i);//get Tree entries

      if (i == 0) runList.push_back(fRun); //first is added
      else {
        if(std::find(runList.begin(), runList.end(), fRun) != runList.end())
          continue;
        else runList.push_back(fRun);          
      }
    }
    //for (int i=0; i<runList.size(); i++) {cout << runList[i];};

    //Start converting ZDC to EE percentile
    cout << "\n--------------------------------------------------" <<
    "\nConverting ZDC information in effective energy ... \n"<<
    "---------------------------------------------------" << endl;
    
    //output file    
    TFile* outputfile = new TFile (Form("ExtractedZDCPercentile_%s.root",period), "recreate");

    //Usefull variables
    TH1F * hESumZDC[runList.size()+1];
    TH1F * hCumulative[runList.size()+1]; 
    Double_t IntTOT[runList.size()+1];
    TDirectoryFile *lDir[runList.size()];

    for (int j=0; j<runList.size()+1; j++){

      if (j!=runList.size()) {
        hESumZDC[j] = new TH1F(Form("hESumZDC_%i",runList[j]),Form("hSumEZDC_%i",runList[j]), binnumb ,0.,range);
        hCumulative[j] = new TH1F(*hESumZDC[j]);
        hCumulative[j]->SetName(Form("hCumulative_%i",runList[j]));
        hCumulative[j]->Reset();
      }
      else {
        hESumZDC[j] = new TH1F("hESumZDC_TOT","hSumEZDC_TOT", binnumb ,0.,range);
        hCumulative[j] = new TH1F(*hESumZDC[j]);
        hCumulative[j]->SetName("hCumulative_TOT");
        hCumulative[j]->Reset();
      }

      if (j==runList.size()) runList[j] = 0;
      cout << 
      "\n--------------------------------------------------" <<
      "\nProcessing run number =" << " run # = " << runList[j] << "\n" <<
      "----------------------------------------------------" << endl;

      //
      for(Int_t i=0; i<T->GetEntries();i++)//loop over Tree entries
      {
        T->GetEvent(i);//get Tree entries

        Float_t EZDC, LogE;
        if (i%1000000 == 0)
        cout<<"\nAnalizing event = "<<i+1<<" / "<< T->GetEntries() << endl;
      
        EZDC = ZPCpp + ZNCpp + ZPApp + ZNApp;
        if (EZDC != 0) LogE = TMath::Log10(TMath::Abs(EZDC)+1);
        else LogE = 0;
        if (j==runList.size()) hESumZDC[j]->Fill(LogE);
        else {
          if (fRun == runList[j]) hESumZDC[j]->Fill(LogE);
          else continue;
        }
      }

      //Fill Cumulative
      IntTOT[j] = hESumZDC[j]->Integral(1,binnumb);// == total number of events
      Double_t val = 0;
    
      for (Int_t k=1; k<binnumb; k++) //loop through the bins
      {
        val+= (hESumZDC[j]->GetBinContent(k))/IntTOT[j];
        hCumulative[j]->SetBinContent(k,val);
      }

      //Write
      if (j==runList.size()) {
        hESumZDC[j]->Write();
        hCumulative[j]->Write();  
      }
      else {
        lDir[j] = new TDirectoryFile(Form("%i",runList[j]),Form("%i",runList[j]));
        lDir[j]->cd();
        hESumZDC[j]->Write();
        hCumulative[j]->Write();  
        outputfile->cd();    
      }
    } 
 }
