#include<stdio.h>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TAxis.h"

const int nMCType = 4;
const char *filename[nMCType];
const char *mcType[nMCType];
const char *anchored[nMCType];

const int nV0cent = 11;
//int V0centr[nV0cent+1] = {0,10,40,70,100};
int V0centr[nV0cent+1] = {0,1,5,10,15,20,30,40,50,60,70,100};

const int nZDCcent = 9;
int ZDCcentr[nZDCcent+1] = {0,20,30,40,50,60,70,80,90,100};

const int nSPDcent = 11;
int SPDcentr[nSPDcent+1] = {0,1,5,10,15,20,30,40,50,60,70,100};

int *centr[3] = {V0centr, ZDCcentr, SPDcentr};
int ncentr[3] = {nV0cent, nZDCcent, nSPDcent};

bool isV0fixed=false;
int V0min=0;
int V0max=100;

bool isZDCfixed=false;
int ZDCmin=0;
int ZDCmax=100;

bool isSPDfixed=false;
int SPDmin=0;
int SPDmax=100;

float cCentr[3];

int differential = -1; // 0=V0, 1=ZDC, 2=SPD

TTree *t;

void init(const char *sel);
bool isSelected();
float GetExpression(const char *formula);

void analize(int iMC=0, const char *sel=0, TString X="SPDtracklets", TString Y="effEnergy",const char *tag="",bool normalizedM=false,bool normalizedX=false, const char *filter=0){
  printf("run in this way\n");
  printf("root 'analize.C(1,\"-SPD 0-10 -ZDC diff\",\"_tag\")'\n");
  printf("first number to select what MC production, _tag used for the name of output files\n");
  printf("-DET diff (DET = SPD, ZDC, V0) means to differentiate on that detector percentile\n");
  printf("-DET XX-YY (DET = SPD, ZDC, V0) means to require that detector percentile in that range\n\n");
  printf("root 'analize.C(1,\"-SPD 0-10 -ZDC diff\",\"SPDtracklets\",\"effEnergy\")' to change the variables you want to plot\n\n\n\n");


  if(iMC < 0 || iMC >= nMCType){
    printf("No MC type = %d (max allowed = %d)\n",iMC,nMCType-1);
    return;
  }

  init(sel);

  TFile *f = new TFile(filename[iMC]);
  TTree *tOr = (TTree *) f->Get("fTree");
  TFile *fout = new TFile(Form("out%d%s.root",iMC,tag),"RECREATE");
  if(filter) t = tOr->CopyTree(filter);
  else t =tOr;

  int nev = t->GetEntries();

  float x,y;

  int nPoints = 1;
  if(differential > -1){
    nPoints = ncentr[differential];
  }

  float *xVal = new float[nPoints];
  float *yVal = new float[nPoints];
  float *exVal = new float[nPoints];
  float *eyVal = new float[nPoints];
  float *mVal = new float[nPoints];
  int *xN = new int[nPoints];

  for(int j=0; j < nPoints; j++){
    xVal[j] = yVal[j] = mVal[j] = 0.;
    exVal[j] = eyVal[j] = 0.;
    xN[j] = 0;
  }

  for(int i=0; i < nev/100; i++){
    t->GetEvent(i);
    if(!isSelected()) continue;

    if(!(i % 1000000)) printf("%d/%d\n",i,nev);

    x = GetExpression(X.Data());
    y = GetExpression(Y.Data());

//    x = t->GetLeaf(X.Data())->GetValue();
//    y = t->GetLeaf(Y.Data())->GetValue();
//    if(X.Contains("effEnergy")) x += 13000;
//    if(Y.Contains("effEnergy")) y += 13000;

    int icentr = 0;
    if(differential > -1){
      while(icentr < ncentr[differential]-1 && cCentr[differential] > centr[differential][icentr+1]){
	icentr++;
      }
    }

    xVal[icentr] += x;
    yVal[icentr] += y;
    exVal[icentr] += x*x;
    eyVal[icentr] += y*y;
    mVal[icentr] += t->GetLeaf("nchEta")->GetValue();
    xN[icentr]++;
  }

  for(int j=0; j < nPoints; j++){
    if(xN[j] > 0){
      xVal[j] /= xN[j];
      yVal[j] /= xN[j];
      exVal[j] /= xN[j];
      eyVal[j] /= xN[j];
      mVal[j] /= xN[j];

      exVal[j] -= xVal[j]*xVal[j];
      eyVal[j] -= yVal[j]*yVal[j];

      exVal[j] = sqrt(exVal[j]/xN[j]);
      eyVal[j] = sqrt(eyVal[j]/xN[j]);

      if(normalizedM){
        yVal[j] /= mVal[j];
        eyVal[j] /= mVal[j];
      }
      if(normalizedX){
        yVal[j] /= xVal[j];
        eyVal[j] /= xVal[j];
      }
    }
    printf("%d) %f %f\n",j,xVal[j],yVal[j]);
  }

  TGraphErrors *g = new TGraphErrors(nPoints,xVal,yVal,exVal,eyVal);
  g->SetMarkerStyle(20);
  g->SetTitle("");
  g->SetMarkerSize(1.8);
  //
  TString trueY = Y;
  if (normalizedM) trueY = Form("#frac{%s}{dN_{ch}/d#eta}",Y.Data());
  if (normalizedX) trueY = Form("#frac{%s}{%s}",Y.Data(),X.Data());
  //
  g->GetXaxis()->SetTitle(X.Data());
  g->GetYaxis()->SetTitle(trueY.Data());
  //
  TCanvas* c = new TCanvas("c","",1000,800);
  c->SetRightMargin(0.09);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);
  //
  TString mctype = "Pythia8 Monash";
  TLatex *xlabel = new TLatex();
	xlabel->SetTextFont(42);
  xlabel-> SetTextSize(0.03);
  xlabel-> SetTextAlign(22);
  xlabel-> SetNDC();
  if (iMC==1) mctype = "Pythia8 Strg. Inj.";
  if (iMC==2) mctype = "Pythia6 Perugia";
  if (iMC==3) mctype = "Phojet";

  //Draw to canvas
  g->Draw("AP");
  xlabel->DrawLatex(0.27, 0.83,mctype.Data());

  //Write to file
  g->Write();
  fout->Close();
}

void init(const char *sel){
  filename[0] = "LHC20i2a-Pythia8_Monash2013.root";
  mcType[0] = "Pythia8_Monash2013";
  anchored[0] = "LHC18i";

  filename[1] = "LHC16d3-EPOS.root";
  mcType[1] = "EPOS_pp13_2015";
  anchored[1] = "LHC15f";

  filename[2] = "LHC17h7a-Pythia6_Perugia2011.root";
  mcType[2] = "pythia6_Perugia2011";
  anchored[2] = "LHC17j";

  filename[3] = "LHC17h7b-Phojet.root";
  mcType[3] = "phojet";
  anchored[3] = "LHC17j";

  if(sel){
    istringstream iss(sel);

    do {
      string subs;

      // Get the word from the istringstream
      iss >> subs;

      if(subs == "-V0"){
	int a,b;
	iss >> subs;
	if(subs == "diff"){
	  printf("V0 differential\n");
	  differential = 0;
	}
	else if(sscanf(subs.c_str(),"%d-%d",&a,&b) == 2){
	  printf("V0 Selection %d-%d\n",a,b);
	  V0min = a;
	  V0max = b;
          isV0fixed = true;
	}
      }

      if(subs == "-ZDC"){
	int a,b;
	iss >> subs;
	if(subs == "diff"){
	  printf("ZDC differential\n");
	  differential = 1;
	}
	else if(sscanf(subs.c_str(),"%d-%d",&a,&b) == 2){
	  printf("ZDC Selection %d-%d\n",a,b);
	  ZDCmin = a;
	  ZDCmax = b;
          isZDCfixed = true;
	}
      }

      if(subs == "-SPD"){
	int a,b;
	iss >> subs;
	if(subs == "diff"){
	  printf("SPD differential\n");
	  differential = 2;
	}
	else if(sscanf(subs.c_str(),"%d-%d",&a,&b) == 2){
	  printf("SPD Selection %d-%d\n",a,b);
	  SPDmin = a;
	  SPDmax = b;
          isSPDfixed=true;
	}
      }

    } while (iss);
  }
}

bool isSelected(){
  if(!t->GetLeaf("inelGT0")->GetValue()) return false;
  if(!t->GetLeaf("SPDtracklets")->GetValue()) return false;

  cCentr[0] = t->GetLeaf("v0mPerc")->GetValue();
  cCentr[1] = t->GetLeaf("zdcPerc")->GetValue();
  cCentr[2] = t->GetLeaf("multSPDcl")->GetValue();

 if(isV0fixed){
   if(cCentr[0] < V0min || cCentr[0] > V0max) return false;
 }

 if(isZDCfixed){
   if(cCentr[1] < ZDCmin || cCentr[1] > ZDCmax) return false;
 }

 if(isSPDfixed){
   if(cCentr[2] < SPDmin || cCentr[2] > SPDmax) return false;
 }

  return true;
}

float GetExpression(const char *formula){
    istringstream iss(formula);

    float result = 0,coef;

    do {
      string subs;

      // Get the word from the istringstream
      iss >> subs;

      if(subs[0] == '-'){
        sscanf(subs.c_str(),"%f",&coef);
        iss >> subs;
        if(subs.find("effEnergy") < 10000) coef *=t->GetLeaf(subs.c_str())->GetValue() + 13000;
        else coef *=t->GetLeaf(subs.c_str())->GetValue();
        result += coef;
      }
      else if(subs[0] == '+'){
        sscanf(subs.c_str(),"%f",&coef);
        iss >> subs;
        if(subs.find("effEnergy") < 10000) coef *=t->GetLeaf(subs.c_str())->GetValue() + 13000;
        else coef *=t->GetLeaf(subs.c_str())->GetValue();
        result += coef;
      }
      else if(subs[0] != '\0'){
        if(subs.find("effEnergy") < 10000) coef =t->GetLeaf(subs.c_str())->GetValue() + 13000;
        else coef =t->GetLeaf(subs.c_str())->GetValue();
        result += coef;
      }
    } while (iss);
  return result;
}

