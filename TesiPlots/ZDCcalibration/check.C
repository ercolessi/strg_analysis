void fill(const char* filename, int count);
void setcal();

float refZDCN = 210; 
float refZDCP = 65; 

TProfile2D *hMult[20], *hZDC[20];
TGraph *hScan[20];
TH1F *hMultRatio[20];
TH1F *hZDCRatio[20];

TProfile *hZDCNC[20];
TProfile *hZDCNA[20];
TProfile *hZDCPC[20];
TProfile *hZDCPA[20];

TProfile *hZDCNpro;
TProfile *hZDCPpro;

TProfile *hMultSPD[20];

TH2F *hSelection[20];

TProfile *hRunMult;
TProfile *hRunZDC;
TH2F *hRunZDCdist;

TH2F *hRunZDCNC;
TH2F *hRunZDCNA;
TH2F *hRunZDCPC;
TH2F *hRunZDCPA;

TH1F *hPeriodMult[20];
TH1F *hPeriodZDC[20];

// adjusting trackletes calibrations per period
float multScale[] = {0.760200, 1.000000, 0.935682};

// ZDC calibrations per period
float cal1[] = {1.210442, 1.280468, 1.192709, 0.940499};
float cal2[] = {1.000000, 1.069985, 1.000000, 0.895500};
float cal3[] = {1.018118, 1.061667, 1.166989, 0.913540};

float *calZDC[] = {cal1, cal2, cal3};

float adjZDC[] = {1., 1., 1.};

float calZNC[300000];
float calZNA[300000];
float calZPC[300000];
float calZPA[300000];

void check(bool printcalib=false){
  setcal();
  
  hRunMult = new TProfile("hRunMult","; ; < SPD Tracklets >",1,0,1);
  hRunZDC = new TProfile("hRunZDC","; ; < ZDC Sum >",1,0,1);
  hRunZDCdist = new TH2F("hRunZDCdist","; ; ZDC Sum",1,0,1,100,0,2000);
  hRunZDCNC = new TH2F("hRunZDCNC","; ; ZDCN-C",1,0,1,100,0,1000);
  hRunZDCNA = new TH2F("hRunZDCNA","; ; ZDCN-A",1,0,1,100,0,1000);
  hRunZDCPC = new TH2F("hRunZDCPC","; ; ZDCP-C",1,0,1,100,0,1000);
  hRunZDCPA = new TH2F("hRunZDCPA","; ; ZDCP-A",1,0,1,100,0,1000);

  hZDCNpro = new TProfile("hZDCNpro",";<ZDC sum> (a.u); <ZDC-N sum> (a.u.)",50,0,100);
  hZDCPpro = new TProfile("hZDCPpro",";<ZDC sum> (a.u); <ZDC-N sum> (a.u.)",50,0,100);

  for(int i=0; i < 100;i++){
    hRunZDCdist->Fill(Form("RUNS"), Form("%d",i*20), 0.);
    hRunZDCNC->Fill(Form("RUNS"), Form("%d",i*10), 0.);
    hRunZDCNA->Fill(Form("RUNS"), Form("%d",i*10), 0.);
    hRunZDCPC->Fill(Form("RUNS"), Form("%d",i*10), 0.);
    hRunZDCPA->Fill(Form("RUNS"), Form("%d",i*10), 0.);
  }
  
  system("ls LHC*.root >lista");
  
  FILE *f = fopen("lista","r");
  char filename[100];
  int nprod = 0;
  float x[100],y[100];

  int marker[] = {20,21,22};
  int color[] = {1, 2, 4};
  
  while(fscanf(f,"%s",filename) == 1){
    hMult[nprod] = new TProfile2D(Form("h%dMult",nprod),Form("%s; SPD (%); VZERO (%);< N tracklets >",filename),10,0,100,10,0,100);
    hZDC[nprod] = new TProfile2D(Form("h%dZDC",nprod),Form("%s; SPD (%); VZERO (%);< ZDC Sum > (a.u.)",filename),10,0,100,10,0,100);

    hMultRatio[nprod] = new TH1F(Form("h%dMultRatio",nprod),Form("%s;Ratio vs LHC15f;",filename),100,0,2);
    hZDCRatio[nprod] = new TH1F(Form("h%dZDCRatio",nprod),Form("%s;Ratio vs LHC15f;",filename),100,0,2);

    hZDCNC[nprod] = new TProfile(Form("h%dZDCNC",nprod),Form("%s;VZERO (%);<ZDCN-C>",filename),10,0,100);
    hZDCNA[nprod] = new TProfile(Form("h%dZDCNA",nprod),Form("%s;VZERO (%);<ZDCN-A>",filename),10,0,100);
    hZDCPC[nprod] = new TProfile(Form("h%dZDCPC",nprod),Form("%s;VZERO (%);<ZDCP-C>",filename),10,0,100);
    hZDCPA[nprod] = new TProfile(Form("h%dZDCPA",nprod),Form("%s;VZERO (%);<ZDCP-A>",filename),10,0,100);

    hMultSPD[nprod] = new TProfile(Form("h%dMultSPD",nprod),Form("%s;SPD (%);<N SPD tracklets>",filename),10,0,100);

    hSelection[nprod] = new TH2F(Form("h%dSelection",nprod),Form("%s; SPD (%); VZERO (%);",filename),10,0,100,10,0,100);

    hPeriodMult[nprod] = new TH1F(Form("h%dPeriodMult",nprod),Form("%s;SPD Tracklets;",filename),50,0,100);
    hPeriodZDC[nprod] = new TH1F(Form("h%dPeriodZDC",nprod),Form("%s;ZDC Sum (a.u.);",filename),200,0,2000);
    
    fill(filename, nprod);
    TCanvas *c = new TCanvas;
    c->Divide(2,1);
    c->cd(1);
    hMult[nprod]->Draw("surf3");
    c->cd(2);
    hZDC[nprod]->Draw("surf3");
    hMult[nprod]->SetMaximum(20);
    hZDC[nprod]->SetMaximum(800);

    for(int i=0; i <10; i++){
      for(int j=0; j <10; j++){
	x[i+j*10] = hMult[nprod]->GetBinContent(i+1,j+1);
	y[i+j*10] = hZDC[nprod]->GetBinContent(i+1,j+1);

	hMultRatio[nprod]->Fill(hMult[nprod]->GetBinContent(i+1,j+1) / hMult[0]->GetBinContent(i+1,j+1));
        hZDCRatio[nprod]->Fill(hZDC[nprod]->GetBinContent(i+1,j+1) / hZDC[0]->GetBinContent(i+1,j+1));
      }
    }
    
    hScan[nprod] = new TGraph(100,x,y);
    hScan[nprod]->SetMarkerStyle(marker[nprod]);
    hScan[nprod]->SetMarkerSize(1.5);
    hScan[nprod]->SetMarkerColor(color[nprod]);
    hScan[nprod]->SetLineColor(color[nprod]);
    hScan[nprod]->SetName(filename);
    nprod++;
  }

  int run;
  TProfile *hncp = hRunZDCNC->ProfileX();
  TProfile *hnap = hRunZDCNA->ProfileX();
  TProfile *hpcp = hRunZDCPC->ProfileX();
  TProfile *hpap = hRunZDCPA->ProfileX();
  for(int i=2; i <= hRunZDCNC->GetNbinsX(); i++){
    if(hncp->GetBinContent(i) > 10 && printcalib) printf("calZNC[%s]=%f;\n",hRunZDCNC->GetXaxis()->GetBinLabel(i),hncp->GetBinContent(i)/refZDCN);
  }
  for(int i=2; i <= hRunZDCNA->GetNbinsX(); i++){
    if(hnap->GetBinContent(i) > 10 && printcalib) printf("calZNA[%s]=%f;\n",hRunZDCNA->GetXaxis()->GetBinLabel(i),hnap->GetBinContent(i)/refZDCN);
  }
  for(int i=2; i <= hRunZDCNC->GetNbinsX(); i++){
    if(hpcp->GetBinContent(i) > 10 && printcalib) printf("calZPC[%s]=%f;\n",hRunZDCNA->GetXaxis()->GetBinLabel(i),hpcp->GetBinContent(i)/refZDCP);
  }
  for(int i=2; i <= hRunZDCNC->GetNbinsX(); i++){
    if(hpap->GetBinContent(i) > 10 && printcalib) printf("calZPA[%s]=%f;\n",hRunZDCPA->GetXaxis()->GetBinLabel(i),hpap->GetBinContent(i)/refZDCP);
  }
  
  TFile *fout = new TFile("perfSel.root","RECREATE");
  hZDCNpro->Write();
  hZDCPpro->Write();
  for(int i=0; i < nprod; i++){
    hMult[i]->Write();
    hZDC[i]->Write();
    hScan[i]->Write();
    hMultRatio[i]->Write();
    hZDCRatio[i]->Write();
    hPeriodMult[i]->Write();
    hPeriodZDC[i]->Write();
  }
  hRunMult->Write();
  hRunZDC->Write();
  hRunZDCdist->Write();
  fout->Close();

  fout = new TFile("cal.root","RECREATE");
  for(int i=0; i < nprod; i++){
    hZDCNC[i]->Write();
    hZDCNA[i]->Write();
    hZDCPC[i]->Write();
    hZDCPA[i]->Write();
    hMultSPD[i]->Write();
    hSelection[i]->Write();
  }
  hRunZDCNC->Write();
  hRunZDCNA->Write();
  hRunZDCPC->Write();
  hRunZDCPA->Write();
  fout->Close();
}

void fill(const char* filename, int count){
  printf("processing period %d -> %s\n",count,filename);
  TFile fin(filename);
  TTree *t = (TTree *) fin.Get("PWGLF_StrVsMult/fTreeEvent");

  int n = t->GetEntries();

  float spd,vzero;
  float zdc, nch, zdcn, zdcp;
  float zdcnc,zdcna,zdcpc,zdcpa;

  int nrun;
  
  TProfile2D *hM = hMult[count];
  TProfile2D *hZ = hZDC[count];
  
  for(int i=0; i < t->GetEntries(); i+=1){
    t->GetEvent(i);

    //    if(t->GetLeaf("fRun")->GetValue() < 225106){
    //      continue;
    //    }
    
    //    if(t->GetLeaf("fSPDtracklets")->GetValue() < 1){
    //      continue;
    //    }

    nrun = t->GetLeaf("fRun")->GetValue();
    
    if(nrun == 225052 || nrun == 225051 || nrun == 225050 || nrun == 225043 || nrun == 225041 || nrun ==  225037 || nrun == 225035 || nrun == 225031 || nrun == 225026){
      continue;
    }

    spd = t->GetLeaf("fCentrality_SPDClusters")->GetValue();
    vzero = t->GetLeaf("fCentrality_V0M")->GetValue();
    zdcnc = t->GetLeaf("fZNCpp")->GetValue();
    zdcna = t->GetLeaf("fZNApp")->GetValue();
    zdcpc = t->GetLeaf("fZPCpp")->GetValue();
    zdcpa = t->GetLeaf("fZPApp")->GetValue();
  
    zdcn = zdcnc / calZNC[nrun]; //calZDC[count][0];
    zdcn += zdcna / calZNA[nrun]; //calZDC[count][1];
    zdcp = zdcpc / calZPC[nrun]; //calZDC[count][2];
    zdcp += zdcpa / calZPA[nrun]; //calZDC[count][3];
    zdc = zdcn + zdcp;
    zdc /= adjZDC[count];

    if(t->GetLeaf("fCentrality_ZDCFired")->GetValue() > 100){
      continue;
    }
    
    nch = t->GetLeaf("fSPDtracklets")->GetValue() / multScale[count];

    hM->Fill(spd, vzero, nch);
    hZ->Fill(spd, vzero, zdc);
    hSelection[count]->Fill(spd, vzero);

    hZDCNpro->Fill(vzero, zdcn);
    hZDCPpro->Fill(vzero, zdcp);

    if(zdc > 20){   
      hZDCNC[count]->Fill(vzero, t->GetLeaf("fZNCpp")->GetValue());
      hZDCNA[count]->Fill(vzero, t->GetLeaf("fZNApp")->GetValue());
      hZDCPC[count]->Fill(vzero, t->GetLeaf("fZPCpp")->GetValue());
      hZDCPA[count]->Fill(vzero, t->GetLeaf("fZPApp")->GetValue());
      
      hMultSPD[count]->Fill(spd, t->GetLeaf("fSPDtracklets")->GetValue());
    }
    
    hRunMult->Fill(Form("%d",int(t->GetLeaf("fRun")->GetValue())), nch);
    hRunZDC->Fill(Form("%d",int(t->GetLeaf("fRun")->GetValue())), zdc);
    if(zdc > 0 && zdc < 2000) hRunZDCdist->Fill(Form("%d",int(t->GetLeaf("fRun")->GetValue())), Form("%d",int(zdc/20)*20),1.);
    if(zdcnc > 0 && zdcnc < 1000) hRunZDCNC->Fill(Form("%d",int(t->GetLeaf("fRun")->GetValue())), Form("%d",int(zdcnc/10)*10),1.);
    if(zdcna > 0 && zdcna < 1000) hRunZDCNA->Fill(Form("%d",int(t->GetLeaf("fRun")->GetValue())), Form("%d",int(zdcna/10)*10),1.);
    if(zdcpc > 0 && zdcpc < 1000) hRunZDCPC->Fill(Form("%d",int(t->GetLeaf("fRun")->GetValue())), Form("%d",int(zdcpc/10)*10),1.);
    if(zdcpa > 0 && zdcpa < 1000) hRunZDCPA->Fill(Form("%d",int(t->GetLeaf("fRun")->GetValue())), Form("%d",int(zdcpa/10)*10),1.);
    hPeriodMult[count]->Fill(nch);
    hPeriodZDC[count]->Fill(zdc);
   }
  
  fin.Close();
}

void setcal(){
  for(int i=0; i < 300000; i++){
    calZNC[i] = 1;
    calZNA[i] = 1;
    calZPC[i] = 1;
    calZPA[i] = 1;
  }

  // put calibrations below
calZNC[225026]=0.618391;
calZNC[225031]=0.604778;
calZNC[225035]=0.608277;
calZNC[225037]=0.618396;
calZNC[225041]=0.619290;
calZNC[225043]=0.623473;
calZNC[225050]=0.614762;
calZNC[225051]=0.610657;
calZNC[225052]=0.609806;
calZNC[225106]=0.972489;
calZNC[225305]=1.012351;
calZNC[225307]=1.012778;
calZNC[225309]=1.009719;
calZNC[225310]=1.015562;
calZNC[225313]=1.014303;
calZNC[225314]=1.012674;
calZNC[225315]=1.012687;
calZNC[225322]=1.013926;
calZNC[225576]=0.948686;
calZNC[225578]=0.943114;
calZNC[225579]=0.945919;
calZNC[225586]=0.939245;
calZNC[225587]=0.941783;
calZNC[225705]=0.998331;
calZNC[225707]=1.005177;
calZNC[225708]=1.007029;
calZNC[225709]=1.002363;
calZNC[225710]=1.008924;
calZNC[225716]=1.008302;
calZNC[225717]=1.009229;
calZNC[225719]=1.004874;
calZNC[226062]=0.930859;
calZNC[226220]=1.000360;
calZNC[226225]=1.000393;
calZNC[226444]=0.929292;
calZNC[226445]=0.927996;
calZNC[226452]=0.926484;
calZNC[226466]=0.922750;
calZNC[226468]=0.926384;
calZNC[226472]=0.918173;
calZNC[226483]=0.928372;
calZNC[226495]=0.960105;
calZNC[226500]=0.979162;
calZNC[274593]=0.804307;
calZNC[274595]=0.820072;
calZNC[274596]=0.850225;
calZNC[274601]=0.820727;
calZNC[274653]=0.807770;
calZNC[274657]=0.816646;
calZNC[274667]=0.822760;
calZNC[274669]=0.860158;
calZNC[274671]=0.861408;
calZNC[288861]=0.814726;
calZNC[288862]=0.789854;
calZNC[288863]=0.774181;
calZNC[288864]=0.787751;
calZNC[288868]=0.763914;
calZNC[288902]=0.818994;
calZNC[288903]=0.693714;
calZNC[288908]=0.681814;
calZNC[288909]=0.677778;
calZNA[225026]=0.680101;
calZNA[225031]=0.679574;
calZNA[225035]=0.676620;
calZNA[225037]=0.681085;
calZNA[225041]=0.681085;
calZNA[225043]=0.686537;
calZNA[225050]=0.677094;
calZNA[225051]=0.677850;
calZNA[225052]=0.688838;
calZNA[225106]=1.026481;
calZNA[225305]=1.039856;
calZNA[225307]=1.041056;
calZNA[225309]=1.042120;
calZNA[225310]=1.045719;
calZNA[225313]=1.047292;
calZNA[225314]=1.041202;
calZNA[225315]=1.046653;
calZNA[225322]=1.045784;
calZNA[225576]=1.040063;
calZNA[225578]=1.039818;
calZNA[225579]=1.043311;
calZNA[225586]=1.048425;
calZNA[225587]=1.047026;
calZNA[225705]=1.027321;
calZNA[225707]=1.038868;
calZNA[225708]=1.040850;
calZNA[225709]=1.038762;
calZNA[225710]=1.049552;
calZNA[225716]=1.050831;
calZNA[225717]=1.055336;
calZNA[225719]=1.051108;
calZNA[226062]=0.979882;
calZNA[226220]=1.043605;
calZNA[226225]=1.047788;
calZNA[226444]=0.951476;
calZNA[226445]=0.954208;
calZNA[226452]=0.953552;
calZNA[226466]=0.950414;
calZNA[226468]=0.952056;
calZNA[226472]=0.950569;
calZNA[226483]=1.032522;
calZNA[226495]=1.081936;
calZNA[226500]=1.117612;
calZNA[274593]=0.897741;
calZNA[274595]=0.926096;
calZNA[274596]=0.930601;
calZNA[274601]=0.887910;
calZNA[274653]=0.877006;
calZNA[274657]=0.883778;
calZNA[274667]=0.883975;
calZNA[274669]=0.930154;
calZNA[274671]=0.929753;
calZNA[288861]=0.857138;
calZNA[288862]=0.823391;
calZNA[288863]=0.828552;
calZNA[288864]=0.816812;
calZNA[288868]=0.809478;
calZNA[288902]=0.893094;
calZNA[288903]=0.895676;
calZNA[288908]=0.860983;
calZNA[288909]=0.862350;
calZPC[225026]=1.151780;
calZPC[225031]=1.131607;
calZPC[225035]=1.112768;
calZPC[225037]=1.148523;
calZPC[225041]=1.142203;
calZPC[225043]=1.159697;
calZPC[225050]=1.127227;
calZPC[225051]=1.131536;
calZPC[225052]=1.145854;
calZPC[225106]=1.003392;
calZPC[225305]=1.019013;
calZPC[225307]=1.024539;
calZPC[225309]=1.025325;
calZPC[225310]=1.032439;
calZPC[225313]=1.031985;
calZPC[225314]=1.037478;
calZPC[225315]=1.029629;
calZPC[225322]=1.030048;
calZPC[225576]=1.025669;
calZPC[225578]=1.023696;
calZPC[225579]=1.015383;
calZPC[225586]=1.020168;
calZPC[225587]=1.022718;
calZPC[225705]=0.991768;
calZPC[225707]=1.012048;
calZPC[225708]=1.014445;
calZPC[225709]=0.999665;
calZPC[225710]=1.023736;
calZPC[225716]=1.022173;
calZPC[225717]=1.032815;
calZPC[225719]=1.022380;
calZPC[226062]=0.992213;
calZPC[226220]=1.020119;
calZPC[226225]=1.021057;
calZPC[226444]=0.989543;
calZPC[226445]=0.998756;
calZPC[226452]=0.998902;
calZPC[226466]=0.995001;
calZPC[226468]=0.992226;
calZPC[226472]=1.000960;
calZPC[226483]=1.001375;
calZPC[226495]=1.071444;
calZPC[226500]=1.128116;
calZPC[274593]=0.848674;
calZPC[274595]=0.831793;
calZPC[274596]=1.029565;
calZPC[274601]=0.823529;
calZPC[274653]=0.845691;
calZPC[274657]=0.875099;
calZPC[274667]=0.891961;
calZPC[274669]=1.087647;
calZPC[274671]=1.125912;
calZPC[288861]=0.937071;
calZPC[288862]=0.898419;
calZPC[288863]=0.905629;
calZPC[288864]=0.932558;
calZPC[288868]=0.898199;
calZPC[288902]=0.937276;
calZPC[288903]=0.748719;
calZPC[288908]=0.734559;
calZPC[288909]=0.725329;
calZPA[225026]=0.859956;
calZPA[225031]=0.858352;
calZPA[225035]=0.850923;
calZPA[225037]=0.856543;
calZPA[225041]=0.859692;
calZPA[225043]=0.871729;
calZPA[225050]=0.871842;
calZPA[225051]=0.863457;
calZPA[225052]=0.865793;
calZPA[225106]=0.760234;
calZPA[225305]=0.779644;
calZPA[225307]=0.790238;
calZPA[225309]=0.781908;
calZPA[225310]=0.805168;
calZPA[225313]=0.800494;
calZPA[225314]=0.811734;
calZPA[225315]=0.803928;
calZPA[225322]=0.807729;
calZPA[225576]=0.813563;
calZPA[225578]=0.820989;
calZPA[225579]=0.819059;
calZPA[225586]=0.820139;
calZPA[225587]=0.822957;
calZPA[225705]=0.801446;
calZPA[225707]=0.825234;
calZPA[225708]=0.821644;
calZPA[225709]=0.822279;
calZPA[225710]=0.838677;
calZPA[225716]=0.838669;
calZPA[225717]=0.838768;
calZPA[225719]=0.832107;
calZPA[226062]=0.816675;
calZPA[226220]=0.827820;
calZPA[226225]=0.837408;
calZPA[226444]=0.804573;
calZPA[226445]=0.803040;
calZPA[226452]=0.801528;
calZPA[226466]=0.802182;
calZPA[226468]=0.802992;
calZPA[226472]=0.811410;
calZPA[226483]=0.792640;
calZPA[226495]=0.878354;
calZPA[226500]=0.935314;
calZPA[274593]=0.729432;
calZPA[274595]=0.759041;
calZPA[274596]=0.810394;
calZPA[274601]=0.754625;
calZPA[274653]=0.759995;
calZPA[274657]=0.764551;
calZPA[274667]=0.763023;
calZPA[274669]=0.822284;
calZPA[274671]=0.820420;
calZPA[288861]=0.754475;
calZPA[288862]=0.717979;
calZPA[288863]=0.683893;
calZPA[288864]=0.712046;
calZPA[288868]=0.675636;
calZPA[288902]=0.800369;
calZPA[288903]=0.814714;
calZPA[288908]=0.774062;
calZPA[288909]=0.773663;
}
