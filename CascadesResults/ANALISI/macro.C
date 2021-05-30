void macro()
{
//=========Macro generated from canvas: EE/
//=========  (Tue Apr 20 10:56:56 2021) by ROOT version 6.20/08
   TCanvas *EE = new TCanvas("EE", "",270,254,1300,1027);
   EE->SetHighLightColor(2);
   EE->Range(-46.66667,-0.44,120,2.493333);
   EE->SetFillColor(0);
   EE->SetBorderSize(2);
   EE->SetGridy();
   EE->SetTickx(1);
   EE->SetTicky(1);
   EE->SetLeftMargin(0.25);
   EE->SetBottomMargin(0.15);
   
   TH1D *g__1 = new TH1D("g__1"," ",12,-5,105);
   g__1->SetBinContent(1,-0.02);
   g__1->SetBinContent(12,2.2);
   g__1->SetMinimum(0);
   g__1->SetMaximum(2.2);
   g__1->SetEntries(2);
   g__1->SetStats(0);
   g__1->SetLineColor(0);
   g__1->GetXaxis()->SetTitle("percentile ZDC (%)");
   g__1->GetXaxis()->SetLabelFont(42);
   g__1->GetXaxis()->SetTitleSize(0.04);
   g__1->GetXaxis()->SetTitleOffset(1.2);
   g__1->GetXaxis()->SetTitleFont(42);
   g__1->GetYaxis()->SetTitle("#left( #frac{ d#it{N}/d#it{y}}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{Sel} / #left( #frac{ dN/dy}{#LT d#it{N}_{ch}/d#it{#eta} #GT _{|#it{#eta}|<0.5}} #right)_{MB}");
   g__1->GetYaxis()->SetLabelFont(42);
   g__1->GetYaxis()->SetTitleSize(0.04);
   g__1->GetYaxis()->SetTitleOffset(1.5);
   g__1->GetYaxis()->SetTitleFont(42);
   g__1->GetZaxis()->SetLabelFont(42);
   g__1->GetZaxis()->SetTitleOffset(1);
   g__1->GetZaxis()->SetTitleFont(42);
   g__1->Draw("");
   
   Double_t NormYieldsvspercentile_Syst_fx3001[9] = {
   10,
   25,
   35,
   45,
   55,
   65,
   75,
   85,
   95};
   Double_t NormYieldsvspercentile_Syst_fy3001[9] = {
   1.140036,
   1.137278,
   1.099444,
   1.11283,
   1.137516,
   1.069624,
   1.19098,
   1.130483,
   1.091707};
   Double_t NormYieldsvspercentile_Syst_felx3001[9] = {
   10,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t NormYieldsvspercentile_Syst_fely3001[9] = {
   0.03968528,
   0.03177994,
   0.03936961,
   0.03902574,
   0.04013729,
   0.04586263,
   0.1117624,
   0.09173524,
   0.06261008};
   Double_t NormYieldsvspercentile_Syst_fehx3001[9] = {
   10,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t NormYieldsvspercentile_Syst_fehy3001[9] = {
   0.03968528,
   0.03177994,
   0.03936961,
   0.03902574,
   0.04013729,
   0.04586263,
   0.1117624,
   0.09173524,
   0.06261008};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(9,NormYieldsvspercentile_Syst_fx3001,NormYieldsvspercentile_Syst_fy3001,NormYieldsvspercentile_Syst_felx3001,NormYieldsvspercentile_Syst_fehx3001,NormYieldsvspercentile_Syst_fely3001,NormYieldsvspercentile_Syst_fehy3001);
   grae->SetName("NormYieldsvspercentile_Syst");
   grae->SetTitle("Graph");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = 1179;
   color = new TColor(ci, 1, 0.8, 1, " ", 0.1);
   grae->SetFillColor(ci);
   grae->SetFillStyle(3000);

   ci = TColor::GetColor("#990099");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#990099");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetMarkerSize(2.4);
   
   TH1F *Graph_NormYieldsvspercentile_Syst3001 = new TH1F("Graph_NormYieldsvspercentile_Syst3001","Graph",100,0,110);
   Graph_NormYieldsvspercentile_Syst3001->SetMinimum(0);
   Graph_NormYieldsvspercentile_Syst3001->SetMaximum(1.5);
   Graph_NormYieldsvspercentile_Syst3001->SetDirectory(0);
   Graph_NormYieldsvspercentile_Syst3001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_NormYieldsvspercentile_Syst3001->SetLineColor(ci);
   Graph_NormYieldsvspercentile_Syst3001->GetXaxis()->SetTitle("percentile ZDC (%)");
   Graph_NormYieldsvspercentile_Syst3001->GetXaxis()->SetRange(0,92);
   Graph_NormYieldsvspercentile_Syst3001->GetXaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Syst3001->GetXaxis()->SetTitleSize(0.04);
   Graph_NormYieldsvspercentile_Syst3001->GetXaxis()->SetTitleOffset(1.2);
   Graph_NormYieldsvspercentile_Syst3001->GetXaxis()->SetTitleFont(42);
   Graph_NormYieldsvspercentile_Syst3001->GetYaxis()->SetTitle("1/n_{ch} <dN/dy> / (1/n_{ch} <dN/dy>)_{MB}");
   Graph_NormYieldsvspercentile_Syst3001->GetYaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Syst3001->GetYaxis()->SetTitleSize(0.04);
   Graph_NormYieldsvspercentile_Syst3001->GetYaxis()->SetTitleOffset(1.2);
   Graph_NormYieldsvspercentile_Syst3001->GetYaxis()->SetTitleFont(42);
   Graph_NormYieldsvspercentile_Syst3001->GetZaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Syst3001->GetZaxis()->SetTitleOffset(1);
   Graph_NormYieldsvspercentile_Syst3001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_NormYieldsvspercentile_Syst3001);
   
   grae->Draw("e2 p ");
   
   Double_t NormYieldsvspercentile_Stat_fx1001[9] = {
   10,
   25,
   35,
   45,
   55,
   65,
   75,
   85,
   95};
   Double_t NormYieldsvspercentile_Stat_fy1001[9] = {
   1.140036,
   1.137278,
   1.099444,
   1.11283,
   1.137516,
   1.069624,
   1.19098,
   1.130483,
   1.091707};
   Double_t NormYieldsvspercentile_Stat_fex1001[9] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t NormYieldsvspercentile_Stat_fey1001[9] = {
   0.01010409,
   0.01687957,
   0.0200597,
   0.02290204,
   0.02618649,
   0.02870787,
   0.03460855,
   0.03664685,
   0.04542779};
   TGraphErrors *gre = new TGraphErrors(9,NormYieldsvspercentile_Stat_fx1001,NormYieldsvspercentile_Stat_fy1001,NormYieldsvspercentile_Stat_fex1001,NormYieldsvspercentile_Stat_fey1001);
   gre->SetName("NormYieldsvspercentile_Stat");
   gre->SetTitle("");
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#990099");
   gre->SetLineColor(ci);
   gre->SetLineWidth(2);

   ci = TColor::GetColor("#990099");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(22);
   gre->SetMarkerSize(2.4);
   
   TH1F *Graph_NormYieldsvspercentile_Stat1001 = new TH1F("Graph_NormYieldsvspercentile_Stat1001","",100,1.5,103.5);
   Graph_NormYieldsvspercentile_Stat1001->SetMinimum(0);
   Graph_NormYieldsvspercentile_Stat1001->SetMaximum(1.5);
   Graph_NormYieldsvspercentile_Stat1001->SetDirectory(0);
   Graph_NormYieldsvspercentile_Stat1001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_NormYieldsvspercentile_Stat1001->SetLineColor(ci);
   Graph_NormYieldsvspercentile_Stat1001->GetXaxis()->SetTitle("percentile ZDC (%)");
   Graph_NormYieldsvspercentile_Stat1001->GetXaxis()->SetRange(0,98);
   Graph_NormYieldsvspercentile_Stat1001->GetXaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Stat1001->GetXaxis()->SetTitleSize(0.04);
   Graph_NormYieldsvspercentile_Stat1001->GetXaxis()->SetTitleOffset(1.2);
   Graph_NormYieldsvspercentile_Stat1001->GetXaxis()->SetTitleFont(42);
   Graph_NormYieldsvspercentile_Stat1001->GetYaxis()->SetTitle("1/n_{ch} <dN/dy> / (1/n_{ch} <dN/dy>)_{MB}");
   Graph_NormYieldsvspercentile_Stat1001->GetYaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Stat1001->GetYaxis()->SetTitleSize(0.04);
   Graph_NormYieldsvspercentile_Stat1001->GetYaxis()->SetTitleOffset(1.2);
   Graph_NormYieldsvspercentile_Stat1001->GetYaxis()->SetTitleFont(42);
   Graph_NormYieldsvspercentile_Stat1001->GetZaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Stat1001->GetZaxis()->SetTitleOffset(1);
   Graph_NormYieldsvspercentile_Stat1001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_NormYieldsvspercentile_Stat1001);
   
   gre->Draw("e p ");
   
   Double_t NormYieldsvspercentile_Syst_fx3002[8] = {
   10,
   30,
   45,
   55,
   65,
   75,
   85,
   95};
   Double_t NormYieldsvspercentile_Syst_fy3002[8] = {
   0.5762461,
   0.6267354,
   0.6884544,
   0.5992229,
   0.5634744,
   0.4963417,
   0.5362466,
   0.4994598};
   Double_t NormYieldsvspercentile_Syst_felx3002[8] = {
   10,
   10,
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t NormYieldsvspercentile_Syst_fely3002[8] = {
   0.09936669,
   0.1106022,
   0.1084069,
   0.07437421,
   0.05462876,
   0.05761811,
   0.06637945,
   0.07179599};
   Double_t NormYieldsvspercentile_Syst_fehx3002[8] = {
   10,
   10,
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t NormYieldsvspercentile_Syst_fehy3002[8] = {
   0.09936669,
   0.1106022,
   0.1084069,
   0.07437421,
   0.05462876,
   0.05761811,
   0.06637945,
   0.07179599};
   grae = new TGraphAsymmErrors(8,NormYieldsvspercentile_Syst_fx3002,NormYieldsvspercentile_Syst_fy3002,NormYieldsvspercentile_Syst_felx3002,NormYieldsvspercentile_Syst_fehx3002,NormYieldsvspercentile_Syst_fely3002,NormYieldsvspercentile_Syst_fehy3002);
   grae->SetName("NormYieldsvspercentile_Syst");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#99cc99");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3000);

   ci = TColor::GetColor("#003300");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#006600");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(23);
   grae->SetMarkerSize(2.4);
   
   TH1F *Graph_NormYieldsvspercentile_Syst3002 = new TH1F("Graph_NormYieldsvspercentile_Syst3002","Graph",100,0,110);
   Graph_NormYieldsvspercentile_Syst3002->SetMinimum(0);
   Graph_NormYieldsvspercentile_Syst3002->SetMaximum(1.5);
   Graph_NormYieldsvspercentile_Syst3002->SetDirectory(0);
   Graph_NormYieldsvspercentile_Syst3002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_NormYieldsvspercentile_Syst3002->SetLineColor(ci);
   Graph_NormYieldsvspercentile_Syst3002->GetXaxis()->SetTitle("percentile ZDC (%)");
   Graph_NormYieldsvspercentile_Syst3002->GetXaxis()->SetRange(0,92);
   Graph_NormYieldsvspercentile_Syst3002->GetXaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Syst3002->GetXaxis()->SetTitleSize(0.04);
   Graph_NormYieldsvspercentile_Syst3002->GetXaxis()->SetTitleOffset(1.2);
   Graph_NormYieldsvspercentile_Syst3002->GetXaxis()->SetTitleFont(42);
   Graph_NormYieldsvspercentile_Syst3002->GetYaxis()->SetTitle("1/n_{ch} <dN/dy> / (1/n_{ch} <dN/dy>)_{MB}");
   Graph_NormYieldsvspercentile_Syst3002->GetYaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Syst3002->GetYaxis()->SetTitleSize(0.04);
   Graph_NormYieldsvspercentile_Syst3002->GetYaxis()->SetTitleOffset(1.2);
   Graph_NormYieldsvspercentile_Syst3002->GetYaxis()->SetTitleFont(42);
   Graph_NormYieldsvspercentile_Syst3002->GetZaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Syst3002->GetZaxis()->SetTitleOffset(1);
   Graph_NormYieldsvspercentile_Syst3002->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_NormYieldsvspercentile_Syst3002);
   
   grae->Draw("e2 p ");
   
   Double_t NormYieldsvspercentile_Stat_fx1002[8] = {
   10,
   30,
   45,
   55,
   65,
   75,
   85,
   95};
   Double_t NormYieldsvspercentile_Stat_fy1002[8] = {
   0.5762461,
   0.6267354,
   0.6884544,
   0.5992229,
   0.5634744,
   0.4963417,
   0.5362466,
   0.4994598};
   Double_t NormYieldsvspercentile_Stat_fex1002[8] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t NormYieldsvspercentile_Stat_fey1002[8] = {
   0.06487372,
   0.0960617,
   0.07179487,
   0.05787257,
   0.05433128,
   0.05178564,
   0.05436411,
   0.04705988};
   gre = new TGraphErrors(8,NormYieldsvspercentile_Stat_fx1002,NormYieldsvspercentile_Stat_fy1002,NormYieldsvspercentile_Stat_fex1002,NormYieldsvspercentile_Stat_fey1002);
   gre->SetName("NormYieldsvspercentile_Stat");
   gre->SetTitle("");
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#006600");
   gre->SetLineColor(ci);
   gre->SetLineWidth(2);

   ci = TColor::GetColor("#006600");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(23);
   gre->SetMarkerSize(2.4);
   
   TH1F *Graph_NormYieldsvspercentile_Stat1002 = new TH1F("Graph_NormYieldsvspercentile_Stat1002","",100,1.5,103.5);
   Graph_NormYieldsvspercentile_Stat1002->SetMinimum(0);
   Graph_NormYieldsvspercentile_Stat1002->SetMaximum(1.5);
   Graph_NormYieldsvspercentile_Stat1002->SetDirectory(0);
   Graph_NormYieldsvspercentile_Stat1002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_NormYieldsvspercentile_Stat1002->SetLineColor(ci);
   Graph_NormYieldsvspercentile_Stat1002->GetXaxis()->SetTitle("percentile ZDC (%)");
   Graph_NormYieldsvspercentile_Stat1002->GetXaxis()->SetRange(0,98);
   Graph_NormYieldsvspercentile_Stat1002->GetXaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Stat1002->GetXaxis()->SetTitleSize(0.04);
   Graph_NormYieldsvspercentile_Stat1002->GetXaxis()->SetTitleOffset(1.2);
   Graph_NormYieldsvspercentile_Stat1002->GetXaxis()->SetTitleFont(42);
   Graph_NormYieldsvspercentile_Stat1002->GetYaxis()->SetTitle("1/n_{ch} <dN/dy> / (1/n_{ch} <dN/dy>)_{MB}");
   Graph_NormYieldsvspercentile_Stat1002->GetYaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Stat1002->GetYaxis()->SetTitleSize(0.04);
   Graph_NormYieldsvspercentile_Stat1002->GetYaxis()->SetTitleOffset(1.2);
   Graph_NormYieldsvspercentile_Stat1002->GetYaxis()->SetTitleFont(42);
   Graph_NormYieldsvspercentile_Stat1002->GetZaxis()->SetLabelFont(42);
   Graph_NormYieldsvspercentile_Stat1002->GetZaxis()->SetTitleOffset(1);
   Graph_NormYieldsvspercentile_Stat1002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_NormYieldsvspercentile_Stat1002);
   
   gre->Draw("e p ");
   
   Double_t hclonehigh_fx3003[9] = {
   10,
   25,
   35,
   45,
   55,
   65,
   75,
   85,
   95};
   Double_t hclonehigh_fy3003[9] = {
   1.140036,
   1.137278,
   1.099444,
   1.11283,
   1.137516,
   1.069624,
   1.19098,
   1.130483,
   1.091707};
   Double_t hclonehigh_felx3003[9] = {
   10,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t hclonehigh_fely3003[9] = {
   0.03968528,
   0.03177994,
   0.03936961,
   0.03902574,
   0.04013729,
   0.04586263,
   0.1117624,
   0.09173524,
   0.06261008};
   Double_t hclonehigh_fehx3003[9] = {
   10,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t hclonehigh_fehy3003[9] = {
   0.03968528,
   0.03177994,
   0.03936961,
   0.03902574,
   0.04013729,
   0.04586263,
   0.1117624,
   0.09173524,
   0.06261008};
   grae = new TGraphAsymmErrors(9,hclonehigh_fx3003,hclonehigh_fy3003,hclonehigh_felx3003,hclonehigh_fehx3003,hclonehigh_fely3003,hclonehigh_fehy3003);
   grae->SetName("hclonehigh");
   grae->SetTitle("Graph");

   ci = 1179;
   color = new TColor(ci, 1, 0.8, 1, " ", 0.1);
   grae->SetFillColor(ci);
   grae->SetFillStyle(3000);

   ci = TColor::GetColor("#990099");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#990099");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetMarkerSize(2.4);
   
   TH1F *Graph_hclonehigh3003 = new TH1F("Graph_hclonehigh3003","Graph",100,0,110);
   Graph_hclonehigh3003->SetMinimum(0);
   Graph_hclonehigh3003->SetMaximum(1.5);
   Graph_hclonehigh3003->SetDirectory(0);
   Graph_hclonehigh3003->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_hclonehigh3003->SetLineColor(ci);
   Graph_hclonehigh3003->GetXaxis()->SetTitle("percentile ZDC (%)");
   Graph_hclonehigh3003->GetXaxis()->SetRange(0,92);
   Graph_hclonehigh3003->GetXaxis()->SetLabelFont(42);
   Graph_hclonehigh3003->GetXaxis()->SetTitleSize(0.04);
   Graph_hclonehigh3003->GetXaxis()->SetTitleOffset(1.2);
   Graph_hclonehigh3003->GetXaxis()->SetTitleFont(42);
   Graph_hclonehigh3003->GetYaxis()->SetTitle("1/n_{ch} <dN/dy> / (1/n_{ch} <dN/dy>)_{MB}");
   Graph_hclonehigh3003->GetYaxis()->SetLabelFont(42);
   Graph_hclonehigh3003->GetYaxis()->SetTitleSize(0.04);
   Graph_hclonehigh3003->GetYaxis()->SetTitleOffset(1.2);
   Graph_hclonehigh3003->GetYaxis()->SetTitleFont(42);
   Graph_hclonehigh3003->GetZaxis()->SetLabelFont(42);
   Graph_hclonehigh3003->GetZaxis()->SetTitleOffset(1);
   Graph_hclonehigh3003->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_hclonehigh3003);
   
   grae->Draw("e2 ");
   
   Double_t hclonelow_fx3004[8] = {
   10,
   30,
   45,
   55,
   65,
   75,
   85,
   95};
   Double_t hclonelow_fy3004[8] = {
   0.5762461,
   0.6267354,
   0.6884544,
   0.5992229,
   0.5634744,
   0.4963417,
   0.5362466,
   0.4994598};
   Double_t hclonelow_felx3004[8] = {
   10,
   10,
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t hclonelow_fely3004[8] = {
   0.09936669,
   0.1106022,
   0.1084069,
   0.07437421,
   0.05462876,
   0.05761811,
   0.06637945,
   0.07179599};
   Double_t hclonelow_fehx3004[8] = {
   10,
   10,
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t hclonelow_fehy3004[8] = {
   0.09936669,
   0.1106022,
   0.1084069,
   0.07437421,
   0.05462876,
   0.05761811,
   0.06637945,
   0.07179599};
   grae = new TGraphAsymmErrors(8,hclonelow_fx3004,hclonelow_fy3004,hclonelow_felx3004,hclonelow_fehx3004,hclonelow_fely3004,hclonelow_fehy3004);
   grae->SetName("hclonelow");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#99cc99");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3000);

   ci = TColor::GetColor("#003300");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#006600");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(23);
   grae->SetMarkerSize(2.4);
   
   TH1F *Graph_hclonelow3004 = new TH1F("Graph_hclonelow3004","Graph",100,0,110);
   Graph_hclonelow3004->SetMinimum(0);
   Graph_hclonelow3004->SetMaximum(1.5);
   Graph_hclonelow3004->SetDirectory(0);
   Graph_hclonelow3004->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_hclonelow3004->SetLineColor(ci);
   Graph_hclonelow3004->GetXaxis()->SetTitle("percentile ZDC (%)");
   Graph_hclonelow3004->GetXaxis()->SetRange(0,92);
   Graph_hclonelow3004->GetXaxis()->SetLabelFont(42);
   Graph_hclonelow3004->GetXaxis()->SetTitleSize(0.04);
   Graph_hclonelow3004->GetXaxis()->SetTitleOffset(1.2);
   Graph_hclonelow3004->GetXaxis()->SetTitleFont(42);
   Graph_hclonelow3004->GetYaxis()->SetTitle("1/n_{ch} <dN/dy> / (1/n_{ch} <dN/dy>)_{MB}");
   Graph_hclonelow3004->GetYaxis()->SetLabelFont(42);
   Graph_hclonelow3004->GetYaxis()->SetTitleSize(0.04);
   Graph_hclonelow3004->GetYaxis()->SetTitleOffset(1.2);
   Graph_hclonelow3004->GetYaxis()->SetTitleFont(42);
   Graph_hclonelow3004->GetZaxis()->SetLabelFont(42);
   Graph_hclonelow3004->GetZaxis()->SetTitleOffset(1);
   Graph_hclonelow3004->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_hclonelow3004);
   
   grae->Draw("e2 ");
   
   TLegend *leg = new TLegend(0.5,0.6238438,0.8898305,0.8129496,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.025);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("NULL","Multiplicity selection:","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("NormYieldsvspercentile_Syst","V0M [0-30] % (high multiplicity)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#990099");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(2.4);
   entry->SetTextFont(42);
   entry=leg->AddEntry("NormYieldsvspercentile_Syst","V0M [70-100] % (low multiplicity)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#006600");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(2.4);
   entry->SetTextFont(42);
   leg->Draw();
   
   TPavesText *pst = new TPavesText(3.795582,1.733621,13.68259,1.959726,5,"br");
   pst->SetBorderSize(1);
   pst->SetLineWidth(0);
   pst->SetTextFont(42);
   pst->SetTextSize(0.02904787);
   TText *pst_LaTex = pst->AddText("stat.");
   pst_LaTex = pst->AddText("syst.");
   pst->Draw();
   TMarker *marker = new TMarker(1.869542,1.908476,2);
   marker->SetMarkerStyle(2);
   marker->SetMarkerSize(2.9);
   marker->Draw();
   marker = new TMarker(1.869542,1.793916,25);
   marker->SetMarkerStyle(25);
   marker->SetMarkerSize(2.8);
   marker->Draw();
   
   TPaveText *pt = new TPaveText(0.01,0.945,0.04617874,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   TText *pt_LaTex = pt->AddText(" ");
   pt->Draw();
   TLatex *   tex = new TLatex(-0.184901,2.044138,"ALICE");
   tex->SetTextFont(42);
   tex->SetTextSize(0.03186022);
   tex->SetLineWidth(2);
   tex->Draw();
   EE->Modified();
   EE->cd();
   EE->SetSelected(EE);
   EE->ToggleToolBar();
}
