#include <string>
#include "setNCUStyle.C"
#include "TFile.h"
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>

using namespace std;
void plotOnlyGraphs3()
{

  TGraphErrors* gr[4];

  TFile *_file0 = TFile::Open("rfull009/rfull009_radius0.4_jetresponse_tcalo.root");
  gr[0] = (TGraphErrors*)(_file0->FindObjectAny("gr_Mean90"));
  gr[0]->SetName("gr0");
  gr[1] = (TGraphErrors*)(_file0->FindObjectAny("gr_RMSMean90"));
  gr[0]->SetName("gr1");
  TFile *_file1 = TFile::Open("rfull012/rfull012_radius0.4_jetresponse_tcalo.root");
  gr[2] = (TGraphErrors*)(_file1->FindObjectAny("gr_Mean90"));
  gr[2]->SetName("gr2");
  gr[3] = (TGraphErrors*)(_file1->FindObjectAny("gr_RMSMean90"));
  gr[3]->SetName("gr3");


  gr[0]->SetMarkerStyle(20);
  gr[0]->SetMarkerColor(1);
  gr[2]->SetMarkerStyle(24);
  gr[2]->SetMarkerColor(4);

  gr[1]->SetMarkerStyle(20);
  gr[1]->SetMarkerColor(1);
  gr[3]->SetMarkerStyle(24);
  gr[3]->SetMarkerColor(4);



  setNCUStyle(true);
  TMultiGraph *mg = new TMultiGraph();
  TMultiGraph *mg2 = new TMultiGraph();

  TCanvas* c1 = new TCanvas("c1","",700,500);

  mg->Add(gr[0]);
  mg->Add(gr[2]);
  mg->Draw("AP");
  mg->GetYaxis()->SetRangeUser(0.7,1.0);   
  mg->GetXaxis()->SetTitle("E_{true} [GeV]");
  mg->GetYaxis()->SetTitleOffset(1.0);
  mg->GetYaxis()->SetTitle("Mean^{90} of E_{jet}/E_{true}");


  TLegend* leg = new TLegend(0.696,0.26,0.946,0.5);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetHeader("Anti-kt Jet");
  leg->AddEntry("gr0", "rfull009","p");
  leg->AddEntry("gr2", "rfull012","p");
  leg->Draw("same");

  c1->Print("rfull009_rfull012_Mean90.pdf");
  c1->Print("rfull009_rfull012_Mean90.gif");
  c1->Clear();


  mg2->Add(gr[1]);
  mg2->Add(gr[3]);
  mg2->Draw("AP");
  mg2->GetYaxis()->SetRangeUser(0.0,0.14);   
  mg2->GetXaxis()->SetTitle("E_{true} [GeV]");
  mg2->GetYaxis()->SetTitleOffset(1.0);
  mg2->GetYaxis()->SetTitle("RMS^{90}/Mean^{90} of E_{jet}/E_{true}");


  leg->Clear();
  leg->SetHeader("Anti-kt Jet");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->AddEntry("gr1", "rfull009","p");
  leg->AddEntry("gr3", "rfull012","p");
  leg->Draw("same");


  c1->Print("rfull009_rfull012_RMSMean90.pdf");
   c1->Print("rfull009_rfull012_RMSMean90.gif");


}
