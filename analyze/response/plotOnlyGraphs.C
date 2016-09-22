#include <string>
#include "setNCUStyle.C"
#include "TFile.h"
#include <TMultiGraph.h>


using namespace std;
void plotOnlyGraphs(float radius)
{

  TGraphErrors* gr[4];

  TFile *_file0 = TFile::Open(Form("rfull009_radius%.1f_jetresponse_tcalo.root",radius));
  gr[0] = (TGraphErrors*)(_file0->FindObjectAny("gr_Mean90"));
  gr[0]->SetName("gr0");
  gr[1] = (TGraphErrors*)(_file0->FindObjectAny("gr_RMS90"));
  gr[0]->SetName("gr1");
  TFile *_file1 = TFile::Open(Form("rfull009_radius%.1f_jetresponse_trawhits2.root",radius));
  gr[2] = (TGraphErrors*)(_file1->FindObjectAny("gr_Mean90"));
  gr[2]->SetName("gr2");
  gr[3] = (TGraphErrors*)(_file1->FindObjectAny("gr_RMS90"));
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
  mg->SetTitle(Form("Jet #Delta R=%.1f",radius));
  mg->GetXaxis()->SetTitle("E_{true} [GeV]");
  mg->GetYaxis()->SetTitleOffset(1.0);
  mg->GetYaxis()->SetTitle("Mean^{90} of E_{jet}/E_{true}");

  
  TLegend* leg = new TLegend(0.148,0.634,0.397,0.877);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetHeader(Form("Jet #Delta R=%.1f",radius));
  leg->AddEntry("gr0", "Anti-kt","p");
  leg->AddEntry("gr2", "Simple","p");
  leg->Draw("same");

  c1->Print(Form("rfull009_radius%.1f_Mean90.pdf",radius));
  c1->Print(Form("rfull009_radius%.1f_Mean90.gif",radius));
  c1->Clear();


  mg2->Add(gr[1]);
  mg2->Add(gr[3]);
  //  mg->GetYaxis()->SetRangeUser();   
  mg2->Draw("AP");
  mg2->SetTitle(Form("Jet #Delta R=%.1f",radius));
  mg2->GetXaxis()->SetTitle("E_{true} [GeV]");
  mg2->GetYaxis()->SetTitleOffset(1.0);
  mg2->GetYaxis()->SetTitle("RMS^{90} of E_{jet}/E_{true}");


  leg->Clear();
  leg->SetHeader(Form("Jet #Delta R=%.1f",radius));
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->AddEntry("gr1", "Anti-kt","p");
  leg->AddEntry("gr3", "Simple","p");
  leg->Draw("same");


  c1->Print(Form("rfull009_radius%.1f_RMS90.pdf",radius));
  c1->Print(Form("rfull009_radius%.1f_RMS90.gif",radius));


}
