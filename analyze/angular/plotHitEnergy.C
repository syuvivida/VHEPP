#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>  
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include "setNCUStyle.C"
#include <TLegend.h>

void plotHitEnergy(string inputDir, float radius=0.4){
  setNCUStyle(true);

  TH1F* hecal;
  TH1F* hhcal;


  TString energy=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; test2=${test%%mumu*}; echo \"${test2##*tev}\"",inputDir.data()));
  cout << "energy = " << energy.Data() << endl;
  
  string inputFile = inputDir + "/radius" + Form("%0.1f",radius)+ "_rawhit.root";
  cout << "opening " << inputFile.data() << endl;


  
  TFile *f = TFile::Open(inputFile.data());
  hecal = (TH1F*)f->FindObjectAny("hecalhit_energy");
  hecal->SetLineWidth(3);
  hecal->SetFillStyle(1001);
  hecal->SetFillColor(4);
  hecal->SetLineColor(4);
  hecal->SetXTitle("EM_BARREL hit energy [GeV]");
  hecal->SetYTitle(Form("Number of hits per %d GeV",(int)hecal->GetBinWidth(1)));
  hecal->SetTitleOffset(1.2,"X");
  hecal->SetTitleOffset(1.4,"Y");
  hhcal = (TH1F*)f->FindObjectAny("hhcalhit_energy");
  hhcal->SetLineWidth(3);
  hhcal->SetFillStyle(1001);
  hhcal->SetFillColor(2);
  hhcal->SetLineColor(2);
  hhcal->SetXTitle("HAD_BARREL hit energy [GeV]");
  hhcal->SetYTitle(Form("Number of hits per %d GeV",(int)hhcal->GetBinWidth(1)));
  hhcal->SetTitleOffset(1.2,"X");
  hhcal->SetTitleOffset(1.4,"Y");


  TLegend* leg = new TLegend(0.444,0.690,0.990,0.903);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  TCanvas* c1 = new TCanvas("c1","",500,500);
  int lastbin=hecal->FindLastBinAbove(0)+50;
  float xmax=hecal->GetBinLowEdge(lastbin);
  hecal->GetXaxis()->SetRangeUser(0,xmax);
  hecal->Draw("hist");
  c1->SetLogy(1);
  
  leg->SetHeader("rfull012");
  leg->AddEntry(hecal,Form("%s-TeV Z'#rightarrow qq",energy.Data()),"f");
  leg->Draw("same");

  c1->Print(Form("rfull012/rfull012_EM_BARREL_hit_energy_%sTeVZp.pdf",energy.Data()));
  c1->Print(Form("rfull012/rfull012_EM_BARREL_hit_energy_%sTeVZp.eps",energy.Data()));
  leg->Clear();


  lastbin=hhcal->FindLastBinAbove(0)+50;
  xmax=hhcal->GetBinLowEdge(lastbin);
  hhcal->GetXaxis()->SetRangeUser(0,xmax);
  hhcal->Draw("hist");
  c1->SetLogy(1);
  
  leg->SetHeader("rfull012");
  leg->AddEntry(hhcal,Form("%s-TeV Z'#rightarrow qq",energy.Data()),"f");
  leg->Draw("same");

  c1->Print(Form("rfull012/rfull012_HAD_BARREL_hit_energy_%sTeVZp.pdf",energy.Data()));
  c1->Print(Form("rfull012/rfull012_HAD_BARREL_hit_energy_%sTeVZp.eps",energy.Data()));


}
  
		 

