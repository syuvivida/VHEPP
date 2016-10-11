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

void plotHitEnergy(string inputDir, int hittype=1, float radius=0.4){
  setNCUStyle(true);
  
  string decversion;
  if(inputDir.find("rfull009")!=std::string::npos)decversion="rfull009";
  else if(inputDir.find("rfull012")!=std::string::npos)decversion="rfull012";

  TH1F* hecal;
  TH1F* hhcal;
  std::string ecalhitname = hittype==1? "EM_BARREL":"EcalBarrelHits";
  std::string hcalhitname = hittype==1? "HAD_BARREL":"HcalBarrelHits";

  TString energy=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; test2=${test%%mumu*}; echo \"${test2##*tev}\"",inputDir.data()));
  cout << "energy = " << energy.Data() << endl;
  
  string inputFile = inputDir + "/radius" + Form("%0.1f",radius)+ (hittype==1? "_rawhit.root": "_rawhit_new.root");
  cout << "opening " << inputFile.data() << endl;


  
  TFile *f = TFile::Open(inputFile.data());
  hecal = (TH1F*)f->FindObjectAny("hecalhit_energy");
  hecal->SetLineWidth(3);
  hecal->SetFillStyle(1001);
  hecal->SetFillColor(4);
  hecal->SetLineColor(4);
  hecal->SetXTitle(Form("%s hit energy [GeV]",ecalhitname.data()));
  hecal->SetYTitle(Form("Number of hits per %.1f GeV",hecal->GetBinWidth(1)));
  hecal->SetTitleOffset(1.2,"X");
  hecal->SetTitleOffset(1.4,"Y");
  hhcal = (TH1F*)f->FindObjectAny("hhcalhit_energy");
  hhcal->SetLineWidth(3);
  hhcal->SetFillStyle(1001);
  hhcal->SetFillColor(2);
  hhcal->SetLineColor(2);
  hhcal->SetXTitle(Form("%s hit energy [GeV]",hcalhitname.data()));
  hhcal->SetYTitle(Form("Number of hits per %.1f GeV",hhcal->GetBinWidth(1)));
  hhcal->SetTitleOffset(1.2,"X");
  hhcal->SetTitleOffset(1.4,"Y");


  TLegend* leg = new TLegend(0.444,0.690,0.990,0.903);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  TCanvas* c1 = new TCanvas("c1","",500,500);
  int lastbin=hecal->FindLastBinAbove(0)+20;
  float xmax=hecal->GetBinLowEdge(lastbin);
  hecal->GetXaxis()->SetRangeUser(0,xmax);
  hecal->Draw("hist");
  c1->SetLogy(1);
  
  leg->SetHeader(decversion.data());
  leg->AddEntry(hecal,Form("%s-TeV Z'#rightarrow qq",energy.Data()),"f");
  leg->Draw("same");

  std::string outputname = decversion + "/" + decversion + "_" + ecalhitname + "hit_energy_" + Form("%s",energy.Data()) + "TeVZp";
  c1->Print(Form("%s.pdf",outputname.data()));
  c1->Print(Form("%s.eps",outputname.data()));
  leg->Clear();


  lastbin=hhcal->FindLastBinAbove(0)+20;
  xmax=hhcal->GetBinLowEdge(lastbin);
  hhcal->GetXaxis()->SetRangeUser(0,xmax);
  hhcal->Draw("hist");
  c1->SetLogy(1);
  
  leg->SetHeader(decversion.data());
  leg->AddEntry(hhcal,Form("%s-TeV Z'#rightarrow qq",energy.Data()),"f");
  leg->Draw("same");

  outputname = decversion + "/" + decversion + "_" + hcalhitname + "hit_energy_" + Form("%s",energy.Data()) + "TeVZp";

  c1->Print(Form("%s.pdf",outputname.data()));
  c1->Print(Form("%s.eps",outputname.data()));

}
  
		 

