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

void plotHitEnergy2(string inputDirectory, int hittype=1, float radius=0.4){
  setNCUStyle(true);
  
  string decversion[2]={"rfull009","rfull012"};

  TString inputDir[2];
  if(inputDirectory.find(decversion[0].data())!=std::string::npos){
    inputDir[0]=Form("%s",inputDirectory.data());   
    inputDir[1]=gSystem->GetFromPipe(Form("file=%s; echo \"${file//%s/%s}\"",inputDir[0].Data(),decversion[0].data(),decversion[1].data()));
    cout << "directory 1 " << inputDir[0] << endl;
    cout << "directory 2 " << inputDir[1] << endl;
  }
  else if(inputDirectory.find(decversion[1].data())!=std::string::npos){
    inputDir[1]=Form("%s",inputDirectory.data());   
    inputDir[0]=gSystem->GetFromPipe(Form("file=%s; echo \"${file//%s/%s}\"",inputDir[1].Data(),decversion[1].data(),decversion[0].data()));
    cout << "directory 1 " << inputDir[0] << endl;
    cout << "directory 2 " << inputDir[1] << endl;
  }


  TH1F* hecal[2];
  TH1F* hhcal[2];
  std::string ecalhitname = hittype==1? "EM_BARREL":"EcalBarrelHits";
  std::string hcalhitname = hittype==1? "HAD_BARREL":"HcalBarrelHits";

  TString energy=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; test2=${test%%mumu*}; echo \"${test2##*tev}\"",inputDirectory.data()));
  cout << "energy = " << energy.Data() << endl;
  

  TFile* f[2];
  for(int ifile=0; ifile<2; ifile++){

    TString inputFile = inputDir[ifile] + "/radius" + Form("%0.1f",radius)+ (hittype==1? "_rawhit.root": "_rawhit_new.root");
    f[ifile] = TFile::Open(inputFile.Data());
    hecal[ifile] = (TH1F*)f[ifile]->FindObjectAny("hecalhit_energy");
    hecal[ifile]->SetLineWidth(3);
    hecal[ifile]->SetLineColor(kBlue-ifile);
    if(ifile==1)
      {
	hecal[ifile]->SetFillColor(kBlue-ifile);
	hecal[ifile]->SetFillStyle(1001);
	hecal[ifile]->SetFillColorAlpha(kBlue-ifile,0.5);
      }
    hecal[ifile]->SetXTitle(Form("%s hit energy [GeV]",ecalhitname.data()));
    hecal[ifile]->SetYTitle(Form("Number of hits per %.1f GeV",hecal[ifile]->GetBinWidth(1)));
    hecal[ifile]->SetTitleOffset(1.2,"X");
    hecal[ifile]->SetTitleOffset(1.4,"Y");
    hhcal[ifile] = (TH1F*)f[ifile]->FindObjectAny("hhcalhit_energy");
    hhcal[ifile]->SetLineWidth(3);
    hhcal[ifile]->SetLineColor(kRed-ifile);
    if(ifile==1)
      {
	hhcal[ifile]->SetFillColor(kRed-ifile);
	hhcal[ifile]->SetFillStyle(1001);
	hhcal[ifile]->SetFillColorAlpha(kRed-ifile,0.5);
      }
    hhcal[ifile]->SetXTitle(Form("%s hit energy [GeV]",hcalhitname.data()));
    hhcal[ifile]->SetYTitle(Form("Number of hits per %.1f GeV",hhcal[ifile]->GetBinWidth(1)));
    hhcal[ifile]->SetTitleOffset(1.2,"X");
    hhcal[ifile]->SetTitleOffset(1.4,"Y");

  }

  TLegend* leg = new TLegend(0.444,0.690,0.990,0.903);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  TCanvas* c1 = new TCanvas("c1","",500,500);
  int lastbin=hecal[0]->FindLastBinAbove(0) > hecal[1]->FindLastBinAbove(0)? 
    hecal[0]->FindLastBinAbove(0)+20:
    hecal[1]->FindLastBinAbove(0)+20;
  float xmax=hecal[0]->GetBinLowEdge(lastbin);
  float ymax=hecal[0]->GetMaximum() > hecal[1]->GetMaximum()?
    1.2*hecal[0]->GetMaximum(): 1.2*hecal[1]->GetMaximum();

  hecal[0]->GetXaxis()->SetRangeUser(0,xmax);
  hecal[0]->SetMaximum(ymax);
  hecal[0]->DrawNormalized("hist");
  hecal[1]->DrawNormalized("histsame");
  c1->SetLogy(1);
  
  leg->SetHeader(Form("%s-TeV Z'#rightarrow qq",energy.Data()));
  for(int ifile=0; ifile<2; ifile++)
    leg->AddEntry(hecal[ifile],decversion[ifile].data(),ifile==0? "l":"f");
  leg->Draw("same");

  std::string outputname = decversion[1] + "/" + decversion[0] + "_" + decversion[1] + "_" + ecalhitname + "hit_energy_" + Form("%s",energy.Data()) + "TeVZp";
  c1->Print(Form("%s.pdf",outputname.data()));
  c1->Print(Form("%s.eps",outputname.data()));


  lastbin=hhcal[0]->FindLastBinAbove(0) > hhcal[1]->FindLastBinAbove(0)? 
    hhcal[0]->FindLastBinAbove(0)+20:
    hhcal[1]->FindLastBinAbove(0)+20;
  xmax=hhcal[0]->GetBinLowEdge(lastbin);
  ymax=hhcal[0]->GetMaximum() > hhcal[1]->GetMaximum()?
    1.2*hhcal[0]->GetMaximum(): 1.2*hhcal[1]->GetMaximum();

  hhcal[0]->GetXaxis()->SetRangeUser(0,xmax);
  hhcal[0]->SetMaximum(ymax);
  hhcal[0]->DrawNormalized("hist");
  hhcal[1]->DrawNormalized("histsame");
  c1->SetLogy(1);
  
  leg->Clear();
  leg->SetHeader(Form("%s-TeV Z'#rightarrow qq",energy.Data()));
  for(int ifile=0; ifile<2; ifile++)
    leg->AddEntry(hhcal[ifile],decversion[ifile].data(),ifile==0? "l":"f");
  leg->Draw("same");


  outputname = decversion[1] + "/" + decversion[0] + "_" + decversion[1] + "_" + hcalhitname + "hit_energy_" + Form("%s",energy.Data()) + "TeVZp";

  c1->Print(Form("%s.pdf",outputname.data()));
  c1->Print(Form("%s.eps",outputname.data()));

}
  
		 

