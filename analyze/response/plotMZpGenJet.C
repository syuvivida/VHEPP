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

void plotMZpGenJet(string inputDirectory){
  setNCUStyle(true);
  
  string decversion[3]={"rfull009","rfull010","rfull012"};

  TString energy=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; test2=${test%%mumu*}; echo \"${test2##*tev}\"",inputDirectory.data()));
  cout << "energy = " << energy.Data() << endl;

  int nfiles=3;
  if(energy=="20")nfiles=2;
  TFile* f[3];
  TH1F* hecal[3];
  for(int ifile=0; ifile<nfiles; ifile++){

    TString inputFile = Form("genjet_tev%smumu_pythia6_zprime%stev_qq%s.root",energy.Data(),energy.Data(),decversion[ifile].data());
    f[ifile] = TFile::Open(inputFile.Data());
    hecal[ifile] = (TH1F*)f[ifile]->FindObjectAny("hmjj");
    hecal[ifile]->SetLineWidth(3);
    hecal[ifile]->SetLineColor(ifile+1);
    hecal[ifile]->SetXTitle("Generator-level M_{jj} [GeV]");
    hecal[ifile]->SetTitleOffset(1.2,"X");
    hecal[ifile]->SetTitleOffset(1.4,"Y");
  }

  TLegend* leg = new TLegend(0.444,0.690,0.990,0.903);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  TCanvas* c1 = new TCanvas("c1","",500,500);
  
  for(int ifile=0; ifile<nfiles; ifile++){
    if(ifile==0)
      hecal[ifile]->DrawNormalized();
    else
      hecal[ifile]->DrawNormalized("same");
  } 
 
  leg->SetHeader(Form("%s-TeV Z'#rightarrow q#bar{q}",energy.Data()));
  for(int ifile=0; ifile<nfiles; ifile++)
    leg->AddEntry(hecal[ifile],decversion[ifile].data(),"l");
  leg->Draw("same");

  std::string outputname = Form("%sTeV_MZpMjj",energy.Data()) ;
  c1->Print(Form("%s.pdf",outputname.data()));
  c1->Print(Form("%s.eps",outputname.data()));
  c1->Print(Form("%s.gif",outputname.data()));

}
  
		 

