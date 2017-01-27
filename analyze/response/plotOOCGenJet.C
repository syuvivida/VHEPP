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
#include <TMultiGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include "setNCUStyle.C"
#include <TLegend.h>

void plotOOCGenJet(){
  setNCUStyle(true);

  const int nfiles=5;
  float energy[nfiles]={1,5,10,20,40};
  const int nRs=2;
  float radius[]={0.4,1.5};

  TFile* f[nfiles];
  TH1F* hecal[nfiles];
  TCanvas* c1 = new TCanvas("c1","",500,500);
  gStyle->SetOptStat(1111111);

  float x[nRs][nfiles];
  float xerr[nRs][nfiles];
  float y[nRs][nfiles];
  float yerr[nRs][nfiles];
  for(int iR=0; iR < nRs; iR++){
    for(int ifile=0; ifile<nfiles; ifile++){

      TString inputFile = Form("genjet_radius%.1f_tev%dmumu_pythia6_zprime%dtev_qqrfull009.root",radius[iR],(int)energy[ifile],(int)energy[ifile]);
      cout << inputFile << endl;
      f[ifile] = TFile::Open(inputFile.Data());
      hecal[ifile] = (TH1F*)f[ifile]->FindObjectAny("hooc");
      hecal[ifile]->SetLineWidth(3);
      //      hecal[ifile]->SetLineColor(ifile+1);
      hecal[ifile]->SetXTitle("Fraction of Out-of-cone Energy");
      hecal[ifile]->SetTitleOffset(1.2,"X");
      hecal[ifile]->SetTitleOffset(1.4,"Y");
      hecal[ifile]->SetTitle(Form("Jet radius=%.1f, %d-TeV Z'#rightarrow q#bar{q}",radius[iR],(int)energy[ifile]));
      hecal[ifile]->Draw();
      c1->Print(Form("OOC_radius%.1f_%dTeVZqq.pdf",radius[iR],(int)energy[ifile]));
      x[iR][ifile]    = energy[ifile]*0.5;
      xerr[iR][ifile] = 0.5;
      y[iR][ifile]   = hecal[ifile]->GetMean();
      yerr[iR][ifile]= hecal[ifile]->GetMeanError();
    }
  }


  TMultiGraph *mg = new TMultiGraph();
  TGraphErrors* gr[nRs];
  for(int iR=0; iR<nRs; iR++){
    gr[iR] = new TGraphErrors(nfiles,x[iR],y[iR],xerr[iR],yerr[iR]);
    gr[iR]->SetName(Form("gr%d",iR));
    gr[iR]->SetTitle("Fraction of Out-of-cone energy");
    gr[iR]->GetXaxis()->SetTitle("E_{jet} [GeV]");
    gr[iR]->GetYaxis()->SetTitle("Mean fraction");
    gr[iR]->SetMarkerStyle(8);
    gr[iR]->SetMarkerSize(1);
    gr[iR]->SetMarkerColor(iR+1);
    gr[iR]->SetLineColor(iR+1);
    gr[iR]->GetXaxis()->SetTitleSize(0.04);
    gr[iR]->GetYaxis()->SetTitleSize(0.04);
    gr[iR]->GetXaxis()->SetNdivisions(5);
    gr[iR]->GetYaxis()->SetTitleOffset(1.2);
    gr[iR]->GetYaxis()->SetDecimals();
    mg->Add(gr[iR]);
  }


  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("E_{jet} [GeV]");
  mg->GetYaxis()->SetTitleOffset(1.4);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitle("Mean Fraction of Out-of-cone Energy");
  
  TLegend* leg = new TLegend(0.696,0.26,0.946,0.5);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  for(int iR=0; iR< nRs; iR++)
    leg->AddEntry(Form("gr%d",iR), Form("#Delta R = %.1f",radius[iR]),"p");
  leg->Draw("same");

  c1->Print("OOC_rfull009_radius_0p4_1p5.pdf");
  c1->Print("OOC_rfull009_radius_0p4_1p5.eps");


}
  
		 

