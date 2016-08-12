#include <string>
#include <map>
#include "setNCUStyle.C"
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>

using namespace std;
void plotJetContent(std::string inputFile)
{
  setNCUStyle(true);
  std::map<unsigned int, std::string> PDGIDs;

  const int nparts=2300;
  const int lastp = nparts-1;
  typedef std::map<unsigned int,std::string>::iterator it_type;

  PDGIDs[22]  = "pho";
  PDGIDs[11]  = "ele";
  PDGIDs[13]  = "muo";
  PDGIDs[211] = "pion";
  PDGIDs[321] = "kaon";
  PDGIDs[130] = "KL";
  PDGIDs[310] = "Ks";
  PDGIDs[2212] = "p";
  PDGIDs[2112] = "n";

  TGraphErrors *gtemp[10];

  int COLOR[2300];
  COLOR[22] = 2;
  COLOR[11] = 9;
  COLOR[13] = kPink-1;
  COLOR[211] = 1;
  COLOR[321] = 4;
  COLOR[130] = kCyan;
  COLOR[310] = kOrange+4;
  COLOR[2212] = kViolet;
  COLOR[2112] = kGreen+1;

  int MARKER[2300];
  MARKER[11] = 24;
  MARKER[13] = 25;
  MARKER[22] = 22;
  MARKER[211] = 20;
  MARKER[321] = 23;
  MARKER[130] = 34;
  MARKER[310] = 33;
  MARKER[2212] = 29;
  MARKER[2112] = 21;

  TFile *f = TFile::Open(inputFile.data());
 
  TCanvas* c1 = new TCanvas("c1","",600,600);
  
  int i=0;

  TMultiGraph *mg = new TMultiGraph();

  for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) {
    
    gtemp[i] = dynamic_cast<TGraphErrors*>(gDirectory->Get(Form("gr_Mean_%s",iterator->second.data())));
    gtemp[i]->SetMarkerColor(COLOR[iterator->first]);
    gtemp[i]->SetLineColor(COLOR[iterator->first]);
    gtemp[i]->SetMarkerStyle(MARKER[iterator->first]);
    gtemp[i]->SetMarkerSize(1.2);
    gtemp[i]->SetName(Form("gr_Mean_%s",iterator->second.data()));
    mg->Add(gtemp[i]);

    std::string tempName =gtemp[i]->GetName();
    i++;
    
  }


  mg->Draw("AP");
  mg->GetXaxis()->SetNdivisions(5);
  mg->GetXaxis()->SetTitle("E_{true} [GeV]");
  mg->GetYaxis()->SetRangeUser(0,0.6);
  mg->GetYaxis()->SetTitle("Mean of E_{particle}/E_{true}");

  TLegend* leg = new TLegend(0.156,0.763,0.931,0.912);
  leg-> SetNColumns(3);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) {
    leg->AddEntry(Form("gr_Mean_%s",iterator->second.data()),iterator->second.data(), "p");  
  }

  leg->Draw("same");
  c1->Print("genJetContent.gif");  
  c1->Print("genJetContent.pdf");

}
