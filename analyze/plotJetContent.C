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

  TFile *f = TFile::Open(inputFile.data());
 
  TCanvas* c1 = new TCanvas("c1","",600,600);
  
  int i=0;

  TLegend* leg = new TLegend(0.156,0.763,0.931,0.912);
  leg-> SetNColumns(3);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) {
    
    gtemp[i] = dynamic_cast<TGraphErrors*>(gDirectory->Get(Form("gr_Mean_%s",iterator->second.data())));
    gtemp[i]->SetMarkerColor(1+i);
    gtemp[i]->SetLineColor(1+i);
    gtemp[i]->SetMarkerSize(1);
    gtemp[i]->SetName(Form("gr_Mean_%s",iterator->second.data()));
    if(i==0)
      gtemp[i]->Draw();
    else
      gtemp[i]->Draw("same");
    std::string tempName =gtemp[i]->GetName();
    leg->AddEntry(tempName.data(), iterator->second.data(), "p");
    i++;
    
  }

  
  leg->Draw();
     
  c1->Print("genJetContent.gif");  
  c1->Print("genJetContent.pdf");

}
