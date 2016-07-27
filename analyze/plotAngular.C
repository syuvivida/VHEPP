#include <string>
#include <map>
#include "setNCUStyle.C"
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>

using namespace std;
void plotAngular(string inputDir, float radius=0.4, int mode=0)
{
  setNCUStyle();
  std::string outputfolder = inputDir + "/angular_plots";
  gSystem->mkdir(outputfolder.data());

  std::string name[]={"dphi","deta","dphi1","dphi2","deta1","deta2","signdR"};
  const unsigned int nhistos = sizeof(name)/sizeof(name[0]);

  TFile *f = TFile::Open(Form("%s/angular_radius%.1f_rawhit.root",inputDir.data(),radius));
 
  TCanvas* c1 = new TCanvas("c1","",600,600);
  
  TLegend* leg = new TLegend(0.128,0.743,0.928,0.894);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  leg->Draw();    

  for(int i = 0; i < nhistos; i++) {
    
    leg->Clear();
    TProfile* hhcal = dynamic_cast<TProfile*>(gDirectory->Get(Form("h_hcal_%s",name[i].data()))); 
    hhcal->SetFillColor(2);
    hhcal->SetFillStyle(1001);
    if(mode!=0)
      hhcal->GetXaxis()->SetNdivisions(5);

    TProfile* htotal = dynamic_cast<TProfile*>(gDirectory->Get(Form("h_total_%s",name[i].data()))); 
    htotal->SetFillColor(4);
    htotal->SetFillStyle(1001);
    if(mode!=0)
      htotal->GetXaxis()->SetNdivisions(5);
    htotal->GetYaxis()->SetTitleOffset(1.2);

    htotal->Draw("hist");
    hhcal->Draw("histsame");
    leg->AddEntry(hhcal,"HCAL","f");
    leg->AddEntry(htotal,"ECAL","f");
    leg->Draw("same");
    c1->Print(Form("%s/h_%s.pdf",outputfolder.data(),name[i].data()));
    c1->Print(Form("%s/h_%s.gif",outputfolder.data(),name[i].data()));

  }

  

     

}
