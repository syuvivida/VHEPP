#include <string>
#include <map>
#include "setNCUStyle.C"
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>

using namespace std;
void plotAngular(string inputDir, float radius=0.4)
{
  setNCUStyle();
  std::string outputfolder = inputDir + "/angular_plots";
  gSystem->mkdir(outputfolder.data());

  std::string name[]={"dphi","deta","dphi1","dphi2","deta1","deta2"};
  const unsigned int nhistos = sizeof(name)/sizeof(name[0]);

  TFile *f = TFile::Open(Form("%s/angular_radius%.1f_rawhit.root",inputDir.data(),radius));
 
  TCanvas* c1 = new TCanvas("c1","",600,600);
  
  TLegend* leg = new TLegend(0.156,0.663,0.931,0.812);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  leg->Draw();    

  for(int i = 0; i < nhistos; i++) {
    
    leg->Clear();
    TProfile* hecal = dynamic_cast<TProfile*>(gDirectory->Get(Form("h_ecal_%s",name[i].data()))); 
    hecal->SetFillColor(4);
    hecal->SetFillStyle(1001);

    TProfile* hhcal = dynamic_cast<TProfile*>(gDirectory->Get(Form("h_total_%s",name[i].data()))); 
    hhcal->SetFillColor(2);
    hhcal->SetFillStyle(1001);
      
    hhcal->Draw("hist");
    hecal->Draw("histsame");
    leg->AddEntry(hecal,"ECAL","f");
    leg->AddEntry(hhcal,"HCAL","f");
    leg->Draw("same");
    c1->Print(Form("%s/h_%s.pdf",outputfolder.data(),name[i].data()));
    c1->Print(Form("%s/h_%s.gif",outputfolder.data(),name[i].data()));

  }

  

     

}
