#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>  
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TProfile.h>
#include "setNCUStyle.C"


using namespace std;

float myDeltaR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta1-eta2;
  float dphi = phi1-phi2;
  if(dphi >= TMath::Pi())
    dphi = dphi-TMath::Pi()*2;
  if(dphi < -TMath::Pi())
    dphi = dphi+TMath::Pi()*2;
  float dr = sqrt(deta*deta+dphi*dphi);
  return dr;
}


void mZp_recoJet(string inputDir, float radius=0.4, int mode=0){

  string decversion;
  if(inputDir.find("rfull009")!=std::string::npos)decversion="rfull009";
  else if(inputDir.find("rfull012")!=std::string::npos)decversion="rfull012";
  else if(inputDir.find("rfull010")!=std::string::npos)decversion="rfull010";

  std::string treeName = mode==0 ? "tcalo":"trawhits2";
  std::string title = mode==0? Form("Anti-kt jet #Delta R = %.1f",radius):
    Form("Simple jet #Delta R = %.1f",radius);
  
  

  string inputFile = inputDir + "/radius" + Form("%0.1f",radius)+ "_response_e2.root";
  cout << "opening " << inputFile.data() << endl;

  TString energy=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; test2=${test%%mumu*}; echo \"${test2##*tev}\"",inputDir.data()));
  cout << "energy = " << energy.Data() << endl;

  float xmin=0,xmax=45000;
  if(energy=="1")
    {
      xmin=0;
      xmax=1200;
    } 
  else if(energy=="5")
    {
      xmin=0;
      xmax=6000;
    }
  else if(energy=="10")
    {
      xmin=0;
      xmax=12000;
    }
  else if(energy=="20")
    {
      xmin=0;
      xmax=22000;
    }
  else if(energy=="40")
    {
      xmin=0;
      xmax=42000;
    }

  TString output=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; echo \"recojet_${test}.root\"",inputDir.data()));
  cout << "writing " << output.Data() << endl;

  TString output2=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; echo \"recojet_${test}\"",inputDir.data()));
  cout << "writing " << output2.Data() << endl;

  TH1F* hmjj = new TH1F("hmjj","",100,xmin,xmax);
  hmjj->SetXTitle("Reconstruction-level M_{jj} [GeV]");
  hmjj->SetTitle(Form("%s-TeV Z'#rightarrow q#bar{q}",energy.Data()));

  TreeReader caloTree(inputFile.data(),"tcalo");
  TreeReader genTree(inputFile.data(),"tGEN_nonu");


  for(Long64_t jEntry=0; jEntry< caloTree.GetEntriesFast() ;jEntry++){

    caloTree.GetEntry(jEntry);
    genTree.GetEntry(jEntry);
    
    Float_t*  gen_je = genTree.GetPtrFloat("je");
    Float_t*  gen_jeta = genTree.GetPtrFloat("jeta");
    Float_t*  gen_jphi = genTree.GetPtrFloat("jphi");
    Float_t*  gen_jpt = genTree.GetPtrFloat("jpt");
    Float_t*  gen_jm = genTree.GetPtrFloat("jmass");
    vector<bool> &gen_jislep = *((vector<bool>*) genTree.GetPtr("jisleptag"));
    Int_t     gen_njets = genTree.GetInt("njets");

    Float_t*  calo_je = caloTree.GetPtrFloat("je");
    Float_t*  calo_jeta = caloTree.GetPtrFloat("jeta");
    Float_t*  calo_jphi = caloTree.GetPtrFloat("jphi");
    Float_t*  calo_jpt = caloTree.GetPtrFloat("jpt");
    Float_t*  calo_jm = caloTree.GetPtrFloat("jmass");
    vector<bool> &calo_jislep = *((vector<bool>*) caloTree.GetPtr("jisleptag"));
    Int_t     calo_njets    = caloTree.GetInt("njets");

    if(calo_njets<2)continue;

    int nGoodJets=0;
    for(int i=0; i< calo_njets; i++){
      
      if(fabs(calo_jeta[i])>1.1)continue;
      if(calo_jislep[i])continue;

      int findGenMatch=-1;
      for(int k=0; k< gen_njets; k++){
	
	
 	float dr = myDeltaR(gen_jeta[k], calo_jeta[i],
			    gen_jphi[k], calo_jphi[i]);

 	if(dr<0.1)
 	  {
 	    findGenMatch=k;
 	    break;
 	  }
      }

      if(findGenMatch<0)continue;

      nGoodJets++;
    }

    if(nGoodJets<2)continue;

    TLorentzVector l4_q[2];

    for(int i=0; i<2; i++){
    
      l4_q[i].SetPtEtaPhiM(
			   calo_jpt[i],
			   calo_jeta[i],
			   calo_jphi[i],
			   calo_jm[i]
			   );

    }
    
    float mjj = (l4_q[0]+l4_q[1]).M();
    hmjj->Fill(mjj);

  } // end loop of entries
     
 
  TFile* outFile = new TFile(output.Data(),"recreate");
  hmjj->Write();
  outFile->Close();
  
  TCanvas* c1 = new TCanvas("c1","",500,500);
  setNCUStyle();
  gStyle->SetOptStat(1111111);
  hmjj->Draw();
  c1->Print(Form("%s.pdf",output2.Data()));
  c1->Print(Form("%s.gif",output2.Data()));
  c1->Print(Form("%s.eps",output2.Data()));

}
  
		 

