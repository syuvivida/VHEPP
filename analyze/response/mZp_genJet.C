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


void mZp_genJet(string inputDir, float radius=0.4, int mode=0){

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
      xmin=750;
      xmax=1020;
    } 
  else if(energy=="5")
    {
      xmin=4330;
      xmax=5060;
    }
  else if(energy=="10")
    {
      xmin=6850;
      xmax=10250;
    }
  else if(energy=="20")
    {
      xmin=17780;
      xmax=20200;
    }
  else if(energy=="40")
    {
      xmin=33700;
      xmax=40500;
    }

  TString output=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; echo \"genjet_radius%.1f_${test}.root\"",inputDir.data(),radius));
  cout << "writing " << output.Data() << endl;

  TString output2=gSystem->GetFromPipe(Form("file=%s; test=${file##*/}; echo \"genjet_radius%.1f_${test}\"",inputDir.data(),radius));
  cout << "writing " << output2.Data() << endl;

  TH1F* hmjj = new TH1F("hmjj","",100,xmin,xmax);
  hmjj->SetXTitle("Generator-level M_{jj} [GeV]");
  hmjj->SetTitle(Form("%s-TeV Z'#rightarrow q#bar{q}",energy.Data()));


  TH1F* hooc = new TH1F("hooc","",120,-0.2,1.0);
  hooc->SetXTitle("Fraction of out-of-cone energy");
  hooc->SetTitle(Form("%s-TeV Z'#rightarrow q#bar{q}",energy.Data()));

  //  TreeReader genTree(inputFile.data(),"tGEN_nonu");
  TreeReader genTree(inputFile.data(),"tGEN");
  TreeReader mcTree(inputFile.data(),"tMC");

  for(Long64_t jEntry=0; jEntry< genTree.GetEntriesFast() ;jEntry++){

    genTree.GetEntry(jEntry);
    mcTree.GetEntry(jEntry);

    Float_t*  gen_je   = genTree.GetPtrFloat("je");
    Float_t*  gen_jeta = genTree.GetPtrFloat("jeta");
    Float_t*  gen_jphi = genTree.GetPtrFloat("jphi");
    Float_t*  gen_jpt  = genTree.GetPtrFloat("jpt");
    Float_t*  gen_jm   = genTree.GetPtrFloat("jmass");

    Int_t     gen_njets = genTree.GetInt("njets");

    Float_t*  mc_je    = mcTree.GetPtrFloat("be");
    Float_t*  mc_jeta  = mcTree.GetPtrFloat("beta");
    Float_t*  mc_jphi  = mcTree.GetPtrFloat("bphi");
    Float_t*  mc_jpt   = mcTree.GetPtrFloat("bpt");
    Float_t*  mc_jm    = mcTree.GetPtrFloat("bmass");


    if(gen_njets<2)continue;

    TLorentzVector l4_j[2];
    TLorentzVector l4_q[2];
    int nGoodJets=0;
    for(int i=0; i<2; i++){

       int findGenMatch=-1;
       for(int k=0; k< 2; k++){
	
	
  	float dr = myDeltaR(gen_jeta[k], mc_jeta[i],
 			    gen_jphi[k], mc_jphi[i]);

  	if(dr<TMath::Pi())
  	  {
  	    findGenMatch=k;
  	    break;
  	  }
       }

       if(findGenMatch<0)continue;

      nGoodJets++;
      cout << i << "\t" << findGenMatch << endl;
      l4_j[i].SetPtEtaPhiM(gen_jpt[i],
			   gen_jeta[i],
			   gen_jphi[i],
			   gen_jm[i]);


      l4_q[i].SetPtEtaPhiM(mc_jpt[findGenMatch],
 			   mc_jeta[findGenMatch],
 			   mc_jphi[findGenMatch],
 			   mc_jm[findGenMatch]);


    }
    
    if(nGoodJets<2)continue;

    float mjj = (l4_j[0]+l4_j[1]).M();
    hmjj->Fill(mjj);

    for(int i=0; i<2; i++)
      hooc->Fill((l4_q[i].E()-l4_j[i].E())/l4_q[i].E());

  } // end loop of entries
     
 
  TFile* outFile = new TFile(output.Data(),"recreate");
  hmjj->Write();
  hooc->Write();
  outFile->Close();
  
  TCanvas* c1 = new TCanvas("c1","",500,500);
  setNCUStyle();
  gStyle->SetOptStat(1111111);
  hmjj->Draw();
  c1->Print(Form("%s.pdf",output2.Data()));
  c1->Print(Form("%s.gif",output2.Data()));
  c1->Print(Form("%s.eps",output2.Data()));

}
  
		 

