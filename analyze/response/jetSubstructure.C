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
#include <TProfile.h>

using namespace std;
const float removal=0.1;

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


void jetSubstructure(string inputDir, float radius=0.4, int mode=0){

  string decversion;
  if(inputDir.find("rfull009")!=std::string::npos)decversion="rfull009";
  else if(inputDir.find("rfull012")!=std::string::npos)decversion="rfull012";
  else if(inputDir.find("rfull010")!=std::string::npos)decversion="rfull010";

  std::string treeName = mode==0 ? "tcalo":"trawhits2";
  std::string title = mode==0? Form("Anti-kt jet #Delta R = %.1f",radius):
    Form("Simple jet #Delta R = %.1f",radius);



  const std::string histonames[]={
    "mass_trim",
    "mass_mmdt",
    "mass_prun",
    "mass_sdb2",
    "mass_sdm1",
    "tau21_b1",
    "c2_b1",
    "d2_b1",
    "n2_b1",
    "tau21_b1_mmdt",
    "c2_b1_mmdt",
    "d2_b1_mmdt",
    "n2_b1_mmdt"
  };



  //  int ntemp_histos=sizeof(histonames)/sizeof(histonames[0]);
  const int nhistos=13;//ntemp_histos;

  std::cout << "Total number of histograms is " << nhistos << std::endl;
  const float xmax[nhistos]={400,400,400,400,400,1,1,20,1,1,1,20,1};

  TH1F* h_sub[nhistos];
  const float xmin=0;
  const int nbins=100;
  for(int ih=0; ih < nhistos; ih++)
    {
      h_sub[ih] = new TH1F(Form("h_%s",histonames[ih].data()), histonames[ih].data(), nbins,xmin,xmax[ih]);
      float binwidth = (xmax[ih]-xmin)/(float)nbins;
      h_sub[ih]->SetXTitle(histonames[ih].data());
      h_sub[ih]->SetYTitle(Form("Number of jets per %.2f",binwidth));
    }

  string inputFile = inputDir + "/radius" + Form("%0.1f",radius)+ "_response_e2.root";
  cout << "opening " << inputFile.data() << endl;
  TreeReader genTree(inputFile.data(),"tGEN_nonu");
  TreeReader caloTree(inputFile.data(),treeName.data());


  for(Long64_t jEntry=0; jEntry< genTree.GetEntriesFast() ;jEntry++){

    genTree.GetEntry(jEntry);
    caloTree.GetEntry(jEntry);

    Float_t*  gen_je = genTree.GetPtrFloat("je");
    Float_t*  gen_jeta = genTree.GetPtrFloat("jeta");
    Float_t*  gen_jphi = genTree.GetPtrFloat("jphi");
    Int_t     gen_njets = genTree.GetInt("njets");

    Float_t*  calo_je = caloTree.GetPtrFloat("je");
    Float_t*  calo_jeta = caloTree.GetPtrFloat("jeta");
    Float_t*  calo_jphi = caloTree.GetPtrFloat("jphi");
    Int_t     calo_njets    = caloTree.GetInt("njets");

    Float_t*  sub[nhistos];

    for(int ih=0; ih < nhistos; ih++)
      sub[ih] = caloTree.GetPtrFloat(Form("j_%s",histonames[ih].data()));


    for(int i=0; i< calo_njets; i++){

      
      if(fabs(calo_jeta[i])>1.1)continue;

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

      for(int ih=0; ih < nhistos; ih++)    
	h_sub[ih]->Fill(sub[ih][i]);


    } // end of loop over calo jets
  } // end loop of entries
     
  

  string outputFile = inputDir + "/radius" + Form("%0.1f",radius)+ "_jetsubstructure_" + treeName + ".root";
  cout << "writing output to " << outputFile.data() << endl;

  TFile* outFile = new TFile(outputFile.data(),"recreate");

  for(int ih=0;ih<nhistos;ih++)
    h_sub[ih]->Write();
  outFile->Close();
  

}
  
		 

