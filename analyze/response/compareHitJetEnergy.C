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

using namespace std;
const float removal=0.1;

float FWHM(TH1F* hist)
{
  int bin1 = hist->FindFirstBinAbove(hist->GetMaximum()/2);
  int bin2 = hist->FindLastBinAbove(hist->GetMaximum()/2);
  float fwhm = hist->GetBinCenter(bin2) - hist->GetBinCenter(bin1);
  return (fwhm/2.36);
  //  return fwhm;
}

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

void compareHitJetEnergy(string inputDir, string title, bool debug=false, float emin=-999, float emax=-999, float radius=0.4){
  setNCUStyle(true);
//   const float xmin=0.6;
//   const float xmax=1.1;
  const float xmin=0;
  const float xmax=0.08;
  TH1F* h_jeratio = new TH1F("h_jeratio", "", 50,xmin,xmax);
  h_jeratio->SetTitle(title.data());
  h_jeratio->SetXTitle("E_{jet}^{hit1}/E_{jet}^{hit2}");  
  
  string inputFile1 = inputDir + "/radius" + Form("%0.1f",radius)+ "_rawhit_new.root";
  string inputFile2 = inputDir + "/radius" + Form("%0.1f",radius)+ "_rawhit.root";
  cout << "opening " << inputFile1.data() << endl;
  cout << "opening " << inputFile2.data() << endl;
  TreeReader genTree(inputFile1.data(),"tGEN_nonu");
  TreeReader caloTree1(inputFile1.data(),"trawhits2");
  TreeReader caloTree2(inputFile2.data(),"trawhits2");

  for(Long64_t jEntry=0; jEntry< genTree.GetEntriesFast() ;jEntry++){

    genTree.GetEntry(jEntry);
    caloTree1.GetEntry(jEntry);
    caloTree2.GetEntry(jEntry);

    if(jEntry>5 && debug)break;

    Float_t*  gen_je    = genTree.GetPtrFloat("je");
    Float_t*  gen_jeta  = genTree.GetPtrFloat("jeta");
    Float_t*  gen_jphi  = genTree.GetPtrFloat("jphi");
    Int_t     gen_njets = genTree.GetInt("njets");

    Float_t*  calo_je1    = caloTree1.GetPtrFloat("je");
    Float_t*  calo_jeta1  = caloTree1.GetPtrFloat("jeta");
    Float_t*  calo_jphi1  = caloTree1.GetPtrFloat("jphi");
    Int_t     calo_njets1 = caloTree1.GetInt("njets");

    Float_t*  calo_je2    = caloTree2.GetPtrFloat("je");
    Float_t*  calo_jeta2  = caloTree2.GetPtrFloat("jeta");
    Float_t*  calo_jphi2  = caloTree2.GetPtrFloat("jphi");
    Int_t     calo_njets2 = caloTree2.GetInt("njets");


    for(int i=0; i< calo_njets1; i++){

      if(fabs(calo_jeta1[i])>1.1)continue;
      if( (gen_je[i] < emin && emin>-1) || (gen_je[i] > emax && emax>-1) )continue;
      int findGenMatch=-1;
      for(int k=0; k< gen_njets; k++){
	
	
	float dr = myDeltaR(gen_jeta[k],calo_jeta1[i],
			    gen_jphi[k],calo_jphi1[i]);

	if(dr<0.1)
	  {
	    findGenMatch=k;
	    break;
	  }
      } // end of loop over genjets

      if(findGenMatch<0)continue;
      if(debug)
	cout << "findGenMatch = " << findGenMatch << endl;
      int findChMatch=-1;
      for(int m=0; m< calo_njets2; m++){
	
	
 	float dr = myDeltaR(gen_jeta[findGenMatch],calo_jeta2[m],
			    gen_jphi[findGenMatch],calo_jphi2[m]);

 	if(dr<0.1)
 	  {
 	    findChMatch=m;
 	    break;
 	  }
      } // end loop over second hit-jets
    
      if(findChMatch < 0 || findGenMatch<0)continue;
      if(debug)
	cout << "findChMatch = " << findChMatch << endl;
      float ratio=calo_je1[i]/calo_je2[findChMatch];
      if(debug)
	cout << "ratio = " << ratio << endl;
      if(ratio<xmax && ratio>xmin)
	{
	  h_jeratio->Fill(ratio);
	}
    } // end of loop over calo jets
  } // end loop of tries
      

  TCanvas* c1 = new TCanvas("c1","",500,500);
  gStyle->SetOptStat(1111);
  h_jeratio->Draw();

}
  
		 

