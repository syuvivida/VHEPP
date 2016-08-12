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

void singleParticleResponse(string inputDir, float radius=0.4){

  const float xmin=0.6;
  const float xmax=1.1;
  TH1F* h_jeratio = new TH1F("h_jeratio", "", 50,xmin,xmax);
  h_jeratio->SetXTitle("E_{jet}/E_{true}");  
  
  string inputFile = inputDir + "/radius" + Form("%0.1f",radius)+ "_of_PanPFA.root";
  cout << "opening " << inputFile.data() << endl;
  TreeReader genTree(inputFile.data(),"tGEN_nonu");
  TreeReader caloTree(inputFile.data(),"tcalo");


  vector<float> jeratio_vec;

  for(Long64_t jEntry=0; jEntry< genTree.GetEntriesFast() ;jEntry++){

    genTree.GetEntry(jEntry);
    caloTree.GetEntry(jEntry);

    Float_t*  gen_je    = genTree.GetPtrFloat("je");
    Float_t*  gen_jeta  = genTree.GetPtrFloat("jeta");
    Float_t*  gen_jphi  = genTree.GetPtrFloat("jphi");
    Int_t     gen_njets = genTree.GetInt("njets");

    Float_t*  calo_je    = caloTree.GetPtrFloat("je");
    Float_t*  calo_jeta  = caloTree.GetPtrFloat("jeta");
    Float_t*  calo_jphi  = caloTree.GetPtrFloat("jphi");
    Int_t     calo_njets = caloTree.GetInt("njets");

    for(unsigned int i=0; i< calo_njets; i++){

      if(fabs(calo_jeta[i])>1.1)continue;

      int findGenMatch=-1;
      for(unsigned int k=0; k< gen_njets; k++){
	
	
	float dr = sqrt(pow(gen_jeta[k]-calo_jeta[i],2)+
			pow(gen_jphi[k]-calo_jphi[i],2));

	if(dr<0.1)
	  {
	    findGenMatch=k;
	    break;
	  }
      }

      if(findGenMatch<0)continue;
      float ratio=calo_je[i]/gen_je[findGenMatch];
      h_jeratio->Fill(ratio);
      jeratio_vec.push_back(ratio);
    } // end of loop over calo jets
  } // end loop of tries
      

  // sort the vector of jeratio first
  std::sort(jeratio_vec.begin(),jeratio_vec.end());

  unsigned int size_temp = jeratio_vec.size();
  unsigned int n_temp_to_be_removed = removal*0.5*size_temp;
  if(size_temp>0){
    // remove the first 5%
    jeratio_vec.erase(jeratio_vec.begin(),jeratio_vec.begin()+n_temp_to_be_removed);
    // remove the last 5%
    jeratio_vec.erase(jeratio_vec.end()-n_temp_to_be_removed,jeratio_vec.end());
  }


  h_jeratio->Draw();
  cout << "RMS = " << h_jeratio->GetRMS() << endl;
  cout << "FWHM= " << FWHM(h_jeratio) << endl;


  // getting mean90 and rms90
  double nsamples = jeratio_vec.size();
  double sum = std::accumulate(jeratio_vec.begin(), jeratio_vec.end(), 0.0);
  double mean90 = sum / nsamples;
    
  double sq_sum = std::inner_product(jeratio_vec.begin(), jeratio_vec.end(), jeratio_vec.begin(), 0.0);
  double RMS90 = std::sqrt(sq_sum / nsamples - mean90 * mean90);

  std::cout << "Mean90 = " << mean90 << " +- " << RMS90/sqrt(nsamples) << std::endl;
  std::cout << "RMS90  = " << RMS90  << " +- " << RMS90/sqrt(2*nsamples) << std::endl;


}
  
		 

