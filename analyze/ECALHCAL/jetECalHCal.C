#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>  
#include <TString.h>
#include <TH2.h>
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

float FWHM(TH1F* hist)
{
  int bin1 = hist->FindFirstBinAbove(hist->GetMaximum()/2);
  int bin2 = hist->FindLastBinAbove(hist->GetMaximum()/2);
  float fwhm = hist->GetBinCenter(bin2) - hist->GetBinCenter(bin1);
  //  return (fwhm/2.36);
  return fwhm;
}

void jetECalHcal(string inputDir, float radius=0.4, bool corr=false){

  vector<float> jeratio_vec;

  const float xmin=0.5;
  const float xmax=1.5;
  TH1F* h_jeratio = new TH1F("h_jeratio", "", 50,xmin,xmax);
  h_jeratio->SetXTitle("E_{jet}/E_{true}");
  if(corr)
    h_jeratio->SetXTitle("E_{jet}^{corr}/E_{true}");

  TH1F* h_jeecalratio = new TH1F("h_jeecalratio", "", 50,0.5,1.5);
  h_jeecalratio->SetXTitle("E_{jet}^{ECAL}/E_{pho}");

  TH1F* h_jehcalratio = new TH1F("h_jehcalratio", "", 50,0,2);
  //  h_jehcalratio->SetXTitle("E_{jet}^{HCAL}/E_{other}");
  h_jehcalratio->SetXTitle("E_{jet}^{HCAL}/(E_{true}-E_{jet}^{ECAL})");
 
  TH2F* h_jeratio_phoratio = new TH2F("h_jeratio_phoratio","",50,0,1,50,0,2);
  h_jeratio_phoratio->SetXTitle("E_{pho}/E_{true}");
  h_jeratio_phoratio->SetYTitle("E_{jet}/E_{true}");
  if(corr)
    h_jeratio_phoratio->SetYTitle("E_{jet}^{corr}/E_{true}");

  TProfile* pf_jeratio_phoratio = new TProfile("pf_jeratio_phoratio","",50,0,1,0,2);
  pf_jeratio_phoratio->SetXTitle("E_{pho}/E_{true}");
  pf_jeratio_phoratio->SetYTitle("<E_{jet}/E_{true}>");
  if(corr)
    pf_jeratio_phoratio->SetYTitle("<E_{jet}^{corr}/E_{true}>");
  


  string inputFile = inputDir + "/radius" + Form("%0.1f",radius)+ "_rawhit.root";
  cout << "opening " << inputFile.data() << endl;
  TreeReader genTree(inputFile.data(),"tGEN_nonu");
  TreeReader caloTree(inputFile.data(),"trawhits2");



  for(Long64_t jEntry=0; jEntry< genTree.GetEntriesFast() ;jEntry++){

    genTree.GetEntry(jEntry);
    caloTree.GetEntry(jEntry);

    Float_t*  gen_je    = genTree.GetPtrFloat("je");
    Float_t*  gen_jeta  = genTree.GetPtrFloat("jeta");
    Float_t*  gen_jphi  = genTree.GetPtrFloat("jphi");
    Float_t*  gen_jepho = genTree.GetPtrFloat("jepho");
    Int_t     gen_njets = genTree.GetInt("njets");

    Float_t*  calo_je       = caloTree.GetPtrFloat("je");
    Float_t*  calo_jeecal   = caloTree.GetPtrFloat("jeecal");
    Float_t*  calo_jehcal   = caloTree.GetPtrFloat("jehcal");
    Float_t*  calo_jeta = caloTree.GetPtrFloat("jeta");
    Float_t*  calo_jphi = caloTree.GetPtrFloat("jphi");
    Int_t     calo_njets= caloTree.GetInt("njets");

    for(unsigned int i=0; i< calo_njets; i++){

      //      cout << jEntry << "\t" << calo_je[i] << "\t" << gen_je[i] << endl;

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
      float ratio = corr? (calo_jeecal[i]+ calo_jehcal[i]/0.7632 )/gen_je[findGenMatch]:
	(calo_jeecal[i]+ calo_jehcal[i])/gen_je[findGenMatch];

      float phoratio = gen_jepho[findGenMatch]/gen_je[findGenMatch];
      if(fabs(calo_je[i]-(calo_jeecal[i]+calo_jehcal[i]))>0.1)
	cout << calo_je[i]-(calo_jeecal[i]+calo_jehcal[i]) << endl;
      h_jeratio->Fill(ratio);

      h_jeecalratio->Fill(calo_jeecal[i]/gen_jepho[findGenMatch]);

      //      h_jehcalratio->Fill(calo_jehcal[i]/(gen_je[findGenMatch]-gen_jepho[findGenMatch]));

      h_jehcalratio->Fill(calo_jehcal[i]/(gen_je[findGenMatch]-calo_jeecal[i]));

      h_jeratio_phoratio->Fill(phoratio,ratio);
      pf_jeratio_phoratio->Fill(phoratio,ratio);

      jeratio_vec.push_back(ratio);

    } // end of loop over calo jets
  } // end loop of tries
      

  std::sort(jeratio_vec.begin(),jeratio_vec.end());
  unsigned int size_temp = jeratio_vec.size();
  const float removal=0.1;
  unsigned int n_temp_to_be_removed = removal*0.5*size_temp;
  // remove the first 5%                         
  jeratio_vec.erase(jeratio_vec.begin(),jeratio_vec.begin()+n_temp_to_be_removed);
  // remove the last 5%
  jeratio_vec.erase(jeratio_vec.end()-n_temp_to_be_removed,jeratio_vec.end());
  double nsamples = jeratio_vec.size();
  double sum = std::accumulate(jeratio_vec.begin(), jeratio_vec.end(), 0.0);
  double mean90 = sum / nsamples;

  double sq_sum = std::inner_product(jeratio_vec.begin(), jeratio_vec.end(), jeratio_vec.begin(), 0.0);
  double RMS90 = std::sqrt(sq_sum / nsamples - mean90 * mean90);


  cout << "Mean90 = " << mean90 << " +- " << RMS90/sqrt(nsamples) << endl;
  cout << "RMS90 = " << RMS90 << " +- " << RMS90/sqrt(2*nsamples) << endl;


  std::string outputFileName = corr? Form("radius%.1f_ECALHCAL_corr.root",radius):
    Form("radius%.1f_ECALHCAL.root",radius);

    TFile* outFile = new TFile(outputFileName.data(),"recreate");
  h_jeratio->Write();
  h_jeecalratio->Write();
  h_jehcalratio->Write();
  h_jeratio_phoratio->Write();
  pf_jeratio_phoratio->Write();
  outFile->Close();

  cout << "Mean = " << h_jeratio->GetMean() 
       << " +- " << h_jeratio->GetMeanError() << endl;

  cout << "RMS = " << h_jeratio->GetRMS() 
       << " +- " << h_jeratio->GetRMSError() << endl;

  cout << "FWHM= " << FWHM(h_jeratio) << endl;

  

}
  
		 

