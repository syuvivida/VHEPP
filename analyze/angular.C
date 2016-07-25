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
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TF1.h>

using namespace std;


void angular(string inputDir, float radius=0.4){

  const int nbins=400;
  const float xmin=-0.2;
  const float xmax= 0.2;


  std::string name[]={"dphi","deta","dphi1","dphi2","deta1","deta2"};
  std::string xtitle[]={"#Delta #phi", "#Delta #eta",
			"#Delta #phi_{q1}",
			"#Delta #phi_{q2}",
			"#Delta #eta_{q1}",
			"#Delta #eta_{q2}"};
  
  std::string ytitle=Form("Entries per %.3f",(xmax-xmin)/(float)nbins);

  TH1F* h_dphi = new TH1F("h_dphi","",nbins,xmin,xmax);
  TProfile* h_pf = new TProfile("h_pf","",nbins,xmin,xmax);
  h_pf->SetYTitle(ytitle.data());

  const unsigned int nhistos = sizeof(name)/sizeof(name[0]);
  TProfile* h_ecal[nhistos]; 
  TProfile* h_total[nhistos];

  for(unsigned int i=0; i<nhistos; i++){

    h_ecal[i]  = (TProfile*)h_pf->Clone(Form("h_ecal_%s",name[i].data()));
    h_ecal[i] -> SetXTitle(xtitle[i].data());
    h_total[i] = (TProfile*)h_pf->Clone(Form("h_total_%s",name[i].data()));
    h_total[i] -> SetXTitle(xtitle[i].data());
  }

  string inputFile = inputDir + "/radius" + Form("%0.1f",radius)+ "_rawhit.root";

  TreeReader caloTree(inputFile.data(),"hittuple");
  //  TreeReader caloTree(inputFile.data(),"fjtuple");


  for(Long64_t jEntry=0; jEntry< caloTree.GetEntriesFast() ;jEntry++){
 
    caloTree.GetEntry(jEntry);


    Float_t*  calo_je = caloTree.GetPtrFloat("hite");
    Float_t*  calo_delta[nhistos];
    for(unsigned int i=0; i < nhistos; i++)
      calo_delta[i]= caloTree.GetPtrFloat(Form("hit%s",name[i].data()));

    Int_t*    calo_dec  = caloTree.GetPtrInt("hitdec");
    Int_t*    calo_index  = caloTree.GetPtrInt("hitindex");
    Float_t*  calo_geneta = caloTree.GetPtrFloat("genMEta");
    Int_t     nhits     = caloTree.GetInt("nhits");


    for(unsigned int ih=0; ih< nhistos; ih++){

      float cale[2][2][nbins]={0};
      float total[2][nbins]={0};
    
      for(unsigned int i=0; i< nhits; i++){

	if(abs(calo_geneta[i])>1.1)continue;
	int bin = h_dphi->FindBin(calo_delta[ih][i]);
	if(bin==0)continue;
	if(bin>nbins)continue;
	bin=bin-1;
      
	if(calo_index[i]<0 || calo_index[i]>1)continue;
	if(calo_dec[i]<0 || calo_dec[i]>1)continue;

	cale[calo_index[i]][calo_dec[i]][bin] += calo_je[i];
	total[calo_index[i]][bin] += calo_je[i];
      } // end of loop over calo hits
    
    for(unsigned int ib=0; ib<nbins; ib++)
      {
	float bincenter = h_dphi->GetBinCenter(ib+1);
	  for(unsigned int ip=0; ip<2; ip++){
	    if(cale[ip][0][ib]>1e-6)
	    h_ecal[ih]->Fill(bincenter,cale[ip][0][ib]);
	  if(total[ip][ib]>1e-6)
	    h_total[ih]->Fill(bincenter,total[ip][ib]);
	
	  } // end loop over mothers
		     
      } // end loop of histogram bins

    } // end loop of histos

  }  // end loop of entries  


  TFile* outFile = new TFile(Form("%s/angular_radius%.1f_rawhit.root",inputDir.data(),radius),"recreate");

  for(unsigned int ih=0; ih<nhistos; ih++){

    h_ecal[ih]->Write();
    h_total[ih]->Write();
  }    
  outFile->Close();

  

}
  
		 

