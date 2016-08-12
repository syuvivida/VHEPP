#include <vector>
#include <map>
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
#include <TGraphErrors.h>
#include <TF1.h>

using namespace std;
std::map<unsigned int, std::string> PDGIDs;

const int nparts=2300;
const int lastp = nparts-1;
typedef std::map<unsigned int,std::string>::iterator it_type;


void jetContent(string inputDir, float radius=0.4){

  PDGIDs[22]  = "pho";
  PDGIDs[11]  = "ele";
  PDGIDs[13]  = "muo";
  PDGIDs[211] = "pion";
  PDGIDs[321] = "kaon";
  PDGIDs[130] = "KL";
  PDGIDs[310] = "Ks";
  PDGIDs[2212] = "p";
  PDGIDs[2112] = "n";


  const float xmin=0, xmax=1;
  TH1F* h_jeratio = new TH1F("h_jeratio", "", 50,xmin,xmax);
  h_jeratio->SetXTitle("E_{jet}/E_{true}");

  const int nparts=2300;
  const int nbins=4;
  TH1F* h_je = new TH1F("h_je","",nbins,2500,22500);

  TH1F* h_jeratiobin[nbins][nparts];
  for(unsigned int i=0; i<nbins; i++){
    for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) {
      h_jeratiobin[i][iterator->first] = (TH1F*)h_jeratio->Clone(Form("h_jeratiobin%d_%s",i,iterator->second.data()));
    }
    h_jeratiobin[i][lastp] = (TH1F*)h_jeratio->Clone(Form("h_jeratiobin%d_others",i));
  }

  string inputFile = inputDir + "/radius" + Form("%0.1f",radius)+ "_jetContent.root";
  cout << "opening " << inputFile.data() << endl;
  TreeReader genTree(inputFile.data(),"tGEN_nonu");


  for(Long64_t jEntry=0; jEntry< genTree.GetEntriesFast() ;jEntry++){

    genTree.GetEntry(jEntry);

    Float_t*  gen_je = genTree.GetPtrFloat("trueje");
    Float_t*  gen_jeratio[nparts];
    unsigned int gen_njets = genTree.GetInt("njets");

    for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) {
      gen_jeratio[iterator->first] = genTree.GetPtrFloat(Form("jeratio_%s",iterator->second.data()));
    }
    gen_jeratio[lastp] = genTree.GetPtrFloat("jeratio_others");

    for(unsigned int i=0; i< gen_njets; i++){

      int bin = h_je->FindBin(gen_je[i]);
      if(bin==0)continue;
      if(bin>nbins)bin=nbins;
      bin=bin-1;

      for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) {
	//	cout << gen_jeratio[iterator->first][i]<< endl;
	h_jeratiobin[bin][iterator->first]->Fill(gen_jeratio[iterator->first][i]);
      } // end of loop over particle types
      h_jeratiobin[bin][lastp]->Fill(gen_jeratio[lastp][i]);

    } // end loop of genjets

  } // end loop of entries 



  TFile* outFile = new TFile(Form("radius%.1f_jetContent.root",radius),"recreate");

  float x[nbins];
  float y[nparts][nbins];
  float yerr[nparts][nbins];
  float xerr[nbins];

  for(unsigned int i=0; i<nbins; i++){

    x[i]  = h_je->GetBinCenter(i+1);
    xerr[i] = h_je->GetBinWidth(i+1)*0.5;


    for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) {
      y[iterator->first][i]    = h_jeratiobin[i][iterator->first]->GetMean();
      yerr[iterator->first][i] = h_jeratiobin[i][iterator->first]->GetMeanError();
      cout << i << ", " << iterator->second << ": " << y[iterator->first][i]  << " +- " <<  yerr[iterator->first][i] << std::endl;
    }
    y[lastp][i] = h_jeratiobin[i][lastp]->GetMean();
    yerr[lastp][i] = h_jeratiobin[i][lastp]->GetMeanError();
  }

  TGraphErrors* gr_Mean[nparts];

  int counter=0;
  for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) {

    unsigned int i = iterator->first;
    gr_Mean[i] = new TGraphErrors(nbins,x,y[iterator->first],xerr,yerr[iterator->first]);
    gr_Mean[i]->Draw("ACP");
    gr_Mean[i]->SetTitle(iterator->second.data());
    gr_Mean[i]->SetName(Form("gr_Mean_%s",iterator->second.data()));
    gr_Mean[i]->GetXaxis()->SetTitle("E_{true} [GeV]");
    gr_Mean[i]->GetYaxis()->SetTitle("Mean of E_{jet}/E_{true}");
    gr_Mean[i]->SetMarkerStyle(20+counter);
    gr_Mean[i]->SetMarkerSize(1);
    gr_Mean[i]->GetXaxis()->SetTitleSize(0.05);
    gr_Mean[i]->GetYaxis()->SetTitleSize(0.05);
    gr_Mean[i]->GetXaxis()->SetNdivisions(5);
    gr_Mean[i]->GetYaxis()->SetTitleOffset(1.2);
    gr_Mean[i]->GetYaxis()->SetDecimals();
    gr_Mean[i]->SetMaximum(0.6);
    gr_Mean[i]->SetMinimum(0);
    gr_Mean[i]->Write();
    counter++;
  }

  gr_Mean[lastp] = new TGraphErrors(nbins,x,y[lastp],xerr,yerr[lastp]);
  gr_Mean[lastp]->Draw("ACP");
  gr_Mean[lastp]->SetTitle("Other particles");
  gr_Mean[lastp]->SetName("gr_Mean_others");
  gr_Mean[lastp]->GetXaxis()->SetTitle("E_{true} [GeV]");
  gr_Mean[lastp]->GetYaxis()->SetTitle("Mean of E_{jet}/E_{true}");
  gr_Mean[lastp]->SetMarkerStyle(20+counter);
  gr_Mean[lastp]->SetMarkerSize(1);
  gr_Mean[lastp]->GetXaxis()->SetTitleSize(0.05);
  gr_Mean[lastp]->GetYaxis()->SetTitleSize(0.05);
  gr_Mean[lastp]->GetXaxis()->SetNdivisions(5);
  gr_Mean[lastp]->GetYaxis()->SetTitleOffset(1.2);
  gr_Mean[lastp]->GetYaxis()->SetDecimals();
  gr_Mean[lastp]->SetMaximum(0.6);
  gr_Mean[lastp]->SetMinimum(0);
  gr_Mean[lastp]->Write();

  outFile->Close();

  

}
  
		 

