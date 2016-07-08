#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
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


float FWHM(TH1F* hist)
{
  int bin1 = hist->FindFirstBinAbove(hist->GetMaximum()/2);
  int bin2 = hist->FindLastBinAbove(hist->GetMaximum()/2);
  float fwhm = hist->GetBinCenter(bin2) - hist->GetBinCenter(bin1);
  return (fwhm/2.36);
  //  return fwhm;
}

void jetResponse(string inputDir, float radius=0.4){


  TH1F* h_jeratio = new TH1F("h_jeratio", "", 50,0.6,1.1);
  h_jeratio->SetXTitle("E_{jet}/E_{true}");

  const int nbins=8;
  TH1F* h_jeratiobin[nbins];
  for(unsigned int i=0; i<nbins; i++){
    h_jeratiobin[i] = (TH1F*)h_jeratio->Clone(Form("h_jeratiobin%d",i));
  }

  TH1F* h_je = new TH1F("h_je","",nbins,0,20000);
  TH1F* h_jebin[nbins];
  for(unsigned int i=0; i<nbins; i++){
    h_jebin[i] = (TH1F*)h_je->Clone(Form("h_jebin%d",i));
  }
  string inputFile = inputDir + "/radius" + Form("%0.1f",radius)+ "_of_PanPFA.root";
  TreeReader genTree(inputFile.data(),"tGEN");
  TreeReader caloTree(inputFile.data(),"tcalo");

  for(Long64_t jEntry=0; jEntry< genTree.GetEntriesFast() ;jEntry++){

    genTree.GetEntry(jEntry);
    caloTree.GetEntry(jEntry);

    Float_t*  gen_je = genTree.GetPtrFloat("je");
    Float_t*  gen_jeta = genTree.GetPtrFloat("jeta");
    Float_t*  gen_jphi = genTree.GetPtrFloat("jphi");
    Float_t*  gen_jislep = genTree.GetPtrFloat("jisleptag");

    Float_t*  calo_je = caloTree.GetPtrFloat("je");
    Float_t*  calo_jeta = caloTree.GetPtrFloat("jeta");
    Float_t*  calo_jphi = caloTree.GetPtrFloat("jphi");
    Float_t*  calo_jislep = caloTree.GetPtrFloat("jisleptag");

    for(unsigned int i=0; i<2; i++){

      if(fabs(calo_jeta[i])>1.1)continue;
      if(calo_jislep[i]>1e-6)continue;

      int findGenMatch=-1;
      for(unsigned int k=0; k<2; k++){
	
	
	float dr = sqrt(pow(gen_jeta[k]-calo_jeta[i],2)+
			pow(gen_jphi[k]-calo_jphi[i],2));

	if(dr<0.1)
	  {
	    findGenMatch=k;
	    break;
	  }
      }

      if(findGenMatch<0)continue;
      h_jeratio->Fill(calo_je[i]/gen_je[findGenMatch]);
      int bin = h_je->FindBin(gen_je[findGenMatch]);
      if(bin>nbins)bin=nbins;
      bin=bin-1;
      h_jeratiobin[bin]->Fill(calo_je[i]/gen_je[findGenMatch]);
      h_jebin[bin]->Fill(gen_je[findGenMatch]);
    }
  }
      
  h_jeratio->Draw();
  cout << "RMS = " << h_jeratio->GetRMS() << endl;
  cout << "FWHM= " << FWHM(h_jeratio) << endl;

  TFile* outFile = new TFile(Form("radius%.1f_jetresponse.root",radius),"recreate");
  h_jeratio->Write();
  float x[nbins];
  float y1[nbins];
  float y2[nbins];
  float y3[nbins];
  float y4[nbins];
  float y5[nbins], y5err[nbins];
  float y6[nbins], y6err[nbins];
  float xerr[nbins];

  for(unsigned int i=0; i<nbins; i++){
    h_jeratiobin[i] ->Write();
    h_jebin[i]->Write();
    xerr[i] = h_je->GetBinWidth(i+1)*0.5;
    x[i]  = h_je->GetBinCenter(i+1);
    y1[i] = h_jeratiobin[i]->GetRMS();
    y2[i] = FWHM(h_jeratiobin[i]);
    y3[i] = h_jeratiobin[i]->GetMean();
    y4[i] = h_jeratiobin[i]->GetBinCenter(h_jeratiobin[i]->GetMaximumBin());
    cout << "RMS = " << y1[i] << endl;
    cout << "FWHM= " << y2[i] << endl;
    h_jeratiobin[i]->Fit("gaus");
    TF1 *myfunc = h_jeratiobin[i]->GetFunction("gaus");
    y5[i]    = myfunc->GetParameter(1);
    y5err[i] = myfunc->GetParError(1);
    y6[i]    = myfunc->GetParameter(2);
    y6err[i] = myfunc->GetParError(2);

  }

  TGraph* gr_RMS = new TGraph(nbins,x,y1);
  gr_RMS->SetName("gr_RMS");
  gr_RMS->SetTitle("");
  gr_RMS->GetXaxis()->SetTitle("E_{true} [GeV]");
  gr_RMS->GetYaxis()->SetTitle("RMS of E_{jet}/E_{true}");
  gr_RMS->Draw("ACP");
  gr_RMS->SetMarkerStyle(8);
  gr_RMS->SetMarkerSize(1);
  gr_RMS->GetXaxis()->SetTitleSize(0.05);
  gr_RMS->GetYaxis()->SetTitleSize(0.05);
  gr_RMS->GetXaxis()->SetNdivisions(5);
  //  gr_RMS->GetXaxis()->SetTitleOffset(1.2);
  gr_RMS->GetYaxis()->SetTitleOffset(1.2);
  gr_RMS->GetYaxis()->SetDecimals();

  TGraph* gr_FWHM = new TGraph(nbins,x,y2);
  gr_FWHM->SetName("gr_FWHM");
  gr_FWHM->SetTitle("");
  gr_FWHM->Draw("ACP");
  gr_FWHM->GetXaxis()->SetTitle("E_{true} [GeV]");
  gr_FWHM->GetYaxis()->SetTitle("(FWHM of E_{jet}/E_{true})/2.36");
  gr_FWHM->SetMarkerStyle(8);
  gr_FWHM->SetMarkerSize(1);
  gr_FWHM->GetXaxis()->SetTitleSize(0.05);
  gr_FWHM->GetYaxis()->SetTitleSize(0.05);
  gr_FWHM->GetXaxis()->SetNdivisions(5);
  //  gr_FWHM->GetXaxis()->SetTitleOffset(1.2);
  gr_FWHM->GetYaxis()->SetTitleOffset(1.2);
  gr_FWHM->GetYaxis()->SetDecimals();

  TGraph* gr_Mean = new TGraph(nbins,x,y3);
  gr_Mean->Draw("ACP");
  gr_Mean->SetTitle("");
  gr_Mean->SetName("gr_Mean");
  gr_Mean->GetXaxis()->SetTitle("E_{true} [GeV]");
  gr_Mean->GetYaxis()->SetTitle("Mean of E_{jet}/E_{true}");
  gr_Mean->SetMarkerStyle(8);
  gr_Mean->SetMarkerSize(1);
  gr_Mean->GetXaxis()->SetTitleSize(0.05);
  gr_Mean->GetYaxis()->SetTitleSize(0.05);
  gr_Mean->GetXaxis()->SetNdivisions(5);
//   gr_Mean->GetXaxis()->SetTitleOffset(1.2);
  gr_Mean->GetYaxis()->SetTitleOffset(1.2);
  gr_Mean->GetYaxis()->SetDecimals();



  TGraph* gr_Peak = new TGraph(nbins,x,y4);
  gr_Peak->Draw("ACP");
  gr_Peak->SetTitle("");
  gr_Peak->SetName("gr_Peak");
  gr_Peak->GetXaxis()->SetTitle("E_{true} [GeV]");
  gr_Peak->GetYaxis()->SetTitle("Most probable E_{jet}/E_{true}");
  gr_Peak->SetMarkerStyle(8);
  gr_Peak->SetMarkerSize(1);
  gr_Peak->GetXaxis()->SetTitleSize(0.05);
  gr_Peak->GetYaxis()->SetTitleSize(0.05);
  gr_Peak->GetXaxis()->SetNdivisions(5);
//   gr_Peak->GetXaxis()->SetTitleOffset(1.2);
  gr_Peak->GetYaxis()->SetTitleOffset(1.2);
  gr_Peak->GetYaxis()->SetDecimals();


  TGraphErrors* gr_GausMean = new TGraphErrors(nbins,x,y5,xerr,y5err);
  gr_GausMean->Draw("ACP");
  gr_GausMean->SetTitle("");
  gr_GausMean->SetName("gr_GausMean");
  gr_GausMean->GetXaxis()->SetTitle("E_{true} [GeV]");
  gr_GausMean->GetYaxis()->SetTitle("Gaussian Mean of E_{jet}/E_{true}");
  gr_GausMean->SetMarkerStyle(8);
  gr_GausMean->SetMarkerSize(1);
  gr_GausMean->GetXaxis()->SetTitleSize(0.05);
  gr_GausMean->GetYaxis()->SetTitleSize(0.05);
  gr_GausMean->GetXaxis()->SetNdivisions(5);
  gr_GausMean->GetYaxis()->SetTitleOffset(1.2);
  gr_GausMean->GetYaxis()->SetDecimals();


  TGraphErrors* gr_GausSigma = new TGraphErrors(nbins,x,y6,xerr,y6err);
  gr_GausSigma->Draw("ACP");
  gr_GausSigma->SetTitle("");
  gr_GausSigma->SetName("gr_GausSigma");
  gr_GausSigma->GetXaxis()->SetTitle("E_{true} [GeV]");
  gr_GausSigma->GetYaxis()->SetTitle("Gaussian #sigma of E_{jet}/E_{true}");
  gr_GausSigma->SetMarkerStyle(8);
  gr_GausSigma->SetMarkerSize(1);
  gr_GausSigma->GetXaxis()->SetTitleSize(0.05);
  gr_GausSigma->GetYaxis()->SetTitleSize(0.05);
  gr_GausSigma->GetXaxis()->SetNdivisions(5);
  gr_GausSigma->GetYaxis()->SetTitleOffset(1.2);
  gr_GausSigma->GetYaxis()->SetDecimals();



  gr_RMS->Write();
  gr_FWHM->Write();
  gr_Peak->Write();
  gr_Mean->Write();

  gr_GausMean->Write();
  gr_GausSigma->Write();

  outFile->Close();

  

}
  
		 

