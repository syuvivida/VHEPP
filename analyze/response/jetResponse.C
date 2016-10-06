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

void jetResponse(string inputDir, float radius=0.4, int mode=0){

  std::string treeName = mode==0 ? "tcalo":"trawhits2";
  std::string title = mode==0? Form("Anti-kt jet #Delta R = %.1f",radius):
    Form("Simple jet #Delta R = %.1f",radius);
  const float xmin=0.6;
  const float xmax=1.1;
  TH1F* h_jeratio = new TH1F("h_jeratio", "", 50,xmin,xmax);
  h_jeratio->SetXTitle("E_{jet}/E_{true}");
  
  const int nbins=5;
  TH1F* h_jeratiobin[nbins];
  for(int i=0; i<nbins; i++)
    h_jeratiobin[i] = (TH1F*)h_jeratio->Clone(Form("h_jeratiobin%d",i));
  
  const float xlows[nbins+1]={0,600,2600,5500,12000,22000};

  float sumEntries[nbins]={0,0,0};
  float numEntries[nbins]={0,0,0};
  float xcenter[nbins]={0,0,0};

  TH1F* h_je = new TH1F("h_je","",nbins,xlows);
  TH1F* h_jebin[nbins];
  for(int i=0; i<nbins; i++){
    h_jebin[i] = (TH1F*)h_je->Clone(Form("h_jebin%d",i));
  }
  string inputFile = inputDir + "/radius" + Form("%0.1f",radius)+ "_rawhit.root";
  cout << "opening " << inputFile.data() << endl;
  TreeReader genTree(inputFile.data(),"tGEN_nonu");
  TreeReader caloTree(inputFile.data(),treeName.data());


  vector<float> jeratio_vec[nbins];

  for(Long64_t jEntry=0; jEntry< genTree.GetEntriesFast() ;jEntry++){

    genTree.GetEntry(jEntry);
    caloTree.GetEntry(jEntry);

    Float_t*  gen_je = genTree.GetPtrFloat("je");
    Float_t*  gen_jeta = genTree.GetPtrFloat("jeta");
    Float_t*  gen_jphi = genTree.GetPtrFloat("jphi");
    vector<bool> &gen_jislep = *((vector<bool>*) genTree.GetPtr("jisleptag"));


    Int_t     gen_njets = genTree.GetInt("njets");

    Float_t*  calo_je = caloTree.GetPtrFloat("je");
    Float_t*  calo_jeta = caloTree.GetPtrFloat("jeta");
    Float_t*  calo_jphi = caloTree.GetPtrFloat("jphi");
    vector<bool> &calo_jislep = *((vector<bool>*) caloTree.GetPtr("jisleptag"));
    Int_t     calo_njets    = caloTree.GetInt("njets");

    for(int i=0; i< calo_njets; i++){

      
      if(fabs(calo_jeta[i])>1.1)continue;
      if(calo_jislep[i])continue;

      int findGenMatch=-1;
      for(int k=0; k< gen_njets; k++){
	
	
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
      int bin = h_je->FindBin(gen_je[findGenMatch]);
      cout << "e = " << gen_je[findGenMatch] << " bin = " << bin << endl;

      if(bin==0)continue;
      if(bin>nbins)bin=nbins;

      bin=bin-1;
      h_jebin[bin]->Fill(gen_je[findGenMatch]);
      h_jeratiobin[bin]->Fill(ratio);

      if(ratio>=xmin && ratio<=xmax)
	{
	  jeratio_vec[bin].push_back(ratio);
	  sumEntries[bin] += gen_je[findGenMatch];
	  numEntries[bin] += 1;
	}

    } // end of loop over calo jets
  } // end loop of entries
     
  
  // sort the vector of jeratio first
  for(int i=0; i<nbins; i++)
    std::sort(jeratio_vec[i].begin(),jeratio_vec[i].end());

  // erase elements in the vectors
  for(int i=0; i<nbins; i++)
    {
      unsigned int size_temp = jeratio_vec[i].size();
      unsigned int n_temp_to_be_removed = removal*0.5*size_temp;
      // remove the first 5%
      jeratio_vec[i].erase(jeratio_vec[i].begin(),jeratio_vec[i].begin()+n_temp_to_be_removed);
      // remove the last 5%
      jeratio_vec[i].erase(jeratio_vec[i].end()-n_temp_to_be_removed,jeratio_vec[i].end());
    }


  h_jeratio->Draw();
  cout << "RMS = " << h_jeratio->GetRMS() << endl;
  cout << "FWHM= " << FWHM(h_jeratio) << endl;

  TFile* outFile = new TFile(Form("rfull009_radius%.1f_jetresponse_%s.root",radius,treeName.data()),"recreate");
  h_jeratio->Write();
  float x[nbins];
  float y1[nbins], y1err[nbins];
  float y2[nbins];
  float y3[nbins], y3err[nbins];
  float y4[nbins];
  float y5[nbins], y5err[nbins];
  float y6[nbins], y6err[nbins];
  float xerr[nbins];
  float yrms90[nbins],yrms90err[nbins];
  float ymean90[nbins],ymean90err[nbins];
  float yrmsmean90[nbins],yrmsmean90err[nbins];

  TH1F* h_Mean = (TH1F*)h_je->Clone("h_Mean");
  h_Mean->Reset();
  TH1F* h_RMS = (TH1F*)h_je->Clone("h_RMS");
  h_RMS->Reset();


  TH1F* h_Mean90 = (TH1F*)h_je->Clone("h_Mean90");
  h_Mean90->Reset();
  TH1F* h_RMS90 = (TH1F*)h_je->Clone("h_RMS90");
  h_RMS90->Reset();



  
  for(int i=0; i<nbins; i++){
    h_jeratiobin[i]  ->Write();
    h_jebin[i]->Write();

    cout << "entries of bin " << i << "= " << h_jeratiobin[i]->GetEntries() << endl;

    xcenter[i] = numEntries[i]>0? sumEntries[i]/numEntries[i]: h_je->GetBinCenter(i+1);
    cout << "xcenter of bin " << i << " = " << xcenter[i] << endl;
    x[i]  = xcenter[i];
    xerr[i] = TMath::Min(fabs(h_je->GetBinLowEdge(i+1)-xcenter[i]),
			 fabs(h_je->GetBinLowEdge(i+2)-xcenter[i]));
    
    bool emptyBin=false;
    if(h_jeratiobin[i]->GetEntries()<50)emptyBin=true;


    y1[i] = emptyBin? -1: h_jeratiobin[i]->GetRMS();
    y1err[i] = emptyBin? 0.001: h_jeratiobin[i]->GetRMSError();

    h_RMS->SetBinContent(i+1,y1[i]);
    h_RMS->SetBinError(i+1,y1err[i]);

    y2[i] = emptyBin? -1: FWHM(h_jeratiobin[i]);

    y3[i] = emptyBin? -1: h_jeratiobin[i]->GetMean();
    y3err[i] = emptyBin? 0.001: h_jeratiobin[i]->GetMeanError();
    h_Mean->SetBinContent(i+1,y3[i]);
    h_Mean->SetBinError(i+1,y3err[i]);

    y4[i] = emptyBin? -1: h_jeratiobin[i]->GetBinCenter(h_jeratiobin[i]->GetMaximumBin());

    cout << "RMS = " << y1[i] << " +- " << y1err[i] << endl;
    cout << "FWHM= " << y2[i] << endl;

    h_jeratiobin[i]->Fit("gaus");
    TF1 *myfunc = h_jeratiobin[i]->GetFunction("gaus");

    y5[i]    = emptyBin? -1: myfunc->GetParameter(1);
    y5err[i] = emptyBin? 0.001: myfunc->GetParError(1);

    y6[i]    = emptyBin? -1: myfunc->GetParameter(2);
    y6err[i] = emptyBin? 0.001: myfunc->GetParError(2);


    // getting mean90 and rms90
    double nsamples = jeratio_vec[i].size();
    double sum = emptyBin? 0: std::accumulate(jeratio_vec[i].begin(), jeratio_vec[i].end(), 0.0);
    double mean90 = emptyBin? -1: sum / nsamples;
    
    double sq_sum = emptyBin? 0: std::inner_product(jeratio_vec[i].begin(), jeratio_vec[i].end(), jeratio_vec[i].begin(), 0.0);
    double RMS90 = emptyBin? -1: std::sqrt(sq_sum / nsamples - mean90 * mean90);

    ymean90[i]    = mean90;
    ymean90err[i] = emptyBin? 0.001: RMS90/sqrt(nsamples);

    h_Mean90->SetBinContent(i+1,mean90);
    h_Mean90->SetBinError(i+1, RMS90/sqrt(nsamples));

    yrms90[i]     = RMS90;
    yrms90err[i]  = emptyBin? 0.001: RMS90/sqrt(2*nsamples);

    h_RMS90->SetBinContent(i+1,RMS90);
    h_RMS90->SetBinError(i+1, RMS90/sqrt(2*nsamples));

    yrmsmean90[i]    = emptyBin? -1: RMS90/mean90;
    yrmsmean90err[i] = emptyBin? 0.001: yrmsmean90[i]*sqrt(pow(yrms90err[i]/yrms90[i],2)+pow(ymean90err[i]/ymean90[i],2));

    cout << "bin " << i << ": RMSMean90 = " << yrmsmean90[i]  << "+-" << yrmsmean90err[i] << endl;
  }

  h_Mean->SetMarkerStyle(8);
  h_Mean->SetMarkerSize(1);
  h_Mean->GetXaxis()->SetTitle("E_{true} [GeV]");
  h_Mean->GetYaxis()->SetTitle("Mean of E_{jet}/E_{true}");
  h_Mean->GetXaxis()->SetNdivisions(5);
  h_Mean->GetYaxis()->SetTitleOffset(1.6);
  h_Mean->GetYaxis()->SetDecimals();
  h_Mean->Write();

  h_RMS->SetMarkerStyle(8);
  h_RMS->SetMarkerSize(1);
  h_RMS->GetXaxis()->SetTitle("E_{true} [GeV]");
  h_RMS->GetYaxis()->SetTitle("RMS of E_{jet}/E_{true}");
  h_RMS->GetXaxis()->SetNdivisions(5);
  h_RMS->GetYaxis()->SetTitleOffset(1.9);
  h_RMS->GetYaxis()->SetDecimals();
  h_RMS->Write();



  h_Mean90->SetMarkerStyle(8);
  h_Mean90->SetMarkerSize(1);
  h_Mean90->GetXaxis()->SetTitle("E_{true} [GeV]");
  h_Mean90->GetYaxis()->SetTitle("Mean^{90} of E_{jet}/E_{true}");
  h_Mean90->GetXaxis()->SetNdivisions(5);
  h_Mean90->GetYaxis()->SetTitleOffset(1.6);
  h_Mean90->GetYaxis()->SetDecimals();
  h_Mean90->Write();

  h_RMS90->SetMarkerStyle(8);
  h_RMS90->SetMarkerSize(1);
  h_RMS90->GetXaxis()->SetTitle("E_{true} [GeV]");
  h_RMS90->GetYaxis()->SetTitle("RMS^{90} of E_{jet}/E_{true}");
  h_RMS90->GetXaxis()->SetNdivisions(5);
  h_RMS90->GetYaxis()->SetTitleOffset(1.9);
  h_RMS90->GetYaxis()->SetDecimals();
  h_RMS90->Write();


  TGraphErrors* gr_RMS90 = new TGraphErrors(nbins,x,yrms90,xerr,yrms90err);
  gr_RMS90->SetName("gr_RMS90");
  gr_RMS90->SetTitle(title.data());
  gr_RMS90->GetXaxis()->SetTitle("E_{true} [GeV]");
  gr_RMS90->GetYaxis()->SetTitle("RMS^{90} of E_{jet}/E_{true}");
  gr_RMS90->Draw("ACP");
  gr_RMS90->SetMarkerStyle(8);
  gr_RMS90->SetMarkerSize(1);
  gr_RMS90->GetXaxis()->SetTitleSize(0.05);
  gr_RMS90->GetYaxis()->SetTitleSize(0.05);
  gr_RMS90->GetXaxis()->SetNdivisions(5);
  gr_RMS90->GetYaxis()->SetTitleOffset(1.2);
  gr_RMS90->GetYaxis()->SetDecimals();



  TGraphErrors* gr_Mean90 = new TGraphErrors(nbins,x,ymean90,xerr,ymean90err);
  gr_Mean90->SetName("gr_Mean90");
  gr_Mean90->SetTitle(title.data());
  gr_Mean90->GetXaxis()->SetTitle("E_{true} [GeV]");
  gr_Mean90->GetYaxis()->SetTitle("Mean^{90} of E_{jet}/E_{true}");
  gr_Mean90->Draw("ACP");
  gr_Mean90->SetMarkerStyle(8);
  gr_Mean90->SetMarkerSize(1);
  gr_Mean90->GetXaxis()->SetTitleSize(0.05);
  gr_Mean90->GetYaxis()->SetTitleSize(0.05);
  gr_Mean90->GetXaxis()->SetNdivisions(5);
  gr_Mean90->GetYaxis()->SetTitleOffset(1.2);
  gr_Mean90->GetYaxis()->SetDecimals();


  TGraphErrors* gr_RMSMean90 = new TGraphErrors(nbins,x,yrmsmean90,xerr,yrmsmean90err);
  gr_RMSMean90->SetName("gr_RMSMean90");
  gr_RMSMean90->SetTitle(title.data());
  gr_RMSMean90->GetXaxis()->SetTitle("E_{true} [GeV]");
  gr_RMSMean90->GetYaxis()->SetTitle("RMS^{90}/Mean^{90} of E_{jet}/E_{true}");
  gr_RMSMean90->Draw("ACP");
  gr_RMSMean90->SetMarkerStyle(8);
  gr_RMSMean90->SetMarkerSize(1);
  gr_RMSMean90->GetXaxis()->SetTitleSize(0.05);
  gr_RMSMean90->GetYaxis()->SetTitleSize(0.05);
  gr_RMSMean90->GetXaxis()->SetNdivisions(5);
  gr_RMSMean90->GetYaxis()->SetTitleOffset(1.2);
  gr_RMSMean90->GetYaxis()->SetDecimals();


  TGraphErrors* gr_RMS = new TGraphErrors(nbins,x,y1,xerr,y1err);
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

  TGraphErrors* gr_Mean = new TGraphErrors(nbins,x,y3,xerr,y3err);
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



  gr_RMS90->Write();
  gr_Mean90->Write();
  gr_RMSMean90->Write();

  gr_RMS->Write();
  gr_FWHM->Write();
  gr_Peak->Write();
  gr_Mean->Write();

  gr_GausMean->Write();
  gr_GausSigma->Write();

  outFile->Close();
  

}
  
		 

