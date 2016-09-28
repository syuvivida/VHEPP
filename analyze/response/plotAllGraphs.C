#include <string>
#include "setNCUStyle.C"
#include "TF1.h"

using namespace std;
void plotAllGraphs()
{
  setNCUStyle(true);
 
  TCanvas* c1 = new TCanvas("c1","",600,600);
  
  gr_RMS->Draw();
  c1->Print("RMS.pdf");
  c1->Print("RMS.gif");
  c1->Clear();

  gr_FWHM->Draw();
  c1->Print("FWHM.pdf");
  c1->Print("FWHM.gif");
  c1->Clear();

  gr_Peak->Draw();
  c1->Print("Peak.pdf");
  c1->Print("Peak.gif");
  c1->Clear();

  gr_Mean->Draw();
  c1->Print("Mean.pdf");
  c1->Print("Mean.gif");
  c1->Clear();

  gr_GausMean->Draw();
  c1->Print("GausMean.pdf");
  c1->Print("GausMean.gif");
  c1->Clear();

  gr_GausSigma->Draw();
  c1->Print("GausSigma.pdf");
  c1->Print("GausSigma.gif");
  c1->Clear();


  gStyle->SetOptFit(11111);

  //  TF1* f1 = new TF1("f1","[0]/sqrt(x)+[1]");
  TF1* f1 = new TF1("f1","sqrt(pow([0]/sqrt(x),2)+pow([1],2))");
  f1->SetNpx(2500);
  h_RMS->Fit("f1");
  c1->Print("RMSFit.pdf");
  c1->Print("RMSFit.gif");
  c1->Clear();

  h_RMS90->Fit("f1");
  c1->Print("RMS90Fit.pdf");
  c1->Print("RMS90Fit.gif");
  c1->Clear();


  TF1* f2 = new TF1("f2","[0]*(1-[1]*exp(-x/[2]))");
  f2->SetNpx(2500);
  f2->SetParameters(1,1,5000);
  gStyle->SetStatX(0.90);
  gStyle->SetStatY(0.37);


  h_Mean->Fit("f2");
  c1->Print("MeanFit.pdf");
  c1->Print("MeanFit.gif");
  c1->Clear();

  h_Mean90->Fit("f2");
  c1->Print("Mean90Fit.pdf");
  c1->Print("Mean90Fit.gif");
  c1->Clear();


}
