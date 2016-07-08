#include <string>
#include "setNCUStyle.C"

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
  

}
