// quick macro to compute pt of Higgs for a X->Yh decay given X mass and Y mass
#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>
#include <TH1.h>
#include <TROOT.h>

using namespace std;

void simulate_randomAngle(const unsigned int nSteps=2300, const unsigned int nTrials=20000, bool debug=false, const double degree=0.016)
{
  const double eachAngle = degree*TMath::DegToRad();
  if(debug)cout << "each deflection in radians = " << eachAngle << endl;
  TRandom3* gRandom3=new TRandom3();
  gRandom3->SetSeed();
  TH1D* h1 = new TH1D("h1","",100,0,TMath::Pi()*2.0);
  
  for(unsigned int i=0; i<nTrials; i++)
    {
      TVector3 firstOriginal(1,0,0);
      TVector3 original = firstOriginal;
      TVector3 thisStep = original;
      TVector3 ortho = original.Orthogonal();
      for(unsigned int j=0; j<nSteps; j++)
	{
	  double angle = gRandom3->Rndm()>0.5? eachAngle: -eachAngle;
	  thisStep.Rotate(angle,ortho);
	  if(debug){
	    cout << "Before rotation: " << endl;
	    thisStep.Print();
	  }
	  // now do the rotation
	  double rotation = -TMath::Pi()+gRandom3->Rndm()*TMath::Pi()*2.0;
	  if(debug)cout << "rotation = " << rotation << endl;
	  thisStep.Rotate(rotation,original);
	  if(debug)
	    {
	      cout << "================================================================" << endl;
	      cout << "original vector: " << endl;
	      original.Print();
	      cout << "thisStep vector: " << endl;
	      thisStep.Print();
	      cout << "Angular separation from last step = " << thisStep.Angle(original) << endl; 
	    }
	  original = thisStep;
	  ortho = original.Orthogonal();
	} // end of loop over steps

      double sep = thisStep.Angle(firstOriginal);
      if(debug)
	cout << "Angular sepration from the very orginal = " << sep << endl;

      h1->Fill(sep);
    } // end of loop over trials

  h1->Draw();
}



