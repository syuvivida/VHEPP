// quick macro to compute pt of Higgs for a X->Yh decay given X mass and Y mass
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>
#include <TH1.h>

using namespace std;

double angular_separation(double M, double m1, double m2, double pMom, double decayAngle,bool debug=true)
{
  double E_mother = sqrt(pMom*pMom + M*M);

  double beta = pMom/E_mother;
  if(debug)cout << "beta = " << beta << endl;

  double gamma = 1/sqrt(1-beta*beta);
  if(debug)cout << "gamma = " << gamma << endl;
  
  // momentum of daughter at the rest frame of mother
  double p= sqrt( (M*M-(m1+m2)*(m1+m2))*
		  (M*M-(m1-m2)*(m1-m2)))*0.5/M;  
  if(debug)cout << "p = " << p << endl;
  
  // energy of daughter 1 at the rest frame of mother
  double E1= sqrt(p*p+m1*m1);
  if(debug)cout << "E1 = " << E1 << endl;

  // energy of daughter 2 at the rest frame of mother
  double E2= sqrt(p*p+m2*m2);
  if(debug)cout << "E2 = " << E2 << endl;

  // check if it's possible for decayAngle 0 to have smaller separation than 90 degree, i.e. if mother particle has too much boost

  bool tooboosted = (beta>p/E2 || beta>p/E1)? true:false;
  if(debug)
  cout << "Mother particle is too boosted? " << tooboosted << endl;
  
  double theta = TMath::DegToRad()*decayAngle;

  // momentum of particle 1 in the lab frame
  double px1 = gamma*(beta*E1+p*TMath::Cos(theta));
  double py1 = p*TMath::Sin(theta);
  double p1 = sqrt(px1*px1+py1*py1);
  if(debug)
  cout << "Momentum of particle 1 in the lab frame is " << p1 << endl;
  
  double tan1= py1/px1;
  double theta1 = TMath::ATan(tan1);
  if(debug)
  cout << "Angle of particle 1 from the mother in the lab frame is "
       << theta1 << endl;

  // momentum of particle 2 in the lab frame
  double px2 = gamma*(beta*E2-p*TMath::Cos(theta));
  double py2 = -p*TMath::Sin(theta);
  double p2 = sqrt(px2*px2+py2*py2);
  if(debug)
  cout << "Momentum of particle 2 in the lab frame is " << p2 << endl;

  double tan2= py2/px2;
  double theta2 = TMath::ATan(tan2);
  if(debug)
  cout << "Angle of particle 2 from the mother in the lab frame is "
       << theta2 << endl;

  double sep = TMath::ACos((px1*px2+py1*py2)/p1/p2);

  if(debug)
  cout << "angular separation = " << sep << " radians" << endl;

  

  return sep;

}

void computeAverage(double M, double m1, double m2, double pMom)
{
  TRandom3* gRandom3=new TRandom3();
  double sep_total=0;  
  TH1F* hangle=new TH1F("hangle","",180,0,180);
  const unsigned int ntrials=20000;
  for(unsigned int i=0; i<ntrials; i++)
    {
      double decayAngle = gRandom3->Rndm()*180;
      hangle->Fill(decayAngle);
      sep_total += abs(angular_separation(M,m1,m2,pMom,decayAngle,false));   
    }
 
  sep_total /= ntrials;

  cout << "average angular separation, assuming a flat decay angle distribution";
  cout << " is " << sep_total << " radians" << endl;
  
  
}



Double_t angularSep(Double_t* x, Double_t* par)
{

  Double_t theta=x[0];
  Double_t M = par[0];
  Double_t m1 = par[1];
  Double_t m2 = par[2];
  Double_t pMom = par[3];

  Double_t sep = angular_separation(M,m1,m2,pMom,theta,false);
  return sep;
}
