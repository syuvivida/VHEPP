// quick macro to compute pt of Higgs for a X->Yh decay given X mass and Y mass
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>
#include <TH1.h>

using namespace std;
const double MW = 80.379;
const double Mu = 0.0022;
const double Md = 0.0047;

void simulate_ZprimeWW(double MZprime, const unsigned int ntrials=20000, bool debug=false)
{

  // momentum of W in the rest frame of Z' (in the lepton collider, it's also the lab frame)
  double pW= sqrt(MZprime*MZprime-4*MW*MW)*0.5;
  TH1F* hdR   = new TH1F("hdR","",100,0,10);
  TH1F* heta  = new TH1F("heta","",100,-2.0,2.0);
  TH1F* hangle = new TH1F("hangle","",100,0,TMath::Pi());
  TH1F* htheta = (TH1F*)hangle->Clone("htheta");
  TH1F* hphi = new TH1F("hphi","",100,-TMath::Pi(),TMath::Pi());
  TH1F* hdecayAngle = (TH1F*)hangle->Clone("hdecayAngle");
  TH1F* hrotation = (TH1F*)hangle->Clone("hrotation");
  TH1F* hopen = (TH1F*)hangle->Clone("hopen");

  TRandom3* gRandom3=new TRandom3();
  for(unsigned int i=0; i<ntrials; i++)
    {
      // first randomize W momentum direction
      
      double eta = -2.0 + gRandom3->Rndm()*4;
      heta->Fill(eta);
      
      double theta = TMath::ATan(TMath::Exp(-eta))*2.0;
      htheta->Fill(theta);
      
      double phi = -TMath::Pi() + gRandom3->Rndm()*2.0*TMath::Pi();
      hphi->Fill(phi);
      
      double pt = pW*TMath::Sin(theta);
      TLorentzVector l4_W;
      l4_W.SetPtEtaPhiM(pt,eta,phi,MW);
      if(debug)
	l4_W.Print();
      
      // get unit vector of W
      TVector3 unit_pW = l4_W.Vect().Unit();
      if(debug)unit_pW.Print();

      // get unit orthogonal vector of W
      TVector3 orth_pW = l4_W.Vect().Orthogonal().Unit();
      if(debug)orth_pW.Print();

      double E_mother = sqrt(pW*pW + MW*MW);

      double beta = pW/E_mother;

      double gamma = 1/sqrt(1-beta*beta);
  
      // momentum of daughter 1 at the rest frame of W

      double p= sqrt( (MW*MW-(Mu+Md)*(Mu+Md))*
		      (MW*MW-(Mu-Md)*(Mu-Md)))*0.5/MW; 
  
      // energy of daughter 1 at the rest frame of W
      double E1= sqrt(p*p+Mu*Mu);

      // energy of daughter 2 at the rest frame of W
      double E2= sqrt(p*p+Md*Md);

      double decayAngle = gRandom3->Rndm()*TMath::Pi();
      hdecayAngle->Fill(decayAngle);

      // momentum of particle 1 in the lab frame
      double px1 = gamma*(beta*E1+p*TMath::Cos(decayAngle));
      double py1 = p*TMath::Sin(decayAngle);
      double p1 = sqrt(px1*px1+py1*py1);

      TVector3 v3_p1 = px1*unit_pW + py1*orth_pW;
      if(debug)v3_p1.Print();

      double tan1= py1/px1;
      if(debug)cout << "tan1 = " << tan1 << endl;
      double theta1 = TMath::ATan(tan1);

      // momentum of particle 2 in the lab frame
      double px2 = gamma*(beta*E2-p*TMath::Cos(decayAngle));
      double py2 = -p*TMath::Sin(decayAngle);
      double p2 = sqrt(px2*px2+py2*py2);

      TVector3 v3_p2 = px2*unit_pW + py2*orth_pW;
      if(debug)v3_p2.Print();
      
      double tan2= py2/px2;
      if(debug)cout << "tan2 = " << tan2 << endl;
      double theta2 = TMath::ATan(tan2);

      if(debug)cout << "Original opening angle = " << v3_p1.Angle(v3_p2) << endl;

      double sep = TMath::ACos((px1*px2+py1*py2)/p1/p2);
      if(debug)cout << "sep = " << sep << endl;


      // now do the rotation
      double rotation = gRandom3->Rndm()*TMath::Pi()*2.0;
      hrotation->Fill(rotation);
      v3_p1.Rotate(rotation,unit_pW);
      v3_p2.Rotate(rotation,unit_pW);
      if(debug)v3_p1.Print();
      if(debug)v3_p2.Print();
      if(debug)cout << "After rotation, opening angle = " << v3_p1.Angle(v3_p2) << endl;
      hopen->Fill(v3_p1.Angle(v3_p2));

      TLorentzVector l4_q1,l4_q2;
      l4_q1.SetPxPyPzE(v3_p1.X(),
		       v3_p1.Y(),
		       v3_p1.Z(),
		       sqrt(Mu*Mu+v3_p1.Mag2()));

      l4_q2.SetPxPyPzE(v3_p2.X(),
		       v3_p2.Y(),
		       v3_p2.Z(),
		       sqrt(Md*Md+v3_p2.Mag2()));
      

      double dR = l4_q1.DeltaR(l4_q2);
      if(debug)cout << "dR = " << dR << endl;
      hdR->Fill(dR);
    }

  //  hdR->Draw();
  
  cout << "average deltaR, assuming a flat decay angle distribution";
  cout << " is " << hdR->GetMean() << " +- " << hdR->GetMeanError() << endl;
  
  
}



