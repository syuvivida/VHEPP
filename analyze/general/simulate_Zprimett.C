// quick macro to compute pt of Higgs for a X->Yh decay given X mass and Y mass
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>
#include <TH1.h>
#include <TROOT.h>
#include <vector>
#include <algorithm> // std::min_element
#include <iterator>  // std::begin, std::end


using namespace std;
const double MW = 80.379;
const double Mu = 0.0022;
const double Md = 0.0047;
const double Mb = 4.18;
const double Mt = 173.0;

void simulate_Zprimett(double MZprime, const unsigned int ntrials=20000, bool extreme=false, bool debug=false)
{

  // momentum of top in the rest frame of Z' (in the lepton collider, it's also the lab frame)
  double ptop= sqrt(MZprime*MZprime-4*Mt*Mt)*0.5;
  TH1F* hdR   = new TH1F("hdR","",100,0,10);
  TH1F* heta  = new TH1F("heta","",100,-2.0,2.0);

  TH1F* hangle = new TH1F("hangle","",100,0,TMath::Pi());
  TH1F* htheta = (TH1F*)hangle->Clone("htheta");
  TH1F* hphi = new TH1F("hphi","",100,-TMath::Pi(),TMath::Pi());
  TH1F* hdecayAngle = (TH1F*)hangle->Clone("hdecayAngle");
  TH1F* hrotation = (TH1F*)hangle->Clone("hrotation");
  TH1F* hopen = (TH1F*)hangle->Clone("hopen");
  TH1F* hsmallest = (TH1F*)hangle->Clone("hsmallest");

  TRandom3* gRandom3=new TRandom3();
  for(unsigned int i=0; i<ntrials; i++)
    {
      // first randomize top momentum direction
      
      double eta = -2.0 + gRandom3->Rndm()*4;
      heta->Fill(eta);
      
      double theta = TMath::ATan(TMath::Exp(-eta))*2.0;
      htheta->Fill(theta);
      
      double phi = -TMath::Pi() + gRandom3->Rndm()*2.0*TMath::Pi();
      hphi->Fill(phi);
      
      double pt = ptop*TMath::Sin(theta);
      TLorentzVector l4_top;
      l4_top.SetPtEtaPhiM(pt,eta,phi,Mt);
      if(debug)
	l4_top.Print();
      
      // get unit vector of top momentum
      TVector3 unit_ptop = l4_top.Vect().Unit();
      if(debug)unit_ptop.Print();

      // get unit orthogonal vector of W
      TVector3 orth_ptop = l4_top.Vect().Orthogonal().Unit();
      if(debug)orth_ptop.Print();

      double E_mother = sqrt(ptop*ptop + Mt*Mt);

      double beta = l4_top.Beta();

      double gamma = l4_top.Gamma();
  
      // momentum of W and b at the rest frame of top

      double p= sqrt( (Mt*Mt-(MW+Mb)*(MW+Mb))*
		      (Mt*Mt-(MW-Mb)*(MW-Mb)))*0.5/Mt; 
  
      // energy of W at the rest frame of top
      double E1= sqrt(p*p+MW*MW);

      // energy of b at the rest frame of top
      double E2= sqrt(p*p+Mb*Mb);

      double decayAngle = extreme? 0.5*TMath::Pi(): gRandom3->Rndm()*TMath::Pi();
      hdecayAngle->Fill(decayAngle);

      // momentum of particle W in the lab frame
      double px1 = gamma*(beta*E1+p*TMath::Cos(decayAngle));
      double py1 = p*TMath::Sin(decayAngle);
      double p1 = sqrt(px1*px1+py1*py1);

      TVector3 v3_p1 = px1*unit_ptop + py1*orth_ptop;
      if(debug)v3_p1.Print();

      // momentum of particle b in the lab frame
      double px2 = gamma*(beta*E2-p*TMath::Cos(decayAngle));
      double py2 = -p*TMath::Sin(decayAngle);
      double p2 = sqrt(px2*px2+py2*py2);

      TVector3 v3_p2 = px2*unit_ptop + py2*orth_ptop;
      if(debug)v3_p2.Print();
      
      if(debug)cout << "Original opening angle = " << v3_p1.Angle(v3_p2) << endl;

      double sep = TMath::ACos((px1*px2+py1*py2)/p1/p2);
      if(debug)cout << "sep = " << sep << endl;

      // now do the rotation
      double rotation = extreme? 0: gRandom3->Rndm()*TMath::Pi()*2.0;
      hrotation->Fill(rotation);
      v3_p1.Rotate(rotation,unit_ptop);
      v3_p2.Rotate(rotation,unit_ptop);
      if(debug)v3_p1.Print();
      if(debug)v3_p2.Print();
      if(debug)cout << "After rotation, opening angle = " << v3_p1.Angle(v3_p2) << endl;
      hopen->Fill(v3_p1.Angle(v3_p2));

      TLorentzVector l4_W,l4_b;
      l4_W.SetPxPyPzE(v3_p1.X(),
		      v3_p1.Y(),
		      v3_p1.Z(),
		      sqrt(MW*MW+v3_p1.Mag2()));
      if(debug)l4_W.Print();

      l4_b.SetPxPyPzE(v3_p2.X(),
		      v3_p2.Y(),
		      v3_p2.Z(),
		      sqrt(Mb*Mb+v3_p2.Mag2()));
      if(debug)l4_b.Print();

      if(debug)cout << "Angular separation between W and b is " << l4_W.Angle(l4_b.Vect()) << endl;
      

      // now check the decay of W
      // get unit vector of W
      TVector3 unit_pW = l4_W.Vect().Unit();
      if(debug)unit_pW.Print();

      // get unit orthogonal vector of W
      TVector3 orth_pW = l4_W.Vect().Orthogonal().Unit();
      if(debug)orth_pW.Print();

      double beta_W = l4_W.Beta();

      double gamma_W = l4_W.Gamma();
  
      // momentum of daughters at the rest frame of W

      double pquark= sqrt( (MW*MW-(Mu+Md)*(Mu+Md))*
			   (MW*MW-(Mu-Md)*(Mu-Md)))*0.5/MW; 
  
      // energy of daughter 1 at the rest frame of W
      double E1_q= sqrt(pquark*pquark+Mu*Mu);

      // energy of daughter 2 at the rest frame of W
      double E2_q= sqrt(pquark*pquark+Md*Md);

      double decayAngle_W = extreme? 0.5*TMath::Pi():gRandom3->Rndm()*TMath::Pi();

      // momentum of particle 1 in the lab frame
      double px1_q = gamma_W*(beta_W*E1_q+pquark*TMath::Cos(decayAngle_W));
      double py1_q = pquark*TMath::Sin(decayAngle_W);
      double p1_q = sqrt(px1_q*px1_q+py1_q*py1_q);

      TVector3 v3_p1_q = px1_q*unit_pW + py1_q*orth_pW;
      if(debug)v3_p1_q.Print();

      // momentum of particle 2 in the lab frame
      double px2_q = gamma_W*(beta_W*E2_q-pquark*TMath::Cos(decayAngle_W));
      double py2_q = -pquark*TMath::Sin(decayAngle_W);
      double p2_q = sqrt(px2_q*px2_q+py2_q*py2_q);

      TVector3 v3_p2_q = px2_q*unit_pW + py2_q*orth_pW;
      if(debug)v3_p2_q.Print();
      
      if(debug)cout << "Original opening angle = " << v3_p1_q.Angle(v3_p2_q) << endl;

      double sep_q = TMath::ACos((px1_q*px2_q+py1_q*py2_q)/p1_q/p2_q);
      if(debug)cout << "sep_q = " << sep_q << endl;


      // now do the rotation
      double rotation_q = extreme? 0: gRandom3->Rndm()*TMath::Pi()*2.0;
      v3_p1_q.Rotate(rotation_q,unit_pW);
      v3_p2_q.Rotate(rotation_q,unit_pW);
      if(debug)v3_p1_q.Print();
      if(debug)v3_p2_q.Print();
      if(debug)cout << "After rotation, opening angle = " << v3_p1_q.Angle(v3_p2_q) << endl;

      TLorentzVector l4_q1,l4_q2;
      l4_q1.SetPxPyPzE(v3_p1_q.X(),
		       v3_p1_q.Y(),
		       v3_p1_q.Z(),
		       sqrt(Mu*Mu+v3_p1_q.Mag2()));

      l4_q2.SetPxPyPzE(v3_p2_q.X(),
		       v3_p2_q.Y(),
		       v3_p2_q.Z(),
		       sqrt(Md*Md+v3_p2_q.Mag2()));
      

      // check the smallest deltaR between b, u, d
      std::vector<double> dRs = {
	l4_b.DeltaR(l4_q1),
	l4_b.DeltaR(l4_q2),
	l4_q1.DeltaR(l4_q2)
      };

      double mindR = -9999;
      auto result = std::min_element(std::begin(dRs), std::end(dRs));
      if (std::end(dRs)!=result)
	mindR = *result;

      if(debug)cout << "dR = " << mindR << endl;
      hdR->Fill(mindR);


      // check the smallest opening angle between b, u, d
      std::vector<double> dthetas = {
	l4_b.Angle(l4_q1.Vect()),
	l4_b.Angle(l4_q2.Vect()),
	l4_q1.Angle(l4_q2.Vect())
      };

      double minAngle = -9999;
      auto result2 = std::min_element(std::begin(dthetas), std::end(dthetas));
      if (std::end(dthetas)!=result2)
	minAngle = *result2;

      if(debug)cout << "smallest angle = " << minAngle << endl;
      hsmallest->Fill(minAngle);



    }

  //  hdR->Draw();
  
  cout << "average smallest deltaR, assuming a flat decay angle distribution";
  cout << " is " << hdR->GetMean() << " +- " << hdR->GetMeanError() << endl;

  cout << "average smallest opening angle, assuming a flat decay angle distribution";
  cout << " is " << hsmallest->GetMean() << " +- " << hsmallest->GetMeanError() << endl;

  
}



