// quick macro to compute pt of Higgs for a X->Yh decay given X mass and Y mass

float momentum_symmetry(float M, float m1, float m2, float pMom, float decayAngle,bool debug=true)
{
  float E_mother = sqrt(pMom*pMom + M*M);

  float beta = pMom/E_mother;
  if(debug)cout << "beta = " << beta << endl;

  float gamma = 1/sqrt(1-beta*beta);
  if(debug)cout << "gamma = " << gamma << endl;
  
  // momentum of daughter at the rest frame of mother
  float p= sqrt( (M*M-(m1+m2)*(m1+m2))*
		  (M*M-(m1-m2)*(m1-m2)))*0.5/M;  
  if(debug)cout << "p = " << p << endl;
  
  // energy of daughter 1 at the rest frame of mother
  float E1= sqrt(p*p+m1*m1);
  if(debug)cout << "E1 = " << E1 << endl;

  // energy of daughter 2 at the rest frame of mother
  float E2= sqrt(p*p+m2*m2);
  if(debug)cout << "E2 = " << E2 << endl;

  // check if it's possible for decayAngle 0 to have smaller separation than 90 degree, i.e. if mother particle has too much boost

  bool tooboosted = (beta>p/E2 || beta>p/E1)? true:false;
  if(debug)
  cout << "Mother particle is too boosted? " << tooboosted << endl;
  
  float theta = TMath::DegToRad()*decayAngle;

  // momentum of particle 1 in the lab frame
  float px1 = gamma*(beta*E1+p*TMath::Cos(theta));
  float py1 = p*TMath::Sin(theta);
  float p1 = sqrt(px1*px1+py1*py1);
  if(debug)
  cout << "Momentum of particle 1 in the lab frame is " << p1 << endl;
  
  float tan1= py1/px1;
  float theta1 = TMath::ATan(tan1);
  if(debug)
  cout << "Angle of particle 1 from the mother in the lab frame is "
       << theta1 << endl;

  // momentum of particle 2 in the lab frame
  float px2 = gamma*(beta*E2-p*TMath::Cos(theta));
  float py2 = -p*TMath::Sin(theta);
  float p2 = sqrt(px2*px2+py2*py2);
  if(debug)
  cout << "Momentum of particle 2 in the lab frame is " << p2 << endl;

  float tan2= py2/px2;
  float theta2 = TMath::ATan(tan2);
  if(debug)
  cout << "Angle of particle 2 from the mother in the lab frame is "
       << theta2 << endl;

  float balence = p1/p2;

  if(debug)
    cout << "momentum ratio = " << balence << endl;

  

  return balence;

}

Double_t momSym(Double_t* x, Double_t* par)
{

  Double_t theta=x[0];
  Double_t M = par[0];
  Double_t m1 = par[1];
  Double_t m2 = par[2];
  Double_t pMom = par[3];

  Double_t balence = momentum_symmetry(M,m1,m2,pMom,theta,false);
  return balence;
}
