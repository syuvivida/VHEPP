#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <cmath>
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"

#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/contrib/ModifiedMassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/DistanceMeasure.hh"

//#ifdef __MAKECINT__
//#pragma link C++ class vector<float>+;
//#endif

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

/*
 
 TO COMPILE:
 
 export ROOTSYS=~/Desktop/root
 export PATH=$ROOTSYS/bin:$PATH
 export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
 export DYLD_LIBRARY_PATH=$ROOTSYS/lib:$DYLD_LIBRARY_PATH 
 
 c++ -o anaResponse `root-config --glibs --cflags` `/Users/ntran/Research/CMS/PhysicsAnalysis/boost2013/fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -lvectorDict -lEnergyCorrelator anaResponse.cpp
 
 TO RUN:     
 
 ./analysis01 ../madevent_3/Events/pietro_test14_unweighted_events.lhe ../Z2j/run_01_unweighted_events.lhe 
 
 */

//! ========================================================================================

////////////////////-----------------------------------------------
// Global variables
ifstream fin;
ifstream finMC;
ifstream finGEN;
ifstream finGEN_charged;
ifstream finGEN_nonu;
ifstream finGEN_response;
ifstream finGEN_resolution;
ifstream finCalo;
ifstream finTrack;
TRandom3* myrndm;
std::vector<fastjet::PseudoJet> ISRPhotons;
TH1F* heta_status3;
TH1F* heta;
TH1F* heta_after1;
TH1F* heta_after2;
TH1F* heta_PF;
TH1F* hy;

int counter;
int evtCtr;

int njets;
double gen_mZp;
double gen_mWW;

std::vector<float> be;
std::vector<float> bp;
std::vector<float> bpt;
std::vector<float> beta;
std::vector<float> bphi;
std::vector<float> bmass;

std::vector<float> je;
std::vector<float> jpt;
std::vector<float> jp;
std::vector<float> jeta;
std::vector<float> jphi;
std::vector<float> jmass;
std::vector<int> jmultiplicity;
std::vector<bool> jisleptag;
std::vector<float> jmass_sd;
std::vector<float> jmass_pr;
std::vector<float> costheta1;

// a bunch of jet substructure variables
std::vector<float> j_tau1_b1;
std::vector<float> j_tau2_b1;
std::vector<float> j_tau3_b1;
std::vector<float> j_tau1_b2;
std::vector<float> j_tau2_b2;
std::vector<float> j_tau3_b2;
std::vector<float> j_tau21_b2;
std::vector<float> j_tau21_b1;
std::vector<float> j_tau32_b2;
std::vector<float> j_tau32_b1;
std::vector<float> j_zlogz;
std::vector<float> j_c1_b0;
std::vector<float> j_c1_b1;
std::vector<float> j_c1_b2;
std::vector<float> j_c2_b1;
std::vector<float> j_c2_b2;
std::vector<float> j_d2_b1;
std::vector<float> j_d2_b2;
std::vector<float> j_d2_a1_b1;
std::vector<float> j_d2_a1_b2;
std::vector<float> j_n2_b1;
std::vector<float> j_n2_b2;
std::vector<float> j_m2_b1;
std::vector<float> j_m2_b2;

std::vector<float> j_tau1_b1_mmdt;
std::vector<float> j_tau2_b1_mmdt;
std::vector<float> j_tau3_b1_mmdt;
std::vector<float> j_tau1_b2_mmdt;
std::vector<float> j_tau2_b2_mmdt;
std::vector<float> j_tau3_b2_mmdt;
std::vector<float> j_tau21_b2_mmdt;
std::vector<float> j_tau21_b1_mmdt;
std::vector<float> j_tau32_b2_mmdt;
std::vector<float> j_tau32_b1_mmdt;
std::vector<float> j_c1_b0_mmdt;
std::vector<float> j_c1_b1_mmdt;
std::vector<float> j_c1_b2_mmdt;
std::vector<float> j_c2_b1_mmdt;
std::vector<float> j_c2_b2_mmdt;
std::vector<float> j_d2_b1_mmdt;
std::vector<float> j_d2_b2_mmdt;
std::vector<float> j_d2_a1_b1_mmdt;
std::vector<float> j_d2_a1_b2_mmdt;
std::vector<float> j_n2_b1_mmdt;
std::vector<float> j_n2_b2_mmdt;
std::vector<float> j_m2_b1_mmdt;
std::vector<float> j_m2_b2_mmdt;

std::vector<float> j_qjetVol;
std::vector<float> j_mass_trim;
std::vector<float> j_mass_mmdt;
std::vector<float> j_mass_prun;
std::vector<float> j_mass_sdb2;
std::vector<float> j_mass_sdm1;

struct genHiggs
{
  TLorentzVector p4;
  float cosTheta1;
  TLorentzVector p4_b1;
  TLorentzVector p4_b2;
};


genHiggs h_gen[2];

float jetRadius;
const float etamax=99999.0;
//const float etamax=1.1;

////////////////////-----------------------------------------------
void readEvent( std::vector< fastjet::PseudoJet > &allParticles );
void readEventMC( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN_charged( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN_nonu( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN_response( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN_resolution( std::vector< fastjet::PseudoJet > &allParticles );
void readEventCalo( std::vector< fastjet::PseudoJet > &allParticles );
void readEventTrack(  std::vector< fastjet::PseudoJet > &allParticles );

void computeCosTheta1();
void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal, std::vector < fastjet::PseudoJet > MCparticles, bool fillJSub=false);
void analyzeMCEvent(std::vector < fastjet::PseudoJet > MCparticles);
void clearVectors(){
    //  INIT
  njets = 0;


  be.clear();
  bpt.clear();
  bp.clear();    
  beta.clear();
  bphi.clear();        
  bmass.clear();

  je.clear();
  jpt.clear();
  jp.clear();    
  jeta.clear();
  jphi.clear();        
  jmass.clear();
  jmass_sd.clear();        
  jmass_pr.clear();        
  jmultiplicity.clear();
  jisleptag.clear();
  costheta1.clear();

  // a bunch of jet substructure variables

  j_tau1_b1.clear();
  j_tau2_b1.clear();
  j_tau3_b1.clear();    
  j_tau1_b2.clear();
  j_tau2_b2.clear();
  j_tau3_b2.clear();    
  j_tau21_b1.clear();
  j_tau21_b2.clear();
  j_tau32_b1.clear();
  j_tau32_b2.clear();
  j_zlogz.clear();
  j_c1_b0.clear();
  j_c1_b1.clear();
  j_c1_b2.clear();
  j_c2_b1.clear();
  j_c2_b2.clear();
  j_d2_b1.clear();
  j_d2_b2.clear();
  j_d2_a1_b1.clear();
  j_d2_a1_b2.clear();
  j_m2_b1.clear();
  j_m2_b2.clear();
  j_n2_b1.clear();
  j_n2_b2.clear();

  j_tau1_b1_mmdt.clear();
  j_tau2_b1_mmdt.clear();
  j_tau3_b1_mmdt.clear();    
  j_tau1_b2_mmdt.clear();
  j_tau2_b2_mmdt.clear();
  j_tau3_b2_mmdt.clear();    
  j_tau21_b1_mmdt.clear();
  j_tau21_b2_mmdt.clear();
  j_tau32_b1_mmdt.clear();
  j_tau32_b2_mmdt.clear();
  j_c1_b0_mmdt.clear();
  j_c1_b1_mmdt.clear();
  j_c1_b2_mmdt.clear();
  j_c2_b1_mmdt.clear();
  j_c2_b2_mmdt.clear();
  j_d2_b1_mmdt.clear();
  j_d2_b2_mmdt.clear();
  j_d2_a1_b1_mmdt.clear();
  j_d2_a1_b2_mmdt.clear();
  j_m2_b1_mmdt.clear();
  j_m2_b2_mmdt.clear();
  j_n2_b1_mmdt.clear();
  j_n2_b2_mmdt.clear();

  j_qjetVol.clear();
  j_mass_trim.clear();
  j_mass_mmdt.clear();
  j_mass_prun.clear();
  j_mass_sdb2.clear();
  j_mass_sdm1.clear();


}


void computeCosTheta1()
{
  for(int i=0;i<2;i++){
    if(h_gen[i].p4.E()<1e-6)return;
    if(h_gen[i].p4_b1.E()<1e-6)return;
    if(h_gen[i].p4_b2.E()<1e-6)return;
  }
  
  for(int i=0; i<2;i++){

    TLorentzVector thisB(h_gen[i].p4_b1);

    thisB.Boost(-h_gen[i].p4.BoostVector());

    float cos=TMath::Cos(thisB.Vect().Angle(h_gen[i].p4.Vect()));

    h_gen[i].cosTheta1=cos;
  }
  return;

}

////////////////////////////////////////////////////////////////////////////////////////
int main (int argc, char **argv) {
    
    // std::cout << "hello world" << std::endl;
    std::string type = "of_PanPFA";   // type "gg" or "qq"
    std::string inputFolder = argv[1];
    jetRadius = (float)atof(argv[2]);

    std::cout << "jetRadius = " << jetRadius << std::endl;
    
    myrndm = new TRandom3();


    // pick only the highest pt particle
    std::vector < fastjet::PseudoJet > allParticles;
    std::vector < fastjet::PseudoJet > allParticlesCalo; // calo clusters matched to PF candidate
    std::vector < fastjet::PseudoJet > allParticlesTrack; // tracks matched to PF candidate
    std::vector < fastjet::PseudoJet > allParticlesMC; //status3
    std::vector < fastjet::PseudoJet > allParticlesGEN; //status1 
    std::vector < fastjet::PseudoJet > allParticlesGEN_charged; //status1 // only charged
    std::vector < fastjet::PseudoJet > allParticlesGEN_nonu; //status1
    std::vector < fastjet::PseudoJet > allParticlesGEN_response; //gen scaled with single particle response
    std::vector < fastjet::PseudoJet > allParticlesGEN_resolution; //gen scaled with single particle resolution

    char fname[150];
    // PF-jets
    sprintf( fname, "%s/%s.dat", inputFolder.c_str(), type.c_str() );
    fin.open(fname);

    char fnameMC[150];
    sprintf( fnameMC, "%s/of_status3.dat", inputFolder.c_str());
    finMC.open(fnameMC);

    char fnameGEN[150];
    sprintf( fnameGEN, "%s/of_status1.dat", inputFolder.c_str());
    finGEN.open(fnameGEN);
    finGEN_charged.open(fnameGEN);
    finGEN_nonu.open(fnameGEN);
    finGEN_response.open(fnameGEN);
    finGEN_resolution.open(fnameGEN);

    char fnameCalo[150];
    sprintf( fnameCalo, "%s/of_CaloClusters.dat", inputFolder.c_str());
    finCalo.open(fnameCalo);

    char fnameTrack[150];
    sprintf( fnameTrack, "%s/of_Tracks.dat", inputFolder.c_str());
    finTrack.open(fnameTrack);

    char outName[192];
    sprintf( outName, "%s/radius%.1f_response_e2.root", inputFolder.c_str(), jetRadius);
    TFile *f = TFile::Open(outName,"RECREATE");

    heta = new TH1F("heta","",200,-2.5,2.5);
    hy   = (TH1F*)heta->Clone("hy");
    heta_status3 = (TH1F*)heta->Clone("heta_status3");
    heta_after1 = (TH1F*)heta->Clone("heta_after1");
    heta_after2 = (TH1F*)heta->Clone("heta_after2");
    heta_PF     = (TH1F*)heta->Clone("heta_PF");


    TTree *tPFA = new TTree("tPFA","Tree with vectors");
    tPFA->Branch("njets"          , &njets      );
    tPFA->Branch("je"             , &je        );
    tPFA->Branch("jpt"            , &jpt        );
    tPFA->Branch("jp"             , &jp        );    
    tPFA->Branch("jeta"           , &jeta       );
    tPFA->Branch("jphi"           , &jphi       );
    tPFA->Branch("jmass"          , &jmass      );    
    tPFA->Branch("jmass_sd"       , &jmass_sd      );    
    tPFA->Branch("jmass_pr"       , &jmass_pr      );    
    tPFA->Branch("jmultiplicity"  , &jmultiplicity      );    
    tPFA->Branch("jisleptag"      , &jisleptag      );    
    tPFA->Branch("costheta1"      , &costheta1      );    



    TTree *tcalo = new TTree("tcalo","Tree with vectors");
    tcalo->Branch("njets"          , &njets      );
    tcalo->Branch("je"             , &je         );
    tcalo->Branch("jpt"            , &jpt        );
    tcalo->Branch("jp"             , &jp        );      
    tcalo->Branch("jeta"           , &jeta       );
    tcalo->Branch("jphi"           , &jphi       );
    tcalo->Branch("jmass"          , &jmass      );    
    tcalo->Branch("jmass_sd"       , &jmass_sd      );    
    tcalo->Branch("jmass_pr"       , &jmass_pr      );    
    tcalo->Branch("jmultiplicity"  , &jmultiplicity      );    
    tcalo->Branch("jisleptag"      , &jisleptag      );    
    tcalo->Branch("costheta1"      , &costheta1      );    

    tcalo->Branch("j_tau1_b1"        , &j_tau1_b1        );
    tcalo->Branch("j_tau2_b1"        , &j_tau2_b1        );
    tcalo->Branch("j_tau3_b1"        , &j_tau3_b1        );    
    tcalo->Branch("j_tau1_b2"        , &j_tau1_b2        );
    tcalo->Branch("j_tau2_b2"        , &j_tau2_b2        );
    tcalo->Branch("j_tau3_b2"        , &j_tau3_b2        );    
    tcalo->Branch("j_tau21_b1"       , &j_tau21_b1       );
    tcalo->Branch("j_tau21_b2"       , &j_tau21_b2       );
    tcalo->Branch("j_tau21_b1"       , &j_tau21_b1       );
    tcalo->Branch("j_tau21_b2"       , &j_tau21_b2       );
    tcalo->Branch("j_tau32_b1"       , &j_tau32_b1       );
    tcalo->Branch("j_tau32_b2"       , &j_tau32_b2       );
    tcalo->Branch("j_zlogz"          , &j_zlogz          );
    tcalo->Branch("j_c1_b0"          , &j_c1_b0          );
    tcalo->Branch("j_c1_b1"          , &j_c1_b1          );
    tcalo->Branch("j_c1_b2"          , &j_c1_b2          );
    tcalo->Branch("j_c2_b1"          , &j_c2_b1          );
    tcalo->Branch("j_c2_b2"          , &j_c2_b2          );
    tcalo->Branch("j_d2_b1"          , &j_d2_b1          );
    tcalo->Branch("j_d2_b2"          , &j_d2_b2          );
    tcalo->Branch("j_d2_a1_b1"       , &j_d2_a1_b1       );
    tcalo->Branch("j_d2_a1_b2"       , &j_d2_a1_b2       );
    tcalo->Branch("j_m2_b1"          , &j_m2_b1          );
    tcalo->Branch("j_m2_b2"          , &j_m2_b2          );
    tcalo->Branch("j_n2_b1"          , &j_n2_b1          );
    tcalo->Branch("j_n2_b2"          , &j_n2_b2          );

    tcalo->Branch("j_tau1_b1_mmdt"        , &j_tau1_b1_mmdt        );
    tcalo->Branch("j_tau2_b1_mmdt"        , &j_tau2_b1_mmdt        );
    tcalo->Branch("j_tau3_b1_mmdt"        , &j_tau3_b1_mmdt        );    
    tcalo->Branch("j_tau1_b2_mmdt"        , &j_tau1_b2_mmdt        );
    tcalo->Branch("j_tau2_b2_mmdt"        , &j_tau2_b2_mmdt        );
    tcalo->Branch("j_tau3_b2_mmdt"        , &j_tau3_b2_mmdt        );    
    tcalo->Branch("j_tau21_b1_mmdt"       , &j_tau21_b1_mmdt       );
    tcalo->Branch("j_tau21_b2_mmdt"       , &j_tau21_b2_mmdt       );
    tcalo->Branch("j_tau21_b1_mmdt"       , &j_tau21_b1_mmdt       );
    tcalo->Branch("j_tau21_b2_mmdt"       , &j_tau21_b2_mmdt       );
    tcalo->Branch("j_tau32_b1_mmdt"       , &j_tau32_b1_mmdt       );
    tcalo->Branch("j_tau32_b2_mmdt"       , &j_tau32_b2_mmdt       );
    tcalo->Branch("j_c1_b0_mmdt"          , &j_c1_b0_mmdt          );
    tcalo->Branch("j_c1_b1_mmdt"          , &j_c1_b1_mmdt          );
    tcalo->Branch("j_c1_b2_mmdt"          , &j_c1_b2_mmdt          );
    tcalo->Branch("j_c2_b1_mmdt"          , &j_c2_b1_mmdt          );
    tcalo->Branch("j_c2_b2_mmdt"          , &j_c2_b2_mmdt          );
    tcalo->Branch("j_d2_b1_mmdt"          , &j_d2_b1_mmdt          );
    tcalo->Branch("j_d2_b2_mmdt"          , &j_d2_b2_mmdt          );
    tcalo->Branch("j_d2_a1_b1_mmdt"       , &j_d2_a1_b1_mmdt       );
    tcalo->Branch("j_d2_a1_b2_mmdt"       , &j_d2_a1_b2_mmdt       );
    tcalo->Branch("j_m2_b1_mmdt"          , &j_m2_b1_mmdt          );
    tcalo->Branch("j_m2_b2_mmdt"          , &j_m2_b2_mmdt          );
    tcalo->Branch("j_n2_b1_mmdt"          , &j_n2_b1_mmdt          );
    tcalo->Branch("j_n2_b2_mmdt"          , &j_n2_b2_mmdt          );

    tcalo->Branch("j_mass_trim"      , &j_mass_trim      );
    tcalo->Branch("j_mass_mmdt"      , &j_mass_mmdt      );
    tcalo->Branch("j_mass_prun"      , &j_mass_prun      );
    tcalo->Branch("j_mass_sdb2"      , &j_mass_sdb2      );
    tcalo->Branch("j_mass_sdm1"      , &j_mass_sdm1      );


    TTree *ttrack = new TTree("ttrack","Tree with vectors");
    ttrack->Branch("njets"          , &njets      );
    ttrack->Branch("je"             , &je         );
    ttrack->Branch("jpt"            , &jpt        );
    ttrack->Branch("jp"             , &jp        );      
    ttrack->Branch("jeta"           , &jeta       );
    ttrack->Branch("jphi"           , &jphi       );
    ttrack->Branch("jmass"          , &jmass      );    
    ttrack->Branch("jmass_sd"       , &jmass_sd      );    
    ttrack->Branch("jmass_pr"       , &jmass_pr      );    
    ttrack->Branch("jmultiplicity"  , &jmultiplicity      );    
    ttrack->Branch("jisleptag"      , &jisleptag      );    
    ttrack->Branch("costheta1"      , &costheta1      );    


    TTree *tGEN = new TTree("tGEN","Tree with vectors");
    tGEN->Branch("njets"          , &njets      );
    tGEN->Branch("je"             , &je         );
    tGEN->Branch("jpt"            , &jpt        );
    tGEN->Branch("jp"             , &jp        );          
    tGEN->Branch("jeta"           , &jeta       );
    tGEN->Branch("jphi"           , &jphi       );
    tGEN->Branch("jmass"          , &jmass      );    
    tGEN->Branch("jmass_sd"       , &jmass_sd      );    
    tGEN->Branch("jmass_pr"       , &jmass_pr      );    
    tGEN->Branch("jmultiplicity"  , &jmultiplicity      );    
    tGEN->Branch("jisleptag"      , &jisleptag      );        
    tGEN->Branch("costheta1"      , &costheta1      );    

    tGEN->Branch("j_tau1_b1"        , &j_tau1_b1        );
    tGEN->Branch("j_tau2_b1"        , &j_tau2_b1        );
    tGEN->Branch("j_tau3_b1"        , &j_tau3_b1        );    
    tGEN->Branch("j_tau1_b2"        , &j_tau1_b2        );
    tGEN->Branch("j_tau2_b2"        , &j_tau2_b2        );
    tGEN->Branch("j_tau3_b2"        , &j_tau3_b2        );    
    tGEN->Branch("j_tau21_b1"       , &j_tau21_b1       );
    tGEN->Branch("j_tau21_b2"       , &j_tau21_b2       );
    tGEN->Branch("j_tau21_b1"       , &j_tau21_b1       );
    tGEN->Branch("j_tau21_b2"       , &j_tau21_b2       );
    tGEN->Branch("j_tau32_b1"       , &j_tau32_b1       );
    tGEN->Branch("j_tau32_b2"       , &j_tau32_b2       );
    tGEN->Branch("j_zlogz"          , &j_zlogz          );
    tGEN->Branch("j_c1_b0"          , &j_c1_b0          );
    tGEN->Branch("j_c1_b1"          , &j_c1_b1          );
    tGEN->Branch("j_c1_b2"          , &j_c1_b2          );
    tGEN->Branch("j_c2_b1"          , &j_c2_b1          );
    tGEN->Branch("j_c2_b2"          , &j_c2_b2          );
    tGEN->Branch("j_d2_b1"          , &j_d2_b1          );
    tGEN->Branch("j_d2_b2"          , &j_d2_b2          );
    tGEN->Branch("j_d2_a1_b1"       , &j_d2_a1_b1       );
    tGEN->Branch("j_d2_a1_b2"       , &j_d2_a1_b2       );
    tGEN->Branch("j_m2_b1"          , &j_m2_b1          );
    tGEN->Branch("j_m2_b2"          , &j_m2_b2          );
    tGEN->Branch("j_n2_b1"          , &j_n2_b1          );
    tGEN->Branch("j_n2_b2"          , &j_n2_b2          );

    tGEN->Branch("j_tau1_b1_mmdt"        , &j_tau1_b1_mmdt        );
    tGEN->Branch("j_tau2_b1_mmdt"        , &j_tau2_b1_mmdt        );
    tGEN->Branch("j_tau3_b1_mmdt"        , &j_tau3_b1_mmdt        );    
    tGEN->Branch("j_tau1_b2_mmdt"        , &j_tau1_b2_mmdt        );
    tGEN->Branch("j_tau2_b2_mmdt"        , &j_tau2_b2_mmdt        );
    tGEN->Branch("j_tau3_b2_mmdt"        , &j_tau3_b2_mmdt        );    
    tGEN->Branch("j_tau21_b1_mmdt"       , &j_tau21_b1_mmdt       );
    tGEN->Branch("j_tau21_b2_mmdt"       , &j_tau21_b2_mmdt       );
    tGEN->Branch("j_tau21_b1_mmdt"       , &j_tau21_b1_mmdt       );
    tGEN->Branch("j_tau21_b2_mmdt"       , &j_tau21_b2_mmdt       );
    tGEN->Branch("j_tau32_b1_mmdt"       , &j_tau32_b1_mmdt       );
    tGEN->Branch("j_tau32_b2_mmdt"       , &j_tau32_b2_mmdt       );
    tGEN->Branch("j_c1_b0_mmdt"          , &j_c1_b0_mmdt          );
    tGEN->Branch("j_c1_b1_mmdt"          , &j_c1_b1_mmdt          );
    tGEN->Branch("j_c1_b2_mmdt"          , &j_c1_b2_mmdt          );
    tGEN->Branch("j_c2_b1_mmdt"          , &j_c2_b1_mmdt          );
    tGEN->Branch("j_c2_b2_mmdt"          , &j_c2_b2_mmdt          );
    tGEN->Branch("j_d2_b1_mmdt"          , &j_d2_b1_mmdt          );
    tGEN->Branch("j_d2_b2_mmdt"          , &j_d2_b2_mmdt          );
    tGEN->Branch("j_d2_a1_b1_mmdt"       , &j_d2_a1_b1_mmdt       );
    tGEN->Branch("j_d2_a1_b2_mmdt"       , &j_d2_a1_b2_mmdt       );
    tGEN->Branch("j_m2_b1_mmdt"          , &j_m2_b1_mmdt          );
    tGEN->Branch("j_m2_b2_mmdt"          , &j_m2_b2_mmdt          );
    tGEN->Branch("j_n2_b1_mmdt"          , &j_n2_b1_mmdt          );
    tGEN->Branch("j_n2_b2_mmdt"          , &j_n2_b2_mmdt          );

    tGEN->Branch("j_mass_trim"      , &j_mass_trim      );
    tGEN->Branch("j_mass_mmdt"      , &j_mass_mmdt      );
    tGEN->Branch("j_mass_prun"      , &j_mass_prun      );
    tGEN->Branch("j_mass_sdb2"      , &j_mass_sdb2      );
    tGEN->Branch("j_mass_sdm1"      , &j_mass_sdm1      );


    TTree *tGEN_charged = new TTree("tGEN_charged","Tree with vectors");
    tGEN_charged->Branch("njets"          , &njets      );
    tGEN_charged->Branch("je"             , &je         );
    tGEN_charged->Branch("jpt"            , &jpt        );
    tGEN_charged->Branch("jp"             , &jp        );          
    tGEN_charged->Branch("jeta"           , &jeta       );
    tGEN_charged->Branch("jphi"           , &jphi       );
    tGEN_charged->Branch("jmass"          , &jmass      );    
    tGEN_charged->Branch("jmass_sd"       , &jmass_sd      );    
    tGEN_charged->Branch("jmass_pr"       , &jmass_pr      );    
    tGEN_charged->Branch("jmultiplicity"  , &jmultiplicity      );    
    tGEN_charged->Branch("jisleptag"      , &jisleptag      );        
    tGEN_charged->Branch("costheta1"      , &costheta1      );    

    TTree *tGEN_nonu = new TTree("tGEN_nonu","Tree with vectors");
    tGEN_nonu->Branch("njets"          , &njets      );
    tGEN_nonu->Branch("je"             , &je         );
    tGEN_nonu->Branch("jpt"            , &jpt        );
    tGEN_nonu->Branch("jp"             , &jp        );          
    tGEN_nonu->Branch("jeta"           , &jeta       );
    tGEN_nonu->Branch("jphi"           , &jphi       );
    tGEN_nonu->Branch("jmass"          , &jmass      );    
    tGEN_nonu->Branch("jmass_sd"       , &jmass_sd      );    
    tGEN_nonu->Branch("jmass_pr"       , &jmass_pr      );    
    tGEN_nonu->Branch("jmultiplicity"  , &jmultiplicity      );    
    tGEN_nonu->Branch("jisleptag"      , &jisleptag      );        
    tGEN_nonu->Branch("costheta1"      , &costheta1      );    

    TTree *tGEN_response = new TTree("tGEN_response","Tree with vectors");
    tGEN_response->Branch("njets"          , &njets      );
    tGEN_response->Branch("je"             , &je         );
    tGEN_response->Branch("jpt"            , &jpt        );
    tGEN_response->Branch("jp"             , &jp        );          
    tGEN_response->Branch("jeta"           , &jeta       );
    tGEN_response->Branch("jphi"           , &jphi       );
    tGEN_response->Branch("jmass"          , &jmass      );    
    tGEN_response->Branch("jmass_sd"       , &jmass_sd      );    
    tGEN_response->Branch("jmass_pr"       , &jmass_pr      );    
    tGEN_response->Branch("jmultiplicity"  , &jmultiplicity      );    
    tGEN_response->Branch("jisleptag"      , &jisleptag      );        
    tGEN_response->Branch("costheta1"      , &costheta1      );    

    TTree *tGEN_resolution = new TTree("tGEN_resolution","Tree with vectors");
    tGEN_resolution->Branch("njets"          , &njets      );
    tGEN_resolution->Branch("je"             , &je         );
    tGEN_resolution->Branch("jpt"            , &jpt        );
    tGEN_resolution->Branch("jp"             , &jp        );          
    tGEN_resolution->Branch("jeta"           , &jeta       );
    tGEN_resolution->Branch("jphi"           , &jphi       );
    tGEN_resolution->Branch("jmass"          , &jmass      );    
    tGEN_resolution->Branch("jmass_sd"       , &jmass_sd      );    
    tGEN_resolution->Branch("jmass_pr"       , &jmass_pr      );    
    tGEN_resolution->Branch("jmultiplicity"  , &jmultiplicity      );    
    tGEN_resolution->Branch("jisleptag"      , &jisleptag      );        
    tGEN_resolution->Branch("costheta1"      , &costheta1      );    


    TTree *tMC = new TTree("tMC","Tree with vectors");
    tMC->Branch("gen_mZp"          , &gen_mZp      );
    tMC->Branch("gen_mWW"          , &gen_mWW      );
    tMC->Branch("be"               , &be        );
    tMC->Branch("bpt"              , &bpt        );
    tMC->Branch("bp"               , &bp        );    
    tMC->Branch("beta"             , &beta       );
    tMC->Branch("bphi"             , &bphi       );
    tMC->Branch("bmass"            , &bmass      );    

    int ctr = 0;
    counter =0;
    while(true){

      for(int i=0; i<2; i++){
	h_gen[i].p4.SetPxPyPzE(0,0,0,0);
	h_gen[i].p4_b1.SetPxPyPzE(0,0,0,0);
	h_gen[i].p4_b2.SetPxPyPzE(0,0,0,0);
	h_gen[i].cosTheta1 = -999;
      }

        readEvent( allParticles );
        readEventGEN( allParticlesGEN ); 
 	readEventGEN_charged( allParticlesGEN_charged ); 
	readEventGEN_nonu( allParticlesGEN_nonu );
	readEventGEN_response( allParticlesGEN_response);
	readEventGEN_resolution( allParticlesGEN_resolution);

        readEventMC( allParticlesMC );
        readEventCalo( allParticlesCalo );
	readEventTrack( allParticlesTrack);

	computeCosTheta1();

        // std::cout << "size of collection = " << allParticles.size() << "," << allParticlesMC.size() << ", " << allParticlesCalo.size() << std::endl;
        if (ctr > 0){
	  
 	    clearVectors();
 	    analyzeEvent( allParticlesGEN, jetRadius, allParticlesMC, true);
	    tGEN->Fill();
	
 	    clearVectors();
 	    analyzeEvent( allParticlesGEN_nonu, jetRadius, allParticlesMC );
	    tGEN_nonu->Fill();
	    
 	    clearVectors();
 	    analyzeEvent( allParticlesGEN_charged, jetRadius, allParticlesMC );
	    tGEN_charged->Fill();
 
 	    clearVectors();
 	    analyzeEvent( allParticlesGEN_response, jetRadius, allParticlesMC );
	    tGEN_response->Fill();

 	    clearVectors();
 	    analyzeEvent( allParticlesGEN_resolution, jetRadius, allParticlesMC );
	    tGEN_resolution->Fill();

            clearVectors();
            analyzeEvent( allParticles, jetRadius, allParticlesMC );
            tPFA->Fill();

            clearVectors();
            analyzeEvent( allParticlesCalo, jetRadius, allParticlesMC, true);
            tcalo->Fill();    

            clearVectors();
	    analyzeEvent( allParticlesTrack, jetRadius, allParticlesMC );
	    ttrack->Fill();

  
            analyzeMCEvent( allParticlesMC );
            tMC->Fill();
        }

        allParticles.clear();
        allParticlesMC.clear();
        allParticlesGEN.clear();
        allParticlesGEN_charged.clear();
        allParticlesGEN_nonu.clear();
        allParticlesGEN_response.clear();
        allParticlesGEN_resolution.clear();
        allParticlesCalo.clear();
        allParticlesTrack.clear();
	
        ctr++;
	std::cout << "ctr = " << ctr << std::endl;

	if(fin.eof()) break;
	//	if(ctr==10)break;
    }

    f->cd();
    heta->Write();
    hy->Write();
    heta_status3->Write();
    heta_after1->Write();
    heta_after2->Write();
    heta_PF->Write();

    tPFA->Write();
    tcalo->Write();
    ttrack->Write();
    tGEN->Write();
    tGEN_charged->Write();
    tGEN_nonu->Write();
    tGEN_response->Write();
    tGEN_resolution->Write();
    tMC->Write();
    f->Close();

}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
// reading PF jets
void readEvent( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = etamax;
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
    
    int ctr = 0;
    while(true){
        
        fin  >> pdgid >> px >> py >> pz >> e;         
        // std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
    
        if (pdgid == -99){
            return;
        }        
    
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        curPseudoJet.set_user_index(pdgid);
	if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
	}


        if(fin.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventMC( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = etamax;
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
    
    int ctr = 0;
    while(true){
        
        finMC  >> pdgid >> px >> py >> pz >> e;         
        // std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
    
        if (pdgid == -99){
	  return;
        }        
    
	if(pdgid==25 && h_gen[0].p4.E()< 1e-6)
	  h_gen[0].p4.SetPxPyPzE(px,py,pz,e);
	else if(pdgid==25 && h_gen[1].p4.E()< 1e-6)
	  h_gen[1].p4.SetPxPyPzE(px,py,pz,e);

	if(pdgid==5 && h_gen[0].p4_b1.E()<1e-6)
	  h_gen[0].p4_b1.SetPxPyPzE(px,py,pz,e);
	else if(pdgid==-5 && h_gen[0].p4_b2.E()<1e-6)
	  h_gen[0].p4_b2.SetPxPyPzE(px,py,pz,e);

	else if(pdgid==5 && h_gen[1].p4_b1.E()<1e-6)
	  h_gen[1].p4_b1.SetPxPyPzE(px,py,pz,e);
	else if(pdgid==-5 && h_gen[1].p4_b2.E()<1e-6)
	  h_gen[1].p4_b2.SetPxPyPzE(px,py,pz,e);


        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
	TLorentzVector temp_l4;
 	temp_l4.SetPxPyPzE(px,py,pz,e);
	if(abs(pdgid)<=5 && abs(pdgid)>=1)
	  heta_status3->Fill(temp_l4.Eta());

        if(finMC.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventGEN( std::vector< fastjet::PseudoJet > &allParticles ){

    
    // hard-coded! 
    float etaMax = etamax; // to check how good the generator-level jets are
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
    int ctr = 0;
    while(true){
        
        finGEN  >> pdgid >> px >> py >> pz >> e;         
        if (pdgid == -99){
            return;
        }        

        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
	curPseudoJet.set_user_index(pdgid);
	if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
	}        
	
	TLorentzVector temp_l4;
 	temp_l4.SetPxPyPzE(px,py,pz,e);
	heta->Fill(temp_l4.Eta());

        if(finGEN.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventGEN_charged( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = etamax; // to check how good the generator-level jets are
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;

    int ctr = 0;
    while(true){
        
        finGEN_charged >> pdgid >> px >> py >> pz >> e;         
        // std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
    
        if (pdgid == -99){
            return;
        }        
 	int pdg=abs(pdgid);
  	if( pdg== 12 || pdg== 14 || pdg== 16 
  	   || pdg== 130 || pdg == 2112){
        if(finGEN_charged.eof()) break;
 	else
 	  continue;
	}
	//	std::cout << pdg << std::endl;
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
        
        if(finGEN_charged.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}



////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventGEN_nonu( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = etamax; // to check how good the generator-level jets are
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;

    int ctr = 0;
    while(true){
        
        finGEN_nonu >> pdgid >> px >> py >> pz >> e;         
        // std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
    
        if (pdgid == -99){
            return;
        }        
 	int pdg=abs(pdgid);
  	if(pdg== 12 || pdg== 14 || pdg== 16){
	  if(finGEN_nonu.eof()) break;
	  else
	    continue;
	}
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
        
        if(finGEN_nonu.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}



////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventGEN_response( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = etamax; // to check how good the generator-level jets are
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;

    int ctr = 0;
    while(true){
        
        finGEN_response >> pdgid >> px >> py >> pz >> e;         
        if (pdgid == -99){
            return;
        }        
 	int pdg=abs(pdgid);
   	if(pdg== 12 || pdg== 14 || pdg== 16 
	   ){
	  if(finGEN_response.eof()) break;
	  else
	    continue;
	}
	//	std::cout << "pdg = " << pdg << std::endl;
        // fill vector of pseudojets

	
	float RMS = 0;
	float response=1;

	if(e<2 && !(pdg== 130 || pdg==310 || pdg == 2112))response=0;


   	if(response<1e-6){
	  if(finGEN_response.eof()) break;
	  else
	    continue;
	}


	switch(pdg){
	case 13:
	  response=2.0/e;
	  RMS=0;
	  break;
	case 22: case 11:
	  response = 1;
	  RMS=sqrt(pow(0.15/sqrt(e),2) + 0.01*0.01)*response;
	  break;
	default:
	  if(e<16)
	    response = 0.7 + 0.077*log(e);
	  else if(e>=16 && e<2000)
	    response = 0.86 + 0.0202*log(e);
	  else if(e>2000)
	    response = 1.0;
	  RMS=sqrt(pow(0.43/sqrt(e),2) + 0.009*0.009)*response;	
	  break;
	}
	

   	if(response<1e-6){
	  if(finGEN_response.eof()) break;
	  else
	    continue;
	}


	float gauss_smear = myrndm->Gaus(response,RMS);
        fastjet::PseudoJet curPseudoJet( px*gauss_smear, 
					 py*gauss_smear, 
					 pz*gauss_smear, 
					 e*gauss_smear);
	

//    	std::cout << "RMS = " << RMS << std::endl;
//   	std::cout << "response = " << response << std::endl;
//     	std::cout << "response = " << gauss_smear << std::endl;

        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
        
        if(finGEN_response.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}

////////////////////////////////////////////////////////////////////////////////////////
void readEventGEN_resolution( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = etamax; // to check how good the generator-level jets are
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;

    int ctr = 0;
    while(true){        
        finGEN_resolution >> pdgid >> px >> py >> pz >> e;             
        if (pdgid == -99){
            return;
        }        
 	int pdg=abs(pdgid);
  	if(pdg== 12 || pdg== 14 || pdg== 16 
	   || pdg== 130 || pdg==310 || pdg == 2112 
	   ){
	  if(finGEN_resolution.eof()) break;
	  else
	    continue;
	}
        // fill vector of pseudojets

	float RMS = 0;
	float response=1;

	if(e<2)response=0;

   	if(response<1e-6){
	  if(finGEN_response.eof()) break;
	  else
	    continue;
	}


	switch(pdg){
	case 13:
	  response=2.0/e;
	  RMS=0;
	  break;
	case 22: case 11:
	  response = 1;
	  RMS=sqrt(pow(0.15/sqrt(e),2) + 0.01*0.01)*response;
	  break;
	default:
	  if(e<16)
	    response = 0.7 + 0.077*log(e);
	  else if(e>=16 && e<2000)
	    response = 0.86 + 0.0202*log(e);
	  else if(e>2000)
	    response = 1.0;
	  RMS=sqrt(pow(0.43/sqrt(e),2) + 0.009*0.009)*response;	
	  break;
	}

	//	std::cout << pdg << "\t" << e << "\t" << response << "\t" << RMS << endl;

   	if(response<1e-6){
	  if(finGEN_response.eof()) break;
	  else
	    continue;
	}

	float gauss_smear = myrndm->Gaus(response,RMS);
//   	std::cout << "e = " << e << std::endl;
//   	std::cout << "RMS = " << RMS << std::endl;
// 	std::cout << "response = " << response << std::endl;
//   	std::cout << "response = " << gauss_smear << std::endl;
        fastjet::PseudoJet curPseudoJet( px*gauss_smear, 
					 py*gauss_smear, 
					 pz*gauss_smear, 
					 e*gauss_smear);

        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
        
        if(finGEN_resolution.eof()) break;
        ctr++;
    }
    
}


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventCalo( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = etamax;
    
    float npart, x, y, z, e, pdgid, isCh, isPU = 0;
    
    int ctr = 0;
    while(true){
        
        finCalo  >> pdgid >> x >> y >> z >> e;         
        // std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
    
        TVector3 caloposition(x,y,z);
        double curPhi = caloposition.Phi();
        double curEta = caloposition.Eta();
        double curPt = sin(caloposition.Theta())*e;
        TLorentzVector cp4;
        cp4.SetPtEtaPhiE(curPt,curEta,curPhi,e);

        if (pdgid == -99){
            return;
        }        
	if(cp4.E()<1.0)
	  {
	    if(finCalo.eof()) break;
	    else continue;
	  }
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( cp4.Px(), cp4.Py(), cp4.Pz(), cp4.E() );
        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
        
        if(finCalo.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventTrack( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = etamax;
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
    
    int ctr = 0;
    while(true){
        
        finTrack  >> pdgid >> px >> py >> pz >> e;         
        // std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
    
        if (pdgid == -99){
            return;
        }        
    
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
        
        if(finTrack.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal, std::vector < fastjet::PseudoJet > MCparticles, bool fillJSub){

  //    std::cout << "analyzing event..." << particles.size() << std::endl;
    double rParam = rVal;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    
    if(counter%8==0)
      for(unsigned int i=0; i < particles.size(); i++)
	heta_after1->Fill(particles[i].eta());
      

    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    
    fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(particles, jetDef, fjAreaDefinition);
    // taking minium E = 0 for pion gun
    std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(25.0));
 
    if(counter%8==0)
      {
	for (unsigned int i = 0; i < out_jets.size() ; i++)
	  {
	    heta_after2->Fill(out_jets[i].eta());
	    hy->Fill(out_jets[i].rapidity());
	  }

	ISRPhotons.clear();	
	for (unsigned int i = 0; i < out_jets.size() ; i++)
	  {
	    if(out_jets[i].constituents().size()!=1)continue;
	    if(out_jets[i].constituents()[0].user_index()!=22)continue;
	    ISRPhotons.push_back(out_jets[i]);	    	    
	  }
	
      }

    else if(counter%8==5)
      for (unsigned int i = 0; i < out_jets.size() ; i++)
	heta_PF->Fill(out_jets[i].eta());




    // groomers/taggers
    fastjet::Pruner pruner1( fastjet::cambridge_algorithm, 0.1, 0.5 );
    fastjet::Filter trimmer1( fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03)) );
    double beta_sd = 1.0;
    double zcut_sd = 0.1;
    double mu_sd   = 1.0;
    fastjet::contrib::SoftDrop soft_drop_mmdt(0.0, zcut_sd, mu_sd);
    fastjet::contrib::SoftDrop soft_drop_sdb2(2.0, zcut_sd, mu_sd);
    fastjet::contrib::SoftDrop soft_drop_sdm1(-1.0, zcut_sd, mu_sd);   

    // n-subjettiness   
    const double RPARAM=jetRadius;
    double beta = 1;      // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
    double R0   = RPARAM; // Characteristic jet radius for normalization              
    double Rcut = RPARAM; // maximum R particles can be from axis to be included in jet   

    // beta = 1                   
    fastjet::contrib::Nsubjettiness nSub1KT_b1(1, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(1));
    fastjet::contrib::Nsubjettiness nSub2KT_b1(2, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(1));
    fastjet::contrib::Nsubjettiness nSub3KT_b1(3, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(1));

    // beta = 2
    fastjet::contrib::Nsubjettiness nSub1KT_b2(1, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(2));
    fastjet::contrib::Nsubjettiness nSub2KT_b2(2, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(2));
    fastjet::contrib::Nsubjettiness nSub3KT_b2(3, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(2));

    // ECF
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C1_b0(1,0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C1_b1(1,1,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C1_b2(1,2,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C2_b1(2,1,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C2_b2(2,2,fastjet::contrib::EnergyCorrelator::pt_R);

    fastjet::contrib::EnergyCorrelatorD2 ECF_D2_b1 (1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorD2 ECF_D2_b2 (2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorGeneralizedD2 ECF_D2_a1_b1 (1.0,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorGeneralizedD2 ECF_D2_a1_b2 (1.0,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorMseries ECF_M2_b1 (2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorMseries ECF_M2_b2 (2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorNseries ECF_N2_b1 (2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorNseries ECF_N2_b2 (2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);


    // take only the two leading jets for Z' and H2 decays
    int numGoodJets=0;

    for (unsigned int i = 0; i < out_jets.size() && numGoodJets<2; i++){

      // check if this jet is matched to an ISR photon
      bool isISRPhoton=false;
      for(unsigned int k =0; k < ISRPhotons.size(); k++){	

	TLorentzVector temp_pho(0,0,0,0);
	temp_pho.SetPxPyPzE(ISRPhotons[k].px(),
			    ISRPhotons[k].py(),
			    ISRPhotons[k].pz(),
			    ISRPhotons[k].e());

	TLorentzVector temp_jet(0,0,0,0);
	temp_jet.SetPxPyPzE(out_jets[i].px(),
			    out_jets[i].py(),
			    out_jets[i].pz(),
			    out_jets[i].e());
			  
	double dr = temp_pho.DeltaR(temp_jet);

	double Eratio = ISRPhotons[k].E()>1e-6? out_jets[i].E()/ISRPhotons[k].E(): -1;
	if(dr < 0.1 && Eratio > 0.7 && Eratio < 1.3)
	  {
	    isISRPhoton=true;
	    break;
	  }
      
      }
      
      if(isISRPhoton)continue;
      numGoodJets++;

      bool isleptag = false;
      double minDR = 9999.;
      for (unsigned int j = 0; j < MCparticles.size(); j++){
	// std::cout << "MCparticles = " << MCparticles[j].pt() << "," << MCparticles[j].user_index() << std::endl;
	double pdgid =  fabs(MCparticles[j].user_index());
	double dr = sqrt( pow(MCparticles[j].phi()-out_jets[i].phi(),2) + pow(MCparticles[j].eta()-out_jets[i].eta(),2) );
	// std::cout << "dr = " << dr << "," << pdgid << std::endl;
	if (minDR > dr && (fabs(pdgid)==11 || fabs(pdgid)==13 || fabs(pdgid)==15) ){ minDR = dr; }
      }
      if (minDR < jetRadius){ isleptag = true; }
      je.push_back( out_jets[i].E() );
      jpt.push_back( out_jets[i].pt() );
      jp.push_back( sqrt(out_jets[i].modp2()) );
      jeta.push_back( out_jets[i].eta() );
      jphi.push_back( out_jets[i].phi() );
      jmass.push_back( out_jets[i].m() );  
      jmass_sd.push_back( soft_drop_mmdt( out_jets.at(i) ).m() );        

      fastjet::Pruner mypruner(fastjet::cambridge_algorithm,0.1,0.5);
      PseudoJet pruned_jet = mypruner(out_jets[i]); 
      jmass_pr.push_back(pruned_jet.m());

      jmultiplicity.push_back( out_jets.at(i).constituents().size() );
      jisleptag.push_back( isleptag );
      

      
      int genHIndex=-1;
      for(unsigned int ip=0; ip<2; ip++)
	{
	  double dr = sqrt(pow(h_gen[i].p4.Eta()-out_jets[i].eta(),2)+
			   pow(h_gen[i].p4.Phi()-out_jets[i].phi(),2));
	  if(dr < jetRadius)
	    {
	      genHIndex=ip;
	      break;
	    }
	}
  
      if(genHIndex<0)
	costheta1.push_back(-999);
      else
	costheta1.push_back(h_gen[i].cosTheta1);


      // a bunch of substructure variables
      if(!fillJSub)continue;
      // N-subjettiness
      fastjet::PseudoJet curjet = out_jets[i];
      j_tau1_b1.push_back(  nSub1KT_b1(curjet) );        
      j_tau2_b1.push_back(  nSub2KT_b1(curjet) );        
      j_tau3_b1.push_back(  nSub3KT_b1(curjet) );        
      j_tau1_b2.push_back(  nSub1KT_b2(curjet) );        
      j_tau2_b2.push_back(  nSub2KT_b2(curjet) );  
      j_tau3_b2.push_back(  nSub3KT_b2(curjet) );  
      j_tau21_b1.push_back( nSub2KT_b1(curjet) / nSub1KT_b1(curjet) );
      j_tau21_b2.push_back( nSub2KT_b2(curjet) / nSub1KT_b2(curjet) );
      j_tau32_b1.push_back( nSub3KT_b1(curjet) / nSub2KT_b1(curjet) );
      j_tau32_b2.push_back( nSub3KT_b2(curjet) / nSub2KT_b2(curjet) );

      // Z log Z
      float zlogz = 0.;
      for (int i = 0; i < curjet.constituents().size(); i++){
	float conste = curjet.constituents().at(i).e()/curjet.e();
	zlogz += conste * log(conste);
      }

      // energy correlator     
      j_zlogz.push_back( zlogz );   
      j_c1_b0.push_back( ECF_C1_b0(curjet) );
      j_c1_b1.push_back( ECF_C1_b1(curjet) );
      j_c1_b2.push_back( ECF_C1_b2(curjet) );
      j_c2_b1.push_back( ECF_C2_b1(curjet) );
      j_c2_b2.push_back( ECF_C2_b2(curjet) );
      j_d2_b1.push_back( ECF_D2_b1(curjet) );
      j_d2_b2.push_back( ECF_D2_b2(curjet) );
      j_d2_a1_b1.push_back( ECF_D2_a1_b1(curjet) );
      j_d2_a1_b2.push_back( ECF_D2_a1_b2(curjet) );
      j_m2_b1.push_back( ECF_M2_b1(curjet) );
      j_m2_b2.push_back( ECF_M2_b2(curjet) );
      j_n2_b1.push_back( ECF_N2_b1(curjet) );
      j_n2_b2.push_back( ECF_N2_b2(curjet) );

      // Groomed variables
      fastjet::PseudoJet mmdtjet=soft_drop_mmdt(curjet);

      j_tau1_b1_mmdt.push_back(  nSub1KT_b1(mmdtjet) );        
      j_tau2_b1_mmdt.push_back(  nSub2KT_b1(mmdtjet) );        
      j_tau3_b1_mmdt.push_back(  nSub3KT_b1(mmdtjet) );        
      j_tau1_b2_mmdt.push_back(  nSub1KT_b2(mmdtjet) );        
      j_tau2_b2_mmdt.push_back(  nSub2KT_b2(mmdtjet) );  
      j_tau3_b2_mmdt.push_back(  nSub3KT_b2(mmdtjet) );  
      j_tau21_b1_mmdt.push_back( nSub2KT_b1(mmdtjet) / nSub1KT_b1(mmdtjet) );
      j_tau21_b2_mmdt.push_back( nSub2KT_b2(mmdtjet) / nSub1KT_b2(mmdtjet) );
      j_tau32_b1_mmdt.push_back( nSub3KT_b1(mmdtjet) / nSub2KT_b1(mmdtjet) );
      j_tau32_b2_mmdt.push_back( nSub3KT_b2(mmdtjet) / nSub2KT_b2(mmdtjet) );
      j_c1_b0_mmdt.push_back( ECF_C1_b0(mmdtjet) );
      j_c1_b1_mmdt.push_back( ECF_C1_b1(mmdtjet) );
      j_c1_b2_mmdt.push_back( ECF_C1_b2(mmdtjet) );
      j_c2_b1_mmdt.push_back( ECF_C2_b1(mmdtjet) );
      j_c2_b2_mmdt.push_back( ECF_C2_b2(mmdtjet) );
      j_d2_b1_mmdt.push_back( ECF_D2_b1(mmdtjet) );
      j_d2_b2_mmdt.push_back( ECF_D2_b2(mmdtjet) );
      j_d2_a1_b1_mmdt.push_back( ECF_D2_a1_b1(mmdtjet) );
      j_d2_a1_b2_mmdt.push_back( ECF_D2_a1_b2(mmdtjet) );
      j_m2_b1_mmdt.push_back( ECF_M2_b1(mmdtjet) );
      j_m2_b2_mmdt.push_back( ECF_M2_b2(mmdtjet) );
      j_n2_b1_mmdt.push_back( ECF_N2_b1(mmdtjet) );
      j_n2_b2_mmdt.push_back( ECF_N2_b2(mmdtjet) );
        
      j_mass_trim.push_back( trimmer1( curjet ).m() );
      j_mass_prun.push_back( pruner1( curjet ).m() );    
      j_mass_mmdt.push_back( mmdtjet.m() );
      j_mass_sdb2.push_back( soft_drop_sdb2( curjet ).m() );
      j_mass_sdm1.push_back( soft_drop_sdm1( curjet ).m() );


    }
    njets = numGoodJets;

    counter++;    

}
////////////////////////////////////////////////////////////////////////////////////////
void analyzeMCEvent(std::vector < fastjet::PseudoJet > MCparticles){

  fastjet::PseudoJet wplus;
  fastjet::PseudoJet wminus;
  fastjet::PseudoJet Zprime;
  for (unsigned int i = 0; i < MCparticles.size(); i++){
    double pdgid =  MCparticles[i].user_index();
    // for Z'->WW
    if (pdgid == 24){         
      wplus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
    }
    if (pdgid == -24){
      wminus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
    }
    if (pdgid == 32){
      Zprime.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
    }
    // for H2->hh
    if (pdgid == 25 && wplus.e()<1e-6){         
      wplus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
    }
    else if(pdgid == 25 && wminus.e()<1e-6){         
      wminus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
    }
    // for Z'->qq
    if (abs(pdgid)<=5 && pdgid>0 && wplus.e()<1e-6){         
      wplus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
    }
    else if(abs(pdgid)<=5 && pdgid<0 && wminus.e()<1e-6){         
      wminus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
    }
    
    if (pdgid == 35){
      Zprime.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
    }
  }
  // fastjet::PseudoJet zprime = wplus + wminus;
  //  std::cout << wplus.px() << ", " << wminus.px() << std::endl;
  gen_mZp = Zprime.m();
  gen_mWW = (wplus+wminus).m();


  be.push_back( wplus.E() );
  bpt.push_back( wplus.pt() );
  bp.push_back( sqrt(wplus.modp2()) );
  beta.push_back( wplus.eta() );
  bphi.push_back( wplus.phi() );
  bmass.push_back( wplus.m() );  

  be.push_back( wminus.E() );
  bpt.push_back( wminus.pt() );
  bp.push_back( sqrt(wminus.modp2()) );
  beta.push_back( wminus.eta() );
  bphi.push_back( wminus.phi() );
  bmass.push_back( wminus.m() );  


}

