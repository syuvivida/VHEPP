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
 
 c++ -o anaRawHit `root-config --glibs --cflags` `/Users/ntran/Research/CMS/PhysicsAnalysis/boost2013/fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -lvectorDict -lEnergyCorrelator anaRawHit.cpp
 
 TO RUN:     
 
 ./analysis01 ../madevent_3/Events/pietro_test14_unweighted_events.lhe ../Z2j/run_01_unweighted_events.lhe 
 
 */

//! ========================================================================================

////////////////////-----------------------------------------------
// Global variables
ifstream fin;
ifstream finMC;
ifstream finGEN;
ifstream finGEN_nonu;
ifstream finGEN_response;
ifstream finCalo;
ifstream finTrack;

ifstream finECalHits;
ifstream finHCalHits;

ifstream finECalHits2;
ifstream finHCalHits2;

TRandom3* myrndm;
std::vector<fastjet::PseudoJet> ISRPhotons;
TH1F* heta_status3;
TH1F* heta;
TH1F* heta_after1;
TH1F* heta_after2;
TH1F* heta_PF;
TH1F* hy;

TH1F* hecalhit_energy;
TH1F* hhcalhit_energy;
TH1F* hecalhit_energy_log;
TH1F* hhcalhit_energy_log;

int counter;
int evtCtr;

const int NJOBS=8;
int njets;
double gen_mZp;
double gen_mWW;

//jet-level ntuple
std::vector<float> je;
std::vector<float> jeecal;
std::vector<float> jehcal;
std::vector<float> jepho;
std::vector<float> jpt;
std::vector<float> jp;
std::vector<float> jeta;
std::vector<float> jphi;
std::vector<float> jmass;
std::vector<int> jmultiplicity;
std::vector<bool> jisleptag;
std::vector<float> jmass_sd;

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

//hit-level ntuple
int nhits;
std::vector<float> genMEta; // mother eta
std::vector<float> hitdphi;   
std::vector<float> hitdeta;   
std::vector<float> hitdphi12;
std::vector<float> hitdeta12;
std::vector<float> hitsigndR;
std::vector<int> hitindex;  
std::vector<int> hitdec;    
std::vector<float> hitpx;     
std::vector<float> hitpy;     
std::vector<float> hitpz;     
std::vector<float> hite;      



float jetRadius=0.4;
int mode=0; // 0: cut on emin directly E> emin GeV, 1: remove the lowest emin%
float emin=0.5;
const float etamax=99999.0;
//const float etamax=1.1;

////////////////////-----------------------------------------------
void readEvent( std::vector< fastjet::PseudoJet > &allParticles );
void readEventMC( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN_nonu( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN_response( std::vector< fastjet::PseudoJet > &allParticles );
void readEventCalo( std::vector< fastjet::PseudoJet > &allParticles );
void readEventRawHits( std::vector< fastjet::PseudoJet > &allParticles );
void readEventRawHits2( std::vector< fastjet::PseudoJet > &allParticles );
void readEventTrack(  std::vector< fastjet::PseudoJet > &allParticles );

void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal, std::vector < fastjet::PseudoJet > MCparticles, bool fillJSub=false);
void analyzeMCEvent(std::vector < fastjet::PseudoJet > MCparticles);
void clearVectors(){
    //  INIT
  njets = 0;
  je.clear();
  jpt.clear();
  jp.clear();    
  jeta.clear();
  jphi.clear();        
  jmass.clear();
  jmass_sd.clear();        
  jmultiplicity.clear();
  jisleptag.clear();
  
  jeecal.clear();
  jehcal.clear();
  jepho.clear();

  nhits = 0;
  genMEta.clear();
  hitdphi.clear();
  hitdeta.clear();
  hitdphi12.clear();
  hitdeta12.clear();
  hitsigndR.clear();
  hitindex.clear();
  hitdec.clear();
  hitpx.clear();
  hitpy.clear();
  hitpz.clear();
  hite.clear();

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
////////////////////////////////////////////////////////////////////////////////////////
int main (int argc, char **argv) {
    
    // std::cout << "hello world" << std::endl;
    std::string type = "of_PanPFA";   // type "gg" or "qq"
    std::string inputFolder = argv[1];
    jetRadius = (float)atof(argv[2]);
    std::cout << "jetRadius = " << jetRadius << std::endl;
    
    mode = (int)atoi(argv[3]);
    emin = (float)atof(argv[4]);

    if(mode==0)
      std::cout << "cutting on cal hit energy E > " << emin << " GeV" << endl;
    else if(mode==1)
      std::cout << "removing hits with the lowest " << emin << "% of energy" << endl;
    else
      {
	std::cout << "Invalid option" << endl;
	return 0;
      }
      

    myrndm = new TRandom3();


    // pick only the highest pt particle
    std::vector < fastjet::PseudoJet > allParticles;
    std::vector < fastjet::PseudoJet > allParticlesCalo; // calo clusters matched to PF candidate
    std::vector < fastjet::PseudoJet > allParticlesTrack; // tracks matched to PF candidate
    std::vector < fastjet::PseudoJet > allParticlesMC; //status3
    std::vector < fastjet::PseudoJet > allParticlesGEN; //status1 
    std::vector < fastjet::PseudoJet > allParticlesGEN_nonu; //status1
    std::vector < fastjet::PseudoJet > allParticlesGEN_response; //gen scaled with single particle response

    std::vector < fastjet::PseudoJet > allParticlesRawHits; // calo hits
    std::vector < fastjet::PseudoJet > allParticlesRawHits2; // calo hits


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
    finGEN_nonu.open(fnameGEN);
    finGEN_response.open(fnameGEN);

    char fnameCalo[150];
    sprintf( fnameCalo, "%s/of_CaloClusters.dat", inputFolder.c_str());
    finCalo.open(fnameCalo);

    sprintf( fnameCalo, "%s/of_ECaloHits_r.dat", inputFolder.c_str());
    finECalHits.open(fnameCalo);
    finECalHits2.open(fnameCalo);

    sprintf( fnameCalo, "%s/of_HCaloHits_r.dat", inputFolder.c_str());
    finHCalHits.open(fnameCalo);
    finHCalHits2.open(fnameCalo);

    char fnameTrack[150];
    sprintf( fnameTrack, "%s/of_Tracks.dat", inputFolder.c_str());
    finTrack.open(fnameTrack);

    char outName[192];
    if(NJOBS==7)
      sprintf( outName, "%s/radius%.1f_rawhit.root", inputFolder.c_str(), jetRadius);
    else
      sprintf( outName, "%s/radius%.1f_rawhit_fastjet.root", inputFolder.c_str(), jetRadius);
    TFile *f = TFile::Open(outName,"RECREATE");

    heta = new TH1F("heta","",200,-2.5,2.5);
    hy   = (TH1F*)heta->Clone("hy");
    heta_status3 = (TH1F*)heta->Clone("heta_status3");
    heta_after1 = (TH1F*)heta->Clone("heta_after1");
    heta_after2 = (TH1F*)heta->Clone("heta_after2");
    heta_PF     = (TH1F*)heta->Clone("heta_PF");

//     hecalhit_energy = new TH1F("hecalhit_energy","",2000,0,200);
//     hhcalhit_energy = new TH1F("hhcalhit_energy","",2000,0,200);
//     hecalhit_energy_log = new TH1F("hecalhit_energy_log","",700,-5,2);
//     hhcalhit_energy_log = new TH1F("hhcalhit_energy_log","",700,-5,2);

    hecalhit_energy = new TH1F("hecalhit_energy","",1000,0,2000);
    hhcalhit_energy = new TH1F("hhcalhit_energy","",1000,0,2000);
    hecalhit_energy_log = new TH1F("hecalhit_energy_log","",900,-5,4);
    hhcalhit_energy_log = new TH1F("hhcalhit_energy_log","",900,-5,4);

    TTree *tPFA = new TTree("tPFA","Tree with vectors");
    tPFA->Branch("njets"          , &njets      );
    tPFA->Branch("je"             , &je        );
    tPFA->Branch("jpt"            , &jpt        );
    tPFA->Branch("jp"             , &jp        );    
    tPFA->Branch("jeta"           , &jeta       );
    tPFA->Branch("jphi"           , &jphi       );
    tPFA->Branch("jmass"          , &jmass      );    
    tPFA->Branch("jmass_sd"       , &jmass_sd      );    
    tPFA->Branch("jmultiplicity"  , &jmultiplicity      );    
    tPFA->Branch("jisleptag"      , &jisleptag      );    

    TTree *tcalo = new TTree("tcalo","Tree with vectors");
    tcalo->Branch("njets"          , &njets      );
    tcalo->Branch("je"             , &je         );
    tcalo->Branch("jpt"            , &jpt        );
    tcalo->Branch("jp"             , &jp        );      
    tcalo->Branch("jeta"           , &jeta       );
    tcalo->Branch("jphi"           , &jphi       );
    tcalo->Branch("jmass"          , &jmass      );    
    tcalo->Branch("jmass_sd"       , &jmass_sd      );    
    tcalo->Branch("jmultiplicity"  , &jmultiplicity      );    
    tcalo->Branch("jisleptag"      , &jisleptag      );    

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

    TTree *trawhits = new TTree("trawhits","Tree with vectors");
    trawhits->Branch("njets"          , &njets      );
    trawhits->Branch("je"             , &je         );
    trawhits->Branch("jpt"            , &jpt        );
    trawhits->Branch("jp"             , &jp        );      
    trawhits->Branch("jeta"           , &jeta       );
    trawhits->Branch("jphi"           , &jphi       );
    trawhits->Branch("jmass"          , &jmass      );    
    trawhits->Branch("jmass_sd"       , &jmass_sd      );    
    trawhits->Branch("jmultiplicity"  , &jmultiplicity      );    
    trawhits->Branch("jisleptag"      , &jisleptag      );    

    trawhits->Branch("j_tau1_b1"        , &j_tau1_b1        );
    trawhits->Branch("j_tau2_b1"        , &j_tau2_b1        );
    trawhits->Branch("j_tau3_b1"        , &j_tau3_b1        );    
    trawhits->Branch("j_tau1_b2"        , &j_tau1_b2        );
    trawhits->Branch("j_tau2_b2"        , &j_tau2_b2        );
    trawhits->Branch("j_tau3_b2"        , &j_tau3_b2        );    
    trawhits->Branch("j_tau21_b1"       , &j_tau21_b1       );
    trawhits->Branch("j_tau21_b2"       , &j_tau21_b2       );
    trawhits->Branch("j_tau21_b1"       , &j_tau21_b1       );
    trawhits->Branch("j_tau21_b2"       , &j_tau21_b2       );
    trawhits->Branch("j_tau32_b1"       , &j_tau32_b1       );
    trawhits->Branch("j_tau32_b2"       , &j_tau32_b2       );
    trawhits->Branch("j_zlogz"          , &j_zlogz          );
    trawhits->Branch("j_c1_b0"          , &j_c1_b0          );
    trawhits->Branch("j_c1_b1"          , &j_c1_b1          );
    trawhits->Branch("j_c1_b2"          , &j_c1_b2          );
    trawhits->Branch("j_c2_b1"          , &j_c2_b1          );
    trawhits->Branch("j_c2_b2"          , &j_c2_b2          );
    trawhits->Branch("j_d2_b1"          , &j_d2_b1          );
    trawhits->Branch("j_d2_b2"          , &j_d2_b2          );
    trawhits->Branch("j_d2_a1_b1"       , &j_d2_a1_b1       );
    trawhits->Branch("j_d2_a1_b2"       , &j_d2_a1_b2       );
    trawhits->Branch("j_m2_b1"          , &j_m2_b1          );
    trawhits->Branch("j_m2_b2"          , &j_m2_b2          );
    trawhits->Branch("j_n2_b1"          , &j_n2_b1          );
    trawhits->Branch("j_n2_b2"          , &j_n2_b2          );

    trawhits->Branch("j_tau1_b1_mmdt"        , &j_tau1_b1_mmdt        );
    trawhits->Branch("j_tau2_b1_mmdt"        , &j_tau2_b1_mmdt        );
    trawhits->Branch("j_tau3_b1_mmdt"        , &j_tau3_b1_mmdt        );    
    trawhits->Branch("j_tau1_b2_mmdt"        , &j_tau1_b2_mmdt        );
    trawhits->Branch("j_tau2_b2_mmdt"        , &j_tau2_b2_mmdt        );
    trawhits->Branch("j_tau3_b2_mmdt"        , &j_tau3_b2_mmdt        );    
    trawhits->Branch("j_tau21_b1_mmdt"       , &j_tau21_b1_mmdt       );
    trawhits->Branch("j_tau21_b2_mmdt"       , &j_tau21_b2_mmdt       );
    trawhits->Branch("j_tau21_b1_mmdt"       , &j_tau21_b1_mmdt       );
    trawhits->Branch("j_tau21_b2_mmdt"       , &j_tau21_b2_mmdt       );
    trawhits->Branch("j_tau32_b1_mmdt"       , &j_tau32_b1_mmdt       );
    trawhits->Branch("j_tau32_b2_mmdt"       , &j_tau32_b2_mmdt       );
    trawhits->Branch("j_c1_b0_mmdt"          , &j_c1_b0_mmdt          );
    trawhits->Branch("j_c1_b1_mmdt"          , &j_c1_b1_mmdt          );
    trawhits->Branch("j_c1_b2_mmdt"          , &j_c1_b2_mmdt          );
    trawhits->Branch("j_c2_b1_mmdt"          , &j_c2_b1_mmdt          );
    trawhits->Branch("j_c2_b2_mmdt"          , &j_c2_b2_mmdt          );
    trawhits->Branch("j_d2_b1_mmdt"          , &j_d2_b1_mmdt          );
    trawhits->Branch("j_d2_b2_mmdt"          , &j_d2_b2_mmdt          );
    trawhits->Branch("j_d2_a1_b1_mmdt"       , &j_d2_a1_b1_mmdt       );
    trawhits->Branch("j_d2_a1_b2_mmdt"       , &j_d2_a1_b2_mmdt       );
    trawhits->Branch("j_m2_b1_mmdt"          , &j_m2_b1_mmdt          );
    trawhits->Branch("j_m2_b2_mmdt"          , &j_m2_b2_mmdt          );
    trawhits->Branch("j_n2_b1_mmdt"          , &j_n2_b1_mmdt          );
    trawhits->Branch("j_n2_b2_mmdt"          , &j_n2_b2_mmdt          );

    trawhits->Branch("j_mass_trim"      , &j_mass_trim      );
    trawhits->Branch("j_mass_mmdt"      , &j_mass_mmdt      );
    trawhits->Branch("j_mass_prun"      , &j_mass_prun      );
    trawhits->Branch("j_mass_sdb2"      , &j_mass_sdb2      );
    trawhits->Branch("j_mass_sdm1"      , &j_mass_sdm1      );


    TTree *trawhits2 = new TTree("trawhits2","Tree with vectors");
    trawhits2->Branch("njets"          , &njets      );
    trawhits2->Branch("je"             , &je         );
    trawhits2->Branch("jeecal"         , &jeecal     );
    trawhits2->Branch("jehcal"         , &jehcal     );
    trawhits2->Branch("jpt"            , &jpt        );
    trawhits2->Branch("jp"             , &jp        );      
    trawhits2->Branch("jeta"           , &jeta       );
    trawhits2->Branch("jphi"           , &jphi       );
    trawhits2->Branch("jmass"          , &jmass      );    
    trawhits2->Branch("jmass_sd"       , &jmass_sd      );    
    trawhits2->Branch("jmultiplicity"  , &jmultiplicity      );    
    trawhits2->Branch("jisleptag"      , &jisleptag      );    

    trawhits2->Branch("j_tau1_b1"        , &j_tau1_b1        );
    trawhits2->Branch("j_tau2_b1"        , &j_tau2_b1        );
    trawhits2->Branch("j_tau3_b1"        , &j_tau3_b1        );    
    trawhits2->Branch("j_tau1_b2"        , &j_tau1_b2        );
    trawhits2->Branch("j_tau2_b2"        , &j_tau2_b2        );
    trawhits2->Branch("j_tau3_b2"        , &j_tau3_b2        );    
    trawhits2->Branch("j_tau21_b1"       , &j_tau21_b1       );
    trawhits2->Branch("j_tau21_b2"       , &j_tau21_b2       );
    trawhits2->Branch("j_tau21_b1"       , &j_tau21_b1       );
    trawhits2->Branch("j_tau21_b2"       , &j_tau21_b2       );
    trawhits2->Branch("j_tau32_b1"       , &j_tau32_b1       );
    trawhits2->Branch("j_tau32_b2"       , &j_tau32_b2       );
    trawhits2->Branch("j_zlogz"          , &j_zlogz          );
    trawhits2->Branch("j_c1_b0"          , &j_c1_b0          );
    trawhits2->Branch("j_c1_b1"          , &j_c1_b1          );
    trawhits2->Branch("j_c1_b2"          , &j_c1_b2          );
    trawhits2->Branch("j_c2_b1"          , &j_c2_b1          );
    trawhits2->Branch("j_c2_b2"          , &j_c2_b2          );
    trawhits2->Branch("j_d2_b1"          , &j_d2_b1          );
    trawhits2->Branch("j_d2_b2"          , &j_d2_b2          );
    trawhits2->Branch("j_d2_a1_b1"       , &j_d2_a1_b1       );
    trawhits2->Branch("j_d2_a1_b2"       , &j_d2_a1_b2       );
    trawhits2->Branch("j_m2_b1"          , &j_m2_b1          );
    trawhits2->Branch("j_m2_b2"          , &j_m2_b2          );
    trawhits2->Branch("j_n2_b1"          , &j_n2_b1          );
    trawhits2->Branch("j_n2_b2"          , &j_n2_b2          );

    trawhits2->Branch("j_tau1_b1_mmdt"        , &j_tau1_b1_mmdt        );
    trawhits2->Branch("j_tau2_b1_mmdt"        , &j_tau2_b1_mmdt        );
    trawhits2->Branch("j_tau3_b1_mmdt"        , &j_tau3_b1_mmdt        );    
    trawhits2->Branch("j_tau1_b2_mmdt"        , &j_tau1_b2_mmdt        );
    trawhits2->Branch("j_tau2_b2_mmdt"        , &j_tau2_b2_mmdt        );
    trawhits2->Branch("j_tau3_b2_mmdt"        , &j_tau3_b2_mmdt        );    
    trawhits2->Branch("j_tau21_b1_mmdt"       , &j_tau21_b1_mmdt       );
    trawhits2->Branch("j_tau21_b2_mmdt"       , &j_tau21_b2_mmdt       );
    trawhits2->Branch("j_tau21_b1_mmdt"       , &j_tau21_b1_mmdt       );
    trawhits2->Branch("j_tau21_b2_mmdt"       , &j_tau21_b2_mmdt       );
    trawhits2->Branch("j_tau32_b1_mmdt"       , &j_tau32_b1_mmdt       );
    trawhits2->Branch("j_tau32_b2_mmdt"       , &j_tau32_b2_mmdt       );
    trawhits2->Branch("j_c1_b0_mmdt"          , &j_c1_b0_mmdt          );
    trawhits2->Branch("j_c1_b1_mmdt"          , &j_c1_b1_mmdt          );
    trawhits2->Branch("j_c1_b2_mmdt"          , &j_c1_b2_mmdt          );
    trawhits2->Branch("j_c2_b1_mmdt"          , &j_c2_b1_mmdt          );
    trawhits2->Branch("j_c2_b2_mmdt"          , &j_c2_b2_mmdt          );
    trawhits2->Branch("j_d2_b1_mmdt"          , &j_d2_b1_mmdt          );
    trawhits2->Branch("j_d2_b2_mmdt"          , &j_d2_b2_mmdt          );
    trawhits2->Branch("j_d2_a1_b1_mmdt"       , &j_d2_a1_b1_mmdt       );
    trawhits2->Branch("j_d2_a1_b2_mmdt"       , &j_d2_a1_b2_mmdt       );
    trawhits2->Branch("j_m2_b1_mmdt"          , &j_m2_b1_mmdt          );
    trawhits2->Branch("j_m2_b2_mmdt"          , &j_m2_b2_mmdt          );
    trawhits2->Branch("j_n2_b1_mmdt"          , &j_n2_b1_mmdt          );
    trawhits2->Branch("j_n2_b2_mmdt"          , &j_n2_b2_mmdt          );

    trawhits2->Branch("j_mass_trim"      , &j_mass_trim      );
    trawhits2->Branch("j_mass_mmdt"      , &j_mass_mmdt      );
    trawhits2->Branch("j_mass_prun"      , &j_mass_prun      );
    trawhits2->Branch("j_mass_sdb2"      , &j_mass_sdb2      );
    trawhits2->Branch("j_mass_sdm1"      , &j_mass_sdm1      );


    TTree* hittuple = new TTree("hittuple","Tree with calo hits");
    hittuple->Branch("nhits", &nhits);
    hittuple->Branch("genMEta",&genMEta);
    hittuple->Branch("hitdphi", &hitdphi);
    hittuple->Branch("hitdeta", &hitdeta);
    hittuple->Branch("hitdphi12", &hitdphi12);
    hittuple->Branch("hitdeta12", &hitdeta12);
    hittuple->Branch("hitsigndR", &hitsigndR);
    hittuple->Branch("hitindex", &hitindex);
    hittuple->Branch("hitdec", &hitdec);
    hittuple->Branch("hitpx", &hitpx);
    hittuple->Branch("hitpy", &hitpy);
    hittuple->Branch("hitpz", &hitpz);
    hittuple->Branch("hite", &hite);

    TTree* fjtuple = new TTree("fjtuple","Tree with calo hits from fastjet");
    fjtuple->Branch("nhits", &nhits);
    fjtuple->Branch("genMEta",&genMEta);
    fjtuple->Branch("hitdphi", &hitdphi);
    fjtuple->Branch("hitdeta", &hitdeta);
    fjtuple->Branch("hitindex", &hitindex);
    fjtuple->Branch("hitdec", &hitdec);
    fjtuple->Branch("hitpx", &hitpx);
    fjtuple->Branch("hitpy", &hitpy);
    fjtuple->Branch("hitpz", &hitpz);
    fjtuple->Branch("hite", &hite);

    TTree *ttrack = new TTree("ttrack","Tree with vectors");
    ttrack->Branch("njets"          , &njets      );
    ttrack->Branch("je"             , &je         );
    ttrack->Branch("jpt"            , &jpt        );
    ttrack->Branch("jp"             , &jp        );      
    ttrack->Branch("jeta"           , &jeta       );
    ttrack->Branch("jphi"           , &jphi       );
    ttrack->Branch("jmass"          , &jmass      );    
    ttrack->Branch("jmass_sd"       , &jmass_sd      );    
    ttrack->Branch("jmultiplicity"  , &jmultiplicity      );    
    ttrack->Branch("jisleptag"      , &jisleptag      );    

    TTree *tGEN = new TTree("tGEN","Tree with vectors");
    tGEN->Branch("njets"          , &njets      );
    tGEN->Branch("je"             , &je         );
    tGEN->Branch("jepho"          , &jepho      );
    tGEN->Branch("jpt"            , &jpt        );
    tGEN->Branch("jp"             , &jp        );          
    tGEN->Branch("jeta"           , &jeta       );
    tGEN->Branch("jphi"           , &jphi       );
    tGEN->Branch("jmass"          , &jmass      );    
    tGEN->Branch("jmass_sd"       , &jmass_sd      );    
    tGEN->Branch("jmultiplicity"  , &jmultiplicity      );    
    tGEN->Branch("jisleptag"      , &jisleptag      );        

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


    TTree *tGEN_nonu = new TTree("tGEN_nonu","Tree with vectors");
    tGEN_nonu->Branch("njets"          , &njets      );
    tGEN_nonu->Branch("je"             , &je         );
    tGEN_nonu->Branch("jepho"          , &jepho      );
    tGEN_nonu->Branch("jpt"            , &jpt        );
    tGEN_nonu->Branch("jp"             , &jp        );          
    tGEN_nonu->Branch("jeta"           , &jeta       );
    tGEN_nonu->Branch("jphi"           , &jphi       );
    tGEN_nonu->Branch("jmass"          , &jmass      );    
    tGEN_nonu->Branch("jmass_sd"       , &jmass_sd      );    
    tGEN_nonu->Branch("jmultiplicity"  , &jmultiplicity      );    
    tGEN_nonu->Branch("jisleptag"      , &jisleptag      );        

    TTree *tGEN_response = new TTree("tGEN_response","Tree with vectors");
    tGEN_response->Branch("njets"          , &njets      );
    tGEN_response->Branch("je"             , &je         );
    tGEN_response->Branch("jepho"          , &jepho      );
    tGEN_response->Branch("jpt"            , &jpt        );
    tGEN_response->Branch("jp"             , &jp        );          
    tGEN_response->Branch("jeta"           , &jeta       );
    tGEN_response->Branch("jphi"           , &jphi       );
    tGEN_response->Branch("jmass"          , &jmass      );    
    tGEN_response->Branch("jmass_sd"       , &jmass_sd      );    
    tGEN_response->Branch("jmultiplicity"  , &jmultiplicity      );    
    tGEN_response->Branch("jisleptag"      , &jisleptag      );        


    TTree *tMC = new TTree("tMC","Tree with vectors");
    tMC->Branch("gen_mZp"          , &gen_mZp      );
    tMC->Branch("gen_mWW"          , &gen_mWW      );

    int ctr = 0;
    counter =0;
    while(true){

        readEvent( allParticles );
        readEventGEN( allParticlesGEN ); 
	readEventGEN_nonu( allParticlesGEN_nonu );
	readEventGEN_response( allParticlesGEN_response);

        readEventMC( allParticlesMC );
        readEventCalo( allParticlesCalo );
	readEventRawHits( allParticlesRawHits );
 	readEventRawHits2( allParticlesRawHits2 );
	readEventTrack( allParticlesTrack);

        // std::cout << "size of collection = " << allParticles.size() << "," << allParticlesMC.size() << ", " << allParticlesCalo.size() << std::endl;
	if (ctr > 0){
  	    clearVectors();
 	    analyzeEvent( allParticlesGEN, jetRadius, allParticlesMC, true );
	    tGEN->Fill();
	
 	    clearVectors();
 	    analyzeEvent( allParticlesGEN_nonu, jetRadius, allParticlesMC );
	    tGEN_nonu->Fill();
	    
 	    clearVectors();
 	    analyzeEvent( allParticlesGEN_response, jetRadius, allParticlesMC );
	    tGEN_response->Fill();

            clearVectors();
            analyzeEvent( allParticles, jetRadius, allParticlesMC );
            tPFA->Fill();

            clearVectors();
            analyzeEvent( allParticlesCalo, jetRadius, allParticlesMC, true );
            tcalo->Fill();    
	    fjtuple->Fill();

            clearVectors();
	    analyzeEvent( allParticlesTrack, jetRadius, allParticlesMC );
	    ttrack->Fill();
	    
	    if(NJOBS==8){
	      clearVectors();
	      analyzeEvent( allParticlesRawHits, jetRadius, allParticlesMC, true );
	      trawhits->Fill();
	    }

	    clearVectors();
 	    analyzeEvent( allParticlesRawHits2, jetRadius, allParticlesMC,true );
 	    trawhits2->Fill();
	    hittuple->Fill();
  
            analyzeMCEvent( allParticlesMC );
            tMC->Fill();
        }

        allParticles.clear();
        allParticlesMC.clear();
        allParticlesGEN.clear();
        allParticlesGEN_nonu.clear();
        allParticlesGEN_response.clear();
        allParticlesCalo.clear();
        allParticlesTrack.clear();
	allParticlesRawHits.clear();
 	allParticlesRawHits2.clear();

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

    hecalhit_energy->Write();
    hhcalhit_energy->Write();
    hecalhit_energy_log->Write();
    hhcalhit_energy_log->Write();

    tPFA->Write();
    tcalo->Write();
    trawhits->Write();
    trawhits2->Write();
    hittuple->Write();
    fjtuple->Write();
    ttrack->Write();
    tGEN->Write();
    tGEN_nonu->Write();
    tGEN_response->Write();
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
    
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
// 	TLorentzVector temp_l4;
//  	temp_l4.SetPxPyPzE(px,py,pz,e);
// 	if(abs(pdgid)<=5 && abs(pdgid)>=1)
// 	  heta_status3->Fill(temp_l4.Eta());

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
	
// 	TLorentzVector temp_l4;
//  	temp_l4.SetPxPyPzE(px,py,pz,e);
// 	heta->Fill(temp_l4.Eta());

        if(finGEN.eof()) break;
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

	if(e<3)response=0;

// 	std::cout << "pdg = " << pdg << std::endl;
//     	std::cout << "e = " << e << std::endl;

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
void readEventRawHits( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = etamax;
    
    float npart, x, y, z, e, pdgid, isCh, isPU = 0;
    
    int ctr = 0;
    //    std::cout << "reading ECAL hits" << std::endl;
    while(true){
        
        finECalHits  >> pdgid >> x >> y >> z >> e;         
    
        TVector3 caloposition(x,y,z);
        double curPhi = caloposition.Phi();
        double curEta = caloposition.Eta();
        double curPt = sin(caloposition.Theta())*e;
        TLorentzVector cp4;
        cp4.SetPtEtaPhiE(curPt,curEta,curPhi,e);

        if (pdgid == -99){
            break;
        }        
 	if(cp4.E()<emin && mode==0)
 	  {
 	    if(finECalHits.eof()) break;
 	    else continue;
 	  }
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( cp4.Px(), cp4.Py(), cp4.Pz(), cp4.E() );
        curPseudoJet.set_user_index(pdgid);
	// 	std::cout << curPseudoJet.px() << ", " << curPseudoJet.py() << ", " << curPseudoJet.pz() << ", " << curPseudoJet.e() << endl;
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
        
        if(finECalHits.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }


    //    std::cout << "The size of ECAL allParticles = " << allParticles.size() << std::endl;
    
    ctr = 0;
    while(true){
      
       finHCalHits  >> pdgid >> x >> y >> z >> e;         
  
       TVector3 caloposition(x,y,z);
       double curPhi = caloposition.Phi();
       double curEta = caloposition.Eta();
       double curPt = sin(caloposition.Theta())*e;
       TLorentzVector cp4;
       cp4.SetPtEtaPhiE(curPt,curEta,curPhi,e);
       
       if (pdgid == -99){
	 break;
       }        
       if(cp4.E()<emin && mode==0)
	 {
	   if(finHCalHits.eof()) break;
	   else continue;
	 }
         // fill vector of pseudojets
       fastjet::PseudoJet curPseudoJet( cp4.Px(), cp4.Py(), cp4.Pz(), cp4.E() );
       curPseudoJet.set_user_index(pdgid);
       // 	std::cout << curPseudoJet.px() << ", " << curPseudoJet.py() << ", " << curPseudoJet.pz() << ", " << curPseudoJet.e() << endl;
       if (fabs(curPseudoJet.eta()) < etaMax){
	 allParticles.push_back( curPseudoJet );
       }
      
       if(finHCalHits.eof()) break;
       ctr++;
       // std::cout << "ctr2 = " << ctr << std::endl;
    }

    const unsigned int size_temp = allParticles.size();
//     std::cout << "Before erasing, the size of ECAL+HCAL allParticles = " << size_temp << std::endl;

    if(mode!=1)return;								      
    std::vector < fastjet::PseudoJet > allParticles_temp = allParticles;
    allParticles=sorted_by_pt(allParticles_temp);
    
    const unsigned int n_temp_to_be_removed = emin*size_temp;
    allParticles.erase(allParticles.end()-n_temp_to_be_removed,allParticles.end());
//     std::cout << "After erasing, the size of ECAL+HCAL allParticles = " << allParticles.size() << std::endl;
    
    return;
}



////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventRawHits2( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = etamax;
    
    float npart, x, y, z, e, pdgid, isCh, isPU = 0;
    
    int ctr = 0;
    while(true){
        
        finECalHits2  >> pdgid >> x >> y >> z >> e;         
	//	std::cout << "pdgid = " << pdgid << ", " << x << ", " << y << ", " << z << std::endl;
    
        TVector3 caloposition(x,y,z);
        double curPhi = caloposition.Phi();
        double curEta = caloposition.Eta();
        double curPt = sin(caloposition.Theta())*e;
        TLorentzVector cp4;
        cp4.SetPtEtaPhiE(curPt,curEta,curPhi,e);

        if (pdgid == -99){
            break;
        }        
// 	if(cp4.E()<1.0)
// 	  {
// 	    if(finECalHits2.eof()) break;
// 	    else continue;
// 	  }
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( cp4.Px(), cp4.Py(), cp4.Pz(), cp4.E() );
        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }

        if (fabs(curPseudoJet.eta()) < 1.5){
	  hecalhit_energy->Fill(cp4.E());
	  hecalhit_energy_log->Fill(TMath::Log10(cp4.E()));
        }
        
        if(finECalHits2.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }

    
    ctr = 0;
    while(true){
        
        finHCalHits2  >> pdgid >> x >> y >> z >> e;         
	//    	std::cout << "pdgid = " << pdgid << ", " << x << ", " << y << ", " << z << std::endl;

        TVector3 caloposition(x,y,z);
        double curPhi = caloposition.Phi();
        double curEta = caloposition.Eta();
        double curPt = sin(caloposition.Theta())*e;
        TLorentzVector cp4;
        cp4.SetPtEtaPhiE(curPt,curEta,curPhi,e);

        if (pdgid == -99){
            return;
        }        
// 	if(cp4.E()<1.0)
// 	  {
// 	    if(finHCalHits2.eof()) break;
// 	    else continue;
// 	  }
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( cp4.Px(), cp4.Py(), cp4.Pz(), cp4.E() );
        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }

        if (fabs(curPseudoJet.eta()) < 1.5){
	  hhcalhit_energy->Fill(cp4.E());
	  hhcalhit_energy_log->Fill(TMath::Log10(cp4.E()));
        }
        
        if(finHCalHits2.eof()) break;
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
    if(counter%NJOBS==0)
      for(unsigned int i=0; i < particles.size(); i++)       
	heta_after1->Fill(particles[i].eta());
      

    // doing my own clustering
    if(counter%NJOBS==NJOBS-1){
      fastjet::PseudoJet w[2];
      fastjet::PseudoJet quarks[2][2];

      for (unsigned int i = 0; i < MCparticles.size(); i++){
        double pdgid =  MCparticles[i].user_index();
            if (pdgid == 24){         
                w[0].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
            else if (pdgid == -24){
                w[1].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
	    // for Z'->qq in mu mu collider
            else if (pdgid>=1 && pdgid <=5 && w[0].e()<1e-6){         
                w[0].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
            else if (pdgid>=-5 && pdgid<=-1 && w[1].e()<1e-6){
                w[1].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
	    
	    
            else if (pdgid == 25 && w[0].e()<1e-6){         
                w[0].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
            else if (pdgid == 25 && w[1].e()<1e-6){
                w[1].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
	    
	    // if they are quarks
	    if(fabs(pdgid)>=1 && fabs(pdgid)<=5)
	      {
		if(quarks[0][0].e()<1e-6)
		  {
		    quarks[0][0].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
		    quarks[0][0].set_user_index(pdgid);
		  }
		else if(quarks[0][1].e()<1e-6)
		  {
		    quarks[0][1].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
		    quarks[0][1].set_user_index(pdgid);
		  }
		else if(quarks[1][0].e()<1e-6)
		  {
		    quarks[1][0].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
		    quarks[1][0].set_user_index(pdgid);
		  }
		else if(quarks[1][1].e()<1e-6)
		  {
		    quarks[1][1].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
		    quarks[1][1].set_user_index(pdgid);
		  }

	      }

      } // end loop of MC particles


      // reset W if it overlaps with ISR photons
	// check if this jet is matched to an ISR photon
      for(int iw=0; iw<2; iw++){
	
	if(w[iw].e()<1e-6)continue;
	//	cout << "W " << iw  << " " << w[iw].px() << " " << w[iw].py() << " " << w[iw].pz() << endl;	
	bool isISRPhoton=false;
	for(unsigned int k =0; k < ISRPhotons.size(); k++){

	  TLorentzVector temp_pho(0,0,0,0);
	  temp_pho.SetPxPyPzE(ISRPhotons[k].px(),
			      ISRPhotons[k].py(),
			      ISRPhotons[k].pz(),
			      ISRPhotons[k].e());

	  TLorentzVector temp_jet(0,0,0,0);
	  temp_jet.SetPxPyPzE(w[iw].px(),
			      w[iw].py(),
			      w[iw].pz(),
			      w[iw].e());
	    
	  double dr = temp_pho.DeltaR(temp_jet);

	  if(dr < jetRadius)
	    {
	      isISRPhoton=true;
	      break;
	    }
      
	}
      
	if(isISRPhoton){
	  std::cout << "W overlaps with an ISR Photon ! " << std::endl;
	  w[iw].reset(0,0,0,0);
	}

      } // end loop of generator-level W
    
      if(w[0].e()<1e-6 && w[1].e()<1e-6)
	{
	  std::cout << "W not found!" << std::endl;
	  counter++;
	  return;
	}

      fastjet::PseudoJet Wjet[2];
      Wjet[0].reset(0,0,0,0);
      Wjet[1].reset(0,0,0,0);
      int nPart[2]={0,0};
      double ecale[2] = {0.,0.};
      double hcale[2] = {0.,0.};
    
      for(unsigned int i=0; i < particles.size(); i++)
	{
	  bool used = false;
	  for(unsigned int iw=0; iw < 2; iw++)
	    {
	      if(w[iw].e()<1e-6)continue;
	      if(used)continue;

//  	      std::cout << w[iw].px() << " " << w[iw].py() << " " << w[iw].pz() << std::endl;
//  	      cout << quarks[iw][0].user_index() << " " << quarks[iw][0].px() << " " << quarks[iw][0].py() << " " << quarks[iw][0].pz() << endl;
//  	      cout << quarks[iw][1].user_index() << " " << quarks[iw][1].px() << " " << quarks[iw][1].py() << " " << quarks[iw][1].pz() << endl;

	      double dphi = TVector2::Phi_mpi_pi(particles[i].phi()-w[iw].phi());
	      double deta = particles[i].eta()-w[iw].eta();

	      double dphiq12 = TVector2::Phi_mpi_pi(quarks[iw][0].phi()-quarks[iw][1].phi());
	      double detaq12 = quarks[iw][0].eta()-quarks[iw][1].eta();

	      double signedDR = (detaq12*deta + dphiq12*dphi)/sqrt(detaq12*detaq12+dphiq12*dphiq12); 
//  	      cout << "deta = " << deta << " dphi = " << dphi << endl;
// 	      cout << "dphi = " << particles[i].phi()-w[iw].phi() << endl;
// 	      cout << "signedDR = " << signedDR << "\t" << (detaq12*deta + dphiq12*dphi) << endl;


	      TLorentzVector a_temp(0,0,0,0);
	      TLorentzVector b_temp(0,0,0,0);
	      a_temp.SetPxPyPzE(particles[i].px(),
				particles[i].py(),
				particles[i].pz(),
				particles[i].e());
	      b_temp.SetPxPyPzE(w[iw].px(),
				w[iw].py(),
				w[iw].pz(),
				w[iw].e());
				
	      double dr = a_temp.DeltaR(b_temp);
	      if(dr > jetRadius)continue;
	      used = true;
	      if((int)particles[i].user_index()==0){
		ecale[iw] += particles[i].e();
	      }
	      else if((int)particles[i].user_index()==1){
		hcale[iw] += particles[i].e();
	      }
	      else
		std::cout << "Warning!  user_index = " << particles[i].user_index() << std::endl;

	      nPart[iw]++;
	      Wjet[iw] += particles[i];
	      nhits++;
	      genMEta.push_back(w[iw].eta());
	      hitdphi.push_back(dphi);
	      hitdeta.push_back(deta);

	      hitdphi12.push_back(dphiq12);
	      hitdeta12.push_back(detaq12);

	      hitsigndR.push_back(signedDR);

	      hitindex.push_back(iw);
	      hitdec.push_back((int)particles[i].user_index());
	      hitpx.push_back(particles[i].px());
	      hitpy.push_back(particles[i].py());
	      hitpz.push_back(particles[i].pz());
	      hite.push_back(particles[i].e());

	    } // end loop of W/Higgs boson
	} // end loop of calo hits


      // filling clustered jet information
      for(unsigned int iw=0; iw < 2; iw++)
	{
	  if(Wjet[iw].e()<1e-6)continue;
	  njets++;
	  je.push_back( Wjet[iw].e() );
	  jeecal.push_back( ecale[iw] );
	  jehcal.push_back( hcale[iw] );
	  jpt.push_back( Wjet[iw].pt() );
	  jp.push_back( sqrt(Wjet[iw].modp2()) );
	  jeta.push_back( Wjet[iw].eta() );
	  jphi.push_back( Wjet[iw].phi() );
	  jmass.push_back( Wjet[iw].m() );  
	  jmass_sd.push_back( Wjet[iw].m() );        
	  jmultiplicity.push_back(nPart[iw]);
	  jisleptag.push_back(false);
	}

    } // if counter%NJOBS!=NJOBS-1
    else{

      int activeAreaRepeats = 1;
      double ghostArea = 0.01;
      double ghostEtaMax = 7.0;
    
      fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
      fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    
      fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(particles, jetDef, fjAreaDefinition);
      std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(25.0));
 
      if(counter%NJOBS==0)
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

      else if(counter%NJOBS==3)
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
//        fastjet::contrib::Nsubjettiness nSub1KT_b2(1, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(2));
//        fastjet::contrib::Nsubjettiness nSub2KT_b2(2, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(2));
//        fastjet::contrib::Nsubjettiness nSub3KT_b2(3, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(2));

       // ECF
//        fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C1_b0(1,0,fastjet::contrib::EnergyCorrelator::pt_R);
//        fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C1_b1(1,1,fastjet::contrib::EnergyCorrelator::pt_R);
//        fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C1_b2(1,2,fastjet::contrib::EnergyCorrelator::pt_R);
       fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C2_b1(2,1,fastjet::contrib::EnergyCorrelator::pt_R);
//        fastjet::contrib::EnergyCorrelatorDoubleRatio ECF_C2_b2(2,2,fastjet::contrib::EnergyCorrelator::pt_R);

//        fastjet::contrib::EnergyCorrelatorD2 ECF_D2_b1 (1.0,fastjet::contrib::EnergyCorrelator::pt_R);
//        fastjet::contrib::EnergyCorrelatorD2 ECF_D2_b2 (2.0,fastjet::contrib::EnergyCorrelator::pt_R);
//        fastjet::contrib::EnergyCorrelatorGeneralizedD2 ECF_D2_a1_b1 (1.0,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
//        fastjet::contrib::EnergyCorrelatorGeneralizedD2 ECF_D2_a1_b2 (1.0,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
//        fastjet::contrib::EnergyCorrelatorMseries ECF_M2_b1 (2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
//        fastjet::contrib::EnergyCorrelatorMseries ECF_M2_b2 (2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
//        fastjet::contrib::EnergyCorrelatorNseries ECF_N2_b1 (2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
//        fastjet::contrib::EnergyCorrelatorNseries ECF_N2_b2 (2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);



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

	// look for constituents of calojets

	if(counter%NJOBS==4)
	  {
	    for(unsigned int jc = 0; jc < out_jets[i].constituents().size(); jc++){	      
	      double dphi = TVector2::Phi_mpi_pi(out_jets[i].constituents()[jc].phi()-out_jets[i].phi());
	      double deta = out_jets[i].constituents()[jc].eta()-out_jets[i].eta();
	      nhits++;
	      genMEta.push_back(out_jets[i].eta());
	      hitdphi.push_back(dphi);
	      hitdeta.push_back(deta);
	      hitindex.push_back(numGoodJets);
	      hitpx.push_back(out_jets[i].constituents()[jc].px());
	      hitpy.push_back(out_jets[i].constituents()[jc].py());
	      hitpz.push_back(out_jets[i].constituents()[jc].pz());
	      hite.push_back(out_jets[i].constituents()[jc].e());
	    }
	  }
	
	numGoodJets++;

	
	if(counter%NJOBS<3){

	  double epho = 0;
	  for (unsigned int ip = 0; ip < out_jets[i].constituents().size() ; ip++)
	    {
	      unsigned int pdg = abs(out_jets[i].constituents()[ip].user_index());
	      if(pdg!=22)continue;
	      epho +=out_jets[i].constituents()[ip].e();      
	    }
	  jepho.push_back(epho);

	} // if this is filling generator-level information

	bool isleptag = false;
	double minDR = 9999.;
	for (unsigned int j = 0; j < MCparticles.size(); j++){
	  // std::cout << "MCparticles = " << MCparticles[j].pt() << "," << MCparticles[j].user_index() << std::endl;
	  double pdgid =  fabs(MCparticles[j].user_index());
	  double dr = sqrt( pow(TVector2::Phi_mpi_pi(MCparticles[j].phi()-out_jets[i].phi()),2) + pow(MCparticles[j].eta()-out_jets[i].eta(),2) );
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
	jmultiplicity.push_back( out_jets.at(i).constituents().size() );
	jisleptag.push_back( isleptag );

	// a bunch of substructure variables
	if(!fillJSub)continue;
// 	// N-subjettiness
 	fastjet::PseudoJet curjet = out_jets[i];
//  	j_tau1_b1.push_back(  nSub1KT_b1(curjet) );        
//  	j_tau2_b1.push_back(  nSub2KT_b1(curjet) );        
//  	j_tau3_b1.push_back(  nSub3KT_b1(curjet) );        
//  	j_tau1_b2.push_back(  nSub1KT_b2(curjet) );        
//  	j_tau2_b2.push_back(  nSub2KT_b2(curjet) );  
//  	j_tau3_b2.push_back(  nSub3KT_b2(curjet) );  
   	j_tau21_b1.push_back( nSub2KT_b1(curjet) / nSub1KT_b1(curjet) );
//   	j_tau21_b2.push_back( nSub2KT_b2(curjet) / nSub1KT_b2(curjet) );
   	j_tau32_b1.push_back( nSub3KT_b1(curjet) / nSub2KT_b1(curjet) );
//   	j_tau32_b2.push_back( nSub3KT_b2(curjet) / nSub2KT_b2(curjet) );

 	// Z log Z
//  	float zlogz = 0.;
//  	for (int i = 0; i < curjet.constituents().size(); i++){
//  	  float conste = curjet.constituents().at(i).e()/curjet.e();
//  	  zlogz += conste * log(conste);
//  	}

//  	// energy correlator     
//  	j_zlogz.push_back( zlogz );   
//   	j_c1_b0.push_back( ECF_C1_b0(curjet) );
//   	j_c1_b1.push_back( ECF_C1_b1(curjet) );
//   	j_c1_b2.push_back( ECF_C1_b2(curjet) );
	j_c2_b1.push_back( ECF_C2_b1(curjet) );
//   	j_c2_b2.push_back( ECF_C2_b2(curjet) );
//   	j_d2_b1.push_back( ECF_D2_b1(curjet) );
//   	j_d2_b2.push_back( ECF_D2_b2(curjet) );
//   	j_d2_a1_b1.push_back( ECF_D2_a1_b1(curjet) );
//   	j_d2_a1_b2.push_back( ECF_D2_a1_b2(curjet) );
//   	j_m2_b1.push_back( ECF_M2_b1(curjet) );
//   	j_m2_b2.push_back( ECF_M2_b2(curjet) );
//   	j_n2_b1.push_back( ECF_N2_b1(curjet) );
//   	j_n2_b2.push_back( ECF_N2_b2(curjet) );

//  	// Groomed variables
//  	fastjet::PseudoJet mmdtjet=soft_drop_mmdt(curjet);

//  	j_tau1_b1_mmdt.push_back(  nSub1KT_b1(mmdtjet) );        
//  	j_tau2_b1_mmdt.push_back(  nSub2KT_b1(mmdtjet) );        
//  	j_tau3_b1_mmdt.push_back(  nSub3KT_b1(mmdtjet) );        
//  	j_tau1_b2_mmdt.push_back(  nSub1KT_b2(mmdtjet) );        
//  	j_tau2_b2_mmdt.push_back(  nSub2KT_b2(mmdtjet) );  
//  	j_tau3_b2_mmdt.push_back(  nSub3KT_b2(mmdtjet) );  
//  	j_tau21_b1_mmdt.push_back( nSub2KT_b1(mmdtjet) / nSub1KT_b1(mmdtjet) );
//  	j_tau21_b2_mmdt.push_back( nSub2KT_b2(mmdtjet) / nSub1KT_b2(mmdtjet) );
//  	j_tau32_b1_mmdt.push_back( nSub3KT_b1(mmdtjet) / nSub2KT_b1(mmdtjet) );
//  	j_tau32_b2_mmdt.push_back( nSub3KT_b2(mmdtjet) / nSub2KT_b2(mmdtjet) );
//  	j_c1_b0_mmdt.push_back( ECF_C1_b0(mmdtjet) );
//  	j_c1_b1_mmdt.push_back( ECF_C1_b1(mmdtjet) );
//  	j_c1_b2_mmdt.push_back( ECF_C1_b2(mmdtjet) );
//  	j_c2_b1_mmdt.push_back( ECF_C2_b1(mmdtjet) );
//  	j_c2_b2_mmdt.push_back( ECF_C2_b2(mmdtjet) );
//  	j_d2_b1_mmdt.push_back( ECF_D2_b1(mmdtjet) );
//  	j_d2_b2_mmdt.push_back( ECF_D2_b2(mmdtjet) );
//  	j_d2_a1_b1_mmdt.push_back( ECF_D2_a1_b1(mmdtjet) );
//  	j_d2_a1_b2_mmdt.push_back( ECF_D2_a1_b2(mmdtjet) );
//  	j_m2_b1_mmdt.push_back( ECF_M2_b1(mmdtjet) );
//  	j_m2_b2_mmdt.push_back( ECF_M2_b2(mmdtjet) );
//  	j_n2_b1_mmdt.push_back( ECF_N2_b1(mmdtjet) );
//  	j_n2_b2_mmdt.push_back( ECF_N2_b2(mmdtjet) );
      
//    	j_mass_trim.push_back( trimmer1( curjet ).m() );
//    	j_mass_prun.push_back( pruner1( curjet ).m() );    
//   	j_mass_mmdt.push_back( mmdtjet.m() );
//    	j_mass_sdb2.push_back( soft_drop_sdb2( curjet ).m() );
//    	j_mass_sdm1.push_back( soft_drop_sdm1( curjet ).m() );
      
      } // loop over out_jets

      njets = numGoodJets;

    } // itcounter%NJOBS!=NJOBS-1
    counter++;    

}
////////////////////////////////////////////////////////////////////////////////////////
void analyzeMCEvent(std::vector < fastjet::PseudoJet > MCparticles){

    fastjet::PseudoJet wplus;
    fastjet::PseudoJet wminus;
    fastjet::PseudoJet Zprime;
    for (unsigned int i = 0; i < MCparticles.size(); i++){
        double pdgid =  MCparticles[i].user_index();
            if (pdgid == 24){         
                wplus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
            if (pdgid == -24){
                wminus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
            if (pdgid == 32){
                Zprime.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
            if (pdgid == 25 && wplus.e()<1e-6){         
 	     wplus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
	    else if(pdgid == 25 && wminus.e()<1e-6){         
	      wminus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }

            if (pdgid == 35){
	      Zprime.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
    }
    // fastjet::PseudoJet zprime = wplus + wminus;
    // std::cout << wplus.px() << ", " << wminus.px() << std::endl;
    gen_mZp = Zprime.m();
    gen_mWW = (wplus+wminus).m();
}

