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
// #include "fastjet/contrib/ModifiedMassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"

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

int counter;
int evtCtr;

int njets;
double gen_mZp;
double gen_mWW;
std::vector<float> je;
std::vector<float> jpt;
std::vector<float> jp;
std::vector<float> jeta;
std::vector<float> jphi;
std::vector<float> jmass;
std::vector<float> jmultiplicity;
std::vector<float> jisleptag;
std::vector<float> jmass_sd;




float jetRadius;
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

void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal, std::vector < fastjet::PseudoJet > MCparticles);
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
    sprintf( fnameCalo, "%s/%s_Calo.dat", inputFolder.c_str(), type.c_str());
    finCalo.open(fnameCalo);

    sprintf( fnameCalo, "%s/of_ECaloHits_r.dat", inputFolder.c_str());
    finECalHits.open(fnameCalo);
    finECalHits2.open(fnameCalo);

    sprintf( fnameCalo, "%s/of_HCaloHits_r.dat", inputFolder.c_str());
    finHCalHits.open(fnameCalo);
    finHCalHits2.open(fnameCalo);

    char fnameTrack[150];
    sprintf( fnameTrack, "%s/%s_Tracks.dat", inputFolder.c_str(), type.c_str());
    finTrack.open(fnameTrack);

    char outName[192];
    sprintf( outName, "%s/radius%.1f_rawhit.root", inputFolder.c_str(), jetRadius);
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

    TTree *trawhits2 = new TTree("trawhits2","Tree with vectors");
    trawhits2->Branch("njets"          , &njets      );
    trawhits2->Branch("je"             , &je         );
    trawhits2->Branch("jpt"            , &jpt        );
    trawhits2->Branch("jp"             , &jp        );      
    trawhits2->Branch("jeta"           , &jeta       );
    trawhits2->Branch("jphi"           , &jphi       );
    trawhits2->Branch("jmass"          , &jmass      );    
    trawhits2->Branch("jmass_sd"       , &jmass_sd      );    
    trawhits2->Branch("jmultiplicity"  , &jmultiplicity      );    
    trawhits2->Branch("jisleptag"      , &jisleptag      );    



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
    tGEN->Branch("jpt"            , &jpt        );
    tGEN->Branch("jp"             , &jp        );          
    tGEN->Branch("jeta"           , &jeta       );
    tGEN->Branch("jphi"           , &jphi       );
    tGEN->Branch("jmass"          , &jmass      );    
    tGEN->Branch("jmass_sd"       , &jmass_sd      );    
    tGEN->Branch("jmultiplicity"  , &jmultiplicity      );    
    tGEN->Branch("jisleptag"      , &jisleptag      );        



    TTree *tGEN_nonu = new TTree("tGEN_nonu","Tree with vectors");
    tGEN_nonu->Branch("njets"          , &njets      );
    tGEN_nonu->Branch("je"             , &je         );
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
 	    analyzeEvent( allParticlesGEN, jetRadius, allParticlesMC );
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
            analyzeEvent( allParticlesCalo, jetRadius, allParticlesMC );
            tcalo->Fill();    

            clearVectors();
	    analyzeEvent( allParticlesTrack, jetRadius, allParticlesMC );
	    ttrack->Fill();

	    clearVectors();
	    analyzeEvent( allParticlesRawHits, jetRadius, allParticlesMC );
	    trawhits->Fill();

	    clearVectors();
 	    analyzeEvent( allParticlesRawHits2, jetRadius, allParticlesMC );
 	    trawhits2->Fill();
  
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
	//	if(ctr==100)break;
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
    trawhits->Write();
    trawhits2->Write();
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
	  RMS=(0.15/sqrt(e) + 0.01)*response;
	  break;
	default:
	  const float A = 1.00647e+00;
	  const float B = -6.45952e+01;
	  const float C = 5.68083e+01;
	  response = A*TMath::Erf((e-B)/C);
	  RMS=(0.38/sqrt(e) + 0.01)*response;	
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
// 	if(cp4.E()<1.0)
// 	  {
// 	    if(finECalHits.eof()) break;
// 	    else continue;
// 	  }
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

    
    ctr = 0;
    //    std::cout << "reading HCAL hits" << std::endl;
    while(true){
        
        finHCalHits  >> pdgid >> x >> y >> z >> e;         
    
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
// 	    if(finHCalHits.eof()) break;
// 	    else continue;
// 	  }
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
        // std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
    
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
        
        if(finECalHits2.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }

    
    ctr = 0;
    while(true){
        
        finHCalHits2  >> pdgid >> x >> y >> z >> e;         
    
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
void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal, std::vector < fastjet::PseudoJet > MCparticles){

  //    std::cout << "analyzing event..." << particles.size() << std::endl;
    double rParam = rVal;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    
    if(counter%8==0)
      for(unsigned int i=0; i < particles.size(); i++)       
	heta_after1->Fill(particles[i].eta());
      

    // doing my own clustering
    if(counter%8==7){
      fastjet::PseudoJet w[2];
      for (unsigned int i = 0; i < MCparticles.size(); i++){
        double pdgid =  MCparticles[i].user_index();
            if (pdgid == 24){         
                w[0].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
            if (pdgid == -24){
                w[1].reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
            }
      }

      if(w[0].e()<1e-6 && w[1].e()<1e-6)
	{
	  std::cout << "W not found!" << std::endl;
	  return;
	}

      fastjet::PseudoJet Wjet[2];
      Wjet[0].reset(0,0,0,0);
      Wjet[1].reset(0,0,0,0);
      int nPart[2]={0,0};
      for(unsigned int i=0; i < particles.size(); i++)
	{
	  for(unsigned int iw=0; iw < 2; iw++)
	    {
	      if(w[iw].e()<1e-6)continue;
	      double dr = sqrt( pow(particles[i].phi()-w[iw].phi(),2) + pow(particles[i].eta()-w[iw].eta(),2) );
	      if(dr > jetRadius)continue;

	      nPart[iw]++;
	      Wjet[iw] += particles[i];
	    }
	}


      for(unsigned int iw=0; iw < 2; iw++)
	{
	  if(Wjet[iw].e()<1e-6)continue;
	  njets++;
	  je.push_back( Wjet[iw].e() );
	  jpt.push_back( Wjet[iw].pt() );
	  jp.push_back( sqrt(Wjet[iw].modp2()) );
	  jeta.push_back( Wjet[iw].eta() );
	  jphi.push_back( Wjet[iw].phi() );
	  jmass.push_back( Wjet[iw].m() );  
	  jmass_sd.push_back( Wjet[iw].m() );        
	  jmultiplicity.push_back(nPart[iw]);
	  jisleptag.push_back(0);
	}

    } // if counter%8==7
    else{

      int activeAreaRepeats = 1;
      double ghostArea = 0.01;
      double ghostEtaMax = 7.0;
    
      fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
      fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    
      fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(particles, jetDef, fjAreaDefinition);
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

      else if(counter%8==3)
	for (unsigned int i = 0; i < out_jets.size() ; i++)
	  heta_PF->Fill(out_jets[i].eta());




      double beta_sd = 1.0;
      double zcut_sd = 0.1;
      double mu_sd   = 1.0; // for mass drop, equivalent to no cut
      fastjet::contrib::SoftDrop soft_drop_mmdt(0.0, zcut_sd, mu_sd);


      // take only the two leading jets for Z' and H2 decays
      int numGoodJets=0;

      for (unsigned int i = 0; i < out_jets.size() && numGoodJets<2; i++){

	// check if this jet is matched to an ISR photon
	bool isISRPhoton=false;
	for(unsigned int k =0; k < ISRPhotons.size(); k++){
	
	  double dr = sqrt( pow(ISRPhotons[k].phi()-out_jets[i].phi(),2) + pow(ISRPhotons[k].eta()-out_jets[i].eta(),2) );
	  double Eratio = ISRPhotons[k].E()>1e-6? out_jets[i].E()/ISRPhotons[k].E(): -1;
	  if(dr < 0.1 && Eratio > 0.7 && Eratio < 1.3)
	    {
	      isISRPhoton=true;
	      break;
	    }
      
	}
      
	if(isISRPhoton)continue;
	numGoodJets++;

	double isleptag = 0;
	double minDR = 9999.;
	for (unsigned int j = 0; j < MCparticles.size(); j++){
	  // std::cout << "MCparticles = " << MCparticles[j].pt() << "," << MCparticles[j].user_index() << std::endl;
	  double pdgid =  fabs(MCparticles[j].user_index());
	  double dr = sqrt( pow(MCparticles[j].phi()-out_jets[i].phi(),2) + pow(MCparticles[j].eta()-out_jets[i].eta(),2) );
	  // std::cout << "dr = " << dr << "," << pdgid << std::endl;
	  if (minDR > dr && (fabs(pdgid)==11 || fabs(pdgid)==13 || fabs(pdgid)==15) ){ minDR = dr; }
	}
	if (minDR < jetRadius){ isleptag = 1.; }
	je.push_back( out_jets[i].E() );
	jpt.push_back( out_jets[i].pt() );
	jp.push_back( sqrt(out_jets[i].modp2()) );
	jeta.push_back( out_jets[i].eta() );
	jphi.push_back( out_jets[i].phi() );
	jmass.push_back( out_jets[i].m() );  
	jmass_sd.push_back( soft_drop_mmdt( out_jets.at(i) ).m() );        
	jmultiplicity.push_back( (float) out_jets.at(i).constituents().size() );
	jisleptag.push_back( isleptag );
      
      }

    njets = numGoodJets;

    } // itcounter!=7
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

