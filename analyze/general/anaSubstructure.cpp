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
 
 c++ -o anaSubstructure `root-config --glibs --cflags` `/Users/ntran/Research/CMS/PhysicsAnalysis/boost2013/fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -lvectorDict -lEnergyCorrelator anaSubstructure.cpp
 
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
ifstream finCalo;
ifstream finTrack;

int evtCtr;

int njets;
double gen_mZp;
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

////////////////////-----------------------------------------------
void readEvent( std::vector< fastjet::PseudoJet > &allParticles );
void readEventMC( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN_charged( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN_nonu( std::vector< fastjet::PseudoJet > &allParticles );
void readEventCalo( std::vector< fastjet::PseudoJet > &allParticles );
void readEventTrack(  std::vector< fastjet::PseudoJet > &allParticles );

void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal, std::vector < fastjet::PseudoJet > MCparticles);
void analyzeMCEvent(std::vector < fastjet::PseudoJet > MCparticles);
void clearVectors(){
    //  INIT
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

    // pick only the highest pt particle
    std::vector < fastjet::PseudoJet > allParticles;
    std::vector < fastjet::PseudoJet > allParticlesCalo; // calo clusters matched to PF candidate
    std::vector < fastjet::PseudoJet > allParticlesTrack; // tracks matched to PF candidate
    std::vector < fastjet::PseudoJet > allParticlesMC; //status3
    std::vector < fastjet::PseudoJet > allParticlesGEN; //status1 //etamax < 1.1
    std::vector < fastjet::PseudoJet > allParticlesGEN_charged; //status1 // only charged
    std::vector < fastjet::PseudoJet > allParticlesGEN_nonu; //status1

    char fname[150];
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

    char fnameCalo[150];
    sprintf( fnameCalo, "%s/%s_Calo.dat", inputFolder.c_str(), type.c_str());
    finCalo.open(fnameCalo);

    char fnameTrack[150];
    sprintf( fnameTrack, "%s/%s_Tracks.dat", inputFolder.c_str(), type.c_str());
    finTrack.open(fnameTrack);

    char outName[192];
    sprintf( outName, "%s/radius%.1f_%s.root", inputFolder.c_str(), jetRadius, type.c_str() );
    TFile *f = TFile::Open(outName,"RECREATE");
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


    TTree *tGEN_charged = new TTree("tGEN_charged","Tree with vectors");
    tGEN_charged->Branch("njets"          , &njets      );
    tGEN_charged->Branch("je"             , &je         );
    tGEN_charged->Branch("jpt"            , &jpt        );
    tGEN_charged->Branch("jp"             , &jp        );          
    tGEN_charged->Branch("jeta"           , &jeta       );
    tGEN_charged->Branch("jphi"           , &jphi       );
    tGEN_charged->Branch("jmass"          , &jmass      );    
    tGEN_charged->Branch("jmass_sd"       , &jmass_sd      );    
    tGEN_charged->Branch("jmultiplicity"  , &jmultiplicity      );    
    tGEN_charged->Branch("jisleptag"      , &jisleptag      );        

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


    TTree *tMC = new TTree("tMC","Tree with vectors");
    tMC->Branch("gen_mZp"          , &gen_mZp      );

    int ctr = 0;
    while(true){

        readEvent( allParticles );
        readEventGEN( allParticlesGEN ); 
 	readEventGEN_charged( allParticlesGEN_charged ); 
	readEventGEN_nonu( allParticlesGEN_nonu );

        readEventMC( allParticlesMC );
        readEventCalo( allParticlesCalo );
	readEventTrack( allParticlesTrack);

        // std::cout << "size of collection = " << allParticles.size() << "," << allParticlesMC.size() << ", " << allParticlesCalo.size() << std::endl;
        if (ctr > 0){
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
	    analyzeEvent( allParticlesGEN, jetRadius, allParticlesMC );
            tGEN->Fill();
	 
 	    clearVectors();
 	    analyzeEvent( allParticlesGEN_charged, jetRadius, allParticlesMC );
	    tGEN_charged->Fill();
 
 	    clearVectors();
 	    analyzeEvent( allParticlesGEN_nonu, jetRadius, allParticlesMC );
	    tGEN_nonu->Fill();
  
            analyzeMCEvent( allParticlesMC );
            tMC->Fill();
        }

        allParticles.clear();
        allParticlesMC.clear();
        allParticlesGEN.clear();
        allParticlesGEN_charged.clear();
        allParticlesGEN_nonu.clear();
        allParticlesCalo.clear();
        allParticlesTrack.clear();
	
        ctr++;
	std::cout << "ctr = " << ctr << std::endl;

        if(fin.eof()) break;
    }

    f->cd();
    tPFA->Write();
    tcalo->Write();
    ttrack->Write();
    tGEN->Write();
    tGEN_charged->Write();
    tGEN_nonu->Write();
    tMC->Write();
    f->Close();

}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEvent( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = 1.1;
    
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
    float etaMax = 1.1;
    
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
        
        if(finMC.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventGEN( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = 1.1; // to check how good the generator-level jets are
    
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
        
        if(finGEN.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventGEN_charged( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = 1.1; // to check how good the generator-level jets are
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;

    int ctr = 0;
    while(true){
        
        finGEN_charged >> pdgid >> px >> py >> pz >> e;         
        // std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
    
        if (pdgid == -99){
            return;
        }        
 	int pdg=abs(pdgid);
  	if(pdg== 22 || pdg== 12 || pdg== 14 || pdg== 16 
  	   || pdg== 130 || pdg== 310 || pdg== 311 || pdg == 2112 ||
  	   pdg== 3122){
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
    float etaMax = 1.1; // to check how good the generator-level jets are
    
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
	//	std::cout << pdg << std::endl;
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
void readEventCalo( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = 1.1;
    
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
    float etaMax = 1.1;
    
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

    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    
    fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(particles, jetDef, fjAreaDefinition);
    // fastjet::ClusterSequenceArea* caloClustering = new fastjet::ClusterSequenceArea(caloclusters, jetDef, fjAreaDefinition);
    // taking minium E = 0 for pion gun
    //    std::vector<fastjet::PseudoJet> out_jets = sorted_by_E(thisClustering->inclusive_jets(0.0));

    std::vector<fastjet::PseudoJet> out_jets = sorted_by_E(thisClustering->inclusive_jets(25.0));
//    std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(1000.0));

    double beta_sd = 1.0;
    double zcut_sd = 0.1;
    double mu_sd   = 1.0; // for mass drop, equivalent to no cut
    fastjet::contrib::SoftDrop soft_drop_mmdt(0.0, zcut_sd, mu_sd);

    njets = out_jets.size();
//     if(njets>2)
//       std::cout << "number of high pT jets = " << njets << std::endl;
    // take only the two leading jets for Z' and H2 decays
    for (unsigned int i = 0; i < out_jets.size() && i < 2; i++){
    // take only the first-leading jet for pion gun
    //    for (unsigned int i = 0; i < out_jets.size() && i < 1; i++){

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
}
////////////////////////////////////////////////////////////////////////////////////////
void analyzeMCEvent(std::vector < fastjet::PseudoJet > MCparticles){

    fastjet::PseudoJet wplus;
    fastjet::PseudoJet wminus;
    fastjet::PseudoJet Zprime;
    for (unsigned int i = 0; i < MCparticles.size(); i++){
        double pdgid =  MCparticles[i].user_index();
//           if (pdgid == 24){         
//               wplus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
//           };
//           if (pdgid == -24){
//               wminus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
//           };
//           if (pdgid == 32){
//               Zprime.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
//           };
          if (pdgid == 25){         
              wplus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
          };
          if (pdgid == 35){
              Zprime.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
          };
    }
    // fastjet::PseudoJet zprime = wplus + wminus;
    // std::cout << wplus.px() << ", " << wminus.px() << std::endl;
    gen_mZp = Zprime.m();

}

