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

// #include "fastjet/contrib/Nsubjettiness.hh"
// #include "fastjet/contrib/Njettiness.hh"
// #include "fastjet/contrib/NjettinessPlugin.hh"
// #include "fastjet/contrib/EnergyCorrelator.hh"
// // #include "fastjet/contrib/ModifiedMassDropTagger.hh"
// #include "fastjet/contrib/SoftDrop.hh"

//#ifdef __MAKECINT__
//#pragma link C++ class vector<float>+;
//#endif

using namespace std;
using namespace fastjet;
//using namespace fastjet::contrib;

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
ifstream finCalo;

int evtCtr;

int njets;
double gen_mZp;
std::vector<float> jpt;
std::vector<float> jp;
std::vector<float> jeta;
std::vector<float> jphi;
std::vector<float> jmass;
std::vector<float> jmultiplicity;
std::vector<float> jisleptag;
std::vector<float> jmass_sd;


////////////////////-----------------------------------------------
void readEvent( std::vector< fastjet::PseudoJet > &allParticles );
void readEventMC( std::vector< fastjet::PseudoJet > &allParticles );
void readEventGEN( std::vector< fastjet::PseudoJet > &allParticles );
void readEventCalo( std::vector< fastjet::PseudoJet > &allParticles );

void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal, std::vector < fastjet::PseudoJet > MCparticles);
void analyzeMCEvent(std::vector < fastjet::PseudoJet > MCparticles);
void clearVectors(){
    //  INIT
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
    std::string type = argv[1];   // type "gg" or "qq"

    std::vector < fastjet::PseudoJet > allParticles;
    std::vector < fastjet::PseudoJet > allParticlesCalo;
    std::vector < fastjet::PseudoJet > allParticlesMC; //status3
    std::vector < fastjet::PseudoJet > allParticlesGEN; //status1

    char fname[150];
    sprintf( fname, "dat/%s.dat", type.c_str() );
    fin.open(fname);

    char fnameMC[150];
    sprintf( fnameMC, "dat/of_status3.dat" );
    finMC.open(fnameMC);

    char fnameGEN[150];
    sprintf( fnameGEN, "dat/of_status1.dat" );
    finGEN.open(fnameGEN);

    char fnameCalo[150];
    sprintf( fnameCalo, "dat/of_CaloHits.dat" );
    finCalo.open(fnameCalo);

    char outName[192];
    sprintf( outName, "dat/%s.root", type.c_str() );
    TFile *f = TFile::Open(outName,"RECREATE");
    TTree *tPFA = new TTree("tPFA","Tree with vectors");
    tPFA->Branch("njets"          , &njets      );
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
    tcalo->Branch("jpt"            , &jpt        );
    tcalo->Branch("jp"             , &jp        );      
    tcalo->Branch("jeta"           , &jeta       );
    tcalo->Branch("jphi"           , &jphi       );
    tcalo->Branch("jmass"          , &jmass      );    
    tcalo->Branch("jmass_sd"       , &jmass_sd      );    
    tcalo->Branch("jmultiplicity"  , &jmultiplicity      );    
    tcalo->Branch("jisleptag"      , &jisleptag      );    

    TTree *tGEN = new TTree("tGEN","Tree with vectors");
    tGEN->Branch("njets"          , &njets      );
    tGEN->Branch("jpt"            , &jpt        );
    tGEN->Branch("jp"             , &jp        );          
    tGEN->Branch("jeta"           , &jeta       );
    tGEN->Branch("jphi"           , &jphi       );
    tGEN->Branch("jmass"          , &jmass      );    
    tGEN->Branch("jmass_sd"       , &jmass_sd      );    
    tGEN->Branch("jmultiplicity"  , &jmultiplicity      );    
    tGEN->Branch("jisleptag"      , &jisleptag      );        

    TTree *tMC = new TTree("tMC","Tree with vectors");
    tMC->Branch("gen_mZp"          , &gen_mZp      );

    int ctr = 0;
    while(true){

        readEvent( allParticles );
        readEventGEN( allParticlesGEN ); 
        readEventMC( allParticlesMC );
        readEventCalo( allParticlesCalo );
        // std::cout << "size of collection = " << allParticles.size() << "," << allParticlesMC.size() << ", " << allParticlesCalo.size() << std::endl;

        if (ctr > 0){
            clearVectors();
            analyzeEvent( allParticles, 0.4, allParticlesMC );
            tPFA->Fill();
            clearVectors();
            analyzeEvent( allParticlesCalo, 0.4, allParticlesMC );
            tcalo->Fill();    
            clearVectors();
            analyzeEvent( allParticlesGEN, 0.4, allParticlesMC );
            tGEN->Fill();  
            analyzeMCEvent( allParticlesMC );
            tMC->Fill();
        }

        allParticles.clear();
        allParticlesMC.clear();
        allParticlesGEN.clear();
        allParticlesCalo.clear();
        
        ctr++;
        std::cout << "ctr = " << ctr << std::endl;

        if(fin.eof()) break;
    }

    f->cd();
    tPFA->Write();
    tcalo->Write();
    tGEN->Write();
    tMC->Write();
    f->Close();

}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEvent( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = 0.4;
    
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
    float etaMax = 0.4;
    
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
    float etaMax = 0.4;
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
    
    int ctr = 0;
    while(true){
        
        finGEN  >> pdgid >> px >> py >> pz >> e;         
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
        
        if(finGEN.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventCalo( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = 0.4;
    
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
void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal, std::vector < fastjet::PseudoJet > MCparticles){

    std::cout << "analyzing event..." << particles.size() << std::endl;
    double rParam = rVal;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    

    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    
    fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(particles, jetDef, fjAreaDefinition);
    // fastjet::ClusterSequenceArea* caloClustering = new fastjet::ClusterSequenceArea(caloclusters, jetDef, fjAreaDefinition);
    std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(25.0));
    // std::vector<fastjet::PseudoJet> calo_jets = sorted_by_pt(caloClustering->inclusive_jets(25.0));

    // double beta_sd = 1.0;
    // double zcut_sd = 0.1;
    // double mu_sd   = 1.0;
    // fastjet::contrib::SoftDrop soft_drop_mmdt(0.0, zcut_sd, mu_sd);

    njets = out_jets.size();
    // std::cout << "number of high pT jets! = " << njets << std::endl;
    for (unsigned int i = 0; i < out_jets.size(); i++){

        double isleptag = 0;
        double minDR = 9999.;
        for (unsigned int j = 0; j < MCparticles.size(); j++){
            // std::cout << "MCparticles = " << MCparticles[j].pt() << "," << MCparticles[j].user_index() << std::endl;
            double pdgid =  fabs(MCparticles[j].user_index());
            double dr = sqrt( pow(MCparticles[j].phi()-out_jets[i].phi(),2) + pow(MCparticles[j].eta()-out_jets[i].eta(),2) );
            // std::cout << "dr = " << dr << "," << pdgid << std::endl;
            if (minDR > dr && (fabs(pdgid) >= 11) && (fabs(pdgid) < 20) ){ minDR = dr; }
        }
        // std::cout<< "minDR = " << minDR << std::endl;
        if (minDR < 0.8){ isleptag = 1.; }

        jpt.push_back( out_jets[i].pt() );
        jp.push_back( sqrt(out_jets[i].modp2()) );
        jeta.push_back( out_jets[i].eta() );
        jphi.push_back( out_jets[i].phi() );
        jmass.push_back( out_jets[i].m() );  
        // jmass_sd.push_back( soft_drop_mmdt( out_jets.at(i) ).m() );        
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
        // if (pdgid == 24){         
        //     wplus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
        // };
        // if (pdgid == -24){
        //     wminus.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
        // };
        // if (pdgid == 32){
        //     Zprime.reset( MCparticles[i].px(), MCparticles[i].py(), MCparticles[i].pz(), MCparticles[i].e() );
        // };
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

