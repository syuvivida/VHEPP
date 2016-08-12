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
#include "fastjet/contrib/SoftDrop.hh"


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
ifstream finGEN_nonu;
ifstream finCalo;

std::vector<fastjet::PseudoJet> ISRPhotons;

int counter;
int evtCtr;

int njets;
double gen_mZp;
double gen_mWW;

std::vector<float> jepho;
std::vector<float> je;
std::vector<float> jpt;
std::vector<float> jp;
std::vector<float> jeta;
std::vector<float> jphi;
std::vector<float> jmass;
std::vector<float> jmultiplicity;
std::vector<float> jmass_sd;
std::vector<float> jmass_pr;



float jetRadius;
const float etamax=99999.0;
//const float etamax=1.1;

////////////////////-----------------------------------------------
void readEventGEN_nonu( std::vector< fastjet::PseudoJet > &allParticles );
void readEventCalo( std::vector< fastjet::PseudoJet > &allParticles );

void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal);
void clearVectors(){
    //  INIT
  njets = 0;
  jepho.clear();
  je.clear();
  jpt.clear();
  jp.clear();    
  jeta.clear();
  jphi.clear();        
  jmass.clear();
  jmass_sd.clear();        
  jmass_pr.clear();        
  jmultiplicity.clear();

}



////////////////////////////////////////////////////////////////////////////////////////
int main (int argc, char **argv) {
    
    // std::cout << "hello world" << std::endl;
    std::string type = "of_PanPFA";   // type "gg" or "qq"
    std::string inputFolder = argv[1];
    jetRadius = (float)atof(argv[2]);

    std::cout << "jetRadius = " << jetRadius << std::endl;
    

    // pick only the highest pt particle
    std::vector < fastjet::PseudoJet > allParticlesCalo; // calo clusters matched to PF candidate
    std::vector < fastjet::PseudoJet > allParticlesGEN_nonu; //status1

    char fnameGEN[150];
    sprintf( fnameGEN, "%s/of_status1.dat", inputFolder.c_str());
    finGEN_nonu.open(fnameGEN);

    char fnameCalo[150];
    sprintf( fnameCalo, "%s/%s_Calo.dat", inputFolder.c_str(), type.c_str());
    finCalo.open(fnameCalo);

    char outName[192];
    sprintf( outName, "%s/radius%.1f_ECALHCAL.root", inputFolder.c_str(), jetRadius);
    TFile *f = TFile::Open(outName,"RECREATE");

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
    tGEN_nonu->Branch("jmass_pr"       , &jmass_pr      );    
    tGEN_nonu->Branch("jmultiplicity"  , &jmultiplicity      );    



    int ctr = 0;
    counter =0;
    while(true){
      
      
      readEventGEN_nonu( allParticlesGEN_nonu );
      readEventCalo( allParticlesCalo );

      // std::cout << "size of collection = " << allParticles.size() << "," << allParticlesMC.size() << ", " << allParticlesCalo.size() << std::endl;
      if (ctr > 0){
	
	clearVectors();
	analyzeEvent( allParticlesGEN_nonu, jetRadius);
	tGEN_nonu->Fill();
	    
	clearVectors();
	analyzeEvent( allParticlesCalo, jetRadius);
	tcalo->Fill();    

      }
      

      allParticlesGEN_nonu.clear();
      allParticlesCalo.clear();

      
      ctr++;
      std::cout << "ctr = " << ctr << std::endl;
      
      if(finGEN_nonu.eof()) break;
      //	if(ctr==100)break;
    }
    
    f->cd();

    tcalo->Write();
    tGEN_nonu->Write();
    f->Close();

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
void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal){

  //    std::cout << "analyzing event..." << particles.size() << std::endl;
    double rParam = rVal;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    

    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    
    fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(particles, jetDef, fjAreaDefinition);
    // taking minium E = 0 for pion gun
    std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(25.0));
 
    if(counter%2==0)
      {
	ISRPhotons.clear();	
	for (unsigned int i = 0; i < out_jets.size() ; i++)
	  {
	    if(out_jets[i].constituents().size()!=1)continue;
	    if(out_jets[i].constituents()[0].user_index()!=22)continue;
	    ISRPhotons.push_back(out_jets[i]);	    	    
	  }
	
      }

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


      if(counter%2==0){

	double epho = 0;
	for (unsigned int ip = 0; ip < out_jets[i].constituents().size() ; ip++)
	  {
	    unsigned int pdg = abs(out_jets[i].constituents()[ip].user_index());
	    if(pdg!=22)continue;
	    epho +=out_jets[i].constituents()[ip].e();	      
	  }
	jepho.push_back(epho);

      }


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
      jmultiplicity.push_back( (float) out_jets.at(i).constituents().size() );

      
      
    }
    njets = numGoodJets;

    counter++;    

}
