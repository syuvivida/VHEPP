#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
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


//! ========================================================================================

////////////////////-----------------------------------------------
// Global variables
ifstream finGEN_nonu;

std::vector<fastjet::PseudoJet> ISRPhotons;
std::map<unsigned int, std::string> PDGIDs;

const int nparts=2300;
const int lastp = nparts-1;
TH1F* hjeratio[nparts]; 


float jetRadius;
const float etamax=99999.0;
//const float etamax=1.1;


vector<float> trueje;
vector<float> jeratio[nparts];
unsigned int njets;

////////////////////-----------------------------------------------
typedef std::map<unsigned int,std::string>::iterator it_type;

void readEventGEN_nonu( std::vector< fastjet::PseudoJet > &allParticles );
void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal);
void clearVectors(){
  njets = 0;
  trueje.clear();
  for(int i=0; i<nparts;i++)jeratio[i].clear();
}


////////////////////////////////////////////////////////////////////////////////////////
int main (int argc, char **argv) {
    

  PDGIDs[22]  = "pho";
  PDGIDs[11]  = "ele";
  PDGIDs[13]  = "muo";
  PDGIDs[211] = "pion";
  PDGIDs[321] = "kaon";
  PDGIDs[130] = "KL";
  PDGIDs[310] = "Ks";
  PDGIDs[2212] = "p";
  PDGIDs[2112] = "n";


  // std::cout << "hello world" << std::endl;
  std::string type = "jetContent";   // type "gg" or "qq"
  std::string inputFolder = argv[1];
  jetRadius = (float)atof(argv[2]);
  
  std::cout << "jetRadius = " << jetRadius << std::endl;
    
  
  for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) {
    hjeratio[iterator->first] = new TH1F(Form("h_jeratio_%s",iterator->second.data()), iterator->second.data(), 50,0,1);      
  }

  hjeratio[lastp] = new TH1F(Form("h_jeratio_others"), "Other particles", 50,0,1);      


  // pick only the highest pt particle

  std::vector < fastjet::PseudoJet > allParticlesGEN_nonu; //status1

  char fnameGEN[150];
  sprintf( fnameGEN, "%s/of_status1.dat", inputFolder.c_str());
  finGEN_nonu.open(fnameGEN);
    

  char outName[192];
  sprintf( outName, "%s/radius%.1f_%s.root", inputFolder.c_str(), jetRadius, type.c_str() );
  TFile *f = TFile::Open(outName,"RECREATE");

  TTree *tGEN_nonu = new TTree("tGEN_nonu","Tree with vectors");
  tGEN_nonu->Branch("njets", &njets);
  tGEN_nonu->Branch("trueje"         , &trueje      );
  for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) 
    tGEN_nonu->Branch(Form("jeratio_%s",iterator->second.data()), &jeratio[iterator->first]);
  tGEN_nonu->Branch("jeratio_others", &jeratio[lastp]);
 



  int ctr = 0;
  while(true){

    readEventGEN_nonu( allParticlesGEN_nonu );
    
    // std::cout << "size of collection = " << allParticles.size() << "," << allParticlesMC.size() << ", " << allParticlesCalo.size() << std::endl;
    if (ctr > 0)
      {
	clearVectors();
	analyzeEvent( allParticlesGEN_nonu, jetRadius);
	tGEN_nonu->Fill();
      }
    allParticlesGEN_nonu.clear();

	
    ctr++;
    std::cout << "ctr = " << ctr << std::endl;

    if(finGEN_nonu.eof()) break;
  }

  for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) 
    hjeratio[iterator->first]->Write();
  
  hjeratio[lastp]->Write();
  tGEN_nonu->Write();
  

  f->cd();
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
	//	std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
	
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
  std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(25.0));
 
  ISRPhotons.clear();	
  for (unsigned int i = 0; i < out_jets.size() ; i++)
    {
      if(out_jets[i].constituents().size()!=1)continue;
      if(out_jets[i].constituents()[0].user_index()!=22)continue;
      ISRPhotons.push_back(out_jets[i]);	    	    
    }
	

  double beta_sd = 1.0;
  double zcut_sd = 0.1;
  double mu_sd   = 1.0; // for mass drop, equivalent to no cut
  fastjet::contrib::SoftDrop soft_drop_mmdt(0.0, zcut_sd, mu_sd);
    
  // take only the two leading jets for Z' and H2 decays
  int numGoodJets=0;

  for (unsigned int i = 0; i < out_jets.size() && numGoodJets<2; i++){

    float je[nparts]={0.};

    // check if this jet is matched to an ISR photon
    bool isISRPhoton=false;
    for(unsigned int k =0; k < ISRPhotons.size(); k++){
	
      TLorentzVector temp_pho(0,0,0,0);
      temp_pho.SetPxPyPzE(ISRPhotons[k].px(),
			  ISRPhotons[k].py(),
			  ISRPhotons[k].pz(),
			  ISRPhotons[k].e());
      //      cout << "phi = " << temp_pho.Phi() << "\t" << ISRPhotons[k].phi() << endl;
      

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

    if(fabs(out_jets[i].eta())>1.1)continue;

    // check the constituents

    for (unsigned int ip = 0; ip < out_jets[i].constituents().size() ; ip++)
      {
 	unsigned int pdg = abs(out_jets[i].constituents()[ip].user_index());
	if(PDGIDs.find(pdg) != PDGIDs.end())
	  je[pdg] += out_jets[i].constituents()[ip].e();
	else
	  je[lastp] += out_jets[i].constituents()[ip].e();
      }
    
    
    float sume=0;
    for(it_type iterator = PDGIDs.begin(); iterator != PDGIDs.end(); iterator++) {

      float r = je[iterator->first]/out_jets[i].e();
      hjeratio[iterator->first]->Fill(r);
      jeratio[iterator->first].push_back(r);
      sume += je[iterator->first];
      //      std::cout << "energy of " << iterator->second << " = " << je[iterator->first] << std::endl;
    }
    //    std::cout << "sum e = " << sume << "\t jet energy = " << out_jets[i].e() << std::endl;

    hjeratio[lastp]->Fill(je[lastp]/out_jets[i].e());
    jeratio[lastp].push_back(je[lastp]/out_jets[i].e());

    trueje.push_back(out_jets[i].e());


    numGoodJets++;   
  } // number of good jets < 2

  njets = numGoodJets;

}
