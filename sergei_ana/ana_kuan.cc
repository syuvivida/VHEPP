#include <iostream>
#include <fstream>
#include <stdlib.h>

// check directory
#include <sys/stat.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include "TMath.h"
#include "TNtupleD.h"
#include "TLorentzVector.h"

#include "time.h"
#include <dirent.h>
#include <string>
#include <vector>
#include <map>
struct stat sb;

#include "lcio.h"
#include <stdio.h>
#include "IO/LCReader.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/ReconstructedParticle.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCCollection.h"
#include "IMPL/LCEventImpl.h"
#include "UTIL/LCTOOLS.h"
#include "UTIL/Operators.h"
#include "UTIL/LCIterator.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/MCParticleImpl.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "EVENT/Cluster.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "Objects/CartesianVector.h"
#include "Objects/Helix.h"
#include "DetectorGeometrySimple.h"
#include "IDDecoder.h"
#include "Pandora/Pandora.h"

#include "LParticle.h"
#include "CParticle.h"
using namespace std;

const double kPI = TMath::Pi();
const double k2PI = 2 * kPI;
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

float EffectiveRadius(PseudoJet Jet, vector<PseudoJet> constituents, double jetR = 0.5)
{
    float Energy = Jet.Et();
    float numerator = 0;
    int size = constituents.size();
    for (int i = 0; i < size; i++)
    {
        if (Jet.delta_R(constituents[i]) > jetR)
            continue;
        numerator += constituents[i].Et() * Jet.delta_R(constituents[i]);
    }
    //cout << Energy << endl;
    //cout << numerator << endl;
    return numerator / Energy;
}

float eccentricity(PseudoJet Jet, vector<PseudoJet> constituents)
{
    unsigned int num = constituents.size();
    double Dphi[num], Deta[num], E[num];
    double etaSum = 0.;
    double phiSum = 0.;
    double eTot = 0.;

    for (unsigned int j = 0; j < num; j++)
    {
        PseudoJet cp = constituents.at(j);
        E[j] = cp.e();
        Dphi[j] = Jet.phi() - cp.phi();
        // if (Dphi[j]>TMath::Pi()) Dphi[j]=2*TMath::Pi()-Dphi[j];
        if (fabs(Dphi[j] - 2. * TMath::Pi()) < fabs(Dphi[j]))
            Dphi[j] -= 2. * TMath::Pi();
        if (fabs(Dphi[j] + 2. * TMath::Pi()) < fabs(Dphi[j]))
            Dphi[j] += 2. * TMath::Pi();
        Deta[j] = Jet.eta() - cp.eta();
        etaSum = etaSum + Deta[j] * E[j];
        phiSum = phiSum + Dphi[j] * E[j];
        eTot = eTot + E[j];
    }
    etaSum = etaSum / eTot;
    phiSum = phiSum / eTot;
    for (unsigned int j = 0; j < num; j++)
    {
        Deta[j] = Deta[j] - etaSum;
        Dphi[j] = Dphi[j] - phiSum;
    }

    double X1 = 0.;
    double X2 = 0;
    for (unsigned int i = 0; i < num; i++)
    {
        X1 += 2. * E[i] * Deta[i] * Dphi[i];
        X2 += E[i] * (Dphi[i] * Dphi[i] - Deta[i] * Deta[i]);
    }

    // variance calculations
    double Theta = .5 * atan(X1 / X2);
    double sinTheta = TMath::Sin(Theta);
    double cosTheta = TMath::Cos(Theta);
    double Theta2 = Theta + 0.5 * TMath::Pi();
    double sinThetaPrime = TMath::Sin(Theta2);
    double cosThetaPrime = TMath::Cos(Theta2);

    double VarX = 0.;
    double VarY = 0.;
    for (unsigned int i = 0; i < num; i++)
    {
        double X = cosTheta * Deta[i] - sinTheta * Dphi[i];
        double Y = sinTheta * Deta[i] + cosTheta * Dphi[i];
        VarX += E[i] * X * X;
        VarY += E[i] * Y * Y;
    }

    double VarianceMax = VarX;
    double VarianceMin = VarY;
    if (VarianceMax < VarianceMin)
    {
        VarianceMax = VarY;
        VarianceMin = VarX;
    }

    double ECC = 1.0 - (VarianceMin / VarianceMax);

    return ECC;
}

double nsubjettiness(PseudoJet Jet, vector<PseudoJet> constituents, int NSubJets, double jetRad = 0.5)
{
    //vector<CParticle> constit = jet.GetConstituents();
    vector<PseudoJet> jetConstit;

    unsigned int num = constituents.size();

    if (num < (unsigned int)NSubJets)
    {
        return -999;
    }
    TLorentzVector Jet_p;
    Jet_p.SetPxPyPzE(Jet.px(), Jet.py(), Jet.pz(), Jet.E());
    num = 0;
    for (unsigned int i = 0; i < constituents.size(); i++)
    {
        PseudoJet c_i = constituents[i];
        if (c_i.delta_R(Jet_p) > jetRad)
            continue;
        jetConstit.push_back(c_i);
        num++;
    }
    if (num < (unsigned int)NSubJets)
    {
        return -999;
    }
    std::vector<std::vector<fastjet::PseudoJet>> kt_subjets_vec; //a vector of vectors of Pseudojets
    fastjet::JetDefinition *m_jetdef = new fastjet::JetDefinition(fastjet::kt_algorithm, 1.5, fastjet::E_scheme, fastjet::Best);
    fastjet::ClusterSequence kt_seq(jetConstit, *m_jetdef);
    delete m_jetdef;

    for (unsigned int i = 0; i < (unsigned int)NSubJets; i++)
    {
        kt_subjets_vec.push_back(fastjet::sorted_by_pt(kt_seq.exclusive_jets((int)NSubJets)));
    }
    double min_dist = 100000.0;
    double sum_pt = 0.0;
    double sum_dist = 0.0;
    //first find the minimum distance.
    for (unsigned int i = 0; i < jetConstit.size(); i++)
    {
        fastjet::PseudoJet theconstit(jetConstit[i]);
        sum_pt += theconstit.perp();
        float min_dist = 1e10;
        for (unsigned int j = 0; j < (unsigned int)NSubJets; j++)
        {
            const std::vector<fastjet::PseudoJet> *kt_subjets = &(kt_subjets_vec[j]);
            float temp_dist = theconstit.perp() * std::sqrt(kt_subjets->at(j).plain\
_distance(theconstit));
            if (temp_dist < min_dist)
                min_dist = temp_dist;
        } //loop over axis (subjets)
        sum_dist += min_dist;
    } //loop over jet constituents

    sum_dist /= (jetRad * sum_pt);

    double nSubJettiness = sum_dist;
    if (sum_dist > 1.0)
        cout << "uh oh" << sum_dist << endl;

    return nSubJettiness;
}

double splittingscale(PseudoJet Jet)
{
    if (!Jet.has_constituents())
        return -5.;

    vector<PseudoJet> jetConstit = Jet.constituents();

    double dR = Jet.associated_cluster_sequence()->jet_def().R();
    fastjet::JetDefinition *m_jetdef = new fastjet::JetDefinition(fastjet::kt_algorithm, dR, fastjet::E_scheme, fastjet::Best);

    fastjet::ClusterSequence kt_seq(jetConstit, *m_jetdef);
    delete m_jetdef;
    return dR * TMath::Sqrt(kt_seq.exclusive_dmerge(1));
}

// find all files inside a directory
std::vector<std::string> open(std::string name = "data.in")
{
    vector<std::string> ntup;
    ifstream myfile;
    myfile.open(name.c_str(), ios::in);

    if (!myfile)
    {
        cerr << " -> Can't open input file:  " << name << endl;
        exit(1);
    }
    else
    {
        cout << "-> Read data file=" << name << endl;
    }

    string temp;
    while (myfile >> temp)
    {
        //the following line trims white space from the beginning of the string
        temp.erase(temp.begin(), std::find_if(temp.begin(), temp.end(), not1(ptr_fun<int, int>(isspace))));
        if (temp.find("#") == 0)
            continue;
        ntup.push_back(temp);
    }
    cout << "-> Number of files=" << ntup.size() << endl;
    myfile.close();

    for (unsigned int i = 0; i < ntup.size(); i++)
    {
        cout << ".. file to analyse=" + ntup[i] << endl;
    }
    return ntup;
}

int main(int argc, char **argv)
{
    std::vector<std::string> files = open("data.in");
    int mtype = 2;

    string outputfile = "root/output.root";
    cout << "\n -> Output file is =" << outputfile << endl;

    TFile *RootFile = new TFile(outputfile.c_str(), "RECREATE", "Histogram file");
    TH1D *h_debug = new TH1D("debug", "events", 10, 0, 10);

    double TTxbins[] = {0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0, 20.0, 30};
    const int TTBins = sizeof(TTxbins) / sizeof(double);

    TH1D *firstLayerTime = new TH1D("firstLayerTime","First Layer Time",100,7.0,8.0);
    firstLayerTime->GetXaxis()->SetTitle("T[ns]");
    firstLayerTime->GetYaxis()->SetTitle("Entries");

    TH1D *lastLayerTime = new TH1D("lastLayerTime","last Layer Time",100,7.0,8.0);
    lastLayerTime->GetXaxis()->SetTitle("T[ns]");
    lastLayerTime->GetYaxis()->SetTitle("Entries");

    TH1D *OwnLayerTime = new TH1D("OwnLayerTime","Time difference between 31 nd 0 layers for 1 GeV",1000,0.0,1.0);

    TH1D *TimeEcalDiff1 = new TH1D("time_ecal_diff1gev", "Time difference between 31 nd 0 layers for 1 GeV", TTBins - 1, TTxbins);
    TH1D *TimeEcalDiff10 = new TH1D("time_ecal_diff10gev", "Time difference between 31 nd 0 layers for 10 GeV", TTBins - 1, TTxbins);
    TProfile *TimeEcalDiffPt = new TProfile("time_ecal_diff vs pt", "Time difference between 31 nd 0 layers vs pT", 20, 0, 20, 's');

    TH1D *binsTT = new TH1D("bins_time_ecal", "bins_time_ecal", TTBins - 1, TTxbins);
    binsTT->Sumw2();
    for (Int_t j = 0; j < TTBins - 1; j++)
    {
        float x = TTxbins[j + 1] - TTxbins[j];
        binsTT->Fill(TTxbins[j] + 0.5 * x, x);
    }

    // read detector geometry for this configuration
    string detector = "./data/rfull009_sifcch7/sifcch7/sifcch7.pandora";
    DetectorGeometrySimple *geom = new DetectorGeometrySimple(detector);
    string caloType = "HAD_BARREL";
    DetectorGeometrySimple::ExtraSubDetectorParameters *xsubdet = geom->getExtraSubDetectorParametersFromType(caloType);
    if (xsubdet == NULL)
    {
        std::cout << "The ExtraSubDetectorParameters for " << caloType << " were not found." << std::endl;
        throw new std::exception;
    }

    string caloType_ecal = "EM_BARREL";
    DetectorGeometrySimple::ExtraSubDetectorParameters *xsubdet_ecal = geom->getExtraSubDetectorParametersFromType(caloType_ecal);
    if (xsubdet_ecal == NULL)
    {
        std::cout << "The ExtraSubDetectorParameters for " << caloType_ecal << " were not found." << std::endl;
        throw new std::exception;
    }

    // Get the decoder.
    IDDecoder *decoder = xsubdet->m_decoder;
    IDDecoder *decoder_ecal = xsubdet_ecal->m_decoder;

    // calculate position
    double layer_size_mm = 27.5 + 5 + 2.5;
    double hcal_inner_R_mm = 2300;

    double layer_size_ecal_mm = 4.0;
    double ecal_inner_R_mm = 2100;

    // calculate weighted average of time
    double XTime = 0;
    double XEnergy = 0;

    const pandora::SubDetectorType subDetectorType = geom->getPandoraSubDetectorType(caloType);

    const char *tupleNames = "px:py:pz:vx:vy:vz:ex:ey:ez:pdg:np:nd:gs:ss";
    TNtupleD *tuple = new TNtupleD("MCParticle", "", tupleNames);

    double minPtConst = 0.2; // min pT on constituents
    // jets
    double Rparam = 0.5;
    const double ETAmax = 1.0;
    double ETmin = 0.5;
    if (mtype == 1)
        ETmin = 200; // increase for boosted
    if (mtype == 31)
        ETmin = 50;
    if (mtype == 32)
        ETmin = 2000;
    if (mtype == 33)
        ETmin = 8000;

    if (mtype == 1)
        cout << "mu+mu- mode for boosted jets" << endl;
    if (mtype == 2)
        cout << "single-particle mode for pgun" << endl;
    if (mtype == 3)
        cout << "pp mode for jet resolutions" << endl;

    cout << "min PT for jets=" << ETmin << endl;
    cout << "eta max for jets=" << ETAmax << endl;
    cout << "R for jets =" << Rparam << endl;

    // fastjet
    Strategy strategy = fastjet::Best;
    JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, strategy);

    int MaxEvents = 100000000;

    int nEvents = 0;

    // loop over all files
    for (unsigned int mfile = 0; mfile < files.size(); mfile++)
    {
        string Rfile = files[mfile];

        cout << " # File=" << Rfile << endl;

        IO::LCReader *lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
        lcReader->open(Rfile.c_str());

        EVENT::LCEvent *evt = 0;

        if (nEvents > MaxEvents)
            break;
        //----------- the event loop -----------
        while ((evt = lcReader->readNextEvent()) != 0)
        {
            if (nEvents == 0)
                UTIL::LCTOOLS::dumpEvent(evt);

            // UTIL::LCTOOLS::dumpEvent( evt ) ;

            nEvents++;
            if ((nEvents < 100 && nEvents % 10 == 0) || (nEvents > 100 && nEvents % 200 == 0))
                cout << " # Events=" << nEvents << endl;

            if (nEvents > MaxEvents)
                break;

            h_debug->Fill(1.0);

            std::string mcpName("MCParticle");
            // get truth
            IMPL::LCCollectionVec *col = (IMPL::LCCollectionVec *)evt->getCollection(mcpName);
            int nMCP = col->getNumberOfElements();

            int neu = 0;

            vector<PseudoJet> avec_truth;     // created by generator
            vector<PseudoJet> avec_truth_sim; // also created by geant

            for (int i = 0; i < nMCP; ++i)
            {
                EVENT::MCParticle *mcp = (EVENT::MCParticle *)col->getElementAt(i);
                double px = mcp->getMomentum()[0];
                double py = mcp->getMomentum()[1];
                double pz = mcp->getMomentum()[2];
                double m = mcp->getMomentum()[3];
                double e = sqrt(px * px + py * py + pz * pz + m * m);

                //cout << "m=" << m << " |p|=" << sqrt(px*px+py*py+pz*pz) << " e=" << e << endl;

                TLorentzVector p;
                p.SetPxPyPzE(px, py, pz, e);

                // find high-pT neutrino
                if ((mcp->getPDG() == 12 || mcp->getPDG() == 14 || mcp->getPDG() == 16) && p.Pt() > 100)
                {
                    neu++;
                }

                // includes Geant4 and backscatter too
                bool isGeant4 = false;
                if (abs(p.PseudoRapidity()) < 2.0 && (mcp->isCreatedInSimulation() == true || mcp->isBackscatter() == true) && p.Pt() > minPtConst)
                {
                    if (mcp->getPDG() != 12 && mcp->getPDG() != 14 && mcp->getPDG() != 16)
                        isGeant4 = true;
                }

                // generator-level
                bool isGen = false;
                if (abs(p.PseudoRapidity()) < 2.0 && mcp->getGeneratorStatus() == 1 && p.Pt() > minPtConst)
                {
                    if (mcp->getPDG() != 12 && mcp->getPDG() != 14 && mcp->getPDG() != 16)
                        isGen = true;
                }

                // only generator level
                if (isGen == true)
                {
                    int pdg = mcp->getPDG();
                    int np = mcp->getParents().size();
                    int nd = mcp->getDaughters().size();
                    int gs = mcp->getGeneratorStatus();
                    int ss = mcp->getSimulatorStatus();

                    avec_truth.push_back(PseudoJet(px, py, pz, e));
                }

                // plus Geant4
                if (isGen == true || isGeant4 == true)
                {
                    int pdg = mcp->getPDG();
                    int np = mcp->getParents().size();
                    int nd = mcp->getDaughters().size();
                    int gs = mcp->getGeneratorStatus();
                    int ss = mcp->getSimulatorStatus();

                    avec_truth_sim.push_back(PseudoJet(px, py, pz, e));
                }
            }
            // assume remnant for mu+mu-
            if (mtype == 1)
            {
                avec_truth.push_back(PseudoJet(0, 0, 5000, 5000));
                avec_truth.push_back(PseudoJet(0, 0, -5000, 5000));
            }
            ClusterSequence clust_seq(avec_truth, jet_def);
            vector<PseudoJet> jets_truth = clust_seq.inclusive_jets(ETmin);
            vector<PseudoJet> sjets_truth = sorted_by_pt(jets_truth);
            vector<LParticle> truthjets;
            int nn = 0;
            for (unsigned int k = 0; k < sjets_truth.size(); k++)
            {
                double eta = sjets_truth[k].pseudorapidity();
                double phi = sjets_truth[k].phi();
                if (phi < 0)
                    phi = phi + k2PI;
                double m = jets_truth[k].m();
                double pt = sjets_truth[k].perp();
                double e = sjets_truth[k].e();

                //TLorentzVector p;
                //p.SetPxPyPzE(px,py,pz,e);
                if (pt < ETmin)
                    continue;
                if (fabs(eta) > ETAmax)
                    continue;
                LParticle p(sjets_truth[k].px(), sjets_truth[k].py(), sjets_truth[k].pz(), sjets_truth[k].e(), 0);

                // fill jet substructure
                p.SetParameter(EffectiveRadius(sjets_truth[k], sjets_truth[k].constituents(), Rparam));  // 0
                p.SetParameter(nsubjettiness(sjets_truth[k], sjets_truth[k].constituents(), 1, Rparam)); // 1
                p.SetParameter(nsubjettiness(sjets_truth[k], sjets_truth[k].constituents(), 2, Rparam)); // 2
                p.SetParameter(nsubjettiness(sjets_truth[k], sjets_truth[k].constituents(), 3, Rparam)); //  3
                p.SetParameter(splittingscale(sjets_truth[k]));                                          // 4
                p.SetParameter(eccentricity(sjets_truth[k], sjets_truth[k].constituents()));             // 5
                truthjets.push_back(p);

                // cout << "jet=" << nn << " E=" << e << " eta=" << eta << " phi=" << phi <<  endl;
                nn++;
                if (nn > 2)
                    break; // take first 3 jets
            }

            float TruthPTjet = 0; // keep pT of partice (for single particle studies)

            for (unsigned int j = 0; j < truthjets.size(); j++)
            {
                LParticle p = (LParticle)truthjets.at(j);
                TLorentzVector L = p.GetP();
                double pt = L.Perp();
                double eta = L.Eta();
                double mass = L.M();
                TruthPTjet = pt;
            }

                vector<float> ECALhits0, ECALhits31;

                double ecalsum_raw = 0;
                vector<PseudoJet> avec_hits_raw;
                vector<PseudoJet> avec_hits_raw_sf;
                vector<LParticle> simhits;
                //ECal hits

                IMPL::LCCollectionVec *col53 = (IMPL::LCCollectionVec *)evt->getCollection("EcalBarrelHits");
                int nCL = col53->getNumberOfElements();

                for (int i = 0; i < nCL; i++)
                {
                    EVENT::SimCalorimeterHit *mcp = (EVENT::SimCalorimeterHit *)col53->getElementAt(i);

                    const float *pos = mcp->getPosition();
                    float x = pos[0];
                    float y = pos[1];
                    float z = pos[2];
                    double e = mcp->getEnergy();
                    double _tmp = std::sqrt(x * x + y * y + z * z);
                    double px = e * x / _tmp;
                    double py = e * y / _tmp;
                    double pz = e * z / _tmp;

                    // Get the two 32-bit chunks of the ID.
                    int cellId0 = mcp->getCellID0();
                    int cellId1 = mcp->getCellID1();
                    // Make a 64-bit id for the IDDecoder.  The type MUST be "long long" and not "long".  (from Tony Johnson)
                    long long cellId = ((long long)cellId1) << 32 | cellId0;
                    int layer = decoder_ecal->getFieldValue("layer", cellId);
                    // 1st layer on middle
                    double layer_pos_cm = 0.1 * (layer_size_ecal_mm * 0.5 + (layer * layer_size_ecal_mm));
                    double Thit = mcp->getTimeCont(0);

                    double avt = 0;
                    double ave = 0;
                    for (int jj = 0; jj < mcp->getNMCContributions(); jj++)
                    {

                        if (layer == 0)
                            ECALhits0.push_back(mcp->getTimeCont(jj));
                        if (layer == 31)
                            ECALhits31.push_back(mcp->getTimeCont(jj));
                        avt = avt + mcp->getEnergyCont(jj) * mcp->getTimeCont(jj);
                        ave = ave + mcp->getEnergyCont(jj);
                    }
                    avt = avt / ave;
                    ecalsum_raw = ecalsum_raw + e;

                    PseudoJet pj(px, py, pz, e);
                    double eta_r = pj.pseudorapidity();
                    double phi_r = pj.phi();
                    double pt_r = pj.pt();
                    //cout << " e=" << e <<  " phi=" << pj.phi() << " eta=" << pj.eta() << endl;
                    avec_hits_raw.push_back(pj);

                    // fill hits
                    LParticle p(px, py, pz, e, layer);
                    p.SetCharge(layer);
                    p.SetType(2); // ECAL
                    p.SetStatus(mcp->getNMCContributions());
                    p.SetParameter(x);
                    p.SetParameter(y);
                    p.SetParameter(z);
                    p.SetParameter(layer_pos_cm);

                    // find fastest hit
                    float timeCont = mcp->getTimeCont(0);
                    EVENT::MCParticle *pmm = mcp->getParticleCont(0);
                    int pdg = pmm->getPDG();
                    int status = pmm->getSimulatorStatus();
                    float rawTime = timeCont;
                    int nCont = mcp->getNMCContributions();
                    for (int jjj = 0; jjj < nCont; jjj++)
                    {
                        if (mcp->getTimeCont(jjj) < rawTime)
                            rawTime = mcp->getTimeCont(jjj);
                        EVENT::MCParticle *pmm = mcp->getParticleCont(jjj);
                        pdg = pmm->getPDG();
                        status = pmm->getSimulatorStatus();
                    }
                    p.SetParameter(rawTime); // fastest hit
                    p.SetParameter(avt);     // average hit time
                    p.SetParameter((double)pdg);
                    p.SetParameter((double)status);

                    simhits.push_back(p);
                    // correct by SF (1st, 2nd, 3rd layer are 1., 0.0184, 0.0092)
                    double ECAL_SF = 0.0184;
                    e = e / ECAL_SF;
                    px = e * x / _tmp;
                    py = e * y / _tmp;
                    pz = e * z / _tmp;
                    PseudoJet pj_sf(px, py, pz, e);
                    avec_hits_raw_sf.push_back(pj_sf);
                }
                sort(ECALhits0.begin(), ECALhits0.end(), less<float>());
                sort(ECALhits31.begin(), ECALhits31.end(), less<float>());

                if (ECALhits31.size() > 0 && ECALhits0.size() > 0)
                {
                    float dtime = ECALhits31[0] - ECALhits0[0];
                    float firsttime = ECALhits0[0];
                    float lasttime = ECALhits31[0];
                    //cout << dtime << endl;
                    if (TruthPTjet > 0.8 && TruthPTjet < 1.2)
                        TimeEcalDiff1->Fill(dtime);
                        firstLayerTime->Fill(firsttime);
                        lastLayerTime->Fill(lasttime);
                        OwnLayerTime->Fill(dtime);
                        
                    if (TruthPTjet > 8 && TruthPTjet < 12)
                        TimeEcalDiff10->Fill(dtime);

                    TimeEcalDiffPt->Fill(dtime, TruthPTjet);
                }
            } //end loop
            lcReader->close();
            delete lcReader;
        }
        RootFile->Write();
        //RootFile->Print();
        RootFile->Close();
        cout << "end" << endl;

        return 0;
    }