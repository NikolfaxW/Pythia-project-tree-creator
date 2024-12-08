//
// Created by nikol on 11/14/2024.
//

#include "functions.h"
#include "particleUnit.h"
#include "MyInfo.h"

#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <deque>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

#include "Pythia8/Pythia.h"
#include "TVector2.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "TCanvas.h"
#include "TH2F.h"


void showProgressBar(int progress, int total) {
    double ratio = static_cast<double>(progress) / total;
    std::cout << "(required number of D_0 found: " << progress << " | " << total << ") ";
    std::cout << int(ratio * 100.0) << "%\r";
    std::cout.flush();
}

Double_t delta_R(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
    Double_t delta_eta = eta1 - eta2;
    Double_t delta_phi = TVector2::Phi_mpi_pi(phi1 - phi2);
//    std::cout << "delta_eta: " << delta_eta << " delta_phi: " << delta_phi << std::endl;
    return TMath::Sqrt(pow(delta_eta, 2) + pow(delta_phi, 2));
}

int getId() {
    std::string fileName = "../status.txt";
    std::ifstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Error: file not opened" << std::endl;
        return -1;
    }
    int id;
    std::string temp;
    file >> temp >> id;
    file.close();
    return id;
}

std::string getStatus() {
    std::string fileName = "../status.txt";
    std::ifstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Error: file not opened" << std::endl;
        return "Error";
    }
    std::string temp;
    file >> temp;
    file.close();
    return temp;
}

bool increaseIdOrChageStatus(int id, std::string status) {
    std::string fileName = "../status.txt";
    std::fstream file(fileName, std::ios::out);
    if (!file.is_open()) {
        std::cerr << "Error: file not opened" << std::endl;
        return false;
    }
    file << status << std::endl << (id + 1) << std::endl;
    file.close();
    return true;
}

bool areAlmostEqual(Float_t a, Float_t b, int ulp = 1) {
    // ULP (Units in the Last Place) is the number of steps allowed
    return std::fabs(a - b) <= std::fabs(std::nextafter(a, b) - a) * ulp;
}

int randomSeed(){
    // Seed for reproducibility
    std::random_device rd; // Use random_device for non-deterministic randomness
    std::mt19937 gen(rd()); // Mersenne Twister random number engine

    // Define normal distribution with mean and standard deviation
    double mean = 450000;
    double stddev = 259807.62;
    std::normal_distribution<> dist(mean, stddev);
    return dist(gen);
}

void mainThreadedSec(const int numThreads,
                std::string seed,
                TTree *&T,
                Float_t &D_0_pT,
                Float_t &Jet_Pt,
                Float_t &rapidity,
                Float_t &z_val,
                const unsigned int &requiredNumberOfD_0,
                unsigned int &numberOfD_0Found,
                Float_t &l11,
                Float_t &l105,
                Float_t &l115,
                Float_t &l12,
                Float_t &l13,
                Float_t &l20,
                Float_t &pl11,
                Float_t &pl105,
                Float_t &pl115,
                Float_t &pl12,
                Float_t &pl13,
                Float_t &pl20,
                TTree *&T2,
                TTree *&T3,
                Float_t &isCharged,
                Float_t &deltaR,
                Float_t &pT_frac) {


    Pythia8::Pythia pythia; //create pythia object


    {
        std::lock_guard<std::mutex> lock(std::mutex);  // Lock the mutex
        pythia.readFile("../config1.cmnd"); //read config file and intialize pythia
        pythia.readString("Random:setSeed = on");  // Enable setting of the seed
        pythia.readString("Random:seed = " + seed);
        pythia.init();

    }


    std::list<particleUnit> j_constituents;
    std::map<TString, fastjet::JetDefinition> jetDefs; //map to store jet definitions

    //parameters for jet finding and to make program less hardcoded
    double R = 0.4; //jet radius
    double pTmin_jet = 5, pThadron = 0.2; //minimum pT for jets and hadrons
    int triggerId = 421; //pdg code of the particle to be found

    double pTMinTrig = 0; //minimum pT for the particle to be found
    double pTMaxTrig = 5.0;
    double mTemp; //This variable are needed to recount momentum after particle mass resets
    Pythia8::Vec4 pTemp; //This variable are needed to recount momentum after particle mass resets


    Double_t R_frac, temp_pT_frac, temp_z_val = -1000, temp_D_0_pT, temp_l11, temp_l105, temp_l115, temp_l12, temp_l13, temp_l20;
    Double_t delta_R_temp;
    Float_t temp_pl11, temp_pl105, temp_pl115, temp_pl12, temp_pl13, temp_pl20;
    std::vector<Double_t>pT_fracs, delta_R_fracs;
    pT_fracs.reserve(100);
    delta_R_fracs.reserve(100);

    //define jet finding algorithms here:
    jetDefs["#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
            fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);


    auto &event = pythia.event; //create a reference to the Pythia8 collision event
    std::vector<fastjet::PseudoJet> fjInputs; //to store particles to be clustered
    std::vector<fastjet::PseudoJet> selectedJets; //to store jets after all cuts



    for (int iEvent = 0; numberOfD_0Found < requiredNumberOfD_0; ++iEvent) { //loop over needed number of events
        if (!pythia.next()) continue; //generate next event, if it is not possible, continue
        int idxD = -1; // to store index of the D_0 particle in the event
        for (int i = pythia.event.size() - 1; i > 0; i--) { //goes through all particles generated in event
            if (pythia.event[i].idAbs() == triggerId &&
                pythia.event[i].pT() >= pTMinTrig) { //finds D_0 particle with required pT cut
                idxD = i; //saved its index in the event
                break;
            }
        }
        if (idxD == -1) //if there is no D_0 particle in the event, skip it to not waste resources
            continue;

        fjInputs.clear(); //clears the vector of particles to be clustered
        for (int i = 0; i < event.size(); ++i) { // saves particles in order to make jets
            auto &p = event[i];
            if (event[i].isFinal() || event[i].idAbs() == triggerId) { //checks if particle is final state or D_0
                pTemp = p.p();
                if (p.idAbs() == 22) //changes mass for any other particle except photons to pion mass
                    mTemp = 0;
                else
                    mTemp = 0.13957;
                fastjet::PseudoJet particleTemp = fastjet::PseudoJet(p.px(), p.py(), p.pz(),
                                                                     sqrt(pTemp.pAbs2() + mTemp *
                                                                                          mTemp)); //recounts 4 momentum and energy after mass reset
                particleTemp.set_user_info(
                        new MyInfo(p.idAbs(), i, event[i].isCharged())); //adds additional info to the particle
                fjInputs.push_back(particleTemp); //saves the particle to the vector to be clustered

            }
        }


        if (fjInputs.empty()) { // Abort to avoid running algorithms through empty vectors
            std::cout << "Error: event with no final state particles" << std::endl;
            continue;
        }


        for (auto jetDef: jetDefs) { //for each jet definition runs jet clustering sequence, then saves it to ROOT tree
            selectedJets.clear(); //empties the vector of jets - prevents from saving jets from previous e2vent
            fastjet::ClusterSequence clusterSequence(fjInputs, jetDef.second); //sets up the cluster sequence
            auto jets = sorted_by_pt(clusterSequence.inclusive_jets(0)); //runs the clustering and sorts jets by pT
            fastjet::Selector eta_selector = fastjet::SelectorEtaRange(-1 + R, 1 - R); // Selects jets with |eta| < 1-R
            selectedJets = eta_selector(jets); //applies the eta selector to the jets

            for (const auto &jet: selectedJets) { //loop through all jets
                Bool_t temp_Has_D_0 = false;
                for (const auto &c: jet.constituents()) {//loop through all jet constituents to check if D_0 is there

                    if (c.user_info<MyInfo>().pdg_id() == triggerId) { // Only for D_0
                        temp_Has_D_0 = true;
                        temp_z_val = (jet.px() * c.px() + jet.py() * c.py()) / jet.pt2();
                        temp_D_0_pT = c.pt(); //saves pT of the D_0 particle
                        break;
                    }
                }
                if (not temp_Has_D_0) continue; // if there is not d_0 particle in the jet, skip it;
                temp_l11 = 0;
                temp_l105 = 0;
                temp_l115 = 0;
                temp_l12 = 0;
                temp_l13 = 0;
                temp_l20 = 0;
                temp_pl11 = -1;
                temp_pl105 = -1;
                temp_pl115 = -1;
                temp_pl12 = -1;
                temp_pl13 = -1;
                temp_pl20 = -1;

                delta_R_fracs.clear();
                pT_fracs.clear();
                for (const auto &c: jet.constituents()) { //loop through all jet constituents to calculate l11, l105, l115, l12, l13, l20
                    if(c.user_info<MyInfo>().isCharged()){
                        std::lock_guard<std::mutex> lock(std::mutex);
                        isCharged = 1;
                        T3->Fill();
                    } else{
                        std::lock_guard<std::mutex> lock(std::mutex);
                        isCharged = -1;
                        T3->Fill();
                        if(c.user_info<MyInfo>().pdg_id() != triggerId) continue;
                    }

                    temp_pT_frac = c.pt() / jet.pt();
                    pT_fracs.push_back(temp_pT_frac);
                    delta_R_temp = delta_R(jet.eta(), jet.phi(), c.eta(), c.phi());
                    R_frac = delta_R_temp / R;
                    delta_R_fracs.push_back(delta_R_temp);
                    temp_l11 += pT_frac * R_frac;
                    temp_l105 += pT_frac * pow(R_frac, 0.5);
                    temp_l115 = pT_frac * pow(R_frac, 1.5);
                    temp_l12 += pT_frac * pow(R_frac, 2);
                    temp_l13 += pT_frac * pow(R_frac, 3);
                    temp_l20 += pow(pT_frac, 2);
                }

                if(temp_l11 > 1 ) temp_pl11 = 1;
                if(temp_l105 > 1 ) temp_pl105 = 1;
                if(temp_l115 > 1 ) temp_pl115 = 1;
                if(temp_l12 > 1 ) temp_pl12 = 1;
                if(temp_l13 > 1 ) temp_pl13 = 1;
                if(temp_l20 > 1 ) temp_pl20 = 1;
                if(areAlmostEqual(temp_pl11,1)|| areAlmostEqual(temp_pl105,1) || areAlmostEqual(temp_pl115,1)
                || areAlmostEqual(temp_pl12,1)||areAlmostEqual(temp_pl13,1) || areAlmostEqual(temp_pl20,1)){
                    std::lock_guard<std::mutex> lock(std::mutex);
                    while(!delta_R_fracs.empty()){
                        pT_frac = pT_fracs.back();
                        pT_fracs.pop_back();
                        deltaR = delta_R_fracs.back();
                        delta_R_fracs.pop_back();
                        T2->Fill();
                    }
                }

                {
                    //! Not safe to fill the tree from multiple threads
                    std::lock_guard<std::mutex> lock(std::mutex);  // Lock the mutex
                    l11 = temp_l11;
                    l105 = temp_l105;
                    l115 = temp_l115;
                    l12 = temp_l12;
                    l13 = temp_l13;
                    l20 = temp_l20;
                    pl11 = temp_pl11;
                    pl105 = temp_pl105;
                    pl115 = temp_pl115;
                    pl12 = temp_pl12;
                    pl13 = temp_pl13;
                    pl20 = temp_pl20;
                    D_0_pT = temp_D_0_pT;
                    z_val = temp_z_val;
                    Jet_Pt = jet.pt();
                    rapidity = jet.rapidity();
                    T->Fill();  // Fill the data vector safely
                    ++numberOfD_0Found;
//                    if (numberOfD_0Found % 10 == 0)
                        showProgressBar(numberOfD_0Found, requiredNumberOfD_0);
                }
            }

        }
    }
}

void mainSec(int id){
    unsigned int requiredNumberOfD_0InOneFile = 100000;
    unsigned int foundNumberOfD_0 = 0; //to store number of D_0 particles found
    unsigned int numThreads = std::thread::hardware_concurrency();
    int seed = randomSeed();
    std::thread **alocatedThreads = new std::thread *[numThreads];
    std::string pathToTheFile = "../results/";
    std::string filename = std::to_string(id) + ','+ std::to_string(requiredNumberOfD_0InOneFile) + ',' + std::to_string(seed) + ".root";
    TTree *T = new TTree("T", "saves Pt, pseudo rapidity and phi of an D_0 jet"); //create a tree to store the data
    TTree *T2 = new TTree("T2", "saves important dubug data(delta R, pT_frac)");
    TTree *T3 = new TTree("T3", "saves important dubug data(charrge)");
    TFile *file = new TFile((pathToTheFile + filename).c_str(),
                            "RECREATE"); //create a file to store the tree as an output

    Float_t D_0_pT = -999, Jet_pT = -999, rapidity = -999, z_val = -999; //variables to store data and used in tree
    Float_t l11 = -1000, l105 = -1000, l115 = -1000, l12 = -1000, l13 = -1000, l20 = -1000;
    Float_t isCharged = -1000, deltaR = -1000, pT_frac = -1000, pl11 = -1000, pl105 = -1000, pl115 = -1000, pl12 = -1000, pl13 = -1000, pl20 = -1000;

    //set tree branches
    T->Branch("l11", &l11, "l11");
    T->Branch("l105", &l105, "l105");
    T->Branch("l115", &l115, "l115");
    T->Branch("l12", &l12, "l12");
    T->Branch("l13", &l13, "l13");
    T->Branch("l20", &l20, "l20");
    T->Branch("pl11", &pl11, "pl11");
    T->Branch("pl105", &pl105, "pl105");
    T->Branch("pl115", &pl115, "pl115");
    T->Branch("pl12", &pl12, "pl12");
    T->Branch("pl13", &pl13, "pl13");
    T->Branch("pl20", &pl20, "pl20");


    T->Branch("z_val", &z_val, "z_val");
    T->Branch("D_0_pT", &D_0_pT, "D_0_pT");
    T->Branch("jet_pT", &Jet_pT, "jet_pT");
    T->Branch("eta", &rapidity, "eta");

    T3->Branch("isCharged", &isCharged, "isCharged");
    T2->Branch("deltaR", &deltaR, "deltaR");
    T2->Branch("pT_frac", &pT_frac, "deltaR");
    for (int i = 0; i < numThreads; i++){
        alocatedThreads[i] = new std::thread(mainThreadedSec, numThreads,
                                             std::to_string(seed),
                                             std::ref(T),
                                             std::ref(D_0_pT),
                                             std::ref(Jet_pT),
                                             std::ref(rapidity),
                                             std::ref(z_val),
                                             std::ref(requiredNumberOfD_0InOneFile),
                                             std::ref(foundNumberOfD_0),
                                             std::ref(l11),
                                             std::ref(l105),
                                             std::ref(l115),
                                             std::ref(l12),
                                             std::ref(l13),
                                             std::ref(l20),
                                             std::ref(pl11),
                                             std::ref(pl105),
                                             std::ref(pl115),
                                             std::ref(pl12),
                                             std::ref(pl13),
                                             std::ref(pl20),
                                             std::ref(T2),
                                             std::ref(T3),
                                             std::ref(isCharged),
                                             std::ref(deltaR),
                                             std::ref(pT_frac));
    }
    for (int i = 0; i < numThreads; i++) {
        alocatedThreads[i]->join();
        std::cout <<"Thread " << i + 1 << " joined" << std::endl;
    }

    T->Print(); //prints the tree structure
    T2->Print();
    T3->Print();
    T->Write();
    T2->Write();
    T3->Write();
    delete file;
    delete T;
}
