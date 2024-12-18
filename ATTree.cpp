//
// Created by Nikol on 12/16/2024.
//

#include "ATTree.h"

#include <random>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <deque>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

#include "particleUnit.h"
#include "MyInfo.h"


#include "Pythia8/Pythia.h"
#include "TVector2.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "TH2F.h"

ATTree::ATTree(){
    angT = new TTree("angT", "saves z_value, pseudo rapidity and phi, jet_pT, D_0_pT and angularities");

    angT->Branch("z_val", &z_val, "z_val/D");
    angT->Branch("D_0_pT", &D_0_pT, "D_0_pT/D");
    angT->Branch("jet_pT", &Jet_pT, "jet_pT/D");
    angT->Branch("eta", &rapidity, "eta/D");

    angT->Branch("l11", &l11, "l11/D");
    angT->Branch("l105", &l105, "l105/D");
    angT->Branch("l115", &l115, "l115/D");
    angT->Branch("l12", &l12, "l12/D");
    angT->Branch("l13", &l13, "l13/D");
    angT->Branch("l20", &l20, "l20/D");
    n = 0;
}

ATTree::~ATTree(){
    delete angT;
}

long long ATTree::getN() {
    return n;
}

Double_t ATTree::delta_R(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
    Double_t delta_eta = eta1 - eta2;
    Double_t delta_phi = TVector2::Phi_mpi_pi(phi1 - phi2);
//    std::cout << "delta_eta: " << delta_eta << " delta_phi: " << delta_phi << std::endl;
    return TMath::Sqrt(pow(delta_eta, 2) + pow(delta_phi, 2));
}

void ATTree::showProgressBar(int progress, int total) {
    double ratio = static_cast<double>(progress) / total;
    std::cout << "(required number of D_0 found: " << progress << " | " << total << ") ";
    std::cout << int(ratio * 100.0) << "%\r";
    std::cout.flush();
}

void ATTree::genSeed() {
    std::random_device rd; // Use random_device for non-deterministic randomness
    std::mt19937 gen(rd()); // Mersenne Twister random number engine

    // Define normal distribution with mean and standard deviation
    double mean = 450000;
    double stddev = 259807.62;
    std::normal_distribution<> dist(mean, stddev);
    seed = std::to_string( static_cast<int>( abs((std::round(dist(gen))))));
    isSeedSet = true;
}

void ATTree::runEvents() {
    if (!isSeedSet) {
        genSeed();
    }
    pythia.readFile("../config1.cmnd"); //read config file and intialize pythia
    pythia.readString("Random:setSeed = on");  // Enable setting of the seed
    pythia.readString("Random:seed = " + seed);
    pythia.init();

    long unsigned int numberOfD_0Found = 0;

    std::list<particleUnit> j_constituents;
    std::map<TString, fastjet::JetDefinition> jetDefs; //map to store jet definitions

    //parameters for jet finding and to make program less hardcoded
    double R = 0.4; //jet radius
    double pTmin_jet = 5, pThadron = 0.2; //minimum pT for jets and hadrons
    int triggerId = 421; //pdg code of the particle to be found

    double pTMinTrig = 0; //minimum pT for the particle to be found
    double pTMaxTrig = 5.0;
    double mTemp; //This variable are needed to recount momentum after particle mass resets
    double firstEtaCut = 4.2;
    Pythia8::Vec4 pTemp; //This variable are needed to recount momentum after particle mass resets


    Double_t R_frac = -1000, pT_frac = - 1000;
    //define jet finding algorithms here:
    jetDefs["#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
            fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);


    auto &event = pythia.event; //create a reference to the Pythia8 collision event
    std::vector<fastjet::PseudoJet> fjInputs; //to store particles to be clustered
    std::vector<fastjet::PseudoJet> selectedJets; //to store jets after all cuts

    for (int iEvent = 0; numberOfD_0Found < requiredNumberOfD_0; ++iEvent) { //loop over needed number of events
        if (!pythia.next()) {
            n++;
            continue; //generate next event, if it is not possible, continue
        }
        int idxD = -1; // to store index of the D_0 particle in the event
        for (int i = pythia.event.size() - 1; i > 0; i--) { //goes through all particles generated in event
            if (pythia.event[i].idAbs() == triggerId &&
                pythia.event[i].pT() >= pTMinTrig && fabs(event[i].eta()) < firstEtaCut) { //finds D_0 particle with required pT cut
                idxD = i; //saved its index in the event
                break;
            }
        }
        if (idxD == -1) //if there is no D_0 particle in the event, skip it to not waste resources
            continue;

        fjInputs.clear(); //clears the vector of particles to be clustered
        for (int i = 0; i < event.size(); ++i) { // saves particles in order to make jets
            auto &p = event[i];
            if ((event[i].isFinal() || event[i].idAbs() == triggerId) && fabs(event[i].eta()) < firstEtaCut) { //checks if particle is final state or D_0
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
                //! it is better to do eta cut here?
                //eta of constituent is bigger than 1
                fjInputs.push_back(particleTemp); //saves the particle to the vector to be clustered

            }
        }


        if (fjInputs.empty()) { // Abort to avoid running algorithms through empty vectors
            std::cout << "Error: event with no final state particles" << std::endl;
            continue;
        }


        for (auto jetDef: jetDefs) { //for each jet definition runs jet clustering sequence, then saves it to ROOT tree
            selectedJets.clear(); //empties the vector of jets - prevents from saving jets from previous event
            fastjet::ClusterSequence clusterSequence(fjInputs, jetDef.second); //sets up the cluster sequence
            auto jets = sorted_by_pt(clusterSequence.inclusive_jets(0)); //runs the clustering and sorts jets by pT
            fastjet::Selector eta_selector = fastjet::SelectorEtaRange(-1 + R, 1 - R); // Selects jets with |eta| < 1-R
            selectedJets = eta_selector(jets); //applies the eta selector to the jets

            for (const auto &jet: selectedJets) { //loop through all jets
                Bool_t temp_Has_D_0 = false;
                for (const auto &c: jet.constituents()) {//loop through all jet constituents to check if D_0 is there

                    if (c.user_info<MyInfo>().pdg_id() == triggerId) { // Only for D_0
                        temp_Has_D_0 = true;
                        z_val = (jet.px() * c.px() + jet.py() * c.py()) / jet.pt2();
                        D_0_pT = c.pt(); //saves pT of the D_0 particle
                        Jet_pT = jet.pt();
                        rapidity = jet.rapidity();
                        break;
                    }
                }
                if (not temp_Has_D_0) continue; // if there is not d_0 particle in the jet, skip it;
                l11 = 0;
                l105 = 0;
                l115 = 0;
                l12 = 0;
                l13 = 0;
                l20 = 0;

                for (const auto &c: jet.constituents()) { //loop through all jet constituents to calculate l11, l105, l115, l12, l13, l20
                    if(not c.user_info<MyInfo>().isCharged()) continue;
                    pT_frac = c.pt() / jet.pt();
                    R_frac = delta_R(jet.eta(), jet.phi(), c.eta(), c.phi())/R;
                    l11 += pT_frac * R_frac;
                    l105 += pT_frac * pow(R_frac, 0.5);
                    l115 = pT_frac * pow(R_frac, 1.5);
                    l12 += pT_frac * pow(R_frac, 2);
                    l13 += pT_frac * pow(R_frac, 3);
                    l20 += pow(pT_frac, 2);
                }

                angT->Fill();
                ++numberOfD_0Found;
                showProgressBar(numberOfD_0Found, requiredNumberOfD_0);

            }

        }
    }
}

void ATTree::saveTree() {
    std::string path = "../results/" +seed + ".root";
    TFile *file = new TFile(path.c_str(), "RECREATE");
    angT->Print();
    angT->Write();
    file->Close();
    delete file;
}
