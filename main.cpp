#include <list>

#include "Pythia8/Pythia.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH2D.h"
#include "TMath.h"
#include "TMarker.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "particleUnit.h"



class MyInfo: public fastjet::UserInfoBase {

public:
    MyInfo(const int & id, const int & i, const int & ch) : _pdg_id(id), _id(i), _ch(ch) {}

    int pdg_id() 	const {return _pdg_id;}
    int id()		const {return _id;}
    int ch()		const {return _ch;}

protected:
    int _pdg_id;
    int _id;
    int _ch;
};


int main() {
    Pythia8::Pythia pythia;
    pythia.readFile("../config1.cmnd");
    pythia.init();

    TFile *file = new TFile("../results/Jet_tree.root", "RECREATE");
    TTree *T = new TTree("T", "saves Pt, pseudo rapidity and phi of an jet");
    Float_t pT, gamma, phi; /*j_iter = 0*/
    std::list<particleUnit> j_constituents;
//    T->Branch("Jet id", &j_iter, "Jet id");//!
    T->Branch("pT", &pT, "pT");
    T->Branch("eta", &gamma, "eta");
    T->Branch("phi", &phi, "phi");
    T->Branch("jet constituents", &j_constituents, "jet constituents");


    std::map<TString, fastjet::JetDefinition> jetDefs;


    double R = 0.4;
    double pTmin_jet = 5, pThadron = 0.2;
    int triggerId = 421; //!D_0
    double pTMinTrig = 0.0;
    double mTemp = 0;

    TString description = "Number of events: " + std::to_string(pythia.mode("Main:numberOfEvents"));

    //define jet finding algorithms here:
    jetDefs["#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
            fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    //CaCambridge-Aachen example
    //        jetDefs["Cambridge-Aachen jets,  #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
    //                fastjet::cambridge_algorithm, R, fastjet::E_scheme, fastjet::Best);

    //till here

    auto &event = pythia.event;
    std::vector<Pythia8::Particle> particles_histogram;
    std::vector<fastjet::PseudoJet> stable_particles;
    std::vector<fastjet::PseudoJet> selectedJets;


    for (int iEvent = 0; iEvent < pythia.mode("Main:numberOfEvents"); ++iEvent) { //choosing final particles only
        particles_histogram.clear();
        stable_particles.clear();
        if (!pythia.next()) continue;

        //coppied from the example
        // Find the index of the D0 meson in an event and (optonally) its pT and eta
        int idxD = -1;
        for (int i = pythia.event.size() - 1; i > 0; i--) {
            if (pythia.event[i].idAbs() == triggerId && pythia.event[i].pT() >= pTMinTrig) {
                idxD = i;
//                pt = pythia.event[i].pT();
//                eta = pythia.event[i].eta();
                break;
            }
        }
        if(idxD == -1) //! works fine
            continue; //! D_0 was not found in the event
//        else{
//            std::cout << "D_0 was found in the event" << std::endl;
//        }



        for (int i = 0; i < event.size(); ++i) { // saves particles in order to make jets
            auto &p = event[i];
            if(not p.isFinal())
                continue;
            if(p.id() == 22) //changes mass for any other particle except photons to pion mass
                mTemp = 0;
            else
                mTemp = 0.13957;
            p.e(sqrt(p.pAbs2() + mTemp * mTemp);
            int pdg_id_abs = std::abs(p.id());
            mTemp = pdg_id_abs == 22 ? 0. : 0.13957;
            stable_particles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
            particles_histogram.push_back(p);
        }




        for (auto jetDef: jetDefs) { //for each jet definition runs jet clustering sequance, then saves it to ROOT tree
            fastjet::ClusterSequence clustSeq(stable_particles, jetDef.second);
            auto jets = sorted_by_pt(clustSeq.inclusive_jets(pTmin_jet));
            fastjet::Selector eta_selector = fastjet::SelectorEtaRange(-1 + R, 1 - R); // Selects jets with |eta| < 1-R
            selectedJets = eta_selector(jets);


            for (auto jet: selectedJets) {
                if (jet.pt() < pTmin_jet) continue; // skipping selectedJets with low rapidity
                j_constituents.clear();
                // For each particle:
//                ++j_iter;
                for (auto c: jet.constituents()) {
                    if (c.pt() < pThadron) continue; //gets rid of bubbles in selectedJets
                    j_constituents.push_back(particleUnit(c.pt(), c.rapidity(), c.phi()));
                }
                pT = jet.pt();
                gamma = jet.rapidity();
                phi = jet.phi();
                T->Fill();
            }
        }
    } //move it to the end in order to split events
    T->Print();
    T->Write();
    file->Close();

    delete file;

    return 0;
}
