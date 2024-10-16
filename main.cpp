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


//simple class to store additional info for fastjet
class MyInfo: public fastjet::PseudoJet::UserInfoBase { //here also was charge attribute

public:
    MyInfo(const int & id, const int & i) : _pdg_id(id), _id(i){}

    int pdg_id() 	const {return _pdg_id;}
    int id()		const {return _id;}

protected:
    int _pdg_id;
    int _id;

};

int main() {
    Pythia8::Pythia pythia;
    pythia.readFile("../config1.cmnd");
    pythia.init();

    TFile *file = new TFile("../results/Jet_tree2.root", "RECREATE");
    TTree *T = new TTree("T", "saves Pt, pseudo rapidity and phi of an D_0 jet");
    Float_t pT, gamma, phi, /*j_iter = 0*/ hasD_0, EVI;
    std::list<particleUnit> j_constituents;
//    T->Branch("Jet id", &j_iter, "Jet id");//!
    T->Branch("EventID", &EVI, "EventID");
    T->Branch("pT", &pT, "pT");
    T->Branch("eta", &gamma, "eta");
    T->Branch("phi", &phi, "phi");
    T->Branch("D_0", &hasD_0, "D_0");
//    T->Branch("jet constituents", &j_constituents, "jet constituents");


    std::map<TString, fastjet::JetDefinition> jetDefs;


    double R = 0.4;
    double pTmin_jet = 5, pThadron = 0.2;
    int triggerId = 421; //!D_0
    double pTMinTrig = 0.0;
    double mTemp; //This variables are needed to recount momentum after particle mass resets
    Pythia8::Vec4 pTemp;


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
    std::vector<fastjet::PseudoJet> fjInputs;
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
//            std::cout << "D_0 was found in the event number : " << iEvent << std::endl;
//        }


        fjInputs.clear();
        for (int i = 0; i < event.size(); ++i) { // saves particles in order to make jets
            auto & p = event[i];
            if(not event[i].isFinal() &&  i != idxD)
                continue;
            pTemp = p.p();
            if(p.idAbs() == 22) //changes mass for any other particle except photons to pion mass
                mTemp = 0;
            else
                mTemp = 0.13957;

//            fastjet::PseudoJet particleTemp = pythia.event[i] //!Does not work
            fastjet::PseudoJet particleTemp = fastjet::PseudoJet(p.px(), p.py(), p.pz(), sqrt(pTemp.pAbs2() + mTemp*mTemp));  //comps energy with new mass
            //idk if it is needed to correct boost or no. Let's consider 'no'.
            particleTemp.set_user_info(new MyInfo(p.idAbs(), i));
            fjInputs.push_back(particleTemp);


//            int pdg_id_abs = std::abs(p.id());
//            stable_particles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
//            particles_histogram.push_back(p);
        }


        if (fjInputs.empty()) { // Abort to avoid running algs through empty vecs
            std::cout << "Error: event with no final state particles" << std::endl;
            continue;
        }


        for (auto jetDef: jetDefs) { //for each jet definition runs jet clustering sequance, then saves it to ROOT tree
            selectedJets.clear();
            fastjet::ClusterSequence clustSeq(fjInputs, jetDef.second); //changed from  stable_particles
            auto jets = sorted_by_pt(clustSeq.inclusive_jets(0)); //!! For some reason sorting requires much lower ptmin than pTmin_jet
            fastjet::Selector eta_selector = fastjet::SelectorEtaRange(-1 + R, 1 - R); // Selects jets with |eta| < 1-R
            selectedJets = eta_selector(jets);


//            for (auto jet: selectedJets) { // writes constituents of jets
//                if (jet.pt() < pTmin_jet) continue; // skipping selectedJets with low rapidity
//                j_constituents.clear();
//                // For each particle:
////                ++j_iter;
//                for (auto c: jet.constituents()) {
//                    if (c.pt() < pThadron) continue; //gets rid of bubbles in selectedJets
//                    j_constituents.push_back(particleUnit(c.pt(), c.rapidity(), c.phi()));
//                }
//                pT = jet.pt();
//                gamma = jet.rapidity();
//                phi = jet.phi();
//                T->Fill();
//            }
            for(const auto& jet: selectedJets){ //!EMPTY

//                if (jet.pt() < pTmin_jet) continue; // skipping selectedJets with low rapidity //SAME PROBLEM HERE
                hasD_0 = 0;
                for(const auto& c : jet.constituents() ){
                    if(c.user_info<MyInfo>().pdg_id() == triggerId){
                        hasD_0 = 1;
                        break;
                    }
                }
                if(hasD_0 == 0) continue;
                EVI = iEvent;
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
