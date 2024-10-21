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
class MyInfo: public fastjet::PseudoJet::UserInfoBase {
public:
    MyInfo(const int & id, const int & i) : _pdg_id(id), _id(i){}  //stores ddg codes of the particle

    int pdg_id() 	const {return _pdg_id;}
    int id()		const {return _id;}

protected:
    int _pdg_id;
    int _id;

};

//!pTmin_jet - minimum pT for jets is not used

int main() {
    Pythia8::Pythia pythia; //create pythia object
    pythia.readFile("../config1.cmnd"); //read config file and intialize pythia
    pythia.init();

    TFile *file = new TFile("../results/Jet_tree2.root", "RECREATE"); //create a file to store the tree as an output
    TTree *T = new TTree("T", "saves Pt, pseudo rapidity and phi of an D_0 jet"); //create a tree to store the data
    Float_t pT, rapidity, phi, hasD_0, EVI; //variables to store data and used in tree
    std::list<particleUnit> j_constituents;
    T->Branch("EventID", &EVI, "EventID"); //set tree branches
    T->Branch("pT", &pT, "pT");
    T->Branch("eta", &rapidity, "eta");
    T->Branch("phi", &phi, "phi");
    T->Branch("D_0", &hasD_0, "D_0");

    std::map<TString, fastjet::JetDefinition> jetDefs; //map to store jet definitions

    //parameters for jet finding and to make program less hardcoded
    double R = 0.4; //jet radius
    double pTmin_jet = 5, pThadron = 0.2; //minimum pT for jets and hadrons
    int triggerId = 421; //!D_0 //pdg code of the particle to be found
    double pTMinTrig = 0.0; //minimum pT for the particle to be found
    double mTemp; //This variable are needed to recount momentum after particle mass resets
    Pythia8::Vec4 pTemp; //This variable are needed to recount momentum after particle mass resets

    //define jet finding algorithms here:
    jetDefs["#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
            fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);


    auto &event = pythia.event; //create a reference to the Pythia8 collision event
    std::vector<fastjet::PseudoJet> fjInputs; //to store particles to be clustered
    std::vector<fastjet::PseudoJet> selectedJets; //to store jets after all cuts


    for (int iEvent = 0; iEvent < pythia.mode("Main:numberOfEvents"); ++iEvent) { //loop over needed number of events
        if (!pythia.next()) continue; //generate next event, if it is not possible, continue
        int idxD = -1; // to store index of the D_0 particle in the event
        for (int i = pythia.event.size() - 1; i > 0; i--) { //goes through all particles generated in event
            if (pythia.event[i].idAbs() == triggerId && pythia.event[i].pT() >= pTMinTrig) { //finds D_0 particle with required pT cut
                idxD = i; //saved its index in the event
                break;
            }
        }
        if(idxD == -1) //if there is no D_0 particle in the event, skip it to not waste resources
            continue;

        fjInputs.clear(); //clears the vector of particles to be clustered
        for (int i = 0; i < event.size(); ++i) { // saves particles in order to make jets
            auto & p = event[i];
            if(!event[i].isFinal() &&  i != idxD) //skips particles that are not final state or not D_0
                continue;
            pTemp = p.p();
            if(p.idAbs() == 22) //changes mass for any other particle except photons to pion mass
                mTemp = 0;
            else
                mTemp = 0.13957;
            fastjet::PseudoJet particleTemp = fastjet::PseudoJet(p.px(), p.py(), p.pz(), sqrt(pTemp.pAbs2() + mTemp*mTemp)); //recounts 4 momentum and energy after mass reset
            particleTemp.set_user_info(new MyInfo(p.idAbs(), i)); //adds additional info to the particle
            fjInputs.push_back(particleTemp); //saves the particle to the vector to be clustered
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

            for(const auto& jet: selectedJets){ //loop through all jets
                hasD_0 = 0;
                for(const auto& c : jet.constituents() ){ //loop through all jet constituents to check if D_0 is there
                    if(c.user_info<MyInfo>().pdg_id() == triggerId){
                        hasD_0 = 1;
                        break;
                    }
                }

                if(hasD_0 == 0) continue; // if there is not d_0 particle in the jet, skip it
                EVI = iEvent;
                pT = jet.pt();
                rapidity = jet.rapidity();
                phi = jet.phi();
                T->Fill(); //saves jet parameters to the tree as a new entry
            }

        }
    } //move it to the end in order to split events
    T->Print(); //prints the tree structure
    T->Write(); //writes the tree to the file
    file->Close();
    delete file;

    return 0;
}
