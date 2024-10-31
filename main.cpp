#include <list>
#include <thread>
#include <mutex>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>

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
class MyInfo : public fastjet::PseudoJet::UserInfoBase {
public:
    MyInfo(const int &id, const int &i) : _pdg_id(id), _id(i) {}  //stores ddg codes of the particle

    int pdg_id() const { return _pdg_id; }

    int id() const { return _id; }

protected:
    int _pdg_id;
    int _id;

};

//!pTmin_jet - minimum pT for jets is not used

void mainSec(int numThreads, std::string  seed, TTree *&T, Float_t &D_0_Pt, Float_t &Jet_Pt, Float_t &rapidity, Float_t &hasD_0,
             Float_t &EVI, Float_t &z_val, unsigned int & requiredNumberOfD_0, unsigned int & numberOfD_0Found) {
    Pythia8::Pythia pythia; //create pythia object

    {
        std::lock_guard<std::mutex> lock(std::mutex);  // Lock the mutex
        pythia.readFile("../config1.cmnd"); //read config file and intialize pythia
    }

    pythia.readString("Random:setSeed = on");  // Enable setting of the seed
    pythia.readString("Random:seed = " + seed);
    pythia.init();

    std::list<particleUnit> j_constituents;
    std::map<TString, fastjet::JetDefinition> jetDefs; //map to store jet definitions

    //parameters for jet finding and to make program less hardcoded
    double R = 0.4; //jet radius
    double pTmin_jet = 5, pThadron = 0.2; //minimum pT for jets and hadrons
    int triggerId = 421; //!D_0 //pdg code of the particle to be found
    double pTMinTrig = 1; //minimum pT for the particle to be found
    double pTMaxTrig = 5.0;
    double mTemp; //This variable are needed to recount momentum after particle mass resets
    Pythia8::Vec4 pTemp; //This variable are needed to recount momentum after particle mass resets

    //define jet finding algorithms here:
    jetDefs["#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
            fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);


    auto &event = pythia.event; //create a reference to the Pythia8 collision event
    std::vector<fastjet::PseudoJet> fjInputs; //to store particles to be clustered
    std::vector<fastjet::PseudoJet> selectedJets; //to store jets after all cuts



    for (int iEvent = 0;
         numberOfD_0Found < requiredNumberOfD_0; ++iEvent) { //loop over needed number of events
        if (!pythia.next()) continue; //generate next event, if it is not possible, continue
        int idxD = -1; // to store index of the D_0 particle in the event
        for (int i = pythia.event.size() - 1; i > 0; i--) { //goes through all particles generated in event
            if (pythia.event[i].idAbs() == triggerId &&
                pythia.event[i].pT() >= pTMinTrig ) { //finds D_0 particle with required pT cut
                idxD = i; //saved its index in the event
                break;
            }
        }
        if (idxD == -1) //if there is no D_0 particle in the event, skip it to not waste resources
            continue;

        fjInputs.clear(); //clears the vector of particles to be clustered
        for (int i = 0; i < event.size(); ++i) { // saves particles in order to make jets
            auto &p = event[i];
            if (!event[i].isFinal() && i != idxD) //skips particles that are not final state or not D_0
                continue;
            pTemp = p.p();
            if (p.idAbs() == 22) //changes mass for any other particle except photons to pion mass
                mTemp = 0;
            else
                mTemp = 0.13957;
            fastjet::PseudoJet particleTemp = fastjet::PseudoJet(p.px(), p.py(), p.pz(), sqrt(pTemp.pAbs2() + mTemp *
                                                                                                              mTemp)); //recounts 4 momentum and energy after mass reset
            particleTemp.set_user_info(new MyInfo(p.idAbs(), i)); //adds additional info to the particle
            fjInputs.push_back(particleTemp); //saves the particle to the vector to be clustered
        }


        if (fjInputs.empty()) { // Abort to avoid running algorithms through empty vectors
            std::cout << "Error: event with no final state particles" << std::endl;
            continue;
        }

        Float_t tempD_0_p_t;
        for (auto jetDef: jetDefs) { //for each jet definition runs jet clustering sequence, then saves it to ROOT tree
            selectedJets.clear(); //empties the vector of jets - prevents from saving jets from previous e2vent
            fastjet::ClusterSequence clusterSequence(fjInputs, jetDef.second); //sets up the cluster sequence
            auto jets = sorted_by_pt(clusterSequence.inclusive_jets(0)); //runs the clustering and sorts jets by pT
            fastjet::Selector eta_selector = fastjet::SelectorEtaRange(-1 + R, 1 - R); // Selects jets with |eta| < 1-R
            selectedJets = eta_selector(jets); //applies the eta selector to the jets

            for (const auto &jet: selectedJets) { //loop through all jets
                int tempHas_D_0 = 0;
                for (const auto &c: jet.constituents()) { //loop through all jet constituents to check if D_0 is there
                    if (c.user_info<MyInfo>().pdg_id() == triggerId) {
                        z_val = (jet.px() * c.px() + jet.py() * jet.py()) / jet.pt2();
                        tempD_0_p_t = c.pt(); //saves pT of the D_0 particle
                        tempHas_D_0 = 1;
                        break;
                    }
                }
                if(tempHas_D_0 == 0) continue; // if there is not d_0 particle in the jet, skip it
                hasD_0 = tempHas_D_0;
                {
                    //! Not safe to fill the tree from multiple threads
                    std::lock_guard<std::mutex> lock(std::mutex);  // Lock the mutex
                    D_0_Pt = tempD_0_p_t;
                    EVI = iEvent;
                    Jet_Pt = jet.pt();
                    rapidity = jet.rapidity();
                    T->Fill();  // Fill the data vector safely
                    ++numberOfD_0Found;
                    if(numberOfD_0Found % 100 == 0) std::cout << numberOfD_0Found << " of D_0 instaces were found,"<< std::endl;
                }
            }

        }
    }


}


int main() {
    unsigned int requiredNumberOfD_0 = 10000000;
    unsigned int foundNumberOfD_0 = 0; //to store number of D_0 particles found
    unsigned int numThreads = std::thread::hardware_concurrency();
    int seed = std::time(0) % (900000000 - numThreads);
//    unsigned int numThreads = 2;
    std::thread **alocatedThreads = new std::thread *[numThreads];
    std::string pathToTheFile = "../results/";
    std::string filename = "jet tree, cut:[D_0 n =" + std::to_string(requiredNumberOfD_0) + ", 1+ GeV].root";
    TTree *T = new TTree("T", "saves Pt, pseudo rapidity and phi of an D_0 jet"); //create a tree to store the data
    TFile *file = new TFile( (pathToTheFile + filename).c_str(), "RECREATE"); //create a file to store the tree as an output

    Float_t D_0_Pt, Jet_Pt, rapidity, hasD_0, EVI, z_val; //variables to store data and used in tree

    T->Branch("EventID", &EVI, "EventID"); //set tree branches
    T->Branch("z_val", &z_val, "z_val");
    T->Branch("D_0_Pt", &D_0_Pt, "D_0_Pt");
    T->Branch("pT", &Jet_Pt, "pT");
    T->Branch("eta", &rapidity, "eta");
    T->Branch("D_0", &hasD_0, "D_0");

    for (int i = 0; i < numThreads; i++) {
        alocatedThreads[i] = new std::thread(mainSec, numThreads, std::to_string(seed + i), std::ref(T), std::ref(D_0_Pt), std::ref(Jet_Pt), std::ref(rapidity), std::ref(hasD_0), std::ref(EVI),
                                             std::ref(z_val), std::ref(requiredNumberOfD_0), std::ref(foundNumberOfD_0));
    }
    for (int i = 0; i < numThreads; i++) {
        alocatedThreads[i]->join();
        std::cout << "Thread " << i + 1 << " joined" << std::endl;
    }

    T->Print(); //prints the tree structure
    T->Write(); //writes the tree to the file
    file->Close();
    delete file;


    return 0;
}
