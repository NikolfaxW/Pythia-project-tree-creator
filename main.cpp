#include <list>
#include <thread>
#include <mutex>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
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
#include "TVector2.h"
#include "TMath.h"

#include "particleUnit.h"


//simple class to store additional info for fastjet
class MyInfo : public fastjet::PseudoJet::UserInfoBase {
public:
    MyInfo(const int &id, const int &i, const bool & is_charged) : _pdg_id(id), _id(i), _is_charged(is_charged) {}  //stores ddg codes of the particle

    int pdg_id() const { return _pdg_id; }

    int id() const { return _id; }

    bool isCharged() const { return _is_charged;}

protected:
    int _pdg_id;
    int _id;
    bool _is_charged;


};

//!pTmin_jet - minimum pT for jets is not used

void showProgressBar(int progress, int total) {
    float ratio = static_cast<float>(progress) / total;
    std::cout << "(required number of D_0 found: " << progress << " | " << total << ") ";
    std::cout << int(ratio * 100.0) << "%\r";
    std::cout.flush();
}

Float_t delta_R(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2) {
    Float_t delta_eta= eta1 - eta2;
    Float_t delta_phi = TVector2::Phi_mpi_pi(phi1 - phi2);
//    std::cout << "delta_eta: " << delta_eta << " delta_phi: " << delta_phi << std::endl;
    return TMath::Sqrt(pow(delta_eta, 2) + pow(delta_phi, 2));
}

void learnD_0JetsOrigin(const unsigned int requiredNumberOfD_0){
    particleDictionarry dict;
    Pythia8::Pythia pythia; //create pythia object
    pythia.readFile("../config1.cmnd"); //read config file and intialize pythia
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
    unsigned int numberOfD_0Found = 0;

    std::ofstream file("../results/D_0-origins-structure.txt");
    if (!file.is_open()){
        std::cout <<  "Error opening file for writing." << std::endl;
        return;
    }
    bool to_print = true;
    for (int iEvent = 0; numberOfD_0Found < requiredNumberOfD_0; ++iEvent) { //loop over needed number of events
        if (!pythia.next()) continue; //generate next event, if it is not possible, continue
        int idxD = -1; // to store index of the D_0 particle in the event
        int i = 0;
        for (i = pythia.event.size() - 1; i > 0; i--) { //goes through all particles generated in event
            if (pythia.event[i].idAbs() == triggerId &&
                pythia.event[i].pT() >= pTMinTrig) { //finds D_0 particle with required pT cut
                idxD = i; //saved its index in the event
                break;
            }
        }
        if (idxD == -1) //if there is no D_0 particle in the event, skip it to not waste resources
        {
            continue;
        }
        ++numberOfD_0Found;
        auto p = pythia.event[i];
        file << "[" << dict.getIdABSName(p.idAbs())  << " ; " << dict.getStatusName(p.status()) << " ; " << i << " ]  =>  ";
        while(i != 0){
//            to_print = p.mother1() != p.mother2();
            i = p.mother1();
            p = pythia.event[i];
            if(to_print) file << "[" << dict.getIdABSName(p.idAbs())  << " ; " << dict.getStatusName(abs(p.status())) << " ; " << i << " ]  =>  ";

        }
        file << "x\n";
    }
    file.close();

}


void mainSec(const int numThreads, std::string  seed, TTree *&T, Float_t &D_0_pT, Float_t &Jet_Pt, Float_t &rapidity,
             Float_t &z_val, const unsigned int & requiredNumberOfD_0, unsigned int & numberOfD_0Found,
             Float_t &l11, Float_t &l105, Float_t &l115, Float_t &l12, Float_t &l13, Float_t &l20, bool &hasUncharged, bool R_frac_mistake, bool &pT_frac_mistake){
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

    Float_t temp_z_val, temp_D_0_pT, temp_l11, temp_l105, temp_l115, temp_l12, temp_l13, temp_l20, R_frac, pT_frac;
    bool temp_hasUncharged, temp_R_frac_mistake, temp_pT_frac_mistake;
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
            if (event[i].isFinal() || event[i].idAbs() == triggerId) { //checks if particle is final state or D_0
                pTemp = p.p();
                if (p.idAbs() == 22) //changes mass for any other particle except photons to pion mass
                    mTemp = 0;
                else
                    mTemp = 0.13957;
                fastjet::PseudoJet particleTemp = fastjet::PseudoJet(p.px(), p.py(), p.pz(), sqrt(pTemp.pAbs2() + mTemp *
                                                                                                                  mTemp)); //recounts 4 momentum and energy after mass reset
                particleTemp.set_user_info(new MyInfo(p.idAbs(), i, event[i].isCharged())); //adds additional info to the particle
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
                        temp_z_val = (jet.px() * c.px() + jet.py() * jet.py()) / jet.pt2();
                        temp_D_0_pT = c.pt(); //saves pT of the D_0 particle
                        break;
                    }
                }
                if(not temp_Has_D_0) continue; // if there is not d_0 particle in the jet, skip it;
                temp_l11 = 0;
                temp_l105 = 0;
                temp_l115 = 0;
                temp_l12 = 0;
                temp_l13 = 0;
                temp_l20 = 0;
                temp_hasUncharged = false;
                temp_R_frac_mistake = false;
                temp_pT_frac_mistake = false;

                for (const auto &c: jet.constituents()) { //loop through all jet constituents to calculate l11, l105, l115, l12, l13, l20

                    if(c.user_info<MyInfo>().pdg_id() != triggerId && !c.user_info<MyInfo>().isCharged()) continue; //skip not D_0 and neutral particles
                    pT_frac = c.pt() / jet.pt();
                    R_frac = delta_R(jet.eta(), jet.phi(), c.eta(), c.phi()) / R;
                    if(R_frac > 1) { temp_R_frac_mistake = true; continue; }
                    if(pT_frac> 1) temp_pT_frac_mistake = true;

                    temp_l11 += pT_frac*R_frac;
                    temp_l105 += pT_frac * pow(R_frac, 0.5);
                    temp_l115 = pT_frac * pow(R_frac, 1.5);
                    temp_l12 += pT_frac * pow(R_frac, 2);
                    temp_l13 += pT_frac * pow(R_frac, 3);
                    temp_l20 += pow(pT_frac,2);
                    if(c.user_info<MyInfo>().isCharged() && c.user_info<MyInfo>().pdg_id() != triggerId) temp_hasUncharged = true;

                }


                {
                    //! Not safe to fill the tree from multiple threads
                    std::lock_guard<std::mutex> lock(std::mutex);  // Lock the mutex
                    hasUncharged = temp_hasUncharged;
                    l11 = temp_l11;
                    l105 = temp_l105;
                    l115 = temp_l115;
                    l12 = temp_l12;
                    l13 = temp_l13;
                    l20 = temp_l20;
                    D_0_pT = temp_D_0_pT;
                    z_val = temp_z_val;
                    Jet_Pt = jet.pt();
                    rapidity = jet.rapidity();
                    R_frac_mistake = temp_R_frac_mistake;
                    pT_frac_mistake = temp_pT_frac_mistake;
                    T->Fill();  // Fill the data vector safely
                    ++numberOfD_0Found;
                    if(numberOfD_0Found % 10 == 0) showProgressBar(numberOfD_0Found, requiredNumberOfD_0);
                }
            }

        }
    }


}
  

int main() {
//    learnD_0JetsOrigin(100);
//
    unsigned int requiredNumberOfD_0 = 200;
    unsigned int foundNumberOfD_0 = 0; //to store number of D_0 particles found
    unsigned int numThreads = std::thread::hardware_concurrency();
    int seed = std::time(0) % (900000000 - numThreads);
    std::thread **alocatedThreads = new std::thread *[numThreads];
    std::string pathToTheFile = "../results/";
    std::string filename = "jet tree, cut [D_0 n =" + std::to_string(requiredNumberOfD_0) + ", 1+ GeV].root";
    TTree *T = new TTree("T", "saves Pt, pseudo rapidity and phi of an D_0 jet"); //create a tree to store the data
    TFile *file = new TFile( (pathToTheFile + filename).c_str(), "RECREATE"); //create a file to store the tree as an output
    Float_t D_0_pT = -999, Jet_pT = -999, rapidity = -999, z_val = -999; //variables to store data and used in tree
    Float_t l11, l105, l115, l12, l13, l20;
    bool hasUncharged, R_frac_mistake, pT_frac_mistake;
    //set tree branches


    T->Branch("R_frac_mistake", &R_frac_mistake, "R_frac_mistake");
    T->Branch("pT_frac_mistake", &pT_frac_mistake, "pT_frac_mistake");
    T->Branch("hasUncharged", &hasUncharged, "hasUncharged");
    T->Branch("l11", &l11, "l11");
    T->Branch("l105", &l105, "l105");
    T->Branch("l115", &l115, "l115");
    T->Branch("l12", &l12, "l12");
    T->Branch("l13", &l13, "l13");
    T->Branch("l20", &l20, "l20");
    T->Branch("z_val", &z_val, "z_val");
    T->Branch("D_0_pT", &D_0_pT, "D_0_pT");
    T->Branch("jet_pT", &Jet_pT, "jet_pT");
    T->Branch("eta", &rapidity, "eta");

    for (int i = 0; i < numThreads; i++)
        alocatedThreads[i] = new std::thread(mainSec, numThreads, std::to_string(seed + i), std::ref(T), std::ref(D_0_pT), std::ref(Jet_pT), std::ref(rapidity),
                                             std::ref(z_val), std::ref(requiredNumberOfD_0), std::ref(foundNumberOfD_0),
                                             std::ref(l11), std::ref(l105), std::ref(l115), std::ref(l12), std::ref(l13), std::ref(l20),
                                             std::ref(hasUncharged), std::ref(R_frac_mistake), std::ref(pT_frac_mistake));

    for (int i = 0; i < numThreads; i++) {
        alocatedThreads[i]->join();
        std::cout << std::endl << "Thread " << i + 1 << " joined" << std::endl;
    }

    T->Print(); //prints the tree structure
    T->Write(); //writes the tree to the file
    file->Close();
    delete file;


    return 0;
}
