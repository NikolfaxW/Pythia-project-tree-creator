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

#include "Pythia8/Pythia.h"
#include "TVector2.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

std::string createBlock( const unsigned int x, const unsigned int y, std::string fill, const unsigned int id, const unsigned int status, const unsigned int evi){
    std::string result = "<rect x = \"" + std::to_string(x) + "\" y = \"" + std::to_string(y) + "\" width = \"150\" height = \"100\" fill = \"" + fill + "\" stroke = \"black\" stroke-width = \"2\"></rect>\n";
    result.append("<text x = \"" + std::to_string(x + 75) + "\" y = \"" + std::to_string(y + 60) + "\" font-size = \"16\" text-anchor = \"middle\" fill = \"black\">" + "Id : " + std::to_string(id) + " ; St : " + std::to_string(status) + " ; Ei : " + std::to_string(evi) + "</text>\n");
    return result;
}

std::string createArrow(const unsigned int x1, const unsigned int y1, const unsigned int x2, const unsigned int y2){
    std::string result = "<line x1 = \"" + std::to_string(x1) + "\" y1 = \"" + std::to_string(y1) + "\" x2 = \"" + std::to_string(x2) + "\" y2 = \"" + std::to_string(y2) + "\" stroke = \"black\" stroke-width = \"2\"></line>\n";
    result.append("<polygon points = \"" + std::to_string(x2) + "," + std::to_string(y2 - 5) + " " + std::to_string(x2) + "," + std::to_string(y2 + 5) + " " + std::to_string(x2 + 10) + "," + std::to_string(y2) + "\" fill = \"black\"></polygon>\n");
    return result;
}

void drawBlocksTest() {
    std::remove("../results/block_diagram.html");
    // Create an HTML file
    std::ofstream htmlFile("../results/block_diagram.html");

    if (!htmlFile) {
        std::cerr << "Error creating file." << std::endl;
        return;
    }

    // Write the HTML content
    htmlFile << R"(
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Block Diagram</title>
           <style>
                body {
                    font-family: 'Arial', sans-serif;
                }
            </style>
</head>
<body>
    <svg viewBox="0 0 500 200" style="width: 20%; height: 20%;">)";
//    <!-- Block 1 -->
    htmlFile << createBlock(50,50, "lightblue", 1, 1, 1);
    htmlFile << createBlock(300,50, "lightgreen", 2, 1, 1);
    htmlFile << createArrow(200, 100, 300, 100);

    htmlFile << R"(
        </svg>
</body>
</html>)";
    htmlFile.close();

    std::cout << "HTML file created: block_diagram.html" << std::endl;
}



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


void mainSec(const int numThreads, std::string seed, TTree *&T, Float_t &D_0_pT, Float_t &Jet_Pt, Float_t &rapidity,
             Float_t &z_val, const unsigned int &requiredNumberOfD_0, unsigned int &numberOfD_0Found,
             Float_t &l11, Float_t &l105, Float_t &l115, Float_t &l12, Float_t &l13, Float_t &l20) {
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
    ;
    Double_t R_frac, pT_frac, temp_z_val, temp_D_0_pT, temp_l11, temp_l105, temp_l115, temp_l12, temp_l13, temp_l20;

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
                        temp_z_val = (jet.px() * c.px() + jet.py() * jet.py()) / jet.pt2();
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


                for (const auto &c: jet.constituents()) { //loop through all jet constituents to calculate l11, l105, l115, l12, l13, l20
                    if (c.user_info<MyInfo>().pdg_id() != triggerId && !c.user_info<MyInfo>().isCharged())
                        continue; //skip not D_0 and neutral particles
                    pT_frac = c.pt() / jet.pt();
                    R_frac = delta_R(jet.eta(), jet.phi(), c.eta(), c.phi()) / R;
                    temp_l11 += pT_frac * R_frac;
                    temp_l105 += pT_frac * pow(R_frac, 0.5);
                    temp_l115 = pT_frac * pow(R_frac, 1.5);
                    temp_l12 += pT_frac * pow(R_frac, 2);
                    temp_l13 += pT_frac * pow(R_frac, 3);
                    temp_l20 += pow(pT_frac, 2);
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
                    D_0_pT = temp_D_0_pT;
                    z_val = temp_z_val;
                    Jet_Pt = jet.pt();
                    rapidity = jet.rapidity();
                    T->Fill();  // Fill the data vector safely
                    ++numberOfD_0Found;
                    if (numberOfD_0Found % 10 == 0) showProgressBar(numberOfD_0Found, requiredNumberOfD_0);
                }
            }

        }
    }


}