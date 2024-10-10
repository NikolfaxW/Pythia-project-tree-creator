// Barbara Trzeciak, 19 May 2023
// HF jets - PYTHIA generation of jets with Dzero mesons
// STAR at RHIC

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"

// ROOT, for histogramming.
#include "TROOT.h"
#include "TH1.h"
#include "TH2D.h"
#include "TVector2.h"
// ROOT, for saving file.
#include "TFile.h"
#include "TNtuple.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/Selector.hh>
#include <fastjet/SharedPtr.hh>


using namespace std;
using namespace fastjet;
using namespace Pythia8;

class MyInfo: public PseudoJet::UserInfoBase {

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

double delta_R(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    return sqrt(deta*deta + dphi*dphi);
}

int main(int argc, char* argv[]) {


    // Select common parameters for FastJet analysis.
    double Rparam = 0.4;
    fastjet::JetAlgorithm 		algorithm = fastjet::antikt_algorithm;
    fastjet::RecombinationScheme	recombScheme = fastjet::E_scheme;
    fastjet::Strategy 			strategy = fastjet::Best;
    fastjet::JetDefinition		*jetDef = nullptr;

    int triggerId = 421;

    double pTMin = 0.2;			// Min particle pT
    double pTMinTrig = 0.0;		// Min triger particle pT
    double pTMinJet   = 0.0;      // Min jet pT.
    double etaMax  = 1.0;    		// Pseudorapidity range of detector.
    bool 	 chJets = 0;
    bool 	 chConst = 1;
    TString confFileName = TString(argv[1]).Remove(TString(argv[1]).Last('.'));

    // Create a ROOT file to store the output
    TString outputFileName;
    if(chJets) outputFileName = TString::Format("out/pythia8_chargedJets_%s_R0%.f_hfpt%.1f_%s.root",confFileName.Data(), Rparam*10, pTMinTrig, argv[2]);
    else  if(chConst) outputFileName = TString::Format("out/pythia8_inclusiveJets_chargedConst_%s_R0%.f_hfpt%.1f_%s.root",confFileName.Data(), Rparam*10, pTMinTrig, argv[2]);
    else outputFileName = TString::Format("out/pythia8_inclusiveJets_%s_R0%.f_hfpt%.1f_%s.root",confFileName.Data(), Rparam*10, pTMinTrig, argv[2]);
    TFile* file = TFile::Open(outputFileName.Data(), "RECREATE");
    // Create a TNtuple to store the event information
    TNtuple ntuple("HFjets", "HF-jets Information", "pt_jet:eta_jet:pt_hf:eta_hf:y_hf:zvalue:lambda_1:lambda_15:lambda_2:lambda_3:lambda_4");
    // TNtuple *ntuple = new TNtuple("HFjets", "HF-jets Information", "pt_jet:eta_jet:pt_hf:eta_hf:zvalue:lambda_1:lambda_15:lambda_2:lambda_3:lambda_4");
    Long64_t maxFileSize = 1e7;  // 0.01 GB (adjust this value as needed)
    // ntuple->SetMaxTreeSize(maxFileSize);
    ntuple.SetMaxTreeSize(maxFileSize);

    // Generator. Shorthand for event.
    Pythia pythia;
    Event& event = pythia.event;

    // Extract settings to be used in the main program.
    pythia.readFile(argv[1]);
    int nEvent = pythia.mode("Main:numberOfEvents");
    int nAbort = pythia.mode("Main:timesAllowErrors");

    pythia.init();

    // Set up FastJet jet finder.
    jetDef = new fastjet::JetDefinition(algorithm, Rparam,
                                        recombScheme, strategy);
    // Fastjet input
    std::vector <fastjet::PseudoJet> fjInputs;

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

        if (!pythia.next()) continue;

        // Find trigger particle
        int idxD = -1;
        double pt, eta;
        for (int i = pythia.event.size() - 1; i > 0; i--) {
            if (pythia.event[i].idAbs() == triggerId && pythia.event[i].pT() >= pTMinTrig) {
                idxD = i;
                pt = pythia.event[i].pT();
                eta = pythia.event[i].eta();
                break;
            }
        }
        if (idxD == -1) {
            //cout << "Error: Could not find trigger particle " << triggerId << endl;
            continue;
        }
        //else cout << "Found trigger particle " << triggerId << "\t id = " << idxD << "\t in event: " << iEvent << "\t pt = " << pt << "\t eta = " << eta << endl;

        fjInputs.resize(0);
        Vec4   pTemp;
        double mTemp;
        int pdg_id, pdg_id_abs, charge;

        for (int i = 0; i < pythia.event.size(); ++i) {
            // Final state only
            if (!pythia.event[i].isFinal())        continue;
            if (chJets) { if(!pythia.event[i].isCharged() && pythia.event[i].idAbs() != triggerId)	continue; }

            pdg_id = pythia.event[i].id();
            pdg_id_abs = pythia.event[i].idAbs();
            charge = pythia.event[i].isCharged();
            // No neutrinos
            //===================!!!! CHECK FOR SECONDARIES !!!!!!
            if (pdg_id_abs == 12 || pdg_id_abs == 14 || pdg_id_abs == 16) continue;
            // if (geantid == 4 || geantid == 5 || geantid == 6 || geantid == 10 || geantid == 16 || geantid == 17 || geantid == 18 || geantid == 20 || geantid == 22 || geantid == 26 || geantid == 28 || geantid == 30 ) continue;


            // |eta| and pt acceptance
            // ==================!!! SHOULD BE PARTICLE PT CUT BE USED ?
            if (abs(pythia.event[i].eta()) > etaMax || pythia.event[i].pT() < pTMin) continue;

            // Create a PseudoJet from the complete Pythia particle.
            fastjet::PseudoJet particleTemp = pythia.event[i];
            // Modify mass to pi0 or gamma mass
            pTemp = pythia.event[i].p();
            mTemp = pythia.event[i].m();
            // ==================!!! CHECK THIS MASS SETTING, for D meson and other neutrals
            mTemp = pdg_id_abs == 22 ? 0. : 0.13957;
            pTemp.e( sqrt(pTemp.pAbs2() + mTemp*mTemp) );
            particleTemp.reset_momentum( pTemp.px(), pTemp.py(), pTemp.pz(), pTemp.e() );
            particleTemp.set_user_info(new MyInfo(pdg_id, i, charge));

            // Store acceptable particles as input to Fastjet.
            fjInputs.push_back(particleTemp);

        }

        if (fjInputs.size() == 0) {
            cout << "Error: event with no final state particles" << endl;
            continue;
        }

        // Run Fastjet algorithm and sort jets in pT order.
        vector <fastjet::PseudoJet> inclusiveJets, sortedJets, selectedJets;
        fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
        inclusiveJets = clustSeq.inclusive_jets(pTMinJet);
        sortedJets    = sorted_by_pt(inclusiveJets);

        fastjet::Selector eta_selector = fastjet::SelectorEtaRange(-1+Rparam, 1-Rparam); // Selects jets with |eta| < 1-R
        selectedJets = eta_selector(sortedJets);

        //cout << "Number of jets: " << sortedJets.size() << "\t selected jets: " << selectedJets.size() << endl;

        for(int jetid = 0; jetid < selectedJets.size(); jetid++){

            vector <fastjet::PseudoJet> constituents = selectedJets[jetid].constituents();
            vector <fastjet::PseudoJet> sortedconstituents = sorted_by_pt(constituents);

            bool d0Jet = false;
            double HFPx = -99.;
            double HFPy = -99.;
            double HFPt = -99.;
            double HFEta = -99.;
            double HFy = -99.;

            double zvalue = -99.;
            double lambda[5] = {0.,0.,0.,0.,0.};
            double lambda_alpha[5] = {1.,1.5,2.,3.,4.};
            double lambda_kappa = 1;

            double JetPx = selectedJets[jetid].px();
            double JetPy = selectedJets[jetid].py();
            double JetPt = selectedJets[jetid].perp();
            double JetEta = selectedJets[jetid].eta();
            double JetPhi = selectedJets[jetid].phi();

            // Loop over jet constituents to get some info
            for (unsigned j = 0; j < sortedconstituents.size(); j++){
                //cout << "id: " << abs(sortedconstituents[j].user_info<MyInfo>().pdg_id()) << endl;

                double consPx = sortedconstituents[j].px();
                double consPy = sortedconstituents[j].py();
                double consPt = sortedconstituents[j].perp();
                double consEta = sortedconstituents[j].eta();
                double consPhi = sortedconstituents[j].phi();

                double Delta_R = delta_R(JetEta,JetPhi,consEta,consPhi);

                if (chConst){
                    if(sortedconstituents[j].user_info<MyInfo>().ch() != 0) {
                        for(int a = 0; a<5; a++) lambda[a]+=pow(consPt / JetPt, lambda_kappa)*pow( Delta_R / Rparam, lambda_alpha[a]);
                    }
                }
                else {
                    for(int a = 0; a<5; a++) lambda[a]+=pow(consPt / JetPt, lambda_kappa)*pow( Delta_R / Rparam, lambda_alpha[a]);
                }
                // Check if trigger particle
                if (abs(sortedconstituents[j].user_info<MyInfo>().pdg_id()) == triggerId){

                    HFPx = sortedconstituents[j].px();
                    HFPy = sortedconstituents[j].py();
                    HFPt = sortedconstituents[j].perp();
                    HFEta = sortedconstituents[j].eta();
                    HFy = sortedconstituents[j].rap();

                    if (HFPt < pTMinTrig) continue;

                    d0Jet = true;

                    zvalue = (JetPx*HFPx + JetPy*HFPy)/(JetPt * JetPt);

                    //LOG_INFO << Form("Jet found with pT: %f, D0 with pT: %f, z: %f", MCJetPt, MCD0Pt, zvalue) << endl;
                }

                // End of jet constituent loop
            }

            if (d0Jet) ntuple.Fill(JetPt, JetEta, HFPt, HFEta, HFy, zvalue, lambda[0],lambda[1], lambda[2], lambda[3], lambda[4]);
            //if (d0Jet) ntuple -> Fill(JetPt, JetEta, HFPt, HFEta, zvalue, lambda[0],lambda[1], lambda[2], lambda[3], lambda[4]);


            // End of jet loop
        }

        // End of event loop.
    }

    //ntuple->Write();
    ntuple.Write();
    file->Close();

    // Statistics. Histograms.
    pythia.stat();

    //cout <<  nJetsF << pTjetsF << RdistF << distJets << pTdiff << nAna << tGen;

    // Done.
    return 0;

}


