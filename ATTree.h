//
// Created by Nikol on 12/16/2024.
//

#ifndef PYTHIA_PROJECT_TREE_CREATOR_ATTREE_H
#define PYTHIA_PROJECT_TREE_CREATOR_ATTREE_H

#include "Pythia8/Pythia.h"
#include "TVector2.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "TCanvas.h"
#include "TH2F.h"
#include "TTree.h"



class ATTree {
    TTree *angT;
    Double_t  D_0_pT = -1000, Jet_pT = -1000, rapidity = -1000, z_val = -1000, l11 = -1000, l105 = -1000, l115 = -1000,
    l12 = -1000, l13 = -1000, l20 = -1000;
    std::string seed;
    bool isSeedSet = false;

    Pythia8::Pythia pythia;

    long unsigned int requiredNumberOfD_0 = 100000;

    Double_t delta_R(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
    void showProgressBar(int progress, int total);

public:
    ATTree();
    ~ATTree();

    void genSeed();
    void runEvents();
    void saveTree();
};


#endif //PYTHIA_PROJECT_TREE_CREATOR_ATTREE_H
