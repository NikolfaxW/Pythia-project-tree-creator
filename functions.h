//
// Created by nikol on 11/14/2024.
//

#ifndef PYTHIA_PROJECT_TREE_CREATOR_FUNCTIONS_H
#define PYTHIA_PROJECT_TREE_CREATOR_FUNCTIONS_H

#include "Rtypes.h"
#include "TTree.h"

void drawBlocksTest();

void showProgressBar(int progress, int total);
Double_t delta_R(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
void learnD_0JetsOrigin(const unsigned int requiredNumberOfD_0);
void mainSec(const int numThreads, std::string  seed, TTree *&T, Float_t &D_0_pT, Float_t &Jet_Pt, Float_t &rapidity,
             Float_t &z_val, const unsigned int & requiredNumberOfD_0, unsigned int & numberOfD_0Found,
             Float_t &l11, Float_t &l105, Float_t &l115, Float_t &l12, Float_t &l13, Float_t &l20);



#endif //PYTHIA_PROJECT_TREE_CREATOR_FUNCTIONS_H
