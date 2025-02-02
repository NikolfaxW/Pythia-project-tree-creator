//
// Created by nikol on 11/14/2024.
//

#ifndef PYTHIA_PROJECT_TREE_CREATOR_FUNCTIONS_H
#define PYTHIA_PROJECT_TREE_CREATOR_FUNCTIONS_H

#include "Rtypes.h"
#include "TTree.h"



void showProgressBar(int progress, int total);
Double_t delta_R(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
long long int getId(); //Ok
std::string getStatus(); //Ok
bool increaseIdOrChageStatus(long long int n0, long long int inc, std::string status); //ok status can be "true" or "false" !!!
int randomSeed();
void mainThreadedSec(const int numThreads,
                std::string seed,
                TTree *&T,
                Float_t &D_0_pT,
                Float_t &Jet_Pt,
                Float_t &rapidity,
                Float_t &z_val,
                const unsigned int &requiredNumberOfD_0,
                unsigned int &numberOfD_0Found,
                Float_t &l11,
                Float_t &l105,
                Float_t &l115,
                Float_t &l12,
                Float_t &l13,
                Float_t &l20,
                Float_t &pl11,
                Float_t &pl105,
                Float_t &pl115,
                Float_t &pl12,
                Float_t &pl13,
                Float_t &pl20,
                TTree *&T2,
                TTree *&T3,
                Float_t &isCharged,
                Float_t &deltaR,
                Float_t &pT_frac);
void mainSec(int id);



#endif //PYTHIA_PROJECT_TREE_CREATOR_FUNCTIONS_H
