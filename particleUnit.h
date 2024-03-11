//
// Created by nikol on 3/11/2024.
//

#ifndef PYTHIA_PROJECT_TREE_CREATOR_PARTICLEUNIT_H
#define PYTHIA_PROJECT_TREE_CREATOR_PARTICLEUNIT_H

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


class particleUnit {
    Float_t pT, gamma, phi;

public:
    particleUnit(Float_t pT, Float_t gamma, Float_t phi) : pT(pT), gamma(gamma), phi(phi){}

};


#endif //PYTHIA_PROJECT_TREE_CREATOR_PARTICLEUNIT_H
