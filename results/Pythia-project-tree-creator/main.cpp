#include <list>
#include <thread>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

#include "TCanvas.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"

#include "functions.h"

//!pTmin_jet - minimum pT for jets is not used


  

int main() {
//    learnD_0JetsOrigin(1);
    learnD_0JetsOriginTest(10);


//    drawBlocksTest();

//    unsigned int requiredNumberOfD_0 = 200;
//    unsigned int foundNumberOfD_0 = 0; //to store number of D_0 particles found
//    unsigned int numThreads = std::thread::hardware_concurrency();
//    int seed = std::time(0) % (900000000 - numThreads);
//    std::thread **alocatedThreads = new std::thread *[numThreads];
//    std::string pathToTheFile = "../results/";
//    std::string filename = "jet tree, cut [D_0 n =" + std::to_string(requiredNumberOfD_0) + ", 1+ GeV].root";
//    TTree *T = new TTree("T", "saves Pt, pseudo rapidity and phi of an D_0 jet"); //create a tree to store the data
//    TFile *file = new TFile( (pathToTheFile + filename).c_str(), "RECREATE"); //create a file to store the tree as an output
//    Float_t D_0_pT = -999, Jet_pT = -999, rapidity = -999, z_val = -999; //variables to store data and used in tree
//    Float_t l11 = -1000, l105 = -1000, l115 = -1000, l12 = -1000, l13 = -1000, l20 = -1000;
//
//    //set tree branches
//    T->Branch("l11", &l11, "l11");
//    T->Branch("l105", &l105, "l105");
//    T->Branch("l115", &l115, "l115");
//    T->Branch("l12", &l12, "l12");
//    T->Branch("l13", &l13, "l13");
//    T->Branch("l20", &l20, "l20");
//    T->Branch("z_val", &z_val, "z_val");
//    T->Branch("D_0_pT", &D_0_pT, "D_0_pT");
//    T->Branch("jet_pT", &Jet_pT, "jet_pT");
//    T->Branch("eta", &rapidity, "eta");
//
//    for (int i = 0; i < numThreads; i++)
//        alocatedThreads[i] = new std::thread(mainSec, numThreads, std::to_string(seed + i), std::ref(T), std::ref(D_0_pT), std::ref(Jet_pT), std::ref(rapidity),
//                                             std::ref(z_val), std::ref(requiredNumberOfD_0), std::ref(foundNumberOfD_0),
//                                             std::ref(l11), std::ref(l105), std::ref(l115), std::ref(l12), std::ref(l13), std::ref(l20));
//
//    for (int i = 0; i < numThreads; i++) {
//        alocatedThreads[i]->join();
//        std::cout << std::endl << "Thread " << i + 1 << " joined" << std::endl;
//    }
//
//    T->Print(); //prints the tree structure
//    T->Write(); //writes the tree to the file
//    file->Close();
//    delete file;


    return 0;
}
