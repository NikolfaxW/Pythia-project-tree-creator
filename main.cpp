#include <list>
#include <thread>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

#include "Pythia8/Pythia.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

#include "TFile.h"

#include "functions.h"
#include "ATTree.h"


//!pTmin_jet - minimum pT for jets is not used


int main() {
    while(getStatus() == "true"){
        ATTree tree;
        tree.runEvents();
        tree.saveTree();
    }

    return 0;
}
