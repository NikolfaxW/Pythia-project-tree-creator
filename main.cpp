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
    std::string status;
    while(true){
        status = getStatus();
        if(status.compare("true") != 0)
            break;
        ATTree tree;
        tree.runEvents();
        increaseIdOrChageStatus(getId(), tree.getN(), getStatus());
        tree.saveTree();
    }

    return 0;
}
