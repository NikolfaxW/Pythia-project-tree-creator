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


//!pTmin_jet - minimum pT for jets is not used


int main() {
    bool DEBUG = false;
    int id;
    if (DEBUG)
        mainSec(0);
    else
        while (true) {
            if (getStatus().compare("false"))
                break;
            while (true) {
                id = getId();
                if (id != -1)
                    break;
            }
            while (true) {
                if (increaseIdOrChageStatus(id, "false"))
                    break;
            }
            mainSec(id);
        }


    return 0;
}
