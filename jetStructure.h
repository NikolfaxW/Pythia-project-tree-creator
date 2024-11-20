//
// Created by Nikol on 11/20/2024.
//

#ifndef PYTHIA_PROJECT_TREE_CREATOR_JETSTRUCTURE_H
#define PYTHIA_PROJECT_TREE_CREATOR_JETSTRUCTURE_H

#include "Pythia8/Pythia.h"

#include <vector>

struct particle{
    Pythia8::Particle& p;
    double depth;
};

class jetStructure {
    std::vector<particle> particles;

    jetStructure();

    void addParticle(Pythia8::Particle& p, double depth);

};


#endif //PYTHIA_PROJECT_TREE_CREATOR_JETSTRUCTURE_H
