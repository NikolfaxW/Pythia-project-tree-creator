//
// Created by Nikol on 11/20/2024.
//

#include "jetStructure.h"

jetStructure::jetStructure() { //TODO
    this->particles.reserve(100);
}

void jetStructure::addParticle(Pythia8::Particle &p, double depth) { //TODO
        this->particles.push_back({p, depth});
}

