//
// Created by nikol on 3/11/2024.
//

#include "particleUnit.h"

particleDictionarry::particleDictionarry(){
    idABSNames[1] = "d";
    idABSNames[-1] = "-d";
    idABSNames[2] = "u";
    idABSNames[-2] = "-u";
    idABSNames[3] = "s";
    idABSNames[-3] = "-s";
    idABSNames[4] = "c";
    idABSNames[-4] = "-c";
    idABSNames[5] = "b";
    idABSNames[-5] = "-b";
    idABSNames[6] = "t";
    idABSNames[-6] = "-t";
    idABSNames[11] = "e";
    idABSNames[12] = "e neutrino";
    idABSNames[13] = "mu";
    idABSNames[14] = "mu neutrino";
    idABSNames[15] = "tau";
    idABSNames[16] = "tau neutrino";
    idABSNames[21] = "g";
    idABSNames[22] = "photon";
    idABSNames[23] = "Z_0";
    idABSNames[24] = "W+";
    idABSNames[25] = "Higgs";
    idABSNames[90] = "Cluster (intermediate hadronization state)";
    idABSNames[111] = "pi_0";
    idABSNames[211] = "pi+";
    idABSNames[221] = "nu";
    idABSNames[331] = "nu'";
    idABSNames[411] = "D+";
    idABSNames[413] = "D_{s1}(2536)^+";
    idABSNames[421] = "D_0";
    idABSNames[423] = "D*_0";
    idABSNames[431] = "D+_s";
    idABSNames[511] = "B_0";
    idABSNames[521] = "B+";
    idABSNames[531] = "B_0_s";
    idABSNames[541] = "B+_c";
    idABSNames[2112] = "N";
    idABSNames[2212] = "P";

    statusNames[1] = "final-state";
    statusNames[2] = "intermediate/unstable";
    statusNames[3] = "partons";
    statusNames[21] = "parton showering, intermediate";
    statusNames[31] = " Initial-state";
    statusNames[31] = " Initial-state";
    statusNames[11] = "whole event";
    statusNames[12] = "incoming beam";
    statusNames[13] = " incoming beam-inside-beam";
    statusNames[14] = "outgoing elastically scattered";
    statusNames[15] = "outgoing diffractively scattered";
    statusNames[21] = "incoming";
    statusNames[22] = "intermediate (intended to have preserved mass)";
    statusNames[23] = "outgoing";
    statusNames[24] = "outgoing, nonperturbatively kicked out in diffraction";
    statusNames[31] = "incoming, subsequent subprocesses";
    statusNames[32] = "intermediate, subsequent subprocesses";
    statusNames[33] = "outgoing, subsequent subprocesses";
    statusNames[34] = " incoming that has already scattered";
    statusNames[41] = "incoming on spacelike main branch";
    statusNames[42] = "incoming copy of recoiler";
    statusNames[43] = "outgoing produced by a branching";
    statusNames[44] = "outgoing shifted by a branching";
    statusNames[45] = "incoming rescattered parton, with changed kinematics owing to ISR in the mother system (cf. status 34)";
    statusNames[46] = "incoming copy of recoiler when this is a rescattered parton (cf. status 42)";
    statusNames[47] = "a W or Z gauge boson produced in the shower evolution";
    statusNames[49] = "a special state in the evolution, where E^2 - p^2 = m^2 is not fulfilled";
    statusNames[51] = "outgoing produced by parton branching";
    statusNames[52] = "outgoing copy of recoiler, with changed momentum";
    statusNames[53] = "copy of recoiler when this is incoming parton, with changed momentum";
    statusNames[54] = "copy of a recoiler, when in the initial state of a different system from the radiator";
    statusNames[55] = "copy of a recoiler, when in the final state of a different system from the radiator";
    statusNames[56] = "a W or Z gauge boson produced in a shower branching (special case of 51)";
    statusNames[59] = "a special state in the evolution, where E^2 - p^2 = m^2 is not fulfilled";
    statusNames[61] = "particles produced by beam-remnant treatment";
    statusNames[61] = "incoming subprocess particle with primordial kT included";
    statusNames[62] = "outgoing subprocess particle with primordial kT included";
    statusNames[63] = "outgoing beam remnant";
    statusNames[64] = "copied particle with new colour according to the colour configuration of the beam remnant";
    statusNames[71] = "copied partons to collect into contiguous colour singlet";
    statusNames[72] = "copied recoiling singlet when ministring collapses to one hadron and momentum has to be reshuffled";
    statusNames[73] = "combination of very nearby partons into one";
    statusNames[74] = "combination of two junction quarks (+ nearby gluons) to a diquark (this index is not changed by recoils or recombinations)";
    statusNames[75] = "gluons split to decouple a junction-antijunction pair";
    statusNames[76] = "partons with momentum shuffled or a new colour to decouple junction-antijunction structures";
    statusNames[77] = "temporary opposing parton when fragmenting first two strings in to junction (should disappear again)";
    statusNames[78] = "temporary combined diquark end when fragmenting last string in to junction (should disappear again)";
    statusNames[79] = "copy of particle with new colour indices after the colour reconnection";
    statusNames[81] = "from ministring into one hadron";
    statusNames[82] = "from ministring into two hadrons";
    statusNames[83] = "from normal string (the difference between the two is technical, whether fragmented off from the top of the string system or from the bottom, useful for debug only)";
    statusNames[84] = "from normal string (the difference between the two is technical, whether fragmented off from the top of the string system or from the bottom, useful for debug only)";
    statusNames[85] = "primary produced hadrons in junction fragmentation of the first two string legs in to the junction, in order of treatment";
    statusNames[86] = "primary produced hadrons in junction fragmentation of the first two string legs in to the junction, in order of treatment";
    statusNames[87] = "primary produced baryon from a junction using standard junction fragmentation of fragmenting two shortest legs into a diquark (status 78)";
    statusNames[88] = "primary produced baryon from a junction using the pearl-on-a-string gluon approximation";
    statusNames[89] = "primary produced baryon from a junction in the ministring framework";
    statusNames[91] = "normal decay products";
    statusNames[92] = "decay products after oscillation B0 ↔ B0bar or B_s0 ↔ B_s0bar";
    statusNames[93] = "decay handled by an external program, normally or with oscillation";
    statusNames[94] = "decay handled by an external program, normally or with oscillation";
    statusNames[95] = "a forced decay handled by an external program, i.e. picked from a subset of the possible channels to enhance the rate of rare signals, normally or with oscillation";
    statusNames[96] = "a forced decay handled by an external program, i.e. picked from a subset of the possible channels to enhance the rate of rare signals, normally or with oscillation";
    statusNames[97] = "decay products from a resonance produced in rescattering";
    statusNames[99] = "particles with momenta shifted by Bose-Einstein effects (not a proper decay, but bookkept as an 1 → 1 such, happening after decays of short-lived resonances but before decays of longer-lived particles)";
    statusNames[101] = "particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles containing a very heavy flavour";
    statusNames[102] = "partons rearranged from the long-lived particle end to prepare for fragmentation from this end";
    statusNames[103] = "intermediate half-R-hadron formed when a colour octet particle (like the gluino) has been fragmented on one side, but not yet on the other";
    statusNames[104] = "an R-hadron";
    statusNames[105] = "partons or particles formed together with the R-hadron during the fragmentation treatment";
    statusNames[106] = "subdivision of an R-hadron into its flavour content, with momentum split accordingly, in preparation of the decay of the heavy new particle, if it is unstable";
    statusNames[107] = "two temporary leftover gluons joined into one in the formation of a gluino-gluon R-hadron";
    statusNames[121] = "deuteron formed by coalescence, and other particles formed simultaneously (like a gamma or pi)";
    statusNames[150] = "reserved for possible future expansion";
    statusNames[151] = "inelastic nondiffractive collisions";
    statusNames[152] = "elastic scattering";
    statusNames[153] = "single diffraction, first particle excited";
    statusNames[154] = "single diffraction, second particle excited";
    statusNames[155] = "double diffraction";
    statusNames[156] = "reserved for future central diffraction";
    statusNames[157] = "low-mass diffraction, forming explicit resonances (NN collisions only)";
    statusNames[158] = "annihilation of incoming quark-antiquark pairs";
    statusNames[159] = "formation of a resonance particle";
    statusNames[181] = "a nucleon hit in the medium, specifically the first (by order in the event record) nucleon hit in the current nucleus";
    statusNames[182] = "a further nucleon hit in a nucleus, subsequent to the first one above";
    statusNames[-84] = "beam remnants after the hard scattering";
    statusNames[-71] = "part of the beam remnant after parton shower evolution";
    statusNames[-63] = "intermediate state";
    statusNames[-12] = "decaying particles";
    statusNames[-11] = "intermediate state particles";

}
particleDictionarry::~particleDictionarry(){
    idABSNames.clear();
    statusNames.clear();
}

void particleDictionarry::printALL(){}
std::string particleDictionarry::getStatusName(int i){
//    if(statusNames.find(i) == statusNames.end())return std::to_string(i);
//    return "(" +std::to_string(i) + ")" + statusNames[i];
    return std::to_string(i);
}
std::string particleDictionarry::getIdABSName(int i){
    if (idABSNames.find(i) == idABSNames.end())return std::to_string(i);
    return "(" + std::to_string(i) + ")"+idABSNames[i];

}
