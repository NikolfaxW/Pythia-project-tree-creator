//
// Created by nikol on 11/14/2024.
//

#ifndef PYTHIA_PROJECT_TREE_CREATOR_MYINFO_H
#define PYTHIA_PROJECT_TREE_CREATOR_MYINFO_H

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"



//simple class to store additional info for fastjet
class MyInfo : public fastjet::PseudoJet::UserInfoBase {
public:
    MyInfo(const int &id, const int &i, const bool & is_charged) : _pdg_id(id), _id(i), _is_charged(is_charged) {}  //stores ddg codes of the particle

    int pdg_id() const { return _pdg_id; }

    int id() const { return _id; }

    bool isCharged() const { return _is_charged;}

protected:
    int _pdg_id;
    int _id;
    bool _is_charged;


};

#endif //PYTHIA_PROJECT_TREE_CREATOR_MYINFO_H
