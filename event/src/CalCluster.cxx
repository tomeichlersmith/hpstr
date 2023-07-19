/**
 * @file CalCluster.cxx
 * @brief Class used to encapsulate calorimeter cluster information.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#include "CalCluster.h"

ClassImp(CalCluster)

CalCluster::CalCluster()
    : TObject() { 
}


//TODO Fix the relation between particles->CalClusters (same as tracks)
CalCluster::~CalCluster() {
    Clear();
}


void CalCluster::Clear(Option_t* /*option*/) {
    TObject::Clear();
    seed_hit_ = 0;
    hits_.clear();
}

void CalCluster::setPosition(const float* position) {
    x_ = position[0];
    y_ = position[1];
    z_ = position[2];
}

void CalCluster::addHit(std::size_t i) { 
    hits_.push_back(i);
}
