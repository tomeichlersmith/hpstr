#include <iostream>

//HPSTR
#include "HpsEvent.h"
#include "Collections.h"
#include "MCParticle.h"
#include "MCTrackerHit.h"
#include "MCEcalHit.h"
#include "MCAnaHistos.h"
#include "Vertex.h"


//ROOT
#include "Processor.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TTree.h"
#include "TFile.h"

class TTree;

struct iDMCandidateEvent {
  MCParticle chi2;
  bool chi2_found;
  MCParticle rad_ele;
  bool rad_ele_found;
  MCParticle ele;
  bool ele_found;
  MCParticle pos;
  bool pos_found;
  Vertex vtx;

  void clear();
};

void iDMCandidateEvent::clear() {
  chi2_found = false;
  rad_ele_found = false;
  ele_found = false;
  pos_found = false;
  vtx.Clear();
  chi2.Clear();
  rad_ele.Clear();
  ele.Clear();
  pos.Clear();
}

/**
 * @brief Insert description here.
 * more details
 */
class iDMPreSelection : public Processor {
  public:
    iDMPreSelection(const std::string& name, Process& process)
      : Processor(name, process) {}
    ~iDMPreSelection() = default;
    virtual bool process(IEvent* ievent);
    virtual void initialize(TTree* tree);
    virtual void finalize() final override {}
    virtual void configure(const ParameterSet& parameters);

  private:
    iDMCandidateEvent event_;

    std::string anaName_{"recoHitAna"}; //!< description
    std::string partColl_{"MCParticle"}; //!< description
    std::string trkrHitColl_{"TrackerHits"}; //!< description
    std::string ecalHitColl_{"EcalHits"}; //!< description
    std::string analysis_{"vertex"}; //!< description

    int debug_{0}; //!< Debug Level
};

void iDMPreSelection::configure(const ParameterSet& parameters) {
  std::cout << "Configuring iDMPreSelection" << std::endl;
  try {
    debug_       = parameters.getInteger("debug");
    anaName_     = parameters.getString("anaName");
    partColl_    = parameters.getString("partColl");
    trkrHitColl_   = parameters.getString("trkrHitColl");
    ecalHitColl_   = parameters.getString("ecalHitColl");
    analysis_    = parameters.getString("analysis");
  } catch (std::runtime_error& error) {
    std::cout << error.what() << std::endl;
  }
}

void iDMPreSelection::initialize(TTree* tree) {
  tree->Branch("iDM", &event_);
}

bool iDMPreSelection::process(IEvent* ievent) {
  event_.clear();

  return true;
}

DECLARE_PROCESSOR(iDMPreSelection);
