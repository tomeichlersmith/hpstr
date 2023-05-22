#include <iostream>

//HPSTR
#include "HpsEvent.h"
#include "Collections.h"
#include "MCParticle.h"
#include "MCTrackerHit.h"
#include "MCEcalHit.h"
#include "MCAnaHistos.h"


//ROOT
#include "Processor.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TTree.h"
#include "TFile.h"

class TTree;

/**
 * @brief Insert description here.
 * more details
 */
class iDM : public Processor {
  public:
    /**
     * @brief Constructor
     * 
     * @param name 
     * @param process 
     */
    iDM(const std::string& name, Process& process)
      : Processor(name, process) {}

    ~iDM() = default;

    /**
     * @brief description
     * 
     * @param ievent 
     * @return true 
     * @return false 
     */
    virtual bool process(IEvent* ievent);

    /**
     * @brief description
     * 
     * @param tree 
     */
    virtual void initialize(TTree* tree);

    /**
     * @brief description
     * 
     */
    virtual void finalize();

    /**
     * @brief description
     * 
     * @param parameters 
     */
    virtual void configure(const ParameterSet& parameters);

  private:

    /** Containers to hold histogrammer info */
    MCAnaHistos* histos{nullptr};
    std::string  histCfgFilename_; //!< description

    std::vector<MCParticle*>   * mcParts_{nullptr}; //!< description
    std::vector<MCTrackerHit*> * mcTrkrHits_{nullptr}; //!< description
    std::vector<MCEcalHit*>  * mcEcalHits_{nullptr}; //!< description

    std::string anaName_{"recoHitAna"}; //!< description
    std::string partColl_{"MCParticle"}; //!< description
    std::string trkrHitColl_{"TrackerHits"}; //!< description
    std::string ecalHitColl_{"EcalHits"}; //!< description
    std::string analysis_{"vertex"}; //!< description

    int debug_{0}; //!< Debug Level
};

void iDM::configure(const ParameterSet& parameters) {
  std::cout << "Configuring iDM" << std::endl;
  try {
    debug_       = parameters.getInteger("debug");
    anaName_     = parameters.getString("anaName");
    partColl_    = parameters.getString("partColl");
    trkrHitColl_   = parameters.getString("trkrHitColl");
    ecalHitColl_   = parameters.getString("ecalHitColl");
    histCfgFilename_ = parameters.getString("histCfg");
    analysis_    = parameters.getString("analysis");
  } catch (std::runtime_error& error) {
    std::cout << error.what() << std::endl;
  }
}

void iDM::initialize(TTree* tree) {
  // init histos
  histos = new MCAnaHistos(anaName_);
  histos->loadHistoConfig(histCfgFilename_);
  histos->DefineHistos();
  histos->Define2DHistos();

  // init TTree
  tree->SetBranchAddress(partColl_.c_str(), &mcParts_);
  if (tree->FindBranch(trkrHitColl_.c_str()))
    tree->SetBranchAddress(trkrHitColl_.c_str(), &mcTrkrHits_);
  else
    std::cout<<"WARNING: No tracker hit collection, will skip FillMCTrackerHits! "<<std::endl;

  if ( tree->FindBranch(ecalHitColl_.c_str()))
    tree->SetBranchAddress(ecalHitColl_.c_str(), &mcEcalHits_);
  else
    std::cout<<"WARNING: No Ecal hit collection, will skip FillMCEcalHits"<<std::endl;
}

bool iDM::process(IEvent* ievent) {
  /**
   * # Primary Electron 
   * has the lowest getID() in the event
   * AND
   * it is the only electron (getPDG()==11) that has a getGenStatus()==3
   *
   * # Gen-Level Electrons
   * The other electrons have (getPDG()==11) and (getGenStatus()==1).
   * The recoil electron has a lower getID() since it is copied into
   * the event before the produced electron.
   *
   * # Gen-Level Positron
   * The positron has (getPDG()==-11) and (getGenStatus()==1)
   *
   * (getGenStatus()==0) particles are Geant4-produced.
   */
  const MCParticle* primary_electron = nullptr,
                  * recoil_electron  = nullptr,
                  * prod_electron    = nullptr,
                  * prod_positron    = nullptr;
  for (const MCParticle* particle : *mcParts_) {
    if (particle->getGenStatus() == 0) continue;
    if (particle->getGenStatus() == 1) {
      if (particle->getPDG() == -11) prod_positron = particle;
      else if (particle->getPDG() == 11) {
        if (recoil_electron == nullptr) recoil_electron = particle;
        else if (prod_electron == nullptr) prod_electron = particle;
        // will need to swap after the loop in case order is fucked
        else {
          // more than outgoing electrons from the generator
        }
      }
    } 
    if (particle->getGenStatus() == 3 and particle->getPDG() == 11) {
      primary_electron = particle;
    }
  }

  if (prod_electron->getID() < recoil_electron->getID()) {
    std::swap(prod_electron, recoil_electron);
  }

  histos->FillMCParticles(mcParts_, analysis_);
  if(mcTrkrHits_)
    histos->FillMCTrackerHits(mcTrkrHits_);
  if(mcEcalHits_)
    histos->FillMCEcalHits(mcEcalHits_);

  return true;
}

void iDM::finalize() {
  histos->saveHistos(outF_, anaName_);
  delete histos;
  histos = nullptr;
}

DECLARE_PROCESSOR(iDM);
