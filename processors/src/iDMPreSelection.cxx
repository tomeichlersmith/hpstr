#include <iostream>
#include <set>

//HPSTR
#include "HpsEvent.h"
#include "Collections.h"
#include "EventHeader.h"
#include "TSData.h"
#include "Vertex.h"
#include "Track.h"
#include "TrackerHit.h"
#include "MCParticle.h"
#include "Particle.h"
#include "Processor.h"
#include "BaseSelector.h"
#include "TrackHistos.h"
#include "MCAnaHistos.h"


//ROOT
#include "Processor.h"
#include "TBranch.h"
#include "TTree.h"
#include "TFile.h"

class TTree;

struct iDMCandidateEvent {
  int ievent;
  Vertex vtx;

  void attach(TTree* t);
  void clear();
};

void iDMCandidateEvent::attach(TTree* t) {
  t->Branch("ievent", &ievent);
  t->Branch("vtx", &vtx);
}

void iDMCandidateEvent::clear() {
  vtx.Clear();
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
    virtual void finalize() final override;
    virtual void configure(const ParameterSet& parameters);

  private:
    std::string selectionCfg_; //!< description
    std::shared_ptr<BaseSelector> vtxSelector; //!< description

    TTree *ana_tree_;
    iDMCandidateEvent event_;
    EventHeader *evth_{}; //!< description
    std::vector<CalCluster> ecal_{}; //!< description
    std::vector<Vertex> vtxs_{}; //!< description
    std::vector<Track> trks_{}; //!< description
    std::vector<TrackerHit> hits_{}; //!< description
    std::vector<MCParticle> mcParts_{}; //!< description

    std::string anaName_{"vtxAna"}; //!< description
    std::string analysis_;
    std::string vtxColl_{"Vertices"}; //!< description
    std::string hitColl_{"RotatedHelicalTrackHits"}; //!< description
    std::string trkColl_{"GBLTracks"}; //!< description
    std::string ecalColl_{"RecoEcalClusters"}; //!< description
    std::string mcColl_{"MCParticle"}; //!< description

    int isData_;
    int debug_{0}; //!< Debug Level
};

void iDMPreSelection::configure(const ParameterSet& parameters) {
  std::cout << "Configuring iDMPreSelection" << std::endl;
  debug_       = parameters.getInteger("debug");
  anaName_ = parameters.getString("anaName",anaName_);
  vtxColl_ = parameters.getString("vtxColl",vtxColl_);
  trkColl_ = parameters.getString("trkColl",trkColl_);
  hitColl_ = parameters.getString("hitColl",hitColl_);
  ecalColl_ = parameters.getString("ecalColl",ecalColl_);
  mcColl_  = parameters.getString("mcColl",mcColl_);
  isData_  = parameters.getInteger("isData",isData_);
  analysis_        = parameters.getString("analysis");
  selectionCfg_   = parameters.getString("vtxSelectionjson",selectionCfg_);
}

template <typename T>
void set_and_check(TTree* tree, const std::string& branch, T& addr) {
  TBranch* br{tree->GetBranch(branch.c_str())};
  if (br == nullptr) {
    throw std::runtime_error("Unable to set "+branch+" on input TTree.");
  }
  br->SetObject(&addr);
}

void iDMPreSelection::initialize(TTree* tree) {
  vtxSelector  = std::make_shared<BaseSelector>(anaName_+"_"+"vtxSelection",selectionCfg_);
  vtxSelector->setDebug(debug_);
  vtxSelector->LoadSelection();

  ana_tree_ = new TTree("iDMCandidates","iDMCandidates");
  event_.attach(ana_tree_);

  //init Reading Tree
  // Get list of branches in tree to help protect accessing them
  int nBr = tree->GetListOfBranches()->GetEntries();
  if (debug_) std::cout << "Tree has " << nBr << " branches" << std::endl;
  std::set<std::string> branches;
  for(int iBr = 0; iBr < nBr; iBr++) {
    TBranch *br = dynamic_cast<TBranch*>(tree->GetListOfBranches()->At(iBr));
    branches.insert(br->GetName());
    if (debug_) std::cout << br->GetName() << std::endl;
  }
  tree->SetBranchAddress("EventHeader", &evth_);
  set_and_check(tree, vtxColl_, vtxs_);
  if (branches.find(hitColl_) != branches.end()) {
    set_and_check(tree, hitColl_, hits_);
  }
  set_and_check(tree, ecalColl_, ecal_);
  if(!isData_ && !mcColl_.empty()) {
    set_and_check(tree, mcColl_, mcParts_);
  }
  //If track collection name is empty take the tracks from the particles. 
  //TODO:: change this
  if (!trkColl_.empty()) {
    set_and_check(tree, trkColl_, trks_);
  }
}

bool iDMPreSelection::process(IEvent* ievent) {
  event_.clear();
  event_.ievent = evth_->getEventNumber();

  if(debug_) {
    std:: cout << "----------------- Event " 
      << evth_->getEventNumber() << " -----------------" << std::endl;
  }
  HpsEvent* hps_evt = (HpsEvent*) ievent;
  double weight = 1.;
  vtxSelector->getCutFlowHisto()->Fill(0.,weight);

  // double check that we have the right trigger
  // this handles the data mixin where other triggers are kept alongside pair1
  if (!vtxSelector->passCutEq("Pair1_eq",(int)evth_->isPair1Trigger(),weight))
    return true;

  // require that we only have one candidate vertex in an event
  if (!vtxSelector->passCutEq("nVtxs_eq", vtxs_.size(), weight))
    return true;

  event_.vtx = vtxs_.at(0);

  ana_tree_->Fill();
  return true;
}

void iDMPreSelection::finalize() {
  outF_->cd();
  ana_tree_->Write();
  vtxSelector->getCutFlowHisto()->Write();
}

DECLARE_PROCESSOR(iDMPreSelection);
