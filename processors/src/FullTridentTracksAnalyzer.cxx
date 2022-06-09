/** 
 * @file FullTridentTrackAnalyzer.cxx
 * @author Tom Eichlersmith, UMN updated to 2022 hpstr
 */

//HPSTR
#include "HpsEvent.h"
#include "Collections.h"
#include "EventHeader.h"
#include "Vertex.h"
#include "Track.h"
#include "TrackerHit.h"
#include "Particle.h"
#include "Processor.h"
#include "BaseSelector.h"
#include "TrackHistos.h"
#include "FlatTupleMaker.h"
#include "AnaHelpers.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TVector3.h"

#include <memory>
#include <iostream>
#include <utility> 

class FullTridentTracksAnalyzer : public Processor {
  std::shared_ptr<BaseSelector> event_selector_;
  std::shared_ptr<HistoManager> histos_;
  std::string selection_cfg_, histo_cfg_;

  std::vector<Particle*>* particles_{};
  std::string particle_coll_{"FinalStateParticles"};

  std::vector<Vertex*> * vtxs_{};
  std::vector<Track*>  * trks_{};
  EventHeader* evth_{nullptr};
  std::string vtxColl_{"Vertices"};
  std::string trkColl_{"GBLTracks"};
  TTree* tree_{nullptr};

  double timeOffset_{-999};
  //In GeV. Default is 2016 value;
  double beamE_{2.3};
  int isData{0};
  std::shared_ptr<AnaHelpers> _ah;
  //Debug level
  int debug_{0};
 public:
  FullTridentTracksAnalyzer(const std::string& name, Process& process);
  ~FullTridentTracksAnalyzer();
  virtual void configure(const ParameterSet& parameters) final override;
  virtual void initialize(TTree* tree) final override;
  virtual bool process(IEvent* ievent) final override;
  virtual void finalize() final override;
};

FullTridentTracksAnalyzer::FullTridentTracksAnalyzer(const std::string& name, Process& process) 
  : Processor(name,process) {}

FullTridentTracksAnalyzer::~FullTridentTracksAnalyzer(){}

void FullTridentTracksAnalyzer::configure(const ParameterSet& parameters) {
  debug_   = parameters.getInteger("debug");
  selection_cfg_   = parameters.getString("event_selection");
  histo_cfg_ = parameters.getString("histo_cfg");

  timeOffset_ = parameters.getDouble("CalTimeOffset");
  beamE_  = parameters.getDouble("beamE");
  isData  = parameters.getInteger("isData");
}

void FullTridentTracksAnalyzer::initialize(TTree* tree) {
  tree_ = tree;
  _ah =  std::make_shared<AnaHelpers>();
  
  event_selector_ = std::make_shared<BaseSelector>("event_selection", selection_cfg_);
  event_selector_->setDebug(debug_);
  event_selector_->LoadSelection();
      
  histos_ = std::make_shared<HistoManager>("event_selection");
  histos_->loadHistoConfig(histo_cfg_);
  histos_->DefineHistos();
  
  //init Reading Tree
  tree_->SetBranchAddress(vtxColl_.c_str(), &vtxs_);
  tree_->SetBranchAddress("EventHeader",&evth_);
  tree_->SetBranchAddress(particle_coll_.c_str(), &particles_);
}

bool FullTridentTracksAnalyzer::process(IEvent* ievent) { 
  /// not sure how to acquire a weight from the event header, so leaving this
  double weight = 1.;
  event_selector_->getCutFlowHisto()->Fill(0.,weight);

  int n_electrons{0}, n_positrons{0};
  /// IF BRANCH NOT SET CORRECTLY, THIS WILL SEG VIO
  for (Particle* p : *particles_) {
    if (p->getPDG() == 11) {
      n_electrons++;
    } else if (p->getPDG() == -11) {
      n_positrons++;
    }
  }

  if (not event_selector_->passCutEq("one_positron", n_positrons, weight)) {
    return true;
  }

  if (not event_selector_->passCutEq("two_electrons", n_electrons, weight)) {
    return true;
  }

  for (Particle* p : *particles_) {
    if (abs(p->getPDG()) != 11) continue;
    histos_->Fill1DHisto("all_elepos_d0_h", p->getTrack().getD0(), weight);
  }
    
  return true;
}

void FullTridentTracksAnalyzer::finalize() {
  //TODO clean this up a little.
  outF_->cd();
  histos_->saveHistos(outF_,histos_->getName());
  outF_->cd(histos_->getName().c_str());
  event_selector_->getCutFlowHisto()->Write();
  outF_->Close();
}

DECLARE_PROCESSOR(FullTridentTracksAnalyzer);
