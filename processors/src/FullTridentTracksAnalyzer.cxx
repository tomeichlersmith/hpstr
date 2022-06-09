/** 
 * @file FullTridentTrackAnalyzer.cxx
 * @author Tom Eichlersmith, UMN updated to 2022 hpstr
 */

//HPSTR
#include "HpsEvent.h"
#include "Collections.h"
#include "EventHeader.h"
#include "Track.h"
#include "Particle.h"
#include "Processor.h"
#include "BaseSelector.h"
#include "TrackHistos.h"
#include "AnaHelpers.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "Math/Vector4D.h"

#include <memory>
#include <iostream>
#include <utility> 

class FullTridentTracksAnalyzer : public Processor {
  std::shared_ptr<BaseSelector> event_selector_;
  std::shared_ptr<TrackHistos> histos_;
  std::string selection_cfg_, histo_cfg_;

  std::vector<Particle*>* particles_{};
  std::string particle_coll_{"FinalStateParticles"};

  EventHeader* evth_{nullptr};

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
  _ah =  std::make_shared<AnaHelpers>();
  
  event_selector_ = std::make_shared<BaseSelector>("event_selection", selection_cfg_);
  event_selector_->setDebug(debug_);
  event_selector_->LoadSelection();
      
  histos_ = std::make_shared<TrackHistos>("full_trident");
  histos_->debugMode(debug_>0);
  histos_->loadHistoConfig(histo_cfg_);
  histos_->DefineHistos({"no_p_sum_cut","p_sum_cut"},"");
  
  //init Reading Tree
  tree->SetBranchAddress("EventHeader",&evth_);
  tree->SetBranchAddress(particle_coll_.c_str(), &particles_);
}

bool FullTridentTracksAnalyzer::process(IEvent* ievent) { 
  /// not sure how to acquire a weight from the event header, so leaving this
  double weight = 1.;
  event_selector_->getCutFlowHisto()->Fill(0.,weight);

  ROOT::Math::PxPyPzEVector total_4momentum{};
  std::vector<ROOT::Math::PxPyPzEVector> electron_4momenta, positron_4momenta;
  std::vector<Particle*> electrons, positrons;
  /// IF BRANCH NOT SET CORRECTLY, THIS WILL SEG VIO
  for (Particle* p : *particles_) {
    auto mom{p->getTrack().getMomentum()};
    if (p->getPDG() == 11) {
      electrons.push_back(p);
      electron_4momenta.emplace_back(mom.at(0), mom.at(1), mom.at(2), p->getEnergy());
      total_4momentum += electron_4momenta.back();
    } else if (p->getPDG() == -11) {
      positrons.push_back(p);
      positron_4momenta.emplace_back(mom.at(0), mom.at(1), mom.at(2), p->getEnergy());
      total_4momentum += positron_4momenta.back();
    }
  }

  if (not event_selector_->passCutEq("one_positron", positrons.size(), weight)) {
    return true;
  }

  Particle* positron = positrons.at(0);

  if (not event_selector_->passCutEq("two_electrons", electrons.size(), weight)) {
    return true;
  }

  auto positron_trk{positron->getTrack()};
  histos_->Fill1DTrack(&positron_trk, weight, "no_p_sum_cut_positron_");
  histos_->Fill2DTrack(&positron_trk, weight, "no_p_sum_cut_positron_");

  for (Particle* p : electrons) {
    auto trk{p->getTrack()};
    histos_->Fill1DTrack(&trk, weight,"no_p_sum_cut_electrons_");
    histos_->Fill2DTrack(&trk, weight,"no_p_sum_cut_electrons_");
  }

  histos_->Fill1DHisto("no_p_sum_cut_Psum_h", total_4momentum.P(), weight);
  histos_->Fill1DHisto("no_p_sum_cut_cluster_Esum_h", total_4momentum.E(), weight);

  if (not event_selector_->passCutLt("max_p_sum", total_4momentum.P(), weight)) {
    return true;
  }

  if (not event_selector_->passCutGt("min_p_sum", total_4momentum.P(), weight)) {
    return true;
  }

  histos_->Fill1DTrack(&positron_trk, weight, "p_sum_cut_positron_");
  histos_->Fill2DTrack(&positron_trk, weight, "p_sum_cut_positron_");

  for (Particle* p : electrons) {
    auto trk{p->getTrack()};
    histos_->Fill1DTrack(&trk, weight,"p_sum_cut_electrons_");
    histos_->Fill2DTrack(&trk, weight,"p_sum_cut_electrons_");
  }

  histos_->Fill1DHisto("p_sum_cut_Psum_h", total_4momentum.P(), weight);
  histos_->Fill1DHisto("p_sum_cut_cluster_Esum_h", total_4momentum.E(), weight);
  
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
