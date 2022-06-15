/** 
 * @file ThreeProngTridentTrackAnalyzer.cxx
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

class ThreeProngTridentTracksAnalyzer : public Processor {
  std::shared_ptr<BaseSelector> event_selector_, cluster_selector_;
  std::shared_ptr<TrackHistos> histos_;

  std::vector<Particle*>* particles_{};
  std::string particle_coll_{"FinalStateParticles"};

  std::vector<CalCluster*>* clusters_{};
  std::string cluster_coll_{"RecoEcalClusters"};

  EventHeader* evth_{nullptr};

  double timeOffset_{-999};
  //In GeV. Default is 2016 value;
  double beamE_{2.3};
  int isData_{0};
  std::shared_ptr<AnaHelpers> _ah;
  //Debug level
  int debug_{0};
 public:
  ThreeProngTridentTracksAnalyzer(const std::string& name, Process& process);
  ~ThreeProngTridentTracksAnalyzer();
  virtual void configure(const ParameterSet& parameters) final override;
  virtual void initialize(TTree* tree) final override;
  virtual bool process(IEvent* ievent) final override;
  virtual void finalize() final override;
};

ThreeProngTridentTracksAnalyzer::ThreeProngTridentTracksAnalyzer(const std::string& name, Process& process) 
  : Processor(name,process) {}

ThreeProngTridentTracksAnalyzer::~ThreeProngTridentTracksAnalyzer(){}

void ThreeProngTridentTracksAnalyzer::configure(const ParameterSet& parameters) {
  debug_   = parameters.getInteger("debug");

  event_selector_ = std::make_shared<BaseSelector>("event_selection", 
      parameters.getString("event_selection"));
  event_selector_->setDebug(debug_);
      
  cluster_selector_ = std::make_shared<BaseSelector>("cluster_selection", 
      parameters.getString("cluster_selection"));
  cluster_selector_->setDebug(debug_);
      
  histos_ = std::make_shared<TrackHistos>("full_trident");
  histos_->debugMode(debug_>0);
  histos_->loadHistoConfig(parameters.getString("histo_cfg"));
  
  cluster_coll_ = parameters.getString("cluster_coll", cluster_coll_);
  particle_coll_ = parameters.getString("particle_coll", particle_coll_);

  timeOffset_ = parameters.getDouble("CalTimeOffset");
  beamE_  = parameters.getDouble("beamE", beamE_);
  isData_ = parameters.getInteger("isData", isData_);
}

void ThreeProngTridentTracksAnalyzer::initialize(TTree* tree) {
  _ah =  std::make_shared<AnaHelpers>();
  event_selector_->LoadSelection();
  cluster_selector_->LoadSelection();
  histos_->DefineHistos({"no_p_sum_cut","p_sum_cut"},"");
  
  //init Reading Tree
  tree->SetBranchAddress("EventHeader",&evth_);
  tree->SetBranchAddress(particle_coll_.c_str(), &particles_);
  tree->SetBranchAddress(cluster_coll_.c_str(), &clusters_);
}

bool ThreeProngTridentTracksAnalyzer::process(IEvent* ievent) { 
  /// not sure how to acquire a weight from the event header, so leaving this
  double weight = 1.;
  event_selector_->getCutFlowHisto()->Fill(0.,weight);

  std::vector<CalCluster*> electron_clusters, positron_clusters;
  /// IF BRANCH NOT SET CORRECTLY, THIS WILL SEG VIO
  for (CalCluster* c : *clusters_) {
    cluster_selector_->getCutFlowHisto()->Fill(0.);
    if (not cluster_selector_->passCutLt("max_energy_frac", c->getEnergy()/beamE_)) continue;
    if (not cluster_selector_->passCutGt("min_energy", c->getEnergy())) continue;
    if (cluster_selector_->passCutLt("electron_max_x", c->getX())) {
      electron_clusters.push_back(c);
    } else if (cluster_selector_->passCutGt("positron_min_x", c->getX())) {
      positron_clusters.push_back(c);
    }
  }

  // make sure at least 1 positron and 2 electrons
  if (not event_selector_->passCutGt("at_least_one_positron", positron_clusters.size(), weight))
    return true;

  if (not event_selector_->passCutGt("at_least_two_electrons", electron_clusters.size(), weight))
    return true;

  // use greater-than here so that the earlier indices in the vector get assigned
  // the larger energy clusters
  static auto energy_sort = [](const CalCluster* lhs, const CalCluster* rhs) {
    return lhs->getEnergy() > rhs->getEnergy();
  };

  std::sort(positron_clusters.begin(), positron_clusters.end(), energy_sort);
  std::sort(electron_clusters.begin(), electron_clusters.end(), energy_sort);

  /**
   * Call the highest energy positron candidate cluster the positron
   * and the two highest energy electron candidate clusters the electrons
   *
   * WARN: This is where we assume that the vectors are at least the correct size.
   */
  CalCluster* positron{positron_clusters.at(0)}, 
              electron0{electron_clusters.at(0)}, 
              electron1{electron_clusters.at(1)};

  /**
   * Deduce maximum time difference between any pair within the three clusters
   * we are focusing on.
   */
  double max_time_diff{0.};
  std::vector<CalCluster*> clusters_of_importance{positron, electron0, electron1};
  for (std::size_t i{0}; i < 3; i++) {
    for (std::size_t j{i}; j < 3; j++) {
      double time_diff{abs(
           clusters_of_importance.at(i)->getTime()
          -clusters_of_importance.at(j)->getTime()
          )};
      if (time_diff > max_time_diff) max_time_diff = time_diff;
    }
  }

  if (not event_selector_->passCutLt("max_cluster_time_diff", max_time_diff, weight))
    return true;

  double cluster_E_sum{0.};
  for (CalCluster* c : clusters_of_importance) {
    cluster_E_sum += c->getEnergy();
  }

  histos_->Fill1DHisto("cluster_E_sum_h", cluster_E_sum, weight);

  /*
  auto positron_trk{positron->getTrack()};
  histos_->Fill1DTrack(&positron_trk, weight, "no_p_sum_cut_positron_");
  histos_->Fill2DTrack(&positron_trk, weight, "no_p_sum_cut_positron_");

  for (Particle* p : electrons) {
    auto trk{p->getTrack()};
    histos_->Fill1DTrack(&trk, weight,"no_p_sum_cut_electrons_");
    histos_->Fill2DTrack(&trk, weight,"no_p_sum_cut_electrons_");
  }

  histos_->Fill1DHisto("no_p_sum_cut_Psum_h", total_4momentum.P(), weight);

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
  */
  
  return true;
}

void ThreeProngTridentTracksAnalyzer::finalize() {
  //TODO clean this up a little.
  outF_->cd();
  histos_->saveHistos(outF_,histos_->getName());
  outF_->cd(histos_->getName().c_str());
  event_selector_->getCutFlowHisto()->Write();
  cluster_selector_->getCutFlowHisto()->Write();
  outF_->Close();
}

DECLARE_PROCESSOR(ThreeProngTridentTracksAnalyzer);
