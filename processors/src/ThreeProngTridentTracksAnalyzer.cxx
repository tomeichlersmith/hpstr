/** 
 * @file ThreeProngTridentTrackAnalyzer.cxx
 * @author Tom Eichlersmith, UMN updated to 2022 hpstr
 *
 * This analysis requires a few different reconstructed objects in order to function properly.
 *  - RecoEcalClusters : needed for doing cuts on clusters and candidate trident clusters
 *  - RecoEcalHits : not explicitly listed but needed to load the seed hit of the ecal clusters
 *  - FinalStateParticles : needed for checking if clusters having a matching track
 *  - KalmanFullTracks : needed for looking at tracks in detail after the cluster selection
 *
 * This analysis requires an update written by Cam copying a unique ID for each cluster in LCIO
 * into the ROOT tuples. This makes finding the FSP for a specific cluster much easier.
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
#include "CalHit.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "Math/Vector4D.h"

#include <memory>
#include <iostream>
#include <utility> 
#include <exception>

class ThreeProngTridentTracksAnalyzer : public Processor {
  std::shared_ptr<BaseSelector> event_selector_, cluster_selector_;
  std::shared_ptr<TrackHistos> histos_;

  std::vector<Particle*>* particles_{};
  std::string particle_coll_{"FinalStateParticles"};

  std::vector<CalCluster*>* clusters_{};
  std::string cluster_coll_{"RecoEcalClusters"};

  std::vector<Track*>* tracks_{};
  std::string track_coll_{"KalmanFullTracks"};

  EventHeader* evth_{nullptr};

  double timeOffset_{-999};
  //In GeV. Default is 2016 value;
  double beamE_{2.3};
  int isData_{0};
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
      
  histos_ = std::make_shared<TrackHistos>("tpt");
  histos_->debugMode(debug_>0);
  histos_->loadHistoConfig(parameters.getString("histo_cfg"));
  
  cluster_coll_ = parameters.getString("cluster_coll", cluster_coll_);
  particle_coll_ = parameters.getString("particle_coll", particle_coll_);
  track_coll_ = parameters.getString("track_coll", track_coll_);

  timeOffset_ = parameters.getDouble("CalTimeOffset");
  beamE_  = parameters.getDouble("beamE", beamE_);
  isData_ = parameters.getInteger("isData", isData_);
}

void ThreeProngTridentTracksAnalyzer::initialize(TTree* tree) {
  event_selector_->LoadSelection();
  cluster_selector_->LoadSelection();

  /**
   * Histos are defined here - any histograms with the string 'follow'
   * within their name will have multiple copies created in order to "follow"
   * the event selections made within this analyzer
   */
  histos_->DefineHistos({"pre_time_cut","pre_fiducial_cut","pre_esum_cut","final_selection"},"follow");
  
  //init Reading Tree
  tree->SetBranchAddress("EventHeader",&evth_);
  tree->SetBranchAddress(particle_coll_.c_str(), &particles_);
  tree->SetBranchAddress(cluster_coll_.c_str(), &clusters_);
  tree->SetBranchAddress(track_coll_.c_str(), &tracks_);
}

bool ThreeProngTridentTracksAnalyzer::process(IEvent* ievent) { 
  /// not sure how to acquire a weight from the event header, so leaving this
  double weight = 1.;
  event_selector_->getCutFlowHisto()->Fill(0.,weight);

  std::vector<CalCluster*> electron_clusters, positron_clusters;
  /// IF BRANCH NOT SET CORRECTLY, THIS WILL SEG VIO
  for (CalCluster* c : *clusters_) {
    cluster_selector_->getCutFlowHisto()->Fill(0.);
    if (not cluster_selector_->passCutLt("max_energy_frac", c->getEnergy()/beamE_,1)) continue;
    if (not cluster_selector_->passCutGt("min_energy", c->getEnergy(),1)) continue;
    double x{c->getPosition().at(0)};
    if (cluster_selector_->passCutLt("electron_max_x", x,1)) {
      electron_clusters.push_back(c);
    } else if (cluster_selector_->passCutGt("positron_min_x", x,1)) {
      positron_clusters.push_back(c);
    }
  }

  histos_->Fill1DHisto("n_positron_candidates_h", positron_clusters.size(), weight);
  histos_->Fill1DHisto("n_electron_candidates_h", electron_clusters.size(), weight);

  // make sure at least 1 positron and 2 electrons
  if (not event_selector_->passCutGt("min_positron_candidates", positron_clusters.size(), weight))
    return true;

  if (not event_selector_->passCutGt("min_electron_candidates", electron_clusters.size(), weight))
    return true;

  if (not event_selector_->passCutLt("max_positron_candidates", positron_clusters.size(), weight))
    return true;

  if (not event_selector_->passCutLt("max_electron_candidates", electron_clusters.size(), weight))
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
            * electron0{electron_clusters.at(0)}, 
            * electron1{electron_clusters.at(1)};

  static auto extract_seed_indices = [](const CalCluster* cluster) {
    auto seed = dynamic_cast<CalHit*>(cluster->getSeed());
    if (not seed) {
      throw std::runtime_error("Cluster seed hit not in ROOT file. RecEcalHits need to be in tuple file.");
    }
    return seed->getCrystalIndices();
  };

  /**
   * Calculate meta-variables of the three trident clusters
   *
   * max_time_diff - maximum abs difference in time between any pair within the three clusters
   * cluster_E_sum - total energy from three clusters
   * cluster_on_edge - count of clusters that are seeded on the edge of the ECal, "fiducial" clusters
   */
  double max_time_diff{0.};
  double cluster_E_sum{0.};
  int clusters_on_edge{0};
  std::vector<CalCluster*> trident_clusters{positron, electron0, electron1};
  for (std::size_t i{0}; i < 3; i++) {
    CalCluster* c{trident_clusters.at(i)};
    // energy
    cluster_E_sum += c->getEnergy();

    // fiducial
    std::vector<int> indices{extract_seed_indices(c)};
    int x{indices.at(0)}, y{indices.at(1)};
    // left or right outer edge
    if (abs(x) == 23) ++clusters_on_edge;
    // top or bottom outer edge
    else if (abs(y) == 5) ++ clusters_on_edge;
    // middle gap outside cutout
    else if (abs(y) == 1 and (x > -2 or x < -10)) ++clusters_on_edge;
    // middle gap inside cutout
    else if (abs(y) == 2 and (x < -1 and x > -11)) ++clusters_on_edge;

    // time difference
    for (std::size_t j{i}; j < 3; j++) {
      double time_diff{abs(c->getTime()
          -trident_clusters.at(j)->getTime()
          )};
      if (time_diff > max_time_diff) max_time_diff = time_diff;
    }
  }

  static auto fill = [&](const std::string& name) {
    for (CalCluster* c : trident_clusters) {
      std::vector<int> indices{extract_seed_indices(c)};
      int x{indices.at(0)}, y{indices.at(1)};
      histos_->Fill2DHisto(name+"_follow_cluster_seed_pos_hh",x,y);
    }

    histos_->Fill2DHisto(name+"_follow_time_diff_hh",
        abs(positron->getTime()  - electron0->getTime()),
        abs(electron1->getTime() - electron0->getTime()));

    histos_->Fill1DHisto(name+"_follow_max_time_diff_h", max_time_diff, weight);
    histos_->Fill1DHisto(name+"_follow_clusters_on_edge_h", clusters_on_edge, weight);
    histos_->Fill2DHisto(name+"_follow_max_time_diff_vs_E_sum_hh", 
        max_time_diff, cluster_E_sum, weight);
    histos_->Fill1DHisto(name+"_follow_cluster_E_sum_h", cluster_E_sum, weight);
    histos_->Fill1DHisto(name+"_follow_electron0_cluster_E_h", electron0->getEnergy(), weight);
    histos_->Fill1DHisto(name+"_follow_electron1_cluster_E_h", electron1->getEnergy(), weight);
    histos_->Fill1DHisto(name+"_follow_positron_cluster_E_h", positron->getEnergy(), weight);
  };

  fill("pre_fiducial_cut");

  if (not event_selector_->passCutEq("no_edge_clusters", clusters_on_edge, weight)) 
    return true;

  fill("pre_time_cut");

  if (not event_selector_->passCutLt("max_cluster_time_diff", max_time_diff, weight))
    return true;

  fill("pre_esum_cut");

  if (not event_selector_->passCutLt("max_cluster_E_sum", cluster_E_sum, weight)
      or not event_selector_->passCutGt("min_cluster_E_sum", cluster_E_sum, weight))
    return true;

  fill("final_selection");

  /**
   * Get matching tracks for the clusters now that we have a final selection
   *
   * We determine matching clusters by looking through the particle collection until
   * a cluster with the same ID is found. The track stored by that particle is then used
   * for track studying.
   */
  Track positron_trk, electron0_trk, electron1_trk;
  for (Particle* p : *particles_) {
    int id{p->getCluster().getID()};
    if (id == positron->getID())
      positron_trk = p->getTrack();
    else if (id == electron0->getID())
      electron0_trk = p->getTrack();
    else if (id == electron1->getID())
      electron1_trk = p->getTrack();
  }

  static auto fill_trk_match_histos = [&](const std::string& name, 
      bool successful_match, CalCluster* c) {
    histos_->Fill1DHisto(name+"_track_N_h", int(successful_match), weight);
    if (successful_match) {
      histos_->Fill2DHisto(name+"_matched_x_y_hh",
          c->getPosition().at(0), c->getPosition().at(1), weight);
    } else {
      histos_->Fill2DHisto(name+"_unmatched_x_y_hh",
          c->getPosition().at(0), c->getPosition().at(1), weight);
    }
  };

  fill_trk_match_histos("electron0", electron0_trk.getID() > 0, electron0);
  fill_trk_match_histos("electron1", electron1_trk.getID() > 0, electron1);
  fill_trk_match_histos("positron", positron_trk.getID() > 0, positron);

  int num_clusters_with_track{0};
  if (positron_trk.getID() > 0) num_clusters_with_track++;

  if (not event_selector_->passCutGt("min_positron_with_track", num_clusters_with_track, weight))
    return true;

  if (electron0_trk.getID() > 0) num_clusters_with_track++;
  if (electron1_trk.getID() > 0) num_clusters_with_track++;

  if (not event_selector_->passCutGt("min_clusters_with_track", num_clusters_with_track, weight))
    return true;

  /**
   * Track Distributions
   *
   * Now that we have selected events that are candidate three prong tridents,
   * we can get the matching tracks and fill up some histograms of their parameters.
   */
  static auto fill_track_histos = [&](Track& trk, const std::string& name) {
    // silent leave if track wasn't found from its cluster
    if (trk.getID() < 0) return;
    // de-reference iterator to pass into function
    histos_->Fill1DTrack(&trk, weight, name+"_");
    histos_->Fill2DTrack(&trk, weight, name+"_");
  };

  fill_track_histos(positron_trk, "positron");
  fill_track_histos(electron0_trk, "electron0");
  fill_track_histos(electron1_trk, "electron1");
  
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
