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
  Particle ele;
  Particle pos;

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

    double timeOffset_;
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
  timeOffset_ = parameters.getDouble("CalTimeOffset",timeOffset_);
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
  int isPair1Trigger = isData_ ? evth_->isPair1Trigger() : 1;
  if (!vtxSelector->passCutEq("Pair1_eq",isPair1Trigger,weight))
    return true;

  // require that we only have one candidate vertex in an event
  // we separate this into a lt and gt requirement so we can see
  // how many events have no vertices and how many have more than one
  if (!vtxSelector->passCutGt("nVtx_gt", vtxs_.size(), weight))
    return true;
  if (!vtxSelector->passCutEq("nVtx_lt", vtxs_.size(), weight))
    return true;

  event_.vtx = vtxs_.at(0);
  bool found_ele{false}, found_pos{false};
  for (std::size_t ipart{0}; ipart < event_.vtx.getParticles().size(); ++ipart) {
    Particle& p{event_.vtx.getParticles()[ipart]};
    if (p.getPDG() == +11) {
      event_.ele = p;
      found_ele = true;
    } else if (p.getPDG() == -11) {
      event_.pos = p;
      found_pos = true;
    }
  }
  if (not found_ele) {
    throw std::runtime_error("Unable to find electron part of vertex.");
  }
  if (not found_pos) {
    throw std::runtime_error("Unable to find positron part of vertex.");
  }

  Track ele_trk{event_.ele.getTrack()},
        pos_trk{event_.pos.getTrack()};

  //Ele Track Time
  if (!vtxSelector->passCutLt("eleTrkTime_lt",fabs(ele_trk.getTrackTime()),weight))
    return true;

  //Pos Track Time
  if (!vtxSelector->passCutLt("posTrkTime_lt",fabs(pos_trk.getTrackTime()),weight))
    return true;

  //Ele Track-cluster match
  if (!vtxSelector->passCutLt("eleTrkCluMatch_lt",event_.ele.getGoodnessOfPID(),weight))
    return true;

  //Pos Track-cluster match
  if (!vtxSelector->passCutLt("posTrkCluMatch_lt",event_.pos.getGoodnessOfPID(),weight))
    return true;

  //Require Positron Cluster exists
  if (!vtxSelector->passCutGt("posClusE_gt",event_.ele.getCluster().getEnergy(),weight))
    return true;

  //Require Positron Cluster does NOT exists
  if (!vtxSelector->passCutLt("posClusE_lt",event_.pos.getCluster().getEnergy(),weight))
    return true;

  double botClusTime = 0.0;
  if(event_.ele.getCluster().getPosition().at(1) < 0.0) 
    botClusTime = event_.ele.getCluster().getTime();
  else 
    botClusTime = event_.pos.getCluster().getTime();

  //Bottom Cluster Time
  if (!vtxSelector->passCutLt("botCluTime_lt", botClusTime, weight))
    return true;
  if (!vtxSelector->passCutGt("botCluTime_gt", botClusTime, weight))
    return true;

  double corr_eleClusterTime = event_.ele.getCluster().getTime() - timeOffset_;
  double corr_posClusterTime = event_.pos.getCluster().getTime() - timeOffset_;
  //Ele Pos Cluster Time Difference
  if (!vtxSelector->passCutLt("eleposCluTimeDiff_lt",fabs(corr_eleClusterTime - corr_posClusterTime),weight))
    return true;

  //Ele Track-Cluster Time Difference
  if (!vtxSelector->passCutLt("eleTrkCluTimeDiff_lt",fabs(ele_trk.getTrackTime() - corr_eleClusterTime),weight))
    return true;

  //Pos Track-Cluster Time Difference
  if (!vtxSelector->passCutLt("posTrkCluTimeDiff_lt",fabs(pos_trk.getTrackTime() - corr_posClusterTime),weight))
    return true;

  //Ele Track Quality - Chi2
  if (!vtxSelector->passCutLt("eleTrkChi2_lt",ele_trk.getChi2(),weight))
    return true;

  //Pos Track Quality - Chi2
  if (!vtxSelector->passCutLt("posTrkChi2_lt",pos_trk.getChi2(),weight))
    return true;

  //Ele Track Quality - Chi2Ndf
  if (!vtxSelector->passCutLt("eleTrkChi2Ndf_lt",ele_trk.getChi2Ndf(),weight))
    return true;

  //Pos Track Quality - Chi2Ndf
  if (!vtxSelector->passCutLt("posTrkChi2Ndf_lt",pos_trk.getChi2Ndf(),weight))
    return true;

  TVector3 ele_mom;
  //ele_mom.SetX(ele.getMomentum()[0]);
  //ele_mom.SetY(ele.getMomentum()[1]);
  //ele_mom.SetZ(ele.getMomentum()[2]);
  ele_mom.SetX(ele_trk.getMomentum()[0]);
  ele_mom.SetY(ele_trk.getMomentum()[1]);
  ele_mom.SetZ(ele_trk.getMomentum()[2]);

  TVector3 pos_mom;
  //pos_mom.SetX(pos.getMomentum()[0]);
  //pos_mom.SetY(pos.getMomentum()[1]);
  //pos_mom.SetZ(pos.getMomentum()[2]);
  pos_mom.SetX(pos_trk.getMomentum()[0]);
  pos_mom.SetY(pos_trk.getMomentum()[1]);
  pos_mom.SetZ(pos_trk.getMomentum()[2]);

  //Beam Electron cut
  if (!vtxSelector->passCutLt("eleMom_lt",ele_mom.Mag(),weight))
    return true;

  //Ele min momentum cut
  if (!vtxSelector->passCutGt("eleMom_gt",ele_mom.Mag(),weight))
    return true;

  //Pos min momentum cut
  if (!vtxSelector->passCutGt("posMom_gt",pos_mom.Mag(),weight))
    return true;

  //Ele nHits
  int ele2dHits = ele_trk.getTrackerHitCount();
  if (!ele_trk.isKalmanTrack())
      ele2dHits*=2;

  if (!vtxSelector->passCutGt("eleN2Dhits_gt",ele2dHits,weight))  {
    return true;
  }

  //Pos nHits
  int pos2dHits = pos_trk.getTrackerHitCount();
  if (!pos_trk.isKalmanTrack())
      pos2dHits*=2;

  if (!vtxSelector->passCutGt("posN2Dhits_gt",pos2dHits,weight))  {
    return true;
  }

  //Less than 4 shared hits for ele/pos track
  if (!vtxSelector->passCutLt("eleNshared_lt",ele_trk.getNShared(),weight)) {
    return true;
  }

  if (!vtxSelector->passCutLt("posNshared_lt",pos_trk.getNShared(),weight)) {
    return true;
  }


  //Vertex Quality
  if (!vtxSelector->passCutLt("chi2unc_lt",event_.vtx.getChi2(),weight))
    return true;

  //Max vtx momentum
  if (!vtxSelector->passCutLt("maxVtxMom_lt",(ele_mom+pos_mom).Mag(),weight))
    return true;

  //Min vtx momentum
  if (!vtxSelector->passCutGt("minVtxMom_gt",(ele_mom+pos_mom).Mag(),weight))
    return true;

  ana_tree_->Fill();
  return true;
}

void iDMPreSelection::finalize() {
  outF_->cd();
  ana_tree_->Write();
  vtxSelector->getCutFlowHisto()->Write();
}

DECLARE_PROCESSOR(iDMPreSelection);
