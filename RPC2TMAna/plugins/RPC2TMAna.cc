// -*- C++ -*-
//
// Package:    rpc2tmTime/RPC2TMAna
// Class:      RPC2TMAna
//
/**\class RPC2TMAna RPC2TMAna.cc rpc2tmTime/RPC2TMAna/plugins/RPC2TMAna.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ahmed Hussein
//         Created:  Mon, 20 Jun 2022 22:31:04 GMT
//
//

// system include files
#include <memory>
#include <iostream> // for priting

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
//#include "DataFormats/TrackReco/interface/Track.h"      // why commented out?
//#include "DataFormats/TrackReco/interface/TrackFwd.h"   // why commented out?

// Add the needed include files
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"

#include "L1Trigger/L1TTwinMux/interface/RPCHitCleaner.h"

#include "DataFormats/L1TMuon/interface/CPPFDigi.h"

// Header files for TFileService.
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Header file for ROOT histogramming.
#include "TH1.h"
#include "TH2.h"

// Header Files for the algo.
#include "L1Trigger/L1TTwinMux/interface/RPCHitCleaner.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class RPC2TMAna : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit RPC2TMAna(const edm::ParameterSet&);
  ~RPC2TMAna();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file

  // RPC part
  edm::EDGetTokenT<RPCDigiCollection> rpcLegacyToken_;
  edm::EDGetTokenT<RPCDigiCollection> rpcTwinMuxToken_;
  //edm::EDGetTokenT<l1t::CPPFDigiCollection> rpcCPPFToken_;

  // DT part
  edm::EDGetTokenT<L1MuDTChambPhContainer> m_tm_phiIn_Token_;
  edm::EDGetTokenT<L1MuDTChambPhContainer> m_tm_phiOut_Token_;
  edm::EDGetTokenT<L1MuDTChambThContainer> m_tm_theta_Token_;

  // Histograms
  TH1D* hist_phiInSize;
  TH1D* hist_phiOutSize;
  TH1D* hist_thetaSize;
  TH1D* hist_clusterSize_RPCTwinMux_before;
  TH1D* hist_clusterSize_RPCTwinMux_after;
  TH2D* hist_clusterSize_bx;
  TH2D* hist_clusterSize_bx_after;

  // Create a vector to store the region of each hit.
  std::vector<int> region_v;

  // Create a vector to store the bx of each hit.
  std::vector<int> bx_v;

  // Create a vector to store the strip of each hit.
  std::vector<int> strip_v;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RPC2TMAna::RPC2TMAna(const edm::ParameterSet& iConfig):
  rpcLegacyToken_(consumes<RPCDigiCollection>(iConfig.getUntrackedParameter<edm::InputTag>("rpcLegacy"))),
  rpcTwinMuxToken_(consumes<MuonDigiCollection<RPCDetId,RPCDigi>>(iConfig.getUntrackedParameter<edm::InputTag>("rpcTwinMux"))),
  //rpcCPPFToken_(consumes<l1t::CPPFDigiCollection>(iConfig.getUntrackedParameter<edm::InputTag>("rpccppf"))),

  m_tm_phiIn_Token_(consumes<L1MuDTChambPhContainer>(iConfig.getUntrackedParameter<edm::InputTag>("inputTagTMphIn"))),
  m_tm_phiOut_Token_(consumes<L1MuDTChambPhContainer>(iConfig.getUntrackedParameter<edm::InputTag>("inputTagTMphOut"))),
  m_tm_theta_Token_(consumes<L1MuDTChambThContainer>(iConfig.getUntrackedParameter<edm::InputTag>("inputTagTMth"))){


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

RPC2TMAna::~RPC2TMAna() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

namespace {
  constexpr int max_rpc_bx = 3;
  constexpr int min_rpc_bx = -3;

  /// Need to shift the index so that index 0 corresponds to min_rpc_bx.
  // Define the class of "BXToStrips"
  class BXToStrips {
  public:
    BXToStrips() : m_strips{} {} //zero initializes

    // A fn that returns true if bx is out of the range [min_rpc_bx, max_rpc_bx].
    static bool outOfRange(int iBX) {return (iBX > max_rpc_bx or iBX < min_rpc_bx);}

    // Make operator such that m_strips[iBX] gives the #strips for iBX.
    int& operator[] (int iBX) {return m_strips[iBX - min_rpc_bx];}

    size_t size() const {return m_strips.size();}

  private:
    // m_strips is an array of type int and of size (max_rpc_bx - min_rpc_bx + 1)
    // m_strips is used to store the #strips for a specific bx.
    // Each bx is assigned an index [-3,3]->[1-6]
    std::array<int, max_rpc_bx - min_rpc_bx + 1> m_strips;
  }; // class
} //namespace

// ------------ method called for each event  ------------
void RPC2TMAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  bool debug = false;

  edm::Handle< RPCDigiCollection > digiCollectionRPCLegacy;
  iEvent.getByToken(rpcLegacyToken_,digiCollectionRPCLegacy);

  edm::Handle< MuonDigiCollection<RPCDetId,RPCDigi> > digiCollectionRPCTwinMux;
  iEvent.getByToken(rpcTwinMuxToken_,digiCollectionRPCTwinMux);

  //edm::Handle<l1t::CPPFDigiCollection> digiCollectionCPPF;
  //iEvent.getByToken(rpcCPPFToken_,digiCollectionCPPF);

  edm::Handle<L1MuDTChambPhContainer> phiInTrigsTM; // https://cmssdt.cern.ch/lxr/source/DataFormats/Common/interface/Handle.h
  iEvent.getByToken(m_tm_phiIn_Token_, phiInTrigsTM);

  edm::Handle<L1MuDTChambPhContainer> phiOutTrigsTM;
  iEvent.getByToken(m_tm_phiOut_Token_, phiOutTrigsTM);

  edm::Handle<L1MuDTChambThContainer> thetaTrigsTM;
  iEvent.getByToken(m_tm_theta_Token_, thetaTrigsTM);

  vector<L1MuDTChambPhDigi> const* v_TwinMuxPhi_in = phiInTrigsTM->getContainer();	//pure DT only tp
  if (debug) std::cout << "container size TwinMuxPhi_in = " << v_TwinMuxPhi_in->size() << std::endl;
  hist_phiInSize->Fill(v_TwinMuxPhi_in->size());

  vector<L1MuDTChambPhDigi> const* v_TwinMuxPhi_out = phiOutTrigsTM->getContainer();	// combination, depend on the TM configuration
  if (debug) std::cout << "container size TwinMuxPhi_out = " << v_TwinMuxPhi_out->size() << std::endl;
  hist_phiOutSize->Fill(v_TwinMuxPhi_out->size());

  vector<L1MuDTChambThDigi> const* v_TwinMuxTheta = thetaTrigsTM->getContainer();
  if (debug) std::cout << "container size TwinMuxTheta = " << v_TwinMuxTheta->size() << std::endl;
  hist_thetaSize->Fill(v_TwinMuxTheta->size());


/////////////////// Taken from RPCHitCleaner.cc //////////////////////

// Some definitions.

// Instance to store the cleaned RPC Digis.
RPCDigiCollection m_outrpcDigis;
// map to store the rpcDetId, bx, strip and the index of the cluster where each digi was clusterized.
// key = struct(rpcDetId, bx, strip)
// value = index of the cluster where each digi was clusterized
std::map<RPCHitCleaner::detId_Ext, int> hits;

// vector for the cluster size of the clusters in each event, ordered by cluster index
std::vector<int> vcluster_size;

// map to store thr RPC DetId for a chamber and the min bx of a digi in the chamber's digi collection.
std::map<RPCDetId, int> bx_hits;

int cluster_size = 0; // assigned zero for each new cluster.
int cluster_id = -1;
int itr = 0;

// map to store values of bx for each cluster_size.
std::map<int, std::vector<int>> clusterSize_bx;


// First: Select digis unpacked from RPC TwinMux digi collection.
edm::Handle<RPCDigiCollection> m_inrpcDigis = digiCollectionRPCTwinMux;

/////// Test
/*int chamber_no = 1;
for(auto chamber = m_inrpcDigis->begin(); chamber != m_inrpcDigis->end(); ++chamber) {
  std::cout << "\n\n\n\n\n\n Chamber #" << chamber_no << ": {";
  for(auto digi = (*chamber).second.first; digi != (*chamber).second.second; ++digi) {
    std::cout << "(" << digi->bx() << "," << digi->strip() << "), ";
  }
  std::cout << "}\n";
  chamber_no++;
}*/
/////// End of test

int bx_n1 = -10000;

// First Loop through the chmabers. RPCDetId specifies a chamber.
for(auto chamber = m_inrpcDigis->begin(); chamber != m_inrpcDigis->end(); ++chamber) {
  //region_v.push_back(rpcDetId.region());
  RPCDetId rpcDetId = (*chamber).first;
  int strip_n1 = -10000;
  //int bx_n1 = -10000;

  // Select Barrel hits only.
  if(rpcDetId.region() != 0) continue;

  // Loop through the digi collection in the specific chamber.
  for(auto digi = (*chamber).second.first; digi != (*chamber).second.second; ++digi) {
    //std::cout << "itr = " << itr;
    // (*hit.second.first) is our digi iterator and (*hit.second.second) is the ending one.
    //bx_v.push_back(digi->bx());
    //strip_v.push_back(digi->strip());

    // Select hits with |bx| <= 3 only.
    if(fabs(digi->bx()) > 3) continue;
    /// Create cluster ids and store their size
    // Check:
    // if the two successive digis have the same bx and in adjacent strips, fill in the same cluster.
    // if any of the two conditions is false, construct a new cluster (cluster_id and cluster_size).
    if(abs(digi->strip() - strip_n1) != 1 || digi->bx() != bx_n1) {
      // Fill the cluster size for the previous cluster in vcluster_size before
      // assigning zero for cluster_size and adding a new \index for the next cluster.
      //std::cout << "itr = " << itr;
      if(itr != 0) {
        //std::cout << "Size: " << cluster_size;
        //if(cluster_size == 0) std::cout << "\n\n\n\n\n WARNING: WE HAVE A ZERO \n\n\n\n\n\n\n\n";
        vcluster_size.push_back(cluster_size); // note: the cluster_id = index for the specific cluster in the vector.
        //hist_clusterSize_RPCTwinMux->Fill(cluster_size);
        clusterSize_bx[cluster_size].push_back(bx_n1);
      }

      //std::cout << "Cluster size = " << vcluster_size[cluster_id] << std::endl;
      cluster_size = 0; // for the new cluster.
      cluster_id++; // assign a new index for the next cluster.
    }
    itr++;
    cluster_size++;

    /// Hits belong to cluster with cluster_id.
    // Store the info of the hit in a tmp variable.
    RPCHitCleaner::detId_Ext tmp{rpcDetId, digi->bx(), digi->strip()};
    // Assign the cluster index to this hit.
    hits[tmp] = cluster_id;
    /// Strip of i-1
    strip_n1 = digi->strip();
    bx_n1 = digi->bx();
  } // End of loop over digis.
} // End of first loop over chambers.
std::cout << "Size: " << cluster_size;
if(cluster_size == 0) std::cout << "\n\n\n\n\n WARNING: WE HAVE A ZERO \n\n\n\n\n\n\n\n";
vcluster_size.push_back(cluster_size); // store size of the last cluster.
// If the event does not have any clusters thats satisfy the conditions, we have cluster_size = 0, the default value.
//hist_clusterSize_RPCTwinMux->Fill(cluster_size);
clusterSize_bx[cluster_size].push_back(bx_n1);
//If the event does not have any clusters thats satisfy the conditions, we have cluster_size = 0, the default value,
// and then the map stores the default value of bx_n1 which is -10000.
//std::cout << "Final Cluster size = " << vcluster_size[cluster_id] << std::endl;

/*/// Another way to form the vcluster_size vector.
std::map<RPCDetId, std::map<RPCHitCleaner::detId_Ext, int>> clusters_all;
//std::map<RPCHitCleaner::detId_Ext, int> clusters;
//int cluster_id = 1;
//int cluster_size = 0;
int itrr = 1;
for(auto chamber = m_inrpcDigis->begin(); chamber != m_inrpcDigis->end(); ++chamber) {
  RPCDetId rpcDetId = (*chamber).first;
  //int strip_n1 = -10000;
  //int bx_n1 = -10000;
  if(rpcDetId.region() != 0) continue;
  for(auto digi = (*chamber).second.first; digi != (*chamber).second.second; ++digi) {
    int strip_n2 = digi->strip();
    int bx_n2 = digi->bx();
    if(fabs(bx_n2) > 3) continue;
    RPCHitCleaner::detId_Ext tmp{rpcDetId, bx_n2, strip_n2};
    //RPCHitCleaner detId_Ext tmp_2(rpcDetId, bx_n1, strip_n1 - 1); // strip to the left
    //RPCHitCleaner detId_Ext tmp_3(rpcDetId, bx_n1, strip_n1 + 1); // strip to the right
    if(itrr == 1) {
      clusters[tmp]++;
      continue;
    }
    // Check if we have recorded the cluster before.
    clusters = clusters_all[rpcDetId]
    for(auto ext = clusters.begin(); ext != clusters.end(); ++ext) {
      ext_strip = ext->first.strip;
      ext_bx = ext->first.bx;
      //strip_before = (ext->first).strip - 1;
      //strip_after = (ext->first).strip + 1;
      //RPCHitCleaner::detId_Ext ext_before{rpcDetId, bx_n1, strip_before}; // strip to the left
      //RPCHitCleaner::detId_Ext ext_after{rpcDetId, bx_n1, strip_after}; // strip to the right
      if(abs(strip_n2 - ext_strip) != 1 || bx_n2 != ext_bx)
      if(tmp == ext->first || tmp == ext_before || tmp == ext_after) clusters[ext->first]++;
      else clusters[tmp]++;
    }
    itr++;
  }
    //cluster_id++;
    //cluster_size = 1;
  }
    /*for(auto digi2 = (*chamber).second.first; digi2 != (*chamber).second.second; ++digi2) {
      int strip_n2 = digi2->strip();
      int bx_n2 = digi2->bx();
      if(fabs(bx_n2) > 3) continue;
      RPCHitCleaner detId_Ext tmp2(rpcDetId, bx_n2, strip_n2);
      RPCHitCleaner detId_Ext tmp3(rpcDetId, bx_n2, strip_n2 - 1);
      RPCHitCleaner detId_Ext tmp4(rpcDetId, bx_n2, strip_n2 + 1);
      if(tmp == tmp3 || tmp == tmp4) {
        cluster_size ++;
      }
    }*/
//  }
//}





////////////////////////////////////// End of second way.*/



// Loop over chambers and store the min bx of a digi in each chamber's digi collection.
// Filling bx_hits.
// Second Loop through the chmabers. RPCDetId specifies a chamber.
for(auto chamber = m_inrpcDigis->begin(); chamber != m_inrpcDigis->end(); ++chamber) {
  RPCDetId rpcDetId = (*chamber).first;
  // Select Barrel hits only.
  if(rpcDetId.region() != 0) continue;
  BXToStrips strips; // ???
  int cluster_n1 = -10;
  bx_hits[rpcDetId] = 10; // a start value

  /// Keep cluster with min bx in a roll.
  // Fitst Inner Loop through the digi collection in the specific chamber.
  for(auto digi = (*chamber).second.first; digi != (*chamber).second.second; ++digi) {
    // Same as
    //if(fabs(digi->bx()) > 3) continue;
    if(BXToStrips::outOfRange(digi->bx())) continue;
    // Store the info of the hit in a tmp variable.
    RPCHitCleaner::detId_Ext tmp{rpcDetId, digi->bx(), digi->strip()};
    // Get the cluster_id of this hit.
    int cluster_id = hits[tmp];
    /// Remove clusters with size>=4
    // if a chamber have all its clusters with size >= 4 then we may have chambers
    // in bx_hits with bx = 10 (start value) or it does not matter because we will
    // neglect digis which belong to clusters of size >=4?
    if(vcluster_size[cluster_id] >= 4) continue;
    // keep cluster with min bx in a roll.
    //if(bx_hits[rpcDetId] > digi->bx())
    if(digi->bx() < bx_hits[rpcDetId])
      bx_hits[rpcDetId] = digi->bx();
  } // End of the first inner loop over digis.

  // Second Inner Loop through the digi collection in the specific chamber.
  for(auto digi = (*chamber).second.first; digi != (*chamber).second.second; ++digi) {
    if(fabs(digi->bx()) > 3) continue;
    // Store the info of the hit in a tmp variable.
    RPCHitCleaner::detId_Ext tmp{rpcDetId, digi->bx(), digi->strip()};
    int cluster_id = hits[tmp];
    /// Remove clusters with size>=4
    if(vcluster_size[cluster_id] >= 4) continue;
    /// Keep only one bx per st/sec/wheel/layer (chamber ?)
    // Keep only the min bx we have constructed in "bx_hits"?
    if(digi->bx() != bx_hits[rpcDetId]) continue;
    /// Count strips in a cluster
    if(cluster_n1 != cluster_id) {
      strips[digi->bx()] = {0}; // restart counter
    }
    strips[digi->bx()]++;
    cluster_n1 = cluster_id;  // assign the current cluster_id to cluster_n1

    // If the cluster_size of the current cluster = 2
    // and
    // the #strips in the bx != 2
    // neglect this digi
    // That means we store all digis which belong to clusters of size < 3,
    // i.e., (1, 2). And for digis which belong to clusters of size = 3, we
    // store only the 2nd digi in the cluster.
    if(vcluster_size[cluster_id] == 3 && strips[digi->bx()] != 2) continue;

    ///Keep clusters with size=2. Calculate and store the mean phi in RPCtoDTTranslator
    RPCDigi digi_out(digi->strip(), digi->bx());
    m_outrpcDigis.insertDigi(rpcDetId, digi_out);

  } // End of the second inner loop over digis.


  /*
  // Print the bx values stored in bx_hits.
  std::cout << "\n\n\nbx = { ";
  for(auto bxHits_it : bx_hits) {
    std::cout << bxHits_it.second << ", ";
  }
  std::cout << "}; \n";*/

} // End of the second loop over chambers.


// Loop through the vcluster_size vector to fill the cluster size for RPCTwinMux clusters.
for(int clu_size : vcluster_size){
  if(clu_size == 0) continue; // To neglect ignored events.
  hist_clusterSize_RPCTwinMux_before->Fill(clu_size);
//  std::cout << clu_size << std::endl;
}
//std::cout << "\n";

// Plot the count of cluster_size(s) after cleaning.
for (auto chamber = m_outrpcDigis.begin(); chamber != m_outrpcDigis.end(); ++chamber) {
  RPCDetId rpcDetId = (*chamber).first;
  // Loop over the digis.
  for(auto digi = (*chamber).second.first; digi != (*chamber).second.second; ++digi) {
    // Store the digi info in tmp.
    RPCHitCleaner::detId_Ext tmp{rpcDetId, digi->bx(), digi->strip()};
    // Get the cluster_id of the digi.
    int cluster_id = hits[tmp];
    // Get the size of the cluster.
    int clu_size = vcluster_size[cluster_id];
    // Fill the cluster size in the histogram.
    if(clu_size == 2)
      hist_clusterSize_RPCTwinMux_after->Fill(clu_size, 0.5);
    else
      hist_clusterSize_RPCTwinMux_after->Fill(clu_size);
  }
}

// Draw the scatter plot of cluster_size vs bx before cleaning.
for(auto cS_bx : clusterSize_bx) {
  if(cS_bx.first == 0) continue;
  for(int BX : cS_bx.second){
    hist_clusterSize_bx->Fill(cS_bx.first, BX);
  }
}

// Draw the scatter plot of cluster_size vs bx after cleaning.
for(auto chamber = m_outrpcDigis.begin(); chamber != m_outrpcDigis.end(); ++chamber) {
  RPCDetId rpcDetId = (*chamber).first;
  for(auto digi = (*chamber).second.first; digi != (*chamber).second.second; ++digi) {
    // Store the digi info in tmp.
    RPCHitCleaner::detId_Ext tmp{rpcDetId, digi->bx(), digi->strip()};
    // Get the cluster_id of the digi.
    int cluster_id = hits[tmp];
    // Get the size of the cluster.
    int clu_size = vcluster_size[cluster_id];
    // Get the bx of the cluster
    int clu_bx = digi->bx();
    // Fill the histogram.
    if(clu_size == 2)
    hist_clusterSize_bx_after->Fill(clu_size, clu_bx, 0.5);
    else
    hist_clusterSize_bx_after->Fill(clu_size, clu_bx);
  }
}



#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void RPC2TMAna::beginJob() {
  // please remove this method if not needed

  // Access the TFileService object.
  edm::Service<TFileService> fs;

  // Book Histograms.
  hist_phiInSize = fs->make<TH1D>("phiIn_size", "phiIn_size", 60, -0.5, 59.5);
  hist_phiOutSize = fs->make<TH1D>("phiOut_size", "phiOut_size", 60, -0.5, 59.5);
  hist_thetaSize = fs->make<TH1D>("theta_size", "theta_size", 50, -0.5, 49.5);
  hist_clusterSize_RPCTwinMux_before = fs->make<TH1D>("clusterSize_RPCTwinMux_before", "clusterSize_RPCTwinMux_before", 60, -0.5, 59.5);
  hist_clusterSize_RPCTwinMux_after = fs->make<TH1D>("clusterSize_RPCTwinMux_after", "clusterSize_RPCTwinMux_after", 5, -0.5, 4.5);
  hist_clusterSize_bx = fs->make<TH2D>("clusterSize_bx", "clusterSize_bx", 60, -0.5, 59.5, 9, -4.5, 4.5);
  hist_clusterSize_bx_after = fs->make<TH2D>("clusterSize_bx_after", "clusterSize_bx_after", 60, -0.5, 59.5, 9, -4.5, 4.5);



}

// ------------ method called once each job just after ending the event loop  ------------
void RPC2TMAna::endJob() {
  // please remove this method if not needed

  // Check overflows
  //std::cout << "\n\n\n\n";
  //std::cout << "#overflows in phi_in = " << hist_phiInSize->GetBinContent(hist_phiInSize->GetNbinsX() + 1) << std::endl;
  //std::cout << "#overflows in phi_out = " << hist_phiOutSize->GetBinContent(hist_phiOutSize->GetNbinsX() + 1) << std::endl;
  //std::cout << "#overflows in theta = " << hist_thetaSize->GetBinContent(hist_thetaSize->GetNbinsX() + 1) << std::endl;
  //std::cout << "#overflows in hist_clusterSize_RPCTwinMux = " << hist_clusterSize_RPCTwinMux->GetBinContent(hist_clusterSize_RPCTwinMux->GetNbinsX() + 1) << std::endl;

  /*
  // Print the region vector.
  std::cout << "\n\n\nregion = { ";
  for(int n : region_v){
    std::cout << n << ", ";
  }
  std::cout << "}; \n";

  // Print the bx vector.
  std::cout << "\n\n\nbx = { ";
  for(int n : bx_v){
    std::cout << n << ", ";
  }
  std::cout << "}; \n";

  // Print the strip vector.
  std::cout << "\n\n\nstrip = { ";
  for(int n : strip_v){
    std::cout << n << ", ";
  }
  std::cout << "}; \n";*/

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void RPC2TMAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RPC2TMAna);
