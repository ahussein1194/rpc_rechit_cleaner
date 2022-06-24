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


// Taken from RPCHitCleaner.cc

//RPCDigiCollection m_inrpcDigis = digiCollectionRPCTwinMux;
typedef  DigiContainerIterator<RPCDetId, RPCDigi> DigiRangeIterator;
DigiRangeIterator m_inrpcDigis = digiCollectionRPCTwinMux;
//std::cout << *digiCollectionRPCTwinMux.first << std::endl;
for(hit = m_inrpcDigis->begin(); hit != m_inrpcDigis->end(); ++hit) {
  RPCDetId rpcDetId = (*hit).first;
  std::cout << "Region is: " << rpcDetId.region() << std::endl;
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
  hist_phiInSize = fs->make<TH1D>("phiIn_size", "phiIn_size", 60, -0.5, 59.5);  // integer
  hist_phiOutSize = fs->make<TH1D>("phiOut_size", "phiOut_size", 60, -0.5, 59.5);  // integer
  hist_thetaSize = fs->make<TH1D>("theta_size", "theta_size", 60, -0.5, 59.5);  // integer
}

// ------------ method called once each job just after ending the event loop  ------------
void RPC2TMAna::endJob() {
  // please remove this method if not needed

  // Check overflows
  //std::cout << "\n\n\n\n";
  //std::cout << "#overflows in phi_in = " << hist_phiInSize->GetBinContent(hist_phiInSize->GetNbinsX() + 1) << std::endl;
  //std::cout << "#overflows in phi_out = " << hist_phiOutSize->GetBinContent(hist_phiOutSize->GetNbinsX() + 1) << std::endl;
  //std::cout << "#overflows in theta = " << hist_thetaSize->GetBinContent(hist_thetaSize->GetNbinsX() + 1) << std::endl;

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
