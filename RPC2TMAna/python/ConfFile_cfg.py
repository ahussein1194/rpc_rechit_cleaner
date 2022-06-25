import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.GlobalTag.globaltag = "123X_dataRun3_Express_v5"		#"122X_dataRun3_Express_v3"

### RPC RawToDigi - from Legacy
process.load("EventFilter.RPCRawToDigi.rpcUnpackingModule_cfi")

### RPC RawToDigi - from TwinMux
process.load("EventFilter.RPCRawToDigi.RPCTwinMuxRawToDigi_cff")

### RPC RawToDigi - from CPPF
#process.load("EventFilter.RPCRawToDigi.RPCCPPFRawToDigi_cff")
#process.rpcCPPFRawToDigi = cms.EDProducer("RPCCPPFUnpacker",
#  inputLabel = cms.InputTag('rawDataCollector'),
#)

# process.load("EventFilter.RPCRawToDigi.RPCCPPFRawToDigi_sqlite_cff") #to load CPPF link maps from the local DB

### RPC RawToDigi - from OMTF
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.omtfStage2Digis = cms.EDProducer("OmtfUnpacker",
  inputLabel = cms.InputTag('rawDataCollector'),
)

#from EventFilter.RPCRawToDigi.rpcTwinMuxRawToDigi_cfi import rpcTwinMuxRawToDigi
#from EventFilter.RPCRawToDigi.RPCCPPFRawToDigi_cfi import rpcCPPFRawToDigi
#from EventFilter.L1TRawToDigi.bmtfDigis_cfi import bmtfDigis
#from EventFilter.L1TRawToDigi.omtfStage2Digis_cfi import omtfStage2Digis
#from EventFilter.L1TRawToDigi.emtfStage2Digis_cfi import emtfStage2Digis
#from EventFilter.L1TRawToDigi.caloLayer1Digis_cfi import caloLayer1Digis
#from EventFilter.L1TRawToDigi.caloStage2Digis_cfi import caloStage2Digis
#from EventFilter.L1TRawToDigi.gmtStage2Digis_cfi import gmtStage2Digis
#from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis
#from EventFilter.L1TXRawToDigi.twinMuxStage2Digis_cfi import twinMuxStage2Digis

### RPC RawToDigi - from TwinMux
process.load("EventFilter.L1TXRawToDigi.twinMuxStage2Digis_cfi")



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(6000) )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(

##'file:/afs/cern.ch/work/m/mileva/offlineAna2022/testLegTM/CMSSW_12_2_1/src/testRPCDigiMerger.root'
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/450/00000/f6ee5ebf-a981-4f54-be78-4a33b2605248.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/450/00000/f6f6e541-8818-4301-885e-b28d272870cb.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/450/00000/f7bea01f-3fd3-4c02-9c52-bfe4b6f1783d.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/450/00000/f7c0d09f-e055-45c9-a931-481fb7be0cfd.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/450/00000/f7c1411f-3c46-4e70-83bf-90fd43999b6c.root'


'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/d222973d-caa9-4bb7-ad7f-e4b67a70d632.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/d366768a-64e0-4ecf-a83d-eb06894224bf.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/d3c53ca8-b714-423d-8d88-d45eb11642d8.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/d4e2ee30-72d4-4f1c-a524-c4c36920637f.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/da0653ee-479a-4493-97df-2ff2f4973b47.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/dbee63c0-cddc-49f4-bad2-933558fea72f.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/dc54b89d-fc3d-481c-aa4f-d3ba886061d3.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/de837b6d-f607-41e0-9473-5c0d36c3f468.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/deb42379-e189-49e6-9b0f-eb82f45dee2b.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/decc28f5-2a06-4904-b7d0-c4b681524613.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/e0ca3c49-8406-46b0-9e13-1e1a951ad1c1.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/e0fbdf34-2abf-435a-a962-1467e49321bf.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/e31dacb1-7cc0-40a9-8bca-2e8818b5abf3.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/e3aa3f9c-a5ef-4455-81d3-f85d9b752cd5.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/e4f95802-01bf-4443-b517-eef4a1db9610.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/e5d83d28-9504-4bff-a25a-df3bed4b79ea.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/e6b35a08-5c52-436b-bf8a-687558ddfd61.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/e7187bc5-503d-4b87-9419-eaf8f373769d.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/e7fc990d-9a2e-4e07-9aa2-6973b0e78ee0.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/e952fd17-3a22-4dc2-a0db-fef2c1f69717.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/e99ceb12-48d2-4479-8453-54020d9e4674.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/eb14eb5b-40f0-4eb8-832c-2294bb644a64.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/eb6af3cd-b6bd-4f52-8271-8bc12e0845fb.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/eba95e20-bf3e-4283-9e30-04d9ee9ad4ee.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/edee7764-d5bf-426c-9aeb-d48b01052ca1.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/ee06e90f-05d1-4eab-a50f-da8278dba162.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/f0a3494c-dec8-4d43-9fc4-3431fcf471ea.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/f2616e61-f8a6-443c-a41a-20ca6549f685.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/f377a859-22ae-4fe5-90a5-f9cab4c56146.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/f56754dd-4688-4ae0-8708-98cd4b0cf94b.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/f7d7fc34-31b4-4907-80a4-599260ce84ec.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/f91905eb-f419-4e2f-a885-9bf2c33bf619.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/f99abebe-fff9-42f6-9b0c-26250ded52b2.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/f9f2e7a2-6618-4762-87a6-d44e5d8be3fa.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/fc2401aa-75c0-4496-adfb-66e4e8fae108.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/fc9542ca-4686-4c4f-9e8d-898f9f3acd94.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/fd200061-4291-4954-b555-f3d646917612.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/fd79993d-6b1e-49e6-9ef1-fd192967bd25.root',
'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/566/00000/fdae2f3d-1e77-45fb-a351-cbd5708c4e52.root'

                )
                            )

process.demo  = cms.EDAnalyzer('RPC2TMAna',
   rpcLegacy  = cms.untracked.InputTag('muonRPCDigis'),
   rpcTwinMux = cms.untracked.InputTag('rpcTwinMuxRawToDigi'),
#   inputDCC   = cms.untracked.InputTag( "dttfDigis" ),		#the container is empty in FEVT and data need to be unpaked
   rpccppf = cms.untracked.InputTag('rpcCPPFRawToDigi'),   #End-Cap part
#twinmux unpacker
   inputTagTMphIn = cms.untracked.InputTag('twinMuxStage2Digis:PhIn'),
   inputTagTMphOut = cms.untracked.InputTag('twinMuxStage2Digis:PhOut'),
   inputTagTMth = cms.untracked.InputTag('twinMuxStage2Digis:ThIn')

                              )

# Enable the TFileService.
process.TFileService = cms.Service("TFileService", fileName = cms.string('hist.root'), closeFileFast = cms.untracked.bool(False))

process.p = cms.Path(process.rpcTwinMuxRawToDigi*process.twinMuxStage2Digis*process.demo)

#process.p = cms.Path(process.rpcCPPFRawToDigi*process.rpcTwinMuxRawToDigi*process.twinMuxStage2Digis*process.demo)
