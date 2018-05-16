import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
import os

from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.options.allowUnscheduled = cms.untracked.bool(True)

from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
jetAlgo="AK4"
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=True, postfix=postfix)

options = VarParsing('analysis')
options.register('isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.int,'Run on real data')
options.register('confFile', 'conf.xml', VarParsing.multiplicity.singleton, VarParsing.varType.string, "Flattree variables configuration")
options.register('bufferSize', 32000, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Buffer size for branches of the flat tree")
##process = cms.Process("FlatTree")
options.parseArguments()

#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load('Configuration.StandardSequences.Services_cff')
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('START70_V7::All')

#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/ttH/testFiles/AOD/ttH_ev.root'
    )
)

# run PAT
#process.p = cms.Path( getattr(process,"patPF2PATSequence"+postfix) )

# top projections in PF2PAT:
##getattr(process,"pfNoPileUpJME"+postfix).enable = True
##getattr(process,"pfNoMuonJME"+postfix).enable = True
##getattr(process,"pfNoElectronJME"+postfix).enable = True
##getattr(process,"pfNoTau"+postfix).enable = False
##getattr(process,"pfNoJet"+postfix).enable = True
# to use tau-cleaned jet collection uncomment the following:
#getattr(process,"pfNoTau"+postfix).enable = True
        
# verbose flags for the PF2PAT modules
##getattr(process,"pfNoMuonJME"+postfix).verbose = False
        
# enable delta beta correction for muon selection in PF2PAT?
##getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = cms.bool(False)

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAodDataTable
process.FlatTree = cms.EDAnalyzer('FlatTreeProducer',

                  dataFormat = cms.string("AOD"),

                  bufferSize        = cms.int32(options.bufferSize),
                  confFile          = cms.string(options.confFile),

                  isData            = cms.bool(options.isData),

                  vertexInput       = cms.InputTag("offlinePrimaryVertices"),
                  electronInput     = cms.InputTag("patElectronsPFlow"),
                  muonInput         = cms.InputTag("patMuonsPFlow"),
                  tauInput          = cms.InputTag("patTausPFlow"),
                  jetInput          = cms.InputTag("patJetsPFlow"),
                  metInput          = cms.InputTag("patMETsPFlow"),
                  rhoInput          = cms.InputTag("fixedGridRhoAll"),
                  genParticlesInput = cms.InputTag("genParticles")
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output.root")
)

process.p = cms.Path(process.FlatTree)

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('patTuple.root'),
#                               outputCommands = cms.untracked.vstring(
#                               'keep *_*PFJets*_*_*',
#                               'keep *_*PFlow*_*_*',
#                               'keep *_*pat*_*_*',
#                               'keep *_*_*_PAT'
#                               )
#)

#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
#process.out.outputCommands = cms.untracked.vstring('drop *',
#                                                   'keep *_*PFlow*_*_PAT'
#                                                   )

# drop all info in the output PAT file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out.outputCommands = cms.untracked.vstring('drop *'
                                                   )

#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
#process.out.outputCommands = cms.untracked.vstring('drop *',
#                                                   'keep recoPFCandidates_particleFlow_*_*',
#                                                   'keep *_selectedPatJets*_*_*',
#                                                   'drop *_selectedPatJets*_caloTowers_*',
#                                                   'keep *_selectedPatElectrons*_*_*',
#                                                   'keep *_selectedPatMuons*_*_*',
#                                                   'keep *_selectedPatTaus*_*_*',
#                                                   'keep *_patMETs*_*_*',
#                                                   'keep *_selectedPatPhotons*_*_*',
#                                                   'keep *_selectedPatTaus*_*_*'
#                                                   )

#process.patDefaultSequence = cms.Sequence( process.goodOfflinePrimaryVerticesPFlow
#                                           *getattr(process,"patPF2PATSequence"+postfix) )

#process.p = cms.Path( process.patDefaultSequence )

# only if needed to save pat output
#process.e = cms.EndPath(process.out)
##process.out.fileName = 'patTuple.root'

