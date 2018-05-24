import FWCore.ParameterSet.Config as cms

#####################
#  Options parsing  #
#####################

from FWCore.ParameterSet.VarParsing import VarParsing
import os, sys

options = VarParsing('analysis')
options.register('isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Run on real data')
options.register('applyMETFilters',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply MET filters')
options.register('applyJEC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply JEC corrections')
options.register('runAK10',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Add AK10 jets')

options.register('runQG',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Run QGTagger')

options.register('fillMCScaleWeight',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Fill PDF weights')
options.register('fillPUInfo',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Fill PU info')
options.register('nPDF', -1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "nPDF")
options.register('confFile', 'conf.xml', VarParsing.multiplicity.singleton, VarParsing.varType.string, "Flattree variables configuration")
options.register('bufferSize', 32000, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Buffer size for branches of the flat tree")
options.parseArguments()

##########################
#  Global configuration  #
##########################

process = cms.Process("FlatTree")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.suppressWarning = cms.untracked.vstring(["JetPtMismatchAtLowPt","NullTransverseMomentum"])

# remove verbose from patTrigger due to missing L1 prescales for some trigger paths
#process.MessageLogger.suppressWarning.append('patTrigger')
#process.MessageLogger.cerr.FwkJob.limit=1
#process.MessageLogger.cerr.ERROR = cms.untracked.PSet( limit = cms.untracked.int32(0) )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

if options.isData:
    process.GlobalTag.globaltag = '94X_dataRun2_v6'    
else:
    process.GlobalTag.globaltag = '94X_mc2017_realistic_v13'

corName="Fall17_17Nov2017_V6_MC"
corTag="JetCorrectorParametersCollection_"+corName
if options.isData:
    corName="Fall17_17Nov2017BCDEF_V6_DATA"
    corTag="JetCorrectorParametersCollection_"+corName
dBFile=corName+".db"

process.load("CondCore.CondDB.CondDB_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
                           DBParameters = cms.PSet(
                           messageLevel = cms.untracked.int32(0)
                           ),
                           timetype = cms.string('runnumber'),
                           toGet = cms.VPSet(
                           cms.PSet(
                                    record = cms.string('JetCorrectionsRecord'),
                                    tag    = cms.string(corTag+"_AK4PF"),
                                    label  = cms.untracked.string('AK4PF')
                                    ),
                           cms.PSet(
                                    record = cms.string('JetCorrectionsRecord'),
                                    tag    = cms.string(corTag+"_AK4PFchs"),
                                    label  = cms.untracked.string('AK4PFchs')
                                    ),
                           cms.PSet(
                                    record = cms.string('JetCorrectionsRecord'),
                                    tag    = cms.string(corTag+"_AK8PF"),
                                    label  = cms.untracked.string('AK8PF')
                                    ),
                           cms.PSet(
                                    record = cms.string('JetCorrectionsRecord'),
                                    tag    = cms.string(corTag+"_AK8PFchs"),
                                    label  = cms.untracked.string('AK8PFchs')
                                    ),
                           ),
                           connect = cms.string("sqlite_file:"+dBFile)
)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')                           

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

corList = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
if options.isData:
    corList = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

# Re-apply JEC to AK4
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

#bTagDiscriminators = [
#    'deepFlavourJetTags:probudsg',
#    'deepFlavourJetTags:probb',
#    'deepFlavourJetTags:probc',
#    'deepFlavourJetTags:probbb',
#    'deepFlavourJetTags:probcc',
#]

updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    labelName = 'UpdatedJEC',
    jetCorrections = ('AK4PFchs', corList, 'None')
)

#updateJetCollection(
#    process,
    #jetSource = cms.InputTag('slimmedJets','','PAT'), #FIXME -- jets reclustering to add DeepFlav caused bug -- fixed by Kirill  !
#    jetSource = cms.InputTag('slimmedJets'),
#    jetCorrections = ('AK4PFchs', corList, 'None'),
#    labelName = 'UpdatedJEC'
#    btagDiscriminators = bTagDiscriminators,
#)

# Re-apply JEC to AK8
updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    labelName = 'UpdatedJECAK8',
    jetCorrections = ('AK8PFchs', corList, 'None')
)

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

#jetsNameAK4="selectedUpdatedPatJetsUpdatedJEC"
jetsNameAK4="updatedPatJetsUpdatedJEC"
#jetsNameAK4="slimmedJets"
jetsNameAK8="selectedUpdatedPatJetsUpdatedJECAK8"
#jetsNameAK10="patJetsReapplyJECAK10"
jetsNameAK10="selectedPatJetsAK10PFCHS"

########################
#  Additional modules  #
########################

from VUBFlatTree.FlatTreeProducer.runTauIdMVA import *
na = TauIDEmbedder(process, cms,
     debug=False,
     toKeep = ["dR0p32017v2"]
)
na.runTauID()

if not options.isData:
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    
    from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
    process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
            particles = "prunedGenParticles"
    )
    
    from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
    process.genJetFlavourInfos = ak4JetFlavourInfos.clone(
            jets = "slimmedGenJets"
    )
    
    from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
    process.matchGenBHadron = matchGenBHadron.clone(
            genParticles = "prunedGenParticles",
            jetFlavourInfos = "genJetFlavourInfos"
    )
    
    from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenCHadron
    process.matchGenCHadron = matchGenCHadron.clone(
            genParticles = "prunedGenParticles",
            jetFlavourInfos = "genJetFlavourInfos"
    )
    
# egamma
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)

# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
my_id_modules = [
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff'
]

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer - added unwillingly by Xavier
###process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
#process.calibratedPatElectrons

#from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
#process = regressionWeights(process)

#process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

# MET
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
#process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
#process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
#process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
#process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

#process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
#process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
#process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

#process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
#process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
#process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

#process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
#    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
#    reverseDecision = cms.bool(False)
#)

#process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
#    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
#    reverseDecision = cms.bool(False)
#)

#####################
# MET Significance  #
#####################
process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")
from RecoMET.METProducers.testInputFiles_cff import recoMETtestInputFiles


#######################
# AK10 collection     #
#######################
if options.runAK10:
    from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
    jetToolbox( process, 'ak10', 'ak10JetSubs', 'out', runOnMC=(not options.isData),
                addPruning=True, addSoftDrop=True , addPrunedSubjets=True, addSoftDropSubjets=True,
                JETCorrPayload='AK3Pachs', subJETCorrPayload='AK10PFchs', JETCorrLevels=['L1FastJet', 'L2Relative', 'L3Absolute'],
                addNsub=True, maxTau=6, addTrimming=True, addFiltering=True,
                addEnergyCorrFunc=True, maxECF=5 )    
                
#######################
# Quark gluon tagging #
#######################
if options.runQG:
    qgDatabaseVersion = 'v2b' # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

    from CondCore.DBCommon.CondDBSetup_cfi import *
    QGPoolDBESSource = cms.ESSource("PoolDBESSource",
          CondDBSetup,
          toGet = cms.VPSet(),
          connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
    )

    for type in ['AK4PFchs','AK4PFchs_antib']:
        QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
              record = cms.string('QGLikelihoodRcd'),
              tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
              label  = cms.untracked.string('QGL_'+type)
        )))

    process.load('RecoJets.JetProducers.QGTagger_cfi')
    process.QGTagger.srcJets          = cms.InputTag(jetsNameAK4) # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
    process.QGTagger.jetsLabel        = cms.string('QGL_AK4PFchs') # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

###########
#  Input  #
###########

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"), # WARNING / FIXME for test only !
    fileNames = cms.untracked.vstring(
         #'/store/data/Run2017D/MuonEG/MINIAOD/17Nov2017-v1/50000/3E5F02AC-33E7-E711-AE42-A0369FC5FBA4.root'
#         'file:0CF65340-0200-E811-ABB7-0025905C53F0.root'
#         '/store/mc/RunIIFall17MiniAOD/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/0CF65340-0200-E811-ABB7-0025905C53F0.root'
         '/store/mc/RunIIFall17MiniAOD/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/50000/08DE33C0-87EB-E711-819D-0242AC1C0500.root'
        )
)

############
#  Output  #
############

process.TFileService = cms.Service("TFileService", fileName = cms.string("output.root"))

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True)	 # needed for ak10 computation (JMEAnalysis/JetToolbox)
)

process.slimmedPatTriggerUnpacked = cms.EDProducer('PATTriggerObjectStandAloneUnpacker',
                                                   patTriggerObjectsStandAlone = cms.InputTag('slimmedPatTrigger'),
                                                   triggerResults = cms.InputTag('TriggerResults::HLT'),
                                                   unpackFilterLabels = cms.bool(True)
)

#############################
#  Flat Tree configuration  #
#############################

process.FlatTree = cms.EDAnalyzer('FlatTreeProducer',

                  dataFormat        = cms.string("MINIAOD"),

                  bufferSize        = cms.int32(options.bufferSize),
                  confFile          = cms.string(options.confFile),

                  isData            = cms.bool(options.isData),
                  applyMETFilters   = cms.bool(options.applyMETFilters),
                  fillMCScaleWeight = cms.bool(options.fillMCScaleWeight),
                  fillPUInfo	    = cms.bool(options.fillPUInfo),
                  nPDF              = cms.int32(options.nPDF),
                  
                  vertexInput              = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  electronInput            = cms.InputTag("slimmedElectrons"),
                  electronPATInput         = cms.InputTag("slimmedElectrons"),
                  #electronPATInput         = cms.InputTag("calibratedPatElectrons"),

                  eleVetoCBIdMap           = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-veto"),
                  eleLooseCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose"),
                  eleMediumCBIdMap         = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium"),
                  eleTightCBIdMap          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight"),

                  ele90NoIsoMVAIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
                  ele80NoIsoMVAIdMap         = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
                  eleLooseNoIsoMVAIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose"),

                  ele90IsoMVAIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp90"),
                  ele80IsoMVAIdMap         = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp80"),
                  eleLooseIsoMVAIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wpLoose"),
                  
                  mvaIsoValuesMap             = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
                  mvaNoIsoValuesMap             = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),

#                  mvaCategoriesMap         = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),

#                  BadMuonFilter              = cms.InputTag("BadPFMuonFilter",""),
#                  BadChargedCandidateFilter  = cms.InputTag("BadChargedCandidateFilter",""),
                  
                  filterTriggerNames       = cms.untracked.vstring(
                  "*"
#                  "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*",
#                  "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*",
#                  "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
#                  "HLT_DoubleMu33NoFiltersNoVtx_v*",
#                  "HLT_DoubleMu38NoFiltersNoVtx_v*",
#                  "HLT_Ele22_eta2p1_WPLoose_Gsf_v*",
#                  "HLT_Ele22_eta2p1_WPTight_Gsf_v*",
#                  "HLT_Ele22_eta2p1_WP75_Gsf_v*",
#                  "HLT_Ele23_WPLoose_Gsf_v*",
#                  "HLT_Ele23_WP75_Gsf_v*",
#                  "HLT_Ele27_eta2p1_WP75_Gsf_v*",
#                  "HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v*",
#                  "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
#                  "HLT_Ele27_eta2p1_WPTight_Gsf_v*",
#                  "HLT_Ele32_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v*",
#                  "HLT_Ele32_eta2p1_WPLoose_Gsf_v*",
#                  "HLT_Ele32_eta2p1_WPTight_Gsf_v*",
#                  "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
#                  "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*",
#                  "HLT_IsoMu17_eta2p1_v*",
#                  "HLT_DoubleIsoMu17_eta2p1_v*",
#                  "HLT_IsoMu18_v*",
#                  "HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v*",
#                  "HLT_IsoMu20_v*",
#                  "HLT_IsoMu20_eta2p1_v*",
#                  "HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v*",
#                  "HLT_IsoMu24_eta2p1_v*",
#                  "HLT_IsoMu27_v*",
#                  "HLT_IsoTkMu20_v*",
#                  "HLT_IsoTkMu20_eta2p1_v*",
#                  "HLT_IsoTkMu24_eta2p1_v*",
#                  "HLT_IsoTkMu27_v*",
#                  "HLT_Mu17_Mu8_v*",
#                  "HLT_Mu17_Mu8_DZ_v*",
#                  "HLT_Mu17_Mu8_SameSign_DZ_v*",
#                  "HLT_Mu20_Mu10_v*",
#                  "HLT_Mu20_Mu10_DZ_v*",
#                  "HLT_Mu20_Mu10_SameSign_DZ_v*",
#                  "HLT_Mu17_TkMu8_DZ_v*",
#                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
#                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
#                  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
#                  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
#                  "HLT_Mu27_TkMu8_v*",
#                  "HLT_Mu30_TkMu11_v*",
#                  "HLT_Mu40_TkMu11_v*",
#                  "HLT_Mu20_v*",
#                  "HLT_TkMu20_v*",
#                  "HLT_Mu24_eta2p1_v*",
#                  "HLT_TkMu24_eta2p1_v*",
#                  "HLT_Mu27_v*",
#                  "HLT_TkMu27_v*",
#                  "HLT_Mu50_v*",
#                  "HLT_Mu55_v*",
#                  "HLT_Mu45_eta2p1_v*",
#                  "HLT_Mu50_eta2p1_v*",
#                  "HLT_PFJet40_v*",
#                  "HLT_PFJet60_v*",
#                  "HLT_PFJet80_v*",
#                  "HLT_PFJet140_v*",
#                  "HLT_PFJet200_v*",
#                  "HLT_PFJet260_v*",
#                  "HLT_PFJet320_v*",
#                  "HLT_PFJet400_v*",
#                  "HLT_PFJet450_v*",
#                  "HLT_PFJet500_v*",
#                  "HLT_Mu8_TrkIsoVVL_v*",
#                  "HLT_Mu17_TrkIsoVVL_v*",
#                  "HLT_Mu24_TrkIsoVVL_v*",
#                  "HLT_Mu34_TrkIsoVVL_v*",
#                  "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
#                  "HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
#                  "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
#                  "HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",
#                  "HLT_BTagMu_DiJet20_Mu5_v*",
#                  "HLT_BTagMu_DiJet40_Mu5_v*",
#                  "HLT_BTagMu_DiJet70_Mu5_v*",
#                  "HLT_BTagMu_DiJet110_Mu5_v*",
#                  "HLT_BTagMu_Jet300_Mu5_v*",
#                  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
#                  "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
#                  "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*",
#                  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
#                  "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*",
#                  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
#                  "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
#                  "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*",
#                  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
#                  "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
#                  "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*",
#                  "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
#                  "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
#                  "HLT_Mu8_v*",
#                  "HLT_Mu17_v*",
#                  "HLT_Mu24_v*",
#                  "HLT_Mu34_v*",
#                  "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*",
#                  "HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v*",
#                  "HLT_Ele18_CaloIdM_TrackIdM_PFJet30_v*",
#                  "HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v*",
#                  "HLT_Ele33_CaloIdM_TrackIdM_PFJet30_v*",
#                  "HLT_Mu300_v*",
#                  "HLT_Mu350_v*",
#                  "HLT_PFHT450_SixJet40_PFBTagCSV0p72_v*",
#                  "HLT_PFHT400_SixJet30_BTagCSV0p55_2PFBTagCSV0p72_v*",
#                  "HLT_PFHT450_SixJet40_PFBTagCSV_v*",
#                  "HLT_PFHT400_SixJet30_BTagCSV0p5_2PFBTagCSV_v*"
                  ),
                  
                  muonInput                = cms.InputTag("slimmedMuons"),
                  #tauInput                 = cms.InputTag("slimmedTaus"),
                  tauInput                 = cms.InputTag("NewTauIDsEmbedded"),
                  jetInput                 = cms.InputTag(jetsNameAK4),
                  jetPuppiInput            = cms.InputTag("slimmedJetsPuppi"),
                  ak8jetInput              = cms.InputTag(jetsNameAK8),
                  ak10jetInput             = cms.InputTag(jetsNameAK10),
                  genJetInput              = cms.InputTag("slimmedGenJets"),
                  jetFlavorMatchTokenInput = cms.InputTag("jetFlavourMatch"),
                  metInput                 = cms.InputTag("slimmedMETs"),
                  metPuppiInput            = cms.InputTag("slimmedMETsPuppi"),
                  metNoHFInput             = cms.InputTag("slimmedMETsNoHF"),
                  metSigInput              = cms.InputTag("METSignificance"),
                  metCovInput              = cms.InputTag("METSignificance","METCovariance"),
                  #rhoInput                 = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                  rhoInput                 = cms.InputTag("fixedGridRhoFastjetAll"),
                  genParticlesInput        = cms.InputTag("prunedGenParticles"),
                  genEventInfoInput        = cms.InputTag("generator"),
                  LHEEventProductInput     = cms.InputTag("externalLHEProducer"),
                  bsInput                  = cms.InputTag("offlineBeamSpot"),
                  pfcandsInput             = cms.InputTag("packedPFCandidates"),
                  hConversionsInput        = cms.InputTag("reducedEgamma","reducedConversions"),
                  puInfoInput		   = cms.InputTag("slimmedAddPileupInfo"),
#                  puInfoInput		   = cms.InputTag("addPileupInfo"),
                  objects                  = cms.InputTag("slimmedPatTriggerUnpacked"),
                  
                  genTTXJets                    = cms.InputTag("slimmedGenJets"),
                  genTTXBHadJetIndex            = cms.InputTag("matchGenBHadron","genBHadJetIndex"),
                  genTTXBHadFlavour             = cms.InputTag("matchGenBHadron","genBHadFlavour"),
                  genTTXBHadFromTopWeakDecay    = cms.InputTag("matchGenBHadron","genBHadFromTopWeakDecay"),
                  genTTXBHadPlusMothers         = cms.InputTag("matchGenBHadron","genBHadPlusMothers"),
                  genTTXBHadPlusMothersIndices  = cms.InputTag("matchGenBHadron","genBHadPlusMothersIndices"),
                  genTTXBHadIndex               = cms.InputTag("matchGenBHadron","genBHadIndex"),
                  genTTXBHadLeptonHadronIndex   = cms.InputTag("matchGenBHadron","genBHadLeptonHadronIndex"),
                  genTTXBHadLeptonViaTau        = cms.InputTag("matchGenBHadron","genBHadLeptonViaTau"),
                  genTTXCHadJetIndex            = cms.InputTag("matchGenCHadron","genCHadJetIndex"),
                  genTTXCHadFlavour             = cms.InputTag("matchGenCHadron","genCHadFlavour"),
                  genTTXCHadFromTopWeakDecay    = cms.InputTag("matchGenCHadron","genCHadFromTopWeakDecay"),
                  genTTXCHadBHadronId           = cms.InputTag("matchGenCHadron","genCHadBHadronId")
)

##########
#  Path  #
##########

process.runQG = cms.Sequence()
if options.runQG:
    process.runQG = cms.Sequence(process.QGTagger)

if not options.isData:    
    process.p = cms.Path(
                       #                     process.calibratedPatElectrons+
                       process.rerunMvaIsolationSequence+
                       process.NewTauIDsEmbedded+
                       process.electronMVAValueMapProducer+
                       process.egmGsfElectronIDSequence+
                       #                     process.regressionApplication+
                       process.METSignificance+
                       process.jecSequence+
                       process.runQG+
                       #                     process.BadChargedCandidateFilter+
                       #                     process.BadPFMuonFilter+
                       process.slimmedPatTriggerUnpacked+
                       process.selectedHadronsAndPartons+
                       process.genJetFlavourInfos+
                       process.matchGenBHadron+
                       process.matchGenCHadron+ 
                       process.FlatTree
    )
 
else:
    process.p = cms.Path(
                       #                     process.calibratedPatElectrons+
                       process.rerunMvaIsolationSequence+
                       process.NewTauIDsEmbedded+
                       process.electronMVAValueMapProducer+
                       process.egmGsfElectronIDSequence+
                       #                     process.regressionApplication+
                       process.METSignificance+
                       process.jecSequence+
                       process.runQG+
                       #                     process.BadChargedCandidateFilter+
                       #                     process.BadPFMuonFilter+
                       process.slimmedPatTriggerUnpacked+
                       process.FlatTree
    )
