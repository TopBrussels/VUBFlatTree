import FWCore.ParameterSet.Config as cms

genPartons = cms.EDProducer("PartonSelector",
    withLeptons = cms.bool(False),
    src = cms.InputTag("prunedGenParticles")
)

jetPartonMatching = cms.EDProducer("JetPartonMatcher",
    jets    = cms.InputTag("slimmedGenJets"),
    partons = cms.InputTag("genPartons"),
    coneSizeToAssociate = cms.double(0.3),
)

genJetFlavour = cms.EDProducer("JetFlavourIdentifier",
    srcByReference    = cms.InputTag("jetPartonMatching"),
    physicsDefinition = cms.bool(False)
)

genJetMatch = cms.EDProducer("GenJetMatcher",        # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = cms.InputTag("slimmedJets"),       # RECO jets (any View<Jet> is ok)
    matched     = cms.InputTag("slimmedGenJets"),    # GEN jets  (must be GenJetCollection)
    mcPdgId     = cms.vint32(),                      # n/a
    mcStatus    = cms.vint32(),                      # n/a
    checkCharge = cms.bool(False),                   # n/a
    maxDeltaR   = cms.double(0.4),                   # Minimum deltaR for the match
    maxDPtRel   = cms.double(3.0),                   # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),          # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(False),         # False = just match input in order; True = pick lowest deltaR pair first
)

# default PAT sequence for jet flavour identification
genJetFlavourAlg = cms.Sequence(genPartons * jetPartonMatching * genJetFlavour)
