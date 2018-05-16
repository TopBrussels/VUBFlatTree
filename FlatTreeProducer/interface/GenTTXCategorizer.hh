#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "VUBFlatTree/FlatTreeProducer/interface/FlatTree.hh"

class GenTTXCategorizer
{
 public:
   
   explicit GenTTXCategorizer();
   ~GenTTXCategorizer() {};

   void Init(FlatTree &tree);
   
   void Run(FlatTree &tree,
	    edm::Handle<reco::GenJetCollection> genTTXJets,
	    edm::Handle<std::vector<int> >genTTXBHadFlavour,
	    edm::Handle<std::vector<int> > genTTXBHadJetIndex,
	    edm::Handle<std::vector<int> > genTTXBHadFromTopWeakDecay,
	    edm::Handle<std::vector<reco::GenParticle> > genTTXBHadPlusMothers,
	    edm::Handle<std::vector<std::vector<int> > > genTTXBHadPlusMothersIndices,
	    edm::Handle<std::vector<int> > genTTXBHadIndex,
	    edm::Handle<std::vector<int> > genTTXBHadLeptonHadronIndex,
	    edm::Handle<std::vector<int> > genTTXBHadLeptonViaTau,
	    edm::Handle<std::vector<int> > genTTXCHadFlavour,
	    edm::Handle<std::vector<int> > genTTXCHadJetIndex,
	    edm::Handle<std::vector<int> > genTTXCHadFromTopWeakDecay,
	    edm::Handle<std::vector<int> > genTTXCHadBHadronId);

 private:
   
   double genJetPtMin_;
   double genJetAbsEtaMax_;
};

