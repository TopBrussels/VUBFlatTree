#ifndef MCTRUTH_H
#define MCTRUTH_H

#include <string>
#include <iostream>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "VUBFlatTree/FlatTreeProducer/interface/FlatTree.hh"

#define DEFVAL -666

class MCTruth
{
 public:
   
   MCTruth() {};
   
   void Init(FlatTree &tree);

   void fillGenParticles(const edm::Event& iEvent,
			 const edm::EventSetup& iSetup,
			 FlatTree& tree,
			 const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);

   void fillGenPV(const edm::Event& iEvent,
		  const edm::EventSetup& iSetup,
		  FlatTree& tree,
		  const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);

   void fillStopNeutralinoMass(const edm::Event& iEvent,
		  const edm::EventSetup& iSetup,
		  FlatTree& tree,
		  const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);

   void fillTopStopDecayChain(const edm::Event& iEvent,
			       const edm::EventSetup& iSetup,
			       FlatTree& tree,
			       const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);

   bool doMatch(const edm::Event& iEvent,
		const edm::EventSetup& iSetup,
		const edm::Handle<std::vector<reco::GenParticle> >& GenParticles,
		reco::GenParticle *genp,
		float &drMin,
		float pt, float eta, float phi, int pdgId);
   
   reco::GenParticle* getUnique(const reco::GenParticle* p,
				bool verbose);
   
   void fillTTHSignalGenParticles(const edm::Event& iEvent,
				  const edm::EventSetup& iSetup,
				  FlatTree& tree,
				  const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);

   void fillTTZSignalGenParticles(const edm::Event& iEvent,
				  const edm::EventSetup& iSetup,
				  FlatTree& tree,
				  const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);

   void fillTTWSignalGenParticles(const edm::Event& iEvent,
				  const edm::EventSetup& iSetup,
				  FlatTree& tree,
				  const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);
   
   void fillTZQSignalGenParticles(const edm::Event& iEvent,
				  const edm::EventSetup& iSetup,
				  FlatTree& tree,
				  const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);

   void fillTHQSignalGenParticles(const edm::Event& iEvent,
				  const edm::EventSetup& iSetup,
				  FlatTree& tree,
				  const edm::Handle<std::vector<reco::GenParticle> >& GenParticles);
   
   void p4toTLV(reco::Particle::LorentzVector vp4,
		TLorentzVector& tlv);
   
   const reco::GenParticle* getMother(const reco::GenParticle&);
};

#endif
