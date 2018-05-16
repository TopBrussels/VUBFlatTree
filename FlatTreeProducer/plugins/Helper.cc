#include "VUBFlatTree/FlatTreeProducer/interface/Helper.hh"
#include <iostream>
#include "TMath.h"
#include "TLorentzVector.h"

namespace
{
   struct ByEta
     {
        bool operator()(const pat::PackedCandidate *c1, const pat::PackedCandidate *c2) const
	  {
	     return c1->eta() < c2->eta();
	  }
        bool operator()(float c1eta, const pat::PackedCandidate *c2) const
	  {
	     return c1eta < c2->eta();
	  }
        bool operator()(const pat::PackedCandidate *c1, float c2eta) const
	  {
	     return c1->eta() < c2eta;
	  }
     };
}

float GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
   float DeltaPhi = TMath::Abs(phi2 - phi1);
   if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}

double getEA(double eta, const std::vector<double>& EA_ETA, const std::vector<double>& EA_VALUE)
{
   if (EA_VALUE.size() != EA_ETA.size()+1)
     throw std::invalid_argument("Isolation::getEA : EA_VALUE.size() != EA_ETA.size()+1");
   double EA = EA_VALUE.back();
   for (unsigned i=0; i<EA_ETA.size(); ++i)
     {
	if (fabs(eta) < EA_ETA.at(i))
	  {
	     EA = EA_VALUE.at(i);
	     break;
	  }
     }
   return EA;
}

// part of miniIso for ttH
float isoSumRaw(const std::vector<const pat::PackedCandidate *> & cands, const reco::Candidate &cand, float dR, float innerR, float threshold, SelfVetoPolicy::SelfVetoPolicy selfVeto, int pdgId)
{
   std::vector<const reco::Candidate *> vetos; vetos.clear();
   
   float dR2 = dR*dR, innerR2 = innerR*innerR;
   
   for( unsigned int i=0,n=cand.numberOfSourceCandidatePtrs();i<n;++i )
     {
        if(selfVeto == SelfVetoPolicy::selfVetoNone) break;
        const reco::CandidatePtr &cp = cand.sourceCandidatePtr(i);
        if( cp.isNonnull() && cp.isAvailable() )
	  {
	     vetos.push_back(&*cp);
	     if (selfVeto == SelfVetoPolicy::selfVetoFirst) break;
	  }
     }
   
   typedef std::vector<const pat::PackedCandidate *>::const_iterator IT;
   IT candsbegin = std::lower_bound(cands.begin(), cands.end(), cand.eta() - dR, ByEta());
   IT candsend = std::upper_bound(candsbegin, cands.end(), cand.eta() + dR, ByEta());
   
   double isosum = 0;
   for( IT icharged=candsbegin;icharged<candsend;++icharged )
     {
	//	if( ! (*icharged)->hasTrackDetails() ) continue;
	
        // pdgId
        if( pdgId > 0 && abs((*icharged)->pdgId()) != pdgId ) continue;
        // threshold
        if( threshold > 0 && (*icharged)->pt() < threshold ) continue;
        // cone
        float mydr2 = reco::deltaR2(**icharged, cand);
        if( mydr2 > dR2 || mydr2 < innerR2 ) continue;
        // veto
        if( std::find(vetos.begin(), vetos.end(), *icharged) != vetos.end() )
	  {
	     continue;
	  }
        // add to sum
        isosum += (*icharged)->pt();
     }
   return isosum;
}

// #############
// # Electrons #
// #############
//
float ElecPfIsoCharged(const pat::Electron& elec,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR)
{
   std::vector<const pat::PackedCandidate *> charged;
   
   for( const pat::PackedCandidate &p : *pfcands )
     {
        if( p.charge() != 0 )
	  {
	     if( fabs(p.pdgId()) == 211 )
	       {
		  //		  if( p.hasTrackDetails() )
		  //		    {		       
		  if (p.fromPV() > 1 && fabs(p.dz()) < 9999. )
		    {
		       charged.push_back(&p);
		    }
//		    }		  
	       }
	  }
     }

   std::sort(charged.begin(), charged.end(), ByEta());
   
   float innerR_Ch = .0;
   if( elec.isEB() ) { innerR_Ch = 0.0; }
   else { innerR_Ch = 0.015; }
   
   float result = isoSumRaw(charged,elec,miniIsoR,innerR_Ch,0.0,SelfVetoPolicy::selfVetoNone);
   
   return result;
}

float ElecPfIsoNeutral(const pat::Electron& elec,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR)
{
   std::vector<const pat::PackedCandidate *> neutral;

   for( const pat::PackedCandidate &p : *pfcands )
     {
        if( p.charge() == 0 )
	  {
	     neutral.push_back(&p);
	  }
     }
   
   std::sort(neutral.begin(), neutral.end(), ByEta());
   
   float innerR_N = .0;
   if( elec.isEB() ) { innerR_N = 0.0; }
   else { innerR_N = 0.08; }
   
   float result1 = isoSumRaw(neutral,elec,miniIsoR,innerR_N,0.0,SelfVetoPolicy::selfVetoNone,22 );
   float result2 = isoSumRaw(neutral,elec,miniIsoR,0.0     ,0.0,SelfVetoPolicy::selfVetoNone,130);
   float result = result1 + result2;

   return result;
}

// #########
// # Muons #
// #########
//
float MuonPfIsoCharged(const pat::Muon& muon,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR)
{   
   std::vector<const pat::PackedCandidate *> charged;
   
   for( const pat::PackedCandidate &p : *pfcands )
     {
        if( p.charge() != 0 )
	  {
	     if( fabs(p.pdgId()) == 211 )
	       {
		  //		  if( p.hasTrackDetails() )
		  //		    {		       
		  if (p.fromPV() > 1 && fabs(p.dz()) < 9999. )
		    {
		       charged.push_back(&p);
		    }
		  //		    }		  
	       }
	  }
     }
   
   std::sort(charged.begin(), charged.end(), ByEta());
   
   return isoSumRaw(charged,muon,miniIsoR,0.0001,0.0,SelfVetoPolicy::selfVetoAll);
}

float MuonPfIsoNeutral(const pat::Muon& muon,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR)
{
   std::vector<const pat::PackedCandidate *> neutral;
   
   for( const pat::PackedCandidate &p : *pfcands )
     {
        if( p.charge() == 0 )
	  {
	     neutral.push_back(&p);
	  }
     }
   
   std::sort(neutral.begin(), neutral.end(), ByEta());
   
   return isoSumRaw(neutral,muon,miniIsoR,0.01,0.5,SelfVetoPolicy::selfVetoAll);
}

// https://twiki.cern.ch/twiki/bin/view/CMS/MiniIsolationSUSY
// https://github.com/manuelfs/CfANtupler/blob/master/minicfa/interface/miniAdHocNTupler.h
double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
		      const reco::Candidate* ptcl,
		      double rho, double EA,
		      double r_iso_min, double r_iso_max, double kt_scale,
		      bool use_pfweight, bool charged_only)
{
   if (ptcl->pt()<5.) return 99999.;
   double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
   if(ptcl->isElectron())
     {
	if (fabs(ptcl->eta())>1.479)
	  {
	     deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;
	  }
     }
   else if(ptcl->isMuon())
     {
	deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
     }
   else
     {
	//deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
     }

   double iso_nh(0.); double iso_ch(0.);
   double iso_ph(0.); double iso_pu(0.);
   double ptThresh(0.5);
   if(ptcl->isElectron()) ptThresh = 0;
   double r_iso = std::max(r_iso_min,std::min(r_iso_max, kt_scale/ptcl->pt()));
   for (const pat::PackedCandidate &pfc : *pfcands)
     {
	if (abs(pfc.pdgId())<7) continue;
	double dr = deltaR(pfc, *ptcl);
	if (dr > r_iso) continue;
	////////////////// NEUTRALS /////////////////////////
	if (pfc.charge()==0)
	  {
	     if (pfc.pt()>ptThresh)
	       {
		  double wpf(1.);
		  if (use_pfweight)
		    {
		       double wpv(0.), wpu(0.);
		       for (const pat::PackedCandidate &jpfc : *pfcands)
			 {
			    double jdr = deltaR(pfc, jpfc);
			    if (pfc.charge()!=0 || jdr<0.00001) continue;
			    double jpt = jpfc.pt();
			    if ( pfc.hasTrackDetails() && pfc.fromPV()>1) wpv *= jpt/jdr;
			    else wpu *= jpt/jdr;
			 }

		       wpv = log(wpv);
		       wpu = log(wpu);
		       wpf = wpv/(wpv+wpu);
		    }

		  /////////// PHOTONS ////////////
		  if (abs(pfc.pdgId())==22)
		    {
		       if(dr < deadcone_ph) continue;
		       iso_ph += wpf*pfc.pt();
		       /////////// NEUTRAL HADRONS ////////////
		    }
		  else if (abs(pfc.pdgId())==130)
		    {
		       if(dr < deadcone_nh) continue;
		       iso_nh += wpf*pfc.pt();
		    }
	       }
	     ////////////////// CHARGED from PV /////////////////////////
	  }
	else if (pfc.hasTrackDetails() && pfc.fromPV()>1)
	  {
	     if (abs(pfc.pdgId())==211)
	       {
		  if(dr < deadcone_ch) continue;
		  iso_ch += pfc.pt();
	       }

	     ////////////////// CHARGED from PU /////////////////////////
	  }
	else
	  {
	     if (pfc.hasTrackDetails() && pfc.pt()>ptThresh)
	       {
		  if(dr < deadcone_pu) continue;
		  iso_pu += pfc.pt();
	       }
	  }
     }

   double iso(0.);
   if (charged_only)
     {
	iso = iso_ch;
     }
   else
     {
	//iso = iso_ph + iso_nh;
	// Formula for delta beta correction
	//if (!use_pfweight) 	iso = std::max(iso_ph+iso_nh - 0.5*iso_pu,0);
	// Formula for EA correction
	//iso += rho*EA;
	iso = iso_ch;
	double dr = 0;
	if(ptcl->pt()<50) dr = 0.2;
	if(ptcl->pt()>=50 && ptcl->pt()<200) dr = 10./ptcl->pt();
	if(ptcl->pt()>=200) dr = 0.05;
	iso += std::max(iso_ph+iso_nh - rho*EA*pow((dr/0.3),2),0.);

     }

   iso = iso/ptcl->pt();

   return iso;
}

double ptRatioElec(const pat::Electron& elec, const pat::Jet& jet)
{
   pat::Jet myCorJet;
   myCorJet.setP4(jet.correctedJet("L1FastJet").p4());

   float           SF          = jet.p4().E() / myCorJet.p4().E();

   auto lepAwareJetp4 = ( myCorJet.p4() - elec.p4() ) * SF + elec.p4();
   
   float ptRatio = elec.pt() / lepAwareJetp4.pt();
   
   //std::cout << "Electron ptRatio:            " << ptRatio                      << std::endl;
   
   return (ptRatio > 0) ? ptRatio : 0.0;
}

float ptRelElec(const pat::Electron& elec,const pat::Jet& jet)
{
   pat::Jet myCorJet;
   myCorJet.setP4(jet.correctedJet("L1FastJet").p4());

   float          SF          = jet.p4().E() / myCorJet.p4().E();

   auto lepAwareJetp4 = ( myCorJet.p4() - elec.p4() ) * SF + elec.p4();

   TLorentzVector elecV = TLorentzVector(elec.px(),elec.py(),elec.pz(),elec.p4().E());
   TLorentzVector jetV = TLorentzVector(lepAwareJetp4.px(),lepAwareJetp4.py(),lepAwareJetp4.pz(),lepAwareJetp4.E());

   float PtRel = elecV.Perp( (jetV - elecV).Vect() );

   //std::cout << "Electron PtRel (lep aware):   " << PtRel                      << std::endl;
   
   return (PtRel > 0) ? PtRel : 0.0;
}

float conePtElec(const pat::Electron& elec,const pat::Jet& jet)
{
   pat::Jet myCorJet;
   myCorJet.setP4(jet.correctedJet("L1FastJet").p4());

   float SF = jet.p4().E() / myCorJet.p4().E();
   
   auto lepAwareJetp4 = ( myCorJet.p4() - elec.p4() ) * SF + elec.p4();
   
   return lepAwareJetp4.pt();
}

double ptRatioMuon(const pat::Muon& muon,const pat::Jet& jet)
{
   pat::Jet myCorJet;
   myCorJet.setP4(jet.correctedJet("L1FastJet").p4());

   float          SF          = jet.p4().E() / myCorJet.p4().E();

   auto lepAwareJetp4 = ( myCorJet.p4() - muon.p4() ) * SF + muon.p4();

   float ptRatio = muon.pt() / lepAwareJetp4.pt();

   //std::cout << "Muon ptRatio (lep aware):     " << ptRatio << std::endl;

   return (ptRatio > 0) ? ptRatio : 0.0;
}

float ptRelMuon(const pat::Muon& muon,const pat::Jet& jet)
{
   pat::Jet myCorJet;
   myCorJet.setP4(jet.correctedJet("L1FastJet").p4());

   float          SF          = jet.p4().E() / myCorJet.p4().E();

   auto lepAwareJetp4 = ( myCorJet.p4() - muon.p4() ) * SF + muon.p4();

   TLorentzVector muonV = TLorentzVector(muon.px(),muon.py(),muon.pz(),muon.p4().E());
   TLorentzVector jetV = TLorentzVector(lepAwareJetp4.px(),lepAwareJetp4.py(),lepAwareJetp4.pz(),lepAwareJetp4.E());
   
   float PtRel = muonV.Perp( (jetV - muonV).Vect() );
   
   //std::cout << "Muon PtRel (lep aware):   " << PtRel                      << std::endl;
   
   return (PtRel > 0) ? PtRel : 0.0;
}

float conePtMuon(const pat::Muon& muon,const pat::Jet& jet)
{
   pat::Jet myCorJet;
   myCorJet.setP4(jet.correctedJet("L1FastJet").p4());
   
   float SF = jet.p4().E() / myCorJet.p4().E();

   auto lepAwareJetp4 = ( myCorJet.p4() - muon.p4() ) * SF + muon.p4();
   
   return lepAwareJetp4.pt();
}

int jetNDauChargedMVASel(const pat::Jet& jet, const reco::Vertex& vtx)
{
   int n = 0;
  
   /*   std::vector<reco::PFCandidatePtr> pfConsts = jet.getPFConstituents();
   for( std::vector<reco::PFCandidatePtr>::const_iterator pfJetConstituent=pfConsts.begin();
	pfJetConstituent!=pfConsts.end();
	++pfJetConstituent )
     {
	const reco::Candidate* icand = pfJetConstituent->get();	
	const pat::PackedCandidate* x = dynamic_cast<const pat::PackedCandidate* >( icand );
    if ( GetDeltaR(x->eta(),x->phi(),jet.eta(),jet.phi()) <= 0.4 && x->charge() != 0 && x->fromPV() > 1 &&
	     qualityTrk(x->pseudoTrack(),vtx) )
	  n++;
     }*/

   if( vtx.isValid() )
     {	
	for( unsigned int id=0,nd=jet.numberOfDaughters();id<nd;++id )
	  {
	     const pat::PackedCandidate *x = dynamic_cast<const pat::PackedCandidate* >( jet.daughter(id) );
	     if( x->hasTrackDetails() )
	       {
		  if ( GetDeltaR(x->eta(),x->phi(),jet.eta(),jet.phi()) <= 0.4 && x->charge() != 0 && x->fromPV() > 1 &&
		       qualityTrk(x->pseudoTrack(),vtx) )
		    n++;
	       }	     
	  }   
     }
   
  return n;
}

bool qualityTrk(const reco::Track trk, const reco::Vertex &vtx) 
{
   bool isgoodtrk = false;
   if(trk.pt()>1 &&
      trk.hitPattern().numberOfValidHits()>=8 &&
      trk.hitPattern().numberOfValidPixelHits()>=2 &&
      trk.normalizedChi2()<5 &&
      std::fabs(trk.dxy(vtx.position()))<0.2 &&
      std::fabs(trk.dz(vtx.position()))<17
     ) isgoodtrk = true;
   return isgoodtrk;
}
