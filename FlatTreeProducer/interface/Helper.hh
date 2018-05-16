#ifndef HELPER_H
#define HELPER_H

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/Track.h"
#include <algorithm>

// EA values for electrons
const static std::vector<double> EL_EA_ETA   {1.0,    1.479,    2.0,    2.2,    2.3,    2.4     };
const static std::vector<double> EL_EA_VALUE {0.1752, 0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687};

// EA values for muons
const static std::vector<double> MU_EA_ETA   {0.8,    1.3,    2.0,    2.3     };
const static std::vector<double> MU_EA_VALUE {0.0735, 0.0619, 0.0465, 0.0433, 0.0577 };

namespace SelfVetoPolicy
{
   enum SelfVetoPolicy
     {
	selfVetoNone=0, selfVetoAll=1, selfVetoFirst=2
     };
}

float GetDeltaR(float,float,float,float);

double ptRatioElec(const pat::Electron& elec,const pat::Jet& jet);
float ptRelElec(const pat::Electron& elec,const pat::Jet& jet);
float conePtElec(const pat::Electron& elec,const pat::Jet& jet);

double ptRatioMuon(const pat::Muon& muon,const pat::Jet& jet);
float ptRelMuon(const pat::Muon& muon,const pat::Jet& jet);
float conePtMuon(const pat::Muon& muon,const pat::Jet& jet);

float ElecPfIsoCharged(const pat::Electron& elec,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR);
float ElecPfIsoNeutral(const pat::Electron& elec,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR);

float MuonPfIsoCharged(const pat::Muon& muon,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR);
float MuonPfIsoNeutral(const pat::Muon& muon,edm::Handle<pat::PackedCandidateCollection> pfcands,float miniIsoR);

float isoSumRaw(const std::vector<const pat::PackedCandidate *> & cands, const reco::Candidate &cand, float dR, float innerR, float threshold, SelfVetoPolicy::SelfVetoPolicy selfVeto, int pdgId=-1);

double getEA(double eta, const std::vector<double>& EA_ETA, const std::vector<double>& EA_VALUE);

double getPFIsolation(edm::Handle<pat::PackedCandidateCollection>,
                      const reco::Candidate* ptcl,
                      double rho, double EA,
                      double r_iso_min, double r_iso_max, double kt_scale,
                      bool use_pfweight, bool charged_only);

bool qualityTrk(const reco::Track trk, const reco::Vertex &vtx);

int jetNDauChargedMVASel(const pat::Jet& jet, const reco::Vertex& vtx);

#endif
