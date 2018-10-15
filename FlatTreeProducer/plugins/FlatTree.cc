#include "VUBFlatTree/FlatTreeProducer/interface/FlatTree.hh"

FlatTree::FlatTree(TTree* _tree)
{
   tree = _tree;
   n_presel_jets_min = 0;
   n_presel_electrons_min = 0;
   n_presel_muons_min = 0;
   n_presel_leptons_min = 0;
   presel_MET_min = 0;
   apply_presel = false;
}

void FlatTree::Init()
{
   n_presel_jets = 0;
   n_presel_btag = 0;
   n_presel_electrons = 0;
   n_presel_muons = 0;
   n_presel_tau = 0;
   
   ev_run = DEFVAL;
   ev_id = DEFVAL;
   ev_lumi = DEFVAL;
   ev_rho = DEFVAL;

   met_px = DEFVAL;
   met_py = DEFVAL;
   met_pt = DEFVAL;
   met_phi = DEFVAL;
   met_sumet = DEFVAL;
   met_sig = DEFVAL;
   
   met_cov00 = DEFVAL;
   met_cov10 = DEFVAL;
   met_cov01 = DEFVAL;
   met_cov11 = DEFVAL;

   metGen_px = DEFVAL;
   metGen_py = DEFVAL;
   metGen_pt = DEFVAL;
   metGen_phi = DEFVAL;
   metGen_sumet = DEFVAL;
   
   metGen_NeutralEMEt = DEFVAL;
   metGen_ChargedEMEt = DEFVAL;
   metGen_NeutralHadEt = DEFVAL;
   metGen_ChargedHadEt = DEFVAL;
   metGen_MuonEt = DEFVAL;
   metGen_InvisibleEt = DEFVAL;
   
   met_uncorrectedPt = DEFVAL;
   met_uncorrectedPhi = DEFVAL;
   met_uncorrectedSumEt = DEFVAL;
   
   met_caloMETPt = DEFVAL;
   met_caloMETPhi = DEFVAL;
   met_caloMETSumEt = DEFVAL;
   
   met_shiftedPx_JetEnUp = DEFVAL;
   met_shiftedPx_JetEnDown = DEFVAL;
   met_shiftedPx_JetResUp = DEFVAL;
   met_shiftedPx_JetResDown = DEFVAL;
   met_shiftedPx_MuonEnUp = DEFVAL;
   met_shiftedPx_MuonEnDown = DEFVAL;
   met_shiftedPx_ElectronEnUp = DEFVAL;
   met_shiftedPx_ElectronEnDown = DEFVAL;
   met_shiftedPx_TauEnUp = DEFVAL;
   met_shiftedPx_TauEnDown = DEFVAL;
   met_shiftedPx_UnclusteredEnUp = DEFVAL;
   met_shiftedPx_UnclusteredEnDown = DEFVAL;
   met_shiftedPx_NoShift = DEFVAL;
   met_shiftedPx_PhotonEnUp = DEFVAL;
   met_shiftedPx_PhotonEnDown = DEFVAL;

   met_shiftedPy_JetEnUp = DEFVAL;
   met_shiftedPy_JetEnDown = DEFVAL;
   met_shiftedPy_JetResUp = DEFVAL;
   met_shiftedPy_JetResDown = DEFVAL;
   met_shiftedPy_MuonEnUp = DEFVAL;
   met_shiftedPy_MuonEnDown = DEFVAL;
   met_shiftedPy_ElectronEnUp = DEFVAL;
   met_shiftedPy_ElectronEnDown = DEFVAL;
   met_shiftedPy_TauEnUp = DEFVAL;
   met_shiftedPy_TauEnDown = DEFVAL;
   met_shiftedPy_UnclusteredEnUp = DEFVAL;
   met_shiftedPy_UnclusteredEnDown = DEFVAL;
   met_shiftedPy_NoShift = DEFVAL;
   met_shiftedPy_PhotonEnUp = DEFVAL;
   met_shiftedPy_PhotonEnDown = DEFVAL;

   met_shiftedPhi_JetEnUp = DEFVAL;
   met_shiftedPhi_JetEnDown = DEFVAL;
   met_shiftedPhi_JetResUp = DEFVAL;
   met_shiftedPhi_JetResDown = DEFVAL;
   met_shiftedPhi_MuonEnUp = DEFVAL;
   met_shiftedPhi_MuonEnDown = DEFVAL;
   met_shiftedPhi_ElectronEnUp = DEFVAL;
   met_shiftedPhi_ElectronEnDown = DEFVAL;
   met_shiftedPhi_TauEnUp = DEFVAL;
   met_shiftedPhi_TauEnDown = DEFVAL;
   met_shiftedPhi_UnclusteredEnUp = DEFVAL;
   met_shiftedPhi_UnclusteredEnDown = DEFVAL;
   met_shiftedPhi_NoShift = DEFVAL;
   met_shiftedPhi_PhotonEnUp = DEFVAL;
   met_shiftedPhi_PhotonEnDown = DEFVAL;
   
   met_shiftedSumEt_JetEnUp = DEFVAL;
   met_shiftedSumEt_JetEnDown = DEFVAL;
   met_shiftedSumEt_JetResUp = DEFVAL;
   met_shiftedSumEt_JetResDown = DEFVAL;
   met_shiftedSumEt_MuonEnUp = DEFVAL;
   met_shiftedSumEt_MuonEnDown = DEFVAL;
   met_shiftedSumEt_ElectronEnUp = DEFVAL;
   met_shiftedSumEt_ElectronEnDown = DEFVAL;
   met_shiftedSumEt_TauEnUp = DEFVAL;
   met_shiftedSumEt_TauEnDown = DEFVAL;
   met_shiftedSumEt_UnclusteredEnUp = DEFVAL;
   met_shiftedSumEt_UnclusteredEnDown = DEFVAL;
   met_shiftedSumEt_NoShift = DEFVAL;
   met_shiftedSumEt_PhotonEnUp = DEFVAL;
   met_shiftedSumEt_PhotonEnDown = DEFVAL;
  
   metNoHF_pt     = DEFVAL;
   metNoHF_phi    = DEFVAL;
   metNoHF_sumet  = DEFVAL;

   metPuppi_pt    = DEFVAL;
   metPuppi_phi   = DEFVAL;
   metPuppi_sumet = DEFVAL;
  
   nvertex = DEFVAL;

   pv_n = DEFVAL;
   pv_x = DEFVAL;
   pv_y = DEFVAL;
   pv_z = DEFVAL;

   pv_xError = DEFVAL;
   pv_yError = DEFVAL;
   pv_zError = DEFVAL;
   
   pv_chi2 = DEFVAL;
   pv_ndof = DEFVAL;
   pv_rho = DEFVAL;
   pv_isFake = DEFVAL;
   
   mc_weight = DEFVAL;
   mc_weight_originalValue = DEFVAL;
   mc_id = DEFVAL;
   mc_f1 = DEFVAL;
   mc_f2 = DEFVAL;
   mc_x1 = DEFVAL;
   mc_x2 = DEFVAL;
   mc_scale = DEFVAL;
   mc_ptHat = DEFVAL;
   
   weight_originalXWGTUP = DEFVAL;
   weight_scale_muF0p5 = DEFVAL;
   weight_scale_muF2   = DEFVAL;
   weight_scale_muR0p5 = DEFVAL;
   weight_scale_muR2   = DEFVAL;
   mc_pdfweights.clear();
   mc_pdfweightIds.clear();

   mc_pu_intime_NumInt = DEFVAL;
   mc_pu_trueNumInt = DEFVAL;
   mc_pu_before_npu = DEFVAL;
   mc_pu_after_npu = DEFVAL;

   mc_pu_Npvi = DEFVAL;
   mc_pu_Nzpositions.clear();
   mc_pu_BunchCrossing.clear();
   for(unsigned int i=0;i<mc_pu_zpositions.size();i++) mc_pu_zpositions[i].clear();
   mc_pu_zpositions.clear();
   for(unsigned int i=0;i<mc_pu_sumpT_lowpT.size();i++) mc_pu_sumpT_lowpT[i].clear();
   mc_pu_sumpT_lowpT.clear();
   for(unsigned int i=0;i<mc_pu_sumpT_highpT.size();i++) mc_pu_sumpT_highpT[i].clear();
   mc_pu_sumpT_highpT.clear();
   for(unsigned int i=0;i<mc_pu_ntrks_lowpT.size();i++) mc_pu_ntrks_lowpT[i].clear();
   mc_pu_ntrks_lowpT.clear();
   for(unsigned int i=0;i<mc_pu_ntrks_highpT.size();i++) mc_pu_ntrks_highpT[i].clear();
   mc_pu_ntrks_highpT.clear();

   trigger_n = 0;
   trigger.clear();
   trigger_name.clear();
   trigger_pass.clear();
   trigger_prescale.clear();
   trigger_HLTprescale.clear();
   trigger_L1prescale.clear();

   triggerobject_n = 0;
   triggerobject_pt.clear();
   triggerobject_eta.clear();
   triggerobject_phi.clear();

   triggerobject_collection.clear();

   triggerobject_filterIds_n.clear();
   triggerobject_filterIds.clear();

   triggerobject_isTriggerL1Mu.clear();
   triggerobject_isTriggerL1NoIsoEG.clear();
   triggerobject_isTriggerL1IsoEG.clear();
   triggerobject_isTriggerL1CenJet.clear();
   triggerobject_isTriggerL1ForJet.clear();
   triggerobject_isTriggerL1TauJet.clear();
   triggerobject_isTriggerL1ETM.clear();
   triggerobject_isTriggerL1ETT.clear();
   triggerobject_isTriggerL1HTT.clear();
   triggerobject_isTriggerL1HTM.clear();
   triggerobject_isTriggerL1JetCounts.clear();
   triggerobject_isTriggerL1HfBitCounts.clear();
   triggerobject_isTriggerL1HfRingEtSums.clear();
   triggerobject_isTriggerL1TechTrig.clear();
   triggerobject_isTriggerL1Castor.clear();
   triggerobject_isTriggerL1BPTX.clear();
   triggerobject_isTriggerL1GtExternal.clear();

   triggerobject_isHLT_TriggerPhoton.clear();
   triggerobject_isHLT_TriggerElectron.clear();
   triggerobject_isHLT_TriggerMuon.clear();
   triggerobject_isHLT_TriggerTau.clear();
   triggerobject_isHLT_TriggerJet.clear();
   triggerobject_isHLT_TriggerBJet.clear();
   triggerobject_isHLT_TriggerMET.clear();
   triggerobject_isHLT_TriggerTET.clear();
   triggerobject_isHLT_TriggerTHT.clear();
   triggerobject_isHLT_TriggerMHT.clear();
   triggerobject_isHLT_TriggerTrack.clear();
   triggerobject_isHLT_TriggerCluster.clear();
   triggerobject_isHLT_TriggerMETSig.clear();
   triggerobject_isHLT_TriggerELongit.clear();
   triggerobject_isHLT_TriggerMHTSig.clear();
   triggerobject_isHLT_TriggerHLongit.clear();

   triggerobject_filterLabels_n.clear();
   triggerobject_filterLabels.clear();

   triggerobject_pathNamesAll_n.clear();
   triggerobject_pathNamesAll.clear();
   triggerobject_pathNamesAll_isBoth.clear();
   triggerobject_pathNamesAll_isL3.clear();
   triggerobject_pathNamesAll_isLF.clear();
   triggerobject_pathNamesAll_isNone.clear();

   el_n = 0;
   el_pt.clear();
   el_eta.clear();
   el_phi.clear();
   el_m.clear();
   el_E.clear();
   el_id.clear();
   el_charge.clear();

   el_isGsfCtfScPixChargeConsistent.clear();
   el_isGsfScPixChargeConsistent.clear();
   el_passConversionVeto.clear();
   
   el_ip3d.clear();
   el_ip3dErr.clear();
   el_ip2d.clear();
   el_ip2dErr.clear();
   el_ip3dBS.clear();
   el_ip3dBSErr.clear();
   el_ip2dBS.clear();
   el_ip2dBSErr.clear();

   el_ecalEnergy.clear();
   el_correctedEcalEnergy.clear();
   el_correctedEcalEnergyError.clear();
   el_trackMomentumError.clear();
   
   el_neutralHadronIso.clear();
   el_chargedHadronIso.clear();
   el_puChargedHadronIso.clear();
   el_ecalIso.clear();
   el_hcalIso.clear();
   el_particleIso.clear();
   el_photonIso.clear();
   el_trackIso.clear();
   
   el_ecalPFClusterIso.clear();
   el_hcalPFClusterIso.clear();

   el_dr03EcalRecHitSumEt.clear();
   el_dr03HcalTowerSumEt.clear();
   el_dr03HcalDepth1TowerSumEt.clear();
   el_dr03HcalDepth2TowerSumEt.clear();
   el_dr03TkSumPt.clear();

   el_dr04EcalRecHitSumEt.clear();
   el_dr04HcalTowerSumEt.clear();
   el_dr04HcalDepth1TowerSumEt.clear();
   el_dr04HcalDepth2TowerSumEt.clear();
   el_dr04TkSumPt.clear();

   el_hcalOverEcal.clear();
   el_hcalOverEcalBc.clear();
   el_hcalDepth1OverEcal.clear();
   el_hcalDepth2OverEcal.clear();
   el_eSeedClusterOverPout.clear();
   el_eSeedClusterOverP.clear();
   el_eEleClusterOverPout.clear();
   el_deltaEtaEleClusterTrackAtCalo.clear();
   el_deltaPhiEleClusterTrackAtCalo.clear();
   
   el_pfIso_sumChargedHadronPt.clear();
   el_pfIso_sumNeutralHadronEt.clear();
   el_pfIso_sumPhotonEt.clear();
   el_pfIso_sumPUPt.clear();
   
   el_miniIso.clear();
   el_miniIsoTTH.clear();

   el_vx.clear();
   el_vy.clear();
   el_vz.clear();

   el_hasGsfTrack.clear();
   el_gsfTrack_d0.clear();
   el_gsfTrack_z0.clear();
   el_gsfTrack_d0Error.clear();
   el_gsfTrack_z0Error.clear();
   el_gsfTrack_PV_dxy.clear();
   el_gsfTrack_PV_dz.clear();
   el_gsfTrack_RP_dxy.clear();
   el_gsfTrack_RP_dz.clear();
   el_gsfTrack_BS_dxy.clear();
   el_gsfTrack_BS_dz.clear();
   el_gsfTrack_dxyError.clear();
   el_gsfTrack_dzError.clear();
   el_gsfTrack_normalizedChi2.clear();
   
   el_superCluster_eta.clear();
   el_superCluster_phi.clear();
   el_superCluster_energy.clear();
   el_superCluster_rawEnergy.clear();
   el_superCluster_preshowerEnergy.clear();
   el_superCluster_etaWidth.clear();
   el_superCluster_phiWidth.clear();
   el_superCluster_preshowerEnergyPlane1.clear();
   el_superCluster_preshowerEnergyPlane2.clear();
   el_superCluster_positionR.clear();
   
   el_basicClustersSize.clear();
   el_e1x5.clear();
   el_e5x5.clear();
   el_e2x5Max.clear();
   el_sigmaEtaEta.clear();
   el_sigmaIetaIeta.clear();
   el_sigmaIphiIphi.clear();
   el_sigmaIetaIphi.clear();
   el_full5x5_sigmaIphiIphi.clear();
   el_full5x5_sigmaEtaEta.clear();
   el_full5x5_sigmaIetaIeta.clear();
   el_full5x5_sigmaIetaIphi.clear();
   el_full5x5_r9.clear();
   el_full5x5_e1x5.clear();
   el_full5x5_e5x5.clear();
   el_full5x5_e2x5Max.clear();
   
   el_numberOfHits.clear();
   el_numberOfValidHits.clear();

   el_expectedMissingOuterHits.clear();
   el_numberOfValidPixelHits.clear();
   el_numberOfLostPixelHits.clear();
   el_trackerLayersWithMeasurement.clear();
   el_pixelLayersWithMeasurement.clear();
   el_numberOfValidStripLayersWithMonoAndStereo.clear();
   el_trackerLayersWithoutMeasurement.clear();
   
   el_hadronicOverEm.clear();
   el_numberOfLostHits.clear();
   el_numberOfLostHitsDefault.clear();

   el_fbrem.clear();
   el_kf_normalizedChi2.clear();
   el_gsf_normalizedChi2.clear();
   el_deltaEtaSuperClusterTrackAtVtx.clear();
   el_deltaPhiSuperClusterTrackAtVtx.clear();
   el_deltaEtaSeedClusterTrackAtCalo.clear();
   el_deltaPhiSeedClusterTrackAtCalo.clear();
   el_full5x5_OneMinusE1x5E5x5.clear();
   el_OneMinusE1x5E5x5.clear();
   el_eSuperClusterOverP.clear();
   el_IoEmIoP.clear();
   el_ooEmooP.clear();
   el_eleEoPout.clear();
   el_PreShowerOverRaw.clear();

   el_mvaIso.clear();
   el_mvaNoIso.clear();
   
   el_vetoCBId.clear();
   el_looseCBId.clear();
   el_mediumCBId.clear();
   el_tightCBId.clear();
  
//   el_vetoStopID.clear();
//   el_mediumStopID.clear();
   
   el_NoIso90MVAId.clear();
   el_NoIso80MVAId.clear();
   el_NoIsoLooseMVAId.clear();

   el_Iso90MVAId.clear();
   el_Iso80MVAId.clear();
   el_IsoLooseMVAId.clear();
   
   el_lepMVA.clear();
   
   el_lepMVA_pt.clear(); 
   el_lepMVA_eta.clear();
   el_lepMVA_miniRelIsoCharged.clear();
   el_lepMVA_miniRelIsoNeutral.clear();
   el_lepMVA_jetPtRatio.clear();
   el_lepMVA_jetPtRelv2.clear();
   el_lepMVA_jetBTagCSV.clear();
   el_lepMVA_sip3d.clear();
   el_lepMVA_dxy.clear();
   el_lepMVA_dz.clear();
   el_lepMVA_mvaId.clear();
   el_lepMVA_jetNDauChargedMVASel.clear();

   el_conept.clear();

   el_hasMCMatch.clear();
   el_gen_pt.clear();
   el_gen_eta.clear();
   el_gen_phi.clear();
   el_gen_m.clear();
   el_gen_status.clear();
   el_gen_id.clear();
   el_gen_charge.clear();
   el_gen_dr.clear();

   el_hasMCMatchPAT.clear();
   el_genPAT_pt.clear();
   el_genPAT_eta.clear();
   el_genPAT_phi.clear();
   el_genPAT_m.clear();
   el_genPAT_status.clear();
   el_genPAT_id.clear();
   el_genPAT_charge.clear();
   
   el_hasMatchedConversion.clear();
   el_expectedMissingInnerHits.clear();

   mu_n = 0;
   mu_pt.clear();
   mu_eta.clear();
   mu_phi.clear();
   mu_m.clear();
   mu_E.clear();
   mu_id.clear();
   mu_charge.clear();

   mu_ip3d.clear();
   mu_ip3dErr.clear();
   mu_ip2d.clear();
   mu_ip2dErr.clear();
   mu_ip3dBS.clear();
   mu_ip3dBSErr.clear();
   mu_ip2dBS.clear();
   mu_ip2dBSErr.clear();
   
   mu_neutralHadronIso.clear();
   mu_chargedHadronIso.clear();
   mu_puChargedHadronIso.clear();
   mu_ecalIso.clear();
   mu_hcalIso.clear();
   mu_photonIso.clear();
   mu_trackIso.clear();

   mu_pfIso03_sumChargedHadronPt.clear();
   mu_pfIso03_sumChargedParticlePt.clear();
   mu_pfIso03_sumNeutralHadronEt.clear();
   mu_pfIso03_sumNeutralHadronEtHighThreshold.clear();
   mu_pfIso03_sumPhotonEt.clear();
   mu_pfIso03_sumPhotonEtHighThreshold.clear();
   mu_pfIso03_sumPUPt.clear();

   mu_pfIso04_sumChargedHadronPt.clear();
   mu_pfIso04_sumChargedParticlePt.clear();
   mu_pfIso04_sumNeutralHadronEt.clear();
   mu_pfIso04_sumNeutralHadronEtHighThreshold.clear();
   mu_pfIso04_sumPhotonEt.clear();
   mu_pfIso04_sumPhotonEtHighThreshold.clear();
   mu_pfIso04_sumPUPt.clear();

   mu_pfMeanIso03_sumChargedHadronPt.clear();
   mu_pfMeanIso03_sumChargedParticlePt.clear();
   mu_pfMeanIso03_sumNeutralHadronEt.clear();
   mu_pfMeanIso03_sumNeutralHadronEtHighThreshold.clear();
   mu_pfMeanIso03_sumPhotonEt.clear();
   mu_pfMeanIso03_sumPhotonEtHighThreshold.clear();
   mu_pfMeanIso03_sumPUPt.clear();

   mu_pfSumIso03_sumChargedHadronPt.clear();
   mu_pfSumIso03_sumChargedParticlePt.clear();
   mu_pfSumIso03_sumNeutralHadronEt.clear();
   mu_pfSumIso03_sumNeutralHadronEtHighThreshold.clear();
   mu_pfSumIso03_sumPhotonEt.clear();
   mu_pfSumIso03_sumPhotonEtHighThreshold.clear();
   mu_pfSumIso03_sumPUPt.clear();

   mu_pfMeanIso04_sumChargedHadronPt.clear();
   mu_pfMeanIso04_sumChargedParticlePt.clear();
   mu_pfMeanIso04_sumNeutralHadronEt.clear();
   mu_pfMeanIso04_sumNeutralHadronEtHighThreshold.clear();
   mu_pfMeanIso04_sumPhotonEt.clear();
   mu_pfMeanIso04_sumPhotonEtHighThreshold.clear();
   mu_pfMeanIso04_sumPUPt.clear();

   mu_pfSumIso04_sumChargedHadronPt.clear();
   mu_pfSumIso04_sumChargedParticlePt.clear();
   mu_pfSumIso04_sumNeutralHadronEt.clear();
   mu_pfSumIso04_sumNeutralHadronEtHighThreshold.clear();
   mu_pfSumIso04_sumPhotonEt.clear();
   mu_pfSumIso04_sumPhotonEtHighThreshold.clear();
   mu_pfSumIso04_sumPUPt.clear();
   
   mu_miniIso.clear();
   mu_miniIsoTTH.clear();

   mu_isGlobalMuon.clear();
   mu_isTrackerMuon.clear();
   mu_isStandAloneMuon.clear();
   mu_isCaloMuon.clear();
   mu_isPFMuon.clear();
   mu_isRPCMuon.clear();

   mu_vx.clear();
   mu_vy.clear();
   mu_vz.clear();
   
   mu_numberOfMatches.clear();
   mu_numberOfMatchedStations.clear();
   
   mu_segmentCompatibility.clear();
   mu_caloCompatibility.clear();

   mu_combinedQuality_chi2LocalPosition.clear();
   mu_combinedQuality_trkKink.clear();
   
   mu_isLooseMuon.clear();
   mu_isMediumMuon.clear();
   mu_isTightMuon.clear();
   mu_isSoftMuon.clear();
   mu_isHighPtMuon.clear();

   mu_isGoodMuon_AllGlobalMuons.clear();
   mu_isGoodMuon_AllStandAloneMuons.clear();
   mu_isGoodMuon_AllTrackerMuons.clear();
   mu_isGoodMuon_TrackerMuonArbitrated.clear();
   mu_isGoodMuon_AllArbitrated.clear();
   mu_isGoodMuon_GlobalMuonPromptTight.clear();
   mu_isGoodMuon_TMLastStationLoose.clear();
   mu_isGoodMuon_TMLastStationTight.clear();
   mu_isGoodMuon_TM2DCompatibilityLoose.clear();
   mu_isGoodMuon_TM2DCompatibilityTight.clear();
   mu_isGoodMuon_TMOneStationLoose.clear();
   mu_isGoodMuon_TMOneStationTight.clear();
   mu_isGoodMuon_TMLastStationOptimizedLowPtLoose.clear();
   mu_isGoodMuon_TMLastStationOptimizedLowPtTight.clear();
   mu_isGoodMuon_GMTkChiCompatibility.clear();
   mu_isGoodMuon_GMStaChiCompatibility.clear();
   mu_isGoodMuon_GMTkKinkTight.clear();
   mu_isGoodMuon_TMLastStationAngLoose.clear();
   mu_isGoodMuon_TMLastStationAngTight.clear();
   mu_isGoodMuon_TMOneStationAngLoose.clear();
   mu_isGoodMuon_TMOneStationAngTight.clear();
   mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtLoose.clear();
   mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight.clear();

   mu_calEnergy_em.clear();
   mu_calEnergy_had.clear();
   mu_calEnergy_ho.clear();
   mu_calEnergy_emS9.clear();
   mu_calEnergy_hadS9.clear();
   mu_calEnergy_hoS9.clear();
   mu_calEnergy_emS25.clear();
   mu_calEnergy_emMax.clear();
   mu_calEnergy_hadMax.clear();
   mu_calEnergy_ecal_time.clear();
   mu_calEnergy_hcal_time.clear();
   mu_calEnergy_ecal_rawId.clear();
   mu_calEnergy_hcal_rawId.clear();

   mu_isolationR03_trackerVetoPt.clear();
   mu_isolationR03_emVetoEt.clear();
   mu_isolationR03_hadVetoEt.clear();
   mu_isolationR03_hoVetoEt.clear();
   mu_isolationR03_sumPt.clear();
   mu_isolationR03_emEt.clear();
   mu_isolationR03_hadEt.clear();
   mu_isolationR03_hoEt.clear();
   mu_isolationR03_nTracks.clear();
   mu_isolationR03_nJets.clear();
   
   mu_isolationR05_trackerVetoPt.clear();
   mu_isolationR05_emVetoEt.clear();
   mu_isolationR05_hadVetoEt.clear();
   mu_isolationR05_hoVetoEt.clear();
   mu_isolationR05_sumPt.clear();
   mu_isolationR05_emEt.clear();
   mu_isolationR05_hadEt.clear();
   mu_isolationR05_hoEt.clear();
   mu_isolationR05_nTracks.clear();
   mu_isolationR05_nJets.clear();

   mu_hasGlobalTrack.clear();
   mu_globalTrack_d0.clear();
   mu_globalTrack_z0.clear();
   mu_globalTrack_d0Error.clear();
   mu_globalTrack_z0Error.clear();
   mu_globalTrack_RP_dxy.clear();
   mu_globalTrack_RP_dz.clear();
   mu_globalTrack_PV_dxy.clear();
   mu_globalTrack_PV_dz.clear();
   mu_globalTrack_BS_dxy.clear();
   mu_globalTrack_BS_dz.clear();
   mu_globalTrack_dxyError.clear();
   mu_globalTrack_dzError.clear();
   mu_globalTrack_normalizedChi2.clear();
   mu_globalTrack_numberOfValidHits.clear();
   mu_globalTrack_numberOfValidMuonHits.clear();
   mu_globalTrack_numberOfLostHits.clear();
   mu_globalTrack_pt.clear();
   mu_globalTrack_eta.clear();
   mu_globalTrack_phi.clear();
   mu_globalTrack_ptError.clear();
   mu_globalTrack_etaError.clear();
   mu_globalTrack_phiError.clear();
   mu_globalTrack_vx.clear();
   mu_globalTrack_vy.clear();
   mu_globalTrack_vz.clear();
   mu_globalTrack_qoverp.clear();
   mu_globalTrack_qoverpError.clear();
   mu_globalTrack_charge.clear();
   mu_globalTrack_trackerLayersWithMeasurement.clear();
   mu_globalTrack_pixelLayersWithMeasurement.clear();
   mu_globalTrack_numberOfValidStripLayersWithMonoAndStereo.clear();
   mu_globalTrack_trackerLayersWithoutMeasurement.clear();
   mu_globalTrack_numberOfValidPixelHits.clear();
   mu_globalTrack_numberOfLostPixelHits.clear();
   mu_globalTrack_numberOfInnerHits.clear();
   mu_globalTrack_numberOfOuterHits.clear();
   mu_globalTrack_validFraction.clear();
   
   mu_bestTrackType.clear();
   mu_hasBestTrack.clear();
   mu_bestTrack_d0.clear();
   mu_bestTrack_z0.clear();
   mu_bestTrack_d0Error.clear();
   mu_bestTrack_z0Error.clear();
   mu_bestTrack_RP_dxy.clear();
   mu_bestTrack_RP_dz.clear();
   mu_bestTrack_PV_dxy.clear();
   mu_bestTrack_PV_dz.clear();
   mu_bestTrack_BS_dxy.clear();
   mu_bestTrack_BS_dz.clear();
   mu_bestTrack_dxyError.clear();
   mu_bestTrack_dzError.clear();
   mu_bestTrack_normalizedChi2.clear();
   mu_bestTrack_numberOfValidHits.clear();
   mu_bestTrack_numberOfLostHits.clear();
   mu_bestTrack_pt.clear();
   mu_bestTrack_eta.clear();
   mu_bestTrack_phi.clear();
   mu_bestTrack_ptError.clear();
   mu_bestTrack_etaError.clear();
   mu_bestTrack_phiError.clear();
   mu_bestTrack_vx.clear();
   mu_bestTrack_vy.clear();
   mu_bestTrack_vz.clear();
   mu_bestTrack_qoverp.clear();
   mu_bestTrack_qoverpError.clear();
   mu_bestTrack_charge.clear();
   mu_bestTrack_trackerLayersWithMeasurement.clear();
   mu_bestTrack_pixelLayersWithMeasurement.clear();
   mu_bestTrack_numberOfValidStripLayersWithMonoAndStereo.clear();
   mu_bestTrack_trackerLayersWithoutMeasurement.clear();
   mu_bestTrack_numberOfValidPixelHits.clear();
   mu_bestTrack_numberOfLostPixelHits.clear();
   mu_bestTrack_numberOfInnerHits.clear();
   mu_bestTrack_numberOfOuterHits.clear();
   mu_bestTrack_validFraction.clear();

   mu_hasInnerTrack.clear();
   mu_innerTrack_d0.clear();
   mu_innerTrack_z0.clear();
   mu_innerTrack_d0Error.clear();
   mu_innerTrack_z0Error.clear();
   mu_innerTrack_RP_dxy.clear();
   mu_innerTrack_RP_dz.clear();
   mu_innerTrack_PV_dxy.clear();
   mu_innerTrack_PV_dz.clear();
   mu_innerTrack_BS_dxy.clear();
   mu_innerTrack_BS_dz.clear();
   mu_innerTrack_dxyError.clear();
   mu_innerTrack_dzError.clear();
   mu_innerTrack_normalizedChi2.clear();
   mu_innerTrack_numberOfValidHits.clear();
   mu_innerTrack_numberOfLostHits.clear();
   mu_innerTrack_pt.clear();
   mu_innerTrack_eta.clear();
   mu_innerTrack_phi.clear();
   mu_innerTrack_ptError.clear();
   mu_innerTrack_etaError.clear();
   mu_innerTrack_phiError.clear();
   mu_innerTrack_vx.clear();
   mu_innerTrack_vy.clear();
   mu_innerTrack_vz.clear();
   mu_innerTrack_qoverp.clear();
   mu_innerTrack_qoverpError.clear();
   mu_innerTrack_charge.clear();
   mu_innerTrack_trackerLayersWithMeasurement.clear();
   mu_innerTrack_pixelLayersWithMeasurement.clear();
   mu_innerTrack_numberOfValidStripLayersWithMonoAndStereo.clear();
   mu_innerTrack_trackerLayersWithoutMeasurement.clear();
   mu_innerTrack_numberOfValidPixelHits.clear();
   mu_innerTrack_numberOfLostPixelHits.clear();
   mu_innerTrack_numberOfInnerHits.clear();
   mu_innerTrack_numberOfOuterHits.clear();
   mu_innerTrack_validFraction.clear();
   
   mu_type.clear();

   mu_lepMVA.clear();

   mu_lepMVA_pt.clear(); 
   mu_lepMVA_eta.clear();
   mu_lepMVA_miniRelIsoCharged.clear();
   mu_lepMVA_miniRelIsoNeutral.clear();
   mu_lepMVA_jetPtRatio.clear();
   mu_lepMVA_jetPtRelv2.clear();
   mu_lepMVA_jetBTagCSV.clear();
   mu_lepMVA_sip3d.clear();
   mu_lepMVA_dxy.clear();
   mu_lepMVA_dz.clear();
   mu_lepMVA_mvaId.clear();
   mu_lepMVA_jetNDauChargedMVASel.clear();

   mu_conept.clear();

   mu_hasMCMatch.clear();
   mu_gen_pt.clear();
   mu_gen_eta.clear();
   mu_gen_phi.clear();
   mu_gen_m.clear();
   mu_gen_status.clear();
   mu_gen_id.clear();
   mu_gen_charge.clear();
   mu_gen_dr.clear();

   mu_hasMCMatchPAT.clear();
   mu_genPAT_pt.clear();
   mu_genPAT_eta.clear();
   mu_genPAT_phi.clear();
   mu_genPAT_m.clear();
   mu_genPAT_status.clear();
   mu_genPAT_id.clear();
   mu_genPAT_charge.clear();
   
   tau_n = 0;
   tau_pt.clear();
   tau_eta.clear();
   tau_phi.clear();
   tau_m.clear();
   tau_E.clear();
   tau_id.clear();
   tau_charge.clear();
   
   tau_hasLeadChargedHadrCand.clear();
   tau_leadingTrackPt.clear();
   tau_leadingTrackCharge.clear();
   tau_leadingTrackDz.clear();
   tau_leadingTrackDxy.clear();
   
   tau_decayMode.clear();
   tau_decayModeFinding.clear();
//   tau_decayModeFindingOldDMs.clear();
   tau_decayModeFindingNewDMs.clear();
   
   tau_puCorrPtSum.clear();
   tau_neutralIsoPtSum.clear();
   tau_chargedIsoPtSum.clear();
   
   tau_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();
   
   tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
   tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
   tau_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();

   tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
   tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT.clear();
   tau_byTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
   tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
   tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
   
   tau_againstMuonLoose3.clear();
   tau_againstMuonTight3.clear();
   
   tau_againstElectronVLooseMVA6.clear();
   tau_againstElectronLooseMVA6.clear();
   tau_againstElectronMediumMVA6.clear();
   tau_againstElectronTightMVA6.clear();
   
   tau_pfEssential_jet_pt.clear();
   tau_pfEssential_jet_eta.clear();
   tau_pfEssential_jet_phi.clear();
   tau_pfEssential_jet_m.clear();

   tau_pfEssential_jetCorr_pt.clear();
   tau_pfEssential_jetCorr_eta.clear();
   tau_pfEssential_jetCorr_phi.clear();
   tau_pfEssential_jetCorr_m.clear();
   
   tau_pfEssential_hasSV.clear();
   tau_pfEssential_sv_x.clear();
   tau_pfEssential_sv_y.clear();
   tau_pfEssential_sv_z.clear();
   
   tau_pfEssential_flightLengthSig.clear();
   tau_pfEssential_dxy.clear();
   tau_pfEssential_dxy_error.clear();
   tau_pfEssential_dxy_Sig.clear();
   
   jet_n = 0;
   jet_pt.clear();
   jet_eta.clear();
   jet_phi.clear();
   jet_m.clear();
   jet_E.clear();

   jet_ntrk.clear();

   jet_JBP.clear();
   jet_JP.clear();
   jet_TCHP.clear();
   jet_TCHE.clear();
   jet_SSVHE.clear();
   jet_SSVHP.clear();
   jet_CMVA.clear();
   jet_CSVv2.clear();
   jet_DeepCSVProbudsg.clear();
   jet_DeepCSVProbb.clear();
   jet_DeepCSVProbc.clear();
   jet_DeepCSVProbbb.clear();
   jet_DeepCSVProbcc.clear();
   jet_cMVAv2.clear();
   jet_CharmCvsL.clear();
   jet_CharmCvsB.clear();
   jet_partonFlavour.clear();
   jet_hadronFlavour.clear();

   jet_neutralHadronEnergy.clear();
   jet_neutralEmEnergy.clear();
   jet_chargedHadronEnergy.clear();
   jet_chargedEmEnergy.clear();
   jet_electronEnergy.clear();
   jet_muonEnergy.clear();
   jet_photonEnergy.clear();

   jet_chargedMultiplicity.clear();
   jet_neutralMultiplicity.clear();
   jet_chargedHadronMultiplicity.clear();
   
   jet_jetArea.clear();
   
   jet_jecFactorUncorrected.clear();
   jet_jecFactorL1FastJet.clear();
   jet_jecFactorL2Relative.clear();
   jet_jecFactorL3Absolute.clear();

   jet_neutralHadronEnergyFraction.clear();
   jet_neutralEmEnergyFraction.clear();
   jet_chargedHadronEnergyFraction.clear();
   jet_muonEnergyFraction.clear();
   jet_chargedEmEnergyFraction.clear();
   
   jet_Unc.clear();
   
   jet_pileupJetId.clear();

//   jet_looseJetID.clear();
   jet_tightJetID.clear();
   jet_tightLepVetoJetID.clear();
   jet_qgtag.clear();   
   
   jet_hasGenJet.clear();   
   jet_genJet_pt.clear();
   jet_genJet_eta.clear();
   jet_genJet_phi.clear();
   jet_genJet_m.clear();
   jet_genJet_E.clear();
   jet_genJet_status.clear();
   jet_genJet_id.clear();

   jet_hasGenParton.clear();   
   jet_genParton_pt.clear();
   jet_genParton_eta.clear();
   jet_genParton_phi.clear();
   jet_genParton_m.clear();
   jet_genParton_E.clear();
   jet_genParton_status.clear();
   jet_genParton_id.clear();
   
   jetPuppi_n = 0;
   jetPuppi_pt.clear();
   jetPuppi_eta.clear();
   jetPuppi_phi.clear();
   jetPuppi_m.clear();
   jetPuppi_E.clear();

   jetPuppi_ntrk.clear();

   jetPuppi_JBP.clear();
   jetPuppi_JP.clear();
   jetPuppi_TCHP.clear();
   jetPuppi_TCHE.clear();
   jetPuppi_SSVHE.clear();
   jetPuppi_SSVHP.clear();
   jetPuppi_CMVA.clear();
   jetPuppi_CSVv2.clear();
   jetPuppi_partonFlavour.clear();
   jetPuppi_hadronFlavour.clear();

   jetPuppi_neutralHadronEnergy.clear();
   jetPuppi_neutralEmEnergy.clear();
   jetPuppi_chargedHadronEnergy.clear();
   jetPuppi_chargedEmEnergy.clear();
   jetPuppi_electronEnergy.clear();
   jetPuppi_muonEnergy.clear();
   jetPuppi_photonEnergy.clear();

   jetPuppi_chargedMultiplicity.clear();
   jetPuppi_neutralMultiplicity.clear();
   jetPuppi_chargedHadronMultiplicity.clear();
   
   jetPuppi_jecFactorUncorrected.clear();
   jetPuppi_jecFactorL1FastJet.clear();
   jetPuppi_jecFactorL2Relative.clear();
   jetPuppi_jecFactorL3Absolute.clear();
   
   jetPuppi_pileupJetId.clear();

   jetPuppi_hasGenJet.clear();   
   jetPuppi_genJet_pt.clear();
   jetPuppi_genJet_eta.clear();
   jetPuppi_genJet_phi.clear();
   jetPuppi_genJet_m.clear();
   jetPuppi_genJet_E.clear();
   jetPuppi_genJet_status.clear();
   jetPuppi_genJet_id.clear();

   jetPuppi_hasGenParton.clear();   
   jetPuppi_genParton_pt.clear();
   jetPuppi_genParton_eta.clear();
   jetPuppi_genParton_phi.clear();
   jetPuppi_genParton_m.clear();
   jetPuppi_genParton_E.clear();
   jetPuppi_genParton_status.clear();
   jetPuppi_genParton_id.clear();
   
   
   //------------------------
   //  ak8 collection
   //------------------------
   
   ak8jet_n = 0;
   ak8jet_pt.clear();
   ak8jet_eta.clear();
   ak8jet_phi.clear();
   ak8jet_m.clear();
   ak8jet_E.clear();
   
   ak8jet_ntrk.clear();
   
   ak8jet_JBP.clear();
   ak8jet_JP.clear();
   ak8jet_TCHP.clear();
   ak8jet_TCHE.clear();
   ak8jet_SSVHE.clear();
   ak8jet_SSVHP.clear();
   ak8jet_CMVA.clear();
   
   ak8jet_CSVv2.clear();
   ak8jet_partonFlavour.clear();
   ak8jet_hadronFlavour.clear();
   
   ak8jet_neutralHadronEnergy.clear();
   ak8jet_neutralEmEnergy.clear();
   ak8jet_chargedHadronEnergy.clear();
   ak8jet_chargedEmEnergy.clear();
   ak8jet_electronEnergy.clear();
   ak8jet_muonEnergy.clear();
   ak8jet_photonEnergy.clear();
   
   ak8jet_chargedMultiplicity.clear();
   ak8jet_neutralMultiplicity.clear();
   ak8jet_chargedHadronMultiplicity.clear();
   
   ak8jet_jetArea.clear();
   
   ak8jet_jecFactorUncorrected.clear();
   ak8jet_jecFactorL1FastJet.clear();
   ak8jet_jecFactorL2Relative.clear();
   ak8jet_jecFactorL3Absolute.clear();
   
   ak8jet_pileupJetId.clear();
   
//   ak8jet_looseJetID.clear();
   ak8jet_tightJetID.clear();
   
   ak8jet_hasGenJet.clear();
   ak8jet_genJet_pt.clear();
   ak8jet_genJet_eta.clear();
   ak8jet_genJet_phi.clear();
   ak8jet_genJet_m.clear();
   ak8jet_genJet_E.clear();
   ak8jet_genJet_status.clear();
   ak8jet_genJet_id.clear();
   
   ak8jet_hasGenParton.clear();
   ak8jet_genParton_pt.clear();
   ak8jet_genParton_eta.clear();
   ak8jet_genParton_phi.clear();
   ak8jet_genParton_m.clear();
   ak8jet_genParton_E.clear();
   ak8jet_genParton_status.clear();
   ak8jet_genParton_id.clear();
   
   ak8jet_tau1.clear();
   ak8jet_tau2.clear();
   ak8jet_tau3.clear();
   ak8jet_softdrop_mass.clear();
   ak8jet_trimmed_mass.clear();
   ak8jet_pruned_mass.clear();
   ak8jet_filtered_mass.clear();
   ak8jet_minMass.clear();
   ak8jet_topMass.clear();
   ak8jet_nSubJets.clear();
  
   //------------------------
   //  ak10 collection
   //------------------------
   
   ak10jet_n = 0;
   ak10jet_pt.clear();
   ak10jet_eta.clear();
   ak10jet_phi.clear();
   ak10jet_m.clear();
   ak10jet_E.clear();
   
   ak10jet_ntrk.clear();
   
   ak10jet_JBP.clear();
   ak10jet_JP.clear();
   ak10jet_TCHP.clear();
   ak10jet_TCHE.clear();
   ak10jet_SSVHE.clear();
   ak10jet_SSVHP.clear();
   ak10jet_CMVA.clear();
   
   ak10jet_CSVv2.clear();
   ak10jet_partonFlavour.clear();
   ak10jet_hadronFlavour.clear();
   
   ak10jet_neutralHadronEnergy.clear();
   ak10jet_neutralEmEnergy.clear();
   ak10jet_chargedHadronEnergy.clear();
   ak10jet_chargedEmEnergy.clear();
   ak10jet_electronEnergy.clear();
   ak10jet_muonEnergy.clear();
   ak10jet_photonEnergy.clear();
   
   ak10jet_chargedMultiplicity.clear();
   ak10jet_neutralMultiplicity.clear();
   ak10jet_chargedHadronMultiplicity.clear();
   
   ak10jet_jetArea.clear();
   
   ak10jet_jecFactorUncorrected.clear();
   ak10jet_jecFactorL1FastJet.clear();
   ak10jet_jecFactorL2Relative.clear();
   ak10jet_jecFactorL3Absolute.clear();
   
   ak10jet_pileupJetId.clear();
   
//   ak10jet_looseJetID.clear();
   ak10jet_tightJetID.clear();
   
   ak10jet_hasGenJet.clear();
   ak10jet_genJet_pt.clear();
   ak10jet_genJet_eta.clear();
   ak10jet_genJet_phi.clear();
   ak10jet_genJet_m.clear();
   ak10jet_genJet_E.clear();
   ak10jet_genJet_status.clear();
   ak10jet_genJet_id.clear();
   
   ak10jet_hasGenParton.clear();
   ak10jet_genParton_pt.clear();
   ak10jet_genParton_eta.clear();
   ak10jet_genParton_phi.clear();
   ak10jet_genParton_m.clear();
   ak10jet_genParton_E.clear();
   ak10jet_genParton_status.clear();
   ak10jet_genParton_id.clear();
   
   ak10jet_tau1.clear();
   ak10jet_tau2.clear();
   ak10jet_tau3.clear();
   ak10jet_softdrop_mass.clear();
   ak10jet_trimmed_mass.clear();
   ak10jet_pruned_mass.clear();
   ak10jet_filtered_mass.clear();
   ak10jet_minMass.clear();
   ak10jet_topMass.clear();
   ak10jet_nSubJets.clear();
   
   //------------------------
   //  GenJet collection
   //------------------------
   
   genJet_n = 0;
   genJet_pt.clear();
   genJet_eta.clear();
   genJet_phi.clear();
   genJet_m.clear();
   genJet_E.clear();
   genJet_emEnergy.clear();
   genJet_hadEnergy.clear();
   genJet_invisibleEnergy.clear();
   genJet_auxiliaryEnergy.clear();
   genJet_flavour.clear();
   
   //------------------------
   //  GenTTXJet collection
   //------------------------
   
   genTTXJet_n = 0;
   genTTXJet_pt.clear();
   genTTXJet_eta.clear();
   genTTXJet_phi.clear();
   genTTXJet_m.clear();
   genTTXJet_E.clear();
   genTTXJet_flavour.clear();
   
   pfcand_n = 0;
   pfcand_pt.clear();
   pfcand_eta.clear();
   pfcand_phi.clear();
   pfcand_E.clear();
   pfcand_charge.clear();
   pfcand_id.clear();
   pfcand_dz.clear();
   pfcand_trackIso.clear();

   pfch_loose_n = 0;
   pfch_loose_sumpt = 0;
   pfch_tight_n = 0;
   pfch_tight_sumpt = 0;

   
   //------------------------
   //  GenInfo stop
   //------------------------
   gen_stop_m.clear();
   gen_neutralino_m.clear();

   genTTX_id = DEFVAL;
}

void FlatTree::CreateBranches(int buffersize = 32000)
{
   if( doWrite("ev_run") ) tree->Branch("ev_run", &ev_run, "ev_run/I", buffersize);
   if( doWrite("ev_id") ) tree->Branch("ev_id", &ev_id, "ev_id/I", buffersize);
   if( doWrite("ev_lumi") ) tree->Branch("ev_lumi", &ev_lumi, "ev_lumi/I", buffersize);
   if( doWrite("ev_rho") ) tree->Branch("ev_rho", &ev_rho, "ev_rho/F", buffersize);
   
   if( doWrite("met_px") ) tree->Branch("met_px", &met_px, "met_px/F", buffersize);
   if( doWrite("met_py") ) tree->Branch("met_py", &met_py, "met_py/F", buffersize);
   if( doWrite("met_pt") ) tree->Branch("met_pt", &met_pt, "met_pt/F", buffersize);
   if( doWrite("met_phi") ) tree->Branch("met_phi", &met_phi, "met_phi/F", buffersize);
   if( doWrite("met_sumet") ) tree->Branch("met_sumet", &met_sumet, "met_sumet/F", buffersize);
   if( doWrite("met_sig") ) tree->Branch("met_sig", &met_sig, "met_sig/D", buffersize);
   
   if( doWrite("met_cov00") ) tree->Branch("met_cov00", &met_cov00, "met_cov00/D", buffersize);
   if( doWrite("met_cov10") ) tree->Branch("met_cov10", &met_cov10, "met_cov10/D", buffersize);
   if( doWrite("met_cov01") ) tree->Branch("met_cov01", &met_cov01, "met_cov01/D", buffersize);
   if( doWrite("met_cov11") ) tree->Branch("met_cov11", &met_cov11, "met_cov11/D", buffersize);

   if( doWrite("metGen_px") ) tree->Branch("metGen_px", &metGen_px, "metGen_px/F", buffersize);
   if( doWrite("metGen_py") ) tree->Branch("metGen_py", &metGen_py, "metGen_py/F", buffersize);
   if( doWrite("metGen_pt") ) tree->Branch("metGen_pt", &metGen_pt, "metGen_pt/F", buffersize);
   if( doWrite("metGen_phi") ) tree->Branch("metGen_phi", &metGen_phi, "metGen_phi/F", buffersize);
   if( doWrite("metGen_sumet") ) tree->Branch("metGen_sumet", &metGen_sumet, "metGen_sumet/F", buffersize);
   
   if( doWrite("metGen_NeutralEMEt") ) tree->Branch("metGen_NeutralEMEt", &metGen_NeutralEMEt, "metGen_NeutralEMEt/F", buffersize);
   if( doWrite("metGen_ChargedEMEt") ) tree->Branch("metGen_ChargedEMEt", &metGen_ChargedEMEt, "metGen_ChargedEMEt/F", buffersize);
   if( doWrite("metGen_NeutralHadEt") ) tree->Branch("metGen_NeutralHadEt", &metGen_NeutralHadEt, "metGen_NeutralHadEt/F", buffersize);
   if( doWrite("metGen_ChargedHadEt") ) tree->Branch("metGen_ChargedHadEt", &metGen_ChargedHadEt, "metGen_ChargedHadEt/F", buffersize);
   if( doWrite("metGen_MuonEt") ) tree->Branch("metGen_MuonEt", &metGen_MuonEt, "metGen_MuonEt/F", buffersize);
   if( doWrite("metGen_InvisibleEt") ) tree->Branch("metGen_InvisibleEt", &metGen_InvisibleEt, "metGen_InvisibleEt/F", buffersize);
   
   if( doWrite("met_uncorrectedPt") ) tree->Branch("met_uncorrectedPt", &met_uncorrectedPt, "met_uncorrectedPt/F", buffersize);
   if( doWrite("met_uncorrectedPhi") ) tree->Branch("met_uncorrectedPhi", &met_uncorrectedPhi, "met_uncorrectedPhi/F", buffersize);
   if( doWrite("met_uncorrectedSumEt") ) tree->Branch("met_uncorrectedSumEt", &met_uncorrectedSumEt, "met_uncorrectedSumEt/F", buffersize);

   if( doWrite("met_caloMETPt") ) tree->Branch("met_caloMETPt", &met_caloMETPt, "met_caloMETPt/F", buffersize);
   if( doWrite("met_caloMETPhi") ) tree->Branch("met_caloMETPhi", &met_caloMETPhi, "met_caloMETPhi/F", buffersize);
   if( doWrite("met_caloMETSumEt") ) tree->Branch("met_caloMETSumEt", &met_caloMETSumEt, "met_caloMETSumEt/F", buffersize);
   
   if( doWrite("met_shiftedPx_JetEnUp") ) tree->Branch("met_shiftedPx_JetEnUp", &met_shiftedPx_JetEnUp, "met_shiftedPx_JetEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_JetEnDown") ) tree->Branch("met_shiftedPx_JetEnDown", &met_shiftedPx_JetEnDown, "met_shiftedPx_JetEnDown/F", buffersize);
   if( doWrite("met_shiftedPx_JetResUp") ) tree->Branch("met_shiftedPx_JetResUp", &met_shiftedPx_JetResUp, "met_shiftedPx_JetResUp/F", buffersize);
   if( doWrite("met_shiftedPx_JetResDown") ) tree->Branch("met_shiftedPx_JetResDown", &met_shiftedPx_JetResDown, "met_shiftedPx_JetResDown/F", buffersize);
   if( doWrite("met_shiftedPx_MuonEnUp") ) tree->Branch("met_shiftedPx_MuonEnUp", &met_shiftedPx_MuonEnUp, "met_shiftedPx_MuonEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_MuonEnDown") ) tree->Branch("met_shiftedPx_MuonEnDown", &met_shiftedPx_MuonEnDown, "met_shiftedPx_MuonEnDown/F", buffersize);
   if( doWrite("met_shiftedPx_ElectronEnUp") ) tree->Branch("met_shiftedPx_ElectronEnUp", &met_shiftedPx_ElectronEnUp, "met_shiftedPx_ElectronEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_ElectronEnDown") ) tree->Branch("met_shiftedPx_ElectronEnDown", &met_shiftedPx_ElectronEnDown, "met_shiftedPx_ElectronEnDown/F", buffersize);
   if( doWrite("met_shiftedPx_TauEnUp") ) tree->Branch("met_shiftedPx_TauEnUp", &met_shiftedPx_TauEnUp, "met_shiftedPx_TauEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_TauEnDown") ) tree->Branch("met_shiftedPx_TauEnDown", &met_shiftedPx_TauEnDown, "met_shiftedPx_TauEnDown/F", buffersize);
   if( doWrite("met_shiftedPx_UnclusteredEnUp") ) tree->Branch("met_shiftedPx_UnclusteredEnUp", &met_shiftedPx_UnclusteredEnUp, "met_shiftedPx_UnclusteredEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_UnclusteredEnDown") ) tree->Branch("met_shiftedPx_UnclusteredEnDown", &met_shiftedPx_UnclusteredEnDown, "met_shiftedPx_UnclusteredEnDown/F", buffersize);
   if( doWrite("met_shiftedPx_NoShift") ) tree->Branch("met_shiftedPx_NoShift", &met_shiftedPx_NoShift, "met_shiftedPx_NoShift/F", buffersize);
   if( doWrite("met_shiftedPx_PhotonEnUp") ) tree->Branch("met_shiftedPx_PhotonEnUp", &met_shiftedPx_PhotonEnUp, "met_shiftedPx_PhotonEnUp/F", buffersize);
   if( doWrite("met_shiftedPx_PhotonEnDown") ) tree->Branch("met_shiftedPx_PhotonEnDown", &met_shiftedPx_PhotonEnDown, "met_shiftedPx_PhotonEnDown/F", buffersize);

   if( doWrite("met_shiftedPy_JetEnUp") ) tree->Branch("met_shiftedPy_JetEnUp", &met_shiftedPy_JetEnUp, "met_shiftedPy_JetEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_JetEnDown") ) tree->Branch("met_shiftedPy_JetEnDown", &met_shiftedPy_JetEnDown, "met_shiftedPy_JetEnDown/F", buffersize);
   if( doWrite("met_shiftedPy_JetResUp") ) tree->Branch("met_shiftedPy_JetResUp", &met_shiftedPy_JetResUp, "met_shiftedPy_JetResUp/F", buffersize);
   if( doWrite("met_shiftedPy_JetResDown") ) tree->Branch("met_shiftedPy_JetResDown", &met_shiftedPy_JetResDown, "met_shiftedPy_JetResDown/F", buffersize);
   if( doWrite("met_shiftedPy_MuonEnUp") ) tree->Branch("met_shiftedPy_MuonEnUp", &met_shiftedPy_MuonEnUp, "met_shiftedPy_MuonEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_MuonEnDown") ) tree->Branch("met_shiftedPy_MuonEnDown", &met_shiftedPy_MuonEnDown, "met_shiftedPy_MuonEnDown/F", buffersize);
   if( doWrite("met_shiftedPy_ElectronEnUp") ) tree->Branch("met_shiftedPy_ElectronEnUp", &met_shiftedPy_ElectronEnUp, "met_shiftedPy_ElectronEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_ElectronEnDown") ) tree->Branch("met_shiftedPy_ElectronEnDown", &met_shiftedPy_ElectronEnDown, "met_shiftedPy_ElectronEnDown/F", buffersize);
   if( doWrite("met_shiftedPy_TauEnUp") ) tree->Branch("met_shiftedPy_TauEnUp", &met_shiftedPy_TauEnUp, "met_shiftedPy_TauEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_TauEnDown") ) tree->Branch("met_shiftedPy_TauEnDown", &met_shiftedPy_TauEnDown, "met_shiftedPy_TauEnDown/F", buffersize);
   if( doWrite("met_shiftedPy_UnclusteredEnUp") ) tree->Branch("met_shiftedPy_UnclusteredEnUp", &met_shiftedPy_UnclusteredEnUp, "met_shiftedPy_UnclusteredEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_UnclusteredEnDown") ) tree->Branch("met_shiftedPy_UnclusteredEnDown", &met_shiftedPy_UnclusteredEnDown, "met_shiftedPy_UnclusteredEnDown/F", buffersize);
   if( doWrite("met_shiftedPy_NoShift") ) tree->Branch("met_shiftedPy_NoShift", &met_shiftedPy_NoShift, "met_shiftedPy_NoShift/F", buffersize);
   if( doWrite("met_shiftedPy_PhotonEnUp") ) tree->Branch("met_shiftedPy_PhotonEnUp", &met_shiftedPy_PhotonEnUp, "met_shiftedPy_PhotonEnUp/F", buffersize);
   if( doWrite("met_shiftedPy_PhotonEnDown") ) tree->Branch("met_shiftedPy_PhotonEnDown", &met_shiftedPy_PhotonEnDown, "met_shiftedPy_PhotonEnDown/F", buffersize);
   
   if( doWrite("met_shiftedPhi_JetEnUp") ) tree->Branch("met_shiftedPhi_JetEnUp", &met_shiftedPhi_JetEnUp, "met_shiftedPhi_JetEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_JetEnDown") ) tree->Branch("met_shiftedPhi_JetEnDown", &met_shiftedPhi_JetEnDown, "met_shiftedPhi_JetEnDown/F", buffersize);
   if( doWrite("met_shiftedPhi_JetResUp") ) tree->Branch("met_shiftedPhi_JetResUp", &met_shiftedPhi_JetResUp, "met_shiftedPhi_JetResUp/F", buffersize);
   if( doWrite("met_shiftedPhi_JetResDown") ) tree->Branch("met_shiftedPhi_JetResDown", &met_shiftedPhi_JetResDown, "met_shiftedPhi_JetResDown/F", buffersize);
   if( doWrite("met_shiftedPhi_MuonEnUp") ) tree->Branch("met_shiftedPhi_MuonEnUp", &met_shiftedPhi_MuonEnUp, "met_shiftedPhi_MuonEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_MuonEnDown") ) tree->Branch("met_shiftedPhi_MuonEnDown", &met_shiftedPhi_MuonEnDown, "met_shiftedPhi_MuonEnDown/F", buffersize);
   if( doWrite("met_shiftedPhi_ElectronEnUp") ) tree->Branch("met_shiftedPhi_ElectronEnUp", &met_shiftedPhi_ElectronEnUp, "met_shiftedPhi_ElectronEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_ElectronEnDown") ) tree->Branch("met_shiftedPhi_ElectronEnDown", &met_shiftedPhi_ElectronEnDown, "met_shiftedPhi_ElectronEnDown/F", buffersize);
   if( doWrite("met_shiftedPhi_TauEnUp") ) tree->Branch("met_shiftedPhi_TauEnUp", &met_shiftedPhi_TauEnUp, "met_shiftedPhi_TauEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_TauEnDown") ) tree->Branch("met_shiftedPhi_TauEnDown", &met_shiftedPhi_TauEnDown, "met_shiftedPhi_TauEnDown/F", buffersize);
   if( doWrite("met_shiftedPhi_UnclusteredEnUp") ) tree->Branch("met_shiftedPhi_UnclusteredEnUp", &met_shiftedPhi_UnclusteredEnUp, "met_shiftedPhi_UnclusteredEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_UnclusteredEnDown") ) tree->Branch("met_shiftedPhi_UnclusteredEnDown", &met_shiftedPhi_UnclusteredEnDown, "met_shiftedPhi_UnclusteredEnDown/F", buffersize);
   if( doWrite("met_shiftedPhi_NoShift") ) tree->Branch("met_shiftedPhi_NoShift", &met_shiftedPhi_NoShift, "met_shiftedPhi_NoShift/F", buffersize);
   if( doWrite("met_shiftedPhi_PhotonEnUp") ) tree->Branch("met_shiftedPhi_PhotonEnUp", &met_shiftedPhi_PhotonEnUp, "met_shiftedPhi_PhotonEnUp/F", buffersize);
   if( doWrite("met_shiftedPhi_PhotonEnDown") ) tree->Branch("met_shiftedPhi_PhotonEnDown", &met_shiftedPhi_PhotonEnDown, "met_shiftedPhi_PhotonEnDown/F", buffersize);

   if( doWrite("met_shiftedSumEt_JetEnUp") ) tree->Branch("met_shiftedSumEt_JetEnUp", &met_shiftedSumEt_JetEnUp, "met_shiftedSumEt_JetEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_JetEnDown") ) tree->Branch("met_shiftedSumEt_JetEnDown", &met_shiftedSumEt_JetEnDown, "met_shiftedSumEt_JetEnDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_JetResUp") ) tree->Branch("met_shiftedSumEt_JetResUp", &met_shiftedSumEt_JetResUp, "met_shiftedSumEt_JetResUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_JetResDown") ) tree->Branch("met_shiftedSumEt_JetResDown", &met_shiftedSumEt_JetResDown, "met_shiftedSumEt_JetResDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_MuonEnUp") ) tree->Branch("met_shiftedSumEt_MuonEnUp", &met_shiftedSumEt_MuonEnUp, "met_shiftedSumEt_MuonEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_MuonEnDown") ) tree->Branch("met_shiftedSumEt_MuonEnDown", &met_shiftedSumEt_MuonEnDown, "met_shiftedSumEt_MuonEnDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_ElectronEnUp") ) tree->Branch("met_shiftedSumEt_ElectronEnUp", &met_shiftedSumEt_ElectronEnUp, "met_shiftedSumEt_ElectronEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_ElectronEnDown") ) tree->Branch("met_shiftedSumEt_ElectronEnDown", &met_shiftedSumEt_ElectronEnDown, "met_shiftedSumEt_ElectronEnDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_TauEnUp") ) tree->Branch("met_shiftedSumEt_TauEnUp", &met_shiftedSumEt_TauEnUp, "met_shiftedSumEt_TauEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_TauEnDown") ) tree->Branch("met_shiftedSumEt_TauEnDown", &met_shiftedSumEt_TauEnDown, "met_shiftedSumEt_TauEnDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_UnclusteredEnUp") ) tree->Branch("met_shiftedSumEt_UnclusteredEnUp", &met_shiftedSumEt_UnclusteredEnUp, "met_shiftedSumEt_UnclusteredEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_UnclusteredEnDown") ) tree->Branch("met_shiftedSumEt_UnclusteredEnDown", &met_shiftedSumEt_UnclusteredEnDown, "met_shiftedSumEt_UnclusteredEnDown/F", buffersize);
   if( doWrite("met_shiftedSumEt_NoShift") ) tree->Branch("met_shiftedSumEt_NoShift", &met_shiftedSumEt_NoShift, "met_shiftedSumEt_NoShift/F", buffersize);
   if( doWrite("met_shiftedSumEt_PhotonEnUp") ) tree->Branch("met_shiftedSumEt_PhotonEnUp", &met_shiftedSumEt_PhotonEnUp, "met_shiftedSumEt_PhotonEnUp/F", buffersize);
   if( doWrite("met_shiftedSumEt_PhotonEnDown") ) tree->Branch("met_shiftedSumEt_PhotonEnDown", &met_shiftedSumEt_PhotonEnDown, "met_shiftedSumEt_PhotonEnDown/F", buffersize);
  
   if( doWrite("metNoHF_pt") )     tree->Branch("metNoHF_pt",    &metNoHF_pt,    "metNoHF_pt/F",    buffersize);
   if( doWrite("metNoHF_phi") )    tree->Branch("metNoHF_phi",   &metNoHF_phi,   "metNoHF_phi/F",   buffersize);
   if( doWrite("metNoHF_sumet") )  tree->Branch("metNoHF_sumet", &metNoHF_sumet, "metNoHF_sumet/F", buffersize);

   if( doWrite("metPuppi_pt") )    tree->Branch("metPuppi_pt",    &metPuppi_pt,    "metPuppi_pt/F",    buffersize);
   if( doWrite("metPuppi_phi") )   tree->Branch("metPuppi_phi",   &metPuppi_phi,   "metPuppi_phi/F",   buffersize);
   if( doWrite("metPuppi_sumet") ) tree->Branch("metPuppi_sumet", &metPuppi_sumet, "metPuppi_sumet/F", buffersize);
   
   if( doWrite("pv_n") ) tree->Branch("pv_n", &pv_n, "pv_n/I", buffersize);
   if( doWrite("pv_x") ) tree->Branch("pv_x", &pv_x, "pv_x/F", buffersize);
   if( doWrite("pv_y") ) tree->Branch("pv_y", &pv_y, "pv_y/F", buffersize);
   if( doWrite("pv_z") ) tree->Branch("pv_z", &pv_z, "pv_z/F", buffersize);

   if( doWrite("pv_xError") ) tree->Branch("pv_xError", &pv_xError, "pv_xError/F", buffersize);
   if( doWrite("pv_yError") ) tree->Branch("pv_yError", &pv_yError, "pv_yError/F", buffersize);
   if( doWrite("pv_zError") ) tree->Branch("pv_zError", &pv_zError, "pv_zError/F", buffersize);
   
   if( doWrite("pv_chi2") ) tree->Branch("pv_chi2", &pv_chi2, "pv_chi2/F", buffersize);
   if( doWrite("pv_ndof") ) tree->Branch("pv_ndof", &pv_ndof, "pv_ndof/I", buffersize);
   if( doWrite("pv_rho") ) tree->Branch("pv_rho", &pv_rho, "pv_rho/F", buffersize);
   if( doWrite("pv_isFake") ) tree->Branch("pv_isFake", &pv_isFake, "pv_isFake/I", buffersize);

   if( doWrite("mc_weight") ) tree->Branch("mc_weight", &mc_weight, "mc_weight/F", buffersize);
   if( doWrite("mc_weight_originalValue") ) tree->Branch("mc_weight_originalValue", &mc_weight_originalValue, "mc_weight_originalValue/F", buffersize);
   if( doWrite("mc_id") ) tree->Branch("mc_id", &mc_id, "mc_id/I", buffersize);
   if( doWrite("mc_f1") ) tree->Branch("mc_f1", &mc_f1, "mc_f1/I", buffersize);
   if( doWrite("mc_f2") ) tree->Branch("mc_f2", &mc_f2, "mc_f2/I", buffersize);
   if( doWrite("mc_x1") ) tree->Branch("mc_x1", &mc_x1, "mc_x1/F", buffersize);
   if( doWrite("mc_x2") ) tree->Branch("mc_x2", &mc_x2, "mc_x2/F", buffersize);
   if( doWrite("mc_scale") ) tree->Branch("mc_scale", &mc_scale, "mc_scale/F", buffersize);
   if( doWrite("mc_ptHat") ) tree->Branch("mc_ptHat", &mc_ptHat, "mc_ptHat/F", buffersize);
   
   if( doWrite("weight_originalXWGTUP") ) tree->Branch("weight_originalXWGTUP", &weight_originalXWGTUP, "weight_originalXWGTUP/F", buffersize);
   if( doWrite("weight_scale_muF0p5") ) tree->Branch("weight_scale_muF0p5", &weight_scale_muF0p5, "weight_scale_muF0p5/F", buffersize);
   if( doWrite("weight_scale_muF2"  ) ) tree->Branch("weight_scale_muF2",   &weight_scale_muF2,   "weight_scale_muF2/F", buffersize);
   if( doWrite("weight_scale_muR0p5") ) tree->Branch("weight_scale_muR0p5", &weight_scale_muR0p5, "weight_scale_muR0p5/F", buffersize);
   if( doWrite("weight_scale_muR2"  ) ) tree->Branch("weight_scale_muR2",   &weight_scale_muR2,   "weight_scale_muR2/F", buffersize);
   if( doWrite("mc_pdfweights") ) tree->Branch("mc_pdfweights", "std::vector<float>", &mc_pdfweights, buffersize);
   if( doWrite("mc_pdfweightIds") ) tree->Branch("mc_pdfweightIds", "std::vector<std::string>", &mc_pdfweightIds, buffersize);
   
   if( doWrite("mc_pu_intime_NumInt") ) tree->Branch("mc_pu_intime_NumInt", &mc_pu_intime_NumInt, "mc_pu_intime_NumInt/I", buffersize);
   if( doWrite("mc_pu_trueNumInt") ) tree->Branch("mc_pu_trueNumInt", &mc_pu_trueNumInt, "mc_pu_trueNumInt/I", buffersize);
   if( doWrite("mc_pu_before_npu") ) tree->Branch("mc_pu_before_npu", &mc_pu_before_npu, "mc_pu_before_npu/I", buffersize);
   if( doWrite("mc_pu_after_npu") ) tree->Branch("mc_pu_after_npu", &mc_pu_after_npu, "mc_pu_after_npu/I", buffersize);
   
   if( doWrite("mc_pu_Npvi") ) tree->Branch("mc_pu_Npvi", &mc_pu_Npvi, "mc_pu_Npvi/I", buffersize);
   if( doWrite("mc_pu_Nzpositions") ) tree->Branch("mc_pu_Nzpositions", "std::vector<int>", &mc_pu_Nzpositions, buffersize);
   if( doWrite("mc_pu_BunchCrossing") ) tree->Branch("mc_pu_BunchCrossing", "std::vector<int>", &mc_pu_BunchCrossing, buffersize );
   if( doWrite("mc_pu_zpositions") ) tree->Branch("mc_pu_zpositions", "std::vector<std::vector<float> >", &mc_pu_zpositions, buffersize);
   if( doWrite("mc_pu_sumpT_lowpT") ) tree->Branch("mc_pu_sumpT_lowpT", "std::vector<std::vector<float> >", &mc_pu_sumpT_lowpT, buffersize);
   if( doWrite("mc_pu_sumpT_highpT") ) tree->Branch("mc_pu_sumpT_highpT", "std::vector<std::vector<float> >", &mc_pu_sumpT_highpT, buffersize);
   if( doWrite("mc_pu_ntrks_lowpT") ) tree->Branch("mc_pu_ntrks_lowpT", "std::vector<std::vector<int> >", &mc_pu_ntrks_lowpT, buffersize);
   if( doWrite("mc_pu_ntrks_highpT") ) tree->Branch("mc_pu_ntrks_highpT", "std::vector<std::vector<int> >", &mc_pu_ntrks_highpT, buffersize);
 
   if( doWrite("trigger_n") ) tree->Branch("trigger_n", &trigger_n, "trigger_n/I", buffersize); 
   if( doWrite("trigger") ) tree->Branch("trigger", "std::vector<int>", &trigger, buffersize);
   if( doWrite("trigger_name") ) tree->Branch("trigger_name", "std::vector<string>", &trigger_name, buffersize);
   if( doWrite("trigger_pass") ) tree->Branch("trigger_pass", "std::vector<bool>", &trigger_pass, buffersize);
   if( doWrite("trigger_prescale") ) tree->Branch("trigger_prescale", "std::vector<int>", &trigger_prescale, buffersize);
   if( doWrite("trigger_HLTprescale") ) tree->Branch("trigger_HLTprescale", "std::vector<int>", &trigger_HLTprescale, buffersize);
   if( doWrite("trigger_L1prescale") ) tree->Branch("trigger_L1prescale", "std::vector<int>", &trigger_L1prescale, buffersize);
   
   if( doWrite("triggerobject_do") )
     {	
	if( doWrite("triggerobject_n") ) tree->Branch("triggerobject_n", &triggerobject_n, "triggerobject_n/I",buffersize);
	if( doWrite("triggerobject_pt") ) tree->Branch("triggerobject_pt", "std::vector<float>", &triggerobject_pt, buffersize);
	if( doWrite("triggerobject_eta") ) tree->Branch("triggerobject_eta", "std::vector<float>", &triggerobject_eta, buffersize);
	if( doWrite("triggerobject_phi") ) tree->Branch("triggerobject_phi", "std::vector<float>", &triggerobject_phi, buffersize);
	
	if( doWrite("triggerobject_collection") ) tree->Branch("triggerobject_collection", "std::vector<std::string>", &triggerobject_collection, buffersize);
	
	if( doWrite("triggerobject_filterIds_n") ) tree->Branch("triggerobject_filterIds_n", "std::vector<int>", &triggerobject_filterIds_n, buffersize);
	if( doWrite("triggerobject_filterIds") ) tree->Branch("triggerobject_filterIds", "std::vector<int>", &triggerobject_filterIds, buffersize);
	
	if( doWrite("triggerobject_isTriggerL1Mu") ) tree->Branch("triggerobject_isTriggerL1Mu", "std::vector<bool>", &triggerobject_isTriggerL1Mu, buffersize);
	if( doWrite("triggerobject_isTriggerL1NoIsoEG") ) tree->Branch("triggerobject_isTriggerL1NoIsoEG", "std::vector<bool>", &triggerobject_isTriggerL1NoIsoEG, buffersize);
	if( doWrite("triggerobject_isTriggerL1IsoEG") ) tree->Branch("triggerobject_isTriggerL1IsoEG", "std::vector<bool>", &triggerobject_isTriggerL1IsoEG, buffersize);
	if( doWrite("triggerobject_isTriggerL1CenJet") ) tree->Branch("triggerobject_isTriggerL1CenJet", "std::vector<bool>", &triggerobject_isTriggerL1CenJet, buffersize);
	if( doWrite("triggerobject_isTriggerL1ForJet") ) tree->Branch("triggerobject_isTriggerL1ForJet", "std::vector<bool>", &triggerobject_isTriggerL1ForJet, buffersize);
	if( doWrite("triggerobject_isTriggerL1TauJet") ) tree->Branch("triggerobject_isTriggerL1TauJet", "std::vector<bool>", &triggerobject_isTriggerL1TauJet, buffersize);
	if( doWrite("triggerobject_isTriggerL1ETM") ) tree->Branch("triggerobject_isTriggerL1ETM", "std::vector<bool>", &triggerobject_isTriggerL1ETM, buffersize);
	if( doWrite("triggerobject_isTriggerL1ETT") ) tree->Branch("triggerobject_isTriggerL1ETT", "std::vector<bool>", &triggerobject_isTriggerL1ETT, buffersize);
	if( doWrite("triggerobject_isTriggerL1HTT") ) tree->Branch("triggerobject_isTriggerL1HTT", "std::vector<bool>", &triggerobject_isTriggerL1HTT, buffersize);
	if( doWrite("triggerobject_isTriggerL1HTM") ) tree->Branch("triggerobject_isTriggerL1HTM", "std::vector<bool>", &triggerobject_isTriggerL1HTM, buffersize);
	if( doWrite("triggerobject_isTriggerL1JetCounts") ) tree->Branch("triggerobject_isTriggerL1JetCounts", "std::vector<bool>", &triggerobject_isTriggerL1JetCounts, buffersize);
	if( doWrite("triggerobject_isTriggerL1HfBitCounts") ) tree->Branch("triggerobject_isTriggerL1HfBitCounts", "std::vector<bool>", &triggerobject_isTriggerL1HfBitCounts, buffersize);
	if( doWrite("triggerobject_isTriggerL1HfRingEtSums") ) tree->Branch("triggerobject_isTriggerL1HfRingEtSums", "std::vector<bool>", &triggerobject_isTriggerL1HfRingEtSums, buffersize);
	if( doWrite("triggerobject_isTriggerL1TechTrig") ) tree->Branch("triggerobject_isTriggerL1TechTrig", "std::vector<bool>", &triggerobject_isTriggerL1TechTrig, buffersize);
	if( doWrite("triggerobject_isTriggerL1Castor") ) tree->Branch("triggerobject_isTriggerL1Castor", "std::vector<bool>", &triggerobject_isTriggerL1Castor, buffersize);
	if( doWrite("triggerobject_isTriggerL1BPTX") ) tree->Branch("triggerobject_isTriggerL1BPTX", "std::vector<bool>", &triggerobject_isTriggerL1BPTX, buffersize);
	if( doWrite("triggerobject_isTriggerL1GtExternal") ) tree->Branch("triggerobject_isTriggerL1GtExternal", "std::vector<bool>", &triggerobject_isTriggerL1GtExternal, buffersize);
	
	if( doWrite("triggerobject_isHLT_TriggerPhoton") ) tree->Branch("triggerobject_isHLT_TriggerPhoton", "std::vector<bool>", &triggerobject_isHLT_TriggerPhoton, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerElectron") ) tree->Branch("triggerobject_isHLT_TriggerElectron", "std::vector<bool>", &triggerobject_isHLT_TriggerElectron, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerMuon") ) tree->Branch("triggerobject_isHLT_TriggerMuon", "std::vector<bool>", &triggerobject_isHLT_TriggerMuon, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerTau") ) tree->Branch("triggerobject_isHLT_TriggerTau", "std::vector<bool>", &triggerobject_isHLT_TriggerTau, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerJet") ) tree->Branch("triggerobject_isHLT_TriggerJet", "std::vector<bool>", &triggerobject_isHLT_TriggerJet, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerBJet") ) tree->Branch("triggerobject_isHLT_TriggerBJet", "std::vector<bool>", &triggerobject_isHLT_TriggerBJet, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerMET") ) tree->Branch("triggerobject_isHLT_TriggerMET", "std::vector<bool>", &triggerobject_isHLT_TriggerMET, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerTET") ) tree->Branch("triggerobject_isHLT_TriggerTET", "std::vector<bool>", &triggerobject_isHLT_TriggerTET, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerTHT") ) tree->Branch("triggerobject_isHLT_TriggerTHT", "std::vector<bool>", &triggerobject_isHLT_TriggerTHT, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerMHT") ) tree->Branch("triggerobject_isHLT_TriggerMHT", "std::vector<bool>", &triggerobject_isHLT_TriggerMHT, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerTrack") ) tree->Branch("triggerobject_isHLT_TriggerTrack", "std::vector<bool>", &triggerobject_isHLT_TriggerTrack, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerCluster") ) tree->Branch("triggerobject_isHLT_TriggerCluster", "std::vector<bool>", &triggerobject_isHLT_TriggerCluster, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerMETSig") ) tree->Branch("triggerobject_isHLT_TriggerMETSig", "std::vector<bool>", &triggerobject_isHLT_TriggerMETSig, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerELongit") ) tree->Branch("triggerobject_isHLT_TriggerELongit", "std::vector<bool>", &triggerobject_isHLT_TriggerELongit, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerMHTSig") ) tree->Branch("triggerobject_isHLT_TriggerMHTSig", "std::vector<bool>", &triggerobject_isHLT_TriggerMHTSig, buffersize);
	if( doWrite("triggerobject_isHLT_TriggerHLongit") ) tree->Branch("triggerobject_isHLT_TriggerHLongit", "std::vector<bool>", &triggerobject_isHLT_TriggerHLongit, buffersize);
	
	if( doWrite("triggerobject_filterLabels_n") ) tree->Branch("triggerobject_filterLabels_n", "std::vector<int>", &triggerobject_filterLabels_n, buffersize);
	if( doWrite("triggerobject_filterLabels") ) tree->Branch("triggerobject_filterLabels", "std::vector<string>", &triggerobject_filterLabels, buffersize);
	
	if( doWrite("triggerobject_pathNamesAll_n") ) tree->Branch("triggerobject_pathNamesAll_n", "std::vector<int>", &triggerobject_pathNamesAll_n, buffersize);
	if( doWrite("triggerobject_pathNamesAll") ) tree->Branch("triggerobject_pathNamesAll", "std::vector<string>", &triggerobject_pathNamesAll, buffersize);
	if( doWrite("triggerobject_pathNamesAll_isBoth") ) tree->Branch("triggerobject_pathNamesAll_isBoth", "std::vector<bool>", &triggerobject_pathNamesAll_isBoth, buffersize);
	if( doWrite("triggerobject_pathNamesAll_isL3") ) tree->Branch("triggerobject_pathNamesAll_isL3", "std::vector<bool>", &triggerobject_pathNamesAll_isL3, buffersize);
	if( doWrite("triggerobject_pathNamesAll_isLF") ) tree->Branch("triggerobject_pathNamesAll_isLF", "std::vector<bool>", &triggerobject_pathNamesAll_isLF, buffersize);
	if( doWrite("triggerobject_pathNamesAll_isNone") ) tree->Branch("triggerobject_pathNamesAll_isNone", "std::vector<bool>", &triggerobject_pathNamesAll_isNone, buffersize);
     }
   
   if( doWrite("nvertex") ) tree->Branch("nvertex", &nvertex, "nvertex/I", buffersize);

   if( doWrite("el_n") ) tree->Branch("el_n", &el_n, "el_n/I", buffersize);
   if( doWrite("el_pt") ) tree->Branch("el_pt", "std::vector<float>", &el_pt, buffersize);
   if( doWrite("el_eta") ) tree->Branch("el_eta", "std::vector<float>", &el_eta, buffersize);
   if( doWrite("el_phi") ) tree->Branch("el_phi", "std::vector<float>", &el_phi, buffersize);
   if( doWrite("el_m") ) tree->Branch("el_m", "std::vector<float>", &el_m, buffersize);
   if( doWrite("el_E") ) tree->Branch("el_E", "std::vector<float>", &el_E, buffersize);
   if( doWrite("el_id") ) tree->Branch("el_id", "std::vector<int>", &el_id, buffersize);
   if( doWrite("el_charge") ) tree->Branch("el_charge", "std::vector<int>", &el_charge, buffersize);
   
   if( doWrite("el_passConversionVeto") ) tree->Branch("el_passConversionVeto", "std::vector<int>", &el_passConversionVeto, buffersize);   
   if( doWrite("el_isGsfCtfScPixChargeConsistent") ) tree->Branch("el_isGsfCtfScPixChargeConsistent", "std::vector<int>", &el_isGsfCtfScPixChargeConsistent, buffersize);
   if( doWrite("el_isGsfScPixChargeConsistent") ) tree->Branch("el_isGsfScPixChargeConsistent", "std::vector<int>", &el_isGsfScPixChargeConsistent, buffersize);

   if( doWrite("el_ecalEnergy") ) tree->Branch("el_ecalEnergy", "std::vector<float>", &el_ecalEnergy, buffersize);
   if( doWrite("el_correctedEcalEnergy") ) tree->Branch("el_correctedEcalEnergy", "std::vector<float>", &el_correctedEcalEnergy, buffersize);
   if( doWrite("el_correctedEcalEnergyError") ) tree->Branch("el_correctedEcalEnergyError", "std::vector<float>", &el_correctedEcalEnergyError, buffersize);
   if( doWrite("el_trackMomentumError") ) tree->Branch("el_trackMomentumError", "std::vector<float>", &el_trackMomentumError, buffersize);
   
   if( doWrite("el_ip3d") ) tree->Branch("el_ip3d", "std::vector<float>", &el_ip3d, buffersize);
   if( doWrite("el_ip3dErr") ) tree->Branch("el_ip3dErr", "std::vector<float>", &el_ip3dErr, buffersize);
   if( doWrite("el_ip2d") ) tree->Branch("el_ip2d", "std::vector<float>", &el_ip2d, buffersize);
   if( doWrite("el_ip2dErr") ) tree->Branch("el_ip2dErr", "std::vector<float>", &el_ip2dErr, buffersize);
   if( doWrite("el_ip3dBS") ) tree->Branch("el_ip3dBS", "std::vector<float>", &el_ip3dBS, buffersize);
   if( doWrite("el_ip3dBSErr") ) tree->Branch("el_ip3dBSErr", "std::vector<float>", &el_ip3dBSErr, buffersize);
   if( doWrite("el_ip2dBS") ) tree->Branch("el_ip2dBS", "std::vector<float>", &el_ip2dBS, buffersize);
   if( doWrite("el_ip2dBSErr") ) tree->Branch("el_ip2dBSErr", "std::vector<float>", &el_ip2dBSErr, buffersize);
   
   if( doWrite("el_neutralHadronIso") ) tree->Branch("el_neutralHadronIso", "std::vector<float>", &el_neutralHadronIso, buffersize);
   if( doWrite("el_chargedHadronIso") ) tree->Branch("el_chargedHadronIso", "std::vector<float>", &el_chargedHadronIso, buffersize);
   if( doWrite("el_puChargedHadronIso") ) tree->Branch("el_puChargedHadronIso", "std::vector<float>", &el_puChargedHadronIso, buffersize);
   if( doWrite("el_ecalIso") ) tree->Branch("el_ecalIso", "std::vector<float>", &el_ecalIso, buffersize);
   if( doWrite("el_hcalIso") ) tree->Branch("el_hcalIso", "std::vector<float>", &el_hcalIso, buffersize);
   if( doWrite("el_particleIso") ) tree->Branch("el_particleIso", "std::vector<float>", &el_particleIso, buffersize);
   if( doWrite("el_photonIso") ) tree->Branch("el_photonIso", "std::vector<float>", &el_photonIso, buffersize);
   if( doWrite("el_trackIso") ) tree->Branch("el_trackIso", "std::vector<float>", &el_trackIso, buffersize);
   
   if( doWrite("el_ecalPFClusterIso") ) tree->Branch("el_ecalPFClusterIso", "std::vector<float>", &el_ecalPFClusterIso, buffersize);
   if( doWrite("el_hcalPFClusterIso") ) tree->Branch("el_hcalPFClusterIso", "std::vector<float>", &el_hcalPFClusterIso, buffersize);

   if( doWrite("el_dr03EcalRecHitSumEt") ) tree->Branch("el_dr03EcalRecHitSumEt", "std::vector<float>", &el_dr03EcalRecHitSumEt, buffersize);
   if( doWrite("el_dr03HcalTowerSumEt") ) tree->Branch("el_dr03HcalTowerSumEt", "std::vector<float>", &el_dr03HcalTowerSumEt, buffersize);
   if( doWrite("el_dr03HcalDepth1TowerSumEt") ) tree->Branch("el_dr03HcalDepth1TowerSumEt", "std::vector<float>", &el_dr03HcalDepth1TowerSumEt, buffersize);
   if( doWrite("el_dr03HcalDepth2TowerSumEt") ) tree->Branch("el_dr03HcalDepth2TowerSumEt", "std::vector<float>", &el_dr03HcalDepth2TowerSumEt, buffersize);
   if( doWrite("el_dr03TkSumPt") ) tree->Branch("el_dr03TkSumPt", "std::vector<float>", &el_dr03TkSumPt, buffersize);

   if( doWrite("el_dr04EcalRecHitSumEt") ) tree->Branch("el_dr04EcalRecHitSumEt", "std::vector<float>", &el_dr04EcalRecHitSumEt, buffersize);
   if( doWrite("el_dr04HcalTowerSumEt") ) tree->Branch("el_dr04HcalTowerSumEt", "std::vector<float>", &el_dr04HcalTowerSumEt, buffersize);
   if( doWrite("el_dr04HcalDepth1TowerSumEt") ) tree->Branch("el_dr04HcalDepth1TowerSumEt", "std::vector<float>", &el_dr04HcalDepth1TowerSumEt, buffersize);
   if( doWrite("el_dr04HcalDepth2TowerSumEt") ) tree->Branch("el_dr04HcalDepth2TowerSumEt", "std::vector<float>", &el_dr04HcalDepth2TowerSumEt, buffersize);
   if( doWrite("el_dr04TkSumPt") ) tree->Branch("el_dr04TkSumPt", "std::vector<float>", &el_dr04TkSumPt, buffersize);
   
   if( doWrite("el_hcalOverEcal") ) tree->Branch("el_hcalOverEcal", "std::vector<float>", &el_hcalOverEcal, buffersize);
   if( doWrite("el_hcalOverEcalBc") ) tree->Branch("el_hcalOverEcalBc", "std::vector<float>", &el_hcalOverEcalBc, buffersize);
   if( doWrite("el_hcalDepth1OverEcal") ) tree->Branch("el_hcalDepth1OverEcal", "std::vector<float>", &el_hcalDepth1OverEcal, buffersize);
   if( doWrite("el_hcalDepth2OverEcal") ) tree->Branch("el_hcalDepth2OverEcal", "std::vector<float>", &el_hcalDepth2OverEcal, buffersize);
   if( doWrite("el_eSeedClusterOverPout") ) tree->Branch("el_eSeedClusterOverPout", "std::vector<float>", &el_eSeedClusterOverPout, buffersize);
   if( doWrite("el_eSeedClusterOverP") ) tree->Branch("el_eSeedClusterOverP", "std::vector<float>", &el_eSeedClusterOverP, buffersize);
   if( doWrite("el_eEleClusterOverPout") ) tree->Branch("el_eEleClusterOverPout", "std::vector<float>", &el_eEleClusterOverPout, buffersize);
   if( doWrite("el_deltaEtaEleClusterTrackAtCalo") ) tree->Branch("el_deltaEtaEleClusterTrackAtCalo", "std::vector<float>", &el_deltaEtaEleClusterTrackAtCalo, buffersize);
   if( doWrite("el_deltaPhiEleClusterTrackAtCalo") ) tree->Branch("el_deltaPhiEleClusterTrackAtCalo", "std::vector<float>", &el_deltaPhiEleClusterTrackAtCalo, buffersize);
   
   if( doWrite("el_pfIso_sumChargedHadronPt") ) tree->Branch("el_pfIso_sumChargedHadronPt", "std::vector<float>", &el_pfIso_sumChargedHadronPt, buffersize);
   if( doWrite("el_pfIso_sumNeutralHadronEt") ) tree->Branch("el_pfIso_sumNeutralHadronEt", "std::vector<float>", &el_pfIso_sumNeutralHadronEt, buffersize);
   if( doWrite("el_pfIso_sumPhotonEt") ) tree->Branch("el_pfIso_sumPhotonEt", "std::vector<float>", &el_pfIso_sumPhotonEt, buffersize);
   if( doWrite("el_pfIso_sumPUPt") ) tree->Branch("el_pfIso_sumPUPt", "std::vector<float>", &el_pfIso_sumPUPt, buffersize);
   
   if( doWrite("el_miniIso") ) tree->Branch("el_miniIso", "std::vector<float>", &el_miniIso, buffersize);
   if( doWrite("el_miniIsoTTH") ) tree->Branch("el_miniIsoTTH", "std::vector<float>", &el_miniIsoTTH, buffersize);
   
   if( doWrite("el_vx") ) tree->Branch("el_vx", "std::vector<float>", &el_vx, buffersize);
   if( doWrite("el_vy") ) tree->Branch("el_vy", "std::vector<float>", &el_vy, buffersize);
   if( doWrite("el_vz") ) tree->Branch("el_vz", "std::vector<float>", &el_vz, buffersize);
   
   if( doWrite("el_hasGsfTrack") ) tree->Branch("el_hasGsfTrack", "std::vector<bool>", &el_hasGsfTrack, buffersize);
   if( doWrite("el_gsfTrack_d0") ) tree->Branch("el_gsfTrack_d0", "std::vector<float>", &el_gsfTrack_d0, buffersize);
   if( doWrite("el_gsfTrack_z0") ) tree->Branch("el_gsfTrack_z0", "std::vector<float>", &el_gsfTrack_z0, buffersize);
   if( doWrite("el_gsfTrack_d0Error") ) tree->Branch("el_gsfTrack_d0Error", "std::vector<float>", &el_gsfTrack_d0Error, buffersize);
   if( doWrite("el_gsfTrack_z0Error") ) tree->Branch("el_gsfTrack_z0Error", "std::vector<float>", &el_gsfTrack_z0Error, buffersize);
   if( doWrite("el_gsfTrack_PV_dxy") ) tree->Branch("el_gsfTrack_PV_dxy", "std::vector<float>", &el_gsfTrack_PV_dxy, buffersize);
   if( doWrite("el_gsfTrack_PV_dz") ) tree->Branch("el_gsfTrack_PV_dz", "std::vector<float>", &el_gsfTrack_PV_dz, buffersize);
   if( doWrite("el_gsfTrack_RP_dxy") ) tree->Branch("el_gsfTrack_RP_dxy", "std::vector<float>", &el_gsfTrack_RP_dxy, buffersize);
   if( doWrite("el_gsfTrack_RP_dz") ) tree->Branch("el_gsfTrack_RP_dz", "std::vector<float>", &el_gsfTrack_RP_dz, buffersize);
   if( doWrite("el_gsfTrack_BS_dxy") ) tree->Branch("el_gsfTrack_BS_dxy", "std::vector<float>", &el_gsfTrack_BS_dxy, buffersize);
   if( doWrite("el_gsfTrack_BS_dz") ) tree->Branch("el_gsfTrack_BS_dz", "std::vector<float>", &el_gsfTrack_BS_dz, buffersize);
   if( doWrite("el_gsfTrack_dxyError") ) tree->Branch("el_gsfTrack_dxyError", "std::vector<float>", &el_gsfTrack_dxyError, buffersize);
   if( doWrite("el_gsfTrack_dzError") ) tree->Branch("el_gsfTrack_dzError", "std::vector<float>", &el_gsfTrack_dzError, buffersize);
   if( doWrite("el_gsfTrack_normalizedChi2") ) tree->Branch("el_gsfTrack_normalizedChi2", "std::vector<float>", &el_gsfTrack_normalizedChi2, buffersize);
   
   if( doWrite("el_superCluster_eta") ) tree->Branch("el_superCluster_eta", "std::vector<float>", &el_superCluster_eta, buffersize);
   if( doWrite("el_superCluster_phi") ) tree->Branch("el_superCluster_phi", "std::vector<float>", &el_superCluster_phi, buffersize);
   if( doWrite("el_superCluster_energy") ) tree->Branch("el_superCluster_energy", "std::vector<float>", &el_superCluster_energy, buffersize);
   if( doWrite("el_superCluster_rawEnergy") ) tree->Branch("el_superCluster_rawEnergy", "std::vector<float>", &el_superCluster_rawEnergy, buffersize);
   if( doWrite("el_superCluster_preshowerEnergy") ) tree->Branch("el_superCluster_preshowerEnergy", "std::vector<float>", &el_superCluster_preshowerEnergy, buffersize);
   if( doWrite("el_superCluster_etaWidth") ) tree->Branch("el_superCluster_etaWidth", "std::vector<float>", &el_superCluster_etaWidth, buffersize);
   if( doWrite("el_superCluster_phiWidth") ) tree->Branch("el_superCluster_phiWidth", "std::vector<float>", &el_superCluster_phiWidth, buffersize);
   if( doWrite("el_superCluster_preshowerEnergyPlane1") ) tree->Branch("el_superCluster_preshowerEnergyPlane1", "std::vector<float>", &el_superCluster_preshowerEnergyPlane1, buffersize);
   if( doWrite("el_superCluster_preshowerEnergyPlane2") ) tree->Branch("el_superCluster_preshowerEnergyPlane2", "std::vector<float>", &el_superCluster_preshowerEnergyPlane2, buffersize);
   if( doWrite("el_superCluster_positionR") ) tree->Branch("el_superCluster_positionR", "std::vector<float>", &el_superCluster_positionR, buffersize);
   
   if( doWrite("el_basicClustersSize") ) tree->Branch("el_basicClustersSize", "std::vector<int>", &el_basicClustersSize, buffersize);
   if( doWrite("el_e1x5") ) tree->Branch("el_e1x5", "std::vector<float>", &el_e1x5, buffersize);
   if( doWrite("el_e5x5") ) tree->Branch("el_e5x5", "std::vector<float>", &el_e5x5, buffersize);
   if( doWrite("el_e2x5Max") ) tree->Branch("el_e2x5Max", "std::vector<float>", &el_e2x5Max, buffersize);
   if( doWrite("el_sigmaEtaEta") ) tree->Branch("el_sigmaEtaEta", "std::vector<float>", &el_sigmaEtaEta, buffersize);
   if( doWrite("el_sigmaIetaIeta") ) tree->Branch("el_sigmaIetaIeta", "std::vector<float>", &el_sigmaIetaIeta, buffersize);
   if( doWrite("el_sigmaIphiIphi") ) tree->Branch("el_sigmaIphiIphi", "std::vector<float>", &el_sigmaIphiIphi, buffersize);
   if( doWrite("el_sigmaIetaIphi") ) tree->Branch("el_sigmaIetaIphi", "std::vector<float>", &el_sigmaIetaIphi, buffersize);
   if( doWrite("el_full5x5_sigmaIphiIphi") ) tree->Branch("el_full5x5_sigmaIphiIphi", "std::vector<float>", &el_full5x5_sigmaIphiIphi, buffersize);
   if( doWrite("el_full5x5_sigmaEtaEta") ) tree->Branch("el_full5x5_sigmaEtaEta", "std::vector<float>", &el_full5x5_sigmaEtaEta, buffersize);
   if( doWrite("el_full5x5_sigmaIetaIeta") ) tree->Branch("el_full5x5_sigmaIetaIeta", "std::vector<float>", &el_full5x5_sigmaIetaIeta, buffersize);
   if( doWrite("el_full5x5_sigmaIetaIphi") ) tree->Branch("el_full5x5_sigmaIetaIphi", "std::vector<float>", &el_full5x5_sigmaIetaIphi, buffersize);
   if( doWrite("el_full5x5_r9") ) tree->Branch("el_full5x5_r9", "std::vector<float>", &el_full5x5_r9, buffersize);
   if( doWrite("el_full5x5_e1x5") ) tree->Branch("el_full5x5_e1x5", "std::vector<float>", &el_full5x5_e1x5, buffersize);
   if( doWrite("el_full5x5_e5x5") ) tree->Branch("el_full5x5_e5x5", "std::vector<float>", &el_full5x5_e5x5, buffersize);
   if( doWrite("el_full5x5_e2x5Max") ) tree->Branch("el_full5x5_e2x5Max", "std::vector<float>", &el_full5x5_e2x5Max, buffersize);

   if( doWrite("el_expectedMissingOuterHits") ) tree->Branch("el_expectedMissingOuterHits", "std::vector<int>", &el_expectedMissingOuterHits, buffersize);
   if( doWrite("el_numberOfValidPixelHits") ) tree->Branch("el_numberOfValidPixelHits", "std::vector<int>", &el_numberOfValidPixelHits, buffersize);
   if( doWrite("el_numberOfLostPixelHits") ) tree->Branch("el_numberOfLostPixelHits", "std::vector<int>", &el_numberOfLostPixelHits, buffersize);
   if( doWrite("el_trackerLayersWithMeasurement") ) tree->Branch("el_trackerLayersWithMeasurement", "std::vector<int>", &el_trackerLayersWithMeasurement, buffersize);
   if( doWrite("el_pixelLayersWithMeasurement") ) tree->Branch("el_pixelLayersWithMeasurement", "std::vector<int>", &el_pixelLayersWithMeasurement, buffersize);
   if( doWrite("el_numberOfValidStripLayersWithMonoAndStereo") ) tree->Branch("el_numberOfValidStripLayersWithMonoAndStereo", "std::vector<int>", &el_numberOfValidStripLayersWithMonoAndStereo, buffersize);
   if( doWrite("el_trackerLayersWithoutMeasurement") ) tree->Branch("el_trackerLayersWithoutMeasurement", "std::vector<int>", &el_trackerLayersWithoutMeasurement, buffersize);
   
   if( doWrite("el_numberOfHits") ) tree->Branch("el_numberOfHits", "std::vector<int>", &el_numberOfHits, buffersize);
   if( doWrite("el_numberOfValidHits") ) tree->Branch("el_numberOfValidHits", "std::vector<int>", &el_numberOfValidHits, buffersize);
   
   if( doWrite("el_hadronicOverEm") ) tree->Branch("el_hadronicOverEm", "std::vector<float>", &el_hadronicOverEm, buffersize);
   if( doWrite("el_numberOfLostHits") ) tree->Branch("el_numberOfLostHits", "std::vector<int>", &el_numberOfLostHits, buffersize);
   if( doWrite("el_numberOfLostHitsDefault") ) tree->Branch("el_numberOfLostHitsDefault", "std::vector<int>", &el_numberOfLostHitsDefault, buffersize);

   if( doWrite("el_fbrem") ) tree->Branch("el_fbrem", "std::vector<float>", &el_fbrem, buffersize);
   if( doWrite("el_kf_normalizedChi2") ) tree->Branch("el_kf_normalizedChi2", "std::vector<float>", &el_kf_normalizedChi2, buffersize);
   if( doWrite("el_gsf_normalizedChi2") ) tree->Branch("el_gsf_normalizedChi2", "std::vector<float>", &el_gsf_normalizedChi2, buffersize);
   if( doWrite("el_deltaEtaSuperClusterTrackAtVtx") ) tree->Branch("el_deltaEtaSuperClusterTrackAtVtx", "std::vector<float>", &el_deltaEtaSuperClusterTrackAtVtx, buffersize);
   if( doWrite("el_deltaPhiSuperClusterTrackAtVtx") ) tree->Branch("el_deltaPhiSuperClusterTrackAtVtx", "std::vector<float>", &el_deltaPhiSuperClusterTrackAtVtx, buffersize);
   if( doWrite("el_deltaEtaSeedClusterTrackAtCalo") ) tree->Branch("el_deltaEtaSeedClusterTrackAtCalo", "std::vector<float>", &el_deltaEtaSeedClusterTrackAtCalo, buffersize);
   if( doWrite("el_deltaPhiSeedClusterTrackAtCalo") ) tree->Branch("el_deltaPhiSeedClusterTrackAtCalo", "std::vector<float>", &el_deltaPhiSeedClusterTrackAtCalo, buffersize);
   if( doWrite("el_full5x5_OneMinusE1x5E5x5") ) tree->Branch("el_full5x5_OneMinusE1x5E5x5", "std::vector<float>", &el_full5x5_OneMinusE1x5E5x5, buffersize);
   if( doWrite("el_OneMinusE1x5E5x5") ) tree->Branch("el_OneMinusE1x5E5x5", "std::vector<float>", &el_OneMinusE1x5E5x5, buffersize);
   if( doWrite("el_eSuperClusterOverP") ) tree->Branch("el_eSuperClusterOverP", "std::vector<float>", &el_eSuperClusterOverP, buffersize);
   if( doWrite("el_IoEmIoP") ) tree->Branch("el_IoEmIoP", "std::vector<float>", &el_IoEmIoP, buffersize);
   if( doWrite("el_ooEmooP") ) tree->Branch("el_ooEmooP", "std::vector<float>", &el_ooEmooP, buffersize);
   if( doWrite("el_eleEoPout") ) tree->Branch("el_eleEoPout", "std::vector<float>", &el_eleEoPout, buffersize);
   if( doWrite("el_PreShowerOverRaw") ) tree->Branch("el_PreShowerOverRaw", "std::vector<float>", &el_PreShowerOverRaw, buffersize);

   if( doWrite("el_mvaIso") ) tree->Branch("el_mvaIso", "std::vector<float>", &el_mvaIso, buffersize);
   if( doWrite("el_mvaNoIso") ) tree->Branch("el_mvaNoIso", "std::vector<float>", &el_mvaNoIso, buffersize);
   
   if( doWrite("el_vetoCBId") ) tree->Branch("el_vetoCBId", "std::vector<bool>", &el_vetoCBId, buffersize);
   if( doWrite("el_looseCBId") ) tree->Branch("el_looseCBId", "std::vector<bool>", &el_looseCBId, buffersize);
   if( doWrite("el_mediumCBId") ) tree->Branch("el_mediumCBId", "std::vector<bool>", &el_mediumCBId, buffersize);
   if( doWrite("el_tightCBId") ) tree->Branch("el_tightCBId", "std::vector<bool>", &el_tightCBId, buffersize);

//   if( doWrite("el_vetoStopID") ) tree->Branch("el_vetoStopID", "std::vector<bool>", &el_vetoStopID, buffersize);
//   if( doWrite("el_mediumStopID") ) tree->Branch("el_mediumStopID", "std::vector<bool>", &el_mediumStopID, buffersize);

   if( doWrite("el_NoIso90MVAId") ) tree->Branch("el_NoIso90MVAId", "std::vector<bool>", &el_NoIso90MVAId, buffersize);
   if( doWrite("el_NoIso80MVAId") ) tree->Branch("el_NoIso80MVAId", "std::vector<bool>", &el_NoIso80MVAId, buffersize);
   if( doWrite("el_NoIsoLooseMVAId") ) tree->Branch("el_NoIsoLooseMVAId", "std::vector<bool>", &el_NoIsoLooseMVAId, buffersize);

   if( doWrite("el_Iso90MVAId") ) tree->Branch("el_Iso90MVAId", "std::vector<bool>", &el_Iso90MVAId, buffersize);
   if( doWrite("el_Iso80MVAId") ) tree->Branch("el_Iso80MVAId", "std::vector<bool>", &el_Iso80MVAId, buffersize);
   if( doWrite("el_IsoLooseMVAId") ) tree->Branch("el_IsoLooseMVAId", "std::vector<bool>", &el_IsoLooseMVAId, buffersize);
   
   if( doWrite("el_lepMVA") ) tree->Branch("el_lepMVA", "std::vector<float>", &el_lepMVA, buffersize);

   if( doWrite("el_lepMVA_pt") ) tree->Branch("el_lepMVA_pt", "std::vector<float>", &el_lepMVA_pt, buffersize);
   if( doWrite("el_lepMVA_eta") ) tree->Branch("el_lepMVA_eta", "std::vector<float>", &el_lepMVA_eta, buffersize);
   if( doWrite("el_lepMVA_miniRelIsoCharged") ) tree->Branch("el_lepMVA_miniRelIsoCharged", "std::vector<float>", &el_lepMVA_miniRelIsoCharged, buffersize);
   if( doWrite("el_lepMVA_miniRelIsoNeutral") ) tree->Branch("el_lepMVA_miniRelIsoNeutral", "std::vector<float>", &el_lepMVA_miniRelIsoNeutral, buffersize);
   if( doWrite("el_lepMVA_jetPtRatio") ) tree->Branch("el_lepMVA_jetPtRatio", "std::vector<float>", &el_lepMVA_jetPtRatio, buffersize);
   if( doWrite("el_lepMVA_jetPtRelv2") ) tree->Branch("el_lepMVA_jetPtRelv2", "std::vector<float>", &el_lepMVA_jetPtRelv2, buffersize);
   if( doWrite("el_lepMVA_jetBTagCSV") ) tree->Branch("el_lepMVA_jetBTagCSV", "std::vector<float>", &el_lepMVA_jetBTagCSV, buffersize);
   if( doWrite("el_lepMVA_sip3d") ) tree->Branch("el_lepMVA_sip3d", "std::vector<float>", &el_lepMVA_sip3d, buffersize);
   if( doWrite("el_lepMVA_dxy") ) tree->Branch("el_lepMVA_dxy", "std::vector<float>", &el_lepMVA_dxy, buffersize);
   if( doWrite("el_lepMVA_dz") ) tree->Branch("el_lepMVA_dz", "std::vector<float>", &el_lepMVA_dz, buffersize);
   if( doWrite("el_lepMVA_mvaId") ) tree->Branch("el_lepMVA_mvaId", "std::vector<float>", &el_lepMVA_mvaId, buffersize);
   if( doWrite("el_lepMVA_jetNDauChargedMVASel") ) tree->Branch("el_lepMVA_jetNDauChargedMVASel", "std::vector<float>", &el_lepMVA_jetNDauChargedMVASel, buffersize);

   if( doWrite("el_conept") ) tree->Branch("el_conept", "std::vector<float>", &el_conept, buffersize);

   if( doWrite("el_hasMCMatch") ) tree->Branch("el_hasMCMatch", "std::vector<int>", &el_hasMCMatch, buffersize);
   if( doWrite("el_gen_pt") ) tree->Branch("el_gen_pt", "std::vector<float>", &el_gen_pt, buffersize);
   if( doWrite("el_gen_eta") ) tree->Branch("el_gen_eta", "std::vector<float>", &el_gen_eta, buffersize);
   if( doWrite("el_gen_phi") ) tree->Branch("el_gen_phi", "std::vector<float>", &el_gen_phi, buffersize);
   if( doWrite("el_gen_m") ) tree->Branch("el_gen_m", "std::vector<float>", &el_gen_m, buffersize);
   if( doWrite("el_gen_status") ) tree->Branch("el_gen_status", "std::vector<int>", &el_gen_status, buffersize);
   if( doWrite("el_gen_id") ) tree->Branch("el_gen_id", "std::vector<int>", &el_gen_id, buffersize);
   if( doWrite("el_gen_charge") ) tree->Branch("el_gen_charge", "std::vector<int>", &el_gen_charge, buffersize);
   if( doWrite("el_gen_dr") ) tree->Branch("el_gen_dr", "std::vector<float>", &el_gen_dr, buffersize);

   if( doWrite("el_hasMCMatchPAT") ) tree->Branch("el_hasMCMatchPAT", "std::vector<int>", &el_hasMCMatchPAT, buffersize);
   if( doWrite("el_genPAT_pt") ) tree->Branch("el_genPAT_pt", "std::vector<float>", &el_genPAT_pt, buffersize);
   if( doWrite("el_genPAT_eta") ) tree->Branch("el_genPAT_eta", "std::vector<float>", &el_genPAT_eta, buffersize);
   if( doWrite("el_genPAT_phi") ) tree->Branch("el_genPAT_phi", "std::vector<float>", &el_genPAT_phi, buffersize);
   if( doWrite("el_genPAT_m") ) tree->Branch("el_genPAT_m", "std::vector<float>", &el_genPAT_m, buffersize);
   if( doWrite("el_genPAT_status") ) tree->Branch("el_genPAT_status", "std::vector<int>", &el_genPAT_status, buffersize);
   if( doWrite("el_genPAT_id") ) tree->Branch("el_genPAT_id", "std::vector<int>", &el_genPAT_id, buffersize);
   if( doWrite("el_genPAT_charge") ) tree->Branch("el_genPAT_charge", "std::vector<int>", &el_genPAT_charge, buffersize);
   
   if( doWrite("el_hasMatchedConversion") ) tree->Branch("el_hasMatchedConversion", "std::vector<bool>", &el_hasMatchedConversion, buffersize);
   if( doWrite("el_expectedMissingInnerHits") ) tree->Branch("el_expectedMissingInnerHits", "std::vector<int>", &el_expectedMissingInnerHits, buffersize);

   if( doWrite("mu_n") ) tree->Branch("mu_n", &mu_n, "mu_n/I", buffersize);
   if( doWrite("mu_pt") ) tree->Branch("mu_pt", "std::vector<float>", &mu_pt, buffersize);
   if( doWrite("mu_eta") ) tree->Branch("mu_eta", "std::vector<float>", &mu_eta, buffersize);
   if( doWrite("mu_phi") ) tree->Branch("mu_phi", "std::vector<float>", &mu_phi, buffersize);
   if( doWrite("mu_m") ) tree->Branch("mu_m", "std::vector<float>", &mu_m, buffersize);
   if( doWrite("mu_E") ) tree->Branch("mu_E", "std::vector<float>", &mu_E, buffersize);
   if( doWrite("mu_id") ) tree->Branch("mu_id", "std::vector<int>", &mu_id, buffersize);
   if( doWrite("mu_charge") ) tree->Branch("mu_charge", "std::vector<int>", &mu_charge, buffersize);

   if( doWrite("mu_ip3d") ) tree->Branch("mu_ip3d", "std::vector<float>", &mu_ip3d, buffersize);
   if( doWrite("mu_ip3dErr") ) tree->Branch("mu_ip3dErr", "std::vector<float>", &mu_ip3dErr, buffersize);
   if( doWrite("mu_ip2d") ) tree->Branch("mu_ip2d", "std::vector<float>", &mu_ip2d, buffersize);
   if( doWrite("mu_ip2dErr") ) tree->Branch("mu_ip2dErr", "std::vector<float>", &mu_ip2dErr, buffersize);
   if( doWrite("mu_ip3dBS") ) tree->Branch("mu_ip3dBS", "std::vector<float>", &mu_ip3dBS, buffersize);
   if( doWrite("mu_ip3dBSErr") ) tree->Branch("mu_ip3dBSErr", "std::vector<float>", &mu_ip3dBSErr, buffersize);
   if( doWrite("mu_ip2dBS") ) tree->Branch("mu_ip2dBS", "std::vector<float>", &mu_ip2dBS, buffersize);
   if( doWrite("mu_ip2dBSErr") ) tree->Branch("mu_ip2dBSErr", "std::vector<float>", &mu_ip2dBSErr, buffersize);
   
   if( doWrite("mu_neutralHadronIso") ) tree->Branch("mu_neutralHadronIso", "std::vector<float>", &mu_neutralHadronIso, buffersize);
   if( doWrite("mu_chargedHadronIso") ) tree->Branch("mu_chargedHadronIso", "std::vector<float>", &mu_chargedHadronIso, buffersize);
   if( doWrite("mu_puChargedHadronIso") ) tree->Branch("mu_puChargedHadronIso", "std::vector<float>", &mu_puChargedHadronIso, buffersize);
   if( doWrite("mu_ecalIso") ) tree->Branch("mu_ecalIso", "std::vector<float>", &mu_ecalIso, buffersize);
   if( doWrite("mu_hcalIso") ) tree->Branch("mu_hcalIso", "std::vector<float>", &mu_hcalIso, buffersize);
   if( doWrite("mu_photonIso") ) tree->Branch("mu_photonIso", "std::vector<float>", &mu_photonIso, buffersize);
   if( doWrite("mu_trackIso") ) tree->Branch("mu_trackIso", "std::vector<float>", &mu_trackIso, buffersize);

   if( doWrite("mu_pfIso03_sumChargedHadronPt") ) tree->Branch("mu_pfIso03_sumChargedHadronPt", "std::vector<float>", &mu_pfIso03_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfIso03_sumChargedParticlePt") ) tree->Branch("mu_pfIso03_sumChargedParticlePt", "std::vector<float>", &mu_pfIso03_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfIso03_sumNeutralHadronEt") ) tree->Branch("mu_pfIso03_sumNeutralHadronEt", "std::vector<float>", &mu_pfIso03_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfIso03_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfIso03_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfIso03_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfIso03_sumPhotonEt") ) tree->Branch("mu_pfIso03_sumPhotonEt", "std::vector<float>", &mu_pfIso03_sumPhotonEt, buffersize);
   if( doWrite("mu_pfIso03_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfIso03_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfIso03_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfIso03_sumPUPt") ) tree->Branch("mu_pfIso03_sumPUPt", "std::vector<float>", &mu_pfIso03_sumPUPt, buffersize);

   if( doWrite("mu_pfIso04_sumChargedHadronPt") ) tree->Branch("mu_pfIso04_sumChargedHadronPt", "std::vector<float>", &mu_pfIso04_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfIso04_sumChargedParticlePt") ) tree->Branch("mu_pfIso04_sumChargedParticlePt", "std::vector<float>", &mu_pfIso04_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfIso04_sumNeutralHadronEt") ) tree->Branch("mu_pfIso04_sumNeutralHadronEt", "std::vector<float>", &mu_pfIso04_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfIso04_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfIso04_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfIso04_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfIso04_sumPhotonEt") ) tree->Branch("mu_pfIso04_sumPhotonEt", "std::vector<float>", &mu_pfIso04_sumPhotonEt, buffersize);
   if( doWrite("mu_pfIso04_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfIso04_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfIso04_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfIso04_sumPUPt") ) tree->Branch("mu_pfIso04_sumPUPt", "std::vector<float>", &mu_pfIso04_sumPUPt, buffersize);

   if( doWrite("mu_pfMeanIso03_sumChargedHadronPt") ) tree->Branch("mu_pfMeanIso03_sumChargedHadronPt", "std::vector<float>", &mu_pfMeanIso03_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfMeanIso03_sumChargedParticlePt") ) tree->Branch("mu_pfMeanIso03_sumChargedParticlePt", "std::vector<float>", &mu_pfMeanIso03_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfMeanIso03_sumNeutralHadronEt") ) tree->Branch("mu_pfMeanIso03_sumNeutralHadronEt", "std::vector<float>", &mu_pfMeanIso03_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfMeanIso03_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfMeanIso03_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfMeanIso03_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfMeanIso03_sumPhotonEt") ) tree->Branch("mu_pfMeanIso03_sumPhotonEt", "std::vector<float>", &mu_pfMeanIso03_sumPhotonEt, buffersize);
   if( doWrite("mu_pfMeanIso03_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfMeanIso03_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfMeanIso03_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfMeanIso03_sumPUPt") ) tree->Branch("mu_pfMeanIso03_sumPUPt", "std::vector<float>", &mu_pfMeanIso03_sumPUPt, buffersize);

   if( doWrite("mu_pfSumIso03_sumChargedHadronPt") ) tree->Branch("mu_pfSumIso03_sumChargedHadronPt", "std::vector<float>", &mu_pfSumIso03_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfSumIso03_sumChargedParticlePt") ) tree->Branch("mu_pfSumIso03_sumChargedParticlePt", "std::vector<float>", &mu_pfSumIso03_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfSumIso03_sumNeutralHadronEt") ) tree->Branch("mu_pfSumIso03_sumNeutralHadronEt", "std::vector<float>", &mu_pfSumIso03_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfSumIso03_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfSumIso03_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfSumIso03_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfSumIso03_sumPhotonEt") ) tree->Branch("mu_pfSumIso03_sumPhotonEt", "std::vector<float>", &mu_pfSumIso03_sumPhotonEt, buffersize);
   if( doWrite("mu_pfSumIso03_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfSumIso03_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfSumIso03_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfSumIso03_sumPUPt") ) tree->Branch("mu_pfSumIso03_sumPUPt", "std::vector<float>", &mu_pfSumIso03_sumPUPt, buffersize);

   if( doWrite("mu_pfMeanIso04_sumChargedHadronPt") ) tree->Branch("mu_pfMeanIso04_sumChargedHadronPt", "std::vector<float>", &mu_pfMeanIso04_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfMeanIso04_sumChargedParticlePt") ) tree->Branch("mu_pfMeanIso04_sumChargedParticlePt", "std::vector<float>", &mu_pfMeanIso04_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfMeanIso04_sumNeutralHadronEt") ) tree->Branch("mu_pfMeanIso04_sumNeutralHadronEt", "std::vector<float>", &mu_pfMeanIso04_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfMeanIso04_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfMeanIso04_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfMeanIso04_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfMeanIso04_sumPhotonEt") ) tree->Branch("mu_pfMeanIso04_sumPhotonEt", "std::vector<float>", &mu_pfMeanIso04_sumPhotonEt, buffersize);
   if( doWrite("mu_pfMeanIso04_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfMeanIso04_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfMeanIso04_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfMeanIso04_sumPUPt") ) tree->Branch("mu_pfMeanIso04_sumPUPt", "std::vector<float>", &mu_pfMeanIso04_sumPUPt, buffersize);

   if( doWrite("mu_pfSumIso04_sumChargedHadronPt") ) tree->Branch("mu_pfSumIso04_sumChargedHadronPt", "std::vector<float>", &mu_pfSumIso04_sumChargedHadronPt, buffersize);
   if( doWrite("mu_pfSumIso04_sumChargedParticlePt") ) tree->Branch("mu_pfSumIso04_sumChargedParticlePt", "std::vector<float>", &mu_pfSumIso04_sumChargedParticlePt, buffersize);
   if( doWrite("mu_pfSumIso04_sumNeutralHadronEt") ) tree->Branch("mu_pfSumIso04_sumNeutralHadronEt", "std::vector<float>", &mu_pfSumIso04_sumNeutralHadronEt, buffersize);
   if( doWrite("mu_pfSumIso04_sumNeutralHadronEtHighThreshold") ) tree->Branch("mu_pfSumIso04_sumNeutralHadronEtHighThreshold", "std::vector<float>", &mu_pfSumIso04_sumNeutralHadronEtHighThreshold, buffersize);
   if( doWrite("mu_pfSumIso04_sumPhotonEt") ) tree->Branch("mu_pfSumIso04_sumPhotonEt", "std::vector<float>", &mu_pfSumIso04_sumPhotonEt, buffersize);
   if( doWrite("mu_pfSumIso04_sumPhotonEtHighThreshold") ) tree->Branch("mu_pfSumIso04_sumPhotonEtHighThreshold", "std::vector<float>", &mu_pfSumIso04_sumPhotonEtHighThreshold, buffersize);
   if( doWrite("mu_pfSumIso04_sumPUPt") ) tree->Branch("mu_pfSumIso04_sumPUPt", "std::vector<float>", &mu_pfSumIso04_sumPUPt, buffersize);
   
   if( doWrite("mu_miniIso") ) tree->Branch("mu_miniIso", "std::vector<float>", &mu_miniIso, buffersize);
   if( doWrite("mu_miniIsoTTH") ) tree->Branch("mu_miniIsoTTH", "std::vector<float>", &mu_miniIsoTTH, buffersize);

   if( doWrite("mu_isGlobalMuon") ) tree->Branch("mu_isGlobalMuon", "std::vector<int>", &mu_isGlobalMuon, buffersize);
   if( doWrite("mu_isTrackerMuon") ) tree->Branch("mu_isTrackerMuon", "std::vector<int>", &mu_isTrackerMuon, buffersize);
   if( doWrite("mu_isStandAloneMuon") ) tree->Branch("mu_isStandAloneMuon", "std::vector<int>", &mu_isStandAloneMuon, buffersize);
   if( doWrite("mu_isCaloMuon") ) tree->Branch("mu_isCaloMuon", "std::vector<int>", &mu_isCaloMuon, buffersize);
   if( doWrite("mu_isPFMuon") ) tree->Branch("mu_isPFMuon", "std::vector<int>", &mu_isPFMuon, buffersize);
   if( doWrite("mu_isRPCMuon") ) tree->Branch("mu_isRPCMuon", "std::vector<int>", &mu_isRPCMuon, buffersize);

   if( doWrite("mu_vx") ) tree->Branch("mu_vx", "std::vector<float>", &mu_vx, buffersize);
   if( doWrite("mu_vy") ) tree->Branch("mu_vy", "std::vector<float>", &mu_vy, buffersize);
   if( doWrite("mu_vz") ) tree->Branch("mu_vz", "std::vector<float>", &mu_vz, buffersize);
   
   if( doWrite("mu_numberOfMatches") ) tree->Branch("mu_numberOfMatches", "std::vector<int>", &mu_numberOfMatches, buffersize);
   if( doWrite("mu_numberOfMatchedStations") ) tree->Branch("mu_numberOfMatchedStations", "std::vector<int>", &mu_numberOfMatchedStations, buffersize);
   
   if( doWrite("mu_segmentCompatibility") ) tree->Branch("mu_segmentCompatibility", "std::vector<float>", &mu_segmentCompatibility, buffersize);
   if( doWrite("mu_caloCompatibility") ) tree->Branch("mu_caloCompatibility", "std::vector<float>", &mu_caloCompatibility, buffersize);
   
   if( doWrite("mu_combinedQuality_chi2LocalPosition") ) tree->Branch("mu_combinedQuality_chi2LocalPosition", "std::vector<float>", &mu_combinedQuality_chi2LocalPosition, buffersize);
   if( doWrite("mu_combinedQuality_trkKink") ) tree->Branch("mu_combinedQuality_trkKink", "std::vector<float>", &mu_combinedQuality_trkKink, buffersize);

   if( doWrite("mu_isLooseMuon") ) tree->Branch("mu_isLooseMuon", "std::vector<bool>", &mu_isLooseMuon, buffersize);
   if( doWrite("mu_isMediumMuon") ) tree->Branch("mu_isMediumMuon", "std::vector<bool>", &mu_isMediumMuon, buffersize);
   if( doWrite("mu_isTightMuon") ) tree->Branch("mu_isTightMuon", "std::vector<bool>", &mu_isTightMuon, buffersize);
   if( doWrite("mu_isSoftMuon") ) tree->Branch("mu_isSoftMuon", "std::vector<bool>", &mu_isSoftMuon, buffersize);
   if( doWrite("mu_isHighPtMuon") ) tree->Branch("mu_isHighPtMuon", "std::vector<bool>", &mu_isHighPtMuon, buffersize);
   
   if( doWrite("mu_isGoodMuon_AllGlobalMuons") ) tree->Branch("mu_isGoodMuon_AllGlobalMuons", "std::vector<bool>", &mu_isGoodMuon_AllGlobalMuons, buffersize);
   if( doWrite("mu_isGoodMuon_AllStandAloneMuons") ) tree->Branch("mu_isGoodMuon_AllStandAloneMuons", "std::vector<bool>", &mu_isGoodMuon_AllStandAloneMuons, buffersize);
   if( doWrite("mu_isGoodMuon_AllTrackerMuons") ) tree->Branch("mu_isGoodMuon_AllTrackerMuons", "std::vector<bool>", &mu_isGoodMuon_AllTrackerMuons, buffersize);
   if( doWrite("mu_isGoodMuon_TrackerMuonArbitrated") ) tree->Branch("mu_isGoodMuon_TrackerMuonArbitrated", "std::vector<bool>", &mu_isGoodMuon_TrackerMuonArbitrated, buffersize);
   if( doWrite("mu_isGoodMuon_AllArbitrated") ) tree->Branch("mu_isGoodMuon_AllArbitrated", "std::vector<bool>", &mu_isGoodMuon_AllArbitrated, buffersize);
   if( doWrite("mu_isGoodMuon_GlobalMuonPromptTight") ) tree->Branch("mu_isGoodMuon_GlobalMuonPromptTight", "std::vector<bool>", &mu_isGoodMuon_GlobalMuonPromptTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationLoose") ) tree->Branch("mu_isGoodMuon_TMLastStationLoose", "std::vector<bool>", &mu_isGoodMuon_TMLastStationLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationTight") ) tree->Branch("mu_isGoodMuon_TMLastStationTight", "std::vector<bool>", &mu_isGoodMuon_TMLastStationTight, buffersize);
   if( doWrite("mu_isGoodMuon_TM2DCompatibilityLoose") ) tree->Branch("mu_isGoodMuon_TM2DCompatibilityLoose", "std::vector<bool>", &mu_isGoodMuon_TM2DCompatibilityLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TM2DCompatibilityTight") ) tree->Branch("mu_isGoodMuon_TM2DCompatibilityTight", "std::vector<bool>", &mu_isGoodMuon_TM2DCompatibilityTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMOneStationLoose") ) tree->Branch("mu_isGoodMuon_TMOneStationLoose", "std::vector<bool>", &mu_isGoodMuon_TMOneStationLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMOneStationTight") ) tree->Branch("mu_isGoodMuon_TMOneStationTight", "std::vector<bool>", &mu_isGoodMuon_TMOneStationTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationOptimizedLowPtLoose") ) tree->Branch("mu_isGoodMuon_TMLastStationOptimizedLowPtLoose", "std::vector<bool>", &mu_isGoodMuon_TMLastStationOptimizedLowPtLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationOptimizedLowPtTight") ) tree->Branch("mu_isGoodMuon_TMLastStationOptimizedLowPtTight", "std::vector<bool>", &mu_isGoodMuon_TMLastStationOptimizedLowPtTight, buffersize);
   if( doWrite("mu_isGoodMuon_GMTkChiCompatibility") ) tree->Branch("mu_isGoodMuon_GMTkChiCompatibility", "std::vector<bool>", &mu_isGoodMuon_GMTkChiCompatibility, buffersize);
   if( doWrite("mu_isGoodMuon_GMStaChiCompatibility") ) tree->Branch("mu_isGoodMuon_GMStaChiCompatibility", "std::vector<bool>", &mu_isGoodMuon_GMStaChiCompatibility, buffersize);
   if( doWrite("mu_isGoodMuon_GMTkKinkTight") ) tree->Branch("mu_isGoodMuon_GMTkKinkTight", "std::vector<bool>", &mu_isGoodMuon_GMTkKinkTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationAngLoose") ) tree->Branch("mu_isGoodMuon_TMLastStationAngLoose", "std::vector<bool>", &mu_isGoodMuon_TMLastStationAngLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationAngTight") ) tree->Branch("mu_isGoodMuon_TMLastStationAngTight", "std::vector<bool>", &mu_isGoodMuon_TMLastStationAngTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMOneStationAngLoose") ) tree->Branch("mu_isGoodMuon_TMOneStationAngLoose", "std::vector<bool>", &mu_isGoodMuon_TMOneStationAngLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMOneStationAngTight") ) tree->Branch("mu_isGoodMuon_TMOneStationAngTight", "std::vector<bool>", &mu_isGoodMuon_TMOneStationAngTight, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtLoose") ) tree->Branch("mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtLoose", "std::vector<bool>", &mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtLoose, buffersize);
   if( doWrite("mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight") ) tree->Branch("mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight", "std::vector<bool>", &mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight, buffersize);

   if( doWrite("mu_calEnergy_em") ) tree->Branch("mu_calEnergy_em", "std::vector<float>", &mu_calEnergy_em, buffersize);
   if( doWrite("mu_calEnergy_had") ) tree->Branch("mu_calEnergy_had", "std::vector<float>", &mu_calEnergy_had, buffersize);
   if( doWrite("mu_calEnergy_ho") ) tree->Branch("mu_calEnergy_ho", "std::vector<float>", &mu_calEnergy_ho, buffersize);
   if( doWrite("mu_calEnergy_emS9") ) tree->Branch("mu_calEnergy_emS9", "std::vector<float>", &mu_calEnergy_emS9, buffersize);
   if( doWrite("mu_calEnergy_hadS9") ) tree->Branch("mu_calEnergy_hadS9", "std::vector<float>", &mu_calEnergy_hadS9, buffersize);
   if( doWrite("mu_calEnergy_hoS9") ) tree->Branch("mu_calEnergy_hoS9", "std::vector<float>", &mu_calEnergy_hoS9, buffersize);
   if( doWrite("mu_calEnergy_emS25") ) tree->Branch("mu_calEnergy_emS25", "std::vector<float>", &mu_calEnergy_emS25, buffersize);
   if( doWrite("mu_calEnergy_emMax") ) tree->Branch("mu_calEnergy_emMax", "std::vector<float>", &mu_calEnergy_emMax, buffersize);
   if( doWrite("mu_calEnergy_hadMax") ) tree->Branch("mu_calEnergy_hadMax", "std::vector<float>", &mu_calEnergy_hadMax, buffersize);
   if( doWrite("mu_calEnergy_ecal_time") ) tree->Branch("mu_calEnergy_ecal_time", "std::vector<float>", &mu_calEnergy_ecal_time, buffersize);
   if( doWrite("mu_calEnergy_hcal_time") ) tree->Branch("mu_calEnergy_hcal_time", "std::vector<float>", &mu_calEnergy_hcal_time, buffersize);
   if( doWrite("mu_calEnergy_ecal_rawId") ) tree->Branch("mu_calEnergy_ecal_rawId", "std::vector<int>", &mu_calEnergy_ecal_rawId, buffersize);
   if( doWrite("mu_calEnergy_hcal_rawId") ) tree->Branch("mu_calEnergy_hcal_rawId", "std::vector<int>", &mu_calEnergy_hcal_rawId, buffersize);

   if( doWrite("mu_isolationR03_trackerVetoPt") ) tree->Branch("mu_isolationR03_trackerVetoPt", "std::vector<float>", &mu_isolationR03_trackerVetoPt, buffersize);
   if( doWrite("mu_isolationR03_emVetoEt") ) tree->Branch("mu_isolationR03_emVetoEt", "std::vector<float>", &mu_isolationR03_emVetoEt, buffersize);
   if( doWrite("mu_isolationR03_hadVetoEt") ) tree->Branch("mu_isolationR03_hadVetoEt", "std::vector<float>", &mu_isolationR03_hadVetoEt, buffersize);
   if( doWrite("mu_isolationR03_hoVetoEt") ) tree->Branch("mu_isolationR03_hoVetoEt", "std::vector<float>", &mu_isolationR03_hoVetoEt, buffersize);
   if( doWrite("mu_isolationR03_sumPt") ) tree->Branch("mu_isolationR03_sumPt", "std::vector<float>", &mu_isolationR03_sumPt, buffersize);
   if( doWrite("mu_isolationR03_emEt") ) tree->Branch("mu_isolationR03_emEt", "std::vector<float>", &mu_isolationR03_emEt, buffersize);
   if( doWrite("mu_isolationR03_hadEt") ) tree->Branch("mu_isolationR03_hadEt", "std::vector<float>", &mu_isolationR03_hadEt, buffersize);
   if( doWrite("mu_isolationR03_hoEt") ) tree->Branch("mu_isolationR03_hoEt", "std::vector<float>", &mu_isolationR03_hoEt, buffersize);
   if( doWrite("mu_isolationR03_nTracks") ) tree->Branch("mu_isolationR03_nTracks", "std::vector<int>", &mu_isolationR03_nTracks, buffersize);
   if( doWrite("mu_isolationR03_nJets") ) tree->Branch("mu_isolationR03_nJets", "std::vector<int>", &mu_isolationR03_nJets, buffersize);

   if( doWrite("mu_isolationR05_trackerVetoPt") ) tree->Branch("mu_isolationR05_trackerVetoPt", "std::vector<float>", &mu_isolationR05_trackerVetoPt, buffersize);
   if( doWrite("mu_isolationR05_emVetoEt") ) tree->Branch("mu_isolationR05_emVetoEt", "std::vector<float>", &mu_isolationR05_emVetoEt, buffersize);
   if( doWrite("mu_isolationR05_hadVetoEt") ) tree->Branch("mu_isolationR05_hadVetoEt", "std::vector<float>", &mu_isolationR05_hadVetoEt, buffersize);
   if( doWrite("mu_isolationR05_hoVetoEt") ) tree->Branch("mu_isolationR05_hoVetoEt", "std::vector<float>", &mu_isolationR05_hoVetoEt, buffersize);
   if( doWrite("mu_isolationR05_sumPt") ) tree->Branch("mu_isolationR05_sumPt", "std::vector<float>", &mu_isolationR05_sumPt, buffersize);
   if( doWrite("mu_isolationR05_emEt") ) tree->Branch("mu_isolationR05_emEt", "std::vector<float>", &mu_isolationR05_emEt, buffersize);
   if( doWrite("mu_isolationR05_hadEt") ) tree->Branch("mu_isolationR05_hadEt", "std::vector<float>", &mu_isolationR05_hadEt, buffersize);
   if( doWrite("mu_isolationR05_hoEt") ) tree->Branch("mu_isolationR05_hoEt", "std::vector<float>", &mu_isolationR05_hoEt, buffersize);
   if( doWrite("mu_isolationR05_nTracks") ) tree->Branch("mu_isolationR05_nTracks", "std::vector<int>", &mu_isolationR05_nTracks, buffersize);
   if( doWrite("mu_isolationR05_nJets") ) tree->Branch("mu_isolationR05_nJets", "std::vector<int>", &mu_isolationR05_nJets, buffersize);

   if( doWrite("mu_hasGlobalTrack") ) tree->Branch("mu_hasGlobalTrack", "std::vector<int>", &mu_hasGlobalTrack, buffersize);
   if( doWrite("mu_globalTrack_d0") ) tree->Branch("mu_globalTrack_d0", "std::vector<float>", &mu_globalTrack_d0, buffersize);
   if( doWrite("mu_globalTrack_z0") ) tree->Branch("mu_globalTrack_z0", "std::vector<float>", &mu_globalTrack_z0, buffersize);
   if( doWrite("mu_globalTrack_d0Error") ) tree->Branch("mu_globalTrack_d0Error", "std::vector<float>", &mu_globalTrack_d0Error, buffersize);
   if( doWrite("mu_globalTrack_z0Error") ) tree->Branch("mu_globalTrack_z0Error", "std::vector<float>", &mu_globalTrack_z0Error, buffersize);
   if( doWrite("mu_globalTrack_PV_dxy") ) tree->Branch("mu_globalTrack_PV_dxy", "std::vector<float>", &mu_globalTrack_PV_dxy, buffersize);
   if( doWrite("mu_globalTrack_PV_dz") ) tree->Branch("mu_globalTrack_PV_dz", "std::vector<float>", &mu_globalTrack_PV_dz, buffersize);
   if( doWrite("mu_globalTrack_RP_dxy") ) tree->Branch("mu_globalTrack_RP_dxy", "std::vector<float>", &mu_globalTrack_RP_dxy, buffersize);
   if( doWrite("mu_globalTrack_RP_dz") ) tree->Branch("mu_globalTrack_RP_dz", "std::vector<float>", &mu_globalTrack_RP_dz, buffersize);
   if( doWrite("mu_globalTrack_BS_dxy") ) tree->Branch("mu_globalTrack_BS_dxy", "std::vector<float>", &mu_globalTrack_BS_dxy, buffersize);
   if( doWrite("mu_globalTrack_BS_dz") ) tree->Branch("mu_globalTrack_BS_dz", "std::vector<float>", &mu_globalTrack_BS_dz, buffersize);
   if( doWrite("mu_globalTrack_dxyError") ) tree->Branch("mu_globalTrack_dxyError", "std::vector<float>", &mu_globalTrack_dxyError, buffersize);
   if( doWrite("mu_globalTrack_dzError") ) tree->Branch("mu_globalTrack_dzError", "std::vector<float>", &mu_globalTrack_dzError, buffersize);
   if( doWrite("mu_globalTrack_normalizedChi2") ) tree->Branch("mu_globalTrack_normalizedChi2", "std::vector<float>", &mu_globalTrack_normalizedChi2, buffersize);
   if( doWrite("mu_globalTrack_numberOfValidHits") ) tree->Branch("mu_globalTrack_numberOfValidHits", "std::vector<int>", &mu_globalTrack_numberOfValidHits, buffersize);
   if( doWrite("mu_globalTrack_numberOfValidMuonHits") ) tree->Branch("mu_globalTrack_numberOfValidMuonHits", "std::vector<int>", &mu_globalTrack_numberOfValidMuonHits, buffersize);
   if( doWrite("mu_globalTrack_numberOfLostHits") ) tree->Branch("mu_globalTrack_numberOfLostHits", "std::vector<int>", &mu_globalTrack_numberOfLostHits, buffersize);
   if( doWrite("mu_globalTrack_pt") ) tree->Branch("mu_globalTrack_pt", "std::vector<float>", &mu_globalTrack_pt, buffersize);
   if( doWrite("mu_globalTrack_eta") ) tree->Branch("mu_globalTrack_eta", "std::vector<float>", &mu_globalTrack_eta, buffersize);
   if( doWrite("mu_globalTrack_phi") ) tree->Branch("mu_globalTrack_phi", "std::vector<float>", &mu_globalTrack_phi, buffersize);
   if( doWrite("mu_globalTrack_ptError") ) tree->Branch("mu_globalTrack_ptError", "std::vector<float>", &mu_globalTrack_ptError, buffersize);
   if( doWrite("mu_globalTrack_etaError") ) tree->Branch("mu_globalTrack_etaError", "std::vector<float>", &mu_globalTrack_etaError, buffersize);
   if( doWrite("mu_globalTrack_phiError") ) tree->Branch("mu_globalTrack_phiError", "std::vector<float>", &mu_globalTrack_phiError, buffersize);
   if( doWrite("mu_globalTrack_vx") ) tree->Branch("mu_globalTrack_vx", "std::vector<float>", &mu_globalTrack_vx, buffersize);
   if( doWrite("mu_globalTrack_vy") ) tree->Branch("mu_globalTrack_vy", "std::vector<float>", &mu_globalTrack_vy, buffersize);
   if( doWrite("mu_globalTrack_vz") ) tree->Branch("mu_globalTrack_vz", "std::vector<float>", &mu_globalTrack_vz, buffersize);
   if( doWrite("mu_globalTrack_qoverp") ) tree->Branch("mu_globalTrack_qoverp", "std::vector<float>", &mu_globalTrack_qoverp, buffersize);
   if( doWrite("mu_globalTrack_qoverpError") ) tree->Branch("mu_globalTrack_qoverpError", "std::vector<float>", &mu_globalTrack_qoverpError, buffersize);
   if( doWrite("mu_globalTrack_charge") ) tree->Branch("mu_globalTrack_charge", "std::vector<int>", &mu_globalTrack_charge, buffersize);
   if( doWrite("mu_globalTrack_trackerLayersWithMeasurement") ) tree->Branch("mu_globalTrack_trackerLayersWithMeasurement", "std::vector<int>", &mu_globalTrack_trackerLayersWithMeasurement, buffersize);
   if( doWrite("mu_globalTrack_pixelLayersWithMeasurement") ) tree->Branch("mu_globalTrack_pixelLayersWithMeasurement", "std::vector<int>", &mu_globalTrack_pixelLayersWithMeasurement, buffersize);
   if( doWrite("mu_globalTrack_numberOfValidStripLayersWithMonoAndStereo") ) tree->Branch("mu_globalTrack_numberOfValidStripLayersWithMonoAndStereo", "std::vector<int>", &mu_globalTrack_numberOfValidStripLayersWithMonoAndStereo, buffersize);
   if( doWrite("mu_globalTrack_trackerLayersWithoutMeasurement") ) tree->Branch("mu_globalTrack_trackerLayersWithoutMeasurement", "std::vector<int>", &mu_globalTrack_trackerLayersWithoutMeasurement, buffersize);
   if( doWrite("mu_globalTrack_numberOfValidPixelHits") ) tree->Branch("mu_globalTrack_numberOfValidPixelHits", "std::vector<int>", &mu_globalTrack_numberOfValidPixelHits, buffersize);
   if( doWrite("mu_globalTrack_numberOfLostPixelHits") ) tree->Branch("mu_globalTrack_numberOfLostPixelHits", "std::vector<int>", &mu_globalTrack_numberOfLostPixelHits, buffersize);
   if( doWrite("mu_globalTrack_numberOfInnerHits") ) tree->Branch("mu_globalTrack_numberOfInnerHits", "std::vector<int>", &mu_globalTrack_numberOfInnerHits, buffersize);
   if( doWrite("mu_globalTrack_numberOfOuterHits") ) tree->Branch("mu_globalTrack_numberOfOuterHits", "std::vector<int>", &mu_globalTrack_numberOfOuterHits, buffersize);
   if( doWrite("mu_globalTrack_validFraction") ) tree->Branch("mu_globalTrack_validFraction", "std::vector<float>", &mu_globalTrack_validFraction, buffersize);
   
   if( doWrite("mu_bestTrackType") ) tree->Branch("mu_bestTrackType", "std::vector<int>", &mu_bestTrackType, buffersize);
   if( doWrite("mu_hasBestTrack") ) tree->Branch("mu_hasBestTrack", "std::vector<int>", &mu_hasBestTrack, buffersize);
   if( doWrite("mu_bestTrack_d0") ) tree->Branch("mu_bestTrack_d0", "std::vector<float>", &mu_bestTrack_d0, buffersize);
   if( doWrite("mu_bestTrack_z0") ) tree->Branch("mu_bestTrack_z0", "std::vector<float>", &mu_bestTrack_z0, buffersize);
   if( doWrite("mu_bestTrack_d0Error") ) tree->Branch("mu_bestTrack_d0Error", "std::vector<float>", &mu_bestTrack_d0Error, buffersize);
   if( doWrite("mu_bestTrack_z0Error") ) tree->Branch("mu_bestTrack_z0Error", "std::vector<float>", &mu_bestTrack_z0Error, buffersize);
   if( doWrite("mu_bestTrack_PV_dxy") ) tree->Branch("mu_bestTrack_PV_dxy", "std::vector<float>", &mu_bestTrack_PV_dxy, buffersize);
   if( doWrite("mu_bestTrack_PV_dz") ) tree->Branch("mu_bestTrack_PV_dz", "std::vector<float>", &mu_bestTrack_PV_dz, buffersize);
   if( doWrite("mu_bestTrack_RP_dxy") ) tree->Branch("mu_bestTrack_RP_dxy", "std::vector<float>", &mu_bestTrack_RP_dxy, buffersize);
   if( doWrite("mu_bestTrack_RP_dz") ) tree->Branch("mu_bestTrack_RP_dz", "std::vector<float>", &mu_bestTrack_RP_dz, buffersize);
   if( doWrite("mu_bestTrack_BS_dxy") ) tree->Branch("mu_bestTrack_BS_dxy", "std::vector<float>", &mu_bestTrack_BS_dxy, buffersize);
   if( doWrite("mu_bestTrack_BS_dz") ) tree->Branch("mu_bestTrack_BS_dz", "std::vector<float>", &mu_bestTrack_BS_dz, buffersize);
   if( doWrite("mu_bestTrack_dxyError") ) tree->Branch("mu_bestTrack_dxyError", "std::vector<float>", &mu_bestTrack_dxyError, buffersize);
   if( doWrite("mu_bestTrack_dzError") ) tree->Branch("mu_bestTrack_dzError", "std::vector<float>", &mu_bestTrack_dzError, buffersize);
   if( doWrite("mu_bestTrack_normalizedChi2") ) tree->Branch("mu_bestTrack_normalizedChi2", "std::vector<float>", &mu_bestTrack_normalizedChi2, buffersize);
   if( doWrite("mu_bestTrack_numberOfValidHits") ) tree->Branch("mu_bestTrack_numberOfValidHits", "std::vector<int>", &mu_bestTrack_numberOfValidHits, buffersize);
   if( doWrite("mu_bestTrack_numberOfLostHits") ) tree->Branch("mu_bestTrack_numberOfLostHits", "std::vector<int>", &mu_bestTrack_numberOfLostHits, buffersize);
   if( doWrite("mu_bestTrack_pt") ) tree->Branch("mu_bestTrack_pt", "std::vector<float>", &mu_bestTrack_pt, buffersize);
   if( doWrite("mu_bestTrack_eta") ) tree->Branch("mu_bestTrack_eta", "std::vector<float>", &mu_bestTrack_eta, buffersize);
   if( doWrite("mu_bestTrack_phi") ) tree->Branch("mu_bestTrack_phi", "std::vector<float>", &mu_bestTrack_phi, buffersize);
   if( doWrite("mu_bestTrack_ptError") ) tree->Branch("mu_bestTrack_ptError", "std::vector<float>", &mu_bestTrack_ptError, buffersize);
   if( doWrite("mu_bestTrack_etaError") ) tree->Branch("mu_bestTrack_etaError", "std::vector<float>", &mu_bestTrack_etaError, buffersize);
   if( doWrite("mu_bestTrack_phiError") ) tree->Branch("mu_bestTrack_phiError", "std::vector<float>", &mu_bestTrack_phiError, buffersize);
   if( doWrite("mu_bestTrack_vx") ) tree->Branch("mu_bestTrack_vx", "std::vector<float>", &mu_bestTrack_vx, buffersize);
   if( doWrite("mu_bestTrack_vy") ) tree->Branch("mu_bestTrack_vy", "std::vector<float>", &mu_bestTrack_vy, buffersize);
   if( doWrite("mu_bestTrack_vz") ) tree->Branch("mu_bestTrack_vz", "std::vector<float>", &mu_bestTrack_vz, buffersize);
   if( doWrite("mu_bestTrack_qoverp") ) tree->Branch("mu_bestTrack_qoverp", "std::vector<float>", &mu_bestTrack_qoverp, buffersize);
   if( doWrite("mu_bestTrack_qoverpError") ) tree->Branch("mu_bestTrack_qoverpError", "std::vector<float>", &mu_bestTrack_qoverpError, buffersize);
   if( doWrite("mu_bestTrack_charge") ) tree->Branch("mu_bestTrack_charge", "std::vector<int>", &mu_bestTrack_charge, buffersize);
   if( doWrite("mu_bestTrack_trackerLayersWithMeasurement") ) tree->Branch("mu_bestTrack_trackerLayersWithMeasurement", "std::vector<int>", &mu_bestTrack_trackerLayersWithMeasurement, buffersize);
   if( doWrite("mu_bestTrack_pixelLayersWithMeasurement") ) tree->Branch("mu_bestTrack_pixelLayersWithMeasurement", "std::vector<int>", &mu_bestTrack_pixelLayersWithMeasurement, buffersize);
   if( doWrite("mu_bestTrack_numberOfValidStripLayersWithMonoAndStereo") ) tree->Branch("mu_bestTrack_numberOfValidStripLayersWithMonoAndStereo", "std::vector<int>", &mu_bestTrack_numberOfValidStripLayersWithMonoAndStereo, buffersize);
   if( doWrite("mu_bestTrack_trackerLayersWithoutMeasurement") ) tree->Branch("mu_bestTrack_trackerLayersWithoutMeasurement", "std::vector<int>", &mu_bestTrack_trackerLayersWithoutMeasurement, buffersize);
   if( doWrite("mu_bestTrack_numberOfValidPixelHits") ) tree->Branch("mu_bestTrack_numberOfValidPixelHits", "std::vector<int>", &mu_bestTrack_numberOfValidPixelHits, buffersize);
   if( doWrite("mu_bestTrack_numberOfLostPixelHits") ) tree->Branch("mu_bestTrack_numberOfLostPixelHits", "std::vector<int>", &mu_bestTrack_numberOfLostPixelHits, buffersize);
   if( doWrite("mu_bestTrack_numberOfInnerHits") ) tree->Branch("mu_bestTrack_numberOfInnerHits", "std::vector<int>", &mu_bestTrack_numberOfInnerHits, buffersize);
   if( doWrite("mu_bestTrack_numberOfOuterHits") ) tree->Branch("mu_bestTrack_numberOfOuterHits", "std::vector<int>", &mu_bestTrack_numberOfOuterHits, buffersize);
   if( doWrite("mu_bestTrack_validFraction") ) tree->Branch("mu_bestTrack_validFraction", "std::vector<float>", &mu_bestTrack_validFraction, buffersize);
   
   if( doWrite("mu_hasInnerTrack") ) tree->Branch("mu_hasInnerTrack", "std::vector<int>", &mu_hasInnerTrack, buffersize);
   if( doWrite("mu_innerTrack_d0") ) tree->Branch("mu_innerTrack_d0", "std::vector<float>", &mu_innerTrack_d0, buffersize);
   if( doWrite("mu_innerTrack_z0") ) tree->Branch("mu_innerTrack_z0", "std::vector<float>", &mu_innerTrack_z0, buffersize);
   if( doWrite("mu_innerTrack_d0Error") ) tree->Branch("mu_innerTrack_d0Error", "std::vector<float>", &mu_innerTrack_d0Error, buffersize);
   if( doWrite("mu_innerTrack_z0Error") ) tree->Branch("mu_innerTrack_z0Error", "std::vector<float>", &mu_innerTrack_z0Error, buffersize);
   if( doWrite("mu_innerTrack_PV_dxy") ) tree->Branch("mu_innerTrack_PV_dxy", "std::vector<float>", &mu_innerTrack_PV_dxy, buffersize);
   if( doWrite("mu_innerTrack_PV_dz") ) tree->Branch("mu_innerTrack_PV_dz", "std::vector<float>", &mu_innerTrack_PV_dz, buffersize);
   if( doWrite("mu_innerTrack_RP_dxy") ) tree->Branch("mu_innerTrack_RP_dxy", "std::vector<float>", &mu_innerTrack_RP_dxy, buffersize);
   if( doWrite("mu_innerTrack_RP_dz") ) tree->Branch("mu_innerTrack_RP_dz", "std::vector<float>", &mu_innerTrack_RP_dz, buffersize);
   if( doWrite("mu_innerTrack_BS_dxy") ) tree->Branch("mu_innerTrack_BS_dxy", "std::vector<float>", &mu_innerTrack_BS_dxy, buffersize);
   if( doWrite("mu_innerTrack_BS_dz") ) tree->Branch("mu_innerTrack_BS_dz", "std::vector<float>", &mu_innerTrack_BS_dz, buffersize);
   if( doWrite("mu_innerTrack_dxyError") ) tree->Branch("mu_innerTrack_dxyError", "std::vector<float>", &mu_innerTrack_dxyError, buffersize);
   if( doWrite("mu_innerTrack_dzError") ) tree->Branch("mu_innerTrack_dzError", "std::vector<float>", &mu_innerTrack_dzError, buffersize);
   if( doWrite("mu_innerTrack_normalizedChi2") ) tree->Branch("mu_innerTrack_normalizedChi2", "std::vector<float>", &mu_innerTrack_normalizedChi2, buffersize);
   if( doWrite("mu_innerTrack_numberOfValidHits") ) tree->Branch("mu_innerTrack_numberOfValidHits", "std::vector<int>", &mu_innerTrack_numberOfValidHits, buffersize);
   if( doWrite("mu_innerTrack_numberOfLostHits") ) tree->Branch("mu_innerTrack_numberOfLostHits", "std::vector<int>", &mu_innerTrack_numberOfLostHits, buffersize);
   if( doWrite("mu_innerTrack_pt") ) tree->Branch("mu_innerTrack_pt", "std::vector<float>", &mu_innerTrack_pt, buffersize);
   if( doWrite("mu_innerTrack_eta") ) tree->Branch("mu_innerTrack_eta", "std::vector<float>", &mu_innerTrack_eta, buffersize);
   if( doWrite("mu_innerTrack_phi") ) tree->Branch("mu_innerTrack_phi", "std::vector<float>", &mu_innerTrack_phi, buffersize);
   if( doWrite("mu_innerTrack_ptError") ) tree->Branch("mu_innerTrack_ptError", "std::vector<float>", &mu_innerTrack_ptError, buffersize);
   if( doWrite("mu_innerTrack_etaError") ) tree->Branch("mu_innerTrack_etaError", "std::vector<float>", &mu_innerTrack_etaError, buffersize);
   if( doWrite("mu_innerTrack_phiError") ) tree->Branch("mu_innerTrack_phiError", "std::vector<float>", &mu_innerTrack_phiError, buffersize);
   if( doWrite("mu_innerTrack_vx") ) tree->Branch("mu_innerTrack_vx", "std::vector<float>", &mu_innerTrack_vx, buffersize);
   if( doWrite("mu_innerTrack_vy") ) tree->Branch("mu_innerTrack_vy", "std::vector<float>", &mu_innerTrack_vy, buffersize);
   if( doWrite("mu_innerTrack_vz") ) tree->Branch("mu_innerTrack_vz", "std::vector<float>", &mu_innerTrack_vz, buffersize);
   if( doWrite("mu_innerTrack_qoverp") ) tree->Branch("mu_innerTrack_qoverp", "std::vector<float>", &mu_innerTrack_qoverp, buffersize);
   if( doWrite("mu_innerTrack_qoverpError") ) tree->Branch("mu_innerTrack_qoverpError", "std::vector<float>", &mu_innerTrack_qoverpError, buffersize);
   if( doWrite("mu_innerTrack_charge") ) tree->Branch("mu_innerTrack_charge", "std::vector<int>", &mu_innerTrack_charge, buffersize);
   if( doWrite("mu_innerTrack_trackerLayersWithMeasurement") ) tree->Branch("mu_innerTrack_trackerLayersWithMeasurement", "std::vector<int>", &mu_innerTrack_trackerLayersWithMeasurement, buffersize);
   if( doWrite("mu_innerTrack_pixelLayersWithMeasurement") ) tree->Branch("mu_innerTrack_pixelLayersWithMeasurement", "std::vector<int>", &mu_innerTrack_pixelLayersWithMeasurement, buffersize);
   if( doWrite("mu_innerTrack_numberOfValidStripLayersWithMonoAndStereo") ) tree->Branch("mu_innerTrack_numberOfValidStripLayersWithMonoAndStereo", "std::vector<int>", &mu_innerTrack_numberOfValidStripLayersWithMonoAndStereo, buffersize);
   if( doWrite("mu_innerTrack_trackerLayersWithoutMeasurement") ) tree->Branch("mu_innerTrack_trackerLayersWithoutMeasurement", "std::vector<int>", &mu_innerTrack_trackerLayersWithoutMeasurement, buffersize);
   if( doWrite("mu_innerTrack_numberOfValidPixelHits") ) tree->Branch("mu_innerTrack_numberOfValidPixelHits", "std::vector<int>", &mu_innerTrack_numberOfValidPixelHits, buffersize);
   if( doWrite("mu_innerTrack_numberOfLostPixelHits") ) tree->Branch("mu_innerTrack_numberOfLostPixelHits", "std::vector<int>", &mu_innerTrack_numberOfLostPixelHits, buffersize);
   if( doWrite("mu_innerTrack_numberOfInnerHits") ) tree->Branch("mu_innerTrack_numberOfInnerHits", "std::vector<int>", &mu_innerTrack_numberOfInnerHits, buffersize);
   if( doWrite("mu_innerTrack_numberOfOuterHits") ) tree->Branch("mu_innerTrack_numberOfOuterHits", "std::vector<int>", &mu_innerTrack_numberOfOuterHits, buffersize);
   if( doWrite("mu_innerTrack_validFraction") ) tree->Branch("mu_innerTrack_validFraction", "std::vector<float>", &mu_innerTrack_validFraction, buffersize);
   
   if( doWrite("mu_type") ) tree->Branch("mu_type", "std::vector<int>", &mu_type, buffersize);
   
   if( doWrite("mu_lepMVA") ) tree->Branch("mu_lepMVA", "std::vector<float>", &mu_lepMVA, buffersize);

   if( doWrite("mu_lepMVA_pt") ) tree->Branch("mu_lepMVA_pt", "std::vector<float>", &mu_lepMVA_pt, buffersize);
   if( doWrite("mu_lepMVA_eta") ) tree->Branch("mu_lepMVA_eta", "std::vector<float>", &mu_lepMVA_eta, buffersize);
   if( doWrite("mu_lepMVA_miniRelIsoCharged") ) tree->Branch("mu_lepMVA_miniRelIsoCharged", "std::vector<float>", &mu_lepMVA_miniRelIsoCharged, buffersize);
   if( doWrite("mu_lepMVA_miniRelIsoNeutral") ) tree->Branch("mu_lepMVA_miniRelIsoNeutral", "std::vector<float>", &mu_lepMVA_miniRelIsoNeutral, buffersize);
   if( doWrite("mu_lepMVA_jetPtRatio") ) tree->Branch("mu_lepMVA_jetPtRatio", "std::vector<float>", &mu_lepMVA_jetPtRatio, buffersize);
   if( doWrite("mu_lepMVA_jetPtRelv2") ) tree->Branch("mu_lepMVA_jetPtRelv2", "std::vector<float>", &mu_lepMVA_jetPtRelv2, buffersize);
   if( doWrite("mu_lepMVA_jetBTagCSV") ) tree->Branch("mu_lepMVA_jetBTagCSV", "std::vector<float>", &mu_lepMVA_jetBTagCSV, buffersize);
   if( doWrite("mu_lepMVA_sip3d") ) tree->Branch("mu_lepMVA_sip3d", "std::vector<float>", &mu_lepMVA_sip3d, buffersize);
   if( doWrite("mu_lepMVA_dxy") ) tree->Branch("mu_lepMVA_dxy", "std::vector<float>", &mu_lepMVA_dxy, buffersize);
   if( doWrite("mu_lepMVA_dz") ) tree->Branch("mu_lepMVA_dz", "std::vector<float>", &mu_lepMVA_dz, buffersize);
   if( doWrite("mu_lepMVA_mvaId") ) tree->Branch("mu_lepMVA_mvaId", "std::vector<float>", &mu_lepMVA_mvaId, buffersize);
   if( doWrite("mu_lepMVA_jetNDauChargedMVASel") ) tree->Branch("mu_lepMVA_jetNDauChargedMVASel", "std::vector<float>", &mu_lepMVA_jetNDauChargedMVASel, buffersize);

   if( doWrite("mu_conept") ) tree->Branch("mu_conept", "std::vector<float>", &mu_conept, buffersize);

   if( doWrite("mu_hasMCMatch") ) tree->Branch("mu_hasMCMatch", "std::vector<int>", &mu_hasMCMatch, buffersize);
   if( doWrite("mu_gen_pt") ) tree->Branch("mu_gen_pt", "std::vector<float>", &mu_gen_pt, buffersize);
   if( doWrite("mu_gen_eta") ) tree->Branch("mu_gen_eta", "std::vector<float>", &mu_gen_eta, buffersize);
   if( doWrite("mu_gen_phi") ) tree->Branch("mu_gen_phi", "std::vector<float>", &mu_gen_phi, buffersize);
   if( doWrite("mu_gen_m") ) tree->Branch("mu_gen_m", "std::vector<float>", &mu_gen_m, buffersize);
   if( doWrite("mu_gen_status") ) tree->Branch("mu_gen_status", "std::vector<int>", &mu_gen_status, buffersize);
   if( doWrite("mu_gen_id") ) tree->Branch("mu_gen_id", "std::vector<int>", &mu_gen_id, buffersize);
   if( doWrite("mu_gen_charge") ) tree->Branch("mu_gen_charge", "std::vector<int>", &mu_gen_charge, buffersize);
   if( doWrite("mu_gen_dr") ) tree->Branch("mu_gen_dr", "std::vector<float>", &mu_gen_dr, buffersize);

   if( doWrite("mu_hasMCMatchPAT") ) tree->Branch("mu_hasMCMatchPAT", "std::vector<int>", &mu_hasMCMatchPAT, buffersize);
   if( doWrite("mu_genPAT_pt") ) tree->Branch("mu_genPAT_pt", "std::vector<float>", &mu_genPAT_pt, buffersize);
   if( doWrite("mu_genPAT_eta") ) tree->Branch("mu_genPAT_eta", "std::vector<float>", &mu_genPAT_eta, buffersize);
   if( doWrite("mu_genPAT_phi") ) tree->Branch("mu_genPAT_phi", "std::vector<float>", &mu_genPAT_phi, buffersize);
   if( doWrite("mu_genPAT_m") ) tree->Branch("mu_genPAT_m", "std::vector<float>", &mu_genPAT_m, buffersize);
   if( doWrite("mu_genPAT_status") ) tree->Branch("mu_genPAT_status", "std::vector<int>", &mu_genPAT_status, buffersize);
   if( doWrite("mu_genPAT_id") ) tree->Branch("mu_genPAT_id", "std::vector<int>", &mu_genPAT_id, buffersize);
   if( doWrite("mu_genPAT_charge") ) tree->Branch("mu_genPAT_charge", "std::vector<int>", &mu_genPAT_charge, buffersize);
   
   if( doWrite("tau_n") ) tree->Branch("tau_n", &tau_n, "tau_n/I", buffersize);
   if( doWrite("tau_pt") ) tree->Branch("tau_pt", "std::vector<float>", &tau_pt, buffersize);
   if( doWrite("tau_eta") ) tree->Branch("tau_eta", "std::vector<float>", &tau_eta, buffersize);
   if( doWrite("tau_phi") ) tree->Branch("tau_phi", "std::vector<float>", &tau_phi, buffersize);
   if( doWrite("tau_m") ) tree->Branch("tau_m", "std::vector<float>", &tau_m, buffersize);
   if( doWrite("tau_E") ) tree->Branch("tau_E", "std::vector<float>", &tau_E, buffersize);
   if( doWrite("tau_id") ) tree->Branch("tau_id", "std::vector<int>", &tau_id, buffersize);
   if( doWrite("tau_charge") ) tree->Branch("tau_charge", "std::vector<int>", &tau_charge, buffersize);
   
   if( doWrite("tau_hasLeadChargedHadrCand") ) tree->Branch("tau_hasLeadChargedHadrCand", "std::vector<bool>", &tau_hasLeadChargedHadrCand, buffersize);
   if( doWrite("tau_leadingTrackPt") ) tree->Branch("tau_leadingTrackPt", "std::vector<float>", &tau_leadingTrackPt, buffersize);
   if( doWrite("tau_leadingTrackCharge") ) tree->Branch("tau_leadingTrackCharge", "std::vector<int>", &tau_leadingTrackCharge, buffersize);
   if( doWrite("tau_leadingTrackDz") ) tree->Branch("tau_leadingTrackDz", "std::vector<float>", &tau_leadingTrackDz, buffersize);
   if( doWrite("tau_leadingTrackDxy") ) tree->Branch("tau_leadingTrackDxy", "std::vector<float>", &tau_leadingTrackDxy, buffersize);
   
   if( doWrite("tau_decayMode") ) tree->Branch("tau_decayMode", "std::vector<int>", &tau_decayMode, buffersize);
   if( doWrite("tau_decayModeFinding") ) tree->Branch("tau_decayModeFinding", "std::vector<float>", &tau_decayModeFinding, buffersize);
//   if( doWrite("tau_decayModeFindingOldDMs") ) tree->Branch("tau_decayModeFindingOldDMs", "std::vector<float>", &tau_decayModeFindingOldDMs, buffersize);
   if( doWrite("tau_decayModeFindingNewDMs") ) tree->Branch("tau_decayModeFindingNewDMs", "std::vector<float>", &tau_decayModeFindingNewDMs, buffersize);
   
   if( doWrite("tau_puCorrPtSum") ) tree->Branch("tau_puCorrPtSum", "std::vector<float>", &tau_puCorrPtSum, buffersize);
   if( doWrite("tau_neutralIsoPtSum") ) tree->Branch("tau_neutralIsoPtSum", "std::vector<float>", &tau_neutralIsoPtSum, buffersize);
   if( doWrite("tau_chargedIsoPtSum") ) tree->Branch("tau_chargedIsoPtSum", "std::vector<float>", &tau_chargedIsoPtSum, buffersize);
   
   if( doWrite("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits") ) tree->Branch("tau_byCombinedIsolationDeltaBetaCorrRaw3Hits", "std::vector<float>", &tau_byCombinedIsolationDeltaBetaCorrRaw3Hits, buffersize);
   
   if( doWrite("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits") ) tree->Branch("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", "std::vector<float>", &tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, buffersize);
   if( doWrite("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits") ) tree->Branch("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", "std::vector<float>", &tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, buffersize);
   if( doWrite("tau_byTightCombinedIsolationDeltaBetaCorr3Hits") ) tree->Branch("tau_byTightCombinedIsolationDeltaBetaCorr3Hits", "std::vector<float>", &tau_byTightCombinedIsolationDeltaBetaCorr3Hits, buffersize);
 
   if( doWrite("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT") ) tree->Branch("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT", "std::vector<float>", &tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT, buffersize);
   if( doWrite("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT") ) tree->Branch("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT", "std::vector<float>", &tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT, buffersize);
   if( doWrite("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT") ) tree->Branch("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT", "std::vector<float>", &tau_byTightIsolationMVArun2v1DBdR03oldDMwLT, buffersize);
   if( doWrite("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT") ) tree->Branch("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT", "std::vector<float>", &tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT, buffersize);
   if( doWrite("tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT") ) tree->Branch("tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT", "std::vector<float>", &tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT, buffersize);
  
   if( doWrite("tau_againstMuonLoose3") ) tree->Branch("tau_againstMuonLoose3", "std::vector<float>", &tau_againstMuonLoose3, buffersize);
   if( doWrite("tau_againstMuonTight3") ) tree->Branch("tau_againstMuonTight3", "std::vector<float>", &tau_againstMuonTight3, buffersize);
   
   if( doWrite("tau_againstElectronVLooseMVA6") ) tree->Branch("tau_againstElectronVLooseMVA6", "std::vector<float>", &tau_againstElectronVLooseMVA6, buffersize);
   if( doWrite("tau_againstElectronLooseMVA6") ) tree->Branch("tau_againstElectronLooseMVA6", "std::vector<float>", &tau_againstElectronLooseMVA6, buffersize);
   if( doWrite("tau_againstElectronMediumMVA6") ) tree->Branch("tau_againstElectronMediumMVA6", "std::vector<float>", &tau_againstElectronMediumMVA6, buffersize);
   if( doWrite("tau_againstElectronTightMVA6") ) tree->Branch("tau_againstElectronTightMVA6", "std::vector<float>", &tau_againstElectronTightMVA6, buffersize);
   
   if( doWrite("tau_pfEssential_jet_pt") ) tree->Branch("tau_pfEssential_jet_pt", "std::vector<float>", &tau_pfEssential_jet_pt, buffersize);
   if( doWrite("tau_pfEssential_jet_eta") ) tree->Branch("tau_pfEssential_jet_eta", "std::vector<float>", &tau_pfEssential_jet_eta, buffersize);
   if( doWrite("tau_pfEssential_jet_phi") ) tree->Branch("tau_pfEssential_jet_phi", "std::vector<float>", &tau_pfEssential_jet_phi, buffersize);
   if( doWrite("tau_pfEssential_jet_m") ) tree->Branch("tau_pfEssential_jet_m", "std::vector<float>", &tau_pfEssential_jet_m, buffersize);

   if( doWrite("tau_pfEssential_jetCorr_pt") ) tree->Branch("tau_pfEssential_jetCorr_pt", "std::vector<float>", &tau_pfEssential_jetCorr_pt, buffersize);
   if( doWrite("tau_pfEssential_jetCorr_eta") ) tree->Branch("tau_pfEssential_jetCorr_eta", "std::vector<float>", &tau_pfEssential_jetCorr_eta, buffersize);
   if( doWrite("tau_pfEssential_jetCorr_phi") ) tree->Branch("tau_pfEssential_jetCorr_phi", "std::vector<float>", &tau_pfEssential_jetCorr_phi, buffersize);
   if( doWrite("tau_pfEssential_jetCorr_m") ) tree->Branch("tau_pfEssential_jetCorr_m", "std::vector<float>", &tau_pfEssential_jetCorr_m, buffersize);
   
   if( doWrite("tau_pfEssential_hasSV") ) tree->Branch("tau_pfEssential_hasSV", "std::vector<bool>", &tau_pfEssential_hasSV, buffersize);
   if( doWrite("tau_pfEssential_sv_x") ) tree->Branch("tau_pfEssential_sv_x", "std::vector<float>", &tau_pfEssential_sv_x, buffersize);
   if( doWrite("tau_pfEssential_sv_y") ) tree->Branch("tau_pfEssential_sv_y", "std::vector<float>", &tau_pfEssential_sv_y, buffersize);
   if( doWrite("tau_pfEssential_sv_z") ) tree->Branch("tau_pfEssential_sv_z", "std::vector<float>", &tau_pfEssential_sv_z, buffersize);
   
   if( doWrite("tau_pfEssential_flightLengthSig") ) tree->Branch("tau_pfEssential_flightLengthSig", "std::vector<float>", &tau_pfEssential_flightLengthSig, buffersize);
   if( doWrite("tau_pfEssential_dxy") ) tree->Branch("tau_pfEssential_dxy", "std::vector<float>", &tau_pfEssential_dxy, buffersize);
   if( doWrite("tau_pfEssential_dxy_error") ) tree->Branch("tau_pfEssential_dxy_error", "std::vector<float>", &tau_pfEssential_dxy_error, buffersize);
   if( doWrite("tau_pfEssential_dxy_Sig") ) tree->Branch("tau_pfEssential_dxy_Sig", "std::vector<float>", &tau_pfEssential_dxy_Sig, buffersize);
   
   if( doWrite("jet_n") ) tree->Branch("jet_n", &jet_n, "jet_n/I", buffersize);
   if( doWrite("jet_pt") ) tree->Branch("jet_pt", "std::vector<float>", &jet_pt, buffersize);
   if( doWrite("jet_eta") ) tree->Branch("jet_eta", "std::vector<float>", &jet_eta, buffersize);
   if( doWrite("jet_phi") ) tree->Branch("jet_phi", "std::vector<float>", &jet_phi, buffersize);
   if( doWrite("jet_m") ) tree->Branch("jet_m", "std::vector<float>", &jet_m, buffersize);
   if( doWrite("jet_E") ) tree->Branch("jet_E", "std::vector<float>", &jet_E, buffersize);

   if( doWrite("jet_ntrk") ) tree->Branch("jet_ntrk", "std::vector<int>", &jet_ntrk, buffersize);

   if( doWrite("jet_JBP") ) tree->Branch("jet_JBP", "std::vector<float>", &jet_JBP, buffersize);
   if( doWrite("jet_JP") ) tree->Branch("jet_JP", "std::vector<float>", &jet_JP, buffersize);
   if( doWrite("jet_TCHP") ) tree->Branch("jet_TCHP", "std::vector<float>", &jet_TCHP, buffersize);
   if( doWrite("jet_TCHE") ) tree->Branch("jet_TCHE", "std::vector<float>", &jet_TCHE, buffersize);
   if( doWrite("jet_SSVHE") ) tree->Branch("jet_SSVHE", "std::vector<float>", &jet_SSVHE, buffersize);
   if( doWrite("jet_SSVHP") ) tree->Branch("jet_SSVHP", "std::vector<float>", &jet_SSVHP, buffersize);
   if( doWrite("jet_CMVA") ) tree->Branch("jet_CMVA", "std::vector<float>", &jet_CMVA, buffersize);
   if( doWrite("jet_CSVv2") ) tree->Branch("jet_CSVv2", "std::vector<float>", &jet_CSVv2, buffersize);
   if( doWrite("jet_DeepCSVProbudsg") ) tree->Branch("jet_DeepCSVProbudsg", "std::vector<float>", &jet_DeepCSVProbudsg, buffersize);
   if( doWrite("jet_DeepCSVProbb") ) tree->Branch("jet_DeepCSVProbb", "std::vector<float>", &jet_DeepCSVProbb, buffersize);
   if( doWrite("jet_DeepCSVProbc") ) tree->Branch("jet_DeepCSVProbc", "std::vector<float>", &jet_DeepCSVProbc, buffersize);
   if( doWrite("jet_DeepCSVProbbb") ) tree->Branch("jet_DeepCSVProbbb", "std::vector<float>", &jet_DeepCSVProbbb, buffersize);
   if( doWrite("jet_DeepCSVProbcc") ) tree->Branch("jet_DeepCSVProbcc", "std::vector<float>", &jet_DeepCSVProbcc, buffersize);
   if( doWrite("jet_cMVAv2") ) tree->Branch("jet_cMVAv2", "std::vector<float>", &jet_cMVAv2, buffersize);
   if( doWrite("jet_CharmCvsL") ) tree->Branch("jet_CharmCvsL", "std::vector<float>", &jet_CharmCvsL, buffersize);
   if( doWrite("jet_CharmCvsB") ) tree->Branch("jet_CharmCvsB", "std::vector<float>", &jet_CharmCvsB, buffersize);
   if( doWrite("jet_partonFlavour") ) tree->Branch("jet_partonFlavour", "std::vector<int>", &jet_partonFlavour, buffersize);
   if( doWrite("jet_hadronFlavour") ) tree->Branch("jet_hadronFlavour", "std::vector<int>", &jet_hadronFlavour, buffersize);

   if( doWrite("jet_neutralHadronEnergy") ) tree->Branch("jet_neutralHadronEnergy", "std::vector<float>", &jet_neutralHadronEnergy, buffersize);
   if( doWrite("jet_neutralEmEnergy") ) tree->Branch("jet_neutralEmEnergy", "std::vector<float>", &jet_neutralEmEnergy, buffersize);
   if( doWrite("jet_chargedHadronEnergy") ) tree->Branch("jet_chargedHadronEnergy", "std::vector<float>", &jet_chargedHadronEnergy, buffersize);
   if( doWrite("jet_chargedEmEnergy") ) tree->Branch("jet_chargedEmEnergy", "std::vector<float>", &jet_chargedEmEnergy, buffersize);
   if( doWrite("jet_electronEnergy") ) tree->Branch("jet_electronEnergy", "std::vector<float>", &jet_electronEnergy, buffersize);
   if( doWrite("jet_muonEnergy") ) tree->Branch("jet_muonEnergy", "std::vector<float>", &jet_muonEnergy, buffersize);
   if( doWrite("jet_photonEnergy") ) tree->Branch("jet_photonEnergy", "std::vector<float>", &jet_photonEnergy, buffersize);

   if( doWrite("jet_chargedMultiplicity") ) tree->Branch("jet_chargedMultiplicity", "std::vector<int>", &jet_chargedMultiplicity, buffersize);
   if( doWrite("jet_neutralMultiplicity") ) tree->Branch("jet_neutralMultiplicity", "std::vector<int>", &jet_neutralMultiplicity, buffersize);
   if( doWrite("jet_chargedHadronMultiplicity") ) tree->Branch("jet_chargedHadronMultiplicity", "std::vector<int>", &jet_chargedHadronMultiplicity, buffersize);
   
   if( doWrite("jet_jetArea") ) tree->Branch("jet_jetArea", "std::vector<float>", &jet_jetArea, buffersize);
   
   if( doWrite("jet_jecFactorUncorrected") ) tree->Branch("jet_jecFactorUncorrected", "std::vector<float>", &jet_jecFactorUncorrected, buffersize);
   if( doWrite("jet_jecFactorL1FastJet") ) tree->Branch("jet_jecFactorL1FastJet", "std::vector<float>", &jet_jecFactorL1FastJet, buffersize);
   if( doWrite("jet_jecFactorL2Relative") ) tree->Branch("jet_jecFactorL2Relative", "std::vector<float>", &jet_jecFactorL2Relative, buffersize);
   if( doWrite("jet_jecFactorL3Absolute") ) tree->Branch("jet_jecFactorL3Absolute", "std::vector<float>", &jet_jecFactorL3Absolute, buffersize);
   
   if( doWrite("jet_neutralHadronEnergyFraction") ) tree->Branch("jet_neutralHadronEnergyFraction", "std::vector<float>", &jet_neutralHadronEnergyFraction, buffersize);
   if( doWrite("jet_neutralEmEnergyFraction") ) tree->Branch("jet_neutralEmEnergyFraction", "std::vector<float>", &jet_neutralEmEnergyFraction, buffersize);
   if( doWrite("jet_chargedHadronEnergyFraction") ) tree->Branch("jet_chargedHadronEnergyFraction", "std::vector<float>", &jet_chargedHadronEnergyFraction, buffersize);
   if( doWrite("jet_muonEnergyFraction") ) tree->Branch("jet_muonEnergyFraction", "std::vector<float>", &jet_muonEnergyFraction, buffersize);
   if( doWrite("jet_chargedEmEnergyFraction") ) tree->Branch("jet_chargedEmEnergyFraction", "std::vector<float>", &jet_chargedEmEnergyFraction, buffersize);
   
   if( doWrite("jet_Unc") ) tree->Branch("jet_Unc", "std::vector<float>", &jet_Unc, buffersize);
   
   if( doWrite("jet_pileupJetId") ) tree->Branch("jet_pileupJetId", "std::vector<float>", &jet_pileupJetId, buffersize);
   
//   if( doWrite("jet_looseJetID") ) tree->Branch("jet_looseJetID", "std::vector<bool>", &jet_looseJetID, buffersize);
   if( doWrite("jet_tightJetID") ) tree->Branch("jet_tightJetID", "std::vector<bool>", &jet_tightJetID, buffersize);
   if( doWrite("jet_tightLepVetoJetID") ) tree->Branch("jet_tightLepVetoJetID", "std::vector<bool>", &jet_tightLepVetoJetID, buffersize);
   if( doWrite("jet_qgtag") ) tree->Branch("jet_qgtag", "std::vector<float>", &jet_qgtag, buffersize);

   if( doWrite("jet_hasGenJet") ) tree->Branch("jet_hasGenJet", "std::vector<bool>", &jet_hasGenJet, buffersize);   
   if( doWrite("jet_genJet_pt") ) tree->Branch("jet_genJet_pt", "std::vector<float>", &jet_genJet_pt, buffersize);
   if( doWrite("jet_genJet_eta") ) tree->Branch("jet_genJet_eta", "std::vector<float>", &jet_genJet_eta, buffersize);
   if( doWrite("jet_genJet_phi") ) tree->Branch("jet_genJet_phi", "std::vector<float>", &jet_genJet_phi, buffersize);
   if( doWrite("jet_genJet_m") ) tree->Branch("jet_genJet_m", "std::vector<float>", &jet_genJet_m, buffersize);
   if( doWrite("jet_genJet_E") ) tree->Branch("jet_genJet_E", "std::vector<float>", &jet_genJet_E, buffersize);
   if( doWrite("jet_genJet_status") ) tree->Branch("jet_genJet_status", "std::vector<int>", &jet_genJet_status, buffersize);
   if( doWrite("jet_genJet_id") ) tree->Branch("jet_genJet_id", "std::vector<int>", &jet_genJet_id, buffersize);

   if( doWrite("jet_hasGenParton") ) tree->Branch("jet_hasGenParton", "std::vector<bool>", &jet_hasGenParton, buffersize);   
   if( doWrite("jet_genParton_pt") ) tree->Branch("jet_genParton_pt", "std::vector<float>", &jet_genParton_pt, buffersize);
   if( doWrite("jet_genParton_eta") ) tree->Branch("jet_genParton_eta", "std::vector<float>", &jet_genParton_eta, buffersize);
   if( doWrite("jet_genParton_phi") ) tree->Branch("jet_genParton_phi", "std::vector<float>", &jet_genParton_phi, buffersize);
   if( doWrite("jet_genParton_m") ) tree->Branch("jet_genParton_m", "std::vector<float>", &jet_genParton_m, buffersize);
   if( doWrite("jet_genParton_E") ) tree->Branch("jet_genParton_E", "std::vector<float>", &jet_genParton_E, buffersize);
   if( doWrite("jet_genParton_status") ) tree->Branch("jet_genParton_status", "std::vector<int>", &jet_genParton_status, buffersize);
   if( doWrite("jet_genParton_id") ) tree->Branch("jet_genParton_id", "std::vector<int>", &jet_genParton_id, buffersize);

   if( doWrite("jetPuppi_do") )
     {	   
	if( doWrite("jetPuppi_n") ) tree->Branch("jetPuppi_n", &jetPuppi_n, "jetPuppi_n/I", buffersize);
	if( doWrite("jetPuppi_pt") ) tree->Branch("jetPuppi_pt", "std::vector<float>", &jetPuppi_pt, buffersize);
	if( doWrite("jetPuppi_eta") ) tree->Branch("jetPuppi_eta", "std::vector<float>", &jetPuppi_eta, buffersize);
	if( doWrite("jetPuppi_phi") ) tree->Branch("jetPuppi_phi", "std::vector<float>", &jetPuppi_phi, buffersize);
	if( doWrite("jetPuppi_m") ) tree->Branch("jetPuppi_m", "std::vector<float>", &jetPuppi_m, buffersize);
	if( doWrite("jetPuppi_E") ) tree->Branch("jetPuppi_E", "std::vector<float>", &jetPuppi_E, buffersize);
	
	if( doWrite("jetPuppi_ntrk") ) tree->Branch("jetPuppi_ntrk", "std::vector<int>", &jetPuppi_ntrk, buffersize);
	
	if( doWrite("jetPuppi_JBP") ) tree->Branch("jetPuppi_JBP", "std::vector<float>", &jetPuppi_JBP, buffersize);
	if( doWrite("jetPuppi_JP") ) tree->Branch("jetPuppi_JP", "std::vector<float>", &jetPuppi_JP, buffersize);
	if( doWrite("jetPuppi_TCHP") ) tree->Branch("jetPuppi_TCHP", "std::vector<float>", &jetPuppi_TCHP, buffersize);
	if( doWrite("jetPuppi_TCHE") ) tree->Branch("jetPuppi_TCHE", "std::vector<float>", &jetPuppi_TCHE, buffersize);
	if( doWrite("jetPuppi_SSVHE") ) tree->Branch("jetPuppi_SSVHE", "std::vector<float>", &jetPuppi_SSVHE, buffersize);
	if( doWrite("jetPuppi_SSVHP") ) tree->Branch("jetPuppi_SSVHP", "std::vector<float>", &jetPuppi_SSVHP, buffersize);
	if( doWrite("jetPuppi_CMVA") ) tree->Branch("jetPuppi_CMVA", "std::vector<float>", &jetPuppi_CMVA, buffersize);
	if( doWrite("jetPuppi_CSVv2") ) tree->Branch("jetPuppi_CSVv2", "std::vector<float>", &jetPuppi_CSVv2, buffersize);
	if( doWrite("jetPuppi_partonFlavour") ) tree->Branch("jetPuppi_partonFlavour", "std::vector<int>", &jetPuppi_partonFlavour, buffersize);
	if( doWrite("jetPuppi_hadronFlavour") ) tree->Branch("jetPuppi_hadronFlavour", "std::vector<int>", &jetPuppi_hadronFlavour, buffersize);
	
	if( doWrite("jetPuppi_neutralHadronEnergy") ) tree->Branch("jetPuppi_neutralHadronEnergy", "std::vector<float>", &jetPuppi_neutralHadronEnergy, buffersize);
	if( doWrite("jetPuppi_neutralEmEnergy") ) tree->Branch("jetPuppi_neutralEmEnergy", "std::vector<float>", &jetPuppi_neutralEmEnergy, buffersize);
	if( doWrite("jetPuppi_chargedHadronEnergy") ) tree->Branch("jetPuppi_chargedHadronEnergy", "std::vector<float>", &jetPuppi_chargedHadronEnergy, buffersize);
	if( doWrite("jetPuppi_chargedEmEnergy") ) tree->Branch("jetPuppi_chargedEmEnergy", "std::vector<float>", &jetPuppi_chargedEmEnergy, buffersize);
	if( doWrite("jetPuppi_electronEnergy") ) tree->Branch("jetPuppi_electronEnergy", "std::vector<float>", &jetPuppi_electronEnergy, buffersize);
	if( doWrite("jetPuppi_muonEnergy") ) tree->Branch("jetPuppi_muonEnergy", "std::vector<float>", &jetPuppi_muonEnergy, buffersize);
	if( doWrite("jetPuppi_photonEnergy") ) tree->Branch("jetPuppi_photonEnergy", "std::vector<float>", &jetPuppi_photonEnergy, buffersize);
	
	if( doWrite("jetPuppi_chargedMultiplicity") ) tree->Branch("jetPuppi_chargedMultiplicity", "std::vector<int>", &jetPuppi_chargedMultiplicity, buffersize);
	if( doWrite("jetPuppi_neutralMultiplicity") ) tree->Branch("jetPuppi_neutralMultiplicity", "std::vector<int>", &jetPuppi_neutralMultiplicity, buffersize);
	if( doWrite("jetPuppi_chargedHadronMultiplicity") ) tree->Branch("jetPuppi_chargedHadronMultiplicity", "std::vector<int>", &jetPuppi_chargedHadronMultiplicity, buffersize);
	
	if( doWrite("jetPuppi_jecFactorUncorrected") ) tree->Branch("jetPuppi_jecFactorUncorrected", "std::vector<float>", &jetPuppi_jecFactorUncorrected, buffersize);
	if( doWrite("jetPuppi_jecFactorL1FastJet") ) tree->Branch("jetPuppi_jecFactorL1FastJet", "std::vector<float>", &jetPuppi_jecFactorL1FastJet, buffersize);
	if( doWrite("jetPuppi_jecFactorL2Relative") ) tree->Branch("jetPuppi_jecFactorL2Relative", "std::vector<float>", &jetPuppi_jecFactorL2Relative, buffersize);
	if( doWrite("jetPuppi_jecFactorL3Absolute") ) tree->Branch("jetPuppi_jecFactorL3Absolute", "std::vector<float>", &jetPuppi_jecFactorL3Absolute, buffersize);
	
	if( doWrite("jetPuppi_pileupJetId") ) tree->Branch("jetPuppi_pileupJetId", "std::vector<float>", &jetPuppi_pileupJetId, buffersize);
	
	if( doWrite("jetPuppi_hasGenJet") ) tree->Branch("jetPuppi_hasGenJet", "std::vector<bool>", &jetPuppi_hasGenJet, buffersize);   
	if( doWrite("jetPuppi_genJet_pt") ) tree->Branch("jetPuppi_genJet_pt", "std::vector<float>", &jetPuppi_genJet_pt, buffersize);
	if( doWrite("jetPuppi_genJet_eta") ) tree->Branch("jetPuppi_genJet_eta", "std::vector<float>", &jetPuppi_genJet_eta, buffersize);
	if( doWrite("jetPuppi_genJet_phi") ) tree->Branch("jetPuppi_genJet_phi", "std::vector<float>", &jetPuppi_genJet_phi, buffersize);
	if( doWrite("jetPuppi_genJet_m") ) tree->Branch("jetPuppi_genJet_m", "std::vector<float>", &jetPuppi_genJet_m, buffersize);
	if( doWrite("jetPuppi_genJet_E") ) tree->Branch("jetPuppi_genJet_E", "std::vector<float>", &jetPuppi_genJet_E, buffersize);
	if( doWrite("jetPuppi_genJet_status") ) tree->Branch("jetPuppi_genJet_status", "std::vector<int>", &jetPuppi_genJet_status, buffersize);
	if( doWrite("jetPuppi_genJet_id") ) tree->Branch("jetPuppi_genJet_id", "std::vector<int>", &jetPuppi_genJet_id, buffersize);       
	
	if( doWrite("jetPuppi_hasGenParton") ) tree->Branch("jetPuppi_hasGenJet", "std::vector<bool>", &jetPuppi_hasGenParton, buffersize);   
	if( doWrite("jetPuppi_genParton_pt") ) tree->Branch("jetPuppi_genParton_pt", "std::vector<float>", &jetPuppi_genParton_pt, buffersize);
	if( doWrite("jetPuppi_genParton_eta") ) tree->Branch("jetPuppi_genParton_eta", "std::vector<float>", &jetPuppi_genParton_eta, buffersize);
	if( doWrite("jetPuppi_genParton_phi") ) tree->Branch("jetPuppi_genParton_phi", "std::vector<float>", &jetPuppi_genParton_phi, buffersize);
	if( doWrite("jetPuppi_genParton_m") ) tree->Branch("jetPuppi_genParton_m", "std::vector<float>", &jetPuppi_genParton_m, buffersize);
	if( doWrite("jetPuppi_genParton_E") ) tree->Branch("jetPuppi_genParton_E", "std::vector<float>", &jetPuppi_genParton_E, buffersize);
	if( doWrite("jetPuppi_genParton_status") ) tree->Branch("jetPuppi_genParton_status", "std::vector<int>", &jetPuppi_genParton_status, buffersize);
	if( doWrite("jetPuppi_genParton_id") ) tree->Branch("jetPuppi_genParton_id", "std::vector<int>", &jetPuppi_genParton_id, buffersize);
     }   
   
   
   //------------------------
   //  ak8 collection
   //------------------------

   if( doWrite("ak8jet_do") )
     {	
	if( doWrite("ak8jet_n") ) tree->Branch("ak8jet_n", &ak8jet_n, "ak8jet_n/I", buffersize);
	if( doWrite("ak8jet_pt") ) tree->Branch("ak8jet_pt", "std::vector<float>", &ak8jet_pt, buffersize);
	if( doWrite("ak8jet_eta") ) tree->Branch("ak8jet_eta", "std::vector<float>", &ak8jet_eta, buffersize);
	if( doWrite("ak8jet_phi") ) tree->Branch("ak8jet_phi", "std::vector<float>", &ak8jet_phi, buffersize);
	if( doWrite("ak8jet_m") ) tree->Branch("ak8jet_m", "std::vector<float>", &ak8jet_m, buffersize);
	if( doWrite("ak8jet_E") ) tree->Branch("ak8jet_E", "std::vector<float>", &ak8jet_E, buffersize);
	
	if( doWrite("ak8jet_ntrk") ) tree->Branch("ak8jet_ntrk", "std::vector<int>", &ak8jet_ntrk, buffersize);
	
	if( doWrite("ak8jet_JBP") ) tree->Branch("ak8jet_JBP", "std::vector<float>", &ak8jet_JBP, buffersize);
	if( doWrite("ak8jet_JP") ) tree->Branch("ak8jet_JP", "std::vector<float>", &ak8jet_JP, buffersize);
	if( doWrite("ak8jet_TCHP") ) tree->Branch("ak8jet_TCHP", "std::vector<float>", &ak8jet_TCHP, buffersize);
	if( doWrite("ak8jet_TCHE") ) tree->Branch("ak8jet_TCHE", "std::vector<float>", &ak8jet_TCHE, buffersize);
	if( doWrite("ak8jet_SSVHE") ) tree->Branch("ak8jet_SSVHE", "std::vector<float>", &ak8jet_SSVHE, buffersize);
	if( doWrite("ak8jet_SSVHP") ) tree->Branch("ak8jet_SSVHP", "std::vector<float>", &ak8jet_SSVHP, buffersize);
	if( doWrite("ak8jet_CMVA") ) tree->Branch("ak8jet_CMVA", "std::vector<float>", &ak8jet_CMVA, buffersize);
	
	if( doWrite("ak8jet_CSVv2") ) tree->Branch("ak8jet_CSVv2", "std::vector<float>", &ak8jet_CSVv2, buffersize);
	if( doWrite("ak8jet_partonFlavour") ) tree->Branch("ak8jet_partonFlavour", "std::vector<int>", &ak8jet_partonFlavour, buffersize);
	if( doWrite("ak8jet_hadronFlavour") ) tree->Branch("ak8jet_hadronFlavour", "std::vector<int>", &ak8jet_hadronFlavour, buffersize);
	
	if( doWrite("ak8jet_neutralHadronEnergy") ) tree->Branch("ak8jet_neutralHadronEnergy", "std::vector<float>", &ak8jet_neutralHadronEnergy, buffersize);
	if( doWrite("ak8jet_neutralEmEnergy") ) tree->Branch("ak8jet_neutralEmEnergy", "std::vector<float>", &ak8jet_neutralEmEnergy, buffersize);
	if( doWrite("ak8jet_chargedHadronEnergy") ) tree->Branch("ak8jet_chargedHadronEnergy", "std::vector<float>", &ak8jet_chargedHadronEnergy, buffersize);
	if( doWrite("ak8jet_chargedEmEnergy") ) tree->Branch("ak8jet_chargedEmEnergy", "std::vector<float>", &ak8jet_chargedEmEnergy, buffersize);
	if( doWrite("ak8jet_electronEnergy") ) tree->Branch("ak8jet_electronEnergy", "std::vector<float>", &ak8jet_electronEnergy, buffersize);
	if( doWrite("ak8jet_muonEnergy") ) tree->Branch("ak8jet_muonEnergy", "std::vector<float>", &ak8jet_muonEnergy, buffersize);
	if( doWrite("ak8jet_photonEnergy") ) tree->Branch("ak8jet_photonEnergy", "std::vector<float>", &ak8jet_photonEnergy, buffersize);
	
	if( doWrite("ak8jet_chargedMultiplicity") ) tree->Branch("ak8jet_chargedMultiplicity", "std::vector<int>", &ak8jet_chargedMultiplicity, buffersize);
	if( doWrite("ak8jet_neutralMultiplicity") ) tree->Branch("ak8jet_neutralMultiplicity", "std::vector<int>", &ak8jet_neutralMultiplicity, buffersize);
	if( doWrite("ak8jet_chargedHadronMultiplicity") ) tree->Branch("ak8jet_chargedHadronMultiplicity", "std::vector<int>", &ak8jet_chargedHadronMultiplicity, buffersize);
	
	if( doWrite("ak8jet_jetArea") ) tree->Branch("ak8jet_jetArea", "std::vector<float>", &ak8jet_jetArea, buffersize);
	
	if( doWrite("ak8jet_jecFactorUncorrected") ) tree->Branch("ak8jet_jecFactorUncorrected", "std::vector<float>", &ak8jet_jecFactorUncorrected, buffersize);
	if( doWrite("ak8jet_jecFactorL1FastJet") ) tree->Branch("ak8jet_jecFactorL1FastJet", "std::vector<float>", &ak8jet_jecFactorL1FastJet, buffersize);
	if( doWrite("ak8jet_jecFactorL2Relative") ) tree->Branch("ak8jet_jecFactorL2Relative", "std::vector<float>", &ak8jet_jecFactorL2Relative, buffersize);
	if( doWrite("ak8jet_jecFactorL3Absolute") ) tree->Branch("ak8jet_jecFactorL3Absolute", "std::vector<float>", &ak8jet_jecFactorL3Absolute, buffersize);
	
	if( doWrite("ak8jet_pileupJetId") ) tree->Branch("ak8jet_pileupJetId", "std::vector<float>", &ak8jet_pileupJetId, buffersize);
	
	//   if( doWrite("ak8jet_looseJetID") ) tree->Branch("ak8jet_looseJetID", "std::vector<bool>", &ak8jet_looseJetID, buffersize);
	if( doWrite("ak8jet_tightJetID") ) tree->Branch("ak8jet_tightJetID", "std::vector<bool>", &ak8jet_tightJetID, buffersize);
	
	if( doWrite("ak8jet_hasGenJet") ) tree->Branch("ak8jet_hasGenJet", "std::vector<bool>", &ak8jet_hasGenJet, buffersize);
	if( doWrite("ak8jet_genJet_pt") ) tree->Branch("ak8jet_genJet_pt", "std::vector<float>", &ak8jet_genJet_pt, buffersize);
	if( doWrite("ak8jet_genJet_eta") ) tree->Branch("ak8jet_genJet_eta", "std::vector<float>", &ak8jet_genJet_eta, buffersize);
	if( doWrite("ak8jet_genJet_phi") ) tree->Branch("ak8jet_genJet_phi", "std::vector<float>", &ak8jet_genJet_phi, buffersize);
	if( doWrite("ak8jet_genJet_m") ) tree->Branch("ak8jet_genJet_m", "std::vector<float>", &ak8jet_genJet_m, buffersize);
	if( doWrite("ak8jet_genJet_E") ) tree->Branch("ak8jet_genJet_E", "std::vector<float>", &ak8jet_genJet_E, buffersize);
	if( doWrite("ak8jet_genJet_status") ) tree->Branch("ak8jet_genJet_status", "std::vector<int>", &ak8jet_genJet_status, buffersize);
	if( doWrite("ak8jet_genJet_id") ) tree->Branch("ak8jet_genJet_id", "std::vector<int>", &ak8jet_genJet_id, buffersize);
	
	if( doWrite("ak8jet_hasGenParton") ) tree->Branch("ak8jet_hasGenParton", "std::vector<bool>", &ak8jet_hasGenParton, buffersize);
	if( doWrite("ak8jet_genParton_pt") ) tree->Branch("ak8jet_genParton_pt", "std::vector<float>", &ak8jet_genParton_pt, buffersize);
	if( doWrite("ak8jet_genParton_eta") ) tree->Branch("ak8jet_genParton_eta", "std::vector<float>", &ak8jet_genParton_eta, buffersize);
	if( doWrite("ak8jet_genParton_phi") ) tree->Branch("ak8jet_genParton_phi", "std::vector<float>", &ak8jet_genParton_phi, buffersize);
	if( doWrite("ak8jet_genParton_m") ) tree->Branch("ak8jet_genParton_m", "std::vector<float>", &ak8jet_genParton_m, buffersize);
	if( doWrite("ak8jet_genParton_E") ) tree->Branch("ak8jet_genParton_E", "std::vector<float>", &ak8jet_genParton_E, buffersize);
	if( doWrite("ak8jet_genParton_status") ) tree->Branch("ak8jet_genParton_status", "std::vector<int>", &ak8jet_genParton_status, buffersize);
	if( doWrite("ak8jet_genParton_id") ) tree->Branch("ak8jet_genParton_id", "std::vector<int>", &ak8jet_genParton_id, buffersize);
	
	if( doWrite("ak8jet_tau1") ) tree->Branch("ak8jet_tau1", "std::vector<float>", &ak8jet_tau1, buffersize);
	if( doWrite("ak8jet_tau2") ) tree->Branch("ak8jet_tau2", "std::vector<float>", &ak8jet_tau2, buffersize);
	if( doWrite("ak8jet_tau3") ) tree->Branch("ak8jet_tau3", "std::vector<float>", &ak8jet_tau3, buffersize);
	if( doWrite("ak8jet_softdrop_mass") ) tree->Branch("ak8jet_softdrop_mass", "std::vector<float>", &ak8jet_softdrop_mass, buffersize);
	if( doWrite("ak8jet_trimmed_mass") ) tree->Branch("ak8jet_trimmed_mass", "std::vector<float>", &ak8jet_trimmed_mass, buffersize);
	if( doWrite("ak8jet_pruned_mass") ) tree->Branch("ak8jet_pruned_mass", "std::vector<float>", &ak8jet_pruned_mass, buffersize);
	if( doWrite("ak8jet_filtered_mass") ) tree->Branch("ak8jet_filtered_mass", "std::vector<float>", &ak8jet_filtered_mass, buffersize);
	if( doWrite("ak8jet_minMass") ) tree->Branch("ak8jet_minMass", "std::vector<float>", &ak8jet_minMass, buffersize);
	if( doWrite("ak8jet_topMass") ) tree->Branch("ak8jet_topMass", "std::vector<float>", &ak8jet_topMass, buffersize);
	if( doWrite("ak8jet_nSubJets") ) tree->Branch("ak8jet_nSubJets", "std::vector<int>", &ak8jet_nSubJets, buffersize);
     }   
   
   //------------------------
   //  ak10 collection
   //------------------------

   if( doWrite("ak10jet_do") )
     {	
	if( doWrite("ak10jet_n") ) tree->Branch("ak10jet_n", &ak10jet_n, "ak10jet_n/I", buffersize);
	if( doWrite("ak10jet_pt") ) tree->Branch("ak10jet_pt", "std::vector<float>", &ak10jet_pt, buffersize);
	if( doWrite("ak10jet_eta") ) tree->Branch("ak10jet_eta", "std::vector<float>", &ak10jet_eta, buffersize);
	if( doWrite("ak10jet_phi") ) tree->Branch("ak10jet_phi", "std::vector<float>", &ak10jet_phi, buffersize);
	if( doWrite("ak10jet_m") ) tree->Branch("ak10jet_m", "std::vector<float>", &ak10jet_m, buffersize);
	if( doWrite("ak10jet_E") ) tree->Branch("ak10jet_E", "std::vector<float>", &ak10jet_E, buffersize);
	
	if( doWrite("ak10jet_ntrk") ) tree->Branch("ak10jet_ntrk", "std::vector<int>", &ak10jet_ntrk, buffersize);
	
	if( doWrite("ak10jet_JBP") ) tree->Branch("ak10jet_JBP", "std::vector<float>", &ak10jet_JBP, buffersize);
	if( doWrite("ak10jet_JP") ) tree->Branch("ak10jet_JP", "std::vector<float>", &ak10jet_JP, buffersize);
	if( doWrite("ak10jet_TCHP") ) tree->Branch("ak10jet_TCHP", "std::vector<float>", &ak10jet_TCHP, buffersize);
	if( doWrite("ak10jet_TCHE") ) tree->Branch("ak10jet_TCHE", "std::vector<float>", &ak10jet_TCHE, buffersize);
	if( doWrite("ak10jet_SSVHE") ) tree->Branch("ak10jet_SSVHE", "std::vector<float>", &ak10jet_SSVHE, buffersize);
	if( doWrite("ak10jet_SSVHP") ) tree->Branch("ak10jet_SSVHP", "std::vector<float>", &ak10jet_SSVHP, buffersize);
	if( doWrite("ak10jet_CMVA") ) tree->Branch("ak10jet_CMVA", "std::vector<float>", &ak10jet_CMVA, buffersize);
	
	if( doWrite("ak10jet_CSVv2") ) tree->Branch("ak10jet_CSVv2", "std::vector<float>", &ak10jet_CSVv2, buffersize);
	if( doWrite("ak10jet_partonFlavour") ) tree->Branch("ak10jet_partonFlavour", "std::vector<int>", &ak10jet_partonFlavour, buffersize);
	if( doWrite("ak10jet_hadronFlavour") ) tree->Branch("ak10jet_hadronFlavour", "std::vector<int>", &ak10jet_hadronFlavour, buffersize);
	
	if( doWrite("ak10jet_neutralHadronEnergy") ) tree->Branch("ak10jet_neutralHadronEnergy", "std::vector<float>", &ak10jet_neutralHadronEnergy, buffersize);
	if( doWrite("ak10jet_neutralEmEnergy") ) tree->Branch("ak10jet_neutralEmEnergy", "std::vector<float>", &ak10jet_neutralEmEnergy, buffersize);
	if( doWrite("ak10jet_chargedHadronEnergy") ) tree->Branch("ak10jet_chargedHadronEnergy", "std::vector<float>", &ak10jet_chargedHadronEnergy, buffersize);
	if( doWrite("ak10jet_chargedEmEnergy") ) tree->Branch("ak10jet_chargedEmEnergy", "std::vector<float>", &ak10jet_chargedEmEnergy, buffersize);
	if( doWrite("ak10jet_electronEnergy") ) tree->Branch("ak10jet_electronEnergy", "std::vector<float>", &ak10jet_electronEnergy, buffersize);
	if( doWrite("ak10jet_muonEnergy") ) tree->Branch("ak10jet_muonEnergy", "std::vector<float>", &ak10jet_muonEnergy, buffersize);
	if( doWrite("ak10jet_photonEnergy") ) tree->Branch("ak10jet_photonEnergy", "std::vector<float>", &ak10jet_photonEnergy, buffersize);
	
	if( doWrite("ak10jet_chargedMultiplicity") ) tree->Branch("ak10jet_chargedMultiplicity", "std::vector<int>", &ak10jet_chargedMultiplicity, buffersize);
	if( doWrite("ak10jet_neutralMultiplicity") ) tree->Branch("ak10jet_neutralMultiplicity", "std::vector<int>", &ak10jet_neutralMultiplicity, buffersize);
	if( doWrite("ak10jet_chargedHadronMultiplicity") ) tree->Branch("ak10jet_chargedHadronMultiplicity", "std::vector<int>", &ak10jet_chargedHadronMultiplicity, buffersize);
	
	if( doWrite("ak10jet_jetArea") ) tree->Branch("ak10jet_jetArea", "std::vector<float>", &ak10jet_jetArea, buffersize);
	
	if( doWrite("ak10jet_jecFactorUncorrected") ) tree->Branch("ak10jet_jecFactorUncorrected", "std::vector<float>", &ak10jet_jecFactorUncorrected, buffersize);
	if( doWrite("ak10jet_jecFactorL1FastJet") ) tree->Branch("ak10jet_jecFactorL1FastJet", "std::vector<float>", &ak10jet_jecFactorL1FastJet, buffersize);
	if( doWrite("ak10jet_jecFactorL2Relative") ) tree->Branch("ak10jet_jecFactorL2Relative", "std::vector<float>", &ak10jet_jecFactorL2Relative, buffersize);
	if( doWrite("ak10jet_jecFactorL3Absolute") ) tree->Branch("ak10jet_jecFactorL3Absolute", "std::vector<float>", &ak10jet_jecFactorL3Absolute, buffersize);
	
	if( doWrite("ak10jet_pileupJetId") ) tree->Branch("ak10jet_pileupJetId", "std::vector<float>", &ak10jet_pileupJetId, buffersize);
	
	//   if( doWrite("ak10jet_looseJetID") ) tree->Branch("ak10jet_looseJetID", "std::vector<bool>", &ak10jet_looseJetID, buffersize);
	if( doWrite("ak10jet_tightJetID") ) tree->Branch("ak10jet_tightJetID", "std::vector<bool>", &ak10jet_tightJetID, buffersize);
	
	if( doWrite("ak10jet_hasGenJet") ) tree->Branch("ak10jet_hasGenJet", "std::vector<bool>", &ak10jet_hasGenJet, buffersize);
	if( doWrite("ak10jet_genJet_pt") ) tree->Branch("ak10jet_genJet_pt", "std::vector<float>", &ak10jet_genJet_pt, buffersize);
	if( doWrite("ak10jet_genJet_eta") ) tree->Branch("ak10jet_genJet_eta", "std::vector<float>", &ak10jet_genJet_eta, buffersize);
	if( doWrite("ak10jet_genJet_phi") ) tree->Branch("ak10jet_genJet_phi", "std::vector<float>", &ak10jet_genJet_phi, buffersize);
	if( doWrite("ak10jet_genJet_m") ) tree->Branch("ak10jet_genJet_m", "std::vector<float>", &ak10jet_genJet_m, buffersize);
	if( doWrite("ak10jet_genJet_E") ) tree->Branch("ak10jet_genJet_E", "std::vector<float>", &ak10jet_genJet_E, buffersize);
	if( doWrite("ak10jet_genJet_status") ) tree->Branch("ak10jet_genJet_status", "std::vector<int>", &ak10jet_genJet_status, buffersize);
	if( doWrite("ak10jet_genJet_id") ) tree->Branch("ak10jet_genJet_id", "std::vector<int>", &ak10jet_genJet_id, buffersize);
	
	if( doWrite("ak10jet_hasGenParton") ) tree->Branch("ak10jet_hasGenParton", "std::vector<bool>", &ak10jet_hasGenParton, buffersize);
	if( doWrite("ak10jet_genParton_pt") ) tree->Branch("ak10jet_genParton_pt", "std::vector<float>", &ak10jet_genParton_pt, buffersize);
	if( doWrite("ak10jet_genParton_eta") ) tree->Branch("ak10jet_genParton_eta", "std::vector<float>", &ak10jet_genParton_eta, buffersize);
	if( doWrite("ak10jet_genParton_phi") ) tree->Branch("ak10jet_genParton_phi", "std::vector<float>", &ak10jet_genParton_phi, buffersize);
	if( doWrite("ak10jet_genParton_m") ) tree->Branch("ak10jet_genParton_m", "std::vector<float>", &ak10jet_genParton_m, buffersize);
	if( doWrite("ak10jet_genParton_E") ) tree->Branch("ak10jet_genParton_E", "std::vector<float>", &ak10jet_genParton_E, buffersize);
	if( doWrite("ak10jet_genParton_status") ) tree->Branch("ak10jet_genParton_status", "std::vector<int>", &ak10jet_genParton_status, buffersize);
	if( doWrite("ak10jet_genParton_id") ) tree->Branch("ak10jet_genParton_id", "std::vector<int>", &ak10jet_genParton_id, buffersize);
	
	if( doWrite("ak10jet_tau1") ) tree->Branch("ak10jet_tau1", "std::vector<float>", &ak10jet_tau1, buffersize);
	if( doWrite("ak10jet_tau2") ) tree->Branch("ak10jet_tau2", "std::vector<float>", &ak10jet_tau2, buffersize);
	if( doWrite("ak10jet_tau3") ) tree->Branch("ak10jet_tau3", "std::vector<float>", &ak10jet_tau3, buffersize);
	if( doWrite("ak10jet_softdrop_mass") ) tree->Branch("ak10jet_softdrop_mass", "std::vector<float>", &ak10jet_softdrop_mass, buffersize);
	if( doWrite("ak10jet_trimmed_mass") ) tree->Branch("ak10jet_trimmed_mass", "std::vector<float>", &ak10jet_trimmed_mass, buffersize);
	if( doWrite("ak10jet_pruned_mass") ) tree->Branch("ak10jet_pruned_mass", "std::vector<float>", &ak10jet_pruned_mass, buffersize);
	if( doWrite("ak10jet_filtered_mass") ) tree->Branch("ak10jet_filtered_mass", "std::vector<float>", &ak10jet_filtered_mass, buffersize);
	if( doWrite("ak10jet_minMass") ) tree->Branch("ak10jet_minMass", "std::vector<float>", &ak10jet_minMass, buffersize);
	if( doWrite("ak10jet_topMass") ) tree->Branch("ak10jet_topMass", "std::vector<float>", &ak10jet_topMass, buffersize);
	if( doWrite("ak10jet_nSubJets") ) tree->Branch("ak10jet_nSubJets", "std::vector<int>", &ak10jet_nSubJets, buffersize);
     }   
   
   //------------------------
   //  genJet collection
   //------------------------

   
   if( doWrite("genJet_n") ) tree->Branch("genJet_n", &genJet_n, "genJet_n/I", buffersize);
   if( doWrite("genJet_pt") ) tree->Branch("genJet_pt", "std::vector<float>", &genJet_pt, buffersize);
   if( doWrite("genJet_eta") ) tree->Branch("genJet_eta", "std::vector<float>", &genJet_eta, buffersize);
   if( doWrite("genJet_phi") ) tree->Branch("genJet_phi", "std::vector<float>", &genJet_phi, buffersize);
   if( doWrite("genJet_m") ) tree->Branch("genJet_m", "std::vector<float>", &genJet_m, buffersize);
   if( doWrite("genJet_E") ) tree->Branch("genJet_E", "std::vector<float>", &genJet_E, buffersize);
   if( doWrite("genJet_emEnergy") ) tree->Branch("genJet_emEnergy", "std::vector<float>", &genJet_emEnergy, buffersize);
   if( doWrite("genJet_hadEnergy") ) tree->Branch("genJet_hadEnergy", "std::vector<float>", &genJet_hadEnergy, buffersize);
   if( doWrite("genJet_invisibleEnergy") ) tree->Branch("genJet_invisibleEnergy", "std::vector<float>", &genJet_invisibleEnergy, buffersize);
   if( doWrite("genJet_auxiliaryEnergy") ) tree->Branch("genJet_auxiliaryEnergy", "std::vector<float>", &genJet_auxiliaryEnergy, buffersize);
   if( doWrite("genJet_flavour") ) tree->Branch("genJet_flavour", "std::vector<int>", &genJet_flavour, buffersize);
    
   //------------------------
   //  genTTXJet collection
   //------------------------
   
   if( doWrite("genTTXJet_n") ) tree->Branch("genTTXJet_n", &genTTXJet_n, "genTTXJet_n/I", buffersize);
   if( doWrite("genTTXJet_pt") ) tree->Branch("genTTXJet_pt", "std::vector<float>", &genTTXJet_pt, buffersize);
   if( doWrite("genTTXJet_eta") ) tree->Branch("genTTXJet_eta", "std::vector<float>", &genTTXJet_eta, buffersize);
   if( doWrite("genTTXJet_phi") ) tree->Branch("genTTXJet_phi", "std::vector<float>", &genTTXJet_phi, buffersize);
   if( doWrite("genTTXJet_m") ) tree->Branch("genTTXJet_m", "std::vector<float>", &genTTXJet_m, buffersize);
   if( doWrite("genTTXJet_E") ) tree->Branch("genTTXJet_E", "std::vector<float>", &genTTXJet_E, buffersize);
   if( doWrite("genTTXJet_flavour") ) tree->Branch("genTTXJet_flavour", "std::vector<int>", &genTTXJet_flavour, buffersize);


   if( doWrite("pfcand_do") )
     {	
	if( doWrite("pfcand_n") ) tree->Branch("pfcand_n", &pfcand_n, "pfcand_n/I", buffersize);
	if( doWrite("pfcand_pt") ) tree->Branch("pfcand_pt", "std::vector<float>", &pfcand_pt, buffersize);
	if( doWrite("pfcand_eta") ) tree->Branch("pfcand_eta", "std::vector<float>", &pfcand_eta, buffersize);
	if( doWrite("pfcand_phi") ) tree->Branch("pfcand_phi", "std::vector<float>", &pfcand_phi, buffersize);
	if( doWrite("pfcand_E") ) tree->Branch("pfcand_E", "std::vector<float>", &pfcand_E, buffersize);
	if( doWrite("pfcand_charge") ) tree->Branch("pfcand_charge", "std::vector<float>", &pfcand_charge, buffersize);
	if( doWrite("pfcand_id") ) tree->Branch("pfcand_id", "std::vector<int>", &pfcand_id, buffersize);
	if( doWrite("pfcand_dz") ) tree->Branch("pfcand_dz", "std::vector<float>", &pfcand_dz, buffersize);
	if( doWrite("pfcand_trackIso") ) tree->Branch("pfcand_trackIso", "std::vector<float>", &pfcand_trackIso, buffersize);
	
	if( doWrite("pfch_loose_n") ) tree->Branch("pfch_loose_n",  &pfch_loose_n, "pfch_loose_n/I",buffersize);
	if( doWrite("pfch_loose_sumpt") ) tree->Branch("pfch_loose_sumpt",  &pfch_loose_sumpt, "pfch_loose_sumpt/F",buffersize);
	if( doWrite("pfch_tight_n") ) tree->Branch("pfch_tight_n",  &pfch_tight_n, "pfch_tight_n/I",buffersize);
	if( doWrite("pfch_tight_sumpt") ) tree->Branch("pfch_tight_sumpt",  &pfch_tight_sumpt, "pfch_tight_sumpt/F",buffersize);
     }      

   if( doWrite("mc_truth_tth") )
     {
	tree->Branch("mc_truth_tth_channel", &mc_truth_tth_channel, "mc_truth_tth_channel/I", buffersize);

	if( doWrite("mc_truth_p4") )
	  {	
	     tree->Branch("mc_truth_h0_p4", "TLorentzVector", &mc_truth_h0_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0W1_p4", "TLorentzVector", &mc_truth_h0W1_p4, buffersize);
	     tree->Branch("mc_truth_h0W2_p4", "TLorentzVector", &mc_truth_h0W2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wl1_p4", "TLorentzVector", &mc_truth_h0Wl1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnu1_p4", "TLorentzVector", &mc_truth_h0Wnu1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtau1_p4", "TLorentzVector", &mc_truth_h0Wtau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnutau1_p4", "TLorentzVector", &mc_truth_h0Wnutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaul1_p4", "TLorentzVector", &mc_truth_h0Wtaul1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunu1_p4", "TLorentzVector", &mc_truth_h0Wtaunu1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunutau1_p4", "TLorentzVector", &mc_truth_h0Wtaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wl2_p4", "TLorentzVector", &mc_truth_h0Wl2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnu2_p4", "TLorentzVector", &mc_truth_h0Wnu2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtau2_p4", "TLorentzVector", &mc_truth_h0Wtau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnutau2_p4", "TLorentzVector", &mc_truth_h0Wnutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaul2_p4", "TLorentzVector", &mc_truth_h0Wtaul2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunu2_p4", "TLorentzVector", &mc_truth_h0Wtaunu2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunutau2_p4", "TLorentzVector", &mc_truth_h0Wtaunutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq11_p4", "TLorentzVector", &mc_truth_h0Wq11_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq21_p4", "TLorentzVector", &mc_truth_h0Wq21_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq12_p4", "TLorentzVector", &mc_truth_h0Wq12_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq22_p4", "TLorentzVector", &mc_truth_h0Wq22_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq11_IS_p4", "TLorentzVector", &mc_truth_h0Wq11_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq21_IS_p4", "TLorentzVector", &mc_truth_h0Wq21_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq12_IS_p4", "TLorentzVector", &mc_truth_h0Wq12_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq22_IS_p4", "TLorentzVector", &mc_truth_h0Wq22_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0Z1_p4", "TLorentzVector", &mc_truth_h0Z1_p4, buffersize);
	     tree->Branch("mc_truth_h0Z2_p4", "TLorentzVector", &mc_truth_h0Z2_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl11_p4", "TLorentzVector", &mc_truth_h0Zl11_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl21_p4", "TLorentzVector", &mc_truth_h0Zl21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl12_p4", "TLorentzVector", &mc_truth_h0Zl12_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl22_p4", "TLorentzVector", &mc_truth_h0Zl22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau11_p4", "TLorentzVector", &mc_truth_h0Ztau11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau21_p4", "TLorentzVector", &mc_truth_h0Ztau21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul11_p4", "TLorentzVector", &mc_truth_h0Ztaul11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul21_p4", "TLorentzVector", &mc_truth_h0Ztaul21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu11_p4", "TLorentzVector", &mc_truth_h0Ztaunu11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu21_p4", "TLorentzVector", &mc_truth_h0Ztaunu21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau11_p4", "TLorentzVector", &mc_truth_h0Ztaunutau11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau21_p4", "TLorentzVector", &mc_truth_h0Ztaunutau21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq11_p4", "TLorentzVector", &mc_truth_h0Zq11_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq21_p4", "TLorentzVector", &mc_truth_h0Zq21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq12_p4", "TLorentzVector", &mc_truth_h0Zq12_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq22_p4", "TLorentzVector", &mc_truth_h0Zq22_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq11_IS_p4", "TLorentzVector", &mc_truth_h0Zq11_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq21_IS_p4", "TLorentzVector", &mc_truth_h0Zq21_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq12_IS_p4", "TLorentzVector", &mc_truth_h0Zq12_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq22_IS_p4", "TLorentzVector", &mc_truth_h0Zq22_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau12_p4", "TLorentzVector", &mc_truth_h0Ztau12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau22_p4", "TLorentzVector", &mc_truth_h0Ztau22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul12_p4", "TLorentzVector", &mc_truth_h0Ztaul12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul22_p4", "TLorentzVector", &mc_truth_h0Ztaul22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu12_p4", "TLorentzVector", &mc_truth_h0Ztaunu12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu22_p4", "TLorentzVector", &mc_truth_h0Ztaunu22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau12_p4", "TLorentzVector", &mc_truth_h0Ztaunutau12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau22_p4", "TLorentzVector", &mc_truth_h0Ztaunutau22_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu11_p4", "TLorentzVector", &mc_truth_h0Znu11_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu21_p4", "TLorentzVector", &mc_truth_h0Znu21_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu12_p4", "TLorentzVector", &mc_truth_h0Znu12_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu22_p4", "TLorentzVector", &mc_truth_h0Znu22_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0tau1_p4", "TLorentzVector", &mc_truth_h0tau1_p4, buffersize);
	     tree->Branch("mc_truth_h0tau2_p4", "TLorentzVector", &mc_truth_h0tau2_p4, buffersize);
	     tree->Branch("mc_truth_h0taul1_p4", "TLorentzVector", &mc_truth_h0taul1_p4, buffersize);
	     tree->Branch("mc_truth_h0taunutau1_p4", "TLorentzVector", &mc_truth_h0taunutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0taunu1_p4", "TLorentzVector", &mc_truth_h0taunu1_p4, buffersize);
	     tree->Branch("mc_truth_h0taul2_p4", "TLorentzVector", &mc_truth_h0taul2_p4, buffersize);
	     tree->Branch("mc_truth_h0taunutau2_p4", "TLorentzVector", &mc_truth_h0taunutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0taunu2_p4", "TLorentzVector", &mc_truth_h0taunu2_p4, buffersize);
	     
	     tree->Branch("mc_truth_t1_p4", "TLorentzVector", &mc_truth_t1_p4, buffersize);
	     tree->Branch("mc_truth_t2_p4", "TLorentzVector", &mc_truth_t2_p4, buffersize);
	     tree->Branch("mc_truth_tb1_p4", "TLorentzVector", &mc_truth_tb1_p4, buffersize);
	     tree->Branch("mc_truth_tb2_p4", "TLorentzVector", &mc_truth_tb2_p4, buffersize);
	     tree->Branch("mc_truth_tb1_IS_p4", "TLorentzVector", &mc_truth_tb1_IS_p4, buffersize);
	     tree->Branch("mc_truth_tb2_IS_p4", "TLorentzVector", &mc_truth_tb2_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW1_p4", "TLorentzVector", &mc_truth_tW1_p4, buffersize);
	     tree->Branch("mc_truth_tWnu1_p4", "TLorentzVector", &mc_truth_tWnu1_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau1_p4", "TLorentzVector", &mc_truth_tWnutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWl1_p4", "TLorentzVector", &mc_truth_tWl1_p4, buffersize);
	     tree->Branch("mc_truth_tWtau1_p4", "TLorentzVector", &mc_truth_tWtau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu1_p4", "TLorentzVector", &mc_truth_tWtaunu1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau1_p4", "TLorentzVector", &mc_truth_tWtaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul1_p4", "TLorentzVector", &mc_truth_tWtaul1_p4, buffersize);
	     tree->Branch("mc_truth_tWq11_p4", "TLorentzVector", &mc_truth_tWq11_p4, buffersize);
	     tree->Branch("mc_truth_tWq21_p4", "TLorentzVector", &mc_truth_tWq21_p4, buffersize);
	     tree->Branch("mc_truth_tWq11_IS_p4", "TLorentzVector", &mc_truth_tWq11_IS_p4, buffersize);
	     tree->Branch("mc_truth_tWq21_IS_p4", "TLorentzVector", &mc_truth_tWq21_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW2_p4", "TLorentzVector", &mc_truth_tW2_p4, buffersize);
	     tree->Branch("mc_truth_tWnu2_p4", "TLorentzVector", &mc_truth_tWnu2_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau2_p4", "TLorentzVector", &mc_truth_tWnutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWl2_p4", "TLorentzVector", &mc_truth_tWl2_p4, buffersize);
	     tree->Branch("mc_truth_tWtau2_p4", "TLorentzVector", &mc_truth_tWtau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu2_p4", "TLorentzVector", &mc_truth_tWtaunu2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau2_p4", "TLorentzVector", &mc_truth_tWtaunutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul2_p4", "TLorentzVector", &mc_truth_tWtaul2_p4, buffersize);
	     tree->Branch("mc_truth_tWq12_p4", "TLorentzVector", &mc_truth_tWq12_p4, buffersize);
	     tree->Branch("mc_truth_tWq22_p4", "TLorentzVector", &mc_truth_tWq22_p4, buffersize);
	     tree->Branch("mc_truth_tWq12_IS_p4", "TLorentzVector", &mc_truth_tWq12_IS_p4, buffersize);
	     tree->Branch("mc_truth_tWq22_IS_p4", "TLorentzVector", &mc_truth_tWq22_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_j1_p4", "TLorentzVector", &mc_truth_j1_p4, buffersize);
	     tree->Branch("mc_truth_j2_p4", "TLorentzVector", &mc_truth_j2_p4, buffersize);
	     tree->Branch("mc_truth_j3_p4", "TLorentzVector", &mc_truth_j3_p4, buffersize);
	  }

	tree->Branch("mc_truth_h0_pt", &mc_truth_h0_pt, "mc_truth_h0_pt/F", buffersize);

	tree->Branch("mc_truth_h0W1_pt", &mc_truth_h0W1_pt, "mc_truth_h0W1_pt/F", buffersize);
	tree->Branch("mc_truth_h0W2_pt", &mc_truth_h0W2_pt, "mc_truth_h0W2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_pt", &mc_truth_h0Wl1_pt, "mc_truth_h0Wl1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_pt", &mc_truth_h0Wnu1_pt, "mc_truth_h0Wnu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_pt", &mc_truth_h0Wtau1_pt, "mc_truth_h0Wtau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_pt", &mc_truth_h0Wnutau1_pt, "mc_truth_h0Wnutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_pt", &mc_truth_h0Wtaul1_pt, "mc_truth_h0Wtaul1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_pt", &mc_truth_h0Wtaunu1_pt, "mc_truth_h0Wtaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_pt", &mc_truth_h0Wtaunutau1_pt, "mc_truth_h0Wtaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_pt", &mc_truth_h0Wl2_pt, "mc_truth_h0Wl2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_pt", &mc_truth_h0Wnu2_pt, "mc_truth_h0Wnu2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_pt", &mc_truth_h0Wtau2_pt, "mc_truth_h0Wtau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_pt", &mc_truth_h0Wnutau2_pt, "mc_truth_h0Wnutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_pt", &mc_truth_h0Wtaul2_pt, "mc_truth_h0Wtaul2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_pt", &mc_truth_h0Wtaunu2_pt, "mc_truth_h0Wtaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_pt", &mc_truth_h0Wtaunutau2_pt, "mc_truth_h0Wtaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_pt", &mc_truth_h0Wq11_pt, "mc_truth_h0Wq11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_pt", &mc_truth_h0Wq21_pt, "mc_truth_h0Wq21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_pt", &mc_truth_h0Wq12_pt, "mc_truth_h0Wq12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_pt", &mc_truth_h0Wq22_pt, "mc_truth_h0Wq22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_pt", &mc_truth_h0Wq11_IS_pt, "mc_truth_h0Wq11_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_pt", &mc_truth_h0Wq21_IS_pt, "mc_truth_h0Wq21_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_pt", &mc_truth_h0Wq12_IS_pt, "mc_truth_h0Wq12_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_pt", &mc_truth_h0Wq22_IS_pt, "mc_truth_h0Wq22_IS_pt/F", buffersize);

	tree->Branch("mc_truth_h0Z1_pt", &mc_truth_h0Z1_pt, "mc_truth_h0Z1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Z2_pt", &mc_truth_h0Z2_pt, "mc_truth_h0Z2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_pt", &mc_truth_h0Zl11_pt, "mc_truth_h0Zl11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_pt", &mc_truth_h0Zl21_pt, "mc_truth_h0Zl21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_pt", &mc_truth_h0Ztau11_pt, "mc_truth_h0Ztau11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_pt", &mc_truth_h0Ztau21_pt, "mc_truth_h0Ztau21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_pt", &mc_truth_h0Ztaul11_pt, "mc_truth_h0Ztaul11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_pt", &mc_truth_h0Ztaul21_pt, "mc_truth_h0Ztaul21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_pt", &mc_truth_h0Ztaunu11_pt, "mc_truth_h0Ztaunu11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_pt", &mc_truth_h0Ztaunu21_pt, "mc_truth_h0Ztaunu21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_pt", &mc_truth_h0Ztaunutau11_pt, "mc_truth_h0Ztaunutau11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_pt", &mc_truth_h0Ztaunutau21_pt, "mc_truth_h0Ztaunutau21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_pt", &mc_truth_h0Zq11_pt, "mc_truth_h0Zq11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_pt", &mc_truth_h0Zq21_pt, "mc_truth_h0Zq21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_pt", &mc_truth_h0Zq11_IS_pt, "mc_truth_h0Zq11_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_pt", &mc_truth_h0Zq21_IS_pt, "mc_truth_h0Zq21_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_pt", &mc_truth_h0Zl12_pt, "mc_truth_h0Zl12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_pt", &mc_truth_h0Zl22_pt, "mc_truth_h0Zl22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_pt", &mc_truth_h0Ztau12_pt, "mc_truth_h0Ztau12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_pt", &mc_truth_h0Ztau22_pt, "mc_truth_h0Ztau22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_pt", &mc_truth_h0Ztaul12_pt, "mc_truth_h0Ztaul12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_pt", &mc_truth_h0Ztaul22_pt, "mc_truth_h0Ztaul22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_pt", &mc_truth_h0Ztaunu12_pt, "mc_truth_h0Ztaunu12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_pt", &mc_truth_h0Ztaunu22_pt, "mc_truth_h0Ztaunu22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_pt", &mc_truth_h0Ztaunutau12_pt, "mc_truth_h0Ztaunutau12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_pt", &mc_truth_h0Ztaunutau22_pt, "mc_truth_h0Ztaunutau22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_pt", &mc_truth_h0Zq12_pt, "mc_truth_h0Zq12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_pt", &mc_truth_h0Zq22_pt, "mc_truth_h0Zq22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_pt", &mc_truth_h0Zq12_IS_pt, "mc_truth_h0Zq12_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_pt", &mc_truth_h0Zq22_IS_pt, "mc_truth_h0Zq22_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_pt", &mc_truth_h0Znu11_pt, "mc_truth_h0Znu11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_pt", &mc_truth_h0Znu21_pt, "mc_truth_h0Znu21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_pt", &mc_truth_h0Znu12_pt, "mc_truth_h0Znu12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_pt", &mc_truth_h0Znu22_pt, "mc_truth_h0Znu22_pt/F", buffersize);

	tree->Branch("mc_truth_h0tau1_pt", &mc_truth_h0tau1_pt, "mc_truth_h0tau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0tau2_pt", &mc_truth_h0tau2_pt, "mc_truth_h0tau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taul1_pt", &mc_truth_h0taul1_pt, "mc_truth_h0taul1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_pt", &mc_truth_h0taunutau1_pt, "mc_truth_h0taunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_pt", &mc_truth_h0taunu1_pt, "mc_truth_h0taunu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taul2_pt", &mc_truth_h0taul2_pt, "mc_truth_h0taul2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_pt", &mc_truth_h0taunutau2_pt, "mc_truth_h0taunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_pt", &mc_truth_h0taunu2_pt, "mc_truth_h0taunu2_pt/F", buffersize);

	tree->Branch("mc_truth_t1_pt", &mc_truth_t1_pt, "mc_truth_t1_pt/F", buffersize);
	tree->Branch("mc_truth_t2_pt", &mc_truth_t2_pt, "mc_truth_t2_pt/F", buffersize);
	tree->Branch("mc_truth_tb1_pt", &mc_truth_tb1_pt, "mc_truth_tb1_pt/F", buffersize);
	tree->Branch("mc_truth_tb2_pt", &mc_truth_tb2_pt, "mc_truth_tb2_pt/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_pt", &mc_truth_tb1_IS_pt, "mc_truth_tb1_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_pt", &mc_truth_tb2_IS_pt, "mc_truth_tb2_IS_pt/F", buffersize);

	tree->Branch("mc_truth_tW1_pt", &mc_truth_tW1_pt, "mc_truth_tW1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu1_pt", &mc_truth_tWnu1_pt, "mc_truth_tWnu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_pt", &mc_truth_tWnutau1_pt, "mc_truth_tWnutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWl1_pt", &mc_truth_tWl1_pt, "mc_truth_tWl1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau1_pt", &mc_truth_tWtau1_pt, "mc_truth_tWtau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_pt", &mc_truth_tWtaunu1_pt, "mc_truth_tWtaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_pt", &mc_truth_tWtaunutau1_pt, "mc_truth_tWtaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_pt", &mc_truth_tWtaul1_pt, "mc_truth_tWtaul1_pt/F", buffersize);
	tree->Branch("mc_truth_tWq11_pt", &mc_truth_tWq11_pt, "mc_truth_tWq11_pt/F", buffersize);
	tree->Branch("mc_truth_tWq21_pt", &mc_truth_tWq21_pt, "mc_truth_tWq21_pt/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_pt", &mc_truth_tWq11_IS_pt, "mc_truth_tWq11_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_pt", &mc_truth_tWq21_IS_pt, "mc_truth_tWq21_IS_pt/F", buffersize);

	tree->Branch("mc_truth_tW2_pt", &mc_truth_tW2_pt, "mc_truth_tW2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu2_pt", &mc_truth_tWnu2_pt, "mc_truth_tWnu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_pt", &mc_truth_tWnutau2_pt, "mc_truth_tWnutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWl2_pt", &mc_truth_tWl2_pt, "mc_truth_tWl2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau2_pt", &mc_truth_tWtau2_pt, "mc_truth_tWtau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_pt", &mc_truth_tWtaunu2_pt, "mc_truth_tWtaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_pt", &mc_truth_tWtaunutau2_pt, "mc_truth_tWtaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_pt", &mc_truth_tWtaul2_pt, "mc_truth_tWtaul2_pt/F", buffersize);
	tree->Branch("mc_truth_tWq12_pt", &mc_truth_tWq12_pt, "mc_truth_tWq12_pt/F", buffersize);
	tree->Branch("mc_truth_tWq22_pt", &mc_truth_tWq22_pt, "mc_truth_tWq22_pt/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_pt", &mc_truth_tWq12_IS_pt, "mc_truth_tWq12_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_pt", &mc_truth_tWq22_IS_pt, "mc_truth_tWq22_IS_pt/F", buffersize);

	tree->Branch("mc_truth_j1_pt", &mc_truth_j1_pt, "mc_truth_j1_pt/F", buffersize);
	tree->Branch("mc_truth_j2_pt", &mc_truth_j2_pt, "mc_truth_j2_pt/F", buffersize);
	tree->Branch("mc_truth_j3_pt", &mc_truth_j3_pt, "mc_truth_j3_pt/F", buffersize);
		
	tree->Branch("mc_truth_h0_eta", &mc_truth_h0_eta, "mc_truth_h0_eta/F", buffersize);

	tree->Branch("mc_truth_h0W1_eta", &mc_truth_h0W1_eta, "mc_truth_h0W1_eta/F", buffersize);
	tree->Branch("mc_truth_h0W2_eta", &mc_truth_h0W2_eta, "mc_truth_h0W2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_eta", &mc_truth_h0Wl1_eta, "mc_truth_h0Wl1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_eta", &mc_truth_h0Wnu1_eta, "mc_truth_h0Wnu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_eta", &mc_truth_h0Wtau1_eta, "mc_truth_h0Wtau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_eta", &mc_truth_h0Wnutau1_eta, "mc_truth_h0Wnutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_eta", &mc_truth_h0Wtaul1_eta, "mc_truth_h0Wtaul1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_eta", &mc_truth_h0Wtaunu1_eta, "mc_truth_h0Wtaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_eta", &mc_truth_h0Wtaunutau1_eta, "mc_truth_h0Wtaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_eta", &mc_truth_h0Wl2_eta, "mc_truth_h0Wl2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_eta", &mc_truth_h0Wnu2_eta, "mc_truth_h0Wnu2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_eta", &mc_truth_h0Wtau2_eta, "mc_truth_h0Wtau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_eta", &mc_truth_h0Wnutau2_eta, "mc_truth_h0Wnutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_eta", &mc_truth_h0Wtaul2_eta, "mc_truth_h0Wtaul2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_eta", &mc_truth_h0Wtaunu2_eta, "mc_truth_h0Wtaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_eta", &mc_truth_h0Wtaunutau2_eta, "mc_truth_h0Wtaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_eta", &mc_truth_h0Wq11_eta, "mc_truth_h0Wq11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_eta", &mc_truth_h0Wq21_eta, "mc_truth_h0Wq21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_eta", &mc_truth_h0Wq12_eta, "mc_truth_h0Wq12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_eta", &mc_truth_h0Wq22_eta, "mc_truth_h0Wq22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_eta", &mc_truth_h0Wq11_IS_eta, "mc_truth_h0Wq11_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_eta", &mc_truth_h0Wq21_IS_eta, "mc_truth_h0Wq21_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_eta", &mc_truth_h0Wq12_IS_eta, "mc_truth_h0Wq12_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_eta", &mc_truth_h0Wq22_IS_eta, "mc_truth_h0Wq22_IS_eta/F", buffersize);

	tree->Branch("mc_truth_h0Z1_eta", &mc_truth_h0Z1_eta, "mc_truth_h0Z1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Z2_eta", &mc_truth_h0Z2_eta, "mc_truth_h0Z2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_eta", &mc_truth_h0Zl11_eta, "mc_truth_h0Zl11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_eta", &mc_truth_h0Zl21_eta, "mc_truth_h0Zl21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_eta", &mc_truth_h0Ztau11_eta, "mc_truth_h0Ztau11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_eta", &mc_truth_h0Ztau21_eta, "mc_truth_h0Ztau21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_eta", &mc_truth_h0Ztaul11_eta, "mc_truth_h0Ztaul11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_eta", &mc_truth_h0Ztaul21_eta, "mc_truth_h0Ztaul21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_eta", &mc_truth_h0Ztaunu11_eta, "mc_truth_h0Ztaunu11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_eta", &mc_truth_h0Ztaunu21_eta, "mc_truth_h0Ztaunu21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_eta", &mc_truth_h0Ztaunutau11_eta, "mc_truth_h0Ztaunutau11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_eta", &mc_truth_h0Ztaunutau21_eta, "mc_truth_h0Ztaunutau21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_eta", &mc_truth_h0Zq11_eta, "mc_truth_h0Zq11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_eta", &mc_truth_h0Zq21_eta, "mc_truth_h0Zq21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_eta", &mc_truth_h0Zq11_IS_eta, "mc_truth_h0Zq11_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_eta", &mc_truth_h0Zq21_IS_eta, "mc_truth_h0Zq21_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_eta", &mc_truth_h0Zl12_eta, "mc_truth_h0Zl12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_eta", &mc_truth_h0Zl22_eta, "mc_truth_h0Zl22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_eta", &mc_truth_h0Ztau12_eta, "mc_truth_h0Ztau12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_eta", &mc_truth_h0Ztau22_eta, "mc_truth_h0Ztau22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_eta", &mc_truth_h0Ztaul12_eta, "mc_truth_h0Ztaul12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_eta", &mc_truth_h0Ztaul22_eta, "mc_truth_h0Ztaul22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_eta", &mc_truth_h0Ztaunu12_eta, "mc_truth_h0Ztaunu12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_eta", &mc_truth_h0Ztaunu22_eta, "mc_truth_h0Ztaunu22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_eta", &mc_truth_h0Ztaunutau12_eta, "mc_truth_h0Ztaunutau12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_eta", &mc_truth_h0Ztaunutau22_eta, "mc_truth_h0Ztaunutau22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_eta", &mc_truth_h0Zq12_eta, "mc_truth_h0Zq12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_eta", &mc_truth_h0Zq22_eta, "mc_truth_h0Zq22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_eta", &mc_truth_h0Zq12_IS_eta, "mc_truth_h0Zq12_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_eta", &mc_truth_h0Zq22_IS_eta, "mc_truth_h0Zq22_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_eta", &mc_truth_h0Znu11_eta, "mc_truth_h0Znu11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_eta", &mc_truth_h0Znu21_eta, "mc_truth_h0Znu21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_eta", &mc_truth_h0Znu12_eta, "mc_truth_h0Znu12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_eta", &mc_truth_h0Znu22_eta, "mc_truth_h0Znu22_eta/F", buffersize);

	tree->Branch("mc_truth_h0tau1_eta", &mc_truth_h0tau1_eta, "mc_truth_h0tau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0tau2_eta", &mc_truth_h0tau2_eta, "mc_truth_h0tau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taul1_eta", &mc_truth_h0taul1_eta, "mc_truth_h0taul1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_eta", &mc_truth_h0taunutau1_eta, "mc_truth_h0taunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_eta", &mc_truth_h0taunu1_eta, "mc_truth_h0taunu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taul2_eta", &mc_truth_h0taul2_eta, "mc_truth_h0taul2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_eta", &mc_truth_h0taunutau2_eta, "mc_truth_h0taunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_eta", &mc_truth_h0taunu2_eta, "mc_truth_h0taunu2_eta/F", buffersize);

	tree->Branch("mc_truth_t1_eta", &mc_truth_t1_eta, "mc_truth_t1_eta/F", buffersize);
	tree->Branch("mc_truth_t2_eta", &mc_truth_t2_eta, "mc_truth_t2_eta/F", buffersize);
	tree->Branch("mc_truth_tb1_eta", &mc_truth_tb1_eta, "mc_truth_tb1_eta/F", buffersize);
	tree->Branch("mc_truth_tb2_eta", &mc_truth_tb2_eta, "mc_truth_tb2_eta/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_eta", &mc_truth_tb1_IS_eta, "mc_truth_tb1_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_eta", &mc_truth_tb2_IS_eta, "mc_truth_tb2_IS_eta/F", buffersize);

	tree->Branch("mc_truth_tW1_eta", &mc_truth_tW1_eta, "mc_truth_tW1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu1_eta", &mc_truth_tWnu1_eta, "mc_truth_tWnu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_eta", &mc_truth_tWnutau1_eta, "mc_truth_tWnutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWl1_eta", &mc_truth_tWl1_eta, "mc_truth_tWl1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau1_eta", &mc_truth_tWtau1_eta, "mc_truth_tWtau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_eta", &mc_truth_tWtaunu1_eta, "mc_truth_tWtaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_eta", &mc_truth_tWtaunutau1_eta, "mc_truth_tWtaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_eta", &mc_truth_tWtaul1_eta, "mc_truth_tWtaul1_eta/F", buffersize);
	tree->Branch("mc_truth_tWq11_eta", &mc_truth_tWq11_eta, "mc_truth_tWq11_eta/F", buffersize);
	tree->Branch("mc_truth_tWq21_eta", &mc_truth_tWq21_eta, "mc_truth_tWq21_eta/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_eta", &mc_truth_tWq11_IS_eta, "mc_truth_tWq11_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_eta", &mc_truth_tWq21_IS_eta, "mc_truth_tWq21_IS_eta/F", buffersize);

	tree->Branch("mc_truth_tW2_eta", &mc_truth_tW2_eta, "mc_truth_tW2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu2_eta", &mc_truth_tWnu2_eta, "mc_truth_tWnu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_eta", &mc_truth_tWnutau2_eta, "mc_truth_tWnutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWl2_eta", &mc_truth_tWl2_eta, "mc_truth_tWl2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau2_eta", &mc_truth_tWtau2_eta, "mc_truth_tWtau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_eta", &mc_truth_tWtaunu2_eta, "mc_truth_tWtaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_eta", &mc_truth_tWtaunutau2_eta, "mc_truth_tWtaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_eta", &mc_truth_tWtaul2_eta, "mc_truth_tWtaul2_eta/F", buffersize);
	tree->Branch("mc_truth_tWq12_eta", &mc_truth_tWq12_eta, "mc_truth_tWq12_eta/F", buffersize);
	tree->Branch("mc_truth_tWq22_eta", &mc_truth_tWq22_eta, "mc_truth_tWq22_eta/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_eta", &mc_truth_tWq12_IS_eta, "mc_truth_tWq12_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_eta", &mc_truth_tWq22_IS_eta, "mc_truth_tWq22_IS_eta/F", buffersize);

	tree->Branch("mc_truth_j1_eta", &mc_truth_j1_eta, "mc_truth_j1_eta/F", buffersize);
	tree->Branch("mc_truth_j2_eta", &mc_truth_j2_eta, "mc_truth_j2_eta/F", buffersize);
	tree->Branch("mc_truth_j3_eta", &mc_truth_j3_eta, "mc_truth_j3_eta/F", buffersize);
		
	tree->Branch("mc_truth_h0_phi", &mc_truth_h0_phi, "mc_truth_h0_phi/F", buffersize);

	tree->Branch("mc_truth_h0W1_phi", &mc_truth_h0W1_phi, "mc_truth_h0W1_phi/F", buffersize);
	tree->Branch("mc_truth_h0W2_phi", &mc_truth_h0W2_phi, "mc_truth_h0W2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_phi", &mc_truth_h0Wl1_phi, "mc_truth_h0Wl1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_phi", &mc_truth_h0Wnu1_phi, "mc_truth_h0Wnu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_phi", &mc_truth_h0Wtau1_phi, "mc_truth_h0Wtau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_phi", &mc_truth_h0Wnutau1_phi, "mc_truth_h0Wnutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_phi", &mc_truth_h0Wtaul1_phi, "mc_truth_h0Wtaul1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_phi", &mc_truth_h0Wtaunu1_phi, "mc_truth_h0Wtaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_phi", &mc_truth_h0Wtaunutau1_phi, "mc_truth_h0Wtaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_phi", &mc_truth_h0Wl2_phi, "mc_truth_h0Wl2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_phi", &mc_truth_h0Wnu2_phi, "mc_truth_h0Wnu2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_phi", &mc_truth_h0Wtau2_phi, "mc_truth_h0Wtau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_phi", &mc_truth_h0Wnutau2_phi, "mc_truth_h0Wnutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_phi", &mc_truth_h0Wtaul2_phi, "mc_truth_h0Wtaul2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_phi", &mc_truth_h0Wtaunu2_phi, "mc_truth_h0Wtaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_phi", &mc_truth_h0Wtaunutau2_phi, "mc_truth_h0Wtaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_phi", &mc_truth_h0Wq11_phi, "mc_truth_h0Wq11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_phi", &mc_truth_h0Wq21_phi, "mc_truth_h0Wq21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_phi", &mc_truth_h0Wq12_phi, "mc_truth_h0Wq12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_phi", &mc_truth_h0Wq22_phi, "mc_truth_h0Wq22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_phi", &mc_truth_h0Wq11_IS_phi, "mc_truth_h0Wq11_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_phi", &mc_truth_h0Wq21_IS_phi, "mc_truth_h0Wq21_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_phi", &mc_truth_h0Wq12_IS_phi, "mc_truth_h0Wq12_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_phi", &mc_truth_h0Wq22_IS_phi, "mc_truth_h0Wq22_IS_phi/F", buffersize);

	tree->Branch("mc_truth_h0Z1_phi", &mc_truth_h0Z1_phi, "mc_truth_h0Z1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Z2_phi", &mc_truth_h0Z2_phi, "mc_truth_h0Z2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_phi", &mc_truth_h0Zl11_phi, "mc_truth_h0Zl11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_phi", &mc_truth_h0Zl21_phi, "mc_truth_h0Zl21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_phi", &mc_truth_h0Ztau11_phi, "mc_truth_h0Ztau11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_phi", &mc_truth_h0Ztau21_phi, "mc_truth_h0Ztau21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_phi", &mc_truth_h0Ztaul11_phi, "mc_truth_h0Ztaul11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_phi", &mc_truth_h0Ztaul21_phi, "mc_truth_h0Ztaul21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_phi", &mc_truth_h0Ztaunu11_phi, "mc_truth_h0Ztaunu11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_phi", &mc_truth_h0Ztaunu21_phi, "mc_truth_h0Ztaunu21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_phi", &mc_truth_h0Ztaunutau11_phi, "mc_truth_h0Ztaunutau11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_phi", &mc_truth_h0Ztaunutau21_phi, "mc_truth_h0Ztaunutau21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_phi", &mc_truth_h0Zq11_phi, "mc_truth_h0Zq11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_phi", &mc_truth_h0Zq21_phi, "mc_truth_h0Zq21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_phi", &mc_truth_h0Zq11_IS_phi, "mc_truth_h0Zq11_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_phi", &mc_truth_h0Zq21_IS_phi, "mc_truth_h0Zq21_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_phi", &mc_truth_h0Zl12_phi, "mc_truth_h0Zl12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_phi", &mc_truth_h0Zl22_phi, "mc_truth_h0Zl22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_phi", &mc_truth_h0Ztau12_phi, "mc_truth_h0Ztau12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_phi", &mc_truth_h0Ztau22_phi, "mc_truth_h0Ztau22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_phi", &mc_truth_h0Ztaul12_phi, "mc_truth_h0Ztaul12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_phi", &mc_truth_h0Ztaul22_phi, "mc_truth_h0Ztaul22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_phi", &mc_truth_h0Ztaunu12_phi, "mc_truth_h0Ztaunu12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_phi", &mc_truth_h0Ztaunu22_phi, "mc_truth_h0Ztaunu22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_phi", &mc_truth_h0Ztaunutau12_phi, "mc_truth_h0Ztaunutau12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_phi", &mc_truth_h0Ztaunutau22_phi, "mc_truth_h0Ztaunutau22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_phi", &mc_truth_h0Zq12_phi, "mc_truth_h0Zq12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_phi", &mc_truth_h0Zq22_phi, "mc_truth_h0Zq22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_phi", &mc_truth_h0Zq12_IS_phi, "mc_truth_h0Zq12_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_phi", &mc_truth_h0Zq22_IS_phi, "mc_truth_h0Zq22_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_phi", &mc_truth_h0Znu11_phi, "mc_truth_h0Znu11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_phi", &mc_truth_h0Znu21_phi, "mc_truth_h0Znu21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_phi", &mc_truth_h0Znu12_phi, "mc_truth_h0Znu12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_phi", &mc_truth_h0Znu22_phi, "mc_truth_h0Znu22_phi/F", buffersize);

	tree->Branch("mc_truth_h0tau1_phi", &mc_truth_h0tau1_phi, "mc_truth_h0tau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0tau2_phi", &mc_truth_h0tau2_phi, "mc_truth_h0tau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taul1_phi", &mc_truth_h0taul1_phi, "mc_truth_h0taul1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_phi", &mc_truth_h0taunutau1_phi, "mc_truth_h0taunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_phi", &mc_truth_h0taunu1_phi, "mc_truth_h0taunu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taul2_phi", &mc_truth_h0taul2_phi, "mc_truth_h0taul2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_phi", &mc_truth_h0taunutau2_phi, "mc_truth_h0taunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_phi", &mc_truth_h0taunu2_phi, "mc_truth_h0taunu2_phi/F", buffersize);

	tree->Branch("mc_truth_t1_phi", &mc_truth_t1_phi, "mc_truth_t1_phi/F", buffersize);
	tree->Branch("mc_truth_t2_phi", &mc_truth_t2_phi, "mc_truth_t2_phi/F", buffersize);
	tree->Branch("mc_truth_tb1_phi", &mc_truth_tb1_phi, "mc_truth_tb1_phi/F", buffersize);
	tree->Branch("mc_truth_tb2_phi", &mc_truth_tb2_phi, "mc_truth_tb2_phi/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_phi", &mc_truth_tb1_IS_phi, "mc_truth_tb1_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_phi", &mc_truth_tb2_IS_phi, "mc_truth_tb2_IS_phi/F", buffersize);

	tree->Branch("mc_truth_tW1_phi", &mc_truth_tW1_phi, "mc_truth_tW1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu1_phi", &mc_truth_tWnu1_phi, "mc_truth_tWnu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_phi", &mc_truth_tWnutau1_phi, "mc_truth_tWnutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWl1_phi", &mc_truth_tWl1_phi, "mc_truth_tWl1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau1_phi", &mc_truth_tWtau1_phi, "mc_truth_tWtau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_phi", &mc_truth_tWtaunu1_phi, "mc_truth_tWtaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_phi", &mc_truth_tWtaunutau1_phi, "mc_truth_tWtaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_phi", &mc_truth_tWtaul1_phi, "mc_truth_tWtaul1_phi/F", buffersize);
	tree->Branch("mc_truth_tWq11_phi", &mc_truth_tWq11_phi, "mc_truth_tWq11_phi/F", buffersize);
	tree->Branch("mc_truth_tWq21_phi", &mc_truth_tWq21_phi, "mc_truth_tWq21_phi/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_phi", &mc_truth_tWq11_IS_phi, "mc_truth_tWq11_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_phi", &mc_truth_tWq21_IS_phi, "mc_truth_tWq21_IS_phi/F", buffersize);

	tree->Branch("mc_truth_tW2_phi", &mc_truth_tW2_phi, "mc_truth_tW2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu2_phi", &mc_truth_tWnu2_phi, "mc_truth_tWnu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_phi", &mc_truth_tWnutau2_phi, "mc_truth_tWnutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWl2_phi", &mc_truth_tWl2_phi, "mc_truth_tWl2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau2_phi", &mc_truth_tWtau2_phi, "mc_truth_tWtau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_phi", &mc_truth_tWtaunu2_phi, "mc_truth_tWtaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_phi", &mc_truth_tWtaunutau2_phi, "mc_truth_tWtaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_phi", &mc_truth_tWtaul2_phi, "mc_truth_tWtaul2_phi/F", buffersize);
	tree->Branch("mc_truth_tWq12_phi", &mc_truth_tWq12_phi, "mc_truth_tWq12_phi/F", buffersize);
	tree->Branch("mc_truth_tWq22_phi", &mc_truth_tWq22_phi, "mc_truth_tWq22_phi/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_phi", &mc_truth_tWq12_IS_phi, "mc_truth_tWq12_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_phi", &mc_truth_tWq22_IS_phi, "mc_truth_tWq22_IS_phi/F", buffersize);

	tree->Branch("mc_truth_j1_phi", &mc_truth_j1_phi, "mc_truth_j1_phi/F", buffersize);
	tree->Branch("mc_truth_j2_phi", &mc_truth_j2_phi, "mc_truth_j2_phi/F", buffersize);
	tree->Branch("mc_truth_j3_phi", &mc_truth_j3_phi, "mc_truth_j3_phi/F", buffersize);
	
	tree->Branch("mc_truth_h0_E", &mc_truth_h0_E, "mc_truth_h0_E/F", buffersize);

	tree->Branch("mc_truth_h0W1_E", &mc_truth_h0W1_E, "mc_truth_h0W1_E/F", buffersize);
	tree->Branch("mc_truth_h0W2_E", &mc_truth_h0W2_E, "mc_truth_h0W2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_E", &mc_truth_h0Wl1_E, "mc_truth_h0Wl1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_E", &mc_truth_h0Wnu1_E, "mc_truth_h0Wnu1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_E", &mc_truth_h0Wtau1_E, "mc_truth_h0Wtau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_E", &mc_truth_h0Wnutau1_E, "mc_truth_h0Wnutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_E", &mc_truth_h0Wtaul1_E, "mc_truth_h0Wtaul1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_E", &mc_truth_h0Wtaunu1_E, "mc_truth_h0Wtaunu1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_E", &mc_truth_h0Wtaunutau1_E, "mc_truth_h0Wtaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_E", &mc_truth_h0Wl2_E, "mc_truth_h0Wl2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_E", &mc_truth_h0Wnu2_E, "mc_truth_h0Wnu2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_E", &mc_truth_h0Wtau2_E, "mc_truth_h0Wtau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_E", &mc_truth_h0Wnutau2_E, "mc_truth_h0Wnutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_E", &mc_truth_h0Wtaul2_E, "mc_truth_h0Wtaul2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_E", &mc_truth_h0Wtaunu2_E, "mc_truth_h0Wtaunu2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_E", &mc_truth_h0Wtaunutau2_E, "mc_truth_h0Wtaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_E", &mc_truth_h0Wq11_E, "mc_truth_h0Wq11_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_E", &mc_truth_h0Wq21_E, "mc_truth_h0Wq21_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_E", &mc_truth_h0Wq12_E, "mc_truth_h0Wq12_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_E", &mc_truth_h0Wq22_E, "mc_truth_h0Wq22_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_E", &mc_truth_h0Wq11_IS_E, "mc_truth_h0Wq11_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_E", &mc_truth_h0Wq21_IS_E, "mc_truth_h0Wq21_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_E", &mc_truth_h0Wq12_IS_E, "mc_truth_h0Wq12_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_E", &mc_truth_h0Wq22_IS_E, "mc_truth_h0Wq22_IS_E/F", buffersize);

	tree->Branch("mc_truth_h0Z1_E", &mc_truth_h0Z1_E, "mc_truth_h0Z1_E/F", buffersize);
	tree->Branch("mc_truth_h0Z2_E", &mc_truth_h0Z2_E, "mc_truth_h0Z2_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_E", &mc_truth_h0Zl11_E, "mc_truth_h0Zl11_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_E", &mc_truth_h0Zl21_E, "mc_truth_h0Zl21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_E", &mc_truth_h0Ztau11_E, "mc_truth_h0Ztau11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_E", &mc_truth_h0Ztau21_E, "mc_truth_h0Ztau21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_E", &mc_truth_h0Ztaul11_E, "mc_truth_h0Ztaul11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_E", &mc_truth_h0Ztaul21_E, "mc_truth_h0Ztaul21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_E", &mc_truth_h0Ztaunu11_E, "mc_truth_h0Ztaunu11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_E", &mc_truth_h0Ztaunu21_E, "mc_truth_h0Ztaunu21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_E", &mc_truth_h0Ztaunutau11_E, "mc_truth_h0Ztaunutau11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_E", &mc_truth_h0Ztaunutau21_E, "mc_truth_h0Ztaunutau21_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_E", &mc_truth_h0Zq11_E, "mc_truth_h0Zq11_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_E", &mc_truth_h0Zq21_E, "mc_truth_h0Zq21_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_E", &mc_truth_h0Zq11_IS_E, "mc_truth_h0Zq11_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_E", &mc_truth_h0Zq21_IS_E, "mc_truth_h0Zq21_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_E", &mc_truth_h0Zl12_E, "mc_truth_h0Zl12_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_E", &mc_truth_h0Zl22_E, "mc_truth_h0Zl22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_E", &mc_truth_h0Ztau12_E, "mc_truth_h0Ztau12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_E", &mc_truth_h0Ztau22_E, "mc_truth_h0Ztau22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_E", &mc_truth_h0Ztaul12_E, "mc_truth_h0Ztaul12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_E", &mc_truth_h0Ztaul22_E, "mc_truth_h0Ztaul22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_E", &mc_truth_h0Ztaunu12_E, "mc_truth_h0Ztaunu12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_E", &mc_truth_h0Ztaunu22_E, "mc_truth_h0Ztaunu22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_E", &mc_truth_h0Ztaunutau12_E, "mc_truth_h0Ztaunutau12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_E", &mc_truth_h0Ztaunutau22_E, "mc_truth_h0Ztaunutau22_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_E", &mc_truth_h0Zq12_E, "mc_truth_h0Zq12_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_E", &mc_truth_h0Zq22_E, "mc_truth_h0Zq22_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_E", &mc_truth_h0Zq12_IS_E, "mc_truth_h0Zq12_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_E", &mc_truth_h0Zq22_IS_E, "mc_truth_h0Zq22_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_E", &mc_truth_h0Znu11_E, "mc_truth_h0Znu11_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_E", &mc_truth_h0Znu21_E, "mc_truth_h0Znu21_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_E", &mc_truth_h0Znu12_E, "mc_truth_h0Znu12_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_E", &mc_truth_h0Znu22_E, "mc_truth_h0Znu22_E/F", buffersize);

	tree->Branch("mc_truth_h0tau1_E", &mc_truth_h0tau1_E, "mc_truth_h0tau1_E/F", buffersize);
	tree->Branch("mc_truth_h0tau2_E", &mc_truth_h0tau2_E, "mc_truth_h0tau2_E/F", buffersize);
	tree->Branch("mc_truth_h0taul1_E", &mc_truth_h0taul1_E, "mc_truth_h0taul1_E/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_E", &mc_truth_h0taunutau1_E, "mc_truth_h0taunutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_E", &mc_truth_h0taunu1_E, "mc_truth_h0taunu1_E/F", buffersize);
	tree->Branch("mc_truth_h0taul2_E", &mc_truth_h0taul2_E, "mc_truth_h0taul2_E/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_E", &mc_truth_h0taunutau2_E, "mc_truth_h0taunutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_E", &mc_truth_h0taunu2_E, "mc_truth_h0taunu2_E/F", buffersize);

	tree->Branch("mc_truth_t1_E", &mc_truth_t1_E, "mc_truth_t1_E/F", buffersize);
	tree->Branch("mc_truth_t2_E", &mc_truth_t2_E, "mc_truth_t2_E/F", buffersize);
	tree->Branch("mc_truth_tb1_E", &mc_truth_tb1_E, "mc_truth_tb1_E/F", buffersize);
	tree->Branch("mc_truth_tb2_E", &mc_truth_tb2_E, "mc_truth_tb2_E/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_E", &mc_truth_tb1_IS_E, "mc_truth_tb1_IS_E/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_E", &mc_truth_tb2_IS_E, "mc_truth_tb2_IS_E/F", buffersize);

	tree->Branch("mc_truth_tW1_E", &mc_truth_tW1_E, "mc_truth_tW1_E/F", buffersize);
	tree->Branch("mc_truth_tWnu1_E", &mc_truth_tWnu1_E, "mc_truth_tWnu1_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_E", &mc_truth_tWnutau1_E, "mc_truth_tWnutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWl1_E", &mc_truth_tWl1_E, "mc_truth_tWl1_E/F", buffersize);
	tree->Branch("mc_truth_tWtau1_E", &mc_truth_tWtau1_E, "mc_truth_tWtau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_E", &mc_truth_tWtaunu1_E, "mc_truth_tWtaunu1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_E", &mc_truth_tWtaunutau1_E, "mc_truth_tWtaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_E", &mc_truth_tWtaul1_E, "mc_truth_tWtaul1_E/F", buffersize);
	tree->Branch("mc_truth_tWq11_E", &mc_truth_tWq11_E, "mc_truth_tWq11_E/F", buffersize);
	tree->Branch("mc_truth_tWq21_E", &mc_truth_tWq21_E, "mc_truth_tWq21_E/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_E", &mc_truth_tWq11_IS_E, "mc_truth_tWq11_IS_E/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_E", &mc_truth_tWq21_IS_E, "mc_truth_tWq21_IS_E/F", buffersize);

	tree->Branch("mc_truth_tW2_E", &mc_truth_tW2_E, "mc_truth_tW2_E/F", buffersize);
	tree->Branch("mc_truth_tWnu2_E", &mc_truth_tWnu2_E, "mc_truth_tWnu2_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_E", &mc_truth_tWnutau2_E, "mc_truth_tWnutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWl2_E", &mc_truth_tWl2_E, "mc_truth_tWl2_E/F", buffersize);
	tree->Branch("mc_truth_tWtau2_E", &mc_truth_tWtau2_E, "mc_truth_tWtau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_E", &mc_truth_tWtaunu2_E, "mc_truth_tWtaunu2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_E", &mc_truth_tWtaunutau2_E, "mc_truth_tWtaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_E", &mc_truth_tWtaul2_E, "mc_truth_tWtaul2_E/F", buffersize);
	tree->Branch("mc_truth_tWq12_E", &mc_truth_tWq12_E, "mc_truth_tWq12_E/F", buffersize);
	tree->Branch("mc_truth_tWq22_E", &mc_truth_tWq22_E, "mc_truth_tWq22_E/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_E", &mc_truth_tWq12_IS_E, "mc_truth_tWq12_IS_E/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_E", &mc_truth_tWq22_IS_E, "mc_truth_tWq22_IS_E/F", buffersize);

	tree->Branch("mc_truth_j1_E", &mc_truth_j1_E, "mc_truth_j1_E/F", buffersize);
	tree->Branch("mc_truth_j2_E", &mc_truth_j2_E, "mc_truth_j2_E/F", buffersize);
	tree->Branch("mc_truth_j3_E", &mc_truth_j3_E, "mc_truth_j3_E/F", buffersize);
		
	tree->Branch("mc_truth_h0_id", &mc_truth_h0_id, "mc_truth_h0_id/I", buffersize);

	tree->Branch("mc_truth_h0W1_id", &mc_truth_h0W1_id, "mc_truth_h0W1_id/I", buffersize);
	tree->Branch("mc_truth_h0W2_id", &mc_truth_h0W2_id, "mc_truth_h0W2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wl1_id", &mc_truth_h0Wl1_id, "mc_truth_h0Wl1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnu1_id", &mc_truth_h0Wnu1_id, "mc_truth_h0Wnu1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtau1_id", &mc_truth_h0Wtau1_id, "mc_truth_h0Wtau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_id", &mc_truth_h0Wnutau1_id, "mc_truth_h0Wnutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_id", &mc_truth_h0Wtaul1_id, "mc_truth_h0Wtaul1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_id", &mc_truth_h0Wtaunu1_id, "mc_truth_h0Wtaunu1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_id", &mc_truth_h0Wtaunutau1_id, "mc_truth_h0Wtaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wl2_id", &mc_truth_h0Wl2_id, "mc_truth_h0Wl2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnu2_id", &mc_truth_h0Wnu2_id, "mc_truth_h0Wnu2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtau2_id", &mc_truth_h0Wtau2_id, "mc_truth_h0Wtau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_id", &mc_truth_h0Wnutau2_id, "mc_truth_h0Wnutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_id", &mc_truth_h0Wtaul2_id, "mc_truth_h0Wtaul2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_id", &mc_truth_h0Wtaunu2_id, "mc_truth_h0Wtaunu2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_id", &mc_truth_h0Wtaunutau2_id, "mc_truth_h0Wtaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_id", &mc_truth_h0Wq11_id, "mc_truth_h0Wq11_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_id", &mc_truth_h0Wq21_id, "mc_truth_h0Wq21_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_id", &mc_truth_h0Wq12_id, "mc_truth_h0Wq12_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_id", &mc_truth_h0Wq22_id, "mc_truth_h0Wq22_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_id", &mc_truth_h0Wq11_IS_id, "mc_truth_h0Wq11_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_id", &mc_truth_h0Wq21_IS_id, "mc_truth_h0Wq21_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_id", &mc_truth_h0Wq12_IS_id, "mc_truth_h0Wq12_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_id", &mc_truth_h0Wq22_IS_id, "mc_truth_h0Wq22_IS_id/I", buffersize);

	tree->Branch("mc_truth_h0Z1_id", &mc_truth_h0Z1_id, "mc_truth_h0Z1_id/I", buffersize);
	tree->Branch("mc_truth_h0Z2_id", &mc_truth_h0Z2_id, "mc_truth_h0Z2_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl11_id", &mc_truth_h0Zl11_id, "mc_truth_h0Zl11_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl21_id", &mc_truth_h0Zl21_id, "mc_truth_h0Zl21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau11_id", &mc_truth_h0Ztau11_id, "mc_truth_h0Ztau11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau21_id", &mc_truth_h0Ztau21_id, "mc_truth_h0Ztau21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_id", &mc_truth_h0Ztaul11_id, "mc_truth_h0Ztaul11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_id", &mc_truth_h0Ztaul21_id, "mc_truth_h0Ztaul21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_id", &mc_truth_h0Ztaunu11_id, "mc_truth_h0Ztaunu11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_id", &mc_truth_h0Ztaunu21_id, "mc_truth_h0Ztaunu21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_id", &mc_truth_h0Ztaunutau11_id, "mc_truth_h0Ztaunutau11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_id", &mc_truth_h0Ztaunutau21_id, "mc_truth_h0Ztaunutau21_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_id", &mc_truth_h0Zq11_id, "mc_truth_h0Zq11_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_id", &mc_truth_h0Zq21_id, "mc_truth_h0Zq21_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_id", &mc_truth_h0Zq11_IS_id, "mc_truth_h0Zq11_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_id", &mc_truth_h0Zq21_IS_id, "mc_truth_h0Zq21_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl12_id", &mc_truth_h0Zl12_id, "mc_truth_h0Zl12_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl22_id", &mc_truth_h0Zl22_id, "mc_truth_h0Zl22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau12_id", &mc_truth_h0Ztau12_id, "mc_truth_h0Ztau12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau22_id", &mc_truth_h0Ztau22_id, "mc_truth_h0Ztau22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_id", &mc_truth_h0Ztaul12_id, "mc_truth_h0Ztaul12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_id", &mc_truth_h0Ztaul22_id, "mc_truth_h0Ztaul22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_id", &mc_truth_h0Ztaunu12_id, "mc_truth_h0Ztaunu12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_id", &mc_truth_h0Ztaunu22_id, "mc_truth_h0Ztaunu22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_id", &mc_truth_h0Ztaunutau12_id, "mc_truth_h0Ztaunutau12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_id", &mc_truth_h0Ztaunutau22_id, "mc_truth_h0Ztaunutau22_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_id", &mc_truth_h0Zq12_id, "mc_truth_h0Zq12_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_id", &mc_truth_h0Zq22_id, "mc_truth_h0Zq22_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_id", &mc_truth_h0Zq12_IS_id, "mc_truth_h0Zq12_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_id", &mc_truth_h0Zq22_IS_id, "mc_truth_h0Zq22_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu11_id", &mc_truth_h0Znu11_id, "mc_truth_h0Znu11_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu21_id", &mc_truth_h0Znu21_id, "mc_truth_h0Znu21_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu12_id", &mc_truth_h0Znu12_id, "mc_truth_h0Znu12_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu22_id", &mc_truth_h0Znu22_id, "mc_truth_h0Znu22_id/I", buffersize);

	tree->Branch("mc_truth_h0tau1_id", &mc_truth_h0tau1_id, "mc_truth_h0tau1_id/I", buffersize);
	tree->Branch("mc_truth_h0tau2_id", &mc_truth_h0tau2_id, "mc_truth_h0tau2_id/I", buffersize);
	tree->Branch("mc_truth_h0taul1_id", &mc_truth_h0taul1_id, "mc_truth_h0taul1_id/I", buffersize);
	tree->Branch("mc_truth_h0taunutau1_id", &mc_truth_h0taunutau1_id, "mc_truth_h0taunutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0taunu1_id", &mc_truth_h0taunu1_id, "mc_truth_h0taunu1_id/I", buffersize);
	tree->Branch("mc_truth_h0taul2_id", &mc_truth_h0taul2_id, "mc_truth_h0taul2_id/I", buffersize);
	tree->Branch("mc_truth_h0taunutau2_id", &mc_truth_h0taunutau2_id, "mc_truth_h0taunutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0taunu2_id", &mc_truth_h0taunu2_id, "mc_truth_h0taunu2_id/I", buffersize);

	tree->Branch("mc_truth_t1_id", &mc_truth_t1_id, "mc_truth_t1_id/I", buffersize);
	tree->Branch("mc_truth_t2_id", &mc_truth_t2_id, "mc_truth_t2_id/I", buffersize);
	tree->Branch("mc_truth_tb1_id", &mc_truth_tb1_id, "mc_truth_tb1_id/I", buffersize);
	tree->Branch("mc_truth_tb2_id", &mc_truth_tb2_id, "mc_truth_tb2_id/I", buffersize);
	tree->Branch("mc_truth_tb1_IS_id", &mc_truth_tb1_IS_id, "mc_truth_tb1_IS_id/I", buffersize);
	tree->Branch("mc_truth_tb2_IS_id", &mc_truth_tb2_IS_id, "mc_truth_tb2_IS_id/I", buffersize);

	tree->Branch("mc_truth_tW1_id", &mc_truth_tW1_id, "mc_truth_tW1_id/I", buffersize);
	tree->Branch("mc_truth_tWnu1_id", &mc_truth_tWnu1_id, "mc_truth_tWnu1_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_id", &mc_truth_tWnutau1_id, "mc_truth_tWnutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWl1_id", &mc_truth_tWl1_id, "mc_truth_tWl1_id/I", buffersize);
	tree->Branch("mc_truth_tWtau1_id", &mc_truth_tWtau1_id, "mc_truth_tWtau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_id", &mc_truth_tWtaunu1_id, "mc_truth_tWtaunu1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_id", &mc_truth_tWtaunutau1_id, "mc_truth_tWtaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_id", &mc_truth_tWtaul1_id, "mc_truth_tWtaul1_id/I", buffersize);
	tree->Branch("mc_truth_tWq11_id", &mc_truth_tWq11_id, "mc_truth_tWq11_id/I", buffersize);
	tree->Branch("mc_truth_tWq21_id", &mc_truth_tWq21_id, "mc_truth_tWq21_id/I", buffersize);
	tree->Branch("mc_truth_tWq11_IS_id", &mc_truth_tWq11_IS_id, "mc_truth_tWq11_IS_id/I", buffersize);
	tree->Branch("mc_truth_tWq21_IS_id", &mc_truth_tWq21_IS_id, "mc_truth_tWq21_IS_id/I", buffersize);

	tree->Branch("mc_truth_tW2_id", &mc_truth_tW2_id, "mc_truth_tW2_id/I", buffersize);
	tree->Branch("mc_truth_tWnu2_id", &mc_truth_tWnu2_id, "mc_truth_tWnu2_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_id", &mc_truth_tWnutau2_id, "mc_truth_tWnutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWl2_id", &mc_truth_tWl2_id, "mc_truth_tWl2_id/I", buffersize);
	tree->Branch("mc_truth_tWtau2_id", &mc_truth_tWtau2_id, "mc_truth_tWtau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_id", &mc_truth_tWtaunu2_id, "mc_truth_tWtaunu2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_id", &mc_truth_tWtaunutau2_id, "mc_truth_tWtaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_id", &mc_truth_tWtaul2_id, "mc_truth_tWtaul2_id/I", buffersize);
	tree->Branch("mc_truth_tWq12_id", &mc_truth_tWq12_id, "mc_truth_tWq12_id/I", buffersize);
	tree->Branch("mc_truth_tWq22_id", &mc_truth_tWq22_id, "mc_truth_tWq22_id/I", buffersize);
	tree->Branch("mc_truth_tWq12_IS_id", &mc_truth_tWq12_IS_id, "mc_truth_tWq12_IS_id/I", buffersize);
	tree->Branch("mc_truth_tWq22_IS_id", &mc_truth_tWq22_IS_id, "mc_truth_tWq22_IS_id/I", buffersize);

	tree->Branch("mc_truth_j1_id", &mc_truth_j1_id, "mc_truth_j1_id/I", buffersize);
	tree->Branch("mc_truth_j2_id", &mc_truth_j2_id, "mc_truth_j2_id/I", buffersize);
	tree->Branch("mc_truth_j3_id", &mc_truth_j3_id, "mc_truth_j3_id/I", buffersize);

	tree->Branch("mc_truth_h0_status", &mc_truth_h0_status, "mc_truth_h0_status/I", buffersize);

	tree->Branch("mc_truth_h0W1_status", &mc_truth_h0W1_status, "mc_truth_h0W1_status/I", buffersize);
	tree->Branch("mc_truth_h0W2_status", &mc_truth_h0W2_status, "mc_truth_h0W2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wl1_status", &mc_truth_h0Wl1_status, "mc_truth_h0Wl1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnu1_status", &mc_truth_h0Wnu1_status, "mc_truth_h0Wnu1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtau1_status", &mc_truth_h0Wtau1_status, "mc_truth_h0Wtau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_status", &mc_truth_h0Wnutau1_status, "mc_truth_h0Wnutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_status", &mc_truth_h0Wtaul1_status, "mc_truth_h0Wtaul1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_status", &mc_truth_h0Wtaunu1_status, "mc_truth_h0Wtaunu1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_status", &mc_truth_h0Wtaunutau1_status, "mc_truth_h0Wtaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wl2_status", &mc_truth_h0Wl2_status, "mc_truth_h0Wl2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnu2_status", &mc_truth_h0Wnu2_status, "mc_truth_h0Wnu2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtau2_status", &mc_truth_h0Wtau2_status, "mc_truth_h0Wtau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_status", &mc_truth_h0Wnutau2_status, "mc_truth_h0Wnutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_status", &mc_truth_h0Wtaul2_status, "mc_truth_h0Wtaul2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_status", &mc_truth_h0Wtaunu2_status, "mc_truth_h0Wtaunu2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_status", &mc_truth_h0Wtaunutau2_status, "mc_truth_h0Wtaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_status", &mc_truth_h0Wq11_status, "mc_truth_h0Wq11_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_status", &mc_truth_h0Wq21_status, "mc_truth_h0Wq21_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_status", &mc_truth_h0Wq12_status, "mc_truth_h0Wq12_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_status", &mc_truth_h0Wq22_status, "mc_truth_h0Wq22_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_status", &mc_truth_h0Wq11_IS_status, "mc_truth_h0Wq11_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_status", &mc_truth_h0Wq21_IS_status, "mc_truth_h0Wq21_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_status", &mc_truth_h0Wq12_IS_status, "mc_truth_h0Wq12_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_status", &mc_truth_h0Wq22_IS_status, "mc_truth_h0Wq22_IS_status/I", buffersize);

	tree->Branch("mc_truth_h0Z1_status", &mc_truth_h0Z1_status, "mc_truth_h0Z1_status/I", buffersize);
	tree->Branch("mc_truth_h0Z2_status", &mc_truth_h0Z2_status, "mc_truth_h0Z2_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl11_status", &mc_truth_h0Zl11_status, "mc_truth_h0Zl11_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl21_status", &mc_truth_h0Zl21_status, "mc_truth_h0Zl21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau11_status", &mc_truth_h0Ztau11_status, "mc_truth_h0Ztau11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau21_status", &mc_truth_h0Ztau21_status, "mc_truth_h0Ztau21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_status", &mc_truth_h0Ztaul11_status, "mc_truth_h0Ztaul11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_status", &mc_truth_h0Ztaul21_status, "mc_truth_h0Ztaul21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_status", &mc_truth_h0Ztaunu11_status, "mc_truth_h0Ztaunu11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_status", &mc_truth_h0Ztaunu21_status, "mc_truth_h0Ztaunu21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_status", &mc_truth_h0Ztaunutau11_status, "mc_truth_h0Ztaunutau11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_status", &mc_truth_h0Ztaunutau21_status, "mc_truth_h0Ztaunutau21_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_status", &mc_truth_h0Zq11_status, "mc_truth_h0Zq11_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_status", &mc_truth_h0Zq21_status, "mc_truth_h0Zq21_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_status", &mc_truth_h0Zq11_IS_status, "mc_truth_h0Zq11_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_status", &mc_truth_h0Zq21_IS_status, "mc_truth_h0Zq21_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl12_status", &mc_truth_h0Zl12_status, "mc_truth_h0Zl12_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl22_status", &mc_truth_h0Zl22_status, "mc_truth_h0Zl22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau12_status", &mc_truth_h0Ztau12_status, "mc_truth_h0Ztau12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau22_status", &mc_truth_h0Ztau22_status, "mc_truth_h0Ztau22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_status", &mc_truth_h0Ztaul12_status, "mc_truth_h0Ztaul12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_status", &mc_truth_h0Ztaul22_status, "mc_truth_h0Ztaul22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_status", &mc_truth_h0Ztaunu12_status, "mc_truth_h0Ztaunu12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_status", &mc_truth_h0Ztaunu22_status, "mc_truth_h0Ztaunu22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_status", &mc_truth_h0Ztaunutau12_status, "mc_truth_h0Ztaunutau12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_status", &mc_truth_h0Ztaunutau22_status, "mc_truth_h0Ztaunutau22_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_status", &mc_truth_h0Zq12_status, "mc_truth_h0Zq12_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_status", &mc_truth_h0Zq22_status, "mc_truth_h0Zq22_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_status", &mc_truth_h0Zq12_IS_status, "mc_truth_h0Zq12_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_status", &mc_truth_h0Zq22_IS_status, "mc_truth_h0Zq22_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu11_status", &mc_truth_h0Znu11_status, "mc_truth_h0Znu11_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu21_status", &mc_truth_h0Znu21_status, "mc_truth_h0Znu21_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu12_status", &mc_truth_h0Znu12_status, "mc_truth_h0Znu12_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu22_status", &mc_truth_h0Znu22_status, "mc_truth_h0Znu22_status/I", buffersize);

	tree->Branch("mc_truth_h0tau1_status", &mc_truth_h0tau1_status, "mc_truth_h0tau1_status/I", buffersize);
	tree->Branch("mc_truth_h0tau2_status", &mc_truth_h0tau2_status, "mc_truth_h0tau2_status/I", buffersize);
	tree->Branch("mc_truth_h0taul1_status", &mc_truth_h0taul1_status, "mc_truth_h0taul1_status/I", buffersize);
	tree->Branch("mc_truth_h0taunutau1_status", &mc_truth_h0taunutau1_status, "mc_truth_h0taunutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0taunu1_status", &mc_truth_h0taunu1_status, "mc_truth_h0taunu1_status/I", buffersize);
	tree->Branch("mc_truth_h0taul2_status", &mc_truth_h0taul2_status, "mc_truth_h0taul2_status/I", buffersize);
	tree->Branch("mc_truth_h0taunutau2_status", &mc_truth_h0taunutau2_status, "mc_truth_h0taunutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0taunu2_status", &mc_truth_h0taunu2_status, "mc_truth_h0taunu2_status/I", buffersize);

	tree->Branch("mc_truth_t1_status", &mc_truth_t1_status, "mc_truth_t1_status/I", buffersize);
	tree->Branch("mc_truth_t2_status", &mc_truth_t2_status, "mc_truth_t2_status/I", buffersize);
	tree->Branch("mc_truth_tb1_status", &mc_truth_tb1_status, "mc_truth_tb1_status/I", buffersize);
	tree->Branch("mc_truth_tb2_status", &mc_truth_tb2_status, "mc_truth_tb2_status/I", buffersize);
	tree->Branch("mc_truth_tb1_IS_status", &mc_truth_tb1_IS_status, "mc_truth_tb1_IS_status/I", buffersize);
	tree->Branch("mc_truth_tb2_IS_status", &mc_truth_tb2_IS_status, "mc_truth_tb2_IS_status/I", buffersize);

	tree->Branch("mc_truth_tW1_status", &mc_truth_tW1_status, "mc_truth_tW1_status/I", buffersize);
	tree->Branch("mc_truth_tWnu1_status", &mc_truth_tWnu1_status, "mc_truth_tWnu1_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_status", &mc_truth_tWnutau1_status, "mc_truth_tWnutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWl1_status", &mc_truth_tWl1_status, "mc_truth_tWl1_status/I", buffersize);
	tree->Branch("mc_truth_tWtau1_status", &mc_truth_tWtau1_status, "mc_truth_tWtau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_status", &mc_truth_tWtaunu1_status, "mc_truth_tWtaunu1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_status", &mc_truth_tWtaunutau1_status, "mc_truth_tWtaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_status", &mc_truth_tWtaul1_status, "mc_truth_tWtaul1_status/I", buffersize);
	tree->Branch("mc_truth_tWq11_status", &mc_truth_tWq11_status, "mc_truth_tWq11_status/I", buffersize);
	tree->Branch("mc_truth_tWq21_status", &mc_truth_tWq21_status, "mc_truth_tWq21_status/I", buffersize);
	tree->Branch("mc_truth_tWq11_IS_status", &mc_truth_tWq11_IS_status, "mc_truth_tWq11_IS_status/I", buffersize);
	tree->Branch("mc_truth_tWq21_IS_status", &mc_truth_tWq21_IS_status, "mc_truth_tWq21_IS_status/I", buffersize);

	tree->Branch("mc_truth_tW2_status", &mc_truth_tW2_status, "mc_truth_tW2_status/I", buffersize);
	tree->Branch("mc_truth_tWnu2_status", &mc_truth_tWnu2_status, "mc_truth_tWnu2_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_status", &mc_truth_tWnutau2_status, "mc_truth_tWnutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWl2_status", &mc_truth_tWl2_status, "mc_truth_tWl2_status/I", buffersize);
	tree->Branch("mc_truth_tWtau2_status", &mc_truth_tWtau2_status, "mc_truth_tWtau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_status", &mc_truth_tWtaunu2_status, "mc_truth_tWtaunu2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_status", &mc_truth_tWtaunutau2_status, "mc_truth_tWtaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_status", &mc_truth_tWtaul2_status, "mc_truth_tWtaul2_status/I", buffersize);
	tree->Branch("mc_truth_tWq12_status", &mc_truth_tWq12_status, "mc_truth_tWq12_status/I", buffersize);
	tree->Branch("mc_truth_tWq22_status", &mc_truth_tWq22_status, "mc_truth_tWq22_status/I", buffersize);
	tree->Branch("mc_truth_tWq12_IS_status", &mc_truth_tWq12_IS_status, "mc_truth_tWq12_IS_status/I", buffersize);
	tree->Branch("mc_truth_tWq22_IS_status", &mc_truth_tWq22_IS_status, "mc_truth_tWq22_IS_status/I", buffersize);

	tree->Branch("mc_truth_j1_status", &mc_truth_j1_status, "mc_truth_j1_status/I", buffersize);
	tree->Branch("mc_truth_j2_status", &mc_truth_j2_status, "mc_truth_j2_status/I", buffersize);
	tree->Branch("mc_truth_j3_status", &mc_truth_j3_status, "mc_truth_j3_status/I", buffersize);
     }

   if( doWrite("mc_truth_ttz") )
     {
	if( doWrite("mc_truth_p4") )
	  {
	     tree->Branch("mc_truth_Z_p4", "TLorentzVector", &mc_truth_Z_p4, buffersize);
	     tree->Branch("mc_truth_Zl1_p4", "TLorentzVector", &mc_truth_Zl1_p4, buffersize);
	     tree->Branch("mc_truth_Zl2_p4", "TLorentzVector", &mc_truth_Zl2_p4, buffersize);
	     tree->Branch("mc_truth_Ztau1_p4", "TLorentzVector", &mc_truth_Ztau1_p4, buffersize);
	     tree->Branch("mc_truth_Ztau2_p4", "TLorentzVector", &mc_truth_Ztau2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaul1_p4", "TLorentzVector", &mc_truth_Ztaul1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaul2_p4", "TLorentzVector", &mc_truth_Ztaul2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunu1_p4", "TLorentzVector", &mc_truth_Ztaunu1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunu2_p4", "TLorentzVector", &mc_truth_Ztaunu2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunutau1_p4", "TLorentzVector", &mc_truth_Ztaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunutau2_p4", "TLorentzVector", &mc_truth_Ztaunutau2_p4, buffersize);	
	     tree->Branch("mc_truth_Zq1_p4", "TLorentzVector", &mc_truth_Zq1_p4, buffersize);
	     tree->Branch("mc_truth_Zq2_p4", "TLorentzVector", &mc_truth_Zq2_p4, buffersize);
	     tree->Branch("mc_truth_Zq1_IS_p4", "TLorentzVector", &mc_truth_Zq1_IS_p4, buffersize);
	     tree->Branch("mc_truth_Zq2_IS_p4", "TLorentzVector", &mc_truth_Zq2_IS_p4, buffersize);
	     tree->Branch("mc_truth_Znu1_p4", "TLorentzVector", &mc_truth_Znu1_p4, buffersize);
	     tree->Branch("mc_truth_Znu2_p4", "TLorentzVector", &mc_truth_Znu2_p4, buffersize);

	     tree->Branch("mc_truth_gammal1_p4", "TLorentzVector", &mc_truth_gammal1_p4, buffersize);
	     tree->Branch("mc_truth_gammal2_p4", "TLorentzVector", &mc_truth_gammal2_p4, buffersize);
	     tree->Branch("mc_truth_gammatau1_p4", "TLorentzVector", &mc_truth_gammatau1_p4, buffersize);
	     tree->Branch("mc_truth_gammatau2_p4", "TLorentzVector", &mc_truth_gammatau2_p4, buffersize);
	     tree->Branch("mc_truth_gammataul1_p4", "TLorentzVector", &mc_truth_gammataul1_p4, buffersize);
	     tree->Branch("mc_truth_gammataul2_p4", "TLorentzVector", &mc_truth_gammataul2_p4, buffersize);
	     tree->Branch("mc_truth_gammataunu1_p4", "TLorentzVector", &mc_truth_gammataunu1_p4, buffersize);
	     tree->Branch("mc_truth_gammataunu2_p4", "TLorentzVector", &mc_truth_gammataunu2_p4, buffersize);
	     tree->Branch("mc_truth_gammataunutau1_p4", "TLorentzVector", &mc_truth_gammataunutau1_p4, buffersize);
	     tree->Branch("mc_truth_gammataunutau2_p4", "TLorentzVector", &mc_truth_gammataunutau2_p4, buffersize);	
	     
	     tree->Branch("mc_truth_t1_p4", "TLorentzVector", &mc_truth_t1_p4, buffersize);
	     tree->Branch("mc_truth_t2_p4", "TLorentzVector", &mc_truth_t2_p4, buffersize);
	     tree->Branch("mc_truth_tb1_p4", "TLorentzVector", &mc_truth_tb1_p4, buffersize);
	     tree->Branch("mc_truth_tb2_p4", "TLorentzVector", &mc_truth_tb2_p4, buffersize);
	     tree->Branch("mc_truth_tb1_IS_p4", "TLorentzVector", &mc_truth_tb1_IS_p4, buffersize);
	     tree->Branch("mc_truth_tb2_IS_p4", "TLorentzVector", &mc_truth_tb2_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW1_p4", "TLorentzVector", &mc_truth_tW1_p4, buffersize);
	     tree->Branch("mc_truth_tWnu1_p4", "TLorentzVector", &mc_truth_tWnu1_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau1_p4", "TLorentzVector", &mc_truth_tWnutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWl1_p4", "TLorentzVector", &mc_truth_tWl1_p4, buffersize);
	     tree->Branch("mc_truth_tWtau1_p4", "TLorentzVector", &mc_truth_tWtau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu1_p4", "TLorentzVector", &mc_truth_tWtaunu1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau1_p4", "TLorentzVector", &mc_truth_tWtaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul1_p4", "TLorentzVector", &mc_truth_tWtaul1_p4, buffersize);
	     tree->Branch("mc_truth_tWq11_p4", "TLorentzVector", &mc_truth_tWq11_p4, buffersize);
	     tree->Branch("mc_truth_tWq21_p4", "TLorentzVector", &mc_truth_tWq21_p4, buffersize);
	     tree->Branch("mc_truth_tWq11_IS_p4", "TLorentzVector", &mc_truth_tWq11_IS_p4, buffersize);
	     tree->Branch("mc_truth_tWq21_IS_p4", "TLorentzVector", &mc_truth_tWq21_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW2_p4", "TLorentzVector", &mc_truth_tW2_p4, buffersize);
	     tree->Branch("mc_truth_tWnu2_p4", "TLorentzVector", &mc_truth_tWnu2_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau2_p4", "TLorentzVector", &mc_truth_tWnutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWl2_p4", "TLorentzVector", &mc_truth_tWl2_p4, buffersize);
	     tree->Branch("mc_truth_tWtau2_p4", "TLorentzVector", &mc_truth_tWtau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu2_p4", "TLorentzVector", &mc_truth_tWtaunu2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau2_p4", "TLorentzVector", &mc_truth_tWtaunutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul2_p4", "TLorentzVector", &mc_truth_tWtaul2_p4, buffersize);
	     tree->Branch("mc_truth_tWq12_p4", "TLorentzVector", &mc_truth_tWq12_p4, buffersize);
	     tree->Branch("mc_truth_tWq22_p4", "TLorentzVector", &mc_truth_tWq22_p4, buffersize);
	     tree->Branch("mc_truth_tWq12_IS_p4", "TLorentzVector", &mc_truth_tWq12_IS_p4, buffersize);
	     tree->Branch("mc_truth_tWq22_IS_p4", "TLorentzVector", &mc_truth_tWq22_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_j1_p4", "TLorentzVector", &mc_truth_j1_p4, buffersize);
	     tree->Branch("mc_truth_j2_p4", "TLorentzVector", &mc_truth_j2_p4, buffersize);
	     tree->Branch("mc_truth_j3_p4", "TLorentzVector", &mc_truth_j3_p4, buffersize);
	  }

	tree->Branch("mc_truth_Z_pt", &mc_truth_Z_pt, "mc_truth_Z_pt/F", buffersize);
	tree->Branch("mc_truth_Zl1_pt", &mc_truth_Zl1_pt, "mc_truth_Zl1_pt/F", buffersize);
	tree->Branch("mc_truth_Zl2_pt", &mc_truth_Zl2_pt, "mc_truth_Zl2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztau1_pt", &mc_truth_Ztau1_pt, "mc_truth_Ztau1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztau2_pt", &mc_truth_Ztau2_pt, "mc_truth_Ztau2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_pt", &mc_truth_Ztaul1_pt, "mc_truth_Ztaul1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_pt", &mc_truth_Ztaul2_pt, "mc_truth_Ztaul2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_pt", &mc_truth_Ztaunu1_pt, "mc_truth_Ztaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_pt", &mc_truth_Ztaunu2_pt, "mc_truth_Ztaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_pt", &mc_truth_Ztaunutau1_pt, "mc_truth_Ztaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_pt", &mc_truth_Ztaunutau2_pt, "mc_truth_Ztaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_Zq1_pt", &mc_truth_Zq1_pt, "mc_truth_Zq1_pt/F", buffersize);
	tree->Branch("mc_truth_Zq2_pt", &mc_truth_Zq2_pt, "mc_truth_Zq2_pt/F", buffersize);
	tree->Branch("mc_truth_Zq1_IS_pt", &mc_truth_Zq1_IS_pt, "mc_truth_Zq1_IS_pt/F", buffersize);
	tree->Branch("mc_truth_Zq2_IS_pt", &mc_truth_Zq2_IS_pt, "mc_truth_Zq2_IS_pt/F", buffersize);
	tree->Branch("mc_truth_Znu1_pt", &mc_truth_Znu1_pt, "mc_truth_Znu1_pt/F", buffersize);
	tree->Branch("mc_truth_Znu2_pt", &mc_truth_Znu2_pt, "mc_truth_Znu2_pt/F", buffersize);

	tree->Branch("mc_truth_gammal1_pt", &mc_truth_gammal1_pt, "mc_truth_gammal1_pt/F", buffersize);
	tree->Branch("mc_truth_gammal2_pt", &mc_truth_gammal2_pt, "mc_truth_gammal2_pt/F", buffersize);
	tree->Branch("mc_truth_gammatau1_pt", &mc_truth_gammatau1_pt, "mc_truth_gammatau1_pt/F", buffersize);
	tree->Branch("mc_truth_gammatau2_pt", &mc_truth_gammatau2_pt, "mc_truth_gammatau2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataul1_pt", &mc_truth_gammataul1_pt, "mc_truth_gammataul1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataul2_pt", &mc_truth_gammataul2_pt, "mc_truth_gammataul2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_pt", &mc_truth_gammataunu1_pt, "mc_truth_gammataunu1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_pt", &mc_truth_gammataunu2_pt, "mc_truth_gammataunu2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_pt", &mc_truth_gammataunutau1_pt, "mc_truth_gammataunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_pt", &mc_truth_gammataunutau2_pt, "mc_truth_gammataunutau2_pt/F", buffersize);
	
	tree->Branch("mc_truth_t1_pt", &mc_truth_t1_pt, "mc_truth_t1_pt/F", buffersize);
	tree->Branch("mc_truth_t2_pt", &mc_truth_t2_pt, "mc_truth_t2_pt/F", buffersize);
	tree->Branch("mc_truth_tb1_pt", &mc_truth_tb1_pt, "mc_truth_tb1_pt/F", buffersize);
	tree->Branch("mc_truth_tb2_pt", &mc_truth_tb2_pt, "mc_truth_tb2_pt/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_pt", &mc_truth_tb1_IS_pt, "mc_truth_tb1_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_pt", &mc_truth_tb2_IS_pt, "mc_truth_tb2_IS_pt/F", buffersize);

	tree->Branch("mc_truth_tW1_pt", &mc_truth_tW1_pt, "mc_truth_tW1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu1_pt", &mc_truth_tWnu1_pt, "mc_truth_tWnu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_pt", &mc_truth_tWnutau1_pt, "mc_truth_tWnutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWl1_pt", &mc_truth_tWl1_pt, "mc_truth_tWl1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau1_pt", &mc_truth_tWtau1_pt, "mc_truth_tWtau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_pt", &mc_truth_tWtaunu1_pt, "mc_truth_tWtaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_pt", &mc_truth_tWtaunutau1_pt, "mc_truth_tWtaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_pt", &mc_truth_tWtaul1_pt, "mc_truth_tWtaul1_pt/F", buffersize);
	tree->Branch("mc_truth_tWq11_pt", &mc_truth_tWq11_pt, "mc_truth_tWq11_pt/F", buffersize);
	tree->Branch("mc_truth_tWq21_pt", &mc_truth_tWq21_pt, "mc_truth_tWq21_pt/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_pt", &mc_truth_tWq11_IS_pt, "mc_truth_tWq11_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_pt", &mc_truth_tWq21_IS_pt, "mc_truth_tWq21_IS_pt/F", buffersize);

	tree->Branch("mc_truth_tW2_pt", &mc_truth_tW2_pt, "mc_truth_tW2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu2_pt", &mc_truth_tWnu2_pt, "mc_truth_tWnu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_pt", &mc_truth_tWnutau2_pt, "mc_truth_tWnutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWl2_pt", &mc_truth_tWl2_pt, "mc_truth_tWl2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau2_pt", &mc_truth_tWtau2_pt, "mc_truth_tWtau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_pt", &mc_truth_tWtaunu2_pt, "mc_truth_tWtaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_pt", &mc_truth_tWtaunutau2_pt, "mc_truth_tWtaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_pt", &mc_truth_tWtaul2_pt, "mc_truth_tWtaul2_pt/F", buffersize);
	tree->Branch("mc_truth_tWq12_pt", &mc_truth_tWq12_pt, "mc_truth_tWq12_pt/F", buffersize);
	tree->Branch("mc_truth_tWq22_pt", &mc_truth_tWq22_pt, "mc_truth_tWq22_pt/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_pt", &mc_truth_tWq12_IS_pt, "mc_truth_tWq12_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_pt", &mc_truth_tWq22_IS_pt, "mc_truth_tWq22_IS_pt/F", buffersize);

	tree->Branch("mc_truth_j1_pt", &mc_truth_j1_pt, "mc_truth_j1_pt/F", buffersize);
	tree->Branch("mc_truth_j2_pt", &mc_truth_j2_pt, "mc_truth_j2_pt/F", buffersize);
	tree->Branch("mc_truth_j3_pt", &mc_truth_j3_pt, "mc_truth_j3_pt/F", buffersize);
		
	tree->Branch("mc_truth_Z_eta", &mc_truth_Z_eta, "mc_truth_Z_eta/F", buffersize);
	tree->Branch("mc_truth_Zl1_eta", &mc_truth_Zl1_eta, "mc_truth_Zl1_eta/F", buffersize);
	tree->Branch("mc_truth_Zl2_eta", &mc_truth_Zl2_eta, "mc_truth_Zl2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztau1_eta", &mc_truth_Ztau1_eta, "mc_truth_Ztau1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztau2_eta", &mc_truth_Ztau2_eta, "mc_truth_Ztau2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_eta", &mc_truth_Ztaul1_eta, "mc_truth_Ztaul1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_eta", &mc_truth_Ztaul2_eta, "mc_truth_Ztaul2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_eta", &mc_truth_Ztaunu1_eta, "mc_truth_Ztaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_eta", &mc_truth_Ztaunu2_eta, "mc_truth_Ztaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_eta", &mc_truth_Ztaunutau1_eta, "mc_truth_Ztaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_eta", &mc_truth_Ztaunutau2_eta, "mc_truth_Ztaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_Zq1_eta", &mc_truth_Zq1_eta, "mc_truth_Zq1_eta/F", buffersize);
	tree->Branch("mc_truth_Zq2_eta", &mc_truth_Zq2_eta, "mc_truth_Zq2_eta/F", buffersize);
	tree->Branch("mc_truth_Zq1_IS_eta", &mc_truth_Zq1_IS_eta, "mc_truth_Zq1_IS_eta/F", buffersize);
	tree->Branch("mc_truth_Zq2_IS_eta", &mc_truth_Zq2_IS_eta, "mc_truth_Zq2_IS_eta/F", buffersize);
	tree->Branch("mc_truth_Znu1_eta", &mc_truth_Znu1_eta, "mc_truth_Znu1_eta/F", buffersize);
	tree->Branch("mc_truth_Znu2_eta", &mc_truth_Znu2_eta, "mc_truth_Znu2_eta/F", buffersize);

	tree->Branch("mc_truth_gammal1_eta", &mc_truth_gammal1_eta, "mc_truth_gammal1_eta/F", buffersize);
	tree->Branch("mc_truth_gammal2_eta", &mc_truth_gammal2_eta, "mc_truth_gammal2_eta/F", buffersize);
	tree->Branch("mc_truth_gammatau1_eta", &mc_truth_gammatau1_eta, "mc_truth_gammatau1_eta/F", buffersize);
	tree->Branch("mc_truth_gammatau2_eta", &mc_truth_gammatau2_eta, "mc_truth_gammatau2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataul1_eta", &mc_truth_gammataul1_eta, "mc_truth_gammataul1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataul2_eta", &mc_truth_gammataul2_eta, "mc_truth_gammataul2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_eta", &mc_truth_gammataunu1_eta, "mc_truth_gammataunu1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_eta", &mc_truth_gammataunu2_eta, "mc_truth_gammataunu2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_eta", &mc_truth_gammataunutau1_eta, "mc_truth_gammataunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_eta", &mc_truth_gammataunutau2_eta, "mc_truth_gammataunutau2_eta/F", buffersize);
	
	tree->Branch("mc_truth_t1_eta", &mc_truth_t1_eta, "mc_truth_t1_eta/F", buffersize);
	tree->Branch("mc_truth_t2_eta", &mc_truth_t2_eta, "mc_truth_t2_eta/F", buffersize);
	tree->Branch("mc_truth_tb1_eta", &mc_truth_tb1_eta, "mc_truth_tb1_eta/F", buffersize);
	tree->Branch("mc_truth_tb2_eta", &mc_truth_tb2_eta, "mc_truth_tb2_eta/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_eta", &mc_truth_tb1_IS_eta, "mc_truth_tb1_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_eta", &mc_truth_tb2_IS_eta, "mc_truth_tb2_IS_eta/F", buffersize);

	tree->Branch("mc_truth_tW1_eta", &mc_truth_tW1_eta, "mc_truth_tW1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu1_eta", &mc_truth_tWnu1_eta, "mc_truth_tWnu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_eta", &mc_truth_tWnutau1_eta, "mc_truth_tWnutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWl1_eta", &mc_truth_tWl1_eta, "mc_truth_tWl1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau1_eta", &mc_truth_tWtau1_eta, "mc_truth_tWtau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_eta", &mc_truth_tWtaunu1_eta, "mc_truth_tWtaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_eta", &mc_truth_tWtaunutau1_eta, "mc_truth_tWtaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_eta", &mc_truth_tWtaul1_eta, "mc_truth_tWtaul1_eta/F", buffersize);
	tree->Branch("mc_truth_tWq11_eta", &mc_truth_tWq11_eta, "mc_truth_tWq11_eta/F", buffersize);
	tree->Branch("mc_truth_tWq21_eta", &mc_truth_tWq21_eta, "mc_truth_tWq21_eta/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_eta", &mc_truth_tWq11_IS_eta, "mc_truth_tWq11_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_eta", &mc_truth_tWq21_IS_eta, "mc_truth_tWq21_IS_eta/F", buffersize);

	tree->Branch("mc_truth_tW2_eta", &mc_truth_tW2_eta, "mc_truth_tW2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu2_eta", &mc_truth_tWnu2_eta, "mc_truth_tWnu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_eta", &mc_truth_tWnutau2_eta, "mc_truth_tWnutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWl2_eta", &mc_truth_tWl2_eta, "mc_truth_tWl2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau2_eta", &mc_truth_tWtau2_eta, "mc_truth_tWtau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_eta", &mc_truth_tWtaunu2_eta, "mc_truth_tWtaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_eta", &mc_truth_tWtaunutau2_eta, "mc_truth_tWtaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_eta", &mc_truth_tWtaul2_eta, "mc_truth_tWtaul2_eta/F", buffersize);
	tree->Branch("mc_truth_tWq12_eta", &mc_truth_tWq12_eta, "mc_truth_tWq12_eta/F", buffersize);
	tree->Branch("mc_truth_tWq22_eta", &mc_truth_tWq22_eta, "mc_truth_tWq22_eta/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_eta", &mc_truth_tWq12_IS_eta, "mc_truth_tWq12_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_eta", &mc_truth_tWq22_IS_eta, "mc_truth_tWq22_IS_eta/F", buffersize);

	tree->Branch("mc_truth_j1_eta", &mc_truth_j1_eta, "mc_truth_j1_eta/F", buffersize);
	tree->Branch("mc_truth_j2_eta", &mc_truth_j2_eta, "mc_truth_j2_eta/F", buffersize);
	tree->Branch("mc_truth_j3_eta", &mc_truth_j3_eta, "mc_truth_j3_eta/F", buffersize);
		
	tree->Branch("mc_truth_Z_phi", &mc_truth_Z_phi, "mc_truth_Z_phi/F", buffersize);
	tree->Branch("mc_truth_Zl1_phi", &mc_truth_Zl1_phi, "mc_truth_Zl1_phi/F", buffersize);
	tree->Branch("mc_truth_Zl2_phi", &mc_truth_Zl2_phi, "mc_truth_Zl2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztau1_phi", &mc_truth_Ztau1_phi, "mc_truth_Ztau1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztau2_phi", &mc_truth_Ztau2_phi, "mc_truth_Ztau2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_phi", &mc_truth_Ztaul1_phi, "mc_truth_Ztaul1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_phi", &mc_truth_Ztaul2_phi, "mc_truth_Ztaul2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_phi", &mc_truth_Ztaunu1_phi, "mc_truth_Ztaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_phi", &mc_truth_Ztaunu2_phi, "mc_truth_Ztaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_phi", &mc_truth_Ztaunutau1_phi, "mc_truth_Ztaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_phi", &mc_truth_Ztaunutau2_phi, "mc_truth_Ztaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_Zq1_phi", &mc_truth_Zq1_phi, "mc_truth_Zq1_phi/F", buffersize);
	tree->Branch("mc_truth_Zq2_phi", &mc_truth_Zq2_phi, "mc_truth_Zq2_phi/F", buffersize);
	tree->Branch("mc_truth_Zq1_IS_phi", &mc_truth_Zq1_IS_phi, "mc_truth_Zq1_IS_phi/F", buffersize);
	tree->Branch("mc_truth_Zq2_IS_phi", &mc_truth_Zq2_IS_phi, "mc_truth_Zq2_IS_phi/F", buffersize);
	tree->Branch("mc_truth_Znu1_phi", &mc_truth_Znu1_phi, "mc_truth_Znu1_phi/F", buffersize);
	tree->Branch("mc_truth_Znu2_phi", &mc_truth_Znu2_phi, "mc_truth_Znu2_phi/F", buffersize);

	tree->Branch("mc_truth_gammal1_phi", &mc_truth_gammal1_phi, "mc_truth_gammal1_phi/F", buffersize);
	tree->Branch("mc_truth_gammal2_phi", &mc_truth_gammal2_phi, "mc_truth_gammal2_phi/F", buffersize);
	tree->Branch("mc_truth_gammatau1_phi", &mc_truth_gammatau1_phi, "mc_truth_gammatau1_phi/F", buffersize);
	tree->Branch("mc_truth_gammatau2_phi", &mc_truth_gammatau2_phi, "mc_truth_gammatau2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataul1_phi", &mc_truth_gammataul1_phi, "mc_truth_gammataul1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataul2_phi", &mc_truth_gammataul2_phi, "mc_truth_gammataul2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_phi", &mc_truth_gammataunu1_phi, "mc_truth_gammataunu1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_phi", &mc_truth_gammataunu2_phi, "mc_truth_gammataunu2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_phi", &mc_truth_gammataunutau1_phi, "mc_truth_gammataunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_phi", &mc_truth_gammataunutau2_phi, "mc_truth_gammataunutau2_phi/F", buffersize);
	
	tree->Branch("mc_truth_t1_phi", &mc_truth_t1_phi, "mc_truth_t1_phi/F", buffersize);
	tree->Branch("mc_truth_t2_phi", &mc_truth_t2_phi, "mc_truth_t2_phi/F", buffersize);
	tree->Branch("mc_truth_tb1_phi", &mc_truth_tb1_phi, "mc_truth_tb1_phi/F", buffersize);
	tree->Branch("mc_truth_tb2_phi", &mc_truth_tb2_phi, "mc_truth_tb2_phi/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_phi", &mc_truth_tb1_IS_phi, "mc_truth_tb1_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_phi", &mc_truth_tb2_IS_phi, "mc_truth_tb2_IS_phi/F", buffersize);

	tree->Branch("mc_truth_tW1_phi", &mc_truth_tW1_phi, "mc_truth_tW1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu1_phi", &mc_truth_tWnu1_phi, "mc_truth_tWnu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_phi", &mc_truth_tWnutau1_phi, "mc_truth_tWnutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWl1_phi", &mc_truth_tWl1_phi, "mc_truth_tWl1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau1_phi", &mc_truth_tWtau1_phi, "mc_truth_tWtau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_phi", &mc_truth_tWtaunu1_phi, "mc_truth_tWtaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_phi", &mc_truth_tWtaunutau1_phi, "mc_truth_tWtaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_phi", &mc_truth_tWtaul1_phi, "mc_truth_tWtaul1_phi/F", buffersize);
	tree->Branch("mc_truth_tWq11_phi", &mc_truth_tWq11_phi, "mc_truth_tWq11_phi/F", buffersize);
	tree->Branch("mc_truth_tWq21_phi", &mc_truth_tWq21_phi, "mc_truth_tWq21_phi/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_phi", &mc_truth_tWq11_IS_phi, "mc_truth_tWq11_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_phi", &mc_truth_tWq21_IS_phi, "mc_truth_tWq21_IS_phi/F", buffersize);

	tree->Branch("mc_truth_tW2_phi", &mc_truth_tW2_phi, "mc_truth_tW2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu2_phi", &mc_truth_tWnu2_phi, "mc_truth_tWnu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_phi", &mc_truth_tWnutau2_phi, "mc_truth_tWnutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWl2_phi", &mc_truth_tWl2_phi, "mc_truth_tWl2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau2_phi", &mc_truth_tWtau2_phi, "mc_truth_tWtau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_phi", &mc_truth_tWtaunu2_phi, "mc_truth_tWtaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_phi", &mc_truth_tWtaunutau2_phi, "mc_truth_tWtaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_phi", &mc_truth_tWtaul2_phi, "mc_truth_tWtaul2_phi/F", buffersize);
	tree->Branch("mc_truth_tWq12_phi", &mc_truth_tWq12_phi, "mc_truth_tWq12_phi/F", buffersize);
	tree->Branch("mc_truth_tWq22_phi", &mc_truth_tWq22_phi, "mc_truth_tWq22_phi/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_phi", &mc_truth_tWq12_IS_phi, "mc_truth_tWq12_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_phi", &mc_truth_tWq22_IS_phi, "mc_truth_tWq22_IS_phi/F", buffersize);

	tree->Branch("mc_truth_j1_phi", &mc_truth_j1_phi, "mc_truth_j1_phi/F", buffersize);
	tree->Branch("mc_truth_j2_phi", &mc_truth_j2_phi, "mc_truth_j2_phi/F", buffersize);
	tree->Branch("mc_truth_j3_phi", &mc_truth_j3_phi, "mc_truth_j3_phi/F", buffersize);
	
	tree->Branch("mc_truth_Z_E", &mc_truth_Z_E, "mc_truth_Z_E/F", buffersize);
	tree->Branch("mc_truth_Zl1_E", &mc_truth_Zl1_E, "mc_truth_Zl1_E/F", buffersize);
	tree->Branch("mc_truth_Zl2_E", &mc_truth_Zl2_E, "mc_truth_Zl2_E/F", buffersize);
	tree->Branch("mc_truth_Ztau1_E", &mc_truth_Ztau1_E, "mc_truth_Ztau1_E/F", buffersize);
	tree->Branch("mc_truth_Ztau2_E", &mc_truth_Ztau2_E, "mc_truth_Ztau2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_E", &mc_truth_Ztaul1_E, "mc_truth_Ztaul1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_E", &mc_truth_Ztaul2_E, "mc_truth_Ztaul2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_E", &mc_truth_Ztaunu1_E, "mc_truth_Ztaunu1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_E", &mc_truth_Ztaunu2_E, "mc_truth_Ztaunu2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_E", &mc_truth_Ztaunutau1_E, "mc_truth_Ztaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_E", &mc_truth_Ztaunutau2_E, "mc_truth_Ztaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_Zq1_E", &mc_truth_Zq1_E, "mc_truth_Zq1_E/F", buffersize);
	tree->Branch("mc_truth_Zq2_E", &mc_truth_Zq2_E, "mc_truth_Zq2_E/F", buffersize);
	tree->Branch("mc_truth_Zq1_IS_E", &mc_truth_Zq1_IS_E, "mc_truth_Zq1_IS_E/F", buffersize);
	tree->Branch("mc_truth_Zq2_IS_E", &mc_truth_Zq2_IS_E, "mc_truth_Zq2_IS_E/F", buffersize);
	tree->Branch("mc_truth_Znu1_E", &mc_truth_Znu1_E, "mc_truth_Znu1_E/F", buffersize);
	tree->Branch("mc_truth_Znu2_E", &mc_truth_Znu2_E, "mc_truth_Znu2_E/F", buffersize);

	tree->Branch("mc_truth_gammal1_E", &mc_truth_gammal1_E, "mc_truth_gammal1_E/F", buffersize);
	tree->Branch("mc_truth_gammal2_E", &mc_truth_gammal2_E, "mc_truth_gammal2_E/F", buffersize);
	tree->Branch("mc_truth_gammatau1_E", &mc_truth_gammatau1_E, "mc_truth_gammatau1_E/F", buffersize);
	tree->Branch("mc_truth_gammatau2_E", &mc_truth_gammatau2_E, "mc_truth_gammatau2_E/F", buffersize);
	tree->Branch("mc_truth_gammataul1_E", &mc_truth_gammataul1_E, "mc_truth_gammataul1_E/F", buffersize);
	tree->Branch("mc_truth_gammataul2_E", &mc_truth_gammataul2_E, "mc_truth_gammataul2_E/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_E", &mc_truth_gammataunu1_E, "mc_truth_gammataunu1_E/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_E", &mc_truth_gammataunu2_E, "mc_truth_gammataunu2_E/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_E", &mc_truth_gammataunutau1_E, "mc_truth_gammataunutau1_E/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_E", &mc_truth_gammataunutau2_E, "mc_truth_gammataunutau2_E/F", buffersize);
	
	tree->Branch("mc_truth_t1_E", &mc_truth_t1_E, "mc_truth_t1_E/F", buffersize);
	tree->Branch("mc_truth_t2_E", &mc_truth_t2_E, "mc_truth_t2_E/F", buffersize);
	tree->Branch("mc_truth_tb1_E", &mc_truth_tb1_E, "mc_truth_tb1_E/F", buffersize);
	tree->Branch("mc_truth_tb2_E", &mc_truth_tb2_E, "mc_truth_tb2_E/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_E", &mc_truth_tb1_IS_E, "mc_truth_tb1_IS_E/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_E", &mc_truth_tb2_IS_E, "mc_truth_tb2_IS_E/F", buffersize);

	tree->Branch("mc_truth_tW1_E", &mc_truth_tW1_E, "mc_truth_tW1_E/F", buffersize);
	tree->Branch("mc_truth_tWnu1_E", &mc_truth_tWnu1_E, "mc_truth_tWnu1_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_E", &mc_truth_tWnutau1_E, "mc_truth_tWnutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWl1_E", &mc_truth_tWl1_E, "mc_truth_tWl1_E/F", buffersize);
	tree->Branch("mc_truth_tWtau1_E", &mc_truth_tWtau1_E, "mc_truth_tWtau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_E", &mc_truth_tWtaunu1_E, "mc_truth_tWtaunu1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_E", &mc_truth_tWtaunutau1_E, "mc_truth_tWtaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_E", &mc_truth_tWtaul1_E, "mc_truth_tWtaul1_E/F", buffersize);
	tree->Branch("mc_truth_tWq11_E", &mc_truth_tWq11_E, "mc_truth_tWq11_E/F", buffersize);
	tree->Branch("mc_truth_tWq21_E", &mc_truth_tWq21_E, "mc_truth_tWq21_E/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_E", &mc_truth_tWq11_IS_E, "mc_truth_tWq11_IS_E/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_E", &mc_truth_tWq21_IS_E, "mc_truth_tWq21_IS_E/F", buffersize);

	tree->Branch("mc_truth_tW2_E", &mc_truth_tW2_E, "mc_truth_tW2_E/F", buffersize);
	tree->Branch("mc_truth_tWnu2_E", &mc_truth_tWnu2_E, "mc_truth_tWnu2_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_E", &mc_truth_tWnutau2_E, "mc_truth_tWnutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWl2_E", &mc_truth_tWl2_E, "mc_truth_tWl2_E/F", buffersize);
	tree->Branch("mc_truth_tWtau2_E", &mc_truth_tWtau2_E, "mc_truth_tWtau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_E", &mc_truth_tWtaunu2_E, "mc_truth_tWtaunu2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_E", &mc_truth_tWtaunutau2_E, "mc_truth_tWtaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_E", &mc_truth_tWtaul2_E, "mc_truth_tWtaul2_E/F", buffersize);
	tree->Branch("mc_truth_tWq12_E", &mc_truth_tWq12_E, "mc_truth_tWq12_E/F", buffersize);
	tree->Branch("mc_truth_tWq22_E", &mc_truth_tWq22_E, "mc_truth_tWq22_E/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_E", &mc_truth_tWq12_IS_E, "mc_truth_tWq12_IS_E/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_E", &mc_truth_tWq22_IS_E, "mc_truth_tWq22_IS_E/F", buffersize);

	tree->Branch("mc_truth_j1_E", &mc_truth_j1_E, "mc_truth_j1_E/F", buffersize);
	tree->Branch("mc_truth_j2_E", &mc_truth_j2_E, "mc_truth_j2_E/F", buffersize);
	tree->Branch("mc_truth_j3_E", &mc_truth_j3_E, "mc_truth_j3_E/F", buffersize);
		
	tree->Branch("mc_truth_Z_id", &mc_truth_Z_id, "mc_truth_Z_id/I", buffersize);
	tree->Branch("mc_truth_Zl1_id", &mc_truth_Zl1_id, "mc_truth_Zl1_id/I", buffersize);
	tree->Branch("mc_truth_Zl2_id", &mc_truth_Zl2_id, "mc_truth_Zl2_id/I", buffersize);
	tree->Branch("mc_truth_Ztau1_id", &mc_truth_Ztau1_id, "mc_truth_Ztau1_id/I", buffersize);
	tree->Branch("mc_truth_Ztau2_id", &mc_truth_Ztau2_id, "mc_truth_Ztau2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaul1_id", &mc_truth_Ztaul1_id, "mc_truth_Ztaul1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaul2_id", &mc_truth_Ztaul2_id, "mc_truth_Ztaul2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunu1_id", &mc_truth_Ztaunu1_id, "mc_truth_Ztaunu1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunu2_id", &mc_truth_Ztaunu2_id, "mc_truth_Ztaunu2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_id", &mc_truth_Ztaunutau1_id, "mc_truth_Ztaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_id", &mc_truth_Ztaunutau2_id, "mc_truth_Ztaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_Zq1_id", &mc_truth_Zq1_id, "mc_truth_Zq1_id/I", buffersize);
	tree->Branch("mc_truth_Zq2_id", &mc_truth_Zq2_id, "mc_truth_Zq2_id/I", buffersize);
	tree->Branch("mc_truth_Zq1_IS_id", &mc_truth_Zq1_IS_id, "mc_truth_Zq1_IS_id/I", buffersize);
	tree->Branch("mc_truth_Zq2_IS_id", &mc_truth_Zq2_IS_id, "mc_truth_Zq2_IS_id/I", buffersize);
	tree->Branch("mc_truth_Znu1_id", &mc_truth_Znu1_id, "mc_truth_Znu1_id/I", buffersize);
	tree->Branch("mc_truth_Znu2_id", &mc_truth_Znu2_id, "mc_truth_Znu2_id/I", buffersize);

	tree->Branch("mc_truth_gammal1_id", &mc_truth_gammal1_id, "mc_truth_gammal1_id/I", buffersize);
	tree->Branch("mc_truth_gammal2_id", &mc_truth_gammal2_id, "mc_truth_gammal2_id/I", buffersize);
	tree->Branch("mc_truth_gammatau1_id", &mc_truth_gammatau1_id, "mc_truth_gammatau1_id/I", buffersize);
	tree->Branch("mc_truth_gammatau2_id", &mc_truth_gammatau2_id, "mc_truth_gammatau2_id/I", buffersize);
	tree->Branch("mc_truth_gammataul1_id", &mc_truth_gammataul1_id, "mc_truth_gammataul1_id/I", buffersize);
	tree->Branch("mc_truth_gammataul2_id", &mc_truth_gammataul2_id, "mc_truth_gammataul2_id/I", buffersize);
	tree->Branch("mc_truth_gammataunu1_id", &mc_truth_gammataunu1_id, "mc_truth_gammataunu1_id/I", buffersize);
	tree->Branch("mc_truth_gammataunu2_id", &mc_truth_gammataunu2_id, "mc_truth_gammataunu2_id/I", buffersize);
	tree->Branch("mc_truth_gammataunutau1_id", &mc_truth_gammataunutau1_id, "mc_truth_gammataunutau1_id/I", buffersize);
	tree->Branch("mc_truth_gammataunutau2_id", &mc_truth_gammataunutau2_id, "mc_truth_gammataunutau2_id/I", buffersize);
	
	tree->Branch("mc_truth_t1_id", &mc_truth_t1_id, "mc_truth_t1_id/I", buffersize);
	tree->Branch("mc_truth_t2_id", &mc_truth_t2_id, "mc_truth_t2_id/I", buffersize);
	tree->Branch("mc_truth_tb1_id", &mc_truth_tb1_id, "mc_truth_tb1_id/I", buffersize);
	tree->Branch("mc_truth_tb2_id", &mc_truth_tb2_id, "mc_truth_tb2_id/I", buffersize);
	tree->Branch("mc_truth_tb1_IS_id", &mc_truth_tb1_IS_id, "mc_truth_tb1_IS_id/I", buffersize);
	tree->Branch("mc_truth_tb2_IS_id", &mc_truth_tb2_IS_id, "mc_truth_tb2_IS_id/I", buffersize);

	tree->Branch("mc_truth_tW1_id", &mc_truth_tW1_id, "mc_truth_tW1_id/I", buffersize);
	tree->Branch("mc_truth_tWnu1_id", &mc_truth_tWnu1_id, "mc_truth_tWnu1_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_id", &mc_truth_tWnutau1_id, "mc_truth_tWnutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWl1_id", &mc_truth_tWl1_id, "mc_truth_tWl1_id/I", buffersize);
	tree->Branch("mc_truth_tWtau1_id", &mc_truth_tWtau1_id, "mc_truth_tWtau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_id", &mc_truth_tWtaunu1_id, "mc_truth_tWtaunu1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_id", &mc_truth_tWtaunutau1_id, "mc_truth_tWtaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_id", &mc_truth_tWtaul1_id, "mc_truth_tWtaul1_id/I", buffersize);
	tree->Branch("mc_truth_tWq11_id", &mc_truth_tWq11_id, "mc_truth_tWq11_id/I", buffersize);
	tree->Branch("mc_truth_tWq21_id", &mc_truth_tWq21_id, "mc_truth_tWq21_id/I", buffersize);
	tree->Branch("mc_truth_tWq11_IS_id", &mc_truth_tWq11_IS_id, "mc_truth_tWq11_IS_id/I", buffersize);
	tree->Branch("mc_truth_tWq21_IS_id", &mc_truth_tWq21_IS_id, "mc_truth_tWq21_IS_id/I", buffersize);

	tree->Branch("mc_truth_tW2_id", &mc_truth_tW2_id, "mc_truth_tW2_id/I", buffersize);
	tree->Branch("mc_truth_tWnu2_id", &mc_truth_tWnu2_id, "mc_truth_tWnu2_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_id", &mc_truth_tWnutau2_id, "mc_truth_tWnutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWl2_id", &mc_truth_tWl2_id, "mc_truth_tWl2_id/I", buffersize);
	tree->Branch("mc_truth_tWtau2_id", &mc_truth_tWtau2_id, "mc_truth_tWtau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_id", &mc_truth_tWtaunu2_id, "mc_truth_tWtaunu2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_id", &mc_truth_tWtaunutau2_id, "mc_truth_tWtaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_id", &mc_truth_tWtaul2_id, "mc_truth_tWtaul2_id/I", buffersize);
	tree->Branch("mc_truth_tWq12_id", &mc_truth_tWq12_id, "mc_truth_tWq12_id/I", buffersize);
	tree->Branch("mc_truth_tWq22_id", &mc_truth_tWq22_id, "mc_truth_tWq22_id/I", buffersize);
	tree->Branch("mc_truth_tWq12_IS_id", &mc_truth_tWq12_IS_id, "mc_truth_tWq12_IS_id/I", buffersize);
	tree->Branch("mc_truth_tWq22_IS_id", &mc_truth_tWq22_IS_id, "mc_truth_tWq22_IS_id/I", buffersize);

	tree->Branch("mc_truth_j1_id", &mc_truth_j1_id, "mc_truth_j1_id/I", buffersize);
	tree->Branch("mc_truth_j2_id", &mc_truth_j2_id, "mc_truth_j2_id/I", buffersize);
	tree->Branch("mc_truth_j3_id", &mc_truth_j3_id, "mc_truth_j3_id/I", buffersize);

	tree->Branch("mc_truth_Z_status", &mc_truth_Z_status, "mc_truth_Z_status/I", buffersize);
	tree->Branch("mc_truth_Zl1_status", &mc_truth_Zl1_status, "mc_truth_Zl1_status/I", buffersize);
	tree->Branch("mc_truth_Zl2_status", &mc_truth_Zl2_status, "mc_truth_Zl2_status/I", buffersize);
	tree->Branch("mc_truth_Ztau1_status", &mc_truth_Ztau1_status, "mc_truth_Ztau1_status/I", buffersize);
	tree->Branch("mc_truth_Ztau2_status", &mc_truth_Ztau2_status, "mc_truth_Ztau2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaul1_status", &mc_truth_Ztaul1_status, "mc_truth_Ztaul1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaul2_status", &mc_truth_Ztaul2_status, "mc_truth_Ztaul2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunu1_status", &mc_truth_Ztaunu1_status, "mc_truth_Ztaunu1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunu2_status", &mc_truth_Ztaunu2_status, "mc_truth_Ztaunu2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_status", &mc_truth_Ztaunutau1_status, "mc_truth_Ztaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_status", &mc_truth_Ztaunutau2_status, "mc_truth_Ztaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_Zq1_status", &mc_truth_Zq1_status, "mc_truth_Zq1_status/I", buffersize);
	tree->Branch("mc_truth_Zq2_status", &mc_truth_Zq2_status, "mc_truth_Zq2_status/I", buffersize);
	tree->Branch("mc_truth_Zq1_IS_status", &mc_truth_Zq1_IS_status, "mc_truth_Zq1_IS_status/I", buffersize);
	tree->Branch("mc_truth_Zq2_IS_status", &mc_truth_Zq2_IS_status, "mc_truth_Zq2_IS_status/I", buffersize);
	tree->Branch("mc_truth_Znu1_status", &mc_truth_Znu1_status, "mc_truth_Znu1_status/I", buffersize);
	tree->Branch("mc_truth_Znu2_status", &mc_truth_Znu2_status, "mc_truth_Znu2_status/I", buffersize);

	tree->Branch("mc_truth_gammal1_status", &mc_truth_gammal1_status, "mc_truth_gammal1_status/I", buffersize);
	tree->Branch("mc_truth_gammal2_status", &mc_truth_gammal2_status, "mc_truth_gammal2_status/I", buffersize);
	tree->Branch("mc_truth_gammatau1_status", &mc_truth_gammatau1_status, "mc_truth_gammatau1_status/I", buffersize);
	tree->Branch("mc_truth_gammatau2_status", &mc_truth_gammatau2_status, "mc_truth_gammatau2_status/I", buffersize);
	tree->Branch("mc_truth_gammataul1_status", &mc_truth_gammataul1_status, "mc_truth_gammataul1_status/I", buffersize);
	tree->Branch("mc_truth_gammataul2_status", &mc_truth_gammataul2_status, "mc_truth_gammataul2_status/I", buffersize);
	tree->Branch("mc_truth_gammataunu1_status", &mc_truth_gammataunu1_status, "mc_truth_gammataunu1_status/I", buffersize);
	tree->Branch("mc_truth_gammataunu2_status", &mc_truth_gammataunu2_status, "mc_truth_gammataunu2_status/I", buffersize);
	tree->Branch("mc_truth_gammataunutau1_status", &mc_truth_gammataunutau1_status, "mc_truth_gammataunutau1_status/I", buffersize);
	tree->Branch("mc_truth_gammataunutau2_status", &mc_truth_gammataunutau2_status, "mc_truth_gammataunutau2_status/I", buffersize);
	
	tree->Branch("mc_truth_t1_status", &mc_truth_t1_status, "mc_truth_t1_status/I", buffersize);
	tree->Branch("mc_truth_t2_status", &mc_truth_t2_status, "mc_truth_t2_status/I", buffersize);
	tree->Branch("mc_truth_tb1_status", &mc_truth_tb1_status, "mc_truth_tb1_status/I", buffersize);
	tree->Branch("mc_truth_tb2_status", &mc_truth_tb2_status, "mc_truth_tb2_status/I", buffersize);
	tree->Branch("mc_truth_tb1_IS_status", &mc_truth_tb1_IS_status, "mc_truth_tb1_IS_status/I", buffersize);
	tree->Branch("mc_truth_tb2_IS_status", &mc_truth_tb2_IS_status, "mc_truth_tb2_IS_status/I", buffersize);

	tree->Branch("mc_truth_tW1_status", &mc_truth_tW1_status, "mc_truth_tW1_status/I", buffersize);
	tree->Branch("mc_truth_tWnu1_status", &mc_truth_tWnu1_status, "mc_truth_tWnu1_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_status", &mc_truth_tWnutau1_status, "mc_truth_tWnutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWl1_status", &mc_truth_tWl1_status, "mc_truth_tWl1_status/I", buffersize);
	tree->Branch("mc_truth_tWtau1_status", &mc_truth_tWtau1_status, "mc_truth_tWtau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_status", &mc_truth_tWtaunu1_status, "mc_truth_tWtaunu1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_status", &mc_truth_tWtaunutau1_status, "mc_truth_tWtaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_status", &mc_truth_tWtaul1_status, "mc_truth_tWtaul1_status/I", buffersize);
	tree->Branch("mc_truth_tWq11_status", &mc_truth_tWq11_status, "mc_truth_tWq11_status/I", buffersize);
	tree->Branch("mc_truth_tWq21_status", &mc_truth_tWq21_status, "mc_truth_tWq21_status/I", buffersize);
	tree->Branch("mc_truth_tWq11_IS_status", &mc_truth_tWq11_IS_status, "mc_truth_tWq11_IS_status/I", buffersize);
	tree->Branch("mc_truth_tWq21_IS_status", &mc_truth_tWq21_IS_status, "mc_truth_tWq21_IS_status/I", buffersize);

	tree->Branch("mc_truth_tW2_status", &mc_truth_tW2_status, "mc_truth_tW2_status/I", buffersize);
	tree->Branch("mc_truth_tWnu2_status", &mc_truth_tWnu2_status, "mc_truth_tWnu2_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_status", &mc_truth_tWnutau2_status, "mc_truth_tWnutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWl2_status", &mc_truth_tWl2_status, "mc_truth_tWl2_status/I", buffersize);
	tree->Branch("mc_truth_tWtau2_status", &mc_truth_tWtau2_status, "mc_truth_tWtau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_status", &mc_truth_tWtaunu2_status, "mc_truth_tWtaunu2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_status", &mc_truth_tWtaunutau2_status, "mc_truth_tWtaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_status", &mc_truth_tWtaul2_status, "mc_truth_tWtaul2_status/I", buffersize);
	tree->Branch("mc_truth_tWq12_status", &mc_truth_tWq12_status, "mc_truth_tWq12_status/I", buffersize);
	tree->Branch("mc_truth_tWq22_status", &mc_truth_tWq22_status, "mc_truth_tWq22_status/I", buffersize);
	tree->Branch("mc_truth_tWq12_IS_status", &mc_truth_tWq12_IS_status, "mc_truth_tWq12_IS_status/I", buffersize);
	tree->Branch("mc_truth_tWq22_IS_status", &mc_truth_tWq22_IS_status, "mc_truth_tWq22_IS_status/I", buffersize);

	tree->Branch("mc_truth_j1_status", &mc_truth_j1_status, "mc_truth_j1_status/I", buffersize);
	tree->Branch("mc_truth_j2_status", &mc_truth_j2_status, "mc_truth_j2_status/I", buffersize);
	tree->Branch("mc_truth_j3_status", &mc_truth_j3_status, "mc_truth_j3_status/I", buffersize);
     }

   if( doWrite("mc_truth_ttw") )
     {
	if( doWrite("mc_truth_p4") )
	  {
	     tree->Branch("mc_truth_W_p4", "TLorentzVector", &mc_truth_W_p4, buffersize);
	     tree->Branch("mc_truth_Wnu_p4", "TLorentzVector", &mc_truth_Wnu_p4, buffersize);
	     tree->Branch("mc_truth_Wnutau_p4", "TLorentzVector", &mc_truth_Wnutau_p4, buffersize);
	     tree->Branch("mc_truth_Wl_p4", "TLorentzVector", &mc_truth_Wl_p4, buffersize);
	     tree->Branch("mc_truth_Wtau_p4", "TLorentzVector", &mc_truth_Wtau_p4, buffersize);
	     tree->Branch("mc_truth_Wtaunu_p4", "TLorentzVector", &mc_truth_Wtaunu_p4, buffersize);
	     tree->Branch("mc_truth_Wtaunutau_p4", "TLorentzVector", &mc_truth_Wtaunutau_p4, buffersize);
	     tree->Branch("mc_truth_Wtaul_p4", "TLorentzVector", &mc_truth_Wtaul_p4, buffersize);
	     tree->Branch("mc_truth_Wq1_p4", "TLorentzVector", &mc_truth_Wq1_p4, buffersize);
	     tree->Branch("mc_truth_Wq2_p4", "TLorentzVector", &mc_truth_Wq2_p4, buffersize);
	     tree->Branch("mc_truth_Wq1_IS_p4", "TLorentzVector", &mc_truth_Wq1_IS_p4, buffersize);
	     tree->Branch("mc_truth_Wq2_IS_p4", "TLorentzVector", &mc_truth_Wq2_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_gammal1_p4", "TLorentzVector", &mc_truth_gammal1_p4, buffersize);
	     tree->Branch("mc_truth_gammal2_p4", "TLorentzVector", &mc_truth_gammal2_p4, buffersize);
	     tree->Branch("mc_truth_gammatau1_p4", "TLorentzVector", &mc_truth_gammatau1_p4, buffersize);
	     tree->Branch("mc_truth_gammatau2_p4", "TLorentzVector", &mc_truth_gammatau2_p4, buffersize);
	     tree->Branch("mc_truth_gammataul1_p4", "TLorentzVector", &mc_truth_gammataul1_p4, buffersize);
	     tree->Branch("mc_truth_gammataul2_p4", "TLorentzVector", &mc_truth_gammataul2_p4, buffersize);
	     tree->Branch("mc_truth_gammataunu1_p4", "TLorentzVector", &mc_truth_gammataunu1_p4, buffersize);
	     tree->Branch("mc_truth_gammataunu2_p4", "TLorentzVector", &mc_truth_gammataunu2_p4, buffersize);
	     tree->Branch("mc_truth_gammataunutau1_p4", "TLorentzVector", &mc_truth_gammataunutau1_p4, buffersize);
	     tree->Branch("mc_truth_gammataunutau2_p4", "TLorentzVector", &mc_truth_gammataunutau2_p4, buffersize);	
	     
	     tree->Branch("mc_truth_t1_p4", "TLorentzVector", &mc_truth_t1_p4, buffersize);
	     tree->Branch("mc_truth_t2_p4", "TLorentzVector", &mc_truth_t2_p4, buffersize);
	     tree->Branch("mc_truth_tb1_p4", "TLorentzVector", &mc_truth_tb1_p4, buffersize);
	     tree->Branch("mc_truth_tb2_p4", "TLorentzVector", &mc_truth_tb2_p4, buffersize);
	     tree->Branch("mc_truth_tb1_IS_p4", "TLorentzVector", &mc_truth_tb1_IS_p4, buffersize);
	     tree->Branch("mc_truth_tb2_IS_p4", "TLorentzVector", &mc_truth_tb2_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW1_p4", "TLorentzVector", &mc_truth_tW1_p4, buffersize);
	     tree->Branch("mc_truth_tWnu1_p4", "TLorentzVector", &mc_truth_tWnu1_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau1_p4", "TLorentzVector", &mc_truth_tWnutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWl1_p4", "TLorentzVector", &mc_truth_tWl1_p4, buffersize);
	     tree->Branch("mc_truth_tWtau1_p4", "TLorentzVector", &mc_truth_tWtau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu1_p4", "TLorentzVector", &mc_truth_tWtaunu1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau1_p4", "TLorentzVector", &mc_truth_tWtaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul1_p4", "TLorentzVector", &mc_truth_tWtaul1_p4, buffersize);
	     tree->Branch("mc_truth_tWq11_p4", "TLorentzVector", &mc_truth_tWq11_p4, buffersize);
	     tree->Branch("mc_truth_tWq21_p4", "TLorentzVector", &mc_truth_tWq21_p4, buffersize);
	     tree->Branch("mc_truth_tWq11_IS_p4", "TLorentzVector", &mc_truth_tWq11_IS_p4, buffersize);
	     tree->Branch("mc_truth_tWq21_IS_p4", "TLorentzVector", &mc_truth_tWq21_IS_p4, buffersize);

	     tree->Branch("mc_truth_tW2_p4", "TLorentzVector", &mc_truth_tW2_p4, buffersize);
	     tree->Branch("mc_truth_tWnu2_p4", "TLorentzVector", &mc_truth_tWnu2_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau2_p4", "TLorentzVector", &mc_truth_tWnutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWl2_p4", "TLorentzVector", &mc_truth_tWl2_p4, buffersize);
	     tree->Branch("mc_truth_tWtau2_p4", "TLorentzVector", &mc_truth_tWtau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu2_p4", "TLorentzVector", &mc_truth_tWtaunu2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau2_p4", "TLorentzVector", &mc_truth_tWtaunutau2_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul2_p4", "TLorentzVector", &mc_truth_tWtaul2_p4, buffersize);
	     tree->Branch("mc_truth_tWq12_p4", "TLorentzVector", &mc_truth_tWq12_p4, buffersize);
	     tree->Branch("mc_truth_tWq22_p4", "TLorentzVector", &mc_truth_tWq22_p4, buffersize);
	     tree->Branch("mc_truth_tWq12_IS_p4", "TLorentzVector", &mc_truth_tWq12_IS_p4, buffersize);
	     tree->Branch("mc_truth_tWq22_IS_p4", "TLorentzVector", &mc_truth_tWq22_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_j1_p4", "TLorentzVector", &mc_truth_j1_p4, buffersize);
	     tree->Branch("mc_truth_j2_p4", "TLorentzVector", &mc_truth_j2_p4, buffersize);
	     tree->Branch("mc_truth_j3_p4", "TLorentzVector", &mc_truth_j3_p4, buffersize);
	  }	

	tree->Branch("mc_truth_W_pt", &mc_truth_W_pt, "mc_truth_W_pt/F", buffersize);
	tree->Branch("mc_truth_Wnu_pt", &mc_truth_Wnu_pt, "mc_truth_Wnu_pt/F", buffersize);
	tree->Branch("mc_truth_Wnutau_pt", &mc_truth_Wnutau_pt, "mc_truth_Wnutau_pt/F", buffersize);
	tree->Branch("mc_truth_Wl_pt", &mc_truth_Wl_pt, "mc_truth_Wl_pt/F", buffersize);
	tree->Branch("mc_truth_Wtau_pt", &mc_truth_Wtau_pt, "mc_truth_Wtau_pt/F", buffersize);
	tree->Branch("mc_truth_Wtaunu_pt", &mc_truth_Wtaunu_pt, "mc_truth_Wtaunu_pt/F", buffersize);
	tree->Branch("mc_truth_Wtaunutau_pt", &mc_truth_Wtaunutau_pt, "mc_truth_Wtaunutau_pt/F", buffersize);
	tree->Branch("mc_truth_Wtaul_pt", &mc_truth_Wtaul_pt, "mc_truth_Wtaul_pt/F", buffersize);
	tree->Branch("mc_truth_Wq1_pt", &mc_truth_Wq1_pt, "mc_truth_Wq1_pt/F", buffersize);
	tree->Branch("mc_truth_Wq2_pt", &mc_truth_Wq2_pt, "mc_truth_Wq2_pt/F", buffersize);
	tree->Branch("mc_truth_Wq1_IS_pt", &mc_truth_Wq1_IS_pt, "mc_truth_Wq1_IS_pt/F", buffersize);
	tree->Branch("mc_truth_Wq2_IS_pt", &mc_truth_Wq2_IS_pt, "mc_truth_Wq2_IS_pt/F", buffersize);

	tree->Branch("mc_truth_gammal1_pt", &mc_truth_gammal1_pt, "mc_truth_gammal1_pt/F", buffersize);
	tree->Branch("mc_truth_gammal2_pt", &mc_truth_gammal2_pt, "mc_truth_gammal2_pt/F", buffersize);
	tree->Branch("mc_truth_gammatau1_pt", &mc_truth_gammatau1_pt, "mc_truth_gammatau1_pt/F", buffersize);
	tree->Branch("mc_truth_gammatau2_pt", &mc_truth_gammatau2_pt, "mc_truth_gammatau2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataul1_pt", &mc_truth_gammataul1_pt, "mc_truth_gammataul1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataul2_pt", &mc_truth_gammataul2_pt, "mc_truth_gammataul2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_pt", &mc_truth_gammataunu1_pt, "mc_truth_gammataunu1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_pt", &mc_truth_gammataunu2_pt, "mc_truth_gammataunu2_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_pt", &mc_truth_gammataunutau1_pt, "mc_truth_gammataunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_pt", &mc_truth_gammataunutau2_pt, "mc_truth_gammataunutau2_pt/F", buffersize);
	
	tree->Branch("mc_truth_t1_pt", &mc_truth_t1_pt, "mc_truth_t1_pt/F", buffersize);
	tree->Branch("mc_truth_t2_pt", &mc_truth_t2_pt, "mc_truth_t2_pt/F", buffersize);
	tree->Branch("mc_truth_tb1_pt", &mc_truth_tb1_pt, "mc_truth_tb1_pt/F", buffersize);
	tree->Branch("mc_truth_tb2_pt", &mc_truth_tb2_pt, "mc_truth_tb2_pt/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_pt", &mc_truth_tb1_IS_pt, "mc_truth_tb1_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_pt", &mc_truth_tb2_IS_pt, "mc_truth_tb2_IS_pt/F", buffersize);

	tree->Branch("mc_truth_tW1_pt", &mc_truth_tW1_pt, "mc_truth_tW1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu1_pt", &mc_truth_tWnu1_pt, "mc_truth_tWnu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_pt", &mc_truth_tWnutau1_pt, "mc_truth_tWnutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWl1_pt", &mc_truth_tWl1_pt, "mc_truth_tWl1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau1_pt", &mc_truth_tWtau1_pt, "mc_truth_tWtau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_pt", &mc_truth_tWtaunu1_pt, "mc_truth_tWtaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_pt", &mc_truth_tWtaunutau1_pt, "mc_truth_tWtaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_pt", &mc_truth_tWtaul1_pt, "mc_truth_tWtaul1_pt/F", buffersize);
	tree->Branch("mc_truth_tWq11_pt", &mc_truth_tWq11_pt, "mc_truth_tWq11_pt/F", buffersize);
	tree->Branch("mc_truth_tWq21_pt", &mc_truth_tWq21_pt, "mc_truth_tWq21_pt/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_pt", &mc_truth_tWq11_IS_pt, "mc_truth_tWq11_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_pt", &mc_truth_tWq21_IS_pt, "mc_truth_tWq21_IS_pt/F", buffersize);

	tree->Branch("mc_truth_tW2_pt", &mc_truth_tW2_pt, "mc_truth_tW2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu2_pt", &mc_truth_tWnu2_pt, "mc_truth_tWnu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_pt", &mc_truth_tWnutau2_pt, "mc_truth_tWnutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWl2_pt", &mc_truth_tWl2_pt, "mc_truth_tWl2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau2_pt", &mc_truth_tWtau2_pt, "mc_truth_tWtau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_pt", &mc_truth_tWtaunu2_pt, "mc_truth_tWtaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_pt", &mc_truth_tWtaunutau2_pt, "mc_truth_tWtaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_pt", &mc_truth_tWtaul2_pt, "mc_truth_tWtaul2_pt/F", buffersize);
	tree->Branch("mc_truth_tWq12_pt", &mc_truth_tWq12_pt, "mc_truth_tWq12_pt/F", buffersize);
	tree->Branch("mc_truth_tWq22_pt", &mc_truth_tWq22_pt, "mc_truth_tWq22_pt/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_pt", &mc_truth_tWq12_IS_pt, "mc_truth_tWq12_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_pt", &mc_truth_tWq22_IS_pt, "mc_truth_tWq22_IS_pt/F", buffersize);

	tree->Branch("mc_truth_j1_pt", &mc_truth_j1_pt, "mc_truth_j1_pt/F", buffersize);
	tree->Branch("mc_truth_j2_pt", &mc_truth_j2_pt, "mc_truth_j2_pt/F", buffersize);
	tree->Branch("mc_truth_j3_pt", &mc_truth_j3_pt, "mc_truth_j3_pt/F", buffersize);
			
	tree->Branch("mc_truth_W_eta", &mc_truth_W_eta, "mc_truth_W_eta/F", buffersize);
	tree->Branch("mc_truth_Wnu_eta", &mc_truth_Wnu_eta, "mc_truth_Wnu_eta/F", buffersize);
	tree->Branch("mc_truth_Wnutau_eta", &mc_truth_Wnutau_eta, "mc_truth_Wnutau_eta/F", buffersize);
	tree->Branch("mc_truth_Wl_eta", &mc_truth_Wl_eta, "mc_truth_Wl_eta/F", buffersize);
	tree->Branch("mc_truth_Wtau_eta", &mc_truth_Wtau_eta, "mc_truth_Wtau_eta/F", buffersize);
	tree->Branch("mc_truth_Wtaunu_eta", &mc_truth_Wtaunu_eta, "mc_truth_Wtaunu_eta/F", buffersize);
	tree->Branch("mc_truth_Wtaunutau_eta", &mc_truth_Wtaunutau_eta, "mc_truth_Wtaunutau_eta/F", buffersize);
	tree->Branch("mc_truth_Wtaul_eta", &mc_truth_Wtaul_eta, "mc_truth_Wtaul_eta/F", buffersize);
	tree->Branch("mc_truth_Wq1_eta", &mc_truth_Wq1_eta, "mc_truth_Wq1_eta/F", buffersize);
	tree->Branch("mc_truth_Wq2_eta", &mc_truth_Wq2_eta, "mc_truth_Wq2_eta/F", buffersize);
	tree->Branch("mc_truth_Wq1_IS_eta", &mc_truth_Wq1_IS_eta, "mc_truth_Wq1_IS_eta/F", buffersize);
	tree->Branch("mc_truth_Wq2_IS_eta", &mc_truth_Wq2_IS_eta, "mc_truth_Wq2_IS_eta/F", buffersize);

	tree->Branch("mc_truth_gammal1_eta", &mc_truth_gammal1_eta, "mc_truth_gammal1_eta/F", buffersize);
	tree->Branch("mc_truth_gammal2_eta", &mc_truth_gammal2_eta, "mc_truth_gammal2_eta/F", buffersize);
	tree->Branch("mc_truth_gammatau1_eta", &mc_truth_gammatau1_eta, "mc_truth_gammatau1_eta/F", buffersize);
	tree->Branch("mc_truth_gammatau2_eta", &mc_truth_gammatau2_eta, "mc_truth_gammatau2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataul1_eta", &mc_truth_gammataul1_eta, "mc_truth_gammataul1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataul2_eta", &mc_truth_gammataul2_eta, "mc_truth_gammataul2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_eta", &mc_truth_gammataunu1_eta, "mc_truth_gammataunu1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_eta", &mc_truth_gammataunu2_eta, "mc_truth_gammataunu2_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_eta", &mc_truth_gammataunutau1_eta, "mc_truth_gammataunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_eta", &mc_truth_gammataunutau2_eta, "mc_truth_gammataunutau2_eta/F", buffersize);
	
	tree->Branch("mc_truth_t1_eta", &mc_truth_t1_eta, "mc_truth_t1_eta/F", buffersize);
	tree->Branch("mc_truth_t2_eta", &mc_truth_t2_eta, "mc_truth_t2_eta/F", buffersize);
	tree->Branch("mc_truth_tb1_eta", &mc_truth_tb1_eta, "mc_truth_tb1_eta/F", buffersize);
	tree->Branch("mc_truth_tb2_eta", &mc_truth_tb2_eta, "mc_truth_tb2_eta/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_eta", &mc_truth_tb1_IS_eta, "mc_truth_tb1_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_eta", &mc_truth_tb2_IS_eta, "mc_truth_tb2_IS_eta/F", buffersize);

	tree->Branch("mc_truth_tW1_eta", &mc_truth_tW1_eta, "mc_truth_tW1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu1_eta", &mc_truth_tWnu1_eta, "mc_truth_tWnu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_eta", &mc_truth_tWnutau1_eta, "mc_truth_tWnutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWl1_eta", &mc_truth_tWl1_eta, "mc_truth_tWl1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau1_eta", &mc_truth_tWtau1_eta, "mc_truth_tWtau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_eta", &mc_truth_tWtaunu1_eta, "mc_truth_tWtaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_eta", &mc_truth_tWtaunutau1_eta, "mc_truth_tWtaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_eta", &mc_truth_tWtaul1_eta, "mc_truth_tWtaul1_eta/F", buffersize);
	tree->Branch("mc_truth_tWq11_eta", &mc_truth_tWq11_eta, "mc_truth_tWq11_eta/F", buffersize);
	tree->Branch("mc_truth_tWq21_eta", &mc_truth_tWq21_eta, "mc_truth_tWq21_eta/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_eta", &mc_truth_tWq11_IS_eta, "mc_truth_tWq11_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_eta", &mc_truth_tWq21_IS_eta, "mc_truth_tWq21_IS_eta/F", buffersize);

	tree->Branch("mc_truth_tW2_eta", &mc_truth_tW2_eta, "mc_truth_tW2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu2_eta", &mc_truth_tWnu2_eta, "mc_truth_tWnu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_eta", &mc_truth_tWnutau2_eta, "mc_truth_tWnutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWl2_eta", &mc_truth_tWl2_eta, "mc_truth_tWl2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau2_eta", &mc_truth_tWtau2_eta, "mc_truth_tWtau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_eta", &mc_truth_tWtaunu2_eta, "mc_truth_tWtaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_eta", &mc_truth_tWtaunutau2_eta, "mc_truth_tWtaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_eta", &mc_truth_tWtaul2_eta, "mc_truth_tWtaul2_eta/F", buffersize);
	tree->Branch("mc_truth_tWq12_eta", &mc_truth_tWq12_eta, "mc_truth_tWq12_eta/F", buffersize);
	tree->Branch("mc_truth_tWq22_eta", &mc_truth_tWq22_eta, "mc_truth_tWq22_eta/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_eta", &mc_truth_tWq12_IS_eta, "mc_truth_tWq12_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_eta", &mc_truth_tWq22_IS_eta, "mc_truth_tWq22_IS_eta/F", buffersize);

	tree->Branch("mc_truth_j1_eta", &mc_truth_j1_eta, "mc_truth_j1_eta/F", buffersize);
	tree->Branch("mc_truth_j2_eta", &mc_truth_j2_eta, "mc_truth_j2_eta/F", buffersize);
	tree->Branch("mc_truth_j3_eta", &mc_truth_j3_eta, "mc_truth_j3_eta/F", buffersize);
		
	tree->Branch("mc_truth_W_phi", &mc_truth_W_phi, "mc_truth_W_phi/F", buffersize);
	tree->Branch("mc_truth_Wnu_phi", &mc_truth_Wnu_phi, "mc_truth_Wnu_phi/F", buffersize);
	tree->Branch("mc_truth_Wnutau_phi", &mc_truth_Wnutau_phi, "mc_truth_Wnutau_phi/F", buffersize);
	tree->Branch("mc_truth_Wl_phi", &mc_truth_Wl_phi, "mc_truth_Wl_phi/F", buffersize);
	tree->Branch("mc_truth_Wtau_phi", &mc_truth_Wtau_phi, "mc_truth_Wtau_phi/F", buffersize);
	tree->Branch("mc_truth_Wtaunu_phi", &mc_truth_Wtaunu_phi, "mc_truth_Wtaunu_phi/F", buffersize);
	tree->Branch("mc_truth_Wtaunutau_phi", &mc_truth_Wtaunutau_phi, "mc_truth_Wtaunutau_phi/F", buffersize);
	tree->Branch("mc_truth_Wtaul_phi", &mc_truth_Wtaul_phi, "mc_truth_Wtaul_phi/F", buffersize);
	tree->Branch("mc_truth_Wq1_phi", &mc_truth_Wq1_phi, "mc_truth_Wq1_phi/F", buffersize);
	tree->Branch("mc_truth_Wq2_phi", &mc_truth_Wq2_phi, "mc_truth_Wq2_phi/F", buffersize);
	tree->Branch("mc_truth_Wq1_IS_phi", &mc_truth_Wq1_IS_phi, "mc_truth_Wq1_IS_phi/F", buffersize);
	tree->Branch("mc_truth_Wq2_IS_phi", &mc_truth_Wq2_IS_phi, "mc_truth_Wq2_IS_phi/F", buffersize);

	tree->Branch("mc_truth_gammal1_phi", &mc_truth_gammal1_phi, "mc_truth_gammal1_phi/F", buffersize);
	tree->Branch("mc_truth_gammal2_phi", &mc_truth_gammal2_phi, "mc_truth_gammal2_phi/F", buffersize);
	tree->Branch("mc_truth_gammatau1_phi", &mc_truth_gammatau1_phi, "mc_truth_gammatau1_phi/F", buffersize);
	tree->Branch("mc_truth_gammatau2_phi", &mc_truth_gammatau2_phi, "mc_truth_gammatau2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataul1_phi", &mc_truth_gammataul1_phi, "mc_truth_gammataul1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataul2_phi", &mc_truth_gammataul2_phi, "mc_truth_gammataul2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_phi", &mc_truth_gammataunu1_phi, "mc_truth_gammataunu1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_phi", &mc_truth_gammataunu2_phi, "mc_truth_gammataunu2_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_phi", &mc_truth_gammataunutau1_phi, "mc_truth_gammataunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_phi", &mc_truth_gammataunutau2_phi, "mc_truth_gammataunutau2_phi/F", buffersize);
	
	tree->Branch("mc_truth_t1_phi", &mc_truth_t1_phi, "mc_truth_t1_phi/F", buffersize);
	tree->Branch("mc_truth_t2_phi", &mc_truth_t2_phi, "mc_truth_t2_phi/F", buffersize);
	tree->Branch("mc_truth_tb1_phi", &mc_truth_tb1_phi, "mc_truth_tb1_phi/F", buffersize);
	tree->Branch("mc_truth_tb2_phi", &mc_truth_tb2_phi, "mc_truth_tb2_phi/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_phi", &mc_truth_tb1_IS_phi, "mc_truth_tb1_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_phi", &mc_truth_tb2_IS_phi, "mc_truth_tb2_IS_phi/F", buffersize);

	tree->Branch("mc_truth_tW1_phi", &mc_truth_tW1_phi, "mc_truth_tW1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu1_phi", &mc_truth_tWnu1_phi, "mc_truth_tWnu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_phi", &mc_truth_tWnutau1_phi, "mc_truth_tWnutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWl1_phi", &mc_truth_tWl1_phi, "mc_truth_tWl1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau1_phi", &mc_truth_tWtau1_phi, "mc_truth_tWtau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_phi", &mc_truth_tWtaunu1_phi, "mc_truth_tWtaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_phi", &mc_truth_tWtaunutau1_phi, "mc_truth_tWtaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_phi", &mc_truth_tWtaul1_phi, "mc_truth_tWtaul1_phi/F", buffersize);
	tree->Branch("mc_truth_tWq11_phi", &mc_truth_tWq11_phi, "mc_truth_tWq11_phi/F", buffersize);
	tree->Branch("mc_truth_tWq21_phi", &mc_truth_tWq21_phi, "mc_truth_tWq21_phi/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_phi", &mc_truth_tWq11_IS_phi, "mc_truth_tWq11_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_phi", &mc_truth_tWq21_IS_phi, "mc_truth_tWq21_IS_phi/F", buffersize);

	tree->Branch("mc_truth_tW2_phi", &mc_truth_tW2_phi, "mc_truth_tW2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu2_phi", &mc_truth_tWnu2_phi, "mc_truth_tWnu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_phi", &mc_truth_tWnutau2_phi, "mc_truth_tWnutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWl2_phi", &mc_truth_tWl2_phi, "mc_truth_tWl2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau2_phi", &mc_truth_tWtau2_phi, "mc_truth_tWtau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_phi", &mc_truth_tWtaunu2_phi, "mc_truth_tWtaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_phi", &mc_truth_tWtaunutau2_phi, "mc_truth_tWtaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_phi", &mc_truth_tWtaul2_phi, "mc_truth_tWtaul2_phi/F", buffersize);
	tree->Branch("mc_truth_tWq12_phi", &mc_truth_tWq12_phi, "mc_truth_tWq12_phi/F", buffersize);
	tree->Branch("mc_truth_tWq22_phi", &mc_truth_tWq22_phi, "mc_truth_tWq22_phi/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_phi", &mc_truth_tWq12_IS_phi, "mc_truth_tWq12_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_phi", &mc_truth_tWq22_IS_phi, "mc_truth_tWq22_IS_phi/F", buffersize);

	tree->Branch("mc_truth_j1_phi", &mc_truth_j1_phi, "mc_truth_j1_phi/F", buffersize);
	tree->Branch("mc_truth_j2_phi", &mc_truth_j2_phi, "mc_truth_j2_phi/F", buffersize);
	tree->Branch("mc_truth_j3_phi", &mc_truth_j3_phi, "mc_truth_j3_phi/F", buffersize);

	tree->Branch("mc_truth_W_E", &mc_truth_W_E, "mc_truth_W_E/F", buffersize);
	tree->Branch("mc_truth_Wnu_E", &mc_truth_Wnu_E, "mc_truth_Wnu_E/F", buffersize);
	tree->Branch("mc_truth_Wnutau_E", &mc_truth_Wnutau_E, "mc_truth_Wnutau_E/F", buffersize);
	tree->Branch("mc_truth_Wl_E", &mc_truth_Wl_E, "mc_truth_Wl_E/F", buffersize);
	tree->Branch("mc_truth_Wtau_E", &mc_truth_Wtau_E, "mc_truth_Wtau_E/F", buffersize);
	tree->Branch("mc_truth_Wtaunu_E", &mc_truth_Wtaunu_E, "mc_truth_Wtaunu_E/F", buffersize);
	tree->Branch("mc_truth_Wtaunutau_E", &mc_truth_Wtaunutau_E, "mc_truth_Wtaunutau_E/F", buffersize);
	tree->Branch("mc_truth_Wtaul_E", &mc_truth_Wtaul_E, "mc_truth_Wtaul_E/F", buffersize);
	tree->Branch("mc_truth_Wq1_E", &mc_truth_Wq1_E, "mc_truth_Wq1_E/F", buffersize);
	tree->Branch("mc_truth_Wq2_E", &mc_truth_Wq2_E, "mc_truth_Wq2_E/F", buffersize);
	tree->Branch("mc_truth_Wq1_IS_E", &mc_truth_Wq1_IS_E, "mc_truth_Wq1_IS_E/F", buffersize);
	tree->Branch("mc_truth_Wq2_IS_E", &mc_truth_Wq2_IS_E, "mc_truth_Wq2_IS_E/F", buffersize);

	tree->Branch("mc_truth_gammal1_E", &mc_truth_gammal1_E, "mc_truth_gammal1_E/F", buffersize);
	tree->Branch("mc_truth_gammal2_E", &mc_truth_gammal2_E, "mc_truth_gammal2_E/F", buffersize);
	tree->Branch("mc_truth_gammatau1_E", &mc_truth_gammatau1_E, "mc_truth_gammatau1_E/F", buffersize);
	tree->Branch("mc_truth_gammatau2_E", &mc_truth_gammatau2_E, "mc_truth_gammatau2_E/F", buffersize);
	tree->Branch("mc_truth_gammataul1_E", &mc_truth_gammataul1_E, "mc_truth_gammataul1_E/F", buffersize);
	tree->Branch("mc_truth_gammataul2_E", &mc_truth_gammataul2_E, "mc_truth_gammataul2_E/F", buffersize);
	tree->Branch("mc_truth_gammataunu1_E", &mc_truth_gammataunu1_E, "mc_truth_gammataunu1_E/F", buffersize);
	tree->Branch("mc_truth_gammataunu2_E", &mc_truth_gammataunu2_E, "mc_truth_gammataunu2_E/F", buffersize);
	tree->Branch("mc_truth_gammataunutau1_E", &mc_truth_gammataunutau1_E, "mc_truth_gammataunutau1_E/F", buffersize);
	tree->Branch("mc_truth_gammataunutau2_E", &mc_truth_gammataunutau2_E, "mc_truth_gammataunutau2_E/F", buffersize);
	
	tree->Branch("mc_truth_t1_E", &mc_truth_t1_E, "mc_truth_t1_E/F", buffersize);
	tree->Branch("mc_truth_t2_E", &mc_truth_t2_E, "mc_truth_t2_E/F", buffersize);
	tree->Branch("mc_truth_tb1_E", &mc_truth_tb1_E, "mc_truth_tb1_E/F", buffersize);
	tree->Branch("mc_truth_tb2_E", &mc_truth_tb2_E, "mc_truth_tb2_E/F", buffersize);
	tree->Branch("mc_truth_tb1_IS_E", &mc_truth_tb1_IS_E, "mc_truth_tb1_IS_E/F", buffersize);
	tree->Branch("mc_truth_tb2_IS_E", &mc_truth_tb2_IS_E, "mc_truth_tb2_IS_E/F", buffersize);

	tree->Branch("mc_truth_tW1_E", &mc_truth_tW1_E, "mc_truth_tW1_E/F", buffersize);
	tree->Branch("mc_truth_tWnu1_E", &mc_truth_tWnu1_E, "mc_truth_tWnu1_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau1_E", &mc_truth_tWnutau1_E, "mc_truth_tWnutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWl1_E", &mc_truth_tWl1_E, "mc_truth_tWl1_E/F", buffersize);
	tree->Branch("mc_truth_tWtau1_E", &mc_truth_tWtau1_E, "mc_truth_tWtau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu1_E", &mc_truth_tWtaunu1_E, "mc_truth_tWtaunu1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_E", &mc_truth_tWtaunutau1_E, "mc_truth_tWtaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul1_E", &mc_truth_tWtaul1_E, "mc_truth_tWtaul1_E/F", buffersize);
	tree->Branch("mc_truth_tWq11_E", &mc_truth_tWq11_E, "mc_truth_tWq11_E/F", buffersize);
	tree->Branch("mc_truth_tWq21_E", &mc_truth_tWq21_E, "mc_truth_tWq21_E/F", buffersize);
	tree->Branch("mc_truth_tWq11_IS_E", &mc_truth_tWq11_IS_E, "mc_truth_tWq11_IS_E/F", buffersize);
	tree->Branch("mc_truth_tWq21_IS_E", &mc_truth_tWq21_IS_E, "mc_truth_tWq21_IS_E/F", buffersize);

	tree->Branch("mc_truth_tW2_E", &mc_truth_tW2_E, "mc_truth_tW2_E/F", buffersize);
	tree->Branch("mc_truth_tWnu2_E", &mc_truth_tWnu2_E, "mc_truth_tWnu2_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau2_E", &mc_truth_tWnutau2_E, "mc_truth_tWnutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWl2_E", &mc_truth_tWl2_E, "mc_truth_tWl2_E/F", buffersize);
	tree->Branch("mc_truth_tWtau2_E", &mc_truth_tWtau2_E, "mc_truth_tWtau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu2_E", &mc_truth_tWtaunu2_E, "mc_truth_tWtaunu2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_E", &mc_truth_tWtaunutau2_E, "mc_truth_tWtaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul2_E", &mc_truth_tWtaul2_E, "mc_truth_tWtaul2_E/F", buffersize);
	tree->Branch("mc_truth_tWq12_E", &mc_truth_tWq12_E, "mc_truth_tWq12_E/F", buffersize);
	tree->Branch("mc_truth_tWq22_E", &mc_truth_tWq22_E, "mc_truth_tWq22_E/F", buffersize);
	tree->Branch("mc_truth_tWq12_IS_E", &mc_truth_tWq12_IS_E, "mc_truth_tWq12_IS_E/F", buffersize);
	tree->Branch("mc_truth_tWq22_IS_E", &mc_truth_tWq22_IS_E, "mc_truth_tWq22_IS_E/F", buffersize);

	tree->Branch("mc_truth_j1_E", &mc_truth_j1_E, "mc_truth_j1_E/F", buffersize);
	tree->Branch("mc_truth_j2_E", &mc_truth_j2_E, "mc_truth_j2_E/F", buffersize);
	tree->Branch("mc_truth_j3_E", &mc_truth_j3_E, "mc_truth_j3_E/F", buffersize);
			
	tree->Branch("mc_truth_W_id", &mc_truth_W_id, "mc_truth_W_id/I", buffersize);
	tree->Branch("mc_truth_Wnu_id", &mc_truth_Wnu_id, "mc_truth_Wnu_id/I", buffersize);
	tree->Branch("mc_truth_Wnutau_id", &mc_truth_Wnutau_id, "mc_truth_Wnutau_id/I", buffersize);
	tree->Branch("mc_truth_Wl_id", &mc_truth_Wl_id, "mc_truth_Wl_id/I", buffersize);
	tree->Branch("mc_truth_Wtau_id", &mc_truth_Wtau_id, "mc_truth_Wtau_id/I", buffersize);
	tree->Branch("mc_truth_Wtaunu_id", &mc_truth_Wtaunu_id, "mc_truth_Wtaunu_id/I", buffersize);
	tree->Branch("mc_truth_Wtaunutau_id", &mc_truth_Wtaunutau_id, "mc_truth_Wtaunutau_id/I", buffersize);
	tree->Branch("mc_truth_Wtaul_id", &mc_truth_Wtaul_id, "mc_truth_Wtaul_id/I", buffersize);
	tree->Branch("mc_truth_Wq1_id", &mc_truth_Wq1_id, "mc_truth_Wq1_id/I", buffersize);
	tree->Branch("mc_truth_Wq2_id", &mc_truth_Wq2_id, "mc_truth_Wq2_id/I", buffersize);
	tree->Branch("mc_truth_Wq1_IS_id", &mc_truth_Wq1_IS_id, "mc_truth_Wq1_IS_id/I", buffersize);
	tree->Branch("mc_truth_Wq2_IS_id", &mc_truth_Wq2_IS_id, "mc_truth_Wq2_IS_id/I", buffersize);

	tree->Branch("mc_truth_gammal1_id", &mc_truth_gammal1_id, "mc_truth_gammal1_id/I", buffersize);
	tree->Branch("mc_truth_gammal2_id", &mc_truth_gammal2_id, "mc_truth_gammal2_id/I", buffersize);
	tree->Branch("mc_truth_gammatau1_id", &mc_truth_gammatau1_id, "mc_truth_gammatau1_id/I", buffersize);
	tree->Branch("mc_truth_gammatau2_id", &mc_truth_gammatau2_id, "mc_truth_gammatau2_id/I", buffersize);
	tree->Branch("mc_truth_gammataul1_id", &mc_truth_gammataul1_id, "mc_truth_gammataul1_id/I", buffersize);
	tree->Branch("mc_truth_gammataul2_id", &mc_truth_gammataul2_id, "mc_truth_gammataul2_id/I", buffersize);
	tree->Branch("mc_truth_gammataunu1_id", &mc_truth_gammataunu1_id, "mc_truth_gammataunu1_id/I", buffersize);
	tree->Branch("mc_truth_gammataunu2_id", &mc_truth_gammataunu2_id, "mc_truth_gammataunu2_id/I", buffersize);
	tree->Branch("mc_truth_gammataunutau1_id", &mc_truth_gammataunutau1_id, "mc_truth_gammataunutau1_id/I", buffersize);
	tree->Branch("mc_truth_gammataunutau2_id", &mc_truth_gammataunutau2_id, "mc_truth_gammataunutau2_id/I", buffersize);
	
	tree->Branch("mc_truth_t1_id", &mc_truth_t1_id, "mc_truth_t1_id/I", buffersize);
	tree->Branch("mc_truth_t2_id", &mc_truth_t2_id, "mc_truth_t2_id/I", buffersize);
	tree->Branch("mc_truth_tb1_id", &mc_truth_tb1_id, "mc_truth_tb1_id/I", buffersize);
	tree->Branch("mc_truth_tb2_id", &mc_truth_tb2_id, "mc_truth_tb2_id/I", buffersize);
	tree->Branch("mc_truth_tb1_IS_id", &mc_truth_tb1_IS_id, "mc_truth_tb1_IS_id/I", buffersize);
	tree->Branch("mc_truth_tb2_IS_id", &mc_truth_tb2_IS_id, "mc_truth_tb2_IS_id/I", buffersize);

	tree->Branch("mc_truth_tW1_id", &mc_truth_tW1_id, "mc_truth_tW1_id/I", buffersize);
	tree->Branch("mc_truth_tWnu1_id", &mc_truth_tWnu1_id, "mc_truth_tWnu1_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_id", &mc_truth_tWnutau1_id, "mc_truth_tWnutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWl1_id", &mc_truth_tWl1_id, "mc_truth_tWl1_id/I", buffersize);
	tree->Branch("mc_truth_tWtau1_id", &mc_truth_tWtau1_id, "mc_truth_tWtau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_id", &mc_truth_tWtaunu1_id, "mc_truth_tWtaunu1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_id", &mc_truth_tWtaunutau1_id, "mc_truth_tWtaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_id", &mc_truth_tWtaul1_id, "mc_truth_tWtaul1_id/I", buffersize);
	tree->Branch("mc_truth_tWq11_id", &mc_truth_tWq11_id, "mc_truth_tWq11_id/I", buffersize);
	tree->Branch("mc_truth_tWq21_id", &mc_truth_tWq21_id, "mc_truth_tWq21_id/I", buffersize);
	tree->Branch("mc_truth_tWq11_IS_id", &mc_truth_tWq11_IS_id, "mc_truth_tWq11_IS_id/I", buffersize);
	tree->Branch("mc_truth_tWq21_IS_id", &mc_truth_tWq21_IS_id, "mc_truth_tWq21_IS_id/I", buffersize);

	tree->Branch("mc_truth_tW2_id", &mc_truth_tW2_id, "mc_truth_tW2_id/I", buffersize);
	tree->Branch("mc_truth_tWnu2_id", &mc_truth_tWnu2_id, "mc_truth_tWnu2_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_id", &mc_truth_tWnutau2_id, "mc_truth_tWnutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWl2_id", &mc_truth_tWl2_id, "mc_truth_tWl2_id/I", buffersize);
	tree->Branch("mc_truth_tWtau2_id", &mc_truth_tWtau2_id, "mc_truth_tWtau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_id", &mc_truth_tWtaunu2_id, "mc_truth_tWtaunu2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_id", &mc_truth_tWtaunutau2_id, "mc_truth_tWtaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_id", &mc_truth_tWtaul2_id, "mc_truth_tWtaul2_id/I", buffersize);
	tree->Branch("mc_truth_tWq12_id", &mc_truth_tWq12_id, "mc_truth_tWq12_id/I", buffersize);
	tree->Branch("mc_truth_tWq22_id", &mc_truth_tWq22_id, "mc_truth_tWq22_id/I", buffersize);
	tree->Branch("mc_truth_tWq12_IS_id", &mc_truth_tWq12_IS_id, "mc_truth_tWq12_IS_id/I", buffersize);
	tree->Branch("mc_truth_tWq22_IS_id", &mc_truth_tWq22_IS_id, "mc_truth_tWq22_IS_id/I", buffersize);

	tree->Branch("mc_truth_j1_id", &mc_truth_j1_id, "mc_truth_j1_id/I", buffersize);
	tree->Branch("mc_truth_j2_id", &mc_truth_j2_id, "mc_truth_j2_id/I", buffersize);
	tree->Branch("mc_truth_j3_id", &mc_truth_j3_id, "mc_truth_j3_id/I", buffersize);

	tree->Branch("mc_truth_W_status", &mc_truth_W_status, "mc_truth_W_status/I", buffersize);
	tree->Branch("mc_truth_Wnu_status", &mc_truth_Wnu_status, "mc_truth_Wnu_status/I", buffersize);
	tree->Branch("mc_truth_Wnutau_status", &mc_truth_Wnutau_status, "mc_truth_Wnutau_status/I", buffersize);
	tree->Branch("mc_truth_Wl_status", &mc_truth_Wl_status, "mc_truth_Wl_status/I", buffersize);
	tree->Branch("mc_truth_Wtau_status", &mc_truth_Wtau_status, "mc_truth_Wtau_status/I", buffersize);
	tree->Branch("mc_truth_Wtaunu_status", &mc_truth_Wtaunu_status, "mc_truth_Wtaunu_status/I", buffersize);
	tree->Branch("mc_truth_Wtaunutau_status", &mc_truth_Wtaunutau_status, "mc_truth_Wtaunutau_status/I", buffersize);
	tree->Branch("mc_truth_Wtaul_status", &mc_truth_Wtaul_status, "mc_truth_Wtaul_status/I", buffersize);
	tree->Branch("mc_truth_Wq1_status", &mc_truth_Wq1_status, "mc_truth_Wq1_status/I", buffersize);
	tree->Branch("mc_truth_Wq2_status", &mc_truth_Wq2_status, "mc_truth_Wq2_status/I", buffersize);
	tree->Branch("mc_truth_Wq1_IS_status", &mc_truth_Wq1_IS_status, "mc_truth_Wq1_IS_status/I", buffersize);
	tree->Branch("mc_truth_Wq2_IS_status", &mc_truth_Wq2_IS_status, "mc_truth_Wq2_IS_status/I", buffersize);

	tree->Branch("mc_truth_gammal1_status", &mc_truth_gammal1_status, "mc_truth_gammal1_status/I", buffersize);
	tree->Branch("mc_truth_gammal2_status", &mc_truth_gammal2_status, "mc_truth_gammal2_status/I", buffersize);
	tree->Branch("mc_truth_gammatau1_status", &mc_truth_gammatau1_status, "mc_truth_gammatau1_status/I", buffersize);
	tree->Branch("mc_truth_gammatau2_status", &mc_truth_gammatau2_status, "mc_truth_gammatau2_status/I", buffersize);
	tree->Branch("mc_truth_gammataul1_status", &mc_truth_gammataul1_status, "mc_truth_gammataul1_status/I", buffersize);
	tree->Branch("mc_truth_gammataul2_status", &mc_truth_gammataul2_status, "mc_truth_gammataul2_status/I", buffersize);
	tree->Branch("mc_truth_gammataunu1_status", &mc_truth_gammataunu1_status, "mc_truth_gammataunu1_status/I", buffersize);
	tree->Branch("mc_truth_gammataunu2_status", &mc_truth_gammataunu2_status, "mc_truth_gammataunu2_status/I", buffersize);
	tree->Branch("mc_truth_gammataunutau1_status", &mc_truth_gammataunutau1_status, "mc_truth_gammataunutau1_status/I", buffersize);
	tree->Branch("mc_truth_gammataunutau2_status", &mc_truth_gammataunutau2_status, "mc_truth_gammataunutau2_status/I", buffersize);
	
	tree->Branch("mc_truth_t1_status", &mc_truth_t1_status, "mc_truth_t1_status/I", buffersize);
	tree->Branch("mc_truth_t2_status", &mc_truth_t2_status, "mc_truth_t2_status/I", buffersize);
	tree->Branch("mc_truth_tb1_status", &mc_truth_tb1_status, "mc_truth_tb1_status/I", buffersize);
	tree->Branch("mc_truth_tb2_status", &mc_truth_tb2_status, "mc_truth_tb2_status/I", buffersize);
	tree->Branch("mc_truth_tb1_IS_status", &mc_truth_tb1_IS_status, "mc_truth_tb1_IS_status/I", buffersize);
	tree->Branch("mc_truth_tb2_IS_status", &mc_truth_tb2_IS_status, "mc_truth_tb2_IS_status/I", buffersize);

	tree->Branch("mc_truth_tW1_status", &mc_truth_tW1_status, "mc_truth_tW1_status/I", buffersize);
	tree->Branch("mc_truth_tWnu1_status", &mc_truth_tWnu1_status, "mc_truth_tWnu1_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau1_status", &mc_truth_tWnutau1_status, "mc_truth_tWnutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWl1_status", &mc_truth_tWl1_status, "mc_truth_tWl1_status/I", buffersize);
	tree->Branch("mc_truth_tWtau1_status", &mc_truth_tWtau1_status, "mc_truth_tWtau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu1_status", &mc_truth_tWtaunu1_status, "mc_truth_tWtaunu1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau1_status", &mc_truth_tWtaunutau1_status, "mc_truth_tWtaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul1_status", &mc_truth_tWtaul1_status, "mc_truth_tWtaul1_status/I", buffersize);
	tree->Branch("mc_truth_tWq11_status", &mc_truth_tWq11_status, "mc_truth_tWq11_status/I", buffersize);
	tree->Branch("mc_truth_tWq21_status", &mc_truth_tWq21_status, "mc_truth_tWq21_status/I", buffersize);
	tree->Branch("mc_truth_tWq11_IS_status", &mc_truth_tWq11_IS_status, "mc_truth_tWq11_IS_status/I", buffersize);
	tree->Branch("mc_truth_tWq21_IS_status", &mc_truth_tWq21_IS_status, "mc_truth_tWq21_IS_status/I", buffersize);

	tree->Branch("mc_truth_tW2_status", &mc_truth_tW2_status, "mc_truth_tW2_status/I", buffersize);
	tree->Branch("mc_truth_tWnu2_status", &mc_truth_tWnu2_status, "mc_truth_tWnu2_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau2_status", &mc_truth_tWnutau2_status, "mc_truth_tWnutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWl2_status", &mc_truth_tWl2_status, "mc_truth_tWl2_status/I", buffersize);
	tree->Branch("mc_truth_tWtau2_status", &mc_truth_tWtau2_status, "mc_truth_tWtau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu2_status", &mc_truth_tWtaunu2_status, "mc_truth_tWtaunu2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau2_status", &mc_truth_tWtaunutau2_status, "mc_truth_tWtaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul2_status", &mc_truth_tWtaul2_status, "mc_truth_tWtaul2_status/I", buffersize);
	tree->Branch("mc_truth_tWq12_status", &mc_truth_tWq12_status, "mc_truth_tWq12_status/I", buffersize);
	tree->Branch("mc_truth_tWq22_status", &mc_truth_tWq22_status, "mc_truth_tWq22_status/I", buffersize);
	tree->Branch("mc_truth_tWq12_IS_status", &mc_truth_tWq12_IS_status, "mc_truth_tWq12_IS_status/I", buffersize);
	tree->Branch("mc_truth_tWq22_IS_status", &mc_truth_tWq22_IS_status, "mc_truth_tWq22_IS_status/I", buffersize);

	tree->Branch("mc_truth_j1_status", &mc_truth_j1_status, "mc_truth_j1_status/I", buffersize);
	tree->Branch("mc_truth_j2_status", &mc_truth_j2_status, "mc_truth_j2_status/I", buffersize);
	tree->Branch("mc_truth_j3_status", &mc_truth_j3_status, "mc_truth_j3_status/I", buffersize);
     }
   
   if( doWrite("mc_truth_tzq") )
     {
	tree->Branch("mc_truth_tzq_channel", &mc_truth_tzq_channel, "mc_truth_tzq_channel/I", buffersize);

	if( doWrite("mc_truth_p4") )
	  {	
	     tree->Branch("mc_truth_Z_p4", "TLorentzVector", &mc_truth_Z_p4, buffersize);
	     tree->Branch("mc_truth_Zl1_p4", "TLorentzVector", &mc_truth_Zl1_p4, buffersize);
	     tree->Branch("mc_truth_Zl2_p4", "TLorentzVector", &mc_truth_Zl2_p4, buffersize);
	     tree->Branch("mc_truth_Ztau1_p4", "TLorentzVector", &mc_truth_Ztau1_p4, buffersize);
	     tree->Branch("mc_truth_Ztau2_p4", "TLorentzVector", &mc_truth_Ztau2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaul1_p4", "TLorentzVector", &mc_truth_Ztaul1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaul2_p4", "TLorentzVector", &mc_truth_Ztaul2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunu1_p4", "TLorentzVector", &mc_truth_Ztaunu1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunu2_p4", "TLorentzVector", &mc_truth_Ztaunu2_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunutau1_p4", "TLorentzVector", &mc_truth_Ztaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_Ztaunutau2_p4", "TLorentzVector", &mc_truth_Ztaunutau2_p4, buffersize);
	     
	     tree->Branch("mc_truth_t_p4", "TLorentzVector", &mc_truth_t_p4, buffersize);
	     tree->Branch("mc_truth_tb_p4", "TLorentzVector", &mc_truth_tb_p4, buffersize);
	     tree->Branch("mc_truth_tb_IS_p4", "TLorentzVector", &mc_truth_tb_IS_p4, buffersize);
	     tree->Branch("mc_truth_tW_p4", "TLorentzVector", &mc_truth_tW_p4, buffersize);
	     tree->Branch("mc_truth_tWnu_p4", "TLorentzVector", &mc_truth_tWnu_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau_p4", "TLorentzVector", &mc_truth_tWnutau_p4, buffersize);
	     tree->Branch("mc_truth_tWl_p4", "TLorentzVector", &mc_truth_tWl_p4, buffersize);
	     tree->Branch("mc_truth_tWtau_p4", "TLorentzVector", &mc_truth_tWtau_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu_p4", "TLorentzVector", &mc_truth_tWtaunu_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau_p4", "TLorentzVector", &mc_truth_tWtaunutau_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul_p4", "TLorentzVector", &mc_truth_tWtaul_p4, buffersize);
	     tree->Branch("mc_truth_tWq1_p4", "TLorentzVector", &mc_truth_tWq1_p4, buffersize);
	     tree->Branch("mc_truth_tWq2_p4", "TLorentzVector", &mc_truth_tWq2_p4, buffersize);
	     tree->Branch("mc_truth_tWq1_IS_p4", "TLorentzVector", &mc_truth_tWq1_IS_p4, buffersize);
	     tree->Branch("mc_truth_tWq2_IS_p4", "TLorentzVector", &mc_truth_tWq2_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_j1_p4", "TLorentzVector", &mc_truth_j1_p4, buffersize);
	     tree->Branch("mc_truth_j2_p4", "TLorentzVector", &mc_truth_j2_p4, buffersize);
	     tree->Branch("mc_truth_j3_p4", "TLorentzVector", &mc_truth_j3_p4, buffersize);
	  }

	tree->Branch("mc_truth_Z_pt", &mc_truth_Z_pt, "mc_truth_Z_pt/F", buffersize);
	tree->Branch("mc_truth_Zl1_pt", &mc_truth_Zl1_pt, "mc_truth_Zl1_pt/F", buffersize);
	tree->Branch("mc_truth_Zl2_pt", &mc_truth_Zl2_pt, "mc_truth_Zl2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztau1_pt", &mc_truth_Ztau1_pt, "mc_truth_Ztau1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztau2_pt", &mc_truth_Ztau2_pt, "mc_truth_Ztau2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_pt", &mc_truth_Ztaul1_pt, "mc_truth_Ztaul1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_pt", &mc_truth_Ztaul2_pt, "mc_truth_Ztaul2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_pt", &mc_truth_Ztaunu1_pt, "mc_truth_Ztaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_pt", &mc_truth_Ztaunu2_pt, "mc_truth_Ztaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_pt", &mc_truth_Ztaunutau1_pt, "mc_truth_Ztaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_pt", &mc_truth_Ztaunutau2_pt, "mc_truth_Ztaunutau2_pt/F", buffersize);

	tree->Branch("mc_truth_t_pt", &mc_truth_t_pt, "mc_truth_t_pt/F", buffersize);
	tree->Branch("mc_truth_tb_pt", &mc_truth_tb_pt, "mc_truth_tb_pt/F", buffersize);
	tree->Branch("mc_truth_tb_IS_pt", &mc_truth_tb_IS_pt, "mc_truth_tb_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tW_pt", &mc_truth_tW_pt, "mc_truth_tW_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu_pt", &mc_truth_tWnu_pt, "mc_truth_tWnu_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau_pt", &mc_truth_tWnutau_pt, "mc_truth_tWnutau_pt/F", buffersize);
	tree->Branch("mc_truth_tWl_pt", &mc_truth_tWl_pt, "mc_truth_tWl_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau_pt", &mc_truth_tWtau_pt, "mc_truth_tWtau_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_pt", &mc_truth_tWtaunu_pt, "mc_truth_tWtaunu_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_pt", &mc_truth_tWtaunutau_pt, "mc_truth_tWtaunutau_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul_pt", &mc_truth_tWtaul_pt, "mc_truth_tWtaul_pt/F", buffersize);
	tree->Branch("mc_truth_tWq1_pt", &mc_truth_tWq1_pt, "mc_truth_tWq1_pt/F", buffersize);
	tree->Branch("mc_truth_tWq2_pt", &mc_truth_tWq2_pt, "mc_truth_tWq2_pt/F", buffersize);
	tree->Branch("mc_truth_tWq1_IS_pt", &mc_truth_tWq1_IS_pt, "mc_truth_tWq1_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tWq2_IS_pt", &mc_truth_tWq2_IS_pt, "mc_truth_tWq2_IS_pt/F", buffersize);

	tree->Branch("mc_truth_j1_pt", &mc_truth_j1_pt, "mc_truth_j1_pt/F", buffersize);
	tree->Branch("mc_truth_j2_pt", &mc_truth_j2_pt, "mc_truth_j2_pt/F", buffersize);
	tree->Branch("mc_truth_j3_pt", &mc_truth_j3_pt, "mc_truth_j3_pt/F", buffersize);	
	
	tree->Branch("mc_truth_Z_eta", &mc_truth_Z_eta, "mc_truth_Z_eta/F", buffersize);
	tree->Branch("mc_truth_Zl1_eta", &mc_truth_Zl1_eta, "mc_truth_Zl1_eta/F", buffersize);
	tree->Branch("mc_truth_Zl2_eta", &mc_truth_Zl2_eta, "mc_truth_Zl2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztau1_eta", &mc_truth_Ztau1_eta, "mc_truth_Ztau1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztau2_eta", &mc_truth_Ztau2_eta, "mc_truth_Ztau2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_eta", &mc_truth_Ztaul1_eta, "mc_truth_Ztaul1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_eta", &mc_truth_Ztaul2_eta, "mc_truth_Ztaul2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_eta", &mc_truth_Ztaunu1_eta, "mc_truth_Ztaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_eta", &mc_truth_Ztaunu2_eta, "mc_truth_Ztaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_eta", &mc_truth_Ztaunutau1_eta, "mc_truth_Ztaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_eta", &mc_truth_Ztaunutau2_eta, "mc_truth_Ztaunutau2_eta/F", buffersize);

	tree->Branch("mc_truth_t_eta", &mc_truth_t_eta, "mc_truth_t_eta/F", buffersize);
	tree->Branch("mc_truth_tb_eta", &mc_truth_tb_eta, "mc_truth_tb_eta/F", buffersize);
	tree->Branch("mc_truth_tb_IS_eta", &mc_truth_tb_IS_eta, "mc_truth_tb_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tW_eta", &mc_truth_tW_eta, "mc_truth_tW_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu_eta", &mc_truth_tWnu_eta, "mc_truth_tWnu_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau_eta", &mc_truth_tWnutau_eta, "mc_truth_tWnutau_eta/F", buffersize);
	tree->Branch("mc_truth_tWl_eta", &mc_truth_tWl_eta, "mc_truth_tWl_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau_eta", &mc_truth_tWtau_eta, "mc_truth_tWtau_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_eta", &mc_truth_tWtaunu_eta, "mc_truth_tWtaunu_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_eta", &mc_truth_tWtaunutau_eta, "mc_truth_tWtaunutau_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul_eta", &mc_truth_tWtaul_eta, "mc_truth_tWtaul_eta/F", buffersize);
	tree->Branch("mc_truth_tWq1_eta", &mc_truth_tWq1_eta, "mc_truth_tWq1_eta/F", buffersize);
	tree->Branch("mc_truth_tWq2_eta", &mc_truth_tWq2_eta, "mc_truth_tWq2_eta/F", buffersize);
	tree->Branch("mc_truth_tWq1_IS_eta", &mc_truth_tWq1_IS_eta, "mc_truth_tWq1_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tWq2_IS_eta", &mc_truth_tWq2_IS_eta, "mc_truth_tWq2_IS_eta/F", buffersize);

	tree->Branch("mc_truth_j1_eta", &mc_truth_j1_eta, "mc_truth_j1_eta/F", buffersize);
	tree->Branch("mc_truth_j2_eta", &mc_truth_j2_eta, "mc_truth_j2_eta/F", buffersize);
	tree->Branch("mc_truth_j3_eta", &mc_truth_j3_eta, "mc_truth_j3_eta/F", buffersize);
	
	tree->Branch("mc_truth_Z_phi", &mc_truth_Z_phi, "mc_truth_Z_phi/F", buffersize);
	tree->Branch("mc_truth_Zl1_phi", &mc_truth_Zl1_phi, "mc_truth_Zl1_phi/F", buffersize);
	tree->Branch("mc_truth_Zl2_phi", &mc_truth_Zl2_phi, "mc_truth_Zl2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztau1_phi", &mc_truth_Ztau1_phi, "mc_truth_Ztau1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztau2_phi", &mc_truth_Ztau2_phi, "mc_truth_Ztau2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_phi", &mc_truth_Ztaul1_phi, "mc_truth_Ztaul1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_phi", &mc_truth_Ztaul2_phi, "mc_truth_Ztaul2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_phi", &mc_truth_Ztaunu1_phi, "mc_truth_Ztaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_phi", &mc_truth_Ztaunu2_phi, "mc_truth_Ztaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_phi", &mc_truth_Ztaunutau1_phi, "mc_truth_Ztaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_phi", &mc_truth_Ztaunutau2_phi, "mc_truth_Ztaunutau2_phi/F", buffersize);

	tree->Branch("mc_truth_t_phi", &mc_truth_t_phi, "mc_truth_t_phi/F", buffersize);
	tree->Branch("mc_truth_tb_phi", &mc_truth_tb_phi, "mc_truth_tb_phi/F", buffersize);
	tree->Branch("mc_truth_tb_IS_phi", &mc_truth_tb_IS_phi, "mc_truth_tb_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tW_phi", &mc_truth_tW_phi, "mc_truth_tW_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu_phi", &mc_truth_tWnu_phi, "mc_truth_tWnu_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau_phi", &mc_truth_tWnutau_phi, "mc_truth_tWnutau_phi/F", buffersize);
	tree->Branch("mc_truth_tWl_phi", &mc_truth_tWl_phi, "mc_truth_tWl_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau_phi", &mc_truth_tWtau_phi, "mc_truth_tWtau_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_phi", &mc_truth_tWtaunu_phi, "mc_truth_tWtaunu_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_phi", &mc_truth_tWtaunutau_phi, "mc_truth_tWtaunutau_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul_phi", &mc_truth_tWtaul_phi, "mc_truth_tWtaul_phi/F", buffersize);
	tree->Branch("mc_truth_tWq1_phi", &mc_truth_tWq1_phi, "mc_truth_tWq1_phi/F", buffersize);
	tree->Branch("mc_truth_tWq2_phi", &mc_truth_tWq2_phi, "mc_truth_tWq2_phi/F", buffersize);
	tree->Branch("mc_truth_tWq1_IS_phi", &mc_truth_tWq1_IS_phi, "mc_truth_tWq1_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tWq2_IS_phi", &mc_truth_tWq2_IS_phi, "mc_truth_tWq2_IS_phi/F", buffersize);

	tree->Branch("mc_truth_j1_phi", &mc_truth_j1_phi, "mc_truth_j1_phi/F", buffersize);
	tree->Branch("mc_truth_j2_phi", &mc_truth_j2_phi, "mc_truth_j2_phi/F", buffersize);
	tree->Branch("mc_truth_j3_phi", &mc_truth_j3_phi, "mc_truth_j3_phi/F", buffersize);
	
	tree->Branch("mc_truth_Z_E", &mc_truth_Z_E, "mc_truth_Z_E/F", buffersize);
	tree->Branch("mc_truth_Zl1_E", &mc_truth_Zl1_E, "mc_truth_Zl1_E/F", buffersize);
	tree->Branch("mc_truth_Zl2_E", &mc_truth_Zl2_E, "mc_truth_Zl2_E/F", buffersize);
	tree->Branch("mc_truth_Ztau1_E", &mc_truth_Ztau1_E, "mc_truth_Ztau1_E/F", buffersize);
	tree->Branch("mc_truth_Ztau2_E", &mc_truth_Ztau2_E, "mc_truth_Ztau2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaul1_E", &mc_truth_Ztaul1_E, "mc_truth_Ztaul1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaul2_E", &mc_truth_Ztaul2_E, "mc_truth_Ztaul2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunu1_E", &mc_truth_Ztaunu1_E, "mc_truth_Ztaunu1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunu2_E", &mc_truth_Ztaunu2_E, "mc_truth_Ztaunu2_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_E", &mc_truth_Ztaunutau1_E, "mc_truth_Ztaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_E", &mc_truth_Ztaunutau2_E, "mc_truth_Ztaunutau2_E/F", buffersize);

	tree->Branch("mc_truth_t_E", &mc_truth_t_E, "mc_truth_t_E/F", buffersize);
	tree->Branch("mc_truth_tb_E", &mc_truth_tb_E, "mc_truth_tb_E/F", buffersize);
	tree->Branch("mc_truth_tb_IS_E", &mc_truth_tb_IS_E, "mc_truth_tb_IS_E/F", buffersize);
	tree->Branch("mc_truth_tW_E", &mc_truth_tW_E, "mc_truth_tW_E/F", buffersize);
	tree->Branch("mc_truth_tWnu_E", &mc_truth_tWnu_E, "mc_truth_tWnu_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau_E", &mc_truth_tWnutau_E, "mc_truth_tWnutau_E/F", buffersize);
	tree->Branch("mc_truth_tWl_E", &mc_truth_tWl_E, "mc_truth_tWl_E/F", buffersize);
	tree->Branch("mc_truth_tWtau_E", &mc_truth_tWtau_E, "mc_truth_tWtau_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_E", &mc_truth_tWtaunu_E, "mc_truth_tWtaunu_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_E", &mc_truth_tWtaunutau_E, "mc_truth_tWtaunutau_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul_E", &mc_truth_tWtaul_E, "mc_truth_tWtaul_E/F", buffersize);
	tree->Branch("mc_truth_tWq1_E", &mc_truth_tWq1_E, "mc_truth_tWq1_E/F", buffersize);
	tree->Branch("mc_truth_tWq2_E", &mc_truth_tWq2_E, "mc_truth_tWq2_E/F", buffersize);
	tree->Branch("mc_truth_tWq1_IS_E", &mc_truth_tWq1_IS_E, "mc_truth_tWq1_IS_E/F", buffersize);
	tree->Branch("mc_truth_tWq2_IS_E", &mc_truth_tWq2_IS_E, "mc_truth_tWq2_IS_E/F", buffersize);

	tree->Branch("mc_truth_j1_E", &mc_truth_j1_E, "mc_truth_j1_E/F", buffersize);
	tree->Branch("mc_truth_j2_E", &mc_truth_j2_E, "mc_truth_j2_E/F", buffersize);
	tree->Branch("mc_truth_j3_E", &mc_truth_j3_E, "mc_truth_j3_E/F", buffersize);
		
	tree->Branch("mc_truth_Z_id", &mc_truth_Z_id, "mc_truth_Z_id/I", buffersize);
	tree->Branch("mc_truth_Zl1_id", &mc_truth_Zl1_id, "mc_truth_Zl1_id/I", buffersize);
	tree->Branch("mc_truth_Zl2_id", &mc_truth_Zl2_id, "mc_truth_Zl2_id/I", buffersize);
	tree->Branch("mc_truth_Ztau1_id", &mc_truth_Ztau1_id, "mc_truth_Ztau1_id/I", buffersize);
	tree->Branch("mc_truth_Ztau2_id", &mc_truth_Ztau2_id, "mc_truth_Ztau2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaul1_id", &mc_truth_Ztaul1_id, "mc_truth_Ztaul1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaul2_id", &mc_truth_Ztaul2_id, "mc_truth_Ztaul2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunu1_id", &mc_truth_Ztaunu1_id, "mc_truth_Ztaunu1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunu2_id", &mc_truth_Ztaunu2_id, "mc_truth_Ztaunu2_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_id", &mc_truth_Ztaunutau1_id, "mc_truth_Ztaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_id", &mc_truth_Ztaunutau2_id, "mc_truth_Ztaunutau2_id/I", buffersize);

	tree->Branch("mc_truth_t_id", &mc_truth_t_id, "mc_truth_t_id/I", buffersize);
	tree->Branch("mc_truth_tb_id", &mc_truth_tb_id, "mc_truth_tb_id/I", buffersize);
	tree->Branch("mc_truth_tb_IS_id", &mc_truth_tb_IS_id, "mc_truth_tb_IS_id/I", buffersize);
	tree->Branch("mc_truth_tW_id", &mc_truth_tW_id, "mc_truth_tW_id/I", buffersize);
	tree->Branch("mc_truth_tWnu_id", &mc_truth_tWnu_id, "mc_truth_tWnu_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau_id", &mc_truth_tWnutau_id, "mc_truth_tWnutau_id/I", buffersize);
	tree->Branch("mc_truth_tWl_id", &mc_truth_tWl_id, "mc_truth_tWl_id/I", buffersize);
	tree->Branch("mc_truth_tWtau_id", &mc_truth_tWtau_id, "mc_truth_tWtau_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu_id", &mc_truth_tWtaunu_id, "mc_truth_tWtaunu_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau_id", &mc_truth_tWtaunutau_id, "mc_truth_tWtaunutau_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul_id", &mc_truth_tWtaul_id, "mc_truth_tWtaul_id/I", buffersize);
	tree->Branch("mc_truth_tWq1_id", &mc_truth_tWq1_id, "mc_truth_tWq1_id/I", buffersize);
	tree->Branch("mc_truth_tWq2_id", &mc_truth_tWq2_id, "mc_truth_tWq2_id/I", buffersize);
	tree->Branch("mc_truth_tWq1_IS_id", &mc_truth_tWq1_IS_id, "mc_truth_tWq1_IS_id/I", buffersize);
	tree->Branch("mc_truth_tWq2_IS_id", &mc_truth_tWq2_IS_id, "mc_truth_tWq2_IS_id/I", buffersize);

	tree->Branch("mc_truth_j1_id", &mc_truth_j1_id, "mc_truth_j1_id/I", buffersize);
	tree->Branch("mc_truth_j2_id", &mc_truth_j2_id, "mc_truth_j2_id/I", buffersize);
	tree->Branch("mc_truth_j3_id", &mc_truth_j3_id, "mc_truth_j3_id/I", buffersize);

	tree->Branch("mc_truth_Z_status", &mc_truth_Z_status, "mc_truth_Z_status/I", buffersize);
	tree->Branch("mc_truth_Zl1_status", &mc_truth_Zl1_status, "mc_truth_Zl1_status/I", buffersize);
	tree->Branch("mc_truth_Zl2_status", &mc_truth_Zl2_status, "mc_truth_Zl2_status/I", buffersize);
	tree->Branch("mc_truth_Ztau1_status", &mc_truth_Ztau1_status, "mc_truth_Ztau1_status/I", buffersize);
	tree->Branch("mc_truth_Ztau2_status", &mc_truth_Ztau2_status, "mc_truth_Ztau2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaul1_status", &mc_truth_Ztaul1_status, "mc_truth_Ztaul1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaul2_status", &mc_truth_Ztaul2_status, "mc_truth_Ztaul2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunu1_status", &mc_truth_Ztaunu1_status, "mc_truth_Ztaunu1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunu2_status", &mc_truth_Ztaunu2_status, "mc_truth_Ztaunu2_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau1_status", &mc_truth_Ztaunutau1_status, "mc_truth_Ztaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_Ztaunutau2_status", &mc_truth_Ztaunutau2_status, "mc_truth_Ztaunutau2_status/I", buffersize);

	tree->Branch("mc_truth_t_status", &mc_truth_t_status, "mc_truth_t_status/I", buffersize);
	tree->Branch("mc_truth_tb_status", &mc_truth_tb_status, "mc_truth_tb_status/I", buffersize);
	tree->Branch("mc_truth_tb_IS_status", &mc_truth_tb_IS_status, "mc_truth_tb_IS_status/I", buffersize);
	tree->Branch("mc_truth_tW_status", &mc_truth_tW_status, "mc_truth_tW_status/I", buffersize);
	tree->Branch("mc_truth_tWnu_status", &mc_truth_tWnu_status, "mc_truth_tWnu_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau_status", &mc_truth_tWnutau_status, "mc_truth_tWnutau_status/I", buffersize);
	tree->Branch("mc_truth_tWl_status", &mc_truth_tWl_status, "mc_truth_tWl_status/I", buffersize);
	tree->Branch("mc_truth_tWtau_status", &mc_truth_tWtau_status, "mc_truth_tWtau_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu_status", &mc_truth_tWtaunu_status, "mc_truth_tWtaunu_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau_status", &mc_truth_tWtaunutau_status, "mc_truth_tWtaunutau_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul_status", &mc_truth_tWtaul_status, "mc_truth_tWtaul_status/I", buffersize);
	tree->Branch("mc_truth_tWq1_status", &mc_truth_tWq1_status, "mc_truth_tWq1_status/I", buffersize);
	tree->Branch("mc_truth_tWq2_status", &mc_truth_tWq2_status, "mc_truth_tWq2_status/I", buffersize);
	tree->Branch("mc_truth_tWq1_IS_status", &mc_truth_tWq1_IS_status, "mc_truth_tWq1_IS_status/I", buffersize);
	tree->Branch("mc_truth_tWq2_IS_status", &mc_truth_tWq2_IS_status, "mc_truth_tWq2_IS_status/I", buffersize);

	tree->Branch("mc_truth_j1_status", &mc_truth_j1_status, "mc_truth_j1_status/I", buffersize);
	tree->Branch("mc_truth_j2_status", &mc_truth_j2_status, "mc_truth_j2_status/I", buffersize);
	tree->Branch("mc_truth_j3_status", &mc_truth_j3_status, "mc_truth_j3_status/I", buffersize);
     }

   if( doWrite("mc_truth_thq") )
     {
	tree->Branch("mc_truth_thq_channel", &mc_truth_thq_channel, "mc_truth_thq_channel/I", buffersize);

	if( doWrite("mc_truth_p4") )
	  {	
	     tree->Branch("mc_truth_h0_p4", "TLorentzVector", &mc_truth_h0_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0W1_p4", "TLorentzVector", &mc_truth_h0W1_p4, buffersize);
	     tree->Branch("mc_truth_h0W2_p4", "TLorentzVector", &mc_truth_h0W2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wl1_p4", "TLorentzVector", &mc_truth_h0Wl1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnu1_p4", "TLorentzVector", &mc_truth_h0Wnu1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtau1_p4", "TLorentzVector", &mc_truth_h0Wtau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnutau1_p4", "TLorentzVector", &mc_truth_h0Wnutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaul1_p4", "TLorentzVector", &mc_truth_h0Wtaul1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunu1_p4", "TLorentzVector", &mc_truth_h0Wtaunu1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunutau1_p4", "TLorentzVector", &mc_truth_h0Wtaunutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0Wl2_p4", "TLorentzVector", &mc_truth_h0Wl2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnu2_p4", "TLorentzVector", &mc_truth_h0Wnu2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtau2_p4", "TLorentzVector", &mc_truth_h0Wtau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wnutau2_p4", "TLorentzVector", &mc_truth_h0Wnutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaul2_p4", "TLorentzVector", &mc_truth_h0Wtaul2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunu2_p4", "TLorentzVector", &mc_truth_h0Wtaunu2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wtaunutau2_p4", "TLorentzVector", &mc_truth_h0Wtaunutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq11_p4", "TLorentzVector", &mc_truth_h0Wq11_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq21_p4", "TLorentzVector", &mc_truth_h0Wq21_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq12_p4", "TLorentzVector", &mc_truth_h0Wq12_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq22_p4", "TLorentzVector", &mc_truth_h0Wq22_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq11_IS_p4", "TLorentzVector", &mc_truth_h0Wq11_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq21_IS_p4", "TLorentzVector", &mc_truth_h0Wq21_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq12_IS_p4", "TLorentzVector", &mc_truth_h0Wq12_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Wq22_IS_p4", "TLorentzVector", &mc_truth_h0Wq22_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0Z1_p4", "TLorentzVector", &mc_truth_h0Z1_p4, buffersize);
	     tree->Branch("mc_truth_h0Z2_p4", "TLorentzVector", &mc_truth_h0Z2_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl11_p4", "TLorentzVector", &mc_truth_h0Zl11_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl21_p4", "TLorentzVector", &mc_truth_h0Zl21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl12_p4", "TLorentzVector", &mc_truth_h0Zl12_p4, buffersize);
	     tree->Branch("mc_truth_h0Zl22_p4", "TLorentzVector", &mc_truth_h0Zl22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau11_p4", "TLorentzVector", &mc_truth_h0Ztau11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau21_p4", "TLorentzVector", &mc_truth_h0Ztau21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul11_p4", "TLorentzVector", &mc_truth_h0Ztaul11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul21_p4", "TLorentzVector", &mc_truth_h0Ztaul21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu11_p4", "TLorentzVector", &mc_truth_h0Ztaunu11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu21_p4", "TLorentzVector", &mc_truth_h0Ztaunu21_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau11_p4", "TLorentzVector", &mc_truth_h0Ztaunutau11_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau21_p4", "TLorentzVector", &mc_truth_h0Ztaunutau21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq11_p4", "TLorentzVector", &mc_truth_h0Zq11_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq21_p4", "TLorentzVector", &mc_truth_h0Zq21_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq12_p4", "TLorentzVector", &mc_truth_h0Zq12_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq22_p4", "TLorentzVector", &mc_truth_h0Zq22_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq11_IS_p4", "TLorentzVector", &mc_truth_h0Zq11_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq21_IS_p4", "TLorentzVector", &mc_truth_h0Zq21_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq12_IS_p4", "TLorentzVector", &mc_truth_h0Zq12_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Zq22_IS_p4", "TLorentzVector", &mc_truth_h0Zq22_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau12_p4", "TLorentzVector", &mc_truth_h0Ztau12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztau22_p4", "TLorentzVector", &mc_truth_h0Ztau22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul12_p4", "TLorentzVector", &mc_truth_h0Ztaul12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaul22_p4", "TLorentzVector", &mc_truth_h0Ztaul22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu12_p4", "TLorentzVector", &mc_truth_h0Ztaunu12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunu22_p4", "TLorentzVector", &mc_truth_h0Ztaunu22_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau12_p4", "TLorentzVector", &mc_truth_h0Ztaunutau12_p4, buffersize);
	     tree->Branch("mc_truth_h0Ztaunutau22_p4", "TLorentzVector", &mc_truth_h0Ztaunutau22_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu11_p4", "TLorentzVector", &mc_truth_h0Znu11_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu21_p4", "TLorentzVector", &mc_truth_h0Znu21_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu12_p4", "TLorentzVector", &mc_truth_h0Znu12_p4, buffersize);
	     tree->Branch("mc_truth_h0Znu22_p4", "TLorentzVector", &mc_truth_h0Znu22_p4, buffersize);

	     tree->Branch("mc_truth_h0tau1_p4", "TLorentzVector", &mc_truth_h0tau1_p4, buffersize);
	     tree->Branch("mc_truth_h0tau2_p4", "TLorentzVector", &mc_truth_h0tau2_p4, buffersize);
	     tree->Branch("mc_truth_h0taul1_p4", "TLorentzVector", &mc_truth_h0taul1_p4, buffersize);
	     tree->Branch("mc_truth_h0taunutau1_p4", "TLorentzVector", &mc_truth_h0taunutau1_p4, buffersize);
	     tree->Branch("mc_truth_h0taunu1_p4", "TLorentzVector", &mc_truth_h0taunu1_p4, buffersize);
	     tree->Branch("mc_truth_h0taul2_p4", "TLorentzVector", &mc_truth_h0taul2_p4, buffersize);
	     tree->Branch("mc_truth_h0taunutau2_p4", "TLorentzVector", &mc_truth_h0taunutau2_p4, buffersize);
	     tree->Branch("mc_truth_h0taunu2_p4", "TLorentzVector", &mc_truth_h0taunu2_p4, buffersize);
	     
	     tree->Branch("mc_truth_h0b1_p4", "TLorentzVector", &mc_truth_h0b1_p4, buffersize);
	     tree->Branch("mc_truth_h0b2_p4", "TLorentzVector", &mc_truth_h0b2_p4, buffersize);
	     tree->Branch("mc_truth_h0b1_IS_p4", "TLorentzVector", &mc_truth_h0b1_IS_p4, buffersize);
	     tree->Branch("mc_truth_h0b2_IS_p4", "TLorentzVector", &mc_truth_h0b2_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_t_p4", "TLorentzVector", &mc_truth_t_p4, buffersize);
	     tree->Branch("mc_truth_tb_p4", "TLorentzVector", &mc_truth_tb_p4, buffersize);
	     tree->Branch("mc_truth_tb_IS_p4", "TLorentzVector", &mc_truth_tb_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_tW_p4", "TLorentzVector", &mc_truth_tW_p4, buffersize);
	     tree->Branch("mc_truth_tWnu_p4", "TLorentzVector", &mc_truth_tWnu_p4, buffersize);
	     tree->Branch("mc_truth_tWnutau_p4", "TLorentzVector", &mc_truth_tWnutau_p4, buffersize);
	     tree->Branch("mc_truth_tWl_p4", "TLorentzVector", &mc_truth_tWl_p4, buffersize);
	     tree->Branch("mc_truth_tWtau_p4", "TLorentzVector", &mc_truth_tWtau_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunu_p4", "TLorentzVector", &mc_truth_tWtaunu_p4, buffersize);
	     tree->Branch("mc_truth_tWtaunutau_p4", "TLorentzVector", &mc_truth_tWtaunutau_p4, buffersize);
	     tree->Branch("mc_truth_tWtaul_p4", "TLorentzVector", &mc_truth_tWtaul_p4, buffersize);
	     tree->Branch("mc_truth_tWq1_p4", "TLorentzVector", &mc_truth_tWq1_p4, buffersize);
	     tree->Branch("mc_truth_tWq2_p4", "TLorentzVector", &mc_truth_tWq2_p4, buffersize);
	     tree->Branch("mc_truth_tWq1_IS_p4", "TLorentzVector", &mc_truth_tWq1_IS_p4, buffersize);
	     tree->Branch("mc_truth_tWq2_IS_p4", "TLorentzVector", &mc_truth_tWq2_IS_p4, buffersize);
	     
	     tree->Branch("mc_truth_j1_p4", "TLorentzVector", &mc_truth_j1_p4, buffersize);
	     tree->Branch("mc_truth_j2_p4", "TLorentzVector", &mc_truth_j2_p4, buffersize);
	     tree->Branch("mc_truth_j3_p4", "TLorentzVector", &mc_truth_j3_p4, buffersize);
	  }

	tree->Branch("mc_truth_h0_pt", &mc_truth_h0_pt, "mc_truth_h0_pt/F", buffersize);

	tree->Branch("mc_truth_h0W1_pt", &mc_truth_h0W1_pt, "mc_truth_h0W1_pt/F", buffersize);
	tree->Branch("mc_truth_h0W2_pt", &mc_truth_h0W2_pt, "mc_truth_h0W2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_pt", &mc_truth_h0Wl1_pt, "mc_truth_h0Wl1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_pt", &mc_truth_h0Wnu1_pt, "mc_truth_h0Wnu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_pt", &mc_truth_h0Wtau1_pt, "mc_truth_h0Wtau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_pt", &mc_truth_h0Wnutau1_pt, "mc_truth_h0Wnutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_pt", &mc_truth_h0Wtaul1_pt, "mc_truth_h0Wtaul1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_pt", &mc_truth_h0Wtaunu1_pt, "mc_truth_h0Wtaunu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_pt", &mc_truth_h0Wtaunutau1_pt, "mc_truth_h0Wtaunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_pt", &mc_truth_h0Wl2_pt, "mc_truth_h0Wl2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_pt", &mc_truth_h0Wnu2_pt, "mc_truth_h0Wnu2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_pt", &mc_truth_h0Wtau2_pt, "mc_truth_h0Wtau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_pt", &mc_truth_h0Wnutau2_pt, "mc_truth_h0Wnutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_pt", &mc_truth_h0Wtaul2_pt, "mc_truth_h0Wtaul2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_pt", &mc_truth_h0Wtaunu2_pt, "mc_truth_h0Wtaunu2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_pt", &mc_truth_h0Wtaunutau2_pt, "mc_truth_h0Wtaunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_pt", &mc_truth_h0Wq11_pt, "mc_truth_h0Wq11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_pt", &mc_truth_h0Wq21_pt, "mc_truth_h0Wq21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_pt", &mc_truth_h0Wq12_pt, "mc_truth_h0Wq12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_pt", &mc_truth_h0Wq22_pt, "mc_truth_h0Wq22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_pt", &mc_truth_h0Wq11_IS_pt, "mc_truth_h0Wq11_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_pt", &mc_truth_h0Wq21_IS_pt, "mc_truth_h0Wq21_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_pt", &mc_truth_h0Wq12_IS_pt, "mc_truth_h0Wq12_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_pt", &mc_truth_h0Wq22_IS_pt, "mc_truth_h0Wq22_IS_pt/F", buffersize);

	tree->Branch("mc_truth_h0Z1_pt", &mc_truth_h0Z1_pt, "mc_truth_h0Z1_pt/F", buffersize);
	tree->Branch("mc_truth_h0Z2_pt", &mc_truth_h0Z2_pt, "mc_truth_h0Z2_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_pt", &mc_truth_h0Zl11_pt, "mc_truth_h0Zl11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_pt", &mc_truth_h0Zl21_pt, "mc_truth_h0Zl21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_pt", &mc_truth_h0Ztau11_pt, "mc_truth_h0Ztau11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_pt", &mc_truth_h0Ztau21_pt, "mc_truth_h0Ztau21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_pt", &mc_truth_h0Ztaul11_pt, "mc_truth_h0Ztaul11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_pt", &mc_truth_h0Ztaul21_pt, "mc_truth_h0Ztaul21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_pt", &mc_truth_h0Ztaunu11_pt, "mc_truth_h0Ztaunu11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_pt", &mc_truth_h0Ztaunu21_pt, "mc_truth_h0Ztaunu21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_pt", &mc_truth_h0Ztaunutau11_pt, "mc_truth_h0Ztaunutau11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_pt", &mc_truth_h0Ztaunutau21_pt, "mc_truth_h0Ztaunutau21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_pt", &mc_truth_h0Zq11_pt, "mc_truth_h0Zq11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_pt", &mc_truth_h0Zq21_pt, "mc_truth_h0Zq21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_pt", &mc_truth_h0Zq11_IS_pt, "mc_truth_h0Zq11_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_pt", &mc_truth_h0Zq21_IS_pt, "mc_truth_h0Zq21_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_pt", &mc_truth_h0Zl12_pt, "mc_truth_h0Zl12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_pt", &mc_truth_h0Zl22_pt, "mc_truth_h0Zl22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_pt", &mc_truth_h0Ztau12_pt, "mc_truth_h0Ztau12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_pt", &mc_truth_h0Ztau22_pt, "mc_truth_h0Ztau22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_pt", &mc_truth_h0Ztaul12_pt, "mc_truth_h0Ztaul12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_pt", &mc_truth_h0Ztaul22_pt, "mc_truth_h0Ztaul22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_pt", &mc_truth_h0Ztaunu12_pt, "mc_truth_h0Ztaunu12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_pt", &mc_truth_h0Ztaunu22_pt, "mc_truth_h0Ztaunu22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_pt", &mc_truth_h0Ztaunutau12_pt, "mc_truth_h0Ztaunutau12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_pt", &mc_truth_h0Ztaunutau22_pt, "mc_truth_h0Ztaunutau22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_pt", &mc_truth_h0Zq12_pt, "mc_truth_h0Zq12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_pt", &mc_truth_h0Zq22_pt, "mc_truth_h0Zq22_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_pt", &mc_truth_h0Zq12_IS_pt, "mc_truth_h0Zq12_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_pt", &mc_truth_h0Zq22_IS_pt, "mc_truth_h0Zq22_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_pt", &mc_truth_h0Znu11_pt, "mc_truth_h0Znu11_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_pt", &mc_truth_h0Znu21_pt, "mc_truth_h0Znu21_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_pt", &mc_truth_h0Znu12_pt, "mc_truth_h0Znu12_pt/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_pt", &mc_truth_h0Znu22_pt, "mc_truth_h0Znu22_pt/F", buffersize);

	tree->Branch("mc_truth_h0tau1_pt", &mc_truth_h0tau1_pt, "mc_truth_h0tau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0tau2_pt", &mc_truth_h0tau2_pt, "mc_truth_h0tau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taul1_pt", &mc_truth_h0taul1_pt, "mc_truth_h0taul1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_pt", &mc_truth_h0taunutau1_pt, "mc_truth_h0taunutau1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_pt", &mc_truth_h0taunu1_pt, "mc_truth_h0taunu1_pt/F", buffersize);
	tree->Branch("mc_truth_h0taul2_pt", &mc_truth_h0taul2_pt, "mc_truth_h0taul2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_pt", &mc_truth_h0taunutau2_pt, "mc_truth_h0taunutau2_pt/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_pt", &mc_truth_h0taunu2_pt, "mc_truth_h0taunu2_pt/F", buffersize);

	tree->Branch("mc_truth_h0b1_pt", &mc_truth_h0b1_pt, "mc_truth_h0b1_pt/F", buffersize);
	tree->Branch("mc_truth_h0b2_pt", &mc_truth_h0b2_pt, "mc_truth_h0b2_pt/F", buffersize);
	tree->Branch("mc_truth_h0b1_IS_pt", &mc_truth_h0b1_IS_pt, "mc_truth_h0b1_IS_pt/F", buffersize);
	tree->Branch("mc_truth_h0b2_IS_pt", &mc_truth_h0b2_IS_pt, "mc_truth_h0b2_IS_pt/F", buffersize);
	
	tree->Branch("mc_truth_t_pt", &mc_truth_t_pt, "mc_truth_t_pt/F", buffersize);
	tree->Branch("mc_truth_tb_pt", &mc_truth_tb_pt, "mc_truth_tb_pt/F", buffersize);
	tree->Branch("mc_truth_tb_IS_pt", &mc_truth_tb_IS_pt, "mc_truth_tb_IS_pt/F", buffersize);

	tree->Branch("mc_truth_tW_pt", &mc_truth_tW_pt, "mc_truth_tW_pt/F", buffersize);
	tree->Branch("mc_truth_tWnu_pt", &mc_truth_tWnu_pt, "mc_truth_tWnu_pt/F", buffersize);
	tree->Branch("mc_truth_tWnutau_pt", &mc_truth_tWnutau_pt, "mc_truth_tWnutau_pt/F", buffersize);
	tree->Branch("mc_truth_tWl_pt", &mc_truth_tWl_pt, "mc_truth_tWl_pt/F", buffersize);
	tree->Branch("mc_truth_tWtau_pt", &mc_truth_tWtau_pt, "mc_truth_tWtau_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_pt", &mc_truth_tWtaunu_pt, "mc_truth_tWtaunu_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_pt", &mc_truth_tWtaunutau_pt, "mc_truth_tWtaunutau_pt/F", buffersize);
	tree->Branch("mc_truth_tWtaul_pt", &mc_truth_tWtaul_pt, "mc_truth_tWtaul_pt/F", buffersize);
	tree->Branch("mc_truth_tWq1_pt", &mc_truth_tWq1_pt, "mc_truth_tWq1_pt/F", buffersize);
	tree->Branch("mc_truth_tWq2_pt", &mc_truth_tWq2_pt, "mc_truth_tWq2_pt/F", buffersize);
	tree->Branch("mc_truth_tWq1_IS_pt", &mc_truth_tWq1_IS_pt, "mc_truth_tWq1_IS_pt/F", buffersize);
	tree->Branch("mc_truth_tWq2_IS_pt", &mc_truth_tWq2_IS_pt, "mc_truth_tWq2_IS_pt/F", buffersize);

	tree->Branch("mc_truth_j1_pt", &mc_truth_j1_pt, "mc_truth_j1_pt/F", buffersize);
	tree->Branch("mc_truth_j2_pt", &mc_truth_j2_pt, "mc_truth_j2_pt/F", buffersize);
	tree->Branch("mc_truth_j3_pt", &mc_truth_j3_pt, "mc_truth_j3_pt/F", buffersize);
		
	tree->Branch("mc_truth_h0_eta", &mc_truth_h0_eta, "mc_truth_h0_eta/F", buffersize);

	tree->Branch("mc_truth_h0W1_eta", &mc_truth_h0W1_eta, "mc_truth_h0W1_eta/F", buffersize);
	tree->Branch("mc_truth_h0W2_eta", &mc_truth_h0W2_eta, "mc_truth_h0W2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_eta", &mc_truth_h0Wl1_eta, "mc_truth_h0Wl1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_eta", &mc_truth_h0Wnu1_eta, "mc_truth_h0Wnu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_eta", &mc_truth_h0Wtau1_eta, "mc_truth_h0Wtau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_eta", &mc_truth_h0Wnutau1_eta, "mc_truth_h0Wnutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_eta", &mc_truth_h0Wtaul1_eta, "mc_truth_h0Wtaul1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_eta", &mc_truth_h0Wtaunu1_eta, "mc_truth_h0Wtaunu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_eta", &mc_truth_h0Wtaunutau1_eta, "mc_truth_h0Wtaunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_eta", &mc_truth_h0Wl2_eta, "mc_truth_h0Wl2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_eta", &mc_truth_h0Wnu2_eta, "mc_truth_h0Wnu2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_eta", &mc_truth_h0Wtau2_eta, "mc_truth_h0Wtau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_eta", &mc_truth_h0Wnutau2_eta, "mc_truth_h0Wnutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_eta", &mc_truth_h0Wtaul2_eta, "mc_truth_h0Wtaul2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_eta", &mc_truth_h0Wtaunu2_eta, "mc_truth_h0Wtaunu2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_eta", &mc_truth_h0Wtaunutau2_eta, "mc_truth_h0Wtaunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_eta", &mc_truth_h0Wq11_eta, "mc_truth_h0Wq11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_eta", &mc_truth_h0Wq21_eta, "mc_truth_h0Wq21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_eta", &mc_truth_h0Wq12_eta, "mc_truth_h0Wq12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_eta", &mc_truth_h0Wq22_eta, "mc_truth_h0Wq22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_eta", &mc_truth_h0Wq11_IS_eta, "mc_truth_h0Wq11_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_eta", &mc_truth_h0Wq21_IS_eta, "mc_truth_h0Wq21_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_eta", &mc_truth_h0Wq12_IS_eta, "mc_truth_h0Wq12_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_eta", &mc_truth_h0Wq22_IS_eta, "mc_truth_h0Wq22_IS_eta/F", buffersize);

	tree->Branch("mc_truth_h0Z1_eta", &mc_truth_h0Z1_eta, "mc_truth_h0Z1_eta/F", buffersize);
	tree->Branch("mc_truth_h0Z2_eta", &mc_truth_h0Z2_eta, "mc_truth_h0Z2_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_eta", &mc_truth_h0Zl11_eta, "mc_truth_h0Zl11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_eta", &mc_truth_h0Zl21_eta, "mc_truth_h0Zl21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_eta", &mc_truth_h0Ztau11_eta, "mc_truth_h0Ztau11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_eta", &mc_truth_h0Ztau21_eta, "mc_truth_h0Ztau21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_eta", &mc_truth_h0Ztaul11_eta, "mc_truth_h0Ztaul11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_eta", &mc_truth_h0Ztaul21_eta, "mc_truth_h0Ztaul21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_eta", &mc_truth_h0Ztaunu11_eta, "mc_truth_h0Ztaunu11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_eta", &mc_truth_h0Ztaunu21_eta, "mc_truth_h0Ztaunu21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_eta", &mc_truth_h0Ztaunutau11_eta, "mc_truth_h0Ztaunutau11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_eta", &mc_truth_h0Ztaunutau21_eta, "mc_truth_h0Ztaunutau21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_eta", &mc_truth_h0Zq11_eta, "mc_truth_h0Zq11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_eta", &mc_truth_h0Zq21_eta, "mc_truth_h0Zq21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_eta", &mc_truth_h0Zq11_IS_eta, "mc_truth_h0Zq11_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_eta", &mc_truth_h0Zq21_IS_eta, "mc_truth_h0Zq21_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_eta", &mc_truth_h0Zl12_eta, "mc_truth_h0Zl12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_eta", &mc_truth_h0Zl22_eta, "mc_truth_h0Zl22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_eta", &mc_truth_h0Ztau12_eta, "mc_truth_h0Ztau12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_eta", &mc_truth_h0Ztau22_eta, "mc_truth_h0Ztau22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_eta", &mc_truth_h0Ztaul12_eta, "mc_truth_h0Ztaul12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_eta", &mc_truth_h0Ztaul22_eta, "mc_truth_h0Ztaul22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_eta", &mc_truth_h0Ztaunu12_eta, "mc_truth_h0Ztaunu12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_eta", &mc_truth_h0Ztaunu22_eta, "mc_truth_h0Ztaunu22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_eta", &mc_truth_h0Ztaunutau12_eta, "mc_truth_h0Ztaunutau12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_eta", &mc_truth_h0Ztaunutau22_eta, "mc_truth_h0Ztaunutau22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_eta", &mc_truth_h0Zq12_eta, "mc_truth_h0Zq12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_eta", &mc_truth_h0Zq22_eta, "mc_truth_h0Zq22_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_eta", &mc_truth_h0Zq12_IS_eta, "mc_truth_h0Zq12_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_eta", &mc_truth_h0Zq22_IS_eta, "mc_truth_h0Zq22_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_eta", &mc_truth_h0Znu11_eta, "mc_truth_h0Znu11_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_eta", &mc_truth_h0Znu21_eta, "mc_truth_h0Znu21_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_eta", &mc_truth_h0Znu12_eta, "mc_truth_h0Znu12_eta/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_eta", &mc_truth_h0Znu22_eta, "mc_truth_h0Znu22_eta/F", buffersize);

	tree->Branch("mc_truth_h0tau1_eta", &mc_truth_h0tau1_eta, "mc_truth_h0tau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0tau2_eta", &mc_truth_h0tau2_eta, "mc_truth_h0tau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taul1_eta", &mc_truth_h0taul1_eta, "mc_truth_h0taul1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_eta", &mc_truth_h0taunutau1_eta, "mc_truth_h0taunutau1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_eta", &mc_truth_h0taunu1_eta, "mc_truth_h0taunu1_eta/F", buffersize);
	tree->Branch("mc_truth_h0taul2_eta", &mc_truth_h0taul2_eta, "mc_truth_h0taul2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_eta", &mc_truth_h0taunutau2_eta, "mc_truth_h0taunutau2_eta/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_eta", &mc_truth_h0taunu2_eta, "mc_truth_h0taunu2_eta/F", buffersize);

	tree->Branch("mc_truth_h0b1_eta", &mc_truth_h0b1_eta, "mc_truth_h0b1_eta/F", buffersize);
	tree->Branch("mc_truth_h0b2_eta", &mc_truth_h0b2_eta, "mc_truth_h0b2_eta/F", buffersize);
	tree->Branch("mc_truth_h0b1_IS_eta", &mc_truth_h0b1_IS_eta, "mc_truth_h0b1_IS_eta/F", buffersize);
	tree->Branch("mc_truth_h0b2_IS_eta", &mc_truth_h0b2_IS_eta, "mc_truth_h0b2_IS_eta/F", buffersize);
	
	tree->Branch("mc_truth_t_eta", &mc_truth_t_eta, "mc_truth_t_eta/F", buffersize);
	tree->Branch("mc_truth_tb_eta", &mc_truth_tb_eta, "mc_truth_tb_eta/F", buffersize);
	tree->Branch("mc_truth_tb_IS_eta", &mc_truth_tb_IS_eta, "mc_truth_tb_IS_eta/F", buffersize);

	tree->Branch("mc_truth_tW_eta", &mc_truth_tW_eta, "mc_truth_tW_eta/F", buffersize);
	tree->Branch("mc_truth_tWnu_eta", &mc_truth_tWnu_eta, "mc_truth_tWnu_eta/F", buffersize);
	tree->Branch("mc_truth_tWnutau_eta", &mc_truth_tWnutau_eta, "mc_truth_tWnutau_eta/F", buffersize);
	tree->Branch("mc_truth_tWl_eta", &mc_truth_tWl_eta, "mc_truth_tWl_eta/F", buffersize);
	tree->Branch("mc_truth_tWtau_eta", &mc_truth_tWtau_eta, "mc_truth_tWtau_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_eta", &mc_truth_tWtaunu_eta, "mc_truth_tWtaunu_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_eta", &mc_truth_tWtaunutau_eta, "mc_truth_tWtaunutau_eta/F", buffersize);
	tree->Branch("mc_truth_tWtaul_eta", &mc_truth_tWtaul_eta, "mc_truth_tWtaul_eta/F", buffersize);
	tree->Branch("mc_truth_tWq1_eta", &mc_truth_tWq1_eta, "mc_truth_tWq1_eta/F", buffersize);
	tree->Branch("mc_truth_tWq2_eta", &mc_truth_tWq2_eta, "mc_truth_tWq2_eta/F", buffersize);
	tree->Branch("mc_truth_tWq1_IS_eta", &mc_truth_tWq1_IS_eta, "mc_truth_tWq1_IS_eta/F", buffersize);
	tree->Branch("mc_truth_tWq2_IS_eta", &mc_truth_tWq2_IS_eta, "mc_truth_tWq2_IS_eta/F", buffersize);

	tree->Branch("mc_truth_j1_eta", &mc_truth_j1_eta, "mc_truth_j1_eta/F", buffersize);
	tree->Branch("mc_truth_j2_eta", &mc_truth_j2_eta, "mc_truth_j2_eta/F", buffersize);
	tree->Branch("mc_truth_j3_eta", &mc_truth_j3_eta, "mc_truth_j3_eta/F", buffersize);
		
	tree->Branch("mc_truth_h0_phi", &mc_truth_h0_phi, "mc_truth_h0_phi/F", buffersize);

	tree->Branch("mc_truth_h0W1_phi", &mc_truth_h0W1_phi, "mc_truth_h0W1_phi/F", buffersize);
	tree->Branch("mc_truth_h0W2_phi", &mc_truth_h0W2_phi, "mc_truth_h0W2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_phi", &mc_truth_h0Wl1_phi, "mc_truth_h0Wl1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_phi", &mc_truth_h0Wnu1_phi, "mc_truth_h0Wnu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_phi", &mc_truth_h0Wtau1_phi, "mc_truth_h0Wtau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_phi", &mc_truth_h0Wnutau1_phi, "mc_truth_h0Wnutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_phi", &mc_truth_h0Wtaul1_phi, "mc_truth_h0Wtaul1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_phi", &mc_truth_h0Wtaunu1_phi, "mc_truth_h0Wtaunu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_phi", &mc_truth_h0Wtaunutau1_phi, "mc_truth_h0Wtaunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_phi", &mc_truth_h0Wl2_phi, "mc_truth_h0Wl2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_phi", &mc_truth_h0Wnu2_phi, "mc_truth_h0Wnu2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_phi", &mc_truth_h0Wtau2_phi, "mc_truth_h0Wtau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_phi", &mc_truth_h0Wnutau2_phi, "mc_truth_h0Wnutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_phi", &mc_truth_h0Wtaul2_phi, "mc_truth_h0Wtaul2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_phi", &mc_truth_h0Wtaunu2_phi, "mc_truth_h0Wtaunu2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_phi", &mc_truth_h0Wtaunutau2_phi, "mc_truth_h0Wtaunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_phi", &mc_truth_h0Wq11_phi, "mc_truth_h0Wq11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_phi", &mc_truth_h0Wq21_phi, "mc_truth_h0Wq21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_phi", &mc_truth_h0Wq12_phi, "mc_truth_h0Wq12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_phi", &mc_truth_h0Wq22_phi, "mc_truth_h0Wq22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_phi", &mc_truth_h0Wq11_IS_phi, "mc_truth_h0Wq11_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_phi", &mc_truth_h0Wq21_IS_phi, "mc_truth_h0Wq21_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_phi", &mc_truth_h0Wq12_IS_phi, "mc_truth_h0Wq12_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_phi", &mc_truth_h0Wq22_IS_phi, "mc_truth_h0Wq22_IS_phi/F", buffersize);

	tree->Branch("mc_truth_h0Z1_phi", &mc_truth_h0Z1_phi, "mc_truth_h0Z1_phi/F", buffersize);
	tree->Branch("mc_truth_h0Z2_phi", &mc_truth_h0Z2_phi, "mc_truth_h0Z2_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_phi", &mc_truth_h0Zl11_phi, "mc_truth_h0Zl11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_phi", &mc_truth_h0Zl21_phi, "mc_truth_h0Zl21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_phi", &mc_truth_h0Ztau11_phi, "mc_truth_h0Ztau11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_phi", &mc_truth_h0Ztau21_phi, "mc_truth_h0Ztau21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_phi", &mc_truth_h0Ztaul11_phi, "mc_truth_h0Ztaul11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_phi", &mc_truth_h0Ztaul21_phi, "mc_truth_h0Ztaul21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_phi", &mc_truth_h0Ztaunu11_phi, "mc_truth_h0Ztaunu11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_phi", &mc_truth_h0Ztaunu21_phi, "mc_truth_h0Ztaunu21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_phi", &mc_truth_h0Ztaunutau11_phi, "mc_truth_h0Ztaunutau11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_phi", &mc_truth_h0Ztaunutau21_phi, "mc_truth_h0Ztaunutau21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_phi", &mc_truth_h0Zq11_phi, "mc_truth_h0Zq11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_phi", &mc_truth_h0Zq21_phi, "mc_truth_h0Zq21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_phi", &mc_truth_h0Zq11_IS_phi, "mc_truth_h0Zq11_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_phi", &mc_truth_h0Zq21_IS_phi, "mc_truth_h0Zq21_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_phi", &mc_truth_h0Zl12_phi, "mc_truth_h0Zl12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_phi", &mc_truth_h0Zl22_phi, "mc_truth_h0Zl22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_phi", &mc_truth_h0Ztau12_phi, "mc_truth_h0Ztau12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_phi", &mc_truth_h0Ztau22_phi, "mc_truth_h0Ztau22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_phi", &mc_truth_h0Ztaul12_phi, "mc_truth_h0Ztaul12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_phi", &mc_truth_h0Ztaul22_phi, "mc_truth_h0Ztaul22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_phi", &mc_truth_h0Ztaunu12_phi, "mc_truth_h0Ztaunu12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_phi", &mc_truth_h0Ztaunu22_phi, "mc_truth_h0Ztaunu22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_phi", &mc_truth_h0Ztaunutau12_phi, "mc_truth_h0Ztaunutau12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_phi", &mc_truth_h0Ztaunutau22_phi, "mc_truth_h0Ztaunutau22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_phi", &mc_truth_h0Zq12_phi, "mc_truth_h0Zq12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_phi", &mc_truth_h0Zq22_phi, "mc_truth_h0Zq22_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_phi", &mc_truth_h0Zq12_IS_phi, "mc_truth_h0Zq12_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_phi", &mc_truth_h0Zq22_IS_phi, "mc_truth_h0Zq22_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_phi", &mc_truth_h0Znu11_phi, "mc_truth_h0Znu11_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_phi", &mc_truth_h0Znu21_phi, "mc_truth_h0Znu21_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_phi", &mc_truth_h0Znu12_phi, "mc_truth_h0Znu12_phi/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_phi", &mc_truth_h0Znu22_phi, "mc_truth_h0Znu22_phi/F", buffersize);

	tree->Branch("mc_truth_h0tau1_phi", &mc_truth_h0tau1_phi, "mc_truth_h0tau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0tau2_phi", &mc_truth_h0tau2_phi, "mc_truth_h0tau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taul1_phi", &mc_truth_h0taul1_phi, "mc_truth_h0taul1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_phi", &mc_truth_h0taunutau1_phi, "mc_truth_h0taunutau1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_phi", &mc_truth_h0taunu1_phi, "mc_truth_h0taunu1_phi/F", buffersize);
	tree->Branch("mc_truth_h0taul2_phi", &mc_truth_h0taul2_phi, "mc_truth_h0taul2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_phi", &mc_truth_h0taunutau2_phi, "mc_truth_h0taunutau2_phi/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_phi", &mc_truth_h0taunu2_phi, "mc_truth_h0taunu2_phi/F", buffersize);

	tree->Branch("mc_truth_h0b1_phi", &mc_truth_h0b1_phi, "mc_truth_h0b1_phi/F", buffersize);
	tree->Branch("mc_truth_h0b2_phi", &mc_truth_h0b2_phi, "mc_truth_h0b2_phi/F", buffersize);
	tree->Branch("mc_truth_h0b1_IS_phi", &mc_truth_h0b1_IS_phi, "mc_truth_h0b1_IS_phi/F", buffersize);
	tree->Branch("mc_truth_h0b2_IS_phi", &mc_truth_h0b2_IS_phi, "mc_truth_h0b2_IS_phi/F", buffersize);
	
	tree->Branch("mc_truth_t_phi", &mc_truth_t_phi, "mc_truth_t_phi/F", buffersize);
	tree->Branch("mc_truth_tb_phi", &mc_truth_tb_phi, "mc_truth_tb_phi/F", buffersize);
	tree->Branch("mc_truth_tb_IS_phi", &mc_truth_tb_IS_phi, "mc_truth_tb_IS_phi/F", buffersize);

	tree->Branch("mc_truth_tW_phi", &mc_truth_tW_phi, "mc_truth_tW_phi/F", buffersize);
	tree->Branch("mc_truth_tWnu_phi", &mc_truth_tWnu_phi, "mc_truth_tWnu_phi/F", buffersize);
	tree->Branch("mc_truth_tWnutau_phi", &mc_truth_tWnutau_phi, "mc_truth_tWnutau_phi/F", buffersize);
	tree->Branch("mc_truth_tWl_phi", &mc_truth_tWl_phi, "mc_truth_tWl_phi/F", buffersize);
	tree->Branch("mc_truth_tWtau_phi", &mc_truth_tWtau_phi, "mc_truth_tWtau_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_phi", &mc_truth_tWtaunu_phi, "mc_truth_tWtaunu_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_phi", &mc_truth_tWtaunutau_phi, "mc_truth_tWtaunutau_phi/F", buffersize);
	tree->Branch("mc_truth_tWtaul_phi", &mc_truth_tWtaul_phi, "mc_truth_tWtaul_phi/F", buffersize);
	tree->Branch("mc_truth_tWq1_phi", &mc_truth_tWq1_phi, "mc_truth_tWq1_phi/F", buffersize);
	tree->Branch("mc_truth_tWq2_phi", &mc_truth_tWq2_phi, "mc_truth_tWq2_phi/F", buffersize);
	tree->Branch("mc_truth_tWq1_IS_phi", &mc_truth_tWq1_IS_phi, "mc_truth_tWq1_IS_phi/F", buffersize);
	tree->Branch("mc_truth_tWq2_IS_phi", &mc_truth_tWq2_IS_phi, "mc_truth_tWq2_IS_phi/F", buffersize);

	tree->Branch("mc_truth_j1_phi", &mc_truth_j1_phi, "mc_truth_j1_phi/F", buffersize);
	tree->Branch("mc_truth_j2_phi", &mc_truth_j2_phi, "mc_truth_j2_phi/F", buffersize);
	tree->Branch("mc_truth_j3_phi", &mc_truth_j3_phi, "mc_truth_j3_phi/F", buffersize);
	
	tree->Branch("mc_truth_h0_E", &mc_truth_h0_E, "mc_truth_h0_E/F", buffersize);

	tree->Branch("mc_truth_h0W1_E", &mc_truth_h0W1_E, "mc_truth_h0W1_E/F", buffersize);
	tree->Branch("mc_truth_h0W2_E", &mc_truth_h0W2_E, "mc_truth_h0W2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wl1_E", &mc_truth_h0Wl1_E, "mc_truth_h0Wl1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnu1_E", &mc_truth_h0Wnu1_E, "mc_truth_h0Wnu1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtau1_E", &mc_truth_h0Wtau1_E, "mc_truth_h0Wtau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_E", &mc_truth_h0Wnutau1_E, "mc_truth_h0Wnutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_E", &mc_truth_h0Wtaul1_E, "mc_truth_h0Wtaul1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_E", &mc_truth_h0Wtaunu1_E, "mc_truth_h0Wtaunu1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_E", &mc_truth_h0Wtaunutau1_E, "mc_truth_h0Wtaunutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0Wl2_E", &mc_truth_h0Wl2_E, "mc_truth_h0Wl2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnu2_E", &mc_truth_h0Wnu2_E, "mc_truth_h0Wnu2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtau2_E", &mc_truth_h0Wtau2_E, "mc_truth_h0Wtau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_E", &mc_truth_h0Wnutau2_E, "mc_truth_h0Wnutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_E", &mc_truth_h0Wtaul2_E, "mc_truth_h0Wtaul2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_E", &mc_truth_h0Wtaunu2_E, "mc_truth_h0Wtaunu2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_E", &mc_truth_h0Wtaunutau2_E, "mc_truth_h0Wtaunutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_E", &mc_truth_h0Wq11_E, "mc_truth_h0Wq11_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_E", &mc_truth_h0Wq21_E, "mc_truth_h0Wq21_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_E", &mc_truth_h0Wq12_E, "mc_truth_h0Wq12_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_E", &mc_truth_h0Wq22_E, "mc_truth_h0Wq22_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_E", &mc_truth_h0Wq11_IS_E, "mc_truth_h0Wq11_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_E", &mc_truth_h0Wq21_IS_E, "mc_truth_h0Wq21_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_E", &mc_truth_h0Wq12_IS_E, "mc_truth_h0Wq12_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_E", &mc_truth_h0Wq22_IS_E, "mc_truth_h0Wq22_IS_E/F", buffersize);

	tree->Branch("mc_truth_h0Z1_E", &mc_truth_h0Z1_E, "mc_truth_h0Z1_E/F", buffersize);
	tree->Branch("mc_truth_h0Z2_E", &mc_truth_h0Z2_E, "mc_truth_h0Z2_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl11_E", &mc_truth_h0Zl11_E, "mc_truth_h0Zl11_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl21_E", &mc_truth_h0Zl21_E, "mc_truth_h0Zl21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau11_E", &mc_truth_h0Ztau11_E, "mc_truth_h0Ztau11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau21_E", &mc_truth_h0Ztau21_E, "mc_truth_h0Ztau21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_E", &mc_truth_h0Ztaul11_E, "mc_truth_h0Ztaul11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_E", &mc_truth_h0Ztaul21_E, "mc_truth_h0Ztaul21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_E", &mc_truth_h0Ztaunu11_E, "mc_truth_h0Ztaunu11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_E", &mc_truth_h0Ztaunu21_E, "mc_truth_h0Ztaunu21_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_E", &mc_truth_h0Ztaunutau11_E, "mc_truth_h0Ztaunutau11_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_E", &mc_truth_h0Ztaunutau21_E, "mc_truth_h0Ztaunutau21_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_E", &mc_truth_h0Zq11_E, "mc_truth_h0Zq11_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_E", &mc_truth_h0Zq21_E, "mc_truth_h0Zq21_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_E", &mc_truth_h0Zq11_IS_E, "mc_truth_h0Zq11_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_E", &mc_truth_h0Zq21_IS_E, "mc_truth_h0Zq21_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl12_E", &mc_truth_h0Zl12_E, "mc_truth_h0Zl12_E/F", buffersize);
	tree->Branch("mc_truth_h0Zl22_E", &mc_truth_h0Zl22_E, "mc_truth_h0Zl22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau12_E", &mc_truth_h0Ztau12_E, "mc_truth_h0Ztau12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztau22_E", &mc_truth_h0Ztau22_E, "mc_truth_h0Ztau22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_E", &mc_truth_h0Ztaul12_E, "mc_truth_h0Ztaul12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_E", &mc_truth_h0Ztaul22_E, "mc_truth_h0Ztaul22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_E", &mc_truth_h0Ztaunu12_E, "mc_truth_h0Ztaunu12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_E", &mc_truth_h0Ztaunu22_E, "mc_truth_h0Ztaunu22_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_E", &mc_truth_h0Ztaunutau12_E, "mc_truth_h0Ztaunutau12_E/F", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_E", &mc_truth_h0Ztaunutau22_E, "mc_truth_h0Ztaunutau22_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_E", &mc_truth_h0Zq12_E, "mc_truth_h0Zq12_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_E", &mc_truth_h0Zq22_E, "mc_truth_h0Zq22_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_E", &mc_truth_h0Zq12_IS_E, "mc_truth_h0Zq12_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_E", &mc_truth_h0Zq22_IS_E, "mc_truth_h0Zq22_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu11_E", &mc_truth_h0Znu11_E, "mc_truth_h0Znu11_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu21_E", &mc_truth_h0Znu21_E, "mc_truth_h0Znu21_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu12_E", &mc_truth_h0Znu12_E, "mc_truth_h0Znu12_E/F", buffersize);
	tree->Branch("mc_truth_h0Znu22_E", &mc_truth_h0Znu22_E, "mc_truth_h0Znu22_E/F", buffersize);

	tree->Branch("mc_truth_h0tau1_E", &mc_truth_h0tau1_E, "mc_truth_h0tau1_E/F", buffersize);
	tree->Branch("mc_truth_h0tau2_E", &mc_truth_h0tau2_E, "mc_truth_h0tau2_E/F", buffersize);
	tree->Branch("mc_truth_h0taul1_E", &mc_truth_h0taul1_E, "mc_truth_h0taul1_E/F", buffersize);
	tree->Branch("mc_truth_h0taunutau1_E", &mc_truth_h0taunutau1_E, "mc_truth_h0taunutau1_E/F", buffersize);
	tree->Branch("mc_truth_h0taunu1_E", &mc_truth_h0taunu1_E, "mc_truth_h0taunu1_E/F", buffersize);
	tree->Branch("mc_truth_h0taul2_E", &mc_truth_h0taul2_E, "mc_truth_h0taul2_E/F", buffersize);
	tree->Branch("mc_truth_h0taunutau2_E", &mc_truth_h0taunutau2_E, "mc_truth_h0taunutau2_E/F", buffersize);
	tree->Branch("mc_truth_h0taunu2_E", &mc_truth_h0taunu2_E, "mc_truth_h0taunu2_E/F", buffersize);

	tree->Branch("mc_truth_h0b1_E", &mc_truth_h0b1_E, "mc_truth_h0b1_E/F", buffersize);
	tree->Branch("mc_truth_h0b2_E", &mc_truth_h0b2_E, "mc_truth_h0b2_E/F", buffersize);
	tree->Branch("mc_truth_h0b1_IS_E", &mc_truth_h0b1_IS_E, "mc_truth_h0b1_IS_E/F", buffersize);
	tree->Branch("mc_truth_h0b2_IS_E", &mc_truth_h0b2_IS_E, "mc_truth_h0b2_IS_E/F", buffersize);
	
	tree->Branch("mc_truth_t_E", &mc_truth_t_E, "mc_truth_t_E/F", buffersize);
	tree->Branch("mc_truth_tb_E", &mc_truth_tb_E, "mc_truth_tb_E/F", buffersize);
	tree->Branch("mc_truth_tb_IS_E", &mc_truth_tb_IS_E, "mc_truth_tb_IS_E/F", buffersize);

	tree->Branch("mc_truth_tW_E", &mc_truth_tW_E, "mc_truth_tW_E/F", buffersize);
	tree->Branch("mc_truth_tWnu_E", &mc_truth_tWnu_E, "mc_truth_tWnu_E/F", buffersize);
	tree->Branch("mc_truth_tWnutau_E", &mc_truth_tWnutau_E, "mc_truth_tWnutau_E/F", buffersize);
	tree->Branch("mc_truth_tWl_E", &mc_truth_tWl_E, "mc_truth_tWl_E/F", buffersize);
	tree->Branch("mc_truth_tWtau_E", &mc_truth_tWtau_E, "mc_truth_tWtau_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunu_E", &mc_truth_tWtaunu_E, "mc_truth_tWtaunu_E/F", buffersize);
	tree->Branch("mc_truth_tWtaunutau_E", &mc_truth_tWtaunutau_E, "mc_truth_tWtaunutau_E/F", buffersize);
	tree->Branch("mc_truth_tWtaul_E", &mc_truth_tWtaul_E, "mc_truth_tWtaul_E/F", buffersize);
	tree->Branch("mc_truth_tWq1_E", &mc_truth_tWq1_E, "mc_truth_tWq1_E/F", buffersize);
	tree->Branch("mc_truth_tWq2_E", &mc_truth_tWq2_E, "mc_truth_tWq2_E/F", buffersize);
	tree->Branch("mc_truth_tWq1_IS_E", &mc_truth_tWq1_IS_E, "mc_truth_tWq1_IS_E/F", buffersize);
	tree->Branch("mc_truth_tWq2_IS_E", &mc_truth_tWq2_IS_E, "mc_truth_tWq2_IS_E/F", buffersize);

	tree->Branch("mc_truth_j1_E", &mc_truth_j1_E, "mc_truth_j1_E/F", buffersize);
	tree->Branch("mc_truth_j2_E", &mc_truth_j2_E, "mc_truth_j2_E/F", buffersize);
	tree->Branch("mc_truth_j3_E", &mc_truth_j3_E, "mc_truth_j3_E/F", buffersize);
				
	tree->Branch("mc_truth_h0_id", &mc_truth_h0_id, "mc_truth_h0_id/I", buffersize);

	tree->Branch("mc_truth_h0W1_id", &mc_truth_h0W1_id, "mc_truth_h0W1_id/I", buffersize);
	tree->Branch("mc_truth_h0W2_id", &mc_truth_h0W2_id, "mc_truth_h0W2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wl1_id", &mc_truth_h0Wl1_id, "mc_truth_h0Wl1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnu1_id", &mc_truth_h0Wnu1_id, "mc_truth_h0Wnu1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtau1_id", &mc_truth_h0Wtau1_id, "mc_truth_h0Wtau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_id", &mc_truth_h0Wnutau1_id, "mc_truth_h0Wnutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_id", &mc_truth_h0Wtaul1_id, "mc_truth_h0Wtaul1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_id", &mc_truth_h0Wtaunu1_id, "mc_truth_h0Wtaunu1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_id", &mc_truth_h0Wtaunutau1_id, "mc_truth_h0Wtaunutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0Wl2_id", &mc_truth_h0Wl2_id, "mc_truth_h0Wl2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnu2_id", &mc_truth_h0Wnu2_id, "mc_truth_h0Wnu2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtau2_id", &mc_truth_h0Wtau2_id, "mc_truth_h0Wtau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_id", &mc_truth_h0Wnutau2_id, "mc_truth_h0Wnutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_id", &mc_truth_h0Wtaul2_id, "mc_truth_h0Wtaul2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_id", &mc_truth_h0Wtaunu2_id, "mc_truth_h0Wtaunu2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_id", &mc_truth_h0Wtaunutau2_id, "mc_truth_h0Wtaunutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_id", &mc_truth_h0Wq11_id, "mc_truth_h0Wq11_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_id", &mc_truth_h0Wq21_id, "mc_truth_h0Wq21_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_id", &mc_truth_h0Wq12_id, "mc_truth_h0Wq12_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_id", &mc_truth_h0Wq22_id, "mc_truth_h0Wq22_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_id", &mc_truth_h0Wq11_IS_id, "mc_truth_h0Wq11_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_id", &mc_truth_h0Wq21_IS_id, "mc_truth_h0Wq21_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_id", &mc_truth_h0Wq12_IS_id, "mc_truth_h0Wq12_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_id", &mc_truth_h0Wq22_IS_id, "mc_truth_h0Wq22_IS_id/I", buffersize);

	tree->Branch("mc_truth_h0Z1_id", &mc_truth_h0Z1_id, "mc_truth_h0Z1_id/I", buffersize);
	tree->Branch("mc_truth_h0Z2_id", &mc_truth_h0Z2_id, "mc_truth_h0Z2_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl11_id", &mc_truth_h0Zl11_id, "mc_truth_h0Zl11_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl21_id", &mc_truth_h0Zl21_id, "mc_truth_h0Zl21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau11_id", &mc_truth_h0Ztau11_id, "mc_truth_h0Ztau11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau21_id", &mc_truth_h0Ztau21_id, "mc_truth_h0Ztau21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_id", &mc_truth_h0Ztaul11_id, "mc_truth_h0Ztaul11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_id", &mc_truth_h0Ztaul21_id, "mc_truth_h0Ztaul21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_id", &mc_truth_h0Ztaunu11_id, "mc_truth_h0Ztaunu11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_id", &mc_truth_h0Ztaunu21_id, "mc_truth_h0Ztaunu21_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_id", &mc_truth_h0Ztaunutau11_id, "mc_truth_h0Ztaunutau11_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_id", &mc_truth_h0Ztaunutau21_id, "mc_truth_h0Ztaunutau21_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_id", &mc_truth_h0Zq11_id, "mc_truth_h0Zq11_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_id", &mc_truth_h0Zq21_id, "mc_truth_h0Zq21_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_id", &mc_truth_h0Zq11_IS_id, "mc_truth_h0Zq11_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_id", &mc_truth_h0Zq21_IS_id, "mc_truth_h0Zq21_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl12_id", &mc_truth_h0Zl12_id, "mc_truth_h0Zl12_id/I", buffersize);
	tree->Branch("mc_truth_h0Zl22_id", &mc_truth_h0Zl22_id, "mc_truth_h0Zl22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau12_id", &mc_truth_h0Ztau12_id, "mc_truth_h0Ztau12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztau22_id", &mc_truth_h0Ztau22_id, "mc_truth_h0Ztau22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_id", &mc_truth_h0Ztaul12_id, "mc_truth_h0Ztaul12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_id", &mc_truth_h0Ztaul22_id, "mc_truth_h0Ztaul22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_id", &mc_truth_h0Ztaunu12_id, "mc_truth_h0Ztaunu12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_id", &mc_truth_h0Ztaunu22_id, "mc_truth_h0Ztaunu22_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_id", &mc_truth_h0Ztaunutau12_id, "mc_truth_h0Ztaunutau12_id/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_id", &mc_truth_h0Ztaunutau22_id, "mc_truth_h0Ztaunutau22_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_id", &mc_truth_h0Zq12_id, "mc_truth_h0Zq12_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_id", &mc_truth_h0Zq22_id, "mc_truth_h0Zq22_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_id", &mc_truth_h0Zq12_IS_id, "mc_truth_h0Zq12_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_id", &mc_truth_h0Zq22_IS_id, "mc_truth_h0Zq22_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu11_id", &mc_truth_h0Znu11_id, "mc_truth_h0Znu11_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu21_id", &mc_truth_h0Znu21_id, "mc_truth_h0Znu21_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu12_id", &mc_truth_h0Znu12_id, "mc_truth_h0Znu12_id/I", buffersize);
	tree->Branch("mc_truth_h0Znu22_id", &mc_truth_h0Znu22_id, "mc_truth_h0Znu22_id/I", buffersize);

	tree->Branch("mc_truth_h0tau1_id", &mc_truth_h0tau1_id, "mc_truth_h0tau1_id/I", buffersize);
	tree->Branch("mc_truth_h0tau2_id", &mc_truth_h0tau2_id, "mc_truth_h0tau2_id/I", buffersize);
	tree->Branch("mc_truth_h0taul1_id", &mc_truth_h0taul1_id, "mc_truth_h0taul1_id/I", buffersize);
	tree->Branch("mc_truth_h0taunutau1_id", &mc_truth_h0taunutau1_id, "mc_truth_h0taunutau1_id/I", buffersize);
	tree->Branch("mc_truth_h0taunu1_id", &mc_truth_h0taunu1_id, "mc_truth_h0taunu1_id/I", buffersize);
	tree->Branch("mc_truth_h0taul2_id", &mc_truth_h0taul2_id, "mc_truth_h0taul2_id/I", buffersize);
	tree->Branch("mc_truth_h0taunutau2_id", &mc_truth_h0taunutau2_id, "mc_truth_h0taunutau2_id/I", buffersize);
	tree->Branch("mc_truth_h0taunu2_id", &mc_truth_h0taunu2_id, "mc_truth_h0taunu2_id/I", buffersize);

	tree->Branch("mc_truth_h0b1_id", &mc_truth_h0b1_id, "mc_truth_h0b1_id/I", buffersize);
	tree->Branch("mc_truth_h0b2_id", &mc_truth_h0b2_id, "mc_truth_h0b2_id/I", buffersize);
	tree->Branch("mc_truth_h0b1_IS_id", &mc_truth_h0b1_IS_id, "mc_truth_h0b1_IS_id/I", buffersize);
	tree->Branch("mc_truth_h0b2_IS_id", &mc_truth_h0b2_IS_id, "mc_truth_h0b2_IS_id/I", buffersize);
	
	tree->Branch("mc_truth_t_id", &mc_truth_t_id, "mc_truth_t_id/I", buffersize);
	tree->Branch("mc_truth_tb_id", &mc_truth_tb_id, "mc_truth_tb_id/I", buffersize);
	tree->Branch("mc_truth_tb_IS_id", &mc_truth_tb_IS_id, "mc_truth_tb_IS_id/I", buffersize);

	tree->Branch("mc_truth_tW_id", &mc_truth_tW_id, "mc_truth_tW_id/I", buffersize);
	tree->Branch("mc_truth_tWnu_id", &mc_truth_tWnu_id, "mc_truth_tWnu_id/I", buffersize);
	tree->Branch("mc_truth_tWnutau_id", &mc_truth_tWnutau_id, "mc_truth_tWnutau_id/I", buffersize);
	tree->Branch("mc_truth_tWl_id", &mc_truth_tWl_id, "mc_truth_tWl_id/I", buffersize);
	tree->Branch("mc_truth_tWtau_id", &mc_truth_tWtau_id, "mc_truth_tWtau_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunu_id", &mc_truth_tWtaunu_id, "mc_truth_tWtaunu_id/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau_id", &mc_truth_tWtaunutau_id, "mc_truth_tWtaunutau_id/I", buffersize);
	tree->Branch("mc_truth_tWtaul_id", &mc_truth_tWtaul_id, "mc_truth_tWtaul_id/I", buffersize);
	tree->Branch("mc_truth_tWq1_id", &mc_truth_tWq1_id, "mc_truth_tWq1_id/I", buffersize);
	tree->Branch("mc_truth_tWq2_id", &mc_truth_tWq2_id, "mc_truth_tWq2_id/I", buffersize);
	tree->Branch("mc_truth_tWq1_IS_id", &mc_truth_tWq1_IS_id, "mc_truth_tWq1_IS_id/I", buffersize);
	tree->Branch("mc_truth_tWq2_IS_id", &mc_truth_tWq2_IS_id, "mc_truth_tWq2_IS_id/I", buffersize);

	tree->Branch("mc_truth_j1_id", &mc_truth_j1_id, "mc_truth_j1_id/I", buffersize);
	tree->Branch("mc_truth_j2_id", &mc_truth_j2_id, "mc_truth_j2_id/I", buffersize);
	tree->Branch("mc_truth_j3_id", &mc_truth_j3_id, "mc_truth_j3_id/I", buffersize);

	tree->Branch("mc_truth_h0_status", &mc_truth_h0_status, "mc_truth_h0_status/I", buffersize);

	tree->Branch("mc_truth_h0W1_status", &mc_truth_h0W1_status, "mc_truth_h0W1_status/I", buffersize);
	tree->Branch("mc_truth_h0W2_status", &mc_truth_h0W2_status, "mc_truth_h0W2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wl1_status", &mc_truth_h0Wl1_status, "mc_truth_h0Wl1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnu1_status", &mc_truth_h0Wnu1_status, "mc_truth_h0Wnu1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtau1_status", &mc_truth_h0Wtau1_status, "mc_truth_h0Wtau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau1_status", &mc_truth_h0Wnutau1_status, "mc_truth_h0Wnutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul1_status", &mc_truth_h0Wtaul1_status, "mc_truth_h0Wtaul1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu1_status", &mc_truth_h0Wtaunu1_status, "mc_truth_h0Wtaunu1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau1_status", &mc_truth_h0Wtaunutau1_status, "mc_truth_h0Wtaunutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0Wl2_status", &mc_truth_h0Wl2_status, "mc_truth_h0Wl2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnu2_status", &mc_truth_h0Wnu2_status, "mc_truth_h0Wnu2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtau2_status", &mc_truth_h0Wtau2_status, "mc_truth_h0Wtau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wnutau2_status", &mc_truth_h0Wnutau2_status, "mc_truth_h0Wnutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaul2_status", &mc_truth_h0Wtaul2_status, "mc_truth_h0Wtaul2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunu2_status", &mc_truth_h0Wtaunu2_status, "mc_truth_h0Wtaunu2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wtaunutau2_status", &mc_truth_h0Wtaunutau2_status, "mc_truth_h0Wtaunutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_status", &mc_truth_h0Wq11_status, "mc_truth_h0Wq11_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_status", &mc_truth_h0Wq21_status, "mc_truth_h0Wq21_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_status", &mc_truth_h0Wq12_status, "mc_truth_h0Wq12_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_status", &mc_truth_h0Wq22_status, "mc_truth_h0Wq22_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq11_IS_status", &mc_truth_h0Wq11_IS_status, "mc_truth_h0Wq11_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq21_IS_status", &mc_truth_h0Wq21_IS_status, "mc_truth_h0Wq21_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq12_IS_status", &mc_truth_h0Wq12_IS_status, "mc_truth_h0Wq12_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Wq22_IS_status", &mc_truth_h0Wq22_IS_status, "mc_truth_h0Wq22_IS_status/I", buffersize);

	tree->Branch("mc_truth_h0Z1_status", &mc_truth_h0Z1_status, "mc_truth_h0Z1_status/I", buffersize);
	tree->Branch("mc_truth_h0Z2_status", &mc_truth_h0Z2_status, "mc_truth_h0Z2_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl11_status", &mc_truth_h0Zl11_status, "mc_truth_h0Zl11_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl21_status", &mc_truth_h0Zl21_status, "mc_truth_h0Zl21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau11_status", &mc_truth_h0Ztau11_status, "mc_truth_h0Ztau11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau21_status", &mc_truth_h0Ztau21_status, "mc_truth_h0Ztau21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul11_status", &mc_truth_h0Ztaul11_status, "mc_truth_h0Ztaul11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul21_status", &mc_truth_h0Ztaul21_status, "mc_truth_h0Ztaul21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu11_status", &mc_truth_h0Ztaunu11_status, "mc_truth_h0Ztaunu11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu21_status", &mc_truth_h0Ztaunu21_status, "mc_truth_h0Ztaunu21_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau11_status", &mc_truth_h0Ztaunutau11_status, "mc_truth_h0Ztaunutau11_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau21_status", &mc_truth_h0Ztaunutau21_status, "mc_truth_h0Ztaunutau21_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_status", &mc_truth_h0Zq11_status, "mc_truth_h0Zq11_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_status", &mc_truth_h0Zq21_status, "mc_truth_h0Zq21_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq11_IS_status", &mc_truth_h0Zq11_IS_status, "mc_truth_h0Zq11_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq21_IS_status", &mc_truth_h0Zq21_IS_status, "mc_truth_h0Zq21_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl12_status", &mc_truth_h0Zl12_status, "mc_truth_h0Zl12_status/I", buffersize);
	tree->Branch("mc_truth_h0Zl22_status", &mc_truth_h0Zl22_status, "mc_truth_h0Zl22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau12_status", &mc_truth_h0Ztau12_status, "mc_truth_h0Ztau12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztau22_status", &mc_truth_h0Ztau22_status, "mc_truth_h0Ztau22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul12_status", &mc_truth_h0Ztaul12_status, "mc_truth_h0Ztaul12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaul22_status", &mc_truth_h0Ztaul22_status, "mc_truth_h0Ztaul22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu12_status", &mc_truth_h0Ztaunu12_status, "mc_truth_h0Ztaunu12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunu22_status", &mc_truth_h0Ztaunu22_status, "mc_truth_h0Ztaunu22_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau12_status", &mc_truth_h0Ztaunutau12_status, "mc_truth_h0Ztaunutau12_status/I", buffersize);
	tree->Branch("mc_truth_h0Ztaunutau22_status", &mc_truth_h0Ztaunutau22_status, "mc_truth_h0Ztaunutau22_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_status", &mc_truth_h0Zq12_status, "mc_truth_h0Zq12_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_status", &mc_truth_h0Zq22_status, "mc_truth_h0Zq22_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq12_IS_status", &mc_truth_h0Zq12_IS_status, "mc_truth_h0Zq12_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Zq22_IS_status", &mc_truth_h0Zq22_IS_status, "mc_truth_h0Zq22_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu11_status", &mc_truth_h0Znu11_status, "mc_truth_h0Znu11_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu21_status", &mc_truth_h0Znu21_status, "mc_truth_h0Znu21_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu12_status", &mc_truth_h0Znu12_status, "mc_truth_h0Znu12_status/I", buffersize);
	tree->Branch("mc_truth_h0Znu22_status", &mc_truth_h0Znu22_status, "mc_truth_h0Znu22_status/I", buffersize);

	tree->Branch("mc_truth_h0tau1_status", &mc_truth_h0tau1_status, "mc_truth_h0tau1_status/I", buffersize);
	tree->Branch("mc_truth_h0tau2_status", &mc_truth_h0tau2_status, "mc_truth_h0tau2_status/I", buffersize);
	tree->Branch("mc_truth_h0taul1_status", &mc_truth_h0taul1_status, "mc_truth_h0taul1_status/I", buffersize);
	tree->Branch("mc_truth_h0taunutau1_status", &mc_truth_h0taunutau1_status, "mc_truth_h0taunutau1_status/I", buffersize);
	tree->Branch("mc_truth_h0taunu1_status", &mc_truth_h0taunu1_status, "mc_truth_h0taunu1_status/I", buffersize);
	tree->Branch("mc_truth_h0taul2_status", &mc_truth_h0taul2_status, "mc_truth_h0taul2_status/I", buffersize);
	tree->Branch("mc_truth_h0taunutau2_status", &mc_truth_h0taunutau2_status, "mc_truth_h0taunutau2_status/I", buffersize);
	tree->Branch("mc_truth_h0taunu2_status", &mc_truth_h0taunu2_status, "mc_truth_h0taunu2_status/I", buffersize);

	tree->Branch("mc_truth_h0b1_status", &mc_truth_h0b1_status, "mc_truth_h0b1_status/I", buffersize);
	tree->Branch("mc_truth_h0b2_status", &mc_truth_h0b2_status, "mc_truth_h0b2_status/I", buffersize);
	tree->Branch("mc_truth_h0b1_IS_status", &mc_truth_h0b1_IS_status, "mc_truth_h0b1_IS_status/I", buffersize);
	tree->Branch("mc_truth_h0b2_IS_status", &mc_truth_h0b2_IS_status, "mc_truth_h0b2_IS_status/I", buffersize);
	
	tree->Branch("mc_truth_t_status", &mc_truth_t_status, "mc_truth_t_status/I", buffersize);
	tree->Branch("mc_truth_tb_status", &mc_truth_tb_status, "mc_truth_tb_status/I", buffersize);
	tree->Branch("mc_truth_tb_IS_status", &mc_truth_tb_IS_status, "mc_truth_tb_IS_status/I", buffersize);

	tree->Branch("mc_truth_tW_status", &mc_truth_tW_status, "mc_truth_tW_status/I", buffersize);
	tree->Branch("mc_truth_tWnu_status", &mc_truth_tWnu_status, "mc_truth_tWnu_status/I", buffersize);
	tree->Branch("mc_truth_tWnutau_status", &mc_truth_tWnutau_status, "mc_truth_tWnutau_status/I", buffersize);
	tree->Branch("mc_truth_tWl_status", &mc_truth_tWl_status, "mc_truth_tWl_status/I", buffersize);
	tree->Branch("mc_truth_tWtau_status", &mc_truth_tWtau_status, "mc_truth_tWtau_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunu_status", &mc_truth_tWtaunu_status, "mc_truth_tWtaunu_status/I", buffersize);
	tree->Branch("mc_truth_tWtaunutau_status", &mc_truth_tWtaunutau_status, "mc_truth_tWtaunutau_status/I", buffersize);
	tree->Branch("mc_truth_tWtaul_status", &mc_truth_tWtaul_status, "mc_truth_tWtaul_status/I", buffersize);
	tree->Branch("mc_truth_tWq1_status", &mc_truth_tWq1_status, "mc_truth_tWq1_status/I", buffersize);
	tree->Branch("mc_truth_tWq2_status", &mc_truth_tWq2_status, "mc_truth_tWq2_status/I", buffersize);
	tree->Branch("mc_truth_tWq1_IS_status", &mc_truth_tWq1_IS_status, "mc_truth_tWq1_IS_status/I", buffersize);
	tree->Branch("mc_truth_tWq2_IS_status", &mc_truth_tWq2_IS_status, "mc_truth_tWq2_IS_status/I", buffersize);

	tree->Branch("mc_truth_j1_status", &mc_truth_j1_status, "mc_truth_j1_status/I", buffersize);
	tree->Branch("mc_truth_j2_status", &mc_truth_j2_status, "mc_truth_j2_status/I", buffersize);
	tree->Branch("mc_truth_j3_status", &mc_truth_j3_status, "mc_truth_j3_status/I", buffersize);
     }
   
   if( doWrite("gen_PVz") ) tree->Branch("gen_PVz", &gen_PVz, "gen_PVz/F", buffersize);
   
   if( doWrite("gen_all") || doWrite("gen_stop"))
     {
	tree->Branch("gen_n", &gen_n, "gen_n/I", buffersize);
	tree->Branch("gen_pt", "std::vector<float>", &gen_pt, buffersize);
	tree->Branch("gen_eta", "std::vector<float>", &gen_eta, buffersize);
	tree->Branch("gen_phi", "std::vector<float>", &gen_phi, buffersize);
	tree->Branch("gen_m", "std::vector<float>", &gen_m, buffersize);
	tree->Branch("gen_E", "std::vector<float>", &gen_E, buffersize);
	tree->Branch("gen_status", "std::vector<int>", &gen_status, buffersize);
	tree->Branch("gen_id", "std::vector<int>", &gen_id, buffersize);
	tree->Branch("gen_charge", "std::vector<int>", &gen_charge, buffersize);
	tree->Branch("gen_index", "std::vector<int>", &gen_index, buffersize);
	tree->Branch("gen_mother_index", "std::vector<int>", &gen_mother_index, buffersize);
	tree->Branch("gen_daughter_n", "std::vector<int>", &gen_daughter_n, buffersize);
	tree->Branch("gen_daughter_index", "std::vector<std::vector<int> >", &gen_daughter_index, buffersize);
     }
   if( doWrite("gen_stop") || doWrite("gen_stop_mass"))
     {
	tree->Branch("gen_stop_m", "std::vector<float>", &gen_stop_m, buffersize);
	tree->Branch("gen_neutralino_m", "std::vector<float>", &gen_neutralino_m, buffersize);
     }
   
   if( doWrite("genTTX_id") ) tree->Branch("genTTX_id", &genTTX_id, "genTTX_id/I", buffersize);
   
   tree->Branch("n_presel_jets",     &n_presel_jets,     "n_presel_jets/I",      buffersize);
   tree->Branch("n_presel_btag",     &n_presel_btag,     "n_presel_btag/I",      buffersize);
   tree->Branch("n_presel_electrons", &n_presel_electrons, "n_presel_electrons/I",  buffersize);
   tree->Branch("n_presel_muons",     &n_presel_muons,     "n_presel_muons/I",      buffersize);
   tree->Branch("n_presel_tau",      &n_presel_tau,      "n_presel_tau/I",       buffersize);
   
}

bool FlatTree::doWrite(const std::string& name)
{
   std::map<std::string,bool>::iterator it = conf.find(name);
   if( it != conf.end() )
     return it->second;
   return 0;
}
