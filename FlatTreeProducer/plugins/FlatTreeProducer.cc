#include <memory>
#include <iostream>
#include <sstream>

#include "TRegexp.h"
#include "TString.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "VUBFlatTree/FlatTreeProducer/interface/tinyxml2.h"

#include "VUBFlatTree/FlatTreeProducer/interface/FlatTree.hh"
#include "VUBFlatTree/FlatTreeProducer/interface/MCTruth.hh"
#include "VUBFlatTree/FlatTreeProducer/interface/GenTTXCategorizer.hh"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TFile.h"
#include "TThreadSlots.h"
#include "TROOT.h"
#include "Compression.h"

using namespace tinyxml2;

class FlatTreeProducer : public edm::EDAnalyzer
{
 public:
   
   explicit FlatTreeProducer(const edm::ParameterSet&);
   ~FlatTreeProducer();
   
   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   bool foundTrigger(const std::string& name) const;
   
 private:
   
   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
   
   virtual void beginRun(const edm::Run&, const edm::EventSetup&);
   virtual void endRun(const edm::Run&, const edm::EventSetup&);
   
   TMVA::Reader* BookLeptonMVAReaderMoriond18(std::string basePath, std::string weightFileName, std::string type);
   
   void KeepEvent();
   bool isFloat(const std::string& s);
   bool isFloat(const boost::any & operand);
   bool isBool(const boost::any & operand);
   bool isInt(const boost::any & operand);
   void fillCutVector(const char* cut_type, std::string& cut_value, std::map<std::string, boost::any>& vec);
   void AddValue(const std::string& name);
   void ReadConfFile(const std::string& confFile);
   int CheckAlgo(const std::map<std::string, boost::any>& jet_algo, const char* name, std::string& algo);
   std::string CheckAlgos();
   
   template <typename T>
     boost::any comp(const std::string& op, T v1, T v2);
   template <typename T>
     void CheckVectorCut(const std::vector<T>& v, const std::string& name);
   
   template <typename T>
     void CompareAlgo(const std::string& algo, T conf_algo_value);
   template <typename T>
     void CheckJet(const std::vector<T>& vJetPt, const std::string& jetpt, const std::vector<T>& vJetEta, const std::string& jeteta, const std::string& algo);
   template <typename T>
     void CheckElectron(const std::vector<T>& vElPt, const std::string& elpt, const std::vector<T>& vElEta, const std::string& eleta);
   template <typename T>
     void CheckMuon(const std::vector<T>& vMuonPt, const std::string& muonpt, const std::vector<T>& vMuonEta, const std::string& muoneta);
   
   FlatTree* ftree;
   const edm::Service<TFileService> fs;
   
   TH1D* hcount;
   TH1D* hweight;
   
   TMVA::Reader* ele_reader;
   TMVA::Reader* mu_reader;
   
   float lepMVA_pt;
   float lepMVA_eta;
   float lepMVA_miniRelIsoCharged;
   float lepMVA_miniRelIsoNeutral;
   float lepMVA_jetPtRatio;
   float lepMVA_jetPtRelv2;
   float lepMVA_jetBTagCSV;
   float lepMVA_sip3d;
   float lepMVA_dxy;
   float lepMVA_dz;
   float lepMVA_mvaId;
   float lepMVA_jetNDauChargedMVASel;
   
   XMLDocument xmlconf;
   
   std::string dataFormat_;
   bool isData_;
   bool applyMETFilters_;
   bool fillMCScaleWeight_;
   bool fillPUInfo_;
   int nPdf_;
   
   HLTConfigProvider hltConfig_;
   HLTPrescaleProvider hltPrescale_;
   
   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
   edm::EDGetTokenT<edm::TriggerResults> triggerBitsPAT_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
   edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
   
   edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
   edm::EDGetTokenT<pat::ElectronCollection> electronPATToken_;
   edm::EDGetTokenT<edm::View<reco::GsfElectron> > electronToken_;
   edm::EDGetTokenT<pat::MuonCollection> muonToken_;
   edm::EDGetTokenT<pat::TauCollection> tauToken_;
   edm::EDGetTokenT<pat::JetCollection> jetToken_;
   edm::EDGetTokenT<edm::View<pat::Jet> > viewJetToken_;
   edm::EDGetTokenT<pat::JetCollection> jetPuppiToken_;
   edm::EDGetTokenT<pat::JetCollection> ak8jetToken_;
   edm::EDGetTokenT<pat::JetCollection> ak10jetToken_;
   edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
   edm::EDGetTokenT<std::vector<pat::MET> > metTokenAOD_;
   edm::EDGetTokenT<pat::METCollection> metTokenPuppi_;
   edm::EDGetTokenT<pat::METCollection> metTokenNoHF_;
   edm::EDGetTokenT<pat::METCollection> metTokenMINIAOD_;
   edm::EDGetTokenT<double> rhoToken_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
   edm::EDGetTokenT<LHEEventProduct> LHEEventProductToken_;
   edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;
   edm::EDGetTokenT<reco::BeamSpot> bsToken_;
   edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;
   edm::EDGetTokenT<reco::ConversionCollection> hConversionsToken_;

   //        edm::EDGetTokenT<bool> badMuonFilterToken_;
   //        edm::EDGetTokenT<bool> badChargedCandidateFilterToken_;

   edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoCBIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseCBIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumCBIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleTightCBIdMapToken_;

   edm::EDGetTokenT<edm::ValueMap<bool> > ele90NoIsoMVAIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > ele80NoIsoMVAIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseNoIsoMVAIdMapToken_;
   
   edm::EDGetTokenT<edm::ValueMap<bool> > ele90IsoMVAIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > ele80IsoMVAIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIsoMVAIdMapToken_;
   
   edm::EDGetTokenT<edm::ValueMap<float> > mvaIsoValuesMapToken_;
   edm::EDGetTokenT<edm::ValueMap<float> > mvaNoIsoValuesMapToken_;
   
   //        edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;
   // 
   //        edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > vetoIdFullInfoMapToken_;
   //        edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > mediumIdFullInfoMapToken_;
   // 
   edm::EDGetTokenT<double> metSigToken_;
   edm::EDGetTokenT<math::Error<2>::type> metCovToken_;
   edm::EDGetTokenT<edm::ValueMap<float> > qgToken_;
   
   std::vector<std::string> filterTriggerNames_;
   
   JetCorrectionUncertainty *jecUnc;
   
   edm::EDGetTokenT<reco::GenJetCollection> genTTXJetsToken_;
   
   edm::EDGetTokenT<std::vector<int> > genTTXBHadJetIndexToken_;
   edm::EDGetTokenT<std::vector<int> > genTTXBHadFlavourToken_;
   edm::EDGetTokenT<std::vector<int> > genTTXBHadFromTopWeakDecayToken_;
   edm::EDGetTokenT<std::vector<reco::GenParticle> > genTTXBHadPlusMothersToken_;
   edm::EDGetTokenT<std::vector<std::vector<int> > > genTTXBHadPlusMothersIndicesToken_;
   edm::EDGetTokenT<std::vector<int> > genTTXBHadIndexToken_;
   edm::EDGetTokenT<std::vector<int> > genTTXBHadLeptonHadronIndexToken_;
   edm::EDGetTokenT<std::vector<int> > genTTXBHadLeptonViaTauToken_;
   
   edm::EDGetTokenT<std::vector<int> > genTTXCHadJetIndexToken_;
   edm::EDGetTokenT<std::vector<int> > genTTXCHadFlavourToken_;
   edm::EDGetTokenT<std::vector<int> > genTTXCHadFromTopWeakDecayToken_;
   edm::EDGetTokenT<std::vector<int> > genTTXCHadBHadronIdToken_;
};

bool FlatTreeProducer::isInt(const boost::any & operand)
{
   return operand.type() == typeid(int);
}
bool FlatTreeProducer::isFloat(const boost::any & operand)
{
   return operand.type() == typeid(float);
}
bool FlatTreeProducer::isBool(const boost::any& operand)
{
   return operand.type() == typeid(bool);
}

template<typename T>
  boost::any FlatTreeProducer::comp(const std::string& op, T v1, T v2)
{
   if( !op.compare("||") )
     {
        return v1 || v2;
     }
   else if( !op.compare("&&") )
     {
        return v1 && v2;
     }
   
   if( v1 == typeid(int) && v2 == typeid(int) )
     {
        try
	  {
	     if( !op.compare("|") )
	       {
		  return v1 | v2;
	       }
	     else if( !op.compare("&") )
	       {
		  return v1 & v2;
	       }
	  }
        catch (...)
	  {
	  }
     }
   return "";
}

void FlatTreeProducer::AddValue(const std::string& name)
{
   if( !name.compare("n_presel_jets") )
     {
        ftree->n_presel_jets += 1;
     }
   else if( !name.compare("n_presel_electrons") )
     {
        ftree->n_presel_electrons += 1;
     }
   else if( !name.compare("n_presel_muons") )
     {
        ftree->n_presel_muons += 1;
     }
}

template <typename T>
  void FlatTreeProducer::CheckVectorCut(const std::vector<T>& v, const std::string& name)
{
   if( ftree->keep_conf.find(name) != ftree->keep_conf.end() )
     {
        std::map<std::string, std::map<std::string, boost::any> > keep_conf = ftree->keep_conf;
        std::map<std::string, boost::any> map_conf = keep_conf[name];
        for( unsigned int i=0;i<v.size();++i )
	  {
	     if( isInt(map_conf["cut_min"]) && isInt(map_conf["cut_max"]) )
            {
	       if( v[i] < boost::any_cast<int>(map_conf["cut_min"]) || v[i] > boost::any_cast<int>(map_conf["cut_max"]) )
		 {
                    AddValue(name);
		 }
            }
	     else if( isFloat(map_conf["cut_min"]) && isFloat(map_conf["cut_max"]) )
	       {
		  if( v[i] < boost::any_cast<float>(map_conf["cut_min"]) || v[i] > boost::any_cast<float>(map_conf["cut_max"]) )
		    {
		       AddValue(name);
		    }
	       }
	     else{std::cout << "Wrong types" << std::endl;
	     }
	  }
     }
}

template <typename T>
  void FlatTreeProducer::CompareAlgo(const std::string& algo, T conf_algo_value)
{
   if( !algo.compare("jet_JBP") )
     {
        for( unsigned int i=0;i<ftree->jet_JBP.size();++i )
	  if( ftree->jet_JBP[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_JP") )
     {
        for( unsigned int i=0;i<ftree->jet_JP.size();++i )
	  if( ftree->jet_JP[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_TCHP") )
     {
        for( unsigned int i=0;i<ftree->jet_TCHP.size();++i )
	  if( ftree->jet_TCHP[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_TCHE") )
     {
        for( unsigned int i=0;i<ftree->jet_TCHE.size();++i )
	  if( ftree->jet_TCHE[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_SSVHP") )
     {
        for( unsigned int i=0;i<ftree->jet_SSVHP.size();++i )
	  if( ftree->jet_SSVHP[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_SSVHE") )
     {
        for( unsigned int i=0;i<ftree->jet_SSVHE.size();++i )
	  if( ftree->jet_SSVHE[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_CMVA") )
     {
        for( unsigned int i=0;i<ftree->jet_CMVA.size();++i )
	  if( ftree->jet_CMVA[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_CSVv2") )
     {
        for( unsigned int i=0;i<ftree->jet_CSVv2.size();++i )
	  if( ftree->jet_CSVv2[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
   else if( !algo.compare("jet_partonFlavour") )
     {
        for( unsigned int i=0;i<ftree->jet_partonFlavour.size();++i )
	  if( ftree->jet_partonFlavour[i] > conf_algo_value )
	    AddValue("n_presel_jets");
     }
}

// FIXME : refactorize CheckJet, CheckElectron and CheckMuon
template <typename T>
  void FlatTreeProducer::CheckJet(const std::vector<T>& vJetPt, const std::string& jetpt, const std::vector<T>& vJetEta, const std::string& jeteta, const std::string& algo)
{
   std::map<std::string, std::map<std::string, boost::any> > keep_conf = ftree->keep_conf;
   
   /*for(auto it = keep_conf.cbegin(); it != keep_conf.cend(); ++it)
    std::cout << it->first << " " << std::endl;
    */
   
   //-- check that vJetPt and vJetEta have the same size
   if(vJetPt.size()!=vJetEta.size()){
      std::cerr<<" CheckJet: pt and eta vectors have different sizes !  - Function not applied" <<std::endl;
      return ;
   }
   
   if( keep_conf.find(jetpt) != keep_conf.end() )
     {
        std::map<std::string, boost::any> conf_pt = keep_conf[jetpt];
        for( unsigned int i=0;i<vJetPt.size();++i )
	  {
	     //std::cout<<"pt = "<<vJetPt[i]<<std::endl;
	     if( isInt(conf_pt["cut_min"]) )
	       {
		  //std::cout<<"pass"<<std::endl;
		  if( vJetPt[i] > boost::any_cast<int>(conf_pt["cut_min"]) )
		    {
		       //std::cout<<"here too .."<<std::endl;
		       if( keep_conf.find(jeteta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[jeteta];
			    //for( unsigned int i=0;i<vJetEta.size();++i )
			    // {
			    if( isInt(conf_eta["cut_max"]) )
			      //@EC@
			      //if( isFloat(conf_eta["cut_max"]) )
			      {
				 if( fabs(vJetEta[i]) < boost::any_cast<int>(conf_eta["cut_max"]) )
				   {
				      if( !algo.empty() )
					{
					   std::map<std::string, boost::any> conf_algo = keep_conf[algo];
					   if( isInt(conf_algo["cut_algo"]) )
					     {
						CompareAlgo(algo, boost::any_cast<int>(conf_algo["cut_algo"]));
					     }
					   else
					     {
						std::cout << "'jet_pt' : 'cut_min' and 'cut_max' are type 'int', 'cut_algo' cannot be a 'float'." << std::endl;
					     }
					}
				      else
					{
					   AddValue("n_presel_jets");
					}
				   }
				 //}
				 else
				   {
				      std::cout << "'jet_pt' : 'cut_min' is a type 'int', 'jet_eta' : 'cut_max' cannot be a 'float'." << std::endl;
				      break;
				   }
			      }
			 }
		       else
			 {
			    std::cout << "'jet_pt' : 'cut_min' set, but not 'jet_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else if( isFloat(conf_pt["cut_min"]) )
	       {
		  if( vJetPt[i] > boost::any_cast<float>(conf_pt["cut_min"]) )
		    {
		       if( keep_conf.find(jeteta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[jeteta];
			    if( isFloat(conf_eta["cut_max"]) )
			      {
				 if( fabs(vJetEta[i]) < boost::any_cast<float>(conf_eta["cut_max"]) )
				   {
				      if( !algo.empty() )
					{
					   std::map<std::string, boost::any> conf_algo = keep_conf[algo];
					   if( isFloat(conf_algo["cut_algo"]) )
					     {
						CompareAlgo(algo, boost::any_cast<float>(conf_algo["cut_algo"]));
					     }
					   else
					     {
						std::cout << "'jet_pt' : 'cut_min' and 'cut_max' are type 'float', 'cut_algo' cannot be an 'int'." << std::endl;
					     }
					}
				      else
					{
					   AddValue("n_presel_jets");
					}
				   }
			      }
			    else
			      {
				 std::cout << "'jet_pt' : 'cut_min' is a type 'float', 'jet_eta' : 'cut_max' cannot be an 'int'." << std::endl;
				 break;
			      }
			 }
		       else
			 {
			    std::cout << "'jet_pt' : 'cut_min' set, but not 'jet_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else{std::cout << "'jet_pt' : wrong types." << std::endl;
	     }
	  }
     }
}

template <typename T>
  void FlatTreeProducer::CheckElectron(const std::vector<T>& vElPt, const std::string& elpt, const std::vector<T>& vElEta, const std::string& eleta)
{
   //-- check that vElPt and vElEta have the same size
   if(vElPt.size()!=vElEta.size()){
      std::cerr<<" CheckElectron: pt and eta vectors have different sizes !  - Function not applied" <<std::endl;
      return ;
   }
   std::map<std::string, std::map<std::string, boost::any> > keep_conf = ftree->keep_conf;
   if( keep_conf.find(elpt) != keep_conf.end() )
     {
        std::map<std::string, boost::any> conf_pt = keep_conf[elpt];
        for( unsigned int i=0;i<vElPt.size();++i )
	  {
	     if( isInt(conf_pt["cut_min"]) )
	       {
		  if( vElPt[i] > boost::any_cast<int>(conf_pt["cut_min"]) )
		    {
		       if( keep_conf.find(eleta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[eleta];
			    //for( unsigned int i=0;i<vElEta.size();++i )
			    // {
			    if( isInt(conf_eta["cut_max"]) )
			      {
				 if( fabs(vElEta[i]) < boost::any_cast<int>(conf_eta["cut_max"]) )
				   {
				      AddValue("n_presel_electrons");
				   }
				 //}
				 else
				   {
				      std::cout << "'el_pt' : 'cut_min' is a type 'int', 'el_eta' : 'cut_max' cannot be a 'float'." << std::endl;
				      break;
				   }
			      }
			 }
		       else
			 {
			    std::cout << "'el_pt' : 'cut_min' set, but not 'el_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else if( isFloat(conf_pt["cut_min"]) )
	       {
		  if( vElPt[i] > boost::any_cast<float>(conf_pt["cut_min"]) )
		    {
		       if( keep_conf.find(eleta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[eleta];
			    //for( unsigned int i=0;i<vElEta.size();++i )
			    //{
			    if( isFloat(conf_eta["cut_max"]) )
			      {
				 if( fabs(vElEta[i]) < boost::any_cast<float>(conf_eta["cut_max"]) )
				   {
				      AddValue("n_presel_electrons");
				   }
			      }
			    else
			      {
				 std::cout << "'el_pt' : 'cut_min' is a type 'float', 'el_eta' : 'cut_max' cannot be an 'int'." << std::endl;
				 break;
			      }
			    //}
			 }
		       else
			 {
			    std::cout << "'el_pt' : 'cut_min' set, but not 'el_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else{std::cout << "'el_pt' : wrong types." << std::endl;
	     }
	  }
     }
}

template <typename T>
  void FlatTreeProducer::CheckMuon(const std::vector<T>& vMuonPt, const std::string& muonpt, const std::vector<T>& vMuonEta, const std::string& muoneta)
{
   //-- check that vMuonPt and vMuonEta have the same size
   if(vMuonPt.size()!=vMuonEta.size()){
      std::cerr<<" CheckMuon: pt and eta vectors have different sizes !  - Function not applied" <<std::endl;
      return ;
   }
   std::map<std::string, std::map<std::string, boost::any> > keep_conf = ftree->keep_conf;
   if( keep_conf.find(muonpt) != keep_conf.end() )
     {
        std::map<std::string, boost::any> conf_pt = keep_conf[muonpt];
        for( unsigned int i=0;i<vMuonPt.size();++i )
	  {
	     if( isInt(conf_pt["cut_min"]) )
	       {
		  if( vMuonPt[i] > boost::any_cast<int>(conf_pt["cut_min"]) )
		    {
		       if( keep_conf.find(muoneta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[muoneta];
			    //for( unsigned int i=0;i<vMuonEta.size();++i )
			    //{
			    if( isInt(conf_eta["cut_max"]) )
			      {
				 if( fabs(vMuonEta[i]) < boost::any_cast<int>(conf_eta["cut_max"]) )
				   {
				      AddValue("n_presel_muons");
				   }
			      }
			    else
			      {
				 std::cout << "'mu_pt' : 'cut_min' is a type 'int', 'mu_eta' : 'cut_max' cannot be a 'float'." << std::endl;
				 break;
			      }
			    //}
			 }
		       else
			 {
			    std::cout << "'mu_pt' : 'cut_min' set, but not 'mu_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else if( isFloat(conf_pt["cut_min"]) )
	       {
		  if( vMuonPt[i] > boost::any_cast<float>(conf_pt["cut_min"]) )
		    {
		       if( keep_conf.find(muoneta) != keep_conf.end() )
			 {
			    std::map<std::string, boost::any> conf_eta = keep_conf[muoneta];
			    //for( unsigned int i=0;i<vMuonEta.size();++i )
			    // {
			    if( isFloat(conf_eta["cut_max"]) )
			      {
				 if( fabs(vMuonEta[i]) < boost::any_cast<float>(conf_eta["cut_max"]) )
				   {
				      AddValue("n_presel_muons");
				   }
			      }
			    else
			      {
				 std::cout << "'mu_pt' : 'cut_min' is a type 'float', 'mu_eta' : 'cut_max' cannot be an 'int'." << std::endl;
				 break;
			      }
			    //}
			 }
		       else
			 {
			    std::cout << "'mu_pt' : 'cut_min' set, but not 'mu_eta' : 'cut_max'" << std::endl;
			    break;
			 }
		    }
	       }
	     else{std::cout << "'mu_pt' : wrong types." << std::endl;
	     }
	  }
     }
}

int FlatTreeProducer::CheckAlgo(const std::map<std::string, boost::any>& jet_algo, const char* name, std::string& algo)
{
   if( jet_algo.count("cut_algo") == 1 )
     {
        algo = name;
        return 1;
     }
   return 0;
}

std::string FlatTreeProducer::CheckAlgos()
{
   int nbAlgo = 0;
   std::string algo("");
   std::map<std::string, std::map<std::string, boost::any> > keep_conf = ftree->keep_conf;
   std::map<std::string, boost::any> jet_JBP = keep_conf["jet_JBP"];
   std::map<std::string, boost::any> jet_JP = keep_conf["jet_JP"];
   std::map<std::string, boost::any> jet_TCHP = keep_conf["jet_TCHP"];
   std::map<std::string, boost::any> jet_TCHE = keep_conf["jet_TCHE"];
   std::map<std::string, boost::any> jet_SSVHP = keep_conf["jet_SSVHP"];
   std::map<std::string, boost::any> jet_SSVHE = keep_conf["jet_SSVHE"];
   std::map<std::string, boost::any> jet_CMVA = keep_conf["jet_CMVA"];
   std::map<std::string, boost::any> jet_CSVv2 = keep_conf["jet_CSVv2"];
   std::map<std::string, boost::any> jet_partonFlavour = keep_conf["jet_partonFlavour"];
   
   nbAlgo += CheckAlgo(jet_JBP, "jet_JBP", algo);
   nbAlgo += CheckAlgo(jet_JP, "jet_JP", algo);
   nbAlgo += CheckAlgo(jet_TCHP, "jet_TCHP", algo);
   nbAlgo += CheckAlgo(jet_TCHE, "jet_TCHE", algo);
   nbAlgo += CheckAlgo(jet_SSVHP, "jet_SSVHP", algo);
   nbAlgo += CheckAlgo(jet_SSVHE, "jet_SSHVE", algo);
   nbAlgo += CheckAlgo(jet_CMVA, "jet_CMVA", algo);
   nbAlgo += CheckAlgo(jet_CSVv2, "jet_CSVv2", algo);
   nbAlgo += CheckAlgo(jet_partonFlavour, "jet_partonFlavour", algo);
   
   if( nbAlgo > 1 )
     {
        std::cout << "Different algorithms are set, please choose only one. Algorithm not considered." << std::endl;
        algo.clear();
        return algo;
     }
   return algo;
}

void FlatTreeProducer::KeepEvent()
{
   // Jets
   std::string algo = this->CheckAlgos();
   // jet_pt > cut_min && jet_eta >
   this->CheckJet(ftree->jet_pt, "jet_pt", ftree->jet_eta, "jet_eta", algo);
   
   // Electron
   // el_pt> XX && fabs(el_eta)> YY && iso < ZZ
   // TODO ISO
   this->CheckElectron(ftree->el_pt, "el_pt", ftree->el_eta, "el_eta");
   
   // Muon idem
   this->CheckMuon(ftree->mu_pt, "mu_pt", ftree->mu_eta, "mu_eta");
}

bool FlatTreeProducer::isFloat(const std::string& s)
{
   return !std::all_of(s.begin(), s.end(), ::isdigit);
}

void FlatTreeProducer::fillCutVector(const char* cut_type, std::string& cut_value, std::map<std::string, boost::any> & vmap)
{
   if( !cut_value.empty() )
     {
        bool t_cut = isFloat(cut_value);
        if( t_cut )
	  {
	     float fcut = atof(cut_value.c_str());
	     vmap[cut_type] = fcut;
	  }
        else if( !t_cut )
	  {
	     int fcut = atoi(cut_value.c_str());
	     vmap[cut_type] = fcut;
	  }
        else
	  {
	     std::cout << "Warning ! Different types of cut (int or float). Cut Skipped." << std::endl;
	  }
     }
   cut_value = "";
}

void FlatTreeProducer::ReadConfFile(const std::string& confFile)
{
   xmlconf.LoadFile(confFile.c_str());
   
   //-----------------------------
   //-- Read preselection section
   //-----------------------------
   XMLElement* tElement = xmlconf.FirstChildElement("preselection")->FirstChildElement("obj");
   
   for( XMLElement* child=tElement;child!=0;child=child->NextSiblingElement() )
     {
        const std::string& obj = child->ToElement()->Attribute("name");
        const std::string& min = child->ToElement()->Attribute("min");
        int imin = atoi(min.c_str());
        if(!obj.compare("jet")) 	ftree->n_presel_jets_min=imin; 
        if(!obj.compare("electron")) 	ftree->n_presel_electrons_min=imin; 
        if(!obj.compare("lepton")) 	ftree->n_presel_leptons_min=imin; 
        if(!obj.compare("muon")) 	ftree->n_presel_muons_min=imin; 
        if(!obj.compare("MET")) 	ftree->presel_MET_min=imin; 	
     }
   tElement = xmlconf.FirstChildElement("preselection")->FirstChildElement("info");
   const std::string& activ = tElement->ToElement()->Attribute("activate");
   ftree->apply_presel = (bool) (atoi(activ.c_str()));
   //std::cout<<activ<<" "<<ftree->apply_presel<<std::endl;
   //if(ftree->apply_presel) std::cout<<"TEST"<<std::endl;
   
   
   //-----------------------------
   //-- Read variables section
   //-----------------------------
   tElement = xmlconf.FirstChildElement("variables")->FirstChildElement("var");
   
   std::string vcutmin("");
   std::string vcutmax("");
   std::string vcutiso("");
   std::string vcutalgo("");
   std::string vcutexpr("");
   for( XMLElement* child=tElement;child!=0;child=child->NextSiblingElement() )
     {
        const std::string& vname = child->ToElement()->Attribute("name");
        const std::string& vsave = child->ToElement()->Attribute("save");
        if (child->ToElement()->Attribute("cut_min"))
	  {
	     vcutmin = child->ToElement()->Attribute("cut_min");
	  }
        if (child->ToElement()->Attribute("cut_max"))
	  {
	     vcutmax = child->ToElement()->Attribute("cut_max");
	  }
        if (child->ToElement()->Attribute("cut_iso"))
	  {
	     vcutiso = child->ToElement()->Attribute("cut_iso");
	  }
        if (child->ToElement()->Attribute("cut_expr"))
	  {
	     ;
	  }
        if (child->ToElement()->Attribute("cut_algo"))
	  {
	     vcutalgo = child->ToElement()->Attribute("cut_algo");
	  }
	
        std::map<std::string, boost::any > vmap;
        bool bsave = atoi(vsave.c_str());
        fillCutVector("cut_min", vcutmin, vmap);
        fillCutVector("cut_max", vcutmax, vmap);
        fillCutVector("cut_iso", vcutiso, vmap);
        fillCutVector("cut_algo", vcutalgo, vmap);
        if (vmap.size() > 0)
	  ftree->keep_conf.insert(std::make_pair(vname, vmap));
        ftree->conf.insert(std::make_pair(vname,bsave));
     }
}

TMVA::Reader* FlatTreeProducer::BookLeptonMVAReaderMoriond18(std::string basePath, std::string weightFileName, std::string type)
{
   TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");
   
   reader->AddVariable("LepGood_pt",                                     &lepMVA_pt);
   reader->AddVariable("LepGood_eta",                                    &lepMVA_eta); 
   reader->AddVariable("LepGood_jetNDauChargedMVASel",                   &lepMVA_jetNDauChargedMVASel);
   reader->AddVariable("LepGood_miniRelIsoCharged",                      &lepMVA_miniRelIsoCharged);
   reader->AddVariable("LepGood_miniRelIsoNeutral",                      &lepMVA_miniRelIsoNeutral);
   reader->AddVariable("LepGood_jetPtRelv2",                             &lepMVA_jetPtRelv2);
   reader->AddVariable("max(LepGood_jetBTagCSV,0)",                      &lepMVA_jetBTagCSV);
   reader->AddVariable("(LepGood_jetBTagCSV>-5)*min(LepGood_jetPtRatiov2,1.5)+(LepGood_jetBTagCSV<-5)/(1+LepGood_relIso04)", &lepMVA_jetPtRatio);
   reader->AddVariable("LepGood_sip3d",                                  &lepMVA_sip3d);
   reader->AddVariable("log(abs(LepGood_dxy))",                          &lepMVA_dxy);
   reader->AddVariable("log(abs(LepGood_dz))",                           &lepMVA_dz);
   if( type == "ele" ) reader->AddVariable("LepGood_mvaIdFall17noIso",   &lepMVA_mvaId);
   else reader->AddVariable("LepGood_segmentCompatibility",              &lepMVA_mvaId);
   
   reader->BookMVA("BDTG method", basePath+"/"+weightFileName);
   
   return reader;
}

FlatTreeProducer::FlatTreeProducer(const edm::ParameterSet& iConfig):
hltPrescale_(iConfig,consumesCollector(),*this)
{
   // ###
   // Temporarily redirecting stdout to avoid huge TMVA loading dump
   // ###
   //    std::cout << "Temporarily redirecting stdout to avoid huge TMVA dump when loading MVA readers..." << std::endl;
   //    std::stringstream tmpBuffer;
   //    std::streambuf* oldStdout = std::cout.rdbuf(tmpBuffer.rdbuf());
   // 
   // ###############
   // #  Load MVAs  #
   // ###############
   //
   const char* cmssw_base = std::getenv("CMSSW_BASE");
   std::string FlatTreeProducerLepMVAPath = std::string(cmssw_base)+"/src/VUBFlatTree/FlatTreeProducer/data/lepMVA";
   mu_reader        = BookLeptonMVAReaderMoriond18(FlatTreeProducerLepMVAPath, "mu_BDTG.weights.xml", "mu");
   ele_reader       = BookLeptonMVAReaderMoriond18(FlatTreeProducerLepMVAPath, "el_BDTG.weights.xml", "ele");
   
   // ###
   // Restore stdout
   // ###
   //    std::cout.rdbuf(oldStdout);
   //    std::cout << "Stdout now restored." << std::endl;
   // 
   // ########################
   // #  Create output tree  #
   // ########################
   //
   TFile& f = fs->file();
   f.SetCompressionAlgorithm(ROOT::kZLIB);
   f.SetCompressionLevel(9);
   ftree = new FlatTree(fs->make<TTree>("tree","tree"));
   
   // #############################################################
   // #  Read parameters from python file and get consume tokens  #
   // #############################################################
   //
   dataFormat_           = iConfig.getParameter<std::string>("dataFormat");
   nPdf_                 = iConfig.getParameter<int>("nPDF");
   fillMCScaleWeight_    = iConfig.getParameter<bool>("fillMCScaleWeight");
   fillPUInfo_           = iConfig.getParameter<bool>("fillPUInfo");
   isData_               = iConfig.getParameter<bool>("isData");
   applyMETFilters_      = iConfig.getParameter<bool>("applyMETFilters");
   triggerBits_          = consumes<edm::TriggerResults>(edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("HLT")));
   triggerBitsPAT_       = consumes<edm::TriggerResults>(edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("PAT")));
   triggerObjects_       = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"));
   triggerPrescales_     = consumes<pat::PackedTriggerPrescales>(edm::InputTag(std::string("patTrigger")));
   vertexToken_          = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexInput"));
   electronPATToken_     = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronPATInput"));
   electronToken_        = consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electronInput"));
   muonToken_            = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonInput"));
   tauToken_             = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tauInput"));
   jetToken_             = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetInput"));
   viewJetToken_         = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetInput"));
   ak8jetToken_          = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8jetInput"));
   ak10jetToken_         = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak10jetInput"));
   jetPuppiToken_        = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetPuppiInput"));
   genJetToken_          = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetInput"));
   metTokenAOD_          = consumes<std::vector<pat::MET> >(iConfig.getParameter<edm::InputTag>("metInput"));
   metTokenMINIAOD_      = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metInput"));
   metTokenNoHF_         = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metNoHFInput"));
   metTokenPuppi_        = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metPuppiInput"));
   rhoToken_             = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInput"));
   genParticlesToken_    = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticlesInput"));
   genEventInfoToken_    = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoInput"));
   LHEEventProductToken_ = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventProductInput"));
   puInfoToken_          = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfoInput"));
   bsToken_              = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bsInput"));
   pfcandsToken_         = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfcandsInput"));
   hConversionsToken_    = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("hConversionsInput"));
   
   //    badMuonFilterToken_   = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadMuonFilter"));
   //    badChargedCandidateFilterToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"));
   // 
   eleVetoCBIdMapToken_    = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoCBIdMap"));
   eleLooseCBIdMapToken_   = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseCBIdMap"));
   eleMediumCBIdMapToken_  = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumCBIdMap"));
   eleTightCBIdMapToken_   = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightCBIdMap"));
   
   ele90NoIsoMVAIdMapToken_ = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ele90NoIsoMVAIdMap"));
   ele80NoIsoMVAIdMapToken_  = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ele80NoIsoMVAIdMap"));
   eleLooseNoIsoMVAIdMapToken_  = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseNoIsoMVAIdMap"));
   
   ele90IsoMVAIdMapToken_ = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ele90IsoMVAIdMap"));
   ele80IsoMVAIdMapToken_  = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ele80IsoMVAIdMap"));
   eleLooseIsoMVAIdMapToken_  = consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIsoMVAIdMap"));
   
   mvaIsoValuesMapToken_      = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaIsoValuesMap"));
   mvaNoIsoValuesMapToken_      = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaNoIsoValuesMap"));
   //    mvaCategoriesMapToken_  = consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap"));
   // 
   //for stop analysis
   //    vetoIdFullInfoMapToken_ = consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleVetoCBIdMap"));
   //    mediumIdFullInfoMapToken_ = consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleMediumCBIdMap"));
   // 
   filterTriggerNames_     = iConfig.getUntrackedParameter<std::vector<std::string> >("filterTriggerNames");
   
   metSigToken_            = consumes<double>(iConfig.getParameter<edm::InputTag>("metSigInput"));
   metCovToken_            = consumes<math::Error<2>::type>(iConfig.getParameter<edm::InputTag>("metCovInput"));
   qgToken_                = consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"));
   
   genTTXJetsToken_                    = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genTTXJets"));
   genTTXBHadJetIndexToken_            = consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genTTXBHadJetIndex"));
   genTTXBHadFlavourToken_             = consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genTTXBHadFlavour"));
   genTTXBHadFromTopWeakDecayToken_    = consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genTTXBHadFromTopWeakDecay"));
   genTTXBHadPlusMothersToken_         = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genTTXBHadPlusMothers"));
   genTTXBHadPlusMothersIndicesToken_  = consumes<std::vector<std::vector<int> > >(iConfig.getParameter<edm::InputTag>("genTTXBHadPlusMothersIndices"));
   genTTXBHadIndexToken_               = consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genTTXBHadIndex"));
   genTTXBHadLeptonHadronIndexToken_   = consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genTTXBHadLeptonHadronIndex"));
   genTTXBHadLeptonViaTauToken_        = consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genTTXBHadLeptonViaTau"));
   genTTXCHadJetIndexToken_            = consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genTTXCHadJetIndex"));
   genTTXCHadFlavourToken_             = consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genTTXCHadFlavour"));
   genTTXCHadFromTopWeakDecayToken_    = consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genTTXCHadFromTopWeakDecay"));
   genTTXCHadBHadronIdToken_           = consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genTTXCHadBHadronId"));      
   
   // #########################
   // #  Read XML config file #
   // #########################
   //
   std::string confFile = iConfig.getParameter<std::string>("confFile");
   // ------
   // Read the config file and extract value used to apply preselection
   // Initially developed by Thibault - Modified by Eric
   // -----
   ReadConfFile(confFile);
   int buffersize = iConfig.getParameter<int>("bufferSize");
   if (buffersize <= 0) buffersize = 32000;
   
   xmlconf.LoadFile("conf.xml");
   XMLElement* tElement = xmlconf.FirstChildElement("variables")->FirstChildElement("var");
   
   for( XMLElement* child=tElement;child!=0;child=child->NextSiblingElement() )
     {
        std::string vname = child->ToElement()->Attribute("name");
        std::string vsave = child->ToElement()->Attribute("save");
        bool bsave = atoi(vsave.c_str());
        if( child->ToElement()->Attribute("mc") )
	  {
	     std::string vmc = child->ToElement()->Attribute("mc");
	     bool mc =  atoi(vmc.c_str());
	     if( isData_ && mc ) bsave = 0; // force the exclusion of mc-related variables when running on data
	  }
	
        ftree->conf.insert(std::make_pair(vname,bsave));
     }
   
   ftree->CreateBranches(buffersize);
   
   // ###############################
   // #  Add count & weight histos  #
   // ###############################
   //
   hcount = fs->make<TH1D>("hcount","hcount",1,0.,1.);
   hweight = fs->make<TH1D>("hweight","hweight",1,0.,1.);
}

FlatTreeProducer::~FlatTreeProducer()
{
}

void FlatTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   hcount->SetBinContent(1,hcount->GetBinContent(1)+1);
   
   ftree->Init();
   
   // Initial-state info
   edm::Handle<GenEventInfoProduct> genEventInfo;
   if( !isData_ ) iEvent.getByToken(genEventInfoToken_,genEventInfo);
   
   // LHE
   edm::Handle<LHEEventProduct> lheEventProduct;
   if( !isData_ && fillMCScaleWeight_ ) iEvent.getByToken(LHEEventProductToken_,lheEventProduct);
   
   // Gen particles
   edm::Handle<reco::GenParticleCollection> genParticlesHandle;
   if( !isData_ ) iEvent.getByToken(genParticlesToken_,genParticlesHandle);
   
   // Beamspot
   edm::Handle<reco::BeamSpot> bsHandle;
   iEvent.getByToken(bsToken_, bsHandle);
   const reco::BeamSpot &beamspot = *bsHandle.product();
   
   // Primary vertex
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vertexToken_,vertices);
   
   // Triggers
   edm::Handle<edm::TriggerResults> triggerBits;
   iEvent.getByToken(triggerBits_,triggerBits);
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   
   //    edm::Handle<bool> badMuonFilter;
   //    edm::Handle<bool> badChargedCandidateFilter;
   // 
   //    iEvent.getByToken(badMuonFilterToken_,badMuonFilter);
   //    iEvent.getByToken(badChargedCandidateFilterToken_,badChargedCandidateFilter);
   // 
   edm::Handle<edm::TriggerResults> triggerBitsPAT;
   iEvent.getByToken(triggerBitsPAT_,triggerBitsPAT);
   edm::TriggerNames namesPAT;   
   if( triggerBitsPAT.isValid() )
     {	
        namesPAT = iEvent.triggerNames(*triggerBitsPAT);
     }   
   
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerObjects_, triggerObjects);
   
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   iEvent.getByToken(triggerPrescales_, triggerPrescales);
   
   // Pile-up
   edm::Handle<std::vector< PileupSummaryInfo> > pileupInfo;
   if( !isData_ && fillPUInfo_ ) iEvent.getByToken(puInfoToken_,pileupInfo);
   
   // Rho info
   edm::Handle<double> rhoPtr;
   iEvent.getByToken(rhoToken_,rhoPtr);
   
   // MET significance
   edm::Handle<double> metSigPtr;
   iEvent.getByToken(metSigToken_,metSigPtr);
   std::vector<edm::Handle<double> > doubleVec;
   iEvent.getManyByType(doubleVec);
   for(unsigned int i=0;i<doubleVec.size();i++)
     {
        if(doubleVec[i].provenance()->moduleLabel() == "METSignificance")
	  metSigPtr = doubleVec[i];
     }
   if( !metSigPtr.isValid() or metSigPtr.failedToGet() )
     std::cerr << " Fail to access METSignificance branch " << std::endl;
   
   // MET covariance
   edm::Handle<math::Error<2>::type> metCovPtr;
   iEvent.getByToken(metCovToken_,metCovPtr);
   ftree->met_cov00 = (*metCovPtr)(0,0);
   ftree->met_cov10 = (*metCovPtr)(1,0);
   ftree->met_cov01 = (*metCovPtr)(0,1);
   ftree->met_cov11 = (*metCovPtr)(1,1);
   
   // Packed candidate collection
   edm::Handle<pat::PackedCandidateCollection> pfcands;
   if( dataFormat_ != "AOD" ) iEvent.getByToken(pfcandsToken_,pfcands);
   
   // Jets
   edm::Handle<pat::JetCollection> jets;
   iEvent.getByToken(jetToken_,jets);
   
   edm::Handle<edm::View<pat::Jet> > view_jets;
   iEvent.getByToken(viewJetToken_,view_jets);
   
   edm::Handle<edm::ValueMap<float> > qgHandle;
   iEvent.getByToken(qgToken_, qgHandle);
   
   // PuppiJets
   edm::Handle<pat::JetCollection> jetsPuppi;
   try
     {
        iEvent.getByToken(jetPuppiToken_,jetsPuppi);
     }
   catch (...) {;
   }
   
   // W-jets
   edm::Handle<pat::JetCollection> ak8jets;
   iEvent.getByToken(ak8jetToken_,ak8jets);
   
   // AK10: W-jets
   edm::Handle<pat::JetCollection> ak10jets;
   iEvent.getByToken(ak10jetToken_,ak10jets);
   
   // GenJets
   edm::Handle<reco::GenJetCollection> genJets;
   if( !isData_ ) iEvent.getByToken(genJetToken_,genJets);
   edm::Handle<reco::JetFlavourMatchingCollection> genJetFlavourMatching;
   if( !isData_ )
     {
        //	iEvent.getByLabel("genJetFlavour",genJetFlavourMatching);
     }
   
   // Muons
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_,muons);
   
   // Electrons
   edm::Handle<edm::View<reco::GsfElectron> > electrons;
   iEvent.getByToken(electronToken_,electrons);
   
   edm::Handle<pat::ElectronCollection> electronsPAT;
   iEvent.getByToken(electronPATToken_,electronsPAT);
   
   edm::Handle<edm::ValueMap<bool> > veto_cbid_decisions;
   edm::Handle<edm::ValueMap<bool> > loose_cbid_decisions;
   edm::Handle<edm::ValueMap<bool> > medium_cbid_decisions;
   edm::Handle<edm::ValueMap<bool> > tight_cbid_decisions;
   
   edm::Handle<edm::ValueMap<bool> > NoIso90_mvaid_decisions;
   edm::Handle<edm::ValueMap<bool> > NoIso80_mvaid_decisions;
   edm::Handle<edm::ValueMap<bool> > NoIsoLoose_mvaid_decisions;
   
   edm::Handle<edm::ValueMap<bool> > Iso90_mvaid_decisions;
   edm::Handle<edm::ValueMap<bool> > Iso80_mvaid_decisions;
   edm::Handle<edm::ValueMap<bool> > IsoLoose_mvaid_decisions;
   
   iEvent.getByToken(eleVetoCBIdMapToken_,veto_cbid_decisions);
   iEvent.getByToken(eleLooseCBIdMapToken_,loose_cbid_decisions);
   iEvent.getByToken(eleMediumCBIdMapToken_,medium_cbid_decisions);
   iEvent.getByToken(eleTightCBIdMapToken_,tight_cbid_decisions);
   
   iEvent.getByToken(ele90NoIsoMVAIdMapToken_,NoIso90_mvaid_decisions);
   iEvent.getByToken(ele80NoIsoMVAIdMapToken_,NoIso80_mvaid_decisions);
   iEvent.getByToken(eleLooseNoIsoMVAIdMapToken_,NoIsoLoose_mvaid_decisions);
   
   iEvent.getByToken(ele90IsoMVAIdMapToken_,Iso90_mvaid_decisions);
   iEvent.getByToken(ele80IsoMVAIdMapToken_,Iso80_mvaid_decisions);
   iEvent.getByToken(eleLooseIsoMVAIdMapToken_,IsoLoose_mvaid_decisions);
   
   edm::Handle<edm::ValueMap<float> > mvaIsoValues;
   edm::Handle<edm::ValueMap<float> > mvaNoIsoValues;
   //    edm::Handle<edm::ValueMap<int> > mvaCategories;
   // 
   iEvent.getByToken(mvaIsoValuesMapToken_,mvaIsoValues);
   iEvent.getByToken(mvaNoIsoValuesMapToken_,mvaNoIsoValues);
   //    iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);
   // 
   //    edm::Handle<edm::ValueMap<vid::CutFlowResult> > veto_id_cutflow_;
   //    edm::Handle<edm::ValueMap<vid::CutFlowResult> > medium_id_cutflow_;
   // 
   //    iEvent.getByToken(vetoIdFullInfoMapToken_, veto_id_cutflow_);
   //    iEvent.getByToken(mediumIdFullInfoMapToken_, medium_id_cutflow_);
   // 
   // Taus
   edm::Handle<pat::TauCollection> taus;
   iEvent.getByToken(tauToken_,taus);
   
   // Conversions info
   edm::Handle<reco::ConversionCollection> hConversions;
   if( dataFormat_ != "AOD" ) iEvent.getByToken(hConversionsToken_,hConversions);
   
   // HadronMatcher
   edm::Handle<reco::GenJetCollection> genTTXJets;
   iEvent.getByToken(genTTXJetsToken_,genTTXJets);
   
   edm::Handle<std::vector<int> >genTTXBHadFlavour;
   iEvent.getByToken(genTTXBHadFlavourToken_,genTTXBHadFlavour);
   
   edm::Handle<std::vector<int> > genTTXBHadJetIndex;
   iEvent.getByToken(genTTXBHadJetIndexToken_,genTTXBHadJetIndex);
   
   edm::Handle<std::vector<int> > genTTXBHadFromTopWeakDecay;
   iEvent.getByToken(genTTXBHadFromTopWeakDecayToken_,genTTXBHadFromTopWeakDecay);
   
   edm::Handle<std::vector<reco::GenParticle> > genTTXBHadPlusMothers;
   iEvent.getByToken(genTTXBHadPlusMothersToken_,genTTXBHadPlusMothers);
   
   edm::Handle<std::vector<std::vector<int> > > genTTXBHadPlusMothersIndices;
   iEvent.getByToken(genTTXBHadPlusMothersIndicesToken_,genTTXBHadPlusMothersIndices);
   
   edm::Handle<std::vector<int> > genTTXBHadIndex;
   iEvent.getByToken(genTTXBHadIndexToken_,genTTXBHadIndex);
   
   edm::Handle<std::vector<int> > genTTXBHadLeptonHadronIndex;
   iEvent.getByToken(genTTXBHadLeptonHadronIndexToken_,genTTXBHadLeptonHadronIndex);
   
   edm::Handle<std::vector<int> > genTTXBHadLeptonViaTau;
   iEvent.getByToken(genTTXBHadLeptonViaTauToken_, genTTXBHadLeptonViaTau);
   
   edm::Handle<std::vector<int> > genTTXCHadFlavour;
   iEvent.getByToken(genTTXCHadFlavourToken_,genTTXCHadFlavour);
   
   edm::Handle<std::vector<int> > genTTXCHadJetIndex;
   iEvent.getByToken(genTTXCHadJetIndexToken_,genTTXCHadJetIndex);
   
   edm::Handle<std::vector<int> > genTTXCHadFromTopWeakDecay;
   iEvent.getByToken(genTTXCHadFromTopWeakDecayToken_,genTTXCHadFromTopWeakDecay);
   
   edm::Handle<std::vector<int> > genTTXCHadBHadronId;
   iEvent.getByToken(genTTXCHadBHadronIdToken_,genTTXCHadBHadronId);
   
   // ###############################################################
   // #    ____                           _     _        __         #
   // #   / ___| ___ _ __   ___ _ __ __ _| |   (_)_ __  / _| ___    #
   // #  | |  _ / _ \ '_ \ / _ \ '__/ _` | |   | | '_ \| |_ / _ \   #
   // #  | |_| |  __/ | | |  __/ | | (_| | |   | | | | |  _| (_) |  #
   // #   \____|\___|_| |_|\___|_|  \__,_|_|   |_|_| |_|_|  \___/   #
   // #                                                             #
   // ###############################################################
   //
   ftree->ev_run = iEvent.id().run();
   ftree->ev_id = iEvent.id().event();
   
   ftree->ev_lumi = iEvent.id().luminosityBlock();
   
   //std::cout << " Event =================================================================== " << std::endl << "No: " << iEvent.id().event() << std::endl ;
   
   // ##########################################################
   // #   ___       _ _   _       _         _        _         #
   // #  |_ _|_ __ (_) |_(_) __ _| |    ___| |_ __ _| |_ ___   #
   // #   | || '_ \| | __| |/ _` | |   / __| __/ _` | __/ _ \  #
   // #   | || | | | | |_| | (_| | |   \__ \ || (_| | ||  __/  #
   // #  |___|_| |_|_|\__|_|\__,_|_|   |___/\__\__,_|\__\___|  #
   // #                                                        #
   // ##########################################################
   //
   float mc_weight = 1.;
   float mc_weight_originalValue = 1.;
   
   if( genEventInfo.isValid() )
     {
        // for some 94x amcatnlo samples, the value of the mc weight is not +-1, but some large (pos/neg) number
		// so we save both the +-1 label and the original number
		// the latter is needed for a proper reweighting of the weight_scale_* variables
		float wGen = genEventInfo->weight();
		mc_weight_originalValue = wGen;
        mc_weight = (wGen > 0 ) ? 1. : -1.;
		
	
        ftree->mc_id = genEventInfo->signalProcessID();
        ftree->mc_f1 = genEventInfo->pdf()->id.first;
        ftree->mc_f2 = genEventInfo->pdf()->id.second;
        ftree->mc_x1 = genEventInfo->pdf()->x.first;
        ftree->mc_x2 = genEventInfo->pdf()->x.second;
        ftree->mc_scale = genEventInfo->pdf()->scalePDF;
        if( genEventInfo->binningValues().size() > 0 ) ftree->mc_ptHat = genEventInfo->binningValues()[0];
     }

   if(! lheEventProduct.failedToGet())
     {
        if( !isData_ && fillMCScaleWeight_ )
	  {
	     ftree->weight_originalXWGTUP = lheEventProduct->originalXWGTUP();
	     
	     if( lheEventProduct->weights().size() > 0 )
	       {
		  ftree->weight_scale_muF0p5 = (genEventInfo->weight())*(lheEventProduct->weights()[2].wgt)/(lheEventProduct->originalXWGTUP()); // muF = 0.5 | muR = 1
		  ftree->weight_scale_muF2   = (genEventInfo->weight())*(lheEventProduct->weights()[1].wgt)/(lheEventProduct->originalXWGTUP()); // muF = 2   | muR = 1
		  ftree->weight_scale_muR0p5 = (genEventInfo->weight())*(lheEventProduct->weights()[6].wgt)/(lheEventProduct->originalXWGTUP()); // muF = 1   | muR = 0.5
		  ftree->weight_scale_muR2   = (genEventInfo->weight())*(lheEventProduct->weights()[3].wgt)/(lheEventProduct->originalXWGTUP()); // muF = 1   | muR = 2
	       }
	     
	     int nPdfAll = lheEventProduct->weights().size();
	     if( nPdf_ < nPdfAll && nPdf_ >= 0 ) nPdfAll = nPdf_;
	     for( int w=0;w<nPdfAll;w++ )
	       {
		  const LHEEventProduct::WGT& wgt = lheEventProduct->weights().at(w);
		  ftree->mc_pdfweights.push_back(wgt.wgt);
		  ftree->mc_pdfweightIds.push_back(wgt.id);
	       }	     
	  }
     }
   
   ftree->mc_weight = mc_weight;
   ftree->mc_weight_originalValue = mc_weight_originalValue;
   
   hweight->SetBinContent(1,hweight->GetBinContent(1)+ftree->mc_weight);
   
   // ####################################
   // #   ____  _ _                      #
   // #  |  _ \(_) | ___   _   _ _ __    #
   // #  | |_) | | |/ _ \ | | | | '_ \   #
   // #  |  __/| | |  __/ | |_| | |_) |  #
   // #  |_|   |_|_|\___|  \__,_| .__/   #
   // #                         |_|      #
   // #                                  #
   // ####################################
   //
   if( !isData_ && fillPUInfo_)
     {
        ftree->mc_pu_Npvi = pileupInfo->size();
        for(std::vector<PileupSummaryInfo>::const_iterator pvi=pileupInfo->begin();
	    pvi!=pileupInfo->end();pvi++)
	  {
	     signed int n_bc = pvi->getBunchCrossing();
	     ftree->mc_pu_BunchCrossing.push_back(n_bc);
	     if( n_bc == 0 )
	       {
		  ftree->mc_pu_intime_NumInt = pvi->getPU_NumInteractions();
		  ftree->mc_pu_trueNumInt = pvi->getTrueNumInteractions();
	       }
	     else if( n_bc == -1 ) ftree->mc_pu_before_npu = pvi->getPU_NumInteractions();
	     else if( n_bc == +1 ) ftree->mc_pu_after_npu  = pvi->getPU_NumInteractions();
	     
	     std::vector<float> mc_pu_zpositions;
	     std::vector<float> mc_pu_sumpT_lowpT;
	     std::vector<float> mc_pu_sumpT_highpT;
	     std::vector<int> mc_pu_ntrks_lowpT;
	     std::vector<int> mc_pu_ntrks_highpT;
	     
	     ftree->mc_pu_Nzpositions.push_back(pvi->getPU_zpositions().size());
	     for( unsigned int ipu=0;ipu<pvi->getPU_zpositions().size();ipu++ )
	       {
		  mc_pu_zpositions.push_back((pvi->getPU_zpositions())[ipu]);
		  mc_pu_sumpT_lowpT.push_back((pvi->getPU_sumpT_lowpT())[ipu]);
		  mc_pu_sumpT_highpT.push_back((pvi->getPU_sumpT_highpT())[ipu]);
		  mc_pu_ntrks_lowpT.push_back((pvi->getPU_ntrks_lowpT())[ipu]);
		  mc_pu_ntrks_highpT.push_back((pvi->getPU_ntrks_highpT())[ipu]);
	       }
	     
	     ftree->mc_pu_zpositions.push_back(mc_pu_zpositions);
	     ftree->mc_pu_sumpT_lowpT.push_back(mc_pu_sumpT_lowpT);
	     ftree->mc_pu_sumpT_highpT.push_back(mc_pu_sumpT_highpT);
	     ftree->mc_pu_ntrks_lowpT.push_back(mc_pu_ntrks_lowpT);
	     ftree->mc_pu_ntrks_highpT.push_back(mc_pu_ntrks_highpT);
	  }
     }
   
   GenTTXCategorizer *genTTX = new GenTTXCategorizer();
   
   if( !isData_ )
     {
	genTTX->Init(*ftree);
     }   
   
   // ##################################################
   // #   __  __  ____     _____           _   _       #
   // #  |  \/  |/ ___|   |_   _| __ _   _| |_| |__    #
   // #  | |\/| | |         | || '__| | | | __| '_ \   #
   // #  | |  | | |___      | || |  | |_| | |_| | | |  #
   // #  |_|  |_|\____|     |_||_|   \__,_|\__|_| |_|  #
   // #                                                #
   // ##################################################
   //
   bool do_mc_truth_tth = ftree->doWrite("mc_truth_tth");
   bool do_mc_truth_ttz = ftree->doWrite("mc_truth_ttz");
   bool do_mc_truth_ttw = ftree->doWrite("mc_truth_ttw");
   bool do_mc_truth_tzq = ftree->doWrite("mc_truth_tzq");
   bool do_mc_truth_thq = ftree->doWrite("mc_truth_thq");
   
   MCTruth *mc_truth = new MCTruth();
   
   bool reqMCTruth = 0;
   if( (
	do_mc_truth_tth ||
	do_mc_truth_tzq ||
	do_mc_truth_ttz ||
	do_mc_truth_ttw ||
	do_mc_truth_thq
       ) &&
       !isData_ )
     {
        mc_truth->Init(*ftree);
        if( do_mc_truth_tth ) mc_truth->fillTTHSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
        if( do_mc_truth_ttz ) mc_truth->fillTTZSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
        if( do_mc_truth_ttw ) mc_truth->fillTTWSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
        if( do_mc_truth_tzq ) mc_truth->fillTZQSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
        if( do_mc_truth_thq ) mc_truth->fillTHQSignalGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
        reqMCTruth = 1;
     }
   
   bool do_gen_all = ftree->doWrite("gen_all");
   bool do_gen_stop_mass = ftree->doWrite("gen_stop_mass");
   bool do_gen_stop = ftree->doWrite("gen_stop");
   
   if( do_gen_all &&
       !isData_ )
     {
        if( !reqMCTruth ) mc_truth->Init(*ftree);
        mc_truth->fillGenParticles(iEvent,iSetup,*ftree,genParticlesHandle);
	reqMCTruth = 1;
	if(do_gen_stop_mass)
	  mc_truth->fillStopNeutralinoMass(iEvent,iSetup,*ftree,genParticlesHandle);
     }
   if( !do_gen_all && do_gen_stop && !isData_){
      if(!reqMCTruth ) mc_truth->Init(*ftree);
      mc_truth->fillTopStopDecayChain(iEvent,iSetup,*ftree,genParticlesHandle);
      reqMCTruth = 1;
      
   }
   if( !isData_ )
     {
        if( !reqMCTruth ) mc_truth->Init(*ftree);
        mc_truth->fillGenPV(iEvent,iSetup,*ftree,genParticlesHandle);
        reqMCTruth = 1;
     }
   
   // #########################################
   // #   _____     _                         #
   // #  |_   _| __(_) __ _  __ _  ___ _ __   #
   // #    | || '__| |/ _` |/ _` |/ _ \ '__|  #
   // #    | || |  | | (_| | (_| |  __/ |     #
   // #    |_||_|  |_|\__, |\__, |\___|_|     #
   // #               |___/ |___/             #
   // #                                       #
   // #########################################
   //
   bool passMETFilters = 1;
   bool pass_HBHENoiseFilter = 1;
   bool pass_HBHENoiseIsoFilter = 1;
   bool pass_EcalDeadCellTriggerPrimitiveFilter = 1;
   bool pass_goodVertices = 1;
   bool pass_eeBadScFilter = 1;
   bool pass_globalTightHalo2016Filter = 1;
   bool pass_BadPFMuonFilter = 1;
   bool pass_BadChargedCandidateFilter = 1;
   bool pass_ecalBadCalibFilter = 1;
   
   //    bool pass_badMuonFilter = *badMuonFilter;
   //    bool pass_badChargedCandidateFilter = *badChargedCandidateFilter;
   // 
   // https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Moriond_2018
   // 
   if( triggerBitsPAT.isValid() )
     {	
        for (unsigned int i = 0, n = triggerBitsPAT->size(); i < n; ++i)
	  {
	     std::string triggerName = namesPAT.triggerName(i);
	     
	     bool isFired = (triggerBitsPAT->accept(i) ? true : false);
	     
	     if( strcmp(triggerName.c_str(),"Flag_HBHENoiseFilter") == 0 )
	       {
		  if( !isFired ) pass_HBHENoiseFilter = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_HBHENoiseIsoFilter") == 0 )
	       {
		  if( !isFired ) pass_HBHENoiseIsoFilter = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_EcalDeadCellTriggerPrimitiveFilter") == 0 )
	       {
		  if( !isFired ) pass_EcalDeadCellTriggerPrimitiveFilter = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_goodVertices") == 0 )
	       {
		  if( !isFired ) pass_goodVertices = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_eeBadScFilter") == 0 && isData_ == 1 )
	       {
		  if( !isFired ) pass_eeBadScFilter = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_globalTightHalo2016Filter") == 0 )
	       {
		  if( !isFired ) pass_globalTightHalo2016Filter = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_BadPFMuonFilter") == 0 )
	       {
		  if( !isFired ) pass_BadPFMuonFilter = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_BadChargedCandidateFilter") == 0 )
	       {
		  if( !isFired ) pass_BadChargedCandidateFilter = 0;
	       }
	     else if( strcmp(triggerName.c_str(),"Flag_ecalBadCalibFilter") == 0 )
	       {
		  if( !isFired ) pass_ecalBadCalibFilter = 0;
	       }
	  }
     }

   passMETFilters = (pass_HBHENoiseFilter &&
            pass_HBHENoiseIsoFilter &&
		     pass_EcalDeadCellTriggerPrimitiveFilter &&
		     pass_goodVertices && 
		     pass_eeBadScFilter &&
		     pass_globalTightHalo2016Filter &&
		     pass_BadPFMuonFilter &&
		     pass_BadChargedCandidateFilter &&
		     pass_ecalBadCalibFilter);
   
   //std::cout << "\n === TRIGGER PATHS === " << std::endl;
   for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
     {
        //std::cout << "[" << i << "] " << (triggerBits->accept(i) ? "1" : "0") << "  " << names.triggerName(i)  << std::endl;
	
        std::string triggerName = names.triggerName(i);
	
        if( !foundTrigger(triggerName) ) continue;
	
        ftree->trigger.push_back(i);
        ftree->trigger_name.push_back(triggerName);
        ftree->trigger_pass.push_back(triggerBits->accept(i) ? true : false);
        ftree->trigger_prescale.push_back(triggerPrescales->getPrescaleForIndex(i));
	
        float HLTprescale = 1.;
        float L1prescale = 1.;
	
        if( isData_ )
	  {
	     std::pair<std::vector<std::pair<std::string,int> >,int> detailedPrescaleInfo =
	       hltPrescale_.prescaleValuesInDetail(iEvent,iSetup,triggerName);
	     
	     HLTprescale = triggerPrescales.isValid() ? detailedPrescaleInfo.second : -1;
	     
	     std::vector <int> l1prescalevals;
	     for( size_t varind=0;varind<detailedPrescaleInfo.first.size();varind++ )
	       {
		  l1prescalevals.push_back(detailedPrescaleInfo.first.at(varind).second);
	       }
	     
	     // find and save minimum l1 prescale of any ORed L1 that seeds the HLT
	     std::vector<int>::iterator result = std::min_element(std::begin(l1prescalevals),
								  std::end(l1prescalevals));
	     size_t minind = std::distance(std::begin(l1prescalevals), result);
	     // sometimes there are no L1s associated with a HLT.
	     // In that case, this branch stores -1 for the l1prescale
	     L1prescale = minind < l1prescalevals.size() ? l1prescalevals.at(minind) : -1;
	     
	     //	     std::cout << HLTprescale << " " << L1prescale << std::endl;
	  }
	
        ftree->trigger_HLTprescale.push_back(HLTprescale);
        ftree->trigger_L1prescale.push_back(L1prescale);
     }
   
   if( ftree->doWrite("triggerobject_do") )
     {	
	//std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
	for (pat::TriggerObjectStandAlone obj : *triggerObjects)
	  {
	     // note: not "const &" since we want to call unpackPathNames
	     obj.unpackPathNames(names);
	     
	     // Trigger object basic informations (pt, eta, phi)
	     //std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
	     
	     ftree->triggerobject_pt.push_back(obj.pt());
	     ftree->triggerobject_eta.push_back(obj.eta());
	     ftree->triggerobject_phi.push_back(obj.phi());
	     
	     // Trigger object collection
	     //std::cout << "\t   Collection: " << obj.collection() << std::endl;
	     ftree->triggerobject_collection.push_back(obj.collection());
	     
	     // Trigger object type IDs
	     ftree->triggerobject_filterIds_n.push_back(obj.filterIds().size());
	     for (unsigned h = 0; h < obj.filterIds().size(); ++h)
	       {
		  ftree->triggerobject_isTriggerL1Mu.push_back(obj.filterIds()[h] == -81 ? true : false);
		  ftree->triggerobject_isTriggerL1NoIsoEG.push_back(obj.filterIds()[h] == -82 ? true : false);
		  ftree->triggerobject_isTriggerL1IsoEG.push_back(obj.filterIds()[h] == -83 ? true : false);
		  ftree->triggerobject_isTriggerL1CenJet.push_back(obj.filterIds()[h] == -84 ? true : false);
		  ftree->triggerobject_isTriggerL1ForJet.push_back(obj.filterIds()[h] == -85 ? true : false);
		  ftree->triggerobject_isTriggerL1TauJet.push_back(obj.filterIds()[h] == -86 ? true : false);
		  ftree->triggerobject_isTriggerL1ETM.push_back(obj.filterIds()[h] == -87 ? true : false);
		  ftree->triggerobject_isTriggerL1ETT.push_back(obj.filterIds()[h] == -88 ? true : false);
		  ftree->triggerobject_isTriggerL1HTT.push_back(obj.filterIds()[h] == -89 ? true : false);
		  ftree->triggerobject_isTriggerL1HTM.push_back(obj.filterIds()[h] == -90 ? true : false);
		  ftree->triggerobject_isTriggerL1JetCounts.push_back(obj.filterIds()[h] == -91 ? true : false);
		  ftree->triggerobject_isTriggerL1HfBitCounts.push_back(obj.filterIds()[h] == -92 ? true : false);
		  ftree->triggerobject_isTriggerL1HfRingEtSums.push_back(obj.filterIds()[h] == -93 ? true : false);
		  ftree->triggerobject_isTriggerL1TechTrig.push_back(obj.filterIds()[h] == -94 ? true : false);
		  ftree->triggerobject_isTriggerL1Castor.push_back(obj.filterIds()[h] == -95 ? true : false);
		  ftree->triggerobject_isTriggerL1BPTX.push_back(obj.filterIds()[h] == -96 ? true : false);
		  ftree->triggerobject_isTriggerL1GtExternal.push_back(obj.filterIds()[h] == -97 ? true : false);
		  
		  ftree->triggerobject_isHLT_TriggerPhoton.push_back(obj.filterIds()[h] == 81 ? true : false);
		  ftree->triggerobject_isHLT_TriggerElectron.push_back(obj.filterIds()[h] == 82 ? true : false);
		  ftree->triggerobject_isHLT_TriggerMuon.push_back(obj.filterIds()[h] == 83 ? true : false);
		  ftree->triggerobject_isHLT_TriggerTau.push_back(obj.filterIds()[h] == 84 ? true : false);
		  ftree->triggerobject_isHLT_TriggerJet.push_back(obj.filterIds()[h] == 85 ? true : false);
		  ftree->triggerobject_isHLT_TriggerBJet.push_back(obj.filterIds()[h] == 86 ? true : false);
		  ftree->triggerobject_isHLT_TriggerMET.push_back(obj.filterIds()[h] == 87 ? true : false);
		  ftree->triggerobject_isHLT_TriggerTET.push_back(obj.filterIds()[h] == 88 ? true : false);
		  ftree->triggerobject_isHLT_TriggerTHT.push_back(obj.filterIds()[h] == 89 ? true : false);
		  ftree->triggerobject_isHLT_TriggerMHT.push_back(obj.filterIds()[h] == 90 ? true : false);
		  ftree->triggerobject_isHLT_TriggerTrack.push_back(obj.filterIds()[h] == 91 ? true : false);
		  ftree->triggerobject_isHLT_TriggerCluster.push_back(obj.filterIds()[h] == 92 ? true : false);
		  ftree->triggerobject_isHLT_TriggerMETSig.push_back(obj.filterIds()[h] == 93 ? true : false);
		  ftree->triggerobject_isHLT_TriggerELongit.push_back(obj.filterIds()[h] == 94 ? true : false);
		  ftree->triggerobject_isHLT_TriggerMHTSig.push_back(obj.filterIds()[h] == 95 ? true : false);
		  ftree->triggerobject_isHLT_TriggerHLongit.push_back(obj.filterIds()[h] == 96 ? true : false);
		  
		  ftree->triggerobject_filterIds.push_back(obj.filterIds()[h]);
	       }
	     
	     // Trigger object filter
	     ftree->triggerobject_filterLabels_n.push_back(obj.filterLabels().size());
	     for (unsigned h = 0; h < obj.filterLabels().size(); ++h)
	       {
		  //std::cout << "FilterLabel: " << obj.filterLabels()[h] << std::endl;
		  ftree->triggerobject_filterLabels.push_back(obj.filterLabels()[h]);
	       }
	     
	     //std::cout << std::endl;
	     std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	     std::vector<std::string> pathNamesLast = obj.pathNames(true);
	     
	     // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
	     // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
	     // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
	     //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
	     ftree->triggerobject_pathNamesAll_n.push_back(pathNamesAll.size());
	     for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h)
	       {
		  bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
		  bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
		  bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
		  bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
		  //std::cout << "   " << pathNamesAll[h];
		  //if (isBoth) std::cout << "(L,3)" << std::endl;
		  //if (isL3 && !isBoth) std::cout << "(*,3)" << std::endl;
		  //if (isLF && !isBoth) std::cout << "(L,*)" << std::endl;
		  //if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)" << std::endl;
		  
		  ftree->triggerobject_pathNamesAll.push_back(pathNamesAll[h]);
		  ftree->triggerobject_pathNamesAll_isBoth.push_back(isBoth);
		  ftree->triggerobject_pathNamesAll_isL3.push_back(isL3);
		  ftree->triggerobject_pathNamesAll_isLF.push_back(isLF);
		  ftree->triggerobject_pathNamesAll_isNone.push_back(isNone);
	       }
	     
	     ftree->triggerobject_n = ftree->triggerobject_pt.size();
	  }
     }
   
   reco::Vertex *primVtx = NULL;
   
   ftree->nvertex = int(vertices->size());
   
   if( ! vertices->empty() )
     {
        const reco::Vertex &PV = vertices->front();
        primVtx = (reco::Vertex*)&PV;
     }
   
   if( primVtx )
     {
        ftree->pv_x = primVtx->position().x();
        ftree->pv_y = primVtx->position().y();
        ftree->pv_z = primVtx->position().z();
	
        ftree->pv_xError = primVtx->xError();
        ftree->pv_yError = primVtx->yError();
        ftree->pv_zError = primVtx->zError();

        ftree->pv_ndof = primVtx->chi2();
        ftree->pv_ndof = primVtx->ndof();
        ftree->pv_rho = primVtx->position().Rho();
        ftree->pv_isFake = primVtx->isFake();
     }
   
   // Rho
   ftree->ev_rho = *rhoPtr;
   
   // ####################################################
   // #   __  __ _         _               _____ _____   #
   // #  |  \/  (_)___ ___(_)_ __   __ _  | ____|_   _|  #
   // #  | |\/| | / __/ __| | '_ \ / _` | |  _|   | |    #
   // #  | |  | | \__ \__ \ | | | | (_| | | |___  | |    #
   // #  |_|  |_|_|___/___/_|_| |_|\__, | |_____| |_|    #
   // #                            |___/                 #
   // #                                                  #
   // ####################################################
   //
   // met significance
   ftree->met_sig = *metSigPtr;
   
   // MET
   edm::Handle<std::vector<pat::MET> > metAOD;
   edm::Handle<pat::METCollection> metMINIAOD;
   if( dataFormat_ != "AOD" )
     {
        iEvent.getByToken(metTokenMINIAOD_,metMINIAOD);
        const pat::MET &metv = metMINIAOD->front();
	
        ftree->met_px    = metv.px();
        ftree->met_py    = metv.py();
        ftree->met_pt    = metv.pt();
        ftree->met_phi   = metv.phi();
        ftree->met_sumet = metv.sumEt();
	
        if( !isData_ )
	  {
	     ftree->metGen_px    = metv.genMET()->px();
	     ftree->metGen_py    = metv.genMET()->py();
	     ftree->metGen_pt    = metv.genMET()->pt();
	     ftree->metGen_phi   = metv.genMET()->phi();
	     ftree->metGen_sumet = metv.genMET()->sumEt();
	     
	     ftree->metGen_NeutralEMEt  = metv.genMET()->NeutralEMEt();
	     ftree->metGen_ChargedEMEt  = metv.genMET()->ChargedEMEt();
	     ftree->metGen_NeutralHadEt = metv.genMET()->NeutralHadEt();
	     ftree->metGen_ChargedHadEt = metv.genMET()->ChargedHadEt();
	     ftree->metGen_MuonEt       = metv.genMET()->MuonEt();
	     ftree->metGen_InvisibleEt  = metv.genMET()->InvisibleEt();
	  }
	
        ftree->met_uncorrectedPt    = metv.uncorPt();
        ftree->met_uncorrectedPhi   = metv.uncorPhi();
        ftree->met_uncorrectedSumEt = metv.uncorSumEt();
	
        ftree->met_caloMETPt    = metv.caloMETPt();
        ftree->met_caloMETPhi   = metv.caloMETPhi();
        ftree->met_caloMETSumEt = metv.caloMETSumEt();

        pat::MET::METCorrectionLevel level = pat::MET::METCorrectionLevel::Type1;

        ftree->met_shiftedPx_JetEnUp = metv.shiftedPx(pat::MET::METUncertainty::JetEnUp,level);
        ftree->met_shiftedPx_JetEnDown = metv.shiftedPx(pat::MET::METUncertainty::JetEnDown,level);
        ftree->met_shiftedPx_JetResUp = metv.shiftedPx(pat::MET::METUncertainty::JetResUp,level);
        ftree->met_shiftedPx_JetResDown = metv.shiftedPx(pat::MET::METUncertainty::JetResDown,level);
        ftree->met_shiftedPx_MuonEnUp = metv.shiftedPx(pat::MET::METUncertainty::MuonEnUp,level);
        ftree->met_shiftedPx_MuonEnDown = metv.shiftedPx(pat::MET::METUncertainty::MuonEnDown,level);
        ftree->met_shiftedPx_ElectronEnUp = metv.shiftedPx(pat::MET::METUncertainty::ElectronEnUp,level);
        ftree->met_shiftedPx_ElectronEnDown = metv.shiftedPx(pat::MET::METUncertainty::ElectronEnDown,level);
        ftree->met_shiftedPx_TauEnUp = metv.shiftedPx(pat::MET::METUncertainty::TauEnUp,level);
        ftree->met_shiftedPx_TauEnDown = metv.shiftedPx(pat::MET::METUncertainty::TauEnDown,level);
        ftree->met_shiftedPx_UnclusteredEnUp = metv.shiftedPx(pat::MET::METUncertainty::UnclusteredEnUp,level);
        ftree->met_shiftedPx_UnclusteredEnDown = metv.shiftedPx(pat::MET::METUncertainty::UnclusteredEnDown,level);
        ftree->met_shiftedPx_NoShift = metv.shiftedPx(pat::MET::METUncertainty::NoShift,level);
        ftree->met_shiftedPx_PhotonEnUp = metv.shiftedPx(pat::MET::METUncertainty::PhotonEnUp,level);
        ftree->met_shiftedPx_PhotonEnDown = metv.shiftedPx(pat::MET::METUncertainty::PhotonEnDown,level);

        ftree->met_shiftedPy_JetEnUp = metv.shiftedPy(pat::MET::METUncertainty::JetEnUp,level);
        ftree->met_shiftedPy_JetEnDown = metv.shiftedPy(pat::MET::METUncertainty::JetEnDown,level);
        ftree->met_shiftedPy_JetResUp = metv.shiftedPy(pat::MET::METUncertainty::JetResUp,level);
        ftree->met_shiftedPy_JetResDown = metv.shiftedPy(pat::MET::METUncertainty::JetResDown,level);
        ftree->met_shiftedPy_MuonEnUp = metv.shiftedPy(pat::MET::METUncertainty::MuonEnUp,level);
        ftree->met_shiftedPy_MuonEnDown = metv.shiftedPy(pat::MET::METUncertainty::MuonEnDown,level);
        ftree->met_shiftedPy_ElectronEnUp = metv.shiftedPy(pat::MET::METUncertainty::ElectronEnUp,level);
        ftree->met_shiftedPy_ElectronEnDown = metv.shiftedPy(pat::MET::METUncertainty::ElectronEnDown,level);
        ftree->met_shiftedPy_TauEnUp = metv.shiftedPy(pat::MET::METUncertainty::TauEnUp,level);
        ftree->met_shiftedPy_TauEnDown = metv.shiftedPy(pat::MET::METUncertainty::TauEnDown,level);
        ftree->met_shiftedPy_UnclusteredEnUp = metv.shiftedPy(pat::MET::METUncertainty::UnclusteredEnUp,level);
        ftree->met_shiftedPy_UnclusteredEnDown = metv.shiftedPy(pat::MET::METUncertainty::UnclusteredEnDown,level);
        ftree->met_shiftedPy_NoShift = metv.shiftedPy(pat::MET::METUncertainty::NoShift,level);
        ftree->met_shiftedPy_PhotonEnUp = metv.shiftedPy(pat::MET::METUncertainty::PhotonEnUp,level);
        ftree->met_shiftedPy_PhotonEnDown = metv.shiftedPy(pat::MET::METUncertainty::PhotonEnDown,level);

        ftree->met_shiftedPhi_JetEnUp = metv.shiftedPhi(pat::MET::METUncertainty::JetEnUp,level);
        ftree->met_shiftedPhi_JetEnDown = metv.shiftedPhi(pat::MET::METUncertainty::JetEnDown,level);
        ftree->met_shiftedPhi_JetResUp = metv.shiftedPhi(pat::MET::METUncertainty::JetResUp,level);
        ftree->met_shiftedPhi_JetResDown = metv.shiftedPhi(pat::MET::METUncertainty::JetResDown,level);
        ftree->met_shiftedPhi_MuonEnUp = metv.shiftedPhi(pat::MET::METUncertainty::MuonEnUp,level);
        ftree->met_shiftedPhi_MuonEnDown = metv.shiftedPhi(pat::MET::METUncertainty::MuonEnDown,level);
        ftree->met_shiftedPhi_ElectronEnUp = metv.shiftedPhi(pat::MET::METUncertainty::ElectronEnUp,level);
        ftree->met_shiftedPhi_ElectronEnDown = metv.shiftedPhi(pat::MET::METUncertainty::ElectronEnDown,level);
        ftree->met_shiftedPhi_TauEnUp = metv.shiftedPhi(pat::MET::METUncertainty::TauEnUp,level);
        ftree->met_shiftedPhi_TauEnDown = metv.shiftedPhi(pat::MET::METUncertainty::TauEnDown,level);
        ftree->met_shiftedPhi_UnclusteredEnUp = metv.shiftedPhi(pat::MET::METUncertainty::UnclusteredEnUp,level);
        ftree->met_shiftedPhi_UnclusteredEnDown = metv.shiftedPhi(pat::MET::METUncertainty::UnclusteredEnDown,level);
        ftree->met_shiftedPhi_NoShift = metv.shiftedPhi(pat::MET::METUncertainty::NoShift,level);
        ftree->met_shiftedPhi_PhotonEnUp = metv.shiftedPhi(pat::MET::METUncertainty::PhotonEnUp,level);
        ftree->met_shiftedPhi_PhotonEnDown = metv.shiftedPhi(pat::MET::METUncertainty::PhotonEnDown,level);

        ftree->met_shiftedSumEt_JetEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::JetEnUp,level);
        ftree->met_shiftedSumEt_JetEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::JetEnDown,level);
        ftree->met_shiftedSumEt_JetResUp = metv.shiftedSumEt(pat::MET::METUncertainty::JetResUp,level);
        ftree->met_shiftedSumEt_JetResDown = metv.shiftedSumEt(pat::MET::METUncertainty::JetResDown,level);
        ftree->met_shiftedSumEt_MuonEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::MuonEnUp,level);
        ftree->met_shiftedSumEt_MuonEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::MuonEnDown,level);
        ftree->met_shiftedSumEt_ElectronEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::ElectronEnUp,level);
        ftree->met_shiftedSumEt_ElectronEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::ElectronEnDown,level);
        ftree->met_shiftedSumEt_TauEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::TauEnUp,level);
        ftree->met_shiftedSumEt_TauEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::TauEnDown,level);
        ftree->met_shiftedSumEt_UnclusteredEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::UnclusteredEnUp,level);
        ftree->met_shiftedSumEt_UnclusteredEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::UnclusteredEnDown,level);
        ftree->met_shiftedSumEt_NoShift = metv.shiftedSumEt(pat::MET::METUncertainty::NoShift,level);
        ftree->met_shiftedSumEt_PhotonEnUp = metv.shiftedSumEt(pat::MET::METUncertainty::PhotonEnUp,level);
        ftree->met_shiftedSumEt_PhotonEnDown = metv.shiftedSumEt(pat::MET::METUncertainty::PhotonEnDown,level);
     }
   else
     {
        iEvent.getByToken(metTokenAOD_,metAOD);
        const pat::MET &metv = metAOD->front();
	
        ftree->met_pt = metv.pt();
        ftree->met_phi = metv.phi();
        ftree->met_sumet = metv.sumEt();
     }

   // MET no HF
   edm::Handle<pat::METCollection> metNoHF;
   try
     {
        iEvent.getByToken(metTokenNoHF_,metNoHF);
     }
   catch (...) {;
   }
   if( metNoHF.isValid() )
     {
        const pat::MET &met = metNoHF->front();
	
        ftree->metNoHF_pt    = met.pt();
        ftree->metNoHF_phi   = met.phi();
        ftree->metNoHF_sumet = met.sumEt();
     }
   
   // MET Puppi
   edm::Handle<pat::METCollection> metPuppi;
   try
     {
        iEvent.getByToken(metTokenPuppi_,metPuppi);
     }
   catch (...) {;
   }
   if( metPuppi.isValid() )
     {
        const pat::MET &metv = metPuppi->front();
	
        ftree->metPuppi_pt    = metv.pt();
        ftree->metPuppi_phi   = metv.phi();
        ftree->metPuppi_sumet = metv.sumEt();
     }
   
   if( !isData_ )
     {
	genTTX->Run(*ftree,
		    genTTXJets,
		    genTTXBHadFlavour,
		    genTTXBHadJetIndex,
		    genTTXBHadFromTopWeakDecay,
		    genTTXBHadPlusMothers,
		    genTTXBHadPlusMothersIndices,
		    genTTXBHadIndex,
		    genTTXBHadLeptonHadronIndex,
		    genTTXBHadLeptonViaTau,
		    genTTXCHadFlavour,
		    genTTXCHadJetIndex,
		    genTTXCHadFromTopWeakDecay,
		    genTTXCHadBHadronId);
     }   
   
   // #################################################
   // #   _____ _           _                         #
   // #  | ____| | ___  ___| |_ _ __ ___  _ __  ___   #
   // #  |  _| | |/ _ \/ __| __| '__/ _ \| '_ \/ __|  #
   // #  | |___| |  __/ (__| |_| | | (_) | | | \__ \  #
   // #  |_____|_|\___|\___|\__|_|  \___/|_| |_|___/  #
   // #                                               #
   // #################################################
   //
   int nElec = electrons->size();
   
   for(int ie=0;ie<nElec;ie++)
     {
        const pat::Electron& elec = electronsPAT->at(ie);
	
        // Skimming electrons with pT < 5 GeV.
        //if (elec.pt() < 5) continue;
	
        ftree->el_pt.push_back(elec.pt());
        ftree->el_eta.push_back(elec.eta());
        ftree->el_phi.push_back(elec.phi());
        ftree->el_m.push_back(elec.mass());
        ftree->el_E.push_back(elec.energy());
        ftree->el_id.push_back(elec.pdgId());
        ftree->el_charge.push_back(elec.charge());
	
        ftree->el_isGsfCtfScPixChargeConsistent.push_back(elec.isGsfCtfScPixChargeConsistent());
        ftree->el_isGsfScPixChargeConsistent.push_back(elec.isGsfScPixChargeConsistent());
        ftree->el_hadronicOverEm.push_back(elec.hadronicOverEm());
	
        // IP
        const reco::GsfTrackRef gsfTrack = elec.gsfTrack();
        bool hasGsfTrack = ( gsfTrack.isNonnull() );
        ftree->el_hasGsfTrack.push_back(hasGsfTrack);
        ftree->el_gsfTrack_d0.push_back((hasGsfTrack) ? gsfTrack->d0() : -666);
        ftree->el_gsfTrack_z0.push_back((hasGsfTrack) ? gsfTrack->dz() : -666);
        ftree->el_gsfTrack_d0Error.push_back((hasGsfTrack) ? gsfTrack->d0Error() : -666);
        ftree->el_gsfTrack_z0Error.push_back((hasGsfTrack) ? gsfTrack->dzError() : -666);
        ftree->el_gsfTrack_PV_dxy.push_back((hasGsfTrack) ? gsfTrack->dxy(primVtx->position()) : -666);
        ftree->el_gsfTrack_PV_dz.push_back((hasGsfTrack) ? gsfTrack->dz(primVtx->position()) : -666);
        ftree->el_gsfTrack_RP_dxy.push_back((hasGsfTrack) ? gsfTrack->dxy(gsfTrack->referencePoint()) : -666);
        ftree->el_gsfTrack_RP_dz.push_back((hasGsfTrack) ? gsfTrack->dz(gsfTrack->referencePoint()) : -666);
        ftree->el_gsfTrack_BS_dxy.push_back((hasGsfTrack) ? gsfTrack->dxy(beamspot.position()) : -666);
        ftree->el_gsfTrack_BS_dz.push_back((hasGsfTrack) ? gsfTrack->dz(beamspot.position()) : -666);
        ftree->el_gsfTrack_dxyError.push_back((hasGsfTrack) ? gsfTrack->dxyError() : -666);
        ftree->el_gsfTrack_dzError.push_back((hasGsfTrack) ? gsfTrack->dzError() : -666);
        ftree->el_gsfTrack_normalizedChi2.push_back((hasGsfTrack) ? gsfTrack->normalizedChi2() : -666);

        ftree->el_ip3d.push_back(elec.dB(pat::Electron::PV3D));
        ftree->el_ip3dErr.push_back(elec.edB(pat::Electron::PV3D));
        ftree->el_ip2d.push_back(elec.dB(pat::Electron::PV2D));
        ftree->el_ip2dErr.push_back(elec.edB(pat::Electron::PV2D));
        ftree->el_ip3dBS.push_back(elec.dB(pat::Electron::BS3D));
        ftree->el_ip3dBSErr.push_back(elec.edB(pat::Electron::BS3D));
        ftree->el_ip2dBS.push_back(elec.dB(pat::Electron::BS2D));
        ftree->el_ip2dBSErr.push_back(elec.edB(pat::Electron::BS2D));

        // Energy cluster
        ftree->el_superCluster_eta.push_back(elec.superCluster()->eta());
        ftree->el_superCluster_phi.push_back(elec.superCluster()->phi());
        ftree->el_superCluster_energy.push_back(elec.superCluster()->energy());
        ftree->el_superCluster_rawEnergy.push_back(elec.superCluster()->rawEnergy());
        ftree->el_superCluster_preshowerEnergy.push_back(elec.superCluster()->preshowerEnergy());
        ftree->el_superCluster_etaWidth.push_back(elec.superCluster()->etaWidth());
        ftree->el_superCluster_phiWidth.push_back(elec.superCluster()->phiWidth());
        ftree->el_superCluster_preshowerEnergyPlane1.push_back(elec.superCluster()->preshowerEnergyPlane1());
        ftree->el_superCluster_preshowerEnergyPlane2.push_back(elec.superCluster()->preshowerEnergyPlane2());
        ftree->el_superCluster_positionR.push_back(elec.superCluster()->position().R());

        ftree->el_basicClustersSize.push_back(elec.basicClustersSize());
        ftree->el_e1x5.push_back(elec.e1x5());
        ftree->el_e5x5.push_back(elec.e5x5());
        ftree->el_e2x5Max.push_back(elec.e2x5Max());
        ftree->el_sigmaEtaEta.push_back(elec.sigmaEtaEta());
        ftree->el_sigmaIetaIeta.push_back(elec.sigmaIetaIeta());
        ftree->el_sigmaIphiIphi.push_back(elec.sigmaIphiIphi());
        ftree->el_sigmaIetaIphi.push_back(elec.sigmaIetaIphi());
        ftree->el_full5x5_sigmaIphiIphi.push_back(elec.full5x5_sigmaIphiIphi());
        ftree->el_full5x5_sigmaEtaEta.push_back(elec.full5x5_sigmaEtaEta());
        ftree->el_full5x5_sigmaIetaIeta.push_back(elec.full5x5_sigmaIetaIeta());
        ftree->el_full5x5_sigmaIetaIphi.push_back(elec.full5x5_sigmaIetaIphi());
        ftree->el_full5x5_r9.push_back(elec.full5x5_r9());
        ftree->el_full5x5_e1x5.push_back(elec.full5x5_e1x5());
        ftree->el_full5x5_e5x5.push_back(elec.full5x5_e5x5());
        ftree->el_full5x5_e2x5Max.push_back(elec.full5x5_e2x5Max());

        double OneMinusE1x5E5x5 = (elec.e5x5() != 0.) ? 1.-(elec.e1x5()/elec.e5x5()) : -1.;
        double full5x5_OneMinusE1x5E5x5 = (elec.full5x5_e5x5() != 0.) ? 1.-(elec.full5x5_e1x5()/elec.full5x5_e5x5()) : -1.;

        ftree->el_full5x5_OneMinusE1x5E5x5.push_back(full5x5_OneMinusE1x5E5x5);
        ftree->el_OneMinusE1x5E5x5.push_back(OneMinusE1x5E5x5);

        double IoEmIoP = (1.0/elec.ecalEnergy())-(1.0/elec.p());
        ftree->el_IoEmIoP.push_back(IoEmIoP);
        ftree->el_eleEoPout.push_back(elec.eEleClusterOverPout());
        double PreShowerOverRaw = elec.superCluster()->preshowerEnergy()/elec.superCluster()->rawEnergy();
        ftree->el_PreShowerOverRaw.push_back(PreShowerOverRaw);
        double ooEmooP = fabs(1.0/elec.ecalEnergy()-elec.eSuperClusterOverP()/elec.ecalEnergy());
        ftree->el_ooEmooP.push_back(ooEmooP);

        // Track hits
        const reco::HitPattern& pattern = gsfTrack->hitPattern();
        ftree->el_numberOfLostHits.push_back(pattern.numberOfLostHits(reco::HitPattern::HitCategory::MISSING_INNER_HITS));
        ftree->el_expectedMissingInnerHits.push_back(pattern.numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS));
        ftree->el_numberOfHits.push_back(pattern.numberOfAllHits(reco::HitPattern::HitCategory::MISSING_INNER_HITS));

        ftree->el_expectedMissingOuterHits.push_back(pattern.numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS));
        ftree->el_numberOfValidPixelHits.push_back(pattern.numberOfValidPixelHits());
        ftree->el_numberOfLostPixelHits.push_back(pattern.numberOfLostPixelHits(reco::HitPattern::TRACK_HITS));
        ftree->el_trackerLayersWithMeasurement.push_back(pattern.trackerLayersWithMeasurement());
        ftree->el_pixelLayersWithMeasurement.push_back(pattern.pixelLayersWithMeasurement());
        ftree->el_numberOfValidStripLayersWithMonoAndStereo.push_back(pattern.numberOfValidStripLayersWithMonoAndStereo());
        ftree->el_trackerLayersWithoutMeasurement.push_back(pattern.trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS));
	
        ftree->el_numberOfValidHits.push_back((hasGsfTrack) ? gsfTrack->numberOfValidHits() : -666);
        ftree->el_numberOfLostHitsDefault.push_back((hasGsfTrack) ? gsfTrack->numberOfLostHits() : -666);

        ftree->el_fbrem.push_back(elec.fbrem());

        ftree->el_deltaEtaSuperClusterTrackAtVtx.push_back(elec.deltaEtaSuperClusterTrackAtVtx());
        ftree->el_deltaPhiSuperClusterTrackAtVtx.push_back(elec.deltaPhiSuperClusterTrackAtVtx());
        ftree->el_deltaEtaSeedClusterTrackAtCalo.push_back(elec.deltaEtaSeedClusterTrackAtCalo());
        ftree->el_deltaPhiSeedClusterTrackAtCalo.push_back(elec.deltaPhiSeedClusterTrackAtCalo());
        ftree->el_eSuperClusterOverP.push_back(elec.eSuperClusterOverP());

        const auto el = electrons->ptrAt(ie);
        ftree->el_mvaIso.push_back((*mvaIsoValues)[el]);
        ftree->el_mvaNoIso.push_back((*mvaNoIsoValues)[el]);
	//        ftree->el_mvaNonTrigCat.push_back((*mvaCategories)[el]);
	// 
        ftree->el_NoIso90MVAId.push_back((*NoIso90_mvaid_decisions)[el]);
        ftree->el_NoIso80MVAId.push_back((*NoIso80_mvaid_decisions)[el]);
        ftree->el_NoIsoLooseMVAId.push_back((*NoIsoLoose_mvaid_decisions)[el]);
	
        ftree->el_Iso90MVAId.push_back((*Iso90_mvaid_decisions)[el]);
        ftree->el_Iso80MVAId.push_back((*Iso80_mvaid_decisions)[el]);
        ftree->el_IsoLooseMVAId.push_back((*IsoLoose_mvaid_decisions)[el]);
	
        ftree->el_vetoCBId.push_back((*veto_cbid_decisions)[el]);
        ftree->el_looseCBId.push_back((*loose_cbid_decisions)[el]);
        ftree->el_mediumCBId.push_back((*medium_cbid_decisions)[el]);
        ftree->el_tightCBId.push_back((*tight_cbid_decisions)[el]);
	
        //for stop analysis
	//        vid::CutFlowResult vetoIdIsoMasked = (*veto_id_cutflow_)[el].getCutFlowResultMasking("GsfEleEffAreaPFIsoCut_0");
	//        ftree->el_vetoStopID.push_back(vetoIdIsoMasked.cutFlowPassed());
	//        vid::CutFlowResult mediumIdIsoMasked = (*medium_id_cutflow_)[el].getCutFlowResultMasking("GsfEleEffAreaPFIsoCut_0");
	//        ftree->el_mediumStopID.push_back(mediumIdIsoMasked.cutFlowPassed());
        //end
	
        ftree->el_ecalEnergy.push_back(elec.ecalEnergy());
        ftree->el_correctedEcalEnergy.push_back(elec.correctedEcalEnergy());
        ftree->el_correctedEcalEnergyError.push_back(elec.correctedEcalEnergyError());
        ftree->el_trackMomentumError.push_back(elec.trackMomentumError());
	
        ftree->el_hcalOverEcal.push_back(elec.hcalOverEcal());
        ftree->el_hcalOverEcalBc.push_back(elec.hcalOverEcalBc());
        ftree->el_hcalDepth1OverEcal.push_back(elec.hcalDepth1OverEcal());
        ftree->el_hcalDepth2OverEcal.push_back(elec.hcalDepth2OverEcal());
        ftree->el_eSeedClusterOverPout.push_back(elec.eSeedClusterOverPout());
        ftree->el_eSeedClusterOverP.push_back(elec.eSeedClusterOverP());
        ftree->el_eEleClusterOverPout.push_back(elec.eEleClusterOverPout());
        ftree->el_deltaEtaEleClusterTrackAtCalo.push_back(elec.deltaEtaEleClusterTrackAtCalo());
        ftree->el_deltaPhiEleClusterTrackAtCalo.push_back(elec.deltaPhiEleClusterTrackAtCalo());
	
        ftree->el_neutralHadronIso.push_back(elec.neutralHadronIso());
        ftree->el_chargedHadronIso.push_back(elec.chargedHadronIso());
        ftree->el_puChargedHadronIso.push_back(elec.puChargedHadronIso());
        ftree->el_ecalIso.push_back(elec.ecalIso());
        ftree->el_hcalIso.push_back(elec.hcalIso());
        ftree->el_particleIso.push_back(elec.particleIso());
        ftree->el_photonIso.push_back(elec.photonIso());
        ftree->el_trackIso.push_back(elec.trackIso());
	
        ftree->el_pfIso_sumChargedHadronPt.push_back(elec.pfIsolationVariables().sumChargedHadronPt);
        ftree->el_pfIso_sumNeutralHadronEt.push_back(elec.pfIsolationVariables().sumNeutralHadronEt);
        ftree->el_pfIso_sumPhotonEt.push_back(elec.pfIsolationVariables().sumPhotonEt);
        ftree->el_pfIso_sumPUPt.push_back(elec.pfIsolationVariables().sumPUPt);

        ftree->el_ecalPFClusterIso.push_back(elec.ecalPFClusterIso());
        ftree->el_hcalPFClusterIso.push_back(elec.hcalPFClusterIso());

        ftree->el_dr03EcalRecHitSumEt.push_back(elec.dr03EcalRecHitSumEt());
        ftree->el_dr03HcalTowerSumEt.push_back(elec.dr03HcalTowerSumEt());
        ftree->el_dr03HcalDepth1TowerSumEt.push_back(elec.dr03HcalDepth1TowerSumEt());
        ftree->el_dr03HcalDepth2TowerSumEt.push_back(elec.dr03HcalDepth2TowerSumEt());
        ftree->el_dr03TkSumPt.push_back(elec.dr03TkSumPt());

        ftree->el_dr04EcalRecHitSumEt.push_back(elec.dr04EcalRecHitSumEt());
        ftree->el_dr04HcalTowerSumEt.push_back(elec.dr04HcalTowerSumEt());
        ftree->el_dr04HcalDepth1TowerSumEt.push_back(elec.dr04HcalDepth1TowerSumEt());
        ftree->el_dr04HcalDepth2TowerSumEt.push_back(elec.dr04HcalDepth2TowerSumEt());
        ftree->el_dr04TkSumPt.push_back(elec.dr04TkSumPt());

        // mini-iso
        float miniIso           = -666;
        float miniIsoTTH        = -666;
        float miniIsoTTHCharged = -666;
        float miniIsoTTHNeutral = -666;
	
        if( dataFormat_ != "AOD" )
	  {
	     // Variables used for stop analysis
	     double EA_miniIso = getEA(elec.eta(), EL_EA_ETA, EL_EA_VALUE);
	     miniIso = getPFIsolation(pfcands,dynamic_cast<const reco::Candidate*>(&elec),*rhoPtr, EA_miniIso, 0.05,0.2,10.,false, false);
	     // ---------------------------------
	     //
	     // Below: Variables used for tt-H analysis
	     //
	     float miniIsoR = 10.0/std::min(std::max(float(elec.pt()),float(50.)),float(200.));
	     
	     float EffArea = 0.;
	     float eta = elec.superCluster()->eta();
	     
	     if(      fabs(eta) > 0      && fabs(eta) < 1.0 )   EffArea = 0.1566;
	     else if( fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffArea = 0.1626;
	     else if( fabs(eta) >= 1.479 && fabs(eta) < 2.0 )   EffArea = 0.1073;
	     else if( fabs(eta) >= 2.0   && fabs(eta) < 2.2 )   EffArea = 0.0854;
	     else if( fabs(eta) >= 2.2   && fabs(eta) < 2.3 )   EffArea = 0.1051;
	     else if( fabs(eta) >= 2.3   && fabs(eta) < 2.4 )   EffArea = 0.1204;
	     else if( fabs(eta) >= 2.4   && fabs(eta) < 2.5 )   EffArea = 0.1524;
	     
	     float correction = ftree->ev_rho*EffArea*(miniIsoR/0.3)*(miniIsoR/0.3);
	     
	     float pfIsoCharged = ElecPfIsoCharged(elec,pfcands,miniIsoR);
	     float pfIsoNeutral = ElecPfIsoNeutral(elec,pfcands,miniIsoR);
	     float pfIsoPUSubtracted = std::max(float(0.0),float(pfIsoNeutral-correction));
	     
	     miniIsoTTH        = (pfIsoCharged + pfIsoPUSubtracted)/elec.pt();
	     miniIsoTTHCharged = pfIsoCharged / elec.pt();
	     miniIsoTTHNeutral = pfIsoPUSubtracted / elec.pt();
	     //miniIsoTTHNeutral = pfIsoNeutral;
	     // ---------------------------------
	  }
	
        //std::cout << elec.pt() << "   " << miniIsoTTH << "   " << miniIsoTTHCharged << "  " << miniIsoTTHNeutral << std::endl;
	
        ftree->el_miniIso.push_back(miniIso);
        ftree->el_miniIsoTTH.push_back(miniIsoTTH);
	
        //ftree->el_miniIsoTTHCharged.pushback(miniIsoTTHCharged); //?
        //ftree->el_miniIsoTTHNeutral.pushback(miniIsoTTHNeutral); //?
	
        ftree->el_vx.push_back(elec.vx());
        ftree->el_vy.push_back(elec.vy());
        ftree->el_vz.push_back(elec.vz());
	
        ftree->el_hasGsfTrack.push_back(hasGsfTrack);

        ftree->el_passConversionVeto.push_back(elec.passConversionVeto());
	
        double el_pt = elec.pt();
        double el_eta = elec.eta();
        double el_phi = elec.phi();
        double el_lepMVA = -666.;
	
        float drmin = 0.5;
        int jcl = -1;
        for(unsigned int ij=0;ij<jets->size();ij++)
	  {
	     if( jets->at(ij).pt() < 10. ) continue;
	     float dr = GetDeltaR(jets->at(ij).eta(),
				  jets->at(ij).phi(),
				  el_eta,
				  el_phi);
	     if( dr < drmin )
	       {
		  drmin = dr;
		  jcl = ij;
	       }
	  }
	
        lepMVA_pt = el_pt; 
        lepMVA_eta = el_eta;
        lepMVA_miniRelIsoNeutral = miniIsoTTHNeutral;
        lepMVA_miniRelIsoCharged = miniIsoTTHCharged;
        lepMVA_jetPtRatio = (jcl >= 0) ? std::min(ptRatioElec(elec,jets->at(jcl)),1.5) : 1.5;
        lepMVA_jetPtRelv2 = (jcl >= 0) ? ptRelElec(elec,jets->at(jcl)) : 0.0;
        float csv = (jcl >= 0) ? jets->at(jcl).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -666;
        lepMVA_jetBTagCSV = std::max(double(csv),0.);
        lepMVA_sip3d = fabs(ftree->el_ip3d.back()/ftree->el_ip3dErr.back());
        lepMVA_dxy = log(fabs(ftree->el_gsfTrack_PV_dxy.back()));
        lepMVA_dz = log(fabs(ftree->el_gsfTrack_PV_dz.back()));
        lepMVA_mvaId = ftree->el_NoIsoLooseMVAId.back();
        lepMVA_jetNDauChargedMVASel = (jcl >= 0) ? jetNDauChargedMVASel(jets->at(jcl), *primVtx) : 0.0;

        el_lepMVA = ele_reader->EvaluateMVA("BDTG method");

        ftree->el_lepMVA.push_back(el_lepMVA);
	
        float conept = (jcl >= 0) ? conePtElec(elec,jets->at(jcl)) : -100;
        ftree->el_conept.push_back( conept );
	
        bool debug = false;
	
        if(debug)// && ( ftree->ev_id == 892573) )
	  {
	     std::cout << "el["              << ie                      << "]: "
	       << " pt= "            << lepMVA_pt
	       << " eta= "           << lepMVA_eta
	       << " miniIsoN= "      << lepMVA_miniRelIsoNeutral
	       << " miniIsoC= "      << lepMVA_miniRelIsoCharged
	       << " jetPtRatio= "    << lepMVA_jetPtRatio
	       << " jetPtRel= "      << lepMVA_jetPtRelv2
                << " ellepMVA= "      << el_lepMVA                      << std::endl;
	  }

        if(false)
	  {
	     std::cout << ftree->ev_id                   << " "
	       << lepMVA_pt                      << " "
	       << lepMVA_eta                     << " "
	       << elec.phi()                     << " "
	       << elec.energy()                  << " "
	       << elec.pdgId()                   << " "
	       << elec.charge()                  << " "
	       << lepMVA_jetNDauChargedMVASel    << " "
	       << miniIsoTTH                     << " "
	       << lepMVA_miniRelIsoCharged       << " "
	       << lepMVA_miniRelIsoNeutral       << " "
	       << lepMVA_jetPtRelv2              << " "
	       << lepMVA_jetBTagCSV              << " "
	       << lepMVA_jetPtRatio              << " "
	       << lepMVA_sip3d                   << " "
	       << fabs(ftree->el_gsfTrack_PV_dxy.back())                     << " "
	       << fabs(ftree->el_gsfTrack_PV_dz.back())                      << " "
	       << ftree->el_NoIsoLooseMVAId.back()                              << " "
	       << el_lepMVA
	       << std::endl;
	  }
	
        if( !isData_ )
	  {
	     // Internal matching
	     reco::GenParticle *genp = new reco::GenParticle();
	     
	     float drmin;
	     bool hasMCMatch = mc_truth->doMatch(iEvent,iSetup,genParticlesHandle,genp,drmin,
						 elec.pt(),elec.eta(),elec.phi(),elec.pdgId());
	     ftree->el_hasMCMatch.push_back(hasMCMatch);
	     if( hasMCMatch )
	       {
		  ftree->el_gen_pt.push_back(genp->pt());
		  ftree->el_gen_eta.push_back(genp->eta());
		  ftree->el_gen_phi.push_back(genp->phi());
		  ftree->el_gen_m.push_back(genp->mass());
		  ftree->el_gen_status.push_back(genp->status());
		  ftree->el_gen_id.push_back(genp->pdgId());
		  ftree->el_gen_charge.push_back(genp->charge());
		  ftree->el_gen_dr.push_back(drmin);
	       }
	     else
	       {
		  ftree->el_gen_pt.push_back(-666);
		  ftree->el_gen_eta.push_back(-666);
		  ftree->el_gen_phi.push_back(-666);
		  ftree->el_gen_m.push_back(-666);
		  ftree->el_gen_status.push_back(-666);
		  ftree->el_gen_id.push_back(-666);
		  ftree->el_gen_charge.push_back(-666);
		  ftree->el_gen_dr.push_back(-666);
	       }
	     delete genp;
	     
	     // PAT matching
	     const reco::GenParticle *genpPAT = elec.genParticle();
	     bool hasMCMatchPAT = (genpPAT != 0);
	     ftree->el_hasMCMatchPAT.push_back(hasMCMatchPAT);
	     if( hasMCMatchPAT )
	       {
		  ftree->el_genPAT_pt.push_back(genpPAT->pt());
		  ftree->el_genPAT_eta.push_back(genpPAT->eta());
		  ftree->el_genPAT_phi.push_back(genpPAT->phi());
		  ftree->el_genPAT_m.push_back(genpPAT->mass());
		  ftree->el_genPAT_status.push_back(genpPAT->status());
		  ftree->el_genPAT_id.push_back(genpPAT->pdgId());
		  ftree->el_genPAT_charge.push_back(genpPAT->charge());
	       }
	     else
	       {
		  ftree->el_genPAT_pt.push_back(-666);
		  ftree->el_genPAT_eta.push_back(-666);
		  ftree->el_genPAT_phi.push_back(-666);
		  ftree->el_genPAT_m.push_back(-666);
		  ftree->el_genPAT_status.push_back(-666);
		  ftree->el_genPAT_id.push_back(-666);
		  ftree->el_genPAT_charge.push_back(-666);
	       }
	  }
	
        ftree->el_lepMVA_pt.push_back(lepMVA_pt);
        ftree->el_lepMVA_eta.push_back(lepMVA_eta);
        ftree->el_lepMVA_miniRelIsoCharged.push_back(lepMVA_miniRelIsoCharged);
        ftree->el_lepMVA_miniRelIsoNeutral.push_back(lepMVA_miniRelIsoNeutral);
        ftree->el_lepMVA_jetPtRatio.push_back(lepMVA_jetPtRatio);
        ftree->el_lepMVA_jetPtRelv2.push_back(lepMVA_jetPtRelv2);
        ftree->el_lepMVA_jetBTagCSV.push_back(lepMVA_jetBTagCSV);
        ftree->el_lepMVA_sip3d.push_back(lepMVA_sip3d);
        ftree->el_lepMVA_dxy.push_back(lepMVA_dxy);
        ftree->el_lepMVA_dz.push_back(lepMVA_dz);
        ftree->el_lepMVA_mvaId.push_back(lepMVA_mvaId);
        ftree->el_lepMVA_jetNDauChargedMVASel.push_back(lepMVA_jetNDauChargedMVASel);
	
        bool allowCkfMatch = true;
        float lxyMin = 2.0;
        float probMin = 1e-6;
        uint nHitsBeforeVtxMax = 0;
	
        if( dataFormat_ != "AOD" )
	  {
	     bool matchConv = 0;
	     matchConv = ConversionTools::hasMatchedConversion(elec,hConversions,beamspot.position(),allowCkfMatch,lxyMin,probMin,nHitsBeforeVtxMax);
	     ftree->el_hasMatchedConversion.push_back(matchConv);
	  }
     }
   ftree->el_n = ftree->el_pt.size();
   
   // ####################################
   // #   __  __                         #
   // #  |  \/  |_   _  ___  _ __  ___   #
   // #  | |\/| | | | |/ _ \| '_ \/ __|  #
   // #  | |  | | |_| | (_) | | | \__ \  #
   // #  |_|  |_|\__,_|\___/|_| |_|___/  #
   // #                                  #
   // ####################################
   //
   int nMuon = muons->size();
   for(int im=0;im<nMuon;im++)
     {
        const pat::Muon& muon = muons->at(im);
	
        // Skimming muons with pT < 5 GeV.
        //if (muon.pt() < 5) continue;
	
        ftree->mu_pt.push_back(muon.pt());
        ftree->mu_eta.push_back(muon.eta());
        ftree->mu_phi.push_back(muon.phi());
        ftree->mu_m.push_back(muon.mass());
        ftree->mu_E.push_back(muon.energy());
        ftree->mu_id.push_back(muon.pdgId());
        ftree->mu_charge.push_back(muon.charge());
	
        // IP
        ftree->mu_ip3d.push_back(muon.dB(pat::Muon::PV3D));
        ftree->mu_ip3dErr.push_back(muon.edB(pat::Muon::PV3D));
        ftree->mu_ip2d.push_back(muon.dB(pat::Muon::PV2D));
        ftree->mu_ip2dErr.push_back(muon.edB(pat::Muon::PV2D));
        ftree->mu_ip3dBS.push_back(muon.dB(pat::Muon::BS3D));
        ftree->mu_ip3dBSErr.push_back(muon.edB(pat::Muon::BS3D));
        ftree->mu_ip2dBS.push_back(muon.dB(pat::Muon::BS2D));
        ftree->mu_ip2dBSErr.push_back(muon.edB(pat::Muon::BS2D));

        const reco::MuonQuality combQuality = muon.combinedQuality();
        ftree->mu_combinedQuality_chi2LocalPosition.push_back(combQuality.chi2LocalPosition);
        ftree->mu_combinedQuality_trkKink.push_back(combQuality.trkKink);
	
        ftree->mu_numberOfMatches.push_back(muon.isMatchesValid() ? muon.numberOfMatches() : -666);
        ftree->mu_numberOfMatchedStations.push_back(muon.numberOfMatchedStations());

        // GlobalTrack
        const reco::TrackRef globalTrack = muon.globalTrack();
        bool hasGlobalTrack = globalTrack.isNonnull();

        ftree->mu_hasGlobalTrack.push_back(hasGlobalTrack);
        ftree->mu_globalTrack_d0.push_back((hasGlobalTrack) ? globalTrack->d0() : -666);
        ftree->mu_globalTrack_z0.push_back((hasGlobalTrack) ? globalTrack->dz() : -666);
        ftree->mu_globalTrack_d0Error.push_back((hasGlobalTrack) ? globalTrack->d0Error() : -666);
        ftree->mu_globalTrack_z0Error.push_back((hasGlobalTrack) ? globalTrack->dzError() : -666);
        ftree->mu_globalTrack_PV_dxy.push_back((hasGlobalTrack) ? globalTrack->dxy(primVtx->position()) : -666);
        ftree->mu_globalTrack_PV_dz.push_back((hasGlobalTrack) ? globalTrack->dz(primVtx->position()) : -666);
        ftree->mu_globalTrack_RP_dxy.push_back((hasGlobalTrack) ? globalTrack->dxy(globalTrack->referencePoint()) : -666);
        ftree->mu_globalTrack_RP_dz.push_back((hasGlobalTrack) ? globalTrack->dz(globalTrack->referencePoint()) : -666);
        ftree->mu_globalTrack_BS_dxy.push_back((hasGlobalTrack) ? globalTrack->dxy(beamspot.position()) : -666);
        ftree->mu_globalTrack_BS_dz.push_back((hasGlobalTrack) ? globalTrack->dz(beamspot.position()) : -666);
        ftree->mu_globalTrack_dxyError.push_back((hasGlobalTrack) ? globalTrack->dxyError() : -666);
        ftree->mu_globalTrack_dzError.push_back((hasGlobalTrack) ? globalTrack->dzError() : -666);
        ftree->mu_globalTrack_normalizedChi2.push_back((hasGlobalTrack) ? globalTrack->normalizedChi2() : -666);
        ftree->mu_globalTrack_numberOfValidHits.push_back((hasGlobalTrack) ? globalTrack->numberOfValidHits() : -666);
        ftree->mu_globalTrack_numberOfValidMuonHits.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfValidMuonHits() : -666);
        ftree->mu_globalTrack_numberOfLostHits.push_back((hasGlobalTrack) ? globalTrack->numberOfLostHits() : -666);
        ftree->mu_globalTrack_pt.push_back((hasGlobalTrack) ? globalTrack->pt() : -666);
        ftree->mu_globalTrack_eta.push_back((hasGlobalTrack) ? globalTrack->eta() : -666);
        ftree->mu_globalTrack_phi.push_back((hasGlobalTrack) ? globalTrack->phi() : -666);
        ftree->mu_globalTrack_ptError.push_back((hasGlobalTrack) ? globalTrack->ptError() : -666);
        ftree->mu_globalTrack_etaError.push_back((hasGlobalTrack) ? globalTrack->etaError() : -666);
        ftree->mu_globalTrack_phiError.push_back((hasGlobalTrack) ? globalTrack->phiError() : -666);
        ftree->mu_globalTrack_vx.push_back((hasGlobalTrack) ? globalTrack->vx() : -666);
        ftree->mu_globalTrack_vy.push_back((hasGlobalTrack) ? globalTrack->vy() : -666);
        ftree->mu_globalTrack_vz.push_back((hasGlobalTrack) ? globalTrack->vz() : -666);
        ftree->mu_globalTrack_qoverp.push_back((hasGlobalTrack) ? globalTrack->qoverp() : -666);
        ftree->mu_globalTrack_qoverpError.push_back((hasGlobalTrack) ? globalTrack->qoverpError() : -666);
        ftree->mu_globalTrack_charge.push_back((hasGlobalTrack) ? globalTrack->charge() : -666);
        ftree->mu_globalTrack_trackerLayersWithMeasurement.push_back((hasGlobalTrack) ? globalTrack->hitPattern().trackerLayersWithMeasurement() : -666);
        ftree->mu_globalTrack_pixelLayersWithMeasurement.push_back((hasGlobalTrack) ? globalTrack->hitPattern().pixelLayersWithMeasurement() : -666);
        ftree->mu_globalTrack_numberOfValidStripLayersWithMonoAndStereo.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfValidStripLayersWithMonoAndStereo() : -666);
        ftree->mu_globalTrack_trackerLayersWithoutMeasurement.push_back((hasGlobalTrack) ? globalTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_globalTrack_numberOfValidPixelHits.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfValidPixelHits() : -666);
        ftree->mu_globalTrack_numberOfLostPixelHits.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_globalTrack_numberOfInnerHits.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) : -666);
        ftree->mu_globalTrack_numberOfOuterHits.push_back((hasGlobalTrack) ? globalTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS) : -666);
        ftree->mu_globalTrack_validFraction.push_back((hasGlobalTrack) ? globalTrack->validFraction() : -666);

        // BestTrack
        ftree->mu_bestTrackType.push_back(muon.muonBestTrackType());
        const reco::TrackRef bestTrack = muon.muonBestTrack();
        bool hasBestTrack = bestTrack.isNonnull();

        ftree->mu_hasBestTrack.push_back(hasBestTrack);
        ftree->mu_bestTrack_d0.push_back((hasBestTrack) ? bestTrack->d0() : -666);
        ftree->mu_bestTrack_z0.push_back((hasBestTrack) ? bestTrack->dz() : -666);
        ftree->mu_bestTrack_d0Error.push_back((hasBestTrack) ? bestTrack->d0Error() : -666);
        ftree->mu_bestTrack_z0Error.push_back((hasBestTrack) ? bestTrack->dzError() : -666);
        ftree->mu_bestTrack_PV_dxy.push_back((hasBestTrack) ? bestTrack->dxy(primVtx->position()) : -666);
        ftree->mu_bestTrack_PV_dz.push_back((hasBestTrack) ? bestTrack->dz(primVtx->position()) : -666);
        ftree->mu_bestTrack_RP_dxy.push_back((hasBestTrack) ? bestTrack->dxy(bestTrack->referencePoint()) : -666);
        ftree->mu_bestTrack_RP_dz.push_back((hasBestTrack) ? bestTrack->dz(bestTrack->referencePoint()) : -666);
        ftree->mu_bestTrack_BS_dxy.push_back((hasBestTrack) ? bestTrack->dxy(beamspot.position()) : -666);
        ftree->mu_bestTrack_BS_dz.push_back((hasBestTrack) ? bestTrack->dz(beamspot.position()) : -666);
        ftree->mu_bestTrack_dxyError.push_back((hasBestTrack) ? bestTrack->dxyError() : -666);
        ftree->mu_bestTrack_dzError.push_back((hasBestTrack) ? bestTrack->dzError() : -666);
        ftree->mu_bestTrack_normalizedChi2.push_back((hasBestTrack) ? bestTrack->normalizedChi2() : -666);
        ftree->mu_bestTrack_numberOfValidHits.push_back((hasBestTrack) ? bestTrack->numberOfValidHits() : -666);
        ftree->mu_bestTrack_numberOfLostHits.push_back((hasBestTrack) ? bestTrack->numberOfLostHits() : -666);
        ftree->mu_bestTrack_pt.push_back((hasBestTrack) ? bestTrack->pt() : -666);
        ftree->mu_bestTrack_eta.push_back((hasBestTrack) ? bestTrack->eta() : -666);
        ftree->mu_bestTrack_phi.push_back((hasBestTrack) ? bestTrack->phi() : -666);
        ftree->mu_bestTrack_ptError.push_back((hasBestTrack) ? bestTrack->ptError() : -666);
        ftree->mu_bestTrack_etaError.push_back((hasBestTrack) ? bestTrack->etaError() : -666);
        ftree->mu_bestTrack_phiError.push_back((hasBestTrack) ? bestTrack->phiError() : -666);
        ftree->mu_bestTrack_vx.push_back((hasBestTrack) ? bestTrack->vx() : -666);
        ftree->mu_bestTrack_vy.push_back((hasBestTrack) ? bestTrack->vy() : -666);
        ftree->mu_bestTrack_vz.push_back((hasBestTrack) ? bestTrack->vz() : -666);
        ftree->mu_bestTrack_qoverp.push_back((hasBestTrack) ? bestTrack->qoverp() : -666);
        ftree->mu_bestTrack_qoverpError.push_back((hasBestTrack) ? bestTrack->qoverpError() : -666);
        ftree->mu_bestTrack_charge.push_back((hasBestTrack) ? bestTrack->charge() : -666);
        ftree->mu_bestTrack_trackerLayersWithMeasurement.push_back((hasBestTrack) ? bestTrack->hitPattern().trackerLayersWithMeasurement() : -666);
        ftree->mu_bestTrack_pixelLayersWithMeasurement.push_back((hasBestTrack) ? bestTrack->hitPattern().pixelLayersWithMeasurement() : -666);
        ftree->mu_bestTrack_numberOfValidStripLayersWithMonoAndStereo.push_back((hasBestTrack) ? bestTrack->hitPattern().numberOfValidStripLayersWithMonoAndStereo() : -666);
        ftree->mu_bestTrack_trackerLayersWithoutMeasurement.push_back((hasBestTrack) ? bestTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_bestTrack_numberOfValidPixelHits.push_back((hasBestTrack) ? bestTrack->hitPattern().numberOfValidPixelHits() : -666);
        ftree->mu_bestTrack_numberOfLostPixelHits.push_back((hasBestTrack) ? bestTrack->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_bestTrack_numberOfInnerHits.push_back((hasBestTrack) ? bestTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) : -666);
        ftree->mu_bestTrack_numberOfOuterHits.push_back((hasBestTrack) ? bestTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS) : -666);
        ftree->mu_bestTrack_validFraction.push_back((hasBestTrack) ? bestTrack->validFraction() : -666);

        // InnerTrack
        const reco::TrackRef innerTrack = muon.innerTrack();
        bool hasInnerTrack = innerTrack.isNonnull();

        ftree->mu_hasInnerTrack.push_back(hasInnerTrack);
        ftree->mu_innerTrack_d0.push_back((hasInnerTrack) ? innerTrack->d0() : -666);
        ftree->mu_innerTrack_z0.push_back((hasInnerTrack) ? innerTrack->dz() : -666);
        ftree->mu_innerTrack_d0Error.push_back((hasInnerTrack) ? innerTrack->d0Error() : -666);
        ftree->mu_innerTrack_z0Error.push_back((hasInnerTrack) ? innerTrack->dzError() : -666);
        ftree->mu_innerTrack_PV_dxy.push_back((hasInnerTrack) ? innerTrack->dxy(primVtx->position()) : -666);
        ftree->mu_innerTrack_PV_dz.push_back((hasInnerTrack) ? innerTrack->dz(primVtx->position()) : -666);
        ftree->mu_innerTrack_RP_dxy.push_back((hasInnerTrack) ? innerTrack->dxy(bestTrack->referencePoint()) : -666);
        ftree->mu_innerTrack_RP_dz.push_back((hasInnerTrack) ? innerTrack->dz(bestTrack->referencePoint()) : -666);
        ftree->mu_innerTrack_BS_dxy.push_back((hasInnerTrack) ? innerTrack->dxy(beamspot.position()) : -666);
        ftree->mu_innerTrack_BS_dz.push_back((hasInnerTrack) ? innerTrack->dz(beamspot.position()) : -666);
        ftree->mu_innerTrack_dxyError.push_back((hasInnerTrack) ? innerTrack->dxyError() : -666);
        ftree->mu_innerTrack_dzError.push_back((hasInnerTrack) ? innerTrack->dzError() : -666);
        ftree->mu_innerTrack_normalizedChi2.push_back((hasInnerTrack) ? innerTrack->normalizedChi2() : -666);
        ftree->mu_innerTrack_numberOfValidHits.push_back((hasInnerTrack) ? innerTrack->numberOfValidHits() : -666);
        ftree->mu_innerTrack_numberOfLostHits.push_back((hasInnerTrack) ? innerTrack->numberOfLostHits() : -666);
        ftree->mu_innerTrack_pt.push_back((hasInnerTrack) ? innerTrack->pt() : -666);
        ftree->mu_innerTrack_eta.push_back((hasInnerTrack) ? innerTrack->eta() : -666);
        ftree->mu_innerTrack_phi.push_back((hasInnerTrack) ? innerTrack->phi() : -666);
        ftree->mu_innerTrack_ptError.push_back((hasInnerTrack) ? innerTrack->ptError() : -666);
        ftree->mu_innerTrack_etaError.push_back((hasInnerTrack) ? innerTrack->etaError() : -666);
        ftree->mu_innerTrack_phiError.push_back((hasInnerTrack) ? innerTrack->phiError() : -666);
        ftree->mu_innerTrack_vx.push_back((hasInnerTrack) ? innerTrack->vx() : -666);
        ftree->mu_innerTrack_vy.push_back((hasInnerTrack) ? innerTrack->vy() : -666);
        ftree->mu_innerTrack_vz.push_back((hasInnerTrack) ? innerTrack->vz() : -666);
        ftree->mu_innerTrack_qoverp.push_back((hasInnerTrack) ? innerTrack->qoverp() : -666);
        ftree->mu_innerTrack_qoverpError.push_back((hasInnerTrack) ? innerTrack->qoverpError() : -666);
        ftree->mu_innerTrack_charge.push_back((hasInnerTrack) ? innerTrack->charge() : -666);
        ftree->mu_innerTrack_trackerLayersWithMeasurement.push_back((hasInnerTrack) ? innerTrack->hitPattern().trackerLayersWithMeasurement() : -666);
        ftree->mu_innerTrack_pixelLayersWithMeasurement.push_back((hasInnerTrack) ? innerTrack->hitPattern().pixelLayersWithMeasurement() : -666);
        ftree->mu_innerTrack_numberOfValidStripLayersWithMonoAndStereo.push_back((hasInnerTrack) ? innerTrack->hitPattern().numberOfValidStripLayersWithMonoAndStereo() : -666);
        ftree->mu_innerTrack_trackerLayersWithoutMeasurement.push_back((hasInnerTrack) ? innerTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_innerTrack_numberOfValidPixelHits.push_back((hasInnerTrack) ? innerTrack->hitPattern().numberOfValidPixelHits() : -666);
        ftree->mu_innerTrack_numberOfLostPixelHits.push_back((hasInnerTrack) ? innerTrack->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS) : -666);
        ftree->mu_innerTrack_numberOfInnerHits.push_back((hasInnerTrack) ? innerTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) : -666);
        ftree->mu_innerTrack_numberOfOuterHits.push_back((hasInnerTrack) ? innerTrack->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS) : -666);
        ftree->mu_innerTrack_validFraction.push_back((hasInnerTrack) ? innerTrack->validFraction() : -666);

        // PF Isolation
        reco::MuonPFIsolation pfR03 = muon.pfIsolationR03();
        ftree->mu_pfIso03_sumChargedHadronPt.push_back(pfR03.sumChargedHadronPt);
        ftree->mu_pfIso03_sumChargedParticlePt.push_back(pfR03.sumChargedParticlePt);
        ftree->mu_pfIso03_sumNeutralHadronEt.push_back(pfR03.sumNeutralHadronEt);
        ftree->mu_pfIso03_sumNeutralHadronEtHighThreshold.push_back(pfR03.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfIso03_sumPhotonEt.push_back(pfR03.sumPhotonEt);
        ftree->mu_pfIso03_sumPhotonEtHighThreshold.push_back(pfR03.sumPhotonEtHighThreshold);
        ftree->mu_pfIso03_sumPUPt.push_back(pfR03.sumPUPt);

        reco::MuonPFIsolation pfR04 = muon.pfIsolationR04();
        ftree->mu_pfIso04_sumChargedHadronPt.push_back(pfR04.sumChargedHadronPt);
        ftree->mu_pfIso04_sumChargedParticlePt.push_back(pfR04.sumChargedParticlePt);
        ftree->mu_pfIso04_sumNeutralHadronEt.push_back(pfR04.sumNeutralHadronEt);
        ftree->mu_pfIso04_sumNeutralHadronEtHighThreshold.push_back(pfR04.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfIso04_sumPhotonEt.push_back(pfR04.sumPhotonEt);
        ftree->mu_pfIso04_sumPhotonEtHighThreshold.push_back(pfR04.sumPhotonEtHighThreshold);
        ftree->mu_pfIso04_sumPUPt.push_back(pfR04.sumPUPt);

        reco::MuonPFIsolation pfMeanR03 = muon.pfMeanDRIsoProfileR03();
        ftree->mu_pfMeanIso03_sumChargedHadronPt.push_back(pfMeanR03.sumChargedHadronPt);
        ftree->mu_pfMeanIso03_sumChargedParticlePt.push_back(pfMeanR03.sumChargedParticlePt);
        ftree->mu_pfMeanIso03_sumNeutralHadronEt.push_back(pfMeanR03.sumNeutralHadronEt);
        ftree->mu_pfMeanIso03_sumNeutralHadronEtHighThreshold.push_back(pfMeanR03.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfMeanIso03_sumPhotonEt.push_back(pfMeanR03.sumPhotonEt);
        ftree->mu_pfMeanIso03_sumPhotonEtHighThreshold.push_back(pfMeanR03.sumPhotonEtHighThreshold);
        ftree->mu_pfMeanIso03_sumPUPt.push_back(pfMeanR03.sumPUPt);

        reco::MuonPFIsolation pfSumR03 = muon.pfSumDRIsoProfileR03();
        ftree->mu_pfSumIso03_sumChargedHadronPt.push_back(pfSumR03.sumChargedHadronPt);
        ftree->mu_pfSumIso03_sumChargedParticlePt.push_back(pfSumR03.sumChargedParticlePt);
        ftree->mu_pfSumIso03_sumNeutralHadronEt.push_back(pfSumR03.sumNeutralHadronEt);
        ftree->mu_pfSumIso03_sumNeutralHadronEtHighThreshold.push_back(pfSumR03.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfSumIso03_sumPhotonEt.push_back(pfSumR03.sumPhotonEt);
        ftree->mu_pfSumIso03_sumPhotonEtHighThreshold.push_back(pfSumR03.sumPhotonEtHighThreshold);
        ftree->mu_pfSumIso03_sumPUPt.push_back(pfSumR03.sumPUPt);

        reco::MuonPFIsolation pfMeanR04 = muon.pfMeanDRIsoProfileR04();
        ftree->mu_pfMeanIso04_sumChargedHadronPt.push_back(pfMeanR04.sumChargedHadronPt);
        ftree->mu_pfMeanIso04_sumChargedParticlePt.push_back(pfMeanR04.sumChargedParticlePt);
        ftree->mu_pfMeanIso04_sumNeutralHadronEt.push_back(pfMeanR04.sumNeutralHadronEt);
        ftree->mu_pfMeanIso04_sumNeutralHadronEtHighThreshold.push_back(pfMeanR04.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfMeanIso04_sumPhotonEt.push_back(pfMeanR04.sumPhotonEt);
        ftree->mu_pfMeanIso04_sumPhotonEtHighThreshold.push_back(pfMeanR04.sumPhotonEtHighThreshold);
        ftree->mu_pfMeanIso04_sumPUPt.push_back(pfMeanR04.sumPUPt);

        reco::MuonPFIsolation pfSumR04 = muon.pfSumDRIsoProfileR04();
        ftree->mu_pfSumIso04_sumChargedHadronPt.push_back(pfSumR04.sumChargedHadronPt);
        ftree->mu_pfSumIso04_sumChargedParticlePt.push_back(pfSumR04.sumChargedParticlePt);
        ftree->mu_pfSumIso04_sumNeutralHadronEt.push_back(pfSumR04.sumNeutralHadronEt);
        ftree->mu_pfSumIso04_sumNeutralHadronEtHighThreshold.push_back(pfSumR04.sumNeutralHadronEtHighThreshold);
        ftree->mu_pfSumIso04_sumPhotonEt.push_back(pfSumR04.sumPhotonEt);
        ftree->mu_pfSumIso04_sumPhotonEtHighThreshold.push_back(pfSumR04.sumPhotonEtHighThreshold);
        ftree->mu_pfSumIso04_sumPUPt.push_back(pfSumR04.sumPUPt);

        ftree->mu_neutralHadronIso.push_back(muon.neutralHadronIso());
        ftree->mu_chargedHadronIso.push_back(muon.chargedHadronIso());
        ftree->mu_puChargedHadronIso.push_back(muon.puChargedHadronIso());
        ftree->mu_ecalIso.push_back(muon.ecalIso());
        ftree->mu_hcalIso.push_back(muon.hcalIso());
        ftree->mu_photonIso.push_back(muon.photonIso());
        ftree->mu_trackIso.push_back(muon.trackIso());

        // ID
        ftree->mu_isGlobalMuon.push_back(muon.isGlobalMuon());
        ftree->mu_isTrackerMuon.push_back(muon.isTrackerMuon());
        ftree->mu_isStandAloneMuon.push_back(muon.isStandAloneMuon());
        ftree->mu_isCaloMuon.push_back(muon.isCaloMuon());
        ftree->mu_isPFMuon.push_back(muon.isPFMuon());
        ftree->mu_isRPCMuon.push_back(muon.isRPCMuon());

        ftree->mu_isLooseMuon.push_back(muon.isLooseMuon());
        ftree->mu_isMediumMuon.push_back(muon.isMediumMuon());

        bool isTightMuon = 0;
        if( primVtx ) isTightMuon = muon.isTightMuon(*primVtx);
        ftree->mu_isTightMuon.push_back(isTightMuon);
        bool isSoftMuon = 0;
        if( primVtx ) isSoftMuon = muon.isSoftMuon(*primVtx);
        ftree->mu_isSoftMuon.push_back(isSoftMuon);
        bool isHighPtMuon = 0;
        if( primVtx ) isHighPtMuon = muon.isHighPtMuon(*primVtx);
        ftree->mu_isHighPtMuon.push_back(isHighPtMuon);

        ftree->mu_type.push_back(muon.type());

        ftree->mu_caloCompatibility.push_back(muon.caloCompatibility());
        ftree->mu_segmentCompatibility.push_back(muon.segmentCompatibility());

        // vertex
        ftree->mu_vx.push_back(muon.vx());
        ftree->mu_vy.push_back(muon.vy());
        ftree->mu_vz.push_back(muon.vz());

        // mini-iso
        float miniIso           = -666;
        float miniIsoTTH        = -666;
        float miniIsoTTHCharged = -666;
        float miniIsoTTHNeutral = -666;
        if( dataFormat_ != "AOD" )
	  {
	     // Variables used for stop analysis
	     double EA_miniIso = getEA(muon.eta(), MU_EA_ETA, MU_EA_VALUE);
	     miniIso = getPFIsolation(pfcands,dynamic_cast<const reco::Candidate*>(&muon),*rhoPtr, EA_miniIso, 0.05,0.2,10.,false, false);
	     // -------------------------------------------
	     //
	     // Below: Variables used for ttH analysis
	     float miniIsoR = 10.0/std::min(std::max(float(muon.pt()),float(50.)),float(200.));
	     float EffArea = 0.;
	     float eta = muon.eta();
	     
	     if(      fabs(eta) > 0    && fabs(eta) < 0.8 ) EffArea = 0.0566;
	     else if( fabs(eta) >= 0.8 && fabs(eta) < 1.3 ) EffArea = 0.0562;
	     else if( fabs(eta) >= 1.3 && fabs(eta) < 2.0 ) EffArea = 0.0363;
	     else if( fabs(eta) >= 2.0 && fabs(eta) < 2.2 ) EffArea = 0.0119;
	     else if( fabs(eta) >= 2.2 && fabs(eta) < 2.5 ) EffArea = 0.0064;
	     
	     float correction = ftree->ev_rho*EffArea*(miniIsoR/0.3)*(miniIsoR/0.3);
	     
	     //miniIso = getPFIsolation(pfcands,dynamic_cast<const reco::Candidate*>(&muon),0.05,0.2,10.,false,false);
	     
	     float pfIsoCharged      = MuonPfIsoCharged(muon,pfcands,miniIsoR);
	     float pfIsoNeutral      = MuonPfIsoNeutral(muon,pfcands,miniIsoR);
	     float pfIsoPUSubtracted = std::max(float(0.0),float(pfIsoNeutral-correction));
	     
	     miniIsoTTH        = (pfIsoCharged + pfIsoPUSubtracted) / muon.pt();
	     miniIsoTTHCharged = pfIsoCharged / muon.pt();
	     miniIsoTTHNeutral = pfIsoPUSubtracted / muon.pt();
	  }
	
        ftree->mu_miniIso.push_back(miniIso);
        ftree->mu_miniIsoTTH.push_back(miniIsoTTH);
	
        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMuonAnalysis#Muon_identification
        ftree->mu_isGoodMuon_AllGlobalMuons.push_back(muon::isGoodMuon(muon,muon::AllGlobalMuons));
        ftree->mu_isGoodMuon_AllStandAloneMuons.push_back(muon::isGoodMuon(muon,muon::AllStandAloneMuons));
        ftree->mu_isGoodMuon_AllTrackerMuons.push_back(muon::isGoodMuon(muon,muon::AllTrackerMuons));
        ftree->mu_isGoodMuon_TrackerMuonArbitrated.push_back(muon::isGoodMuon(muon,muon::TrackerMuonArbitrated));
        ftree->mu_isGoodMuon_AllArbitrated.push_back(muon::isGoodMuon(muon,muon::AllArbitrated));
        ftree->mu_isGoodMuon_GlobalMuonPromptTight.push_back(muon::isGoodMuon(muon,muon::GlobalMuonPromptTight));
        ftree->mu_isGoodMuon_TMLastStationLoose.push_back(muon::isGoodMuon(muon,muon::TMLastStationLoose));
        ftree->mu_isGoodMuon_TMLastStationTight.push_back(muon::isGoodMuon(muon,muon::TMLastStationTight));
        ftree->mu_isGoodMuon_TM2DCompatibilityLoose.push_back(muon::isGoodMuon(muon,muon::TM2DCompatibilityLoose));
        ftree->mu_isGoodMuon_TM2DCompatibilityTight.push_back(muon::isGoodMuon(muon,muon::TM2DCompatibilityTight));
        ftree->mu_isGoodMuon_TMOneStationLoose.push_back(muon::isGoodMuon(muon,muon::TMOneStationLoose));
        ftree->mu_isGoodMuon_TMOneStationTight.push_back(muon::isGoodMuon(muon,muon::TMOneStationTight));
        ftree->mu_isGoodMuon_TMLastStationOptimizedLowPtLoose.push_back(muon::isGoodMuon(muon,muon::TMLastStationOptimizedLowPtLoose));
        ftree->mu_isGoodMuon_TMLastStationOptimizedLowPtTight.push_back(muon::isGoodMuon(muon,muon::TMLastStationOptimizedLowPtTight));
        ftree->mu_isGoodMuon_GMTkChiCompatibility.push_back(muon::isGoodMuon(muon,muon::GMTkChiCompatibility));
        ftree->mu_isGoodMuon_GMStaChiCompatibility.push_back(muon::isGoodMuon(muon,muon::GMStaChiCompatibility));
        ftree->mu_isGoodMuon_GMTkKinkTight.push_back(muon::isGoodMuon(muon,muon::GMTkKinkTight));
        ftree->mu_isGoodMuon_TMLastStationAngLoose.push_back(muon::isGoodMuon(muon,muon::TMLastStationAngLoose));
        ftree->mu_isGoodMuon_TMLastStationAngTight.push_back(muon::isGoodMuon(muon,muon::TMLastStationAngTight));
        ftree->mu_isGoodMuon_TMOneStationAngLoose.push_back(muon::isGoodMuon(muon,muon::TMOneStationAngLoose));
        ftree->mu_isGoodMuon_TMOneStationAngTight.push_back(muon::isGoodMuon(muon,muon::TMOneStationAngTight));
        ftree->mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtLoose.push_back(muon::isGoodMuon(muon,muon::TMLastStationOptimizedBarrelLowPtLoose));
        ftree->mu_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight.push_back(muon::isGoodMuon(muon,muon::TMLastStationOptimizedBarrelLowPtTight));

        bool energyIsValid = muon.isEnergyValid();
        ftree->mu_calEnergy_em.push_back(energyIsValid ? muon.calEnergy().em : -666.);
        ftree->mu_calEnergy_had.push_back(energyIsValid ? muon.calEnergy().had : -666.);
        ftree->mu_calEnergy_ho.push_back(energyIsValid ? muon.calEnergy().ho : -666.);
        ftree->mu_calEnergy_emS9.push_back(energyIsValid ? muon.calEnergy().emS9 : -666.);
        ftree->mu_calEnergy_hadS9.push_back(energyIsValid ? muon.calEnergy().hadS9 : -666.);
        ftree->mu_calEnergy_hoS9.push_back(energyIsValid ? muon.calEnergy().hoS9 : -666.);
        ftree->mu_calEnergy_emS25.push_back(energyIsValid ? muon.calEnergy().emS25 : -666.);
        ftree->mu_calEnergy_emMax.push_back(energyIsValid ? muon.calEnergy().emMax : -666.);
        ftree->mu_calEnergy_hadMax.push_back(energyIsValid ? muon.calEnergy().hadMax : -666.);
        ftree->mu_calEnergy_ecal_time.push_back(energyIsValid ? muon.calEnergy().ecal_time : -666.);
        ftree->mu_calEnergy_hcal_time.push_back(energyIsValid ? muon.calEnergy().hcal_time : -666.);
        ftree->mu_calEnergy_ecal_rawId.push_back(energyIsValid ? muon.calEnergy().ecal_id.rawId() : -666.);
        ftree->mu_calEnergy_hcal_rawId.push_back(energyIsValid ? muon.calEnergy().hcal_id.rawId() : -666.);

        bool isoIsValid = muon.isIsolationValid();

        ftree->mu_isolationR03_trackerVetoPt.push_back(isoIsValid ? muon.isolationR03().trackerVetoPt : -666.);
        ftree->mu_isolationR03_emVetoEt.push_back(isoIsValid ? muon.isolationR03().emVetoEt : -666.);
        ftree->mu_isolationR03_hadVetoEt.push_back(isoIsValid ? muon.isolationR03().hadVetoEt : -666.);
        ftree->mu_isolationR03_hoVetoEt.push_back(isoIsValid ? muon.isolationR03().hoVetoEt : -666.);
        ftree->mu_isolationR03_sumPt.push_back(isoIsValid ? muon.isolationR03().sumPt : -666.);
        ftree->mu_isolationR03_emEt.push_back(isoIsValid ? muon.isolationR03().emEt : -666.);
        ftree->mu_isolationR03_hadEt.push_back(isoIsValid ? muon.isolationR03().hadEt : -666.);
        ftree->mu_isolationR03_hoEt.push_back(isoIsValid ? muon.isolationR03().hoEt : -666.);
        ftree->mu_isolationR03_nTracks.push_back(isoIsValid ? muon.isolationR03().nTracks : -666.);
        ftree->mu_isolationR03_nJets.push_back(isoIsValid ? muon.isolationR03().nJets : -666.);

        ftree->mu_isolationR05_trackerVetoPt.push_back(isoIsValid ? muon.isolationR05().trackerVetoPt : -666.);
        ftree->mu_isolationR05_emVetoEt.push_back(isoIsValid ? muon.isolationR05().emVetoEt : -666.);
        ftree->mu_isolationR05_hadVetoEt.push_back(isoIsValid ? muon.isolationR05().hadVetoEt : -666.);
        ftree->mu_isolationR05_hoVetoEt.push_back(isoIsValid ? muon.isolationR05().hoVetoEt : -666.);
        ftree->mu_isolationR05_sumPt.push_back(isoIsValid ? muon.isolationR05().sumPt : -666.);
        ftree->mu_isolationR05_emEt.push_back(isoIsValid ? muon.isolationR05().emEt : -666.);
        ftree->mu_isolationR05_hadEt.push_back(isoIsValid ? muon.isolationR05().hadEt : -666.);
        ftree->mu_isolationR05_hoEt.push_back(isoIsValid ? muon.isolationR05().hoEt : -666.);
        ftree->mu_isolationR05_nTracks.push_back(isoIsValid ? muon.isolationR05().nTracks : -666.);
        ftree->mu_isolationR05_nJets.push_back(isoIsValid ? muon.isolationR05().nJets : -666.);

        // ttH lepton MVA
        double mu_pt = muon.pt();
        double mu_eta = muon.eta();
        double mu_phi = muon.phi();
        double mu_lepMVA = -666.;
	
        float drmin = 0.5;
        int jcl = -1;
        for(unsigned int ij=0;ij<jets->size();ij++)
	  {
	     if( jets->at(ij).pt() < 10. ) continue;
	     float dr = GetDeltaR(jets->at(ij).eta(),
				  jets->at(ij).phi(),
				  mu_eta,
				  mu_phi);
	     if( dr < drmin )
	       {
		  drmin = dr;
		  jcl = ij;
	       }
	  }
	
        lepMVA_pt                                   = mu_pt;
        lepMVA_eta                                  = mu_eta;
        lepMVA_miniRelIsoNeutral                    = miniIsoTTHNeutral; //CAREFUL! WAS CHANGED TO MATCH GEOFF DEFINITION...
        lepMVA_miniRelIsoCharged                    = miniIsoTTHCharged;
        lepMVA_jetPtRatio                           = (jcl >= 0) ? ptRatioMuon(muon,jets->at(jcl)) : 1.5;
        lepMVA_jetPtRelv2                           = (jcl >= 0) ? ptRelMuon(muon,jets->at(jcl)) : 0.0;
        float csv                                   = (jcl >= 0) ? jets->at(jcl).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -666;
        lepMVA_jetBTagCSV                           = std::max(double(csv),0.);
        lepMVA_sip3d                                = fabs(ftree->mu_ip3d.back()/ftree->mu_ip3dErr.back());
        lepMVA_dxy                                  = log(fabs(ftree->mu_innerTrack_PV_dxy.back()));
        lepMVA_dz                                   = log(fabs(ftree->mu_innerTrack_PV_dz.back()));
        lepMVA_mvaId                                = ftree->mu_segmentCompatibility.back();
        lepMVA_jetNDauChargedMVASel                 = (jcl >= 0) ? jetNDauChargedMVASel(jets->at(jcl), *primVtx) : 0.0; //?? correct default value
	
        mu_lepMVA = mu_reader->EvaluateMVA("BDTG method");
	
        ftree->mu_lepMVA.push_back(mu_lepMVA);
        ftree->mu_lepMVA_pt.push_back(lepMVA_pt); 
        ftree->mu_lepMVA_eta.push_back(lepMVA_eta);
        ftree->mu_lepMVA_miniRelIsoCharged.push_back(lepMVA_miniRelIsoCharged);
        ftree->mu_lepMVA_miniRelIsoNeutral.push_back(lepMVA_miniRelIsoNeutral);
        ftree->mu_lepMVA_jetPtRatio.push_back(lepMVA_jetPtRatio);
        ftree->mu_lepMVA_jetPtRelv2.push_back(lepMVA_jetPtRelv2);
        ftree->mu_lepMVA_jetBTagCSV.push_back(lepMVA_jetBTagCSV);
        ftree->mu_lepMVA_sip3d.push_back(lepMVA_sip3d);
        ftree->mu_lepMVA_dxy.push_back(lepMVA_dxy);
        ftree->mu_lepMVA_dz.push_back(lepMVA_dz);
        ftree->mu_lepMVA_mvaId.push_back(lepMVA_mvaId);
        ftree->mu_lepMVA_jetNDauChargedMVASel.push_back(lepMVA_jetNDauChargedMVASel);
	
        float conept = (jcl >= 0) ? conePtMuon(muon,jets->at(jcl)) : -100;
        ftree->mu_conept.push_back( conept );
	
        bool debug = false;

        if(debug)
	  {
	     std::cout << "mu["              << im                      << "]: "
	       << " pt= "            << lepMVA_pt
	       << " eta= "           << lepMVA_eta
	       << " miniIsoN= "      << lepMVA_miniRelIsoNeutral
	       << " miniIsoC= "      << lepMVA_miniRelIsoCharged
	       << " jetPtRatio= "    << lepMVA_jetPtRatio
	       << " jetPtRel= "      << lepMVA_jetPtRelv2 
	       << " mulepMVA= "      << mu_lepMVA                      << std::endl;
	  }
	
        if(false)
	  {
	     std::cout << ftree->ev_id                   << " "
	       << lepMVA_pt                      << " "
	       << lepMVA_eta                     << " "
	       << muon.phi()                     << " "
	       << muon.energy()                  << " "
	       << muon.pdgId()                   << " "
	       << muon.charge()                  << " "
	       << lepMVA_jetNDauChargedMVASel    << " "
	       << miniIsoTTH                     << " "
	       << lepMVA_miniRelIsoCharged       << " "
	       << lepMVA_miniRelIsoNeutral       << " "
	       << lepMVA_jetPtRelv2              << " "
	       << lepMVA_jetBTagCSV              << " "
	       << lepMVA_jetPtRatio              << " "
	       << fabs(ftree->mu_ip3d.back()/ftree->mu_ip3dErr.back())         << " "
	       << fabs(ftree->mu_innerTrack_PV_dxy.back())                     << " "
	       << fabs(ftree->mu_innerTrack_PV_dz.back())                      << " "
	       << lepMVA_mvaId                   << " "
	       << mu_lepMVA            << " "
	       << std::endl;
	  }
	
        if( !isData_ )
	  {
	     // Internal matching
	     reco::GenParticle *genp = new reco::GenParticle();
	     
	     float drmin;
	     bool hasMCMatch = mc_truth->doMatch(iEvent,iSetup,genParticlesHandle,genp,drmin,
						 muon.pt(),muon.eta(),muon.phi(),muon.pdgId());
	     ftree->mu_hasMCMatch.push_back(hasMCMatch);
	     if( hasMCMatch )
	       {
		  ftree->mu_gen_pt.push_back(genp->pt());
		  ftree->mu_gen_eta.push_back(genp->eta());
		  ftree->mu_gen_phi.push_back(genp->phi());
		  ftree->mu_gen_m.push_back(genp->mass());
		  ftree->mu_gen_status.push_back(genp->status());
		  ftree->mu_gen_id.push_back(genp->pdgId());
		  ftree->mu_gen_charge.push_back(genp->charge());
		  ftree->mu_gen_dr.push_back(drmin);
	       }
	     else
	       {
		  ftree->mu_gen_pt.push_back(-666);
		  ftree->mu_gen_eta.push_back(-666);
		  ftree->mu_gen_phi.push_back(-666);
		  ftree->mu_gen_m.push_back(-666);
		  ftree->mu_gen_status.push_back(-666);
		  ftree->mu_gen_id.push_back(-666);
		  ftree->mu_gen_charge.push_back(-666);
		  ftree->mu_gen_dr.push_back(-666);
	       }
	     delete genp;
	     
	     // PAT matching
	     const reco::GenParticle *genpPAT = muon.genParticle();
	     bool hasMCMatchPAT = (genpPAT != 0);
	     ftree->mu_hasMCMatchPAT.push_back(hasMCMatchPAT);
	     if( hasMCMatchPAT )
	       {
		  ftree->mu_genPAT_pt.push_back(genpPAT->pt());
		  ftree->mu_genPAT_eta.push_back(genpPAT->eta());
		  ftree->mu_genPAT_phi.push_back(genpPAT->phi());
		  ftree->mu_genPAT_m.push_back(genpPAT->mass());
		  ftree->mu_genPAT_status.push_back(genpPAT->status());
		  ftree->mu_genPAT_id.push_back(genpPAT->pdgId());
		  ftree->mu_genPAT_charge.push_back(genpPAT->charge());
	       }
	     else
	       {
		  ftree->mu_genPAT_pt.push_back(-666);
		  ftree->mu_genPAT_eta.push_back(-666);
		  ftree->mu_genPAT_phi.push_back(-666);
		  ftree->mu_genPAT_m.push_back(-666);
		  ftree->mu_genPAT_status.push_back(-666);
		  ftree->mu_genPAT_id.push_back(-666);
		  ftree->mu_genPAT_charge.push_back(-666);
	       }
	  }
     }
   ftree->mu_n = ftree->mu_pt.size();
   
   // ########################
   // #  _                   #
   // # | |_ __ _ _   _ ___  #
   // # | __/ _` | | | / __| #
   // # | || (_| | |_| \__ \ #
   // #  \__\__,_|\__,_|___/ #
   // #                      #
   // ########################
   // 
   int nTau = taus->size();
   for(int it=0;it<nTau;it++)
     {
        const pat::Tau& tau = taus->at(it);

        // Skimming taus with pT < 5 GeV. (should do nothing for miniAOD where pT > 18 GeV is applied)
        //if (tau.pt() < 5) continue;
	
        ftree->tau_pt.push_back(tau.pt());
        ftree->tau_eta.push_back(tau.eta());
        ftree->tau_phi.push_back(tau.phi());
        ftree->tau_m.push_back(tau.mass());
        ftree->tau_E.push_back(tau.energy());
        ftree->tau_id.push_back(tau.pdgId());
        ftree->tau_charge.push_back(tau.charge());

        // https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Taus
        //
        float tau_leadingTrackPt = -666;
        int tau_leadingTrackCharge = -666;
        float tau_leadingTrackDz = -666;
        float tau_leadingTrackDxy = -666;

        ftree->tau_hasLeadChargedHadrCand.push_back(tau.leadChargedHadrCand().isNonnull());

	if( tau.leadChargedHadrCand().isNonnull() )
	  {
	     pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
	     tau_leadingTrackPt = packedLeadTauCand->pt();
	     tau_leadingTrackCharge = packedLeadTauCand->charge();
	     tau_leadingTrackDz = packedLeadTauCand->dz();
	     tau_leadingTrackDxy = packedLeadTauCand->dxy();
	  }
	
	ftree->tau_leadingTrackPt.push_back(tau_leadingTrackPt);       
	ftree->tau_leadingTrackCharge.push_back(tau_leadingTrackCharge);
	ftree->tau_leadingTrackDz.push_back(tau_leadingTrackDz);
	ftree->tau_leadingTrackDxy.push_back(tau_leadingTrackDxy);
	
	ftree->tau_decayMode.push_back(tau.decayMode());       
	ftree->tau_decayModeFinding.push_back(tau.tauID("decayModeFinding"));
	//       ftree->tau_decayModeFindingOldDMs.push_back(tau.tauID("decayModeFindingOldDMs"));
	ftree->tau_decayModeFindingNewDMs.push_back(tau.tauID("decayModeFindingNewDMs"));
	
	ftree->tau_puCorrPtSum.push_back(tau.tauID("puCorrPtSum"));
	ftree->tau_neutralIsoPtSum.push_back(tau.tauID("neutralIsoPtSum"));
	ftree->tau_chargedIsoPtSum.push_back(tau.tauID("chargedIsoPtSum"));
	ftree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
	
	ftree->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
	ftree->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
	ftree->tau_byTightCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
	
	ftree->tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau.tauID("byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
	ftree->tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau.tauID("byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
	ftree->tau_byTightIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau.tauID("byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
	ftree->tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau.tauID("byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
	ftree->tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT.push_back(tau.tauID("byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
	
	ftree->tau_againstMuonLoose3.push_back(tau.tauID("againstMuonLoose3"));
	ftree->tau_againstMuonTight3.push_back(tau.tauID("againstMuonTight3"));
	
	ftree->tau_againstElectronVLooseMVA6.push_back(tau.tauID("againstElectronVLooseMVA6"));
	ftree->tau_againstElectronLooseMVA6.push_back(tau.tauID("againstElectronLooseMVA6"));
	ftree->tau_againstElectronMediumMVA6.push_back(tau.tauID("againstElectronMediumMVA6"));
	ftree->tau_againstElectronTightMVA6.push_back(tau.tauID("againstElectronTightMVA6"));
	
	ftree->tau_pfEssential_jet_pt.push_back(tau.pfEssential().p4Jet_.pt());
	ftree->tau_pfEssential_jet_eta.push_back(tau.pfEssential().p4Jet_.eta());
	ftree->tau_pfEssential_jet_phi.push_back(tau.pfEssential().p4Jet_.phi());
	ftree->tau_pfEssential_jet_m.push_back(tau.pfEssential().p4Jet_.mass());
	
	ftree->tau_pfEssential_jetCorr_pt.push_back(tau.pfEssential().p4CorrJet_.pt());
	ftree->tau_pfEssential_jetCorr_eta.push_back(tau.pfEssential().p4CorrJet_.eta());
	ftree->tau_pfEssential_jetCorr_phi.push_back(tau.pfEssential().p4CorrJet_.phi());
	ftree->tau_pfEssential_jetCorr_m.push_back(tau.pfEssential().p4CorrJet_.mass());
	
        float tau_pfEssential_sv_x = -666;
        float tau_pfEssential_sv_y = -666;
        float tau_pfEssential_sv_z = -666;
	
        ftree->tau_pfEssential_hasSV.push_back(tau.pfEssential().sv_.isNonnull());

        if( tau.pfEssential().sv_.isNonnull() )
	  {
	     tau_pfEssential_sv_x = tau.pfEssential().sv_->x();
            tau_pfEssential_sv_y = tau.pfEssential().sv_->y();
	     tau_pfEssential_sv_z = tau.pfEssential().sv_->z();
	  }
	
        ftree->tau_pfEssential_sv_x.push_back(tau_pfEssential_sv_x);
        ftree->tau_pfEssential_sv_y.push_back(tau_pfEssential_sv_y);
        ftree->tau_pfEssential_sv_z.push_back(tau_pfEssential_sv_z);
	
        ftree->tau_pfEssential_flightLengthSig.push_back(tau.pfEssential().flightLengthSig_);
        ftree->tau_pfEssential_dxy.push_back(tau.pfEssential().dxy_);
        ftree->tau_pfEssential_dxy_error.push_back(tau.pfEssential().dxy_error_);
        ftree->tau_pfEssential_dxy_Sig.push_back(tau.pfEssential().dxy_Sig_);
       
        /*	ftree->tau_pfEssential_flightLengthSig.push_back(tau.pfEssential().flightLengthSig);
	 ftree->tau_pfEssential_dxy.push_back(tau.pfEssential().dxy);
	 ftree->tau_pfEssential_dxy_error.push_back(tau.pfEssential().dxy_error);
	 ftree->tau_pfEssential_dxy_Sig.push_back(tau.pfEssential().dxy_Sig);*/
     }

   ftree->tau_n = ftree->tau_pt.size();
   
   // ##########################
   // #       _      _         #
   // #      | | ___| |_ ___   #
   // #   _  | |/ _ \ __/ __|  #
   // #  | |_| |  __/ |_\__ \  #
   // #   \___/ \___|\__|___/  #
   // #                        #
   // ##########################
   //
   // Jets
   //
   int nJet = jets->size();
   ftree->jet_n = nJet;
   for(int ij=0;ij<nJet;ij++)
     {
        const pat::Jet& jet = jets->at(ij);
	
        ftree->jet_pt.push_back(jet.pt());
        ftree->jet_eta.push_back(jet.eta());
        ftree->jet_phi.push_back(jet.phi());
        ftree->jet_m.push_back(jet.mass());
        ftree->jet_E.push_back(jet.energy());
        ftree->jet_JBP.push_back(jet.bDiscriminator("pfJetBProbabilityBJetTags"));
        ftree->jet_JP.push_back(jet.bDiscriminator("pfJetProbabilityBJetTags"));
        ftree->jet_TCHP.push_back(jet.bDiscriminator("pfTrackCountingHighPurBJetTags"));
        ftree->jet_TCHE.push_back(jet.bDiscriminator("pfTrackCountingHighEffBJetTags"));
        ftree->jet_SSVHE.push_back(jet.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"));
        ftree->jet_SSVHP.push_back(jet.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
        ftree->jet_CMVA.push_back(jet.bDiscriminator("pfCombinedMVABJetTags"));

        ftree->jet_chargedMultiplicity.push_back(jet.chargedMultiplicity());
        ftree->jet_neutralMultiplicity.push_back(jet.neutralMultiplicity());
        ftree->jet_chargedHadronMultiplicity.push_back(jet.chargedHadronMultiplicity());

        ftree->jet_jecFactorUncorrected.push_back(jet.jecFactor("Uncorrected"));
        ftree->jet_jecFactorL1FastJet.push_back(jet.jecFactor("L1FastJet"));
        ftree->jet_jecFactorL2Relative.push_back(jet.jecFactor("L2Relative"));
        ftree->jet_jecFactorL3Absolute.push_back(jet.jecFactor("L3Absolute"));
	
	ftree->jet_jetArea.push_back(jet.jetArea());
	
        jecUnc->setJetEta(fabs(jet.eta()));
        jecUnc->setJetPt(jet.pt());
	
        ftree->jet_Unc.push_back(jecUnc->getUncertainty(true));
	
        ftree->jet_ntrk.push_back(jet.associatedTracks().size());
        //	std::cout << jet.hasTagInfo("pfInclusiveSecondaryVertexFinderTagInfos") << std::endl;

        float CSVIVF = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        ftree->jet_CSVv2.push_back(CSVIVF);
	
        ftree->jet_DeepCSVProbudsg.push_back(jet.bDiscriminator("pfDeepCSVJetTags:probudsg"));
        ftree->jet_DeepCSVProbb.push_back(jet.bDiscriminator("pfDeepCSVJetTags:probb"));
        ftree->jet_DeepCSVProbc.push_back(jet.bDiscriminator("pfDeepCSVJetTags:probc"));
        ftree->jet_DeepCSVProbbb.push_back(jet.bDiscriminator("pfDeepCSVJetTags:probbb"));
        ftree->jet_DeepCSVProbcc.push_back(jet.bDiscriminator("pfDeepCSVJetTags:probcc"));
		
		ftree->jet_DeepFlavourProbuds.push_back(jet.bDiscriminator("pfDeepFlavourJetTags:probuds"));
		ftree->jet_DeepFlavourProbg.push_back(jet.bDiscriminator("pfDeepFlavourJetTags:probg"));
        ftree->jet_DeepFlavourProbb.push_back(jet.bDiscriminator("pfDeepFlavourJetTags:probb"));
        ftree->jet_DeepFlavourProbbb.push_back(jet.bDiscriminator("pfDeepFlavourJetTags:probbb"));
        ftree->jet_DeepFlavourProblepb.push_back(jet.bDiscriminator("pfDeepFlavourJetTags:problepb"));	
		ftree->jet_DeepFlavourProbc.push_back(jet.bDiscriminator("pfDeepFlavourJetTags:probc"));
			
	
	/*
	 cout<<"jet.bDiscriminator(pfCombinedInclusiveSecondaryVertexV2BJetTags) = "<<jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<<endl;
	 cout<<"jet.bDiscriminator(deepFlavourJetTags:probudsg) = "<<jet.bDiscriminator("deepFlavourJetTags:probudsg")<<endl;
	 cout<<"jet.bDiscriminator(deepFlavourJetTags:probb) = "<<jet.bDiscriminator("deepFlavourJetTags:probb")<<endl;
	 cout<<"jet.bDiscriminator(deepFlavourJetTags:probc) = "<<jet.bDiscriminator("deepFlavourJetTags:probc")<<endl;
	 cout<<"jet.bDiscriminator(deepFlavourJetTags:probbb) = "<<jet.bDiscriminator("deepFlavourJetTags:probbb")<<endl;
	 cout<<"jet.bDiscriminator(deepFlavourJetTags:probcc) = "<<jet.bDiscriminator("deepFlavourJetTags:probcc")<<endl;
	 */
	
        ftree->jet_cMVAv2.push_back(jet.bDiscriminator("pfCombinedMVAV2BJetTags"));
        ftree->jet_CharmCvsL.push_back(jet.bDiscriminator("pfCombinedCvsLJetTags"));
        ftree->jet_CharmCvsB.push_back(jet.bDiscriminator("pfCombinedCvsBJetTags"));

        ftree->jet_partonFlavour.push_back(jet.partonFlavour());
        ftree->jet_hadronFlavour.push_back(jet.hadronFlavour());
	
        ftree->jet_neutralHadronEnergy.push_back(jet.neutralHadronEnergy());
        ftree->jet_neutralEmEnergy.push_back(jet.neutralEmEnergy());
        ftree->jet_chargedHadronEnergy.push_back(jet.chargedHadronEnergy());
        ftree->jet_chargedEmEnergy.push_back(jet.chargedEmEnergy());
        ftree->jet_electronEnergy.push_back(jet.electronEnergy());
        ftree->jet_muonEnergy.push_back(jet.muonEnergy());
        ftree->jet_photonEnergy.push_back(jet.photonEnergy());
	
        if( jet.hasUserFloat("pileupJetId:fullDiscriminant") )
	  ftree->jet_pileupJetId.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
        else
	  ftree->jet_pileupJetId.push_back(-666.);
	
        // Jet ID
	// 
	// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
	// 
        float NHF = jet.neutralHadronEnergyFraction();
        float NEMF = jet.neutralEmEnergyFraction();
        float CHF = jet.chargedHadronEnergyFraction();
        float MUF = jet.muonEnergyFraction();
        float CEMF = jet.chargedEmEnergyFraction();
        float NEM = jet.neutralMultiplicity();
        float CHM = jet.chargedMultiplicity();
        float NumConst = NEM + CHM;
        float eta = jet.eta();
	//        bool looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8) && ((fabs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(eta)>2.4);
	// 
        bool tightJetID = 1;
        bool tightLepVetoJetID = 1;
	
        if( fabs(eta) < 2.7 ) tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8);
        if( fabs(eta) <= 2.4 ) tightLepVetoJetID = (CHF>0 && CHM>0 && CEMF<0.80);
	
        if( fabs(eta) < 2.7 ) tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1);
        if( fabs(eta) <= 2.4 ) tightJetID = (CHF>0 && CHM>0);
        if( fabs(eta) >= 2.7 && fabs(eta) < 3.0 ) tightJetID = (NEMF>0.02 && NEMF<0.99 && NEM>2);
        if( fabs(eta) >= 3.0 ) tightJetID = (NEMF<0.90 && NHF>0.02 && NEM>10);
	
        ftree->jet_neutralHadronEnergyFraction.push_back(jet.neutralHadronEnergyFraction());
        ftree->jet_neutralEmEnergyFraction.push_back(jet.neutralEmEnergyFraction());
        ftree->jet_chargedHadronEnergyFraction.push_back(jet.chargedHadronEnergyFraction());
        ftree->jet_muonEnergyFraction.push_back(jet.muonEnergyFraction());
        ftree->jet_chargedEmEnergyFraction.push_back(jet.chargedEmEnergyFraction());
	
	//        ftree->jet_looseJetID.push_back(looseJetID);
        ftree->jet_tightJetID.push_back(tightJetID);
        ftree->jet_tightLepVetoJetID.push_back(tightLepVetoJetID);
	
        //Quark-gluon tagging
        const auto jetRef = view_jets->ptrAt(ij);
        if( ! qgHandle.failedToGet() )
	  ftree->jet_qgtag.push_back((*qgHandle)[jetRef]);
        else
	  ftree->jet_qgtag.push_back(-666.);
	
        const reco::GenJet* genJet = jet.genJet();
        bool hasGenInfo = (genJet);
        ftree->jet_hasGenJet.push_back(hasGenInfo);
	
        float gen_jet_pt = -666;
        float gen_jet_eta = -666;
        float gen_jet_phi = -666;
        float gen_jet_m = -666;
        float gen_jet_E = -666;
        int gen_jet_status = -666;
        int gen_jet_id = -666;
	
        if( hasGenInfo )
	  {
	     gen_jet_pt = genJet->pt();
	     gen_jet_eta = genJet->eta();
	     gen_jet_phi = genJet->phi();
	     gen_jet_m = genJet->mass();
	     gen_jet_E = genJet->energy();
	     gen_jet_status = genJet->status();
	     gen_jet_id = genJet->pdgId();
	  }
	
        ftree->jet_genJet_pt.push_back(gen_jet_pt);
        ftree->jet_genJet_eta.push_back(gen_jet_eta);
        ftree->jet_genJet_phi.push_back(gen_jet_phi);
        ftree->jet_genJet_m.push_back(gen_jet_m);
        ftree->jet_genJet_E.push_back(gen_jet_E);
	
        ftree->jet_genJet_status.push_back(gen_jet_status);
        ftree->jet_genJet_id.push_back(gen_jet_id);
	
        const reco::GenParticle* genParton = (!isData_) ? jet.genParton() : 0;
        bool hasGenPartonInfo = (genParton);
        ftree->jet_hasGenParton.push_back(hasGenPartonInfo);

        float gen_parton_pt = -666;
        float gen_parton_eta = -666;
        float gen_parton_phi = -666;
        float gen_parton_m = -666;
        float gen_parton_E = -666;
        int gen_parton_status = -666;
        int gen_parton_id = -666;
	
        if( hasGenPartonInfo )
	  {
	     gen_parton_pt = genParton->pt();
	     gen_parton_eta = genParton->eta();
	     gen_parton_phi = genParton->phi();
	     gen_parton_m = genParton->mass();
	     gen_parton_E = genParton->energy();
	     gen_parton_status = genParton->status();
	     gen_parton_id = genParton->pdgId();
	  }
	
        ftree->jet_genParton_pt.push_back(gen_parton_pt);
        ftree->jet_genParton_eta.push_back(gen_parton_eta);
        ftree->jet_genParton_phi.push_back(gen_parton_phi);
        ftree->jet_genParton_m.push_back(gen_parton_m);
        ftree->jet_genParton_E.push_back(gen_parton_E);
	
        ftree->jet_genParton_status.push_back(gen_parton_status);
        ftree->jet_genParton_id.push_back(gen_parton_id);
     }
   
   // Puppi Jets
   //
   
   bool jetPuppi_do = ftree->doWrite("jetPuppi_do");
   
   if( jetsPuppi.isValid() && jetPuppi_do )
     {
        int nJetPuppi = jetsPuppi->size();
        ftree->jetPuppi_n = nJetPuppi;
        for(int ij=0;ij<nJetPuppi;ij++)
	  {
	     const pat::Jet& jet = jetsPuppi->at(ij);
	     
	     ftree->jetPuppi_pt.push_back(jet.pt());
	     ftree->jetPuppi_eta.push_back(jet.eta());
	     ftree->jetPuppi_phi.push_back(jet.phi());
	     ftree->jetPuppi_m.push_back(jet.mass());
	     ftree->jetPuppi_E.push_back(jet.energy());
	     ftree->jetPuppi_JBP.push_back(jet.bDiscriminator("pfJetBProbabilityBJetTags"));
	     ftree->jetPuppi_JP.push_back(jet.bDiscriminator("pfJetProbabilityBJetTags"));
	     ftree->jetPuppi_TCHP.push_back(jet.bDiscriminator("pfTrackCountingHighPurBJetTags"));
	     ftree->jetPuppi_TCHE.push_back(jet.bDiscriminator("pfTrackCountingHighEffBJetTags"));
	     ftree->jetPuppi_SSVHE.push_back(jet.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"));
	     ftree->jetPuppi_SSVHP.push_back(jet.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
	     ftree->jetPuppi_CMVA.push_back(jet.bDiscriminator("pfCombinedMVABJetTags"));
	     
	     ftree->jetPuppi_chargedMultiplicity.push_back(jet.chargedMultiplicity());
	     ftree->jetPuppi_neutralMultiplicity.push_back(jet.neutralMultiplicity());
	     ftree->jetPuppi_chargedHadronMultiplicity.push_back(jet.chargedHadronMultiplicity());
	     
	     ftree->jetPuppi_jecFactorUncorrected.push_back(jet.jecFactor("Uncorrected"));
	     ftree->jetPuppi_jecFactorL1FastJet.push_back(-666.);
	     ftree->jetPuppi_jecFactorL2Relative.push_back(jet.jecFactor("L2Relative"));
	     ftree->jetPuppi_jecFactorL3Absolute.push_back(jet.jecFactor("L3Absolute"));
	     
	     ftree->jetPuppi_ntrk.push_back(jet.associatedTracks().size());
	     
	     float CSVIVF = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	     ftree->jetPuppi_CSVv2.push_back(CSVIVF);
	     
	     ftree->jetPuppi_partonFlavour.push_back(jet.partonFlavour());
	     ftree->jetPuppi_hadronFlavour.push_back(jet.hadronFlavour());
	     
	     ftree->jetPuppi_neutralHadronEnergy.push_back(jet.neutralHadronEnergy());
	     ftree->jetPuppi_neutralEmEnergy.push_back(jet.neutralEmEnergy());
	     ftree->jetPuppi_chargedHadronEnergy.push_back(jet.chargedHadronEnergy());
	     ftree->jetPuppi_chargedEmEnergy.push_back(jet.chargedEmEnergy());
	     ftree->jetPuppi_electronEnergy.push_back(jet.electronEnergy());
	     ftree->jetPuppi_muonEnergy.push_back(jet.muonEnergy());
	     ftree->jetPuppi_photonEnergy.push_back(jet.photonEnergy());
	     
	     if( jet.hasUserFloat("pileupJetId:fullDiscriminant") )
	       ftree->jetPuppi_pileupJetId.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
	     else
	       ftree->jetPuppi_pileupJetId.push_back(-666.);
	     
	     const reco::GenJet* genJet = jet.genJet();
	     bool hasGenInfo = (genJet);
	     ftree->jetPuppi_hasGenJet.push_back(hasGenInfo);
	     
	     float gen_jet_pt = -666;
	     float gen_jet_eta = -666;
	     float gen_jet_phi = -666;
	     float gen_jet_m = -666;
	     float gen_jet_E = -666;
	     int gen_jet_status = -666;
	     int gen_jet_id = -666;
	     
	     if( hasGenInfo )
	       {
		  gen_jet_pt = genJet->pt();
		  gen_jet_eta = genJet->eta();
		  gen_jet_phi = genJet->phi();
		  gen_jet_m = genJet->mass();
		  gen_jet_E = genJet->energy();
		  gen_jet_status = genJet->status();
		  gen_jet_id = genJet->pdgId();
	       }
	     
	     ftree->jetPuppi_genJet_pt.push_back(gen_jet_pt);
	     ftree->jetPuppi_genJet_eta.push_back(gen_jet_eta);
	     ftree->jetPuppi_genJet_phi.push_back(gen_jet_phi);
	     ftree->jetPuppi_genJet_m.push_back(gen_jet_m);
	     ftree->jetPuppi_genJet_E.push_back(gen_jet_E);
	     
	     ftree->jetPuppi_genJet_status.push_back(gen_jet_status);
	     ftree->jetPuppi_genJet_id.push_back(gen_jet_id);
	     
	     const reco::GenParticle* genParton = (!isData_) ? jet.genParton() : 0;
	     bool hasGenPartonInfo = (genParton);
	     ftree->jetPuppi_hasGenParton.push_back(hasGenPartonInfo);
	     
	     float gen_parton_pt = -666;
	     float gen_parton_eta = -666;
	     float gen_parton_phi = -666;
	     float gen_parton_m = -666;
	     float gen_parton_E = -666;
	     int gen_parton_status = -666;
	     int gen_parton_id = -666;
	     
	     if( hasGenPartonInfo )
	       {
		  gen_parton_pt = genParton->pt();
		  gen_parton_eta = genParton->eta();
		  gen_parton_phi = genParton->phi();
		  gen_parton_m = genParton->mass();
		  gen_parton_E = genParton->energy();
		  gen_parton_status = genParton->status();
		  gen_parton_id = genParton->pdgId();
	       }
	     
	     ftree->jetPuppi_genParton_pt.push_back(gen_parton_pt);
	     ftree->jetPuppi_genParton_eta.push_back(gen_parton_eta);
	     ftree->jetPuppi_genParton_phi.push_back(gen_parton_phi);
	     ftree->jetPuppi_genParton_m.push_back(gen_parton_m);
	     ftree->jetPuppi_genParton_E.push_back(gen_parton_E);
	     
	     ftree->jetPuppi_genParton_status.push_back(gen_parton_status);
	     ftree->jetPuppi_genParton_id.push_back(gen_parton_id);
	  }
     }
   
   // ak8 jets (W-jets)
    
   bool ak8jet_do = ftree->doWrite("ak8jet_do");
   
   if( ak8jets.isValid() && ak8jet_do )
     {
        int nak8Jet = ak8jets->size();
        ftree->ak8jet_n = nak8Jet;
        for(int ij=0;ij<nak8Jet;ij++)
	  {
	     const pat::Jet& jet = ak8jets->at(ij);
	     
	     ftree->ak8jet_pt.push_back(jet.pt());
	     ftree->ak8jet_eta.push_back(jet.eta());
	     ftree->ak8jet_phi.push_back(jet.phi());
	     ftree->ak8jet_m.push_back(jet.mass());
	     ftree->ak8jet_E.push_back(jet.energy());
	     
	     ftree->ak8jet_JBP.push_back(jet.bDiscriminator("pfJetBProbabilityBJetTags"));
	     ftree->ak8jet_JP.push_back(jet.bDiscriminator("pfJetProbabilityBJetTags"));
	     ftree->ak8jet_TCHP.push_back(jet.bDiscriminator("pfTrackCountingHighPurBJetTags"));
	     ftree->ak8jet_TCHE.push_back(jet.bDiscriminator("pfTrackCountingHighEffBJetTags"));
	     ftree->ak8jet_SSVHE.push_back(jet.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"));
	     ftree->ak8jet_SSVHP.push_back(jet.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
	     ftree->ak8jet_CMVA.push_back(jet.bDiscriminator("pfCombinedMVABJetTags"));
	     
	     ftree->ak8jet_chargedMultiplicity.push_back(jet.chargedMultiplicity());
	     ftree->ak8jet_neutralMultiplicity.push_back(jet.neutralMultiplicity());
	     ftree->ak8jet_chargedHadronMultiplicity.push_back(jet.chargedHadronMultiplicity());
	     
	     ftree->ak8jet_jecFactorUncorrected.push_back(jet.jecFactor("Uncorrected"));
	     ftree->ak8jet_jecFactorL1FastJet.push_back(jet.jecFactor("L1FastJet"));
	     ftree->ak8jet_jecFactorL2Relative.push_back(jet.jecFactor("L2Relative"));
	     ftree->ak8jet_jecFactorL3Absolute.push_back(jet.jecFactor("L3Absolute"));
	     
	     ftree->ak8jet_ntrk.push_back(jet.associatedTracks().size());
	     
	     float CSVIVF = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	     ftree->ak8jet_CSVv2.push_back(CSVIVF);
	     
	     ftree->ak8jet_partonFlavour.push_back(jet.partonFlavour());
	     ftree->ak8jet_hadronFlavour.push_back(jet.hadronFlavour());
	     
	     ftree->ak8jet_neutralHadronEnergy.push_back(jet.neutralHadronEnergy());
	     ftree->ak8jet_neutralEmEnergy.push_back(jet.neutralEmEnergy());
	     ftree->ak8jet_chargedHadronEnergy.push_back(jet.chargedHadronEnergy());
	     ftree->ak8jet_chargedEmEnergy.push_back(jet.chargedEmEnergy());
	     ftree->ak8jet_electronEnergy.push_back(jet.electronEnergy());
	     ftree->ak8jet_muonEnergy.push_back(jet.muonEnergy());
	     ftree->ak8jet_photonEnergy.push_back(jet.photonEnergy());
	     
	     if( jet.hasUserFloat("pileupJetId:fullDiscriminant") )
	       ftree->ak8jet_pileupJetId.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
	     else
	       ftree->ak8jet_pileupJetId.push_back(-666.);
	     
	     ftree->ak8jet_jetArea.push_back(jet.jetArea());
	     
	     // Jet ID
	     float NHF = jet.neutralHadronEnergyFraction();
	     float NEMF = jet.neutralEmEnergyFraction();
	     float CHF = jet.chargedHadronEnergyFraction();
	     float MUF = jet.muonEnergyFraction();
	     float CEMF = jet.chargedEmEnergyFraction();
	     float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	     float CHM = jet.chargedMultiplicity();
	     float eta = jet.eta();
	     //            bool looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8) && ((fabs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(eta)>2.4);
	     bool tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((fabs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || fabs(eta)>2.4);
	     
	     //            ftree->ak8jet_looseJetID.push_back(looseJetID);
	     ftree->ak8jet_tightJetID.push_back(tightJetID);
	     
	     const reco::GenJet* genJet = jet.genJet();
	     bool hasGenInfo = (genJet);
	     ftree->ak8jet_hasGenJet.push_back(hasGenInfo);
	     
	     float gen_jet_pt = -666;
	     float gen_jet_eta = -666;
	     float gen_jet_phi = -666;
	     float gen_jet_m = -666;
	     float gen_jet_E = -666;
	     int gen_jet_status = -666;
	     int gen_jet_id = -666;
	     
	     if( hasGenInfo )
	       {
		  gen_jet_pt = genJet->pt();
		  gen_jet_eta = genJet->eta();
		  gen_jet_phi = genJet->phi();
		  gen_jet_m = genJet->mass();
		  gen_jet_E = genJet->energy();
		  gen_jet_status = genJet->status();
		  gen_jet_id = genJet->pdgId();
	       }
	     
	     ftree->ak8jet_genJet_pt.push_back(gen_jet_pt);
	     ftree->ak8jet_genJet_eta.push_back(gen_jet_eta);
	     ftree->ak8jet_genJet_phi.push_back(gen_jet_phi);
	     ftree->ak8jet_genJet_m.push_back(gen_jet_m);
	     ftree->ak8jet_genJet_E.push_back(gen_jet_E);
	     
	     ftree->ak8jet_genJet_status.push_back(gen_jet_status);
	     ftree->ak8jet_genJet_id.push_back(gen_jet_id);
	     
	     const reco::GenParticle* genParton = (!isData_) ? jet.genParton() : 0;
	     bool hasGenPartonInfo = (genParton);
	     ftree->ak8jet_hasGenParton.push_back(hasGenPartonInfo);
	     
	     float gen_parton_pt = -666;
	     float gen_parton_eta = -666;
	     float gen_parton_phi = -666;
	     float gen_parton_m = -666;
	     float gen_parton_E = -666;
	     int gen_parton_status = -666;
	     int gen_parton_id = -666;
	     
	     if( hasGenPartonInfo )
	       {
		  gen_parton_pt = genParton->pt();
		  gen_parton_eta = genParton->eta();
		  gen_parton_phi = genParton->phi();
		  gen_parton_m = genParton->mass();
		  gen_parton_E = genParton->energy();
		  gen_parton_status = genParton->status();
		  gen_parton_id = genParton->pdgId();
	       }
	     
	     ftree->ak8jet_genParton_pt.push_back(gen_parton_pt);
	     ftree->ak8jet_genParton_eta.push_back(gen_parton_eta);
	     ftree->ak8jet_genParton_phi.push_back(gen_parton_phi);
	     ftree->ak8jet_genParton_m.push_back(gen_parton_m);
	     ftree->ak8jet_genParton_E.push_back(gen_parton_E);
	     
	     ftree->ak8jet_genParton_status.push_back(gen_parton_status);
	     ftree->ak8jet_genParton_id.push_back(gen_parton_id);
	     
	     // access to the W-tagging variables
	     //	     ftree->ak8jet_tau1.push_back(jet.userFloat("NjettinessAK8:tau1")); //
	     //	     ftree->ak8jet_tau2.push_back(jet.userFloat("NjettinessAK8:tau2")); // Access the n-subjettiness variables
	     //	     ftree->ak8jet_tau3.push_back(jet.userFloat("NjettinessAK8:tau3")); //
	     // not available since 7_6_4
	     ftree->ak8jet_tau1.push_back(-666.); //
	     ftree->ak8jet_tau2.push_back(-666.); // Access the n-subjettiness variables
	     ftree->ak8jet_tau3.push_back(-666.); //
	     
	     // not available since 7_6_4
	     //	     ftree->ak8jet_softdrop_mass.push_back(jet.userFloat("ak8PFJetsCHSSoftDropMass")); // access to filtered mass
	     //	     ftree->ak8jet_trimmed_mass.push_back(jet.userFloat("ak8PFJetsCHSTrimmedMass"));   // access to trimmed mass
	     //	     ftree->ak8jet_pruned_mass.push_back(jet.userFloat("ak8PFJetsCHSPrunedMass"));     // access to pruned mass
	     //	     ftree->ak8jet_filtered_mass.push_back(jet.userFloat("ak8PFJetsCHSFilteredMass")); // access to filtered mass
	     
	     ftree->ak8jet_softdrop_mass.push_back(-666.); // access to filtered mass
	     ftree->ak8jet_trimmed_mass.push_back(-666.);   // access to trimmed mass
	     ftree->ak8jet_pruned_mass.push_back(-666.);     // access to pruned mass
	     ftree->ak8jet_filtered_mass.push_back(-666.); // access to filtered mass
	     
	     // access to the top-tagging variables
	     reco::CATopJetTagInfo const * tagInfo =  dynamic_cast<reco::CATopJetTagInfo const *>( jet.tagInfo("caTop"));
	     if ( tagInfo != 0 )
	       {
		  ftree->ak8jet_minMass.push_back(tagInfo->properties().minMass);
		  ftree->ak8jet_topMass.push_back(tagInfo->properties().topMass);
		  ftree->ak8jet_nSubJets.push_back(tagInfo->properties().nSubJets);
	       }
	     else
	       {
		  ftree->ak8jet_minMass.push_back(-666);
		  ftree->ak8jet_topMass.push_back(-666);
		  ftree->ak8jet_nSubJets.push_back(-666);
	       }
	  }
     }
   
   // ak10 jets
    
   bool ak10jet_do = ftree->doWrite("ak10jet_do");
   
   if( ak10jets.isValid() && ak10jet_do )
     {
        int nak10Jet = ak10jets->size();
        ftree->ak10jet_n = nak10Jet;
        for(int ij=0;ij<nak10Jet;ij++)
	  {
	     const pat::Jet& jet = ak10jets->at(ij);
	     
	     ftree->ak10jet_pt.push_back(jet.pt());
	     ftree->ak10jet_eta.push_back(jet.eta());
	     ftree->ak10jet_phi.push_back(jet.phi());
	     ftree->ak10jet_m.push_back(jet.mass());
	     ftree->ak10jet_E.push_back(jet.energy());
	     
	     ftree->ak10jet_JBP.push_back(jet.bDiscriminator("pfJetBProbabilityBJetTags"));
	     ftree->ak10jet_JP.push_back(jet.bDiscriminator("pfJetProbabilityBJetTags"));
	     ftree->ak10jet_TCHP.push_back(jet.bDiscriminator("pfTrackCountingHighPurBJetTags"));
	     ftree->ak10jet_TCHE.push_back(jet.bDiscriminator("pfTrackCountingHighEffBJetTags"));
	     ftree->ak10jet_SSVHE.push_back(jet.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"));
	     ftree->ak10jet_SSVHP.push_back(jet.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
	     ftree->ak10jet_CMVA.push_back(jet.bDiscriminator("pfCombinedMVABJetTags"));
	     
	     ftree->ak10jet_chargedMultiplicity.push_back(jet.chargedMultiplicity());
	     ftree->ak10jet_neutralMultiplicity.push_back(jet.neutralMultiplicity());
	     ftree->ak10jet_chargedHadronMultiplicity.push_back(jet.chargedHadronMultiplicity());
	     
	     ftree->ak10jet_jecFactorUncorrected.push_back(jet.jecFactor("Uncorrected"));
	     ftree->ak10jet_jecFactorL1FastJet.push_back(jet.jecFactor("L1FastJet"));
	     ftree->ak10jet_jecFactorL2Relative.push_back(jet.jecFactor("L2Relative"));
	     ftree->ak10jet_jecFactorL3Absolute.push_back(-666.);
	     
	     ftree->ak10jet_ntrk.push_back(jet.associatedTracks().size());
	     
	     float CSVIVF = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	     ftree->ak10jet_CSVv2.push_back(CSVIVF);
	     
	     ftree->ak10jet_partonFlavour.push_back(jet.partonFlavour());
	     ftree->ak10jet_hadronFlavour.push_back(jet.hadronFlavour());
	     
	     ftree->ak10jet_neutralHadronEnergy.push_back(jet.neutralHadronEnergy());
	     ftree->ak10jet_neutralEmEnergy.push_back(jet.neutralEmEnergy());
	     ftree->ak10jet_chargedHadronEnergy.push_back(jet.chargedHadronEnergy());
	     ftree->ak10jet_chargedEmEnergy.push_back(jet.chargedEmEnergy());
	     ftree->ak10jet_electronEnergy.push_back(jet.electronEnergy());
	     ftree->ak10jet_muonEnergy.push_back(jet.muonEnergy());
	     ftree->ak10jet_photonEnergy.push_back(jet.photonEnergy());

	     if( jet.hasUserFloat("pileupJetId:fullDiscriminant") )
	       ftree->ak10jet_pileupJetId.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
	     else
	       ftree->ak10jet_pileupJetId.push_back(-666.);
	     
	     ftree->ak10jet_jetArea.push_back(jet.jetArea());
	     
	     // Jet ID
	     float NHF = jet.neutralHadronEnergyFraction();
	     float NEMF = jet.neutralEmEnergyFraction();
	     float CHF = jet.chargedHadronEnergyFraction();
	     float MUF = jet.muonEnergyFraction();
	     float CEMF = jet.chargedEmEnergyFraction();
	     float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	     float CHM = jet.chargedMultiplicity();
	     float eta = jet.eta();
	     //            bool looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8) && ((fabs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(eta)>2.4);
	     bool tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((fabs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || fabs(eta)>2.4);
	     
	     //            ftree->ak10jet_looseJetID.push_back(looseJetID);
	     ftree->ak10jet_tightJetID.push_back(tightJetID);
	     
	     const reco::GenJet* genJet = jet.genJet();
	     bool hasGenInfo = (genJet);
	     ftree->ak10jet_hasGenJet.push_back(hasGenInfo);
	     
	     float gen_jet_pt = -666;
	     float gen_jet_eta = -666;
	     float gen_jet_phi = -666;
	     float gen_jet_m = -666;
	     float gen_jet_E = -666;
	     int gen_jet_status = -666;
	     int gen_jet_id = -666;
	     
	     if( hasGenInfo )
	       {
		  gen_jet_pt = genJet->pt();
		  gen_jet_eta = genJet->eta();
		  gen_jet_phi = genJet->phi();
		  gen_jet_m = genJet->mass();
		  gen_jet_E = genJet->energy();
		  gen_jet_status = genJet->status();
		  gen_jet_id = genJet->pdgId();
	       }
	     
	     ftree->ak10jet_genJet_pt.push_back(gen_jet_pt);
	     ftree->ak10jet_genJet_eta.push_back(gen_jet_eta);
	     ftree->ak10jet_genJet_phi.push_back(gen_jet_phi);
	     ftree->ak10jet_genJet_m.push_back(gen_jet_m);
	     ftree->ak10jet_genJet_E.push_back(gen_jet_E);
	     
	     ftree->ak10jet_genJet_status.push_back(gen_jet_status);
	     ftree->ak10jet_genJet_id.push_back(gen_jet_id);
	     
	     const reco::GenParticle* genParton = (!isData_) ? jet.genParton() : 0;
	     bool hasGenPartonInfo = (genParton);
	     ftree->ak10jet_hasGenParton.push_back(hasGenPartonInfo);
	     
	     float gen_parton_pt = -666;
	     float gen_parton_eta = -666;
	     float gen_parton_phi = -666;
	     float gen_parton_m = -666;
	     float gen_parton_E = -666;
	     int gen_parton_status = -666;
	     int gen_parton_id = -666;
	     
	     if( hasGenPartonInfo )
	       {
		  gen_parton_pt = genParton->pt();
		  gen_parton_eta = genParton->eta();
		  gen_parton_phi = genParton->phi();
		  gen_parton_m = genParton->mass();
		  gen_parton_E = genParton->energy();
		  gen_parton_status = genParton->status();
		  gen_parton_id = genParton->pdgId();
	       }
	     
	     ftree->ak10jet_genParton_pt.push_back(gen_parton_pt);
	     ftree->ak10jet_genParton_eta.push_back(gen_parton_eta);
	     ftree->ak10jet_genParton_phi.push_back(gen_parton_phi);
	     ftree->ak10jet_genParton_m.push_back(gen_parton_m);
	     ftree->ak10jet_genParton_E.push_back(gen_parton_E);
	     
	     ftree->ak10jet_genParton_status.push_back(gen_parton_status);
	     ftree->ak10jet_genParton_id.push_back(gen_parton_id);
	     
	     // access to the W-tagging variables
	     ftree->ak10jet_tau1.push_back(jet.userFloat("NjettinessAK8:tau1")); //
	     ftree->ak10jet_tau2.push_back(jet.userFloat("NjettinessAK8:tau2")); // Access the n-subjettiness variables
	     ftree->ak10jet_tau3.push_back(jet.userFloat("NjettinessAK8:tau3")); //
	     
	     ftree->ak10jet_softdrop_mass.push_back(jet.userFloat("ak10PFJetsCHSSoftDropMass")); // access to filtered mass
	     ftree->ak10jet_trimmed_mass.push_back(jet.userFloat("ak10PFJetsCHSTrimmedMass"));   // access to trimmed mass
	     ftree->ak10jet_pruned_mass.push_back(jet.userFloat("ak10PFJetsCHSPrunedMass"));     // access to pruned mass
	     ftree->ak10jet_filtered_mass.push_back(jet.userFloat("ak10PFJetsCHSFilteredMass")); // access to filtered mass
	     
	     // access to the top-tagging variables
	     reco::CATopJetTagInfo const * tagInfo =  dynamic_cast<reco::CATopJetTagInfo const *>( jet.tagInfo("caTop"));
	     if ( tagInfo != 0 )
	       {
		  ftree->ak10jet_minMass.push_back(tagInfo->properties().minMass);
		  ftree->ak10jet_topMass.push_back(tagInfo->properties().topMass);
		  ftree->ak10jet_nSubJets.push_back(tagInfo->properties().nSubJets);
	       }
	     else
	       {
		  ftree->ak10jet_minMass.push_back(-666);
		  ftree->ak10jet_topMass.push_back(-666);
		  ftree->ak10jet_nSubJets.push_back(-666);
	       }
	  }
     }
   
   // GenJets
   //
   if( !isData_ )
     {
        int nGenJet = genJets->size();
        ftree->genJet_n = nGenJet;
        for(int ij=0;ij<nGenJet;ij++)
	  {
	     const reco::GenJet& genJet = genJets->at(ij);
	     
	     ftree->genJet_pt.push_back(genJet.pt());
	     ftree->genJet_eta.push_back(genJet.eta());
	     ftree->genJet_phi.push_back(genJet.phi());
	     ftree->genJet_m.push_back(genJet.mass());
	     ftree->genJet_E.push_back(genJet.energy());

	     ftree->genJet_emEnergy.push_back(genJet.emEnergy());
	     ftree->genJet_hadEnergy.push_back(genJet.hadEnergy());
	     ftree->genJet_invisibleEnergy.push_back(genJet.invisibleEnergy());
	     ftree->genJet_auxiliaryEnergy.push_back(genJet.auxiliaryEnergy());
	     
	     //	     RefToBase<reco::Jet> jetRef(RefToBaseProd<reco::Jet>(genJets),ij);
	     //	     int genJet_flavour = (*genJetFlavourMatching)[jetRef].getFlavour();
	     //	     ftree->genJet_flavour.push_back(genJet_flavour);
	     ftree->genJet_flavour.push_back(-666); // FIXME
	  }
     }
   
   // ##########################
   //   PF candidates
   // ##########################
   //
   
   bool pfcand_do = ftree->doWrite("pfcand_do");
   if( pfcand_do )
     {	   
	int nPfcand = pfcands->size();
	ftree->pfcand_n = nPfcand;
	bool do_sel_pfc = ftree->doWrite("sel_pfcand");
	
	ftree->pfch_loose_n = 0;
	ftree->pfch_loose_sumpt = 0;
	ftree->pfch_tight_n = 0;
	ftree->pfch_tight_sumpt = 0;
	
	for( const pat::PackedCandidate &pfc : *pfcands )
	  {
	     if( pfc.hasTrackDetails() )
	       {	     	
		  //make a selection based on TOP-15-017
		  if(pfc.charge()!=0 && pfc.pt()>0.5 && fabs(pfc.eta())<2.1)
		    {
		       if(fabs(pfc.dz())<1 && fabs(pfc.dz()/pfc.dzError())<10 && fabs(pfc.dxy())<3 && fabs(pfc.dxy()/pfc.dxyError())<10)
			 {
			    ftree->pfch_loose_n++;
			    ftree->pfch_loose_sumpt+=pfc.pt();
			 }
		       if(fabs(pfc.dz())<0.1 && fabs(pfc.dz()/pfc.dzError())<5 && fabs(pfc.dxy())<0.5 && fabs(pfc.dxy()/pfc.dxyError())<5)
			 {
			    ftree->pfch_tight_n++;
			    ftree->pfch_tight_sumpt+=pfc.pt();
			 }
		    }
		  
		  //compute track Iso
		  double trackIso = 0;
		  // run over all pf candidates
		  for( const pat::PackedCandidate &pfc2 : *pfcands )
		    {
		       if(&pfc==&pfc2) continue ; // do not count particle itself
		       if(pfc2.charge()==0) continue;
		       if(fabs(pfc2.dz())>0.1) continue;
		       if(pfc2.pt()<0) continue;
		       if(deltaR(pfc,pfc2)<0.3) trackIso+=pfc2.pt();
		    }
		  
		  if(do_sel_pfc)
		    {
		       if(pfc.charge()==0) continue;
		       if(abs(pfc.dz())>=0.1) continue;
		       if(pfc.pt()<=10) continue;
		       if(fabs(pfc.eta())>=2.4) continue;
		       if( (trackIso<6 && pfc.pt()>=60 ) || (trackIso/pfc.pt()<0.1 && pfc.pt()<60) )
			 {
			    ftree->pfcand_pt.push_back(pfc.pt());
			    ftree->pfcand_eta.push_back(pfc.eta());
			    ftree->pfcand_phi.push_back(pfc.phi());
			    ftree->pfcand_E.push_back(pfc.energy());
			    ftree->pfcand_charge.push_back(pfc.charge());
			    ftree->pfcand_id.push_back(pfc.pdgId());
			    ftree->pfcand_dz.push_back(pfc.dz());
			    ftree->pfcand_trackIso.push_back(trackIso);
			 }
		    }
		  else
		    {
		       ftree->pfcand_pt.push_back(pfc.pt());
		       ftree->pfcand_eta.push_back(pfc.eta());
		       ftree->pfcand_phi.push_back(pfc.phi());
		       ftree->pfcand_E.push_back(pfc.energy());
		       ftree->pfcand_charge.push_back(pfc.charge());
		       ftree->pfcand_id.push_back(pfc.pdgId());
		       ftree->pfcand_dz.push_back(pfc.dz());
		       ftree->pfcand_trackIso.push_back(trackIso);
		    }
	       }
	  }   
     }
   
   this->KeepEvent();
   if( (applyMETFilters_ && passMETFilters) || !applyMETFilters_ ){
      //std::cout<<"here we are !"<<std::endl;
      if(ftree->apply_presel){
	 //std::cout<<"apply preselection !"<<std::endl;
	 //std::cout<<"MET "<<ftree->met_pt<<" "<<ftree->presel_MET_min<<std::endl;
	 if(ftree->n_presel_electrons >= ftree->n_presel_electrons_min &&
	    ftree->n_presel_muons >= ftree->n_presel_muons_min &&
	    (ftree->n_presel_electrons+ftree->n_presel_muons) >= ftree->n_presel_leptons_min &&
	    ftree->n_presel_jets >= ftree->n_presel_jets_min &&
	    ftree->met_pt >= ftree->presel_MET_min){
	    //std::cout<<"pass "<<std::endl;
	    ftree->tree->Fill();
	 }
      }
      else{
	 ftree->tree->Fill();
      }
      
   }
   /*
    std::cout<<"nof jets passing the preselection: "<<ftree->n_presel_jets<<std::endl;
    std::cout<<"nof muons passing the preselection: "<<ftree->n_presel_muons<<std::endl;
    std::cout<<"nof electrons passing the preselection: "<<ftree->n_presel_electrons<<std::endl;
       */
   
   delete mc_truth;
   delete genTTX;
}

// ------------ method called when starting to processes a run  ------------
void FlatTreeProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
   bool changed = true;
   if(!hltConfig_.init(iRun, iSetup, "HLT", changed))
     std::cout << "Warning, didn't find HLTConfigProvider with label "
     << "HLT" << " in run " << iRun.run() << std::endl;
   
   if(!hltPrescale_.init(iRun, iSetup, "HLT", changed))
     std::cout << "Warning, didn't find HLTPrescaleProvider with label "
     << "HLT" << " in run " << iRun.run() << std::endl;
   
   /*   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    * 
    iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);
    //   iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl);
    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
    * 
    jecUnc = new JetCorrectionUncertainty(JetCorPar);*/
   
   const char* cmssw_base = std::getenv("CMSSW_BASE");
   std::string JECUncertaintyPath = std::string(cmssw_base)+"/src/VUBFlatTree/FlatTreeProducer/data/jecFiles/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_Uncertainty_AK4PFchs.txt";
   jecUnc = new JetCorrectionUncertainty(JECUncertaintyPath.c_str());
}

// ------------ method called when ending the processing of a run  ------------
void FlatTreeProducer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
   delete jecUnc;
}

void FlatTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

bool FlatTreeProducer::foundTrigger(const std::string& name) const
{
   for( unsigned int i=0;i<filterTriggerNames_.size();++i )
     {
        TString pattern(filterTriggerNames_[i]);
        pattern.ToLower();
        TRegexp reg(Form("%s",pattern.Data()),true);
        TString sname(name);
        sname.ToLower();
        if( sname.Index(reg) >= 0 ) return true;
     }

   return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(FlatTreeProducer);

