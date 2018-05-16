#include "VUBFlatTree/FlatTreeProducer/interface/MCTruth.hh"

void MCTruth::fillTopStopDecayChain(const edm::Event& iEvent,
				    const edm::EventSetup& iSetup,
				    FlatTree& tree,
				    const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   //@MJ@ TODO new method to fill only useful gen information for stop analysis
   // mening stop/top decay chain, radiation
   // not fully tested yet!!!! 
   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;
   
   int gen_n = 0;
   
   std::vector<float> gen_pt;
   std::vector<float> gen_eta;
   std::vector<float> gen_phi;
   std::vector<float> gen_m;
   std::vector<float> gen_E;
   std::vector<int> gen_id;
   std::vector<int> gen_status;
   std::vector<int> gen_charge;
   std::vector<int> gen_index;
   std::vector<int> gen_mother_index;
   std::vector<int> gen_daughter_n;
   std::vector<std::vector<int> > gen_daughter_index;
   
   std::vector<float> gen_stop_m;
   std::vector<float> gen_neutralino_m;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	float ptGen = mcp->pt();
	float etaGen = mcp->eta();
	float phiGen = mcp->phi();
	float mGen = mcp->mass();
	float EGen = mcp->energy();
	int idGen = mcp->pdgId();
	int statusGen = mcp->status();
	int chargeGen = mcp->charge();
	int indexGen = gen_n;
        float mStopGen = idGen == 1000006 ? mGen : -1;
        float mNeutralinoGen = idGen == 1000022 ? mGen : -1;
        //std::cout << "mStopGen = " << mStopGen << ", mNeutralinoGen = " << mNeutralinoGen << std::endl;

	const reco::GenParticle* mom = getMother(*mcp);

	reco::GenParticleCollection genParticlesCollection_m = *GenParticles;
	reco::GenParticleCollection::const_iterator genParticleSrc_m;
	
	int mother_index = 0;
	for(genParticleSrc_m = genParticlesCollection_m.begin();
	    genParticleSrc_m != genParticlesCollection_m.end();
	    genParticleSrc_m++)
	  {
	     reco::GenParticle *mcp_m = &(const_cast<reco::GenParticle&>(*genParticleSrc_m));
	     if( fabs(mcp_m->pt()-mom->pt()) < 10E-6 && fabs(mcp_m->eta()-mom->eta()) < 10E-6 )
	       {
		  break;
	       }		       
	     mother_index++;
	  }		  
	
	int daughter_n = 0;
	std::vector<int> daughter_index;
	
	const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
	for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
	  {
	     if( idr->isAvailable() ) 
	       {		       
		  const reco::GenParticleRef& genParticle = (*idr);
		  const reco::GenParticle *d = genParticle.get();

		  reco::GenParticleCollection genParticlesCollection_s = *GenParticles;
		  reco::GenParticleCollection::const_iterator genParticleSrc_s;
		  
		  int index = 0;
		  for(genParticleSrc_s = genParticlesCollection_s.begin();
		      genParticleSrc_s != genParticlesCollection_s.end();
		      genParticleSrc_s++)
		    {
		       reco::GenParticle *mcp_s = &(const_cast<reco::GenParticle&>(*genParticleSrc_s));
		       if( fabs(mcp_s->pt()-(*d).pt()) < 10E-6 && fabs(mcp_s->eta()-(*d).eta()) < 10E-6 )
			 {
			    break;
			 }		       
		       index++;
		    }		  
		  
		  daughter_index.push_back(index);
		  daughter_n++;
	       }
	  }
 

	gen_pt.push_back(ptGen);
	gen_eta.push_back(etaGen);
	gen_phi.push_back(phiGen);
	gen_m.push_back(mGen);
	gen_E.push_back(EGen);
	gen_id.push_back(idGen);
	gen_status.push_back(statusGen);
	gen_charge.push_back(chargeGen);
	gen_index.push_back(indexGen);
	gen_mother_index.push_back(mother_index);
	gen_daughter_n.push_back(daughter_n);
	gen_daughter_index.push_back(daughter_index);


        if(mStopGen != -1)
            gen_stop_m.push_back(mStopGen);
        if(mNeutralinoGen != -1)
            gen_neutralino_m.push_back(mNeutralinoGen);
        
        gen_n++;
     }
   
      
   int gen_n_slimmed = 0;
   std::vector<float> gen_pt_slimmed;
   std::vector<float> gen_eta_slimmed;
   std::vector<float> gen_phi_slimmed;
   std::vector<float> gen_m_slimmed;
   std::vector<float> gen_E_slimmed;
   std::vector<int> gen_id_slimmed;
   std::vector<int> gen_status_slimmed;
   std::vector<int> gen_charge_slimmed;
   std::vector<int> gen_index_slimmed;
   std::vector<int> gen_mother_index_slimmed;
   std::vector<int> gen_daughter_n_slimmed;
   std::vector<std::vector<int> > gen_daughter_index_slimmed;
      
   bool slimGenInfo = true;
   if(slimGenInfo == true)
   {
       for(uint32_t p = 0; p < gen_id.size(); p++)
       {

           //std::cout << "gen id size " << gen_id.size() << std::endl;
           bool stopTopDecay = (abs(gen_id.at(p)) == 6 || abs(gen_id.at(p)) == 1000006 || abs(gen_id.at(p)) == 5) ? true: false;
           bool charginoAndNeutralino = (abs(gen_id.at(p)) == 1000024 || abs(gen_id.at(p))==1000022) ? true: false;

           uint32_t mi = gen_mother_index.at(p);
           //std::cout << "mother index: " << mi << std::endl;
           uint32_t mmi = gen_mother_index.at(mi);
           //std::cout << "mother mother index: " << mmi << std::endl;
           bool WFromStopTop = (abs(gen_id.at(p)) == 24  && (abs(gen_id.at(mi)) == 6 || abs(gen_id.at(mi)) == 1000006)) ? true : false; 
           bool productsOfW = (abs(gen_id.at(mi)) == 24  && (abs(gen_id.at(mmi)) == 6 || abs(gen_id.at(mmi)) == 1000006)) ? true : false;
           bool leptonFromWNotFromStopTop = ( (abs(gen_id.at(p)) > 10 && abs(gen_id.at(p)) < 19)  && abs(gen_id.at(mi)) == 24  && (abs(gen_id.at(mmi)) != 6 || abs(gen_id.at(mmi)) != 1000006)) ? true : false;
           bool gluonRadiation = ((abs(gen_id.at(p)) == 9 || abs(gen_id.at(p)) == 21) && (abs(gen_id.at(mmi)) == 2212 || abs(gen_id.at(mmi)) == 6 || abs(gen_id.at(mmi)) == 1000006 )) ? true : false;
           bool gammaRadiation = (abs(gen_id.at(mi)) == 24  && gen_id.at(p) == 22) ? true : false;
       
          if(stopTopDecay || charginoAndNeutralino || WFromStopTop || productsOfW || leptonFromWNotFromStopTop || gluonRadiation || gammaRadiation)
          {

              gen_pt_slimmed.push_back(gen_pt.at(p));
              gen_eta_slimmed.push_back(gen_eta.at(p));
              gen_phi_slimmed.push_back(gen_phi.at(p));
              gen_m_slimmed.push_back(gen_m.at(p));
              gen_E_slimmed.push_back(gen_E.at(p));
              gen_id_slimmed.push_back(gen_id.at(p));
              gen_status_slimmed.push_back(gen_status.at(p));
              gen_charge_slimmed.push_back(gen_charge.at(p));
              gen_index_slimmed.push_back(gen_index.at(p));
              gen_mother_index_slimmed.push_back(gen_mother_index.at(p));
              gen_daughter_n_slimmed.push_back(gen_daughter_n.at(p));
              gen_daughter_index_slimmed.push_back(gen_daughter_index.at(p));

              gen_n_slimmed++;
       
         }
 
       }
   }
   else
   {
       gen_pt_slimmed = gen_pt;
       gen_eta_slimmed = gen_eta;
       gen_phi_slimmed = gen_phi;
       gen_m_slimmed = gen_m;
       gen_E_slimmed = gen_E;
       gen_id_slimmed = gen_id;
       gen_status_slimmed = gen_status;
       gen_charge_slimmed = gen_charge;
       gen_index_slimmed = gen_index;
       gen_mother_index_slimmed = gen_mother_index;
       gen_daughter_n_slimmed = gen_daughter_n;
       gen_daughter_index_slimmed = gen_daughter_index;
   }
   tree.gen_n = gen_n_slimmed;
   tree.gen_pt = gen_pt_slimmed;
   tree.gen_eta = gen_eta_slimmed;
   tree.gen_phi = gen_phi_slimmed;
   tree.gen_m = gen_m_slimmed;
   tree.gen_E = gen_E_slimmed;
   tree.gen_status = gen_status_slimmed;
   tree.gen_id = gen_id_slimmed;
   tree.gen_charge = gen_charge_slimmed;
   tree.gen_index = gen_index_slimmed;
   tree.gen_mother_index = gen_mother_index_slimmed;
   tree.gen_daughter_n = gen_daughter_n_slimmed;
   tree.gen_daughter_index = gen_daughter_index_slimmed;

   tree.gen_stop_m = gen_stop_m;
   tree.gen_neutralino_m = gen_neutralino_m;

}

void MCTruth::fillStopNeutralinoMass(const edm::Event& iEvent,
			       const edm::EventSetup& iSetup,
			       FlatTree& tree,
			       const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;
   std::vector<float> gen_stop_m;
   std::vector<float> gen_neutralino_m;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	float mGen = mcp->mass();
	int idGen = mcp->pdgId();
        float mStopGen = idGen == 1000006 ? mGen : -1;
        float mNeutralinoGen = idGen == 1000022 ? mGen : -1;
        //std::cout << "mStopGen = " << mStopGen << ", mNeutralinoGen = " << mNeutralinoGen << std::endl;

        if(mStopGen != -1)
            gen_stop_m.push_back(mStopGen);
        if(mNeutralinoGen != -1)
            gen_neutralino_m.push_back(mNeutralinoGen);
        
     }
   
   tree.gen_stop_m = gen_stop_m;
   tree.gen_neutralino_m = gen_neutralino_m;
}

void MCTruth::fillGenParticles(const edm::Event& iEvent,
			       const edm::EventSetup& iSetup,
			       FlatTree& tree,
			       const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;

   int gen_n = 0;
   
   std::vector<float> gen_pt;
   std::vector<float> gen_eta;
   std::vector<float> gen_phi;
   std::vector<float> gen_m;
   std::vector<float> gen_E;
   std::vector<int> gen_id;
   std::vector<int> gen_status;
   std::vector<int> gen_charge;
   std::vector<int> gen_index;
   std::vector<int> gen_mother_index;
   std::vector<int> gen_daughter_n;
   std::vector<std::vector<int> > gen_daughter_index;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	float ptGen = mcp->pt();
	float etaGen = mcp->eta();
	float phiGen = mcp->phi();
	float mGen = mcp->mass();
	float EGen = mcp->energy();
	int idGen = mcp->pdgId();
	int statusGen = mcp->status();
	int chargeGen = mcp->charge();
	int indexGen = gen_n;

	const reco::GenParticle* mom = getMother(*mcp);

	reco::GenParticleCollection genParticlesCollection_m = *GenParticles;
	reco::GenParticleCollection::const_iterator genParticleSrc_m;
	
	int mother_index = 0;
	for(genParticleSrc_m = genParticlesCollection_m.begin();
	    genParticleSrc_m != genParticlesCollection_m.end();
	    genParticleSrc_m++)
	  {
	     reco::GenParticle *mcp_m = &(const_cast<reco::GenParticle&>(*genParticleSrc_m));
	     if( fabs(mcp_m->pt()-mom->pt()) < 10E-6 && fabs(mcp_m->eta()-mom->eta()) < 10E-6 )
	       {
		  break;
	       }		       
	     mother_index++;
	  }		  
	
	int daughter_n = 0;
	std::vector<int> daughter_index;
	
	const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
	for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
	  {
	     if( idr->isAvailable() ) 
	       {		       
		  const reco::GenParticleRef& genParticle = (*idr);
		  const reco::GenParticle *d = genParticle.get();

		  reco::GenParticleCollection genParticlesCollection_s = *GenParticles;
		  reco::GenParticleCollection::const_iterator genParticleSrc_s;
		  
		  int index = 0;
		  for(genParticleSrc_s = genParticlesCollection_s.begin();
		      genParticleSrc_s != genParticlesCollection_s.end();
		      genParticleSrc_s++)
		    {
		       reco::GenParticle *mcp_s = &(const_cast<reco::GenParticle&>(*genParticleSrc_s));
		       if( fabs(mcp_s->pt()-(*d).pt()) < 10E-6 && fabs(mcp_s->eta()-(*d).eta()) < 10E-6 )
			 {
			    break;
			 }		       
		       index++;
		    }		  
		  
		  daughter_index.push_back(index);
		  daughter_n++;
	       }
	  }	

	gen_pt.push_back(ptGen);
	gen_eta.push_back(etaGen);
	gen_phi.push_back(phiGen);
	gen_m.push_back(mGen);
	gen_E.push_back(EGen);
	gen_id.push_back(idGen);
	gen_status.push_back(statusGen);
	gen_charge.push_back(chargeGen);
	gen_index.push_back(indexGen);
	gen_mother_index.push_back(mother_index);
	gen_daughter_n.push_back(daughter_n);
	gen_daughter_index.push_back(daughter_index);
	
	gen_n++;
     }
   
   tree.gen_n = gen_n;
   tree.gen_pt = gen_pt;
   tree.gen_eta = gen_eta;
   tree.gen_phi = gen_phi;
   tree.gen_m = gen_m;
   tree.gen_E = gen_E;
   tree.gen_status = gen_status;
   tree.gen_id = gen_id;
   tree.gen_charge = gen_charge;
   tree.gen_index = gen_index;
   tree.gen_mother_index = gen_mother_index;
   tree.gen_daughter_n = gen_daughter_n;
   tree.gen_daughter_index = gen_daughter_index;
}

void MCTruth::fillGenPV(const edm::Event& iEvent,
			const edm::EventSetup& iSetup,
			FlatTree& tree,
			const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;

   float gen_PVz = -666;
   for( size_t i=0;i<GenParticles->size();++i )
     {	
	const reco::GenParticle & genIt = (*GenParticles)[i];
	
	int status = genIt.status();
	if( (status>=21 && status<=29) || status==3 ) 
	  {
	     gen_PVz = genIt.vz();
	  }	
     }
   tree.gen_PVz = gen_PVz;
}

bool MCTruth::doMatch(const edm::Event& iEvent,
		      const edm::EventSetup& iSetup,
		      const edm::Handle<std::vector<reco::GenParticle> >& GenParticles,
		      reco::GenParticle *genp,
		      float &drMin,
		      float pt, float eta, float phi, int pdgId)
{
   bool foundMatch = 0;
   
   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;

   float drmin = 0.2;
   float ptRatMin = 0.5;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	float ptGen = mcp->pt();
	float etaGen = mcp->eta();
	float phiGen = mcp->phi();
	int idGen = mcp->pdgId();
	int statusGen = mcp->status();

	if( statusGen != 1 && statusGen != 3 ) continue;
	if( abs(pdgId) != abs(idGen) ) continue;
	
	float dr = GetDeltaR(eta,phi,etaGen,phiGen);
	float ptRat = (pt > 0.) ? fabs(pt-ptGen)/pt : 10E+10;
	
	if( dr < drmin && ptRat < ptRatMin )
	  {
	     drmin = dr;
	     foundMatch = 1;
	     genp = mcp;
	  }	
     }
   
   drMin = drmin;
   
   return foundMatch;
}

void MCTruth::Init(FlatTree &tree)
{
   tree.mc_truth_tth_channel = DEFVAL;

   // TLV
   
   tree.mc_truth_h0_p4.Clear();

   tree.mc_truth_h0W1_p4.Clear();
   tree.mc_truth_h0W2_p4.Clear();
   tree.mc_truth_h0Wl1_p4.Clear();
   tree.mc_truth_h0Wnu1_p4.Clear();
   tree.mc_truth_h0Wtau1_p4.Clear();
   tree.mc_truth_h0Wnutau1_p4.Clear();
   tree.mc_truth_h0Wtaul1_p4.Clear();
   tree.mc_truth_h0Wtaunu1_p4.Clear();
   tree.mc_truth_h0Wtaunutau1_p4.Clear();
   tree.mc_truth_h0Wl2_p4.Clear();
   tree.mc_truth_h0Wnu2_p4.Clear();
   tree.mc_truth_h0Wtau2_p4.Clear();
   tree.mc_truth_h0Wnutau2_p4.Clear();
   tree.mc_truth_h0Wtaul2_p4.Clear();
   tree.mc_truth_h0Wtaunu2_p4.Clear();
   tree.mc_truth_h0Wtaunutau2_p4.Clear();
   tree.mc_truth_h0Wq11_p4.Clear();
   tree.mc_truth_h0Wq21_p4.Clear();
   tree.mc_truth_h0Wq12_p4.Clear();
   tree.mc_truth_h0Wq22_p4.Clear();
   tree.mc_truth_h0Wq11_IS_p4.Clear();
   tree.mc_truth_h0Wq21_IS_p4.Clear();
   tree.mc_truth_h0Wq12_IS_p4.Clear();
   tree.mc_truth_h0Wq22_IS_p4.Clear();
   
   tree.mc_truth_h0Z1_p4.Clear();
   tree.mc_truth_h0Z2_p4.Clear();
   tree.mc_truth_h0Zl11_p4.Clear();
   tree.mc_truth_h0Zl21_p4.Clear();
   tree.mc_truth_h0Ztau11_p4.Clear();
   tree.mc_truth_h0Ztau21_p4.Clear();
   tree.mc_truth_h0Ztaul11_p4.Clear();
   tree.mc_truth_h0Ztaul21_p4.Clear();
   tree.mc_truth_h0Ztaunu11_p4.Clear();
   tree.mc_truth_h0Ztaunu21_p4.Clear();
   tree.mc_truth_h0Ztaunutau11_p4.Clear();
   tree.mc_truth_h0Ztaunutau21_p4.Clear();
   tree.mc_truth_h0Zq11_p4.Clear();
   tree.mc_truth_h0Zq21_p4.Clear();
   tree.mc_truth_h0Zq11_IS_p4.Clear();
   tree.mc_truth_h0Zq21_IS_p4.Clear();
   tree.mc_truth_h0Zl12_p4.Clear();
   tree.mc_truth_h0Zl22_p4.Clear();
   tree.mc_truth_h0Ztau12_p4.Clear();
   tree.mc_truth_h0Ztau22_p4.Clear();
   tree.mc_truth_h0Ztaul12_p4.Clear();
   tree.mc_truth_h0Ztaul22_p4.Clear();
   tree.mc_truth_h0Ztaunu12_p4.Clear();
   tree.mc_truth_h0Ztaunu22_p4.Clear();
   tree.mc_truth_h0Ztaunutau12_p4.Clear();
   tree.mc_truth_h0Ztaunutau22_p4.Clear();
   tree.mc_truth_h0Zq12_p4.Clear();
   tree.mc_truth_h0Zq22_p4.Clear();
   tree.mc_truth_h0Zq12_IS_p4.Clear();
   tree.mc_truth_h0Zq22_IS_p4.Clear();
   tree.mc_truth_h0Znu11_p4.Clear();
   tree.mc_truth_h0Znu21_p4.Clear();
   tree.mc_truth_h0Znu12_p4.Clear();
   tree.mc_truth_h0Znu22_p4.Clear();
   
   tree.mc_truth_h0tau1_p4.Clear();
   tree.mc_truth_h0tau2_p4.Clear();
   tree.mc_truth_h0taul1_p4.Clear();
   tree.mc_truth_h0taunutau1_p4.Clear();
   tree.mc_truth_h0taunu1_p4.Clear();
   tree.mc_truth_h0taul2_p4.Clear();
   tree.mc_truth_h0taunutau2_p4.Clear();
   tree.mc_truth_h0taunu2_p4.Clear();

   tree.mc_truth_h0b1_p4.Clear();
   tree.mc_truth_h0b2_p4.Clear();
   tree.mc_truth_h0b1_IS_p4.Clear();
   tree.mc_truth_h0b2_IS_p4.Clear();
   
   tree.mc_truth_t1_p4.Clear();
   tree.mc_truth_t2_p4.Clear();
   tree.mc_truth_tb1_p4.Clear();
   tree.mc_truth_tb2_p4.Clear();
   tree.mc_truth_tb1_IS_p4.Clear();
   tree.mc_truth_tb2_IS_p4.Clear();
   
   tree.mc_truth_tW1_p4.Clear();
   tree.mc_truth_tWnu1_p4.Clear();
   tree.mc_truth_tWnutau1_p4.Clear();
   tree.mc_truth_tWl1_p4.Clear();
   tree.mc_truth_tWtau1_p4.Clear();
   tree.mc_truth_tWtaunu1_p4.Clear();
   tree.mc_truth_tWtaunutau1_p4.Clear();
   tree.mc_truth_tWtaul1_p4.Clear();
   tree.mc_truth_tWq11_p4.Clear();
   tree.mc_truth_tWq21_p4.Clear();
   tree.mc_truth_tWq11_IS_p4.Clear();
   tree.mc_truth_tWq21_IS_p4.Clear();

   tree.mc_truth_tW2_p4.Clear();
   tree.mc_truth_tWnu2_p4.Clear();
   tree.mc_truth_tWnutau2_p4.Clear();
   tree.mc_truth_tWl2_p4.Clear();
   tree.mc_truth_tWtau2_p4.Clear();
   tree.mc_truth_tWtaunu2_p4.Clear();
   tree.mc_truth_tWtaunutau2_p4.Clear();
   tree.mc_truth_tWtaul2_p4.Clear();
   tree.mc_truth_tWq12_p4.Clear();
   tree.mc_truth_tWq22_p4.Clear();
   tree.mc_truth_tWq12_IS_p4.Clear();
   tree.mc_truth_tWq22_IS_p4.Clear();

   tree.mc_truth_j1_p4.Clear();
   tree.mc_truth_j2_p4.Clear();
   tree.mc_truth_j3_p4.Clear();
   
   tree.mc_truth_gammal1_p4.Clear();
   tree.mc_truth_gammal2_p4.Clear();
   tree.mc_truth_gammatau1_p4.Clear();
   tree.mc_truth_gammatau2_p4.Clear();
   tree.mc_truth_gammataul1_p4.Clear();
   tree.mc_truth_gammataul2_p4.Clear();
   tree.mc_truth_gammataunu1_p4.Clear();
   tree.mc_truth_gammataunu2_p4.Clear();
   tree.mc_truth_gammataunutau1_p4.Clear();
   tree.mc_truth_gammataunutau2_p4.Clear();

   // pt

   tree.mc_truth_h0_pt = DEFVAL;

   tree.mc_truth_h0W1_pt = DEFVAL;
   tree.mc_truth_h0W2_pt = DEFVAL;
   tree.mc_truth_h0Wl1_pt = DEFVAL;
   tree.mc_truth_h0Wnu1_pt = DEFVAL;
   tree.mc_truth_h0Wtau1_pt = DEFVAL;
   tree.mc_truth_h0Wnutau1_pt = DEFVAL;
   tree.mc_truth_h0Wtaul1_pt = DEFVAL;
   tree.mc_truth_h0Wtaunu1_pt = DEFVAL;
   tree.mc_truth_h0Wtaunutau1_pt = DEFVAL;
   tree.mc_truth_h0Wl2_pt = DEFVAL;
   tree.mc_truth_h0Wnu2_pt = DEFVAL;
   tree.mc_truth_h0Wtau2_pt = DEFVAL;
   tree.mc_truth_h0Wnutau2_pt = DEFVAL;
   tree.mc_truth_h0Wtaul2_pt = DEFVAL;
   tree.mc_truth_h0Wtaunu2_pt = DEFVAL;
   tree.mc_truth_h0Wtaunutau2_pt = DEFVAL;
   tree.mc_truth_h0Wq11_pt = DEFVAL;
   tree.mc_truth_h0Wq21_pt = DEFVAL;
   tree.mc_truth_h0Wq12_pt = DEFVAL;
   tree.mc_truth_h0Wq22_pt = DEFVAL;
   tree.mc_truth_h0Wq11_IS_pt = DEFVAL;
   tree.mc_truth_h0Wq21_IS_pt = DEFVAL;
   tree.mc_truth_h0Wq12_IS_pt = DEFVAL;
   tree.mc_truth_h0Wq22_IS_pt = DEFVAL;
   
   tree.mc_truth_h0Z1_pt = DEFVAL;
   tree.mc_truth_h0Z2_pt = DEFVAL;
   tree.mc_truth_h0Zl11_pt = DEFVAL;
   tree.mc_truth_h0Zl21_pt = DEFVAL;
   tree.mc_truth_h0Ztau11_pt = DEFVAL;
   tree.mc_truth_h0Ztau21_pt = DEFVAL;
   tree.mc_truth_h0Ztaul11_pt = DEFVAL;
   tree.mc_truth_h0Ztaul21_pt = DEFVAL;
   tree.mc_truth_h0Ztaunu11_pt = DEFVAL;
   tree.mc_truth_h0Ztaunu21_pt = DEFVAL;
   tree.mc_truth_h0Ztaunutau11_pt = DEFVAL;
   tree.mc_truth_h0Ztaunutau21_pt = DEFVAL;
   tree.mc_truth_h0Zq11_pt = DEFVAL;
   tree.mc_truth_h0Zq21_pt = DEFVAL;
   tree.mc_truth_h0Zq11_IS_pt = DEFVAL;
   tree.mc_truth_h0Zq21_IS_pt = DEFVAL;
   tree.mc_truth_h0Zl12_pt = DEFVAL;
   tree.mc_truth_h0Zl22_pt = DEFVAL;
   tree.mc_truth_h0Ztau12_pt = DEFVAL;
   tree.mc_truth_h0Ztau22_pt = DEFVAL;
   tree.mc_truth_h0Ztaul12_pt = DEFVAL;
   tree.mc_truth_h0Ztaul22_pt = DEFVAL;
   tree.mc_truth_h0Ztaunu12_pt = DEFVAL;
   tree.mc_truth_h0Ztaunu22_pt = DEFVAL;
   tree.mc_truth_h0Ztaunutau12_pt = DEFVAL;
   tree.mc_truth_h0Ztaunutau22_pt = DEFVAL;
   tree.mc_truth_h0Zq12_pt = DEFVAL;
   tree.mc_truth_h0Zq22_pt = DEFVAL;
   tree.mc_truth_h0Zq12_IS_pt = DEFVAL;
   tree.mc_truth_h0Zq22_IS_pt = DEFVAL;
   tree.mc_truth_h0Znu11_pt = DEFVAL;
   tree.mc_truth_h0Znu21_pt = DEFVAL;
   tree.mc_truth_h0Znu12_pt = DEFVAL;
   tree.mc_truth_h0Znu22_pt = DEFVAL;
   
   tree.mc_truth_h0tau1_pt = DEFVAL;
   tree.mc_truth_h0tau2_pt = DEFVAL;
   tree.mc_truth_h0taul1_pt = DEFVAL;
   tree.mc_truth_h0taunutau1_pt = DEFVAL;
   tree.mc_truth_h0taunu1_pt = DEFVAL;
   tree.mc_truth_h0taul2_pt = DEFVAL;
   tree.mc_truth_h0taunutau2_pt = DEFVAL;
   tree.mc_truth_h0taunu2_pt = DEFVAL;

   tree.mc_truth_h0b1_pt = DEFVAL;
   tree.mc_truth_h0b2_pt = DEFVAL;
   tree.mc_truth_h0b1_IS_pt = DEFVAL;
   tree.mc_truth_h0b2_IS_pt = DEFVAL;
   
   tree.mc_truth_t1_pt = DEFVAL;
   tree.mc_truth_t2_pt = DEFVAL;
   tree.mc_truth_tb1_pt = DEFVAL;
   tree.mc_truth_tb2_pt = DEFVAL;
   tree.mc_truth_tb1_IS_pt = DEFVAL;
   tree.mc_truth_tb2_IS_pt = DEFVAL;
   
   tree.mc_truth_tW1_pt = DEFVAL;
   tree.mc_truth_tWnu1_pt = DEFVAL;
   tree.mc_truth_tWnutau1_pt = DEFVAL;
   tree.mc_truth_tWl1_pt = DEFVAL;
   tree.mc_truth_tWtau1_pt = DEFVAL;
   tree.mc_truth_tWtaunu1_pt = DEFVAL;
   tree.mc_truth_tWtaunutau1_pt = DEFVAL;
   tree.mc_truth_tWtaul1_pt = DEFVAL;
   tree.mc_truth_tWq11_pt = DEFVAL;
   tree.mc_truth_tWq21_pt = DEFVAL;
   tree.mc_truth_tWq11_IS_pt = DEFVAL;
   tree.mc_truth_tWq21_IS_pt = DEFVAL;

   tree.mc_truth_tW2_pt = DEFVAL;
   tree.mc_truth_tWnu2_pt = DEFVAL;
   tree.mc_truth_tWnutau2_pt = DEFVAL;
   tree.mc_truth_tWl2_pt = DEFVAL;
   tree.mc_truth_tWtau2_pt = DEFVAL;
   tree.mc_truth_tWtaunu2_pt = DEFVAL;
   tree.mc_truth_tWtaunutau2_pt = DEFVAL;
   tree.mc_truth_tWtaul2_pt = DEFVAL;
   tree.mc_truth_tWq12_pt = DEFVAL;
   tree.mc_truth_tWq22_pt = DEFVAL;
   tree.mc_truth_tWq12_IS_pt = DEFVAL;
   tree.mc_truth_tWq22_IS_pt = DEFVAL;

   tree.mc_truth_j1_pt = DEFVAL;
   tree.mc_truth_j2_pt = DEFVAL;
   tree.mc_truth_j3_pt = DEFVAL;
   
   tree.mc_truth_gammal1_pt = DEFVAL;
   tree.mc_truth_gammal2_pt = DEFVAL;
   tree.mc_truth_gammatau1_pt = DEFVAL;
   tree.mc_truth_gammatau2_pt = DEFVAL;
   tree.mc_truth_gammataul1_pt = DEFVAL;
   tree.mc_truth_gammataul2_pt = DEFVAL;
   tree.mc_truth_gammataunu1_pt = DEFVAL;
   tree.mc_truth_gammataunu2_pt = DEFVAL;
   tree.mc_truth_gammataunutau1_pt = DEFVAL;
   tree.mc_truth_gammataunutau2_pt = DEFVAL;

   // eta

   tree.mc_truth_h0_eta = DEFVAL;

   tree.mc_truth_h0W1_eta = DEFVAL;
   tree.mc_truth_h0W2_eta = DEFVAL;
   tree.mc_truth_h0Wl1_eta = DEFVAL;
   tree.mc_truth_h0Wnu1_eta = DEFVAL;
   tree.mc_truth_h0Wtau1_eta = DEFVAL;
   tree.mc_truth_h0Wnutau1_eta = DEFVAL;
   tree.mc_truth_h0Wtaul1_eta = DEFVAL;
   tree.mc_truth_h0Wtaunu1_eta = DEFVAL;
   tree.mc_truth_h0Wtaunutau1_eta = DEFVAL;
   tree.mc_truth_h0Wl2_eta = DEFVAL;
   tree.mc_truth_h0Wnu2_eta = DEFVAL;
   tree.mc_truth_h0Wtau2_eta = DEFVAL;
   tree.mc_truth_h0Wnutau2_eta = DEFVAL;
   tree.mc_truth_h0Wtaul2_eta = DEFVAL;
   tree.mc_truth_h0Wtaunu2_eta = DEFVAL;
   tree.mc_truth_h0Wtaunutau2_eta = DEFVAL;
   tree.mc_truth_h0Wq11_eta = DEFVAL;
   tree.mc_truth_h0Wq21_eta = DEFVAL;
   tree.mc_truth_h0Wq12_eta = DEFVAL;
   tree.mc_truth_h0Wq22_eta = DEFVAL;
   tree.mc_truth_h0Wq11_IS_eta = DEFVAL;
   tree.mc_truth_h0Wq21_IS_eta = DEFVAL;
   tree.mc_truth_h0Wq12_IS_eta = DEFVAL;
   tree.mc_truth_h0Wq22_IS_eta = DEFVAL;
   
   tree.mc_truth_h0Z1_eta = DEFVAL;
   tree.mc_truth_h0Z2_eta = DEFVAL;
   tree.mc_truth_h0Zl11_eta = DEFVAL;
   tree.mc_truth_h0Zl21_eta = DEFVAL;
   tree.mc_truth_h0Ztau11_eta = DEFVAL;
   tree.mc_truth_h0Ztau21_eta = DEFVAL;
   tree.mc_truth_h0Ztaul11_eta = DEFVAL;
   tree.mc_truth_h0Ztaul21_eta = DEFVAL;
   tree.mc_truth_h0Ztaunu11_eta = DEFVAL;
   tree.mc_truth_h0Ztaunu21_eta = DEFVAL;
   tree.mc_truth_h0Ztaunutau11_eta = DEFVAL;
   tree.mc_truth_h0Ztaunutau21_eta = DEFVAL;
   tree.mc_truth_h0Zq11_eta = DEFVAL;
   tree.mc_truth_h0Zq21_eta = DEFVAL;
   tree.mc_truth_h0Zq11_IS_eta = DEFVAL;
   tree.mc_truth_h0Zq21_IS_eta = DEFVAL;
   tree.mc_truth_h0Zl12_eta = DEFVAL;
   tree.mc_truth_h0Zl22_eta = DEFVAL;
   tree.mc_truth_h0Ztau12_eta = DEFVAL;
   tree.mc_truth_h0Ztau22_eta = DEFVAL;
   tree.mc_truth_h0Ztaul12_eta = DEFVAL;
   tree.mc_truth_h0Ztaul22_eta = DEFVAL;
   tree.mc_truth_h0Ztaunu12_eta = DEFVAL;
   tree.mc_truth_h0Ztaunu22_eta = DEFVAL;
   tree.mc_truth_h0Ztaunutau12_eta = DEFVAL;
   tree.mc_truth_h0Ztaunutau22_eta = DEFVAL;
   tree.mc_truth_h0Zq12_eta = DEFVAL;
   tree.mc_truth_h0Zq22_eta = DEFVAL;
   tree.mc_truth_h0Zq12_IS_eta = DEFVAL;
   tree.mc_truth_h0Zq22_IS_eta = DEFVAL;
   tree.mc_truth_h0Znu11_eta = DEFVAL;
   tree.mc_truth_h0Znu21_eta = DEFVAL;
   tree.mc_truth_h0Znu12_eta = DEFVAL;
   tree.mc_truth_h0Znu22_eta = DEFVAL;
   
   tree.mc_truth_h0tau1_eta = DEFVAL;
   tree.mc_truth_h0tau2_eta = DEFVAL;
   tree.mc_truth_h0taul1_eta = DEFVAL;
   tree.mc_truth_h0taunutau1_eta = DEFVAL;
   tree.mc_truth_h0taunu1_eta = DEFVAL;
   tree.mc_truth_h0taul2_eta = DEFVAL;
   tree.mc_truth_h0taunutau2_eta = DEFVAL;
   tree.mc_truth_h0taunu2_eta = DEFVAL;

   tree.mc_truth_h0b1_eta = DEFVAL;
   tree.mc_truth_h0b2_eta = DEFVAL;
   tree.mc_truth_h0b1_IS_eta = DEFVAL;
   tree.mc_truth_h0b2_IS_eta = DEFVAL;
   
   tree.mc_truth_t1_eta = DEFVAL;
   tree.mc_truth_t2_eta = DEFVAL;
   tree.mc_truth_tb1_eta = DEFVAL;
   tree.mc_truth_tb2_eta = DEFVAL;
   tree.mc_truth_tb1_IS_eta = DEFVAL;
   tree.mc_truth_tb2_IS_eta = DEFVAL;
   
   tree.mc_truth_tW1_eta = DEFVAL;
   tree.mc_truth_tWnu1_eta = DEFVAL;
   tree.mc_truth_tWnutau1_eta = DEFVAL;
   tree.mc_truth_tWl1_eta = DEFVAL;
   tree.mc_truth_tWtau1_eta = DEFVAL;
   tree.mc_truth_tWtaunu1_eta = DEFVAL;
   tree.mc_truth_tWtaunutau1_eta = DEFVAL;
   tree.mc_truth_tWtaul1_eta = DEFVAL;
   tree.mc_truth_tWq11_eta = DEFVAL;
   tree.mc_truth_tWq21_eta = DEFVAL;
   tree.mc_truth_tWq11_IS_eta = DEFVAL;
   tree.mc_truth_tWq21_IS_eta = DEFVAL;

   tree.mc_truth_tW2_eta = DEFVAL;
   tree.mc_truth_tWnu2_eta = DEFVAL;
   tree.mc_truth_tWnutau2_eta = DEFVAL;
   tree.mc_truth_tWl2_eta = DEFVAL;
   tree.mc_truth_tWtau2_eta = DEFVAL;
   tree.mc_truth_tWtaunu2_eta = DEFVAL;
   tree.mc_truth_tWtaunutau2_eta = DEFVAL;
   tree.mc_truth_tWtaul2_eta = DEFVAL;
   tree.mc_truth_tWq12_eta = DEFVAL;
   tree.mc_truth_tWq22_eta = DEFVAL;
   tree.mc_truth_tWq12_IS_eta = DEFVAL;
   tree.mc_truth_tWq22_IS_eta = DEFVAL;

   tree.mc_truth_j1_eta = DEFVAL;
   tree.mc_truth_j2_eta = DEFVAL;
   tree.mc_truth_j3_eta = DEFVAL;
   
   tree.mc_truth_gammal1_eta = DEFVAL;
   tree.mc_truth_gammal2_eta = DEFVAL;
   tree.mc_truth_gammatau1_eta = DEFVAL;
   tree.mc_truth_gammatau2_eta = DEFVAL;
   tree.mc_truth_gammataul1_eta = DEFVAL;
   tree.mc_truth_gammataul2_eta = DEFVAL;
   tree.mc_truth_gammataunu1_eta = DEFVAL;
   tree.mc_truth_gammataunu2_eta = DEFVAL;
   tree.mc_truth_gammataunutau1_eta = DEFVAL;
   tree.mc_truth_gammataunutau2_eta = DEFVAL;
   
   // phi
   
   tree.mc_truth_h0_phi = DEFVAL;

   tree.mc_truth_h0W1_phi = DEFVAL;
   tree.mc_truth_h0W2_phi = DEFVAL;
   tree.mc_truth_h0Wl1_phi = DEFVAL;
   tree.mc_truth_h0Wnu1_phi = DEFVAL;
   tree.mc_truth_h0Wtau1_phi = DEFVAL;
   tree.mc_truth_h0Wnutau1_phi = DEFVAL;
   tree.mc_truth_h0Wtaul1_phi = DEFVAL;
   tree.mc_truth_h0Wtaunu1_phi = DEFVAL;
   tree.mc_truth_h0Wtaunutau1_phi = DEFVAL;
   tree.mc_truth_h0Wl2_phi = DEFVAL;
   tree.mc_truth_h0Wnu2_phi = DEFVAL;
   tree.mc_truth_h0Wtau2_phi = DEFVAL;
   tree.mc_truth_h0Wnutau2_phi = DEFVAL;
   tree.mc_truth_h0Wtaul2_phi = DEFVAL;
   tree.mc_truth_h0Wtaunu2_phi = DEFVAL;
   tree.mc_truth_h0Wtaunutau2_phi = DEFVAL;
   tree.mc_truth_h0Wq11_phi = DEFVAL;
   tree.mc_truth_h0Wq21_phi = DEFVAL;
   tree.mc_truth_h0Wq12_phi = DEFVAL;
   tree.mc_truth_h0Wq22_phi = DEFVAL;
   tree.mc_truth_h0Wq11_IS_phi = DEFVAL;
   tree.mc_truth_h0Wq21_IS_phi = DEFVAL;
   tree.mc_truth_h0Wq12_IS_phi = DEFVAL;
   tree.mc_truth_h0Wq22_IS_phi = DEFVAL;
   
   tree.mc_truth_h0Z1_phi = DEFVAL;
   tree.mc_truth_h0Z2_phi = DEFVAL;
   tree.mc_truth_h0Zl11_phi = DEFVAL;
   tree.mc_truth_h0Zl21_phi = DEFVAL;
   tree.mc_truth_h0Ztau11_phi = DEFVAL;
   tree.mc_truth_h0Ztau21_phi = DEFVAL;
   tree.mc_truth_h0Ztaul11_phi = DEFVAL;
   tree.mc_truth_h0Ztaul21_phi = DEFVAL;
   tree.mc_truth_h0Ztaunu11_phi = DEFVAL;
   tree.mc_truth_h0Ztaunu21_phi = DEFVAL;
   tree.mc_truth_h0Ztaunutau11_phi = DEFVAL;
   tree.mc_truth_h0Ztaunutau21_phi = DEFVAL;
   tree.mc_truth_h0Zq11_phi = DEFVAL;
   tree.mc_truth_h0Zq21_phi = DEFVAL;
   tree.mc_truth_h0Zq11_IS_phi = DEFVAL;
   tree.mc_truth_h0Zq21_IS_phi = DEFVAL;
   tree.mc_truth_h0Zl12_phi = DEFVAL;
   tree.mc_truth_h0Zl22_phi = DEFVAL;
   tree.mc_truth_h0Ztau12_phi = DEFVAL;
   tree.mc_truth_h0Ztau22_phi = DEFVAL;
   tree.mc_truth_h0Ztaul12_phi = DEFVAL;
   tree.mc_truth_h0Ztaul22_phi = DEFVAL;
   tree.mc_truth_h0Ztaunu12_phi = DEFVAL;
   tree.mc_truth_h0Ztaunu22_phi = DEFVAL;
   tree.mc_truth_h0Ztaunutau12_phi = DEFVAL;
   tree.mc_truth_h0Ztaunutau22_phi = DEFVAL;
   tree.mc_truth_h0Zq12_phi = DEFVAL;
   tree.mc_truth_h0Zq22_phi = DEFVAL;
   tree.mc_truth_h0Zq12_IS_phi = DEFVAL;
   tree.mc_truth_h0Zq22_IS_phi = DEFVAL;
   tree.mc_truth_h0Znu11_phi = DEFVAL;
   tree.mc_truth_h0Znu21_phi = DEFVAL;
   tree.mc_truth_h0Znu12_phi = DEFVAL;
   tree.mc_truth_h0Znu22_phi = DEFVAL;
   
   tree.mc_truth_h0tau1_phi = DEFVAL;
   tree.mc_truth_h0tau2_phi = DEFVAL;
   tree.mc_truth_h0taul1_phi = DEFVAL;
   tree.mc_truth_h0taunutau1_phi = DEFVAL;
   tree.mc_truth_h0taunu1_phi = DEFVAL;
   tree.mc_truth_h0taul2_phi = DEFVAL;
   tree.mc_truth_h0taunutau2_phi = DEFVAL;
   tree.mc_truth_h0taunu2_phi = DEFVAL;

   tree.mc_truth_h0b1_phi = DEFVAL;
   tree.mc_truth_h0b2_phi = DEFVAL;
   tree.mc_truth_h0b1_IS_phi = DEFVAL;
   tree.mc_truth_h0b2_IS_phi = DEFVAL;
   
   tree.mc_truth_t1_phi = DEFVAL;
   tree.mc_truth_t2_phi = DEFVAL;
   tree.mc_truth_tb1_phi = DEFVAL;
   tree.mc_truth_tb2_phi = DEFVAL;
   tree.mc_truth_tb1_IS_phi = DEFVAL;
   tree.mc_truth_tb2_IS_phi = DEFVAL;
   
   tree.mc_truth_tW1_phi = DEFVAL;
   tree.mc_truth_tWnu1_phi = DEFVAL;
   tree.mc_truth_tWnutau1_phi = DEFVAL;
   tree.mc_truth_tWl1_phi = DEFVAL;
   tree.mc_truth_tWtau1_phi = DEFVAL;
   tree.mc_truth_tWtaunu1_phi = DEFVAL;
   tree.mc_truth_tWtaunutau1_phi = DEFVAL;
   tree.mc_truth_tWtaul1_phi = DEFVAL;
   tree.mc_truth_tWq11_phi = DEFVAL;
   tree.mc_truth_tWq21_phi = DEFVAL;
   tree.mc_truth_tWq11_IS_phi = DEFVAL;
   tree.mc_truth_tWq21_IS_phi = DEFVAL;

   tree.mc_truth_tW2_phi = DEFVAL;
   tree.mc_truth_tWnu2_phi = DEFVAL;
   tree.mc_truth_tWnutau2_phi = DEFVAL;
   tree.mc_truth_tWl2_phi = DEFVAL;
   tree.mc_truth_tWtau2_phi = DEFVAL;
   tree.mc_truth_tWtaunu2_phi = DEFVAL;
   tree.mc_truth_tWtaunutau2_phi = DEFVAL;
   tree.mc_truth_tWtaul2_phi = DEFVAL;
   tree.mc_truth_tWq12_phi = DEFVAL;
   tree.mc_truth_tWq22_phi = DEFVAL;
   tree.mc_truth_tWq12_IS_phi = DEFVAL;
   tree.mc_truth_tWq22_IS_phi = DEFVAL;

   tree.mc_truth_j1_phi = DEFVAL;
   tree.mc_truth_j2_phi = DEFVAL;
   tree.mc_truth_j3_phi = DEFVAL;
   
   tree.mc_truth_gammal1_phi = DEFVAL;
   tree.mc_truth_gammal2_phi = DEFVAL;
   tree.mc_truth_gammatau1_phi = DEFVAL;
   tree.mc_truth_gammatau2_phi = DEFVAL;
   tree.mc_truth_gammataul1_phi = DEFVAL;
   tree.mc_truth_gammataul2_phi = DEFVAL;
   tree.mc_truth_gammataunu1_phi = DEFVAL;
   tree.mc_truth_gammataunu2_phi = DEFVAL;
   tree.mc_truth_gammataunutau1_phi = DEFVAL;
   tree.mc_truth_gammataunutau2_phi = DEFVAL;
  
   // E
   
   tree.mc_truth_h0_E = DEFVAL;

   tree.mc_truth_h0W1_E = DEFVAL;
   tree.mc_truth_h0W2_E = DEFVAL;
   tree.mc_truth_h0Wl1_E = DEFVAL;
   tree.mc_truth_h0Wnu1_E = DEFVAL;
   tree.mc_truth_h0Wtau1_E = DEFVAL;
   tree.mc_truth_h0Wnutau1_E = DEFVAL;
   tree.mc_truth_h0Wtaul1_E = DEFVAL;
   tree.mc_truth_h0Wtaunu1_E = DEFVAL;
   tree.mc_truth_h0Wtaunutau1_E = DEFVAL;
   tree.mc_truth_h0Wl2_E = DEFVAL;
   tree.mc_truth_h0Wnu2_E = DEFVAL;
   tree.mc_truth_h0Wtau2_E = DEFVAL;
   tree.mc_truth_h0Wnutau2_E = DEFVAL;
   tree.mc_truth_h0Wtaul2_E = DEFVAL;
   tree.mc_truth_h0Wtaunu2_E = DEFVAL;
   tree.mc_truth_h0Wtaunutau2_E = DEFVAL;
   tree.mc_truth_h0Wq11_E = DEFVAL;
   tree.mc_truth_h0Wq21_E = DEFVAL;
   tree.mc_truth_h0Wq12_E = DEFVAL;
   tree.mc_truth_h0Wq22_E = DEFVAL;
   tree.mc_truth_h0Wq11_IS_E = DEFVAL;
   tree.mc_truth_h0Wq21_IS_E = DEFVAL;
   tree.mc_truth_h0Wq12_IS_E = DEFVAL;
   tree.mc_truth_h0Wq22_IS_E = DEFVAL;
  
   tree.mc_truth_h0Z1_E = DEFVAL;
   tree.mc_truth_h0Z2_E = DEFVAL;
   tree.mc_truth_h0Zl11_E = DEFVAL;
   tree.mc_truth_h0Zl21_E = DEFVAL;
   tree.mc_truth_h0Ztau11_E = DEFVAL;
   tree.mc_truth_h0Ztau21_E = DEFVAL;
   tree.mc_truth_h0Ztaul11_E = DEFVAL;
   tree.mc_truth_h0Ztaul21_E = DEFVAL;
   tree.mc_truth_h0Ztaunu11_E = DEFVAL;
   tree.mc_truth_h0Ztaunu21_E = DEFVAL;
   tree.mc_truth_h0Ztaunutau11_E = DEFVAL;
   tree.mc_truth_h0Ztaunutau21_E = DEFVAL;
   tree.mc_truth_h0Zq11_E = DEFVAL;
   tree.mc_truth_h0Zq21_E = DEFVAL;
   tree.mc_truth_h0Zq11_IS_E = DEFVAL;
   tree.mc_truth_h0Zq21_IS_E = DEFVAL;
   tree.mc_truth_h0Zl12_E = DEFVAL;
   tree.mc_truth_h0Zl22_E = DEFVAL;
   tree.mc_truth_h0Ztau12_E = DEFVAL;
   tree.mc_truth_h0Ztau22_E = DEFVAL;
   tree.mc_truth_h0Ztaul12_E = DEFVAL;
   tree.mc_truth_h0Ztaul22_E = DEFVAL;
   tree.mc_truth_h0Ztaunu12_E = DEFVAL;
   tree.mc_truth_h0Ztaunu22_E = DEFVAL;
   tree.mc_truth_h0Ztaunutau12_E = DEFVAL;
   tree.mc_truth_h0Ztaunutau22_E = DEFVAL;
   tree.mc_truth_h0Zq12_E = DEFVAL;
   tree.mc_truth_h0Zq22_E = DEFVAL;
   tree.mc_truth_h0Zq12_IS_E = DEFVAL;
   tree.mc_truth_h0Zq22_IS_E = DEFVAL;
   tree.mc_truth_h0Znu11_E = DEFVAL;
   tree.mc_truth_h0Znu21_E = DEFVAL;
   tree.mc_truth_h0Znu12_E = DEFVAL;
   tree.mc_truth_h0Znu22_E = DEFVAL;
   
   tree.mc_truth_h0tau1_E = DEFVAL;
   tree.mc_truth_h0tau2_E = DEFVAL;
   tree.mc_truth_h0taul1_E = DEFVAL;
   tree.mc_truth_h0taunutau1_E = DEFVAL;
   tree.mc_truth_h0taunu1_E = DEFVAL;
   tree.mc_truth_h0taul2_E = DEFVAL;
   tree.mc_truth_h0taunutau2_E = DEFVAL;
   tree.mc_truth_h0taunu2_E = DEFVAL;

   tree.mc_truth_h0b1_E = DEFVAL;
   tree.mc_truth_h0b2_E = DEFVAL;
   tree.mc_truth_h0b1_IS_E = DEFVAL;
   tree.mc_truth_h0b2_IS_E = DEFVAL;
   
   tree.mc_truth_t1_E = DEFVAL;
   tree.mc_truth_t2_E = DEFVAL;
   tree.mc_truth_tb1_E = DEFVAL;
   tree.mc_truth_tb2_E = DEFVAL;
   tree.mc_truth_tb1_IS_E = DEFVAL;
   tree.mc_truth_tb2_IS_E = DEFVAL;
   
   tree.mc_truth_tW1_E = DEFVAL;
   tree.mc_truth_tWnu1_E = DEFVAL;
   tree.mc_truth_tWnutau1_E = DEFVAL;
   tree.mc_truth_tWl1_E = DEFVAL;
   tree.mc_truth_tWtau1_E = DEFVAL;
   tree.mc_truth_tWtaunu1_E = DEFVAL;
   tree.mc_truth_tWtaunutau1_E = DEFVAL;
   tree.mc_truth_tWtaul1_E = DEFVAL;
   tree.mc_truth_tWq11_E = DEFVAL;
   tree.mc_truth_tWq21_E = DEFVAL;
   tree.mc_truth_tWq11_IS_E = DEFVAL;
   tree.mc_truth_tWq21_IS_E = DEFVAL;

   tree.mc_truth_tW2_E = DEFVAL;
   tree.mc_truth_tWnu2_E = DEFVAL;
   tree.mc_truth_tWnutau2_E = DEFVAL;
   tree.mc_truth_tWl2_E = DEFVAL;
   tree.mc_truth_tWtau2_E = DEFVAL;
   tree.mc_truth_tWtaunu2_E = DEFVAL;
   tree.mc_truth_tWtaunutau2_E = DEFVAL;
   tree.mc_truth_tWtaul2_E = DEFVAL;
   tree.mc_truth_tWq12_E = DEFVAL;
   tree.mc_truth_tWq22_E = DEFVAL;
   tree.mc_truth_tWq12_IS_E = DEFVAL;
   tree.mc_truth_tWq22_IS_E = DEFVAL;

   tree.mc_truth_j1_E = DEFVAL;
   tree.mc_truth_j2_E = DEFVAL;
   tree.mc_truth_j3_E = DEFVAL;
   
   tree.mc_truth_gammal1_E = DEFVAL;
   tree.mc_truth_gammal2_E = DEFVAL;
   tree.mc_truth_gammatau1_E = DEFVAL;
   tree.mc_truth_gammatau2_E = DEFVAL;
   tree.mc_truth_gammataul1_E = DEFVAL;
   tree.mc_truth_gammataul2_E = DEFVAL;
   tree.mc_truth_gammataunu1_E = DEFVAL;
   tree.mc_truth_gammataunu2_E = DEFVAL;
   tree.mc_truth_gammataunutau1_E = DEFVAL;
   tree.mc_truth_gammataunutau2_E = DEFVAL;
  
   // pdgId
   
   tree.mc_truth_h0_id = DEFVAL;

   tree.mc_truth_h0W1_id = DEFVAL;
   tree.mc_truth_h0W2_id = DEFVAL;
   tree.mc_truth_h0Wl1_id = DEFVAL;
   tree.mc_truth_h0Wnu1_id = DEFVAL;
   tree.mc_truth_h0Wtau1_id = DEFVAL;
   tree.mc_truth_h0Wnutau1_id = DEFVAL;
   tree.mc_truth_h0Wtaul1_id = DEFVAL;
   tree.mc_truth_h0Wtaunu1_id = DEFVAL;
   tree.mc_truth_h0Wtaunutau1_id = DEFVAL;
   tree.mc_truth_h0Wl2_id = DEFVAL;
   tree.mc_truth_h0Wnu2_id = DEFVAL;
   tree.mc_truth_h0Wtau2_id = DEFVAL;
   tree.mc_truth_h0Wnutau2_id = DEFVAL;
   tree.mc_truth_h0Wtaul2_id = DEFVAL;
   tree.mc_truth_h0Wtaunu2_id = DEFVAL;
   tree.mc_truth_h0Wtaunutau2_id = DEFVAL;
   tree.mc_truth_h0Wq11_id = DEFVAL;
   tree.mc_truth_h0Wq21_id = DEFVAL;
   tree.mc_truth_h0Wq12_id = DEFVAL;
   tree.mc_truth_h0Wq22_id = DEFVAL;
   tree.mc_truth_h0Wq11_IS_id = DEFVAL;
   tree.mc_truth_h0Wq21_IS_id = DEFVAL;
   tree.mc_truth_h0Wq12_IS_id = DEFVAL;
   tree.mc_truth_h0Wq22_IS_id = DEFVAL;
   
   tree.mc_truth_h0Z1_id = DEFVAL;
   tree.mc_truth_h0Z2_id = DEFVAL;
   tree.mc_truth_h0Zl11_id = DEFVAL;
   tree.mc_truth_h0Zl21_id = DEFVAL;
   tree.mc_truth_h0Ztau11_id = DEFVAL;
   tree.mc_truth_h0Ztau21_id = DEFVAL;
   tree.mc_truth_h0Ztaul11_id = DEFVAL;
   tree.mc_truth_h0Ztaul21_id = DEFVAL;
   tree.mc_truth_h0Ztaunu11_id = DEFVAL;
   tree.mc_truth_h0Ztaunu21_id = DEFVAL;
   tree.mc_truth_h0Ztaunutau11_id = DEFVAL;
   tree.mc_truth_h0Ztaunutau21_id = DEFVAL;
   tree.mc_truth_h0Zq11_id = DEFVAL;
   tree.mc_truth_h0Zq21_id = DEFVAL;
   tree.mc_truth_h0Zq11_IS_id = DEFVAL;
   tree.mc_truth_h0Zq21_IS_id = DEFVAL;
   tree.mc_truth_h0Zl12_id = DEFVAL;
   tree.mc_truth_h0Zl22_id = DEFVAL;
   tree.mc_truth_h0Ztau12_id = DEFVAL;
   tree.mc_truth_h0Ztau22_id = DEFVAL;
   tree.mc_truth_h0Ztaul12_id = DEFVAL;
   tree.mc_truth_h0Ztaul22_id = DEFVAL;
   tree.mc_truth_h0Ztaunu12_id = DEFVAL;
   tree.mc_truth_h0Ztaunu22_id = DEFVAL;
   tree.mc_truth_h0Ztaunutau12_id = DEFVAL;
   tree.mc_truth_h0Ztaunutau22_id = DEFVAL;
   tree.mc_truth_h0Zq12_id = DEFVAL;
   tree.mc_truth_h0Zq22_id = DEFVAL;
   tree.mc_truth_h0Zq12_IS_id = DEFVAL;
   tree.mc_truth_h0Zq22_IS_id = DEFVAL;
   tree.mc_truth_h0Znu11_id = DEFVAL;
   tree.mc_truth_h0Znu21_id = DEFVAL;
   tree.mc_truth_h0Znu12_id = DEFVAL;
   tree.mc_truth_h0Znu22_id = DEFVAL;
   
   tree.mc_truth_h0tau1_id = DEFVAL;
   tree.mc_truth_h0tau2_id = DEFVAL;
   tree.mc_truth_h0taul1_id = DEFVAL;
   tree.mc_truth_h0taunutau1_id = DEFVAL;
   tree.mc_truth_h0taunu1_id = DEFVAL;
   tree.mc_truth_h0taul2_id = DEFVAL;
   tree.mc_truth_h0taunutau2_id = DEFVAL;
   tree.mc_truth_h0taunu2_id = DEFVAL;

   tree.mc_truth_h0b1_id = DEFVAL;
   tree.mc_truth_h0b2_id = DEFVAL;
   tree.mc_truth_h0b1_IS_id = DEFVAL;
   tree.mc_truth_h0b2_IS_id = DEFVAL;
   
   tree.mc_truth_t1_id = DEFVAL;
   tree.mc_truth_t2_id = DEFVAL;
   tree.mc_truth_tb1_id = DEFVAL;
   tree.mc_truth_tb2_id = DEFVAL;
   tree.mc_truth_tb1_IS_id = DEFVAL;
   tree.mc_truth_tb2_IS_id = DEFVAL;
   
   tree.mc_truth_tW1_id = DEFVAL;
   tree.mc_truth_tWnu1_id = DEFVAL;
   tree.mc_truth_tWnutau1_id = DEFVAL;
   tree.mc_truth_tWl1_id = DEFVAL;
   tree.mc_truth_tWtau1_id = DEFVAL;
   tree.mc_truth_tWtaunu1_id = DEFVAL;
   tree.mc_truth_tWtaunutau1_id = DEFVAL;
   tree.mc_truth_tWtaul1_id = DEFVAL;
   tree.mc_truth_tWq11_id = DEFVAL;
   tree.mc_truth_tWq21_id = DEFVAL;
   tree.mc_truth_tWq11_IS_id = DEFVAL;
   tree.mc_truth_tWq21_IS_id = DEFVAL;

   tree.mc_truth_tW2_id = DEFVAL;
   tree.mc_truth_tWnu2_id = DEFVAL;
   tree.mc_truth_tWnutau2_id = DEFVAL;
   tree.mc_truth_tWl2_id = DEFVAL;
   tree.mc_truth_tWtau2_id = DEFVAL;
   tree.mc_truth_tWtaunu2_id = DEFVAL;
   tree.mc_truth_tWtaunutau2_id = DEFVAL;
   tree.mc_truth_tWtaul2_id = DEFVAL;
   tree.mc_truth_tWq12_id = DEFVAL;
   tree.mc_truth_tWq22_id = DEFVAL;
   tree.mc_truth_tWq12_IS_id = DEFVAL;
   tree.mc_truth_tWq22_IS_id = DEFVAL;

   tree.mc_truth_j1_id = DEFVAL;
   tree.mc_truth_j2_id = DEFVAL;
   tree.mc_truth_j3_id = DEFVAL;
   
   tree.mc_truth_gammal1_id = DEFVAL;
   tree.mc_truth_gammal2_id = DEFVAL;
   tree.mc_truth_gammatau1_id = DEFVAL;
   tree.mc_truth_gammatau2_id = DEFVAL;
   tree.mc_truth_gammataul1_id = DEFVAL;
   tree.mc_truth_gammataul2_id = DEFVAL;
   tree.mc_truth_gammataunu1_id = DEFVAL;
   tree.mc_truth_gammataunu2_id = DEFVAL;
   tree.mc_truth_gammataunutau1_id = DEFVAL;
   tree.mc_truth_gammataunutau2_id = DEFVAL;
  
   // status

   tree.mc_truth_h0_status = DEFVAL;

   tree.mc_truth_h0W1_status = DEFVAL;
   tree.mc_truth_h0W2_status = DEFVAL;
   tree.mc_truth_h0Wl1_status = DEFVAL;
   tree.mc_truth_h0Wnu1_status = DEFVAL;
   tree.mc_truth_h0Wtau1_status = DEFVAL;
   tree.mc_truth_h0Wnutau1_status = DEFVAL;
   tree.mc_truth_h0Wtaul1_status = DEFVAL;
   tree.mc_truth_h0Wtaunu1_status = DEFVAL;
   tree.mc_truth_h0Wtaunutau1_status = DEFVAL;
   tree.mc_truth_h0Wl2_status = DEFVAL;
   tree.mc_truth_h0Wnu2_status = DEFVAL;
   tree.mc_truth_h0Wtau2_status = DEFVAL;
   tree.mc_truth_h0Wnutau2_status = DEFVAL;
   tree.mc_truth_h0Wtaul2_status = DEFVAL;
   tree.mc_truth_h0Wtaunu2_status = DEFVAL;
   tree.mc_truth_h0Wtaunutau2_status = DEFVAL;
   tree.mc_truth_h0Wq11_status = DEFVAL;
   tree.mc_truth_h0Wq21_status = DEFVAL;
   tree.mc_truth_h0Wq12_status = DEFVAL;
   tree.mc_truth_h0Wq22_status = DEFVAL;
   tree.mc_truth_h0Wq11_IS_status = DEFVAL;
   tree.mc_truth_h0Wq21_IS_status = DEFVAL;
   tree.mc_truth_h0Wq12_IS_status = DEFVAL;
   tree.mc_truth_h0Wq22_IS_status = DEFVAL;
   
   tree.mc_truth_h0Z1_status = DEFVAL;
   tree.mc_truth_h0Z2_status = DEFVAL;
   tree.mc_truth_h0Zl11_status = DEFVAL;
   tree.mc_truth_h0Zl21_status = DEFVAL;
   tree.mc_truth_h0Ztau11_status = DEFVAL;
   tree.mc_truth_h0Ztau21_status = DEFVAL;
   tree.mc_truth_h0Ztaul11_status = DEFVAL;
   tree.mc_truth_h0Ztaul21_status = DEFVAL;
   tree.mc_truth_h0Ztaunu11_status = DEFVAL;
   tree.mc_truth_h0Ztaunu21_status = DEFVAL;
   tree.mc_truth_h0Ztaunutau11_status = DEFVAL;
   tree.mc_truth_h0Ztaunutau21_status = DEFVAL;
   tree.mc_truth_h0Zq11_status = DEFVAL;
   tree.mc_truth_h0Zq21_status = DEFVAL;
   tree.mc_truth_h0Zq11_IS_status = DEFVAL;
   tree.mc_truth_h0Zq21_IS_status = DEFVAL;
   tree.mc_truth_h0Zl12_status = DEFVAL;
   tree.mc_truth_h0Zl22_status = DEFVAL;
   tree.mc_truth_h0Ztau12_status = DEFVAL;
   tree.mc_truth_h0Ztau22_status = DEFVAL;
   tree.mc_truth_h0Ztaul12_status = DEFVAL;
   tree.mc_truth_h0Ztaul22_status = DEFVAL;
   tree.mc_truth_h0Ztaunu12_status = DEFVAL;
   tree.mc_truth_h0Ztaunu22_status = DEFVAL;
   tree.mc_truth_h0Ztaunutau12_status = DEFVAL;
   tree.mc_truth_h0Ztaunutau22_status = DEFVAL;
   tree.mc_truth_h0Zq12_status = DEFVAL;
   tree.mc_truth_h0Zq22_status = DEFVAL;
   tree.mc_truth_h0Zq12_IS_status = DEFVAL;
   tree.mc_truth_h0Zq22_IS_status = DEFVAL;
   tree.mc_truth_h0Znu11_status = DEFVAL;
   tree.mc_truth_h0Znu21_status = DEFVAL;
   tree.mc_truth_h0Znu12_status = DEFVAL;
   tree.mc_truth_h0Znu22_status = DEFVAL;
   
   tree.mc_truth_h0tau1_status = DEFVAL;
   tree.mc_truth_h0tau2_status = DEFVAL;
   tree.mc_truth_h0taul1_status = DEFVAL;
   tree.mc_truth_h0taunutau1_status = DEFVAL;
   tree.mc_truth_h0taunu1_status = DEFVAL;
   tree.mc_truth_h0taul2_status = DEFVAL;
   tree.mc_truth_h0taunutau2_status = DEFVAL;
   tree.mc_truth_h0taunu2_status = DEFVAL;

   tree.mc_truth_h0b1_status = DEFVAL;
   tree.mc_truth_h0b2_status = DEFVAL;
   tree.mc_truth_h0b1_IS_status = DEFVAL;
   tree.mc_truth_h0b2_IS_status = DEFVAL;
   
   tree.mc_truth_t1_status = DEFVAL;
   tree.mc_truth_t2_status = DEFVAL;
   tree.mc_truth_tb1_status = DEFVAL;
   tree.mc_truth_tb2_status = DEFVAL;
   tree.mc_truth_tb1_IS_status = DEFVAL;
   tree.mc_truth_tb2_IS_status = DEFVAL;

   tree.mc_truth_tW1_status = DEFVAL;
   tree.mc_truth_tWnu1_status = DEFVAL;
   tree.mc_truth_tWnutau1_status = DEFVAL;
   tree.mc_truth_tWl1_status = DEFVAL;
   tree.mc_truth_tWtau1_status = DEFVAL;
   tree.mc_truth_tWtaunu1_status = DEFVAL;
   tree.mc_truth_tWtaunutau1_status = DEFVAL;
   tree.mc_truth_tWtaul1_status = DEFVAL;
   tree.mc_truth_tWq11_status = DEFVAL;
   tree.mc_truth_tWq21_status = DEFVAL;
   tree.mc_truth_tWq11_IS_status = DEFVAL;
   tree.mc_truth_tWq21_IS_status = DEFVAL;

   tree.mc_truth_tW2_status = DEFVAL;
   tree.mc_truth_tWnu2_status = DEFVAL;
   tree.mc_truth_tWnutau2_status = DEFVAL;
   tree.mc_truth_tWl2_status = DEFVAL;
   tree.mc_truth_tWtau2_status = DEFVAL;
   tree.mc_truth_tWtaunu2_status = DEFVAL;
   tree.mc_truth_tWtaunutau2_status = DEFVAL;
   tree.mc_truth_tWtaul2_status = DEFVAL;
   tree.mc_truth_tWq12_status = DEFVAL;
   tree.mc_truth_tWq22_status = DEFVAL;
   tree.mc_truth_tWq12_IS_status = DEFVAL;
   tree.mc_truth_tWq22_IS_status = DEFVAL;

   tree.mc_truth_j1_status = DEFVAL;
   tree.mc_truth_j2_status = DEFVAL;
   tree.mc_truth_j3_status = DEFVAL;
   
   tree.mc_truth_gammal1_status = DEFVAL;
   tree.mc_truth_gammal2_status = DEFVAL;
   tree.mc_truth_gammatau1_status = DEFVAL;
   tree.mc_truth_gammatau2_status = DEFVAL;
   tree.mc_truth_gammataul1_status = DEFVAL;
   tree.mc_truth_gammataul2_status = DEFVAL;
   tree.mc_truth_gammataunu1_status = DEFVAL;
   tree.mc_truth_gammataunu2_status = DEFVAL;
   tree.mc_truth_gammataunutau1_status = DEFVAL;
   tree.mc_truth_gammataunutau2_status = DEFVAL;
  

   // tZq
   tree.mc_truth_tzq_channel = DEFVAL;

   // tHq
   tree.mc_truth_thq_channel = DEFVAL;
   
   // TLV

   tree.mc_truth_Z_p4.Clear();
   tree.mc_truth_Zl1_p4.Clear();
   tree.mc_truth_Zl2_p4.Clear();
   tree.mc_truth_Ztau1_p4.Clear();
   tree.mc_truth_Ztau2_p4.Clear();
   tree.mc_truth_Ztaul1_p4.Clear();
   tree.mc_truth_Ztaul2_p4.Clear();
   tree.mc_truth_Ztaunu1_p4.Clear();
   tree.mc_truth_Ztaunu2_p4.Clear();
   tree.mc_truth_Ztaunutau1_p4.Clear();
   tree.mc_truth_Ztaunutau2_p4.Clear();
   tree.mc_truth_Zq1_p4.Clear();
   tree.mc_truth_Zq2_p4.Clear();
   tree.mc_truth_Zq1_IS_p4.Clear();
   tree.mc_truth_Zq2_IS_p4.Clear();

   tree.mc_truth_W_p4.Clear();
   tree.mc_truth_Wl_p4.Clear();
   tree.mc_truth_Wnu_p4.Clear();
   tree.mc_truth_Wtau_p4.Clear();
   tree.mc_truth_Wtaunu_p4.Clear();
   tree.mc_truth_Wtaunutau_p4.Clear();
   tree.mc_truth_Wtaul_p4.Clear();
   tree.mc_truth_Wnutau_p4.Clear();
   tree.mc_truth_Wq1_p4.Clear();
   tree.mc_truth_Wq2_p4.Clear();
   tree.mc_truth_Wq1_IS_p4.Clear();
   tree.mc_truth_Wq2_IS_p4.Clear();
   
   tree.mc_truth_t_p4.Clear();
   tree.mc_truth_tb_p4.Clear();
   tree.mc_truth_tb_IS_p4.Clear();
   tree.mc_truth_tW_p4.Clear();
   tree.mc_truth_tWnu_p4.Clear();
   tree.mc_truth_tWnutau_p4.Clear();
   tree.mc_truth_tWl_p4.Clear();
   tree.mc_truth_tWtau_p4.Clear();
   tree.mc_truth_tWtaunu_p4.Clear();
   tree.mc_truth_tWtaunutau_p4.Clear();
   tree.mc_truth_tWtaul_p4.Clear();
   tree.mc_truth_tWq1_p4.Clear();
   tree.mc_truth_tWq2_p4.Clear();
   tree.mc_truth_tWq1_IS_p4.Clear();
   tree.mc_truth_tWq2_IS_p4.Clear();

   // pt

   tree.mc_truth_Z_pt = DEFVAL;
   tree.mc_truth_Zl1_pt = DEFVAL;
   tree.mc_truth_Zl2_pt = DEFVAL;
   tree.mc_truth_Ztau1_pt = DEFVAL;
   tree.mc_truth_Ztau2_pt = DEFVAL;
   tree.mc_truth_Ztaul1_pt = DEFVAL;
   tree.mc_truth_Ztaul2_pt = DEFVAL;
   tree.mc_truth_Ztaunu1_pt = DEFVAL;
   tree.mc_truth_Ztaunu2_pt = DEFVAL;
   tree.mc_truth_Ztaunutau1_pt = DEFVAL;
   tree.mc_truth_Ztaunutau2_pt = DEFVAL;
   tree.mc_truth_Zq1_pt = DEFVAL;
   tree.mc_truth_Zq2_pt = DEFVAL;
   tree.mc_truth_Zq1_IS_pt = DEFVAL;
   tree.mc_truth_Zq2_IS_pt = DEFVAL;

   tree.mc_truth_W_pt = DEFVAL;
   tree.mc_truth_Wl_pt = DEFVAL;
   tree.mc_truth_Wnu_pt = DEFVAL;
   tree.mc_truth_Wtau_pt = DEFVAL;
   tree.mc_truth_Wtaunu_pt = DEFVAL;
   tree.mc_truth_Wtaunutau_pt = DEFVAL;
   tree.mc_truth_Wtaul_pt = DEFVAL;
   tree.mc_truth_Wnutau_pt = DEFVAL;
   tree.mc_truth_Wq1_pt = DEFVAL;
   tree.mc_truth_Wq2_pt = DEFVAL;
   tree.mc_truth_Wq1_IS_pt = DEFVAL;
   tree.mc_truth_Wq2_IS_pt = DEFVAL;
   
   tree.mc_truth_t_pt = DEFVAL;
   tree.mc_truth_tb_pt = DEFVAL;
   tree.mc_truth_tb_IS_pt = DEFVAL;
   tree.mc_truth_tW_pt = DEFVAL;
   tree.mc_truth_tWnu_pt = DEFVAL;
   tree.mc_truth_tWnutau_pt = DEFVAL;
   tree.mc_truth_tWl_pt = DEFVAL;
   tree.mc_truth_tWtau_pt = DEFVAL;
   tree.mc_truth_tWtaunu_pt = DEFVAL;
   tree.mc_truth_tWtaunutau_pt = DEFVAL;
   tree.mc_truth_tWtaul_pt = DEFVAL;
   tree.mc_truth_tWq1_pt = DEFVAL;
   tree.mc_truth_tWq2_pt = DEFVAL;
   tree.mc_truth_tWq1_IS_pt = DEFVAL;
   tree.mc_truth_tWq2_IS_pt = DEFVAL;
   
   // eta

   tree.mc_truth_Z_eta = DEFVAL;
   tree.mc_truth_Zl1_eta = DEFVAL;
   tree.mc_truth_Zl2_eta = DEFVAL;
   tree.mc_truth_Ztau1_eta = DEFVAL;
   tree.mc_truth_Ztau2_eta = DEFVAL;
   tree.mc_truth_Ztaul1_eta = DEFVAL;
   tree.mc_truth_Ztaul2_eta = DEFVAL;
   tree.mc_truth_Ztaunu1_eta = DEFVAL;
   tree.mc_truth_Ztaunu2_eta = DEFVAL;
   tree.mc_truth_Ztaunutau1_eta = DEFVAL;
   tree.mc_truth_Ztaunutau2_eta = DEFVAL;
   tree.mc_truth_Zq1_eta = DEFVAL;
   tree.mc_truth_Zq2_eta = DEFVAL;
   tree.mc_truth_Zq1_IS_eta = DEFVAL;
   tree.mc_truth_Zq2_IS_eta = DEFVAL;

   tree.mc_truth_W_eta = DEFVAL;
   tree.mc_truth_Wl_eta = DEFVAL;
   tree.mc_truth_Wnu_eta = DEFVAL;
   tree.mc_truth_Wtau_eta = DEFVAL;
   tree.mc_truth_Wtaunu_eta = DEFVAL;
   tree.mc_truth_Wtaunutau_eta = DEFVAL;
   tree.mc_truth_Wtaul_eta = DEFVAL;
   tree.mc_truth_Wnutau_eta = DEFVAL;
   tree.mc_truth_Wq1_eta = DEFVAL;
   tree.mc_truth_Wq2_eta = DEFVAL;
   tree.mc_truth_Wq1_IS_eta = DEFVAL;
   tree.mc_truth_Wq2_IS_eta = DEFVAL;
   
   tree.mc_truth_t_eta = DEFVAL;
   tree.mc_truth_tb_eta = DEFVAL;
   tree.mc_truth_tb_IS_eta = DEFVAL;
   tree.mc_truth_tW_eta = DEFVAL;
   tree.mc_truth_tWnu_eta = DEFVAL;
   tree.mc_truth_tWnutau_eta = DEFVAL;
   tree.mc_truth_tWl_eta = DEFVAL;
   tree.mc_truth_tWtau_eta = DEFVAL;
   tree.mc_truth_tWtaunu_eta = DEFVAL;
   tree.mc_truth_tWtaunutau_eta = DEFVAL;
   tree.mc_truth_tWtaul_eta = DEFVAL;
   tree.mc_truth_tWq1_eta = DEFVAL;
   tree.mc_truth_tWq2_eta = DEFVAL;
   tree.mc_truth_tWq1_IS_eta = DEFVAL;
   tree.mc_truth_tWq2_IS_eta = DEFVAL;
   
   // phi

   tree.mc_truth_Z_phi = DEFVAL;
   tree.mc_truth_Zl1_phi = DEFVAL;
   tree.mc_truth_Zl2_phi = DEFVAL;
   tree.mc_truth_Ztau1_phi = DEFVAL;
   tree.mc_truth_Ztau2_phi = DEFVAL;
   tree.mc_truth_Ztaul1_phi = DEFVAL;
   tree.mc_truth_Ztaul2_phi = DEFVAL;
   tree.mc_truth_Ztaunu1_phi = DEFVAL;
   tree.mc_truth_Ztaunu2_phi = DEFVAL;
   tree.mc_truth_Ztaunutau1_phi = DEFVAL;
   tree.mc_truth_Ztaunutau2_phi = DEFVAL;
   tree.mc_truth_Zq1_phi = DEFVAL;
   tree.mc_truth_Zq2_phi = DEFVAL;
   tree.mc_truth_Zq1_IS_phi = DEFVAL;
   tree.mc_truth_Zq2_IS_phi = DEFVAL;

   tree.mc_truth_W_phi = DEFVAL;
   tree.mc_truth_Wl_phi = DEFVAL;
   tree.mc_truth_Wnu_phi = DEFVAL;
   tree.mc_truth_Wtau_phi = DEFVAL;
   tree.mc_truth_Wtaunu_phi = DEFVAL;
   tree.mc_truth_Wtaunutau_phi = DEFVAL;
   tree.mc_truth_Wtaul_phi = DEFVAL;
   tree.mc_truth_Wnutau_phi = DEFVAL;
   tree.mc_truth_Wq1_phi = DEFVAL;
   tree.mc_truth_Wq2_phi = DEFVAL;
   tree.mc_truth_Wq1_IS_phi = DEFVAL;
   tree.mc_truth_Wq2_IS_phi = DEFVAL;
   
   tree.mc_truth_t_phi = DEFVAL;
   tree.mc_truth_tb_phi = DEFVAL;
   tree.mc_truth_tb_IS_phi = DEFVAL;
   tree.mc_truth_tW_phi = DEFVAL;
   tree.mc_truth_tWnu_phi = DEFVAL;
   tree.mc_truth_tWnutau_phi = DEFVAL;
   tree.mc_truth_tWl_phi = DEFVAL;
   tree.mc_truth_tWtau_phi = DEFVAL;
   tree.mc_truth_tWtaunu_phi = DEFVAL;
   tree.mc_truth_tWtaunutau_phi = DEFVAL;
   tree.mc_truth_tWtaul_phi = DEFVAL;
   tree.mc_truth_tWq1_phi = DEFVAL;
   tree.mc_truth_tWq2_phi = DEFVAL;
   tree.mc_truth_tWq1_IS_phi = DEFVAL;
   tree.mc_truth_tWq2_IS_phi = DEFVAL;
   
   // E

   tree.mc_truth_Z_E = DEFVAL;
   tree.mc_truth_Zl1_E = DEFVAL;
   tree.mc_truth_Zl2_E = DEFVAL;
   tree.mc_truth_Ztau1_E = DEFVAL;
   tree.mc_truth_Ztau2_E = DEFVAL;
   tree.mc_truth_Ztaul1_E = DEFVAL;
   tree.mc_truth_Ztaul2_E = DEFVAL;
   tree.mc_truth_Ztaunu1_E = DEFVAL;
   tree.mc_truth_Ztaunu2_E = DEFVAL;
   tree.mc_truth_Ztaunutau1_E = DEFVAL;
   tree.mc_truth_Ztaunutau2_E = DEFVAL;
   tree.mc_truth_Zq1_E = DEFVAL;
   tree.mc_truth_Zq2_E = DEFVAL;
   tree.mc_truth_Zq1_IS_E = DEFVAL;
   tree.mc_truth_Zq2_IS_E = DEFVAL;

   tree.mc_truth_W_E = DEFVAL;
   tree.mc_truth_Wl_E = DEFVAL;
   tree.mc_truth_Wnu_E = DEFVAL;
   tree.mc_truth_Wtau_E = DEFVAL;
   tree.mc_truth_Wtaunu_E = DEFVAL;
   tree.mc_truth_Wtaunutau_E = DEFVAL;
   tree.mc_truth_Wtaul_E = DEFVAL;
   tree.mc_truth_Wnutau_E = DEFVAL;
   tree.mc_truth_Wq1_E = DEFVAL;
   tree.mc_truth_Wq2_E = DEFVAL;
   tree.mc_truth_Wq1_IS_E = DEFVAL;
   tree.mc_truth_Wq2_IS_E = DEFVAL;
   
   tree.mc_truth_t_E = DEFVAL;
   tree.mc_truth_tb_E = DEFVAL;
   tree.mc_truth_tb_IS_E = DEFVAL;
   tree.mc_truth_tW_E = DEFVAL;
   tree.mc_truth_tWnu_E = DEFVAL;
   tree.mc_truth_tWnutau_E = DEFVAL;
   tree.mc_truth_tWl_E = DEFVAL;
   tree.mc_truth_tWtau_E = DEFVAL;
   tree.mc_truth_tWtaunu_E = DEFVAL;
   tree.mc_truth_tWtaunutau_E = DEFVAL;
   tree.mc_truth_tWtaul_E = DEFVAL;
   tree.mc_truth_tWq1_E = DEFVAL;
   tree.mc_truth_tWq2_E = DEFVAL;
   tree.mc_truth_tWq1_IS_E = DEFVAL;
   tree.mc_truth_tWq2_IS_E = DEFVAL;
   
   // pdgId

   tree.mc_truth_Z_id = DEFVAL;
   tree.mc_truth_Zl1_id = DEFVAL;
   tree.mc_truth_Zl2_id = DEFVAL;
   tree.mc_truth_Ztau1_id = DEFVAL;
   tree.mc_truth_Ztau2_id = DEFVAL;
   tree.mc_truth_Ztaul1_id = DEFVAL;
   tree.mc_truth_Ztaul2_id = DEFVAL;
   tree.mc_truth_Ztaunu1_id = DEFVAL;
   tree.mc_truth_Ztaunu2_id = DEFVAL;
   tree.mc_truth_Ztaunutau1_id = DEFVAL;
   tree.mc_truth_Ztaunutau2_id = DEFVAL;
   tree.mc_truth_Zq1_id = DEFVAL;
   tree.mc_truth_Zq2_id = DEFVAL;
   tree.mc_truth_Zq1_IS_id = DEFVAL;
   tree.mc_truth_Zq2_IS_id = DEFVAL;

   tree.mc_truth_W_id = DEFVAL;
   tree.mc_truth_Wl_id = DEFVAL;
   tree.mc_truth_Wnu_id = DEFVAL;
   tree.mc_truth_Wtau_id = DEFVAL;
   tree.mc_truth_Wtaunu_id = DEFVAL;
   tree.mc_truth_Wtaunutau_id = DEFVAL;
   tree.mc_truth_Wtaul_id = DEFVAL;
   tree.mc_truth_Wnutau_id = DEFVAL;
   tree.mc_truth_Wq1_id = DEFVAL;
   tree.mc_truth_Wq2_id = DEFVAL;
   tree.mc_truth_Wq1_IS_id = DEFVAL;
   tree.mc_truth_Wq2_IS_id = DEFVAL;
   
   tree.mc_truth_t_id = DEFVAL;
   tree.mc_truth_tb_id = DEFVAL;
   tree.mc_truth_tb_IS_id = DEFVAL;
   tree.mc_truth_tW_id = DEFVAL;
   tree.mc_truth_tWnu_id = DEFVAL;
   tree.mc_truth_tWnutau_id = DEFVAL;
   tree.mc_truth_tWl_id = DEFVAL;
   tree.mc_truth_tWtau_id = DEFVAL;
   tree.mc_truth_tWtaunu_id = DEFVAL;
   tree.mc_truth_tWtaunutau_id = DEFVAL;
   tree.mc_truth_tWtaul_id = DEFVAL;
   tree.mc_truth_tWq1_id = DEFVAL;
   tree.mc_truth_tWq2_id = DEFVAL;
   tree.mc_truth_tWq1_IS_id = DEFVAL;
   tree.mc_truth_tWq2_IS_id = DEFVAL;

   // status
   
   tree.mc_truth_Z_status = DEFVAL;
   tree.mc_truth_Zl1_status = DEFVAL;
   tree.mc_truth_Zl2_status = DEFVAL;
   tree.mc_truth_Ztau1_status = DEFVAL;
   tree.mc_truth_Ztau2_status = DEFVAL;
   tree.mc_truth_Ztaul1_status = DEFVAL;
   tree.mc_truth_Ztaul2_status = DEFVAL;
   tree.mc_truth_Ztaunu1_status = DEFVAL;
   tree.mc_truth_Ztaunu2_status = DEFVAL;
   tree.mc_truth_Ztaunutau1_status = DEFVAL;
   tree.mc_truth_Ztaunutau2_status = DEFVAL;
   tree.mc_truth_Zq1_status = DEFVAL;
   tree.mc_truth_Zq2_status = DEFVAL;
   tree.mc_truth_Zq1_IS_status = DEFVAL;
   tree.mc_truth_Zq2_IS_status = DEFVAL;

   tree.mc_truth_W_status = DEFVAL;
   tree.mc_truth_Wl_status = DEFVAL;
   tree.mc_truth_Wnu_status = DEFVAL;
   tree.mc_truth_Wtau_status = DEFVAL;
   tree.mc_truth_Wtaunu_status = DEFVAL;
   tree.mc_truth_Wtaunutau_status = DEFVAL;
   tree.mc_truth_Wtaul_status = DEFVAL;
   tree.mc_truth_Wnutau_status = DEFVAL;
   tree.mc_truth_Wq1_status = DEFVAL;
   tree.mc_truth_Wq2_status = DEFVAL;
   tree.mc_truth_Wq1_IS_status = DEFVAL;
   tree.mc_truth_Wq2_IS_status = DEFVAL;
   
   tree.mc_truth_t_status = DEFVAL;
   tree.mc_truth_tb_status = DEFVAL;
   tree.mc_truth_tb_IS_status = DEFVAL;
   tree.mc_truth_tW_status = DEFVAL;
   tree.mc_truth_tWnu_status = DEFVAL;
   tree.mc_truth_tWnutau_status = DEFVAL;
   tree.mc_truth_tWl_status = DEFVAL;
   tree.mc_truth_tWtau_status = DEFVAL;
   tree.mc_truth_tWtaunu_status = DEFVAL;
   tree.mc_truth_tWtaunutau_status = DEFVAL;
   tree.mc_truth_tWtaul_status = DEFVAL;
   tree.mc_truth_tWq1_status = DEFVAL;
   tree.mc_truth_tWq2_status = DEFVAL;
   tree.mc_truth_tWq1_IS_status = DEFVAL;
   tree.mc_truth_tWq2_IS_status = DEFVAL;
   
   // gen
   tree.gen_pt.clear();
   tree.gen_eta.clear();
   tree.gen_phi.clear();
   tree.gen_m.clear();
   tree.gen_status.clear();
   tree.gen_id.clear();
   tree.gen_charge.clear();
   tree.gen_index.clear();
   tree.gen_mother_index.clear();
   tree.gen_daughter_n.clear();
   tree.gen_daughter_index.clear();
}

void MCTruth::fillTTHSignalGenParticles(const edm::Event& iEvent,
					const edm::EventSetup& iSetup,
					FlatTree& tree,
					const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   reco::GenParticle *h0 = 0;
   
   reco::GenParticle *h0W1 = 0;
   reco::GenParticle *h0W2 = 0;
   reco::GenParticle *h0Wl1 = 0;
   reco::GenParticle *h0Wnu1 = 0;
   reco::GenParticle *h0Wtau1 = 0;
   reco::GenParticle *h0Wnutau1 = 0;
   reco::GenParticle *h0Wtaul1 = 0;
   reco::GenParticle *h0Wtaunu1 = 0;
   reco::GenParticle *h0Wtaunutau1 = 0;
   reco::GenParticle *h0Wl2 = 0;
   reco::GenParticle *h0Wnu2 = 0;
   reco::GenParticle *h0Wtau2 = 0;
   reco::GenParticle *h0Wnutau2 = 0;
   reco::GenParticle *h0Wtaul2 = 0;
   reco::GenParticle *h0Wtaunu2 = 0;
   reco::GenParticle *h0Wtaunutau2 = 0;
   reco::GenParticle *h0Wq11 = 0;
   reco::GenParticle *h0Wq21 = 0;
   reco::GenParticle *h0Wq12 = 0;
   reco::GenParticle *h0Wq22 = 0;
   reco::GenParticle *h0Wq11_IS = 0;
   reco::GenParticle *h0Wq21_IS = 0;
   reco::GenParticle *h0Wq12_IS = 0;
   reco::GenParticle *h0Wq22_IS = 0;

   reco::GenParticle *h0Z1 = 0;
   reco::GenParticle *h0Z2 = 0;
   reco::GenParticle *h0Zl11 = 0;
   reco::GenParticle *h0Zl21 = 0;
   reco::GenParticle *h0Ztau11 = 0;
   reco::GenParticle *h0Ztau21 = 0;
   reco::GenParticle *h0Ztaul11 = 0;
   reco::GenParticle *h0Ztaul21 = 0;
   reco::GenParticle *h0Ztaunu11 = 0;
   reco::GenParticle *h0Ztaunu21 = 0;
   reco::GenParticle *h0Ztaunutau11 = 0;
   reco::GenParticle *h0Ztaunutau21 = 0;
   reco::GenParticle *h0Zq11 = 0;
   reco::GenParticle *h0Zq21 = 0;
   reco::GenParticle *h0Zq11_IS = 0;
   reco::GenParticle *h0Zq21_IS = 0;
   reco::GenParticle *h0Zl12 = 0;
   reco::GenParticle *h0Zl22 = 0;
   reco::GenParticle *h0Ztau12 = 0;
   reco::GenParticle *h0Ztau22 = 0;
   reco::GenParticle *h0Ztaul12 = 0;
   reco::GenParticle *h0Ztaul22 = 0;
   reco::GenParticle *h0Ztaunu12 = 0;
   reco::GenParticle *h0Ztaunu22 = 0;
   reco::GenParticle *h0Ztaunutau12 = 0;
   reco::GenParticle *h0Ztaunutau22 = 0;
   reco::GenParticle *h0Zq12 = 0;
   reco::GenParticle *h0Zq22 = 0;
   reco::GenParticle *h0Zq12_IS = 0;
   reco::GenParticle *h0Zq22_IS = 0;
   reco::GenParticle *h0Znu11 = 0;
   reco::GenParticle *h0Znu21 = 0;
   reco::GenParticle *h0Znu12 = 0;
   reco::GenParticle *h0Znu22 = 0;
   
   reco::GenParticle *h0tau1 = 0;
   reco::GenParticle *h0tau2 = 0;
   reco::GenParticle *h0taul1 = 0;
   reco::GenParticle *h0taunutau1 = 0;
   reco::GenParticle *h0taunu1 = 0;
   reco::GenParticle *h0taul2 = 0;
   reco::GenParticle *h0taunutau2 = 0;
   reco::GenParticle *h0taunu2 = 0;
   
   reco::GenParticle *t1 = 0;
   reco::GenParticle *t2 = 0;   

   reco::GenParticle *tb1 = 0;
   reco::GenParticle *tb2 = 0;
   reco::GenParticle *tb1_IS = 0;
   reco::GenParticle *tb2_IS = 0;
   
   reco::GenParticle *tW1 = 0;
   reco::GenParticle *tWnu1 = 0;
   reco::GenParticle *tWnutau1 = 0;
   reco::GenParticle *tWl1 = 0;
   reco::GenParticle *tWtau1 = 0;
   reco::GenParticle *tWtaunu1 = 0;
   reco::GenParticle *tWtaunutau1 = 0;
   reco::GenParticle *tWtaul1 = 0;
   reco::GenParticle *tWq11 = 0;
   reco::GenParticle *tWq21 = 0;
   reco::GenParticle *tWq11_IS = 0;
   reco::GenParticle *tWq21_IS = 0;
   
   reco::GenParticle *tW2 = 0;
   reco::GenParticle *tWnu2 = 0;
   reco::GenParticle *tWnutau2 = 0;
   reco::GenParticle *tWl2 = 0;
   reco::GenParticle *tWtau2 = 0;
   reco::GenParticle *tWtaunu2 = 0;
   reco::GenParticle *tWtaunutau2 = 0;
   reco::GenParticle *tWtaul2 = 0;
   reco::GenParticle *tWq12 = 0;
   reco::GenParticle *tWq22 = 0;
   reco::GenParticle *tWq12_IS = 0;
   reco::GenParticle *tWq22_IS = 0;

   reco::GenParticle *j1 = 0;
   reco::GenParticle *j2 = 0;
   reco::GenParticle *j3 = 0;
   
   int chan = -666;

   // 0   = (t->bW,W->lnu)(t->bW,W->lnu)
   // 1   = (t->bW,W->qq)(t->bW,W->qq)
   // 2   = (t->bW,W->qq)(t->bW,W->lnu)
   // 3   = (t->bW,W->lnu)(t->bW,W->qq)
   // 4   = (t->bW,W->tauLnu)(t->bW,W->tauLnu)
   // 5   = (t->bW,W->qq)(t->bW,W->tauLnu)
   // 6   = (t->bW,W->tauLnu)(t->bW,W->qq)
   // 7   = (t->bW,W->tauHnu)(t->bW,W->tauHnu)
   // 8   = (t->bW,W->qq)(t->bW,W->tauHnu)
   // 9   = (t->bW,W->tauHnu)(t->bW,W->qq)
   // 10  = (t->bW,W->tauHnu)(t->bW,W->tauLnu)
   // 11  = (t->bW,W->tauLnu)(t->bW,W->tauHnu)
   // 12  = (t->bW,W->lnu)(t->bW,W->tauLnu)
   // 13  = (t->bW,W->lnu)(t->bW,W->tauHnu)
   // 14  = (t->bW,W->tauLnu)(t->bW,W->lnu)
   // 15  = (t->bW,W->tauHnu)(t->bW,W->lnu)
   
   // (H->WW:W->lnu,W->lnu)                 +0
   // (H->WW:W->tauLtaunu,W->tauLtaunu)     +20
   // (H->WW:W->tauHtaunu,W->tauHtaunu)     +40
   // (H->WW:W->tauHtaunu,W->tauLtaunu)     +60
   // (H->WW:W->tauLtaunu,W->tauHtaunu)     +80
   // (H->WW:W->qq,W->qq)                   +100
   // (H->WW:W->lnu,W->qq)                  +120
   // (H->WW:W->tauLnutau,W->qq)            +140
   // (H->WW:W->tauHnutau,W->qq)            +160
   // (H->WW:W->qq,W->lnu)                  +180
   // (H->WW:W->qq,W->tauLtaunu)            +200
   // (H->WW:W->qq,W->tauHtaunu)            +220
   // (H->WW:W->lnu,W->tauLtaunu)           +240
   // (H->WW:W->lnu,W->tauHtaunu)           +260
   // (H->WW:W->tauLtaunu,W->lnu)           +280
   // (H->WW:W->tauHtaunu,W->lnu)           +300
   
   // +1000
   // (H->ZZ:Z->ll,Z->ll)                   +0
   // (H->ZZ:Z->tauLtauL,Z->tauLtauL)       +20
   // (H->ZZ:Z->tauLtauL,Z->tauHtauH)       +40
   // (H->ZZ:Z->tauLtauL,Z->tauLtauH)       +60
   // (H->ZZ:Z->tauLtauL,Z->tauHtauL)       +80
   // (H->ZZ:Z->tauHtauH,Z->tauLtauL)       +100
   // (H->ZZ:Z->tauHtauH,Z->tauHtauH)       +120
   // (H->ZZ:Z->tauHtauH,Z->tauLtauH)       +140
   // (H->ZZ:Z->tauHtauH,Z->tauHtauL)       +160
   // (H->ZZ:Z->tauLtauH,Z->tauLtauL)       +180
   // (H->ZZ:Z->tauLtauH,Z->tauHtauH)       +200
   // (H->ZZ:Z->tauLtauH,Z->tauLtauH)       +220
   // (H->ZZ:Z->tauLtauH,Z->tauHtauL)       +240
   // (H->ZZ:Z->tauHtauL,Z->tauLtauL)       +260
   // (H->ZZ:Z->tauHtauL,Z->tauHtauH)       +280
   // (H->ZZ:Z->tauHtauL,Z->tauLtauH)       +300
   // (H->ZZ:Z->tauHtauL,Z->tauHtauL)       +320
   // (H->ZZ:Z->qq,Z->qq)                   +340
   // (H->ZZ:Z->ll,Z->qq)                   +360
   // (H->ZZ:Z->tauHtauH,Z->qq)             +380
   // (H->ZZ:Z->tauLtauL,Z->qq)             +400
   // (H->ZZ:Z->tauHtauL,Z->qq)             +420
   // (H->ZZ:Z->tauLtauH,Z->qq)             +440
   // (H->ZZ:Z->qq,Z->ll)                   +460
   // (H->ZZ:Z->qq,Z->tauHtauH)             +480
   // (H->ZZ:Z->qq,Z->tauLtauL)             +500
   // (H->ZZ:Z->qq,Z->tauHtauL)             +520
   // (H->ZZ:Z->qq,Z->tauLtauH)             +540
   // (H->ZZ:Z->ll,Z->nunu)                 +560
   // (H->ZZ:Z->tauHtauH,Z->nunu)           +580
   // (H->ZZ:Z->tauLtauL,Z->nunu)           +600
   // (H->ZZ:Z->tauHtauL,Z->nunu)           +620
   // (H->ZZ:Z->tauLtauH,Z->nunu)           +640
   // (H->ZZ:Z->qq,Z->nunu)                 +660
   // (H->ZZ:Z->nunu,Z->qq)                 +680
   // (H->ZZ:Z->nunu,Z->ll)                 +700
   // (H->ZZ:Z->nunu,Z->tauHtauH)           +720
   // (H->ZZ:Z->nunu,Z->tauLtauL)           +740
   // (H->ZZ:Z->nunu,Z->tauHtauL)           +760
   // (H->ZZ:Z->nunu,Z->tauLtauH)           +780
   // (H->ZZ:Z->nunu,Z->nunu)               +800
   // (H->ZZ:Z->tauHtauH,Z->ll)             +820
   // (H->ZZ:Z->tauLtauL,Z->ll)             +840
   // (H->ZZ:Z->tauHtauL,Z->ll)             +860
   // (H->ZZ:Z->tauLtauH,Z->ll)             +880
   // (H->ZZ:Z->ll,Z->tauHtauH)             +900
   // (H->ZZ:Z->ll,Z->tauLtauL)             +920
   // (H->ZZ:Z->ll,Z->tauHtauL)             +940
   // (H->ZZ:Z->ll,Z->tauLtauH)             +960
   
   // +2000
   // (H->tautau:tauL,tauL)     +2020
   // (H->tautau:tauH,tauH)     +2040
   // (H->tautau:tauL,tauH)     +2060
   // (H->tautau:tauH,tauL)     +2080

   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;

   int ipart = 0;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	int barcode = ipart; // in CMSSW barcode is the index of genParticle in the event
	// https://twiki.cern.ch/twiki/bin/view/CMS/GenParticles2HepMCConverter
	ipart++;
	
	// Additional partons (up to three)
	if( (fabs(mcp->pdgId()) <= 6 || fabs(mcp->pdgId()) == 21) &&
	    mcp->status() == 23 && barcode == 8 )
	  {
	     if( !j1 )
	       j1 = mcp;
	     else if( !j2 )
	       j2 = mcp;
	     else if( !j3 )
	       j3 = mcp;
	  }	

	// Higgs decays
	if( fabs(mcp->pdgId()) == 25 )
	  {
//	     if( iEvent.id().event() == 23539 )
//	       {
//		  std::cout << "found " << mcp->status() << " " << mcp->numberOfDaughters() << std::endl;
//	       }		  
	     
	     if( (mcp->status() == 62) ||
		 (mcp->status() == 3) )
	       {
		  h0 = const_cast<reco::GenParticle*>(mcp);

		  const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
		  for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
		    {
		       if( idr->isAvailable() ) 
			 {		       
			    const reco::GenParticleRef& genParticle = (*idr);
			    const reco::GenParticle *d = genParticle.get();
			    reco::GenParticle *pf = getUnique(d,0);

			    // h0 -> bbbar
			    if( fabs(pf->pdgId()) == 5 ) chan = 10000;
			    // h0 -> ee/mumu
			    if( fabs(pf->pdgId()) == 11 || fabs(pf->pdgId()) == 13 ) chan = 10001;
			    // h0 -> gg
			    if( fabs(pf->pdgId()) == 21 ) chan = 10002;
			    // h0 -> gammagamma
			    if( fabs(pf->pdgId()) == 22 ) chan = 10003;
			    // h0 -> qqbar (non-b)
			    if( fabs(pf->pdgId()) == 6 || fabs(pf->pdgId()) <= 4 ) chan = 10004;
			    
			    // h0 -> WW
			    if( fabs(pf->pdgId()) == 24 )
			      {
				 if( h0W1 && !h0W2 ) {h0W2 = pf;}
				 if( !h0W1 ) {h0W1 = pf;}				 

				 if( h0W1 && !h0W2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0W1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0W1_idr = daughterRefs.begin(); 
					  h0W1_idr!= daughterRefs.end(); ++h0W1_idr) 
					{
					   if( h0W1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0W1_idr);
						const reco::GenParticle *h0W1_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0W1_d,0);
						
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu
						  {
						     h0Wnu1 = pff;
						  }
						if( fabs(pff->pdgId()) == 16 ) // nutau
						  {
						     h0Wnutau1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     h0Wl1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     h0Wtau1 = pff;
						     
						     const reco::GenParticleRefVector& daughterRefs = h0Wtau1->daughterRefVector();
						     for(reco::GenParticleRefVector::const_iterator h0Wtau1_idr = daughterRefs.begin();
							 h0Wtau1_idr!= daughterRefs.end(); ++h0Wtau1_idr) 
						       {
							  if( h0Wtau1_idr->isAvailable() ) 
							    {		       
							       const reco::GenParticleRef& genParticle = (*h0Wtau1_idr);
							       const reco::GenParticle *h0Wtau1_d = genParticle.get();
							       reco::GenParticle *pfff = getUnique(h0Wtau1_d,0);
							       
							       if( fabs(pfff->pdgId()) == 12 ||
								   fabs(pfff->pdgId()) == 14 ) // nu
								 {
								    h0Wtaunu1 = pfff;
								 }
							       if( fabs(pfff->pdgId()) == 16 ) // nutau
								 {
								    h0Wtaunutau1 = pfff;
								 }		
							       if( fabs(pfff->pdgId()) == 11 ||
								   fabs(pfff->pdgId()) == 13 ) // l
								 {
								    h0Wtaul1 = pfff;  
								 }
							    }
						       }
						  }						
						if( fabs(pff->pdgId()) <= 6 )
						  {
						     if( h0Wq11 && !h0Wq21 ) h0Wq21 = pff;
						     if( !h0Wq11 ) h0Wq11 = pff;
						     if( h0Wq11_IS && !h0Wq21_IS ) h0Wq21_IS = const_cast<reco::GenParticle*>(h0W1_d);
						     if( !h0Wq11_IS ) h0Wq11_IS = const_cast<reco::GenParticle*>(h0W1_d);
						  }
					     }
					}				      
				   }
				 if( h0W2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0W2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0W2_idr = daughterRefs.begin(); 
					  h0W2_idr!= daughterRefs.end(); ++h0W2_idr) 
					{
					   if( h0W2_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0W2_idr);
						const reco::GenParticle *h0W2_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0W2_d,0);
						
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu
						  {
						     h0Wnu2 = pff;
						  }
						if( fabs(pff->pdgId()) == 16 ) // nutau
						  {
						     h0Wnutau2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     h0Wl2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     h0Wtau2 = pff;
						     
						     const reco::GenParticleRefVector& daughterRefs = h0Wtau2->daughterRefVector();
						     for(reco::GenParticleRefVector::const_iterator h0Wtau2_idr = daughterRefs.begin();
							 h0Wtau2_idr!= daughterRefs.end(); ++h0Wtau2_idr) 
						       {
							  if( h0Wtau2_idr->isAvailable() ) 
							    {		       
							       const reco::GenParticleRef& genParticle = (*h0Wtau2_idr);
							       const reco::GenParticle *h0Wtau2_d = genParticle.get();
							       reco::GenParticle *pfff = getUnique(h0Wtau2_d,0);
							       
							       if( fabs(pfff->pdgId()) == 12 ||
								   fabs(pfff->pdgId()) == 14 ) // nu
								 {
								    h0Wtaunu2 = pfff;
								 }
							       if( fabs(pfff->pdgId()) == 16 ) // nutau
								 {
								    h0Wtaunutau2 = pfff;
								 }		
							       if( fabs(pfff->pdgId()) == 11 ||
								   fabs(pfff->pdgId()) == 13 ) // l
								 {
								    h0Wtaul2 = pfff;  
								 }
							    }
						       }
						  }						
						if( fabs(pff->pdgId()) <= 6 )
						  {
						     if( h0Wq12 && !h0Wq22 ) h0Wq22 = pff;
						     if( !h0Wq12 ) h0Wq12 = pff;
						     if( h0Wq12_IS && !h0Wq22_IS ) h0Wq22_IS = const_cast<reco::GenParticle*>(h0W2_d);
						     if( !h0Wq12_IS ) h0Wq12_IS = const_cast<reco::GenParticle*>(h0W2_d);
						  }
					     }
					}				      				      
				   }				 
			      }

			    // h0 -> ZZ
			    if( fabs(pf->pdgId()) == 23 )
			      {
				 if( h0Z1 && !h0Z2 ) {h0Z2 = pf;}
				 if( !h0Z1 ) {h0Z1 = pf;}

				 if( h0Z1 && !h0Z2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0Z1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0Z1_idr = daughterRefs.begin(); 
					  h0Z1_idr!= daughterRefs.end(); ++h0Z1_idr) 
					{
					   if( h0Z1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0Z1_idr);
						const reco::GenParticle *h0Z1_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0Z1_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     if( h0Zl11 && !h0Zl21 ) {h0Zl21 = pff;}
						     if( !h0Zl11 ) {h0Zl11 = pff;}
						  }				
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     if( h0Ztau11 && !h0Ztau21 )
						       {
							  h0Ztau21 = pff;
							  
							  const reco::GenParticleRefVector& daughterRefs = h0Ztau21->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau21_idr = daughterRefs.begin(); 
							      h0Ztau21_idr!= daughterRefs.end(); ++h0Ztau21_idr) 
							    {
							       if( h0Ztau21_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau21_idr);
								    const reco::GenParticle *h0Ztau21_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau21_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu21 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau21 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul21 = pfff;
								      }
								 }
							    }							  
						       }
						     if( !h0Ztau11 )
						       {
							  h0Ztau11 = pff;

							  const reco::GenParticleRefVector& daughterRefs = h0Ztau11->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau11_idr = daughterRefs.begin(); 
							      h0Ztau11_idr!= daughterRefs.end(); ++h0Ztau11_idr) 
							    {
							       if( h0Ztau11_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau11_idr);
								    const reco::GenParticle *h0Ztau11_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau11_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu11 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau11 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul11 = pfff;
								      }
								 }
							    }							  
						       }						     
						  }
						if( fabs(pff->pdgId()) <= 6 ) // q
						  {
						     if( h0Zq11 && !h0Zq21 ) h0Zq21 = pff;
						     if( !h0Zq11 ) h0Zq11 = pff;
						     if( h0Zq11_IS && !h0Zq21_IS ) h0Zq21_IS = const_cast<reco::GenParticle*>(h0Z1_d);
						     if( !h0Zq11_IS ) h0Zq11_IS = const_cast<reco::GenParticle*>(h0Z1_d);
						  }				
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ||
						    fabs(pff->pdgId()) == 16 ) // nu
						  {
						     if( h0Znu11 && !h0Znu21 ) {h0Znu21 = pff;}
						     if( !h0Znu11 ) {h0Znu11 = pff;}
						  }						
					     }					   
					}				      
				   }
				 if( h0Z2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0Z2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0Z2_idr = daughterRefs.begin();
					  h0Z2_idr!= daughterRefs.end(); ++h0Z2_idr) 
					{
					   if( h0Z2_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0Z2_idr);
						const reco::GenParticle *h0Z2_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0Z2_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     if( h0Zl12 && !h0Zl22 ) {h0Zl22 = pff;}
						     if( !h0Zl12 ) {h0Zl12 = pff;}
						  }				
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     if( h0Ztau12 && !h0Ztau22 )
						       {
							  h0Ztau22 = pff;
							  
							  const reco::GenParticleRefVector& daughterRefs = h0Ztau22->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau22_idr = daughterRefs.begin(); 
							      h0Ztau22_idr!= daughterRefs.end(); ++h0Ztau22_idr) 
							    {
							       if( h0Ztau22_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau22_idr);
								    const reco::GenParticle *h0Ztau22_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau22_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu22 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau22 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul22 = pfff;
								      }
								 }
							    }							  
						       }
						     if( !h0Ztau12 )
						       {
							  h0Ztau12 = pff;

							  const reco::GenParticleRefVector& daughterRefs = h0Ztau12->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau12_idr = daughterRefs.begin(); 
							      h0Ztau12_idr!= daughterRefs.end(); ++h0Ztau12_idr)
							    {
							       if( h0Ztau12_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau12_idr);
								    const reco::GenParticle *h0Ztau12_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau12_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu12 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau12 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul12 = pfff;
								      }
								 }
							    }							  
						       }						     
						  }
						if( fabs(pff->pdgId()) <= 6 ) // q
						  {
						     if( h0Zq12 && !h0Zq22 ) h0Zq22 = pff;
						     if( !h0Zq12 ) h0Zq12 = pff;
						     if( h0Zq12_IS && !h0Zq22_IS ) h0Zq22_IS = const_cast<reco::GenParticle*>(h0Z2_d);
						     if( !h0Zq12_IS ) h0Zq12_IS = const_cast<reco::GenParticle*>(h0Z2_d);
						  }				
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ||
						    fabs(pff->pdgId()) == 16 ) // nu
						  {
						     if( h0Znu12 && !h0Znu22 ) {h0Znu22 = pff;}
						     if( !h0Znu12 ) {h0Znu12 = pff;}
						  }						
					     }					   
					}				      				      
				   }				 
			      }			    
			    
			    // h0 -> tautau
			    if( fabs(pf->pdgId()) == 15 && pf->status() == 2 ) // tau to decay
			      {
				 if( h0tau1 && !h0tau2 ) {h0tau2 = pf;}
				 if( !h0tau1 ) {h0tau1 = pf;}

				 if( h0tau1 && !h0tau2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0tau1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0tau1_idr = daughterRefs.begin();
					  h0tau1_idr!= daughterRefs.end(); ++h0tau1_idr) 
					{
					   if( h0tau1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0tau1_idr);
						const reco::GenParticle *h0tau1_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0tau1_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ||
						    fabs(pff->pdgId()) == 15 ) // l
						  {
						     h0taul1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 16 ) // nu_tau
						  {
						     h0taunutau1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu_e or nu_mu
						  {
						     h0taunu1 = pff;
						  }						
					     }
					}				      				      
				   }
				 if( h0tau2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0tau2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0tau2_idr = daughterRefs.begin();
					  h0tau2_idr!= daughterRefs.end(); ++h0tau2_idr) 
					{
					   if( h0tau2_idr->isAvailable() )
					     {		       
						const reco::GenParticleRef& genParticle = (*h0tau2_idr);
						const reco::GenParticle *h0tau2_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0tau2_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ||
						    fabs(pff->pdgId()) == 15 ) // l
						  {
						     h0taul2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 16 ) // nu_tau
						  {
						     h0taunutau2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu_e or nu_mu
						  {
						     h0taunu2 = pff;
						  }						
					     }
					}				      				      
				   }				 
			      }			    
			 }		       
		    }		  
	       }	     
	  }	

	// top decays
	if( fabs(mcp->pdgId()) == 6
	    && ( (mcp->status() == 62) || 
		 (mcp->status() == 3)
	       ) )
	  {
	     if( t1 && !t2 ) {t2 = const_cast<reco::GenParticle*>(mcp);}
	     if( !t1 ) {t1 = const_cast<reco::GenParticle*>(mcp);}

	     const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
	     for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
	       {
		  if( idr->isAvailable() ) 
		    {		       
		       const reco::GenParticleRef& genParticle = (*idr);
		       const reco::GenParticle *d = genParticle.get();
		       reco::GenParticle *pf = getUnique(d,0);

//		       if( pf->status() != 3 && pf->status() != 62 ) continue;
		       
		       if( fabs(pf->pdgId()) == 5 || fabs(pf->pdgId()) == 3 || fabs(pf->pdgId()) == 1 ) // b or s or d
			 {
			    if( tb1 && !tb2 ) tb2 = pf;
			    if( !tb1 ) tb1 = pf;
			    if( tb1_IS && !tb2_IS ) tb2_IS = const_cast<reco::GenParticle*>(d);
			    if( !tb1_IS ) tb1_IS = const_cast<reco::GenParticle*>(d);
			 }		       
		       
		       if( fabs(pf->pdgId()) == 24 ) // W
			 {
			    if( tW1 && !tW2 )
			      {
				 tW2 = pf;
				 const reco::GenParticleRefVector& tW2_daughterRefs = tW2->daughterRefVector();
				 for(reco::GenParticleRefVector::const_iterator tW2_idr = tW2_daughterRefs.begin();
				     tW2_idr!= tW2_daughterRefs.end(); ++tW2_idr) 
				   {
				      if( tW2_idr->isAvailable() ) 
					{		       
					   const reco::GenParticleRef& tW2_genParticle = (*tW2_idr);
					   const reco::GenParticle *tW2_d = tW2_genParticle.get();
					   reco::GenParticle *pff = getUnique(tW2_d,0);
					   
					   if( fabs(pff->pdgId()) == 12 ||
					       fabs(pff->pdgId()) == 14 ) // nu
					     {
						tWnu2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 16 ) // nu_tau
					     {
						tWnutau2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 11 ||
					       fabs(pff->pdgId()) == 13 ) // l
					     {
						tWl2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 15 ) // tau
					     {
						tWtau2 = pff;
						
						const reco::GenParticleRefVector& tWtau2_daughterRefs = tWtau2->daughterRefVector();
						for(reco::GenParticleRefVector::const_iterator tWtau2_idr = tWtau2_daughterRefs.begin();
						    tWtau2_idr!= tWtau2_daughterRefs.end(); ++tWtau2_idr) 
						  {
						     if( tWtau2_idr->isAvailable() ) 
						       {		       
							  const reco::GenParticleRef& tWtau2_genParticle = (*tWtau2_idr);
							  const reco::GenParticle *tWtau2_d = tWtau2_genParticle.get();
							  reco::GenParticle *pfff = getUnique(tWtau2_d,0);
							  
							  if( fabs(pfff->pdgId()) == 12 ||
							      fabs(pfff->pdgId()) == 14 ) // nu
							    {
							       tWtaunu2 = pfff;
							    }		
							  if( fabs(pfff->pdgId()) == 16 ) // nu_tau
							    {
							       tWtaunutau2 = pfff;
							    }			
							  if( fabs(pfff->pdgId()) == 11 ||
							      fabs(pfff->pdgId()) == 13 ) // l
							    {
							       tWtaul2 = pfff;
							    }							  
						       }
						  }						
					     }
					   if( fabs(pff->pdgId()) <= 6 ) // q
					     {
						if( tWq12 && !tWq22 ) tWq22 = pff;
						if( !tWq12 ) tWq12 = pff;
						if( tWq12_IS && !tWq22_IS ) tWq22_IS = const_cast<reco::GenParticle*>(tW2_d);
						if( !tWq12_IS ) tWq12_IS = const_cast<reco::GenParticle*>(tW2_d);
					     }					   					   
					}
				   }				
			      }
			    
			    if( !tW1 )
			      {
				 tW1 = pf;
				 const reco::GenParticleRefVector& tW1_daughterRefs = tW1->daughterRefVector();
				 for(reco::GenParticleRefVector::const_iterator tW1_idr = tW1_daughterRefs.begin();
				     tW1_idr!= tW1_daughterRefs.end(); ++tW1_idr) 
				   {
				      if( tW1_idr->isAvailable() ) 
					{		       
					   const reco::GenParticleRef& tW1_genParticle = (*tW1_idr);
					   const reco::GenParticle *tW1_d = tW1_genParticle.get();
					   reco::GenParticle *pff = getUnique(tW1_d,0);
					   
					   if( fabs(pff->pdgId()) == 12 ||
					       fabs(pff->pdgId()) == 14 ) // nu
					     {
						tWnu1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 16 ) // nu_tau
					     {
						tWnutau1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 11 ||
					       fabs(pff->pdgId()) == 13 ) // l
					     {
						tWl1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 15 ) // tau
					     {
						tWtau1 = pff;
						
						const reco::GenParticleRefVector& tWtau1_daughterRefs = tWtau1->daughterRefVector();
						for(reco::GenParticleRefVector::const_iterator tWtau1_idr = tWtau1_daughterRefs.begin();
						    tWtau1_idr!= tWtau1_daughterRefs.end(); ++tWtau1_idr) 
						  {
						     if( tWtau1_idr->isAvailable() ) 
						       {		       
							  const reco::GenParticleRef& tWtau1_genParticle = (*tWtau1_idr);
							  const reco::GenParticle *tWtau1_d = tWtau1_genParticle.get();
							  reco::GenParticle *pfff = getUnique(tWtau1_d,0);
							  
							  if( fabs(pfff->pdgId()) == 12 ||
							      fabs(pfff->pdgId()) == 14 ) // nu
							    {
							       tWtaunu1 = pfff;
							    }		
							  if( fabs(pfff->pdgId()) == 16 ) // nu_tau
							    {
							       tWtaunutau1 = pfff;
							    }			
							  if( fabs(pfff->pdgId()) == 11 ||
							      fabs(pfff->pdgId()) == 13 ) // l
							    {
							       tWtaul1 = pfff;
							    }							  
						       }
						  }						
					     }
					   if( fabs(pff->pdgId()) <= 6 ) // q
					     {
						if( tWq11 && !tWq21 ) tWq21 = pff;
						if( !tWq11 ) tWq11 = pff;
						if( tWq11_IS && !tWq21_IS ) tWq21_IS = const_cast<reco::GenParticle*>(tW1_d);
						if( !tWq11_IS ) tWq11_IS = const_cast<reco::GenParticle*>(tW1_d);
					     }					   					   
					}
				   }				
			      }			    
			 }		       
		    }		  
	       }	     
	  }
     }   

   bool doCheck = 0;

   if( h0 && t1 && t2 && tb1 && tb2 && tW1 && tW2 )
     {	
	int tchan = -666;
	if( tWl1 && tWl2 )   tchan = 0;
	if( tWq11 && tWq12 ) tchan = 1;
	if( tWq11 && tWl2 )  tchan = 2;
	if( tWl1 && tWq12 )  tchan = 3;
	if( tWtaul1 && tWtaul2 )  tchan = 4;
	if( tWq11 && tWtaul2 )  tchan = 5;
	if( tWtaul1 && tWq12 )  tchan = 6;
	if( tWtaunutau1 && tWtaunutau2 && !tWtaul1 && !tWtaul2 )  tchan = 7;
	if( tWq11 && tWtaunutau2 && !tWtaul2 )  tchan = 8;
	if( tWtaunutau1 && !tWtaul1 && tWq12 )  tchan = 9;
	if( tWtaunutau1 && tWtaul2 && !tWtaul1 )  tchan = 10;
	if( tWtaul1 && tWtaunutau2 && !tWtaul2 )  tchan = 11;
	if( tWl1 && tWtaul2 )  tchan = 12;
	if( tWl1 && !tWtaul2 && tWtaunutau2 )  tchan = 13;
	if( tWl2 && tWtaul1 )  tchan = 14;
	if( tWl2 && !tWtaul1 && tWtaunutau1 )  tchan = 15;
	
	if( tchan < 0 && doCheck )
	  {	     
	     std::cout << "Failed to identify top-quark decay chain" << std::endl;
	     
	     std::cout << "t1 = " << bool(t1) << std::endl;
	     std::cout << "t1->W = " << bool(tW1) << std::endl;
	     std::cout << "t1->W->l = " << bool(tWl1) << std::endl;
	     std::cout << "t1->W->nu = " << bool(tWnu1) << std::endl;
	     std::cout << "t1->W->nutau = " << bool(tWnutau1) << std::endl;
	     std::cout << "t1->W->tau = " << bool(tWtau1) << std::endl;
	     std::cout << "t1->W->tau->l = " << bool(tWtaul1) << std::endl;
	     std::cout << "t1->W->tau->nu = " << bool(tWtaunu1) << std::endl;
	     std::cout << "t1->W->tau->nutau = " << bool(tWtaunutau1) << std::endl;
	     std::cout << "t1->W->q = " << bool(tWq11) << std::endl;
	             
	     std::cout << "t2 = " << bool(t2) << std::endl;
	     std::cout << "t2->W = " << bool(tW2) << std::endl;
	     std::cout << "t2->W->l = " << bool(tWl2) << std::endl;
	     std::cout << "t2->W->nu = " << bool(tWnu2) << std::endl;
	     std::cout << "t2->W->nutau = " << bool(tWnutau2) << std::endl;
	     std::cout << "t2->W->tau = " << bool(tWtau2) << std::endl;
	     std::cout << "t2->W->tau->l = " << bool(tWtaul2) << std::endl;
	     std::cout << "t2->W->tau->nu = " << bool(tWtaunu2) << std::endl;
	     std::cout << "t2->W->tau->nutau = " << bool(tWtaunutau2) << std::endl;
	     std::cout << "t2->W->q = " << bool(tWq12) << std::endl;
	             
	     exit(1);
	  }
	
	if( h0W1 && h0W2 )
	  {	              
	     int chan0 = 0;
	     if( h0Wl1 && h0Wl2 ) chan = chan0 + 0 + tchan;
	     if( h0Wtaul1 && h0Wtaul2 ) chan = chan0 + 20 + tchan;
	     if( h0Wtaunutau1 && !h0Wtaul1 && h0Wtaunutau2 && !h0Wtaul2 ) chan = chan0 + 40 + tchan;
	     if( h0Wtaunutau1 && !h0Wtaul1 && h0Wtaul2 ) chan = chan0 + 60 + tchan;
	     if( h0Wtaul1 && h0Wtaunutau2 && !h0Wtaul2 ) chan = chan0 + 80 + tchan;
	     if( h0Wq11 && h0Wq12 ) chan = chan0 + 100 + tchan;
	     if( h0Wl1 && h0Wq12 ) chan = chan0 + 120 + tchan;
	     if( h0Wtaul1 && h0Wq12 ) chan = chan0 + 140 + tchan;
	     if( h0Wtaunutau1 && !h0Wtaul1 && h0Wq12 ) chan = chan0 + 160 + tchan;
	     if( h0Wq11 && h0Wl2 ) chan = chan0 + 180 + tchan;
	     if( h0Wq11 && h0Wtaul2 ) chan = chan0 + 200 + tchan;
	     if( h0Wq11 && h0Wtaunutau2 && !h0Wtaul2 ) chan = chan0 + 220 + tchan;
	     if( h0Wl1 && h0Wtaul2 ) chan = chan0 + 240 + tchan;
	     if( h0Wl1 && !h0Wtaul2 && h0Wtaunutau2 ) chan = chan0 + 260 + tchan;
	     if( h0Wl2 && h0Wtaul1 ) chan = chan0 + 280 + tchan;
	     if( h0Wl2 && !h0Wtaul1 && h0Wtaunutau1 ) chan = chan0 + 300 + tchan;
	  }

	if( h0Z1 && h0Z2 )
	  {
	     int chan0 = 1000;
	     if( h0Zl11 && h0Zl12 ) chan = chan0 + 0 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 20 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 40 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 60 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Ztaul22 && !h0Ztaul12 && h0Ztaunutau12 ) chan = chan0 + 80 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 100 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 120 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 140 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && !h0Ztaul12 && h0Ztaul22 && h0Ztaunutau12 ) chan = chan0 + 160 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 180 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 200 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 220 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && !h0Ztaul12 && h0Ztaul22 && h0Ztaunutau12 ) chan = chan0 + 240 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 260 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 280 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 300 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && !h0Ztaul12 && h0Ztaul22 && h0Ztaunutau12 ) chan = chan0 + 320 + tchan;
	     if( h0Zq11 && h0Zq12 ) chan = chan0 + 340 + tchan;
	     if( h0Zl11 && h0Zq12 ) chan = chan0 + 360 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && !h0Ztaul21 && h0Ztaunutau21 && h0Zq12 ) chan = chan0 + 380 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Zq12 ) chan = chan0 + 400 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Zq12 ) chan = chan0 + 420 + tchan;
	     if( h0Ztaul11 && h0Ztaunutau21 && !h0Ztaul21 && h0Zq12 ) chan = chan0 + 440 + tchan;
	     if( h0Zq11 && h0Zl12 ) chan = chan0 + 460 + tchan;
	     if( h0Zq11 && !h0Ztaul12 && h0Ztaunutau12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 480 + tchan;
	     if( h0Zq11 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 500 + tchan;
	     if( h0Zq11 && !h0Ztaul12 && h0Ztaunutau12 && h0Ztaul22 ) chan = chan0 + 520 + tchan;
	     if( h0Zq11 && h0Ztaul12 && h0Ztaunutau22 && !h0Ztaul22 ) chan = chan0 + 540 + tchan;
	     if( h0Zl11 && h0Znu12 ) chan = chan0 + 560 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && !h0Ztaul21 && h0Ztaunutau21 && h0Znu12 ) chan = chan0 + 580 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Znu12 ) chan = chan0 + 600 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Znu12 ) chan = chan0 + 620 + tchan;
	     if( h0Ztaul11 && h0Ztaunutau21 && !h0Ztaul21 && h0Znu12 ) chan = chan0 + 640 + tchan;
	     if( h0Zq11 && h0Znu12 ) chan = chan0 + 660 + tchan;
	     if( h0Znu11 && h0Zq12 ) chan = chan0 + 680 + tchan;
	     if( h0Znu11 && h0Zl12 ) chan = chan0 + 700 + tchan;
	     if( h0Znu11 && !h0Ztaul12 && h0Ztaunutau12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 720 + tchan;
	     if( h0Znu11 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 740 + tchan;
	     if( h0Znu11 && !h0Ztaul12 && h0Ztaunutau12 && h0Ztaul22 ) chan = chan0 + 760 + tchan;
	     if( h0Znu11 && h0Ztaul12 && h0Ztaunutau22 && !h0Ztaul22 ) chan = chan0 + 780 + tchan;
	     if( h0Znu11 && h0Znu12 ) chan = chan0 + 800 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && !h0Ztaul21 && h0Ztaunutau21 && h0Zl12 ) chan = chan0 + 820 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Zl12 ) chan = chan0 + 840 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Zl12 ) chan = chan0 + 860 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && h0Zl12 ) chan = chan0 + 880 + tchan;
	     if( h0Zl11 && !h0Ztaul12 && h0Ztaunutau12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 900 + tchan;
	     if( h0Zl11 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 920 + tchan;
	     if( h0Zl11 && !h0Ztaul12 && h0Ztaunutau12 && h0Ztaul22 ) chan = chan0 + 940 + tchan;
	     if( h0Zl11 && h0Ztaul12 && h0Ztaunutau22 && !h0Ztaul22 ) chan = chan0 + 960 + tchan;	     
	  }	
	
	if( h0tau1 && h0tau2 )
	  {	     
	     int chan0 = 2000;
	     if( h0taul1 && h0taul2 && h0taunu1 && h0taunu2 ) chan = chan0 + 20 + tchan;
	     if( !h0taul1 && !h0taul2 && h0taunutau1 && h0taunutau2 ) chan = chan0 + 40 + tchan;
	     if( h0taul1 && !h0taul2 && h0taunutau2 && h0taunu1 ) chan = chan0 + 60 + tchan;
	     if( h0taul2 && !h0taul1 && h0taunutau1 && h0taunu2 ) chan = chan0 + 80 + tchan;
	  }	

	if( chan < 0 && doCheck )
	  {
	     std::cout << "Unknown channel found" << std::endl;
	     std::cout << "chan = " << chan << std::endl;

	     std::cout << "j1 = " << bool(j1) << std::endl;
	     std::cout << "j2 = " << bool(j2) << std::endl;
	     std::cout << "j3 = " << bool(j3) << std::endl;
	     
	     std::cout << "h0 = " << bool(h0) << std::endl;
	     std::cout << "h0->W1 = " << bool(h0W1) << std::endl;
	     std::cout << "h0->W1->l = " << bool(h0Wl1) << std::endl;
	     std::cout << "h0->W1->nu = " << bool(h0Wnu1) << std::endl;
	     std::cout << "h0->W1->tau = " << bool(h0Wtau1) << std::endl;
	     std::cout << "h0->W1->nutau = " << bool(h0Wnutau1) << std::endl;
	     std::cout << "h0->W1->tau->l = " << bool(h0Wtaul1) << std::endl;
	     std::cout << "h0->W1->tau->nu = " << bool(h0Wtaunu1) << std::endl;
	     std::cout << "h0->W1->tau->nutau = " << bool(h0Wtaunutau1) << std::endl;
	     std::cout << "h0->W1->q1 = " << bool(h0Wq11) << std::endl;
	     std::cout << "h0->W1->q2 = " << bool(h0Wq21) << std::endl;
	     std::cout << "h0->W2 = " << bool(h0W2) << std::endl;
	     std::cout << "h0->W2->l = " << bool(h0Wl2) << std::endl;
	     std::cout << "h0->W2->nu = " << bool(h0Wnu2) << std::endl;
	     std::cout << "h0->W2->tau = " << bool(h0Wtau2) << std::endl;
	     std::cout << "h0->W2->nutau = " << bool(h0Wnutau2) << std::endl;
	     std::cout << "h0->W2->tau->l = " << bool(h0Wtaul2) << std::endl;
	     std::cout << "h0->W2->tau->nu = " << bool(h0Wtaunu2) << std::endl;
	     std::cout << "h0->W2->tau->nutau = " << bool(h0Wtaunutau2) << std::endl;
	     std::cout << "h0->W2->q1 = " << bool(h0Wq12) << std::endl;
	     std::cout << "h0->W2->q2 = " << bool(h0Wq22) << std::endl;
	     
	     std::cout << "h0->Z1 = " << bool(h0Z1) << std::endl;
	     std::cout << "h0->Z1->l1 = " << bool(h0Zl11) << std::endl;
	     std::cout << "h0->Z1->l2 = " << bool(h0Zl21) << std::endl;
	     std::cout << "h0->Z1->tau1 = " << bool(h0Ztau11) << std::endl;
	     std::cout << "h0->Z1->tau1->l = " << bool(h0Ztaul11) << std::endl;
	     std::cout << "h0->Z1->tau1->nu = " << bool(h0Ztaunu11) << std::endl;
	     std::cout << "h0->Z1->tau1->nutau = " << bool(h0Ztaunutau11) << std::endl;
	     std::cout << "h0->Z1->tau2 = " << bool(h0Ztau21) << std::endl;
	     std::cout << "h0->Z1->tau2->l = " << bool(h0Ztaul21) << std::endl;
	     std::cout << "h0->Z1->tau2->nu = " << bool(h0Ztaunu21) << std::endl;
	     std::cout << "h0->Z1->tau2->nutau = " << bool(h0Ztaunutau21) << std::endl;
	     std::cout << "h0->Z1->q1 = " << bool(h0Zq11) << std::endl;
	     std::cout << "h0->Z1->q2 = " << bool(h0Zq21) << std::endl;
	     std::cout << "h0->Z1->nu1 = " << bool(h0Znu11) << std::endl;
	     std::cout << "h0->Z1->nu2 = " << bool(h0Znu21) << std::endl;
	     
	     std::cout << "h0->Z2 = " << bool(h0Z2) << std::endl;
	     std::cout << "h0->Z2->l1 = " << bool(h0Zl12) << std::endl;
	     std::cout << "h0->Z2->l2 = " << bool(h0Zl22) << std::endl;
	     std::cout << "h0->Z2->tau1 = " << bool(h0Ztau12) << std::endl;
	     std::cout << "h0->Z2->tau1->l = " << bool(h0Ztaul12) << std::endl;
	     std::cout << "h0->Z2->tau1->nu = " << bool(h0Ztaunu12) << std::endl;
	     std::cout << "h0->Z2->tau1->nutau = " << bool(h0Ztaunutau12) << std::endl;
	     std::cout << "h0->Z2->tau2 = " << bool(h0Ztau22) << std::endl;
	     std::cout << "h0->Z2->tau2->l = " << bool(h0Ztaul22) << std::endl;
	     std::cout << "h0->Z2->tau2->nu = " << bool(h0Ztaunu22) << std::endl;
	     std::cout << "h0->Z2->tau2->nutau = " << bool(h0Ztaunutau22) << std::endl;
	     std::cout << "h0->Z2->q1 = " << bool(h0Zq12) << std::endl;
	     std::cout << "h0->Z2->q2 = " << bool(h0Zq22) << std::endl;
	     std::cout << "h0->Z2->nu1 = " << bool(h0Znu12) << std::endl;
	     std::cout << "h0->Z2->nu2 = " << bool(h0Znu22) << std::endl;
	     
	     std::cout << "h0->tau1 = " << bool(h0tau1) << std::endl;
	     std::cout << "h0->tau1->l = " << bool(h0taul1) << std::endl;
	     std::cout << "h0->tau1->nu_tau = " << bool(h0taunutau1) << std::endl;
	     std::cout << "h0->tau1->nu = " << bool(h0taunu1) << std::endl;
	             
	     std::cout << "h0->tau2 = " << bool(h0tau2) << std::endl;
	     std::cout << "h0->tau2->l = " << bool(h0taul2) << std::endl;
	     std::cout << "h0->tau2->nu_tau = " << bool(h0taunutau2) << std::endl;
	     std::cout << "h0->tau2->nu = " << bool(h0taunu2) << std::endl;
	     
	     exit(1);
	  }
     }

   tree.mc_truth_tth_channel = chan;
   
   // TLV

   if( h0 ) p4toTLV(h0->p4(),tree.mc_truth_h0_p4);      
   
   if( h0W1 ) p4toTLV(h0W1->p4(),tree.mc_truth_h0W1_p4);
   if( h0W2 ) p4toTLV(h0W2->p4(),tree.mc_truth_h0W2_p4);
   if( h0Wl1 ) p4toTLV(h0Wl1->p4(),tree.mc_truth_h0Wl1_p4);
   if( h0Wnu1 ) p4toTLV(h0Wnu1->p4(),tree.mc_truth_h0Wnu1_p4);
   if( h0Wtau1 ) p4toTLV(h0Wtau1->p4(),tree.mc_truth_h0Wtau1_p4);
   if( h0Wnutau1 ) p4toTLV(h0Wnutau1->p4(),tree.mc_truth_h0Wnutau1_p4);
   if( h0Wtaul1 ) p4toTLV(h0Wtaul1->p4(),tree.mc_truth_h0Wtaul1_p4);
   if( h0Wtaunu1 ) p4toTLV(h0Wtaunu1->p4(),tree.mc_truth_h0Wtaunu1_p4);
   if( h0Wtaunutau1 ) p4toTLV(h0Wtaunutau1->p4(),tree.mc_truth_h0Wtaunutau1_p4);
   if( h0Wl2 ) p4toTLV(h0Wl2->p4(),tree.mc_truth_h0Wl2_p4);
   if( h0Wnu2 ) p4toTLV(h0Wnu2->p4(),tree.mc_truth_h0Wnu2_p4);
   if( h0Wtau2 ) p4toTLV(h0Wtau2->p4(),tree.mc_truth_h0Wtau2_p4);
   if( h0Wnutau2 ) p4toTLV(h0Wnutau2->p4(),tree.mc_truth_h0Wnutau2_p4);
   if( h0Wtaul2 ) p4toTLV(h0Wtaul2->p4(),tree.mc_truth_h0Wtaul2_p4);
   if( h0Wtaunu2 ) p4toTLV(h0Wtaunu2->p4(),tree.mc_truth_h0Wtaunu2_p4);
   if( h0Wtaunutau2 ) p4toTLV(h0Wtaunutau2->p4(),tree.mc_truth_h0Wtaunutau2_p4);
   if( h0Wq11 ) p4toTLV(h0Wq11->p4(),tree.mc_truth_h0Wq11_p4);
   if( h0Wq21 ) p4toTLV(h0Wq21->p4(),tree.mc_truth_h0Wq21_p4);
   if( h0Wq12 ) p4toTLV(h0Wq12->p4(),tree.mc_truth_h0Wq12_p4);
   if( h0Wq22 ) p4toTLV(h0Wq22->p4(),tree.mc_truth_h0Wq22_p4);
   if( h0Wq11_IS ) p4toTLV(h0Wq11_IS->p4(),tree.mc_truth_h0Wq11_IS_p4);
   if( h0Wq21_IS ) p4toTLV(h0Wq21_IS->p4(),tree.mc_truth_h0Wq21_IS_p4);
   if( h0Wq12_IS ) p4toTLV(h0Wq12_IS->p4(),tree.mc_truth_h0Wq12_IS_p4);
   if( h0Wq22_IS ) p4toTLV(h0Wq22_IS->p4(),tree.mc_truth_h0Wq22_IS_p4);
   
   if( h0Z1 ) p4toTLV(h0Z1->p4(),tree.mc_truth_h0Z1_p4);
   if( h0Z2 ) p4toTLV(h0Z2->p4(),tree.mc_truth_h0Z2_p4);
   if( h0Zl11 ) p4toTLV(h0Zl11->p4(),tree.mc_truth_h0Zl11_p4);
   if( h0Zl21 ) p4toTLV(h0Zl21->p4(),tree.mc_truth_h0Zl21_p4);
   if( h0Ztau11 ) p4toTLV(h0Ztau11->p4(),tree.mc_truth_h0Ztau11_p4);
   if( h0Ztau21 ) p4toTLV(h0Ztau21->p4(),tree.mc_truth_h0Ztau21_p4);
   if( h0Ztaul11 ) p4toTLV(h0Ztaul11->p4(),tree.mc_truth_h0Ztaul11_p4);
   if( h0Ztaul21 ) p4toTLV(h0Ztaul21->p4(),tree.mc_truth_h0Ztaul21_p4);
   if( h0Ztaunu11 ) p4toTLV(h0Ztaunu11->p4(),tree.mc_truth_h0Ztaunu11_p4);
   if( h0Ztaunu21 ) p4toTLV(h0Ztaunu21->p4(),tree.mc_truth_h0Ztaunu21_p4);
   if( h0Ztaunutau11 ) p4toTLV(h0Ztaunutau11->p4(),tree.mc_truth_h0Ztaunutau11_p4);
   if( h0Ztaunutau21 ) p4toTLV(h0Ztaunutau21->p4(),tree.mc_truth_h0Ztaunutau21_p4);
   if( h0Zq11 ) p4toTLV(h0Zq11->p4(),tree.mc_truth_h0Zq11_p4);
   if( h0Zq21 ) p4toTLV(h0Zq21->p4(),tree.mc_truth_h0Zq21_p4);
   if( h0Zq11_IS ) p4toTLV(h0Zq11_IS->p4(),tree.mc_truth_h0Zq11_IS_p4);
   if( h0Zq21_IS ) p4toTLV(h0Zq21_IS->p4(),tree.mc_truth_h0Zq21_IS_p4);
   if( h0Zl12 ) p4toTLV(h0Zl12->p4(),tree.mc_truth_h0Zl12_p4);
   if( h0Zl22 ) p4toTLV(h0Zl22->p4(),tree.mc_truth_h0Zl22_p4);
   if( h0Ztau12 ) p4toTLV(h0Ztau12->p4(),tree.mc_truth_h0Ztau12_p4);
   if( h0Ztau22 ) p4toTLV(h0Ztau22->p4(),tree.mc_truth_h0Ztau22_p4);
   if( h0Ztaul12 ) p4toTLV(h0Ztaul12->p4(),tree.mc_truth_h0Ztaul12_p4);
   if( h0Ztaul22 ) p4toTLV(h0Ztaul22->p4(),tree.mc_truth_h0Ztaul22_p4);
   if( h0Ztaunu12 ) p4toTLV(h0Ztaunu12->p4(),tree.mc_truth_h0Ztaunu12_p4);
   if( h0Ztaunu22 ) p4toTLV(h0Ztaunu22->p4(),tree.mc_truth_h0Ztaunu22_p4);
   if( h0Ztaunutau12 ) p4toTLV(h0Ztaunutau12->p4(),tree.mc_truth_h0Ztaunutau12_p4);
   if( h0Ztaunutau22 ) p4toTLV(h0Ztaunutau22->p4(),tree.mc_truth_h0Ztaunutau22_p4);
   if( h0Zq12 ) p4toTLV(h0Zq12->p4(),tree.mc_truth_h0Zq12_p4);
   if( h0Zq22 ) p4toTLV(h0Zq22->p4(),tree.mc_truth_h0Zq22_p4);
   if( h0Zq12_IS ) p4toTLV(h0Zq12_IS->p4(),tree.mc_truth_h0Zq12_IS_p4);
   if( h0Zq22_IS ) p4toTLV(h0Zq22_IS->p4(),tree.mc_truth_h0Zq22_IS_p4);
   if( h0Znu11 ) p4toTLV(h0Znu11->p4(),tree.mc_truth_h0Znu11_p4);
   if( h0Znu21 ) p4toTLV(h0Znu21->p4(),tree.mc_truth_h0Znu21_p4);
   if( h0Znu12 ) p4toTLV(h0Znu12->p4(),tree.mc_truth_h0Znu12_p4);
   if( h0Znu22 ) p4toTLV(h0Znu22->p4(),tree.mc_truth_h0Znu22_p4);
   
   if( h0tau1 ) p4toTLV(h0tau1->p4(),tree.mc_truth_h0tau1_p4);
   if( h0tau2 ) p4toTLV(h0tau2->p4(),tree.mc_truth_h0tau2_p4);
   if( h0taul1 ) p4toTLV(h0taul1->p4(),tree.mc_truth_h0taul1_p4);
   if( h0taunutau1 ) p4toTLV(h0taunutau1->p4(),tree.mc_truth_h0taunutau1_p4);
   if( h0taunu1 ) p4toTLV(h0taunu1->p4(),tree.mc_truth_h0taunu1_p4);
   if( h0taul2 ) p4toTLV(h0taul2->p4(),tree.mc_truth_h0taul2_p4);
   if( h0taunutau2 ) p4toTLV(h0taunutau2->p4(),tree.mc_truth_h0taunutau2_p4);
   if( h0taunu2 ) p4toTLV(h0taunu2->p4(),tree.mc_truth_h0taunu2_p4);
   
   if( t1 ) p4toTLV(t1->p4(),tree.mc_truth_t1_p4);
   if( t2 ) p4toTLV(t2->p4(),tree.mc_truth_t2_p4);
   if( tb1 ) p4toTLV(tb1->p4(),tree.mc_truth_tb1_p4);
   if( tb2 ) p4toTLV(tb2->p4(),tree.mc_truth_tb2_p4);
   if( tb1_IS ) p4toTLV(tb1_IS->p4(),tree.mc_truth_tb1_IS_p4);
   if( tb2_IS ) p4toTLV(tb2_IS->p4(),tree.mc_truth_tb2_IS_p4);
   
   if( tW1 ) p4toTLV(tW1->p4(),tree.mc_truth_tW1_p4);
   if( tWnu1 ) p4toTLV(tWnu1->p4(),tree.mc_truth_tWnu1_p4);
   if( tWnutau1 ) p4toTLV(tWnutau1->p4(),tree.mc_truth_tWnutau1_p4);
   if( tWl1 ) p4toTLV(tWl1->p4(),tree.mc_truth_tWl1_p4);
   if( tWtau1 ) p4toTLV(tWtau1->p4(),tree.mc_truth_tWtau1_p4);
   if( tWtaunu1 ) p4toTLV(tWtaunu1->p4(),tree.mc_truth_tWtaunu1_p4);
   if( tWtaunutau1 ) p4toTLV(tWtaunutau1->p4(),tree.mc_truth_tWtaunutau1_p4);
   if( tWtaul1 ) p4toTLV(tWtaul1->p4(),tree.mc_truth_tWtaul1_p4);
   if( tWq11 ) p4toTLV(tWq11->p4(),tree.mc_truth_tWq11_p4);
   if( tWq21 ) p4toTLV(tWq21->p4(),tree.mc_truth_tWq21_p4);
   if( tWq11_IS ) p4toTLV(tWq11_IS->p4(),tree.mc_truth_tWq11_IS_p4);
   if( tWq21_IS ) p4toTLV(tWq21_IS->p4(),tree.mc_truth_tWq21_IS_p4);

   if( tW2 ) p4toTLV(tW2->p4(),tree.mc_truth_tW2_p4);
   if( tWnu2 ) p4toTLV(tWnu2->p4(),tree.mc_truth_tWnu2_p4);
   if( tWnutau2 ) p4toTLV(tWnutau2->p4(),tree.mc_truth_tWnutau2_p4);
   if( tWl2 ) p4toTLV(tWl2->p4(),tree.mc_truth_tWl2_p4);
   if( tWtau2 ) p4toTLV(tWtau2->p4(),tree.mc_truth_tWtau2_p4);
   if( tWtaunu2 ) p4toTLV(tWtaunu2->p4(),tree.mc_truth_tWtaunu2_p4);
   if( tWtaunutau2 ) p4toTLV(tWtaunutau2->p4(),tree.mc_truth_tWtaunutau2_p4);
   if( tWtaul2 ) p4toTLV(tWtaul2->p4(),tree.mc_truth_tWtaul2_p4);
   if( tWq12 ) p4toTLV(tWq12->p4(),tree.mc_truth_tWq12_p4);
   if( tWq22 ) p4toTLV(tWq22->p4(),tree.mc_truth_tWq22_p4);
   if( tWq12_IS ) p4toTLV(tWq12_IS->p4(),tree.mc_truth_tWq12_IS_p4);
   if( tWq22_IS ) p4toTLV(tWq22_IS->p4(),tree.mc_truth_tWq22_IS_p4);

   if( j1 ) p4toTLV(j1->p4(),tree.mc_truth_j1_p4);
   if( j2 ) p4toTLV(j2->p4(),tree.mc_truth_j2_p4);
   if( j3 ) p4toTLV(j3->p4(),tree.mc_truth_j3_p4);
   
   // pt

   if( h0 ) tree.mc_truth_h0_pt = h0->p4().pt();

   if( h0W1 ) tree.mc_truth_h0W1_pt = h0W1->p4().pt();
   if( h0W2 ) tree.mc_truth_h0W2_pt = h0W2->p4().pt();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_pt = h0Wl1->p4().pt();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_pt = h0Wnu1->p4().pt();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_pt = h0Wtau1->p4().pt();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_pt = h0Wnutau1->p4().pt();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_pt = h0Wtaul1->p4().pt();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_pt = h0Wtaunu1->p4().pt();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_pt = h0Wtaunutau1->p4().pt();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_pt = h0Wl2->p4().pt();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_pt = h0Wnu2->p4().pt();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_pt = h0Wtau2->p4().pt();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_pt = h0Wnutau2->p4().pt();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_pt = h0Wtaul2->p4().pt();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_pt = h0Wtaunu2->p4().pt();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_pt = h0Wtaunutau2->p4().pt();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_pt = h0Wq11->p4().pt();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_pt = h0Wq21->p4().pt();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_pt = h0Wq12->p4().pt();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_pt = h0Wq22->p4().pt();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_pt = h0Wq11_IS->p4().pt();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_pt = h0Wq21_IS->p4().pt();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_pt = h0Wq12_IS->p4().pt();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_pt = h0Wq22_IS->p4().pt();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_pt = h0Z1->p4().pt();
   if( h0Z2 ) tree.mc_truth_h0Z2_pt = h0Z2->p4().pt();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_pt = h0Zl11->p4().pt();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_pt = h0Zl21->p4().pt();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_pt = h0Zl12->p4().pt();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_pt = h0Zl22->p4().pt();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_pt = h0Ztau11->p4().pt();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_pt = h0Ztau21->p4().pt();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_pt = h0Ztaul11->p4().pt();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_pt = h0Ztaul21->p4().pt();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_pt = h0Ztaunu11->p4().pt();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_pt = h0Ztaunu21->p4().pt();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_pt = h0Ztaunutau11->p4().pt();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_pt = h0Ztaunutau21->p4().pt();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_pt = h0Zq11->p4().pt();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_pt = h0Zq21->p4().pt();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_pt = h0Zq12->p4().pt();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_pt = h0Zq22->p4().pt();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_pt = h0Zq11_IS->p4().pt();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_pt = h0Zq21_IS->p4().pt();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_pt = h0Zq12_IS->p4().pt();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_pt = h0Zq22_IS->p4().pt();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_pt = h0Ztau12->p4().pt();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_pt = h0Ztau22->p4().pt();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_pt = h0Ztaul12->p4().pt();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_pt = h0Ztaul22->p4().pt();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_pt = h0Ztaunu12->p4().pt();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_pt = h0Ztaunu22->p4().pt();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_pt = h0Ztaunutau12->p4().pt();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_pt = h0Ztaunutau22->p4().pt();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_pt = h0Znu11->p4().pt();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_pt = h0Znu21->p4().pt();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_pt = h0Znu12->p4().pt();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_pt = h0Znu22->p4().pt();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_pt = h0tau1->p4().pt();
   if( h0tau2 ) tree.mc_truth_h0tau2_pt = h0tau2->p4().pt();
   if( h0taul1 ) tree.mc_truth_h0taul1_pt = h0taul1->p4().pt();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_pt = h0taunutau1->p4().pt();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_pt = h0taunu1->p4().pt();
   if( h0taul2 ) tree.mc_truth_h0taul2_pt = h0taul2->p4().pt();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_pt = h0taunutau2->p4().pt();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_pt = h0taunu2->p4().pt();
   
   if( t1 ) tree.mc_truth_t1_pt = t1->p4().pt();
   if( t2 ) tree.mc_truth_t2_pt = t2->p4().pt();
   if( tb1 ) tree.mc_truth_tb1_pt = tb1->p4().pt();
   if( tb2 ) tree.mc_truth_tb2_pt = tb2->p4().pt();
   if( tb1_IS ) tree.mc_truth_tb1_IS_pt = tb1_IS->p4().pt();
   if( tb2_IS ) tree.mc_truth_tb2_IS_pt = tb2_IS->p4().pt();
   
   if( tW1 ) tree.mc_truth_tW1_pt = tW1->p4().pt();
   if( tWnu1 ) tree.mc_truth_tWnu1_pt = tWnu1->p4().pt();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_pt = tWnutau1->p4().pt();
   if( tWl1 ) tree.mc_truth_tWl1_pt = tWl1->p4().pt();
   if( tWtau1 ) tree.mc_truth_tWtau1_pt = tWtau1->p4().pt();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_pt = tWtaunu1->p4().pt();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_pt = tWtaunutau1->p4().pt();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_pt = tWtaul1->p4().pt();
   if( tWq11 ) tree.mc_truth_tWq11_pt = tWq11->p4().pt();
   if( tWq21 ) tree.mc_truth_tWq21_pt = tWq21->p4().pt();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_pt = tWq11_IS->p4().pt();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_pt = tWq21_IS->p4().pt();
   
   if( tW2 ) tree.mc_truth_tW2_pt = tW2->p4().pt();
   if( tWnu2 ) tree.mc_truth_tWnu2_pt = tWnu2->p4().pt();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_pt = tWnutau2->p4().pt();
   if( tWl2 ) tree.mc_truth_tWl2_pt = tWl2->p4().pt();
   if( tWtau2 ) tree.mc_truth_tWtau2_pt = tWtau2->p4().pt();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_pt = tWtaunu2->p4().pt();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_pt = tWtaunutau2->p4().pt();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_pt = tWtaul2->p4().pt();
   if( tWq12 ) tree.mc_truth_tWq12_pt = tWq12->p4().pt();
   if( tWq22 ) tree.mc_truth_tWq22_pt = tWq22->p4().pt();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_pt = tWq12_IS->p4().pt();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_pt = tWq22_IS->p4().pt();
   
   if( j1 ) tree.mc_truth_j1_pt = j1->p4().pt();
   if( j2 ) tree.mc_truth_j2_pt = j2->p4().pt();
   if( j3 ) tree.mc_truth_j3_pt = j3->p4().pt();

   // eta

   if( h0 ) tree.mc_truth_h0_eta = h0->p4().eta();

   if( h0W1 ) tree.mc_truth_h0W1_eta = h0W1->p4().eta();
   if( h0W2 ) tree.mc_truth_h0W2_eta = h0W2->p4().eta();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_eta = h0Wl1->p4().eta();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_eta = h0Wnu1->p4().eta();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_eta = h0Wtau1->p4().eta();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_eta = h0Wnutau1->p4().eta();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_eta = h0Wtaul1->p4().eta();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_eta = h0Wtaunu1->p4().eta();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_eta = h0Wtaunutau1->p4().eta();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_eta = h0Wl2->p4().eta();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_eta = h0Wnu2->p4().eta();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_eta = h0Wtau2->p4().eta();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_eta = h0Wnutau2->p4().eta();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_eta = h0Wtaul2->p4().eta();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_eta = h0Wtaunu2->p4().eta();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_eta = h0Wtaunutau2->p4().eta();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_eta = h0Wq11->p4().eta();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_eta = h0Wq21->p4().eta();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_eta = h0Wq12->p4().eta();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_eta = h0Wq22->p4().eta();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_eta = h0Wq11_IS->p4().eta();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_eta = h0Wq21_IS->p4().eta();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_eta = h0Wq12_IS->p4().eta();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_eta = h0Wq22_IS->p4().eta();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_eta = h0Z1->p4().eta();
   if( h0Z2 ) tree.mc_truth_h0Z2_eta = h0Z2->p4().eta();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_eta = h0Zl11->p4().eta();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_eta = h0Zl21->p4().eta();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_eta = h0Zl12->p4().eta();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_eta = h0Zl22->p4().eta();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_eta = h0Ztau11->p4().eta();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_eta = h0Ztau21->p4().eta();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_eta = h0Ztaul11->p4().eta();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_eta = h0Ztaul21->p4().eta();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_eta = h0Ztaunu11->p4().eta();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_eta = h0Ztaunu21->p4().eta();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_eta = h0Ztaunutau11->p4().eta();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_eta = h0Ztaunutau21->p4().eta();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_eta = h0Zq11->p4().eta();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_eta = h0Zq21->p4().eta();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_eta = h0Zq12->p4().eta();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_eta = h0Zq22->p4().eta();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_eta = h0Zq11_IS->p4().eta();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_eta = h0Zq21_IS->p4().eta();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_eta = h0Zq12_IS->p4().eta();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_eta = h0Zq22_IS->p4().eta();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_eta = h0Ztau12->p4().eta();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_eta = h0Ztau22->p4().eta();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_eta = h0Ztaul12->p4().eta();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_eta = h0Ztaul22->p4().eta();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_eta = h0Ztaunu12->p4().eta();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_eta = h0Ztaunu22->p4().eta();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_eta = h0Ztaunutau12->p4().eta();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_eta = h0Ztaunutau22->p4().eta();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_eta = h0Znu11->p4().eta();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_eta = h0Znu21->p4().eta();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_eta = h0Znu12->p4().eta();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_eta = h0Znu22->p4().eta();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_eta = h0tau1->p4().eta();
   if( h0tau2 ) tree.mc_truth_h0tau2_eta = h0tau2->p4().eta();
   if( h0taul1 ) tree.mc_truth_h0taul1_eta = h0taul1->p4().eta();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_eta = h0taunutau1->p4().eta();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_eta = h0taunu1->p4().eta();
   if( h0taul2 ) tree.mc_truth_h0taul2_eta = h0taul2->p4().eta();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_eta = h0taunutau2->p4().eta();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_eta = h0taunu2->p4().eta();
   
   if( t1 ) tree.mc_truth_t1_eta = t1->p4().eta();
   if( t2 ) tree.mc_truth_t2_eta = t2->p4().eta();
   if( tb1 ) tree.mc_truth_tb1_eta = tb1->p4().eta();
   if( tb2 ) tree.mc_truth_tb2_eta = tb2->p4().eta();
   if( tb1_IS ) tree.mc_truth_tb1_IS_eta = tb1_IS->p4().eta();
   if( tb2_IS ) tree.mc_truth_tb2_IS_eta = tb2_IS->p4().eta();
   
   if( tW1 ) tree.mc_truth_tW1_eta = tW1->p4().eta();
   if( tWnu1 ) tree.mc_truth_tWnu1_eta = tWnu1->p4().eta();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_eta = tWnutau1->p4().eta();
   if( tWl1 ) tree.mc_truth_tWl1_eta = tWl1->p4().eta();
   if( tWtau1 ) tree.mc_truth_tWtau1_eta = tWtau1->p4().eta();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_eta = tWtaunu1->p4().eta();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_eta = tWtaunutau1->p4().eta();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_eta = tWtaul1->p4().eta();
   if( tWq11 ) tree.mc_truth_tWq11_eta = tWq11->p4().eta();
   if( tWq21 ) tree.mc_truth_tWq21_eta = tWq21->p4().eta();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_eta = tWq11_IS->p4().eta();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_eta = tWq21_IS->p4().eta();
   
   if( tW2 ) tree.mc_truth_tW2_eta = tW2->p4().eta();
   if( tWnu2 ) tree.mc_truth_tWnu2_eta = tWnu2->p4().eta();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_eta = tWnutau2->p4().eta();
   if( tWl2 ) tree.mc_truth_tWl2_eta = tWl2->p4().eta();
   if( tWtau2 ) tree.mc_truth_tWtau2_eta = tWtau2->p4().eta();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_eta = tWtaunu2->p4().eta();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_eta = tWtaunutau2->p4().eta();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_eta = tWtaul2->p4().eta();
   if( tWq12 ) tree.mc_truth_tWq12_eta = tWq12->p4().eta();
   if( tWq22 ) tree.mc_truth_tWq22_eta = tWq22->p4().eta();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_eta = tWq12_IS->p4().eta();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_eta = tWq22_IS->p4().eta();
   
   if( j1 ) tree.mc_truth_j1_eta = j1->p4().eta();
   if( j2 ) tree.mc_truth_j2_eta = j2->p4().eta();
   if( j3 ) tree.mc_truth_j3_eta = j3->p4().eta();

   // phi

   if( h0 ) tree.mc_truth_h0_phi = h0->p4().phi();

   if( h0W1 ) tree.mc_truth_h0W1_phi = h0W1->p4().phi();
   if( h0W2 ) tree.mc_truth_h0W2_phi = h0W2->p4().phi();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_phi = h0Wl1->p4().phi();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_phi = h0Wnu1->p4().phi();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_phi = h0Wtau1->p4().phi();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_phi = h0Wnutau1->p4().phi();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_phi = h0Wtaul1->p4().phi();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_phi = h0Wtaunu1->p4().phi();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_phi = h0Wtaunutau1->p4().phi();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_phi = h0Wl2->p4().phi();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_phi = h0Wnu2->p4().phi();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_phi = h0Wtau2->p4().phi();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_phi = h0Wnutau2->p4().phi();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_phi = h0Wtaul2->p4().phi();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_phi = h0Wtaunu2->p4().phi();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_phi = h0Wtaunutau2->p4().phi();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_phi = h0Wq11->p4().phi();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_phi = h0Wq21->p4().phi();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_phi = h0Wq12->p4().phi();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_phi = h0Wq22->p4().phi();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_phi = h0Wq11_IS->p4().phi();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_phi = h0Wq21_IS->p4().phi();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_phi = h0Wq12_IS->p4().phi();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_phi = h0Wq22_IS->p4().phi();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_phi = h0Z1->p4().phi();
   if( h0Z2 ) tree.mc_truth_h0Z2_phi = h0Z2->p4().phi();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_phi = h0Zl11->p4().phi();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_phi = h0Zl21->p4().phi();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_phi = h0Zl12->p4().phi();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_phi = h0Zl22->p4().phi();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_phi = h0Ztau11->p4().phi();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_phi = h0Ztau21->p4().phi();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_phi = h0Ztaul11->p4().phi();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_phi = h0Ztaul21->p4().phi();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_phi = h0Ztaunu11->p4().phi();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_phi = h0Ztaunu21->p4().phi();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_phi = h0Ztaunutau11->p4().phi();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_phi = h0Ztaunutau21->p4().phi();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_phi = h0Zq11->p4().phi();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_phi = h0Zq21->p4().phi();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_phi = h0Zq12->p4().phi();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_phi = h0Zq22->p4().phi();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_phi = h0Zq11_IS->p4().phi();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_phi = h0Zq21_IS->p4().phi();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_phi = h0Zq12_IS->p4().phi();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_phi = h0Zq22_IS->p4().phi();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_phi = h0Ztau12->p4().phi();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_phi = h0Ztau22->p4().phi();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_phi = h0Ztaul12->p4().phi();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_phi = h0Ztaul22->p4().phi();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_phi = h0Ztaunu12->p4().phi();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_phi = h0Ztaunu22->p4().phi();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_phi = h0Ztaunutau12->p4().phi();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_phi = h0Ztaunutau22->p4().phi();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_phi = h0Znu11->p4().phi();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_phi = h0Znu21->p4().phi();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_phi = h0Znu12->p4().phi();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_phi = h0Znu22->p4().phi();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_phi = h0tau1->p4().phi();
   if( h0tau2 ) tree.mc_truth_h0tau2_phi = h0tau2->p4().phi();
   if( h0taul1 ) tree.mc_truth_h0taul1_phi = h0taul1->p4().phi();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_phi = h0taunutau1->p4().phi();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_phi = h0taunu1->p4().phi();
   if( h0taul2 ) tree.mc_truth_h0taul2_phi = h0taul2->p4().phi();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_phi = h0taunutau2->p4().phi();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_phi = h0taunu2->p4().phi();
   
   if( t1 ) tree.mc_truth_t1_phi = t1->p4().phi();
   if( t2 ) tree.mc_truth_t2_phi = t2->p4().phi();
   if( tb1 ) tree.mc_truth_tb1_phi = tb1->p4().phi();
   if( tb2 ) tree.mc_truth_tb2_phi = tb2->p4().phi();
   if( tb1_IS ) tree.mc_truth_tb1_IS_phi = tb1_IS->p4().phi();
   if( tb2_IS ) tree.mc_truth_tb2_IS_phi = tb2_IS->p4().phi();
   
   if( tW1 ) tree.mc_truth_tW1_phi = tW1->p4().phi();
   if( tWnu1 ) tree.mc_truth_tWnu1_phi = tWnu1->p4().phi();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_phi = tWnutau1->p4().phi();
   if( tWl1 ) tree.mc_truth_tWl1_phi = tWl1->p4().phi();
   if( tWtau1 ) tree.mc_truth_tWtau1_phi = tWtau1->p4().phi();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_phi = tWtaunu1->p4().phi();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_phi = tWtaunutau1->p4().phi();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_phi = tWtaul1->p4().phi();
   if( tWq11 ) tree.mc_truth_tWq11_phi = tWq11->p4().phi();
   if( tWq21 ) tree.mc_truth_tWq21_phi = tWq21->p4().phi();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_phi = tWq11_IS->p4().phi();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_phi = tWq21_IS->p4().phi();
   
   if( tW2 ) tree.mc_truth_tW2_phi = tW2->p4().phi();
   if( tWnu2 ) tree.mc_truth_tWnu2_phi = tWnu2->p4().phi();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_phi = tWnutau2->p4().phi();
   if( tWl2 ) tree.mc_truth_tWl2_phi = tWl2->p4().phi();
   if( tWtau2 ) tree.mc_truth_tWtau2_phi = tWtau2->p4().phi();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_phi = tWtaunu2->p4().phi();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_phi = tWtaunutau2->p4().phi();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_phi = tWtaul2->p4().phi();
   if( tWq12 ) tree.mc_truth_tWq12_phi = tWq12->p4().phi();
   if( tWq22 ) tree.mc_truth_tWq22_phi = tWq22->p4().phi();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_phi = tWq12_IS->p4().phi();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_phi = tWq22_IS->p4().phi();
   
   if( j1 ) tree.mc_truth_j1_phi = j1->p4().phi();
   if( j2 ) tree.mc_truth_j2_phi = j2->p4().phi();
   if( j3 ) tree.mc_truth_j3_phi = j3->p4().phi();

   // E

   if( h0 ) tree.mc_truth_h0_E = h0->p4().E();

   if( h0W1 ) tree.mc_truth_h0W1_E = h0W1->p4().E();
   if( h0W2 ) tree.mc_truth_h0W2_E = h0W2->p4().E();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_E = h0Wl1->p4().E();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_E = h0Wnu1->p4().E();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_E = h0Wtau1->p4().E();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_E = h0Wnutau1->p4().E();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_E = h0Wtaul1->p4().E();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_E = h0Wtaunu1->p4().E();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_E = h0Wtaunutau1->p4().E();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_E = h0Wl2->p4().E();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_E = h0Wnu2->p4().E();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_E = h0Wtau2->p4().E();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_E = h0Wnutau2->p4().E();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_E = h0Wtaul2->p4().E();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_E = h0Wtaunu2->p4().E();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_E = h0Wtaunutau2->p4().E();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_E = h0Wq11->p4().E();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_E = h0Wq21->p4().E();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_E = h0Wq12->p4().E();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_E = h0Wq22->p4().E();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_E = h0Wq11_IS->p4().E();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_E = h0Wq21_IS->p4().E();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_E = h0Wq12_IS->p4().E();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_E = h0Wq22_IS->p4().E();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_E = h0Z1->p4().E();
   if( h0Z2 ) tree.mc_truth_h0Z2_E = h0Z2->p4().E();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_E = h0Zl11->p4().E();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_E = h0Zl21->p4().E();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_E = h0Zl12->p4().E();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_E = h0Zl22->p4().E();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_E = h0Ztau11->p4().E();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_E = h0Ztau21->p4().E();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_E = h0Ztaul11->p4().E();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_E = h0Ztaul21->p4().E();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_E = h0Ztaunu11->p4().E();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_E = h0Ztaunu21->p4().E();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_E = h0Ztaunutau11->p4().E();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_E = h0Ztaunutau21->p4().E();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_E = h0Zq11->p4().E();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_E = h0Zq21->p4().E();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_E = h0Zq12->p4().E();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_E = h0Zq22->p4().E();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_E = h0Zq11_IS->p4().E();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_E = h0Zq21_IS->p4().E();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_E = h0Zq12_IS->p4().E();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_E = h0Zq22_IS->p4().E();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_E = h0Ztau12->p4().E();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_E = h0Ztau22->p4().E();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_E = h0Ztaul12->p4().E();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_E = h0Ztaul22->p4().E();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_E = h0Ztaunu12->p4().E();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_E = h0Ztaunu22->p4().E();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_E = h0Ztaunutau12->p4().E();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_E = h0Ztaunutau22->p4().E();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_E = h0Znu11->p4().E();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_E = h0Znu21->p4().E();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_E = h0Znu12->p4().E();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_E = h0Znu22->p4().E();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_E = h0tau1->p4().E();
   if( h0tau2 ) tree.mc_truth_h0tau2_E = h0tau2->p4().E();
   if( h0taul1 ) tree.mc_truth_h0taul1_E = h0taul1->p4().E();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_E = h0taunutau1->p4().E();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_E = h0taunu1->p4().E();
   if( h0taul2 ) tree.mc_truth_h0taul2_E = h0taul2->p4().E();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_E = h0taunutau2->p4().E();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_E = h0taunu2->p4().E();
   
   if( t1 ) tree.mc_truth_t1_E = t1->p4().E();
   if( t2 ) tree.mc_truth_t2_E = t2->p4().E();
   if( tb1 ) tree.mc_truth_tb1_E = tb1->p4().E();
   if( tb2 ) tree.mc_truth_tb2_E = tb2->p4().E();
   if( tb1_IS ) tree.mc_truth_tb1_IS_E = tb1_IS->p4().E();
   if( tb2_IS ) tree.mc_truth_tb2_IS_E = tb2_IS->p4().E();
   
   if( tW1 ) tree.mc_truth_tW1_E = tW1->p4().E();
   if( tWnu1 ) tree.mc_truth_tWnu1_E = tWnu1->p4().E();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_E = tWnutau1->p4().E();
   if( tWl1 ) tree.mc_truth_tWl1_E = tWl1->p4().E();
   if( tWtau1 ) tree.mc_truth_tWtau1_E = tWtau1->p4().E();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_E = tWtaunu1->p4().E();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_E = tWtaunutau1->p4().E();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_E = tWtaul1->p4().E();
   if( tWq11 ) tree.mc_truth_tWq11_E = tWq11->p4().E();
   if( tWq21 ) tree.mc_truth_tWq21_E = tWq21->p4().E();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_E = tWq11_IS->p4().E();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_E = tWq21_IS->p4().E();
   
   if( tW2 ) tree.mc_truth_tW2_E = tW2->p4().E();
   if( tWnu2 ) tree.mc_truth_tWnu2_E = tWnu2->p4().E();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_E = tWnutau2->p4().E();
   if( tWl2 ) tree.mc_truth_tWl2_E = tWl2->p4().E();
   if( tWtau2 ) tree.mc_truth_tWtau2_E = tWtau2->p4().E();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_E = tWtaunu2->p4().E();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_E = tWtaunutau2->p4().E();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_E = tWtaul2->p4().E();
   if( tWq12 ) tree.mc_truth_tWq12_E = tWq12->p4().E();
   if( tWq22 ) tree.mc_truth_tWq22_E = tWq22->p4().E();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_E = tWq12_IS->p4().E();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_E = tWq22_IS->p4().E();
   
   if( j1 ) tree.mc_truth_j1_E = j1->p4().E();
   if( j2 ) tree.mc_truth_j2_E = j2->p4().E();
   if( j3 ) tree.mc_truth_j3_E = j3->p4().E();
   
   // pdgId

   if( h0 ) tree.mc_truth_h0_id = h0->pdgId();

   if( h0W1 ) tree.mc_truth_h0W1_id = h0W1->pdgId();
   if( h0W2 ) tree.mc_truth_h0W2_id = h0W2->pdgId();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_id = h0Wl1->pdgId();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_id = h0Wnu1->pdgId();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_id = h0Wtau1->pdgId();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_id = h0Wnutau1->pdgId();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_id = h0Wtaul1->pdgId();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_id = h0Wtaunu1->pdgId();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_id = h0Wtaunutau1->pdgId();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_id = h0Wl2->pdgId();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_id = h0Wnu2->pdgId();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_id = h0Wtau2->pdgId();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_id = h0Wnutau2->pdgId();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_id = h0Wtaul2->pdgId();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_id = h0Wtaunu2->pdgId();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_id = h0Wtaunutau2->pdgId();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_id = h0Wq11->pdgId();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_id = h0Wq21->pdgId();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_id = h0Wq12->pdgId();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_id = h0Wq22->pdgId();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_id = h0Wq11_IS->pdgId();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_id = h0Wq21_IS->pdgId();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_id = h0Wq12_IS->pdgId();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_id = h0Wq22_IS->pdgId();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_id = h0Z1->pdgId();
   if( h0Z2 ) tree.mc_truth_h0Z2_id = h0Z2->pdgId();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_id = h0Zl11->pdgId();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_id = h0Zl21->pdgId();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_id = h0Zl12->pdgId();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_id = h0Zl22->pdgId();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_id = h0Ztau11->pdgId();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_id = h0Ztau21->pdgId();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_id = h0Ztaul11->pdgId();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_id = h0Ztaul21->pdgId();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_id = h0Ztaunu11->pdgId();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_id = h0Ztaunu21->pdgId();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_id = h0Ztaunutau11->pdgId();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_id = h0Ztaunutau21->pdgId();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_id = h0Zq11->pdgId();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_id = h0Zq21->pdgId();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_id = h0Zq12->pdgId();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_id = h0Zq22->pdgId();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_id = h0Zq11_IS->pdgId();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_id = h0Zq21_IS->pdgId();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_id = h0Zq12_IS->pdgId();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_id = h0Zq22_IS->pdgId();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_id = h0Ztau12->pdgId();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_id = h0Ztau22->pdgId();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_id = h0Ztaul12->pdgId();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_id = h0Ztaul22->pdgId();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_id = h0Ztaunu12->pdgId();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_id = h0Ztaunu22->pdgId();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_id = h0Ztaunutau12->pdgId();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_id = h0Ztaunutau22->pdgId();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_id = h0Znu11->pdgId();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_id = h0Znu21->pdgId();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_id = h0Znu12->pdgId();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_id = h0Znu22->pdgId();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_id = h0tau1->pdgId();
   if( h0tau2 ) tree.mc_truth_h0tau2_id = h0tau2->pdgId();
   if( h0taul1 ) tree.mc_truth_h0taul1_id = h0taul1->pdgId();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_id = h0taunutau1->pdgId();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_id = h0taunu1->pdgId();
   if( h0taul2 ) tree.mc_truth_h0taul2_id = h0taul2->pdgId();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_id = h0taunutau2->pdgId();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_id = h0taunu2->pdgId();
   
   if( t1 ) tree.mc_truth_t1_id = t1->pdgId();
   if( t2 ) tree.mc_truth_t2_id = t2->pdgId();
   if( tb1 ) tree.mc_truth_tb1_id = tb1->pdgId();
   if( tb2 ) tree.mc_truth_tb2_id = tb2->pdgId();
   if( tb1_IS ) tree.mc_truth_tb1_IS_id = tb1_IS->pdgId();
   if( tb2_IS ) tree.mc_truth_tb2_IS_id = tb2_IS->pdgId();
   
   if( tW1 ) tree.mc_truth_tW1_id = tW1->pdgId();
   if( tWnu1 ) tree.mc_truth_tWnu1_id = tWnu1->pdgId();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_id = tWnutau1->pdgId();
   if( tWl1 ) tree.mc_truth_tWl1_id = tWl1->pdgId();
   if( tWtau1 ) tree.mc_truth_tWtau1_id = tWtau1->pdgId();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_id = tWtaunu1->pdgId();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_id = tWtaunutau1->pdgId();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_id = tWtaul1->pdgId();
   if( tWq11 ) tree.mc_truth_tWq11_id = tWq11->pdgId();
   if( tWq21 ) tree.mc_truth_tWq21_id = tWq21->pdgId();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_id = tWq11_IS->pdgId();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_id = tWq21_IS->pdgId();
   
   if( tW2 ) tree.mc_truth_tW2_id = tW2->pdgId();
   if( tWnu2 ) tree.mc_truth_tWnu2_id = tWnu2->pdgId();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_id = tWnutau2->pdgId();
   if( tWl2 ) tree.mc_truth_tWl2_id = tWl2->pdgId();
   if( tWtau2 ) tree.mc_truth_tWtau2_id = tWtau2->pdgId();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_id = tWtaunu2->pdgId();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_id = tWtaunutau2->pdgId();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_id = tWtaul2->pdgId();
   if( tWq12 ) tree.mc_truth_tWq12_id = tWq12->pdgId();
   if( tWq22 ) tree.mc_truth_tWq22_id = tWq22->pdgId();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_id = tWq12_IS->pdgId();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_id = tWq22_IS->pdgId();
   
   if( j1 ) tree.mc_truth_j1_id = j1->pdgId();
   if( j2 ) tree.mc_truth_j2_id = j2->pdgId();
   if( j3 ) tree.mc_truth_j3_id = j3->pdgId();
   
   // status

   if( h0 ) tree.mc_truth_h0_status = h0->status();

   if( h0W1 ) tree.mc_truth_h0W1_status = h0W1->status();
   if( h0W2 ) tree.mc_truth_h0W2_status = h0W2->status();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_status = h0Wl1->status();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_status = h0Wnu1->status();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_status = h0Wtau1->status();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_status = h0Wnutau1->status();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_status = h0Wtaul1->status();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_status = h0Wtaunu1->status();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_status = h0Wtaunutau1->status();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_status = h0Wl2->status();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_status = h0Wnu2->status();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_status = h0Wtau2->status();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_status = h0Wnutau2->status();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_status = h0Wtaul2->status();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_status = h0Wtaunu2->status();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_status = h0Wtaunutau2->status();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_status = h0Wq11->status();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_status = h0Wq21->status();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_status = h0Wq12->status();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_status = h0Wq22->status();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_status = h0Wq11_IS->status();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_status = h0Wq21_IS->status();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_status = h0Wq12_IS->status();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_status = h0Wq22_IS->status();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_status = h0Z1->status();
   if( h0Z2 ) tree.mc_truth_h0Z2_status = h0Z2->status();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_status = h0Zl11->status();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_status = h0Zl21->status();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_status = h0Ztau11->status();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_status = h0Ztau21->status();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_status = h0Ztaul11->status();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_status = h0Ztaul21->status();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_status = h0Ztaunu11->status();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_status = h0Ztaunu21->status();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_status = h0Ztaunutau11->status();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_status = h0Ztaunutau21->status();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_status = h0Zq11->status();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_status = h0Zq21->status();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_status = h0Zq11_IS->status();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_status = h0Zq21_IS->status();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_status = h0Zl12->status();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_status = h0Zl22->status();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_status = h0Ztau12->status();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_status = h0Ztau22->status();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_status = h0Ztaul12->status();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_status = h0Ztaul22->status();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_status = h0Ztaunu12->status();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_status = h0Ztaunu22->status();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_status = h0Ztaunutau12->status();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_status = h0Ztaunutau22->status();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_status = h0Zq12->status();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_status = h0Zq22->status();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_status = h0Zq12_IS->status();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_status = h0Zq22_IS->status();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_status = h0Znu11->status();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_status = h0Znu21->status();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_status = h0Znu12->status();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_status = h0Znu22->status();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_status = h0tau1->status();
   if( h0tau2 ) tree.mc_truth_h0tau2_status = h0tau2->status();
   if( h0taul1 ) tree.mc_truth_h0taul1_status = h0taul1->status();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_status = h0taunutau1->status();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_status = h0taunu1->status();
   if( h0taul2 ) tree.mc_truth_h0taul2_status = h0taul2->status();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_status = h0taunutau2->status();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_status = h0taunu2->status();
   
   if( t1 ) tree.mc_truth_t1_status = t1->status();
   if( t2 ) tree.mc_truth_t2_status = t2->status();
   if( tb1 ) tree.mc_truth_tb1_status = tb1->status();
   if( tb2 ) tree.mc_truth_tb2_status = tb2->status();
   if( tb1_IS ) tree.mc_truth_tb1_IS_status = tb1_IS->status();
   if( tb2_IS ) tree.mc_truth_tb2_IS_status = tb2_IS->status();
   
   if( tW1 ) tree.mc_truth_tW1_status = tW1->status();
   if( tWnu1 ) tree.mc_truth_tWnu1_status = tWnu1->status();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_status = tWnutau1->status();
   if( tWl1 ) tree.mc_truth_tWl1_status = tWl1->status();
   if( tWtau1 ) tree.mc_truth_tWtau1_status = tWtau1->status();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_status = tWtaunu1->status();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_status = tWtaunutau1->status();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_status = tWtaul1->status();
   if( tWq11 ) tree.mc_truth_tWq11_status = tWq11->status();
   if( tWq21 ) tree.mc_truth_tWq21_status = tWq21->status();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_status = tWq11_IS->status();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_status = tWq21_IS->status();

   if( tW2 ) tree.mc_truth_tW2_status = tW2->status();
   if( tWnu2 ) tree.mc_truth_tWnu2_status = tWnu2->status();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_status = tWnutau2->status();
   if( tWl2 ) tree.mc_truth_tWl2_status = tWl2->status();
   if( tWtau2 ) tree.mc_truth_tWtau2_status = tWtau2->status();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_status = tWtaunu2->status();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_status = tWtaunutau2->status();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_status = tWtaul2->status();
   if( tWq12 ) tree.mc_truth_tWq12_status = tWq12->status();
   if( tWq22 ) tree.mc_truth_tWq22_status = tWq22->status();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_status = tWq12_IS->status();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_status = tWq22_IS->status();

   if( j1 ) tree.mc_truth_j1_status = j1->status();
   if( j2 ) tree.mc_truth_j2_status = j2->status();
   if( j3 ) tree.mc_truth_j3_status = j3->status();
}

void MCTruth::fillTTZSignalGenParticles(const edm::Event& iEvent,
					const edm::EventSetup& iSetup,
					FlatTree& tree,
					const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   int chan = -666;
   
   reco::GenParticle *Z = 0;
   
   reco::GenParticle *Zl1 = 0;
   reco::GenParticle *Zl2 = 0;
   reco::GenParticle *Ztau1 = 0;
   reco::GenParticle *Ztau2 = 0;
   reco::GenParticle *Ztaul1 = 0;
   reco::GenParticle *Ztaul2 = 0;
   reco::GenParticle *Ztaunu1 = 0;
   reco::GenParticle *Ztaunu2 = 0;
   reco::GenParticle *Ztaunutau1 = 0;
   reco::GenParticle *Ztaunutau2 = 0;
   reco::GenParticle *Zq1 = 0;
   reco::GenParticle *Zq2 = 0;
   reco::GenParticle *Zq1_IS = 0;
   reco::GenParticle *Zq2_IS = 0;
   reco::GenParticle *Znu1 = 0;
   reco::GenParticle *Znu2 = 0;

   reco::GenParticle *gammal1 = 0;
   reco::GenParticle *gammal2 = 0;
   reco::GenParticle *gammatau1 = 0;
   reco::GenParticle *gammatau2 = 0;
   reco::GenParticle *gammataul1 = 0;
   reco::GenParticle *gammataul2 = 0;
   reco::GenParticle *gammataunu1 = 0;
   reco::GenParticle *gammataunu2 = 0;
   reco::GenParticle *gammataunutau1 = 0;
   reco::GenParticle *gammataunutau2 = 0;
   
   reco::GenParticle *t1 = 0;
   reco::GenParticle *t2 = 0;   

   reco::GenParticle *tb1 = 0;
   reco::GenParticle *tb2 = 0;
   reco::GenParticle *tb1_IS = 0;
   reco::GenParticle *tb2_IS = 0;
   
   reco::GenParticle *tW1 = 0;
   reco::GenParticle *tWnu1 = 0;
   reco::GenParticle *tWnutau1 = 0;
   reco::GenParticle *tWl1 = 0;
   reco::GenParticle *tWtau1 = 0;
   reco::GenParticle *tWtaunu1 = 0;
   reco::GenParticle *tWtaunutau1 = 0;
   reco::GenParticle *tWtaul1 = 0;
   reco::GenParticle *tWq11 = 0;
   reco::GenParticle *tWq21 = 0;
   reco::GenParticle *tWq11_IS = 0;
   reco::GenParticle *tWq21_IS = 0;
   
   reco::GenParticle *tW2 = 0;
   reco::GenParticle *tWnu2 = 0;
   reco::GenParticle *tWnutau2 = 0;
   reco::GenParticle *tWl2 = 0;
   reco::GenParticle *tWtau2 = 0;
   reco::GenParticle *tWtaunu2 = 0;
   reco::GenParticle *tWtaunutau2 = 0;
   reco::GenParticle *tWtaul2 = 0;
   reco::GenParticle *tWq12 = 0;
   reco::GenParticle *tWq22 = 0;
   reco::GenParticle *tWq12_IS = 0;
   reco::GenParticle *tWq22_IS = 0;

   reco::GenParticle *j1 = 0;
   reco::GenParticle *j2 = 0;
   reco::GenParticle *j3 = 0;

   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;
   
   int ipart = 0;
   
    
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	int barcode = ipart; // in CMSSW barcode is the index of genParticle in the event
	// https://twiki.cern.ch/twiki/bin/view/CMS/GenParticles2HepMCConverter
	ipart++;
	
	// Additional partons (up to three)
	if( (fabs(mcp->pdgId()) <= 6 || fabs(mcp->pdgId()) == 21) &&
	    mcp->status() == 23 && barcode == 8 )
	  {
	     if( !j1 )
	       j1 = mcp;
	     else if( !j2 )
	       j2 = mcp;
	     else if( !j3 )
	       j3 = mcp;
	  }	
	  
            
	if(( (fabs(mcp->pdgId()) == 11 ||
	     fabs(mcp->pdgId()) == 13) &&
	     mcp->status() == 1
	     )      
	     ||
	    (fabs(mcp->pdgId()) == 15 &&
		mcp->status() == 2
		) 
	    )
	  {  
	      
	  //std::cout <<"lepton "<< mcp->pdgId()<<" " << mcp->status() <<" "<<  mcp->p4().pt()<< " " << barcode << std::endl;
	  //std::cout <<"lepton "<< getMother(*mcp)->pdgId() << " " << getMother(*mcp)->status() << std::endl;
	  //std::cout << "mother form mother "<<(getMother(*getMother(*mcp)))->pdgId() << std::endl;
	   
	 if  ((fabs(getMother(*mcp)->pdgId()) == 21 || (fabs(getMother(*mcp)->pdgId()) >= 1 && fabs(getMother(*mcp)->pdgId()) <= 6)) && getMother(*mcp)->status() == 21) 
	   { 
	   
	       if( fabs(mcp->pdgId()) == 11 ||
		 fabs(mcp->pdgId()) == 13 ) // l
	       {  
	       
	          //std::cout <<"youpiii "<< fabs(getMother(*mcp)->pdgId()) <<" "<< getMother(*mcp)->status() <<  std::endl;//}
		  
	          if( gammal1 && !gammal2 ) {gammal2 = mcp;}
		  if( !gammal1 ) {gammal1 = mcp;}
		  
	       }			    
	     if( fabs(mcp->pdgId()) == 15 ) // tau
	       {
		  if( gammatau1 )
		    {
		       gammatau2 = mcp;
		             
		       const reco::GenParticleRefVector& daughterRefs = gammatau2->daughterRefVector();
		       for(reco::GenParticleRefVector::const_iterator gammatau2_idr = daughterRefs.begin();
			   gammatau2_idr!= daughterRefs.end(); ++gammatau2_idr)
			 {
			    if( gammatau2_idr->isAvailable() ) 
			      {		       
				 const reco::GenParticleRef& genParticle = (*gammatau2_idr);
				 const reco::GenParticle *gammatau2_d = genParticle.get();
				 reco::GenParticle *pfff = getUnique(gammatau2_d,0);
				 
				 if( fabs(pfff->pdgId()) == 12 ||
				     fabs(pfff->pdgId()) == 14 ) // nu
				   {
				      gammataunu2 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 16 ) // nutau
				   {
				      gammataunutau2 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 11 ||
				     fabs(pfff->pdgId()) == 13 ) // l
				   {
				      gammataul2 = pfff;
				   }
			      }
			 }							  
		    }
		  if( !gammatau1 )
		    {
		       gammatau1 = mcp;
		       
		       const reco::GenParticleRefVector& daughterRefs = gammatau1->daughterRefVector();
		       for(reco::GenParticleRefVector::const_iterator gammatau1_idr = daughterRefs.begin(); 
			   gammatau1_idr!= daughterRefs.end(); ++gammatau1_idr)
			 {
			    if( gammatau1_idr->isAvailable() ) 
			      {		       
				 const reco::GenParticleRef& genParticle = (*gammatau1_idr);
				 const reco::GenParticle *gammatau1_d = genParticle.get();
				 reco::GenParticle *pfff = getUnique(gammatau1_d,0);
				 
				 if( fabs(pfff->pdgId()) == 12 ||
				     fabs(pfff->pdgId()) == 14 ) // nu
				   {
				      gammataunu1 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 16 ) // nutau
				   {
				      gammataunutau1 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 11 ||
				     fabs(pfff->pdgId()) == 13 ) // l
				   {
				      gammataul1 = pfff;
				   }
			      }
			 }							  
		    }
	       }	     
	  }		  
       }
	
	// Z
	if( fabs(mcp->pdgId()) == 23 &&
	     (mcp->status() == 62 || mcp->status() == 3) 
	    )
	  {  
	     //std::cout << "Z " <<fabs(mcp->pdgId()) << " " << mcp->status() <<" " << mcp->mass()<< std::endl;
	  
	     if( !Z ) {Z = mcp;}
	     
	     if( Z )
	       {
		  const reco::GenParticleRefVector& daughterRefs = Z->daughterRefVector();
		  for(reco::GenParticleRefVector::const_iterator Z_idr = daughterRefs.begin(); 
		      Z_idr!= daughterRefs.end(); ++Z_idr)
		    {
		       if( Z_idr->isAvailable() ) 
			 {		       
			    const reco::GenParticleRef& genParticle = (*Z_idr);
			    const reco::GenParticle *Z_d = genParticle.get();
			    reco::GenParticle *pff = getUnique(Z_d,0);
			      
			      
			    if( fabs(pff->pdgId()) == 11 ||
				fabs(pff->pdgId()) == 13 ) // l
			      {
				 if( Zl1 && !Zl2 ) {Zl2 = pff;}
				 if( !Zl1 ) {Zl1 = pff;}
			      }			    
			    if( fabs(pff->pdgId()) == 15 ) // tau
			      {
				 if( Ztau1 )
				   {
				      Ztau2 = pff;
				      
				      const reco::GenParticleRefVector& daughterRefs = Ztau2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator Ztau2_idr = daughterRefs.begin(); 
					  Ztau2_idr!= daughterRefs.end(); ++Ztau2_idr)
					{
					   if( Ztau2_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*Ztau2_idr);
						const reco::GenParticle *Ztau2_d = genParticle.get();
						reco::GenParticle *pfff = getUnique(Ztau2_d,0);
						
						if( fabs(pfff->pdgId()) == 12 ||
						    fabs(pfff->pdgId()) == 14 ) // nu
						  {
						     Ztaunu2 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 16 ) // nutau
						  {
						     Ztaunutau2 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 11 ||
						    fabs(pfff->pdgId()) == 13 ) // l
						  {
						     Ztaul2 = pfff;
						  }
					     }
					}							  
				   }
				 if( !Ztau1 )
				   {
				      Ztau1 = pff;
				      
				      const reco::GenParticleRefVector& daughterRefs = Ztau1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator Ztau1_idr = daughterRefs.begin(); 
					  Ztau1_idr!= daughterRefs.end(); ++Ztau1_idr)
					{
					   if( Ztau1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*Ztau1_idr);
						const reco::GenParticle *Ztau1_d = genParticle.get();
						reco::GenParticle *pfff = getUnique(Ztau1_d,0);
						
						if( fabs(pfff->pdgId()) == 12 ||
						    fabs(pfff->pdgId()) == 14 ) // nu
						  {
						     Ztaunu1 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 16 ) // nutau
						  {
						     Ztaunutau1 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 11 ||
						    fabs(pfff->pdgId()) == 13 ) // l
						  {
						     Ztaul1 = pfff;
						  }
					     }					   
					}				      
				   }				 
			      }			    
			    if( fabs(pff->pdgId()) <= 6 ) // q
			      {
				 if( Zq1 && !Zq2 ) Zq2 = pff;
				 if( !Zq1 ) Zq1 = pff;
				 if( Zq1_IS && !Zq2_IS ) Zq2_IS = const_cast<reco::GenParticle*>(Z_d);
				 if( !Zq1_IS ) Zq1_IS = const_cast<reco::GenParticle*>(Z_d);
			      }
			    if( fabs(pff->pdgId()) == 12 ||
				fabs(pff->pdgId()) == 14 ||
				fabs(pff->pdgId()) == 16 
			      ) // nu
			      {
				 if( Znu1 && !Znu2 ) {Znu2 = pff;}
				 if( !Znu1 ) Znu1 = pff;
			      }
			 }					   
		    }				      
	       }
	  }	
	
	
	
	// top decays
	if( fabs(mcp->pdgId()) == 6
	    && ( (mcp->status() == 62) || 
		 (mcp->status() == 3)
	       ) 
	       )
	  { 
	     //std::cout <<"top "<< mcp->status() <<" "<<mcp->pdgId() << std::endl;
	     	  
	     if( t1 && !t2 ) {t2 = const_cast<reco::GenParticle*>(mcp);}
	     if( !t1 ) {t1 = const_cast<reco::GenParticle*>(mcp);}

	     const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
	     for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
	       {
		  if( idr->isAvailable() ) 
		    {		       
		       const reco::GenParticleRef& genParticle = (*idr);
		       const reco::GenParticle *d = genParticle.get();
		       reco::GenParticle *pf = getUnique(d,0);

//		       if( pf->status() != 3 && pf->status() != 62 ) continue;
		       
		       if( fabs(pf->pdgId()) == 5 || fabs(pf->pdgId()) == 3 || fabs(pf->pdgId()) == 1 ) // b or s or d
			 {
			    if( tb1 && !tb2 ) tb2 = pf;
			    if( !tb1 ) tb1 = pf;			    
			    if( tb1_IS && !tb2_IS ) tb2_IS = const_cast<reco::GenParticle*>(d);
			    if( !tb1_IS ) tb1_IS = const_cast<reco::GenParticle*>(d);
			 }
		       
		       if( fabs(pf->pdgId()) == 24 ) // W
			 {
			    if( tW1 && !tW2 )
			      {
				 tW2 = pf;
				 const reco::GenParticleRefVector& tW2_daughterRefs = tW2->daughterRefVector();
				 for(reco::GenParticleRefVector::const_iterator tW2_idr = tW2_daughterRefs.begin();
				     tW2_idr!= tW2_daughterRefs.end(); ++tW2_idr) 
				   {
				      if( tW2_idr->isAvailable() ) 
					{		       
					   const reco::GenParticleRef& tW2_genParticle = (*tW2_idr);
					   const reco::GenParticle *tW2_d = tW2_genParticle.get();
					   reco::GenParticle *pff = getUnique(tW2_d,0);
					   
					   if( fabs(pff->pdgId()) == 12 ||
					       fabs(pff->pdgId()) == 14 ) // nu
					     {
						tWnu2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 16 ) // nu_tau
					     {
						tWnutau2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 11 ||
					       fabs(pff->pdgId()) == 13 ) // l
					     {
						tWl2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 15 ) // tau
					     {
						tWtau2 = pff;
						
						const reco::GenParticleRefVector& tWtau2_daughterRefs = tWtau2->daughterRefVector();
						for(reco::GenParticleRefVector::const_iterator tWtau2_idr = tWtau2_daughterRefs.begin();
						    tWtau2_idr!= tWtau2_daughterRefs.end(); ++tWtau2_idr) 
						  {
						     if( tWtau2_idr->isAvailable() ) 
						       {		       
							  const reco::GenParticleRef& tWtau2_genParticle = (*tWtau2_idr);
							  const reco::GenParticle *tWtau2_d = tWtau2_genParticle.get();
							  reco::GenParticle *pfff = getUnique(tWtau2_d,0);
							  
							  if( fabs(pfff->pdgId()) == 12 ||
							      fabs(pfff->pdgId()) == 14 ) // nu
							    {
							       tWtaunu2 = pfff;
							    }		
							  if( fabs(pfff->pdgId()) == 16 ) // nu_tau
							    {
							       tWtaunutau2 = pfff;
							    }			
							  if( fabs(pfff->pdgId()) == 11 ||
							      fabs(pfff->pdgId()) == 13 ) // l
							    {
							       tWtaul2 = pfff;
							    }							  
						       }
						  }						
					     }
					   if( fabs(pff->pdgId()) <= 6 ) // q
					     {
						if( tWq12 && !tWq22 ) tWq22 = pff;
						if( !tWq12 ) tWq12 = pff;
						if( tWq12_IS && !tWq22_IS ) tWq22_IS = const_cast<reco::GenParticle*>(tW2_d);
						if( !tWq12_IS ) tWq12_IS = const_cast<reco::GenParticle*>(tW2_d);
					     }					   					   
					}
				   }				
			      }
			    
			    if( !tW1 )
			      {
				 tW1 = pf;
				 const reco::GenParticleRefVector& tW1_daughterRefs = tW1->daughterRefVector();
				 for(reco::GenParticleRefVector::const_iterator tW1_idr = tW1_daughterRefs.begin();
				     tW1_idr!= tW1_daughterRefs.end(); ++tW1_idr) 
				   {
				      if( tW1_idr->isAvailable() ) 
					{		       
					   const reco::GenParticleRef& tW1_genParticle = (*tW1_idr);
					   const reco::GenParticle *tW1_d = tW1_genParticle.get();
					   reco::GenParticle *pff = getUnique(tW1_d,0);
					   
					   if( fabs(pff->pdgId()) == 12 ||
					       fabs(pff->pdgId()) == 14 ) // nu
					     {
						tWnu1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 16 ) // nu_tau
					     {
						tWnutau1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 11 ||
					       fabs(pff->pdgId()) == 13 ) // l
					     {
						tWl1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 15 ) // tau
					     {
						tWtau1 = pff;
						
						const reco::GenParticleRefVector& tWtau1_daughterRefs = tWtau1->daughterRefVector();
						for(reco::GenParticleRefVector::const_iterator tWtau1_idr = tWtau1_daughterRefs.begin();
						    tWtau1_idr!= tWtau1_daughterRefs.end(); ++tWtau1_idr) 
						  {
						     if( tWtau1_idr->isAvailable() ) 
						       {		       
							  const reco::GenParticleRef& tWtau1_genParticle = (*tWtau1_idr);
							  const reco::GenParticle *tWtau1_d = tWtau1_genParticle.get();
							  reco::GenParticle *pfff = getUnique(tWtau1_d,0);
							  
							  if( fabs(pfff->pdgId()) == 12 ||
							      fabs(pfff->pdgId()) == 14 ) // nu
							    {
							       tWtaunu1 = pfff;
							    }		
							  if( fabs(pfff->pdgId()) == 16 ) // nu_tau
							    {
							       tWtaunutau1 = pfff;
							    }			
							  if( fabs(pfff->pdgId()) == 11 ||
							      fabs(pfff->pdgId()) == 13 ) // l
							    {
							       tWtaul1 = pfff;
							    }							  
						       }
						  }						
					     }
					   if( fabs(pff->pdgId()) <= 6 ) // q
					     {
						if( tWq11 && !tWq21 ) tWq21 = pff;
						if( !tWq11 ) tWq11 = pff;
						if( tWq11_IS && !tWq21_IS ) tWq21_IS = const_cast<reco::GenParticle*>(tW1_d);
						if( !tWq11_IS ) tWq11_IS = const_cast<reco::GenParticle*>(tW1_d);
					     }
					}
				   }				
			      }			    
			 }		       
		    }		  
	       }	     
	  }
     }   

   bool doCheck = 0;
   
   if( doCheck )
     {	
	if( Zl1 || Ztaul1 || Ztaunutau1 || Zq1 || Znu1 )
	  {
	     int chan0 = 0;
	     if( Zl1 ) chan = chan0 + 0;
	     if( Ztaul1 && Ztaul2 ) chan = chan0 + 20;
	     if( Ztaunutau1 && ! Ztaul1 && Ztaunutau2 && ! Ztaul2 ) chan = chan0 + 40;
	     if( Ztaul1 && Ztaunutau2 && ! Ztaul2 ) chan = chan0 + 60;
	     if( Ztaul2 && Ztaunutau1 && ! Ztaul1 ) chan = chan0 + 80;
	     if( Zq1 && Zq2 ) chan = chan0 + 100;
	     if( Znu1 && Znu2 ) chan = chan0 + 120;
	  }	
	
	if( chan < 0 )
	  {
	     std::cout << "Unknown channel found" << std::endl;
	     std::cout << "chan = " << chan << std::endl;

	     std::cout << "j1 = " << bool(j1) << std::endl;
	     std::cout << "j2 = " << bool(j2) << std::endl;
	     std::cout << "j3 = " << bool(j3) << std::endl;

	     std::cout << "gamma->l1 = " << bool(gammal1) << std::endl;
	     std::cout << "gamma->l2 = " << bool(gammal2) << std::endl;
	     std::cout << "gamma->tau1 = " << bool(gammatau1) << std::endl;
	     std::cout << "gamma->tau1->l = " << bool(gammataul1) << std::endl;
	     std::cout << "gamma->tau1->nu = " << bool(gammataunu1) << std::endl;
	     std::cout << "gamma->tau1->nutau = " << bool(gammataunutau1) << std::endl;
	     std::cout << "gamma->tau2 = " << bool(gammatau2) << std::endl;
	     std::cout << "gamma->tau2->l = " << bool(gammataul2) << std::endl;
	     std::cout << "gamma->tau2->nu = " << bool(gammataunu2) << std::endl;
	     std::cout << "gamma->tau2->nutau = " << bool(gammataunutau2) << std::endl;
	     
	     std::cout << "Z = " << bool(Z) << std::endl;
	     std::cout << "Z->l1 = " << bool(Zl1) << std::endl;
	     std::cout << "Z->l2 = " << bool(Zl2) << std::endl;
	     std::cout << "Z->tau1 = " << bool(Ztau1) << std::endl;
	     std::cout << "Z->tau1->l = " << bool(Ztaul1) << std::endl;
	     std::cout << "Z->tau1->nu = " << bool(Ztaunu1) << std::endl;
	     std::cout << "Z->tau1->nutau = " << bool(Ztaunutau1) << std::endl;
	     std::cout << "Z->tau2 = " << bool(Ztau2) << std::endl;
	     std::cout << "Z->tau2->l = " << bool(Ztaul2) << std::endl;
	     std::cout << "Z->tau2->nu = " << bool(Ztaunu2) << std::endl;
	     std::cout << "Z->tau2->nutau = " << bool(Ztaunutau2) << std::endl;
	     std::cout << "Z->q1 = " << bool(Zq1) << std::endl;
	     std::cout << "Z->q2 = " << bool(Zq2) << std::endl;
	     std::cout << "Z->nu1 = " << bool(Znu1) << std::endl;
	     std::cout << "Z->nu2 = " << bool(Znu2) << std::endl;
	     
	     exit(1);
	  }
     }
   
   // TLV

   if( Z ) p4toTLV(Z->p4(),tree.mc_truth_Z_p4);
   if( Zl1 ) p4toTLV(Zl1->p4(),tree.mc_truth_Zl1_p4);
   if( Zl2 ) p4toTLV(Zl2->p4(),tree.mc_truth_Zl2_p4);
   if( Ztau1 ) p4toTLV(Ztau1->p4(),tree.mc_truth_Ztau1_p4);
   if( Ztau2 ) p4toTLV(Ztau2->p4(),tree.mc_truth_Ztau2_p4);
   if( Ztaul1 ) p4toTLV(Ztaul1->p4(),tree.mc_truth_Ztaul1_p4);
   if( Ztaul2 ) p4toTLV(Ztaul2->p4(),tree.mc_truth_Ztaul2_p4);
   if( Ztaunu1 ) p4toTLV(Ztaunu1->p4(),tree.mc_truth_Ztaunu1_p4);
   if( Ztaunu2 ) p4toTLV(Ztaunu2->p4(),tree.mc_truth_Ztaunu2_p4);
   if( Ztaunutau1 ) p4toTLV(Ztaunutau1->p4(),tree.mc_truth_Ztaunutau1_p4);
   if( Ztaunutau2 ) p4toTLV(Ztaunutau2->p4(),tree.mc_truth_Ztaunutau2_p4);
   if( Zq1 ) p4toTLV(Zq1->p4(),tree.mc_truth_Zq1_p4);
   if( Zq2 ) p4toTLV(Zq2->p4(),tree.mc_truth_Zq2_p4);
   if( Zq1_IS ) p4toTLV(Zq1_IS->p4(),tree.mc_truth_Zq1_IS_p4);
   if( Zq2_IS ) p4toTLV(Zq2_IS->p4(),tree.mc_truth_Zq2_IS_p4);
   if( Znu1 ) p4toTLV(Znu1->p4(),tree.mc_truth_Znu1_p4);
   if( Znu2 ) p4toTLV(Znu2->p4(),tree.mc_truth_Znu2_p4);

   if( gammal1 ) p4toTLV(gammal1->p4(),tree.mc_truth_gammal1_p4);
   if( gammal2 ) p4toTLV(gammal2->p4(),tree.mc_truth_gammal2_p4);
   if( gammatau1 ) p4toTLV(gammatau1->p4(),tree.mc_truth_gammatau1_p4);
   if( gammatau2 ) p4toTLV(gammatau2->p4(),tree.mc_truth_gammatau2_p4);
   if( gammataul1 ) p4toTLV(gammataul1->p4(),tree.mc_truth_gammataul1_p4);
   if( gammataul2 ) p4toTLV(gammataul2->p4(),tree.mc_truth_gammataul2_p4);
   if( gammataunu1 ) p4toTLV(gammataunu1->p4(),tree.mc_truth_gammataunu1_p4);
   if( gammataunu2 ) p4toTLV(gammataunu2->p4(),tree.mc_truth_gammataunu2_p4);
   if( gammataunutau1 ) p4toTLV(gammataunutau1->p4(),tree.mc_truth_gammataunutau1_p4);
   if( gammataunutau2 ) p4toTLV(gammataunutau2->p4(),tree.mc_truth_gammataunutau2_p4);
   
   if( t1 ) p4toTLV(t1->p4(),tree.mc_truth_t1_p4);
   if( t2 ) p4toTLV(t2->p4(),tree.mc_truth_t2_p4);
   if( tb1 ) p4toTLV(tb1->p4(),tree.mc_truth_tb1_p4);
   if( tb2 ) p4toTLV(tb2->p4(),tree.mc_truth_tb2_p4);
   if( tb1_IS ) p4toTLV(tb1_IS->p4(),tree.mc_truth_tb1_IS_p4);
   if( tb2_IS ) p4toTLV(tb2_IS->p4(),tree.mc_truth_tb2_IS_p4);
   
   if( tW1 ) p4toTLV(tW1->p4(),tree.mc_truth_tW1_p4);
   if( tWnu1 ) p4toTLV(tWnu1->p4(),tree.mc_truth_tWnu1_p4);
   if( tWnutau1 ) p4toTLV(tWnutau1->p4(),tree.mc_truth_tWnutau1_p4);
   if( tWl1 ) p4toTLV(tWl1->p4(),tree.mc_truth_tWl1_p4);
   if( tWtau1 ) p4toTLV(tWtau1->p4(),tree.mc_truth_tWtau1_p4);
   if( tWtaunu1 ) p4toTLV(tWtaunu1->p4(),tree.mc_truth_tWtaunu1_p4);
   if( tWtaunutau1 ) p4toTLV(tWtaunutau1->p4(),tree.mc_truth_tWtaunutau1_p4);
   if( tWtaul1 ) p4toTLV(tWtaul1->p4(),tree.mc_truth_tWtaul1_p4);
   if( tWq11 ) p4toTLV(tWq11->p4(),tree.mc_truth_tWq11_p4);
   if( tWq21 ) p4toTLV(tWq21->p4(),tree.mc_truth_tWq21_p4);
   if( tWq11_IS ) p4toTLV(tWq11_IS->p4(),tree.mc_truth_tWq11_IS_p4);
   if( tWq21_IS ) p4toTLV(tWq21_IS->p4(),tree.mc_truth_tWq21_IS_p4);

   if( tW2 ) p4toTLV(tW2->p4(),tree.mc_truth_tW2_p4);
   if( tWnu2 ) p4toTLV(tWnu2->p4(),tree.mc_truth_tWnu2_p4);
   if( tWnutau2 ) p4toTLV(tWnutau2->p4(),tree.mc_truth_tWnutau2_p4);
   if( tWl2 ) p4toTLV(tWl2->p4(),tree.mc_truth_tWl2_p4);
   if( tWtau2 ) p4toTLV(tWtau2->p4(),tree.mc_truth_tWtau2_p4);
   if( tWtaunu2 ) p4toTLV(tWtaunu2->p4(),tree.mc_truth_tWtaunu2_p4);
   if( tWtaunutau2 ) p4toTLV(tWtaunutau2->p4(),tree.mc_truth_tWtaunutau2_p4);
   if( tWtaul2 ) p4toTLV(tWtaul2->p4(),tree.mc_truth_tWtaul2_p4);
   if( tWq12 ) p4toTLV(tWq12->p4(),tree.mc_truth_tWq12_p4);
   if( tWq22 ) p4toTLV(tWq22->p4(),tree.mc_truth_tWq22_p4);
   if( tWq12_IS ) p4toTLV(tWq12_IS->p4(),tree.mc_truth_tWq12_IS_p4);
   if( tWq22_IS ) p4toTLV(tWq22_IS->p4(),tree.mc_truth_tWq22_IS_p4);

   if( j1 ) p4toTLV(j1->p4(),tree.mc_truth_j1_p4);
   if( j2 ) p4toTLV(j2->p4(),tree.mc_truth_j2_p4);
   if( j3 ) p4toTLV(j3->p4(),tree.mc_truth_j3_p4);

   // pt

   if( Z ) tree.mc_truth_Z_pt = Z->p4().pt();
   if( Zl1 ) tree.mc_truth_Zl1_pt = Zl1->p4().pt();
   if( Zl2 ) tree.mc_truth_Zl2_pt = Zl2->p4().pt();
   if( Ztau1 ) tree.mc_truth_Ztau1_pt = Ztau1->p4().pt();
   if( Ztau2 ) tree.mc_truth_Ztau2_pt = Ztau2->p4().pt();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_pt = Ztaul1->p4().pt();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_pt = Ztaul2->p4().pt();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_pt = Ztaunu1->p4().pt();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_pt = Ztaunu2->p4().pt();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_pt = Ztaunutau1->p4().pt();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_pt = Ztaunutau2->p4().pt();
   if( Zq1 ) tree.mc_truth_Zq1_pt = Zq1->p4().pt();
   if( Zq2 ) tree.mc_truth_Zq2_pt = Zq2->p4().pt();
   if( Zq1_IS ) tree.mc_truth_Zq1_IS_pt = Zq1_IS->p4().pt();
   if( Zq2_IS ) tree.mc_truth_Zq2_IS_pt = Zq2_IS->p4().pt();
   if( Znu1 ) tree.mc_truth_Znu1_pt = Znu1->p4().pt();
   if( Znu2 ) tree.mc_truth_Znu2_pt = Znu2->p4().pt();

   if( gammal1 ) tree.mc_truth_gammal1_pt = gammal1->p4().pt();
   if( gammal2 ) tree.mc_truth_gammal2_pt = gammal2->p4().pt();
   if( gammatau1 ) tree.mc_truth_gammatau1_pt = gammatau1->p4().pt();
   if( gammatau2 ) tree.mc_truth_gammatau2_pt = gammatau2->p4().pt();
   if( gammataul1 ) tree.mc_truth_gammataul1_pt = gammataul1->p4().pt();
   if( gammataul2 ) tree.mc_truth_gammataul2_pt = gammataul2->p4().pt();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_pt = gammataunu1->p4().pt();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_pt = gammataunu2->p4().pt();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_pt = gammataunutau1->p4().pt();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_pt = gammataunutau2->p4().pt();
   
   if( t1 ) tree.mc_truth_t1_pt = t1->p4().pt();
   if( t2 ) tree.mc_truth_t2_pt = t2->p4().pt();
   if( tb1 ) tree.mc_truth_tb1_pt = tb1->p4().pt();
   if( tb2 ) tree.mc_truth_tb2_pt = tb2->p4().pt();
   if( tb1_IS ) tree.mc_truth_tb1_IS_pt = tb1_IS->p4().pt();
   if( tb2_IS ) tree.mc_truth_tb2_IS_pt = tb2_IS->p4().pt();
   
   if( tW1 ) tree.mc_truth_tW1_pt = tW1->p4().pt();
   if( tWnu1 ) tree.mc_truth_tWnu1_pt = tWnu1->p4().pt();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_pt = tWnutau1->p4().pt();
   if( tWl1 ) tree.mc_truth_tWl1_pt = tWl1->p4().pt();
   if( tWtau1 ) tree.mc_truth_tWtau1_pt = tWtau1->p4().pt();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_pt = tWtaunu1->p4().pt();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_pt = tWtaunutau1->p4().pt();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_pt = tWtaul1->p4().pt();
   if( tWq11 ) tree.mc_truth_tWq11_pt = tWq11->p4().pt();
   if( tWq21 ) tree.mc_truth_tWq21_pt = tWq21->p4().pt();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_pt = tWq11_IS->p4().pt();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_pt = tWq21_IS->p4().pt();
   
   if( tW2 ) tree.mc_truth_tW2_pt = tW2->p4().pt();
   if( tWnu2 ) tree.mc_truth_tWnu2_pt = tWnu2->p4().pt();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_pt = tWnutau2->p4().pt();
   if( tWl2 ) tree.mc_truth_tWl2_pt = tWl2->p4().pt();
   if( tWtau2 ) tree.mc_truth_tWtau2_pt = tWtau2->p4().pt();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_pt = tWtaunu2->p4().pt();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_pt = tWtaunutau2->p4().pt();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_pt = tWtaul2->p4().pt();
   if( tWq12 ) tree.mc_truth_tWq12_pt = tWq12->p4().pt();
   if( tWq22 ) tree.mc_truth_tWq22_pt = tWq22->p4().pt();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_pt = tWq12_IS->p4().pt();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_pt = tWq22_IS->p4().pt();
   
   if( j1 ) tree.mc_truth_j1_pt = j1->p4().pt();
   if( j2 ) tree.mc_truth_j2_pt = j2->p4().pt();
   if( j3 ) tree.mc_truth_j3_pt = j3->p4().pt();

   // eta

   if( Z ) tree.mc_truth_Z_eta = Z->p4().eta();
   if( Zl1 ) tree.mc_truth_Zl1_eta = Zl1->p4().eta();
   if( Zl2 ) tree.mc_truth_Zl2_eta = Zl2->p4().eta();
   if( Ztau1 ) tree.mc_truth_Ztau1_eta = Ztau1->p4().eta();
   if( Ztau2 ) tree.mc_truth_Ztau2_eta = Ztau2->p4().eta();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_eta = Ztaul1->p4().eta();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_eta = Ztaul2->p4().eta();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_eta = Ztaunu1->p4().eta();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_eta = Ztaunu2->p4().eta();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_eta = Ztaunutau1->p4().eta();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_eta = Ztaunutau2->p4().eta();
   if( Zq1 ) tree.mc_truth_Zq1_eta = Zq1->p4().eta();
   if( Zq2 ) tree.mc_truth_Zq2_eta = Zq2->p4().eta();
   if( Zq1_IS ) tree.mc_truth_Zq1_IS_eta = Zq1_IS->p4().eta();
   if( Zq2_IS ) tree.mc_truth_Zq2_IS_eta = Zq2_IS->p4().eta();
   if( Znu1 ) tree.mc_truth_Znu1_eta = Znu1->p4().eta();
   if( Znu2 ) tree.mc_truth_Znu2_eta = Znu2->p4().eta();

   if( gammal1 ) tree.mc_truth_gammal1_eta = gammal1->p4().eta();
   if( gammal2 ) tree.mc_truth_gammal2_eta = gammal2->p4().eta();
   if( gammatau1 ) tree.mc_truth_gammatau1_eta = gammatau1->p4().eta();
   if( gammatau2 ) tree.mc_truth_gammatau2_eta = gammatau2->p4().eta();
   if( gammataul1 ) tree.mc_truth_gammataul1_eta = gammataul1->p4().eta();
   if( gammataul2 ) tree.mc_truth_gammataul2_eta = gammataul2->p4().eta();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_eta = gammataunu1->p4().eta();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_eta = gammataunu2->p4().eta();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_eta = gammataunutau1->p4().eta();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_eta = gammataunutau2->p4().eta();
   
   if( t1 ) tree.mc_truth_t1_eta = t1->p4().eta();
   if( t2 ) tree.mc_truth_t2_eta = t2->p4().eta();
   if( tb1 ) tree.mc_truth_tb1_eta = tb1->p4().eta();
   if( tb2 ) tree.mc_truth_tb2_eta = tb2->p4().eta();
   if( tb1_IS ) tree.mc_truth_tb1_IS_eta = tb1_IS->p4().eta();
   if( tb2_IS ) tree.mc_truth_tb2_IS_eta = tb2_IS->p4().eta();
   
   if( tW1 ) tree.mc_truth_tW1_eta = tW1->p4().eta();
   if( tWnu1 ) tree.mc_truth_tWnu1_eta = tWnu1->p4().eta();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_eta = tWnutau1->p4().eta();
   if( tWl1 ) tree.mc_truth_tWl1_eta = tWl1->p4().eta();
   if( tWtau1 ) tree.mc_truth_tWtau1_eta = tWtau1->p4().eta();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_eta = tWtaunu1->p4().eta();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_eta = tWtaunutau1->p4().eta();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_eta = tWtaul1->p4().eta();
   if( tWq11 ) tree.mc_truth_tWq11_eta = tWq11->p4().eta();
   if( tWq21 ) tree.mc_truth_tWq21_eta = tWq21->p4().eta();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_eta = tWq11_IS->p4().eta();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_eta = tWq21_IS->p4().eta();
   
   if( tW2 ) tree.mc_truth_tW2_eta = tW2->p4().eta();
   if( tWnu2 ) tree.mc_truth_tWnu2_eta = tWnu2->p4().eta();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_eta = tWnutau2->p4().eta();
   if( tWl2 ) tree.mc_truth_tWl2_eta = tWl2->p4().eta();
   if( tWtau2 ) tree.mc_truth_tWtau2_eta = tWtau2->p4().eta();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_eta = tWtaunu2->p4().eta();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_eta = tWtaunutau2->p4().eta();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_eta = tWtaul2->p4().eta();
   if( tWq12 ) tree.mc_truth_tWq12_eta = tWq12->p4().eta();
   if( tWq22 ) tree.mc_truth_tWq22_eta = tWq22->p4().eta();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_eta = tWq12_IS->p4().eta();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_eta = tWq22_IS->p4().eta();
   
   if( j1 ) tree.mc_truth_j1_eta = j1->p4().eta();
   if( j2 ) tree.mc_truth_j2_eta = j2->p4().eta();
   if( j3 ) tree.mc_truth_j3_eta = j3->p4().eta();

   // phi

   if( Z ) tree.mc_truth_Z_phi = Z->p4().phi();
   if( Zl1 ) tree.mc_truth_Zl1_phi = Zl1->p4().phi();
   if( Zl2 ) tree.mc_truth_Zl2_phi = Zl2->p4().phi();
   if( Ztau1 ) tree.mc_truth_Ztau1_phi = Ztau1->p4().phi();
   if( Ztau2 ) tree.mc_truth_Ztau2_phi = Ztau2->p4().phi();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_phi = Ztaul1->p4().phi();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_phi = Ztaul2->p4().phi();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_phi = Ztaunu1->p4().phi();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_phi = Ztaunu2->p4().phi();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_phi = Ztaunutau1->p4().phi();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_phi = Ztaunutau2->p4().phi();
   if( Zq1 ) tree.mc_truth_Zq1_phi = Zq1->p4().phi();
   if( Zq2 ) tree.mc_truth_Zq2_phi = Zq2->p4().phi();
   if( Zq1_IS ) tree.mc_truth_Zq1_IS_phi = Zq1_IS->p4().phi();
   if( Zq2_IS ) tree.mc_truth_Zq2_IS_phi = Zq2_IS->p4().phi();
   if( Znu1 ) tree.mc_truth_Znu1_phi = Znu1->p4().phi();
   if( Znu2 ) tree.mc_truth_Znu2_phi = Znu2->p4().phi();

   if( gammal1 ) tree.mc_truth_gammal1_phi = gammal1->p4().phi();
   if( gammal2 ) tree.mc_truth_gammal2_phi = gammal2->p4().phi();
   if( gammatau1 ) tree.mc_truth_gammatau1_phi = gammatau1->p4().phi();
   if( gammatau2 ) tree.mc_truth_gammatau2_phi = gammatau2->p4().phi();
   if( gammataul1 ) tree.mc_truth_gammataul1_phi = gammataul1->p4().phi();
   if( gammataul2 ) tree.mc_truth_gammataul2_phi = gammataul2->p4().phi();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_phi = gammataunu1->p4().phi();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_phi = gammataunu2->p4().phi();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_phi = gammataunutau1->p4().phi();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_phi = gammataunutau2->p4().phi();
   
   if( t1 ) tree.mc_truth_t1_phi = t1->p4().phi();
   if( t2 ) tree.mc_truth_t2_phi = t2->p4().phi();
   if( tb1 ) tree.mc_truth_tb1_phi = tb1->p4().phi();
   if( tb2 ) tree.mc_truth_tb2_phi = tb2->p4().phi();
   if( tb1_IS ) tree.mc_truth_tb1_IS_phi = tb1_IS->p4().phi();
   if( tb2_IS ) tree.mc_truth_tb2_IS_phi = tb2_IS->p4().phi();
   
   if( tW1 ) tree.mc_truth_tW1_phi = tW1->p4().phi();
   if( tWnu1 ) tree.mc_truth_tWnu1_phi = tWnu1->p4().phi();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_phi = tWnutau1->p4().phi();
   if( tWl1 ) tree.mc_truth_tWl1_phi = tWl1->p4().phi();
   if( tWtau1 ) tree.mc_truth_tWtau1_phi = tWtau1->p4().phi();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_phi = tWtaunu1->p4().phi();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_phi = tWtaunutau1->p4().phi();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_phi = tWtaul1->p4().phi();
   if( tWq11 ) tree.mc_truth_tWq11_phi = tWq11->p4().phi();
   if( tWq21 ) tree.mc_truth_tWq21_phi = tWq21->p4().phi();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_phi = tWq11_IS->p4().phi();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_phi = tWq21_IS->p4().phi();
   
   if( tW2 ) tree.mc_truth_tW2_phi = tW2->p4().phi();
   if( tWnu2 ) tree.mc_truth_tWnu2_phi = tWnu2->p4().phi();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_phi = tWnutau2->p4().phi();
   if( tWl2 ) tree.mc_truth_tWl2_phi = tWl2->p4().phi();
   if( tWtau2 ) tree.mc_truth_tWtau2_phi = tWtau2->p4().phi();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_phi = tWtaunu2->p4().phi();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_phi = tWtaunutau2->p4().phi();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_phi = tWtaul2->p4().phi();
   if( tWq12 ) tree.mc_truth_tWq12_phi = tWq12->p4().phi();
   if( tWq22 ) tree.mc_truth_tWq22_phi = tWq22->p4().phi();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_phi = tWq12_IS->p4().phi();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_phi = tWq22_IS->p4().phi();
   
   if( j1 ) tree.mc_truth_j1_phi = j1->p4().phi();
   if( j2 ) tree.mc_truth_j2_phi = j2->p4().phi();
   if( j3 ) tree.mc_truth_j3_phi = j3->p4().phi();

   // E

   if( Z ) tree.mc_truth_Z_E = Z->p4().E();
   if( Zl1 ) tree.mc_truth_Zl1_E = Zl1->p4().E();
   if( Zl2 ) tree.mc_truth_Zl2_E = Zl2->p4().E();
   if( Ztau1 ) tree.mc_truth_Ztau1_E = Ztau1->p4().E();
   if( Ztau2 ) tree.mc_truth_Ztau2_E = Ztau2->p4().E();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_E = Ztaul1->p4().E();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_E = Ztaul2->p4().E();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_E = Ztaunu1->p4().E();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_E = Ztaunu2->p4().E();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_E = Ztaunutau1->p4().E();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_E = Ztaunutau2->p4().E();
   if( Zq1 ) tree.mc_truth_Zq1_E = Zq1->p4().E();
   if( Zq2 ) tree.mc_truth_Zq2_E = Zq2->p4().E();
   if( Zq1_IS ) tree.mc_truth_Zq1_IS_E = Zq1_IS->p4().E();
   if( Zq2_IS ) tree.mc_truth_Zq2_IS_E = Zq2_IS->p4().E();
   if( Znu1 ) tree.mc_truth_Znu1_E = Znu1->p4().E();
   if( Znu2 ) tree.mc_truth_Znu2_E = Znu2->p4().E();

   if( gammal1 ) tree.mc_truth_gammal1_E = gammal1->p4().E();
   if( gammal2 ) tree.mc_truth_gammal2_E = gammal2->p4().E();
   if( gammatau1 ) tree.mc_truth_gammatau1_E = gammatau1->p4().E();
   if( gammatau2 ) tree.mc_truth_gammatau2_E = gammatau2->p4().E();
   if( gammataul1 ) tree.mc_truth_gammataul1_E = gammataul1->p4().E();
   if( gammataul2 ) tree.mc_truth_gammataul2_E = gammataul2->p4().E();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_E = gammataunu1->p4().E();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_E = gammataunu2->p4().E();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_E = gammataunutau1->p4().E();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_E = gammataunutau2->p4().E();
   
   if( t1 ) tree.mc_truth_t1_E = t1->p4().E();
   if( t2 ) tree.mc_truth_t2_E = t2->p4().E();
   if( tb1 ) tree.mc_truth_tb1_E = tb1->p4().E();
   if( tb2 ) tree.mc_truth_tb2_E = tb2->p4().E();
   if( tb1_IS ) tree.mc_truth_tb1_IS_E = tb1_IS->p4().E();
   if( tb2_IS ) tree.mc_truth_tb2_IS_E = tb2_IS->p4().E();
   
   if( tW1 ) tree.mc_truth_tW1_E = tW1->p4().E();
   if( tWnu1 ) tree.mc_truth_tWnu1_E = tWnu1->p4().E();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_E = tWnutau1->p4().E();
   if( tWl1 ) tree.mc_truth_tWl1_E = tWl1->p4().E();
   if( tWtau1 ) tree.mc_truth_tWtau1_E = tWtau1->p4().E();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_E = tWtaunu1->p4().E();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_E = tWtaunutau1->p4().E();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_E = tWtaul1->p4().E();
   if( tWq11 ) tree.mc_truth_tWq11_E = tWq11->p4().E();
   if( tWq21 ) tree.mc_truth_tWq21_E = tWq21->p4().E();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_E = tWq11_IS->p4().E();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_E = tWq21_IS->p4().E();
   
   if( tW2 ) tree.mc_truth_tW2_E = tW2->p4().E();
   if( tWnu2 ) tree.mc_truth_tWnu2_E = tWnu2->p4().E();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_E = tWnutau2->p4().E();
   if( tWl2 ) tree.mc_truth_tWl2_E = tWl2->p4().E();
   if( tWtau2 ) tree.mc_truth_tWtau2_E = tWtau2->p4().E();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_E = tWtaunu2->p4().E();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_E = tWtaunutau2->p4().E();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_E = tWtaul2->p4().E();
   if( tWq12 ) tree.mc_truth_tWq12_E = tWq12->p4().E();
   if( tWq22 ) tree.mc_truth_tWq22_E = tWq22->p4().E();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_E = tWq12_IS->p4().E();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_E = tWq22_IS->p4().E();
   
   if( j1 ) tree.mc_truth_j1_E = j1->p4().E();
   if( j2 ) tree.mc_truth_j2_E = j2->p4().E();
   if( j3 ) tree.mc_truth_j3_E = j3->p4().E();
   
   // pdgId

   if( Z ) tree.mc_truth_Z_id = Z->pdgId();
   if( Zl1 ) tree.mc_truth_Zl1_id = Zl1->pdgId();
   if( Zl2 ) tree.mc_truth_Zl2_id = Zl2->pdgId();
   if( Ztau1 ) tree.mc_truth_Ztau1_id = Ztau1->pdgId();
   if( Ztau2 ) tree.mc_truth_Ztau2_id = Ztau2->pdgId();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_id = Ztaul1->pdgId();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_id = Ztaul2->pdgId();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_id = Ztaunu1->pdgId();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_id = Ztaunu2->pdgId();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_id = Ztaunutau1->pdgId();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_id = Ztaunutau2->pdgId();
   if( Zq1 ) tree.mc_truth_Zq1_id = Zq1->pdgId();
   if( Zq2 ) tree.mc_truth_Zq2_id = Zq2->pdgId();
   if( Zq1_IS ) tree.mc_truth_Zq1_IS_id = Zq1_IS->pdgId();
   if( Zq2_IS ) tree.mc_truth_Zq2_IS_id = Zq2_IS->pdgId();
   if( Znu1 ) tree.mc_truth_Znu1_id = Znu1->pdgId();
   if( Znu2 ) tree.mc_truth_Znu2_id = Znu2->pdgId();

   if( gammal1 ) tree.mc_truth_gammal1_id = gammal1->pdgId();
   if( gammal2 ) tree.mc_truth_gammal2_id = gammal2->pdgId();
   if( gammatau1 ) tree.mc_truth_gammatau1_id = gammatau1->pdgId();
   if( gammatau2 ) tree.mc_truth_gammatau2_id = gammatau2->pdgId();
   if( gammataul1 ) tree.mc_truth_gammataul1_id = gammataul1->pdgId();
   if( gammataul2 ) tree.mc_truth_gammataul2_id = gammataul2->pdgId();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_id = gammataunu1->pdgId();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_id = gammataunu2->pdgId();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_id = gammataunutau1->pdgId();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_id = gammataunutau2->pdgId();
   
   if( t1 ) tree.mc_truth_t1_id = t1->pdgId();
   if( t2 ) tree.mc_truth_t2_id = t2->pdgId();
   if( tb1 ) tree.mc_truth_tb1_id = tb1->pdgId();
   if( tb2 ) tree.mc_truth_tb2_id = tb2->pdgId();
   if( tb1_IS ) tree.mc_truth_tb1_IS_id = tb1_IS->pdgId();
   if( tb2_IS ) tree.mc_truth_tb2_IS_id = tb2_IS->pdgId();
   
   if( tW1 ) tree.mc_truth_tW1_id = tW1->pdgId();
   if( tWnu1 ) tree.mc_truth_tWnu1_id = tWnu1->pdgId();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_id = tWnutau1->pdgId();
   if( tWl1 ) tree.mc_truth_tWl1_id = tWl1->pdgId();
   if( tWtau1 ) tree.mc_truth_tWtau1_id = tWtau1->pdgId();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_id = tWtaunu1->pdgId();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_id = tWtaunutau1->pdgId();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_id = tWtaul1->pdgId();
   if( tWq11 ) tree.mc_truth_tWq11_id = tWq11->pdgId();
   if( tWq21 ) tree.mc_truth_tWq21_id = tWq21->pdgId();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_id = tWq11_IS->pdgId();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_id = tWq21_IS->pdgId();
   
   if( tW2 ) tree.mc_truth_tW2_id = tW2->pdgId();
   if( tWnu2 ) tree.mc_truth_tWnu2_id = tWnu2->pdgId();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_id = tWnutau2->pdgId();
   if( tWl2 ) tree.mc_truth_tWl2_id = tWl2->pdgId();
   if( tWtau2 ) tree.mc_truth_tWtau2_id = tWtau2->pdgId();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_id = tWtaunu2->pdgId();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_id = tWtaunutau2->pdgId();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_id = tWtaul2->pdgId();
   if( tWq12 ) tree.mc_truth_tWq12_id = tWq12->pdgId();
   if( tWq22 ) tree.mc_truth_tWq22_id = tWq22->pdgId();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_id = tWq12_IS->pdgId();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_id = tWq22_IS->pdgId();
   
   if( j1 ) tree.mc_truth_j1_id = j1->pdgId();
   if( j2 ) tree.mc_truth_j2_id = j2->pdgId();
   if( j3 ) tree.mc_truth_j3_id = j3->pdgId();
   
   // status

   if( Z ) tree.mc_truth_Z_status = Z->status();
   if( Zl1 ) tree.mc_truth_Zl1_status = Zl1->status();
   if( Zl2 ) tree.mc_truth_Zl2_status = Zl2->status();
   if( Ztau1 ) tree.mc_truth_Ztau1_status = Ztau1->status();
   if( Ztau2 ) tree.mc_truth_Ztau2_status = Ztau2->status();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_status = Ztaul1->status();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_status = Ztaul2->status();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_status = Ztaunu1->status();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_status = Ztaunu2->status();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_status = Ztaunutau1->status();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_status = Ztaunutau2->status();
   if( Zq1 ) tree.mc_truth_Zq1_status = Zq1->status();
   if( Zq2 ) tree.mc_truth_Zq2_status = Zq2->status();
   if( Zq1_IS ) tree.mc_truth_Zq1_IS_status = Zq1_IS->status();
   if( Zq2_IS ) tree.mc_truth_Zq2_IS_status = Zq2_IS->status();
   if( Znu1 ) tree.mc_truth_Znu1_status = Znu1->status();
   if( Znu2 ) tree.mc_truth_Znu2_status = Znu2->status();
   
   if( gammal1 ) tree.mc_truth_gammal1_status = gammal1->status();
   if( gammal2 ) tree.mc_truth_gammal2_status = gammal2->status();
   if( gammatau1 ) tree.mc_truth_gammatau1_status = gammatau1->status();
   if( gammatau2 ) tree.mc_truth_gammatau2_status = gammatau2->status();
   if( gammataul1 ) tree.mc_truth_gammataul1_status = gammataul1->status();
   if( gammataul2 ) tree.mc_truth_gammataul2_status = gammataul2->status();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_status = gammataunu1->status();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_status = gammataunu2->status();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_status = gammataunutau1->status();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_status = gammataunutau2->status();
   
   if( t1 ) tree.mc_truth_t1_status = t1->status();
   if( t2 ) tree.mc_truth_t2_status = t2->status();
   if( tb1 ) tree.mc_truth_tb1_status = tb1->status();
   if( tb2 ) tree.mc_truth_tb2_status = tb2->status();
   if( tb1_IS ) tree.mc_truth_tb1_IS_status = tb1_IS->status();
   if( tb2_IS ) tree.mc_truth_tb2_IS_status = tb2_IS->status();
   
   if( tW1 ) tree.mc_truth_tW1_status = tW1->status();
   if( tWnu1 ) tree.mc_truth_tWnu1_status = tWnu1->status();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_status = tWnutau1->status();
   if( tWl1 ) tree.mc_truth_tWl1_status = tWl1->status();
   if( tWtau1 ) tree.mc_truth_tWtau1_status = tWtau1->status();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_status = tWtaunu1->status();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_status = tWtaunutau1->status();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_status = tWtaul1->status();
   if( tWq11 ) tree.mc_truth_tWq11_status = tWq11->status();
   if( tWq21 ) tree.mc_truth_tWq21_status = tWq21->status();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_status = tWq11_IS->status();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_status = tWq21_IS->status();

   if( tW2 ) tree.mc_truth_tW2_status = tW2->status();
   if( tWnu2 ) tree.mc_truth_tWnu2_status = tWnu2->status();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_status = tWnutau2->status();
   if( tWl2 ) tree.mc_truth_tWl2_status = tWl2->status();
   if( tWtau2 ) tree.mc_truth_tWtau2_status = tWtau2->status();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_status = tWtaunu2->status();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_status = tWtaunutau2->status();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_status = tWtaul2->status();
   if( tWq12 ) tree.mc_truth_tWq12_status = tWq12->status();
   if( tWq22 ) tree.mc_truth_tWq22_status = tWq22->status();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_status = tWq12_IS->status();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_status = tWq22_IS->status();

   if( j1 ) tree.mc_truth_j1_status = j1->status();
   if( j2 ) tree.mc_truth_j2_status = j2->status();
   if( j3 ) tree.mc_truth_j3_status = j3->status();
   	
}

void MCTruth::fillTTWSignalGenParticles(const edm::Event& iEvent,
					const edm::EventSetup& iSetup,
					FlatTree& tree,
					const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   int chan = -666;
   
   reco::GenParticle *W = 0;
   reco::GenParticle *Wnu = 0;
   reco::GenParticle *Wnutau = 0;
   reco::GenParticle *Wl = 0;
   reco::GenParticle *Wtau = 0;
   reco::GenParticle *Wtaunu = 0;
   reco::GenParticle *Wtaunutau = 0;
   reco::GenParticle *Wtaul = 0;
   reco::GenParticle *Wq1 = 0;
   reco::GenParticle *Wq2 = 0;
   reco::GenParticle *Wq1_IS = 0;
   reco::GenParticle *Wq2_IS = 0;

   /*
   reco::GenParticle *gammal1 = 0;
   reco::GenParticle *gammal2 = 0;
   reco::GenParticle *gammatau1 = 0;
   reco::GenParticle *gammatau2 = 0;
   reco::GenParticle *gammataul1 = 0;
   reco::GenParticle *gammataul2 = 0;
   reco::GenParticle *gammataunu1 = 0;
   reco::GenParticle *gammataunu2 = 0;
   reco::GenParticle *gammataunutau1 = 0;
   reco::GenParticle *gammataunutau2 = 0;*/

   reco::GenParticle *t1 = 0;
   reco::GenParticle *t2 = 0;   

   reco::GenParticle *tb1 = 0;
   reco::GenParticle *tb2 = 0;
   reco::GenParticle *tb1_IS = 0;
   reco::GenParticle *tb2_IS = 0;
   
   reco::GenParticle *tW1 = 0;
   reco::GenParticle *tWnu1 = 0;
   reco::GenParticle *tWnutau1 = 0;
   reco::GenParticle *tWl1 = 0;
   reco::GenParticle *tWtau1 = 0;
   reco::GenParticle *tWtaunu1 = 0;
   reco::GenParticle *tWtaunutau1 = 0;
   reco::GenParticle *tWtaul1 = 0;
   reco::GenParticle *tWq11 = 0;
   reco::GenParticle *tWq21 = 0;
   reco::GenParticle *tWq11_IS = 0;
   reco::GenParticle *tWq21_IS = 0;
   
   reco::GenParticle *tW2 = 0;
   reco::GenParticle *tWnu2 = 0;
   reco::GenParticle *tWnutau2 = 0;
   reco::GenParticle *tWl2 = 0;
   reco::GenParticle *tWtau2 = 0;
   reco::GenParticle *tWtaunu2 = 0;
   reco::GenParticle *tWtaunutau2 = 0;
   reco::GenParticle *tWtaul2 = 0;
   reco::GenParticle *tWq12 = 0;
   reco::GenParticle *tWq22 = 0;
   reco::GenParticle *tWq12_IS = 0;
   reco::GenParticle *tWq22_IS = 0;

   reco::GenParticle *j1 = 0;
   reco::GenParticle *j2 = 0;
   reco::GenParticle *j3 = 0;

   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;
   
   int ipart = 0;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	int barcode = ipart; // in CMSSW barcode is the index of genParticle in the event
	// https://twiki.cern.ch/twiki/bin/view/CMS/GenParticles2HepMCConverter
	ipart++;
	
	// Additional partons (up to three)
	if( (fabs(mcp->pdgId()) <= 6 || fabs(mcp->pdgId()) == 21) &&
	    mcp->status() == 23 && barcode == 8 )
	  {
	     if( !j1 )
	       j1 = mcp;
	     else if( !j2 )
	       j2 = mcp;
	     else if( !j3 )
	       j3 = mcp;
	  }	

	/*if( ((fabs(mcp->pdgId()) == 11 ||
	      fabs(mcp->pdgId()) == 13) &&
	     mcp->status() == 3) ||
	    (fabs(mcp->pdgId()) == 15 &&
		mcp->status() == 2) )
	  {
	     if( fabs(mcp->pdgId()) == 11 ||
		 fabs(mcp->pdgId()) == 13 ) // l
	       {
		  if( gammal1 && !gammal2 ) {gammal2 = mcp;}
		  if( !gammal1 ) {gammal1 = mcp;}
	       }			    
	     if( fabs(mcp->pdgId()) == 15 ) // tau
	       {
		  if( gammatau1 )
		    {
		       gammatau2 = mcp;
				      
		       const reco::GenParticleRefVector& daughterRefs = gammatau2->daughterRefVector();
		       for(reco::GenParticleRefVector::const_iterator gammatau2_idr = daughterRefs.begin(); 
			   gammatau2_idr!= daughterRefs.end(); ++gammatau2_idr)
			 {
			    if( gammatau2_idr->isAvailable() ) 
			      {		       
				 const reco::GenParticleRef& genParticle = (*gammatau2_idr);
				 const reco::GenParticle *gammatau2_d = genParticle.get();
				 reco::GenParticle *pfff = getUnique(gammatau2_d,0);
				 
				 if( fabs(pfff->pdgId()) == 12 ||
				     fabs(pfff->pdgId()) == 14 ) // nu
				   {
				      gammataunu2 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 16 ) // nutau
				   {
				      gammataunutau2 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 11 ||
				     fabs(pfff->pdgId()) == 13 ) // l
				   {
				      gammataul2 = pfff;
				   }
			      }
			 }							  
		    }
		  if( !gammatau1 )
		    {
		       gammatau1 = mcp;
		       
		       const reco::GenParticleRefVector& daughterRefs = gammatau1->daughterRefVector();
		       for(reco::GenParticleRefVector::const_iterator gammatau1_idr = daughterRefs.begin(); 
			   gammatau1_idr!= daughterRefs.end(); ++gammatau1_idr)
			 {
			    if( gammatau1_idr->isAvailable() ) 
			      {		       
				 const reco::GenParticleRef& genParticle = (*gammatau1_idr);
				 const reco::GenParticle *gammatau1_d = genParticle.get();
				 reco::GenParticle *pfff = getUnique(gammatau1_d,0);
				 
				 if( fabs(pfff->pdgId()) == 12 ||
				     fabs(pfff->pdgId()) == 14 ) // nu
				   {
				      gammataunu1 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 16 ) // nutau
				   {
				      gammataunutau1 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 11 ||
				     fabs(pfff->pdgId()) == 13 ) // l
				   {
				      gammataul1 = pfff;
				   }
			      }
			 }							  
		    }
	       }	     
	  }*/		  
	
	// W
	if( fabs(mcp->pdgId()) == 24 &&
	    (mcp->status() == 62 || mcp->status() == 3) )
	  {
	     if( !W ) {W = mcp;}
	     
	     if( W )
	       {
		  const reco::GenParticleRefVector& daughterRefs = W->daughterRefVector();
		  for(reco::GenParticleRefVector::const_iterator W_idr = daughterRefs.begin(); 
		      W_idr!= daughterRefs.end(); ++W_idr)
		    {
		       if( W_idr->isAvailable() )
			 {		       
			    const reco::GenParticleRef& genParticle = (*W_idr);
			    const reco::GenParticle *W_d = genParticle.get();
			    reco::GenParticle *pff = getUnique(W_d,0);
			    
			    if( fabs(pff->pdgId()) == 12 ||
				fabs(pff->pdgId()) == 14 ) // nu
			      {
				 Wnu = pff;
			      }		
			    if( fabs(pff->pdgId()) == 16 ) // nu_tau
			      {
				 Wnutau = pff;
			      }
			    if( fabs(pff->pdgId()) == 11 ||
				fabs(pff->pdgId()) == 13 ) // l
			      {
				 Wl = pff;
			      }		
			    if( fabs(pff->pdgId()) == 15 ) // tau
			      {
				 Wtau = pff;
				 
				 const reco::GenParticleRefVector& Wtau_daughterRefs = Wtau->daughterRefVector();
				 for(reco::GenParticleRefVector::const_iterator Wtau_idr = Wtau_daughterRefs.begin();
				     Wtau_idr!= Wtau_daughterRefs.end(); ++Wtau_idr) 
				   {
				      if( Wtau_idr->isAvailable() )
					{		       
					   const reco::GenParticleRef& Wtau_genParticle = (*Wtau_idr);
					   const reco::GenParticle *Wtau_d = Wtau_genParticle.get();
					   reco::GenParticle *pfff = getUnique(Wtau_d,0);
					   
					   if( fabs(pfff->pdgId()) == 12 ||
					       fabs(pfff->pdgId()) == 14 ) // nu
					     {
						Wtaunu = pfff;
					     }		
					   if( fabs(pfff->pdgId()) == 16 ) // nu_tau
					     {
						Wtaunutau = pfff;
					     }			
					   if( fabs(pfff->pdgId()) == 11 ||
					       fabs(pfff->pdgId()) == 13 ) // l
					     {
						Wtaul = pfff;
					     }							  
					}
				   }						
			      }
			    if( fabs(pff->pdgId()) <= 6 ) // q
			      {
				 if( Wq1 && !Wq2 ) Wq2 = pff;
				 if( !Wq1 ) Wq1 = pff;
				 if( Wq1_IS && !Wq2_IS ) Wq2_IS = const_cast<reco::GenParticle*>(W_d);
				 if( !Wq1_IS ) Wq1_IS = const_cast<reco::GenParticle*>(W_d);
			      }					   					   
			 }
		    }
	       }
	  }	

	// top decays
	if( fabs(mcp->pdgId()) == 6
	    && ( (mcp->status() == 62) || 
		 (mcp->status() == 3)
	       ) )
	  {
	     if( t1 && !t2 ) {t2 = const_cast<reco::GenParticle*>(mcp);}
	     if( !t1 ) {t1 = const_cast<reco::GenParticle*>(mcp);}

	     const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
	     for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
	       {
		  if( idr->isAvailable() ) 
		    {		       
		       const reco::GenParticleRef& genParticle = (*idr);
		       const reco::GenParticle *d = genParticle.get();
		       reco::GenParticle *pf = getUnique(d,0);

//		       if( pf->status() != 3 && pf->status() != 62 ) continue;
		       
		       if( fabs(pf->pdgId()) == 5 || fabs(pf->pdgId()) == 3 || fabs(pf->pdgId()) == 1 ) // b or s or d
			 {
			    if( tb1 && !tb2 ) tb2 = pf;
			    if( !tb1 ) tb1 = pf;			    
			    if( tb1_IS && !tb2_IS ) tb2_IS = const_cast<reco::GenParticle*>(d);
			    if( !tb1_IS ) tb1_IS = const_cast<reco::GenParticle*>(d);
			 }		       
		       
		       if( fabs(pf->pdgId()) == 24 ) // W
			 {
			    if( tW1 && !tW2 )
			      {
				 tW2 = pf;
				 const reco::GenParticleRefVector& tW2_daughterRefs = tW2->daughterRefVector();
				 for(reco::GenParticleRefVector::const_iterator tW2_idr = tW2_daughterRefs.begin();
				     tW2_idr!= tW2_daughterRefs.end(); ++tW2_idr) 
				   {
				      if( tW2_idr->isAvailable() ) 
					{		       
					   const reco::GenParticleRef& tW2_genParticle = (*tW2_idr);
					   const reco::GenParticle *tW2_d = tW2_genParticle.get();
					   reco::GenParticle *pff = getUnique(tW2_d,0);
					   
					   if( fabs(pff->pdgId()) == 12 ||
					       fabs(pff->pdgId()) == 14 ) // nu
					     {
						tWnu2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 16 ) // nu_tau
					     {
						tWnutau2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 11 ||
					       fabs(pff->pdgId()) == 13 ) // l
					     {
						tWl2 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 15 ) // tau
					     {
						tWtau2 = pff;
						
						const reco::GenParticleRefVector& tWtau2_daughterRefs = tWtau2->daughterRefVector();
						for(reco::GenParticleRefVector::const_iterator tWtau2_idr = tWtau2_daughterRefs.begin();
						    tWtau2_idr!= tWtau2_daughterRefs.end(); ++tWtau2_idr) 
						  {
						     if( tWtau2_idr->isAvailable() ) 
						       {		       
							  const reco::GenParticleRef& tWtau2_genParticle = (*tWtau2_idr);
							  const reco::GenParticle *tWtau2_d = tWtau2_genParticle.get();
							  reco::GenParticle *pfff = getUnique(tWtau2_d,0);
							  
							  if( fabs(pfff->pdgId()) == 12 ||
							      fabs(pfff->pdgId()) == 14 ) // nu
							    {
							       tWtaunu2 = pfff;
							    }		
							  if( fabs(pfff->pdgId()) == 16 ) // nu_tau
							    {
							       tWtaunutau2 = pfff;
							    }			
							  if( fabs(pfff->pdgId()) == 11 ||
							      fabs(pfff->pdgId()) == 13 ) // l
							    {
							       tWtaul2 = pfff;
							    }							  
						       }
						  }						
					     }
					   if( fabs(pff->pdgId()) <= 6 ) // q
					     {
						if( tWq12 && !tWq22 ) tWq22 = pff;
						if( !tWq12 ) tWq12 = pff;
						if( tWq12_IS && !tWq22_IS ) tWq22_IS = const_cast<reco::GenParticle*>(tW2_d);
						if( !tWq12_IS ) tWq12_IS = const_cast<reco::GenParticle*>(tW2_d);
					     }					   					   
					}
				   }				
			      }
			    
			    if( !tW1 )
			      {
				 tW1 = pf;
				 const reco::GenParticleRefVector& tW1_daughterRefs = tW1->daughterRefVector();
				 for(reco::GenParticleRefVector::const_iterator tW1_idr = tW1_daughterRefs.begin();
				     tW1_idr!= tW1_daughterRefs.end(); ++tW1_idr) 
				   {
				      if( tW1_idr->isAvailable() ) 
					{		       
					   const reco::GenParticleRef& tW1_genParticle = (*tW1_idr);
					   const reco::GenParticle *tW1_d = tW1_genParticle.get();
					   reco::GenParticle *pff = getUnique(tW1_d,0);
					   
					   if( fabs(pff->pdgId()) == 12 ||
					       fabs(pff->pdgId()) == 14 ) // nu
					     {
						tWnu1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 16 ) // nu_tau
					     {
						tWnutau1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 11 ||
					       fabs(pff->pdgId()) == 13 ) // l
					     {
						tWl1 = pff;
					     }		
					   if( fabs(pff->pdgId()) == 15 ) // tau
					     {
						tWtau1 = pff;
						
						const reco::GenParticleRefVector& tWtau1_daughterRefs = tWtau1->daughterRefVector();
						for(reco::GenParticleRefVector::const_iterator tWtau1_idr = tWtau1_daughterRefs.begin();
						    tWtau1_idr!= tWtau1_daughterRefs.end(); ++tWtau1_idr) 
						  {
						     if( tWtau1_idr->isAvailable() ) 
						       {		       
							  const reco::GenParticleRef& tWtau1_genParticle = (*tWtau1_idr);
							  const reco::GenParticle *tWtau1_d = tWtau1_genParticle.get();
							  reco::GenParticle *pfff = getUnique(tWtau1_d,0);
							  
							  if( fabs(pfff->pdgId()) == 12 ||
							      fabs(pfff->pdgId()) == 14 ) // nu
							    {
							       tWtaunu1 = pfff;
							    }		
							  if( fabs(pfff->pdgId()) == 16 ) // nu_tau
							    {
							       tWtaunutau1 = pfff;
							    }			
							  if( fabs(pfff->pdgId()) == 11 ||
							      fabs(pfff->pdgId()) == 13 ) // l
							    {
							       tWtaul1 = pfff;
							    }							  
						       }
						  }						
					     }
					   if( fabs(pff->pdgId()) <= 6 ) // q
					     {
						if( tWq11 && !tWq21 ) tWq21 = pff;
						if( !tWq11 ) tWq11 = pff;
						if( tWq11_IS && !tWq21_IS ) tWq21_IS = const_cast<reco::GenParticle*>(tW1_d);
						if( !tWq11_IS ) tWq11_IS = const_cast<reco::GenParticle*>(tW1_d);
					     }					   					   
					}
				   }				
			      }			    
			 }		       
		    }		  
	       }	     
	  }
     }   

   bool doCheck = 0;
   
   if( doCheck )
     {	
	if( Wl || Wtau || Wq1 )
	  {
	     int chan0 = 0;
	     if( Wl ) chan = chan0 + 0;
	     if( Wtaul ) chan = chan0 + 20;
	     if( Wtaunutau && ! Wtaul ) chan = chan0 + 40;
	     if( Wq1 && Wq2 ) chan = chan0 + 60;
	  }	
	
	if( chan < 0 )
	  {
	     std::cout << "Unknown channel found" << std::endl;
	     std::cout << "chan = " << chan << std::endl;

	     std::cout << "j1 = " << bool(j1) << std::endl;
	     std::cout << "j2 = " << bool(j2) << std::endl;
	     std::cout << "j3 = " << bool(j3) << std::endl;

	     /*
	     std::cout << "gamma->l1 = " << bool(gammal1) << std::endl;
	     std::cout << "gamma->l2 = " << bool(gammal2) << std::endl;
	     std::cout << "gamma->tau1 = " << bool(gammatau1) << std::endl;
	     std::cout << "gamma->tau1->l = " << bool(gammataul1) << std::endl;
	     std::cout << "gamma->tau1->nu = " << bool(gammataunu1) << std::endl;
	     std::cout << "gamma->tau1->nutau = " << bool(gammataunutau1) << std::endl;
	     std::cout << "gamma->tau2 = " << bool(gammatau2) << std::endl;
	     std::cout << "gamma->tau2->l = " << bool(gammataul2) << std::endl;
	     std::cout << "gamma->tau2->nu = " << bool(gammataunu2) << std::endl;
	     std::cout << "gamma->tau2->nutau = " << bool(gammataunutau2) << std::endl;*/

	     std::cout << "W = " << bool(W) << std::endl;
	     std::cout << "W->l = " << bool(Wl) << std::endl;
	     std::cout << "W->nu = " << bool(Wnu) << std::endl;
	     std::cout << "W->nutau = " << bool(Wnutau) << std::endl;
	     std::cout << "W->tau = " << bool(Wtau) << std::endl;
	     std::cout << "W->tau->l = " << bool(Wtaul) << std::endl;
	     std::cout << "W->tau->nu = " << bool(Wtaunu) << std::endl;
	     std::cout << "W->tau->nutau = " << bool(Wtaunutau) << std::endl;
	     std::cout << "W->q1 = " << bool(Wq1) << std::endl;
	     std::cout << "W->q2 = " << bool(Wq2) << std::endl;

	     exit(1);
	  }
     }
   
   // TLV

   /*
   if( gammal1 ) p4toTLV(gammal1->p4(),tree.mc_truth_gammal1_p4);
   if( gammal2 ) p4toTLV(gammal2->p4(),tree.mc_truth_gammal2_p4);
   if( gammatau1 ) p4toTLV(gammatau1->p4(),tree.mc_truth_gammatau1_p4);
   if( gammatau2 ) p4toTLV(gammatau2->p4(),tree.mc_truth_gammatau2_p4);
   if( gammataul1 ) p4toTLV(gammataul1->p4(),tree.mc_truth_gammataul1_p4);
   if( gammataul2 ) p4toTLV(gammataul2->p4(),tree.mc_truth_gammataul2_p4);
   if( gammataunu1 ) p4toTLV(gammataunu1->p4(),tree.mc_truth_gammataunu1_p4);
   if( gammataunu2 ) p4toTLV(gammataunu2->p4(),tree.mc_truth_gammataunu2_p4);
   if( gammataunutau1 ) p4toTLV(gammataunutau1->p4(),tree.mc_truth_gammataunutau1_p4);
   if( gammataunutau2 ) p4toTLV(gammataunutau2->p4(),tree.mc_truth_gammataunutau2_p4);*/

   if( W ) p4toTLV(W->p4(),tree.mc_truth_W_p4);
   if( Wnu ) p4toTLV(Wnu->p4(),tree.mc_truth_Wnu_p4);
   if( Wnutau ) p4toTLV(Wnutau->p4(),tree.mc_truth_Wnutau_p4);
   if( Wl ) p4toTLV(Wl->p4(),tree.mc_truth_Wl_p4);
   if( Wtau ) p4toTLV(Wtau->p4(),tree.mc_truth_Wtau_p4);
   if( Wtaunu ) p4toTLV(Wtaunu->p4(),tree.mc_truth_Wtaunu_p4);
   if( Wtaunutau ) p4toTLV(Wtaunutau->p4(),tree.mc_truth_Wtaunutau_p4);
   if( Wtaul ) p4toTLV(Wtaul->p4(),tree.mc_truth_Wtaul_p4);
   if( Wq1 ) p4toTLV(Wq1->p4(),tree.mc_truth_Wq1_p4);
   if( Wq2 ) p4toTLV(Wq2->p4(),tree.mc_truth_Wq2_p4);
   if( Wq1_IS ) p4toTLV(Wq1_IS->p4(),tree.mc_truth_Wq1_IS_p4);
   if( Wq2_IS ) p4toTLV(Wq2_IS->p4(),tree.mc_truth_Wq2_IS_p4);
   
   if( t1 ) p4toTLV(t1->p4(),tree.mc_truth_t1_p4);
   if( t2 ) p4toTLV(t2->p4(),tree.mc_truth_t2_p4);
   if( tb1 ) p4toTLV(tb1->p4(),tree.mc_truth_tb1_p4);
   if( tb2 ) p4toTLV(tb2->p4(),tree.mc_truth_tb2_p4);
   if( tb1_IS ) p4toTLV(tb1_IS->p4(),tree.mc_truth_tb1_IS_p4);
   if( tb2_IS ) p4toTLV(tb2_IS->p4(),tree.mc_truth_tb2_IS_p4);
   
   if( tW1 ) p4toTLV(tW1->p4(),tree.mc_truth_tW1_p4);
   if( tWnu1 ) p4toTLV(tWnu1->p4(),tree.mc_truth_tWnu1_p4);
   if( tWnutau1 ) p4toTLV(tWnutau1->p4(),tree.mc_truth_tWnutau1_p4);
   if( tWl1 ) p4toTLV(tWl1->p4(),tree.mc_truth_tWl1_p4);
   if( tWtau1 ) p4toTLV(tWtau1->p4(),tree.mc_truth_tWtau1_p4);
   if( tWtaunu1 ) p4toTLV(tWtaunu1->p4(),tree.mc_truth_tWtaunu1_p4);
   if( tWtaunutau1 ) p4toTLV(tWtaunutau1->p4(),tree.mc_truth_tWtaunutau1_p4);
   if( tWtaul1 ) p4toTLV(tWtaul1->p4(),tree.mc_truth_tWtaul1_p4);
   if( tWq11 ) p4toTLV(tWq11->p4(),tree.mc_truth_tWq11_p4);
   if( tWq21 ) p4toTLV(tWq21->p4(),tree.mc_truth_tWq21_p4);
   if( tWq11_IS ) p4toTLV(tWq11_IS->p4(),tree.mc_truth_tWq11_IS_p4);
   if( tWq21_IS ) p4toTLV(tWq21_IS->p4(),tree.mc_truth_tWq21_IS_p4);

   if( tW2 ) p4toTLV(tW2->p4(),tree.mc_truth_tW2_p4);
   if( tWnu2 ) p4toTLV(tWnu2->p4(),tree.mc_truth_tWnu2_p4);
   if( tWnutau2 ) p4toTLV(tWnutau2->p4(),tree.mc_truth_tWnutau2_p4);
   if( tWl2 ) p4toTLV(tWl2->p4(),tree.mc_truth_tWl2_p4);
   if( tWtau2 ) p4toTLV(tWtau2->p4(),tree.mc_truth_tWtau2_p4);
   if( tWtaunu2 ) p4toTLV(tWtaunu2->p4(),tree.mc_truth_tWtaunu2_p4);
   if( tWtaunutau2 ) p4toTLV(tWtaunutau2->p4(),tree.mc_truth_tWtaunutau2_p4);
   if( tWtaul2 ) p4toTLV(tWtaul2->p4(),tree.mc_truth_tWtaul2_p4);
   if( tWq12 ) p4toTLV(tWq12->p4(),tree.mc_truth_tWq12_p4);
   if( tWq22 ) p4toTLV(tWq22->p4(),tree.mc_truth_tWq22_p4);
   if( tWq12_IS ) p4toTLV(tWq12_IS->p4(),tree.mc_truth_tWq12_IS_p4);
   if( tWq22_IS ) p4toTLV(tWq22_IS->p4(),tree.mc_truth_tWq22_IS_p4);

   if( j1 ) p4toTLV(j1->p4(),tree.mc_truth_j1_p4);
   if( j2 ) p4toTLV(j2->p4(),tree.mc_truth_j2_p4);
   if( j3 ) p4toTLV(j3->p4(),tree.mc_truth_j3_p4);

   // pt
   
   /*
   if( gammal1 ) tree.mc_truth_gammal1_pt = gammal1->p4().pt();
   if( gammal2 ) tree.mc_truth_gammal2_pt = gammal2->p4().pt();
   if( gammatau1 ) tree.mc_truth_gammatau1_pt = gammatau1->p4().pt();
   if( gammatau2 ) tree.mc_truth_gammatau2_pt = gammatau2->p4().pt();
   if( gammataul1 ) tree.mc_truth_gammataul1_pt = gammataul1->p4().pt();
   if( gammataul2 ) tree.mc_truth_gammataul2_pt = gammataul2->p4().pt();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_pt = gammataunu1->p4().pt();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_pt = gammataunu2->p4().pt();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_pt = gammataunutau1->p4().pt();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_pt = gammataunutau2->p4().pt();*/

   if( W ) tree.mc_truth_W_pt = W->p4().pt();
   if( Wnu ) tree.mc_truth_Wnu_pt = Wnu->p4().pt();
   if( Wnutau ) tree.mc_truth_Wnutau_pt = Wnutau->p4().pt();
   if( Wl ) tree.mc_truth_Wl_pt = Wl->p4().pt();
   if( Wtau ) tree.mc_truth_Wtau_pt = Wtau->p4().pt();
   if( Wtaunu ) tree.mc_truth_Wtaunu_pt = Wtaunu->p4().pt();
   if( Wtaunutau ) tree.mc_truth_Wtaunutau_pt = Wtaunutau->p4().pt();
   if( Wtaul ) tree.mc_truth_Wtaul_pt = Wtaul->p4().pt();
   if( Wq1 ) tree.mc_truth_Wq1_pt = Wq1->p4().pt();
   if( Wq2 ) tree.mc_truth_Wq2_pt = Wq2->p4().pt();
   if( Wq1_IS ) tree.mc_truth_Wq1_IS_pt = Wq1_IS->p4().pt();
   if( Wq2_IS ) tree.mc_truth_Wq2_IS_pt = Wq2_IS->p4().pt();
   
   if( t1 ) tree.mc_truth_t1_pt = t1->p4().pt();
   if( t2 ) tree.mc_truth_t2_pt = t2->p4().pt();
   if( tb1 ) tree.mc_truth_tb1_pt = tb1->p4().pt();
   if( tb2 ) tree.mc_truth_tb2_pt = tb2->p4().pt();
   if( tb1_IS ) tree.mc_truth_tb1_IS_pt = tb1_IS->p4().pt();
   if( tb2_IS ) tree.mc_truth_tb2_IS_pt = tb2_IS->p4().pt();
   
   if( tW1 ) tree.mc_truth_tW1_pt = tW1->p4().pt();
   if( tWnu1 ) tree.mc_truth_tWnu1_pt = tWnu1->p4().pt();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_pt = tWnutau1->p4().pt();
   if( tWl1 ) tree.mc_truth_tWl1_pt = tWl1->p4().pt();
   if( tWtau1 ) tree.mc_truth_tWtau1_pt = tWtau1->p4().pt();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_pt = tWtaunu1->p4().pt();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_pt = tWtaunutau1->p4().pt();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_pt = tWtaul1->p4().pt();
   if( tWq11 ) tree.mc_truth_tWq11_pt = tWq11->p4().pt();
   if( tWq21 ) tree.mc_truth_tWq21_pt = tWq21->p4().pt();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_pt = tWq11_IS->p4().pt();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_pt = tWq21_IS->p4().pt();
   
   if( tW2 ) tree.mc_truth_tW2_pt = tW2->p4().pt();
   if( tWnu2 ) tree.mc_truth_tWnu2_pt = tWnu2->p4().pt();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_pt = tWnutau2->p4().pt();
   if( tWl2 ) tree.mc_truth_tWl2_pt = tWl2->p4().pt();
   if( tWtau2 ) tree.mc_truth_tWtau2_pt = tWtau2->p4().pt();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_pt = tWtaunu2->p4().pt();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_pt = tWtaunutau2->p4().pt();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_pt = tWtaul2->p4().pt();
   if( tWq12 ) tree.mc_truth_tWq12_pt = tWq12->p4().pt();
   if( tWq22 ) tree.mc_truth_tWq22_pt = tWq22->p4().pt();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_pt = tWq12_IS->p4().pt();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_pt = tWq22_IS->p4().pt();
   
   if( j1 ) tree.mc_truth_j1_pt = j1->p4().pt();
   if( j2 ) tree.mc_truth_j2_pt = j2->p4().pt();
   if( j3 ) tree.mc_truth_j3_pt = j3->p4().pt();

   // eta
   
   /*
   if( gammal1 ) tree.mc_truth_gammal1_eta = gammal1->p4().eta();
   if( gammal2 ) tree.mc_truth_gammal2_eta = gammal2->p4().eta();
   if( gammatau1 ) tree.mc_truth_gammatau1_eta = gammatau1->p4().eta();
   if( gammatau2 ) tree.mc_truth_gammatau2_eta = gammatau2->p4().eta();
   if( gammataul1 ) tree.mc_truth_gammataul1_eta = gammataul1->p4().eta();
   if( gammataul2 ) tree.mc_truth_gammataul2_eta = gammataul2->p4().eta();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_eta = gammataunu1->p4().eta();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_eta = gammataunu2->p4().eta();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_eta = gammataunutau1->p4().eta();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_eta = gammataunutau2->p4().eta();*/

   if( W ) tree.mc_truth_W_eta = W->p4().eta();
   if( Wnu ) tree.mc_truth_Wnu_eta = Wnu->p4().eta();
   if( Wnutau ) tree.mc_truth_Wnutau_eta = Wnutau->p4().eta();
   if( Wl ) tree.mc_truth_Wl_eta = Wl->p4().eta();
   if( Wtau ) tree.mc_truth_Wtau_eta = Wtau->p4().eta();
   if( Wtaunu ) tree.mc_truth_Wtaunu_eta = Wtaunu->p4().eta();
   if( Wtaunutau ) tree.mc_truth_Wtaunutau_eta = Wtaunutau->p4().eta();
   if( Wtaul ) tree.mc_truth_Wtaul_eta = Wtaul->p4().eta();
   if( Wq1 ) tree.mc_truth_Wq1_eta = Wq1->p4().eta();
   if( Wq2 ) tree.mc_truth_Wq2_eta = Wq2->p4().eta();
   if( Wq1_IS ) tree.mc_truth_Wq1_IS_eta = Wq1_IS->p4().eta();
   if( Wq2_IS ) tree.mc_truth_Wq2_IS_eta = Wq2_IS->p4().eta();
   
   if( t1 ) tree.mc_truth_t1_eta = t1->p4().eta();
   if( t2 ) tree.mc_truth_t2_eta = t2->p4().eta();
   if( tb1 ) tree.mc_truth_tb1_eta = tb1->p4().eta();
   if( tb2 ) tree.mc_truth_tb2_eta = tb2->p4().eta();
   if( tb1_IS ) tree.mc_truth_tb1_IS_eta = tb1_IS->p4().eta();
   if( tb2_IS ) tree.mc_truth_tb2_IS_eta = tb2_IS->p4().eta();
   
   if( tW1 ) tree.mc_truth_tW1_eta = tW1->p4().eta();
   if( tWnu1 ) tree.mc_truth_tWnu1_eta = tWnu1->p4().eta();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_eta = tWnutau1->p4().eta();
   if( tWl1 ) tree.mc_truth_tWl1_eta = tWl1->p4().eta();
   if( tWtau1 ) tree.mc_truth_tWtau1_eta = tWtau1->p4().eta();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_eta = tWtaunu1->p4().eta();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_eta = tWtaunutau1->p4().eta();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_eta = tWtaul1->p4().eta();
   if( tWq11 ) tree.mc_truth_tWq11_eta = tWq11->p4().eta();
   if( tWq21 ) tree.mc_truth_tWq21_eta = tWq21->p4().eta();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_eta = tWq11_IS->p4().eta();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_eta = tWq21_IS->p4().eta();
   
   if( tW2 ) tree.mc_truth_tW2_eta = tW2->p4().eta();
   if( tWnu2 ) tree.mc_truth_tWnu2_eta = tWnu2->p4().eta();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_eta = tWnutau2->p4().eta();
   if( tWl2 ) tree.mc_truth_tWl2_eta = tWl2->p4().eta();
   if( tWtau2 ) tree.mc_truth_tWtau2_eta = tWtau2->p4().eta();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_eta = tWtaunu2->p4().eta();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_eta = tWtaunutau2->p4().eta();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_eta = tWtaul2->p4().eta();
   if( tWq12 ) tree.mc_truth_tWq12_eta = tWq12->p4().eta();
   if( tWq22 ) tree.mc_truth_tWq22_eta = tWq22->p4().eta();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_eta = tWq12_IS->p4().eta();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_eta = tWq22_IS->p4().eta();
   
   if( j1 ) tree.mc_truth_j1_eta = j1->p4().eta();
   if( j2 ) tree.mc_truth_j2_eta = j2->p4().eta();
   if( j3 ) tree.mc_truth_j3_eta = j3->p4().eta();

   // phi
   
   /*
   if( gammal1 ) tree.mc_truth_gammal1_phi = gammal1->p4().phi();
   if( gammal2 ) tree.mc_truth_gammal2_phi = gammal2->p4().phi();
   if( gammatau1 ) tree.mc_truth_gammatau1_phi = gammatau1->p4().phi();
   if( gammatau2 ) tree.mc_truth_gammatau2_phi = gammatau2->p4().phi();
   if( gammataul1 ) tree.mc_truth_gammataul1_phi = gammataul1->p4().phi();
   if( gammataul2 ) tree.mc_truth_gammataul2_phi = gammataul2->p4().phi();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_phi = gammataunu1->p4().phi();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_phi = gammataunu2->p4().phi();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_phi = gammataunutau1->p4().phi();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_phi = gammataunutau2->p4().phi();*/

   if( W ) tree.mc_truth_W_phi = W->p4().phi();
   if( Wnu ) tree.mc_truth_Wnu_phi = Wnu->p4().phi();
   if( Wnutau ) tree.mc_truth_Wnutau_phi = Wnutau->p4().phi();
   if( Wl ) tree.mc_truth_Wl_phi = Wl->p4().phi();
   if( Wtau ) tree.mc_truth_Wtau_phi = Wtau->p4().phi();
   if( Wtaunu ) tree.mc_truth_Wtaunu_phi = Wtaunu->p4().phi();
   if( Wtaunutau ) tree.mc_truth_Wtaunutau_phi = Wtaunutau->p4().phi();
   if( Wtaul ) tree.mc_truth_Wtaul_phi = Wtaul->p4().phi();
   if( Wq1 ) tree.mc_truth_Wq1_phi = Wq1->p4().phi();
   if( Wq2 ) tree.mc_truth_Wq2_phi = Wq2->p4().phi();
   if( Wq1_IS ) tree.mc_truth_Wq1_IS_phi = Wq1_IS->p4().phi();
   if( Wq2_IS ) tree.mc_truth_Wq2_IS_phi = Wq2_IS->p4().phi();
   
   if( t1 ) tree.mc_truth_t1_phi = t1->p4().phi();
   if( t2 ) tree.mc_truth_t2_phi = t2->p4().phi();
   if( tb1 ) tree.mc_truth_tb1_phi = tb1->p4().phi();
   if( tb2 ) tree.mc_truth_tb2_phi = tb2->p4().phi();
   if( tb1_IS ) tree.mc_truth_tb1_IS_phi = tb1_IS->p4().phi();
   if( tb2_IS ) tree.mc_truth_tb2_IS_phi = tb2_IS->p4().phi();
   
   if( tW1 ) tree.mc_truth_tW1_phi = tW1->p4().phi();
   if( tWnu1 ) tree.mc_truth_tWnu1_phi = tWnu1->p4().phi();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_phi = tWnutau1->p4().phi();
   if( tWl1 ) tree.mc_truth_tWl1_phi = tWl1->p4().phi();
   if( tWtau1 ) tree.mc_truth_tWtau1_phi = tWtau1->p4().phi();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_phi = tWtaunu1->p4().phi();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_phi = tWtaunutau1->p4().phi();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_phi = tWtaul1->p4().phi();
   if( tWq11 ) tree.mc_truth_tWq11_phi = tWq11->p4().phi();
   if( tWq21 ) tree.mc_truth_tWq21_phi = tWq21->p4().phi();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_phi = tWq11_IS->p4().phi();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_phi = tWq21_IS->p4().phi();
   
   if( tW2 ) tree.mc_truth_tW2_phi = tW2->p4().phi();
   if( tWnu2 ) tree.mc_truth_tWnu2_phi = tWnu2->p4().phi();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_phi = tWnutau2->p4().phi();
   if( tWl2 ) tree.mc_truth_tWl2_phi = tWl2->p4().phi();
   if( tWtau2 ) tree.mc_truth_tWtau2_phi = tWtau2->p4().phi();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_phi = tWtaunu2->p4().phi();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_phi = tWtaunutau2->p4().phi();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_phi = tWtaul2->p4().phi();
   if( tWq12 ) tree.mc_truth_tWq12_phi = tWq12->p4().phi();
   if( tWq22 ) tree.mc_truth_tWq22_phi = tWq22->p4().phi();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_phi = tWq12_IS->p4().phi();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_phi = tWq22_IS->p4().phi();
   
   if( j1 ) tree.mc_truth_j1_phi = j1->p4().phi();
   if( j2 ) tree.mc_truth_j2_phi = j2->p4().phi();
   if( j3 ) tree.mc_truth_j3_phi = j3->p4().phi();

   // E
   
   /*
   if( gammal1 ) tree.mc_truth_gammal1_E = gammal1->p4().E();
   if( gammal2 ) tree.mc_truth_gammal2_E = gammal2->p4().E();
   if( gammatau1 ) tree.mc_truth_gammatau1_E = gammatau1->p4().E();
   if( gammatau2 ) tree.mc_truth_gammatau2_E = gammatau2->p4().E();
   if( gammataul1 ) tree.mc_truth_gammataul1_E = gammataul1->p4().E();
   if( gammataul2 ) tree.mc_truth_gammataul2_E = gammataul2->p4().E();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_E = gammataunu1->p4().E();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_E = gammataunu2->p4().E();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_E = gammataunutau1->p4().E();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_E = gammataunutau2->p4().E();*/

   if( W ) tree.mc_truth_W_E = W->p4().E();
   if( Wnu ) tree.mc_truth_Wnu_E = Wnu->p4().E();
   if( Wnutau ) tree.mc_truth_Wnutau_E = Wnutau->p4().E();
   if( Wl ) tree.mc_truth_Wl_E = Wl->p4().E();
   if( Wtau ) tree.mc_truth_Wtau_E = Wtau->p4().E();
   if( Wtaunu ) tree.mc_truth_Wtaunu_E = Wtaunu->p4().E();
   if( Wtaunutau ) tree.mc_truth_Wtaunutau_E = Wtaunutau->p4().E();
   if( Wtaul ) tree.mc_truth_Wtaul_E = Wtaul->p4().E();
   if( Wq1 ) tree.mc_truth_Wq1_E = Wq1->p4().E();
   if( Wq2 ) tree.mc_truth_Wq2_E = Wq2->p4().E();
   if( Wq1_IS ) tree.mc_truth_Wq1_IS_E = Wq1_IS->p4().E();
   if( Wq2_IS ) tree.mc_truth_Wq2_IS_E = Wq2_IS->p4().E();
   
   if( t1 ) tree.mc_truth_t1_E = t1->p4().E();
   if( t2 ) tree.mc_truth_t2_E = t2->p4().E();
   if( tb1 ) tree.mc_truth_tb1_E = tb1->p4().E();
   if( tb2 ) tree.mc_truth_tb2_E = tb2->p4().E();
   if( tb1_IS ) tree.mc_truth_tb1_IS_E = tb1_IS->p4().E();
   if( tb2_IS ) tree.mc_truth_tb2_IS_E = tb2_IS->p4().E();
   
   if( tW1 ) tree.mc_truth_tW1_E = tW1->p4().E();
   if( tWnu1 ) tree.mc_truth_tWnu1_E = tWnu1->p4().E();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_E = tWnutau1->p4().E();
   if( tWl1 ) tree.mc_truth_tWl1_E = tWl1->p4().E();
   if( tWtau1 ) tree.mc_truth_tWtau1_E = tWtau1->p4().E();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_E = tWtaunu1->p4().E();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_E = tWtaunutau1->p4().E();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_E = tWtaul1->p4().E();
   if( tWq11 ) tree.mc_truth_tWq11_E = tWq11->p4().E();
   if( tWq21 ) tree.mc_truth_tWq21_E = tWq21->p4().E();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_E = tWq11_IS->p4().E();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_E = tWq21_IS->p4().E();
   
   if( tW2 ) tree.mc_truth_tW2_E = tW2->p4().E();
   if( tWnu2 ) tree.mc_truth_tWnu2_E = tWnu2->p4().E();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_E = tWnutau2->p4().E();
   if( tWl2 ) tree.mc_truth_tWl2_E = tWl2->p4().E();
   if( tWtau2 ) tree.mc_truth_tWtau2_E = tWtau2->p4().E();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_E = tWtaunu2->p4().E();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_E = tWtaunutau2->p4().E();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_E = tWtaul2->p4().E();
   if( tWq12 ) tree.mc_truth_tWq12_E = tWq12->p4().E();
   if( tWq22 ) tree.mc_truth_tWq22_E = tWq22->p4().E();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_E = tWq12_IS->p4().E();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_E = tWq22_IS->p4().E();
   
   if( j1 ) tree.mc_truth_j1_E = j1->p4().E();
   if( j2 ) tree.mc_truth_j2_E = j2->p4().E();
   if( j3 ) tree.mc_truth_j3_E = j3->p4().E();
   
   // pdgId

   /*
   if( gammal1 ) tree.mc_truth_gammal1_id = gammal1->pdgId();
   if( gammal2 ) tree.mc_truth_gammal2_id = gammal2->pdgId();
   if( gammatau1 ) tree.mc_truth_gammatau1_id = gammatau1->pdgId();
   if( gammatau2 ) tree.mc_truth_gammatau2_id = gammatau2->pdgId();
   if( gammataul1 ) tree.mc_truth_gammataul1_id = gammataul1->pdgId();
   if( gammataul2 ) tree.mc_truth_gammataul2_id = gammataul2->pdgId();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_id = gammataunu1->pdgId();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_id = gammataunu2->pdgId();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_id = gammataunutau1->pdgId();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_id = gammataunutau2->pdgId();*/

   if( W ) tree.mc_truth_W_id = W->pdgId();
   if( Wnu ) tree.mc_truth_Wnu_id = Wnu->pdgId();
   if( Wnutau ) tree.mc_truth_Wnutau_id = Wnutau->pdgId();
   if( Wl ) tree.mc_truth_Wl_id = Wl->pdgId();
   if( Wtau ) tree.mc_truth_Wtau_id = Wtau->pdgId();
   if( Wtaunu ) tree.mc_truth_Wtaunu_id = Wtaunu->pdgId();
   if( Wtaunutau ) tree.mc_truth_Wtaunutau_id = Wtaunutau->pdgId();
   if( Wtaul ) tree.mc_truth_Wtaul_id = Wtaul->pdgId();
   if( Wq1 ) tree.mc_truth_Wq1_id = Wq1->pdgId();
   if( Wq2 ) tree.mc_truth_Wq2_id = Wq2->pdgId();
   if( Wq1_IS ) tree.mc_truth_Wq1_IS_id = Wq1_IS->pdgId();
   if( Wq2_IS ) tree.mc_truth_Wq2_IS_id = Wq2_IS->pdgId();
   
   if( t1 ) tree.mc_truth_t1_id = t1->pdgId();
   if( t2 ) tree.mc_truth_t2_id = t2->pdgId();
   if( tb1 ) tree.mc_truth_tb1_id = tb1->pdgId();
   if( tb2 ) tree.mc_truth_tb2_id = tb2->pdgId();
   if( tb1_IS ) tree.mc_truth_tb1_IS_id = tb1_IS->pdgId();
   if( tb2_IS ) tree.mc_truth_tb2_IS_id = tb2_IS->pdgId();
   
   if( tW1 ) tree.mc_truth_tW1_id = tW1->pdgId();
   if( tWnu1 ) tree.mc_truth_tWnu1_id = tWnu1->pdgId();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_id = tWnutau1->pdgId();
   if( tWl1 ) tree.mc_truth_tWl1_id = tWl1->pdgId();
   if( tWtau1 ) tree.mc_truth_tWtau1_id = tWtau1->pdgId();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_id = tWtaunu1->pdgId();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_id = tWtaunutau1->pdgId();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_id = tWtaul1->pdgId();
   if( tWq11 ) tree.mc_truth_tWq11_id = tWq11->pdgId();
   if( tWq21 ) tree.mc_truth_tWq21_id = tWq21->pdgId();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_id = tWq11_IS->pdgId();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_id = tWq21_IS->pdgId();
   
   if( tW2 ) tree.mc_truth_tW2_id = tW2->pdgId();
   if( tWnu2 ) tree.mc_truth_tWnu2_id = tWnu2->pdgId();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_id = tWnutau2->pdgId();
   if( tWl2 ) tree.mc_truth_tWl2_id = tWl2->pdgId();
   if( tWtau2 ) tree.mc_truth_tWtau2_id = tWtau2->pdgId();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_id = tWtaunu2->pdgId();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_id = tWtaunutau2->pdgId();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_id = tWtaul2->pdgId();
   if( tWq12 ) tree.mc_truth_tWq12_id = tWq12->pdgId();
   if( tWq22 ) tree.mc_truth_tWq22_id = tWq22->pdgId();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_id = tWq12_IS->pdgId();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_id = tWq22_IS->pdgId();
   
   if( j1 ) tree.mc_truth_j1_id = j1->pdgId();
   if( j2 ) tree.mc_truth_j2_id = j2->pdgId();
   if( j3 ) tree.mc_truth_j3_id = j3->pdgId();
   
   // status

   /*
   if( gammal1 ) tree.mc_truth_gammal1_status = gammal1->status();
   if( gammal2 ) tree.mc_truth_gammal2_status = gammal2->status();
   if( gammatau1 ) tree.mc_truth_gammatau1_status = gammatau1->status();
   if( gammatau2 ) tree.mc_truth_gammatau2_status = gammatau2->status();
   if( gammataul1 ) tree.mc_truth_gammataul1_status = gammataul1->status();
   if( gammataul2 ) tree.mc_truth_gammataul2_status = gammataul2->status();
   if( gammataunu1 ) tree.mc_truth_gammataunu1_status = gammataunu1->status();
   if( gammataunu2 ) tree.mc_truth_gammataunu2_status = gammataunu2->status();
   if( gammataunutau1 ) tree.mc_truth_gammataunutau1_status = gammataunutau1->status();
   if( gammataunutau2 ) tree.mc_truth_gammataunutau2_status = gammataunutau2->status();*/

   if( W ) tree.mc_truth_W_status = W->status();
   if( Wnu ) tree.mc_truth_Wnu_status = Wnu->status();
   if( Wnutau ) tree.mc_truth_Wnutau_status = Wnutau->status();
   if( Wl ) tree.mc_truth_Wl_status = Wl->status();
   if( Wtau ) tree.mc_truth_Wtau_status = Wtau->status();
   if( Wtaunu ) tree.mc_truth_Wtaunu_status = Wtaunu->status();
   if( Wtaunutau ) tree.mc_truth_Wtaunutau_status = Wtaunutau->status();
   if( Wtaul ) tree.mc_truth_Wtaul_status = Wtaul->status();
   if( Wq1 ) tree.mc_truth_Wq1_status = Wq1->status();
   if( Wq2 ) tree.mc_truth_Wq2_status = Wq2->status();
   if( Wq1_IS ) tree.mc_truth_Wq1_IS_status = Wq1_IS->status();
   if( Wq2_IS ) tree.mc_truth_Wq2_IS_status = Wq2_IS->status();
   
   if( t1 ) tree.mc_truth_t1_status = t1->status();
   if( t2 ) tree.mc_truth_t2_status = t2->status();
   if( tb1 ) tree.mc_truth_tb1_status = tb1->status();
   if( tb2 ) tree.mc_truth_tb2_status = tb2->status();
   if( tb1_IS ) tree.mc_truth_tb1_IS_status = tb1_IS->status();
   if( tb2_IS ) tree.mc_truth_tb2_IS_status = tb2_IS->status();
   
   if( tW1 ) tree.mc_truth_tW1_status = tW1->status();
   if( tWnu1 ) tree.mc_truth_tWnu1_status = tWnu1->status();
   if( tWnutau1 ) tree.mc_truth_tWnutau1_status = tWnutau1->status();
   if( tWl1 ) tree.mc_truth_tWl1_status = tWl1->status();
   if( tWtau1 ) tree.mc_truth_tWtau1_status = tWtau1->status();
   if( tWtaunu1 ) tree.mc_truth_tWtaunu1_status = tWtaunu1->status();
   if( tWtaunutau1 ) tree.mc_truth_tWtaunutau1_status = tWtaunutau1->status();
   if( tWtaul1 ) tree.mc_truth_tWtaul1_status = tWtaul1->status();
   if( tWq11 ) tree.mc_truth_tWq11_status = tWq11->status();
   if( tWq21 ) tree.mc_truth_tWq21_status = tWq21->status();
   if( tWq11_IS ) tree.mc_truth_tWq11_IS_status = tWq11_IS->status();
   if( tWq21_IS ) tree.mc_truth_tWq21_IS_status = tWq21_IS->status();

   if( tW2 ) tree.mc_truth_tW2_status = tW2->status();
   if( tWnu2 ) tree.mc_truth_tWnu2_status = tWnu2->status();
   if( tWnutau2 ) tree.mc_truth_tWnutau2_status = tWnutau2->status();
   if( tWl2 ) tree.mc_truth_tWl2_status = tWl2->status();
   if( tWtau2 ) tree.mc_truth_tWtau2_status = tWtau2->status();
   if( tWtaunu2 ) tree.mc_truth_tWtaunu2_status = tWtaunu2->status();
   if( tWtaunutau2 ) tree.mc_truth_tWtaunutau2_status = tWtaunutau2->status();
   if( tWtaul2 ) tree.mc_truth_tWtaul2_status = tWtaul2->status();
   if( tWq12 ) tree.mc_truth_tWq12_status = tWq12->status();
   if( tWq22 ) tree.mc_truth_tWq22_status = tWq22->status();
   if( tWq12_IS ) tree.mc_truth_tWq12_IS_status = tWq12_IS->status();
   if( tWq22_IS ) tree.mc_truth_tWq22_IS_status = tWq22_IS->status();

   if( j1 ) tree.mc_truth_j1_status = j1->status();
   if( j2 ) tree.mc_truth_j2_status = j2->status();
   if( j3 ) tree.mc_truth_j3_status = j3->status();
}

// tZq MC analyzer
void MCTruth::fillTZQSignalGenParticles(const edm::Event& iEvent,
					const edm::EventSetup& iSetup,
					FlatTree& tree,
					const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   reco::GenParticle *Z = 0;
   
   reco::GenParticle *Zl1 = 0;
   reco::GenParticle *Zl2 = 0;
   reco::GenParticle *Ztau1 = 0;
   reco::GenParticle *Ztau2 = 0;
   reco::GenParticle *Ztaul1 = 0;
   reco::GenParticle *Ztaul2 = 0;
   reco::GenParticle *Ztaunu1 = 0;
   reco::GenParticle *Ztaunu2 = 0;
   reco::GenParticle *Ztaunutau1 = 0;
   reco::GenParticle *Ztaunutau2 = 0;
   
   reco::GenParticle *t = 0;
   reco::GenParticle *tb = 0;
   reco::GenParticle *tb_IS = 0;
   reco::GenParticle *tW = 0;
   reco::GenParticle *tWnu = 0;
   reco::GenParticle *tWnutau = 0;
   reco::GenParticle *tWl = 0;
   reco::GenParticle *tWtau = 0;
   reco::GenParticle *tWtaunu = 0;
   reco::GenParticle *tWtaunutau = 0;
   reco::GenParticle *tWtaul = 0;
   reco::GenParticle *tWq1 = 0;
   reco::GenParticle *tWq2 = 0;
   reco::GenParticle *tWq1_IS = 0;
   reco::GenParticle *tWq2_IS = 0;

   reco::GenParticle *j1 = 0;
   reco::GenParticle *j2 = 0;
   reco::GenParticle *j3 = 0;
   
   int chan = -666;

   // 0   = (t->bW,W->lnu)
   // 1   = (t->bW,W->qq)
   // 2   = (t->bW,W->tauLnu)
   // 3   = (t->bW,W->tauHnu)
   
   // (Z->ll)             +0
   // (Z->tauLtauL)       +20
   // (Z->tauHtauH)       +40
   // (Z->tauLtauH)       +60
   // (Z->tauHtauL)       +80
         
   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;
   
   int ipart = 0;

//   std::cout << "event" << std::endl;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	int barcode = ipart; // in CMSSW barcode is the index of genParticle in the event
	// https://twiki.cern.ch/twiki/bin/view/CMS/GenParticles2HepMCConverter
	ipart++;
	
//	if( fabs(mcp->pdgId()) <= 6 )
//	  {	     
//	     std::cout << "pdgId=" << mcp->pdgId() << " pt=" <<
//	       mcp->pt() << " status=" << mcp->status() << 
//	       " motherPdg=" << mcp->mother()->pdgId() << std::endl;
//	  }	
	
	// Additional partons (up to three)
	if( (fabs(mcp->pdgId()) <= 6 || fabs(mcp->pdgId()) == 21) &&
	    mcp->status() == 23 && barcode == 8 )
	  {
	     if( !j1 )
	       j1 = mcp;
	     else if( !j2 )
	       j2 = mcp;
	     else if( !j3 )
	       j3 = mcp;
	  }	

//	if( fabs(mcp->pdgId()) == 15 )
//	  {	     
//	     std::cout << mcp->pdgId() << " " << mcp->status() << " " << mcp->numberOfDaughters() << std::endl;
//	     
//	     const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
//	     for(reco::GenParticleRefVector::const_iterator mcp_idr = daughterRefs.begin(); 
//		 mcp_idr!= daughterRefs.end(); ++mcp_idr)
//	       {
//		  if( mcp_idr->isAvailable() ) 
//		    {
//		       const reco::GenParticleRef& genParticle = (*mcp_idr);
//		       const reco::GenParticle *mcp_d = genParticle.get();
//		       std::cout << "daughter " << mcp_d->pdgId() << " " << mcp_d->status() << std::endl;
//		    }
//	       }	     	     
//	  }

	if( ((fabs(mcp->pdgId()) == 11 ||
	      fabs(mcp->pdgId()) == 13) &&
	     mcp->status() == 3) ||
	    (fabs(mcp->pdgId()) == 15 &&
		mcp->status() == 2) )
	  {
	     if( fabs(mcp->pdgId()) == 11 ||
		 fabs(mcp->pdgId()) == 13 ) // l
	       {
		  if( Zl1 && !Zl2 ) {Zl2 = mcp;}
		  if( !Zl1 ) {Zl1 = mcp;}
	       }			    
	     if( fabs(mcp->pdgId()) == 15 ) // tau
	       {
		  if( Ztau1 )
		    {
		       Ztau2 = mcp;
				      
		       const reco::GenParticleRefVector& daughterRefs = Ztau2->daughterRefVector();
		       for(reco::GenParticleRefVector::const_iterator Ztau2_idr = daughterRefs.begin(); 
			   Ztau2_idr!= daughterRefs.end(); ++Ztau2_idr)
			 {
			    if( Ztau2_idr->isAvailable() ) 
			      {		       
				 const reco::GenParticleRef& genParticle = (*Ztau2_idr);
				 const reco::GenParticle *Ztau2_d = genParticle.get();
				 reco::GenParticle *pfff = getUnique(Ztau2_d,0);
				 
				 if( fabs(pfff->pdgId()) == 12 ||
				     fabs(pfff->pdgId()) == 14 ) // nu
				   {
				      Ztaunu2 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 16 ) // nutau
				   {
				      Ztaunutau2 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 11 ||
				     fabs(pfff->pdgId()) == 13 ) // l
				   {
				      Ztaul2 = pfff;
				   }
			      }
			 }							  
		    }
		  if( !Ztau1 )
		    {
		       Ztau1 = mcp;
		       
		       const reco::GenParticleRefVector& daughterRefs = Ztau1->daughterRefVector();
		       for(reco::GenParticleRefVector::const_iterator Ztau1_idr = daughterRefs.begin(); 
			   Ztau1_idr!= daughterRefs.end(); ++Ztau1_idr)
			 {
			    if( Ztau1_idr->isAvailable() ) 
			      {		       
				 const reco::GenParticleRef& genParticle = (*Ztau1_idr);
				 const reco::GenParticle *Ztau1_d = genParticle.get();
				 reco::GenParticle *pfff = getUnique(Ztau1_d,0);
				 
				 if( fabs(pfff->pdgId()) == 12 ||
				     fabs(pfff->pdgId()) == 14 ) // nu
				   {
				      Ztaunu1 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 16 ) // nutau
				   {
				      Ztaunutau1 = pfff;
				   }
				 if( fabs(pfff->pdgId()) == 11 ||
				     fabs(pfff->pdgId()) == 13 ) // l
				   {
				      Ztaul1 = pfff;
				   }
			      }
			 }							  
		    }
	       }	     
	  }		  
	
	// Z
	if( fabs(mcp->pdgId()) == 23 )
	  {
	     if( !Z ) {Z = mcp;}
	     
	     if( Z )
	       {
		  const reco::GenParticleRefVector& daughterRefs = Z->daughterRefVector();
		  for(reco::GenParticleRefVector::const_iterator Z_idr = daughterRefs.begin(); 
		      Z_idr!= daughterRefs.end(); ++Z_idr)
		    {
		       if( Z_idr->isAvailable() ) 
			 {		       
			    const reco::GenParticleRef& genParticle = (*Z_idr);
			    const reco::GenParticle *Z_d = genParticle.get();
			    reco::GenParticle *pff = getUnique(Z_d,0);
			    
			    if( fabs(pff->pdgId()) == 11 ||
				fabs(pff->pdgId()) == 13 ) // l
			      {
				 if( Zl1 && !Zl2 ) {Zl2 = pff;}
				 if( !Zl1 ) {Zl1 = pff;}
			      }			    
			    if( fabs(pff->pdgId()) == 15 ) // tau
			      {
				 if( Ztau1 )
				   {
				      Ztau2 = pff;
				      
				      const reco::GenParticleRefVector& daughterRefs = Ztau2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator Ztau2_idr = daughterRefs.begin(); 
					  Ztau2_idr!= daughterRefs.end(); ++Ztau2_idr)
					{
					   if( Ztau2_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*Ztau2_idr);
						const reco::GenParticle *Ztau2_d = genParticle.get();
						reco::GenParticle *pfff = getUnique(Ztau2_d,0);
						
						if( fabs(pfff->pdgId()) == 12 ||
						    fabs(pfff->pdgId()) == 14 ) // nu
						  {
						     Ztaunu2 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 16 ) // nutau
						  {
						     Ztaunutau2 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 11 ||
						    fabs(pfff->pdgId()) == 13 ) // l
						  {
						     Ztaul2 = pfff;
						  }
					     }
					}							  
				   }
				 if( !Ztau1 )
				   {
				      Ztau1 = pff;
				      
				      const reco::GenParticleRefVector& daughterRefs = Ztau1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator Ztau1_idr = daughterRefs.begin(); 
					  Ztau1_idr!= daughterRefs.end(); ++Ztau1_idr)
					{
					   if( Ztau1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*Ztau1_idr);
						const reco::GenParticle *Ztau1_d = genParticle.get();
						reco::GenParticle *pfff = getUnique(Ztau1_d,0);
						
						if( fabs(pfff->pdgId()) == 12 ||
						    fabs(pfff->pdgId()) == 14 ) // nu
						  {
						     Ztaunu1 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 16 ) // nutau
						  {
						     Ztaunutau1 = pfff;
						  }
						if( fabs(pfff->pdgId()) == 11 ||
						    fabs(pfff->pdgId()) == 13 ) // l
						  {
						     Ztaul1 = pfff;
						  }
					     }
					}							  
				   }						     
			      }
			 }					   
		    }				      
	       }
	  }	
	
	// top decays
	if( fabs(mcp->pdgId()) == 6
	    && ( (mcp->status() == 62) || 
		 (mcp->status() == 3)
	       ) )
	  {
	     if( !t ) {t = const_cast<reco::GenParticle*>(mcp);}

	     const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
	     for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
	       {
		  if( idr->isAvailable() ) 
		    {		       
		       const reco::GenParticleRef& genParticle = (*idr);
		       const reco::GenParticle *d = genParticle.get();
		       reco::GenParticle *pf = getUnique(d,0);
		       
//		       if( pf->status() != 3 && pf->status() != 62 ) continue;
		       
		       if( fabs(pf->pdgId()) == 5 || fabs(pf->pdgId()) == 3 || fabs(pf->pdgId()) == 1 ) // b or s or d
			 {
			    if( !tb ) tb = pf;
			    if( !tb_IS ) tb_IS = const_cast<reco::GenParticle*>(d);
			 }		       
		       
		       if( fabs(pf->pdgId()) == 24 ) // W
			 {
			    if( !tW )
			      {
				 tW = pf;
				 const reco::GenParticleRefVector& tW_daughterRefs = tW->daughterRefVector();
				 for(reco::GenParticleRefVector::const_iterator tW_idr = tW_daughterRefs.begin();
				     tW_idr!= tW_daughterRefs.end(); ++tW_idr)
				   {
				      if( tW_idr->isAvailable() ) 
					{		       
					   const reco::GenParticleRef& tW_genParticle = (*tW_idr);
					   const reco::GenParticle *tW_d = tW_genParticle.get();
					   reco::GenParticle *pff = getUnique(tW_d,0);
					   
					   if( fabs(pff->pdgId()) == 12 ||
					       fabs(pff->pdgId()) == 14 ) // nu
					     {
						tWnu = pff;
					     }		
					   if( fabs(pff->pdgId()) == 16 ) // nu_tau
					     {
						tWnutau = pff;
					     }		
					   if( fabs(pff->pdgId()) == 11 ||
					       fabs(pff->pdgId()) == 13 ) // l
					     {
						tWl = pff;
					     }		
					   if( fabs(pff->pdgId()) == 15 ) // tau
					     {
						tWtau = pff;
						
						const reco::GenParticleRefVector& tWtau_daughterRefs = tWtau->daughterRefVector();
						for(reco::GenParticleRefVector::const_iterator tWtau_idr = tWtau_daughterRefs.begin();
						    tWtau_idr!= tWtau_daughterRefs.end(); ++tWtau_idr)
						  {
						     if( tWtau_idr->isAvailable() ) 
						       {		       
							  const reco::GenParticleRef& tWtau_genParticle = (*tWtau_idr);
							  const reco::GenParticle *tWtau_d = tWtau_genParticle.get();
							  reco::GenParticle *pfff = getUnique(tWtau_d,0);
							  
							  if( fabs(pfff->pdgId()) == 12 ||
							      fabs(pfff->pdgId()) == 14 ) // nu
							    {
							       tWtaunu = pfff;
							    }		
							  if( fabs(pfff->pdgId()) == 16 ) // nu_tau
							    {
							       tWtaunutau = pfff;
							    }			
							  if( fabs(pfff->pdgId()) == 11 ||
							      fabs(pfff->pdgId()) == 13 ) // l
							    {
							       tWtaul = pfff;
							    }							  
						       }
						  }						
					     }
					   if( fabs(pff->pdgId()) <= 6 ) // q
					     {
						if( tWq1 && !tWq2 ) {tWq2 = pff;}
						if( !tWq1 ) {tWq1 = pff;}
						if( tWq1_IS && !tWq2_IS ) tWq2_IS = const_cast<reco::GenParticle*>(tW_d);
						if( !tWq1_IS ) tWq1_IS = const_cast<reco::GenParticle*>(tW_d);
					     }
					}
				   }				
			      }			    
			 }		       
		    }
	       }	     
	  }
     }

   bool doCheck = 0;
   
   if( t && tb && tW )
     {	
	int tchan = -666;
	if( tWl )   tchan = 0;
	if( tWq1 && tWq2 ) tchan = 1;
	if( tWtaul )  tchan = 2;
	if( tWtaunutau && !tWtaul )  tchan = 3;
	
	if( tchan < 0 && doCheck )
	  {	     
	     std::cout << "Failed to identify top-quark decay chain" << std::endl;
	     
	     std::cout << "t = " << bool(t) << std::endl;
	     std::cout << "t->W = " << bool(tW) << std::endl;
	     std::cout << "t->W->l = " << bool(tWl) << std::endl;
	     std::cout << "t->W->nu = " << bool(tWnu) << std::endl;
	     std::cout << "t->W->nutau = " << bool(tWnutau) << std::endl;
	     std::cout << "t->W->tau = " << bool(tWtau) << std::endl;
	     std::cout << "t->W->tau->l = " << bool(tWtaul) << std::endl;
	     std::cout << "t->W->tau->nu = " << bool(tWtaunu) << std::endl;
	     std::cout << "t->W->tau->nutau = " << bool(tWtaunutau) << std::endl;
	     std::cout << "t->W->q = " << bool(tWq1) << std::endl;

	     exit(1);
	  }
	
	if( Zl1 || Ztaul1 || Ztaunutau1 )
	  {
	     int chan0 = 0;
	     if( Zl1 ) chan = chan0 + 0 + tchan;
	     if( Ztaul1 && Ztaul2 ) chan = chan0 + 20 + tchan;
	     if( Ztaunutau1 && ! Ztaul1 && Ztaunutau2 && ! Ztaul2 ) chan = chan0 + 40 + tchan;
	     if( Ztaul1 && Ztaunutau2 && ! Ztaul2 ) chan = chan0 + 60 + tchan;
	     if( Ztaul2 && Ztaunutau1 && ! Ztaul1 ) chan = chan0 + 80 + tchan;
	  }	

	if( chan < 0 && doCheck )
	  {
	     std::cout << "Unknown channel found" << std::endl;
	     std::cout << "chan = " << chan << std::endl;

	     std::cout << "j1 = " << bool(j1) << std::endl;
	     std::cout << "j2 = " << bool(j2) << std::endl;
	     std::cout << "j3 = " << bool(j3) << std::endl;
	     
	     std::cout << "Z = " << bool(Z) << std::endl;
	     std::cout << "Z->l1 = " << bool(Zl1) << std::endl;
	     std::cout << "Z->l2 = " << bool(Zl2) << std::endl;
	     std::cout << "Z->tau1 = " << bool(Ztau1) << std::endl;
	     std::cout << "Z->tau1->l = " << bool(Ztaul1) << std::endl;
	     std::cout << "Z->tau1->nu = " << bool(Ztaunu1) << std::endl;
	     std::cout << "Z->tau1->nutau = " << bool(Ztaunutau1) << std::endl;
	     std::cout << "Z->tau2 = " << bool(Ztau2) << std::endl;
	     std::cout << "Z->tau2->l = " << bool(Ztaul2) << std::endl;
	     std::cout << "Z->tau2->nu = " << bool(Ztaunu2) << std::endl;
	     std::cout << "Z->tau2->nutau = " << bool(Ztaunutau2) << std::endl;
	     
	     exit(1);
	  }
     }

   tree.mc_truth_tzq_channel = chan;

   // TLV

   if( Z ) p4toTLV(Z->p4(),tree.mc_truth_Z_p4);
   if( Zl1 ) p4toTLV(Zl1->p4(),tree.mc_truth_Zl1_p4);
   if( Zl2 ) p4toTLV(Zl2->p4(),tree.mc_truth_Zl2_p4);
   if( Ztau1 ) p4toTLV(Ztau1->p4(),tree.mc_truth_Ztau1_p4);
   if( Ztau2 ) p4toTLV(Ztau2->p4(),tree.mc_truth_Ztau2_p4);
   if( Ztaul1 ) p4toTLV(Ztaul1->p4(),tree.mc_truth_Ztaul1_p4);
   if( Ztaul2 ) p4toTLV(Ztaul2->p4(),tree.mc_truth_Ztaul2_p4);
   if( Ztaunu1 ) p4toTLV(Ztaunu1->p4(),tree.mc_truth_Ztaunu1_p4);
   if( Ztaunu2 ) p4toTLV(Ztaunu2->p4(),tree.mc_truth_Ztaunu2_p4);
   if( Ztaunutau1 ) p4toTLV(Ztaunutau1->p4(),tree.mc_truth_Ztaunutau1_p4);
   if( Ztaunutau2 ) p4toTLV(Ztaunutau2->p4(),tree.mc_truth_Ztaunutau2_p4);
   
   if( t ) p4toTLV(t->p4(),tree.mc_truth_t_p4);
   if( tb ) p4toTLV(tb->p4(),tree.mc_truth_tb_p4);
   if( tb_IS ) p4toTLV(tb_IS->p4(),tree.mc_truth_tb_IS_p4);
   if( tW ) p4toTLV(tW->p4(),tree.mc_truth_tW_p4);
   if( tWnu ) p4toTLV(tWnu->p4(),tree.mc_truth_tWnu_p4);
   if( tWnutau ) p4toTLV(tWnutau->p4(),tree.mc_truth_tWnutau_p4);
   if( tWl ) p4toTLV(tWl->p4(),tree.mc_truth_tWl_p4);
   if( tWtau ) p4toTLV(tWtau->p4(),tree.mc_truth_tWtau_p4);
   if( tWtaunu ) p4toTLV(tWtaunu->p4(),tree.mc_truth_tWtaunu_p4);
   if( tWtaunutau ) p4toTLV(tWtaunutau->p4(),tree.mc_truth_tWtaunutau_p4);
   if( tWtaul ) p4toTLV(tWtaul->p4(),tree.mc_truth_tWtaul_p4);
   if( tWq1 ) p4toTLV(tWq1->p4(),tree.mc_truth_tWq1_p4);
   if( tWq2 ) p4toTLV(tWq2->p4(),tree.mc_truth_tWq2_p4);
   if( tWq1_IS ) p4toTLV(tWq1_IS->p4(),tree.mc_truth_tWq1_IS_p4);
   if( tWq2_IS ) p4toTLV(tWq2_IS->p4(),tree.mc_truth_tWq2_IS_p4);

   if( j1 ) p4toTLV(j1->p4(),tree.mc_truth_j1_p4);
   if( j2 ) p4toTLV(j2->p4(),tree.mc_truth_j2_p4);
   if( j3 ) p4toTLV(j3->p4(),tree.mc_truth_j3_p4);

   // pt

   if( Z ) tree.mc_truth_Z_pt = Z->p4().pt();
   if( Zl1 ) tree.mc_truth_Zl1_pt = Zl1->p4().pt();
   if( Zl2 ) tree.mc_truth_Zl2_pt = Zl2->p4().pt();
   if( Ztau1 ) tree.mc_truth_Ztau1_pt = Ztau1->p4().pt();
   if( Ztau2 ) tree.mc_truth_Ztau2_pt = Ztau2->p4().pt();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_pt = Ztaul1->p4().pt();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_pt = Ztaul2->p4().pt();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_pt = Ztaunu1->p4().pt();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_pt = Ztaunu2->p4().pt();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_pt = Ztaunutau1->p4().pt();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_pt = Ztaunutau2->p4().pt();
   
   if( t ) tree.mc_truth_t_pt = t->p4().pt();
   if( tb ) tree.mc_truth_tb_pt = tb->p4().pt();
   if( tb_IS ) tree.mc_truth_tb_IS_pt = tb_IS->p4().pt();
   if( tW ) tree.mc_truth_tW_pt = tW->p4().pt();
   if( tWnu ) tree.mc_truth_tWnu_pt = tWnu->p4().pt();
   if( tWnutau ) tree.mc_truth_tWnutau_pt = tWnutau->p4().pt();
   if( tWl ) tree.mc_truth_tWl_pt = tWl->p4().pt();
   if( tWtau ) tree.mc_truth_tWtau_pt = tWtau->p4().pt();
   if( tWtaunu ) tree.mc_truth_tWtaunu_pt = tWtaunu->p4().pt();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_pt = tWtaunutau->p4().pt();
   if( tWtaul ) tree.mc_truth_tWtaul_pt = tWtaul->p4().pt();
   if( tWq1 ) tree.mc_truth_tWq1_pt = tWq1->p4().pt();
   if( tWq2 ) tree.mc_truth_tWq2_pt = tWq2->p4().pt();
   if( tWq1_IS ) tree.mc_truth_tWq1_IS_pt = tWq1_IS->p4().pt();
   if( tWq2_IS ) tree.mc_truth_tWq2_IS_pt = tWq2_IS->p4().pt();
   
   if( j1 ) tree.mc_truth_j1_pt = j1->p4().pt();
   if( j2 ) tree.mc_truth_j2_pt = j2->p4().pt();
   if( j3 ) tree.mc_truth_j3_pt = j3->p4().pt();

   // eta

   if( Z ) tree.mc_truth_Z_eta = Z->p4().eta();
   if( Zl1 ) tree.mc_truth_Zl1_eta = Zl1->p4().eta();
   if( Zl2 ) tree.mc_truth_Zl2_eta = Zl2->p4().eta();
   if( Ztau1 ) tree.mc_truth_Ztau1_eta = Ztau1->p4().eta();
   if( Ztau2 ) tree.mc_truth_Ztau2_eta = Ztau2->p4().eta();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_eta = Ztaul1->p4().eta();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_eta = Ztaul2->p4().eta();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_eta = Ztaunu1->p4().eta();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_eta = Ztaunu2->p4().eta();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_eta = Ztaunutau1->p4().eta();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_eta = Ztaunutau2->p4().eta();
   
   if( t ) tree.mc_truth_t_eta = t->p4().eta();
   if( tb ) tree.mc_truth_tb_eta = tb->p4().eta();
   if( tb_IS ) tree.mc_truth_tb_IS_eta = tb_IS->p4().eta();
   if( tW ) tree.mc_truth_tW_eta = tW->p4().eta();
   if( tWnu ) tree.mc_truth_tWnu_eta = tWnu->p4().eta();
   if( tWnutau ) tree.mc_truth_tWnutau_eta = tWnutau->p4().eta();
   if( tWl ) tree.mc_truth_tWl_eta = tWl->p4().eta();
   if( tWtau ) tree.mc_truth_tWtau_eta = tWtau->p4().eta();
   if( tWtaunu ) tree.mc_truth_tWtaunu_eta = tWtaunu->p4().eta();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_eta = tWtaunutau->p4().eta();
   if( tWtaul ) tree.mc_truth_tWtaul_eta = tWtaul->p4().eta();
   if( tWq1 ) tree.mc_truth_tWq1_eta = tWq1->p4().eta();
   if( tWq2 ) tree.mc_truth_tWq2_eta = tWq2->p4().eta();
   if( tWq1_IS ) tree.mc_truth_tWq1_IS_eta = tWq1_IS->p4().eta();
   if( tWq2_IS ) tree.mc_truth_tWq2_IS_eta = tWq2_IS->p4().eta();
   
   if( j1 ) tree.mc_truth_j1_eta = j1->p4().eta();
   if( j2 ) tree.mc_truth_j2_eta = j2->p4().eta();
   if( j3 ) tree.mc_truth_j3_eta = j3->p4().eta();
   
   // phi

   if( Z ) tree.mc_truth_Z_phi = Z->p4().phi();
   if( Zl1 ) tree.mc_truth_Zl1_phi = Zl1->p4().phi();
   if( Zl2 ) tree.mc_truth_Zl2_phi = Zl2->p4().phi();
   if( Ztau1 ) tree.mc_truth_Ztau1_phi = Ztau1->p4().phi();
   if( Ztau2 ) tree.mc_truth_Ztau2_phi = Ztau2->p4().phi();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_phi = Ztaul1->p4().phi();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_phi = Ztaul2->p4().phi();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_phi = Ztaunu1->p4().phi();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_phi = Ztaunu2->p4().phi();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_phi = Ztaunutau1->p4().phi();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_phi = Ztaunutau2->p4().phi();
   
   if( t ) tree.mc_truth_t_phi = t->p4().phi();
   if( tb ) tree.mc_truth_tb_phi = tb->p4().phi();
   if( tb_IS ) tree.mc_truth_tb_IS_phi = tb_IS->p4().phi();
   if( tW ) tree.mc_truth_tW_phi = tW->p4().phi();
   if( tWnu ) tree.mc_truth_tWnu_phi = tWnu->p4().phi();
   if( tWnutau ) tree.mc_truth_tWnutau_phi = tWnutau->p4().phi();
   if( tWl ) tree.mc_truth_tWl_phi = tWl->p4().phi();
   if( tWtau ) tree.mc_truth_tWtau_phi = tWtau->p4().phi();
   if( tWtaunu ) tree.mc_truth_tWtaunu_phi = tWtaunu->p4().phi();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_phi = tWtaunutau->p4().phi();
   if( tWtaul ) tree.mc_truth_tWtaul_phi = tWtaul->p4().phi();
   if( tWq1 ) tree.mc_truth_tWq1_phi = tWq1->p4().phi();
   if( tWq2 ) tree.mc_truth_tWq2_phi = tWq2->p4().phi();
   
   if( j1 ) tree.mc_truth_j1_phi = j1->p4().phi();
   if( j2 ) tree.mc_truth_j2_phi = j2->p4().phi();
   if( j3 ) tree.mc_truth_j3_phi = j3->p4().phi();
    
   // E

   if( Z ) tree.mc_truth_Z_E = Z->p4().E();
   if( Zl1 ) tree.mc_truth_Zl1_E = Zl1->p4().E();
   if( Zl2 ) tree.mc_truth_Zl2_E = Zl2->p4().E();
   if( Ztau1 ) tree.mc_truth_Ztau1_E = Ztau1->p4().E();
   if( Ztau2 ) tree.mc_truth_Ztau2_E = Ztau2->p4().E();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_E = Ztaul1->p4().E();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_E = Ztaul2->p4().E();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_E = Ztaunu1->p4().E();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_E = Ztaunu2->p4().E();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_E = Ztaunutau1->p4().E();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_E = Ztaunutau2->p4().E();
   
   if( t ) tree.mc_truth_t_E = t->p4().E();
   if( tb ) tree.mc_truth_tb_E = tb->p4().E();
   if( tb_IS ) tree.mc_truth_tb_IS_E = tb_IS->p4().E();
   if( tW ) tree.mc_truth_tW_E = tW->p4().E();
   if( tWnu ) tree.mc_truth_tWnu_E = tWnu->p4().E();
   if( tWnutau ) tree.mc_truth_tWnutau_E = tWnutau->p4().E();
   if( tWl ) tree.mc_truth_tWl_E = tWl->p4().E();
   if( tWtau ) tree.mc_truth_tWtau_E = tWtau->p4().E();
   if( tWtaunu ) tree.mc_truth_tWtaunu_E = tWtaunu->p4().E();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_E = tWtaunutau->p4().E();
   if( tWtaul ) tree.mc_truth_tWtaul_E = tWtaul->p4().E();
   if( tWq1 ) tree.mc_truth_tWq1_E = tWq1->p4().E();
   if( tWq2 ) tree.mc_truth_tWq2_E = tWq2->p4().E();
   if( tWq1_IS ) tree.mc_truth_tWq1_IS_E = tWq1_IS->p4().E();
   if( tWq2_IS ) tree.mc_truth_tWq2_IS_E = tWq2_IS->p4().E();
   
   if( j1 ) tree.mc_truth_j1_E = j1->p4().E();
   if( j2 ) tree.mc_truth_j2_E = j2->p4().E();
   if( j3 ) tree.mc_truth_j3_E = j3->p4().E();
  
   // pdgId

   if( Z ) tree.mc_truth_Z_id = Z->pdgId();
   if( Zl1 ) tree.mc_truth_Zl1_id = Zl1->pdgId();
   if( Zl2 ) tree.mc_truth_Zl2_id = Zl2->pdgId();
   if( Ztau1 ) tree.mc_truth_Ztau1_id = Ztau1->pdgId();
   if( Ztau2 ) tree.mc_truth_Ztau2_id = Ztau2->pdgId();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_id = Ztaul1->pdgId();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_id = Ztaul2->pdgId();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_id = Ztaunu1->pdgId();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_id = Ztaunu2->pdgId();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_id = Ztaunutau1->pdgId();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_id = Ztaunutau2->pdgId();
   
   if( t ) tree.mc_truth_t_id = t->pdgId();
   if( tb ) tree.mc_truth_tb_id = tb->pdgId();
   if( tb_IS ) tree.mc_truth_tb_IS_id = tb_IS->pdgId();
   if( tW ) tree.mc_truth_tW_id = tW->pdgId();
   if( tWnu ) tree.mc_truth_tWnu_id = tWnu->pdgId();
   if( tWnutau ) tree.mc_truth_tWnutau_id = tWnutau->pdgId();
   if( tWl ) tree.mc_truth_tWl_id = tWl->pdgId();
   if( tWtau ) tree.mc_truth_tWtau_id = tWtau->pdgId();
   if( tWtaunu ) tree.mc_truth_tWtaunu_id = tWtaunu->pdgId();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_id = tWtaunutau->pdgId();
   if( tWtaul ) tree.mc_truth_tWtaul_id = tWtaul->pdgId();
   if( tWq1 ) tree.mc_truth_tWq1_id = tWq1->pdgId();
   if( tWq2 ) tree.mc_truth_tWq2_id = tWq2->pdgId();
   if( tWq1_IS ) tree.mc_truth_tWq1_IS_id = tWq1_IS->pdgId();
   if( tWq2_IS ) tree.mc_truth_tWq2_IS_id = tWq2_IS->pdgId();
   
   if( j1 ) tree.mc_truth_j1_id = j1->pdgId();
   if( j2 ) tree.mc_truth_j2_id = j2->pdgId();
   if( j3 ) tree.mc_truth_j3_id = j3->pdgId();

   // status
   
   if( Z ) tree.mc_truth_Z_status = Z->status();
   if( Zl1 ) tree.mc_truth_Zl1_status = Zl1->status();
   if( Zl2 ) tree.mc_truth_Zl2_status = Zl2->status();
   if( Ztau1 ) tree.mc_truth_Ztau1_status = Ztau1->status();
   if( Ztau2 ) tree.mc_truth_Ztau2_status = Ztau2->status();
   if( Ztaul1 ) tree.mc_truth_Ztaul1_status = Ztaul1->status();
   if( Ztaul2 ) tree.mc_truth_Ztaul2_status = Ztaul2->status();
   if( Ztaunu1 ) tree.mc_truth_Ztaunu1_status = Ztaunu1->status();
   if( Ztaunu2 ) tree.mc_truth_Ztaunu2_status = Ztaunu2->status();
   if( Ztaunutau1 ) tree.mc_truth_Ztaunutau1_status = Ztaunutau1->status();
   if( Ztaunutau2 ) tree.mc_truth_Ztaunutau2_status = Ztaunutau2->status();
   
   if( t ) tree.mc_truth_t_status = t->status();
   if( tb ) tree.mc_truth_tb_status = tb->status();
   if( tb_IS ) tree.mc_truth_tb_IS_status = tb_IS->status();
   if( tW ) tree.mc_truth_tW_status = tW->status();
   if( tWnu ) tree.mc_truth_tWnu_status = tWnu->status();
   if( tWnutau ) tree.mc_truth_tWnutau_status = tWnutau->status();
   if( tWl ) tree.mc_truth_tWl_status = tWl->status();
   if( tWtau ) tree.mc_truth_tWtau_status = tWtau->status();
   if( tWtaunu ) tree.mc_truth_tWtaunu_status = tWtaunu->status();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_status = tWtaunutau->status();
   if( tWtaul ) tree.mc_truth_tWtaul_status = tWtaul->status();
   if( tWq1 ) tree.mc_truth_tWq1_status = tWq1->status();
   if( tWq2 ) tree.mc_truth_tWq2_status = tWq2->status();
   if( tWq1_IS ) tree.mc_truth_tWq1_IS_status = tWq1_IS->status();
   if( tWq2_IS ) tree.mc_truth_tWq2_IS_status = tWq2_IS->status();
   
   if( j1 ) tree.mc_truth_j1_status = j1->status();
   if( j2 ) tree.mc_truth_j2_status = j2->status();
   if( j3 ) tree.mc_truth_j3_status = j3->status();
}

void MCTruth::fillTHQSignalGenParticles(const edm::Event& iEvent,
					const edm::EventSetup& iSetup,
					FlatTree& tree,
					const edm::Handle<std::vector<reco::GenParticle> >& GenParticles)
{
   reco::GenParticle *h0 = 0;
   
   reco::GenParticle *h0W1 = 0;
   reco::GenParticle *h0W2 = 0;
   reco::GenParticle *h0Wl1 = 0;
   reco::GenParticle *h0Wnu1 = 0;
   reco::GenParticle *h0Wtau1 = 0;
   reco::GenParticle *h0Wnutau1 = 0;
   reco::GenParticle *h0Wtaul1 = 0;
   reco::GenParticle *h0Wtaunu1 = 0;
   reco::GenParticle *h0Wtaunutau1 = 0;
   reco::GenParticle *h0Wl2 = 0;
   reco::GenParticle *h0Wnu2 = 0;
   reco::GenParticle *h0Wtau2 = 0;
   reco::GenParticle *h0Wnutau2 = 0;
   reco::GenParticle *h0Wtaul2 = 0;
   reco::GenParticle *h0Wtaunu2 = 0;
   reco::GenParticle *h0Wtaunutau2 = 0;
   reco::GenParticle *h0Wq11 = 0;
   reco::GenParticle *h0Wq21 = 0;
   reco::GenParticle *h0Wq12 = 0;
   reco::GenParticle *h0Wq22 = 0;
   reco::GenParticle *h0Wq11_IS = 0;
   reco::GenParticle *h0Wq21_IS = 0;
   reco::GenParticle *h0Wq12_IS = 0;
   reco::GenParticle *h0Wq22_IS = 0;

   reco::GenParticle *h0Z1 = 0;
   reco::GenParticle *h0Z2 = 0;
   reco::GenParticle *h0Zl11 = 0;
   reco::GenParticle *h0Zl21 = 0;
   reco::GenParticle *h0Ztau11 = 0;
   reco::GenParticle *h0Ztau21 = 0;
   reco::GenParticle *h0Ztaul11 = 0;
   reco::GenParticle *h0Ztaul21 = 0;
   reco::GenParticle *h0Ztaunu11 = 0;
   reco::GenParticle *h0Ztaunu21 = 0;
   reco::GenParticle *h0Ztaunutau11 = 0;
   reco::GenParticle *h0Ztaunutau21 = 0;
   reco::GenParticle *h0Zq11 = 0;
   reco::GenParticle *h0Zq21 = 0;
   reco::GenParticle *h0Zq11_IS = 0;
   reco::GenParticle *h0Zq21_IS = 0;
   reco::GenParticle *h0Zl12 = 0;
   reco::GenParticle *h0Zl22 = 0;
   reco::GenParticle *h0Ztau12 = 0;
   reco::GenParticle *h0Ztau22 = 0;
   reco::GenParticle *h0Ztaul12 = 0;
   reco::GenParticle *h0Ztaul22 = 0;
   reco::GenParticle *h0Ztaunu12 = 0;
   reco::GenParticle *h0Ztaunu22 = 0;
   reco::GenParticle *h0Ztaunutau12 = 0;
   reco::GenParticle *h0Ztaunutau22 = 0;
   reco::GenParticle *h0Zq12 = 0;
   reco::GenParticle *h0Zq22 = 0;
   reco::GenParticle *h0Zq12_IS = 0;
   reco::GenParticle *h0Zq22_IS = 0;
   reco::GenParticle *h0Znu11 = 0;
   reco::GenParticle *h0Znu21 = 0;
   reco::GenParticle *h0Znu12 = 0;
   reco::GenParticle *h0Znu22 = 0;
   
   reco::GenParticle *h0tau1 = 0;
   reco::GenParticle *h0tau2 = 0;
   reco::GenParticle *h0taul1 = 0;
   reco::GenParticle *h0taunutau1 = 0;
   reco::GenParticle *h0taunu1 = 0;
   reco::GenParticle *h0taul2 = 0;
   reco::GenParticle *h0taunutau2 = 0;
   reco::GenParticle *h0taunu2 = 0;

   reco::GenParticle *h0b1 = 0;
   reco::GenParticle *h0b2 = 0;
   reco::GenParticle *h0b1_IS = 0;
   reco::GenParticle *h0b2_IS = 0;
   
   reco::GenParticle *t = 0;

   reco::GenParticle *tb = 0;
   reco::GenParticle *tb_IS = 0;
   
   reco::GenParticle *tW = 0;
   reco::GenParticle *tWnu = 0;
   reco::GenParticle *tWnutau = 0;
   reco::GenParticle *tWl = 0;
   reco::GenParticle *tWtau = 0;
   reco::GenParticle *tWtaunu = 0;
   reco::GenParticle *tWtaunutau = 0;
   reco::GenParticle *tWtaul = 0;
   reco::GenParticle *tWq1 = 0;
   reco::GenParticle *tWq2 = 0;
   reco::GenParticle *tWq1_IS = 0;
   reco::GenParticle *tWq2_IS = 0;

   reco::GenParticle *j1 = 0;
   reco::GenParticle *j2 = 0;
   reco::GenParticle *j3 = 0;
   
   int chan = -666;

   // 0   = (t->bW,W->lnu)
   // 1   = (t->bW,W->qq)
   // 2   = (t->bW,W->tauLnu)
   // 3   = (t->bW,W->tauHnu)
   
   // (H->WW:W->lnu,W->lnu)                 +0
   // (H->WW:W->tauLtaunu,W->tauLtaunu)     +20
   // (H->WW:W->tauHtaunu,W->tauHtaunu)     +40
   // (H->WW:W->tauHtaunu,W->tauLtaunu)     +60
   // (H->WW:W->tauLtaunu,W->tauHtaunu)     +80
   // (H->WW:W->qq,W->qq)                   +100
   // (H->WW:W->lnu,W->qq)                  +120
   // (H->WW:W->tauLnutau,W->qq)            +140
   // (H->WW:W->tauHnutau,W->qq)            +160
   // (H->WW:W->qq,W->lnu)                  +180
   // (H->WW:W->qq,W->tauLtaunu)            +200
   // (H->WW:W->qq,W->tauHtaunu)            +220
   // (H->WW:W->lnu,W->tauLtaunu)           +240
   // (H->WW:W->lnu,W->tauHtaunu)           +260
   // (H->WW:W->tauLtaunu,W->lnu)           +280
   // (H->WW:W->tauHtaunu,W->lnu)           +300
   
   // +1000
   // (H->ZZ:Z->ll,Z->ll)                   +0
   // (H->ZZ:Z->tauLtauL,Z->tauLtauL)       +20
   // (H->ZZ:Z->tauLtauL,Z->tauHtauH)       +40
   // (H->ZZ:Z->tauLtauL,Z->tauLtauH)       +60
   // (H->ZZ:Z->tauLtauL,Z->tauHtauL)       +80
   // (H->ZZ:Z->tauHtauH,Z->tauLtauL)       +100
   // (H->ZZ:Z->tauHtauH,Z->tauHtauH)       +120
   // (H->ZZ:Z->tauHtauH,Z->tauLtauH)       +140
   // (H->ZZ:Z->tauHtauH,Z->tauHtauL)       +160
   // (H->ZZ:Z->tauLtauH,Z->tauLtauL)       +180
   // (H->ZZ:Z->tauLtauH,Z->tauHtauH)       +200
   // (H->ZZ:Z->tauLtauH,Z->tauLtauH)       +220
   // (H->ZZ:Z->tauLtauH,Z->tauHtauL)       +240
   // (H->ZZ:Z->tauHtauL,Z->tauLtauL)       +260
   // (H->ZZ:Z->tauHtauL,Z->tauHtauH)       +280
   // (H->ZZ:Z->tauHtauL,Z->tauLtauH)       +300
   // (H->ZZ:Z->tauHtauL,Z->tauHtauL)       +320
   // (H->ZZ:Z->qq,Z->qq)                   +340
   // (H->ZZ:Z->ll,Z->qq)                   +360
   // (H->ZZ:Z->tauHtauH,Z->qq)             +380
   // (H->ZZ:Z->tauLtauL,Z->qq)             +400
   // (H->ZZ:Z->tauHtauL,Z->qq)             +420
   // (H->ZZ:Z->tauLtauH,Z->qq)             +440
   // (H->ZZ:Z->qq,Z->ll)                   +460
   // (H->ZZ:Z->qq,Z->tauHtauH)             +480
   // (H->ZZ:Z->qq,Z->tauLtauL)             +500
   // (H->ZZ:Z->qq,Z->tauHtauL)             +520
   // (H->ZZ:Z->qq,Z->tauLtauH)             +540
   // (H->ZZ:Z->ll,Z->nunu)                 +560
   // (H->ZZ:Z->tauHtauH,Z->nunu)           +580
   // (H->ZZ:Z->tauLtauL,Z->nunu)           +600
   // (H->ZZ:Z->tauHtauL,Z->nunu)           +620
   // (H->ZZ:Z->tauLtauH,Z->nunu)           +640
   // (H->ZZ:Z->qq,Z->nunu)                 +660
   // (H->ZZ:Z->nunu,Z->qq)                 +680
   // (H->ZZ:Z->nunu,Z->ll)                 +700
   // (H->ZZ:Z->nunu,Z->tauHtauH)           +720
   // (H->ZZ:Z->nunu,Z->tauLtauL)           +740
   // (H->ZZ:Z->nunu,Z->tauHtauL)           +760
   // (H->ZZ:Z->nunu,Z->tauLtauH)           +780
   // (H->ZZ:Z->nunu,Z->nunu)               +800
   // (H->ZZ:Z->tauHtauH,Z->ll)             +820
   // (H->ZZ:Z->tauLtauL,Z->ll)             +840
   // (H->ZZ:Z->tauHtauL,Z->ll)             +860
   // (H->ZZ:Z->tauLtauH,Z->ll)             +880
   // (H->ZZ:Z->ll,Z->tauHtauH)             +900
   // (H->ZZ:Z->ll,Z->tauLtauL)             +920
   // (H->ZZ:Z->ll,Z->tauHtauL)             +940
   // (H->ZZ:Z->ll,Z->tauLtauH)             +960
   
   // +2000
   // (H->tautau:tauL,tauL)     +2020
   // (H->tautau:tauH,tauH)     +2040
   // (H->tautau:tauL,tauH)     +2060
   // (H->tautau:tauH,tauL)     +2080

   // +3000
   // (H->bbbar)     +3020

   reco::GenParticleCollection genParticlesCollection = *GenParticles;
   reco::GenParticleCollection::const_iterator genParticleSrc;
   
   int ipart = 0;
   
   for(genParticleSrc = genParticlesCollection.begin();
       genParticleSrc != genParticlesCollection.end(); 
       genParticleSrc++)
     {
	reco::GenParticle *mcp = &(const_cast<reco::GenParticle&>(*genParticleSrc));

	int barcode = ipart; // in CMSSW barcode is the index of genParticle in the event
	// https://twiki.cern.ch/twiki/bin/view/CMS/GenParticles2HepMCConverter
	ipart++;
	
	// Additional partons (up to three)
	if( (fabs(mcp->pdgId()) <= 6 || fabs(mcp->pdgId()) == 21) &&
	    mcp->status() == 23 && barcode == 8 )
	  {
	     if( !j1 )
	       j1 = mcp;
	     else if( !j2 )
	       j2 = mcp;
	     else if( !j3 )
	       j3 = mcp;
	  }	

	// Higgs decays
	if( fabs(mcp->pdgId()) == 25 )
	  {
//	     if( iEvent.id().event() == 23539 )
//	       {
//		  std::cout << "found " << mcp->status() << " " << mcp->numberOfDaughters() << std::endl;
//	       }		  
	     	     	     
	     if( (mcp->status() == 62) || // higgs produced in initial state (MG)
		 (mcp->status() == 52) || // higgs comes from top decay (MG)
		 (mcp->status() == 3) )
	       {
		  h0 = const_cast<reco::GenParticle*>(mcp);
		  
		  const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
		  for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
		    {
		       if( idr->isAvailable() ) 
			 {
			    const reco::GenParticleRef& genParticle = (*idr);
			    const reco::GenParticle *d = genParticle.get();
			    reco::GenParticle *pf = getUnique(d,0);
			    reco::GenParticle *di = const_cast<reco::GenParticle*>(d);
			    
			    // h0 -> bbbar
			    if( fabs(pf->pdgId()) == 5 ) chan = 10000;
			    // h0 -> ee/mumu
			    if( fabs(pf->pdgId()) == 11 || fabs(pf->pdgId()) == 13 ) chan = 10001;
			    // h0 -> gg
			    if( fabs(pf->pdgId()) == 21 ) chan = 10002;
			    // h0 -> gammagamma
			    if( fabs(pf->pdgId()) == 22 ) chan = 10003;
			    // h0 -> qqbar (non-b)
			    if( fabs(pf->pdgId()) == 6 || fabs(pf->pdgId()) <= 4 ) chan = 10004;
			    
			    // h0 -> WW
			    if( fabs(pf->pdgId()) == 24 )
			      {
				 if( h0W1 && !h0W2 ) {h0W2 = pf;}
				 if( !h0W1 ) {h0W1 = pf;}				 

				 if( h0W1 && !h0W2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0W1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0W1_idr = daughterRefs.begin(); 
					  h0W1_idr!= daughterRefs.end(); ++h0W1_idr) 
					{
					   if( h0W1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0W1_idr);
						const reco::GenParticle *h0W1_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0W1_d,0);
						
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu
						  {
						     h0Wnu1 = pff;
						  }
						if( fabs(pff->pdgId()) == 16 ) // nutau
						  {
						     h0Wnutau1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     h0Wl1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     h0Wtau1 = pff;
						     
						     const reco::GenParticleRefVector& daughterRefs = h0Wtau1->daughterRefVector();
						     for(reco::GenParticleRefVector::const_iterator h0Wtau1_idr = daughterRefs.begin();
							 h0Wtau1_idr!= daughterRefs.end(); ++h0Wtau1_idr) 
						       {
							  if( h0Wtau1_idr->isAvailable() ) 
							    {		       
							       const reco::GenParticleRef& genParticle = (*h0Wtau1_idr);
							       const reco::GenParticle *h0Wtau1_d = genParticle.get();
							       reco::GenParticle *pfff = getUnique(h0Wtau1_d,0);
							       
							       if( fabs(pfff->pdgId()) == 12 ||
								   fabs(pfff->pdgId()) == 14 ) // nu
								 {
								    h0Wtaunu1 = pfff;
								 }
							       if( fabs(pfff->pdgId()) == 16 ) // nutau
								 {
								    h0Wtaunutau1 = pfff;
								 }		
							       if( fabs(pfff->pdgId()) == 11 ||
								   fabs(pfff->pdgId()) == 13 ) // l
								 {
								    h0Wtaul1 = pfff;  
								 }
							    }
						       }
						  }						
						if( fabs(pff->pdgId()) <= 6 )
						  {
						     if( h0Wq11 && !h0Wq21 ) h0Wq21 = pff;
						     if( !h0Wq11 ) h0Wq11 = pff;
						     if( h0Wq11_IS && !h0Wq21_IS ) h0Wq21_IS = const_cast<reco::GenParticle*>(h0W1_d);
						     if( !h0Wq11_IS ) h0Wq11_IS = const_cast<reco::GenParticle*>(h0W1_d);
						  }
					     }
					}				      
				   }
				 if( h0W2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0W2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0W2_idr = daughterRefs.begin(); 
					  h0W2_idr!= daughterRefs.end(); ++h0W2_idr) 
					{
					   if( h0W2_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0W2_idr);
						const reco::GenParticle *h0W2_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0W2_d,0);
						
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu
						  {
						     h0Wnu2 = pff;
						  }
						if( fabs(pff->pdgId()) == 16 ) // nutau
						  {
						     h0Wnutau2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     h0Wl2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     h0Wtau2 = pff;
						     
						     const reco::GenParticleRefVector& daughterRefs = h0Wtau2->daughterRefVector();
						     for(reco::GenParticleRefVector::const_iterator h0Wtau2_idr = daughterRefs.begin();
							 h0Wtau2_idr!= daughterRefs.end(); ++h0Wtau2_idr) 
						       {
							  if( h0Wtau2_idr->isAvailable() ) 
							    {		       
							       const reco::GenParticleRef& genParticle = (*h0Wtau2_idr);
							       const reco::GenParticle *h0Wtau2_d = genParticle.get();
							       reco::GenParticle *pfff = getUnique(h0Wtau2_d,0);
							       
							       if( fabs(pfff->pdgId()) == 12 ||
								   fabs(pfff->pdgId()) == 14 ) // nu
								 {
								    h0Wtaunu2 = pfff;
								 }
							       if( fabs(pfff->pdgId()) == 16 ) // nutau
								 {
								    h0Wtaunutau2 = pfff;
								 }		
							       if( fabs(pfff->pdgId()) == 11 ||
								   fabs(pfff->pdgId()) == 13 ) // l
								 {
								    h0Wtaul2 = pfff;  
								 }
							    }
						       }
						  }						
						if( fabs(pff->pdgId()) <= 6 )
						  {
						     if( h0Wq12 && !h0Wq22 ) {h0Wq22 = pff;}
						     if( !h0Wq12 ) {h0Wq12 = pff;}
						     if( h0Wq12_IS && !h0Wq22_IS ) h0Wq22_IS = const_cast<reco::GenParticle*>(h0W2_d);
						     if( !h0Wq12_IS ) h0Wq12_IS = const_cast<reco::GenParticle*>(h0W2_d);
						  }
					     }
					}				      				      
				   }				 
			      }

			    // h0 -> ZZ
			    if( fabs(pf->pdgId()) == 23 )
			      {
				 if( h0Z1 && !h0Z2 ) {h0Z2 = pf;}
				 if( !h0Z1 ) {h0Z1 = pf;}

				 if( h0Z1 && !h0Z2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0Z1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0Z1_idr = daughterRefs.begin(); 
					  h0Z1_idr!= daughterRefs.end(); ++h0Z1_idr) 
					{
					   if( h0Z1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0Z1_idr);
						const reco::GenParticle *h0Z1_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0Z1_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     if( h0Zl11 && !h0Zl21 ) {h0Zl21 = pff;}
						     if( !h0Zl11 ) {h0Zl11 = pff;}
						  }				
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     if( h0Ztau11 && !h0Ztau21 )
						       {
							  h0Ztau21 = pff;
							  
							  const reco::GenParticleRefVector& daughterRefs = h0Ztau21->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau21_idr = daughterRefs.begin(); 
							      h0Ztau21_idr!= daughterRefs.end(); ++h0Ztau21_idr) 
							    {
							       if( h0Ztau21_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau21_idr);
								    const reco::GenParticle *h0Ztau21_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau21_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu21 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau21 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul21 = pfff;
								      }
								 }
							    }							  
						       }
						     if( !h0Ztau11 )
						       {
							  h0Ztau11 = pff;

							  const reco::GenParticleRefVector& daughterRefs = h0Ztau11->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau11_idr = daughterRefs.begin(); 
							      h0Ztau11_idr!= daughterRefs.end(); ++h0Ztau11_idr) 
							    {
							       if( h0Ztau11_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau11_idr);
								    const reco::GenParticle *h0Ztau11_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau11_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu11 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau11 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul11 = pfff;
								      }
								 }
							    }							  
						       }						     
						  }
						if( fabs(pff->pdgId()) <= 6 ) // q
						  {
						     if( h0Zq11 && !h0Zq21 ) h0Zq21 = pff;
						     if( !h0Zq11 ) h0Zq11 = pff;
						     if( h0Zq11_IS && !h0Zq21_IS ) h0Zq21_IS = const_cast<reco::GenParticle*>(h0Z1_d);
						     if( !h0Zq11_IS ) h0Zq11_IS = const_cast<reco::GenParticle*>(h0Z1_d);
						  }				
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ||
						    fabs(pff->pdgId()) == 16 ) // nu
						  {
						     if( h0Znu11 && !h0Znu21 ) {h0Znu21 = pff;}
						     if( !h0Znu11 ) {h0Znu11 = pff;}
						  }						
					     }					   
					}				      
				   }
				 if( h0Z2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0Z2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0Z2_idr = daughterRefs.begin();
					  h0Z2_idr!= daughterRefs.end(); ++h0Z2_idr) 
					{
					   if( h0Z2_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0Z2_idr);
						const reco::GenParticle *h0Z2_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0Z2_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ) // l
						  {
						     if( h0Zl12 && !h0Zl22 ) {h0Zl22 = pff;}
						     if( !h0Zl12 ) {h0Zl12 = pff;}
						  }				
						if( fabs(pff->pdgId()) == 15 ) // tau
						  {
						     if( h0Ztau12 && !h0Ztau22 )
						       {
							  h0Ztau22 = pff;
							  
							  const reco::GenParticleRefVector& daughterRefs = h0Ztau22->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau22_idr = daughterRefs.begin(); 
							      h0Ztau22_idr!= daughterRefs.end(); ++h0Ztau22_idr) 
							    {
							       if( h0Ztau22_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau22_idr);
								    const reco::GenParticle *h0Ztau22_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau22_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu22 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau22 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul22 = pfff;
								      }
								 }
							    }							  
						       }
						     if( !h0Ztau12 )
						       {
							  h0Ztau12 = pff;

							  const reco::GenParticleRefVector& daughterRefs = h0Ztau12->daughterRefVector();
							  for(reco::GenParticleRefVector::const_iterator h0Ztau12_idr = daughterRefs.begin(); 
							      h0Ztau12_idr!= daughterRefs.end(); ++h0Ztau12_idr)
							    {
							       if( h0Ztau12_idr->isAvailable() ) 
								 {		       
								    const reco::GenParticleRef& genParticle = (*h0Ztau12_idr);
								    const reco::GenParticle *h0Ztau12_d = genParticle.get();
								    reco::GenParticle *pfff = getUnique(h0Ztau12_d,0);
								    
								    if( fabs(pfff->pdgId()) == 12 ||
									fabs(pfff->pdgId()) == 14 ) // nu
								      {
									 h0Ztaunu12 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 16 ) // nutau
								      {
									 h0Ztaunutau12 = pfff;
								      }
								    if( fabs(pfff->pdgId()) == 11 ||
									fabs(pfff->pdgId()) == 13 ) // l
								      {
									 h0Ztaul12 = pfff;
								      }
								 }
							    }							  
						       }						     
						  }
						if( fabs(pff->pdgId()) <= 6 ) // q
						  {
						     if( h0Zq12 && !h0Zq22 ) h0Zq22 = pff;
						     if( !h0Zq12 ) h0Zq12 = pff;
						     if( h0Zq12_IS && !h0Zq22_IS ) h0Zq22_IS = const_cast<reco::GenParticle*>(h0Z2_d);
						     if( !h0Zq12_IS ) h0Zq12_IS = const_cast<reco::GenParticle*>(h0Z2_d);
						  }				
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ||
						    fabs(pff->pdgId()) == 16 ) // nu
						  {
						     if( h0Znu12 && !h0Znu22 ) {h0Znu22 = pff;}
						     if( !h0Znu12 ) {h0Znu12 = pff;}
						  }						
					     }					   
					}				      				      
				   }				 
			      }			    
			    
			    // h0 -> tautau
			    if( fabs(pf->pdgId()) == 15 && pf->status() == 2 ) // tau to decay
			      {
				 if( h0tau1 && !h0tau2 ) {h0tau2 = pf;}
				 if( !h0tau1 ) {h0tau1 = pf;}

				 if( h0tau1 && !h0tau2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0tau1->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0tau1_idr = daughterRefs.begin();
					  h0tau1_idr!= daughterRefs.end(); ++h0tau1_idr) 
					{
					   if( h0tau1_idr->isAvailable() ) 
					     {		       
						const reco::GenParticleRef& genParticle = (*h0tau1_idr);
						const reco::GenParticle *h0tau1_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0tau1_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ||
						    fabs(pff->pdgId()) == 15 ) // l
						  {
						     h0taul1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 16 ) // nu_tau
						  {
						     h0taunutau1 = pff;
						  }		
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu_e or nu_mu
						  {
						     h0taunu1 = pff;
						  }						
					     }
					}				      				      
				   }
				 if( h0tau2 )
				   {
				      const reco::GenParticleRefVector& daughterRefs = h0tau2->daughterRefVector();
				      for(reco::GenParticleRefVector::const_iterator h0tau2_idr = daughterRefs.begin();
					  h0tau2_idr!= daughterRefs.end(); ++h0tau2_idr) 
					{
					   if( h0tau2_idr->isAvailable() )
					     {		       
						const reco::GenParticleRef& genParticle = (*h0tau2_idr);
						const reco::GenParticle *h0tau2_d = genParticle.get();
						reco::GenParticle *pff = getUnique(h0tau2_d,0);
						
						if( fabs(pff->pdgId()) == 11 ||
						    fabs(pff->pdgId()) == 13 ||
						    fabs(pff->pdgId()) == 15 ) // l
						  {
						     h0taul2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 16 ) // nu_tau
						  {
						     h0taunutau2 = pff;
						  }		
						if( fabs(pff->pdgId()) == 12 ||
						    fabs(pff->pdgId()) == 14 ) // nu_e or nu_mu
						  {
						     h0taunu2 = pff;
						  }						
					     }
					}				      				      
				   }				 
			      }			    
			    
			    // h0 -> bbbar (FS)
			    if( fabs(pf->pdgId()) == 5 )
			      {
				 if( h0b1 && !h0b2 ) {h0b2 = pf;}
				 if( !h0b1 ) {h0b1 = pf;}
			      }			    
			    // h0 -> bbbar (IS)
			    if( fabs(di->pdgId()) == 5 && di->status() == 23 )
			      {
				 if( h0b1_IS && !h0b2_IS ) 
				   {
				      h0b2_IS = di;
				   }
				 if( !h0b1_IS ) 
				   {
				      h0b1_IS = di;
				   }
			      }			    
			 }		       
		    }		  
	       }	     
	  }	

	// top decays
	if( fabs(mcp->pdgId()) == 6
	    && ( (mcp->status() == 62) || 
		 (mcp->status() == 3)
	       ) )
	  {
	     t = const_cast<reco::GenParticle*>(mcp);	     	
	     
	     const reco::GenParticleRefVector& daughterRefs = mcp->daughterRefVector();
	     for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) 
	       {
		  if( idr->isAvailable() ) 
		    {		       
		       const reco::GenParticleRef& genParticle = (*idr);
		       const reco::GenParticle *d = genParticle.get();
		       reco::GenParticle *pf = getUnique(d,0);
//		       reco::GenParticle *di = const_cast<reco::GenParticle*>(d);

//		       if( pf->status() != 3 && pf->status() != 62 ) continue;
		       
		       if( fabs(pf->pdgId()) == 5 || fabs(pf->pdgId()) == 3 || fabs(pf->pdgId()) == 1 ) // b or s or d
			 {
			    tb = pf;
			    if( !tb_IS ) tb_IS = const_cast<reco::GenParticle*>(d);
			 }		       
/*		       if( (fabs(di->pdgId()) == 5 || fabs(di->pdgId()) == 3 || fabs(di->pdgId()) == 1) &&
			   di->status() == 23 ) // b or s or d
			 {
			    tb_IS = di;
			 }*/

		       if( fabs(pf->pdgId()) == 24 ) // W
			 {
			    tW = pf;
			    const reco::GenParticleRefVector& tW_daughterRefs = tW->daughterRefVector();
			    for(reco::GenParticleRefVector::const_iterator tW_idr = tW_daughterRefs.begin();
				tW_idr!= tW_daughterRefs.end(); ++tW_idr) 
			      {
				 if( tW_idr->isAvailable() ) 
				   {		       
				      const reco::GenParticleRef& tW_genParticle = (*tW_idr);
				      const reco::GenParticle *tW_d = tW_genParticle.get();
				      reco::GenParticle *pff = getUnique(tW_d,0);
					   
				      if( fabs(pff->pdgId()) == 12 ||
					  fabs(pff->pdgId()) == 14 ) // nu
					{
					   tWnu = pff;
					}		
				      if( fabs(pff->pdgId()) == 16 ) // nu_tau
					{
					   tWnutau = pff;
					}		
				      if( fabs(pff->pdgId()) == 11 ||
					  fabs(pff->pdgId()) == 13 ) // l
					{
					   tWl = pff;
					}		
				      if( fabs(pff->pdgId()) == 15 ) // tau
					{
					   tWtau = pff;
						
					   const reco::GenParticleRefVector& tWtau_daughterRefs = tWtau->daughterRefVector();
					   for(reco::GenParticleRefVector::const_iterator tWtau_idr = tWtau_daughterRefs.begin();
					       tWtau_idr!= tWtau_daughterRefs.end(); ++tWtau_idr)
					     {
						if( tWtau_idr->isAvailable() ) 
						  {		       
						     const reco::GenParticleRef& tWtau_genParticle = (*tWtau_idr);
						     const reco::GenParticle *tWtau_d = tWtau_genParticle.get();
						     reco::GenParticle *pfff = getUnique(tWtau_d,0);
							  
						     if( fabs(pfff->pdgId()) == 12 ||
							 fabs(pfff->pdgId()) == 14 ) // nu
						       {
							  tWtaunu = pfff;
						       }		
						     if( fabs(pfff->pdgId()) == 16 ) // nu_tau
						       {
							  tWtaunutau = pfff;
						       }			
						     if( fabs(pfff->pdgId()) == 11 ||
							 fabs(pfff->pdgId()) == 13 ) // l
						       {
							  tWtaul = pfff;
						       }							  
						  }						
					     }
					}
				      
				      if( fabs(pff->pdgId()) <= 6 ) // q
					{
					   if( tWq1 && !tWq2 ) tWq2 = pff;
					   if( !tWq1 ) tWq1 = pff;
					   if( tWq1_IS && !tWq2_IS ) tWq2_IS = const_cast<reco::GenParticle*>(tW_d);
					   if( !tWq1_IS ) tWq1_IS = const_cast<reco::GenParticle*>(tW_d);
					}					   					   
				   }				
			      }			    
			 }		       
		    }		  
	       }
	  }	
     }   
   
   bool doCheck = 0;

   if( h0 && t && tb && tW )
     {	
	int tchan = -666;
	if( tWl )   tchan = 0;
	if( tWq1 ) tchan = 1;
	if( tWtaul )  tchan = 2;
	if( tWtaunutau )  tchan = 3;
	
	if( tchan < 0 && doCheck )
	  {	     
	     std::cout << "Failed to identify top-quark decay chain" << std::endl;
	     
	     std::cout << "t = " << bool(t) << std::endl;
	     std::cout << "t->W = " << bool(tW) << std::endl;
	     std::cout << "t->W->l = " << bool(tWl) << std::endl;
	     std::cout << "t->W->nu = " << bool(tWnu) << std::endl;
	     std::cout << "t->W->nutau = " << bool(tWnutau) << std::endl;
	     std::cout << "t->W->tau = " << bool(tWtau) << std::endl;
	     std::cout << "t->W->tau->l = " << bool(tWtaul) << std::endl;
	     std::cout << "t->W->tau->nu = " << bool(tWtaunu) << std::endl;
	     std::cout << "t->W->tau->nutau = " << bool(tWtaunutau) << std::endl;
	     std::cout << "t->W->q = " << bool(tWq1) << std::endl;
	             
	     exit(1);
	  }
	
	if( h0W1 && h0W2 )
	  {	              
	     int chan0 = 0;
	     if( h0Wl1 && h0Wl2 ) chan = chan0 + 0 + tchan;
	     if( h0Wtaul1 && h0Wtaul2 ) chan = chan0 + 20 + tchan;
	     if( h0Wtaunutau1 && !h0Wtaul1 && h0Wtaunutau2 && !h0Wtaul2 ) chan = chan0 + 40 + tchan;
	     if( h0Wtaunutau1 && !h0Wtaul1 && h0Wtaul2 ) chan = chan0 + 60 + tchan;
	     if( h0Wtaul1 && h0Wtaunutau2 && !h0Wtaul2 ) chan = chan0 + 80 + tchan;
	     if( h0Wq11 && h0Wq12 ) chan = chan0 + 100 + tchan;
	     if( h0Wl1 && h0Wq12 ) chan = chan0 + 120 + tchan;
	     if( h0Wtaul1 && h0Wq12 ) chan = chan0 + 140 + tchan;
	     if( h0Wtaunutau1 && !h0Wtaul1 && h0Wq12 ) chan = chan0 + 160 + tchan;
	     if( h0Wq11 && h0Wl2 ) chan = chan0 + 180 + tchan;
	     if( h0Wq11 && h0Wtaul2 ) chan = chan0 + 200 + tchan;
	     if( h0Wq11 && h0Wtaunutau2 && !h0Wtaul2 ) chan = chan0 + 220 + tchan;
	     if( h0Wl1 && h0Wtaul2 ) chan = chan0 + 240 + tchan;
	     if( h0Wl1 && !h0Wtaul2 && h0Wtaunutau2 ) chan = chan0 + 260 + tchan;
	     if( h0Wl2 && h0Wtaul1 ) chan = chan0 + 280 + tchan;
	     if( h0Wl2 && !h0Wtaul1 && h0Wtaunutau1 ) chan = chan0 + 300 + tchan;
	  }

	if( h0Z1 && h0Z2 )
	  {
	     int chan0 = 1000;
	     if( h0Zl11 && h0Zl12 ) chan = chan0 + 0 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 20 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 40 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 60 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Ztaul22 && !h0Ztaul12 && h0Ztaunutau12 ) chan = chan0 + 80 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 100 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 120 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 140 + tchan;
	     if( !h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau11 && h0Ztaunutau21 && !h0Ztaul12 && h0Ztaul22 && h0Ztaunutau12 ) chan = chan0 + 160 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 180 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 200 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 220 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && !h0Ztaul12 && h0Ztaul22 && h0Ztaunutau12 ) chan = chan0 + 240 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 260 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && !h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau12 && h0Ztaunutau22 ) chan = chan0 + 280 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Ztaul12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 300 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && !h0Ztaul12 && h0Ztaul22 && h0Ztaunutau12 ) chan = chan0 + 320 + tchan;
	     if( h0Zq11 && h0Zq12 ) chan = chan0 + 340 + tchan;
	     if( h0Zl11 && h0Zq12 ) chan = chan0 + 360 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && !h0Ztaul21 && h0Ztaunutau21 && h0Zq12 ) chan = chan0 + 380 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Zq12 ) chan = chan0 + 400 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Zq12 ) chan = chan0 + 420 + tchan;
	     if( h0Ztaul11 && h0Ztaunutau21 && !h0Ztaul21 && h0Zq12 ) chan = chan0 + 440 + tchan;
	     if( h0Zq11 && h0Zl12 ) chan = chan0 + 460 + tchan;
	     if( h0Zq11 && !h0Ztaul12 && h0Ztaunutau12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 480 + tchan;
	     if( h0Zq11 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 500 + tchan;
	     if( h0Zq11 && !h0Ztaul12 && h0Ztaunutau12 && h0Ztaul22 ) chan = chan0 + 520 + tchan;
	     if( h0Zq11 && h0Ztaul12 && h0Ztaunutau22 && !h0Ztaul22 ) chan = chan0 + 540 + tchan;
	     if( h0Zl11 && h0Znu12 ) chan = chan0 + 560 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && !h0Ztaul21 && h0Ztaunutau21 && h0Znu12 ) chan = chan0 + 580 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Znu12 ) chan = chan0 + 600 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Znu12 ) chan = chan0 + 620 + tchan;
	     if( h0Ztaul11 && h0Ztaunutau21 && !h0Ztaul21 && h0Znu12 ) chan = chan0 + 640 + tchan;
	     if( h0Zq11 && h0Znu12 ) chan = chan0 + 660 + tchan;
	     if( h0Znu11 && h0Zq12 ) chan = chan0 + 680 + tchan;
	     if( h0Znu11 && h0Zl12 ) chan = chan0 + 700 + tchan;
	     if( h0Znu11 && !h0Ztaul12 && h0Ztaunutau12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 720 + tchan;
	     if( h0Znu11 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 740 + tchan;
	     if( h0Znu11 && !h0Ztaul12 && h0Ztaunutau12 && h0Ztaul22 ) chan = chan0 + 760 + tchan;
	     if( h0Znu11 && h0Ztaul12 && h0Ztaunutau22 && !h0Ztaul22 ) chan = chan0 + 780 + tchan;
	     if( h0Znu11 && h0Znu12 ) chan = chan0 + 800 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && !h0Ztaul21 && h0Ztaunutau21 && h0Zl12 ) chan = chan0 + 820 + tchan;
	     if( h0Ztaul11 && h0Ztaul21 && h0Zl12 ) chan = chan0 + 840 + tchan;
	     if( !h0Ztaul11 && h0Ztaunutau11 && h0Ztaul21 && h0Zl12 ) chan = chan0 + 860 + tchan;
	     if( h0Ztaul11 && !h0Ztaul21 && h0Ztaunutau21 && h0Zl12 ) chan = chan0 + 880 + tchan;
	     if( h0Zl11 && !h0Ztaul12 && h0Ztaunutau12 && !h0Ztaul22 && h0Ztaunutau22 ) chan = chan0 + 900 + tchan;
	     if( h0Zl11 && h0Ztaul12 && h0Ztaul22 ) chan = chan0 + 920 + tchan;
	     if( h0Zl11 && !h0Ztaul12 && h0Ztaunutau12 && h0Ztaul22 ) chan = chan0 + 940 + tchan;
	     if( h0Zl11 && h0Ztaul12 && h0Ztaunutau22 && !h0Ztaul22 ) chan = chan0 + 960 + tchan;	     
	  }	
	
	if( h0tau1 && h0tau2 )
	  {	     
	     int chan0 = 2000;
	     if( h0taul1 && h0taul2 && h0taunu1 && h0taunu2 ) chan = chan0 + 20 + tchan;
	     if( !h0taul1 && !h0taul2 && h0taunutau1 && h0taunutau2 ) chan = chan0 + 40 + tchan;
	     if( h0taul1 && !h0taul2 && h0taunutau2 && h0taunu1 ) chan = chan0 + 60 + tchan;
	     if( h0taul2 && !h0taul1 && h0taunutau1 && h0taunu2 ) chan = chan0 + 80 + tchan;
	  }	

	if( h0b1 && h0b2 )
	  {	     
	     int chan0 = 3000;
	     if( h0b1 && h0b2 ) chan = chan0 + 20 + tchan;
	  }	
	
	if( chan < 0 && doCheck )
	  {
	     std::cout << "Unknown channel found" << std::endl;
	     std::cout << "chan = " << chan << std::endl;

	     std::cout << "j1 = " << bool(j1) << std::endl;
	     std::cout << "j2 = " << bool(j2) << std::endl;
	     std::cout << "j3 = " << bool(j3) << std::endl;
	     
	     std::cout << "h0 = " << bool(h0) << std::endl;
	     std::cout << "h0->W1 = " << bool(h0W1) << std::endl;
	     std::cout << "h0->W1->l = " << bool(h0Wl1) << std::endl;
	     std::cout << "h0->W1->nu = " << bool(h0Wnu1) << std::endl;
	     std::cout << "h0->W1->tau = " << bool(h0Wtau1) << std::endl;
	     std::cout << "h0->W1->nutau = " << bool(h0Wnutau1) << std::endl;
	     std::cout << "h0->W1->tau->l = " << bool(h0Wtaul1) << std::endl;
	     std::cout << "h0->W1->tau->nu = " << bool(h0Wtaunu1) << std::endl;
	     std::cout << "h0->W1->tau->nutau = " << bool(h0Wtaunutau1) << std::endl;
	     std::cout << "h0->W1->q1 = " << bool(h0Wq11) << std::endl;
	     std::cout << "h0->W1->q2 = " << bool(h0Wq21) << std::endl;
	     std::cout << "h0->W2 = " << bool(h0W2) << std::endl;
	     std::cout << "h0->W2->l = " << bool(h0Wl2) << std::endl;
	     std::cout << "h0->W2->nu = " << bool(h0Wnu2) << std::endl;
	     std::cout << "h0->W2->tau = " << bool(h0Wtau2) << std::endl;
	     std::cout << "h0->W2->nutau = " << bool(h0Wnutau2) << std::endl;
	     std::cout << "h0->W2->tau->l = " << bool(h0Wtaul2) << std::endl;
	     std::cout << "h0->W2->tau->nu = " << bool(h0Wtaunu2) << std::endl;
	     std::cout << "h0->W2->tau->nutau = " << bool(h0Wtaunutau2) << std::endl;
	     std::cout << "h0->W2->q1 = " << bool(h0Wq12) << std::endl;
	     std::cout << "h0->W2->q2 = " << bool(h0Wq22) << std::endl;
	     
	     std::cout << "h0->Z1 = " << bool(h0Z1) << std::endl;
	     std::cout << "h0->Z1->l1 = " << bool(h0Zl11) << std::endl;
	     std::cout << "h0->Z1->l2 = " << bool(h0Zl21) << std::endl;
	     std::cout << "h0->Z1->tau1 = " << bool(h0Ztau11) << std::endl;
	     std::cout << "h0->Z1->tau1->l = " << bool(h0Ztaul11) << std::endl;
	     std::cout << "h0->Z1->tau1->nu = " << bool(h0Ztaunu11) << std::endl;
	     std::cout << "h0->Z1->tau1->nutau = " << bool(h0Ztaunutau11) << std::endl;
	     std::cout << "h0->Z1->tau2 = " << bool(h0Ztau21) << std::endl;
	     std::cout << "h0->Z1->tau2->l = " << bool(h0Ztaul21) << std::endl;
	     std::cout << "h0->Z1->tau2->nu = " << bool(h0Ztaunu21) << std::endl;
	     std::cout << "h0->Z1->tau2->nutau = " << bool(h0Ztaunutau21) << std::endl;
	     std::cout << "h0->Z1->q1 = " << bool(h0Zq11) << std::endl;
	     std::cout << "h0->Z1->q2 = " << bool(h0Zq21) << std::endl;
	     std::cout << "h0->Z1->nu1 = " << bool(h0Znu11) << std::endl;
	     std::cout << "h0->Z1->nu2 = " << bool(h0Znu21) << std::endl;
	     
	     std::cout << "h0->Z2 = " << bool(h0Z2) << std::endl;
	     std::cout << "h0->Z2->l1 = " << bool(h0Zl12) << std::endl;
	     std::cout << "h0->Z2->l2 = " << bool(h0Zl22) << std::endl;
	     std::cout << "h0->Z2->tau1 = " << bool(h0Ztau12) << std::endl;
	     std::cout << "h0->Z2->tau1->l = " << bool(h0Ztaul12) << std::endl;
	     std::cout << "h0->Z2->tau1->nu = " << bool(h0Ztaunu12) << std::endl;
	     std::cout << "h0->Z2->tau1->nutau = " << bool(h0Ztaunutau12) << std::endl;
	     std::cout << "h0->Z2->tau2 = " << bool(h0Ztau22) << std::endl;
	     std::cout << "h0->Z2->tau2->l = " << bool(h0Ztaul22) << std::endl;
	     std::cout << "h0->Z2->tau2->nu = " << bool(h0Ztaunu22) << std::endl;
	     std::cout << "h0->Z2->tau2->nutau = " << bool(h0Ztaunutau22) << std::endl;
	     std::cout << "h0->Z2->q1 = " << bool(h0Zq12) << std::endl;
	     std::cout << "h0->Z2->q2 = " << bool(h0Zq22) << std::endl;
	     std::cout << "h0->Z2->nu1 = " << bool(h0Znu12) << std::endl;
	     std::cout << "h0->Z2->nu2 = " << bool(h0Znu22) << std::endl;
	     
	     std::cout << "h0->tau1 = " << bool(h0tau1) << std::endl;
	     std::cout << "h0->tau1->l = " << bool(h0taul1) << std::endl;
	     std::cout << "h0->tau1->nu_tau = " << bool(h0taunutau1) << std::endl;
	     std::cout << "h0->tau1->nu = " << bool(h0taunu1) << std::endl;
	             
	     std::cout << "h0->tau2 = " << bool(h0tau2) << std::endl;
	     std::cout << "h0->tau2->l = " << bool(h0taul2) << std::endl;
	     std::cout << "h0->tau2->nu_tau = " << bool(h0taunutau2) << std::endl;
	     std::cout << "h0->tau2->nu = " << bool(h0taunu2) << std::endl;

	     std::cout << "h0->b1 = " << bool(h0b1) << std::endl;
	     std::cout << "h0->b2 = " << bool(h0b2) << std::endl;
	     
	     exit(1);
	  }
     }

   tree.mc_truth_thq_channel = chan;

   // TLV

   if( h0 ) p4toTLV(h0->p4(),tree.mc_truth_h0_p4);      
   
   if( h0W1 ) p4toTLV(h0W1->p4(),tree.mc_truth_h0W1_p4);
   if( h0W2 ) p4toTLV(h0W2->p4(),tree.mc_truth_h0W2_p4);
   if( h0Wl1 ) p4toTLV(h0Wl1->p4(),tree.mc_truth_h0Wl1_p4);
   if( h0Wnu1 ) p4toTLV(h0Wnu1->p4(),tree.mc_truth_h0Wnu1_p4);
   if( h0Wtau1 ) p4toTLV(h0Wtau1->p4(),tree.mc_truth_h0Wtau1_p4);
   if( h0Wnutau1 ) p4toTLV(h0Wnutau1->p4(),tree.mc_truth_h0Wnutau1_p4);
   if( h0Wtaul1 ) p4toTLV(h0Wtaul1->p4(),tree.mc_truth_h0Wtaul1_p4);
   if( h0Wtaunu1 ) p4toTLV(h0Wtaunu1->p4(),tree.mc_truth_h0Wtaunu1_p4);
   if( h0Wtaunutau1 ) p4toTLV(h0Wtaunutau1->p4(),tree.mc_truth_h0Wtaunutau1_p4);
   if( h0Wl2 ) p4toTLV(h0Wl2->p4(),tree.mc_truth_h0Wl2_p4);
   if( h0Wnu2 ) p4toTLV(h0Wnu2->p4(),tree.mc_truth_h0Wnu2_p4);
   if( h0Wtau2 ) p4toTLV(h0Wtau2->p4(),tree.mc_truth_h0Wtau2_p4);
   if( h0Wnutau2 ) p4toTLV(h0Wnutau2->p4(),tree.mc_truth_h0Wnutau2_p4);
   if( h0Wtaul2 ) p4toTLV(h0Wtaul2->p4(),tree.mc_truth_h0Wtaul2_p4);
   if( h0Wtaunu2 ) p4toTLV(h0Wtaunu2->p4(),tree.mc_truth_h0Wtaunu2_p4);
   if( h0Wtaunutau2 ) p4toTLV(h0Wtaunutau2->p4(),tree.mc_truth_h0Wtaunutau2_p4);
   if( h0Wq11 ) p4toTLV(h0Wq11->p4(),tree.mc_truth_h0Wq11_p4);
   if( h0Wq21 ) p4toTLV(h0Wq21->p4(),tree.mc_truth_h0Wq21_p4);
   if( h0Wq12 ) p4toTLV(h0Wq12->p4(),tree.mc_truth_h0Wq12_p4);
   if( h0Wq22 ) p4toTLV(h0Wq22->p4(),tree.mc_truth_h0Wq22_p4);
   if( h0Wq11_IS ) p4toTLV(h0Wq11_IS->p4(),tree.mc_truth_h0Wq11_IS_p4);
   if( h0Wq21_IS ) p4toTLV(h0Wq21_IS->p4(),tree.mc_truth_h0Wq21_IS_p4);
   if( h0Wq12_IS ) p4toTLV(h0Wq12_IS->p4(),tree.mc_truth_h0Wq12_IS_p4);
   if( h0Wq22_IS ) p4toTLV(h0Wq22_IS->p4(),tree.mc_truth_h0Wq22_IS_p4);
   
   if( h0Z1 ) p4toTLV(h0Z1->p4(),tree.mc_truth_h0Z1_p4);
   if( h0Z2 ) p4toTLV(h0Z2->p4(),tree.mc_truth_h0Z2_p4);
   if( h0Zl11 ) p4toTLV(h0Zl11->p4(),tree.mc_truth_h0Zl11_p4);
   if( h0Zl21 ) p4toTLV(h0Zl21->p4(),tree.mc_truth_h0Zl21_p4);
   if( h0Ztau11 ) p4toTLV(h0Ztau11->p4(),tree.mc_truth_h0Ztau11_p4);
   if( h0Ztau21 ) p4toTLV(h0Ztau21->p4(),tree.mc_truth_h0Ztau21_p4);
   if( h0Ztaul11 ) p4toTLV(h0Ztaul11->p4(),tree.mc_truth_h0Ztaul11_p4);
   if( h0Ztaul21 ) p4toTLV(h0Ztaul21->p4(),tree.mc_truth_h0Ztaul21_p4);
   if( h0Ztaunu11 ) p4toTLV(h0Ztaunu11->p4(),tree.mc_truth_h0Ztaunu11_p4);
   if( h0Ztaunu21 ) p4toTLV(h0Ztaunu21->p4(),tree.mc_truth_h0Ztaunu21_p4);
   if( h0Ztaunutau11 ) p4toTLV(h0Ztaunutau11->p4(),tree.mc_truth_h0Ztaunutau11_p4);
   if( h0Ztaunutau21 ) p4toTLV(h0Ztaunutau21->p4(),tree.mc_truth_h0Ztaunutau21_p4);
   if( h0Zq11 ) p4toTLV(h0Zq11->p4(),tree.mc_truth_h0Zq11_p4);
   if( h0Zq21 ) p4toTLV(h0Zq21->p4(),tree.mc_truth_h0Zq21_p4);
   if( h0Zq11_IS ) p4toTLV(h0Zq11_IS->p4(),tree.mc_truth_h0Zq11_IS_p4);
   if( h0Zq21_IS ) p4toTLV(h0Zq21_IS->p4(),tree.mc_truth_h0Zq21_IS_p4);
   if( h0Zl12 ) p4toTLV(h0Zl12->p4(),tree.mc_truth_h0Zl12_p4);
   if( h0Zl22 ) p4toTLV(h0Zl22->p4(),tree.mc_truth_h0Zl22_p4);
   if( h0Ztau12 ) p4toTLV(h0Ztau12->p4(),tree.mc_truth_h0Ztau12_p4);
   if( h0Ztau22 ) p4toTLV(h0Ztau22->p4(),tree.mc_truth_h0Ztau22_p4);
   if( h0Ztaul12 ) p4toTLV(h0Ztaul12->p4(),tree.mc_truth_h0Ztaul12_p4);
   if( h0Ztaul22 ) p4toTLV(h0Ztaul22->p4(),tree.mc_truth_h0Ztaul22_p4);
   if( h0Ztaunu12 ) p4toTLV(h0Ztaunu12->p4(),tree.mc_truth_h0Ztaunu12_p4);
   if( h0Ztaunu22 ) p4toTLV(h0Ztaunu22->p4(),tree.mc_truth_h0Ztaunu22_p4);
   if( h0Ztaunutau12 ) p4toTLV(h0Ztaunutau12->p4(),tree.mc_truth_h0Ztaunutau12_p4);
   if( h0Ztaunutau22 ) p4toTLV(h0Ztaunutau22->p4(),tree.mc_truth_h0Ztaunutau22_p4);
   if( h0Zq12 ) p4toTLV(h0Zq12->p4(),tree.mc_truth_h0Zq12_p4);
   if( h0Zq22 ) p4toTLV(h0Zq22->p4(),tree.mc_truth_h0Zq22_p4);
   if( h0Zq12_IS ) p4toTLV(h0Zq12_IS->p4(),tree.mc_truth_h0Zq12_IS_p4);
   if( h0Zq22_IS ) p4toTLV(h0Zq22_IS->p4(),tree.mc_truth_h0Zq22_IS_p4);
   if( h0Znu11 ) p4toTLV(h0Znu11->p4(),tree.mc_truth_h0Znu11_p4);
   if( h0Znu21 ) p4toTLV(h0Znu21->p4(),tree.mc_truth_h0Znu21_p4);
   if( h0Znu12 ) p4toTLV(h0Znu12->p4(),tree.mc_truth_h0Znu12_p4);
   if( h0Znu22 ) p4toTLV(h0Znu22->p4(),tree.mc_truth_h0Znu22_p4);
   
   if( h0tau1 ) p4toTLV(h0tau1->p4(),tree.mc_truth_h0tau1_p4);
   if( h0tau2 ) p4toTLV(h0tau2->p4(),tree.mc_truth_h0tau2_p4);
   if( h0taul1 ) p4toTLV(h0taul1->p4(),tree.mc_truth_h0taul1_p4);
   if( h0taunutau1 ) p4toTLV(h0taunutau1->p4(),tree.mc_truth_h0taunutau1_p4);
   if( h0taunu1 ) p4toTLV(h0taunu1->p4(),tree.mc_truth_h0taunu1_p4);
   if( h0taul2 ) p4toTLV(h0taul2->p4(),tree.mc_truth_h0taul2_p4);
   if( h0taunutau2 ) p4toTLV(h0taunutau2->p4(),tree.mc_truth_h0taunutau2_p4);
   if( h0taunu2 ) p4toTLV(h0taunu2->p4(),tree.mc_truth_h0taunu2_p4);

   if( h0b1 ) p4toTLV(h0b1->p4(),tree.mc_truth_h0b1_p4);
   if( h0b2 ) p4toTLV(h0b2->p4(),tree.mc_truth_h0b2_p4);
   if( h0b1_IS ) p4toTLV(h0b1_IS->p4(),tree.mc_truth_h0b1_IS_p4);
   if( h0b2_IS ) p4toTLV(h0b2_IS->p4(),tree.mc_truth_h0b2_IS_p4);
   
   if( t ) p4toTLV(t->p4(),tree.mc_truth_t_p4);
   if( tb ) p4toTLV(tb->p4(),tree.mc_truth_tb_p4);
   if( tb_IS ) p4toTLV(tb_IS->p4(),tree.mc_truth_tb_IS_p4);
   
   if( tW ) p4toTLV(tW->p4(),tree.mc_truth_tW_p4);
   if( tWnu ) p4toTLV(tWnu->p4(),tree.mc_truth_tWnu_p4);
   if( tWnutau ) p4toTLV(tWnutau->p4(),tree.mc_truth_tWnutau_p4);
   if( tWl ) p4toTLV(tWl->p4(),tree.mc_truth_tWl_p4);
   if( tWtau ) p4toTLV(tWtau->p4(),tree.mc_truth_tWtau_p4);
   if( tWtaunu ) p4toTLV(tWtaunu->p4(),tree.mc_truth_tWtaunu_p4);
   if( tWtaunutau ) p4toTLV(tWtaunutau->p4(),tree.mc_truth_tWtaunutau_p4);
   if( tWtaul ) p4toTLV(tWtaul->p4(),tree.mc_truth_tWtaul_p4);
   if( tWq1 ) p4toTLV(tWq1->p4(),tree.mc_truth_tWq1_p4);
   if( tWq2 ) p4toTLV(tWq2->p4(),tree.mc_truth_tWq2_p4);
   if( tWq1_IS ) p4toTLV(tWq1_IS->p4(),tree.mc_truth_tWq1_IS_p4);
   if( tWq2_IS ) p4toTLV(tWq2_IS->p4(),tree.mc_truth_tWq2_IS_p4);

   if( j1 ) p4toTLV(j1->p4(),tree.mc_truth_j1_p4);
   if( j2 ) p4toTLV(j2->p4(),tree.mc_truth_j2_p4);
   if( j3 ) p4toTLV(j3->p4(),tree.mc_truth_j3_p4);

   // pt

   if( h0 ) tree.mc_truth_h0_pt = h0->p4().pt();

   if( h0W1 ) tree.mc_truth_h0W1_pt = h0W1->p4().pt();
   if( h0W2 ) tree.mc_truth_h0W2_pt = h0W2->p4().pt();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_pt = h0Wl1->p4().pt();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_pt = h0Wnu1->p4().pt();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_pt = h0Wtau1->p4().pt();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_pt = h0Wnutau1->p4().pt();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_pt = h0Wtaul1->p4().pt();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_pt = h0Wtaunu1->p4().pt();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_pt = h0Wtaunutau1->p4().pt();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_pt = h0Wl2->p4().pt();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_pt = h0Wnu2->p4().pt();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_pt = h0Wtau2->p4().pt();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_pt = h0Wnutau2->p4().pt();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_pt = h0Wtaul2->p4().pt();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_pt = h0Wtaunu2->p4().pt();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_pt = h0Wtaunutau2->p4().pt();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_pt = h0Wq11->p4().pt();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_pt = h0Wq21->p4().pt();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_pt = h0Wq12->p4().pt();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_pt = h0Wq22->p4().pt();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_pt = h0Wq11_IS->p4().pt();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_pt = h0Wq21_IS->p4().pt();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_pt = h0Wq12_IS->p4().pt();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_pt = h0Wq22_IS->p4().pt();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_pt = h0Z1->p4().pt();
   if( h0Z2 ) tree.mc_truth_h0Z2_pt = h0Z2->p4().pt();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_pt = h0Zl11->p4().pt();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_pt = h0Zl21->p4().pt();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_pt = h0Zl12->p4().pt();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_pt = h0Zl22->p4().pt();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_pt = h0Ztau11->p4().pt();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_pt = h0Ztau21->p4().pt();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_pt = h0Ztaul11->p4().pt();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_pt = h0Ztaul21->p4().pt();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_pt = h0Ztaunu11->p4().pt();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_pt = h0Ztaunu21->p4().pt();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_pt = h0Ztaunutau11->p4().pt();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_pt = h0Ztaunutau21->p4().pt();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_pt = h0Zq11->p4().pt();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_pt = h0Zq21->p4().pt();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_pt = h0Zq12->p4().pt();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_pt = h0Zq22->p4().pt();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_pt = h0Zq11_IS->p4().pt();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_pt = h0Zq21_IS->p4().pt();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_pt = h0Zq12_IS->p4().pt();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_pt = h0Zq22_IS->p4().pt();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_pt = h0Ztau12->p4().pt();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_pt = h0Ztau22->p4().pt();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_pt = h0Ztaul12->p4().pt();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_pt = h0Ztaul22->p4().pt();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_pt = h0Ztaunu12->p4().pt();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_pt = h0Ztaunu22->p4().pt();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_pt = h0Ztaunutau12->p4().pt();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_pt = h0Ztaunutau22->p4().pt();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_pt = h0Znu11->p4().pt();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_pt = h0Znu21->p4().pt();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_pt = h0Znu12->p4().pt();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_pt = h0Znu22->p4().pt();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_pt = h0tau1->p4().pt();
   if( h0tau2 ) tree.mc_truth_h0tau2_pt = h0tau2->p4().pt();
   if( h0taul1 ) tree.mc_truth_h0taul1_pt = h0taul1->p4().pt();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_pt = h0taunutau1->p4().pt();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_pt = h0taunu1->p4().pt();
   if( h0taul2 ) tree.mc_truth_h0taul2_pt = h0taul2->p4().pt();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_pt = h0taunutau2->p4().pt();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_pt = h0taunu2->p4().pt();

   if( h0b1 ) tree.mc_truth_h0b1_pt = h0b1->p4().pt();
   if( h0b2 ) tree.mc_truth_h0b2_pt = h0b2->p4().pt();
   if( h0b1_IS ) tree.mc_truth_h0b1_IS_pt = h0b1_IS->p4().pt();
   if( h0b2_IS ) tree.mc_truth_h0b2_IS_pt = h0b2_IS->p4().pt();
   
   if( t ) tree.mc_truth_t_pt = t->p4().pt();
   if( tb ) tree.mc_truth_tb_pt = tb->p4().pt();
   if( tb_IS ) tree.mc_truth_tb_IS_pt = tb_IS->p4().pt();
   
   if( tW ) tree.mc_truth_tW_pt = tW->p4().pt();
   if( tWnu ) tree.mc_truth_tWnu_pt = tWnu->p4().pt();
   if( tWnutau ) tree.mc_truth_tWnutau_pt = tWnutau->p4().pt();
   if( tWl ) tree.mc_truth_tWl_pt = tWl->p4().pt();
   if( tWtau ) tree.mc_truth_tWtau_pt = tWtau->p4().pt();
   if( tWtaunu ) tree.mc_truth_tWtaunu_pt = tWtaunu->p4().pt();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_pt = tWtaunutau->p4().pt();
   if( tWtaul ) tree.mc_truth_tWtaul_pt = tWtaul->p4().pt();
   if( tWq1 ) tree.mc_truth_tWq1_pt = tWq1->p4().pt();
   if( tWq2 ) tree.mc_truth_tWq2_pt = tWq2->p4().pt();
   if( tWq1_IS ) tree.mc_truth_tWq1_IS_pt = tWq1_IS->p4().pt();
   if( tWq2_IS ) tree.mc_truth_tWq2_IS_pt = tWq2_IS->p4().pt();
   
   if( j1 ) tree.mc_truth_j1_pt = j1->p4().pt();
   if( j2 ) tree.mc_truth_j2_pt = j2->p4().pt();
   if( j3 ) tree.mc_truth_j3_pt = j3->p4().pt();
   
   // eta

   if( h0 ) tree.mc_truth_h0_eta = h0->p4().eta();

   if( h0W1 ) tree.mc_truth_h0W1_eta = h0W1->p4().eta();
   if( h0W2 ) tree.mc_truth_h0W2_eta = h0W2->p4().eta();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_eta = h0Wl1->p4().eta();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_eta = h0Wnu1->p4().eta();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_eta = h0Wtau1->p4().eta();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_eta = h0Wnutau1->p4().eta();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_eta = h0Wtaul1->p4().eta();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_eta = h0Wtaunu1->p4().eta();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_eta = h0Wtaunutau1->p4().eta();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_eta = h0Wl2->p4().eta();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_eta = h0Wnu2->p4().eta();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_eta = h0Wtau2->p4().eta();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_eta = h0Wnutau2->p4().eta();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_eta = h0Wtaul2->p4().eta();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_eta = h0Wtaunu2->p4().eta();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_eta = h0Wtaunutau2->p4().eta();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_eta = h0Wq11->p4().eta();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_eta = h0Wq21->p4().eta();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_eta = h0Wq12->p4().eta();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_eta = h0Wq22->p4().eta();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_eta = h0Wq11_IS->p4().eta();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_eta = h0Wq21_IS->p4().eta();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_eta = h0Wq12_IS->p4().eta();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_eta = h0Wq22_IS->p4().eta();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_eta = h0Z1->p4().eta();
   if( h0Z2 ) tree.mc_truth_h0Z2_eta = h0Z2->p4().eta();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_eta = h0Zl11->p4().eta();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_eta = h0Zl21->p4().eta();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_eta = h0Zl12->p4().eta();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_eta = h0Zl22->p4().eta();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_eta = h0Ztau11->p4().eta();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_eta = h0Ztau21->p4().eta();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_eta = h0Ztaul11->p4().eta();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_eta = h0Ztaul21->p4().eta();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_eta = h0Ztaunu11->p4().eta();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_eta = h0Ztaunu21->p4().eta();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_eta = h0Ztaunutau11->p4().eta();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_eta = h0Ztaunutau21->p4().eta();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_eta = h0Zq11->p4().eta();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_eta = h0Zq21->p4().eta();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_eta = h0Zq12->p4().eta();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_eta = h0Zq22->p4().eta();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_eta = h0Zq11_IS->p4().eta();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_eta = h0Zq21_IS->p4().eta();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_eta = h0Zq12_IS->p4().eta();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_eta = h0Zq22_IS->p4().eta();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_eta = h0Ztau12->p4().eta();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_eta = h0Ztau22->p4().eta();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_eta = h0Ztaul12->p4().eta();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_eta = h0Ztaul22->p4().eta();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_eta = h0Ztaunu12->p4().eta();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_eta = h0Ztaunu22->p4().eta();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_eta = h0Ztaunutau12->p4().eta();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_eta = h0Ztaunutau22->p4().eta();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_eta = h0Znu11->p4().eta();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_eta = h0Znu21->p4().eta();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_eta = h0Znu12->p4().eta();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_eta = h0Znu22->p4().eta();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_eta = h0tau1->p4().eta();
   if( h0tau2 ) tree.mc_truth_h0tau2_eta = h0tau2->p4().eta();
   if( h0taul1 ) tree.mc_truth_h0taul1_eta = h0taul1->p4().eta();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_eta = h0taunutau1->p4().eta();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_eta = h0taunu1->p4().eta();
   if( h0taul2 ) tree.mc_truth_h0taul2_eta = h0taul2->p4().eta();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_eta = h0taunutau2->p4().eta();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_eta = h0taunu2->p4().eta();

   if( h0b1 ) tree.mc_truth_h0b1_eta = h0b1->p4().eta();
   if( h0b2 ) tree.mc_truth_h0b2_eta = h0b2->p4().eta();
   if( h0b1_IS ) tree.mc_truth_h0b1_IS_eta = h0b1_IS->p4().eta();
   if( h0b2_IS ) tree.mc_truth_h0b2_IS_eta = h0b2_IS->p4().eta();
   
   if( t ) tree.mc_truth_t_eta = t->p4().eta();
   if( tb ) tree.mc_truth_tb_eta = tb->p4().eta();
   if( tb_IS ) tree.mc_truth_tb_IS_eta = tb_IS->p4().eta();
   
   if( tW ) tree.mc_truth_tW_eta = tW->p4().eta();
   if( tWnu ) tree.mc_truth_tWnu_eta = tWnu->p4().eta();
   if( tWnutau ) tree.mc_truth_tWnutau_eta = tWnutau->p4().eta();
   if( tWl ) tree.mc_truth_tWl_eta = tWl->p4().eta();
   if( tWtau ) tree.mc_truth_tWtau_eta = tWtau->p4().eta();
   if( tWtaunu ) tree.mc_truth_tWtaunu_eta = tWtaunu->p4().eta();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_eta = tWtaunutau->p4().eta();
   if( tWtaul ) tree.mc_truth_tWtaul_eta = tWtaul->p4().eta();
   if( tWq1 ) tree.mc_truth_tWq1_eta = tWq1->p4().eta();
   if( tWq2 ) tree.mc_truth_tWq2_eta = tWq2->p4().eta();
   if( tWq1_IS ) tree.mc_truth_tWq1_IS_eta = tWq1_IS->p4().eta();
   if( tWq2_IS ) tree.mc_truth_tWq2_IS_eta = tWq2_IS->p4().eta();
   
   if( j1 ) tree.mc_truth_j1_eta = j1->p4().eta();
   if( j2 ) tree.mc_truth_j2_eta = j2->p4().eta();
   if( j3 ) tree.mc_truth_j3_eta = j3->p4().eta();

   // phi

   if( h0 ) tree.mc_truth_h0_phi = h0->p4().phi();

   if( h0W1 ) tree.mc_truth_h0W1_phi = h0W1->p4().phi();
   if( h0W2 ) tree.mc_truth_h0W2_phi = h0W2->p4().phi();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_phi = h0Wl1->p4().phi();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_phi = h0Wnu1->p4().phi();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_phi = h0Wtau1->p4().phi();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_phi = h0Wnutau1->p4().phi();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_phi = h0Wtaul1->p4().phi();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_phi = h0Wtaunu1->p4().phi();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_phi = h0Wtaunutau1->p4().phi();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_phi = h0Wl2->p4().phi();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_phi = h0Wnu2->p4().phi();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_phi = h0Wtau2->p4().phi();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_phi = h0Wnutau2->p4().phi();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_phi = h0Wtaul2->p4().phi();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_phi = h0Wtaunu2->p4().phi();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_phi = h0Wtaunutau2->p4().phi();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_phi = h0Wq11->p4().phi();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_phi = h0Wq21->p4().phi();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_phi = h0Wq12->p4().phi();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_phi = h0Wq22->p4().phi();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_phi = h0Wq11_IS->p4().phi();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_phi = h0Wq21_IS->p4().phi();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_phi = h0Wq12_IS->p4().phi();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_phi = h0Wq22_IS->p4().phi();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_phi = h0Z1->p4().phi();
   if( h0Z2 ) tree.mc_truth_h0Z2_phi = h0Z2->p4().phi();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_phi = h0Zl11->p4().phi();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_phi = h0Zl21->p4().phi();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_phi = h0Zl12->p4().phi();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_phi = h0Zl22->p4().phi();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_phi = h0Ztau11->p4().phi();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_phi = h0Ztau21->p4().phi();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_phi = h0Ztaul11->p4().phi();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_phi = h0Ztaul21->p4().phi();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_phi = h0Ztaunu11->p4().phi();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_phi = h0Ztaunu21->p4().phi();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_phi = h0Ztaunutau11->p4().phi();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_phi = h0Ztaunutau21->p4().phi();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_phi = h0Zq11->p4().phi();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_phi = h0Zq21->p4().phi();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_phi = h0Zq12->p4().phi();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_phi = h0Zq22->p4().phi();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_phi = h0Zq11_IS->p4().phi();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_phi = h0Zq21_IS->p4().phi();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_phi = h0Zq12_IS->p4().phi();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_phi = h0Zq22_IS->p4().phi();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_phi = h0Ztau12->p4().phi();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_phi = h0Ztau22->p4().phi();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_phi = h0Ztaul12->p4().phi();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_phi = h0Ztaul22->p4().phi();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_phi = h0Ztaunu12->p4().phi();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_phi = h0Ztaunu22->p4().phi();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_phi = h0Ztaunutau12->p4().phi();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_phi = h0Ztaunutau22->p4().phi();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_phi = h0Znu11->p4().phi();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_phi = h0Znu21->p4().phi();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_phi = h0Znu12->p4().phi();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_phi = h0Znu22->p4().phi();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_phi = h0tau1->p4().phi();
   if( h0tau2 ) tree.mc_truth_h0tau2_phi = h0tau2->p4().phi();
   if( h0taul1 ) tree.mc_truth_h0taul1_phi = h0taul1->p4().phi();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_phi = h0taunutau1->p4().phi();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_phi = h0taunu1->p4().phi();
   if( h0taul2 ) tree.mc_truth_h0taul2_phi = h0taul2->p4().phi();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_phi = h0taunutau2->p4().phi();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_phi = h0taunu2->p4().phi();

   if( h0b1 ) tree.mc_truth_h0b1_phi = h0b1->p4().phi();
   if( h0b2 ) tree.mc_truth_h0b2_phi = h0b2->p4().phi();
   if( h0b1_IS ) tree.mc_truth_h0b1_IS_phi = h0b1_IS->p4().phi();
   if( h0b2_IS ) tree.mc_truth_h0b2_IS_phi = h0b2_IS->p4().phi();
   
   if( t ) tree.mc_truth_t_phi = t->p4().phi();
   if( tb ) tree.mc_truth_tb_phi = tb->p4().phi();
   if( tb_IS ) tree.mc_truth_tb_IS_phi = tb_IS->p4().phi();
   
   if( tW ) tree.mc_truth_tW_phi = tW->p4().phi();
   if( tWnu ) tree.mc_truth_tWnu_phi = tWnu->p4().phi();
   if( tWnutau ) tree.mc_truth_tWnutau_phi = tWnutau->p4().phi();
   if( tWl ) tree.mc_truth_tWl_phi = tWl->p4().phi();
   if( tWtau ) tree.mc_truth_tWtau_phi = tWtau->p4().phi();
   if( tWtaunu ) tree.mc_truth_tWtaunu_phi = tWtaunu->p4().phi();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_phi = tWtaunutau->p4().phi();
   if( tWtaul ) tree.mc_truth_tWtaul_phi = tWtaul->p4().phi();
   if( tWq1 ) tree.mc_truth_tWq1_phi = tWq1->p4().phi();
   if( tWq2 ) tree.mc_truth_tWq2_phi = tWq2->p4().phi();
   if( tWq1_IS ) tree.mc_truth_tWq1_IS_phi = tWq1_IS->p4().phi();
   if( tWq2_IS ) tree.mc_truth_tWq2_IS_phi = tWq2_IS->p4().phi();
   
   if( j1 ) tree.mc_truth_j1_phi = j1->p4().phi();
   if( j2 ) tree.mc_truth_j2_phi = j2->p4().phi();
   if( j3 ) tree.mc_truth_j3_phi = j3->p4().phi();

   // E

   if( h0 ) tree.mc_truth_h0_E = h0->p4().E();

   if( h0W1 ) tree.mc_truth_h0W1_E = h0W1->p4().E();
   if( h0W2 ) tree.mc_truth_h0W2_E = h0W2->p4().E();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_E = h0Wl1->p4().E();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_E = h0Wnu1->p4().E();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_E = h0Wtau1->p4().E();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_E = h0Wnutau1->p4().E();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_E = h0Wtaul1->p4().E();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_E = h0Wtaunu1->p4().E();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_E = h0Wtaunutau1->p4().E();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_E = h0Wl2->p4().E();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_E = h0Wnu2->p4().E();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_E = h0Wtau2->p4().E();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_E = h0Wnutau2->p4().E();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_E = h0Wtaul2->p4().E();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_E = h0Wtaunu2->p4().E();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_E = h0Wtaunutau2->p4().E();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_E = h0Wq11->p4().E();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_E = h0Wq21->p4().E();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_E = h0Wq12->p4().E();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_E = h0Wq22->p4().E();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_E = h0Wq11_IS->p4().E();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_E = h0Wq21_IS->p4().E();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_E = h0Wq12_IS->p4().E();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_E = h0Wq22_IS->p4().E();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_E = h0Z1->p4().E();
   if( h0Z2 ) tree.mc_truth_h0Z2_E = h0Z2->p4().E();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_E = h0Zl11->p4().E();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_E = h0Zl21->p4().E();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_E = h0Zl12->p4().E();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_E = h0Zl22->p4().E();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_E = h0Ztau11->p4().E();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_E = h0Ztau21->p4().E();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_E = h0Ztaul11->p4().E();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_E = h0Ztaul21->p4().E();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_E = h0Ztaunu11->p4().E();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_E = h0Ztaunu21->p4().E();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_E = h0Ztaunutau11->p4().E();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_E = h0Ztaunutau21->p4().E();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_E = h0Zq11->p4().E();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_E = h0Zq21->p4().E();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_E = h0Zq12->p4().E();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_E = h0Zq22->p4().E();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_E = h0Zq11_IS->p4().E();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_E = h0Zq21_IS->p4().E();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_E = h0Zq12_IS->p4().E();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_E = h0Zq22_IS->p4().E();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_E = h0Ztau12->p4().E();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_E = h0Ztau22->p4().E();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_E = h0Ztaul12->p4().E();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_E = h0Ztaul22->p4().E();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_E = h0Ztaunu12->p4().E();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_E = h0Ztaunu22->p4().E();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_E = h0Ztaunutau12->p4().E();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_E = h0Ztaunutau22->p4().E();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_E = h0Znu11->p4().E();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_E = h0Znu21->p4().E();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_E = h0Znu12->p4().E();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_E = h0Znu22->p4().E();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_E = h0tau1->p4().E();
   if( h0tau2 ) tree.mc_truth_h0tau2_E = h0tau2->p4().E();
   if( h0taul1 ) tree.mc_truth_h0taul1_E = h0taul1->p4().E();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_E = h0taunutau1->p4().E();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_E = h0taunu1->p4().E();
   if( h0taul2 ) tree.mc_truth_h0taul2_E = h0taul2->p4().E();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_E = h0taunutau2->p4().E();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_E = h0taunu2->p4().E();

   if( h0b1 ) tree.mc_truth_h0b1_E = h0b1->p4().E();
   if( h0b2 ) tree.mc_truth_h0b2_E = h0b2->p4().E();
   if( h0b1_IS ) tree.mc_truth_h0b1_IS_E = h0b1_IS->p4().E();
   if( h0b2_IS ) tree.mc_truth_h0b2_IS_E = h0b2_IS->p4().E();
   
   if( t ) tree.mc_truth_t_E = t->p4().E();
   if( tb ) tree.mc_truth_tb_E = tb->p4().E();
   if( tb_IS ) tree.mc_truth_tb_IS_E = tb_IS->p4().E();
   
   if( tW ) tree.mc_truth_tW_E = tW->p4().E();
   if( tWnu ) tree.mc_truth_tWnu_E = tWnu->p4().E();
   if( tWnutau ) tree.mc_truth_tWnutau_E = tWnutau->p4().E();
   if( tWl ) tree.mc_truth_tWl_E = tWl->p4().E();
   if( tWtau ) tree.mc_truth_tWtau_E = tWtau->p4().E();
   if( tWtaunu ) tree.mc_truth_tWtaunu_E = tWtaunu->p4().E();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_E = tWtaunutau->p4().E();
   if( tWtaul ) tree.mc_truth_tWtaul_E = tWtaul->p4().E();
   if( tWq1 ) tree.mc_truth_tWq1_E = tWq1->p4().E();
   if( tWq2 ) tree.mc_truth_tWq2_E = tWq2->p4().E();
   if( tWq1_IS ) tree.mc_truth_tWq1_IS_E = tWq1_IS->p4().E();
   if( tWq2_IS ) tree.mc_truth_tWq2_IS_E = tWq2_IS->p4().E();
   
   if( j1 ) tree.mc_truth_j1_E = j1->p4().E();
   if( j2 ) tree.mc_truth_j2_E = j2->p4().E();
   if( j3 ) tree.mc_truth_j3_E = j3->p4().E();
   
   // pdgId

   if( h0 ) tree.mc_truth_h0_id = h0->pdgId();

   if( h0W1 ) tree.mc_truth_h0W1_id = h0W1->pdgId();
   if( h0W2 ) tree.mc_truth_h0W2_id = h0W2->pdgId();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_id = h0Wl1->pdgId();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_id = h0Wnu1->pdgId();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_id = h0Wtau1->pdgId();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_id = h0Wnutau1->pdgId();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_id = h0Wtaul1->pdgId();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_id = h0Wtaunu1->pdgId();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_id = h0Wtaunutau1->pdgId();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_id = h0Wl2->pdgId();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_id = h0Wnu2->pdgId();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_id = h0Wtau2->pdgId();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_id = h0Wnutau2->pdgId();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_id = h0Wtaul2->pdgId();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_id = h0Wtaunu2->pdgId();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_id = h0Wtaunutau2->pdgId();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_id = h0Wq11->pdgId();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_id = h0Wq21->pdgId();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_id = h0Wq12->pdgId();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_id = h0Wq22->pdgId();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_id = h0Wq11_IS->pdgId();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_id = h0Wq21_IS->pdgId();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_id = h0Wq12_IS->pdgId();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_id = h0Wq22_IS->pdgId();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_id = h0Z1->pdgId();
   if( h0Z2 ) tree.mc_truth_h0Z2_id = h0Z2->pdgId();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_id = h0Zl11->pdgId();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_id = h0Zl21->pdgId();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_id = h0Zl12->pdgId();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_id = h0Zl22->pdgId();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_id = h0Ztau11->pdgId();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_id = h0Ztau21->pdgId();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_id = h0Ztaul11->pdgId();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_id = h0Ztaul21->pdgId();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_id = h0Ztaunu11->pdgId();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_id = h0Ztaunu21->pdgId();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_id = h0Ztaunutau11->pdgId();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_id = h0Ztaunutau21->pdgId();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_id = h0Zq11->pdgId();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_id = h0Zq21->pdgId();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_id = h0Zq12->pdgId();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_id = h0Zq22->pdgId();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_id = h0Zq11_IS->pdgId();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_id = h0Zq21_IS->pdgId();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_id = h0Zq12_IS->pdgId();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_id = h0Zq22_IS->pdgId();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_id = h0Ztau12->pdgId();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_id = h0Ztau22->pdgId();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_id = h0Ztaul12->pdgId();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_id = h0Ztaul22->pdgId();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_id = h0Ztaunu12->pdgId();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_id = h0Ztaunu22->pdgId();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_id = h0Ztaunutau12->pdgId();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_id = h0Ztaunutau22->pdgId();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_id = h0Znu11->pdgId();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_id = h0Znu21->pdgId();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_id = h0Znu12->pdgId();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_id = h0Znu22->pdgId();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_id = h0tau1->pdgId();
   if( h0tau2 ) tree.mc_truth_h0tau2_id = h0tau2->pdgId();
   if( h0taul1 ) tree.mc_truth_h0taul1_id = h0taul1->pdgId();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_id = h0taunutau1->pdgId();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_id = h0taunu1->pdgId();
   if( h0taul2 ) tree.mc_truth_h0taul2_id = h0taul2->pdgId();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_id = h0taunutau2->pdgId();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_id = h0taunu2->pdgId();

   if( h0b1 ) tree.mc_truth_h0b1_id = h0b1->pdgId();
   if( h0b2 ) tree.mc_truth_h0b2_id = h0b2->pdgId();
   if( h0b1_IS ) tree.mc_truth_h0b1_IS_id = h0b1_IS->pdgId();
   if( h0b2_IS ) tree.mc_truth_h0b2_IS_id = h0b2_IS->pdgId();
   
   if( t ) tree.mc_truth_t_id = t->pdgId();
   if( tb ) tree.mc_truth_tb_id = tb->pdgId();
   if( tb_IS ) tree.mc_truth_tb_IS_id = tb_IS->pdgId();
   
   if( tW ) tree.mc_truth_tW_id = tW->pdgId();
   if( tWnu ) tree.mc_truth_tWnu_id = tWnu->pdgId();
   if( tWnutau ) tree.mc_truth_tWnutau_id = tWnutau->pdgId();
   if( tWl ) tree.mc_truth_tWl_id = tWl->pdgId();
   if( tWtau ) tree.mc_truth_tWtau_id = tWtau->pdgId();
   if( tWtaunu ) tree.mc_truth_tWtaunu_id = tWtaunu->pdgId();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_id = tWtaunutau->pdgId();
   if( tWtaul ) tree.mc_truth_tWtaul_id = tWtaul->pdgId();
   if( tWq1 ) tree.mc_truth_tWq1_id = tWq1->pdgId();
   if( tWq2 ) tree.mc_truth_tWq2_id = tWq2->pdgId();
   if( tWq1_IS ) tree.mc_truth_tWq1_IS_id = tWq1_IS->pdgId();
   if( tWq2_IS ) tree.mc_truth_tWq2_IS_id = tWq2_IS->pdgId();
   
   if( j1 ) tree.mc_truth_j1_id = j1->pdgId();
   if( j2 ) tree.mc_truth_j2_id = j2->pdgId();
   if( j3 ) tree.mc_truth_j3_id = j3->pdgId();
   
   // status

   if( h0 ) tree.mc_truth_h0_status = h0->status();

   if( h0W1 ) tree.mc_truth_h0W1_status = h0W1->status();
   if( h0W2 ) tree.mc_truth_h0W2_status = h0W2->status();
   if( h0Wl1 ) tree.mc_truth_h0Wl1_status = h0Wl1->status();
   if( h0Wnu1 ) tree.mc_truth_h0Wnu1_status = h0Wnu1->status();
   if( h0Wtau1 ) tree.mc_truth_h0Wtau1_status = h0Wtau1->status();
   if( h0Wnutau1 ) tree.mc_truth_h0Wnutau1_status = h0Wnutau1->status();
   if( h0Wtaul1 ) tree.mc_truth_h0Wtaul1_status = h0Wtaul1->status();
   if( h0Wtaunu1 ) tree.mc_truth_h0Wtaunu1_status = h0Wtaunu1->status();
   if( h0Wtaunutau1 ) tree.mc_truth_h0Wtaunutau1_status = h0Wtaunutau1->status();
   if( h0Wl2 ) tree.mc_truth_h0Wl2_status = h0Wl2->status();
   if( h0Wnu2 ) tree.mc_truth_h0Wnu2_status = h0Wnu2->status();
   if( h0Wtau2 ) tree.mc_truth_h0Wtau2_status = h0Wtau2->status();
   if( h0Wnutau2 ) tree.mc_truth_h0Wnutau2_status = h0Wnutau2->status();
   if( h0Wtaul2 ) tree.mc_truth_h0Wtaul2_status = h0Wtaul2->status();
   if( h0Wtaunu2 ) tree.mc_truth_h0Wtaunu2_status = h0Wtaunu2->status();
   if( h0Wtaunutau2 ) tree.mc_truth_h0Wtaunutau2_status = h0Wtaunutau2->status();
   if( h0Wq11 ) tree.mc_truth_h0Wq11_status = h0Wq11->status();
   if( h0Wq21 ) tree.mc_truth_h0Wq21_status = h0Wq21->status();
   if( h0Wq12 ) tree.mc_truth_h0Wq12_status = h0Wq12->status();
   if( h0Wq22 ) tree.mc_truth_h0Wq22_status = h0Wq22->status();
   if( h0Wq11_IS ) tree.mc_truth_h0Wq11_IS_status = h0Wq11_IS->status();
   if( h0Wq21_IS ) tree.mc_truth_h0Wq21_IS_status = h0Wq21_IS->status();
   if( h0Wq12_IS ) tree.mc_truth_h0Wq12_IS_status = h0Wq12_IS->status();
   if( h0Wq22_IS ) tree.mc_truth_h0Wq22_IS_status = h0Wq22_IS->status();
   
   if( h0Z1 ) tree.mc_truth_h0Z1_status = h0Z1->status();
   if( h0Z2 ) tree.mc_truth_h0Z2_status = h0Z2->status();
   if( h0Zl11 ) tree.mc_truth_h0Zl11_status = h0Zl11->status();
   if( h0Zl21 ) tree.mc_truth_h0Zl21_status = h0Zl21->status();
   if( h0Ztau11 ) tree.mc_truth_h0Ztau11_status = h0Ztau11->status();
   if( h0Ztau21 ) tree.mc_truth_h0Ztau21_status = h0Ztau21->status();
   if( h0Ztaul11 ) tree.mc_truth_h0Ztaul11_status = h0Ztaul11->status();
   if( h0Ztaul21 ) tree.mc_truth_h0Ztaul21_status = h0Ztaul21->status();
   if( h0Ztaunu11 ) tree.mc_truth_h0Ztaunu11_status = h0Ztaunu11->status();
   if( h0Ztaunu21 ) tree.mc_truth_h0Ztaunu21_status = h0Ztaunu21->status();
   if( h0Ztaunutau11 ) tree.mc_truth_h0Ztaunutau11_status = h0Ztaunutau11->status();
   if( h0Ztaunutau21 ) tree.mc_truth_h0Ztaunutau21_status = h0Ztaunutau21->status();
   if( h0Zq11 ) tree.mc_truth_h0Zq11_status = h0Zq11->status();
   if( h0Zq21 ) tree.mc_truth_h0Zq21_status = h0Zq21->status();
   if( h0Zq11_IS ) tree.mc_truth_h0Zq11_IS_status = h0Zq11_IS->status();
   if( h0Zq21_IS ) tree.mc_truth_h0Zq21_IS_status = h0Zq21_IS->status();
   if( h0Zl12 ) tree.mc_truth_h0Zl12_status = h0Zl12->status();
   if( h0Zl22 ) tree.mc_truth_h0Zl22_status = h0Zl22->status();
   if( h0Ztau12 ) tree.mc_truth_h0Ztau12_status = h0Ztau12->status();
   if( h0Ztau22 ) tree.mc_truth_h0Ztau22_status = h0Ztau22->status();
   if( h0Ztaul12 ) tree.mc_truth_h0Ztaul12_status = h0Ztaul12->status();
   if( h0Ztaul22 ) tree.mc_truth_h0Ztaul22_status = h0Ztaul22->status();
   if( h0Ztaunu12 ) tree.mc_truth_h0Ztaunu12_status = h0Ztaunu12->status();
   if( h0Ztaunu22 ) tree.mc_truth_h0Ztaunu22_status = h0Ztaunu22->status();
   if( h0Ztaunutau12 ) tree.mc_truth_h0Ztaunutau12_status = h0Ztaunutau12->status();
   if( h0Ztaunutau22 ) tree.mc_truth_h0Ztaunutau22_status = h0Ztaunutau22->status();
   if( h0Zq12 ) tree.mc_truth_h0Zq12_status = h0Zq12->status();
   if( h0Zq22 ) tree.mc_truth_h0Zq22_status = h0Zq22->status();
   if( h0Zq12_IS ) tree.mc_truth_h0Zq12_IS_status = h0Zq12_IS->status();
   if( h0Zq22_IS ) tree.mc_truth_h0Zq22_IS_status = h0Zq22_IS->status();
   if( h0Znu11 ) tree.mc_truth_h0Znu11_status = h0Znu11->status();
   if( h0Znu21 ) tree.mc_truth_h0Znu21_status = h0Znu21->status();
   if( h0Znu12 ) tree.mc_truth_h0Znu12_status = h0Znu12->status();
   if( h0Znu22 ) tree.mc_truth_h0Znu22_status = h0Znu22->status();
   
   if( h0tau1 ) tree.mc_truth_h0tau1_status = h0tau1->status();
   if( h0tau2 ) tree.mc_truth_h0tau2_status = h0tau2->status();
   if( h0taul1 ) tree.mc_truth_h0taul1_status = h0taul1->status();
   if( h0taunutau1 ) tree.mc_truth_h0taunutau1_status = h0taunutau1->status();
   if( h0taunu1 ) tree.mc_truth_h0taunu1_status = h0taunu1->status();
   if( h0taul2 ) tree.mc_truth_h0taul2_status = h0taul2->status();
   if( h0taunutau2 ) tree.mc_truth_h0taunutau2_status = h0taunutau2->status();
   if( h0taunu2 ) tree.mc_truth_h0taunu2_status = h0taunu2->status();

   if( h0b1 ) tree.mc_truth_h0b1_status = h0b1->status();
   if( h0b2 ) tree.mc_truth_h0b2_status = h0b2->status();
   if( h0b1_IS ) tree.mc_truth_h0b1_IS_status = h0b1_IS->status();
   if( h0b2_IS ) tree.mc_truth_h0b2_IS_status = h0b2_IS->status();
   
   if( t ) tree.mc_truth_t_status = t->status();
   if( tb ) tree.mc_truth_tb_status = tb->status();
   if( tb_IS ) tree.mc_truth_tb_IS_status = tb_IS->status();
   
   if( tW ) tree.mc_truth_tW_status = tW->status();
   if( tWnu ) tree.mc_truth_tWnu_status = tWnu->status();
   if( tWnutau ) tree.mc_truth_tWnutau_status = tWnutau->status();
   if( tWl ) tree.mc_truth_tWl_status = tWl->status();
   if( tWtau ) tree.mc_truth_tWtau_status = tWtau->status();
   if( tWtaunu ) tree.mc_truth_tWtaunu_status = tWtaunu->status();
   if( tWtaunutau ) tree.mc_truth_tWtaunutau_status = tWtaunutau->status();
   if( tWtaul ) tree.mc_truth_tWtaul_status = tWtaul->status();
   if( tWq1 ) tree.mc_truth_tWq1_status = tWq1->status();
   if( tWq2 ) tree.mc_truth_tWq2_status = tWq2->status();
   if( tWq1_IS ) tree.mc_truth_tWq1_IS_status = tWq1_IS->status();
   if( tWq2_IS ) tree.mc_truth_tWq2_IS_status = tWq2_IS->status();

   if( j1 ) tree.mc_truth_j1_status = j1->status();
   if( j2 ) tree.mc_truth_j2_status = j2->status();
   if( j3 ) tree.mc_truth_j3_status = j3->status();
}

reco::GenParticle* MCTruth::getUnique(const reco::GenParticle* p,bool verbose)
{
   reco::GenParticle *pcur = const_cast<reco::GenParticle*>(p);
   
   if( verbose )
     {	
	std::cout << "---------b--------" << std::endl;
	std::cout << "INITIAL = " << pcur->pdgId() << " " << pcur->status() << std::endl;
     }
   
   while( 1 )
     {
	bool foundDupl = false;

	const reco::GenParticleRefVector& daughterRefs = pcur->daughterRefVector();
	for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr)
	  {	     
	     if( idr->isAvailable() )
	       {		  
		  const reco::GenParticleRef& genParticle = (*idr);
		  const reco::GenParticle *d = genParticle.get();

		  if( d )
		    {			       
////		       if( fabs(d->pdgId()) != 15 && d->status() == 2 ) continue;
//		       std::cout << d->pdgId() << " " << d->status() << std::endl;
		       if( verbose )
			 {		  
			    std::cout << "current: " << d->pdgId() << " " << d->status() << std::endl;
			    std::cout << "pcur: " << pcur->pdgId() << " " << pcur->status() << std::endl;
			 }
		       
		       if( d->pdgId() == pcur->pdgId() )
			 {
			    pcur = const_cast<reco::GenParticle*>(d);
			    foundDupl = true;
		       
			    if( verbose )
			      {		       
				 std::cout << "Found duplicate, switch to it" << std::endl;
				 std::cout << "Number of children = " << pcur->numberOfDaughters() << std::endl;
			      }		  
			 }
		    }
		  else break; // the world is fcked up in this case
	       }
	     else break;
	  }
	
	if( !foundDupl ) break;
     }   
   
   if( verbose )
     {   
	std::cout << "FINAL: id=" << pcur->pdgId() << " status=" << pcur->status() << 
	  " daughters=" << pcur->numberOfDaughters() << std::endl;
	
	const reco::GenParticleRefVector& daughterRefs = pcur->daughterRefVector();
	int ip = 0;
	for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr)
	  {	     
	     if( idr->isAvailable() )
	       {		  
		  const reco::GenParticleRef& genParticle = (*idr);
		  const reco::GenParticle *d = genParticle.get();
		  std::cout << "daughter #" << ip << ": id=" << d->pdgId() << std::endl;
		  ip++;
	       }
	  }	
	std::cout << "---------e--------" << std::endl;
     }
      
   return pcur;
}

void MCTruth::p4toTLV(reco::Particle::LorentzVector vp4,TLorentzVector& tlv)
{
   return tlv.SetPxPyPzE(vp4.px(),vp4.py(),vp4.pz(),vp4.energy());
}

const reco::GenParticle* MCTruth::getMother(const reco::GenParticle &part)
{
   const reco::GenParticle *mom = &part;
   while( mom->numberOfMothers() > 0 )
     {	     
	for( unsigned int j=0;j<mom->numberOfMothers();++j )
	  {		  
	     mom = dynamic_cast<const reco::GenParticle*>(mom->mother(j));
	     if( mom->pdgId() != part.pdgId() )
	       {		       
		  return mom;
	       }
	  }	     
     }
   
   return mom;
}
