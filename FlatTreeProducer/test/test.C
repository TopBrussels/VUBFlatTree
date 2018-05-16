{
   TFile *f = TFile::Open("output.root");
   
   TLorentzVector *v1 = 0;
   int mc_truth_tzq_channel;
   
   TTree *tr = (TTree*)f->Get("FlatTree/tree");
   tr->SetBranchAddress("mc_truth_tzq_Z_p4",&v1);
   tr->SetBranchAddress("mc_truth_tzq_channel",&mc_truth_tzq_channel);
   
   int nev = tr->GetEntries();
//   std::cout << nev << std::endl;
   for(int i=0;i<nev;i++)
     {
	tr->GetEntry(i);
//	std::cout << mc_truth_tzq_channel << std::endl;
	
	std::cout << mc_truth_tzq_channel << " " << v1->Pt() << std::endl;
     }
   
   gApplication->Terminate();
}
