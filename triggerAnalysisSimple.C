using namespace std;

TH1D *denom = new TH1D("denom","denom",20,0,120);
TH1D *num = new TH1D("num","num",20,0,120);



std::map<unsigned long long, int> runLumiEvtToEntryMap;



unsigned long long keyFromRunLumiEvent(Int_t run,
				       Int_t lumi,
				       ULong64_t event);




void triggerAnalysisSimple(std::string triggerFile = "openHLT.root",
                           std::string inputFile = "pythiahydjet.root",
                           std::string outputFile = "out.root"){


    std::cout << "running triggerAnalysis()" << std::endl;
    std::cout << "inputFile   = " << inputFile.c_str()  << std::endl;
    std::cout << "outputFile  = " << outputFile.c_str() << std::endl;	


    std::cout << "### input file ###" << std::endl; 
    TFile* input = TFile::Open(inputFile.c_str(), "READ");


    TTree* treeggHiNtuplizer = 0;
    TTree* treeJet = 0;
    TTree* treeHiEvt = 0;
    TTree* treeTrig = 0;
        
    TFile* fileTmp = 0;
    TFile* fileTrig = 0;
    
    std::cout << "### HLT bit analysis file ###" << std::endl;
    fileTrig = TFile::Open(triggerFile.c_str(), "READ");
    fileTrig->cd();
    
    std::string treeTrigPath = "hltanalysis/HltTree";
    treeTrig = (TTree*)fileTrig->Get(treeTrigPath.c_str());
    treeTrig->SetBranchStatus("*",0);     // disable all branches

    // specify explicitly which branches to use
    treeTrig->SetBranchStatus("Event", 1);
    treeTrig->SetBranchStatus("LumiBlock", 1);
    treeTrig->SetBranchStatus("Run", 1);

    //Int_t hlt_event;
    ULong64_t       hlt_event;
    Int_t           hlt_lumi;
    Int_t           hlt_run;
    std::cout << "Setting Event, lumi, and run branchAdresses...";
    treeTrig->SetBranchAddress("Event", &hlt_event);
    treeTrig->SetBranchAddress("LumiBlock", &hlt_lumi);
    treeTrig->SetBranchAddress("Run", &hlt_run);
    std::cout << "done" << std::endl;


    std::cout << "get number of entries from HLT file..." << std::endl;

    Long64_t entriesHLT = treeTrig->GetEntries();
    std::cout << "HLT entries = " << entriesHLT << std::endl;
    
    
    
    
    std::cout << "get number of entries from reco file.." << std::endl;
    // read the first file only to get the HiForest info
    std::string inputPath = inputFile.c_str();
    fileTmp = TFile::Open(inputPath.c_str(), "READ");
    fileTmp->cd();


    std::string treePath = "ggHiNtuplizer/EventTree";

    // read one tree only to get the number of entries
    treeggHiNtuplizer = (TTree*)fileTmp->Get(treePath.c_str());
    Long64_t entriesTmp = treeggHiNtuplizer->GetEntries();
    std::cout << "reco entries = " << entriesTmp << std::endl;
 //   treeggHiNtuplizer->Delete();

 //   treeggHiNtuplizer = (TTree*)fileTmp->Get(treePath.c_str());
    treeggHiNtuplizer->SetBranchStatus("*",0);     // disable all branches
    
    std::map<unsigned long long, int> runLumiEvtToEntryMap;
    
    treeJet = (TTree*)fileTmp->Get("akCs4PFJetAnalyzer/t");
    treeJet->SetBranchStatus("*",0);     // disable all branches
    treeJet->SetBranchStatus("jtpt",1);   // enable event information
    treeJet->SetBranchStatus("jteta",1);
    treeJet->SetBranchStatus("nref",1);

    const unsigned int maxJets = 10000;
    
    Float_t jtpt[maxJets];
    Float_t jteta[maxJets];
    Int_t nref;

    treeJet->SetBranchAddress("jtpt",&jtpt);
    treeJet->SetBranchAddress("jteta",&jteta);
    treeJet->SetBranchAddress("nref",&nref);
    
    
    treeHiEvt = (TTree*)fileTmp->Get("hiEvtAnalyzer/HiTree");
    treeHiEvt->SetBranchStatus("*",0);     // disable all branches
    treeHiEvt->SetBranchStatus("run",1);   // enable event information
    treeHiEvt->SetBranchStatus("evt",1);
    treeHiEvt->SetBranchStatus("lumi",1);
    treeHiEvt->SetBranchStatus("vz",1);
    treeHiEvt->SetBranchStatus("hiBin",1);

    Float_t vz;
    Int_t hiBin;
    UInt_t run;
    UInt_t lumi;
    ULong64_t evt;
    Float_t weight;



    treeHiEvt->SetBranchAddress("vz",&vz);
    treeHiEvt->SetBranchAddress("hiBin",&hiBin);
    treeHiEvt->SetBranchAddress("run",&run);
    treeHiEvt->SetBranchAddress("lumi",&lumi);
    treeHiEvt->SetBranchAddress("evt",&evt);
    treeHiEvt->SetBranchAddress("weight",&weight);




    // loop through HLT and create a key for each event

    for(Long64_t i_entry = 0; i_entry < entriesHLT; i_entry++){

       treeTrig->GetEntry(i_entry);
       unsigned long long key = keyFromRunLumiEvent(hlt_run, hlt_lumi, hlt_event);
       runLumiEvtToEntryMap[key] = i_entry;

    }




    // loop through reco objects
    for (Long64_t j_entry = 0; j_entry < entriesTmp; ++j_entry){

        treeggHiNtuplizer->GetEntry(j_entry);
        treeHiEvt->GetEntry(j_entry);
        treeJet->GetEntry(j_entry);
    
	
	// event cuts
 	if(fabs(vz)>15.0) continue;
        if(hiBin>180) continue;	
	
	
	Float_t maxPt_denom = 0;
	// fill denom
	for(Int_t i_jet = 0; i_jet < nref; i_jet++){

	    // jet kinematic cuts
	    if(fabs(jteta[i_jet]) < 1.6) continue;
	    
	    if(jtpt[i_jet] > maxPt_denom) maxPt_denom = jtpt[i_jet];


	}

	if(maxPt_denom>0) denom->Fill(maxPt_denom,weight);

	auto j_entry_status = treeJet->GetEntry(j_entry);
	if(j_entry_status < 0) {
	    std::cout << Form("bad entry %lld in jet tree",j_entry) << std::endl;
	    continue;
	}
	
	unsigned long long key = keyFromRunLumiEvent(run,lumi,evt);

	long long i_entry = -1;
        if(runLumiEvtToEntryMap.count(key) == 0) continue;
	else i_entry = runLumiEvtToEntryMap.at(key);

	//std::cout << "i_entry = " << i_entry << std::endl;
	


	Float_t maxPt_num = 0;
	// fill denom
	for(Int_t i_jet = 0; i_jet < nref; i_jet++){

	    // jet kinematic cuts
	    if(fabs(jteta[i_jet]) < 1.6) continue;
	    
	    if(jtpt[i_jet] > maxPt_num) maxPt_num = jtpt[i_jet];


	}
	if(maxPt_num>0) num->Fill(maxPt_num,weight);
		


    }


auto wf = TFile::Open("out.root","recreate");

num->Write();
denom->Write();

wf->Close();



}


unsigned long long keyFromRunLumiEvent(Int_t run,
                                       Int_t lumi,
                                       ULong64_t event)
{
  const unsigned long long runMult = 1;
  const unsigned long long lumiMult = 1000000;
  const unsigned long long evtMult = 10000000000;
  const unsigned long long evtLimit = 10000000000;

  unsigned long long key = 0;
  if(event >= evtLimit){
    std::cout << "RUNLUMIEVENTKEY WARNING : \'" << event
              << "\' is greated that event limit 10^10. returning key 0"
              << std::endl;
    return key;
  }

  key += runMult* static_cast<unsigned long long>(run);
  key += lumiMult*static_cast<unsigned long long>(lumi);
  key += evtMult*event;

  //std::cout << "key = " << key << std::endl;
  
  return key;
  
}
