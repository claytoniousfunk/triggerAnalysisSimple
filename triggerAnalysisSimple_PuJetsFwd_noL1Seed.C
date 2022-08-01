using namespace std;

Float_t ptMin = 0.;
Float_t ptMax = 300.;
Int_t Nbins = 60;

TH1D *denom = new TH1D("denom","denom",Nbins,ptMin,ptMax);

TH1D *num_40 = new TH1D("num_40","num_40",Nbins,ptMin,ptMax);
TH1D *r_40 = new TH1D("r_40","r_40",Nbins,ptMin,ptMax);

TH1D *num_60 = new TH1D("num_60","num_60",Nbins,ptMin,ptMax);
TH1D *r_60 = new TH1D("r_60","r_60",Nbins,ptMin,ptMax);

TH1D *num_80 = new TH1D("num_80","num_80",Nbins,ptMin,ptMax);
TH1D *r_80 = new TH1D("r_80","r_80",Nbins,ptMin,ptMax);

TH1D *num_100 = new TH1D("num_100","num_100",Nbins,ptMin,ptMax);
TH1D *r_100 = new TH1D("r_100","r_100",Nbins,ptMin,ptMax);

TH1D *num_120 = new TH1D("num_120","num_120",Nbins,ptMin,ptMax);
TH1D *r_120 = new TH1D("r_120","r_120",Nbins,ptMin,ptMax);



std::map<unsigned long long, int> runLumiEvtToEntryMap;



unsigned long long keyFromRunLumiEvent(Int_t run,
				       Int_t lumi,
				       ULong64_t event);




void triggerAnalysisSimple_PuJetsFwd_noL1Seed(std::string triggerFile = "openHLT_CB2_merge.root",
                           std::string inputFile = "HiForestMiniAOD_CB_merge.root",
                           std::string outputFile = "rootFiles/out_PuJetsEta5p1.root"){


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
    treeTrig->SetBranchStatus("HLT_HIPuAK4CaloJet40Fwd_noL1Seed_v", 1);
    treeTrig->SetBranchStatus("HLT_HIPuAK4CaloJet60Fwd_noL1Seed_v", 1);
    treeTrig->SetBranchStatus("HLT_HIPuAK4CaloJet80Fwd_noL1Seed_v", 1);
    treeTrig->SetBranchStatus("HLT_HIPuAK4CaloJet100Fwd_noL1Seed_v", 1);
    treeTrig->SetBranchStatus("HLT_HIPuAK4CaloJet120Fwd_noL1Seed_v", 1);

    //Int_t hlt_event;
    ULong64_t       hlt_event;
    Int_t           hlt_lumi;
    Int_t           hlt_run;
    Bool_t          triggerDecision_40;
    Bool_t          triggerDecision_60;
    Bool_t          triggerDecision_80;
    Bool_t          triggerDecision_100;
    Bool_t          triggerDecision_120;
    std::cout << "Setting Event, lumi, and run branchAdresses...";
    treeTrig->SetBranchAddress("Event", &hlt_event);
    treeTrig->SetBranchAddress("LumiBlock", &hlt_lumi);
    treeTrig->SetBranchAddress("Run", &hlt_run);
    treeTrig->SetBranchAddress("HLT_HIPuAK4CaloJet40Fwd_noL1Seed_v",&triggerDecision_40);
    treeTrig->SetBranchAddress("HLT_HIPuAK4CaloJet60Fwd_noL1Seed_v",&triggerDecision_60);
    treeTrig->SetBranchAddress("HLT_HIPuAK4CaloJet80Fwd_noL1Seed_v",&triggerDecision_80);
    treeTrig->SetBranchAddress("HLT_HIPuAK4CaloJet100Fwd_noL1Seed_v",&triggerDecision_100);
    treeTrig->SetBranchAddress("HLT_HIPuAK4CaloJet120Fwd_noL1Seed_v",&triggerDecision_120);
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
	
        
        auto j_entry_status = treeJet->GetEntry(j_entry);
        
        if(j_entry_status < 0) {
            std::cout << Form("bad entry %lld in jet tree",j_entry) << std::endl;
            continue;
        }
	
        unsigned long long key = keyFromRunLumiEvent(run,lumi,evt);

        long long i_entry = -1;
        
        if(runLumiEvtToEntryMap.count(key) == 0) continue; // skip reco event if there is no HLT event match
        else i_entry = runLumiEvtToEntryMap.at(key);

        //std::cout << "i_entry = " << i_entry << std::endl;
        
        // now fill the denominator
        Float_t maxPt_denom = 0;
        Float_t maxEta_denom = 0;
        
        for(Int_t i_jet = 0; i_jet < nref; i_jet++){

	   	    
            if(jtpt[i_jet] > maxPt_denom) { // find the leading jetPt in the event, regardless of trigger.
                maxPt_denom = jtpt[i_jet];
                maxEta_denom = jteta[i_jet];
            }


	    }

        if(fabs(maxEta_denom)>4.7 || fabs(maxEta_denom)<3.2) continue; // skip event if the leading jet is outside eta range

	    if(maxPt_denom > 0) {
            
            denom->Fill(maxPt_denom,weight);
            //std::cout << "maxPt = " << maxPt_denom <<std::endl;
        }





        treeTrig->GetEntry(i_entry); // get trigger decision from HLT emulation

        
        
        if(triggerDecision_40) {// only fill the numerator if the trigger is on.

            Float_t maxPt_num = 0;
            
            // now fill the numerator
            for(Int_t i_jet = 0; i_jet < nref; i_jet++){

                // no eta cut needed since already applied after the first jet loop.  

                if(jtpt[i_jet] > maxPt_num) { // find the leading jetPt in events with trigger on.
                    maxPt_num = jtpt[i_jet];
                }


            }

            if(maxPt_num > 0){

                num_40->Fill(maxPt_num,weight);
                //std::cout << "maxPt_num = " << maxPt_num << std::endl <<  std::endl;
            } 

        }

        
        
        if(triggerDecision_60) {// only fill the numerator if the trigger is on.

            Float_t maxPt_num = 0;
            
            // now fill the numerator
            for(Int_t i_jet = 0; i_jet < nref; i_jet++){

                // no eta cut needed since already applied after the first jet loop.  

                if(jtpt[i_jet] > maxPt_num) { // find the leading jetPt in events with trigger on.
                    maxPt_num = jtpt[i_jet];
                }


            }

            if(maxPt_num > 0){

                num_60->Fill(maxPt_num,weight);
                //std::cout << "maxPt_num = " << maxPt_num << std::endl <<  std::endl;
            } 

        }

        if(triggerDecision_80) {// only fill the numerator if the trigger is on.

            Float_t maxPt_num = 0;
            
            // now fill the numerator
            for(Int_t i_jet = 0; i_jet < nref; i_jet++){

                // no eta cut needed since already applied after the first jet loop.  

                if(jtpt[i_jet] > maxPt_num) { // find the leading jetPt in events with trigger on.
                    maxPt_num = jtpt[i_jet];
                }


            }

            if(maxPt_num > 0){

                num_80->Fill(maxPt_num,weight);
                //std::cout << "maxPt_num = " << maxPt_num << std::endl <<  std::endl;
            } 

        }

        if(triggerDecision_100) {// only fill the numerator if the trigger is on.

            Float_t maxPt_num = 0;
            
            // now fill the numerator
            for(Int_t i_jet = 0; i_jet < nref; i_jet++){

                // no eta cut needed since already applied after the first jet loop.  

                if(jtpt[i_jet] > maxPt_num) { // find the leading jetPt in events with trigger on.
                    maxPt_num = jtpt[i_jet];
                }


            }

            if(maxPt_num > 0){

                num_100->Fill(maxPt_num,weight);
                //std::cout << "maxPt_num = " << maxPt_num << std::endl <<  std::endl;
            } 

        }

        if(triggerDecision_120) {// only fill the numerator if the trigger is on.

            Float_t maxPt_num = 0;
            
            // now fill the numerator
            for(Int_t i_jet = 0; i_jet < nref; i_jet++){

                // no eta cut needed since already applied after the first jet loop.  

                if(jtpt[i_jet] > maxPt_num) { // find the leading jetPt in events with trigger on.
                    maxPt_num = jtpt[i_jet];
                }


            }

            if(maxPt_num > 0){

                num_120->Fill(maxPt_num,weight);
                //std::cout << "maxPt_num = " << maxPt_num << std::endl <<  std::endl;
            } 

        }




    }

    r_40->Divide(num_40,denom,1,1,"");
    r_60->Divide(num_60,denom,1,1,"");
    r_80->Divide(num_80,denom,1,1,"");
    r_100->Divide(num_100,denom,1,1,"");
    r_120->Divide(num_120,denom,1,1,"");

    r_40->SetLineColor(kRed);
    r_60->SetLineColor(kBlue);
    r_80->SetLineColor(kGreen);
    r_100->SetLineColor(kMagenta);
    r_120->SetLineColor(kCyan+2);

    r_40->SetMarkerColor(kRed);
    r_60->SetMarkerColor(kBlue);
    r_80->SetMarkerColor(kGreen);
    r_100->SetMarkerColor(kMagenta);
    r_120->SetMarkerColor(kCyan+2);

    r_40->SetMarkerStyle(24);
    r_60->SetMarkerStyle(25);
    r_80->SetMarkerStyle(26);
    r_100->SetMarkerStyle(32);
    r_120->SetMarkerStyle(28);

    TCanvas *c1 = new TCanvas("c1","c1",700,600);
    c1->cd();
    TPad *p1 = new TPad("p1","p1",0,0,1,1);
    p1->SetLeftMargin(0.13);
    p1->SetBottomMargin(0.14);
    p1->Draw();
    p1->cd();
    r_40->SetTitle("");
    r_40->SetStats(0);
    r_40->GetXaxis()->SetTitleSize(0.05);
    r_40->GetYaxis()->SetTitleSize(0.05);
    r_40->GetXaxis()->SetTitle("leading caloJet #font[52]{p}_{T} [GeV]");
    r_40->GetYaxis()->SetTitle("Trigger efficiency");
    r_40->GetYaxis()->SetRangeUser(0.0,2.0);
    TLegend *leg = new TLegend(0.48,0.72,0.88,0.88);
    leg->AddEntry(r_40,"HLT_HIPuAK4CaloJet40Fwd_noL1Seed_v","p");
    leg->AddEntry(r_60,"HLT_HIPuAK4CaloJet60Fwd_noL1Seed_v","p");
    leg->AddEntry(r_80,"HLT_HIPuAK4CaloJet80Fwd_noL1Seed_v","p");
    leg->AddEntry(r_100,"HLT_HIPuAK4CaloJet100Fwd_noL1Seed_v","p");
    leg->AddEntry(r_120,"HLT_HIPuAK4CaloJet120Fwd_noL1Seed_v","p");
    leg->SetTextSize(0.021);
    leg->SetBorderSize(0);
    r_40->Draw();
    leg->Draw();
    r_60->Draw("same");
    r_80->Draw("same");
    r_100->Draw("same");
    r_120->Draw("same");




    TF1 *fit_40 = new TF1("fit_40","0.5+0.5*tanh([0]*(x - [1]))",0,500);
    fit_40->SetParameter(0,1.8e-01);
    fit_40->SetParameter(1,40);
    fit_40->SetLineStyle(7);
    fit_40->SetLineColor(kRed);

    TF1 *fit_60 = new TF1("fit_60","0.5+0.5*tanh([0]*(x - [1]))",0,500);
    fit_60->SetParameter(0,7.04514e-02);
    fit_60->SetLineStyle(7);
    fit_60->SetLineColor(kBlue);

    TF1 *fit_80 = new TF1("fit_80","0.5+0.5*tanh([0]*(x - [1]))",0,500);
    fit_80->SetLineStyle(7);
    fit_80->SetLineColor(kGreen);

    TF1 *fit_100 = new TF1("fit_100","0.5+0.5*tanh([0]*(x - [1]))",0,500);
    fit_100->SetLineStyle(7);
    fit_100->SetLineColor(kMagenta);

    TF1 *fit_120 = new TF1("fit_120","0.5+0.5*tanh([0]*(x - [1]))",0,500);
    fit_120->SetParameter(0,7.04514e-02);
    fit_120->SetParameter(1,120);
    fit_120->SetLineStyle(7);
    fit_120->SetLineColor(kCyan+2);

    r_40->Fit("fit_40","R N");
    r_60->Fit("fit_60","R N");
    r_80->Fit("fit_80","R N");
    r_100->Fit("fit_100","R N");
    r_120->Fit("fit_120","R N");


     // fit_40->Draw("same");
     // fit_60->Draw("same");
     // fit_80->Draw("same");
     // fit_100->Draw("same");
     // fit_120->Draw("same");














    TLatex *la = new TLatex();
    la->SetTextFont(42);
    la->SetTextSize(0.03);

    la->DrawLatexNDC(0.18,0.85,"PYTHIA+HYDJET 0-90%");
    la->DrawLatexNDC(0.18,0.79,"Run 3 MC");
    la->DrawLatexNDC(0.18,0.73,"3.2 < |#eta^{jet}| < 4.7, no L1 seed");

 


    auto wf = TFile::Open("out.root","recreate");

        denom->Write();
        
        num_40->Write();
        num_60->Write();
        num_80->Write();
        num_100->Write();
        num_120->Write();

        r_40->Write();
        r_60->Write();
        r_80->Write();
        r_100->Write();
        r_120->Write();

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
