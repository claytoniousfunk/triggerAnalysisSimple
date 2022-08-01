

const unsigned int maxJets = 10000;

class jet{

public:
    jet(){};
    ~jet(){};

    void setupTreeForReading(TTree *t);
    void setupTreeForWriting(TTree *t);

    // Declaration of leaf types
    Float_t jtpt[maxJets];
    Float_t jteta[maxJets];
    Int_t njet;

    // list of branches
    TBranch *b_jtpt;
    TBranch *b_jteta;
    TBranch *b_njet;



};

void jet::setupTreeForReading(TTree *t){

    b_jtpt = 0;
    b_jteta = 0;
    b_njet = 0;

    // set branch address and branch pointers
    if(t->GetBranch("nref")) t->SetBranchAddress("nref",&njet,&b_njet);
    if(t->GetBranch("jtpt")) t->SetBranchAddress("jtpt",&jtpt,&b_jtpt);
    if(t->GetBranch("jteta")) t->SetBranchAddress("jteta",&jteta,&b_jteta);
    





}

void jet::setupTreeForWriting(TTree *t){

    t->Branch("nref",&njet,"njet/F");
    t->Branch("jtpt",&jtpt,"jtpt[nref]/F");
    t->Branch("jteta",&jteta,"jteta[nref]/F");
   

}
