#ifndef SlimReader_h
#define SlimReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include "SimpleAnalysis/AnalysisRunner.h"

using std::vector;

// Fixed size dimensions of array or collections stored in the TTree if any.

#define NUMTYPES 6

class SlimReaderSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           EventNumber;
   Int_t           mcChannel;
   Int_t           susyChannel;
   vector<float>   *mcWeights;
   Float_t         sumet;
   Float_t         met_pt;
   Float_t         met_phi;
   Float_t         genMET;
   Float_t         genHT;

   Int_t           pdf_id1;
   Float_t         pdf_x1;
   Float_t         pdf_pdf1;
   Int_t           pdf_id2;
   Float_t         pdf_x2;
   Float_t         pdf_pdf2;
   Float_t         pdf_scale;
   
   vector<float>   *obj_pt[NUMTYPES];
   vector<float>   *obj_eta[NUMTYPES];
   vector<float>   *obj_phi[NUMTYPES];
   vector<int>     *obj_charge[NUMTYPES];
   vector<int>     *obj_id[NUMTYPES];
   vector<float>   *jet_m;
   vector<float>   *fatjet_m;

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_mcChannel;   //!
   TBranch        *b_susyChannel;   //!
   TBranch        *b_mcWeights;   //!
   TBranch        *b_sumet;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_gen_met;   //!
   TBranch        *b_gen_ht;   //!

   TBranch        *b_pdf_id1;   //!
   TBranch        *b_pdf_x1;    //!
   TBranch        *b_pdf_pdf1;  //!
   TBranch        *b_pdf_id2;   //!
   TBranch        *b_pdf_x2;    //!
   TBranch        *b_pdf_pdf2;  //!
   TBranch        *b_pdf_scale; //!

   TBranch        *b_obj_pt[NUMTYPES];   //!
   TBranch        *b_obj_eta[NUMTYPES];   //!
   TBranch        *b_obj_phi[NUMTYPES];   //!
   TBranch        *b_obj_charge[NUMTYPES];   //!
   TBranch        *b_obj_id[NUMTYPES];   //!
   TBranch        *b_jet_m;   //!
   TBranch        *b_fatjet_m;   //!

   AnalysisRunner* _runner;

 SlimReaderSelector(AnalysisRunner *runner) : fChain(0),_runner(runner) { }
   virtual ~SlimReaderSelector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(SlimReaderSelector,0);
};

class SlimReader : public Reader {
 public:
 SlimReader(std::vector<AnalysisClass*>& analysisList) : Reader(analysisList) { };
  
 protected:
  void processFilesInternal(const std::vector<std::string>& inputNames);
};

#endif

#ifdef SlimReader_cxx
void SlimReaderSelector::Init(TTree *tree)
{
  // Set object pointer
  mcWeights = 0;
  for(int ii=0;ii<NUMTYPES;ii++) {
   obj_pt[ii] = 0;
   obj_eta[ii] = 0;
   obj_phi[ii] = 0;
   obj_charge[ii] = 0;
   obj_id[ii] = 0;
  }
  jet_m = 0;
  fatjet_m = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("Event", &EventNumber, &b_EventNumber);
  fChain->SetBranchAddress("mcChannel", &mcChannel, &b_mcChannel);
  fChain->SetBranchAddress("susyChannel", &susyChannel, &b_susyChannel);
  fChain->SetBranchAddress("mcWeights", &mcWeights, &b_mcWeights);
  fChain->SetBranchAddress("sumet", &sumet, &b_sumet);
  fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
  fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
  fChain->SetBranchAddress("genMET", &genMET, &b_gen_met);
  fChain->SetBranchAddress("genHT", &genHT, &b_gen_ht);
  fChain->SetBranchAddress("pdf_id1", &pdf_id1, &b_pdf_id1);
  fChain->SetBranchAddress("pdf_x1", &pdf_x1, &b_pdf_x1);
  fChain->SetBranchAddress("pdf_pdf1", &pdf_pdf1, &b_pdf_pdf1);
  fChain->SetBranchAddress("pdf_id2", &pdf_id2, &b_pdf_id2);
  fChain->SetBranchAddress("pdf_x2", &pdf_x2, &b_pdf_x2);
  fChain->SetBranchAddress("pdf_pdf2", &pdf_pdf2, &b_pdf_pdf2);
  fChain->SetBranchAddress("pdf_scale", &pdf_scale, &b_pdf_scale);

  int idx=0;
  std::vector<std::string> bNames({"el","mu","tau","ph","jet","fatjet"});
  for(const auto& bName : bNames) {
    fChain->SetBranchAddress((bName+"_pt").c_str(),     &obj_pt[idx], &b_obj_pt[idx]);
    fChain->SetBranchAddress((bName+"_eta").c_str(),    &obj_eta[idx], &b_obj_eta[idx]);
    fChain->SetBranchAddress((bName+"_phi").c_str(),    &obj_phi[idx], &b_obj_phi[idx]);
    fChain->SetBranchAddress((bName+"_charge").c_str(), &obj_charge[idx], &b_obj_charge[idx]);
    fChain->SetBranchAddress((bName+"_id").c_str(),     &obj_id[idx], &b_obj_id[idx]);
    idx++;
  }
  fChain->SetBranchAddress("jet_m", &jet_m, &b_jet_m);
  fChain->SetBranchAddress("fatjet_m", &fatjet_m, &b_fatjet_m);
}

Bool_t SlimReaderSelector::Notify()
{
   return kTRUE;
}

#endif // #ifdef SlimReader_cxx

