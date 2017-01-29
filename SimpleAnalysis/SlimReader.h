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
   Float_t         met_pt;
   Float_t         met_phi;

   vector<float>   *obj_pt[NUMTYPES];
   vector<float>   *obj_eta[NUMTYPES];
   vector<float>   *obj_phi[NUMTYPES];
   vector<int>     *obj_charge[NUMTYPES];
   vector<int>     *obj_id[NUMTYPES];
   vector<float>   *jet_m;
   vector<float>   *fatjet_m;

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_phi;   //!
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
  for(int ii=0;ii<NUMTYPES;ii++) {
   obj_pt[ii] = 0;
   obj_eta[ii] = 0;
   obj_phi[ii] = 0;
   obj_charge[ii] = 0;
   obj_id[ii] = 0;
  }
  jet_m = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);
  fChain->SetBranchAddress("Event", &EventNumber, &b_EventNumber);
  fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
  fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
  int idx=0;
  std::vector<std::string> bNames({"el","mu","tau","photon","jet","fatjet"});
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

