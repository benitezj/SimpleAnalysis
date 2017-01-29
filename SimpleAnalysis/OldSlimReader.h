//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 22 05:42:19 2015 by ROOT version 5.34/14
// from TTree truth/OldSlimmed truth tree
// found on file: slim.root
//////////////////////////////////////////////////////////

#ifndef OldSlimReader_h
#define OldSlimReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include "SimpleAnalysis/AnalysisRunner.h"

using std::vector;

// Fixed size dimensions of array or collections stored in the TTree if any.

class OldSlimReaderSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           EventNumber;
   Int_t           susy_part_id1;
   Int_t           susy_part_id2;
   Int_t           mcevt_pdf_id1;
   Int_t           mcevt_pdf_id2;
   Double_t        mcevt_pdf_x1;
   Double_t        mcevt_pdf_x2;
   Double_t        mcevt_pdf1;
   Double_t        mcevt_pdf2;
   Double_t        mcevt_pdf_scale;
   Float_t         MET_Truth_NonInt_etx;
   Float_t         MET_Truth_NonInt_ety;

   Int_t           llv_n;
   vector<int>     *llv_id;
   vector<float>   *llv_pvx;
   vector<float>   *llv_pvy;
   vector<float>   *llv_pvz;
   vector<float>   *llv_evx;
   vector<float>   *llv_evy;
   vector<float>   *llv_evz;
   Int_t           ph_n;
   vector<float>   *ph_pt;
   vector<float>   *ph_eta;
   vector<float>   *ph_phi;
   Int_t           el_n;
   vector<int>     *el_charge;
   vector<float>   *el_pt;
   vector<float>   *el_eta;
   vector<float>   *el_phi;
   vector<int>     *el_origin;
   Int_t           mu_n;
   vector<int>     *mu_charge;
   vector<float>   *mu_pt;
   vector<float>   *mu_eta;
   vector<float>   *mu_phi;
   vector<int>     *mu_origin;
   Int_t           tau_n;
   vector<int>     *tau_charge;
   vector<float>   *tau_pt;
   vector<float>   *tau_eta;
   vector<float>   *tau_phi;
   vector<int>     *tau_origin;
   Int_t           jet_AntiKt4TruthJets_n;
   vector<float>   *jet_AntiKt4TruthJets_E;
   vector<float>   *jet_AntiKt4TruthJets_pt;
   vector<float>   *jet_AntiKt4TruthJets_m;
   vector<float>   *jet_AntiKt4TruthJets_eta;
   vector<float>   *jet_AntiKt4TruthJets_phi;
   vector<int>     *jet_AntiKt4TruthJets_flavor;
   Int_t           jet_AntiKt10TruthTrimmed_n;
   vector<float>   *jet_AntiKt10TruthTrimmed_E;
   vector<float>   *jet_AntiKt10TruthTrimmed_pt;
   vector<float>   *jet_AntiKt10TruthTrimmed_m;
   vector<float>   *jet_AntiKt10TruthTrimmed_eta;
   vector<float>   *jet_AntiKt10TruthTrimmed_phi;
   vector<int>     *jet_AntiKt10TruthTrimmed_flavor;

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_susy_part_id1;   //!
   TBranch        *b_susy_part_id2;   //!
   TBranch        *b_mcevt_pdf_id1;   //!
   TBranch        *b_mcevt_pdf_id2;   //!
   TBranch        *b_mcevt_pdf_x1;   //!
   TBranch        *b_mcevt_pdf_x2;   //!
   TBranch        *b_mcevt_pdf1;   //!
   TBranch        *b_mcevt_pdf2;   //!
   TBranch        *b_mcevt_pdf_scale;   //!
   TBranch        *b_MET_Truth_NonInt_etx;   //!
   TBranch        *b_MET_Truth_NonInt_ety;   //!
   TBranch        *b_llv_n;   //!
   TBranch        *b_llv_id;   //!
   TBranch        *b_llv_pvx;   //!
   TBranch        *b_llv_pvy;   //!
   TBranch        *b_llv_pvz;   //!
   TBranch        *b_llv_evx;   //!
   TBranch        *b_llv_evy;   //!
   TBranch        *b_llv_evz;   //!
   TBranch        *b_ph_n;   //!
   TBranch        *b_ph_pt;   //!
   TBranch        *b_ph_eta;   //!
   TBranch        *b_ph_phi;   //!
   TBranch        *b_el_n;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_origin;   //!
   TBranch        *b_mu_n;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_origin;   //!
   TBranch        *b_tau_n;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_tau_origin;   //!
   TBranch        *b_jet_AntiKt4TruthJets_n;   //!
   TBranch        *b_jet_AntiKt4TruthJets_E;   //!
   TBranch        *b_jet_AntiKt4TruthJets_pt;   //!
   TBranch        *b_jet_AntiKt4TruthJets_m;   //!
   TBranch        *b_jet_AntiKt4TruthJets_eta;   //!
   TBranch        *b_jet_AntiKt4TruthJets_phi;   //!
   TBranch        *b_jet_AntiKt4TruthJets_flavor;   //!
   TBranch        *b_jet_AntiKt10TruthTrimmed_n;   //!
   TBranch        *b_jet_AntiKt10TruthTrimmed_E;   //!
   TBranch        *b_jet_AntiKt10TruthTrimmed_pt;   //!
   TBranch        *b_jet_AntiKt10TruthTrimmed_m;   //!
   TBranch        *b_jet_AntiKt10TruthTrimmed_eta;   //!
   TBranch        *b_jet_AntiKt10TruthTrimmed_phi;   //!
   TBranch        *b_jet_AntiKt10TruthTrimmed_flavor;   //!

   AnalysisRunner* _runner;

   int xAODTruth; //!

 OldSlimReaderSelector(AnalysisRunner *runner) : fChain(0),_runner(runner) { }
   virtual ~OldSlimReaderSelector() { }
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

   ClassDef(OldSlimReaderSelector,0);
};

class OldSlimReader : public Reader {
 public:
 OldSlimReader(std::vector<AnalysisClass*>& analysisList) : Reader(analysisList) { };
  
 protected:
  void processFilesInternal(const std::vector<std::string>& inputNames);
};

#endif

#ifdef OldSlimReader_cxx
void OldSlimReaderSelector::Init(TTree *tree)
{
   // Set object pointer
   llv_id = 0;
   llv_pvx = 0;
   llv_pvy = 0;
   llv_pvz = 0;
   llv_evx = 0;
   llv_evy = 0;
   llv_evz = 0;
   ph_pt = 0;
   ph_eta = 0;
   ph_phi = 0;
   el_charge = 0;
   el_pt = 0;
   el_eta = 0;
   el_phi = 0;
   el_origin = 0;
   mu_charge = 0;
   mu_pt = 0;
   mu_eta = 0;
   mu_phi = 0;
   mu_origin = 0;
   tau_charge = 0;
   tau_pt = 0;
   tau_eta = 0;
   tau_phi = 0;
   tau_origin = 0;
   jet_AntiKt4TruthJets_E = 0;
   jet_AntiKt4TruthJets_pt = 0;
   jet_AntiKt4TruthJets_m = 0;
   jet_AntiKt4TruthJets_eta = 0;
   jet_AntiKt4TruthJets_phi = 0;
   jet_AntiKt4TruthJets_flavor = 0;
   jet_AntiKt10TruthTrimmed_E = 0;
   jet_AntiKt10TruthTrimmed_pt = 0;
   jet_AntiKt10TruthTrimmed_m = 0;
   jet_AntiKt10TruthTrimmed_eta = 0;
   jet_AntiKt10TruthTrimmed_phi = 0;
   jet_AntiKt10TruthTrimmed_flavor = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);
   xAODTruth=-1;
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("susy_part_id1", &susy_part_id1, &b_susy_part_id1);
   fChain->SetBranchAddress("susy_part_id2", &susy_part_id2, &b_susy_part_id2);
   fChain->SetBranchAddress("mcevt_pdf_id1", &mcevt_pdf_id1, &b_mcevt_pdf_id1);
   fChain->SetBranchAddress("mcevt_pdf_id2", &mcevt_pdf_id2, &b_mcevt_pdf_id2);
   fChain->SetBranchAddress("mcevt_pdf_x1", &mcevt_pdf_x1, &b_mcevt_pdf_x1);
   fChain->SetBranchAddress("mcevt_pdf_x2", &mcevt_pdf_x2, &b_mcevt_pdf_x2);
   fChain->SetBranchAddress("mcevt_pdf1", &mcevt_pdf1, &b_mcevt_pdf1);
   fChain->SetBranchAddress("mcevt_pdf2", &mcevt_pdf2, &b_mcevt_pdf2);
   fChain->SetBranchAddress("mcevt_pdf_scale", &mcevt_pdf_scale, &b_mcevt_pdf_scale);
   fChain->SetBranchAddress("MET_Truth_NonInt_etx", &MET_Truth_NonInt_etx, &b_MET_Truth_NonInt_etx);
   fChain->SetBranchAddress("MET_Truth_NonInt_ety", &MET_Truth_NonInt_ety, &b_MET_Truth_NonInt_ety);
   fChain->SetBranchAddress("llv_n", &llv_n, &b_llv_n);
   fChain->SetBranchAddress("llv_id", &llv_id, &b_llv_id);
   fChain->SetBranchAddress("llv_pvx", &llv_pvx, &b_llv_pvx);
   fChain->SetBranchAddress("llv_pvy", &llv_pvy, &b_llv_pvy);
   fChain->SetBranchAddress("llv_pvz", &llv_pvz, &b_llv_pvz);
   fChain->SetBranchAddress("llv_evx", &llv_evx, &b_llv_evx);
   fChain->SetBranchAddress("llv_evy", &llv_evy, &b_llv_evy);
   fChain->SetBranchAddress("llv_evz", &llv_evz, &b_llv_evz);
   fChain->SetBranchAddress("ph_n", &ph_n, &b_ph_n);
   fChain->SetBranchAddress("ph_pt", &ph_pt, &b_ph_pt);
   fChain->SetBranchAddress("ph_eta", &ph_eta, &b_ph_eta);
   fChain->SetBranchAddress("ph_phi", &ph_phi, &b_ph_phi);
   fChain->SetBranchAddress("el_n", &el_n, &b_el_n);
   fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
   fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_origin", &el_origin, &b_el_origin);
   fChain->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
   fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_origin", &mu_origin, &b_mu_origin);
   fChain->SetBranchAddress("tau_n", &tau_n, &b_tau_n);
   fChain->SetBranchAddress("tau_charge", &tau_charge, &b_tau_charge);
   fChain->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_eta", &tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("tau_phi", &tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("tau_origin", &tau_origin, &b_tau_origin);
   fChain->SetBranchAddress("jet_AntiKt4TruthJets_n", &jet_AntiKt4TruthJets_n, &b_jet_AntiKt4TruthJets_n);
   fChain->SetBranchAddress("jet_AntiKt4TruthJets_E", &jet_AntiKt4TruthJets_E, &b_jet_AntiKt4TruthJets_E);
   fChain->SetBranchAddress("jet_AntiKt4TruthJets_pt", &jet_AntiKt4TruthJets_pt, &b_jet_AntiKt4TruthJets_pt);
   fChain->SetBranchAddress("jet_AntiKt4TruthJets_m", &jet_AntiKt4TruthJets_m, &b_jet_AntiKt4TruthJets_m);
   fChain->SetBranchAddress("jet_AntiKt4TruthJets_eta", &jet_AntiKt4TruthJets_eta, &b_jet_AntiKt4TruthJets_eta);
   fChain->SetBranchAddress("jet_AntiKt4TruthJets_phi", &jet_AntiKt4TruthJets_phi, &b_jet_AntiKt4TruthJets_phi);
   fChain->SetBranchAddress("jet_AntiKt4TruthJets_flavor", &jet_AntiKt4TruthJets_flavor, &b_jet_AntiKt4TruthJets_flavor);
   fChain->SetBranchAddress("jet_AntiKt10TruthTrimmed_n", &jet_AntiKt10TruthTrimmed_n, &b_jet_AntiKt10TruthTrimmed_n);
   fChain->SetBranchAddress("jet_AntiKt10TruthTrimmed_E", &jet_AntiKt10TruthTrimmed_E, &b_jet_AntiKt10TruthTrimmed_E);
   fChain->SetBranchAddress("jet_AntiKt10TruthTrimmed_pt", &jet_AntiKt10TruthTrimmed_pt, &b_jet_AntiKt10TruthTrimmed_pt);
   fChain->SetBranchAddress("jet_AntiKt10TruthTrimmed_m", &jet_AntiKt10TruthTrimmed_m, &b_jet_AntiKt10TruthTrimmed_m);
   fChain->SetBranchAddress("jet_AntiKt10TruthTrimmed_eta", &jet_AntiKt10TruthTrimmed_eta, &b_jet_AntiKt10TruthTrimmed_eta);
   fChain->SetBranchAddress("jet_AntiKt10TruthTrimmed_phi", &jet_AntiKt10TruthTrimmed_phi, &b_jet_AntiKt10TruthTrimmed_phi);
   fChain->SetBranchAddress("jet_AntiKt10TruthTrimmed_flavor", &jet_AntiKt10TruthTrimmed_flavor, &b_jet_AntiKt10TruthTrimmed_flavor);
}

Bool_t OldSlimReaderSelector::Notify()
{
   return kTRUE;
}

#endif // #ifdef OldSlimReader_cxx

