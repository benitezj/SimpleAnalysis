#define OldSlimReader_cxx

#include <cmath>

#include "SimpleAnalysis/OldSlimReader.h"
#include <TH2.h>
#include <TStyle.h>

#include "SimpleAnalysis/TruthEvent.h"
#include "SimpleAnalysis/OutputHandler.h"

void OldSlimReaderSelector::Begin(TTree * /*tree*/)
{
  //  TString option = GetOption();
 
}

void OldSlimReaderSelector::SlaveBegin(TTree * /*tree*/)
{
  //TString option = GetOption();

}

Bool_t OldSlimReaderSelector::Process(Long64_t entry)
{
  if ((entry%10000)==0) std::cout<<"At event "<<entry<<std::endl;
  GetEntry(entry);
  if (xAODTruth==-1) {
    xAODTruth=0;
    if (std::string(fChain->GetTree()->GetTitle())=="OldSlimmed xTruth1 tree") {
      std::cout<<"Reading xAOD Truth"<<std::endl;
      xAODTruth=1;
    }
  }
  
  TruthEvent* event=new TruthEvent(0,MET_Truth_NonInt_etx/1000.,MET_Truth_NonInt_ety/1000.); //FIXME: need sumet as well
  
  TLorentzVector tlv(0.,0.,0.,0.);
  for(int idx=0; idx<el_n; idx++) {
    int pdgID=fabs(el_origin->at(idx));
    int iso=(pdgID==25 || pdgID==24 || pdgID==23 || pdgID==15 || (pdgID>1000000 && pdgID<3000000));
    if (xAODTruth) iso = (pdgID==2); //origin is actually from MCTruthClassifier
    tlv.SetPtEtaPhiM(el_pt->at(idx)/1000.,el_eta->at(idx),el_phi->at(idx),0.510998910/1000.);
    event->addElectron(tlv,el_charge->at(idx),iso?EIsoGood:0,idx);
  }
  for(int idx=0; idx<mu_n; idx++) {
    int pdgID=fabs(mu_origin->at(idx));
    int iso=(pdgID==25 || pdgID==24 || pdgID==23 || pdgID==15 || (pdgID>1000000 && pdgID<3000000));
    if (xAODTruth) iso = (pdgID==6); //origin is actually from MCTruthClassifier
    tlv.SetPtEtaPhiM(mu_pt->at(idx)/1000.,mu_eta->at(idx),mu_phi->at(idx),105.6583715/1000.);
    event->addMuon(tlv,mu_charge->at(idx),iso?MuIsoGood:0,idx);
  }
  for(int idx=0; idx<tau_n; idx++) {
    int pdgID=fabs(tau_origin->at(idx));
    int iso=(pdgID==25 || pdgID==24 || pdgID==23 || pdgID==15 || (pdgID>1000000 && pdgID<3000000));
    if (xAODTruth) iso = (pdgID==10); //origin is actually from MCTruthClassifier

    tlv.SetPtEtaPhiM(tau_pt->at(idx)/1000.,tau_eta->at(idx),tau_phi->at(idx),1776.82/1000.);
    event->addTau(tlv,tau_charge->at(idx),iso?TauIsoGood:0,idx);
  }
  for(int idx=0; idx<jet_AntiKt4TruthJets_n; idx++) {
    tlv.SetPtEtaPhiM(jet_AntiKt4TruthJets_pt->at(idx)/1000.,jet_AntiKt4TruthJets_eta->at(idx),jet_AntiKt4TruthJets_phi->at(idx),jet_AntiKt4TruthJets_m->at(idx)/1000.);
    event->addJet(tlv,(abs(jet_AntiKt4TruthJets_flavor->at(idx))==5)?GoodBJet:GoodJet,idx);
  }
  //TODO: Add photons and fat jets

  double weight=1;
  _runner->processEvent(event,weight,entry);

  delete event;
  return kTRUE;
}

void OldSlimReaderSelector::SlaveTerminate()
{
}

void OldSlimReaderSelector::Terminate()
{
}

void OldSlimReader::processFilesInternal(const std::vector<std::string>& inputFileNames) {
  TChain *tree=new TChain("truth");
  for(const auto& fileName : inputFileNames) {
    if (tree->Add(fileName.c_str())<1) 
      throw std::runtime_error("Could not find truth tree in input file !");
  }
  
  OldSlimReaderSelector* selector = new OldSlimReaderSelector(_analysisRunner);
  tree->Process(selector);
  delete selector;
  delete tree;
}
