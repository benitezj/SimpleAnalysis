#define SlimReader_cxx

#include <cmath>

#include "SimpleAnalysis/SlimReader.h"
#include <TH2.h>
#include <TStyle.h>

#include "SimpleAnalysis/TruthEvent.h"
#include "SimpleAnalysis/OutputHandler.h"

void SlimReaderSelector::Begin(TTree * /*tree*/)
{
  //  TString option = GetOption();
 
}

void SlimReaderSelector::SlaveBegin(TTree * /*tree*/)
{
  //TString option = GetOption();

}

Bool_t SlimReaderSelector::Process(Long64_t entry)
{
  if ((entry%10000)==0) std::cout<<"At event "<<entry<<std::endl;
  GetEntry(entry);

 
  TruthEvent* event=new TruthEvent(sumet,met_pt*cos(met_phi),met_pt*sin(met_phi));
  event->setChannelInfo(mcChannel,susyChannel);
  event->setGenMET(genMET);
  event->setGenHT(genHT);
  event->setPDFInfo(pdf_id1,pdf_x1,pdf_pdf1,pdf_id2,pdf_x2,pdf_pdf2,pdf_scale);

  TLorentzVector tlv(0.,0.,0.,0.);
  for(unsigned int idx=0; idx<obj_pt[0]->size(); idx++) {
    tlv.SetPtEtaPhiM(obj_pt[0]->at(idx),
		     obj_eta[0]->at(idx),
		     obj_phi[0]->at(idx),0.510998910/1000.);
    event->addElectron(tlv,obj_charge[0]->at(idx),obj_id[0]->at(idx),idx);
  }
  for(unsigned int idx=0; idx<obj_pt[1]->size(); idx++) {
    tlv.SetPtEtaPhiM(obj_pt[1]->at(idx),
		     obj_eta[1]->at(idx),
		     obj_phi[1]->at(idx),105.6583715/1000.);
    event->addMuon(tlv,obj_charge[1]->at(idx),obj_id[1]->at(idx),idx);
  }
  for(unsigned int idx=0; idx<obj_pt[2]->size(); idx++) {
    tlv.SetPtEtaPhiM(obj_pt[2]->at(idx),
		     obj_eta[2]->at(idx),
		     obj_phi[2]->at(idx),1776.82/1000.);
    event->addTau(tlv,obj_charge[2]->at(idx),obj_id[2]->at(idx),idx);
  }
  for(unsigned int idx=0; idx<obj_pt[3]->size(); idx++) {
    tlv.SetPtEtaPhiM(obj_pt[3]->at(idx),
		     obj_eta[3]->at(idx),
		     obj_phi[3]->at(idx),0); 
    event->addPhoton(tlv,obj_id[3]->at(idx),idx);
  }
  for(unsigned int idx=0; idx<obj_pt[4]->size(); idx++) {
    tlv.SetPtEtaPhiM(obj_pt[4]->at(idx),
		     obj_eta[4]->at(idx),
		     obj_phi[4]->at(idx),jet_m->at(idx)); 
    event->addJet(tlv,obj_id[4]->at(idx),idx);
  }
  for(unsigned int idx=0; idx<obj_pt[5]->size(); idx++) {
    tlv.SetPtEtaPhiM(obj_pt[5]->at(idx),
		     obj_eta[5]->at(idx),
		     obj_phi[5]->at(idx),fatjet_m->at(idx)); 
    event->addFatJet(tlv,obj_id[5]->at(idx),idx);
  }

  event->setMCWeights(*mcWeights);

  _runner->processEvent(event,EventNumber);

  delete event;
  return kTRUE;
}

void SlimReaderSelector::SlaveTerminate()
{
}

void SlimReaderSelector::Terminate()
{
}

void SlimReader::processFilesInternal(const std::vector<std::string>& inputFileNames) {
  TChain *tree=new TChain("ntuple");
  for(const auto& fileName : inputFileNames) {
    if (tree->Add(fileName.c_str())<1) 
      throw std::runtime_error("Could not find slimmed tree in input file !");
  }
  
  SlimReaderSelector* selector = new SlimReaderSelector(_analysisRunner);
  tree->Process(selector);
  delete selector;
  delete tree;
}
