#define D3PDReader_cxx

#include <cmath>

#include "SimpleAnalysis/D3PDReader.h"
#include <TH2.h>
#include <TStyle.h>

#include "SimpleAnalysis/TruthEvent.h"
#include "SimpleAnalysis/OutputHandler.h"
#include "SUSYTools/SUSYCrossSection.h"

void D3PDReaderSelector::Begin(TTree * /*tree*/)
{
  //  TString option = GetOption();
 
}

void D3PDReaderSelector::SlaveBegin(TTree * /*tree*/)
{
  //TString option = GetOption();

}
int D3PDReaderSelector::parentID(vector<int> & parents,int id)
{
  int parent=0;
  if (parents.size()>0) {
    parent=std::abs(mc_pdgId->at(parents[0]));
    if (parent==id) {
      std::cout<<"found duplicate"<<std::endl;
    }
  }
  return parent;
}

bool D3PDReaderSelector::isIso(vector<int> & parents,int id) {
  int pdgID=parentID(parents,id);
  if (pdgID==25 || pdgID==24 || pdgID==23 || pdgID==15 || (pdgID>1000000 && pdgID<3000000))
    return true;
  return false;
}

Bool_t D3PDReaderSelector::Process(Long64_t entry)
{
  if ((entry%10000)==0) std::cout<<"At event "<<entry<<std::endl;
  GetEntry(entry);

  int mcChannel = mc_channel_number;

  int susy_part_id1 = 0;
  int susy_part_id2 = 0;
  int susy_process  = 0;
  int firstsp=-1;
  int secondsp=-1;
  for(int ii=0;ii<mc_n;ii++) {
    int pdg=std::abs(mc_pdgId->at(ii));
    if (pdg>1000000) {
      int parpdg=std::abs(parentID(mc_parent_index->at(ii),pdg));
      if (parpdg<1000000) {
        if (firstsp<0) {
          firstsp=ii;
        } else {
          secondsp=ii;
          break;
        }
      }
      
    }
  }
  if (secondsp<0) {
    std::cout<<"WARNING did not find susy particles"<<std::endl; 
  } else {
    if (mc_child_index->at(firstsp).size()==1) {
      if (mc_pdgId->at(mc_child_index->at(firstsp)[0]) != mc_pdgId->at(firstsp)) {
        std::cout<<"propagator:"<<(mc_pdgId->at(mc_child_index->at(firstsp)[0]),mc_pdgId->at(firstsp))<<std::endl;
        firstsp=mc_child_index->at(firstsp)[0];
      }
    }
    susy_part_id1 = mc_pdgId->at(firstsp);
    susy_part_id2 = mc_pdgId->at(secondsp);
  }
  if ((abs(susy_part_id1)>1000000) && (abs(susy_part_id1)>1000000)) //only consider BSM particles
    susy_process = SUSY::finalState(susy_part_id1,susy_part_id2);

  float sumet=0;
  //FIXME: approximation since no sumet was stored
  for( const auto jetpt : *jet_AntiKt4TruthJets_pt) sumet += jetpt; 
  float genHT  = sumet/1000.;
  float genMET = sqrt(MET_Truth_NonInt_etx*MET_Truth_NonInt_etx+MET_Truth_NonInt_ety*MET_Truth_NonInt_ety)/1000.;
  for( const auto mupt : *mu_pt) sumet += mupt;

  TruthEvent* event=new TruthEvent(sumet/1000.,MET_Truth_NonInt_etx/1000.,MET_Truth_NonInt_ety/1000.);
  event->setChannelInfo(mcChannel,susy_process);
  event->setGenMET(genMET);
  event->setGenHT(genHT);
  event->setPDFInfo(mcevt_pdf_id1->at(0),mcevt_pdf_x1->at(0),mcevt_pdf1->at(0),
		    mcevt_pdf_id2->at(0),mcevt_pdf_x2->at(0),mcevt_pdf2->at(0),
		    mcevt_pdf_scale->at(0));

  TLorentzVector tlv(0.,0.,0.,0.);
  for(int idx=0; idx<el_n; idx++) {
    tlv.SetPtEtaPhiM(el_pt->at(idx)/1000, el_eta->at(idx), el_phi->at(idx), el_m->at(idx)/1000.);
    event->addElectron(tlv,el_charge->at(idx),isIso(el_parent_index->at(idx),11)?EIsoGood:0,idx);
  }
  for(int idx=0; idx<mu_n; idx++) {
    tlv.SetPtEtaPhiM(mu_pt->at(idx)/1000, mu_eta->at(idx), mu_phi->at(idx), mu_m->at(idx)/1000.);
    event->addMuon(tlv,mu_charge->at(idx),isIso(mu_parent_index->at(idx),13)?MuIsoGood:0,idx);
  }
  int numTau=0;
  for(int idx=0; idx<tau_n; idx++) {
    bool leptonic=false;
    int numCharged=0;
    for(unsigned int jj=0;jj<tau_decay_index->at(idx).size();jj++) {
        int tidx=tau_decay_index->at(idx)[jj];
        int pdg=std::abs(mc_pdgId->at(tidx));
        if (pdg>10 && pdg<15) leptonic=true;
        if (pdg==211||pdg==321) numCharged++;
        else if (pdg>100&&pdg!=310&&pdg!=130) 
          std::cout<<tidx<<" found : "<<pdg<<" "<<mc_pt->at(tidx)<<std::endl;
    }
    if (leptonic and numCharged!=0) {
      std::cout<<"ERROR: "<<entry<<" "<<numCharged<<std::endl;
    }
    if (not leptonic) {
      tlv.SetPtEtaPhiM(tau_pt->at(idx)/1000, tau_eta->at(idx), tau_phi->at(idx),tau_m->at(idx)/1000.);
      int tauId=isIso(tau_parent_index->at(idx),15)?TauIsoGood:0;
      if (numCharged!=1) tauId|=TauThreeProng;
      else tauId|=TauOneProng;
      event->addTau(tlv,tau_charge->at(idx),tauId,numTau++);
    }
  }
  for(int idx=0; idx<ph_n; idx++) {
    tlv.SetPtEtaPhiM(ph_pt->at(idx)/1000., ph_eta->at(idx), ph_phi->at(idx),0); 
    event->addPhoton(tlv,0,idx);
  }
  for(int idx=0; idx<jet_AntiKt4TruthJets_n; idx++) {
    tlv.SetPtEtaPhiM(jet_AntiKt4TruthJets_pt->at(idx)/1000., jet_AntiKt4TruthJets_eta->at(idx),
		     jet_AntiKt4TruthJets_phi->at(idx),jet_AntiKt4TruthJets_m->at(idx)/1000.);
    int flavor = abs(jet_AntiKt4TruthJets_flavor_partonFlavor->at(idx));
    int id=GoodJet;
    if (flavor==4)       id |= TrueCJet;
    else if (flavor==5)  id |= TrueBJet|GoodBJet;
    else if (flavor==15) id |= TrueTau;
    else                 id |= TrueLightJet;
    event->addJet(tlv,id,idx);
  }
  for(int idx=0; idx<jet_AntiKt10TruthTrimmedPtFrac5SmallR30_n; idx++) {
    tlv.SetPtEtaPhiM(jet_AntiKt10TruthTrimmedPtFrac5SmallR30_pt->at(idx)/1000., 
		     jet_AntiKt10TruthTrimmedPtFrac5SmallR30_eta->at(idx),
		     jet_AntiKt10TruthTrimmedPtFrac5SmallR30_phi->at(idx),
		     jet_AntiKt10TruthTrimmedPtFrac5SmallR30_m->at(idx)/1000.);
    int flavor = abs(jet_AntiKt10TruthTrimmedPtFrac5SmallR30_flavor_partonFlavor->at(idx));
    int id=GoodJet;
    if (flavor==4)       id |= TrueCJet;
    else if (flavor==5)  id |= TrueBJet|GoodBJet;
    else if (flavor==15) id |= TrueTau;
    else                 id |= TrueLightJet;
    event->addFatJet(tlv,id,idx);
  }
  std::vector<float> weights;
  weights.push_back(mc_event_weight);
  event->setMCWeights(weights);

  _runner->processEvent(event,EventNumber);

  delete event;
  return kTRUE;
}

void D3PDReaderSelector::SlaveTerminate()
{
}

void D3PDReaderSelector::Terminate()
{
}

void D3PDReader::processFilesInternal(const std::vector<std::string>& inputFileNames) {
  TChain *tree=new TChain("truth");
  for(const auto& fileName : inputFileNames) {
    if (tree->Add(fileName.c_str())<1) 
      throw std::runtime_error("Could not find truth tree in input file !");
  }
  
  D3PDReaderSelector* selector = new D3PDReaderSelector(_analysisRunner);
  tree->Process(selector);
  delete selector;
  delete tree;
}
