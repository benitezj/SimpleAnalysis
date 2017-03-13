#include <cmath>
#include <iostream>
#include <algorithm>

#include "SimpleAnalysis/xAODTruthReader.h"
#include <TH2.h>
#include <TStyle.h>
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODJet/JetContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "SUSYTools/SUSYCrossSection.h"
#include "xAODCore/AuxContainerBase.h"

#include "SimpleAnalysis/TruthEvent.h"
#include "SimpleAnalysis/OutputHandler.h"

using std::vector;

#define WARN_ONCE(warning)          \
  do {                              \
    static bool first=true;         \
    if (first) std::cout<<warning<<std::endl; \
    first=false;                              \
  } while(0)


xAODTruthReader::xAODTruthReader(std::vector<AnalysisClass*>& analysisList) : Reader(analysisList)
{
  _event=new xAOD::TEvent(xAOD::TEvent::kClassAccess);
  _susytools= new ST::SUSYObjDef_xAOD("mySUSYTools");
  _mctool=new MCTruthClassifier("myTruthFinder");
}

int xAODTruthReader::getTruthOrigin(const xAOD::TruthParticle *part) {
  if (part->isAvailable<unsigned int>("classifierParticleOrigin")) {
    return part->auxdata<unsigned int>("classifierParticleOrigin");
  }
  const ElementLink < xAOD::TruthParticleContainer > origPart = part->auxdata< ElementLink< xAOD::TruthParticleContainer > >("originalTruthParticle" );
  if (origPart.isValid()) {
    const auto result = _mctool->particleTruthClassifier(*origPart);
    return result.second;
  }
  return 0;
}

int xAODTruthReader::getTruthType(const xAOD::TruthParticle *part) {
  if (part->isAvailable<unsigned int>("classifierParticleType")) {
    return part->auxdata<unsigned int>("classifierParticleType");
  }
  const ElementLink < xAOD::TruthParticleContainer > origPart = part->auxdata< ElementLink< xAOD::TruthParticleContainer > >("originalTruthParticle" );
  if (origPart.isValid()) {
    const auto result = _mctool->particleTruthClassifier(*origPart);
    return result.first;
  }
  return 0;
}

xAOD::TruthParticleContainer*
xAODTruthReader::findTruthParticles(xAOD::TStore *store,
				  const xAOD::TruthParticleContainer* truthparticles,
				    std::vector<int> pdgIds, int status) {
  xAOD::TruthParticleContainer* truth = new xAOD::TruthParticleContainer;
  xAOD::AuxContainerBase* truthAux = new xAOD::AuxContainerBase();
  truth->setStore( truthAux );
  SG::AuxElement::Decorator< ElementLink<xAOD::TruthParticleContainer> > linkDecorator("originalTruthParticle");
  int idx=0;
  for ( xAOD::TruthParticleContainer::const_iterator it = truthparticles->begin();
	it != truthparticles->end(); ++it ) {

    if ( std::find(pdgIds.begin(), pdgIds.end(), abs((*it)->pdgId()) ) != pdgIds.end()
	 && (*it)->status()==status && (*it)->barcode()<200000) {
      const auto result = _mctool->particleTruthClassifier(*it);
      MCTruthPartClassifier::ParticleOutCome outcome=_mctool->getParticleOutCome();
      if (outcome==MCTruthPartClassifier::DecaytoElectron ||
	  outcome==MCTruthPartClassifier::DecaytoMuon) continue;
      xAOD::TruthParticle* part = new xAOD::TruthParticle();
      truth->push_back(part);
      *part = **it;
      part->auxdecor<unsigned int>("classifierParticleType") = result.first;
      ElementLink<xAOD::TruthParticleContainer> eltp(*truthparticles,idx);
      linkDecorator(*part) = eltp;
      if (status==2) { //tau
	int numPart=0;
	if (outcome==MCTruthPartClassifier::OneProng) numPart=1;
	if (outcome==MCTruthPartClassifier::ThreeProng) numPart=3;
	if (outcome==MCTruthPartClassifier::FiveProng) numPart=5;
	part->auxdecor<char>("IsHadronicTau") = (char)1;
	part->auxdecor<size_t>("numCharged") = numPart;
      }
    }
    idx++;
  }

  std::string pid;
  for(auto ipdg : pdgIds)  pid += std::to_string(ipdg);

  if (!store->record( truth, (std::string("Good")+pid).c_str()).isSuccess())
    throw std::runtime_error("Could not record truth particles");
  if (!store->record( truthAux, (std::string("Good")+pid+"Aux.").c_str()).isSuccess())
    throw std::runtime_error("Could not record truth particles Aux");
  return truth;
}

static SG::AuxElement::Accessor<float> acc_filtHT("GenFiltHT");
static SG::AuxElement::Accessor<float> acc_filtMET("GenFiltMET");
static SG::AuxElement::ConstAccessor<int> acc_HadronConeExclTruthLabelID("HadronConeExclTruthLabelID");

bool xAODTruthReader::processEvent(xAOD::TEvent *xaodEvent,xAOD::TStore *store) {

  const xAOD::EventInfo* eventInfo = 0;
  if ( !xaodEvent->retrieve( eventInfo, "EventInfo").isSuccess() ) {
    throw std::runtime_error("Cannot read EventInfo");
  }

  int eventNumber = eventInfo->eventNumber();
  int mcChannel   = eventInfo->mcChannelNumber();
  if (mcChannel==0) mcChannel = eventInfo->runNumber();
  int susy_part_id1 = 0;
  int susy_part_id2 = 0;
  int susy_process  = 0;

  const xAOD::TruthParticleContainer* truthparticles = 0;
  if ( xaodEvent->contains<xAOD::TruthParticleContainer>("TruthParticles")) {
    if ( !xaodEvent->retrieve( truthparticles, "TruthParticles").isSuccess() ) {
      throw std::runtime_error("Could not retrieve truth particles with key TruthParticles");
    }
  } else {
    if ( !xaodEvent->retrieve( truthparticles, "TruthBSM").isSuccess() ) {
      throw std::runtime_error("Could not retrieve truth particles with key TruthBSM");
    }
  }
  _susytools->FindSusyHardProc(truthparticles,susy_part_id1,susy_part_id2);
  if (susy_part_id2==0 && truthparticles->size()>1) susy_part_id2=truthparticles->at(1)->pdgId();
  if (susy_part_id1==0 && truthparticles->size()) susy_part_id1=truthparticles->at(0)->pdgId();
  if ((abs(susy_part_id1)>1000000) && (abs(susy_part_id1)>1000000)) //only consider BSM particles
    susy_process = SUSY::finalState(susy_part_id1,susy_part_id2);

  const xAOD::MissingETContainer* metCont = 0;
  if ( !xaodEvent->retrieve(metCont, "MET_Truth").isSuccess() ){
    throw std::runtime_error("Could not retrieve truth met with key MET_Truth");
  }
  const xAOD::MissingET* met = (*metCont)["NonInt"];

  TruthEvent* event=new TruthEvent(met->sumet()/1000.,met->mpx()/1000.,met->mpy()/1000.);
  event->setChannelInfo(mcChannel,susy_process);

  TLorentzVector tlv(0.,0.,0.,0.);

  int idx=0;
  const xAOD::TruthParticleContainer* truthelectrons = 0;
  if ( xaodEvent->contains<xAOD::TruthParticleContainer>("TruthElectrons")) {
    if ( !xaodEvent->retrieve( truthelectrons, "TruthElectrons").isSuccess() ) {
      throw std::runtime_error("Could not retrieve truth particles with key TruthElectrons");
    }
  } else {
    truthelectrons=findTruthParticles(store,truthparticles,{11});
  }
  for ( xAOD::TruthParticleContainer::const_iterator it = truthelectrons->begin();
	it != truthelectrons->end(); ++it ){
    const auto electron = *it;
    int iso = getTruthType(electron)==MCTruthPartClassifier::IsoElectron;
    tlv.SetPtEtaPhiM(electron->pt()/1000.,electron->eta(),electron->phi(),electron->m()/1000.);
    event->addElectron(tlv,electron->charge(),iso?EIsoGood:0,idx++);
  }

  idx=0;
  const xAOD::TruthParticleContainer* truthmuons = 0;
  if ( xaodEvent->contains<xAOD::TruthParticleContainer>("TruthMuons")) {
    if ( !xaodEvent->retrieve( truthmuons, "TruthMuons").isSuccess() ) {
      throw std::runtime_error("Could not retrieve truth particles with key TruthMuons");
    }
  } else {
    truthmuons=findTruthParticles(store,truthparticles,{13});
  }
  for ( xAOD::TruthParticleContainer::const_iterator it = truthmuons->begin();
	it != truthmuons->end(); ++it ){
    const auto muon = *it;
    int iso = getTruthType(muon)==MCTruthPartClassifier::IsoMuon;
    tlv.SetPtEtaPhiM(muon->pt()/1000.,muon->eta(),muon->phi(),muon->m()/1000.);
    event->addMuon(tlv,muon->charge(),iso?MuIsoGood:0,idx++);
  }

  idx=0;
  const xAOD::TruthParticleContainer* truthtaus = 0;
  if ( xaodEvent->contains<xAOD::TruthParticleContainer>("TruthTaus")) {
    if ( !xaodEvent->retrieve( truthtaus, "TruthTaus").isSuccess() ) {
      throw std::runtime_error("Could not retrieve truth particles with key TruthTaus");
    }
  } else {
    truthtaus=findTruthParticles(store,truthparticles,{15},2);
  }
  for ( xAOD::TruthParticleContainer::const_iterator it = truthtaus->begin();
	it != truthtaus->end(); ++it ){
    const auto tau = *it;
    if (tau->auxdata<char>("IsHadronicTau")) {
      int iso = getTruthType(tau)==MCTruthPartClassifier::IsoTau;
      tlv.SetPtEtaPhiM(tau->pt()/1000.,tau->eta(),tau->phi(),tau->m()/1000.);
      int tauId=iso?TauIsoGood:0;
      if (tau->auxdata<unsigned long>("numCharged")==3) tauId|=TauThreeProng;
      else tauId|=TauOneProng;
      event->addTau(tlv,tau->charge(),tauId,idx++);
    }
  }

  idx=0;
  const xAOD::TruthParticleContainer* truthphotons = 0;
  std::string photonName="TruthPhotons";
  if ( !xaodEvent->contains<xAOD::TruthParticleContainer>(photonName)) photonName="Truth3Photons";

  if ( xaodEvent->contains<xAOD::TruthParticleContainer>(photonName)) {
    if ( !xaodEvent->retrieve( truthphotons, photonName).isSuccess() ) {
      throw std::runtime_error("Could not retrieve truth particles with key TruthPhotons");
    }
  } else {
    truthphotons=findTruthParticles(store,truthparticles,{22});
  }

  for ( xAOD::TruthParticleContainer::const_iterator it = truthphotons->begin();
	it != truthphotons->end(); ++it ){
    const auto photon = *it;
    int iso = getTruthType(photon)==MCTruthPartClassifier::IsoPhoton;
    tlv.SetPtEtaPhiM(photon->pt()/1000.,photon->eta(),photon->phi(),0);
    event->addPhoton(tlv,iso?PhotonIsoGood:0,idx++);
  }


  //Generator Filter HT (e.g. for ttbar/singleTop samples)
  float gen_ht=0.;
  if ( acc_filtHT.isAvailable(*(eventInfo)) ){
    gen_ht = eventInfo->auxdata<float>("GenFiltHT");
  }
  else{
    WARN_ONCE("Warning : No GenFiltHT decoration available. Setting HT to 0 for now...");
  }
  event->setGenHT( gen_ht/1000. );


  //Generator Filter MET (e.g. for ttbar/singleTop samples)
  float gen_met=0.;
  if ( acc_filtMET.isAvailable(*(eventInfo)) ){
    gen_met = eventInfo->auxdata<float>("GenFiltMET");
  }
  else{ //recompute from particle containers!
    idx=0;
    const xAOD::TruthParticleContainer* truthneutrinos = 0;
    std::string neutrinoName="TruthNeutrinos";
    if ( !xaodEvent->contains<xAOD::TruthParticleContainer>(neutrinoName)) neutrinoName="TruthNeutrinos";

    if ( xaodEvent->contains<xAOD::TruthParticleContainer>(neutrinoName)) {
      if ( !xaodEvent->retrieve( truthneutrinos, neutrinoName).isSuccess() ) {
	throw std::runtime_error("Could not retrieve truth particles with key TruthNeutrinos");
      }
    } else {
      truthneutrinos=findTruthParticles(store,truthparticles,{12,14,16});
    }


    tlv.SetPtEtaPhiM(0.,0.,0.,0.);
    for ( xAOD::TruthParticleContainer::const_iterator it = truthneutrinos->begin();
	  it != truthneutrinos->end(); ++it ){
      const auto nu = *it;
      int iPartOrig = getTruthOrigin(nu);

      switch (iPartOrig) {
      case MCTruthPartClassifier::PhotonConv:
      case MCTruthPartClassifier::DalitzDec:
      case MCTruthPartClassifier::ElMagProc:
      case MCTruthPartClassifier::Mu:
      case MCTruthPartClassifier::TauLep:
      case MCTruthPartClassifier::LightMeson:
      case MCTruthPartClassifier::StrangeMeson:
      case MCTruthPartClassifier::CharmedMeson:
      case MCTruthPartClassifier::BottomMeson:
      case MCTruthPartClassifier::CCbarMeson:
      case MCTruthPartClassifier::JPsi:
      case MCTruthPartClassifier::BBbarMeson:
      case MCTruthPartClassifier::LightBaryon:
      case MCTruthPartClassifier::StrangeBaryon:
      case MCTruthPartClassifier::CharmedBaryon:
      case MCTruthPartClassifier::BottomBaryon:
      case MCTruthPartClassifier::PionDecay:
      case MCTruthPartClassifier::KaonDecay:

      case MCTruthPartClassifier::NonDefined:
	continue;
      default:
	break;
      }
      tlv += nu->p4();

    }
    gen_met = tlv.Pt();
  }
  event->setGenMET( gen_met/1000. );



  idx=0;
  const xAOD::JetContainer* truthjets = 0;
  if ( !xaodEvent->retrieve( truthjets, "AntiKt4TruthJets").isSuccess() ) {
    throw std::runtime_error("Could not retrieve truth particles with key TruthJets");
  }
  for ( xAOD::JetContainer::const_iterator it = truthjets->begin();
	it != truthjets->end(); ++it ){
    const auto jet = *it;
    tlv.SetPtEtaPhiM(jet->pt()/1000.,jet->eta(),jet->phi(),jet->m()/1000.);
    int flavor=0;
    if (jet->isAvailable<int>("HadronConeExclusionID"))
      flavor=jet->auxdata<int>("HadronConeExclusionID");
    else if (jet->isAvailable<int>("ConeTruthLabelID"))
      flavor=jet->auxdata<int>("ConeTruthLabelID");
    else if (jet->isAvailable<int>("PartonTruthLabelID"))
      flavor=abs(jet->auxdata<int>("PartonTruthLabelID"));
    else if (jet->isAvailable<int>("GhostBHadronsFinalCount")) {
      if (jet->auxdata<int>("GhostBHadronsFinalCount")) {
	flavor=5;
      } else if (jet->auxdata<int>("GhostCHadronsFinalCount")) {
	flavor=4;
      } else flavor=1;
    }
    int id=(flavor==5)?GoodBJet:GoodJet;
    if (flavor==4)       id |= TrueCJet;
    else if (flavor==5)  id |= TrueBJet;
    else if (flavor==15) id |= TrueTau;
    else                 id |= TrueLightJet;
    int truth_id = (acc_HadronConeExclTruthLabelID.isAvailable(*(jet))) ? acc_HadronConeExclTruthLabelID(*(jet)) : 1;
    event->addJet(tlv,id,truth_id,idx);
    //event->addJet(tlv,id,flavor,idx);
  }

  idx=0;
  const xAOD::JetContainer* truthfatjets = 0;
  std::string fatjetName="AntiKt10TruthTrimmedPtFrac5SmallR20Jets";
  if ( !xaodEvent->contains<xAOD::JetContainer>(fatjetName)) fatjetName="TrimmedAntiKt10TruthJets";

  if ( xaodEvent->contains<xAOD::JetContainer>(fatjetName)) {
    if ( !xaodEvent->retrieve( truthfatjets, fatjetName).isSuccess() ) {
      throw std::runtime_error("Could not retrieve truth particles with key "+fatjetName);
    }
    for ( xAOD::JetContainer::const_iterator it = truthfatjets->begin();
  	it != truthfatjets->end(); ++it ){
      const auto jet = *it;
      tlv.SetPtEtaPhiM(jet->pt()/1000.,jet->eta(),jet->phi(),jet->m()/1000.);
      int flavor=0; //FIXME check if there are more recent labels for fat jets
      if (jet->isAvailable<int>("PartonTruthLabelID"))
        flavor=abs(jet->auxdata<int>("PartonTruthLabelID"));
      else {
	if (jet->isAvailable<int>("GhostBHadronsFinalCount")) {
	  if (jet->auxdata<int>("GhostBHadronsFinalCount")) {
	    flavor=5;
	  } else if (jet->auxdata<int>("GhostCHadronsFinalCount")) {
	    flavor=4;
	  } else flavor=1;
	}
      }
      event->addFatJet(tlv,(flavor==5)?GoodBJet:GoodJet,idx);
    }
  }

  //Get LHE3 weights
  const xAOD::TruthEventContainer* truthEvtCont;
  if( !xaodEvent->retrieve( truthEvtCont, "TruthEvents").isSuccess() )
    throw std::runtime_error("Could not retrieve truth event container with key TruthEvents");
  const xAOD::TruthEvent *truthevent = (*truthEvtCont)[0];
  const std::vector<float> weights  = truthevent->weights();

  event->setMCWeights(weights);

  _analysisRunner->processEvent(event,eventNumber);

  delete event;
  return true;
}

void xAODTruthReader::processFilesInternal(const std::vector<std::string>& inputFileNames) {
  xAOD::TStore transientStorage;
  transientStorage.setActive();
  TFile *inFile=0;
  for(const auto& inName : inputFileNames) {
    delete inFile;
    std::cout<<"Now reading: "<<inName<<std::endl;
    inFile=TFile::Open(inName.c_str());
    if ( ! _event->readFrom(inFile).isSuccess() ) {
      throw std::runtime_error("Could not connect TEvent to file !");
    }
    Long64_t numEntries=_event->getEntries();
    for(Long64_t index = 0; index<numEntries; index++) {
      Long64_t entry = _event->getEntry(index);
      if (entry<0) break;
      if (index%10000==0)
	std::cout<<"at: "<<index<<std::endl;
      processEvent(_event,&transientStorage);
      transientStorage.clear();
    }
  }
}

xAODTruthReader::~xAODTruthReader() {
  return;
  delete _mctool;
  delete _susytools;
  delete _event;
}
