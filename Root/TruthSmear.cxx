#include "SimpleAnalysis/TruthEvent.h"
#include "SimpleAnalysis/TruthSmear.h"
#include <RootCore/Packages.h>

#ifdef ROOTCORE_PACKAGE_UpgradePerformanceFunctions
#include <UpgradePerformanceFunctions/UpgradePerformanceFunctions.h>
#endif

TruthSmear::TruthSmear(std::vector<std::string>& 
#ifdef ROOTCORE_PACKAGE_UpgradePerformanceFunctions
options
#endif
) :
smearElectrons(true), smearMuons(true), smearTaus(true), smearPhotons(true), smearJets(true), smearMET(true), addPileupJets(false), useHGTD0(false),
useTrackConfirm(true), puEffScheme("PU"), puEff(0.02) {
#ifdef ROOTCORE_PACKAGE_UpgradePerformanceFunctions

  std::string mu="None";
  int seed=12345;
  for(const auto& option : options){
    if (option=="help") {
      std::cout<<"Options for smearing:"<<std::endl;
      std::cout<<" noElectrons   - do not smear electrons"<<std::endl;
      std::cout<<" noMuons       - do not smear muons"<<std::endl;
      std::cout<<" noTaus        - do not smear taus"<<std::endl;
      std::cout<<" noPhotons     - do not smear photons"<<std::endl;
      std::cout<<" noJets        - do not smear jets"<<std::endl;
      std::cout<<" noMET         - do not smear MET"<<std::endl;
      std::cout<<" addPileupJets - add pileup jets"<<std::endl;
      std::cout<<" mu=<value>    - choice pile-up level (required)"<<std::endl;
      std::cout<<"                 Only mu=200 is allowed for now"<<std::endl;
      std::cout<<" seed=<value>  - set seed value (default: "<<seed<<std::endl;
    }
    if (option=="noElectrons")    smearElectrons=false;
    if (option=="noMuon")         smearMuons=false;
    if (option=="noTaus")         smearTaus=false;
    if (option=="noPhotons")      smearPhotons=false;
    if (option=="noJets")         smearJets=false;
    if (option=="noMET")          smearMET=false;
    if (option=="addPileupJets")  addPileupJets=true;
    if (option=="useHGTD0")       useHGTD0=true;
    if (option=="noTrackConfirm") useTrackConfirm=false;
    if (option.find("PUeff=")==0) {
      puEff=std::stof(option.substr(6));
    }
    if (option.find("HSeff=")==0) {
      puEffScheme="HS";
      puEff=std::stof(option.substr(6));
    }
    if (option.find("mu=")==0) {
      mu=option.substr(3);
    }
    if (option.find("seed=")==0) {
      seed=stoi(option.substr(5));
    }
  }
  std::cout<<"Smearing with mu="<<mu<<" and seed="<<seed<<std::endl;

  if (mu!="200") throw std::runtime_error("Unsupported pile-up level. Only mu=200 currently supported");
  m_upgrade = new UpgradePerformanceFunctions();
  m_upgrade->setLayout(UpgradePerformanceFunctions::gold);
  m_upgrade->setAvgMu(stoi(mu));
  m_upgrade->setElectronWorkingPoint(UpgradePerformanceFunctions::looseElectron);
  m_upgrade->setElectronRandomSeed(seed);
  m_upgrade->setMuonWorkingPoint(UpgradePerformanceFunctions::tightMuon);
  m_upgrade->setPhotonWorkingPoint(UpgradePerformanceFunctions::tightPhoton);
  m_upgrade->setTauRandomSeed(seed);
  m_upgrade->setJetRandomSeed(seed);
  m_upgrade->setMETRandomSeed(seed);
  m_upgrade->loadMETHistograms("UpgradePerformanceFunctions/sumetPU_mu200_ttbar_gold.root");
  m_upgrade->setPileupRandomSeed(seed);
  if (useTrackConfirm) m_upgrade->setPileupUseTrackConf(true);
  else m_upgrade->setPileupUseTrackConf(false);
  m_upgrade->setPileupJetPtThresholdMeV(30000.);
  if     (puEffScheme == "PU") m_upgrade->setPileupEfficiencyScheme(UpgradePerformanceFunctions::PU);
  else if (puEffScheme == "HS") m_upgrade->setPileupEfficiencyScheme(UpgradePerformanceFunctions::HS);
  m_upgrade->setPileupEff(puEff);
  if (useHGTD0) m_upgrade->setUseHGTD0(true);
  m_upgrade->setPileupTemplatesPath("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/UpgradePerformanceFunctions/");
  m_upgrade->initPhotonFakeHistograms("UpgradePerformanceFunctions/PhotonFakes.root");
  m_upgrade->setFlavourTaggingCalibrationFilename("UpgradePerformanceFunctions/flavor_tags_v1.1.root");
  m_random.SetSeed(seed);

#else
  throw std::runtime_error("Compiled without smearing support - add UpgradePerformanceFunctions");
#endif
}

TruthEvent *TruthSmear::smearEvent(AnalysisEvent *
#ifdef ROOTCORE_PACKAGE_UpgradePerformanceFunctions
event
#endif
) {
#ifdef ROOTCORE_PACKAGE_UpgradePerformanceFunctions
  auto electrons  = event->getElectrons(1.,4.2); // Filter on pT, eta and "ID"
  auto muons      = event->getMuons(1.,4.2);
  auto taus       = event->getTaus(20.,2.5);
  auto photons    = event->getPhotons(5.,4.2);
  auto jets       = event->getJets(10.,5.2);
  auto fatjets    = event->getFatJets(10.,5.2);
  auto met        = event->getMET();
  auto sumet      = event->getSumET();
  
  double met_x = met.Px();
  double met_y = met.Py();
  if (smearMET) {
    auto smearedMET = m_upgrade->getMETSmeared(sumet*1000., met_x*1000., met_y*1000.);
    met_x = smearedMET.first/1000.;
    met_y = smearedMET.second/1000.;
  }
  TruthEvent* smeared=new TruthEvent(sumet,met_x,met_y);
  smeared->setChannelInfo(event->getMCNumber(),event->getSUSYChannel());
  smeared->setGenMET(event->getGenMET());
  smeared->setGenHT(event->getGenHT());
  smeared->setMCWeights(event->getMCWeights());

  int idx=0;
  for(const auto& electron : electrons) {
    if (smearElectrons){
      float eff = m_upgrade->getElectronEfficiency(electron.Pt()*1000., electron.Eta());
      if (m_random.Uniform(1.0)<eff) {
	float electron_e = m_upgrade->getElectronSmearedEnergy(electron.E()*1000., electron.Eta())/1000.; 
	TLorentzVector eLV;
	eLV.SetPtEtaPhiM(electron.Pt()*electron_e/electron.E(),electron.Eta(),electron.Phi(),0.000510998910);
	int echarge = electron.charge();
	if (m_random.Uniform(1.0)< m_upgrade->getElectronChargeFlipProb(electron.Pt()*1000., electron.Eta())) echarge*=-1;
	smeared->addElectron(eLV,echarge,electron.id(),idx);
	if (m_random.Uniform(1.0)<m_upgrade->getElectronToPhotonFakeRate(electron.Pt()*1000., electron.Eta())) 
	  smeared->addPhoton(eLV,PhotonIsoGood,-2);
      }
      idx++;
    }
    else {
      smeared->addElectron(electron,electron.charge(),electron.id(),idx++);
    }
  }
  idx=0;
  for(const auto& muon : muons) {
    if (smearMuons){
      float eff = m_upgrade->getMuonEfficiency(muon.Pt()*1000., muon.Eta());
      if (m_random.Uniform(1.0)<eff) {
	float muonUnsmearedPt = muon.Pt()*1000.;
	float qoverpt = muon.charge() / muonUnsmearedPt;
	float muonQOverPtResolution = m_upgrade->getMuonQOverPtResolution(muonUnsmearedPt, muon.Eta());
	qoverpt += m_random.Gaus(0., muonQOverPtResolution);
	float muonSmearedPt = fabs(1./qoverpt)/1000.;
	int muonCharge = 1;
	if (qoverpt<0) muonCharge = -1;
	TLorentzVector mLV;
	mLV.SetPtEtaPhiM(muonSmearedPt, muon.Eta(), muon.Phi(), 0.1056583715);
      	smeared->addMuon(mLV,muonCharge,muon.id(),idx);
      }
      idx++;
    }
    else {
      smeared->addMuon(muon,muon.charge(),muon.id(),idx++);
    }
  }
  idx=0;
  for(const auto& tau : taus) {
    if (smearTaus){
      short prong=1;
      if (tau.pass(TauThreeProng)) prong=3;
      float eff = m_upgrade->getTauEfficiency(tau.Pt()*1000., tau.Eta(), prong);
      if (m_random.Uniform(1.0)<eff) {
	float tau_E = m_upgrade->getTauSmearedEnergy(tau.E()*1000., tau.Eta(), prong)/1000.; 
	TLorentzVector tauLV;
	tauLV.SetPtEtaPhiM(tau.Pt()*tau_E/tau.E(),tau.Eta(),tau.Phi(),1.777682);
	smeared->addTau(tauLV,tau.charge(),tau.id(),idx);
      }
      idx++;
    }
    else {
      smeared->addTau(tau,tau.charge(),tau.id(),idx++);
    }
  }
  idx=0;
  for(const auto& photon : photons) {
    if (smearPhotons){
      float eff = m_upgrade->getPhotonEfficiency(photon.Pt()*1000./*, photon.Eta()*/);
      if (m_random.Uniform(1.0)<eff) {
	TLorentzVector pLV;
	pLV.SetPtEtaPhiM(photon.Pt()*1000.,photon.Eta(),photon.Phi(),0.);
	pLV = m_upgrade->getPhotonSmearedVector(&pLV);
	pLV.SetPtEtaPhiM(pLV.Pt()/1000.,pLV.Eta(),pLV.Phi(),0.);
	smeared->addPhoton(pLV,photon.id(),idx);
      }
      idx++;
    }
    else {
      smeared->addPhoton(photon,photon.id(),idx++);
    }
  }
  idx=0;
  for(const auto& jet : jets) {
    if (smearJets){
      float jetpt  = jet.Pt();
      if (jetpt<1500) jetpt = m_upgrade->getJetSmearedEnergy(jet.Pt()*1000.,jet.Eta(),true)/1000.; // FIXME: can only smear jets below 1500 GeV
      float jeteta = jet.Eta();
      float jetphi = jet.Phi();
      float jetE   = jet.E()*jetpt/jet.Pt();
      char jetType = 'L';
      if (jet.pass(TrueBJet)) jetType = 'B';
      if (jet.pass(TrueCJet)) jetType = 'C';
      
      float tagEff70 = m_upgrade->getFlavourTagEfficiency(jetpt*1000., jeteta, jetType, "mv1", 70, m_upgrade->getPileupTrackConfSetting());
      float tagEff85 = m_upgrade->getFlavourTagEfficiency(jetpt*1000., jeteta, jetType, "mv1", 85, m_upgrade->getPileupTrackConfSetting());
      float tag=m_random.Uniform(1.0);
      int jetid=GoodJet;
      if (tag<tagEff85) jetid|=BTag85MV2c20; //FIXME: check if this should set other working points too
      if (tag<tagEff70) jetid = GoodBJet;
      if (jet.pass(TrueLightJet)) jetid|=TrueLightJet;
      if (jet.pass(TrueCJet))     jetid|=TrueCJet;
      if (jet.pass(TrueBJet))     jetid|=TrueBJet;
      if (jet.pass(TrueTau))      jetid|=TrueTau;

      if (addPileupJets) {
	if ( (jetpt*1000) < m_upgrade->getPileupJetPtThresholdMeV()) jetpt=0;
	else {
	  float trackEff = m_upgrade->getTrackJetConfirmEff(jetpt*1000.,jet.Eta(), "HS");
	  std::cout << "Spyros CHECK1: tceff = " << trackEff << std::endl;
	  float hsProb = m_random.Uniform(1.0);
	  if (hsProb > trackEff) {
	    jetpt=0; // FIXME: should couple this to JVT flag
	    std::cout << "Spyros CHECK1b: hsProb = " << hsProb << " , jet pt = " << jetpt << std::endl;
	  }  
	}
      }
      if (jetpt) {
	TLorentzVector j;
	j.SetPtEtaPhiE(jetpt,jeteta,jetphi,jetE);
	smeared->addJet(j,jetid,idx);
	// Add jets faking electrons
	if (m_random.Uniform(1.0)<m_upgrade->getElectronFakeRate(jet.Pt()*1000, jet.Eta())) {
	  float electron_e = m_upgrade->getElectronFakeRescaledEnergy(jet.E()*1000., jet.Eta())/1000.;
	  TLorentzVector eLV;
	  eLV.SetPtEtaPhiM(jet.Pt()*electron_e/jet.E(),jet.Eta(),jet.Phi(),0.000510998910);
	  int echarge = 1;
	  if (m_random.Uniform(1.0)<0.5) echarge = -1;
	  smeared->addElectron(eLV,echarge,EIsoGood,-1);
	}
	// Add jets faking tau
	short tau_prong = 1; // FIXME: how to get prong distribution? Assuming 1-prong is likely overestimate
	if (jetpt<20000 || fabs(jet.Eta())>2.5) continue;
	if (m_random.Uniform(1.0)<m_upgrade->getTauFakeRate(jetpt*1000, jet.Eta(),tau_prong)) {
	  float tau_et = jetpt; // FIXME: no tau smearing exists yet
	  TLorentzVector tauLV;
	  tauLV.SetPtEtaPhiM(tau_et,jet.Eta(),jet.Phi(),1.777682);
	  int taucharge = 1;
	  if (m_random.Uniform(1.0)<0.5) taucharge = -1;
	  smeared->addTau(tauLV,taucharge,TauIsoGood,-1);
	}
	// Add jets faking photon
	if (m_random.Uniform(1.0)<m_upgrade->getPhotonFakeRate(jet.Pt()*1000/*, jet.Eta()*/)) {
	  float photon_et = m_upgrade->getPhotonFakeRescaledET(jet.Pt()*1000/*., jet.Eta()*/)/1000.;
	  TLorentzVector pLV;
	  pLV.SetPtEtaPhiM(photon_et,jet.Eta(),jet.Phi(),0.0);
	  smeared->addPhoton(pLV,PhotonIsoGood,-1);
	}
      }
      idx++;

    }
    else {
      smeared->addJet(jet,jet.id(),idx++);
    }
  }
  if (addPileupJets) {
    for (const auto& pujet : m_upgrade->getPileupJets()) {
      float trackEff = m_upgrade->getTrackJetConfirmEff(pujet.Pt(), pujet.Eta(), "PU");
      std::cout << "Spyros CHECK2: tceff = " << trackEff << std::endl;
      float puProb = m_random.Uniform(1.0);
            
      if (puProb > trackEff) continue; // FIXME: should couple this to JVT flag
      float tagEff70 = m_upgrade->getFlavourTagEfficiency(pujet.Pt(), pujet.Eta(), 'P', "mv1", 70, m_upgrade->getPileupTrackConfSetting());
      float tagEff85 = m_upgrade->getFlavourTagEfficiency(pujet.Pt(), pujet.Eta(), 'P', "mv1", 85, m_upgrade->getPileupTrackConfSetting());
      float tag=m_random.Uniform(1.0);
      int jetid=GoodJet;
      if (tag<tagEff85) jetid|=BTag85MV2c20; //FIXME: check if this should set other working points too
      if (tag<tagEff70) jetid=GoodBJet;
      smeared->addJet(pujet.Px()/1000., pujet.Py()/1000., pujet.Pz()/1000., pujet.E()/1000., jetid, -1);
      // Add jets faking photon
      if (m_random.Uniform(1.0)<m_upgrade->getPhotonPileupFakeRate(pujet.Pt()/*, pujet.Eta()*/)) {
	float photon_et = m_upgrade->getPhotonPileupFakeRescaledET(pujet.Pt()/*., pujet.Eta()*/)/1000.;
	TLorentzVector pLV;
	pLV.SetPtEtaPhiM(photon_et,pujet.Eta(),pujet.Phi(),0.0);
	smeared->addPhoton(pLV,PhotonIsoGood,-3);
      }
    }
  }
  for(const auto& jet : fatjets) {
    smeared->addFatJet(jet,jet.id(),idx++); //FIXME: for now there is no smearing for fat jets
  }

  return smeared;
#else
  return 0;
#endif

}
