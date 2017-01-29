#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(OneLepton2012)

void OneLepton2012::Init()
{
  addRegions({"h1L3j","h1L3j_bin1","h1L3j_bin2","h1L3j_bin3","h1L3j_bin4"});
  addRegions({"h1L5j","h1L5j_bin1","h1L5j_bin2","h1L5j_bin3","h1L5j_bin4"});
  addRegions({"h1L6j","h1L6j_bin1","h1L6j_bin2","h1L6j_bin3","h1L6j_bin4"});
  addRegions({"s1L3j","s1L3j_bin1","s1L3j_bin2","s1L3j_bin3","s1L3j_bin4"});
  addRegions({"s1L5j","s1L5j_bin1","s1L5j_bin2","s1L5j_bin3","s1L5j_bin4"});
  addRegions({"s1L3jinc","s1L3jinc_bin1","s1L3jinc_bin2","s1L3jinc_bin3","s1L3jinc_bin4"});
  addRegion("s2M2j");
}

void OneLepton2012::ProcessEvent(AnalysisEvent *event)
{
  auto hardElectrons  = event->getElectrons(10,2.47);
  auto softElectrons  = event->getElectrons(7,2.47);
  auto hardMuons      = event->getMuons(10,2.4);
  auto softMuons      = event->getMuons(6,2.4);
  auto preJets       = event->getJets(20.,2.8); //FIXME: paper lists 2.5
  auto metVec     = event->getMET();
  float met = metVec.Pt();

  // overlap removal and signal lepton/jet definitions
  auto goodHardJets       = overlapRemoval(preJets,hardElectrons,0.2);
  auto signalHardJets     = filterObjects(goodHardJets,30,2.5);  //FIXME: paper list 25
  auto goodSoftJets       = overlapRemoval(preJets,softElectrons,0.2);
  auto signalSoftJets     = filterObjects(goodSoftJets,25,2.5); 

  auto bjets = filterObjects(signalSoftJets,25.,2.5,true);

  hardMuons = overlapRemoval(hardMuons,goodHardJets,0.4);
  hardElectrons = overlapRemoval(hardElectrons,goodHardJets,0.4);
  softMuons = overlapRemoval(softMuons,goodHardJets,0.4);
  softElectrons = overlapRemoval(softElectrons,goodHardJets,0.4);

  // electron crack vetos
  bool hard_electron_veto=false;
  for(const auto& elec : hardElectrons ) {
    float eta=fabs(elec.Eta());
    if (eta>=1.37 && eta<=1.52 && elec.Pt()>25) hard_electron_veto=true; //FIXME: paper lists 10 GeV
  }

  bool soft_electron_veto=false;
  for(const auto& elec : softElectrons ) {
    float eta=fabs(elec.Eta());
    if (eta>=1.37 && eta<=1.52) soft_electron_veto=true;
  }

  // Hard Lepton signal regions
  auto hardLeptons=hardElectrons+hardMuons;
  if (!hard_electron_veto && 
      hardLeptons.size()>=1 && hardLeptons[0].Pt()>25 &&
      (hardLeptons.size()==1 || hardLeptons[1].Pt()<25)) { //FIXME: paper lists 10 GeV
    float mt=calcMT(hardLeptons[0],metVec);
    float meffIncl=met+sumObjectsPt(hardLeptons,100,25)+sumObjectsPt(signalHardJets); //FIXME: paper lists 10 GeV
    float meffExcl=met+sumObjectsPt(hardLeptons,100,25)+sumObjectsPt(signalHardJets,3); //FIXME: paper lists 10 GeV
    if ( signalHardJets.size()>=3 && signalHardJets[1].Pt()>80 && met>300 && 
	 mt>150 && met/meffExcl>0.3 && meffIncl>800 ) {
      if ( met>500 && meffIncl>1400) accept("h1L3j");
      if ( signalHardJets.size()<5 || signalHardJets[4].Pt()<40 ) {
	if ( meffIncl>1400 )      accept("h1L3j_bin4");
	else if ( meffIncl>1200 ) accept("h1L3j_bin3");
	else if ( meffIncl>1000 ) accept("h1L3j_bin2");
	else                      accept("h1L3j_bin1");
      }
    }
    if ( signalHardJets.size()>=5 && 
	 signalHardJets[0].Pt()>80 && signalHardJets[1].Pt()>50 && signalHardJets[4].Pt()>40 &&
	 met>300 && mt>150 && meffIncl>800 ) {
      if ( mt>200 && meffIncl>1400 ) accept("h1L5j");	
      if ( signalHardJets.size()<6 || signalHardJets[5].Pt()<40 ) {
	if ( meffIncl>1400 )      accept("h1L5j_bin4");
	else if ( meffIncl>1200 ) accept("h1L5j_bin3");
	else if ( meffIncl>1000 ) accept("h1L5j_bin2");
	else                      accept("h1L5j_bin1");
      }
    }
    if ( signalHardJets.size()>=6 && 
	 signalHardJets[0].Pt()>80 && signalHardJets[1].Pt()>50 && signalHardJets[5].Pt()>40 &&
	 met>250 && mt>150 && meffIncl>600 ) {
      if ( met>350 ) accept("h1L6j");	
      if ( met>450 )      accept("h1L6j_bin3");
      else if ( met>350 ) accept("h1L6j_bin2");
      else                accept("h1L6j_bin1");
    }
  } 

  // Soft one lepton signal regions
  auto softLeptons=softElectrons+softMuons;
  sortObjectsByPt(softLeptons);
  if (!soft_electron_veto && 
      softLeptons.size()==1 && softLeptons[0].Pt()<25 ) {
    auto softLeptonsSR3 = overlapRemoval(softLeptons,goodSoftJets,1.0);  
    float mt=calcMT(softLeptons[0],metVec);
    float meffIncl=met+sumObjectsPt(softLeptons)+sumObjectsPt(signalSoftJets);
    if ( softLeptonsSR3.size()==1 && 
	 signalSoftJets.size()>=3 && signalSoftJets.size()<5 && signalSoftJets[0].Pt()>180 &&
	 met>400 && mt>100 && met/meffIncl>0.1) {
      if ( met/meffIncl>0.3) accept("s1L3j");
      if ( met/meffIncl>0.4)      accept("s1L3j_bin4");
      else if ( met/meffIncl>0.3) accept("s1L3j_bin3");
      else if ( met/meffIncl>0.2) accept("s1L3j_bin2");
      else                        accept("s1L3j_bin1");
    }
    if ( signalSoftJets.size()>=5 && signalSoftJets[0].Pt()>180 &&
	 met>300 && mt>100 && met/meffIncl>0.1) {
      if ( met/meffIncl>0.3) accept("s1L5j");
      if ( met/meffIncl>0.4)      accept("s1L5j_bin4");
      else if ( met/meffIncl>0.3) accept("s1L5j_bin3");
      else if ( met/meffIncl>0.2) accept("s1L5j_bin2");
      else                        accept("s1L5j_bin1");
    }
    if ( signalSoftJets.size()>=3 && signalSoftJets[0].Pt()>130 && signalSoftJets[1].Pt()>100 &&
	 bjets.size()==0 && met>180 && mt>120 && met/meffIncl>0.1) {
      if ( met/meffIncl>0.3) accept("s1L3jinc");
      if ( met/meffIncl>0.4)      accept("s1L3jinc_bin4");
      else if ( met/meffIncl>0.3) accept("s1L3jinc_bin3");
      else if ( met/meffIncl>0.2) accept("s1L3jinc_bin2");
      else                        accept("s1L3jinc_bin1");
    }
  }
  // Soft dimuon signal region 
  //FIXME: not validated as it was not in pMSSM code
  if (!soft_electron_veto && softMuons.size()==2 && softLeptons.size()==2 &&
      softMuons[1].Pt()<25) {
    float mass=(softMuons[0]+softMuons[1]).M();
    float mt=calcMT(softLeptons[1],metVec);
    float meffIncl=met+sumObjectsPt(softLeptons)+sumObjectsPt(signalSoftJets);
    //check that muon has deltaR>1 to all jets
    AnalysisObjects softSecondMuon;
    softSecondMuon.push_back(softMuons[1]);
    auto softMuonDR1 = overlapRemoval(softSecondMuon,goodSoftJets,1.0);  
    if ( softMuonDR1.size()==1 && mass>15 && mass<60 &&
	 signalSoftJets.size()>=2 && signalSoftJets[0].Pt()>80 && bjets.size()==0 && 
	 met>180 && mt>40 && met/meffIncl>0.3)
      accept("s2M2j");
  }
  //FIXME: missing hard dilepton (razor)
  return;
}
