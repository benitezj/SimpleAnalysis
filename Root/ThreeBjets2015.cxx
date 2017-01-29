#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(ThreeBjets2015)

void ThreeBjets2015::Init()
{
  addRegions({"SR_Gbb_A","SR_Gbb_B","SR_Gbb_C",
	               "SR_Gtt_0L_A","SR_Gtt_0L_B","SR_Gtt_0L_C",
	               "SR_Gtt_1L_A","SR_Gtt_1L_B"});
  addRegions({"CR_Gbb_A","CR_Gbb_B","CR_Gbb_C",
	               "CR_Gtt_0L_A","CR_Gtt_0L_B","CR_Gtt_0L_C",
	               "CR_Gtt_1L_A","CR_Gtt_1L_B"});
  addHistogram("MET",100,0,2000);
}

void ThreeBjets2015::ProcessEvent(AnalysisEvent *event)
{
    
  auto electrons  = event->getElectrons(20, 2.47, ELooseLH);
  auto muons      = event->getMuons(20, 2.5, MuMedium);
  auto candJets   = event->getJets(20., 2.8);
  auto metVec     = event->getMET();
  double met      = metVec.Et();

  fill("MET",met);
  if (countObjects(candJets, 20, 2.8, NOT(LooseBadJet))!=0) return;
  // No bad muon veto implemented

  candJets = filterObjects(candJets, 20, 2.8, JVT50Jet);
  auto fatJets = reclusterJets(candJets, 1.0, 300, 0.2, 0.05); 
  fatJets = filterObjects(fatJets, 300, 2.0);

  // Overlap removal
  auto radiusCalcJet  = [] (const AnalysisObject& , const AnalysisObject& muon) { return std::min(0.4, 0.04 + 10/muon.Pt()); };
  auto radiusCalcMuon = [] (const AnalysisObject& muon, const AnalysisObject& ) { return std::min(0.4, 0.04 + 10/muon.Pt()); };
  electrons  = overlapRemoval(electrons, muons, 0.01);
  candJets   = overlapRemoval(candJets, electrons, 0.2, NOT(BTag85MV2c20));
  electrons  = overlapRemoval(electrons, candJets, 0.4);
  candJets   = overlapRemoval(candJets, muons, radiusCalcJet, LessThan3Tracks); 
  muons      = overlapRemoval(muons, candJets, radiusCalcMuon);                  

  // No cosmic muon veto implemented

  auto signalJets      = filterObjects(candJets, 30);
  auto signalElectrons = filterObjects(electrons, 20, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoBoosted);
  auto signalMuons     = filterObjects(muons, 20, 2.5, MuD0Sigma3|MuZ05mm|MuIsoBoosted|MuNotCosmic);
  auto signalLeptons   = signalElectrons + signalMuons;

  // get b-jets - need at least two 
  auto bjets = filterObjects(signalJets, 30., 2.5, BTag85MV2c20);
  if (bjets.size()<2 || signalJets.size()<4) return;

  int nTop=0;
  for(const auto& jet : fatJets) {
    if (jet.M()>100) nTop++;
  }

  float meff4j    = met + sumObjectsPt(signalJets,4);
  float meffIncl  = met + sumObjectsPt(signalJets) + sumObjectsPt(signalLeptons);
  float mTbmin   = calcMTmin(bjets, metVec,3);
  float dphiMin4  = minDphi(metVec, signalJets, 4);
  float mT = 0.;
  if(signalLeptons.size()>0) mT=calcMT(signalLeptons[0], metVec);

  int nBjets = bjets.size();
  int nJets = signalJets.size();
  int nLeptons = signalLeptons.size();
  int nCandLeptons = electrons.size() + muons.size();

  // Gbb signal and control regions
  if (nBjets>=3) {
    if (dphiMin4>0.4 && nCandLeptons==0) {
      if (bjets[2].Pt()>90 && signalJets[3].Pt()>90 && met>350 && meff4j>1600) accept("SR_Gbb_A");
      if (bjets[2].Pt()>90 && signalJets[3].Pt()>90 && met>450 && meff4j>1400) accept("SR_Gbb_B");
      if (                                             met>500 && meff4j>1400) accept("SR_Gbb_C");
    }
    if (nLeptons==1 && mT<150) {
      if (bjets[2].Pt()>90 && signalJets[3].Pt()>90 && met>250 && meff4j>1400) accept("CR_Gbb_A");
      if (bjets[2].Pt()>90 && signalJets[3].Pt()>90 && met>300 && meff4j>1000) accept("CR_Gbb_B");
      if (                                             met>400 && meff4j>1200) accept("CR_Gbb_C");
    }
  }
  // Gtt-0L signal and control regions
  if (dphiMin4>0.4 && nJets>=8 && mTbmin>80 && nCandLeptons==0) {
    if (met>400 && meffIncl>1700 && nBjets>=3 && nTop>=1) accept("SR_Gtt_0L_A");
    if (met>350 && meffIncl>1250 && nBjets>=4 && nTop>=1) accept("SR_Gtt_0L_B");
    if (met>350 && meffIncl>1250 && nBjets>=4           ) accept("SR_Gtt_0L_C");
  }
  if (nLeptons==1 && nJets>=7 && mT<150) {
    if (met>250 && meffIncl>1350 && nBjets>=3 && nTop>=1) accept("CR_Gtt_0L_A");
    if (met>200 && meffIncl>1000 && nBjets>=4 && nTop>=1) accept("CR_Gtt_0L_B");
    if (met>200 && meffIncl>1000 && nBjets>=4           ) accept("CR_Gtt_0L_C");
  }
  // Gtt-1L signal and control regions
  if (nLeptons>=1 && nBjets>=3) {
    if (mT>150 && nJets>=6 && mTbmin>160) {
      if (met>200 && meffIncl>1100 && nTop>=1) accept("SR_Gtt_1L_A");
      if (met>300 && meffIncl>900)             accept("SR_Gtt_1L_B");
    }
    if (mT<150 && nJets==6) {
      if (met>200 && meffIncl>1100 && nTop>=1) accept("CR_Gtt_1L_A");
      if (met>300 && meffIncl>900)             accept("CR_Gtt_1L_B");
    }
  }
  ntupVar("met",met);
  ntupVar("signalJets",signalJets);
  return;
}
