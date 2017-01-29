#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(TwoBjets2015)

void TwoBjets2015::Init()
{
  addRegions({"SRA250","SRA350","SRA450","SRB"});
  addRegions({"CRzA","CRttA","CRstA","CRwA","CRzB","CRttB"});
}

void TwoBjets2015::ProcessEvent(AnalysisEvent *event)
{
    
  auto electrons  = event->getElectrons(10, 2.47, ELooseLH);
  auto muons      = event->getMuons(10, 2.5, MuMedium);
  auto candJets   = event->getJets(20., 2.8);
  auto metVec     = event->getMET();
  double met      = metVec.Et();
    
  if (countObjects(candJets, 20, 2.8, NOT(LooseBadJet))!=0) return;

  candJets = filterObjects(candJets, 20, 2.8, JVT50Jet);

  // Overlap removal
  // not doing overlap removal between electrons and muons with identical track
  candJets   = overlapRemoval(candJets,electrons,0.2,NOT(BTag85MV2c20));
  electrons  = overlapRemoval(electrons,candJets,0.4);
  candJets   = overlapRemoval(candJets, muons, 0.4, LessThan3Tracks); // CHECK: CONF is not clear that jet is removed
  muons      = overlapRemoval(muons, candJets, 0.4);

  auto baselineLeptons   = electrons + muons;

  auto signalJets      = filterObjects(candJets, 35);
  auto signalElectrons = filterObjects(electrons, 20, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoGradientLoose);
  auto signalMuons     = filterObjects(muons, 20, 2.5, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);
  auto signalLeptons   = signalElectrons + signalMuons;
  auto signalBjets     = filterObjects(signalJets, 35., 2.5, BTag77MV2c20);

  int nJets  = signalJets.size();
  int nBjets = signalBjets.size();
  if (nJets<2) return;
  if (nJets>=4 && signalJets[3].Pt()>50) return; 

  float meff2j    = met + sumObjectsPt(signalJets,2);
  float dphiMin1  = minDphi(metVec, signalJets, 1);
  float dphiMin4  = minDphi(metVec, signalJets, 4);
  float mbb = 0;
  float mCT = 0;
  if (nBjets>=2) {
    mbb = (signalBjets[0] + signalBjets[1]).M();
    mCT = calcMCT(signalBjets[0], signalBjets[1]);
  }
  bool bjetsLeading = (signalJets[0].pass(BTag77MV2c20) && signalJets[1].pass(BTag77MV2c20));
  bool bjetsSublead = (nJets>=3 && signalJets[0].pass(NOT(BTag77MV2c20)) && signalJets[1].pass(BTag77MV2c20) &&
		       (signalJets[2].pass(BTag77MV2c20) || (nJets>=4 && signalJets[3].pass(BTag77MV2c20))));

  // Signal region definitions
  if (baselineLeptons.size()==0 && nJets>=2 && dphiMin4>0.4 && nBjets==2) {
    if (bjetsLeading && met>250 && signalJets[0].Pt()>130 && signalJets[1].Pt()>50 &&
	met/meff2j>0.25 && mbb>200) {
      if (mCT>250) accept("SRA250");
      if (mCT>350) accept("SRA350");
      if (mCT>450) accept("SRA450");
    }
    if (bjetsSublead && met>400 && signalJets[0].Pt()>300 && signalJets[1].Pt()>50 &&
	dphiMin1>2.5 && met/meff2j>0.25) accept("SRB");
  }
  // 1-lepton Control region definitions
  if (signalLeptons.size()==1 && baselineLeptons.size()==1 && signalLeptons[0].Pt()>26 && met>100)  {
    float mT = calcMT(signalLeptons[0], metVec);
    if (bjetsLeading && mCT>150 && nBjets==2) {
      float mblmin =std::min((signalBjets[0] + signalLeptons[0]).M(), (signalBjets[1] + signalLeptons[0]).M());
      if (mbb<200 && signalJets[0].Pt()>130) accept("CRttA");
      if (mbb>200 && mblmin>170)             accept("CRstA");
    } 
    float mbj=(signalJets[0] + signalJets[1]).M();
    if (signalJets[0].pass(BTag77MV2c20) && signalJets[0].Pt()>130 && mT>30 && mbj>200 && mCT>150) accept("CRwA");
    if (bjetsSublead &&	met>200 && signalJets[0].Pt()>130 && dphiMin1>2.5 && nBjets==2) accept("CRttB");
  }

  // 2-lepton Control region definitions
  if (signalLeptons.size()==2 && baselineLeptons.size()==2 &&
      signalLeptons[0].type()==signalLeptons[1].type() &&
      signalLeptons[0].charge()!=signalLeptons[1].charge() &&
      signalLeptons[0].Pt()>26 && nBjets==2) {
    float mll = (signalLeptons[0] + signalLeptons[1]).M();
    float metCorr = (metVec + signalLeptons[0] + signalLeptons[1]).Pt();
    if (mll<106 && mll>76) {
      if (bjetsLeading && met<100 && metCorr>100 && mCT>150) accept("CRzA");
      if (signalJets[0].Pt()>50 && bjetsSublead && met<70 && metCorr>100 && dphiMin1>2.0) accept("CRzB");
    }
  }

  return;
}
