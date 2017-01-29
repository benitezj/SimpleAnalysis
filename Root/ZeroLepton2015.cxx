#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(ZeroLepton2015)

void ZeroLepton2015::Init()
{
  addRegions({"SR2jl","SR2jm","SR2jt","SR4jt","SR5j","SR6jm","SR6jt"});
  addRegions({"CRW2jl","CRW2jm","CRW2jt","CRW4jt","CRW5j","CRW6jm","CRW6jt"});
  addRegions({"CRT2jl","CRT2jm","CRT2jt","CRT4jt","CRT5j","CRT6jm","CRT6jt"});
}

void ZeroLepton2015::ProcessEvent(AnalysisEvent *event)
{
  auto electrons  = event->getElectrons(10, 2.47, ELooseLH);
  auto muons      = event->getMuons(10, 2.7, MuMedium);
  auto jets       = event->getJets(20., 2.8);
  auto metVec     = event->getMET();
  double met      = metVec.Et();
    
  // Reject events with bad jets
  if (countObjects(jets, 20, 2.8, NOT(LooseBadJet))!=0) return;
  if (jets.size()>0 && jets[0].Pt()>100 && jets[0].pass(NOT(TightBadJet))) return;
  if (jets.size()>1 && jets[1].Pt()>100 && jets[1].pass(NOT(TightBadJet))) return;

  // Standard SUSY overlap removal
  jets       = overlapRemoval(jets, electrons, 0.2);
  electrons  = overlapRemoval(electrons, jets, 0.4);
  muons      = overlapRemoval(muons, jets, 0.4);

  // not doing overlap removal between electrons and muons with identical track
  // FIXME: not doing overlap removal between electrons within DeltaR=0.5

  auto signalMuons     = filterObjects(muons, 25, 2.7, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);
  auto signalElectrons = filterObjects(electrons, 25, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoGradientLoose);
  auto goodJets        = filterObjects(jets, 50);
  auto bjets           = filterObjects(goodJets, 50, 2.5, BTag77MV2c20);

  auto leptons=signalElectrons + signalMuons;
  if (leptons.size()==1 && leptons[0].Pt()>50) { 
    // For CR treat lepton as a jet, but only count if above 50 GeV
    goodJets = goodJets + leptons;
  }

  // Define Meff, delta phi 
  float meff[7];
  int nJets=goodJets.size();
  for(int nJet=2; nJet<=6; nJet++) 
    meff[nJet] = sumObjectsPt(goodJets, nJet) + met;
  
  float meffIncl    = sumObjectsPt(jets) + met;
  float metSig      = met/sqrt(meffIncl - met);
  float dphiMin3    = minDphi(metVec, jets, 3);
  float dphiMinRest = minDphi(metVec, jets);
  float Ap          = aplanarity(jets);

  // Signal regions
  if (muons.size()+electrons.size()==0 && met>200 && nJets>1 && dphiMin3>0.4 && goodJets[0].Pt()>200) {
    if (nJets >= 2) { //2-jet SRs
	if (goodJets[1].Pt()>200 && dphiMin3>0.8 && metSig>=15 && meffIncl>1200)  accept("SR2jl");
	if (goodJets[0].Pt()>300 && goodJets[1].Pt()>50 && metSig>=15 && meffIncl>1600) accept("SR2jm");
	if (goodJets[1].Pt()>200 && dphiMin3>0.8 && metSig>=20 && meffIncl>2000)  accept("SR2jt");
    }
    if (dphiMinRest>0.2 && Ap>0.04) {
	if (nJets >= 4 ) { //4-jet SR
	  if (goodJets[3].Pt()>100 && met/meff[4]>0.2 && meffIncl>2200)             accept("SR4jt");
	}
	if (nJets >= 5 ) { //5-jet SR
	  if (goodJets[3].Pt()>100 && met/meff[5]>0.25 && meffIncl>1600)            accept("SR5j");
	}
	if (nJets >= 6 ) { //6-jet SR
	  if (goodJets[3].Pt()>100 && met/meff[6]>0.25 && meffIncl>1600)            accept("SR6jm");
	  if (goodJets[3].Pt()>100 && met/meff[6]>0.2  && meffIncl>2000)            accept("SR6jt");
	}
    }
  } //end signal regions

  // Control regions
  if (leptons.size()==1) {
    float mt=calcMT(leptons[0], metVec);
    std::string type="CRW";
    if (bjets.size()>0) type="CRT";
    if (mt>30 && mt<100 && met>200 && nJets>1 && goodJets[0].Pt()>200) {
	if (nJets >= 2) { //2-jet CRs
	  if (goodJets[1].Pt()>200 && meffIncl>1200)                   accept(type+"2jl");
	  if (goodJets[0].Pt()>300 && goodJets[1].Pt()>50 && meffIncl>1600) accept(type+"2jm");
	  if (goodJets[1].Pt()>200 && meffIncl>2000)                   accept(type+"2jt");
	}
	if (Ap>0.04) {
	  if (nJets >= 4 ) { //4-jet CR
	    if (goodJets[3].Pt()>100 && meffIncl>2200)                 accept(type+"4jt");
	  }
	  if (nJets >= 5 ) { //5-jet CR
	    if (goodJets[3].Pt()>100 && meffIncl>1600)                 accept(type+"5j");
	  }
	  if (nJets >= 6 ) { //6-jet SR
	    if (goodJets[3].Pt()>100 && meffIncl>1600)                 accept(type+"6jm");
	    if (goodJets[3].Pt()>100 && meffIncl>2000)                 accept(type+"6jt");
	  }
	}
    }
  }

  return;
}
