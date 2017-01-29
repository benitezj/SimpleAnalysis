#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(OneLepton2015)

void OneLepton2015::Init()
{
  addRegions({"SRh4jhx","SRh4jlx","SRh5j","SRh6j","SRs2j","SRs5j"});
  addRegions({"CRWh4jhx","CRWh4jlx","CRWh5j","CRWh6j","CRWs2j","CRWs5j"});
  addRegions({"CRTh4jhx","CRTh4jlx","CRTh5j","CRTh6j","CRTs2j","CRTs5j"});
}

void OneLepton2015::ProcessEvent(AnalysisEvent *event)
{
  auto hardElectrons  = event->getElectrons(10, 2.47, ELooseLH);
  auto softElectrons  = event->getElectrons(7, 2.47, ELooseLH);
  auto hardMuons      = event->getMuons(10,2.4, MuMedium);
  auto softMuons      = event->getMuons(6,2.4, MuMedium);
  auto preJets       = event->getJets(20., 4.5); 
  auto metVec     = event->getMET();
  float met = metVec.Pt();

  // Reject events with bad jets
  if (countObjects(preJets, 20, 4.5, NOT(LooseBadJet))!=0) return;

  // hard lepton overlap removal and signal lepton/jet definitions
  auto hardJets            = overlapRemoval(preJets, hardElectrons, 0.2);
       hardElectrons       = overlapRemoval(hardElectrons, hardJets, 0.4);
       hardElectrons       = overlapRemoval(hardElectrons, hardMuons, 0.01);
       hardJets            = overlapRemoval(hardJets, hardMuons, 0.4, LessThan3Tracks);
       hardMuons           = overlapRemoval(hardMuons, hardJets, 0.4);
       hardJets            = filterObjects(hardJets, 25, 2.8, JVT50Jet);
  auto hardBJets           = filterObjects(hardJets, 25, 2.8, BTag77MV2c20);
  auto hardSignalElectrons = filterObjects(hardElectrons, 35, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoGradientLoose);
  auto hardSignalMuons     = filterObjects(hardMuons, 35, 2.4, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);

  // Hard Lepton signal and control regions
  auto hardLeptons=hardElectrons+hardMuons;
  auto hardSignalLeptons=hardSignalElectrons+hardSignalMuons;
  if (hardLeptons.size()==1 && hardSignalLeptons.size()==1) {
    float mt   = calcMT(hardSignalLeptons[0], metVec);
    float ht   = hardSignalLeptons[0].Pt() + sumObjectsPt(hardJets,999,30);
    float meff = ht + met;
    float Ap   = aplanarity(hardJets);

    if (hardJets.size()>=4 && 
	hardJets[0].Pt()>325 && hardJets[3].Pt()>30 &&
	met>200 && mt>425 && met/meff>0.3 && meff>1800)            accept("SRh4jhx"); 
    if (hardJets.size()>=4 && 
	hardJets[0].Pt()>325 && hardJets[3].Pt()>150 &&
	met>200 && mt>125 && meff>2000 && Ap>0.04)                 accept("SRh4jlx"); 
    if (hardJets.size()>=5 && 
	hardJets[0].Pt()>225 && hardJets[3].Pt()>50 &&
	met>250 && mt>275 && met/meff>0.1 && meff>1800 && Ap>0.04) accept("SRh5j"); 
    if (hardJets.size()>=6 && 
	hardJets[0].Pt()>125 && hardJets[3].Pt()>30 &&
	met>250 && mt>225 && met/meff>0.2 && meff>1000 && Ap>0.04) accept("SRh6j"); 

    std::string type="CRW";
    if (hardBJets.size()>0) type="CRT";
    if (hardJets.size()>=4 && 
	hardJets[0].Pt()>325 && hardJets[3].Pt()>30 &&
	met>200 && mt<125 && mt>60  && met/meff<0.3 && meff>1800)   accept(type+"h4jhx"); 
    if (hardJets.size()>=4 && 
	hardJets[0].Pt()>325 && hardJets[3].Pt()>30 && 
	Ap<0.04 && met>200 && mt<125 && 
	((hardBJets.size()>=1 && mt>95 && meff>1500) || 
	 (hardBJets.size()==0 && mt>60 && meff>1900)))              accept(type+"h4jlx"); 
    if (hardJets.size()>=5 && 
	hardJets[0].Pt()>225 && hardJets[3].Pt()>30 && 
	met>250 && mt<125 && mt>60 && met/meff>0.1 && meff>1500 && Ap<0.04) accept(type+"h5j"); 
    if (hardJets.size()>=6 && 
	hardJets[0].Pt()>125 && hardJets[3].Pt()>30 &&
	met>250 && mt<125 && mt>60 && met/meff>0.2 && meff>1000 && Ap<0.04) accept(type+"h6j"); 
  }

  // Soft lepton overlap removal and signal lepton/jet definitions
  auto softJets            = overlapRemoval(preJets, softElectrons, 0.2);
       softElectrons       = overlapRemoval(softElectrons, softJets, 0.4);
       softElectrons       = overlapRemoval(softElectrons, softMuons, 0.01);
       softJets            = overlapRemoval(softJets, softMuons, 0.4, LessThan3Tracks);
       softMuons           = overlapRemoval(softMuons, softJets, 0.4);
       softJets            = filterObjects(softJets, 25, 2.8, JVT50Jet);
  auto softBJets           = filterObjects(softJets, 25, 2.5, BTag77MV2c20);
  auto softSignalElectrons = filterObjects(softElectrons, 7, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoGradientLoose);
  auto softSignalMuons     = filterObjects(softMuons, 6, 2.4, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);

  // Soft Lepton signal and control regions
  auto softLeptons=softElectrons+softMuons;
  auto softSignalLeptons=softSignalElectrons+softSignalMuons;
  if (softLeptons.size()==1 && softSignalLeptons.size()==1 && softSignalLeptons[0].Pt()<35) {
    float mt   = calcMT(softSignalLeptons[0], metVec);
    float ht   = softSignalLeptons[0].Pt() + sumObjectsPt(softJets,999,30);
    float meff = ht + met;
    float Ap   = aplanarity(softJets);

    if (softJets.size()>=2 && 
	softJets[0].Pt()>180 && softJets[1].Pt()>30 && // CHECK: int note says 25 GeV, not 30
	met>530 && mt>100 && met/meff>0.38) accept("SRs2j"); 
    if (softJets.size()>=5 && 
	softJets[2].Pt()>200 && softJets[4].Pt()>30 && // CHECK: int note says 25 GeV, not 30
	met>375 && ht>1100 && Ap>0.02)      accept("SRs5j"); 

    std::string type="CRW";
    if (softBJets.size()>0) type="CRT";
    if (softJets.size()>=2 && 
	softJets[0].Pt()>180 && softJets[1].Pt()>30 && // CHECK: int note says 25 GeV, paper say 30 GeV
	met<500 && met>250 && mt<80 && mt>40 && met/meff>0.38) accept(type+"s2j"); 
    if (softJets.size()>=5 && 
	softJets[0].Pt()>150 && softJets[1].Pt()>100 && softJets[4].Pt()>30 && // CHECK: int note says 25 GeV, paper says 30 GeV
	met>375 && met>250 && ht<1000 && ht>500 && Ap>0.02)    accept(type+"s5j"); 
  }

  return;
}
