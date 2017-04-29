#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(OneLepton2016)

void OneLepton2016::Init()
{
  addRegions({"SRs2j1bV","SRs2j2bV","SRs2j3bV","SRs2j4bV"});
  addRegions({"SRh4jlx1bV","SRh4jlx2bV","SRh4jlx3bV"});
  addRegions({"SRh4jhx1bV","SRh4jhx2bV","SRh4jhx3bV"});
  addRegions({"SRh6j1bV","SRh6j2bV","SRh6j3bV","SRh6j4bV"});
  addRegions({"SRs2j1bT","SRs2j2bT","SRs2j3bT","SRs2j4bT"});
  addRegions({"SRh4jlx1bT","SRh4jlx2bT","SRh4jlx3bT"});
  addRegions({"SRh4jhx1bT","SRh4jhx2bT","SRh4jhx3bT"});
  addRegions({"SRh6j1bT","SRh6j2bT","SRh6j3bT","SRh6j4bT"});
  addRegions({"SRh9j1","SRh9j2"});
  addRegions({"SRs2j","SRh4jlxGG","SRh4jlxSS","SRh4jhx","SRh6jGG","SRh6jSS","SRh9j"});
}

void OneLepton2016::ProcessEvent(AnalysisEvent *event)
{
  auto hardElectrons  = event->getElectrons(7, 2.47, ELooseLH);
  auto softElectrons  = event->getElectrons(7, 2.47, ELooseLH);
  auto hardMuons      = event->getMuons(6,2.5, MuMedium);
  auto softMuons      = event->getMuons(6,2.5, MuMedium);
  auto preJets       = event->getJets(20., 2.8); 
  auto metVec     = event->getMET();
  float met = metVec.Pt();

  // Reject events with bad jets
  if (countObjects(preJets, 20, 4.5, NOT(LooseBadJet))!=0) return;

  // hard lepton overlap removal and signal lepton/jet definitions
  auto hardJets            = overlapRemoval(preJets, hardElectrons, 0.2, NOT(BTag77MV2c10));
       hardElectrons       = overlapRemoval(hardElectrons, hardJets, 0.4);
       hardElectrons       = overlapRemoval(hardElectrons, hardMuons, 0.01);
       hardJets            = overlapRemoval(hardJets, hardMuons, 0.4, LessThan3Tracks);
       hardMuons           = overlapRemoval(hardMuons, hardJets, 0.4);
       hardJets            = filterObjects(hardJets, 30, 2.8, JVT50Jet);
  auto hardBJets           = filterObjects(hardJets, 30, 2.5, BTag77MV2c10);
  auto hardSignalElectrons = filterObjects(hardElectrons, 35, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoGradientLoose);
  auto hardSignalMuons     = filterObjects(hardMuons, 35, 2.5, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);

  // Hard Lepton signal and control regions
  auto hardLeptons=hardElectrons+hardMuons;
  auto hardSignalLeptons=hardSignalElectrons+hardSignalMuons;
  if (hardLeptons.size()==1 && hardSignalLeptons.size()==1) {
    float mt   = calcMT(hardSignalLeptons[0], metVec);
    float ht   = sumObjectsPt(hardJets,999,30);
    float meff = ht + met+ hardSignalLeptons[0].Pt();
    float Ap   = aplanarity(hardJets);

    // b-Tag signal regions
    if (hardBJets.size()>0){
      if (hardJets.size()>=4 && hardJets.size()<=5 &&
	  met>300 && mt>450 && Ap>0.01 && met/meff>0.25 && meff>1000)
	{
	  if (meff<1500) accept("SRh4jhx1bT");
	  else if (meff<2000) accept("SRh4jhx2bT");
	  else accept("SRh4jhx3bT");
	}

      if (hardJets.size()>=4 && hardJets.size()<=5 &&
	  met>250 && mt>150 && mt<450 && Ap>0.05 && meff>1300)
	{
	  if (meff<1650) accept("SRh4jlx1bT");
	  else if (meff<2000) accept("SRh4jlx2bT");
	  else accept("SRh4jlx3bT");
	}
      
      if (hardJets.size()>=6 &&
	  met>350 && mt>175 && Ap>0.06 && meff>700)
	{
	  if (meff<1233) accept("SRh6j1bT");
	  else if (meff<1766) accept("SRh6j2bT");
	  else if (meff<2300) accept("SRh6j3bT");
	  else accept("SRh6j4bT");
	}
    }

    // b-Veto signal regions
    else{
      if (hardJets.size()>=4 && hardJets.size()<=5 &&
	  met>300 && mt>450 && Ap>0.01 && met/meff>0.25 && meff>1000)
	{
	  if (meff<1500) accept("SRh4jhx1bV");
	  else if (meff<2000) accept("SRh4jhx2bV");
	  else accept("SRh4jhx3bV");
	}
      
      if (hardJets.size()>=4 && hardJets.size()<=5 &&
	  met>250 && mt>150 && mt<450 && Ap>0.05 && meff>1300)
	{
	  if (meff<1650) accept("SRh4jlx1bV");
	  else if (meff<2000) accept("SRh4jlx2bV");
	  else accept("SRh4jlx3bV");
	}

      if (hardJets.size()>=6 &&
	  met>350 && mt>175 && Ap>0.06 && meff>700)
	{
	  if (meff<1233) accept("SRh6j1bV");
	  else if (meff<1766) accept("SRh6j2bV");
	  else if (meff<2300) accept("SRh6j3bV");
	  else accept("SRh6j4bV");
	}
      
    }
    // b-inclusive for SR9j and discovery SRs
    if (hardJets.size()>=9 &&
        met>200 && mt>175 && Ap>0.07 && met/sqrt(ht)>=8 && meff>1000)
      {
	if (meff<1500) accept("SRh9j1");
	else accept("SRh9j2");
	
	if (meff>1500) accept("SRh9j");
      }
    // discovery
    if (hardJets.size()>=4 && hardJets.size()<=5 &&
        met>300 && mt>450 && Ap>0.01 && met/meff>0.25 && meff>1000)
      {
	if (meff>1500) accept("SRh4jhx");
      }
    if (hardJets.size()>=4 && hardJets.size()<=5 &&
        met>250 && mt>150 && mt<450 && Ap>0.05 && meff>1300)
      {
	if (meff>1650) accept("SRh4jlxGG");
	if (meff>1300) accept("SRh4jlxSS");
      }
    if (hardJets.size()>=6 &&
        met>350 && mt>175 && Ap>0.06 && meff>700)
      {
	if (meff>2300) accept("SRh6jGG");
	if (meff>1233) accept("SRh6jSS");
      }
    
  }

  // Soft lepton overlap removal and signal lepton/jet definitions
  auto softJets            = overlapRemoval(preJets, softElectrons, 0.2, NOT(BTag77MV2c10));
       softElectrons       = overlapRemoval(softElectrons, softJets, 0.4);
       softElectrons       = overlapRemoval(softElectrons, softMuons, 0.01);
       softJets            = overlapRemoval(softJets, softMuons, 0.4, LessThan3Tracks);
       softMuons           = overlapRemoval(softMuons, softJets, 0.4);
       softJets            = filterObjects(softJets, 30, 2.8, JVT50Jet);
  auto softBJets           = filterObjects(softJets, 30, 2.5, BTag77MV2c10);
  auto softSignalElectrons = filterObjects(softElectrons, 7, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoGradientLoose);
  auto softSignalMuons     = filterObjects(softMuons, 6, 2.5, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);

  // Soft Lepton signal and control regions
  auto softLeptons=softElectrons+softMuons;
  auto softSignalLeptons=softSignalElectrons+softSignalMuons;
  if (softLeptons.size()==1 && softSignalLeptons.size()==1 && softSignalLeptons[0].Pt()< (35>(5*softJets.size())?5*(softJets.size()):35) ) {
    float mt   = calcMT(softSignalLeptons[0], metVec);
    float ht   = sumObjectsPt(softJets,999,30);
    float meff = ht + met + softSignalLeptons[0].Pt();
    //float Ap   = aplanarity(softJets);

    // b-Tag signal region
    if (softBJets.size()>0){
      if (softJets.size()>=2 && 
	  met>430 && mt>100 && met/meff>0.25 && meff>700)
	{
	  if (meff<1100) accept("SRs2j1bT");
	  else if (meff<1500) accept("SRs2j2bT");
	  else if (meff<1900) accept("SRs2j3bT");
	  else accept("SRs2j4bT");
	}
    }
    // b-Veto signal region
    else{
      if (softJets.size()>=2 && 
	  met>430 && mt>100 && met/meff>0.25 && meff>700)
	{
	  if (meff<1100) accept("SRs2j1bV");
	  else if (meff<1500) accept("SRs2j2bV");
	  else if (meff<1900) accept("SRs2j3bV");
	  else accept("SRs2j4bV");
	}
    }
    // b-inclusive discovery SR
    if (softJets.size()>=2 && 
	met>430 && mt>100 && met/meff>0.25 && meff>700)
      {
	if (meff>1100) accept("SRs2j");
      }
  }

  return;
}
