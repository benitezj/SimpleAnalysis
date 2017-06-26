#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(tHAnalysis)

void tHAnalysis::Init()
{
  // Define signal/control regions
  addRegions({"SR1","SR2"});   // Use same name as in paper, but avoid "-" 
  addRegion("CR_B");            // and similar non-alphanumeric characters

  // Book 1/2D histograms
  addHistogram("MET_nocuts",100,0,1000);
  addHistogram("Njets_nocuts",10,-0.5,9.5);
  addHistogram("Nbjets_nocuts",6,-0.5,5.5);
  
  addHistogram("MET",100,0,1000);
  addHistogram("Njets",10,-0.5,9.5);
  addHistogram("Nbjets",6,-0.5,5.5);
  addHistogram("pTL",100,0,500);
}

void tHAnalysis::ProcessEvent(AnalysisEvent *event)
{
  // Retrieve basic object lists
  // PLEASE NOTE UNITS ARE ALL IN GEV, AND NOT ATLAS STANDARD MEV!!!
  auto electrons  = event->getElectrons(20, 2.47, ELooseLH); // Filter on pT, eta and "ID"
  auto muons      = event->getMuons(20, 2.5, MuMedium);
  auto candJets   = event->getJets(25., 4.5);
  auto metVec     = event->getMET();
  double met      = metVec.Et();

  // Overlap removal - including with object Pt-dependent radius calculation
  electrons  = overlapRemoval(electrons, muons, 0.01);
  candJets   = overlapRemoval(candJets, electrons, 0.2, NOT(BTag85MV2c20));
  electrons  = overlapRemoval(electrons, candJets, 0.4);
  candJets   = overlapRemoval(candJets, muons, 0.4, LessThan3Tracks); 
  muons      = overlapRemoval(muons, candJets, 0.4);                  

  // Basic filtering by pT, eta and "ID"
  auto signalJets      = filterObjects(candJets, 25);
  auto signalElectrons = filterObjects(electrons, 20, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoBoosted);
  auto signalMuons     = filterObjects(muons, 20, 2.5, MuD0Sigma3|MuZ05mm|MuIsoBoosted|MuNotCosmic);
  auto bjets = filterObjects(signalJets, 25., 2.5, BTag85MV2c20);

  // Lists of objects can be merget by simple addition
  auto signalLeptons   = signalElectrons + signalMuons;

  // Object counting
  int numSignalLeptons = signalLeptons.size();               // Object lists are essentially std::vectors so .size() works
  int numSignalJets    = signalJets.size();
  int nBjets = bjets.size();

  // Fill in histograms without cuts
  fill("MET_nocuts",met);
  fill("Njets_nocuts",numSignalJets);
  fill("Nbjets_nocuts",nBjets);

  // Preselection
  if (numSignalLeptons != 1) return;

  // Fill histogram after cuts
  fill("MET",met);
  fill("Njets",numSignalJets);
  fill("Nbjets",nBjets);
  fill("pTL",signalLeptons.at(0).Pt());

  // Flag signal regions that pass selections 
  //if (met>200 && mT>100 && nBjets>2) accept("SR1");
   
  //Fill in optional ntuple variables
  ntupVar("met",met);               // can be simple variables
  ntupVar("nBjets",nBjets);
  ntupVar("signalJets",signalJets); // or even a list of objects (or single object)

  return;
}
