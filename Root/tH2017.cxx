#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(tH2017)

void tH2017::Init()
{
  addRegions({"SR1lep0b","SR1lep1b","SR1lep2b","SR1lep3b","SR1lep4b"});
  
  addHistogram("fwdJetEta0", 25,0,4.); 
  addHistogram("fwdJetEta1", 25,0,4.); 
  addHistogram("fwdJetEta2", 25,0,4.); 
  addHistogram("fwdJetEta3", 25,0,4.); 
  addHistogram("fwdJetEta4", 25,0,4.);

  addHistogram("fwdBJetEta0", 25,0,4.); 
  addHistogram("fwdBJetEta1", 25,0,4.); 
  addHistogram("fwdBJetEta2", 25,0,4.); 
  addHistogram("fwdBJetEta3", 25,0,4.); 
  addHistogram("fwdBJetEta4", 25,0,4.);

  addHistogram("MET_nocuts",100,0,1000);
  addHistogram("Njets_nocuts",10,-0.5,9.5);
  addHistogram("Nbjets_nocuts",6,-0.5,5.5);
  
  addHistogram("MET",100,0,1000);
  addHistogram("Njets",10,-0.5,9.5);
  addHistogram("Nbjets",6,-0.5,5.5);
  addHistogram("pTL",100,0,500);
}

void tH2017::ProcessEvent(AnalysisEvent *event)
{
  auto electrons  = event->getElectrons(25., 2.47, ELooseLH);
  auto muons      = event->getMuons(25., 2.7, MuIsoLoose);
  auto candJets   = event->getJets(25., 3.8, JVT50Jet);
  auto metVec     = event->getMET();
  double met = metVec.Et();
  
  // Choose events with exactly 1 lepton, MET>20GeV, at least 2 jets
  if(candJets.size()<2||met<20) return; 
  //if (countObjects(candJets, 20, 4.5, NOT(LooseBadJet))!=0) return;

  // Overlap removal - including with object Pt-dependent radius calculation
  electrons  = overlapRemoval(electrons, muons, 0.01);
  candJets   = overlapRemoval(candJets, electrons, 0.2, NOT(BTag85MV2c20));
  electrons  = overlapRemoval(electrons, candJets, 0.4);
  candJets   = overlapRemoval(candJets, muons, 0.4, LessThan3Tracks); 
  muons      = overlapRemoval(muons, candJets, 0.4);   

  // Basic filtering by pT, eta and "ID"
  auto signalJets      = filterObjects(candJets, 25.);
  auto signalElectrons = filterObjects(electrons, 20., 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoBoosted);
  auto signalMuons     = filterObjects(muons, 20., 2.5, MuD0Sigma3|MuZ05mm|MuIsoBoosted|MuNotCosmic);
  auto bjets = filterObjects(signalJets, 25., 3.8, BTag85MV2c20);

  // Lists of objects can be merged by simple addition
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
  if (met < 20) return; 

  // Fill histogram after cuts
  fill("MET",met);
  fill("Njets",numSignalJets);
  fill("Nbjets",nBjets);
  fill("pTL",signalLeptons.at(0).Pt());

  // Find forward jets/bjets
  // Find forward jets  
  double fwd_eta[3]; for(int j=0;j<3;j++) fwd_eta[j]=-999; 
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(fabs(signalJets[iJet].Eta())>fwd_eta[0]){
      fwd_eta[0]=fabs(signalJets[iJet].Eta());
    }
  }

  double fwd_bjet_eta[3]; for(int j=0;j<3;j++) fwd_bjet_eta[j]=-999;  
  for(int iJet=0;iJet<nBjets;iJet++){ 
    if(fabs(bjets[iJet].Eta())>fwd_bjet_eta[0]){
      fwd_bjet_eta[0]=fabs(bjets[iJet].Eta());
    }
  }

  // Flag signal regions
  if (numSignalLeptons==1) {
 
    // b-Tag signal regions
    if (nBjets==0){
      //accept("SR1lep0b");
      fill("fwdJetEta0", fwd_eta[0]);
      fill("fwdBJetEta0", fwd_bjet_eta[0]);
    }
    if (nBjets==1){
      //accept("SR1lep1b");
      fill("fwdJetEta1", fwd_eta[0]);
      fill("fwdBJetEta1", fwd_bjet_eta[0]);
    }
    if (nBjets==2){
      //accept("SR1lep2b");
      fill("fwdJetEta2", fwd_eta[0]);
      fill("fwdBJetEta2", fwd_bjet_eta[0]);
    }
    if (nBjets==3){
      //accept("SR1lep3b");
      fill("fwdJetEta3", fwd_eta[0]);
      fill("fwdBJetEta3", fwd_bjet_eta[0]);
    }
    if (nBjets==4){
      //accept("SR1lep4b");
      fill("fwdJetEta4", fwd_eta[0]);
      fill("fwdBJetEta4", fwd_bjet_eta[0]);
    }
  }

  //Fill in optional ntuple variables
  // ntupVar("met",met);    
  // ntupVar("nBjets",nBjets); 
  // ntupVar("ele", signalElectrons);
  // ntupVar("jet", signalJets); 
  // ntupVar("bjet", bjets); 
  // if(fwd_eta[0]!=0) ntupVar("fwdJet", fwdJet); 
  return;
}
