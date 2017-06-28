#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(tH2017)

void tH2017::Init()
{
  addRegions({"SR1lep0b","SR1lep1b","SR1lep2b","SR1lep3b","SR1lep4b"});
  
  //the most forward jets for each signal region
  addHistogram("1stfwdJetEta0", 25,0,4.); 
  addHistogram("1stfwdJetEta1", 25,0,4.); 
  addHistogram("1stfwdJetEta2", 25,0,4.); 
  addHistogram("1stfwdJetEta3", 25,0,4.); 
  addHistogram("1stfwdJetEta4", 25,0,4.);

  addHistogram("2ndfwdJetEta0", 25,0,4.); 
  addHistogram("2ndfwdJetEta1", 25,0,4.); 
  addHistogram("2ndfwdJetEta2", 25,0,4.); 
  addHistogram("2ndfwdJetEta3", 25,0,4.); 
  addHistogram("2ndfwdJetEta4", 25,0,4.);

  addHistogram("3rdfwdJetEta0", 25,0,4.); 
  addHistogram("3rdfwdJetEta1", 25,0,4.); 
  addHistogram("3rdfwdJetEta2", 25,0,4.); 
  addHistogram("3rdfwdJetEta3", 25,0,4.); 
  addHistogram("3rdfwdJetEta4", 25,0,4.);

  //the most forward bjets for each signal region
  addHistogram("1stfwdBJetEta0", 25,0,4.); 
  addHistogram("1stfwdBJetEta1", 25,0,4.); 
  addHistogram("1stfwdBJetEta2", 25,0,4.); 
  addHistogram("1stfwdBJetEta3", 25,0,4.); 
  addHistogram("1stfwdBJetEta4", 25,0,4.);

  addHistogram("2ndfwdBJetEta0", 25,0,4.); 
  addHistogram("2ndfwdBJetEta1", 25,0,4.); 
  addHistogram("2ndfwdBJetEta2", 25,0,4.); 
  addHistogram("2ndfwdBJetEta3", 25,0,4.); 
  addHistogram("2ndfwdBJetEta4", 25,0,4.);

  addHistogram("3rdfwdBJetEta0", 25,0,4.); 
  addHistogram("3rdfwdBJetEta1", 25,0,4.); 
  addHistogram("3rdfwdBJetEta2", 25,0,4.); 
  addHistogram("3rdfwdBJetEta3", 25,0,4.); 
  addHistogram("3rdfwdBJetEta4", 25,0,4.);

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

  // Find 1st,2nd,3rd most forward jets
  double fwd_eta[3]; for(int j=0;j<3;j++) fwd_eta[j]=-999; 
  int fwdJet_idx[3]; for(int j=0;j<3;j++) fwdJet_idx[j]=-1;
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(fabs(signalJets[iJet].Eta())>fwd_eta[0]){
      fwdJet_idx[0]=iJet;
      fwd_eta[0]=fabs(signalJets[iJet].Eta());
    }
  }
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(iJet==fwdJet_idx[0]) continue; 
    if(fabs(signalJets[iJet].Eta())>fwd_eta[1]){
      fwdJet_idx[1]=iJet;
      fwd_eta[1]=fabs(signalJets[iJet].Eta());
    }
  }
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(iJet==fwdJet_idx[0]||iJet==fwdJet_idx[1]) continue; 
    if(fabs(signalJets[iJet].Eta())>fwd_eta[2]){
      fwdJet_idx[2]=iJet;
      fwd_eta[2]=fabs(signalJets[iJet].Eta());
    }
  }

  //Find 1st,2nd,3rd most forward bjets
  double fwd_bjet_eta[3]; for(int j=0;j<3;j++) fwd_bjet_eta[j]=-999;  
  int fwdBJet_idx[3]; for(int j=0;j<3;j++) fwdBJet_idx[j]=-1;
  for(int iJet=0;iJet<nBjets;iJet++){ 
    if(fabs(bjets[iJet].Eta())>fwd_bjet_eta[0]){
      fwdBJet_idx[0]=iJet;
      fwd_bjet_eta[0]=fabs(bjets[iJet].Eta());
    }
  }
  for(int iJet=0;iJet<nBjets;iJet++){ 
    if(iJet==fwdBJet_idx[0]) continue; 
    if(fabs(bjets[iJet].Eta())>fwd_bjet_eta[1]){
      fwdBJet_idx[1]=iJet;
      fwd_bjet_eta[1]=fabs(bjets[iJet].Eta());
    }
  }
  for(int iJet=0;iJet<nBjets;iJet++){ 
    if(iJet==fwdBJet_idx[0]||iJet==fwdBJet_idx[1]) continue; 
    if(fabs(bjets[iJet].Eta())>fwd_bjet_eta[2]){
      fwdBJet_idx[2]=iJet;
      fwd_bjet_eta[2]=fabs(bjets[iJet].Eta());
    }
  }

  // Flag signal regions
  if (numSignalLeptons==1) {
 
    // b-Tag signal regions
    if (nBjets==0){
      //accept("SR1lep0b");
      fill("1stfwdJetEta0", fwd_eta[0]);
      fill("2ndfwdJetEta0", fwd_eta[1]);
      fill("3rdfwdJetEta0", fwd_eta[2]);

      fill("1stfwdBJetEta0", fwd_bjet_eta[0]);
      fill("2ndfwdBJetEta0", fwd_bjet_eta[1]);
      fill("3rdfwdBJetEta0", fwd_bjet_eta[2]);
    }
    if (nBjets==1){
      //accept("SR1lep1b");
      fill("1stfwdJetEta1", fwd_eta[0]);
      fill("2ndfwdJetEta1", fwd_eta[1]);
      fill("3rdfwdJetEta1", fwd_eta[2]);

      fill("1stfwdBJetEta1", fwd_bjet_eta[0]);
      fill("2ndfwdBJetEta1", fwd_bjet_eta[1]);
      fill("3rdfwdBJetEta1", fwd_bjet_eta[2]);
    }
    if (nBjets==2){
      //accept("SR1lep2b");
      fill("1stfwdJetEta2", fwd_eta[0]);
      fill("2ndfwdJetEta2", fwd_eta[1]);
      fill("3rdfwdJetEta2", fwd_eta[2]);

      fill("1stfwdBJetEta2", fwd_bjet_eta[0]);
      fill("2ndfwdBJetEta2", fwd_bjet_eta[1]);
      fill("3rdfwdBJetEta2", fwd_bjet_eta[2]);
    }
    if (nBjets==3){
      //accept("SR1lep3b");
      fill("1stfwdJetEta3", fwd_eta[0]);
      fill("2ndfwdJetEta3", fwd_eta[1]);
      fill("3rdfwdJetEta3", fwd_eta[2]);

      fill("1stfwdBJetEta3", fwd_bjet_eta[0]);
      fill("2ndfwdBJetEta3", fwd_bjet_eta[1]);
      fill("3rdfwdBJetEta3", fwd_bjet_eta[2]);
    }
    if (nBjets==4){
      //accept("SR1lep4b");
      fill("1stfwdJetEta4", fwd_eta[0]);
      fill("2ndfwdJetEta4", fwd_eta[1]);
      fill("3rdfwdJetEta4", fwd_eta[2]);

      fill("1stfwdBJetEta4", fwd_bjet_eta[0]);
      fill("2ndfwdBJetEta4", fwd_bjet_eta[1]);
      fill("3rdfwdBJetEta4", fwd_bjet_eta[2]);
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
