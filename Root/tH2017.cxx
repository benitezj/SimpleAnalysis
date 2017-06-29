#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(tH2017)

void tH2017::Init()
{
  //addRegions({"SR1lep0b","SR1lep1b","SR1lep2b","SR1lep3b","SR1lep4b"});
  
  //the most forward jet per event for each signal region
  addHistogram("fwdJet1Eta_SRB0", 25,0,4.); 
  addHistogram("fwdJet1Eta_SRB1", 25,0,4.); 
  addHistogram("fwdJet1Eta_SRB2", 25,0,4.); 
  addHistogram("fwdJet1Eta_SRB3", 25,0,4.); 
  addHistogram("fwdJet1Eta_SRB4", 25,0,4.);

  addHistogram("fwdJet2Eta_SRB0", 25,0,4.); 
  addHistogram("fwdJet2Eta_SRB1", 25,0,4.); 
  addHistogram("fwdJet2Eta_SRB2", 25,0,4.); 
  addHistogram("fwdJet2Eta_SRB3", 25,0,4.); 
  addHistogram("fwdJet2Eta_SRB4", 25,0,4.);

  addHistogram("fwdJet3Eta_SRB0", 25,0,4.); 
  addHistogram("fwdJet3Eta_SRB1", 25,0,4.); 
  addHistogram("fwdJet3Eta_SRB2", 25,0,4.); 
  addHistogram("fwdJet3Eta_SRB3", 25,0,4.); 
  addHistogram("fwdJet3Eta_SRB4", 25,0,4.);

  addHistogram("fwdJet4Eta_SRB0", 25,0,4.); 
  addHistogram("fwdJet4Eta_SRB1", 25,0,4.); 
  addHistogram("fwdJet4Eta_SRB2", 25,0,4.); 
  addHistogram("fwdJet4Eta_SRB3", 25,0,4.); 
  addHistogram("fwdJet4Eta_SRB4", 25,0,4.);

  // //the most forward bjets for each signal region
  // addHistogram("1stfwdBJetEta0", 25,0,4.); 
  // addHistogram("1stfwdBJetEta1", 25,0,4.); 
  // addHistogram("1stfwdBJetEta2", 25,0,4.); 
  // addHistogram("1stfwdBJetEta3", 25,0,4.); 
  // addHistogram("1stfwdBJetEta4", 25,0,4.);

  // addHistogram("2ndfwdBJetEta0", 25,0,4.); 
  // addHistogram("2ndfwdBJetEta1", 25,0,4.); 
  // addHistogram("2ndfwdBJetEta2", 25,0,4.); 
  // addHistogram("2ndfwdBJetEta3", 25,0,4.); 
  // addHistogram("2ndfwdBJetEta4", 25,0,4.);

  // addHistogram("3rdfwdBJetEta0", 25,0,4.); 
  // addHistogram("3rdfwdBJetEta1", 25,0,4.); 
  // addHistogram("3rdfwdBJetEta2", 25,0,4.); 
  // addHistogram("3rdfwdBJetEta3", 25,0,4.); 
  // addHistogram("3rdfwdBJetEta4", 25,0,4.);

  // addHistogram("4thfwdBJetEta0", 25,0,4.); 
  // addHistogram("4thfwdBJetEta1", 25,0,4.); 
  // addHistogram("4thfwdBJetEta2", 25,0,4.); 
  // addHistogram("4thfwdBJetEta3", 25,0,4.); 
  // addHistogram("4thfwdBJetEta4", 25,0,4.);


  addHistogram("NLep_nocuts",10,-0.5,9.5);
  addHistogram("MET_nocuts",100,0,1000);
  addHistogram("Njets_nocuts",10,-0.5,9.5);
  
  addHistogram("MET",100,0,1000);
  addHistogram("Njets",10,-0.5,9.5);
  addHistogram("Nbjets",6,-0.5,5.5);
  addHistogram("pTL",100,0,500);
  addHistogram("fwdBJet1Eta", 25,0,4.); 
}

void tH2017::ProcessEvent(AnalysisEvent *event)
{
  auto electrons  = event->getElectrons(25., 2.47,ELooseBLLH && EIsoGradient); //add vertex  (ED0Sigma5|EZ05mm ?)
  auto muons      = event->getMuons(25., 2.7, MuLoose && MuIsoGradient);//add vertex (MuD0Sigma3|MuZ05mm ?)
  auto jets   = event->getJets(25., 3.8, JVT50Jet); // what about NOT(LooseBadJet)
  auto metVec     = event->getMET();
  double met = metVec.Et();
  

  // Overlap removal - including with object Pt-dependent radius calculation
  electrons  = overlapRemoval(electrons, muons, 0.01);
  jets   = overlapRemoval(jets, electrons, 0.2, NOT(BTag85MV2c20));//check this cut
  electrons  = overlapRemoval(electrons, jets, 0.4);
  jets   = overlapRemoval(jets, muons, 0.4, LessThan3Tracks); //check this cut (b-jets?)
  muons      = overlapRemoval(muons, jets, 0.4);   

  // Lists of objects can be merged by simple addition
  auto leptons   = electrons + muons;

  //
  auto bjets = filterObjects(jets, 25., 3.8, BTag85MV2c20);


  // Object counting
  int numSignalLeptons = leptons.size();  // Object lists are essentially std::vectors so .size() works
  int numSignalJets    = jets.size();
  int nBjets = bjets.size();
  
  // Fill in histograms without cuts
  fill("NLep_nocuts",numSignalLeptons);
  fill("MET_nocuts",met);
  fill("Njets_nocuts",numSignalJets);

  // Preselection
  if (numSignalLeptons != 1) return;
  if (numSignalJets < 2) return; 
  if (met < 20) return; 

  // Fill histogram after cuts
  fill("MET",met);
  fill("Njets",numSignalJets);
  fill("Nbjets",nBjets);
  fill("pTL",leptons.at(0).Pt());


  
  // Find 1st,2nd,3rd most forward jets
  double fwd_eta[4]; for(int j=0;j<4;j++) fwd_eta[j]=-999; 
  int fwdJet_idx[4]; for(int j=0;j<4;j++) fwdJet_idx[j]=-1;
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(fabs(jets[iJet].Eta())>fwd_eta[0]){
      fwdJet_idx[0]=iJet;
      fwd_eta[0]=fabs(jets[iJet].Eta());
    }
  }
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(iJet==fwdJet_idx[0]) continue; 
    if(fabs(jets[iJet].Eta())>fwd_eta[1]){
      fwdJet_idx[1]=iJet;
      fwd_eta[1]=fabs(jets[iJet].Eta());
    }
  }
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(iJet==fwdJet_idx[0]||iJet==fwdJet_idx[1]) continue; 
    if(fabs(jets[iJet].Eta())>fwd_eta[2]){
      fwdJet_idx[2]=iJet;
      fwd_eta[2]=fabs(jets[iJet].Eta());
    }
  }
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(iJet==fwdJet_idx[0]||iJet==fwdJet_idx[1]||iJet==fwdJet_idx[2]) continue; 
    if(fabs(jets[iJet].Eta())>fwd_eta[3]){
      fwdJet_idx[3]=iJet;
      fwd_eta[3]=fabs(jets[iJet].Eta());
    }
  }

  //Find most forward bjet per event
  double fwd_bjet_eta=-999;
  for(int iJet=0;iJet<nBjets;iJet++){ 
    if(fabs(bjets[iJet].Eta())>fwd_bjet_eta){
      fwd_bjet_eta=fabs(bjets[iJet].Eta());
    }
  }
  fill("fwdBJet1Eta", fwd_bjet_eta); 

  // b-Tag signal regions
  if (nBjets==0){
    fill("fwdJet1Eta_SRB0", fwd_eta[0]);
    fill("fwdJet2Eta_SRB0", fwd_eta[1]);
    fill("fwdJet3Eta_SRB0", fwd_eta[2]);
    fill("fwdJet4Eta_SRB0", fwd_eta[3]);
  }
  if (nBjets==1){
    fill("fwdJet1Eta_SRB1", fwd_eta[0]);
    fill("fwdJet2Eta_SRB1", fwd_eta[1]);
    fill("fwdJet3Eta_SRB1", fwd_eta[2]);
    fill("fwdJet4Eta_SRB1", fwd_eta[3]);

  }
  if (nBjets==2){
    fill("fwdJet1Eta_SRB2", fwd_eta[0]);
    fill("fwdJet2Eta_SRB2", fwd_eta[1]);
    fill("fwdJet3Eta_SRB2", fwd_eta[2]);
    fill("fwdJet4Eta_SRB2", fwd_eta[3]);

  }
  if (nBjets==3){
    fill("fwdJet1Eta_SRB3", fwd_eta[0]);
    fill("fwdJet2Eta_SRB3", fwd_eta[1]);
    fill("fwdJet3Eta_SRB3", fwd_eta[2]);
    fill("fwdJet4Eta_SRB3", fwd_eta[3]);
      
  }
  if (nBjets==4){
    fill("fwdJet1Eta_SRB4", fwd_eta[0]);
    fill("fwdJet2Eta_SRB4", fwd_eta[1]);
    fill("fwdJet3Eta_SRB4", fwd_eta[2]);
    fill("fwdJet4Eta_SRB4", fwd_eta[3]);

  }


  return;
}




  // //Find 1st,2nd,3rd most forward bjets
  // double fwd_bjet_eta[4]; for(int j=0;j<4;j++) fwd_bjet_eta[j]=-999;  
  // int fwdBJet_idx[4]; for(int j=0;j<4;j++) fwdBJet_idx[j]=-1;
  // for(int iJet=0;iJet<nBjets;iJet++){ 
  //   if(fabs(bjets[iJet].Eta())>fwd_bjet_eta[0]){
  //     fwdBJet_idx[0]=iJet;
  //     fwd_bjet_eta[0]=fabs(bjets[iJet].Eta());
  //   }
  // }
  // for(int iJet=0;iJet<nBjets;iJet++){ 
  //   if(iJet==fwdBJet_idx[0]) continue; 
  //   if(fabs(bjets[iJet].Eta())>fwd_bjet_eta[1]){
  //     fwdBJet_idx[1]=iJet;
  //     fwd_bjet_eta[1]=fabs(bjets[iJet].Eta());
  //   }
  // }
  // for(int iJet=0;iJet<nBjets;iJet++){ 
  //   if(iJet==fwdBJet_idx[0]||iJet==fwdBJet_idx[1]) continue; 
  //   if(fabs(bjets[iJet].Eta())>fwd_bjet_eta[2]){
  //     fwdBJet_idx[2]=iJet;
  //     fwd_bjet_eta[2]=fabs(bjets[iJet].Eta());
  //   }
  // }

 
   // fill("1stfwdBJetEta4", fwd_bjet_eta[0]);
    // fill("2ndfwdBJetEta4", fwd_bjet_eta[1]);
    // fill("3rdfwdBJetEta4", fwd_bjet_eta[2]);

  //Fill in optional ntuple variables
  // ntupVar("met",met);    
  // ntupVar("nBjets",nBjets); 
  // ntupVar("ele", signalElectrons);
  // ntupVar("jet", jets); 
  // ntupVar("bjet", bjets); 
  // if(fwd_eta[0]!=0) ntupVar("fwdJet", fwdJet); 
