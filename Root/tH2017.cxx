#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

DefineAnalysis(tH2017)

void tH2017::Init()
{
  //addRegions({"SR1lep0b","SR1lep1b","SR1lep2b","SR1lep3b","SR1lep4b"});

  addHistogram("events",10,-0.5,9.5);
  
  //btagged regions
  for(int i=0;i<5;i++){
    std::string SR=std::to_string(i); 
    for(int j=0;j<5;j++){ //1st,2nd,3rd,4th most forward jets
      std::string jet=std::to_string(j);
      std::string etaHist="fwdJet"+jet+"Eta_SRB"+SR;
      std::string ptHist="fwdJet"+jet+"Pt_SRB"+SR;
      addHistogram(etaHist,25,0,4);
      addHistogram(ptHist,100,0,500);
    }
    addHistogram("fwdBJet1Eta_SRB"+SR,25,0,4);
    addHistogram("fwdBJet1Pt_SRB"+SR,100,0,500);
    
    addHistogram("lep_pt_SRB"+SR,100,0,500); 
    addHistogram("lep_eta_SRB"+SR,100,0,500); 
    
    addHistogram("leadJetPt_SRB"+SR,100,0,500); 
    addHistogram("leadJetEta_SRB"+SR,25,0,4); 
    addHistogram("leadBJetPt_SRB"+SR,100,0,500); 
    addHistogram("leadBJetEta_SRB"+SR,25,0,4); 
    
    addHistogram("Hbjets_m_SRB"+SR,100,0,500); //combined mass of 2 jets closest to mass of Higgs
   }

  //no cuts
  addHistogram("NLep_nocuts",10,-0.5,9.5);
  addHistogram("MET_nocuts",100,0,500);
  addHistogram("Njets_nocuts",10,-0.5,9.5);
  addHistogram("jet_pt_nocuts",100,0,500);
  addHistogram("lep_pt_nocuts",100,0,500); 
  
  //after preselection
  addHistogram("MET",100,0,500);
  addHistogram("Njets",10,-0.5,9.5);
  addHistogram("Nbjets",6,-0.5,5.5);
  addHistogram("sumAllPt",100,0,2000); 
 
  addHistogram("lep_pt",100,0,500);
  addHistogram("lep_eta",25,0,4);   

  addHistogram("Hbjets_m",100,0,500); //combined mass of 2 jets closest to mass of Higgs 
}

void tH2017::ProcessEvent(AnalysisEvent *event)
{
  _output->setEventWeight(1);
  fill("events",0); //Number of processed events

  double eventweight=event->getMCWeights()[0];
  _output->setEventWeight(eventweight); //all histograms after this point get filled with this weight
  fill("events",1); //Sum of initial weights


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
  if(leptons.size()>0) fill("lep_pt_nocuts",leptons.at(0).Pt());
  fill("MET_nocuts",met);
  fill("Njets_nocuts",numSignalJets);
  for(int iJet=0;iJet<numSignalJets;iJet++) fill("jet_pt_nocuts", jets.at(iJet).Pt()); 
  

  // Preselection
  if (numSignalLeptons != 1) return;
  if (numSignalJets < 2) return; 
  if (met < 20) return; 

  // Fill histogram after cuts
  fill("MET",met);
  fill("Njets",numSignalJets);
  fill("Nbjets",nBjets);
  fill("lep_pt",leptons.at(0).Pt());
  fill("lep_eta",leptons.at(0).Eta());
  
  //Sum pT of all objects
  double pTSum=leptons.at(0).Pt();
  for(int i=0;i<numSignalJets;i++) pTSum+=jets.at(i).Pt();
  fill("sumAllPt",pTSum); 

  //Find leading pT jet
  double highPt=-1; 
  double highPtJet_idx=-1; 
  for(int iJet=0;iJet<numSignalJets;iJet++){
    if(jets.at(iJet).Pt()>highPt){ 
      highPt=jets.at(iJet).Pt();
      highPtJet_idx=iJet;
    }
  }
  //Find leading pT bjet
  double highPtB=-1;
  double highPtBJet_idx=-1;
  for(int iJet=0;iJet<nBjets;iJet++){
    if(bjets.at(iJet).Pt()>highPtB){ 
      highPtB=bjets.at(iJet).Pt();
      highPtBJet_idx=iJet;
    }
  }
  
  // Find 1st,2nd,3rd most forward jets per event
  double fwd_eta[4]; for(int j=0;j<4;j++) fwd_eta[j]=-999; 
  int fwdJet_idx[4]; for(int j=0;j<4;j++) fwdJet_idx[j]=-1;
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(fabs(jets.at(iJet).Eta())>fwd_eta[0]){
      fwdJet_idx[0]=iJet;
      fwd_eta[0]=fabs(jets.at(iJet).Eta());
    }
  }
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(iJet==fwdJet_idx[0]) continue; 
    if(fabs(jets.at(iJet).Eta())>fwd_eta[1]){
      fwdJet_idx[1]=iJet;
      fwd_eta[1]=fabs(jets.at(iJet).Eta());
    }
  }
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(iJet==fwdJet_idx[0]||iJet==fwdJet_idx[1]) continue; 
    if(fabs(jets.at(iJet).Eta())>fwd_eta[2]){
      fwdJet_idx[2]=iJet;
      fwd_eta[2]=fabs(jets.at(iJet).Eta());
    }
  }
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(iJet==fwdJet_idx[0]||iJet==fwdJet_idx[1]||iJet==fwdJet_idx[2]) continue; 
    if(fabs(jets.at(iJet).Eta())>fwd_eta[3]){
      fwdJet_idx[3]=iJet;
      fwd_eta[3]=fabs(jets.at(iJet).Eta());
    }
  }

  //Find most forward bjet per event
  double fwdBJet_eta=-999;
  int fwdBJet_idx=-1;
  for(int iJet=0;iJet<nBjets;iJet++){ 
    if(fabs(bjets.at(iJet).Eta())>fwdBJet_eta){
      fwdBJet_eta=fabs(bjets.at(iJet).Eta());
      fwdBJet_idx=iJet;
    }
  }

  //Find 2 jets whose combined mass is closest to the Higgs
  double mHiggs=125; 
  int jet1=-1;
  int jet2=-1;
  double mDiff=1000; 
  for(int J1=0;J1<numSignalJets;J1++){
    for(int J2=0;J2<numSignalJets;J2++){
      if(J2==J1) continue; 
      double mH=(jets.at(J1)+jets.at(J2)).M();
      if(fabs(mH-mHiggs)<mDiff){
  	mDiff=fabs(mH-mHiggs);
  	jet1=J1;
  	jet2=J2;
      }
    }
  }
  double mCombined=-1;
  if(jet1!=-1&&jet2!=-1) mCombined=(jets.at(jet1)+jets.at(jet2)).M();
  fill("Hbjets_m",mCombined);
  

  //fill histos corresponding to b-Tag signal regions
  if(nBjets>4) return; 
  std::string SR=std::to_string(nBjets); 
  
  for(int i=0;i<4;i++){//first 4 forward jets
    std::string jetPosition=std::to_string(i);
    fill("fwdJet"+jetPosition+"Eta_SRB"+SR, fwd_eta[i]);
    if(fwdJet_idx[i]!=-1) fill("fwdJet"+jetPosition+"Pt_SRB"+SR, jets.at(fwdJet_idx[i]).Pt());
  }

  fill("fwdBJet1Eta_SRB"+SR, fwdBJet_eta); 
  if(fwdBJet_idx!=-1) fill("fwdBJet1Pt_SRB"+SR, bjets.at(fwdBJet_idx).Pt()); 
  
  fill("lep_pt_SRB"+SR,leptons.at(0).Pt());
  fill("lep_eta_SRB"+SR,leptons.at(0).Eta()); 
  fill("Hbjets_m_SRB"+SR,mCombined); 
  fill("leadJetPt_SRB"+SR, jets.at(highPtJet_idx).Pt());  
  fill("leadJetEta_SRB"+SR, jets.at(highPtJet_idx).Eta());
  
  if(nBjets>0){
    fill("leadBJetPt_SRB"+SR, bjets.at(highPtBJet_idx).Pt());
    fill("leadBJetEta_SRB"+SR, bjets.at(highPtBJet_idx).Eta());
  }

  // // b-Tag signal regions
  // if (nBjets==0){
  //   fill("fwdJet1Eta_SRB0", fwd_eta[0]);
  //   fill("fwdJet2Eta_SRB0", fwd_eta[1]);
  //   fill("fwdJet3Eta_SRB0", fwd_eta[2]);
  //   fill("fwdJet4Eta_SRB0", fwd_eta[3]);

  //   fill("fwdBJet1Eta_SRB0", fwdBJet_eta); 
  // }
  // if (nBjets==1){
  //   fill("fwdJet1Eta_SRB1", fwd_eta[0]);
  //   fill("fwdJet2Eta_SRB1", fwd_eta[1]);
  //   fill("fwdJet3Eta_SRB1", fwd_eta[2]);
  //   fill("fwdJet4Eta_SRB1", fwd_eta[3]);
    
  //   fill("fwdBJet1Eta_SRB1", fwdBJet_eta); 
  // }
  // if (nBjets==2){
  //   fill("fwdJet1Eta_SRB2", fwd_eta[0]);
  //   fill("fwdJet2Eta_SRB2", fwd_eta[1]);
  //   fill("fwdJet3Eta_SRB2", fwd_eta[2]);
  //   fill("fwdJet4Eta_SRB2", fwd_eta[3]);

  //   fill("fwdBJet1Eta_SRB2", fwdBJet_eta); 
  // }
  // if (nBjets==3){
  //   fill("fwdJet1Eta_SRB3", fwd_eta[0]);
  //   fill("fwdJet2Eta_SRB3", fwd_eta[1]);
  //   fill("fwdJet3Eta_SRB3", fwd_eta[2]);
  //   fill("fwdJet4Eta_SRB3", fwd_eta[3]);
      
  //   fill("fwdBJet1Eta_SRB3", fwdBJet_eta); 
  // }
  // if (nBjets==4){
  //   fill("fwdJet1Eta_SRB4", fwd_eta[0]);
  //   fill("fwdJet2Eta_SRB4", fwd_eta[1]);
  //   fill("fwdJet3Eta_SRB4", fwd_eta[2]);
  //   fill("fwdJet4Eta_SRB4", fwd_eta[3]);

  //   fill("fwdBJet1Eta_SRB4", fwdBJet_eta); 
  // }


  return;
}




  // //Find 1st,2nd,3rd most forward bjets
  // double fwdBJet_eta[4]; for(int j=0;j<4;j++) fwdBJet_eta[j]=-999;  
  // int fwdBJet_idx[4]; for(int j=0;j<4;j++) fwdBJet_idx[j]=-1;
  // for(int iJet=0;iJet<nBjets;iJet++){ 
  //   if(fabs(bjets.at(iJet).Eta())>fwdBJet_eta[0]){
  //     fwdBJet_idx[0]=iJet;
  //     fwdBJet_eta[0]=fabs(bjets.at(iJet).Eta());
  //   }
  // }
  // for(int iJet=0;iJet<nBjets;iJet++){ 
  //   if(iJet==fwdBJet_idx[0]) continue; 
  //   if(fabs(bjets.at(iJet).Eta())>fwdBJet_eta[1]){
  //     fwdBJet_idx[1]=iJet;
  //     fwdBJet_eta[1]=fabs(bjets.at(iJet).Eta());
  //   }
  // }
  // for(int iJet=0;iJet<nBjets;iJet++){ 
  //   if(iJet==fwdBJet_idx[0]||iJet==fwdBJet_idx[1]) continue; 
  //   if(fabs(bjets.at(iJet).Eta())>fwdBJet_eta[2]){
  //     fwdBJet_idx[2]=iJet;
  //     fwdBJet_eta[2]=fabs(bjets.at(iJet).Eta());
  //   }
  // }

 
   // fill("1stfwdBJetEta4", fwdBJet_eta[0]);
    // fill("2ndfwdBJetEta4", fwdBJet_eta[1]);
    // fill("3rdfwdBJetEta4", fwdBJet_eta[2]);

  //Fill in optional ntuple variables
  // ntupVar("met",met);    
  // ntupVar("nBjets",nBjets); 
  // ntupVar("ele", signalElectrons);
  // ntupVar("jet", jets); 
  // ntupVar("bjet", bjets); 
  // if(fwd_eta[0]!=0) ntupVar("fwdJet", fwdJet); 
