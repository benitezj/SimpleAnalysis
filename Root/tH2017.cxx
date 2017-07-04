#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

double dR_fn (float eta1, float eta2, float phi1, float phi2){
  double deta= fabs(eta1 - eta2);      double dphi= fabs(phi1 - phi2);
  if (dphi > 3.14 ) dphi = 2*3.14 - dphi;
  return sqrt((dphi*dphi)+(deta*deta));
}

std::vector<TLorentzVector> findForwardJets(auto jetVec){
  std::vector<TLorentzVector> forwardJets; 
  double fwd_eta[4]; for(int j=0;j<4;j++) fwd_eta[j]=-999; 
  int fwdJet_idx[4]; for(int j=0;j<4;j++) fwdJet_idx[j]=-1;
  int numSignalJets=jetVec.size(); 
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(fabs(jetVec.at(iJet).Eta())>fwd_eta[0]){
      fwdJet_idx[0]=iJet;
      fwd_eta[0]=fabs(jetVec.at(iJet).Eta());
    }
  }
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(iJet==fwdJet_idx[0]) continue; 
    if(fabs(jetVec.at(iJet).Eta())>fwd_eta[1]){
      fwdJet_idx[1]=iJet;
      fwd_eta[1]=fabs(jetVec.at(iJet).Eta());
    }
  }
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(iJet==fwdJet_idx[0]||iJet==fwdJet_idx[1]) continue; 
    if(fabs(jetVec.at(iJet).Eta())>fwd_eta[2]){
      fwdJet_idx[2]=iJet;
      fwd_eta[2]=fabs(jetVec.at(iJet).Eta());
    }
  }
  for(int iJet=0;iJet<numSignalJets;iJet++){ 
    if(iJet==fwdJet_idx[0]||iJet==fwdJet_idx[1]||iJet==fwdJet_idx[2]) continue; 
    if(fabs(jetVec.at(iJet).Eta())>fwd_eta[3]){
      fwdJet_idx[3]=iJet;
      fwd_eta[3]=fabs(jetVec.at(iJet).Eta());
    }
  }
  for(int i=0;i<4;i++) if(fwdJet_idx[i]!=-1) forwardJets.push_back(jetVec.at(fwdJet_idx[i]));
  return forwardJets; 
}

std::vector<TLorentzVector> findLeadingJets(auto jetVec){
  std::vector<TLorentzVector> leadingJets;
  double highPt[2]={-1,-1}; 
  double highPtJet_idx[2]={-1,-1}; 
  for(int iJet=0;iJet<jetVec.size();iJet++){
    if(jetVec.at(iJet).Pt()>highPt[0]){ 
      highPt[0]=jetVec.at(iJet).Pt();
      highPtJet_idx[0]=iJet;
    }
  }
  for(int iJet=0;iJet<jetVec.size();iJet++){
    if(iJet==highPtJet_idx[0]) continue; 
    if(jetVec.at(iJet).Pt()>highPt[1]){ 
      highPt[1]=jetVec.at(iJet).Pt();
      highPtJet_idx[1]=iJet;
    }
  }
  if(highPtJet_idx[0]!=-1) leadingJets.push_back(jetVec.at(highPtJet_idx[0]));
  if(highPtJet_idx[1]!=-1) leadingJets.push_back(jetVec.at(highPtJet_idx[1]));
  return leadingJets; 
}

TLorentzVector findHiggs_method1(auto bjets){
  TLorentzVector higgs1;
  double mHiggs=125; 
  int nBjets = bjets.size();
  int bjet1=-1;
  int bjet2=-1;
  double mDiff=1000; 
  for(int J1=0;J1<nBjets;J1++){
    for(int J2=0;J2<nBjets;J2++){
      if(J2==J1) continue; 
      double mH=(bjets.at(J1)+bjets.at(J2)).M();
      if(fabs(mH-mHiggs)<mDiff){
  	mDiff=fabs(mH-mHiggs);
  	bjet1=J1;
  	bjet2=J2;
      }
    }
  } 
  if(bjet1!=-1) higgs1 = bjets.at(bjet1)+bjets.at(bjet2);
  return higgs1; 
}

TLorentzVector findHiggs_method2(auto bjets, auto leadingBjets){
  TLorentzVector higgs2;
  if(bjets.size()>1) higgs2=leadingBjets.at(0)+leadingBjets.at(1); 
  return higgs2; 
}

TLorentzVector findHiggs_method3(auto bjets){
  TLorentzVector higgs3; 
  double smallDR=100; 
  int nBjets = bjets.size();
  int bjet1=-1;
  int bjet2=-1;
  for(int J1=0;J1<nBjets;J1++){
    for(int J2=0;J2<nBjets;J2++){
      if(J2==J1) continue; 
      double DR_bjets=dR_fn(bjets.at(J1).Eta(),bjets.at(J2).Eta(),bjets.at(J1).Phi(),bjets.at(J2).Phi());
      if(smallDR>DR_bjets){
	smallDR=DR_bjets;
  	bjet1=J1;
  	bjet2=J2;
      }
    }
  }
  if(nBjets>1) higgs3 = bjets.at(bjet1)+bjets.at(bjet2);
  return higgs3;
}

TLorentzVector findTop(auto bjets, auto leptons, auto metVec){
  TLorentzVector top;
  double DR=100; 
  int nBjets=bjets.size(); 
  int top_bjet_idx=-1; 
  for(int iJet=0;iJet<nBjets;iJet++){
    double DR_bjetLep=dR_fn(bjets.at(iJet).Eta(),leptons.at(0).Eta(),bjets.at(iJet).Phi(),leptons.at(0).Phi());
    if(DR>DR_bjetLep){
      DR=DR_bjetLep;
      top_bjet_idx=iJet; 
    }
  }
  if(top_bjet_idx!=-1) top=bjets.at(top_bjet_idx)+leptons.at(0)+metVec; 
  return top; 
}



DefineAnalysis(tH2017)

void tH2017::Init()
{
  //addRegions({"SR1lep0b","SR1lep1b","SR1lep2b","SR1lep3b","SR1lep4b"});

  addHistogram("events",10,-0.5,9.5);
  
  //btagged regions
  for(int i=0;i<5;i++){
    std::string SR=std::to_string(i); 
    for(int j=0;j<4;j++){ //1st,2nd,3rd,4th most forward jets
      std::string jet=std::to_string(j+1);
      addHistogram("fwdJet"+jet+"Eta_SRB"+SR,25,0,4);
      addHistogram("fwdJet"+jet+"Pt_SRB"+SR,100,0,500);
    }
    addHistogram("fwdBJet1Eta_SRB"+SR,25,0,4);
    addHistogram("fwdBJet1Pt_SRB"+SR,100,0,500);
    addHistogram("leadJetPt_SRB"+SR,100,0,500); 
    addHistogram("leadJetEta_SRB"+SR,25,0,4); 
    addHistogram("leadBJetPt_SRB"+SR,100,0,500); 
    addHistogram("leadBJetEta_SRB"+SR,25,0,4); 
    addHistogram("deltaEta_fwdJets_SRB"+SR,50,0,2.5); //|eta(fwdjet)-eta(fwdbjet)|
    addHistogram("deltaEta_leadJets_SRB"+SR,50,0,2.5); //|eta(leadjet)-eta(leadbjet)|
    addHistogram("dR_leadBJets_SRB"+SR,100,0,3);//dR(leading bjet, subleading bjet)

    addHistogram("MET_SRB"+SR,100,0,500);
    addHistogram("Njets_SRB"+SR,10,-0.5,9.5);
    addHistogram("NAntiBjets_SRB"+SR,10,-0.5,9.5);
    addHistogram("sumAllPt_SRB"+SR,100,0,2000); 

    addHistogram("lep_pt_SRB"+SR,100,0,500); 
    addHistogram("lep_eta_SRB"+SR,25,0,4); 

    addHistogram("top_m_SRB"+SR,100,0,500); //reco mass of top quark
    addHistogram("top_pt_SRB"+SR,100,0,500); 
    addHistogram("top_eta_SRB"+SR,25,0,4); 
    
    for(int j=0;j<3;j++){ //higgs reconstruction using 3 different methods
      std::string h_idx = std::to_string(j+1); 
      addHistogram("Higgs_m"+h_idx+"_SRB"+SR,100,0,500);
      addHistogram("Higgs_pt"+h_idx+"_SRB"+SR,100,0,500);
      addHistogram("Higgs_eta"+h_idx+"_SRB"+SR,25,0,4);
      addHistogram("deltaEta_Higgs"+h_idx+"FwdJet_SRB"+SR,100,0,8);//|eta(higgs)-eta(forwardJet)|
    }
   }

  //no cuts
  addHistogram("NLep_nocuts",10,-0.5,9.5);
  addHistogram("MET_nocuts",100,0,500);
  addHistogram("Njets_nocuts",10,-0.5,9.5);
  addHistogram("jet_pt_nocuts",100,0,500);
  addHistogram("lep_pt_nocuts",100,0,500); 
  
  //no pt cut on jets/leptons
  addHistogram("jet_pt_noPtCut",100,0,500); 
  addHistogram("lep_pt_noPtCut",100,0,500); 
  
  //after preselection
  addHistogram("MET",100,0,500);
  addHistogram("Njets",10,-0.5,9.5);
  addHistogram("NAntiBjets",10,-0.5,9.5);
  addHistogram("Nbjets",6,-0.5,5.5);
  addHistogram("sumAllPt",100,0,2000); 
 
  addHistogram("lep_pt",100,0,500);
  addHistogram("lep_eta",25,0,4);   

  addHistogram("fwdJet1Eta",25,0,4);
  addHistogram("fwdJet1Pt",100,0,500);
  addHistogram("fwdBJet1Eta",25,0,4);
  addHistogram("fwdBJet1Pt",100,0,500);						
  
  addHistogram("leadJetPt",100,0,500); 
  addHistogram("leadJetEta",25,0,4); 
  addHistogram("leadBJetPt",100,0,500); 
  addHistogram("leadBJetEta",25,0,4); 
  
  addHistogram("deltaEta_fwdJets",100,0,8); //|eta(fwdjet)-eta(fwdbjet)|
  addHistogram("deltaEta_leadJets",100,0,8); //|eta(leadjet)-eta(leadbjet)|
  addHistogram("dR_leadBJets",100,0,3);//dR(leading bjet, subleading bjet)

  addHistogram("top_m",100,0,500); //reco mass of top quark
  addHistogram("top_pt",100,0,500); 
  addHistogram("top_eta",25,0,4);

  for(int j=0;j<3;j++){ //Higgs reconstruction using 3 different methods
    std::string h_idx = std::to_string(j+1); 
    addHistogram("Higgs_m"+h_idx,100,0,500); 
    addHistogram("Higgs_pt"+h_idx,100,0,500);
    addHistogram("Higgs_eta"+h_idx,25,0,4);
    addHistogram("deltaEta_Higgs"+h_idx+"FwdJet",100,0,8);//|eta(higgs)-eta(forwardJet)|
  }
 
}

void tH2017::ProcessEvent(AnalysisEvent *event)
{
  _output->setEventWeight(1);
  fill("events",0); //Number of processed events

  double eventweight=event->getMCWeights()[0];
  _output->setEventWeight(eventweight); //all histograms after this point get filled with this weight
  fill("events",1); //Sum of initial weights


  auto electrons_noPtCut = event->getElectrons(0., 2.47, ELooseBLLH && EIsoGradient);
  auto muons_noPtCut = event->getMuons(0., 2.7, MuLoose && MuIsoGradient);
  auto jets_noPtCut = event->getJets(0., 3.8); 

  auto electrons  = event->getElectrons(25., 2.47,ELooseBLLH && EIsoGradient); //add vertex  (ED0Sigma5|EZ05mm ?)
  auto muons      = event->getMuons(25., 2.7, MuLoose && MuIsoGradient);//add vertex (MuD0Sigma3|MuZ05mm ?)
  auto jets   = event->getJets(25., 3.8, JVT50Jet); // what about NOT(LooseBadJet)
  auto metVec     = event->getMET();
  double met = metVec.Et();
  

  // Overlap removal - including with object Pt-dependent radius calculation
  electrons_noPtCut  = overlapRemoval(electrons_noPtCut, muons_noPtCut, 0.01);
  jets_noPtCut   = overlapRemoval(jets_noPtCut, electrons_noPtCut, 0.2, NOT(BTag85MV2c20));//check this cut
  electrons_noPtCut  = overlapRemoval(electrons_noPtCut, jets_noPtCut, 0.4);
  jets_noPtCut   = overlapRemoval(jets_noPtCut, muons_noPtCut, 0.4, LessThan3Tracks); //check this cut (b-jets?)
  muons_noPtCut      = overlapRemoval(muons_noPtCut, jets_noPtCut, 0.4);   

  electrons  = overlapRemoval(electrons, muons, 0.01);
  jets   = overlapRemoval(jets, electrons, 0.2, NOT(BTag85MV2c20));//check this cut
  electrons  = overlapRemoval(electrons, jets, 0.4);
  jets   = overlapRemoval(jets, muons, 0.4, LessThan3Tracks); //check this cut (b-jets?)
  muons      = overlapRemoval(muons, jets, 0.4);   


  // Lists of objects can be merged by simple addition
  auto leptons   = electrons + muons;
  auto leptons_noPtCut = electrons_noPtCut + muons_noPtCut; 

  //
  auto bjets = filterObjects(jets, 25., 3.8, BTag85MV2c20);
  auto antiBjets = filterObjects(jets, 25., 3.8, NOT(BTag85MV2c20));


  // Object counting
  int numSignalLeptons = leptons.size();  // Object lists are essentially std::vectors so .size() works
  int numSignalJets    = jets.size();
  int nBjets = bjets.size();
  
  // Fill in histograms without cuts
  fill("NLep_nocuts",numSignalLeptons);
  fill("MET_nocuts",met);
  fill("Njets_nocuts",numSignalJets);

  if(numSignalLeptons>0) fill("lep_pt_nocuts", leptons.at(0).Pt()); 
  if(leptons_noPtCut.size()>0) fill("lep_pt_noPtCut",leptons_noPtCut.at(0).Pt());
  
  for(int iJet=0;iJet<numSignalJets;iJet++) fill("jet_pt_nocuts", jets.at(iJet).Pt()); 
  for(unsigned int iJet=0;iJet<jets_noPtCut.size();iJet++) fill("jet_pt_noPtCut", jets_noPtCut.at(iJet).Pt()); 
  
  // Preselection
  if (numSignalLeptons != 1) return;
  if (numSignalJets < 2) return; 
  if (met < 20) return; 

  // Fill histogram after cuts
  fill("MET",met);
  fill("Njets",numSignalJets);
  fill("NAntiBjets",antiBjets.size()); 
  fill("Nbjets",nBjets);
  fill("lep_pt",leptons.at(0).Pt());
  fill("lep_eta",fabs(leptons.at(0).Eta()));
  
  //Sum pT of all objects
  double pTSum=leptons.at(0).Pt();
  for(int i=0;i<numSignalJets;i++) pTSum+=jets.at(i).Pt();
  fill("sumAllPt",pTSum); 

  //Leading/forward jets and bjets
  std::vector<TLorentzVector> leadingJets=findLeadingJets(jets); 
  std::vector<TLorentzVector> forwardJets=findForwardJets(jets);
  std::vector<TLorentzVector> leadingBjets=findLeadingJets(bjets); 
  std::vector<TLorentzVector> forwardBjets=findForwardJets(bjets); 
  std::vector<TLorentzVector> leadingAntiBjets=findLeadingJets(antiBjets); 
  std::vector<TLorentzVector> forwardAntiBjets=findForwardJets(antiBjets); 
   
  fill("leadJetPt", leadingJets.at(0).Pt());
  fill("leadJetEta", fabs(leadingJets.at(0).Eta()));
  fill("fwdJet1Eta", fabs(forwardJets.at(0).Eta()));
  fill("fwdJet1Pt", forwardJets.at(0).Pt());
  
  if(nBjets>0){
    fill("leadBJetPt", leadingBjets.at(0).Pt());
    fill("leadBJetEta", fabs(leadingBjets.at(0).Eta()));
    fill("fwdBJet1Eta", fabs(forwardBjets.at(0).Eta())); 
    fill("fwdBJet1Pt", forwardBjets.at(0).Pt()); 
  }

  //Delta R between 2 leading bjets
  if(nBjets>1) fill("dR_leadBJets", dR_fn(leadingBjets.at(0).Eta(), leadingBjets.at(1).Eta(), leadingBjets.at(0).Phi(), leadingBjets.at(1).Phi()));
  
  //DeltaEta between leading bjet and leading antibjet
  if(antiBjets.size()>0&&nBjets>0) fill("deltaEta_leadJets", fabs(leadingBjets.at(0).Eta()-leadingAntiBjets.at(0).Eta()));
  
  //Reconstructing Higgs
  TLorentzVector higgs1=findHiggs_method1(bjets); 
  TLorentzVector higgs2=findHiggs_method2(bjets,leadingBjets);
  TLorentzVector higgs3=findHiggs_method3(bjets);

  if(nBjets>1){
    fill("Higgs_m1",higgs1.M());
    fill("Higgs_m2",higgs2.M());
    fill("Higgs_m3",higgs3.M());

    fill("Higgs_pt1",higgs1.Pt()); 
    fill("Higgs_pt2",higgs2.Pt());
    fill("Higgs_pt3",higgs3.Pt());
    
    fill("Higgs_eta1",fabs(higgs1.Eta())); 
    fill("Higgs_eta2",fabs(higgs2.Eta()));
    fill("Higgs_eta3",fabs(higgs3.Eta()));

    if(antiBjets.size()>0){
      fill("deltaEta_Higgs1FwdJet",fabs(higgs1.Eta()-forwardAntiBjets.at(0).Eta())); 
      fill("deltaEta_Higgs2FwdJet",fabs(higgs2.Eta()-forwardAntiBjets.at(0).Eta()));
      fill("deltaEta_Higgs3FwdJet",fabs(higgs3.Eta()-forwardAntiBjets.at(0).Eta()));
    }
  }
  
  //Reconstructed top quark
  TLorentzVector top=findTop(bjets,leptons,metVec); 
  if(nBjets>0){
    fill("top_m", top.M()); 
    fill("top_pt", top.Pt());
    fill("top_eta", top.Eta()); 
  }
  
  //______________Fill histos corresponding to b-Tag signal regions____________
  if(nBjets>4) return; 
  std::string SR=std::to_string(nBjets); 
  
  for(int i=0;i<forwardJets.size();i++){//first 4 forward jets   
    fill("fwdJet"+std::to_string(i+1)+"Pt_SRB"+SR, forwardJets.at(i).Pt());
    fill("fwdJet"+std::to_string(i+1)+"Eta_SRB"+SR, fabs(forwardJets.at(i).Eta()));
    
  }
  if(nBjets>0){
    fill("fwdBJet1Eta_SRB"+SR,forwardBjets.at(0).Eta()); 
    fill("fwdBJet1Pt_SRB"+SR, forwardBjets.at(0).Pt()); 
  }  

  fill("MET_SRB"+SR,met);
  fill("sumAllPt_SRB"+SR,pTSum); 
  fill("lep_pt_SRB"+SR,leptons.at(0).Pt());
  fill("lep_eta_SRB"+SR,fabs(leptons.at(0).Eta())); 

  fill("NAntiBjets_SRB"+SR,antiBjets.size());  
  fill("Njets_SRB"+SR,numSignalJets);
  fill("leadJetPt_SRB"+SR, leadingJets.at(0).Pt());  
  fill("leadJetEta_SRB"+SR, leadingJets.at(0).Eta());
  
    
  if(nBjets>0){
    fill("leadBJetPt_SRB"+SR, leadingBjets.at(0).Pt());
    fill("leadBJetEta_SRB"+SR, fabs(leadingBjets.at(0).Eta()));
    fill("deltaEta_fwdJets_SRB"+SR,fabs(forwardBjets.at(0).Eta()-forwardJets.at(0).Eta())); 
    fill("deltaEta_leadJets_SRB"+SR, fabs(leadingBjets.at(0).Eta()-leadingJets.at(0).Eta()));   

    fill("top_m_SRB"+SR,top.M()); 
    fill("top_pt_SRB"+SR, top.Pt()); 
    fill("top_eta_SRB"+SR, fabs(top.Eta())); 
  }

  if(nBjets>1){
    fill("Higgs_m1_SRB"+SR,higgs1.M()); 
    fill("Higgs_pt1_SRB"+SR,higgs1.Pt()); 
    fill("Higgs_eta1_SRB"+SR,fabs(higgs1.Eta())); 
    fill("Higgs_m2_SRB"+SR,higgs2.M()); 
    fill("Higgs_pt2_SRB"+SR,higgs2.Pt()); 
    fill("Higgs_eta2_SRB"+SR,fabs(higgs2.Eta())); 
    fill("Higgs_m3_SRB"+SR,higgs3.M()); 
    fill("Higgs_pt3_SRB"+SR,higgs3.Pt()); 
    fill("Higgs_eta3_SRB"+SR,fabs(higgs3.Eta())); 

    if(antiBjets.size()>0){
      fill("deltaEta_Higgs1FwdJet_SRB"+SR,fabs(higgs1.Eta()-forwardAntiBjets.at(0).Eta())); 
      fill("deltaEta_Higgs2FwdJet_SRB"+SR,fabs(higgs2.Eta()-forwardAntiBjets.at(0).Eta()));
      fill("deltaEta_Higgs3FwdJet_SRB"+SR,fabs(higgs3.Eta()-forwardAntiBjets.at(0).Eta()));
    }
  }

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
