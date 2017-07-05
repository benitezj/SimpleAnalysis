#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

double dR_fn (float eta1, float eta2, float phi1, float phi2){
  double deta= fabs(eta1 - eta2);      double dphi= fabs(phi1 - phi2);
  if (dphi > 3.14 ) dphi = 2*3.14 - dphi;
  return sqrt((dphi*dphi)+(deta*deta));
}

std::vector<TLorentzVector> findForwardJets(AnalysisObjects jetVec){
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

std::vector<TLorentzVector> findLeadingJets(AnalysisObjects jetVec){
  std::vector<TLorentzVector> leadingJets;
  double highPt[2]={-1,-1}; 
  double highPtJet_idx[2]={-1,-1}; 
  for(unsigned int iJet=0;iJet<jetVec.size();iJet++){
    if(jetVec.at(iJet).Pt()>highPt[0]){ 
      highPt[0]=jetVec.at(iJet).Pt();
      highPtJet_idx[0]=iJet;
    }
  }
  for(unsigned int iJet=0;iJet<jetVec.size();iJet++){
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

TLorentzVector findHiggs_method1(AnalysisObjects bjets){
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

TLorentzVector findHiggs_method2(AnalysisObjects bjets, std::vector<TLorentzVector> leadingBjets){
  TLorentzVector higgs2;
  if(bjets.size()>1) higgs2=leadingBjets.at(0)+leadingBjets.at(1); 
  return higgs2; 
}

TLorentzVector findHiggs_method3(AnalysisObjects bjets){
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

TLorentzVector findTop(AnalysisObjects bjets, AnalysisObjects leptons, AnalysisObject metVec){
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

  //no cuts
  addHistogram("NLep_nocuts",10,-0.5,9.5);
  addHistogram("MET_nocuts",100,0,500);
  addHistogram("Njets_nocuts",10,-0.5,9.5);
  addHistogram("jet_pt_nocuts",100,0,500);
  addHistogram("lep_pt_nocuts",100,0,500); 
  
  addHistogram("Njets_HGTD_nocuts",10,-0.5,9.5); 
  addHistogram("Nbjets_HGTD_nocuts",6,-0.5,5.5); 
  
  //no pt cut on jets/leptons
  addHistogram("jet_pt_noPtCut",100,0,500); 
  addHistogram("lep_pt_noPtCut",100,0,500); 
  
  //after preselection
  addHistogram("MET",100,0,500);
  addHistogram("Njets",10,-0.5,9.5);
  addHistogram("NLightjets",10,-0.5,9.5);
  addHistogram("Nbjets",6,-0.5,5.5);
  addHistogram("sumAllPt",100,0,2000); 
 
  addHistogram("lep_pt",100,0,500);
  addHistogram("lep_eta",25,0,4);   

  addHistogram("leadJetPt_allJets",100,0,500);
  addHistogram("leadJetEta_allJets",25,0,4); 

  addHistogram("fwdJetPt",100,0,500);
  addHistogram("fwdJetEta",25,0,4);
  addHistogram("leadJetPt",100,0,500); 
  addHistogram("leadJetEta",25,0,4); 

  addHistogram("fwdBJetPt",100,0,500);
  addHistogram("fwdBJetEta",25,0,4);
  addHistogram("leadBJetPt",100,0,500); 
  addHistogram("leadBJetEta",25,0,4); 
  
  addHistogram("deltaEta_jfwd_bfwd",100,0,8); //|eta(fwdjet)-eta(fwdbjet)|
  addHistogram("deltaEta_jfwd_b1",100,0,8);
  addHistogram("deltaEta_j1_b1",100,0,8); //|eta(leadjet)-eta(leadbjet)|
  addHistogram("dR_b1_b2",100,0,3);//dR(leading bjet, subleading bjet)

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

  //btagged regions
  for(int i=2;i<5;i++){

    std::string SR=std::to_string(i); 

    addHistogram("MET_SRB"+SR,100,0,500);
    addHistogram("Njets_SRB"+SR,10,-0.5,9.5);
    addHistogram("NLightjets_SRB"+SR,10,-0.5,9.5);
    addHistogram("sumAllPt_SRB"+SR,100,0,2000); 

    addHistogram("lep_pt_SRB"+SR,100,0,500); 
    addHistogram("lep_eta_SRB"+SR,25,0,4); 

    addHistogram("leadJetPt_allJets_SRB"+SR,100,0,500); 
    addHistogram("leadJetEta_allJets_SRB"+SR,25,0,4);

    addHistogram("fwdJetPt_SRB"+SR,100,0,500);
    addHistogram("fwdJetEta_SRB"+SR,25,0,4);
    addHistogram("leadJetPt_SRB"+SR,100,0,500); 
    addHistogram("leadJetEta_SRB"+SR,25,0,4);
    
    addHistogram("fwdJetPtH1_SRMbb_SRB"+SR,100,0,500);
    addHistogram("fwdJetEtaH1_SRMbb_SRB"+SR,25,0,4);
    addHistogram("fwdJetPtH2_SRMbb_SRB"+SR,100,0,500);
    addHistogram("fwdJetEtaH2_SRMbb_SRB"+SR,25,0,4);
    addHistogram("fwdJetPtH3_SRMbb_SRB"+SR,100,0,500);
    addHistogram("fwdJetEtaH3_SRMbb_SRB"+SR,25,0,4);

    addHistogram("fwdBJetPt_SRB"+SR,100,0,500);
    addHistogram("fwdBJetEta_SRB"+SR,25,0,4);     
    addHistogram("leadBJetPt_SRB"+SR,100,0,500); 
    addHistogram("leadBJetEta_SRB"+SR,25,0,4);  
    
    addHistogram("deltaEta_jfwd_bfwd_SRB"+SR,100,0,8); //|eta(fwdjet)-eta(fwdbjet)|
    addHistogram("deltaEta_j1_b1_SRB"+SR,100,0,8); //|eta(leadjet)-eta(leadbjet)|
    addHistogram("deltaEta_jfwd_b1_SRB"+SR,100,0,8);
    addHistogram("dR_b1_b2_SRB"+SR,100,0,3);//dR(leading bjet, subleading bjet)

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
  
  // Jets in HGTD acceptance
  int njets_hgtd(0), nbjets_hgtd(0);
  for (auto j : jets) {
    if ( fabs(j.Eta()) >= 2.4 && fabs(j.Eta()) <= 4.3) {
      ++njets_hgtd;
      if (j.pass(GoodBJet)) ++nbjets_hgtd;
    }
  }
  fill("Njets_HGTD_nocuts", njets_hgtd);
  fill("Nbjets_HGTD_nocuts", nbjets_hgtd);

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
  fill("NLightjets",antiBjets.size()); 
  fill("Nbjets",nBjets);
  fill("lep_pt",leptons.at(0).Pt());
  fill("lep_eta",fabs(leptons.at(0).Eta()));  

  //Leading/forward jets and bjets
  std::vector<TLorentzVector> leadingJets=findLeadingJets(jets); 
  std::vector<TLorentzVector> forwardJets=findForwardJets(jets);
  std::vector<TLorentzVector> leadingBjets=findLeadingJets(bjets); 
  std::vector<TLorentzVector> forwardBjets=findForwardJets(bjets); 
  std::vector<TLorentzVector> leadingLightjets=findLeadingJets(antiBjets); 
  std::vector<TLorentzVector> forwardLightjets=findForwardJets(antiBjets); 
   
  fill("leadJetPt_allJets", leadingJets.at(0).Pt());
  fill("leadJetEta_allJets", fabs(leadingJets.at(0).Eta()));
  
  ////______________bjet cut__________________________
  if( nBjets<2 || nBjets>4 || antiBjets.size()==0)
    return;   
  std::string SR=std::to_string(nBjets); 
  
  ///Sum pT of all objects 
  double pTSum=leptons.at(0).Pt();
  for(int i=0;i<numSignalJets;i++) pTSum+=jets.at(i).Pt();
  fill("sumAllPt",pTSum); //(for >b-tags inclusive)

  ///lead anti b-jet
  fill("fwdJetPt",   forwardLightjets.at(0).Pt()); 
  fill("fwdJetEta",  fabs(forwardLightjets.at(0).Eta())); 
  fill("leadJetPt",  leadingLightjets.at(0).Pt());
  fill("leadJetEta", fabs(leadingLightjets.at(0).Eta()));

  ///lead b-jet
  fill("fwdBJetPt",   forwardBjets.at(0).Pt()); 
  fill("fwdBJetEta",  fabs(forwardBjets.at(0).Eta())); 
  fill("leadBJetPt",  leadingBjets.at(0).Pt());
  fill("leadBJetEta", fabs(leadingBjets.at(0).Eta()));

  //Delta R between 2 leading bjets
  fill("dR_b1_b2", dR_fn(leadingBjets.at(0).Eta(), leadingBjets.at(1).Eta(), leadingBjets.at(0).Phi(), leadingBjets.at(1).Phi()));
  
  //DeltaEta 
  fill("deltaEta_jfwd_bfwd", fabs(forwardLightjets.at(0).Eta()-forwardBjets.at(0).Eta())); 
  fill("deltaEta_j1_b1",     fabs(leadingLightjets.at(0).Eta()-leadingBjets.at(0).Eta())); 
  fill("deltaEta_jfwd_b1",   fabs(forwardLightjets.at(0).Eta()-leadingBjets.at(0).Eta())); 

  //Reconstructing Higgs
  TLorentzVector higgs1=findHiggs_method1(bjets); 
  TLorentzVector higgs2=findHiggs_method2(bjets,leadingBjets);
  TLorentzVector higgs3=findHiggs_method3(bjets);

  fill("Higgs_m1",higgs1.M());
  fill("Higgs_m2",higgs2.M());
  fill("Higgs_m3",higgs3.M());

  fill("Higgs_pt1",higgs1.Pt()); 
  fill("Higgs_pt2",higgs2.Pt());
  fill("Higgs_pt3",higgs3.Pt());
    
  fill("Higgs_eta1",fabs(higgs1.Eta())); 
  fill("Higgs_eta2",fabs(higgs2.Eta()));
  fill("Higgs_eta3",fabs(higgs3.Eta()));

  fill("deltaEta_Higgs1FwdJet",fabs(higgs1.Eta()-forwardLightjets.at(0).Eta())); 
  fill("deltaEta_Higgs2FwdJet",fabs(higgs2.Eta()-forwardLightjets.at(0).Eta()));
  fill("deltaEta_Higgs3FwdJet",fabs(higgs3.Eta()-forwardLightjets.at(0).Eta()));
  
  //Reconstructed top quark
  TLorentzVector top=findTop(bjets,leptons,metVec); 
  fill("top_m", top.M()); 
  fill("top_pt", top.Pt());
  fill("top_eta", top.Eta());   


  ///////////////______________Fill histos corresponding to b-Tag signal regions____________
  fill("MET_SRB"+SR,met);
  fill("Njets_SRB"+SR,numSignalJets);
  fill("NLightjets_SRB"+SR,antiBjets.size());  
  fill("sumAllPt_SRB"+SR,pTSum); 

  fill("lep_pt_SRB"+SR,leptons.at(0).Pt());
  fill("lep_eta_SRB"+SR,fabs(leptons.at(0).Eta())); 

  fill("leadJetPt_allJets_SRB"+SR, leadingJets.at(0).Pt());
  fill("leadJetEta_allJets_SRB"+SR, fabs(leadingJets.at(0).Eta()));

  fill("fwdJetPt_SRB"+SR,    forwardLightjets.at(0).Pt());
  fill("fwdJetEta_SRB"+SR,   fabs(forwardLightjets.at(0).Eta())); 
  fill("leadJetPt_SRB"+SR,   leadingLightjets.at(0).Pt());
  fill("leadJetEta_SRB"+SR,  fabs(leadingLightjets.at(0).Eta()));  

  fill("fwdBJetPt_SRB"+SR,   forwardBjets.at(0).Pt()); 
  fill("fwdBJetEta_SRB"+SR,  fabs(forwardBjets.at(0).Eta()));   			     
  fill("leadBJetPt_SRB"+SR,  leadingBjets.at(0).Pt());
  fill("leadBJetEta_SRB"+SR, fabs(leadingBjets.at(0).Eta()));
  
  fill("dR_b1_b2_SRB"+SR, dR_fn(leadingBjets.at(0).Eta(), leadingBjets.at(1).Eta(), leadingBjets.at(0).Phi(), leadingBjets.at(1).Phi()));

  fill("deltaEta_jfwd_bfwd_SRB"+SR, fabs(forwardLightjets.at(0).Eta()-forwardBjets.at(0).Eta())); 
  fill("deltaEta_j1_b1_SRB"+SR,     fabs(leadingLightjets.at(0).Eta()-leadingBjets.at(0).Eta())); 
  fill("deltaEta_jfwd_b1_SRB"+SR,   fabs(forwardLightjets.at(0).Eta()-leadingBjets.at(0).Eta())); 

  fill("top_m_SRB"+SR,      top.M()); 
  fill("top_pt_SRB"+SR,     top.Pt()); 
  fill("top_eta_SRB"+SR,    fabs(top.Eta()));   

  fill("Higgs_m1_SRB"+SR,   higgs1.M()); 
  fill("Higgs_pt1_SRB"+SR,  higgs1.Pt()); 
  fill("Higgs_eta1_SRB"+SR, fabs(higgs1.Eta())); 
  fill("Higgs_m2_SRB"+SR,   higgs2.M()); 
  fill("Higgs_pt2_SRB"+SR,  higgs2.Pt()); 
  fill("Higgs_eta2_SRB"+SR, fabs(higgs2.Eta())); 
  fill("Higgs_m3_SRB"+SR,   higgs3.M()); 
  fill("Higgs_pt3_SRB"+SR,  higgs3.Pt()); 
  fill("Higgs_eta3_SRB"+SR, fabs(higgs3.Eta())); 


  ////Higgs mass signal regions
  if(80<higgs1.M()&&higgs1.M()<130){
    fill("fwdJetPtH1_SRMbb_SRB"+SR, forwardLightjets.at(0).Pt());
    fill("fwdJetEtaH1_SRMbb_SRB"+SR, forwardLightjets.at(0).Eta());
    fill("deltaEta_Higgs1FwdJet_SRB"+SR,fabs(higgs1.Eta()-forwardLightjets.at(0).Eta())); 
  }
  if(80<higgs2.M()&&higgs2.M()<130){
    fill("fwdJetPtH2_SRMbb_SRB"+SR, forwardLightjets.at(0).Pt());
    fill("fwdJetEtaH2_SRMbb_SRB"+SR, forwardLightjets.at(0).Eta());
    fill("deltaEta_Higgs2FwdJet_SRB"+SR,fabs(higgs2.Eta()-forwardLightjets.at(0).Eta()));
  }
  if(80<higgs3.M()&&higgs3.M()<130){
    fill("fwdJetPtH3_SRMbb_SRB"+SR, forwardLightjets.at(0).Pt());
    fill("fwdJetEtaH3_SRMbb_SRB"+SR, forwardLightjets.at(0).Eta());
    fill("deltaEta_Higgs3FwdJet_SRB"+SR,fabs(higgs3.Eta()-forwardLightjets.at(0).Eta()));
  }

  return;
}

