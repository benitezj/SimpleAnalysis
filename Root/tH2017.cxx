#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

double dR_fn (float eta1, float eta2, float phi1, float phi2){
  double deta= fabs(eta1 - eta2);      double dphi= fabs(phi1 - phi2);
  if (dphi > 3.14 ) dphi = 2*3.14 - dphi;
  return sqrt((dphi*dphi)+(deta*deta));
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
    addHistogram("deltaEta_Higgs_fwdJet_SRB"+SR,100,0,8);//|eta(higgs)-eta(forwardJet)|

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
  addHistogram("deltaEta_Higgs_fwdJet",100,0,8);//|eta(higgs)-eta(forwardJet)|

  addHistogram("top_m",100,0,500); //reco mass of top quark
  addHistogram("top_pt",100,0,500); 
  addHistogram("top_eta",25,0,4);

  for(int j=0;j<3;j++){ //Higgs reconstruction using 3 different methods
    std::string h_idx = std::to_string(j+1); 
    addHistogram("Higgs_m"+h_idx,100,0,500); 
    addHistogram("Higgs_pt"+h_idx,100,0,500);
    addHistogram("Higgs_eta"+h_idx,25,0,4);
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
  if(numSignalLeptons>0) fill("lep_pt_nocuts", leptons.at(0).Pt()); 
  if(leptons_noPtCut.size()>0) fill("lep_pt_noPtCut",leptons_noPtCut.at(0).Pt());
  fill("MET_nocuts",met);
  fill("Njets_nocuts",numSignalJets);
  for(int iJet=0;iJet<numSignalJets;iJet++) fill("jet_pt_nocuts", jets.at(iJet).Pt()); 
  for(unsigned int iJet=0;iJet<jets_noPtCut.size();iJet++) fill("jet_pt_noPtCut", jets_noPtCut.at(iJet).Pt()); 
  
  // Preselection
  if (numSignalLeptons != 1) return;
  if (numSignalJets < 2) return; 
  if (met < 20) return; 

  // Fill histogram after cuts
  fill("MET",met);
  fill("Njets",numSignalJets);
  fill("Nbjets",nBjets);
  fill("lep_pt",leptons.at(0).Pt());
  fill("lep_eta",fabs(leptons.at(0).Eta()));
  
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
  fill("leadJetPt", jets.at(highPtJet_idx).Pt());
  fill("leadJetEta", fabs(jets.at(highPtJet_idx).Eta()));

  //Find 1st/2nd leading bjets 
  double highPtB[2]={-1,-1};
  double highPtBJet_idx[2]={-1,-1};
  for(int iJet=0;iJet<nBjets;iJet++){
    if(bjets.at(iJet).Pt()>highPtB[0]){ 
      highPtB[0]=bjets.at(iJet).Pt();
      highPtBJet_idx[0]=iJet;
    }
  }
  for(int iJet=0;iJet<nBjets;iJet++){
    if(iJet==highPtBJet_idx[0]) continue; 
    if(bjets.at(iJet).Pt()>highPtB[1]){ 
      highPtB[1]=bjets.at(iJet).Pt();
      highPtBJet_idx[1]=iJet;
    }
  }
  if(highPtBJet_idx[0]!=-1){
    fill("leadBJetPt", bjets.at(highPtBJet_idx[0]).Pt());
    fill("leadBJetEta", fabs(bjets.at(highPtBJet_idx[0]).Eta()));
  }

  //Delta R between 2 leading bjets
  if(nBjets>1) fill("dR_leadBJets", dR_fn(bjets.at(highPtBJet_idx[0]).Eta(), bjets.at(highPtBJet_idx[1]).Eta(), bjets.at(highPtBJet_idx[0]).Phi(),  bjets.at(highPtBJet_idx[1]).Phi()));
  
  //Find deltaEta btwn leading bjet and leading antibjet
  int anti_idx=-1;
  double highPt_anti=-1; 
  for(int iJet=0;iJet<antiBjets.size();iJet++){
    if(antiBjets.at(iJet).Pt()>highPt_anti){
      highPt_anti = antiBjets.at(iJet).Pt();
      anti_idx=iJet; 
    }
  }
  if(anti_idx!=-1&&highPtBJet_idx[0]!=-1) fill("deltaEta_leadJets", fabs(bjets.at(highPtBJet_idx[0]).Eta()-antiBjets.at(anti_idx).Eta()));

   

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
  fill("fwdJet1Eta", fabs(jets.at(fwdJet_idx[0]).Eta()));
  fill("fwdJet1Pt", jets.at(fwdJet_idx[0]).Pt()); 

  //Find most forward bjet per event
  double fwdBJet_eta=-999;
  int fwdBJet_idx=-1;
  for(int iJet=0;iJet<nBjets;iJet++){ 
    if(fabs(bjets.at(iJet).Eta())>fwdBJet_eta){
      fwdBJet_eta=fabs(bjets.at(iJet).Eta());
      fwdBJet_idx=iJet;
    }
  }
  if(fwdBJet_idx!=-1){
    fill("fwdBJet1Eta",fabs(bjets.at(fwdBJet_idx).Eta()));
    fill("fwdBJet1Pt",fabs(bjets.at(fwdBJet_idx).Pt()));
    fill("deltaEta_fwdJets",fabs(bjets.at(fwdBJet_idx).Eta()-jets.at(fwdJet_idx[0]).Eta())); 
  }

  /////--------Reconstructing Higgs---------------------------------

  //First method of building Higgs candidate --find 2 jets whose combined mass is closest to the Higgs
  double mHiggs=125; 
  int BJET1=-1;
  int BJET2=-1;
  double mDiff=1000; 
  for(int J1=0;J1<nBjets;J1++){
    for(int J2=0;J2<nBjets;J2++){
      if(J2==J1) continue; 
      double mH=(bjets.at(J1)+bjets.at(J2)).M();
      if(fabs(mH-mHiggs)<mDiff){
  	mDiff=fabs(mH-mHiggs);
  	BJET1=J1;
  	BJET2=J2;
      }
    }
  }

  double higgs_m1=-1;
  double higgs_pt1=-1;
  double higgs_eta1=-999;

  if(nBjets>1&&BJET1!=-1){ 
    higgs_m1=(bjets.at(BJET1)+bjets.at(BJET2)).M();
    higgs_pt1=(bjets.at(BJET1)+bjets.at(BJET2)).Pt();
    higgs_eta1=(bjets.at(BJET1)+bjets.at(BJET2)).Eta();
  }
  fill("Higgs_m1",higgs_m1);
  fill("Higgs_pt1",higgs_pt1); 
  fill("Higgs_eta1",fabs(higgs_eta1)); 
  
  //Second method of building Higgs candidate-- use leading bjets
  double higgs_m2=-1;
  double higgs_pt2=-1;
  double higgs_eta2=-999;
  if(nBjets>1){
    higgs_m2=(bjets.at(highPtBJet_idx[0])+bjets.at(highPtBJet_idx[1])).M();
    higgs_pt2=(bjets.at(highPtBJet_idx[0])+bjets.at(highPtBJet_idx[1])).Pt();
    higgs_eta2=(bjets.at(highPtBJet_idx[0])+bjets.at(highPtBJet_idx[1])).Eta();
  }
  fill("Higgs_m2",higgs_m2);
  fill("Higgs_pt2",higgs_pt2);
  fill("Higgs_eta2",fabs(higgs_eta2)); 

  //Third method of building Higgs candidate --use 2 bjets with smallest dR
  double higgs_m3=-1;
  double higgs_pt3=-1;
  double higgs_eta3=-999;  
  double smallDR=100; 
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
 
  if(nBjets>1){
    higgs_m3=(bjets.at(bjet1)+bjets.at(bjet2)).M();
    higgs_pt3=(bjets.at(bjet1)+bjets.at(bjet2)).Pt();
    higgs_eta3=(bjets.at(bjet1)+bjets.at(bjet2)).Eta();
  }
  
  fill("Higgs_m3",higgs_m3);
  fill("Higgs_pt3",higgs_pt3);
  fill("Higgs_eta3",fabs(higgs_eta3)); 

  //find mass of reconstructed top quark
  //----->first find bjet closest to lepton
  double DR=100; 
  int top_bjet_idx=-1; 
  for(int iJet=0;iJet<nBjets;iJet++){
    double DR_bjetLep=dR_fn(bjets.at(iJet).Eta(),leptons.at(0).Eta(),bjets.at(iJet).Phi(),leptons.at(0).Phi());
    if(DR>DR_bjetLep){
      DR=DR_bjetLep;
      top_bjet_idx=iJet; 
    }
  }
  double topM=-1;
  double topPt=-1;
  double topEta=-999;
  if(top_bjet_idx!=-1){
    topM=(bjets.at(top_bjet_idx)+leptons.at(0)+metVec).M();
    topPt=(bjets.at(top_bjet_idx)+leptons.at(0)+metVec).Pt();
    topEta=(bjets.at(top_bjet_idx)+leptons.at(0)+metVec).Eta();
    fill("top_m",topM);
    fill("top_pt",topPt); 
    fill("top_eta",fabs(topEta)); 
  }

  //fill histos corresponding to b-Tag signal regions
  if(nBjets>4) return; 
  std::string SR=std::to_string(nBjets); 
  
  for(int i=0;i<4;i++){//first 4 forward jets
    std::string jetPosition=std::to_string(i+1);
    if(fwdJet_idx[i]!=-1){
      fill("fwdJet"+jetPosition+"Pt_SRB"+SR, jets.at(fwdJet_idx[i]).Pt());
      fill("fwdJet"+jetPosition+"Eta_SRB"+SR, fabs(jets.at(fwdJet_idx[i]).Eta()));
    }
  }
  fill("fwdBJet1Eta_SRB"+SR, fwdBJet_eta); 
  if(fwdBJet_idx!=-1) fill("fwdBJet1Pt_SRB"+SR, bjets.at(fwdBJet_idx).Pt()); 
  
  fill("lep_pt_SRB"+SR,leptons.at(0).Pt());
  fill("lep_eta_SRB"+SR,fabs(leptons.at(0).Eta())); 
  fill("leadJetPt_SRB"+SR, jets.at(highPtJet_idx).Pt());  
  fill("leadJetEta_SRB"+SR, fabs(jets.at(highPtJet_idx).Eta()));
  fill("Higgs_m1_SRB"+SR,higgs_m1); 
  fill("Higgs_pt1_SRB"+SR,higgs_pt1); 
  fill("Higgs_eta1_SRB"+SR,fabs(higgs_eta1)); 
    
  if(nBjets>0){
    fill("leadBJetPt_SRB"+SR, bjets.at(highPtBJet_idx[0]).Pt());
    fill("leadBJetEta_SRB"+SR, fabs(bjets.at(highPtBJet_idx[0]).Eta()));
    fill("deltaEta_fwdJets_SRB"+SR,fabs(bjets.at(fwdBJet_idx).Eta()-jets.at(fwdJet_idx[0]).Eta())); 
    fill("deltaEta_leadJets_SRB"+SR, fabs(bjets.at(highPtBJet_idx[0]).Eta()-jets.at(highPtJet_idx).Eta()));
    fill("Higgs_m2_SRB"+SR,higgs_m2); 
    fill("Higgs_pt2_SRB"+SR,higgs_pt2); 
    fill("Higgs_eta2_SRB"+SR,fabs(higgs_eta2)); 
    fill("Higgs_m3_SRB"+SR,higgs_m3); 
    fill("Higgs_pt3_SRB"+SR,higgs_pt3); 
    fill("Higgs_eta3_SRB"+SR,fabs(higgs_eta3)); 

    fill("top_m_SRB"+SR,topM); 
    fill("top_pt_SRB"+SR, topPt); 
    fill("top_eta_SRB"+SR, fabs(topEta)); 
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
