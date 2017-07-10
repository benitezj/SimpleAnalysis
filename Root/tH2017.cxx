#include "SimpleAnalysis/AnalysisClass.h"
#include <string>

struct eta_sort {
  bool operator()(const TLorentzVector& v1, const TLorentzVector& v2) const {
    return fabs(v1.Eta()) > fabs(v2.Eta());
  }
};

struct mH_sort {
  bool operator()(const std::pair<AnalysisObject,AnalysisObject> &p1, 
               const std::pair<AnalysisObject,AnalysisObject> &p2) const {    
    return fabs(125.-(p1.first+p1.second).M()) < fabs(125.-(p2.first+p2.second).M());
  }
};

struct dR_sort {
  bool operator()(const std::pair<AnalysisObject,AnalysisObject> &p1, 
               const std::pair<AnalysisObject,AnalysisObject> &p2) const {    
    return p1.first.DeltaR(p1.second) < p2.first.DeltaR(p2.second);
  }
};

void sortObjectsByEta(AnalysisObjects& cands) {
  std::sort(cands.begin(), cands.end(), eta_sort());
}

void sortPairsClosestToHiggs(std::vector<std::pair<AnalysisObject,AnalysisObject>>& cands) {
  std::sort(cands.begin(), cands.end(), mH_sort());
}

void sortPairsClosestdR(std::vector<std::pair<AnalysisObject,AnalysisObject>>& cands) {
  std::sort(cands.begin(), cands.end(), dR_sort());
}


void findHiggs(const AnalysisObjects& bjets, TLorentzVector& h1, TLorentzVector& h2, TLorentzVector &h3){
  // Make vector of b-jet pairs
  std::vector<std::pair<AnalysisObject, AnalysisObject>> bPairs;
  for (unsigned int i=0; i < bjets.size(); ++i) {
    for (unsigned int j=i+1; j < bjets.size(); ++j) {
      bPairs.push_back({bjets.at(i),bjets.at(j)});
    } 
  } 
  
  // Method 1 - closest to Higgs mass
  sortPairsClosestToHiggs(bPairs);
  h1 = bPairs.at(0).first+bPairs.at(0).second;
  
  // Metohd 2 - leading b-jets
  h2 = bjets.at(0) + bjets.at(1);
  
  // Method 3 - closest dR
  sortPairsClosestdR(bPairs);
  h3 = bPairs.at(0).first+bPairs.at(0).second;
}


TLorentzVector findTop(AnalysisObjects& bjets, AnalysisObjects& leptons, AnalysisObject& metVec){
  TLorentzVector top;
  double DR=100; 
  int nBjets=bjets.size(); 
  int top_bjet_idx=-1; 
  for(int iJet=0;iJet<nBjets;iJet++){
    double DR_bjetLep=bjets.at(iJet).DeltaR(leptons.at(0));
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
  addHistogram("Njets_HGTD_cuts",10,-0.5,9.5); 
  addHistogram("Nbjets_HGTD_cuts",6,-0.5,5.5); 
  
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
  auto jets_HGTD = event->getJets(25., 4.5, JVT50Jet); // what about NOT(LooseBadJet)
  auto jets_HGTD_nocuts = event->getJets(0., 4.5, JVT50Jet); // what about NOT(LooseBadJet)
  auto metVec     = event->getMET();
  double met = metVec.Et();
  
  // sorting objects by pt
  sortObjectsByPt(electrons);
  sortObjectsByPt(muons);
  sortObjectsByPt(jets);
  
  // Jets in HGTD acceptance
  int njets_hgtd(0), nbjets_hgtd(0);
  for (auto j : jets_HGTD_nocuts) {
    if ( fabs(j.Eta()) >= 2.4 && fabs(j.Eta()) <= 4.3) {
      ++njets_hgtd;
      if (j.pass(GoodBJet)) ++nbjets_hgtd;
    }
  }
  fill("Njets_HGTD_nocuts", njets_hgtd);
  fill("Nbjets_HGTD_nocuts", nbjets_hgtd);
  njets_hgtd = 0; 
  nbjets_hgtd = 0;
  for (auto j : jets_HGTD) {
    if ( fabs(j.Eta()) >= 2.4 && fabs(j.Eta()) <= 4.3) {
      ++njets_hgtd;
      if (j.pass(GoodBJet)) ++nbjets_hgtd;
    }
  }
  fill("Njets_HGTD_cuts", njets_hgtd);
  fill("Nbjets_HGTD_cuts", nbjets_hgtd);
  

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

  // eta-sorted jets
  auto forwardJets = jets;
  sortObjectsByEta(forwardJets);
  auto forwardBjets = filterObjects(forwardJets, 25., 3.8, BTag85MV2c20);
  auto forwardLightjets = filterObjects(forwardJets, 25., 3.8, NOT(BTag85MV2c20));
  
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
  fill("leadJetPt_allJets", jets.at(0).Pt());
  fill("leadJetEta_allJets", fabs(jets.at(0).Eta()));
  
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
  fill("leadJetPt",  antiBjets.at(0).Pt());
  fill("leadJetEta", fabs(antiBjets.at(0).Eta()));

  ///lead b-jet
  fill("fwdBJetPt",   forwardBjets.at(0).Pt()); 
  fill("fwdBJetEta",  fabs(forwardBjets.at(0).Eta())); 
  fill("leadBJetPt",  bjets.at(0).Pt());
  fill("leadBJetEta", fabs(bjets.at(0).Eta()));

  //Delta R between 2 leading bjets
  fill("dR_b1_b2", bjets.at(0).DeltaR(bjets.at(1)));
  
  //DeltaEta 
  fill("deltaEta_jfwd_bfwd", fabs(forwardLightjets.at(0).Eta()-forwardBjets.at(0).Eta())); 
  fill("deltaEta_j1_b1",     fabs(antiBjets.at(0).Eta()-bjets.at(0).Eta())); 
  fill("deltaEta_jfwd_b1",   fabs(forwardLightjets.at(0).Eta()-bjets.at(0).Eta())); 

  //Reconstructing Higgs
  TLorentzVector higgs1, higgs2, higgs3;
  findHiggs(bjets, higgs1, higgs2, higgs3); 

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

  fill("leadJetPt_allJets_SRB"+SR, jets.at(0).Pt());
  fill("leadJetEta_allJets_SRB"+SR, fabs(jets.at(0).Eta()));

  fill("fwdJetPt_SRB"+SR,    forwardLightjets.at(0).Pt());
  fill("fwdJetEta_SRB"+SR,   fabs(forwardLightjets.at(0).Eta())); 
  fill("leadJetPt_SRB"+SR,   antiBjets.at(0).Pt());
  fill("leadJetEta_SRB"+SR,  fabs(antiBjets.at(0).Eta()));  

  fill("fwdBJetPt_SRB"+SR,   forwardBjets.at(0).Pt()); 
  fill("fwdBJetEta_SRB"+SR,  fabs(forwardBjets.at(0).Eta()));   			     
  fill("leadBJetPt_SRB"+SR,  bjets.at(0).Pt());
  fill("leadBJetEta_SRB"+SR, fabs(bjets.at(0).Eta()));
  
  fill("dR_b1_b2_SRB"+SR, bjets.at(0).DeltaR(bjets.at(1)));

  fill("deltaEta_jfwd_bfwd_SRB"+SR, fabs(forwardLightjets.at(0).Eta()-forwardBjets.at(0).Eta())); 
  fill("deltaEta_j1_b1_SRB"+SR,     fabs(antiBjets.at(0).Eta()-bjets.at(0).Eta())); 
  fill("deltaEta_jfwd_b1_SRB"+SR,   fabs(forwardLightjets.at(0).Eta()-bjets.at(0).Eta())); 

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

