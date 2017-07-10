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


void findHiggs(const AnalysisObjects& bjets, TLorentzVector& h2, TLorentzVector &h3){
  // Make vector of b-jet pairs
  std::vector<std::pair<AnalysisObject, AnalysisObject>> bPairs;
  for (unsigned int i=0; i < bjets.size(); ++i) {
    for (unsigned int j=i+1; j < bjets.size(); ++j) {
      bPairs.push_back({bjets.at(i),bjets.at(j)});
    } 
  } 
  
  // Method 1 - closest to Higgs mass
  //sortPairsClosestToHiggs(bPairs);
  //h1 = bPairs.at(0).first+bPairs.at(0).second;
  
  // Metohd 2 - leading b-jets
  h2 = bjets.at(0) + bjets.at(1);
  
  // Method 3 - closest dR
  sortPairsClosestdR(bPairs);
  h3 = bPairs.at(0).first+bPairs.at(0).second;
  
}

TLorentzVector findTop(const AnalysisObjects& bjets, const AnalysisObjects& leptons, const AnalysisObject& metVec){
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


std::vector<TLorentzVector> findbjets_notH(AnalysisObjects bjets){
  std::vector<TLorentzVector> bjets_notH; 
  double smallDR=100; 
  int nBjets = bjets.size();
  if(nBjets<2) return bjets_notH; 
  int bjet1=-1;
  int bjet2=-1;
  for(int J1=0;J1<nBjets;J1++){
    for(int J2=0;J2<nBjets;J2++){
      if(J2==J1) continue; 
      double DR_bjets=bjets.at(J1).DeltaR(bjets.at(J2));
      if(smallDR>DR_bjets){
	smallDR=DR_bjets;
  	bjet1=J1;
  	bjet2=J2;
      }
    }
  }
  for(int i=0;i<nBjets;i++) if(i!=bjet1&&i!=bjet2) bjets_notH.push_back(bjets.at(i)); 
  return bjets_notH; 
}

std::vector<TLorentzVector> findHiggsRecoil(std::vector<TLorentzVector> bjets_notH, AnalysisObjects leptons, AnalysisObject metVec){
  std::vector<TLorentzVector> higgsRecoil; 
  TLorentzVector HR1;
  TLorentzVector HR2;
  HR1=leptons.at(0)+metVec; 
  HR2=leptons.at(0); 
  for(unsigned int i=0;i<bjets_notH.size();i++){ 
    HR1=HR1+bjets_notH.at(i);
    HR2=HR2+bjets_notH.at(i); 
  }
  higgsRecoil.push_back(HR1); 
  higgsRecoil.push_back(HR2); 
  return higgsRecoil;
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
  addHistogram("HT",100,0,2000); 
 
  addHistogram("lep_pt",100,0,500);
  addHistogram("lep_eta",25,0,4);   

  addHistogram("j1_pt_allJets",100,0,500);
  addHistogram("j1_eta_allJets",25,0,4); 

  addHistogram("jfwd_pt",100,0,500);
  addHistogram("jfwd_eta",25,0,4);
  addHistogram("j1_pt",100,0,500); 
  addHistogram("j1_eta",25,0,4); 

  addHistogram("bfwd_pt",100,0,500);
  addHistogram("bfwd_eta",25,0,4);
  addHistogram("b1_pt",100,0,500); 
  addHistogram("b1_eta",25,0,4); 
  
  addHistogram("dEta_jfwd_bfwd",100,0,8); //|eta(fwdjet)-eta(fwdbjet)|
  addHistogram("dEta_jfwd_b1",100,0,8);
  addHistogram("dEta_j1_b1",100,0,8); //|eta(leadjet)-eta(leadbjet)|
  addHistogram("dR_b1_b2",100,0,8);//dR(leading bjet, subleading bjet)

  addHistogram("top_m",100,0,500); //reco mass of top quark
  addHistogram("top_pt",100,0,500); 
  addHistogram("top_eta",25,0,4);

  for(int j=0;j<3;j++){ //Higgs reconstruction using 3 different methods
    std::string h_idx = std::to_string(j+1); 
    addHistogram("H"+h_idx+"_m",100,0,500); 
    addHistogram("H"+h_idx+"_pt",100,0,500);
    addHistogram("H"+h_idx+"_eta",25,0,4);
    addHistogram("dEta_H"+h_idx+"_jfwd",100,0,8);//|eta(higgs)-eta(forwardJet)|
  }

  //btagged regions
  for(int i=2;i<5;i++){

    std::string SR=std::to_string(i); 

    addHistogram("MET_SRB"+SR,100,0,500);
    addHistogram("Njets_SRB"+SR,10,-0.5,9.5);
    addHistogram("NLightjets_SRB"+SR,10,-0.5,9.5);
    addHistogram("HT_SRB"+SR,100,0,2000); 

    addHistogram("lep_pt_SRB"+SR,100,0,500); 
    addHistogram("lep_eta_SRB"+SR,25,0,4); 

    addHistogram("j1_pt_allJets_SRB"+SR,100,0,500); 
    addHistogram("j1_eta_allJets_SRB"+SR,25,0,4);

    addHistogram("jfwd_pt_SRB"+SR,100,0,500);
    addHistogram("jfwd_eta_SRB"+SR,25,0,4);
    addHistogram("j1_pt_SRB"+SR,100,0,500); 
    addHistogram("j1_eta_SRB"+SR,25,0,4);
    
    addHistogram("bfwd_pt_SRB"+SR,100,0,500);
    addHistogram("bfwd_eta_SRB"+SR,25,0,4);     
    addHistogram("b1_pt_SRB"+SR,100,0,500); 
    addHistogram("b1_eta_SRB"+SR,25,0,4);  
    
    addHistogram("dEta_jfwd_bfwd_SRB"+SR,100,0,8); //|eta(fwdjet)-eta(fwdbjet)|
    addHistogram("dEta_j1_b1_SRB"+SR,100,0,8); //|eta(leadjet)-eta(leadbjet)|
    addHistogram("dEta_jfwd_b1_SRB"+SR,100,0,8);

    addHistogram("dR_b1_b2_SRB"+SR,100,0,8);//dR(leading bjet, subleading bjet)

    addHistogram("top_m_SRB"+SR,100,0,500); //reco mass of top quark
    addHistogram("top_pt_SRB"+SR,100,0,500); 
    addHistogram("top_eta_SRB"+SR,25,0,4); 
    
    for(int j=0;j<3;j++){ //higgs reconstruction using 3 different methods
      std::string h_idx = std::to_string(j+1); 
      addHistogram("H"+h_idx+"_m_SRB"+SR,100,0,500);
      addHistogram("H"+h_idx+"_pt_SRB"+SR,100,0,500);
      addHistogram("H"+h_idx+"_eta_SRB"+SR,25,0,4);
    }

    addHistogram("dEta_H2_jfwd_SRB"+SR,100,0,8);//|eta(higgs)-eta(forwardJet)|
    addHistogram("dEta_H3_jfwd_SRB"+SR,100,0,8);//|eta(higgs)-eta(forwardJet)|

    ////Histograms after mbb cut
    addHistogram("jfwd_pt_SRMbbH2_SRB"+SR,100,0,500);
    addHistogram("jfwd_pt_SRMbbH3_SRB"+SR,100,0,500);

    addHistogram("jfwd_eta_SRMbbH2_SRB"+SR,25,0,4);
    addHistogram("jfwd_eta_SRMbbH3_SRB"+SR,25,0,4);

    addHistogram("dEta_jfwd_b1_SRMbbH2_SRB"+SR,100,0,8);
    addHistogram("dEta_jfwd_b1_SRMbbH3_SRB"+SR,100,0,8);

    addHistogram("dEta_H2_jfwd_SRMbbH2_SRB"+SR,100,0,8);
    addHistogram("dEta_H3_jfwd_SRMbbH3_SRB"+SR,100,0,8);
    
    addHistogram("H2_pt_SRMbbH2_SRB"+SR,100,0,500);
    addHistogram("H3_pt_SRMbbH3_SRB"+SR,100,0,500);

    addHistogram("dR_lep_H3_SRMbbH3_SRB"+SR,100,0,8);
    addHistogram("dR_lep_b1_SRMbbH3_SRB"+SR,100,0,8);
    addHistogram("dR_lep_b1NotH_SRMbbH3_SRB"+SR,100,0,8);
    addHistogram("dR_H3_b1NotH_SRMbbH3_SRB"+SR,100,0,8);
    addHistogram("dR_R1_H3_SRMbbH3_SRB"+SR,100,0,8);//R1=recoil of higgs w/ MET
    addHistogram("dR_R2_H3_SRMbbH3_SRB"+SR,100,0,8);//R2=recoil of higgs w/out MET
    
    addHistogram("dEta_lep_H3_SRMbbH3_SRB"+SR,100,0,8);
    addHistogram("dEta_lep_b1_SRMbbH3_SRB"+SR, 100,0,8);
    addHistogram("dEta_lep_b1NotH_SRMbbH3_SRB"+SR,100,0,8);
    addHistogram("dEta_H3_b1NotH_SRMbbH3_SRB"+SR,100,0,8);
    addHistogram("dEta_R1_H3_SRMbbH3_SRB"+SR,100,0,8);
    addHistogram("dEta_R2_H3_SRMbbH3_SRB"+SR,100,0,8);
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
  auto muons_noPtCut     = event->getMuons(0., 2.7, MuLoose && MuIsoGradient);
  auto jets_noPtCut      = event->getJets(0., 3.8); 

  auto electrons = event->getElectrons(25., 2.47,ELooseBLLH && EIsoGradient); //add vertex  (ED0Sigma5|EZ05mm ?)
  auto muons     = event->getMuons(25., 2.7, MuLoose && MuIsoGradient);//add vertex (MuD0Sigma3|MuZ05mm ?)
  auto jets      = event->getJets(25., 3.8, JVT50Jet); // what about NOT(LooseBadJet)
  auto metVec    = event->getMET();
  double met     = metVec.Et();

  // sorting objects by pt
  sortObjectsByPt(electrons);
  sortObjectsByPt(muons);
  sortObjectsByPt(jets);
  
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
  auto bjets_notH = findbjets_notH(bjets); //bjets not associated with higgs3 - we should avoid looping again here over the jets

  // eta-sorted jets
  auto forwardJets = jets;
  sortObjectsByEta(forwardJets);
  auto forwardBjets = filterObjects(forwardJets, 25., 3.8, BTag85MV2c20);
  auto forwardLightjets = filterObjects(forwardJets, 25., 3.8, NOT(BTag85MV2c20));

  // Object counting
  int numSignalLeptons = leptons.size();  // Object lists are essentially std::vectors so .size() works
  int numSignalJets    = jets.size();
  int nBjets =           bjets.size();
  
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
  if (met < 20 || met>100) return; 
  if (antiBjets.size()==0 || antiBjets.size()>2) return;

  // Fill histogram after cuts
  fill("MET",met);
  fill("Njets",numSignalJets);
  fill("NLightjets",antiBjets.size()); 
  fill("Nbjets",nBjets);
  fill("lep_pt",leptons.at(0).Pt());
  fill("lep_eta",fabs(leptons.at(0).Eta()));  

  ///Leading jet
  fill("j1_pt_allJets", jets.at(0).Pt());
  fill("j1_eta_allJets", fabs(jets.at(0).Eta()));

  ///fwd light jet
  fill("jfwd_pt",   forwardLightjets.at(0).Pt()); 
  fill("jfwd_eta",  fabs(forwardLightjets.at(0).Eta())); 

  ///lead light jet
  fill("j1_pt",  antiBjets.at(0).Pt());
  fill("j1_eta", fabs(antiBjets.at(0).Eta()));

  ///HT
  double pTSum=leptons.at(0).Pt();
  for(int i=0;i<numSignalJets;i++) pTSum+=jets.at(i).Pt();
  fill("HT",pTSum);


  
  ////______________bjet cut_____: 
  // Note some quantities below are not defined without this cut.
  ///Note if the histogram is not defined the code will crash and not write the output
  if( nBjets<2 || nBjets>4 ) return;   
  std::string SR=std::to_string(nBjets); 
  
  ///fwd b-jet
  fill("bfwd_pt",   forwardBjets.at(0).Pt()); 
  fill("bfwd_eta",  fabs(forwardBjets.at(0).Eta())); 

  ///lead b-jet
  fill("b1_pt",  bjets.at(0).Pt());
  fill("b1_eta", fabs(bjets.at(0).Eta()));

  //Delta R between 2 leading bjets
  fill("dR_b1_b2", bjets.at(0).DeltaR(bjets.at(1)));
  
  //DeltaEta 
  fill("dEta_jfwd_bfwd", fabs(forwardLightjets.at(0).Eta()-forwardBjets.at(0).Eta())); 
  fill("dEta_j1_b1",     fabs(antiBjets.at(0).Eta()-bjets.at(0).Eta())); 
  fill("dEta_jfwd_b1",   fabs(forwardLightjets.at(0).Eta()-bjets.at(0).Eta())); 

  //Reconstructed Higgs
  TLorentzVector higgs2, higgs3;
  findHiggs(bjets, higgs2, higgs3); 

  fill("H2_m",higgs2.M());
  fill("H3_m",higgs3.M());

  fill("H2_pt",higgs2.Pt());
  fill("H3_pt",higgs3.Pt());
    
  fill("H2_eta",fabs(higgs2.Eta()));
  fill("H3_eta",fabs(higgs3.Eta()));

  fill("dEta_H2_jfwd",fabs(higgs2.Eta()-forwardLightjets.at(0).Eta()));
  fill("dEta_H3_jfwd",fabs(higgs3.Eta()-forwardLightjets.at(0).Eta()));

  //Higgs recoil
  TLorentzVector R1=findHiggsRecoil(bjets_notH,leptons,metVec).at(0); //non-Higgs (H3) bjets+lepton+MET
  TLorentzVector R2=findHiggsRecoil(bjets_notH,leptons,metVec).at(1); //non-Higgs (H3) bjets+lepton
  
  //Reconstructed top quark
  TLorentzVector top=findTop(bjets,leptons,metVec); 
  fill("top_m", top.M()); 
  fill("top_pt", top.Pt());
  fill("top_eta", top.Eta());   


  ///////////////______________Fill histos corresponding to b-Tag signal regions____________
  fill("MET_SRB"+SR,            met);
  fill("Njets_SRB"+SR,          numSignalJets);
  fill("NLightjets_SRB"+SR,     antiBjets.size());  
  fill("HT_SRB"+SR,             pTSum); 

  fill("top_m_SRB"+SR,          top.M()); 
  fill("top_pt_SRB"+SR,         top.Pt()); 
  fill("top_eta_SRB"+SR,        fabs(top.Eta()));   

  fill("lep_pt_SRB"+SR,         leptons.at(0).Pt());
  fill("lep_eta_SRB"+SR,        fabs(leptons.at(0).Eta())); 

  fill("j1_pt_allJets_SRB"+SR,  jets.at(0).Pt());
  fill("j1_eta_allJets_SRB"+SR, fabs(jets.at(0).Eta()));

  fill("j1_pt_SRB"+SR,          antiBjets.at(0).Pt());
  fill("j1_eta_SRB"+SR,         fabs(antiBjets.at(0).Eta()));  

  fill("b1_pt_SRB"+SR,          bjets.at(0).Pt());
  fill("b1_eta_SRB"+SR,         fabs(bjets.at(0).Eta()));
  
  fill("jfwd_pt_SRB"+SR,        forwardLightjets.at(0).Pt());
  fill("jfwd_eta_SRB"+SR,       fabs(forwardLightjets.at(0).Eta())); 

  fill("bfwd_pt_SRB"+SR,        forwardBjets.at(0).Pt()); 
  fill("bfwd_eta_SRB"+SR,       fabs(forwardBjets.at(0).Eta()));   			     

  fill("dR_b1_b2_SRB"+SR,       bjets.at(0).DeltaR(bjets.at(1)));
  fill("dEta_jfwd_bfwd_SRB"+SR, fabs(forwardLightjets.at(0).Eta()-forwardBjets.at(0).Eta())); 
  fill("dEta_j1_b1_SRB"+SR,     fabs(antiBjets.at(0).Eta()-bjets.at(0).Eta())); 
  fill("dEta_jfwd_b1_SRB"+SR,   fabs(forwardLightjets.at(0).Eta()-bjets.at(0).Eta())); 

  fill("H2_m_SRB"+SR,           higgs2.M()); 
  fill("H2_pt_SRB"+SR,          higgs2.Pt()); 
  fill("H2_eta_SRB"+SR,         fabs(higgs2.Eta())); 
			        
  fill("H3_m_SRB"+SR,           higgs3.M()); 
  fill("H3_pt_SRB"+SR,          higgs3.Pt()); 
  fill("H3_eta_SRB"+SR,         fabs(higgs3.Eta())); 

  fill("dEta_H2_jfwd_SRB"+SR,   fabs(higgs2.Eta()-forwardLightjets.at(0).Eta()));
  fill("dEta_H3_jfwd_SRB"+SR,   fabs(higgs3.Eta()-forwardLightjets.at(0).Eta()));


  ////Histos with mbb cut 
  if(80<higgs2.M()&&higgs2.M()<130){
    fill("H2_pt_SRMbbH2_SRB"+SR,         higgs2.Pt());   
    fill("jfwd_pt_SRMbbH2_SRB"+SR,       forwardLightjets.at(0).Pt());
    fill("jfwd_eta_SRMbbH2_SRB"+SR,      forwardLightjets.at(0).Eta());
    fill("dEta_jfwd_b1_SRMbbH2_SRB"+SR,  fabs(forwardLightjets.at(0).Eta()-bjets.at(0).Eta())); 
    fill("dEta_H2_jfwd_SRMbbH2_SRB"+SR,  fabs(higgs2.Eta()-forwardLightjets.at(0).Eta()));
  }
  if(80<higgs3.M()&&higgs3.M()<130){
    fill("H3_pt_SRMbbH3_SRB"+SR,         higgs3.Pt());
    fill("jfwd_pt_SRMbbH3_SRB"+SR,       forwardLightjets.at(0).Pt());
    fill("jfwd_eta_SRMbbH3_SRB"+SR,      forwardLightjets.at(0).Eta());
    fill("dEta_jfwd_b1_SRMbbH3_SRB"+SR,  fabs(forwardLightjets.at(0).Eta()-bjets.at(0).Eta())); 
    fill("dEta_H3_jfwd_SRMbbH3_SRB"+SR,  fabs(higgs3.Eta()-forwardLightjets.at(0).Eta()));
                    
    fill("dR_lep_H3_SRMbbH3_SRB"+SR,     leptons.at(0).DeltaR(higgs3));
    fill("dR_lep_b1_SRMbbH3_SRB"+SR,     leptons.at(0).DeltaR(bjets.at(0)));
    fill("dEta_lep_H3_SRMbbH3_SRB"+SR,     fabs(leptons.at(0).Eta()-higgs3.Eta()));
    fill("dEta_lep_b1_SRMbbH3_SRB"+SR,     fabs(leptons.at(0).Eta()-bjets.at(0).Eta()));

    if(bjets_notH.size()>0){
      fill("dR_lep_b1NotH_SRMbbH3_SRB"+SR,   leptons.at(0).DeltaR(bjets_notH.at(0)));
      fill("dR_H3_b1NotH_SRMbbH3_SRB"+SR,    higgs3.DeltaR(bjets_notH.at(0)));
      fill("dR_R1_H3_SRMbbH3_SRB"+SR,        R1.DeltaR(higgs3));
      fill("dR_R2_H3_SRMbbH3_SRB"+SR,        R2.DeltaR(higgs3));
      fill("dEta_lep_b1NotH_SRMbbH3_SRB"+SR, fabs(leptons.at(0).Eta()-bjets_notH.at(0).Eta()));
      fill("dEta_H3_b1NotH_SRMbbH3_SRB"+SR,  fabs(higgs3.Eta()-bjets_notH.at(0).Eta()));
      fill("dEta_R1_H3_SRMbbH3_SRB"+SR,      fabs(R1.Eta()-higgs3.Eta()));
      fill("dEta_R2_H3_SRMbbH3_SRB"+SR,      fabs(R2.Eta()-higgs3.Eta()));          
    }          



    ///fill ntuple
    ntupVar("met",met);
    ntupVar("top_m",top.M()); 
    ntupVar("h_pt",higgs3.Pt());
    ntupVar("lep_pt",leptons.at(0).Pt());
    ntupVar("b1_pt",bjets.at(0).Pt());
  }                  
                    
  return;            
}                   
                    
