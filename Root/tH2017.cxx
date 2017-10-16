#include "SimpleAnalysis/AnalysisClass.h"
#include <string>
#include "xAODTruth/TruthParticleContainer.h"

float deltaR(float eta1, float phi1, float eta2, float phi2){
  float deta= fabs(eta1 - eta2);      
  float dphi= fabs(phi1 - phi2);
  if (dphi > TMath::Pi() ) dphi = 2*(TMath::Pi()) - dphi;
  return sqrt((dphi*dphi)+(deta*deta));
}

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

float getMinDRPairs(std::vector<std::pair<AnalysisObject,AnalysisObject>>& cands) {
  std::sort(cands.begin(), cands.end(), dR_sort());
  return (cands.size() > 0) ? cands.at(0).first.DeltaR(cands.at(0).second) : -1.;
}

float getSingleElectronTriggerEfficiency(float ptMeV, float eta) {
  float minPt = 22000.;
  float minPtHighEta = 35000;
  float maxEta = 4.0;
  float eff = 0.95;
  float effHighEta = 0.90;

  // HGTD forward trigger 5 GeV improvement
  minPtHighEta = 30000.;

  if ( ptMeV > 35000. && fabs(eta) < 2.5 ) return 1.0;
  if ( ptMeV > minPt && fabs(eta) < 2.5 )
    return eff;
  if ( ptMeV > minPtHighEta && fabs(eta) < maxEta )
    return effHighEta;
  return 0.0;
}

float muonEtaTriggerEfficiency(float eta) {
    // rpc L0 efficiency data  22 bins for 0<eta<1.1
    const float eta_bin = 0.05;
    const float eff_gold[22] = {0.790656, 0.930483, 0.98033, 0.992508, 0.974555, 0.981241, 0.985142, 0.947444, 0.960144, 0.98223, 0.983938, 0.984972, 0.972907, 0.982902, 0.919753, 0.899409, 0.970952, 0.960322, 0.946016, 0.868755, 0.619748,0};
    //=======    
    float eff = 0.98*0.98; //TGC*MDT efficiency
    if (fabs(eta)>2.65) return 0.;
    
    // RPC efficiencies for gold layout
    if (fabs(eta)<=1.05) {
      int ibin=fabs(eta)/eta_bin;
      eff=eff_gold[ibin]; //RPC recovery with BI RPC chambers
    }
    
    return eff;
}

float getSingleMuonTriggerEfficiency(float etMeV, float eta) {
  //single-mu trigger efficiency w.r.t. reconstruction efficiency (tight=true)
  //using 2012 values from K. Nagano
  float minPt=20000.;
  if (etMeV > minPt) return muonEtaTriggerEfficiency(eta);
  return 0.;
}


float minmaxdEta(const AnalysisObjects& ljets, const AnalysisObjects& bjets, const bool min) {
  std::vector<float> dist;
  for (auto ljet : ljets) {
    for (auto bjet : bjets) {
      dist.push_back(fabs(ljet.Eta()-bjet.Eta()));
    }
  }
  if (min) return *std::min_element(dist.begin(), dist.end());
  else return *std::max_element(dist.begin(), dist.end());
}

void findHiggs(const AnalysisObjects& bjets, TLorentzVector& h2, TLorentzVector &h3, AnalysisObjects & hbjets){
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
  
  // Method 2 - leading b-jets
  h2 = bjets.at(0) + bjets.at(1);
  
  // Method 3 - closest dR
  sortPairsClosestdR(bPairs);
  h3 = bPairs.at(0).first+bPairs.at(0).second;
  hbjets.push_back(bPairs.at(0).first); 
  hbjets.push_back(bPairs.at(0).second); 
  
}

std::vector<double> calculateFoxWMoments(const AnalysisObjects& jets, const AnalysisObjects& leptons){
  std::vector<double> FoxWMoments;
  AnalysisObjects allObj = jets; 
  if (leptons.size() > 0) allObj.push_back(leptons.at(0));
  double fw0=-1;
  double fw1=-1;
  double fw2=-1;
  double fw3=-1;
  
  // for normalization
  double pi = 0.0;
  
  if (allObj.size() > 0) {
    for (unsigned int i=0; i < allObj.size(); ++i) {
      pi += sqrt(pow(allObj.at(i).Px(),2)+ pow(allObj.at(i).Py(),2) + pow(allObj.at(i).Pz(),2));
      for (unsigned int j=0; j < allObj.size(); ++j) {
        // 3-momenta product
        double p1       = sqrt(pow(allObj.at(i).Px(),2)+ pow(allObj.at(i).Py(),2) + pow(allObj.at(i).Pz(),2));
        double p2       = sqrt(pow(allObj.at(j).Px(),2)+ pow(allObj.at(j).Py(),2) + pow(allObj.at(j).Pz(),2));
        double weight   = p1*p2;
        double cosOmega = cos(allObj.at(i).Theta())*cos(allObj.at(j).Theta()) + 
        		       sin(allObj.at(i).Theta())*sin(allObj.at(j).Theta())*cos(allObj.at(i).Phi()-allObj.at(j).Phi());
        fw0 += weight*1;
        fw1 += weight*cosOmega; 
        fw2 += weight*0.5*(3.*pow(cosOmega,2) - 1.); 
        fw3 += weight*0.5*(5.*pow(cosOmega,3) - 3.*cosOmega); 
      } 
    }
  
    // Normalization
    fw0 /= pow(pi,2);
    fw1 /= pow(pi,2);
    fw2 /= pow(pi,2);
    fw3 /= pow(pi,2);
  }
    
  FoxWMoments = {fw0,fw1,fw2,fw3};
  return FoxWMoments; 
}

void findTop(const AnalysisObjects& bjets, const AnalysisObjects& leptons, const AnalysisObject& metVec, AnalysisObject& top1, AnalysisObject &top2){
  double DR=100; 
  int nBjets=bjets.size(); 
  if(nBjets==0) return; 
  int top_bjet_idx=-1; 
  for(int iJet=0;iJet<nBjets;iJet++){
    double DR_bjetLep=bjets.at(iJet).DeltaR(leptons.at(0));
    if(DR>DR_bjetLep){
      DR=DR_bjetLep;
      top_bjet_idx=iJet; 
    }
  }
  //Method 1 - closest bjet to lepton
  if(top_bjet_idx!=-1) top1=bjets.at(top_bjet_idx)+leptons.at(0)+metVec; 
  //Method 2 - leading bjet
  top2=bjets.at(0)+leptons.at(0)+metVec; 
}

void findHiggsRecoil(AnalysisObjects& bjets_notH, AnalysisObjects& leptons, AnalysisObject& metVec, AnalysisObjects& higgsRecoil){
  auto HR1=leptons.at(0)+metVec; 
  auto HR2=leptons.at(0); 
  for(unsigned int i=0;i<bjets_notH.size();i++){ 
    HR1=HR1+bjets_notH.at(i);
    HR2=HR2+bjets_notH.at(i); 
  }
  higgsRecoil.push_back(HR1); 
  higgsRecoil.push_back(HR2); 
}

void printTruthTree(const xAOD::TruthParticleContainer* Truth, unsigned int idx){
  if(!Truth || idx >= Truth->size()){
    std::cout<<"printTruthTree Truth = NULL or idx> Truth size"<<std::endl;
    return;
  }
  
  std::cout<<" TRUTH TREE: "<<Truth->at(idx)->pdgId()<<" status="<<Truth->at(idx)->status()<<" nChildren="<<Truth->at(idx)->nChildren()<<std::endl;
  for(unsigned int i=0;i< Truth->at(idx)->nChildren();i++){
    const xAOD::TruthParticle* part1 = Truth->at(idx)->child(i);
    if(part1==NULL){std::cout<<"printTruthTree part1 = NULL at i="<<i<<"  of nChildren="<<Truth->at(idx)->nChildren()<<std::endl; return;  }
    ///First generation
    std::cout<<"    id="<<part1->pdgId()<<"  status="<<part1->status()<<"  (pt,eta)="<<part1->p4().Pt()<<","<<part1->p4().Eta()<<")"<<std::endl;
    for(unsigned int j=0;j< part1->nChildren();j++){
      const xAOD::TruthParticle* part2 = part1->child(j);
      if(part2==NULL){std::cout<<"printTruthTree part2 = NULL at i="<<j<<"  of nChildren="<<part1->nChildren()<<std::endl; return;  }
      ///Second generation
      std::cout<<"        id="<<part2->pdgId()<<"  status="<<part2->status()<<"  (pt,eta)="<<part2->p4().Pt()<<","<<part2->p4().Eta()<<")"<<std::endl;
      for(unsigned int k=0;k< part2->nChildren();k++){
        const xAOD::TruthParticle* part3 = part2->child(k);
	if(part3==NULL){std::cout<<"printTruthTree part3 = NULL at i="<<j<<"  of nChildren="<<part2->nChildren()<<std::endl; return;  }
        ///Third generation
        std::cout<<"            id="<<part3->pdgId()<<"  status="<<part3->status()<<"  (pt,eta)="<<part3->p4().Pt()<<","<<part3->p4().Eta()<<")"<<std::endl;
	for(unsigned int l=0;l< part3->nChildren();l++){
	  const xAOD::TruthParticle* part4 = part3->child(l);
	  if(part4==NULL){std::cout<<"printTruthTree part4 = NULL at i="<<l<<"  of nChildren="<<part3->nChildren()<<std::endl; return;  }
	  ///Fourth generation
	  std::cout<<"                id="<<part4->pdgId()<<"  status="<<part4->status()<<"  (pt,eta)="<<part4->p4().Pt()<<","<<part4->p4().Eta()<<")"<<std::endl;
	  for(unsigned int m=0;m< part4->nChildren();m++){
	    const xAOD::TruthParticle* part5 = part4->child(m);
	    if(part5==NULL){std::cout<<"printTruthTree part5 = NULL at i="<<m<<"  of nChildren="<<part4->nChildren()<<std::endl; return;  }
	    ///Fith generation
	    std::cout<<"                    id="<<part5->pdgId()<<"  status="<<part5->status()<<"  (pt,eta)="<<part5->p4().Pt()<<","<<part5->p4().Eta()<<")"<<std::endl;
	  }
	}
      }
    }
  }

}


bool findHbb(const xAOD::TruthParticle* parent, const xAOD::TruthParticle* &Hb, const xAOD::TruthParticle* &Hbbar){
  if(!parent){
    std::cout<<"findHbb parent = NULL"<<std::endl;
    return false;
  }
  
  //check if there is a b-bbar pair
  for(unsigned int i=0;i<parent->nChildren();i++){
    const xAOD::TruthParticle* dau = parent->child(i);
    if(dau==NULL){std::cout<<"findHbb dau = NULL at i="<<i<<"  of nChildren="<<parent->nChildren()<<std::endl; return false;  }
    if(dau->pdgId()==-5) Hb = dau;
    if(dau->pdgId()==5)  Hbbar = dau;
  }
  if(Hb!=NULL && Hbbar!=NULL) return true; //exit the recursive loop

  //continue to check the daughters
  for(unsigned int i=0;i<parent->nChildren();i++){
    if(findHbb(parent->child(i),Hb,Hbbar)) return true;
  }

  //std::cout<<"findHbb(): failed to find H->bb "<<std::endl;
  return false;
}
  
bool findTopWb(const xAOD::TruthParticle* parent, const xAOD::TruthParticle* &TopW, const xAOD::TruthParticle* &TopB){
  if(!parent){
    std::cout<<"findTopWb parent = NULL"<<std::endl;
    return false;
  }
  
  //check if there is a b-bbar pair
  for(unsigned int i=0;i<parent->nChildren();i++){
    const xAOD::TruthParticle* dau = parent->child(i);
    if(dau==NULL){std::cout<<"findTopWb dau = NULL at i="<<i<<"  of nChildren="<<parent->nChildren()<<std::endl; return false;  }
    if(abs(dau->pdgId())==24) TopW = dau;
    if(abs(dau->pdgId())==5)  TopB = dau;
  }
  if(TopW!=NULL && TopB!=NULL) return true; //exit the recursive loop

  //continue to check the daughters
  for(unsigned int i=0;i<parent->nChildren();i++){
    if(findTopWb(parent->child(i),TopW,TopB)) return true;
  }

  //std::cout<<"findTopWb(): failed to find Top->Wb "<<std::endl;
  return false;
}
  

bool findWl(const xAOD::TruthParticle* parent, const xAOD::TruthParticle* &Wl){
  if(!parent){
    std::cout<<"findWl parent = NULL"<<std::endl;
    return false;
  }
  
  //check if there is leptonic decay
  for(unsigned int i=0;i<parent->nChildren();i++){
    const xAOD::TruthParticle* dau = parent->child(i);
    if(dau==NULL){std::cout<<"findWl dau = NULL at i="<<i<<"  of nChildren="<<parent->nChildren()<<std::endl; return false;  }
    if((abs(dau->pdgId())==11 || abs(dau->pdgId())==13) && dau->status()==1 ) Wl = dau;
  }
  if(Wl!=NULL) return true; //exit the recursive loop

  //continue to check the daughters
  for(unsigned int i=0;i<parent->nChildren();i++){
    if(findWl(parent->child(i),Wl)) return true;
  }

  //std::cout<<"findWl(): failed to find truth W->l decay "<<std::endl;
  return false;
}
  



DefineAnalysis(tH2017)

void tH2017::Init()
{
  //addRegions({"SR1lep0b","SR1lep1b","SR1lep2b","SR1lep3b","SR1lep4b"});

  addHistogram("events",14,-0.5,13.5);

  //no cuts
  addHistogram("NLep_nocuts",10,-0.5,9.5);
  addHistogram("MET_nocuts",100,0,500);
  addHistogram("Njets_nocuts",10,-0.5,9.5);
  addHistogram("jet_pt_nocuts",100,0,500);
  addHistogram("lep_pt_nocuts",100,0,500); 
  addHistogram("sphericity_nocuts",50,0,1);

  addHistogram("Njets_HGTD_nocuts",10,-0.5,9.5); 
  addHistogram("Nbjets_HGTD_nocuts",6,-0.5,5.5); 
  
  //no pt cut on jets/leptons
  addHistogram("jet_pt_noPtCut",100,0,500); 
  addHistogram("lep_pt_noPtCut",100,0,500); 
  
  //lepton cut but no b-tag cuts
  addHistogram("alljettruth_pt_noBcuts",100,0,500); 
  addHistogram("alljettruth_eta_noBcuts",80,0,4);   
  addHistogram("alljettruthB_pt_noBcuts",100,0,500); 
  addHistogram("alljettruthB_eta_noBcuts",80,0,4);   
  addHistogram("alljettruthC_pt_noBcuts",100,0,500); 
  addHistogram("alljettruthC_eta_noBcuts",80,0,4);   
  addHistogram("alljettruthTau_pt_noBcuts",100,0,500); 
  addHistogram("alljettruthTau_eta_noBcuts",80,0,4);   
  addHistogram("alljettruthL_pt_noBcuts",100,0,500); 
  addHistogram("alljettruthL_eta_noBcuts",80,0,4);   


  //after preselection
  addHistogram("MET",100,0,500);
  addHistogram("Njets",10,-0.5,9.5);
  addHistogram("NLightjets",10,-0.5,9.5);
  addHistogram("Nbjets",6,-0.5,5.5);
  addHistogram("HT",100,0,2000); 
  addHistogram("sphericity",50,0,1);
 
  addHistogram("lep_pt",100,0,500);
  addHistogram("lep_eta",80,0,4);   

  addHistogram("j1_pt_allJets",100,0,500);
  addHistogram("j1_eta_allJets",80,0,4); 

  addHistogram("jfwd_pt",100,0,500);
  addHistogram("jfwd_eta",80,0,4);
  addHistogram("j1_pt",100,0,500); 
  addHistogram("j1_eta",80,0,4); 

  addHistogram("bfwd_pt",100,0,500);
  addHistogram("bfwd_eta",80,0,4);
  addHistogram("b1_pt",100,0,500); 
  addHistogram("b1_eta",80,0,4); 
  
  addHistogram("dEta_jfwd_bfwd",100,0,8); //|eta(fwdjet)-eta(fwdbjet)|
  addHistogram("dEta_jfwd_b1",100,0,8);
  addHistogram("dEta_j1_b1",100,0,8); //|eta(leadjet)-eta(leadbjet)|
  addHistogram("dR_b1_b2",100,0,8);//dR(leading bjet, subleading bjet)
  addHistogram("mindEta_ljets_bjets",100,0,8);
  addHistogram("maxdEta_ljets_bjets",100,0,8);

  addHistogram("top1_m",100,0,500); //reco mass of top quark
  addHistogram("top1_pt",100,0,500); 
  addHistogram("top1_eta",80,0,4);

  addHistogram("top2_m",100,0,500); 
  addHistogram("top2_pt",100,0,500); 
  addHistogram("top2_eta",80,0,4);

  addHistogram("top3_m",100,0,500); 
  addHistogram("top3_pt",100,0,500); 
  addHistogram("top3_eta",80,0,4);

  addHistogram("top4_m",100,0,500); 
  addHistogram("top4_pt",100,0,500); 
  addHistogram("top4_eta",80,0,4);

  addHistogram("H2_m",100,0,500); 
  addHistogram("H2_pt",100,0,500);
  addHistogram("H2_eta",80,0,4);

  addHistogram("H3_m",100,0,500); 
  addHistogram("H3_pt",100,0,500);
  addHistogram("H3_eta",80,0,4);

  addHistogram("dEta_H2_jfwd",100,0,8);//|eta(higgs)-eta(forwardJet)|
  addHistogram("dEta_H3_jfwd",100,0,8);//|eta(higgs)-eta(forwardJet)|

  addHistogram("FoxW0",50,0,1);
  addHistogram("FoxW1",50,0,1);
  addHistogram("FoxW2",50,0,1);
  addHistogram("FoxW3",50,0,1);

  //btagged regions
  for(int i=2;i<5;i++){

    std::string SR=std::to_string(i);  
    std::vector<std::string> SRjfwd_opt = {"", "_SRjfwd"};

    for(auto SRjfwd: SRjfwd_opt){

      addHistogram("MET"+SRjfwd+"_SRB"+SR,100,0,500);
      addHistogram("Njets"+SRjfwd+"_SRB"+SR,10,-0.5,9.5);
      addHistogram("NLightjets"+SRjfwd+"_SRB"+SR,10,-0.5,9.5);
      addHistogram("HT"+SRjfwd+"_SRB"+SR,100,0,2000); 

      addHistogram("lep_pt"+SRjfwd+"_SRB"+SR,100,0,500); 
      addHistogram("lep_eta"+SRjfwd+"_SRB"+SR,80,0,4); 

      addHistogram("j1_pt_allJets"+SRjfwd+"_SRB"+SR,100,0,500); 
      addHistogram("j1_eta_allJets"+SRjfwd+"_SRB"+SR,80,0,4);

      addHistogram("jfwd_pt"+SRjfwd+"_SRB"+SR,100,0,500);
      addHistogram("jfwd_eta"+SRjfwd+"_SRB"+SR,80,0,4);
      addHistogram("j1_pt"+SRjfwd+"_SRB"+SR,100,0,500); 
      addHistogram("j1_eta"+SRjfwd+"_SRB"+SR,80,0,4);
    
      addHistogram("bfwd_pt"+SRjfwd+"_SRB"+SR,100,0,500);
      addHistogram("bfwd_eta"+SRjfwd+"_SRB"+SR,80,0,4);     
      addHistogram("b1_pt"+SRjfwd+"_SRB"+SR,100,0,500); 
      addHistogram("b1_eta"+SRjfwd+"_SRB"+SR,80,0,4);  
    
      addHistogram("dEta_jfwd_bfwd"+SRjfwd+"_SRB"+SR,100,0,8); //|eta(fwdjet)-eta(fwdbjet)|
      addHistogram("dEta_j1_b1"+SRjfwd+"_SRB"+SR,100,0,8); //|eta(leadjet)-eta(leadbjet)|
      addHistogram("dEta_jfwd_b1"+SRjfwd+"_SRB"+SR,100,0,8);

      addHistogram("dR_b1_b2"+SRjfwd+"_SRB"+SR,100,0,8);//dR(leading bjet, subleading bjet)

      addHistogram("top1_m"+SRjfwd+"_SRB"+SR,100,0,500); //reco mass of top quark
      addHistogram("top1_pt"+SRjfwd+"_SRB"+SR,100,0,500); 
      addHistogram("top1_eta"+SRjfwd+"_SRB"+SR,80,0,4); 

      addHistogram("top2_m"+SRjfwd+"_SRB"+SR,100,0,500); //reco mass of top quark
      addHistogram("top2_pt"+SRjfwd+"_SRB"+SR,100,0,500); 
      addHistogram("top2_eta"+SRjfwd+"_SRB"+SR,80,0,4); 

      addHistogram("top3_m"+SRjfwd+"_SRB"+SR,100,0,500); //reco mass of top quark
      addHistogram("top3_pt"+SRjfwd+"_SRB"+SR,100,0,500); 
      addHistogram("top3_eta"+SRjfwd+"_SRB"+SR,80,0,4); 

      addHistogram("top4_m"+SRjfwd+"_SRB"+SR,100,0,500); //reco mass of top quark
      addHistogram("top4_pt"+SRjfwd+"_SRB"+SR,100,0,500); 
      addHistogram("top4_eta"+SRjfwd+"_SRB"+SR,80,0,4); 
    
      addHistogram("H2_m"+SRjfwd+"_SRB"+SR,100,0,500);
      addHistogram("H2_pt"+SRjfwd+"_SRB"+SR,100,0,500);
      addHistogram("H2_eta"+SRjfwd+"_SRB"+SR,80,0,4);

      addHistogram("H3_m"+SRjfwd+"_SRB"+SR,100,0,500);
      addHistogram("H3_pt"+SRjfwd+"_SRB"+SR,100,0,500);
      addHistogram("H3_eta"+SRjfwd+"_SRB"+SR,80,0,4);

      addHistogram("dEta_H2_jfwd"+SRjfwd+"_SRB"+SR,100,0,8);//|eta(higgs)-eta(forwardJet)|
      addHistogram("dEta_H3_jfwd"+SRjfwd+"_SRB"+SR,100,0,8);//|eta(higgs)-eta(forwardJet)|

      addHistogram("mindEta_ljets_bjets"+SRjfwd+"_SRB"+SR,100,0,8);
      addHistogram("maxdEta_ljets_bjets"+SRjfwd+"_SRB"+SR,100,0,8);

      addHistogram("FoxW0"+SRjfwd+"_SRB"+SR,50,0,1);
      addHistogram("FoxW1"+SRjfwd+"_SRB"+SR,50,0,1);
      addHistogram("FoxW2"+SRjfwd+"_SRB"+SR,50,0,1);
      addHistogram("FoxW3"+SRjfwd+"_SRB"+SR,50,0,1);

      
      ///histogram for the reconstruction efficiencies
      // tHqb;  t -> W b -> l v b;  h-> b1 b2, 
      // bin 1 counts total events 
      // bin 2=t, 3=H, 4=q, 5=b;   6=top_l,7=top_b;   8=H_b1, 9=H_b2 
      addHistogram("effReco"+SRjfwd+"_SRB"+SR,9,0.5,9.5);
      
    }

    ////Histograms after mbb cut
    addHistogram("jfwd_pt_SRMbbH2_SRB"+SR,100,0,500);
    addHistogram("jfwd_pt_SRMbbH3_SRB"+SR,100,0,500);

    addHistogram("jfwd_eta_SRMbbH2_SRB"+SR,80,0,4);
    addHistogram("jfwd_eta_SRMbbH3_SRB"+SR,80,0,4);

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

    addHistogram("FoxW0_SRMbbH3_SRB"+SR,50,0,1);
    addHistogram("FoxW1_SRMbbH3_SRB"+SR,50,0,1);
    addHistogram("FoxW2_SRMbbH3_SRB"+SR,50,0,1);
    addHistogram("FoxW3_SRMbbH3_SRB"+SR,50,0,1);


    //b-tagged jet composition 
    addHistogram("allbjettruth_pt_SRB"+SR,100,0,500); 
    addHistogram("allbjettruth_eta_SRB"+SR,80,0,4);   
    addHistogram("allbjettruthB_pt_SRB"+SR,100,0,500); 
    addHistogram("allbjettruthB_eta_SRB"+SR,80,0,4);   
    addHistogram("allbjettruthC_pt_SRB"+SR,100,0,500); 
    addHistogram("allbjettruthC_eta_SRB"+SR,80,0,4);   
    addHistogram("allbjettruthTau_pt_SRB"+SR,100,0,500); 
    addHistogram("allbjettruthTau_eta_SRB"+SR,80,0,4);   
    addHistogram("allbjettruthL_pt_SRB"+SR,100,0,500); 
    addHistogram("allbjettruthL_eta_SRB"+SR,80,0,4);   

    //light jet composition
    addHistogram("lightjettruth_pt_SRB"+SR,100,0,500); 
    addHistogram("lightjettruth_eta_SRB"+SR,80,0,4);   
    addHistogram("lightjettruthB_pt_SRB"+SR,100,0,500); 
    addHistogram("lightjettruthB_eta_SRB"+SR,80,0,4);   
    addHistogram("lightjettruthC_pt_SRB"+SR,100,0,500); 
    addHistogram("lightjettruthC_eta_SRB"+SR,80,0,4);   
    addHistogram("lightjettruthTau_pt_SRB"+SR,100,0,500); 
    addHistogram("lightjettruthTau_eta_SRB"+SR,80,0,4);   
    addHistogram("lightjettruthL_pt_SRB"+SR,100,0,500); 
    addHistogram("lightjettruthL_eta_SRB"+SR,80,0,4);   

  }
}

void tH2017::ProcessEvent(AnalysisEvent *event)
{

  _output->setEventWeight(1);
  fill("events",0); //Number of processed events

  double eventweight=event->getMCWeights()[0];
  _output->setEventWeight(eventweight); //all histograms after this point get filled with this weight
  fill("events",1); //Sum of initial weights

  auto electrons_noPtCut = event->getElectrons(0., 4.0, ELooseBLLH && EIsoGradient);
  auto muons_noPtCut     = event->getMuons(0., 2.7, MuLoose && MuIsoGradient);
  auto jets_noPtCut      = event->getJets(0., 3.8); 

  auto electrons = event->getElectrons(25., 4.0, ELooseBLLH && EIsoGradient); //add vertex  (ED0Sigma5|EZ05mm ?)
  auto muons     = event->getMuons(25., 2.7, MuLoose && MuIsoGradient);//add vertex (MuD0Sigma3|MuZ05mm ?)
  auto jets      = event->getJets(25., 3.8, JVT50Jet); // what about NOT(LooseBadJet)
  auto metVec    = event->getMET();
  double met     = metVec.Et();

  // sorting objects by pt
  sortObjectsByPt(electrons);
  sortObjectsByPt(muons);
  sortObjectsByPt(jets);
    

  // Overlap removal - including with object Pt-dependent radius calculation
  int btagger= BTag70MV2c20; //BTag85MV2c20;
  electrons_noPtCut  = overlapRemoval(electrons_noPtCut, muons_noPtCut, 0.01);
  jets_noPtCut   = overlapRemoval(jets_noPtCut, electrons_noPtCut, 0.2, NOT(btagger));//check this cut
  electrons_noPtCut  = overlapRemoval(electrons_noPtCut, jets_noPtCut, 0.4);
  jets_noPtCut   = overlapRemoval(jets_noPtCut, muons_noPtCut, 0.4, LessThan3Tracks); //check this cut (b-jets?)
  muons_noPtCut      = overlapRemoval(muons_noPtCut, jets_noPtCut, 0.4);   

  electrons  = overlapRemoval(electrons, muons, 0.01);
  jets   = overlapRemoval(jets, electrons, 0.2, NOT(btagger));//check this cut
  electrons  = overlapRemoval(electrons, jets, 0.4);
  jets   = overlapRemoval(jets, muons, 0.4, LessThan3Tracks); //check this cut (b-jets?)
  muons      = overlapRemoval(muons, jets, 0.4);   

  // Lists of objects can be merged by simple addition
  auto leptons   = electrons + muons;
  auto leptons_noPtCut = electrons_noPtCut + muons_noPtCut; 

  //
  auto bjets = filterObjects(jets, 25., 3.8, btagger);
  auto antiBjets = filterObjects(jets, 25., 3.8, NOT(btagger));

  // eta-sorted jets
  auto forwardJets = jets;
  sortObjectsByEta(forwardJets);
  auto forwardBjets = filterObjects(forwardJets, 25., 3.8, btagger);
  auto forwardLightjets = filterObjects(forwardJets, 25., 3.8, NOT(btagger));

  // Object counting
  int numSignalLeptons = leptons.size();  // Object lists are essentially std::vectors so .size() works
  int numSignalJets    = jets.size();
  int nBjets           = bjets.size();
    
  // Jets in HGTD acceptance
  int njets_hgtd(0), nbjets_hgtd(0);
  for (auto j : jets) {
    if ( fabs(j.Eta()) >= 2.4 && fabs(j.Eta()) <= 3.8) {
      ++njets_hgtd;
      if (j.pass(GoodBJet)) ++nbjets_hgtd;
    }
  }
  fill("Njets_HGTD_nocuts", njets_hgtd);
  fill("Nbjets_HGTD_nocuts", nbjets_hgtd);

  // Fill in histograms without cuts
  fill("NLep_nocuts",numSignalLeptons);
  fill("MET_nocuts",met);
  fill("Njets_nocuts",numSignalJets);

  // Sphericity 
  float sph = sphericity(jets);
  fill("sphericity_nocuts", sph); 

  if(numSignalLeptons>0) fill("lep_pt_nocuts", leptons.at(0).Pt()); 
  if(leptons_noPtCut.size()>0) fill("lep_pt_noPtCut",leptons_noPtCut.at(0).Pt());
  
  for(int iJet=0;iJet<numSignalJets;iJet++) fill("jet_pt_nocuts", jets.at(iJet).Pt()); 
  for(unsigned int iJet=0;iJet<jets_noPtCut.size();iJet++) fill("jet_pt_noPtCut", jets_noPtCut.at(iJet).Pt()); 
  
  ///HT
  double pTSum=0.0;
  if (leptons.size() > 0) pTSum+=leptons.at(0).Pt();
  for(int i=0;i<numSignalJets;i++) pTSum+=jets.at(i).Pt();


  /////////////////////////////////////////////
  //truth particles only available in TRUTH1
  const xAOD::TruthParticleContainer* Truth = event->getTruthParticles();  
  const xAOD::TruthParticle* truthFwdJet = Truth?Truth->at(7):NULL;
  const xAOD::TruthParticle* truthB =      Truth?Truth->at(6):NULL;
  const xAOD::TruthParticle* truthHiggs =  Truth?Truth->at(4):NULL;  // printTruthTree(Truth,4);
  const xAOD::TruthParticle* truthHB = NULL;
  const xAOD::TruthParticle* truthHBbar = NULL;
  if(truthHiggs) findHbb(truthHiggs,truthHB,truthHBbar);  
  const xAOD::TruthParticle* truthHB1 = NULL;
  const xAOD::TruthParticle* truthHB2 = NULL;
  if(truthHB && truthHBbar) {
    if(truthHB->pt() > truthHBbar->pt()){ truthHB1 = truthHB;  truthHB2 = truthHBbar;}
    if(truthHB->pt() < truthHBbar->pt()){ truthHB1 = truthHBbar;  truthHB2 = truthHB;}
  }
  const xAOD::TruthParticle* truthTop =  Truth?Truth->at(5):NULL; //printTruthTree(Truth,5);
  const xAOD::TruthParticle* truthTopB = NULL;
  const xAOD::TruthParticle* truthTopW = NULL;
  if(truthTop) findTopWb(truthTop,truthTopW,truthTopB); 
  const xAOD::TruthParticle* truthTopWL = NULL;
  if(truthTopW) findWl(truthTopW,truthTopWL);

  const  xAOD::TruthParticle* truthFwdB = truthB; //most forward b quark
  if(truthFwdB && truthTopB && fabs(truthTopB->eta()) > fabs(truthFwdB->eta())) truthFwdB = truthTopB;
  if(truthFwdB && truthHB1 && fabs(truthHB1->eta()) > fabs(truthFwdB->eta())) truthFwdB = truthHB1;
  if(truthFwdB && truthHB2 && fabs(truthHB2->eta()) > fabs(truthFwdB->eta())) truthFwdB = truthHB2;

  ///fill tuple truth info
  if(truthTopWL)ntupVar("truth_lep_pt",truthTopWL->pt());
  if(truthTopWL)ntupVar("truth_lep_eta",truthTopWL->eta());
  if(truthFwdJet)ntupVar("truth_jfwd_pt",truthFwdJet->pt());
  if(truthFwdJet)ntupVar("truth_jfwd_eta",truthFwdJet->eta());
  if(truthHB)ntupVar("truth_b_pt",truthB->pt());
  if(truthHB)ntupVar("truth_b_eta",truthB->eta());
  if(truthHB1)ntupVar("truth_Hb1_pt",truthHB1->pt());
  if(truthHB1)ntupVar("truth_Hb1_eta",truthHB1->eta());
  if(truthHB2)ntupVar("truth_Hb2_pt",truthHB2->pt());
  if(truthHB2)ntupVar("truth_Hb2_eta",truthHB2->eta());
  if(truthTopB)ntupVar("truth_Topb_pt",truthTopB->pt());
  if(truthTopB)ntupVar("truth_Topb_eta",truthTopB->eta());
  if(truthFwdB)ntupVar("truth_Fwdb_pt",truthFwdB->pt());
  if(truthFwdB)ntupVar("truth_Fwdb_eta",truthFwdB->eta());


  //inclusive jet histograms without b-tagging
  if (numSignalLeptons ==1 && pTSum > 300 && numSignalJets >= 3){
    for (auto j : jets) {
      fill("alljettruth_eta_noBcuts",fabs(j.Eta())); 
      fill("alljettruth_pt_noBcuts",fabs(j.Pt())); 

      if (j.pass(TrueBJet)){
	fill("alljettruthB_eta_noBcuts",fabs(j.Eta())); 
	fill("alljettruthB_pt_noBcuts",fabs(j.Pt())); 
      }else if (j.pass(TrueCJet))  {
	fill("alljettruthC_eta_noBcuts",fabs(j.Eta())); 
	fill("alljettruthC_pt_noBcuts",fabs(j.Pt())); 
      }else if (j.pass(TrueTau)) {
	fill("alljettruthTau_eta_noBcuts",fabs(j.Eta())); 
	fill("alljettruthTau_pt_noBcuts",fabs(j.Pt())); 
      }else {
	fill("alljettruthL_eta_noBcuts",fabs(j.Eta())); 
	fill("alljettruthL_pt_noBcuts",fabs(j.Pt())); 
      }     
    }
  }
 

  // Preselection - decision
  bool passPreSel = true;
  if (numSignalLeptons != 1) passPreSel = false;
  fill("events",2);
  if (numSignalJets < 4 || numSignalJets > 5) passPreSel = false; 
  fill("events",3);
  if (antiBjets.size() < 1 || antiBjets.size() > 3) passPreSel = false;
  fill("events",4);
  if (nBjets < 2 ) passPreSel = false;
  fill("events",5);
  if (pTSum < 300 ) passPreSel = false;
  fill("events",6);
  if (antiBjets.size() > 0) if (fabs(forwardLightjets.at(0).Eta()) < 2.4) passPreSel = false;

  // Trigger
  float trigWeight = 1.0;
  if ( electrons.size() == 1 ) {
    trigWeight = getSingleElectronTriggerEfficiency(electrons.at(0).Pt()*1000., electrons.at(0).Eta());
  }
  else if ( muons.size() == 1 ) {
    trigWeight = getSingleMuonTriggerEfficiency(muons.at(0).Pt()*1000., muons.at(0).Eta());
  }

  // Apply trigger weight
  _output->setEventWeight(eventweight*trigWeight);

  // Calculate FW moments
  std::vector<double> fwMoments = calculateFoxWMoments(jets, leptons); 


  //-------------------------
  // Extra variables for BDT
  //--------
  
  // mindR bjets
  std::vector<std::pair<AnalysisObject, AnalysisObject>> bPairs;
  for (unsigned int i=0; i < bjets.size(); ++i) {
    for (unsigned int j=i+1; j < bjets.size(); ++j) {
      bPairs.push_back({bjets.at(i),bjets.at(j)});
    } 
  } 
  float mindRbb = getMinDRPairs(bPairs);
  
  // mindRLepB
  std::vector<std::pair<AnalysisObject, AnalysisObject>> blPairs;
  float mindRbl = -1.0;
  if (numSignalLeptons > 0) {
    for (unsigned int i=0; i < bjets.size(); ++i) {
      blPairs.push_back({leptons.at(0),bjets.at(i)});
    }
    mindRbl = getMinDRPairs(blPairs);
  }
  //-------------------------
  
  ///fill ntuple
  ntupVar("numSignalLeptons", numSignalLeptons);
  ntupVar("numSignalJets", numSignalJets);
  ntupVar("Nbjets", (float)nBjets);
  ntupVar("Nljets", (float)antiBjets.size());
  ntupVar("pTSum", pTSum);
  ntupVar("met", met); 
  ntupVar("FW3", fwMoments.at(3)); 
  ntupVar("FW1", fwMoments.at(1)); 
  ntupVar("mindRbb", mindRbb);
  ntupVar("mindRbl", mindRbl);
  if (numSignalLeptons > 0) {
    ntupVar("lep_pt",leptons.at(0).Pt());
    ntupVar("lep_eta",leptons.at(0).Eta());    
  }
  if ( forwardLightjets.size() > 0 ) {
    ntupVar("mostFwdLJet_eta", forwardLightjets.at(0).Eta());
    ntupVar("mostFwdLJet_pt", forwardLightjets.at(0).Pt());    
    if (bjets.size() > 0) ntupVar("dEta_mostFwdLJet_LeadBjet", fabs(forwardLightjets.at(0).Eta() - bjets.at(0).Eta()));
  }
  if ( bjets.size() > 0 ) {
    ntupVar("bjet1_pt", bjets.at(0).Pt()); 
  }
  ntupVar("sphericity", sph);

  // Preselection application
  if (!passPreSel) return;

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

  // HT
  fill("HT",pTSum);

  //Shape variables
  fill("FoxW0", fwMoments.at(0)); 
  fill("FoxW1", fwMoments.at(1)); 
  fill("FoxW2", fwMoments.at(2)); 
  fill("FoxW3", fwMoments.at(3)); 

  fill("sphericity",sph); 

  ////______________bjet cut_____: 
  // Note some quantities below are not defined without this cut.
  ///Note if the histogram is not defined the code will crash and not write the output
  //if( nBjets>4 ) return; // this might be unnecessary now
  std::string SR=std::to_string(nBjets); 
  

  for (auto j : bjets) {
    fill("allbjettruth_eta_SRB"+SR,fabs(j.Eta())); 
    fill("allbjettruth_pt_SRB"+SR,fabs(j.Pt()));    
    if (j.pass(TrueBJet)){
      fill("allbjettruthB_eta_SRB"+SR,fabs(j.Eta())); 
      fill("allbjettruthB_pt_SRB"+SR,fabs(j.Pt())); 
    }else if (j.pass(TrueCJet))  {
      fill("allbjettruthC_eta_SRB"+SR,fabs(j.Eta())); 
      fill("allbjettruthC_pt_SRB"+SR,fabs(j.Pt())); 
    }else if (j.pass(TrueTau)) {
      fill("allbjettruthTau_eta_SRB"+SR,fabs(j.Eta())); 
      fill("allbjettruthTau_pt_SRB"+SR,fabs(j.Pt())); 
    }else {
      fill("allbjettruthL_eta_SRB"+SR,fabs(j.Eta())); 
      fill("allbjettruthL_pt_SRB"+SR,fabs(j.Pt())); 
    }     
  }


  fill("lightjettruth_eta_SRB"+SR,fabs(forwardLightjets.at(0).Eta())); 
  fill("lightjettruth_pt_SRB"+SR,fabs(forwardLightjets.at(0).Pt())); 
  if (forwardLightjets.at(0).pass(TrueBJet)){
    fill("lightjettruthB_eta_SRB"+SR,fabs(forwardLightjets.at(0).Eta())); 
    fill("lightjettruthB_pt_SRB"+SR,fabs(forwardLightjets.at(0).Pt())); 
  }else if (forwardLightjets.at(0).pass(TrueCJet))  {
    fill("lightjettruthC_eta_SRB"+SR,fabs(forwardLightjets.at(0).Eta())); 
    fill("lightjettruthC_pt_SRB"+SR,fabs(forwardLightjets.at(0).Pt())); 
  }else if (forwardLightjets.at(0).pass(TrueTau)) {
    fill("lightjettruthTau_eta_SRB"+SR,fabs(forwardLightjets.at(0).Eta())); 
    fill("lightjettruthTau_pt_SRB"+SR,fabs(forwardLightjets.at(0).Pt())); 
  }else {
    fill("lightjettruthL_eta_SRB"+SR,fabs(forwardLightjets.at(0).Eta())); 
    fill("lightjettruthL_pt_SRB"+SR,fabs(forwardLightjets.at(0).Pt())); 
  } 

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
  fill("mindEta_ljets_bjets", minmaxdEta(antiBjets, bjets, true)); 
  fill("maxdEta_ljets_bjets", minmaxdEta(antiBjets, bjets, false)); 
  				  
  //Reconstructed Higgs
  TLorentzVector higgs2, higgs3;
  AnalysisObjects hbjets;
  findHiggs(bjets, higgs2, higgs3, hbjets); 

  AnalysisObjects bjets_notH;  //bjets not associated with higgs3
  for(auto i:bjets) if(hbjets.at(0).DeltaR(i)>0.1&&hbjets.at(1).DeltaR(i)>0.1) bjets_notH.push_back(i); 

  fill("H2_m",higgs2.M());
  fill("H3_m",higgs3.M());

  fill("H2_pt",higgs2.Pt());
  fill("H3_pt",higgs3.Pt());
    
  fill("H2_eta",fabs(higgs2.Eta()));
  fill("H3_eta",fabs(higgs3.Eta()));

  fill("dEta_H2_jfwd",fabs(higgs2.Eta()-forwardLightjets.at(0).Eta()));
  fill("dEta_H3_jfwd",fabs(higgs3.Eta()-forwardLightjets.at(0).Eta()));

  //Higgs recoil
  AnalysisObjects higgsRecoil; 
  findHiggsRecoil(bjets_notH,leptons,metVec,higgsRecoil); //non-Higgs (H3) bjets+lepton+MET
  auto R1 = higgsRecoil.at(0); 
  auto R2 = higgsRecoil.at(1); 
  
  //Reconstructed top quark
  AnalysisObject top1(0.,0.,0.,0.,0,0,COMBINED,0); 
  AnalysisObject top2(0.,0.,0.,0.,0,0,COMBINED,0); 
  AnalysisObject top3(0.,0.,0.,0.,0,0,COMBINED,0); 
  AnalysisObject top4(0.,0.,0.,0.,0,0,COMBINED,0); 
  
  findTop(bjets,leptons,metVec,top1,top2); 
  findTop(bjets_notH,leptons,metVec,top3,top4); 
  
  fill("top1_m", top1.M()); 
  fill("top1_pt", top1.Pt());
  fill("top1_eta", top1.Eta());   

  fill("top2_m", top2.M()); 
  fill("top2_pt", top2.Pt());
  fill("top2_eta", top2.Eta());   

  if(bjets_notH.size()>1){
    fill("top3_m", top3.M()); 
    fill("top3_pt", top3.Pt());
    fill("top3_eta", top3.Eta());   

    fill("top4_m", top4.M()); 
    fill("top4_pt", top4.Pt());
    fill("top4_eta", top4.Eta());
  }   



  ///////////////______________Fill histos corresponding to b-Tag signal regions____________
  if(nBjets==3) fill("events",7);
  if(nBjets==4) fill("events",8);

  //////////  Reco efficiencies
  fill("effReco_SRB"+SR,1);//count total events in this category



  //////check the top
  if(truthTopWL) 
    if(deltaR(truthTopWL->eta(),truthTopWL->phi(),leptons.at(0).Eta(),leptons.at(0).Phi())<0.2) 
      fill("effReco_SRB"+SR,6);
  if(truthTopB)
    if(deltaR(truthTopB->eta(),truthTopB->phi(),bjets.at(0).Eta(),bjets.at(0).Phi())<0.3) ///(use leading b-jet for now)
      fill("effReco_SRB"+SR,7);
  if(truthTopWL && truthTopB)
    if(deltaR(truthTopWL->eta(),truthTopWL->phi(),leptons.at(0).Eta(),leptons.at(0).Phi())<0.2
       && deltaR(truthTopB->eta(),truthTopB->phi(),bjets.at(0).Eta(),bjets.at(0).Phi())<0.3) 
      fill("effReco_SRB"+SR,2);

  //////check the H 
  if(truthHB && truthHBbar) 
    if(deltaR(truthHB->eta(),truthHB->phi(),hbjets.at(0).Eta(),hbjets.at(0).Phi())<0.3
       || deltaR(truthHBbar->eta(),truthHBbar->phi(),hbjets.at(0).Eta(),hbjets.at(0).Phi())<0.3) 
      fill("effReco_SRB"+SR,8);
  if(truthHB && truthHBbar) 
    if(deltaR(truthHB->eta(),truthHB->phi(),hbjets.at(1).Eta(),hbjets.at(1).Phi())<0.3
       || deltaR(truthHBbar->eta(),truthHBbar->phi(),hbjets.at(1).Eta(),hbjets.at(1).Phi())<0.3)
      fill("effReco_SRB"+SR,9);
  if(truthHB && truthHBbar) 
    if((deltaR(truthHB->eta(),truthHB->phi(),hbjets.at(0).Eta(),hbjets.at(0).Phi())<0.3
  	&& deltaR(truthHBbar->eta(),truthHBbar->phi(),hbjets.at(1).Eta(),hbjets.at(1).Phi())<0.3)
       ||(deltaR(truthHB->eta(),truthHB->phi(),hbjets.at(1).Eta(),hbjets.at(1).Phi())<0.3
  	  && deltaR(truthHBbar->eta(),truthHBbar->phi(),hbjets.at(0).Eta(),hbjets.at(0).Phi())<0.3) ) 
      fill("effReco_SRB"+SR,3);
  
  //check the q
  if(truthFwdJet) 
    if(deltaR(truthFwdJet->eta(),truthFwdJet->phi(),forwardLightjets.at(0).Eta(),forwardLightjets.at(0).Phi())<0.3) 
      fill("effReco_SRB"+SR,4);
  
  //check the b
  if(truthB) 
    if(deltaR(truthB->eta(),truthB->phi(),bjets.at(0).Eta(),bjets.at(0).Phi())<0.3) ///(use leading b-jet for now)
      fill("effReco_SRB"+SR,5);


  fill("MET_SRB"+SR,            met);
  fill("Njets_SRB"+SR,          numSignalJets);
  fill("NLightjets_SRB"+SR,     antiBjets.size());  
  fill("HT_SRB"+SR,             pTSum); 

  fill("lep_pt_SRB"+SR,         leptons.at(0).Pt());
  fill("lep_eta_SRB"+SR,        fabs(leptons.at(0).Eta())); 

  fill("j1_pt_allJets_SRB"+SR,  jets.at(0).Pt());
  fill("j1_eta_allJets_SRB"+SR, fabs(jets.at(0).Eta()));

  fill("jfwd_pt_SRB"+SR,        forwardLightjets.at(0).Pt());
  fill("jfwd_eta_SRB"+SR,       fabs(forwardLightjets.at(0).Eta())); 
  fill("j1_pt_SRB"+SR,          antiBjets.at(0).Pt());
  fill("j1_eta_SRB"+SR,         fabs(antiBjets.at(0).Eta()));  

  fill("bfwd_pt_SRB"+SR,        forwardBjets.at(0).Pt()); 
  fill("bfwd_eta_SRB"+SR,       fabs(forwardBjets.at(0).Eta()));   			     
  fill("b1_pt_SRB"+SR,          bjets.at(0).Pt());
  fill("b1_eta_SRB"+SR,         fabs(bjets.at(0).Eta()));
  
  fill("dEta_jfwd_bfwd_SRB"+SR, fabs(forwardLightjets.at(0).Eta()-forwardBjets.at(0).Eta())); 
  fill("dEta_j1_b1_SRB"+SR,     fabs(antiBjets.at(0).Eta()-bjets.at(0).Eta())); 
  fill("dEta_jfwd_b1_SRB"+SR,   fabs(forwardLightjets.at(0).Eta()-bjets.at(0).Eta())); 
  
  fill("dR_b1_b2_SRB"+SR,       bjets.at(0).DeltaR(bjets.at(1)));

  fill("top1_m_SRB"+SR,          top1.M()); 
  fill("top1_pt_SRB"+SR,         top1.Pt()); 
  fill("top1_eta_SRB"+SR,        fabs(top1.Eta()));   

  fill("top2_m_SRB"+SR,          top2.M()); 
  fill("top2_pt_SRB"+SR,         top2.Pt()); 
  fill("top2_eta_SRB"+SR,        fabs(top2.Eta()));   

  if(bjets_notH.size()>1){
    fill("top3_m_SRB"+SR,          top3.M()); 
    fill("top3_pt_SRB"+SR,         top3.Pt()); 
    fill("top3_eta_SRB"+SR,        fabs(top3.Eta()));   

    fill("top4_m_SRB"+SR,          top4.M()); 
    fill("top4_pt_SRB"+SR,         top4.Pt()); 
    fill("top4_eta_SRB"+SR,        fabs(top4.Eta()));   
  }

  fill("H2_m_SRB"+SR,           higgs2.M()); 
  fill("H2_pt_SRB"+SR,          higgs2.Pt()); 
  fill("H2_eta_SRB"+SR,         fabs(higgs2.Eta())); 
			        
  fill("H3_m_SRB"+SR,           higgs3.M()); 
  fill("H3_pt_SRB"+SR,          higgs3.Pt()); 
  fill("H3_eta_SRB"+SR,         fabs(higgs3.Eta())); 

  fill("dEta_H2_jfwd_SRB"+SR,   fabs(higgs2.Eta()-forwardLightjets.at(0).Eta()));
  fill("dEta_H3_jfwd_SRB"+SR,   fabs(higgs3.Eta()-forwardLightjets.at(0).Eta()));
  fill("mindEta_ljets_bjets_SRB"+SR, minmaxdEta(antiBjets, bjets, true)); 
  fill("maxdEta_ljets_bjets_SRB"+SR, minmaxdEta(antiBjets, bjets, false)); 

  fill("FoxW0_SRB"+SR, fwMoments.at(0)); 
  fill("FoxW1_SRB"+SR, fwMoments.at(1)); 
  fill("FoxW2_SRB"+SR, fwMoments.at(2)); 
  fill("FoxW3_SRB"+SR, fwMoments.at(3)); 

  ////Histos with fwdJet cut  
  if(2.5<fabs(forwardLightjets.at(0).Eta())&&fabs(forwardLightjets.at(0).Eta())<3.8){
    if(nBjets==3) fill("events",12);
    if(nBjets==4) fill("events",13);

    fill("MET_SRjfwd_SRB"+SR,            met);
    fill("Njets_SRjfwd_SRB"+SR,          numSignalJets);
    fill("NLightjets_SRjfwd_SRB"+SR,     antiBjets.size());  
    fill("HT_SRjfwd_SRB"+SR,             pTSum); 

    fill("lep_pt_SRjfwd_SRB"+SR,         leptons.at(0).Pt());
    fill("lep_eta_SRjfwd_SRB"+SR,        fabs(leptons.at(0).Eta())); 

    fill("j1_pt_allJets_SRjfwd_SRB"+SR,  jets.at(0).Pt());
    fill("j1_eta_allJets_SRjfwd_SRB"+SR, fabs(jets.at(0).Eta()));

    fill("jfwd_pt_SRjfwd_SRB"+SR,        forwardLightjets.at(0).Pt());
    fill("jfwd_eta_SRjfwd_SRB"+SR,       fabs(forwardLightjets.at(0).Eta()));   
    fill("j1_pt_SRjfwd_SRB"+SR,          antiBjets.at(0).Pt());
    fill("j1_eta_SRjfwd_SRB"+SR,         fabs(antiBjets.at(0).Eta()));  

    fill("bfwd_pt_SRjfwd_SRB"+SR,        forwardBjets.at(0).Pt()); 
    fill("bfwd_eta_SRjfwd_SRB"+SR,       fabs(forwardBjets.at(0).Eta()));   			     
    fill("b1_pt_SRjfwd_SRB"+SR,          bjets.at(0).Pt());
    fill("b1_eta_SRjfwd_SRB"+SR,         fabs(bjets.at(0).Eta()));
  
    fill("dEta_jfwd_bfwd_SRjfwd_SRB"+SR, fabs(forwardLightjets.at(0).Eta()-forwardBjets.at(0).Eta())); 
    fill("dEta_j1_b1_SRjfwd_SRB"+SR,     fabs(antiBjets.at(0).Eta()-bjets.at(0).Eta())); 
    fill("dEta_jfwd_b1_SRjfwd_SRB"+SR,   fabs(forwardLightjets.at(0).Eta()-bjets.at(0).Eta())); 
  
    fill("dR_b1_b2_SRjfwd_SRB"+SR,       bjets.at(0).DeltaR(bjets.at(1)));

    fill("top1_m_SRjfwd_SRB"+SR,         top1.M()); 
    fill("top1_pt_SRjfwd_SRB"+SR,        top1.Pt()); 
    fill("top1_eta_SRjfwd_SRB"+SR,       fabs(top1.Eta()));   

    fill("top2_m_SRjfwd_SRB"+SR,         top2.M()); 
    fill("top2_pt_SRjfwd_SRB"+SR,        top2.Pt()); 
    fill("top2_eta_SRjfwd_SRB"+SR,       fabs(top2.Eta()));   

    if(bjets_notH.size()>1){
      fill("top3_m_SRjfwd_SRB"+SR,       top3.M()); 
      fill("top3_pt_SRjfwd_SRB"+SR,      top3.Pt()); 
      fill("top3_eta_SRjfwd_SRB"+SR,     fabs(top3.Eta()));   

      fill("top4_m_SRjfwd_SRB"+SR,       top4.M()); 
      fill("top4_pt_SRjfwd_SRB"+SR,      top4.Pt()); 
      fill("top4_eta_SRjfwd_SRB"+SR,     fabs(top4.Eta()));   
    }

    fill("H2_m_SRjfwd_SRB"+SR,           higgs2.M()); 
    fill("H2_pt_SRjfwd_SRB"+SR,          higgs2.Pt()); 
    fill("H2_eta_SRjfwd_SRB"+SR,         fabs(higgs2.Eta())); 
			        
    fill("H3_m_SRjfwd_SRB"+SR,           higgs3.M()); 
    fill("H3_pt_SRjfwd_SRB"+SR,          higgs3.Pt()); 
    fill("H3_eta_SRjfwd_SRB"+SR,         fabs(higgs3.Eta())); 

    fill("dEta_H2_jfwd_SRjfwd_SRB"+SR,   fabs(higgs2.Eta()-forwardLightjets.at(0).Eta()));
    fill("dEta_H3_jfwd_SRjfwd_SRB"+SR,   fabs(higgs3.Eta()-forwardLightjets.at(0).Eta()));
    fill("mindEta_ljets_bjets_SRjfwd_SRB"+SR, minmaxdEta(antiBjets, bjets, true)); 
    fill("maxdEta_ljets_bjets_SRjfwd_SRB"+SR, minmaxdEta(antiBjets, bjets, false)); 

    fill("FoxW0_SRjfwd_SRB"+SR, fwMoments.at(0)); 
    fill("FoxW1_SRjfwd_SRB"+SR, fwMoments.at(1)); 
    fill("FoxW2_SRjfwd_SRB"+SR, fwMoments.at(2)); 
    fill("FoxW3_SRjfwd_SRB"+SR, fwMoments.at(3)); 
  }

  ////Histos with mbb cut 
  if(30<higgs2.M()&&higgs2.M()<130){
    fill("events",9);
    if(nBjets==3) fill("events",10);
    if(nBjets==4) fill("events",11);
    fill("H2_pt_SRMbbH2_SRB"+SR,         higgs2.Pt());   
    fill("jfwd_pt_SRMbbH2_SRB"+SR,       forwardLightjets.at(0).Pt());
    fill("jfwd_eta_SRMbbH2_SRB"+SR,      forwardLightjets.at(0).Eta());
    fill("dEta_jfwd_b1_SRMbbH2_SRB"+SR,  fabs(forwardLightjets.at(0).Eta()-bjets.at(0).Eta())); 
    fill("dEta_H2_jfwd_SRMbbH2_SRB"+SR,  fabs(higgs2.Eta()-forwardLightjets.at(0).Eta()));
  }
  if(30<higgs3.M()&&higgs3.M()<130){ 
    fill("H3_pt_SRMbbH3_SRB"+SR,         higgs3.Pt());
    fill("jfwd_pt_SRMbbH3_SRB"+SR,       forwardLightjets.at(0).Pt());
    fill("jfwd_eta_SRMbbH3_SRB"+SR,      forwardLightjets.at(0).Eta());
    fill("dEta_jfwd_b1_SRMbbH3_SRB"+SR,  fabs(forwardLightjets.at(0).Eta()-bjets.at(0).Eta())); 
    fill("dEta_H3_jfwd_SRMbbH3_SRB"+SR,  fabs(higgs3.Eta()-forwardLightjets.at(0).Eta()));
                    
    fill("dR_lep_H3_SRMbbH3_SRB"+SR,     leptons.at(0).DeltaR(higgs3));
    fill("dR_lep_b1_SRMbbH3_SRB"+SR,     leptons.at(0).DeltaR(bjets.at(0)));
    fill("dEta_lep_H3_SRMbbH3_SRB"+SR,   fabs(leptons.at(0).Eta()-higgs3.Eta()));
    fill("dEta_lep_b1_SRMbbH3_SRB"+SR,   fabs(leptons.at(0).Eta()-bjets.at(0).Eta()));

    fill("dR_R1_H3_SRMbbH3_SRB"+SR,      R1.DeltaR(higgs3));
    fill("dR_R2_H3_SRMbbH3_SRB"+SR,      R2.DeltaR(higgs3));
    fill("dEta_R1_H3_SRMbbH3_SRB"+SR,    fabs(R1.Eta()-higgs3.Eta()));
    fill("dEta_R2_H3_SRMbbH3_SRB"+SR,    fabs(R2.Eta()-higgs3.Eta()));          

    if(bjets_notH.size()>0){
      fill("dR_lep_b1NotH_SRMbbH3_SRB"+SR,   leptons.at(0).DeltaR(bjets_notH.at(0)));
      fill("dR_H3_b1NotH_SRMbbH3_SRB"+SR,    higgs3.DeltaR(bjets_notH.at(0)));
      fill("dEta_lep_b1NotH_SRMbbH3_SRB"+SR, fabs(leptons.at(0).Eta()-bjets_notH.at(0).Eta()));
      fill("dEta_H3_b1NotH_SRMbbH3_SRB"+SR,  fabs(higgs3.Eta()-bjets_notH.at(0).Eta()));
    }     
    
    fill("FoxW1_SRMbbH3_SRB"+SR, fwMoments.at(0)); 
    fill("FoxW2_SRMbbH3_SRB"+SR, fwMoments.at(1)); 
    fill("FoxW3_SRMbbH3_SRB"+SR, fwMoments.at(2)); 
  }

  
  // Additional ntuple variables
  /*if(bjets_notH.size()>0){
    ntupVar("dR_lep_b1NotH",leptons.at(0).DeltaR(bjets_notH.at(0)));
    ntupVar("dR_H3_b1NotH",higgs3.DeltaR(bjets_notH.at(0)));
    ntupVar("dEta_lep_b1NotH",fabs(leptons.at(0).Eta()-bjets_notH.at(0).Eta()));
    ntupVar("dEta_H3_b1NotH",fabs(higgs3.Eta()-bjets_notH.at(0).Eta()));
  }*/                    
                 
  return;            
}                   
                    
