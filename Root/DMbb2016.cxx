#include "SimpleAnalysis/AnalysisClass.h"

#include "TMath.h"

DefineAnalysis(DMbb2016)
//
// DM+bb analysis (Run2 2015+2016 data)
// Support Doc: ATL-COM-PHYS-2016-1619
//              https://cds.cern.ch/record/2231810
// @author : M.Tripiana
//
void DMbb2016::Init()
{
  addRegions({"SRM_presel","SRM","SRA","SRB"});
  addRegions({"CRZM","CRWM","CRTM"});
  addRegions({"CRZA","CRWA","CRTA"});
  addRegions({"CRZB","CRWB","CRTB"});

  addRegions({"VRZM","VRWM","VRM"});
  addRegions({"VRA1","VRA2"});
  addRegions({"VRB1","VRB2"});

  // Book 1/2D histograms
  addHistogram("MET",40,0,2000);
  addHistogram("GenMET",40,0,2000);
  addHistogram("j1pt",20,0,1000);
  addHistogram("dphibb",20,0,3.5);

  addHistogram("SRMflow",20,0,20);
  addHistogram("SRAflow",20,0,20);
  addHistogram("SRBflow",20,0,20);

}

void DMbb2016::ProcessEvent(AnalysisEvent *event)
{
    
  auto electrons  = event->getElectrons(7, 2.47, ELooseLH);
  auto muons      = event->getMuons(6, 2.7, MuMedium);
  auto candJets   = event->getJets(20., 2.8);
  auto metVec     = event->getMET();
  double met      = metVec.Et();
    
  // Fill in histogram
  fill("SRMflow",0.5);
  fill("SRAflow",0.5);
  fill("SRBflow",0.5);


  if (countObjects(candJets, 20, 2.8, NOT(LooseBadJet))!=0) return;

  candJets = filterObjects(candJets, 20, 2.8, JVT50Jet);

  // ntupVar("j_pt1",candJets[0].Pt());
  // ntupVar("j_pt2",candJets[1].Pt());
  // ntupVar("j_eta1",candJets[0].Eta());
  // ntupVar("j_eta2",candJets[1].Eta());
 
  // Overlap removal
  // not doing overlap removal between electrons and muons with identical track
  candJets   = overlapRemoval(candJets,electrons,0.2,NOT(BTag85MV2c10));
  electrons  = overlapRemoval(electrons,candJets,0.4);
  candJets   = overlapRemoval(candJets, muons, 0.4, LessThan3Tracks); // CHECK: CONF is not clear that jet is removed
  muons      = overlapRemoval(muons, candJets, 0.4);

  auto baselineLeptons   = electrons + muons;

  auto signalJets      = filterObjects(candJets, 20);
  auto signalElectrons = filterObjects(electrons, 20, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoLooseTrack);
  auto signalMuons     = filterObjects(muons, 20, 2.5, MuD0Sigma3|MuZ05mm|MuIsoLooseTrack);
  auto signalLeptons   = signalElectrons + signalMuons;
  auto signalBjets     = filterObjects(signalJets, 20., 2.5, BTag60MV2c10);


  //MET with invisible leptons
  auto metVecCorr = metVec;
  for(const auto& slep : signalLeptons)
    metVecCorr += slep;
  float metCorr = metVecCorr.Pt();


  //Fill in optional ntuple variables
  ntupVar("met",met);
  ntupVar("gen_met", event->getGenMET());
  
  ntupVar("j_N",(int)candJets.size());
  ntupVar("e_N",(int)electrons.size());
  ntupVar("m_N",(int)muons.size());
  ntupVar("jets",candJets);

  //Common Selection 
  int nJets  = signalJets.size();
  int nBjets = signalBjets.size();


  //fill after OR
  fill("MET",met);
  fill("j1pt", nJets ? signalJets[0].Pt() : -999);
  fill("dphibb", nBjets>1 ? fabs(signalBjets[0].DeltaPhi(signalBjets[1])) : -999);

  if (nJets<2)  return;
  if (nJets>=4) return; 
  if (nBjets<1) return; 


  //  float dphiMin1  = minDphi(metVec, signalJets, 1);
  float dphiMin3  = minDphi(metVec, signalJets, 3);

  if(dphiMin3 < 0.4) return;
      
  if(signalBjets[0].Pt() <= 50.) return;

  fill("SRMflow",1.5);
  fill("SRAflow",1.5);
  fill("SRBflow",1.5);

  //dRmin
  float dRmin  = minDR(signalJets); 
  
  //Imbalance
  float Imb = (nBjets>1 ? (signalBjets[0].Pt()-signalBjets[1].Pt()) / (signalBjets[0].Pt()+signalBjets[1].Pt()) : -999);
  
  //dphibb
  float dphibb =  (nBjets>1 ? fabs(signalBjets[0].DeltaPhi(signalBjets[1])) : -999);

  //j1/HT
  float HT = sumObjectsPt(signalJets);
  float j1HT = signalJets[0].Pt() / HT;

  //float x1
  float x1 = dphiMin3 - dphibb;

  //float y1
  float y1 = TMath::Pi() - dphiMin3 - dphibb;

  //third jet veto 
  bool j3veto = (nJets==2 || signalJets[2].Pt() < 60.);

  //MT
  float mT = ((int)signalLeptons.size() ? calcMT(signalLeptons[0], metVec) : -999.);

  
  // Signal region definitions
  if (baselineLeptons.size()==0 ) { 

    fill("SRMflow",2.5);
    fill("SRAflow",2.5);
    fill("SRBflow",2.5);
    
    if(nBjets>=2){
      
      fill("SRMflow",3.5);
      fill("SRAflow",3.5);
      fill("SRBflow",3.5);

      if(met > 100.){      
	
	fill("SRMflow",4.5);
	fill("SRAflow",4.5);
	fill("SRBflow",4.5);
	
	accept("SRM_presel");


	if( nBjets==2 && signalJets[0].Pt() > 150. && j3veto && met > 180. && dRmin>2.8 && Imb > 0.5 )  accept("SRM");

	if( nBjets>=2 && signalJets[0].Pt() > 150. && j3veto && signalBjets[0].Pt() >150. && met > 180. && dphibb > 0.75 && j1HT > 0.4 && dphiMin3 > 2.2)  accept("SRA");

	if( nBjets>=2 && signalJets[0].Pt() > 150. && j3veto && signalBjets[0].Pt() >150. && met > 180. && Imb > 0.6 && dphibb > 0.5 && j1HT > 0.75 && x1 < 0 && fabs(y1) < 0.5)  accept("SRB");

	if( nBjets==2 && signalJets[0].Pt() > 150. && j3veto && met > 180. && dRmin>2.5 && (Imb > 0.3 && Imb < 0.5) )  accept("VRM");

	if( nBjets>=2 && signalJets[0].Pt() > 150. && j3veto && signalBjets[0].Pt() >150. && met > 180. && dphibb < 0.5 )  accept("VRA1");

	if( nBjets>=2 && signalJets[0].Pt() > 150. && j3veto && signalBjets[0].Pt() >150. && met > 180. && j1HT < 0.6 && dphiMin3 < 1.5 )  accept("VRA2");

	if( nBjets>=2 && signalJets[0].Pt() > 150. && j3veto && signalBjets[0].Pt() >150. && met > 180. && dphibb < 0.5 && j1HT > 0.75 )  accept("VRB1");

	if( nBjets>=2 && signalJets[0].Pt() > 150. && j3veto && signalBjets[0].Pt() >150. && met > 180. && x1 < 0 && fabs(y1) > 0.5 )  accept("VRB2");

      }
    }
  }


  // 1-lepton Control/Validation region definitions
  if (signalLeptons.size()==1 && baselineLeptons.size()==1 && signalLeptons[0].Pt()>30){
    
    if(signalJets[0].Pt()>150. && j3veto && nBjets==1 && met > 180. && dRmin>2.5 && (mT>30. && mT<100.) && Imb>0.5) accept("CRWM");

    //if(signalJets[0].Pt()>150. && j3veto && nBjets==2 && met > 150. && dRmin>2.5 && mT>30. ) accept("CRWT");

    if(signalJets[0].Pt()>150. && j3veto && nBjets==2 && met > 180. && (dRmin>1.25 && dRmin<2.5) && (mT > 30. && mT < 100.) && Imb > 0.5 ) accept("VRWM");

    if(signalJets[0].Pt()>150. && j3veto && nBjets==1 && signalBjets[0].Pt() >150. && met > 180. && j1HT > 0.8 && (mT>30. && mT<100.) && dphiMin3 > 2.2) accept("CRWA");

    if(signalJets[0].Pt()>150. && j3veto && nBjets>=2 && signalBjets[0].Pt() >150. && met > 180. && dphibb > 0.75 && j1HT > 0.4 && mT > 30. && dphiMin3 > 1.2) accept("CRTA");

    if(signalJets[0].Pt()>150. && j3veto && nBjets==1 && signalBjets[0].Pt() >150. && met > 180. && j1HT > 0.75 && (mT>30. && mT<100.) && dphibb > 0.5 && x1 < 1.5) accept("CRWB");

    if(signalJets[0].Pt()>150. && j3veto && nBjets>=2 && signalBjets[0].Pt() >150. && met > 180. && j1HT > 0.75 && mT > 30. && dphibb > 0.5 && Imb > 0.6  && x1 < 0 && fabs(y1) < 0.5) accept("CRTB");

  }

  // 2-lepton Control region definitions
  if (signalLeptons.size()==2 && baselineLeptons.size()==2 &&
      signalLeptons[0].type()==signalLeptons[1].type() &&
      signalLeptons[0].charge()!=signalLeptons[1].charge() &&
      signalLeptons[0].Pt()>30 && signalLeptons[0].Pt()>25 ) {
 
    float mll = (signalLeptons[0] + signalLeptons[1]).M();
    
    if( nBjets==1 && signalJets[0].Pt() > 150. && j3veto && met < 80. && metCorr > 180. && (mll>71. && mll<111.) && dRmin > 2.5 && Imb > 0.5 )  accept("CRZM");

    if( nBjets>=2 && signalJets[0].Pt() > 150. && j3veto && signalBjets[0].Pt() >150. && met < 60. && metCorr > 150. && (mll>71. && mll<111.) && dphibb > 0.75 )  accept("CRZA");

    if( nBjets>=2 && signalJets[0].Pt() > 150. && j3veto && signalBjets[0].Pt() >150. && met < 60. && metCorr > 180. && (mll>71. && mll<111.) && dphibb > 0.5 && j1HT > 0.5 )  accept("CRZB");

  }

  return;
}
