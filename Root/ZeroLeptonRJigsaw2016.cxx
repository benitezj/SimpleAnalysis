#include "SimpleAnalysis/AnalysisClass.h"
#ifdef ROOTCORE_PACKAGE_Ext_RestFrames

DefineAnalysis(ZeroLeptonRJigsaw2016)

void ZeroLeptonRJigsaw2016::Init()
{

  // RJR based SR's
  addRegions({"SR_RJR_S1a","SR_RJR_S1b","SR_RJR_S2a","SR_RJR_S2b","SR_RJR_S3a","SR_RJR_S3b","SR_RJR_S4"});
  addRegions({"SR_RJR_G1a","SR_RJR_G1b","SR_RJR_G2a","SR_RJR_G2b","SR_RJR_G3a","SR_RJR_G3b","SR_RJR_G4"});
  addRegions({"SR_RJR_C1","SR_RJR_C2","SR_RJR_C3","SR_RJR_C4","SR_RJR_C5"});


  //
  // Set up RestFrames trees through m_RF_Helper (RestFramesHelper class)
  //

  //
  // nominal squark and gluino tree
  //
  LabRecoFrame* LAB = m_RF_helper.addLabFrame("LAB");
  DecayRecoFrame* PP = m_RF_helper.addDecayFrame("PP");
  DecayRecoFrame* Pa = m_RF_helper.addDecayFrame("Pa");
  DecayRecoFrame* Pb = m_RF_helper.addDecayFrame("Pb");
  DecayRecoFrame* Ca = m_RF_helper.addDecayFrame("Ca");
  DecayRecoFrame* Cb = m_RF_helper.addDecayFrame("Cb");
  SelfAssemblingRecoFrame* SAV1a = m_RF_helper.addSAFrame("SAV1a");
  SelfAssemblingRecoFrame* SAV1b = m_RF_helper.addSAFrame("SAV1b");
  SelfAssemblingRecoFrame* SAV2a = m_RF_helper.addSAFrame("SAV2a");
  SelfAssemblingRecoFrame* SAV2b = m_RF_helper.addSAFrame("SAV2b");
  VisibleRecoFrame* V1a = m_RF_helper.addVisibleFrame("V1a");
  VisibleRecoFrame* V1b = m_RF_helper.addVisibleFrame("V1b");
  VisibleRecoFrame* V2a = m_RF_helper.addVisibleFrame("V2a");
  VisibleRecoFrame* V2b = m_RF_helper.addVisibleFrame("V2b");
  InvisibleRecoFrame* Ia = m_RF_helper.addInvisibleFrame("Ia");
  InvisibleRecoFrame* Ib = m_RF_helper.addInvisibleFrame("Ib");

  LAB->SetChildFrame(*PP);
  PP->AddChildFrame(*Pa);
  PP->AddChildFrame(*Pb);
  Pa->AddChildFrame(*SAV1a);
  Pb->AddChildFrame(*SAV1b);
  Pa->AddChildFrame(*Ca);
  Pb->AddChildFrame(*Cb);
  SAV1a->AddChildFrame(*V1a);
  SAV1b->AddChildFrame(*V1b);
  Ca->AddChildFrame(*SAV2a);
  Cb->AddChildFrame(*SAV2b);
  Ca->AddChildFrame(*Ia);
  Cb->AddChildFrame(*Ib);
  SAV2a->AddChildFrame(*V2a);
  SAV2b->AddChildFrame(*V2b);

  LAB->InitializeTree();

  InvisibleGroup* INV = m_RF_helper.addInvisibleGroup("INV");
  InvisibleJigsaw* InvMassJigsaw = m_RF_helper.addInvisibleJigsaw("InvMassJigsaw", kSetMass);
  InvisibleJigsaw* InvRapidityJigsaw = m_RF_helper.addInvisibleJigsaw("InvRapidityJigsaw", kSetRapidity);
  InvisibleJigsaw* InvCBJigsaw = m_RF_helper.addInvisibleJigsaw("InvCBJigsaw", kContraBoost);

  INV->AddFrame(*Ia);
  INV->AddFrame(*Ib);

  INV->AddJigsaw(*InvMassJigsaw);

  INV->AddJigsaw(*InvRapidityJigsaw);
  InvRapidityJigsaw->AddVisibleFrames(LAB->GetListVisibleFrames());

  INV->AddJigsaw(*InvCBJigsaw);
  InvCBJigsaw->AddVisibleFrames(Pa->GetListVisibleFrames(), 0);
  InvCBJigsaw->AddVisibleFrames(Pb->GetListVisibleFrames(), 1);
  InvCBJigsaw->AddInvisibleFrame(*Ia, 0);
  InvCBJigsaw->AddInvisibleFrame(*Ib, 1);

  CombinatoricGroup* VIS = m_RF_helper.addCombinatoricGroup("VIS");
  MinMassesCombJigsaw* CombPPJigsaw = m_RF_helper.addCombinatoricJigsaw("CombPPJigsaw", kMinMasses);
  MinMassesCombJigsaw* CombPaJigsaw = m_RF_helper.addCombinatoricJigsaw("CombPaJigsaw", kMinMasses);
  MinMassesCombJigsaw* CombPbJigsaw = m_RF_helper.addCombinatoricJigsaw("CombPbJigsaw", kMinMasses);

  VIS->AddFrame(*V1a);
  VIS->SetNElementsForFrame(*V1a,0,false);
  VIS->AddFrame(*V1b);
  VIS->SetNElementsForFrame(*V1b,0,false);
  VIS->AddFrame(*V2a);
  VIS->SetNElementsForFrame(*V2a,1,false);
  VIS->AddFrame(*V2b);
  VIS->SetNElementsForFrame(*V2b,1,false);

  VIS->AddJigsaw(*CombPPJigsaw);
  CombPPJigsaw->AddFrame(*V1a,0);
  CombPPJigsaw->AddFrame(*V1b,1);
  CombPPJigsaw->AddFrame(*V2a,0);
  CombPPJigsaw->AddFrame(*V2b,1);

  VIS->AddJigsaw(*CombPaJigsaw);
  CombPaJigsaw->AddFrame(*V1a,0);
  CombPaJigsaw->AddFrame(*V2a,1);

  VIS->AddJigsaw(*CombPbJigsaw);
  CombPbJigsaw->AddFrame(*V1b,0);
  CombPbJigsaw->AddFrame(*V2b,1);

  LAB->InitializeAnalysis();

  //
  // self-assembling tree for QCD rejection variables
  //
  LabRecoFrame* LAB_bkg = m_RF_helper.addLabFrame("LAB_bkg");
  SelfAssemblingRecoFrame* S_bkg = m_RF_helper.addSAFrame("S_bkg");
  VisibleRecoFrame* V_bkg = m_RF_helper.addVisibleFrame("V_bkg");
  InvisibleRecoFrame* I_bkg = m_RF_helper.addInvisibleFrame("I_bkg");

  LAB_bkg->SetChildFrame(*S_bkg);
  S_bkg->AddChildFrame(*V_bkg);
  S_bkg->AddChildFrame(*I_bkg);

  LAB_bkg->InitializeTree();

  InvisibleGroup* INV_bkg = m_RF_helper.addInvisibleGroup("INV_bkg");
  InvisibleJigsaw* InvMass_bkg = m_RF_helper.addInvisibleJigsaw("InvMass_bkg", kSetMass);
  InvisibleJigsaw* InvRapidity_bkg = m_RF_helper.addInvisibleJigsaw("InvRapidity_bkg", kSetRapidity);

  INV_bkg->AddFrame(*I_bkg);

  INV_bkg->AddJigsaw(*InvMass_bkg);

  INV_bkg->AddJigsaw(*InvRapidity_bkg);
  InvRapidity_bkg->AddVisibleFrames(LAB_bkg->GetListVisibleFrames());

  CombinatoricGroup* VIS_bkg = m_RF_helper.addCombinatoricGroup("VIS_bkg");

  VIS_bkg->AddFrame(*V_bkg);
  VIS_bkg->SetNElementsForFrame(*V_bkg,1,false);

  LAB_bkg->InitializeAnalysis();

  //
  // ISR-assisted analysis tree for compressed SR's
  //
  LabRecoFrame* LAB_ISR = m_RF_helper.addLabFrame("LAB_ISR");
  DecayRecoFrame* CM_ISR = m_RF_helper.addDecayFrame("CM_ISR");
  DecayRecoFrame* S_ISR = m_RF_helper.addDecayFrame("S_ISR");
  VisibleRecoFrame* ISR_ISR = m_RF_helper.addVisibleFrame("ISR_ISR");
  VisibleRecoFrame* V_ISR = m_RF_helper.addVisibleFrame("V_ISR");
  InvisibleRecoFrame* I_ISR = m_RF_helper.addInvisibleFrame("I_ISR");

  LAB_ISR->SetChildFrame(*CM_ISR);
  CM_ISR->AddChildFrame(*ISR_ISR);
  CM_ISR->AddChildFrame(*S_ISR);
  S_ISR->AddChildFrame(*V_ISR);
  S_ISR->AddChildFrame(*I_ISR);

  LAB_ISR->InitializeTree();

  InvisibleGroup* INV_ISR = m_RF_helper.addInvisibleGroup("INV_ISR");
  InvisibleJigsaw* InvMass_ISR = m_RF_helper.addInvisibleJigsaw("InvMass_ISR", kSetMass);

  INV_ISR->AddFrame(*I_ISR);

  INV_ISR->AddJigsaw(*InvMass_ISR);

  CombinatoricGroup* VIS_ISR = m_RF_helper.addCombinatoricGroup("VIS_ISR");
  MinMassesCombJigsaw* SplitVis_ISR = m_RF_helper.addCombinatoricJigsaw("SplitVis_ISR", kMinMasses);

  VIS_ISR->AddFrame(*ISR_ISR);
  VIS_ISR->SetNElementsForFrame(*ISR_ISR,1,false);
  VIS_ISR->AddFrame(*V_ISR);
  VIS_ISR->SetNElementsForFrame(*V_ISR,0,false);

  VIS_ISR->AddJigsaw(*SplitVis_ISR);
  // "0" group (ISR)
  SplitVis_ISR->AddFrame(*ISR_ISR, 0);
  // "1" group (V + I)
  SplitVis_ISR->AddFrame(*V_ISR,1);
  SplitVis_ISR->AddFrame(*I_ISR,1);

  LAB_ISR->InitializeAnalysis();
}

void ZeroLeptonRJigsaw2016::ProcessEvent(AnalysisEvent *event)
{
  auto electrons  = event->getElectrons(7, 2.47, ELooseLH);
  auto muons      = event->getMuons(7, 2.7, MuMedium);
  auto jets       = event->getJets(50., 2.8);
  auto metVec     = event->getMET();
  double met      = metVec.Et();

  // Reject events with bad jets
  if (countObjects(jets, 50, 2.8, NOT(LooseBadJet))!=0) return;
  //if (jets.size()>0 && jets[0].Pt()>100 && jets[0].pass(NOT(TightBadJet))) return;
  //if (jets.size()>1 && jets[1].Pt()>100 && jets[1].pass(NOT(TightBadJet))) return;

  // Standard SUSY overlap removal
  jets       = overlapRemoval(jets, electrons, 0.2);
  electrons  = overlapRemoval(electrons, jets, 0.4);
  muons      = overlapRemoval(muons, jets, 0.4);

  // not doing overlap removal between electrons and muons with identical track
  // FIXME: not doing overlap removal between electrons within DeltaR=0.5

  // auto signalMuons     = filterObjects(muons, 25, 2.7, MuD0Sigma3|MuZ05mm|MuIsoGradientLoose);
  // auto signalElectrons = filterObjects(electrons, 25, 2.47, ETightLH|ED0Sigma5|EZ05mm|EIsoGradientLoose);
  // auto goodJets        = filterObjects(jets, 50);
  // auto bjets           = filterObjects(goodJets, 50, 2.5, BTag77MV2c20);

  // baseline lepton veto
  auto leptons=electrons + muons;
  if (leptons.size() > 0) {
    return;
  }

  // preselection
  float meffIncl    = sumObjectsPt(jets) + met;
  if(met < 250) return;
  if(jets.size() < 2) return;
  if(jets[0].Pt() < 200.) return;
  if(jets[1].Pt() < 50.) return;
  if(meffIncl < 800) return;

  int Njets = jets.size();

  // RestFrames analysis for different RJR SR's and event trees:

  // nominal tree
  LabRecoFrame* LAB = m_RF_helper.getLabFrame("LAB");
  InvisibleGroup* INV = m_RF_helper.getInvisibleGroup("INV");
  CombinatoricGroup* VIS = m_RF_helper.getCombinatoricGroup("VIS");

  LAB->ClearEvent();
  std::vector<RFKey> jetID;
  for(int i = 0; i < Njets; i++){
    jetID.push_back(VIS->AddLabFrameFourVector(jets[i]));
  }
  INV->SetLabFrameThreeVector(metVec.Vect());
  LAB->AnalyzeEvent();

  // QCD rejection tree
  LabRecoFrame* LAB_bkg = m_RF_helper.getLabFrame("LAB_bkg");
  InvisibleGroup* INV_bkg = m_RF_helper.getInvisibleGroup("INV_bkg");
  CombinatoricGroup* VIS_bkg = m_RF_helper.getCombinatoricGroup("VIS_bkg");

  LAB_bkg->ClearEvent();
  std::vector<RFKey> jetID_bkg;
  for(int i = 0; i < Njets; i++){
    TLorentzVector jet = jets[i];
    jet.SetPtEtaPhiM(jets[i].Pt(),0.0,jets[i].Phi(),jets[i].M());
    jetID_bkg.push_back(VIS_bkg->AddLabFrameFourVector(jet));
  }
  INV_bkg->SetLabFrameThreeVector(metVec.Vect());
  LAB_bkg->AnalyzeEvent();

  // ISR-assisted tree
  LabRecoFrame* LAB_ISR = m_RF_helper.getLabFrame("LAB_ISR");
  InvisibleGroup* INV_ISR = m_RF_helper.getInvisibleGroup("INV_ISR");
  CombinatoricGroup* VIS_ISR = m_RF_helper.getCombinatoricGroup("VIS_ISR");

  LAB_ISR->ClearEvent();
  std::vector<RFKey> jetID_ISR;
  for(int i = 0; i < Njets; i++){
    TLorentzVector jet = jets[i];
    jet.SetPtEtaPhiM(jets[i].Pt(),0.0,jets[i].Phi(),jets[i].M());
    jetID_ISR.push_back(VIS_ISR->AddLabFrameFourVector(jet));
  }
  INV_ISR->SetLabFrameThreeVector(metVec.Vect());
  LAB_ISR->AnalyzeEvent();


  // Event analysis from reconstructed trees

  // QCD rejection variables
  InvisibleRecoFrame* I_bkg = m_RF_helper.getInvisibleFrame("I_bkg");

  TLorentzVector Psib = I_bkg->GetSiblingFrame().GetFourVector(*LAB_bkg);
  TLorentzVector Pmet = I_bkg->GetFourVector(*LAB_bkg);

  double Rsib = std::max(0.,Psib.Vect().Dot(Pmet.Vect().Unit()));
  Rsib = Rsib / (Pmet.Pt() + Rsib);

  TVector3 boostQCD = (Pmet+Psib).BoostVector();
  Psib.Boost(-boostQCD);
  double cosQCD = -1.*Psib.Vect().Unit().Dot(boostQCD.Unit());
  cosQCD = (1.-cosQCD)/2.;

  // deltaQCD definition
  double deltaQCD = (cosQCD-Rsib)/(cosQCD+Rsib);

  // nominal tree variables
  DecayRecoFrame* PP = m_RF_helper.getDecayFrame("PP");
  VisibleRecoFrame* V1a = m_RF_helper.getVisibleFrame("V1a");
  VisibleRecoFrame* V1b = m_RF_helper.getVisibleFrame("V1b");
  VisibleRecoFrame* V2a = m_RF_helper.getVisibleFrame("V2a");
  VisibleRecoFrame* V2b = m_RF_helper.getVisibleFrame("V2b");
  InvisibleRecoFrame* Ia = m_RF_helper.getInvisibleFrame("Ia");
  InvisibleRecoFrame* Ib = m_RF_helper.getInvisibleFrame("Ib");

  int NJ1a = VIS->GetNElementsInFrame(*V1a);
  int NJ1b = VIS->GetNElementsInFrame(*V1b);
  int NJ2a = VIS->GetNElementsInFrame(*V2a);
  int NJ2b = VIS->GetNElementsInFrame(*V2b);
  int NJa = NJ1a+NJ2a;
  int NJb = NJ1b+NJ2b;

  TLorentzVector vP_V1aPP = V1a->GetFourVector(*PP);
  TLorentzVector vP_V2aPP = V2a->GetFourVector(*PP);
  TLorentzVector vP_V1bPP = V1b->GetFourVector(*PP);
  TLorentzVector vP_V2bPP = V2b->GetFourVector(*PP);
  TLorentzVector vP_IaPP  = Ia->GetFourVector(*PP);
  TLorentzVector vP_IbPP  = Ib->GetFourVector(*PP);

  // 3D H variables
  double H2PP = (vP_V1aPP + vP_V2aPP + vP_V1bPP + vP_V2bPP).P() + (vP_IaPP+vP_IbPP).P();
  double H3PP = (vP_V1aPP + vP_V2aPP).P() + (vP_V1bPP + vP_V2bPP).P() + (vP_IaPP + vP_IbPP).P();
  double H5PP = vP_V1aPP.P() + vP_V2aPP.P() + vP_V1bPP.P() + vP_V2bPP.P() + (vP_IaPP + vP_IbPP).P();

  double eta12 = std::max(fabs(jets[0].Eta()),fabs(jets[1].Eta()));

  double pTPP_Va  = PP->GetTransverseMomentum(V1a->GetFourVector()+V2a->GetFourVector());
  double pTPP_V1a = V1a->GetTransverseMomentum(*PP);
  double pTPP_V2a = V2a->GetTransverseMomentum(*PP);
  double pTPP_Vb  = PP->GetTransverseMomentum(V1b->GetFourVector()+V2b->GetFourVector());
  double pTPP_V1b = V1b->GetTransverseMomentum(*PP);
  double pTPP_V2b = V2b->GetTransverseMomentum(*PP);
  double pTPP_Ia = Ia->GetTransverseMomentum(*PP);
  double pTPP_Ib = Ib->GetTransverseMomentum(*PP);

  double pPP_Va  = (V1a->GetFourVector(*PP)+V2a->GetFourVector(*PP)).P();
  double pPP_V1a = V1a->GetMomentum(*PP);
  double pPP_V2a = V2a->GetMomentum(*PP);
  double pPP_Vb  = (V1b->GetFourVector(*PP)+V2b->GetFourVector(*PP)).P();
  double pPP_V1b = V1b->GetMomentum(*PP);
  double pPP_V2b = V2b->GetMomentum(*PP);
  // m_pPP_Ia = Ia->GetMomentum(*PP);
  // m_pPP_Ib = Ib->GetMomentum(*PP);

  double pTPP_jet1a = 0.;
  double pTPP_jet2a = 0.;
  double pTPP_jet1b = 0.;
  double pTPP_jet2b = 0.;
  double eta_jet1a = -5.;
  double eta_jet2a = -5.;
  double eta_jet1b = -5.;
  double eta_jet2b = -5.;

  for(int j = 0; j < Njets; j++){
    RestFrame const& frame = VIS->GetFrame(jetID[j]);
    double pT = PP->GetTransverseMomentum(VIS->GetLabFrameFourVector(jetID[j]));
    //double p  = PP->GetFourVector(VIS->GetLabFrameFourVector(jetID[j])).P();
    double eta = VIS->GetLabFrameFourVector(jetID[j]).Eta();

    if(V1a->IsSame(frame) || V2a->IsSame(frame)){
      if(pT > pTPP_jet1a){
	pTPP_jet2a = pTPP_jet1a;
	eta_jet2a  = eta_jet1a;
	pTPP_jet1a = pT;
	eta_jet1a  = eta;
      } else {
	if(pT > pTPP_jet2a){
	  pTPP_jet2a = pT;
	  eta_jet2a  = eta;
	}
      }
    }
    if(V1b->IsSame(frame) || V2b->IsSame(frame)){
      if(pT > pTPP_jet1b){
	pTPP_jet2b = pTPP_jet1b;
	eta_jet2b  = eta_jet1b;
	pTPP_jet1b = pT;
	eta_jet1b  = eta;
      } else {
	if(pT > pTPP_jet2b){
	  pTPP_jet2b = pT;
	  eta_jet2b  = eta;
	}
      }
    }
  }

  double eta12ab = std::max(std::max(fabs(eta_jet1a),fabs(eta_jet2a)),
			    std::max(fabs(eta_jet1b),fabs(eta_jet2b)));

  double pTPP_jet2 = 0;

  if(pTPP_jet1a > pTPP_jet1b){
    pTPP_jet2 = std::max(pTPP_jet1b,pTPP_jet2a);
  } else {
    pTPP_jet2 = std::max(pTPP_jet1a,pTPP_jet2b);
  }

  double HT3PP = pTPP_Va + pTPP_Vb + H2PP/2.;
  double HT5PP = pTPP_V1a + pTPP_V1b +
    pTPP_V2a + pTPP_V2b + H2PP/2.;

  /// squark ratios
  double R_H2PP_H3PP = H2PP/H3PP;
  double R_pTj2_HT3PP = pTPP_jet2 / HT3PP;

  /// gluino ratios
  double R_HT5PP_H5PP = HT5PP/H5PP;
  double R_H2PP_H5PP = H2PP/H5PP;
  double minR_pTj2i_HT3PPi = std::min(pTPP_jet2a/(pTPP_V1a+pTPP_V2a+pTPP_Ia),
				      pTPP_jet2b/(pTPP_V1b+pTPP_V2b+pTPP_Ib));
  double maxR_H1PPi_H2PPi = std::max(pPP_Va/(pPP_V1a+pPP_V2a),pPP_Vb/(pPP_V1b+pPP_V2b));

  TVector3 vP_PP = PP->GetFourVector(*LAB).Vect();
  double pTCM = vP_PP.Pt();
  double pZCM = fabs(vP_PP.Pz());

  //double RPZ_HT3PP = pZCM / (pZCM + HT3PP);
  double RPZ_HT5PP = pZCM / (pZCM + HT5PP);
  double RPT_HT3PP = pTCM / (pTCM + HT3PP);
  double RPT_HT5PP = pTCM / (pTCM + HT5PP);

  // squark SR's
  // SR-S1
  if(R_H2PP_H3PP > 0.55 && R_H2PP_H3PP < 0.9 && R_pTj2_HT3PP > 0.16 && eta12 < 0.8 && deltaQCD > 0.1 && RPT_HT3PP < 0.08){
    if(HT3PP > 1000.)
      accept("SR_RJR_S1a");
    if(HT3PP > 1200.)
      accept("SR_RJR_S1b");
  }
  // SR-S2

  if(R_H2PP_H3PP > 0.5 && R_H2PP_H3PP < 0.95 && R_pTj2_HT3PP > 0.14 && eta12 < 1.1 && deltaQCD > 0.05 && RPT_HT3PP < 0.08){
    if(HT3PP > 1400.)
      accept("SR_RJR_S2a");
    if(HT3PP > 1600.)
      accept("SR_RJR_S2b");
  }
  // SR-S3
  if(R_H2PP_H3PP > 0.45 && R_H2PP_H3PP < 0.98 && R_pTj2_HT3PP > 0.13 && eta12 < 1.4 && deltaQCD > 0.025 && RPT_HT3PP < 0.08){
    if(HT3PP > 1800.)
      accept("SR_RJR_S3a");
    if(HT3PP > 2100.)
      accept("SR_RJR_S3b");
  }
  // SR-S4
  if(R_pTj2_HT3PP > 0.13 && deltaQCD > 0. && RPT_HT3PP < 0.08){
    if(HT3PP > 2400.)
      accept("SR_RJR_S4");
  }

  // gluino SR's
  if(NJa > 1 && NJb > 1){
    // SR-G1
    if(R_H2PP_H5PP > 0.45 && R_HT5PP_H5PP > 0.7 && maxR_H1PPi_H2PPi < 0.96 && minR_pTj2i_HT3PPi > 0.12 && eta12ab < 1.4 && deltaQCD > 0.05 && RPZ_HT5PP < 0.5 && RPT_HT5PP < 0.08){
      if(HT5PP > 1200.)
	accept("SR_RJR_G1a");
      if(HT5PP > 1400.)
	accept("SR_RJR_G1b");
    }
    // SR-G2
    if(R_H2PP_H5PP > 0.3 && R_HT5PP_H5PP > 0.7 && maxR_H1PPi_H2PPi < 0.97 && minR_pTj2i_HT3PPi > 0.1 && eta12ab < 2.0 && deltaQCD > 0.025 && RPZ_HT5PP < 0.55 && RPT_HT5PP < 0.08){
      if(HT5PP > 1600.)
	accept("SR_RJR_G2a");
      if(HT5PP > 2000.)
	accept("SR_RJR_G2b");
    }
    // SR-G3
    if(R_H2PP_H5PP > 0.2 && R_HT5PP_H5PP > 0.65 && maxR_H1PPi_H2PPi < 0.98 && minR_pTj2i_HT3PPi > 0.08 && deltaQCD > 0.0 && RPZ_HT5PP < 0.6 && RPT_HT5PP < 0.08){
      if(HT5PP > 2400.)
	accept("SR_RJR_G3a");
      if(HT5PP > 2800.)
	accept("SR_RJR_G3b");
    }
    // SR-G4
    if(R_HT5PP_H5PP > 0.65 && maxR_H1PPi_H2PPi < 0.98 && minR_pTj2i_HT3PPi > 0.07 && deltaQCD > 0.0 && RPZ_HT5PP < 0.65 && RPT_HT5PP < 0.08){
      if(HT5PP > 3000.)
	accept("SR_RJR_G4");
    }
  }

  // ISR-assisted analysis

  DecayRecoFrame* CM_ISR = m_RF_helper.getDecayFrame("CM_ISR");
  DecayRecoFrame* S_ISR = m_RF_helper.getDecayFrame("S_ISR");
  VisibleRecoFrame* V_ISR = m_RF_helper.getVisibleFrame("V_ISR");
  VisibleRecoFrame* ISR_ISR = m_RF_helper.getVisibleFrame("ISR_ISR");
  InvisibleRecoFrame* I_ISR = m_RF_helper.getInvisibleFrame("I_ISR");

  int NV = VIS_ISR->GetNElementsInFrame(*V_ISR);
  //int NISR = VIS_ISR->GetNElementsInFrame(*ISR_ISR);

  if(NV < 1) return;

  TVector3 vP_ISR = ISR_ISR->GetFourVector(*CM_ISR).Vect();
  TVector3 vP_I   = I_ISR->GetFourVector(*CM_ISR).Vect();

  double PTISR    = vP_ISR.Mag();
  double MS       = S_ISR->GetMass();
  double RISR     = fabs(vP_I.Dot(vP_ISR.Unit())) / PTISR;
  double dphiISRI = fabs(vP_ISR.DeltaPhi(vP_I));

  int iN = 0;
  double etajV1 = -5.;
  double etajV2 = -5.;
  double etajV3 = -5.;
  for(int i = 0; i < Njets; i++){
    if(VIS_ISR->GetFrame(jetID_ISR[i]) == *V_ISR){
      if(iN == 0){
	etajV1 = jets[i].Eta();
      }
      if(iN == 1){
	etajV2 = jets[i].Eta();
      }
      if(iN == 2){
	etajV3 = jets[i].Eta();
      }
      iN++;
    }
  }

  double etaV1 = fabs(etajV1);
  double etaV2 = std::max(fabs(etajV2), etaV1);
  double etaV3 = std::max(fabs(etajV3), etaV2);

  // compressed SR's
  // SR-C1
  if(RISR > 0.95 && PTISR > 1000. && dphiISRI > acos(-1.)*0.95 && MS >= 0. && NV > 0 && etaV1 <= 2.8){
    accept("SR_RJR_C1");
  }
  // SR-C2
  if(RISR > 0.9 && PTISR > 1000. && dphiISRI > acos(-1.)*0.97 && MS >= 100. && NV > 0 && etaV1 <= 1.2){
    accept("SR_RJR_C2");
  }
  // SR-C3
  if(RISR > 0.8 && PTISR > 800. && dphiISRI > acos(-1.)*0.98 && MS >= 200. && NV > 1 && etaV2 <= 1.4){
    accept("SR_RJR_C3");
  }
  // SR-C4
  if(RISR > 0.7 && PTISR > 700. && dphiISRI > acos(-1.)*0.95 && MS >= 450. && NV > 1 && etaV2 <= 1.4){
    accept("SR_RJR_C4");
  }
  // SR-C5
  if(RISR > 0.7 && PTISR > 700. && dphiISRI > acos(-1.)*0.95 && MS >= 450. && NV > 2 && etaV3 <= 1.4){
    accept("SR_RJR_C5");
  }

  return;
}
#endif
