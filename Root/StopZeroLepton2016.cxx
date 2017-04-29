#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(StopZeroLepton2016)

void StopZeroLepton2016::Init()
{
  addRegions({"SRA_TT","SRA_TW","SRA_T0"});
  addRegions({"SRB_TT","SRB_TW","SRB_T0"});
  addRegions({"SRC1","SRC2","SRC3","SRC4","SRC5"});
  addRegions({"SRD_low","SRD_high"});
  addRegions({"SRE"});

  LabRecoFrame*       LAB = m_RF_helper.addLabFrame("LAB");
  DecayRecoFrame*     CM  = m_RF_helper.addDecayFrame("CM");
  DecayRecoFrame*     S   = m_RF_helper.addDecayFrame("S");
  VisibleRecoFrame*   ISR = m_RF_helper.addVisibleFrame("ISR");
  VisibleRecoFrame*   V   = m_RF_helper.addVisibleFrame("V");
  InvisibleRecoFrame* I   = m_RF_helper.addInvisibleFrame("I");

  LAB->SetChildFrame(*CM);
  CM->AddChildFrame(*ISR);
  CM->AddChildFrame(*S);
  S->AddChildFrame(*V);
  S->AddChildFrame(*I);

  LAB->InitializeTree(); 

  InvisibleGroup*  INV = m_RF_helper.addInvisibleGroup("INV");
  INV->AddFrame(*I);
  CombinatoricGroup* VIS = m_RF_helper.addCombinatoricGroup("VIS");
  VIS->AddFrame(*ISR);
  VIS->SetNElementsForFrame(*ISR,1,false);
  VIS->AddFrame(*V);
  VIS->SetNElementsForFrame(*V,0,false);

  // set the invisible system mass to zero
  InvisibleJigsaw*   InvMass = m_RF_helper.addInvisibleJigsaw("InvMass",kSetMass);
  INV->AddJigsaw(*InvMass);

  // define the rule for partitioning objects between "ISR" and "V"
  MinMassesCombJigsaw* SplitVis = m_RF_helper.addCombinatoricJigsaw("CombPPJigsaw", kMinMasses);
  VIS->AddJigsaw(*SplitVis);
  // "0" group (ISR)
  SplitVis->AddFrame(*ISR, 0);
  // "1" group (V + I)
  SplitVis->AddFrame(*V,1);
  SplitVis->AddFrame(*I,1);

  LAB->InitializeAnalysis();
}

void StopZeroLepton2016::ProcessEvent(AnalysisEvent *event)
{
  auto baseElectrons  = event->getElectrons(7, 2.47, EVeryLooseLH);
  auto baseMuons      = event->getMuons(6,2.7, MuLoose);
  auto baseJets       = event->getJets(20., 2.8); 
  //auto baseTaus       = event->getTaus(20., 2.5, TauLoose); 
  auto metVec         = event->getMET();
  float Met = metVec.Pt();
  
  // Reject events with bad jets
  // CHECK: right place?
  if (countObjects(baseJets, 20, 4.5, NOT(LooseBadJet))!=0) return;

  // overlap removal
  auto radiusCalcLepton = [] (const AnalysisObject& lepton, const AnalysisObject& ) { return std::min(0.4, 0.04 + 10/lepton.Pt()); };
  auto muJetSpecial = [] (const AnalysisObject& jet, const AnalysisObject& muon) { 
    if (jet.pass(NOT(BTag77MV2c20)) && (jet.pass(LessThan3Tracks) || muon.Pt()/jet.Pt()>0.7)) return 0.2;
    else return 0.;
  };
  baseMuons     = overlapRemoval(baseMuons, baseElectrons, 0.01, NOT(MuCaloTaggedOnly));
  baseElectrons = overlapRemoval(baseElectrons, baseMuons, 0.01);

  baseJets      = overlapRemoval(baseJets, baseElectrons, 0.2, NOT(BTag77MV2c20));
  baseElectrons = overlapRemoval(baseElectrons, baseJets, 0.2);

  baseJets      = overlapRemoval(baseJets, baseMuons, muJetSpecial, NOT(BTag77MV2c20));
  baseMuons     = overlapRemoval(baseMuons,baseJets, 0.2);

  baseMuons     = overlapRemoval(baseMuons,baseJets, radiusCalcLepton);
  baseElectrons = overlapRemoval(baseElectrons,baseJets, radiusCalcLepton);
  
  //baseTaus      = overlapRemoval(baseTaus, baseElectrons, 0.1);

  auto baseLeptons     = baseElectrons+baseMuons;
  auto signalJets      = filterObjects(baseJets, 20, 2.5, JVT50Jet);
  auto signalBJets     = filterObjects(signalJets, 20, 2.5, BTag77MV2c20);
  auto nonBJets        = filterObjects(signalJets, 20, 2.5, NOT(BTag77MV2c20));
  auto signalElectrons = filterObjects(baseElectrons, 20, 2.47, ELooseLH|ED0Sigma5|EZ05mm|EIsoLooseTrack);
  auto signalMuons     = filterObjects(baseMuons, 20, 2.7, MuD0Sigma3|MuZ05mm|MuIsoLooseTrack);
  auto signalLeptons   = signalElectrons+signalMuons;
  
  auto fatJetsR8 = reclusterJets(signalJets, 0.8, 150, 0.2, 0.05);  
  auto fatJetsR12 = reclusterJets(signalJets, 1.2, 150, 0.2, 0.05);  

  if (Met<250 || baseLeptons.size()>0 || signalJets.size()<4 || signalBJets.size()<1 || signalJets[3].Pt()<40)  return;

  float DRBB = 0;
  int NBJets = signalBJets.size();

  float AntiKt8M_0 = 0;
  float AntiKt8M_1 = 0;
  float AntiKt12M_0 = 0;
  float AntiKt12M_1 = 0;
  float MtTauCand = -1 ;
  float MtBMin = 0 ; 
  float MtBMax = 0 ; 

  if (fatJetsR8.size()>0)  AntiKt8M_0 = fatJetsR8[0].M() ;
  if (fatJetsR8.size()>1)  AntiKt8M_1 = fatJetsR8[1].M() ;
  if (fatJetsR12.size()>0) AntiKt12M_0 = fatJetsR12[0].M() ; 
  if (fatJetsR12.size()>1) AntiKt12M_1 = fatJetsR12[1].M() ; 
  if (NBJets>1) DRBB = signalBJets[1].DeltaR(signalBJets[0]);
  
  double dPhi_min = 1000.;
  double dPhi_max = 0.;
  if (signalBJets.size()>=2)  {
    for (const auto& jet : signalBJets) {
      double dphi = fabs(metVec.DeltaPhi(jet));
      if (dphi<dPhi_min) {
	dPhi_min=dphi;
	MtBMin=calcMT(jet,metVec);
      }
      if (dphi>dPhi_max) {
	dPhi_max=dphi;
	MtBMax=calcMT(jet,metVec);
      }
    }
  }
  float realWMass = 80.385;
  float realTopMass = 173.210;
    
  //Chi2 method
  double Chi2min = DBL_MAX;
  int W1j1_low = -1,W1j2_low = -1,W2j1_low = -1,W2j2_low = -1,b1_low = -1,b2_low = -1;
    
  double m_mt2Chi2 = 0;
    
  if (signalJets.size()>=4 && signalBJets.size()>=2 && nonBJets.size()>=2)
    {
      for(int W1j1=0; W1j1<(int)nonBJets.size(); W1j1++) {// <------------------This lines has to be replaced
	for(int W2j1=0;W2j1<(int)nonBJets.size(); W2j1++) {// <------------------This lines has to be replaced
	  if (W2j1==W1j1) continue;// <------------------This lines has to be added
	  //		for(int W1j1=0; W1j1<(int)ljets.size()-1; W1j1++) {
	  //		for(int W2j1=W1j1+1;W2j1<(int)ljets.size(); W2j1++) {
	  for(int b1=0;b1<(int)signalBJets.size();b1++){
	    for(int b2=0;b2<(int)signalBJets.size();b2++){
	      if(b2==b1) continue;	
	      double chi21, chi22, mW1, mW2, mt1, mt2;
	      
	      if(W2j1>W1j1){
		
		mW1 = nonBJets[W1j1].M();
		mW2 = nonBJets[W2j1].M();
		mt1 = (nonBJets[W1j1]+signalBJets[b1]).M();
	        mt2 = (nonBJets[W2j1]+signalBJets[b2]).M();
		
		chi21 = (mW1-realWMass)*(mW1-realWMass)/realWMass + (mt1-realTopMass)*(mt1-realTopMass)/realTopMass;
		chi22 = (mW2-realWMass)*(mW2-realWMass)/realWMass + (mt2-realTopMass)*(mt2-realTopMass)/realTopMass;
		
		if(Chi2min > (chi21 + chi22)){
		  Chi2min = chi21 + chi22;
		  if(chi21 < chi22){
		    W1j1_low = W1j1;
		    W1j2_low = -1;
		    W2j1_low = W2j1;
		    W2j2_low = -1;
		    b1_low = b1;
		    b2_low = b2;
		  }
		  else{
		    W2j1_low = W1j1;
		    W2j2_low = -1;
		    W1j1_low = W2j1;
		    W1j2_low = -1;
		    b2_low = b1;
		    b1_low = b2;
		  }
		}
	      }
	     
	      if (nonBJets.size()<3)
		continue;
		  
	      for(int W1j2=W1j1+1;W1j2 < (int)nonBJets.size(); W1j2++) {
		if(W1j2==W2j1) continue;
		    
		//try bll,bl top candidates
		mW1 = (nonBJets[W1j1] + nonBJets[W1j2]).M();
		mW2 = (nonBJets[W2j1]).M();
		mt1 = (nonBJets[W1j1] + nonBJets[W1j2] + signalBJets[b1]).M();
		mt2 = (nonBJets[W2j1] + signalBJets[b2]).M();
		chi21 = (mW1-realWMass)*(mW1-realWMass)/realWMass + (mt1-realTopMass)*(mt1-realTopMass)/realTopMass;
		chi22 = (mW2-realWMass)*(mW2-realWMass)/realWMass + (mt2-realTopMass)*(mt2-realTopMass)/realTopMass;
		if(Chi2min > (chi21 + chi22)){
		  Chi2min = chi21 + chi22;
		  if(chi21 < chi22){
		    W1j1_low = W1j1;
		    W1j2_low = W1j2;
		    W2j1_low = W2j1;
		    W2j2_low = -1;
		    b1_low = b1;
		    b2_low = b2;
		  }
		  else{
		    W2j1_low = W1j1;
		    W2j2_low = W1j2;
		    W1j1_low = W2j1;
		    W1j2_low = -1;
		    b2_low = b1;
		    b1_low = b2;
		  }
		}
		if(nonBJets.size() < 4)continue;
		//try bll, bll top candidates
		for(int W2j2=W2j1+1;W2j2<(int)nonBJets.size(); W2j2++){
		  if((W2j2==W1j1) || (W2j2==W1j2)) continue;
		  if(W2j1<W1j1) continue;  //runtime reasons, we don't want combinations checked twice <--------------------This line should be added
		  mW1 = (nonBJets[W1j1] + nonBJets[W1j2]).M();
		  mW2 = (nonBJets[W2j1] + nonBJets[W2j2]).M();
		  mt1 = (nonBJets[W1j1] + nonBJets[W1j2] + signalBJets[b1]).M();
		  mt2 = (nonBJets[W2j1] + nonBJets[W2j2] + signalBJets[b2]).M();
		  chi21 = (mW1-realWMass)*(mW1-realWMass)/realWMass + (mt1-realTopMass)*(mt1-realTopMass)/realTopMass;
		  chi22 = (mW2-realWMass)*(mW2-realWMass)/realWMass + (mt2-realTopMass)*(mt2-realTopMass)/realTopMass;
		  if(Chi2min > (chi21 + chi22)){
		    Chi2min = chi21 + chi22;
		    if(chi21 < chi22){
		      W1j1_low = W1j1;
		      W1j2_low = W1j2;
		      W2j1_low = W2j1;
		      W2j2_low = W2j2;
		      b1_low = b1;
		      b2_low = b2;
		    }
		    else{
		      W2j1_low = W1j1;
		      W2j2_low = W1j2;
		      W1j1_low = W2j1;
		      W1j2_low = W2j2;
		      b2_low = b1;
		      b1_low = b2;
		    }
		  }
		}
	      }		  
	    }
	  }
	}
      }
	
      AnalysisObject WCand0=nonBJets[W1j1_low];
      if (W1j2_low != -1) WCand0 +=nonBJets[W1j2_low];
      AnalysisObject topCand0 = WCand0 + signalBJets[b1_low];
	
      AnalysisObject WCand1 = nonBJets[W2j1_low];
      if(W2j2_low != -1) WCand1 += nonBJets[W2j2_low];
      AnalysisObject topCand1 = WCand1 + signalBJets[b2_low];

      AnalysisObject top0Chi2 = topCand0;
      AnalysisObject top1Chi2 = topCand1;
      // if (BJets.size()>=2)
      {
	double Energy0=sqrt(173.210*173.210+top0Chi2.Pt()*top0Chi2.Pt());
	double Energy1=sqrt(173.210*173.210+top1Chi2.Pt()*top1Chi2.Pt());
	AnalysisObject top0(top0Chi2.Px(),top0Chi2.Py(),0,Energy0,0,0,COMBINED,0);
	AnalysisObject top1(top1Chi2.Px(),top1Chi2.Py(),0,Energy1,0,0,COMBINED,0);
	m_mt2Chi2 = calcMT2(top0,top1,metVec);
      }
      
    }


  float MT2Chi2 = m_mt2Chi2;

  // RestFrames stuff

  double CA_PTISR=0;
  double CA_MS=0;
  double CA_NbV=0;
  double CA_NjV=0;
  double CA_RISR=0;
  double CA_dphiISRI=0;
  double CA_pTjV4=0;
  double CA_pTbV1=0;
 
  LabRecoFrame* LAB = m_RF_helper.getLabFrame("LAB");
  InvisibleGroup* INV = m_RF_helper.getInvisibleGroup("INV");
  CombinatoricGroup* VIS = m_RF_helper.getCombinatoricGroup("VIS");

  LAB->ClearEvent();
  std::vector<RFKey> jetID;

  // use transverse view of jet 4-vectors
  for(const auto & jet : signalJets) 
    jetID.push_back(VIS->AddLabFrameFourVector(jet.transFourVect()));

  INV->SetLabFrameThreeVector(metVec.Vect());
  if(!LAB->AnalyzeEvent()) std::cout << "Something went wrong..." << std::endl;

  DecayRecoFrame*     CM  = m_RF_helper.getDecayFrame("CM");
  DecayRecoFrame*     S   = m_RF_helper.getDecayFrame("S");
  VisibleRecoFrame*   ISR = m_RF_helper.getVisibleFrame("ISR");
  VisibleRecoFrame*   V   = m_RF_helper.getVisibleFrame("V");
  InvisibleRecoFrame* I   = m_RF_helper.getInvisibleFrame("I");

  int m_NjV(0);
  int m_NbV(0);
  int m_NbISR(0);
  double m_pTjV4(0.);
  double m_pTbV1(0);
  double m_PTISR(0.);
  double m_MS(0.);
  double m_RISR(0.);
  double m_dphiISRI(0.);

  for(int i = 0; i < int(signalJets.size()); i++){
    if (VIS->GetFrame(jetID[i]) == *V){ // sparticle group
      m_NjV++;
      if (m_NjV == 4) m_pTjV4 = signalJets[i].Pt();
      if (signalJets[i].pass(BTag77MV2c20) && fabs(signalJets[i].Eta())<2.5) {
	m_NbV++;
	if (m_NbV == 1) m_pTbV1 = signalJets[i].Pt();
      }
    } else {
      if (signalJets[i].pass(BTag77MV2c20) && fabs(signalJets[i].Eta())<2.5)
	m_NbISR++;
    }
  }    

  // need at least one jet associated with MET-side of event
  if(m_NjV >= 1)
    {
      TVector3 vP_ISR = ISR->GetFourVector(*CM).Vect();
      TVector3 vP_I   = I->GetFourVector(*CM).Vect();
	
      m_PTISR = vP_ISR.Mag();
      m_RISR = fabs(vP_I.Dot(vP_ISR.Unit())) / m_PTISR;

      m_MS = S->GetMass();
 
      m_dphiISRI = fabs(vP_ISR.DeltaPhi(vP_I));

      CA_PTISR=m_PTISR;
      CA_MS=m_MS;
      CA_NbV=m_NbV;
      CA_NjV=m_NjV;
      CA_RISR=m_RISR;
      CA_dphiISRI=m_dphiISRI;
      CA_pTjV4=m_pTjV4;
      CA_pTbV1=m_pTbV1;
    }

  double Ht = sumObjectsPt(signalJets); 
  double HtSig = Met/sqrt(Ht); 
  double sumTagJetPt = sumObjectsPt(signalBJets,2);

  if ( AntiKt12M_0>120. && AntiKt12M_1 > 120. && Met > 400. && AntiKt8M_0>60. && MtBMin > 200. && MtTauCand < 0 && DRBB > 1. && MT2Chi2>400. && NBJets>=2 )
    accept("SRA_TT"); 
  
  if ( AntiKt12M_0>120. && AntiKt12M_1 < 120. && AntiKt12M_1 > 60. && Met > 500. && AntiKt8M_0>60. && MtBMin > 200. && MtTauCand < 0 && MT2Chi2>400. && NBJets>=2)
    accept("SRA_TW");
  
  if (AntiKt12M_0>120. && AntiKt12M_1 < 60. &&  AntiKt8M_0 > 60. && Met>550. && MtBMin > 200. && MtTauCand < 0 && MT2Chi2>500. && NBJets>=2 )
    accept("SRA_T0");

  if (AntiKt12M_0>120. && AntiKt12M_1 > 120. && DRBB > 1.2 && MtBMax > 200. && MtBMin > 200. && MtTauCand < 0 && NBJets>=2)
    accept("SRB_TT");

  if (AntiKt12M_0>120. && AntiKt12M_1 < 120. && AntiKt12M_1 > 60. && DRBB > 1.2 && MtBMax > 200. && MtBMin > 200. && MtTauCand < 0 && NBJets>=2)
    accept("SRB_TW");

  if (AntiKt12M_0>120. && AntiKt12M_1 < 60. &&  MtBMin > 200. && DRBB > 1.2 && MtBMax > 200. && Met > 250. && MtTauCand < 0 && NBJets>=2)
    accept("SRB_T0");

  if (CA_NbV >= 1 && CA_NjV >= 5 && CA_pTbV1 > 40 && CA_MS > 300 && CA_dphiISRI > 3.00 && CA_PTISR > 400 && CA_pTjV4 > 50 && CA_RISR >= 0.30 && CA_RISR <= 0.4)
    accept("SRC1");

  if (CA_NbV >= 1 && CA_NjV >= 5 && CA_pTbV1 > 40 && CA_MS > 300 && CA_dphiISRI > 3.00 && CA_PTISR > 400 && CA_pTjV4 > 50 && CA_RISR >= 0.4 && CA_RISR <= 0.5)
    accept("SRC2");

  if (CA_NbV >= 1 && CA_NjV >= 5 && CA_pTbV1 > 40 && CA_MS > 300 && CA_dphiISRI > 3.00 && CA_PTISR > 400 && CA_pTjV4 > 50 && CA_RISR >= 0.5 && CA_RISR <= 0.6)
    accept("SRC3");

  if (CA_NbV >= 1 && CA_NjV >= 5 && CA_pTbV1 > 40 && CA_MS > 300 && CA_dphiISRI > 3.00 && CA_PTISR > 400 && CA_pTjV4 > 50 && CA_RISR >= 0.6 && CA_RISR <= 0.7)
    accept("SRC4");

  if (CA_NbV >= 1 && CA_NjV >= 5 && CA_pTbV1 > 40 && CA_MS > 300 && CA_dphiISRI > 3.00 && CA_PTISR > 400 && CA_pTjV4 > 50 && CA_RISR >= 0.7 && CA_RISR <= 0.8)
    accept("SRC5");

  if (signalJets.size()>= 5 && NBJets>=2 && Met > 250 && MtBMin > 250 && MtBMax > 300 && DRBB > 0.8 && signalJets[1].Pt() > 150 && signalJets[3].Pt()>100 && signalJets[4].Pt()>60 && sumTagJetPt>300 && MtTauCand < 0)
    accept("SRD_low");

  if (signalJets.size()>= 5 && signalJets[1].Pt()>150 && signalJets[3].Pt()>80 && signalJets[4].Pt()>60 && MtBMin>350 && MtBMax>450 && NBJets>=2 && Met > 250 && DRBB > 0.8 && sumTagJetPt>400 && MtTauCand < 0)
    accept("SRD_high");

  if (Met > 550 && AntiKt8M_0 > 120 && AntiKt8M_1 > 80 && Ht > 800 && HtSig > 18 && MtBMin > 200 && NBJets>=2)
    accept("SRE");

    return;
}
