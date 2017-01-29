#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(ZeroLepton2012)

void ZeroLepton2012::Init()
{
  addRegions({"SR2jl","SR2jm","SR2jt","SR2jW","SR3jt",
		       "SR4jvl","SR4jl","SR4jm","SR4jt","SR4jW",
		       "SR5jt","SR6jl","SR6jm","SR6jt","SR6jvt"});
}

//Recursively go through jets to find close-by pairs with mass consistent with W
static int FindResolvedW(AnalysisObjects &nonWjets) {
  if (nonWjets.size()<2) return 0;
  float minDR=1000.;
  unsigned int imin=0;
  unsigned int jmin=0;
  int foundPairs=0;
  for(unsigned int ii=0; ii<nonWjets.size()-1;ii++) {
    for(unsigned int jj=ii+1; jj<nonWjets.size();jj++) {
      float dR=nonWjets[ii].DeltaR(nonWjets[jj]);
      if (dR<minDR) {
	minDR=dR;
	imin=ii;
	jmin=jj;
      }
    }
  }
  auto sumJets=nonWjets[imin]+nonWjets[jmin];
  if (sumJets.M()>60 && sumJets.M()<100) {
    foundPairs++;
    nonWjets.erase(nonWjets.begin()+jmin);
  }
  nonWjets.erase(nonWjets.begin()+imin);
  foundPairs+=FindResolvedW(nonWjets);
  return foundPairs;
}

void ZeroLepton2012::ProcessEvent(AnalysisEvent *event)
{
    
    auto electrons  = event->getElectrons(10,2.47);
    auto muons      = event->getMuons(10,2.5);
    auto jets       = event->getJets(20.,2.8);
    auto metVec     = event->getMET();
    double met      = metVec.Et();
    
    // Standard SUSY overlap removal
    jets       = overlapRemoval(jets,electrons,0.2);
    electrons  = overlapRemoval(electrons,jets,0.4);
    muons      = overlapRemoval(muons,jets,0.4);

    // Count electrons and muons
    int nb_baseline_muons     = muons.size();
    int nb_baseline_electrons = electrons.size();

    // count jets
    int nb_jets_40  = countObjects(jets,40.);
    int nb_jets_60  = countObjects(jets,60.);

    // Preselection
    if (nb_baseline_muons+nb_baseline_electrons != 0) return;
    if (met<160) return;
    if (nb_jets_60<2) return;
    if (jets[0].Pt()<130.) return;

    // Define Meff, delta phi 
    float meff[7];
    for(int nJets=2; nJets <=6; nJets++) meff[nJets]=sumObjectsPt(jets,nJets,60)+met;
    float meffIncl=sumObjectsPt(jets,999,60)+met;
    float dphi_min3 = minDphi(metVec,jets,3,40);
    float dphi_minRest = minDphi(metVec,jets,999,40);

    //look for (un)resolved W jets
    AnalysisObjects nonWjets;
    int numUnresolvedW=0;
    for (const auto& jet : jets) {
      if (jet.Pt()>40.) {
	if (jet.M()>60 && jet.M()<100) numUnresolvedW++;
	else nonWjets.push_back(jet);
      }
    }
    int numResolvedW=FindResolvedW(nonWjets);
    float met_sig = met/sqrt(meffIncl-met);

    if (dphi_min3<0.4) return;

    if (nb_jets_60 >= 2) { //2-jet SRs
      if (meffIncl>800 && met_sig>=8) accept("SR2jl");
      if (meffIncl>1200 && met_sig>=15) accept("SR2jm");
      if (meffIncl>1600 && met_sig>=15) accept("SR2jt");
      if (meffIncl>1800 && met/meff[2]>0.25 && numUnresolvedW>=2 ) accept("SR2jW");
    }
    if (nb_jets_60 >= 3) { //3-jet SRs
      if (meffIncl>2200 && met/meff[3]>0.3) accept("SR3jt");
    }
    if (nb_jets_60 >= 4 && dphi_minRest>0.2) { //4-jet SRs
      if (meffIncl>700 && met_sig>=10) accept("SR4jvl");
      if (meffIncl>1000 && met_sig>=10) accept("SR4jl");
      if (meffIncl>1300 && met/meff[4]>0.4) accept("SR4jm");
      if (meffIncl>2200 && met/meff[4]>0.25) accept("SR4jt");
    }
    if (nb_jets_60 >= 2 && nb_jets_40 >=4 && dphi_minRest>0.2) { //4-jet SR with W
      float meff4=meff[2]+jets[2].Pt()+jets[3].Pt();
      if (meffIncl>1100 && met/meff4>0.35 && numUnresolvedW>=1 && numResolvedW>=1) {
	accept("SR4jW");
      }
    }
    if (nb_jets_60 >= 5 && dphi_minRest>0.2) { //5-jet SRs
      if (meffIncl>1200 && met/meff[5]>0.2) accept("SR5jt");
    }
    if (nb_jets_60 >= 6 && dphi_minRest>0.2) { //6-jet SRs
      if (meffIncl>900 && met/meff[6]>0.2) accept("SR6jl");
      if (meffIncl>1200 && met/meff[6]>0.2) accept("SR6jm");
      if (meffIncl>1500 && met/meff[6]>0.25) accept("SR6jt");
      if (meffIncl>1700 && met/meff[6]>0.15) accept("SR6jvt");
    }
    return;
}
