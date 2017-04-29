#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(PairedDijets2016)

static int getRoundedWindowLoc(double doubleing_loc){
  double complement = fmod(doubleing_loc,5);
  if (complement > 2.500) return doubleing_loc - complement + 5;
  else {
    return doubleing_loc - complement;
  }
}

static std::pair<int, int> getWindowEdges (int mstop){
  int low_edge = getRoundedWindowLoc(0.89 * mstop + 9);
  int up_edge = getRoundedWindowLoc(1.02 * mstop + 6);
  if (mstop > 2000) up_edge = 99999;
  if (low_edge > up_edge){
    std::cerr << "Invalid window @ "<<mstop<<" GeV: low edge "<<low_edge<<", up edge "<<up_edge<<std::endl;
    low_edge = 0;
    up_edge = 0;
  }
  return std::make_pair(low_edge, up_edge);
}

static bool mass_window(double mass_average, int mass) {
  return (mass_average > getWindowEdges(mass).first) && (mass_average < getWindowEdges(mass).second);
}

static int SRranges[3][3]={ {80, 100, 10},
			    {125, 800, 25},
			    {875, 1750, 125} };

void PairedDijets2016::Init()
{
  addRegions({"Inclusive","Btagged"});
  for(int ii=0; ii<3; ii++) {
    for(int mass = SRranges[ii][0]; mass<=SRranges[ii][1]; mass+=SRranges[ii][2]) {
      addRegion("Inclusive_m"+std::to_string(mass));
      addRegion("Btagged_m"+std::to_string(mass));
    }
  }
}

void PairedDijets2016::ProcessEvent(AnalysisEvent *event)
{
    
  auto candJets   = event->getJets(20., 2.8);
    
  if (countObjects(candJets, 20, 2.8, NOT(TightBadJet))!=0) return;

  candJets = filterObjects(candJets, 20, 2.4);

  auto signalJets      = filterObjects(candJets, 120.);
  int nJets  = signalJets.size();
  if (nJets<4) return;

  //jets pairing
  std::vector<int> vi;
  for (int i=0; i<4; i++)  vi.push_back(i);

  int index1=-1;
  int index2=-1;
  int index3=-1;
  int index4=-1;
  double dr = 0;
  double drmin= 99;

  while(next_permutation(vi.begin(),vi.end())){//12-34 13-24 14-23
    TLorentzVector v1 = signalJets.at(vi[0]);
    TLorentzVector v2 = signalJets.at(vi[1]);
    TLorentzVector v3 = signalJets.at(vi[2]);
    TLorentzVector v4 = signalJets.at(vi[3]);
    dr = fabs(v1.DeltaR(v2) - 1.0)  + fabs(v3.DeltaR(v4) - 1.0);
    if (dr < drmin){
      index1 = vi[0];
      index2 = vi[1];
      index3 = vi[2];
      index4 = vi[3];
      drmin = dr;
    }
  }
  TLorentzVector stop1 = signalJets[index1] + signalJets[index2];
  TLorentzVector stop2 = signalJets[index3] + signalJets[index4];
  TLorentzVector stop12 = stop1 + stop2;

  bool assigned = false;
  if( (signalJets[index1].pass(BTag77MV2c20) || signalJets[index2].pass(BTag77MV2c20))
      && (signalJets[index3].pass(BTag77MV2c20) || signalJets[index4].pass(BTag77MV2c20)) ) assigned = true;

  //variables definition
  double mass_asymmetry = fabs((stop1.M()-stop2.M()) / (stop1.M() + stop2.M()));
  double mass_average = (stop1.M() + stop2.M())/2.;

  TVector3 boostX = - stop12.BoostVector();
  TLorentzVector stop1inXFrame( stop1 );
  stop1inXFrame.Boost( boostX );
  double costhetastar = fabs(stop1inXFrame.CosTheta());

  //  bool drmin_pass = ((mass_average/1000.) <=225  && drmin<(-0.002*((mass_average/1000.)-255)+0.72))    || ((mass_average/1000.) >225  && drmin<(0.0013*((mass_average/1000.)-255)+0.72));
  bool drmin_pass = ((mass_average<=225) && (drmin<(-0.002*(mass_average-255)+0.72))) || 
                    ((mass_average>225)  && (drmin<(0.0013*(mass_average-255)+0.72)));

  if (mass_asymmetry<0.05 && costhetastar<0.3 && drmin_pass) {
    accept("Inclusive");
    if (assigned) accept("Btagged");
    for(int ii=0; ii<3; ii++) {
      for(int mass = SRranges[ii][0]; mass<=SRranges[ii][1]; mass+=SRranges[ii][2]) {
	if (mass_window(mass_average, mass)) {
	  accept("Inclusive_m"+std::to_string(mass));
	  if (assigned) accept("Btagged_m"+std::to_string(mass));
	}
      }
    }
  }

  return;
}


