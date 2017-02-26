#ifndef ANALYSISCLASS_H
#define ANALYSISCLASS_H

#include "SimpleAnalysis/OutputHandler.h"
#include "TLorentzVector.h"
#include <vector>
#include <algorithm>
#include <functional>

enum AnalysisObjectType { ELECTRON,MUON,TAU,PHOTON,JET,FATJET,MET,COMBINED};

enum AnalysisCommonID { NotBit=1<<31 };

#define NOT(x) (NotBit|x)

enum AnalysisElectronID { EVeryLooseLH=1<<0,
			  ELooseLH=1<<1,
			  EMediumLH=1<<2,
			  ETightLH=1<<3,
			  ELooseBLLH=1<<4,
			  EIsoGradientLoose=1<<8,
			  EIsoBoosted=1<<9, // from 3-bjet paper
			  EIsoFixedCutTight=1<<10,
			  EIsoLooseTrack=1<<11,
			  EIsoLoose=1<<12,
			  EIsoGradient=1<<13,
			  EIsoFixedCutLoose=1<<14,
			  EIsoFixedCutTightTrackOnly=1<<15,
			  ED0Sigma5=1<<16,
			  EZ05mm=1<<17,
			  EGood=EVeryLooseLH|ELooseLH|EMediumLH|ETightLH|ELooseBLLH|ED0Sigma5|EZ05mm,
			  EIsoGood=EGood|EIsoGradientLoose|EIsoBoosted|EIsoFixedCutTight|EIsoLooseTrack|EIsoLoose|EIsoGradient|EIsoFixedCutLoose|EIsoFixedCutTightTrackOnly
};

enum AnalysisMuonID {  MuLoose=1<<0,
		       MuMedium=1<<1,
		       MuTight=1<<2,
		       MuVeryLoose=1<<3,
		       MuIsoGradientLoose=1<<8,
		       MuIsoBoosted=1<<9, // from 3-bjet paper
		       MuIsoFixedCutTightTrackOnly=1<<10,
		       MuIsoLooseTrack=1<<11,
		       MuIsoLoose=1<<12,
		       MuIsoGradient=1<<13,
		       MuIsoFixedCutLoose=1<<14,
		       MuD0Sigma3=1<<16,
		       MuZ05mm=1<<17,
		       MuNotCosmic=1<<18,
		       MuQoPSignificance=1<<19,
		       MuCaloTaggedOnly=1<<20,
		       MuGood=MuLoose|MuMedium|MuTight|MuVeryLoose|MuD0Sigma3|MuZ05mm|MuNotCosmic|MuQoPSignificance,
		       MuIsoGood=MuGood|MuIsoGradientLoose|MuIsoBoosted|MuIsoFixedCutTightTrackOnly|MuIsoLooseTrack
};

enum AnalysisTauID { TauLoose=1<<0,
		     TauMedium=1<<1,
		     TauTight=1<<2,
		     TauOneProng=1<<10,
		     TauThreeProng=1<<11,
		     TauGood=TauLoose|TauMedium|TauTight,
		     TauIsoGood=TauGood};

enum AnalysisPhotonID { PhotonLoose=1<<0,
			PhotonTight=1<<1,
			PhotonIsoFixedCutLoose=1<<8,
			PhotonIsoFixedCutTight=1<<9,
			PhotonGood=PhotonLoose|PhotonTight,
			PhotonIsoGood=PhotonGood|PhotonIsoFixedCutLoose|PhotonIsoFixedCutTight};

enum AnalysisJetID { LooseBadJet=1<<8,
		     TightBadJet=1<<9,
		     JVT50Jet=1<<10,
		     LessThan3Tracks=1<<11, // most signal jets should fail this
		     GoodJet=LooseBadJet|TightBadJet|JVT50Jet,
		     BTag85MV2c20=1<<0,
		     BTag80MV2c20=1<<1,
		     BTag77MV2c20=1<<2,
		     BTag70MV2c20=1<<3,
		     BTag85MV2c10=1<<4,
		     BTag77MV2c10=1<<5,
		     BTag70MV2c10=1<<6,
		     BTag60MV2c10=1<<7,
		     GoodBJet=BTag85MV2c20|BTag80MV2c20|BTag77MV2c20|BTag70MV2c20|
		              BTag85MV2c10|BTag77MV2c10|BTag70MV2c10|BTag60MV2c10|GoodJet,
};

enum AnalysisFatJetID { LooseFatJet=1<<8,
			GoodFatJet=LooseFatJet};

class AnalysisObject : public TLorentzVector
{
 public:
 AnalysisObject(double Px, double Py, double Pz, double E, int charge,
		int id, AnalysisObjectType type, int orgIndex) :
  TLorentzVector(Px,Py,Pz,E), _charge(charge), _id(id), _type(type), _orgIndex(orgIndex) {};
 AnalysisObject(TLorentzVector tlv, int charge,
		int id, AnalysisObjectType type, int orgIndex) :
  TLorentzVector(tlv), _charge(charge), _id(id), _type(type), _orgIndex(orgIndex) {};
  virtual bool pass(int id) const { 
    if (id&NotBit)
      return (id&_id)==0;
    else 
      return (id&_id)==id; 
  } ;
  virtual int charge() const { return _charge; };
  virtual int type() const { return _type; };
  virtual int id() const { return _id; }; // not supposed to be used directly except to store in ntuples
 private:
  int _charge; //not used for jets or photons
  int _id; 
  AnalysisObjectType _type;
  int _orgIndex;
};

AnalysisObject operator+(const AnalysisObject& lhs, const AnalysisObject& rhs);

typedef std::vector<AnalysisObject> AnalysisObjects;

AnalysisObjects operator+(const AnalysisObjects& lhs, const AnalysisObjects& rhs);

class AnalysisEvent
{
public:
  virtual AnalysisObjects getElectrons(float ptCut,float etaCut,int isolation=1)=0;
  virtual AnalysisObjects getMuons(float ptCut,float etaCut,int isolation=1)=0;
  virtual AnalysisObjects getTaus(float ptCut,float etaCut,int isolation=1)=0;
  virtual AnalysisObjects getPhotons(float ptCut,float etaCut,int isolation=1)=0;
  virtual AnalysisObjects getJets(float ptCut,float etaCut,int btag=0)=0;
  virtual AnalysisObjects getFatJets(float ptCut,float etaCut,int btag=0)=0;
  virtual AnalysisObject  getMET()=0;
  virtual float           getSumET()=0;
  virtual ~AnalysisEvent() {};
};

class AnalysisClass;

std::vector<AnalysisClass*> *getAnalysisList(); //for automatically tracking which analyses have been defined

class AnalysisClass
{
public:
 AnalysisClass(const std::string &name) : _name(name),_output(0) { getAnalysisList()->push_back(this); };
 AnalysisClass() {};
  virtual void Init() {};
  virtual void ProcessEvent(AnalysisEvent *event)=0;
  virtual void Final() {};
  virtual ~AnalysisClass() {};
  virtual const std::string& name() {return _name; };
  virtual void setOutput(OutputHandler *output) { _output=output; };
  OutputHandler* getOutput() { return _output; };
  virtual void addRegion(const std::string &label) { _output->addEntry(label); };
  virtual void addRegions(const std::vector<std::string> &labels) { _output->addEntries(labels); };
  virtual void addHistogram(const std::string &label, int bins, float min, float max) { _output->addHistogram(label, bins, min, max); };
  virtual void addHistogram(const std::string &label, int binsX, float minX, float maxX,
			    int binsY,float minY,float maxY) { _output->addHistogram(label, binsX, minX, maxX, binsY, minY, maxY); };
  virtual void accept(const std::string &name,double weight=1) { _output->pass(name, weight); };
  virtual void fill(const std::string &name,double x) { _output->fillHistogram(name, x); };
  virtual void fill(const std::string &name,double x,double y) { _output->fillHistogram(name, x, y); };
  void ntupVar(const std::string &label,int value)   { _output->ntupVar(label,value); };
  void ntupVar(const std::string &label,float value) { _output->ntupVar(label,value); };
  void ntupVar(const std::string &label,double value) { _output->ntupVar(label,(float) value); };
  void ntupVar(const std::string &label,AnalysisObject &object, bool saveMass=false, bool saveType=false);
  void ntupVar(const std::string &label,AnalysisObjects &objects, bool saveMass=false, bool saveType=false);

  // Helper functions
  static AnalysisObjects overlapRemoval(const AnalysisObjects &cands,
					const AnalysisObjects &others,
					float deltaR, int passId=0);

  static AnalysisObjects overlapRemoval(const AnalysisObjects &cands,
					const AnalysisObjects &others,
					std::function<float(const AnalysisObject&,const AnalysisObject&)> radiusFunc, 
					int passId=0);

  static AnalysisObjects filterObjects(const AnalysisObjects& cands,
				       float ptCut,float etaCut=100.,int id=0); 

  static AnalysisObjects filterCrack(const AnalysisObjects& cands,float minEta=1.37,float maxEta=1.52);

  static int countObjects(const AnalysisObjects& cands,
			  float ptCut,float etaCut=100.,int id=0);

  static float sumObjectsPt(const AnalysisObjects& cands,
			    unsigned int maxNum=10000,float ptCut=0);

  static void sortObjectsByPt(AnalysisObjects& cands);
  static float minDphi(const AnalysisObject &met, const AnalysisObjects& cands, unsigned int maxNum=10000, float ptCut=0);
  
  static float calcMCT(const AnalysisObject& o1, const AnalysisObject& o2);
  static float calcMT(const AnalysisObject &lepton, const AnalysisObject &met);
  static float calcMTmin(const AnalysisObjects& cands, const AnalysisObject& met, int maxNum=10000);
  static float calcMT2(const AnalysisObject &o1, const AnalysisObject &o2, const AnalysisObject &met);
  static float calcAMT2(const AnalysisObject &o1, const AnalysisObject &o2, const AnalysisObject &met, float m1, float m2);

  static float aplanarity(const AnalysisObjects& jets);

  static AnalysisObjects reclusterJets(const AnalysisObjects &jets, float radius, float ptmin, float rclus=-1, float ptfrac=-1);

protected:
  std::string _name;
  OutputHandler *_output;

};


#define DefineAnalysis(ANALYSISNAME) \
  class ANALYSISNAME : public AnalysisClass { \
public: \
  ANALYSISNAME() : AnalysisClass(#ANALYSISNAME) {}; \
  void Init(); \
  void ProcessEvent(AnalysisEvent *event); \
  };  \
  static const AnalysisClass * ANALYSISNAME_instance __attribute__((used)) = new ANALYSISNAME(); 
  

#endif
