#include "xAODRootAccess/TEvent.h"
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "MCTruthClassifier/MCTruthClassifier.h"
#include "SimpleAnalysis/AnalysisClass.h"
#include "SimpleAnalysis/AnalysisRunner.h"

using std::vector;

class xAODTruthReader : public Reader {

 public:
 xAODTruthReader(std::vector<AnalysisClass*>& analysisList);
  ~xAODTruthReader();

 protected:
  bool processEvent(xAOD::TEvent *event,xAOD::TStore *store);
  void processFilesInternal(const std::vector<std::string>& inputNames);
  int getTruthType(const xAOD::TruthParticle *part);
  xAOD::TruthParticleContainer* findTruthParticles(xAOD::TStore *store,
						   const xAOD::TruthParticleContainer* truthparticles,
						   int pdgId, int status=1);

 private:
  ST::SUSYObjDef_xAOD *_susytools;
  MCTruthClassifier* _mctool;
  xAOD::TEvent* _event;
};
