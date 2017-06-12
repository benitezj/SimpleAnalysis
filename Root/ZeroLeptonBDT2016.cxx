#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(ZeroLeptonBDT2016)

void ZeroLeptonBDT2016::Init()
{
  // BDT based SR's
  addRegions({"SR_BDT_GGdirect","SR_BDT_GGonestep"});

  // Initialize BDT readers
  const std::string xmlFiles[2][2] = {
    // SRBDT-GGdirect
    {"ZeroLepton2016-SRBDT-GGdirect_weights1.xml","ZeroLepton2016-SRBDT-GGdirect_weights2.xml"},
    // SRBDT-GGonestep
    {"ZeroLepton2016-SRBDT-GGonestep_weights1.xml","ZeroLepton2016-SRBDT-GGonestep_weights2.xml"},
  };
  // Input variables for each BDTs {"label","definition"}
  const std::vector< std::vector< std::string > >  variableDefs {
    {
      // Input variables for SRBDT-GGdirect
      { "nJet","meff","METoverMeff3",
	"jetPt0","jetPt1","jetPt2","jetPt3",
	"jetEta0","jetEta1","jetEta2","jetEta3",
	"Sp","Ap","dPhiR"
      },
	// Input variables for SRBDT-GGonestep
      { "nJet","meff","METoverMeff3",
        "jetPt0","jetPt1","jetPt2","jetPt3","jetPt4","jetPt5",
	"jetEta0","jetEta1","jetEta2","jetEta3","jetEta4","jetEta5",
	"Sp","Ap","dPhiR"
      },
    }
  };

  const std::string methodname = "BDTG";
  for( unsigned int i=0 ; i<sizeof(xmlFiles)/sizeof(xmlFiles[0]) ; i++ ){
    BDT::BDTReader* reader = new BDT::BDTReader(methodname,variableDefs[i],xmlFiles[i][0],xmlFiles[i][1]);
    m_BDTReaders.push_back(reader);
  }

}

void ZeroLeptonBDT2016::ProcessEvent(AnalysisEvent *event)
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

  // Meff based analysis regions and selections
  float meff[7];
  for(int nJet=2; nJet<=6; nJet++)
    meff[nJet] = sumObjectsPt(jets, nJet) + met;
  float dphiMin3    = minDphi(metVec, jets, 3);
  float dphiMinRest = minDphi(metVec, jets);
  float Ap          = aplanarity(jets);
  float Sp          = sphericity(jets);


  // BDT-based SR's
  // SR-BDT-GGdirect
  if( Njets >= 4 && meffIncl > 1400. && met > 300. && 
      jets[0].Pt() > 200. && jets[3].Pt() > 50. && 
      dphiMin3 > 0.4 && met/meff[4] > 0.2 ){
    // get BDT score
    BDT::BDTReader * reader = m_BDTReaders[0];
    reader->initInEventLoop();
    reader->setValue("meff"   ,meffIncl);
    reader->setValue("METoverMeff3", met/meff[4]);
    reader->setValue("jetPt0" ,jets[0].Pt());
    reader->setValue("jetPt1" ,jets[1].Pt());
    reader->setValue("jetPt2" ,jets[2].Pt());
    reader->setValue("jetPt3" ,jets[3].Pt());
    reader->setValue("jetEta0",jets[0].Eta());
    reader->setValue("jetEta1",jets[1].Eta());
    reader->setValue("jetEta2",jets[2].Eta());
    reader->setValue("jetEta3",jets[3].Eta());
    reader->setValue("dPhiR"  ,dphiMinRest  );
    reader->setValue("Ap"     ,Ap           );
    // add for A0
    reader->setValue("Sp"     ,Sp           );
    reader->setValue("nJet"   ,Njets        );
    float BDT = reader->getBDT();

    if(BDT>0.90) accept("SR_BDT_GGdirect");
  }
 
  // SR-BDT-GGonestep
  if( Njets >= 4 && meffIncl > 1400. && met > 300. && 
      jets[0].Pt() > 200. && jets[3].Pt() > 50. && 
      dphiMin3 > 0.4 && dphiMinRest > 0.2 && met/meff[4] > 0.2 ){
    // get BDT score
    BDT::BDTReader * reader = m_BDTReaders[1];
    reader->initInEventLoop();
    reader->setValue("meff"   ,meffIncl);
    reader->setValue("METoverMeff3", met/meff[4]);
    reader->setValue("jetPt0" ,jets[0].Pt() );
    reader->setValue("jetPt1" ,jets[1].Pt() );
    reader->setValue("jetPt2" ,jets[2].Pt() );
    reader->setValue("jetPt3" ,jets[3].Pt() );
    reader->setValue("jetPt4" ,(jets.size()>4) ? jets[4].Pt() : 0. );
    reader->setValue("jetPt5" ,(jets.size()>5) ? jets[5].Pt() : 0. );
    reader->setValue("jetEta0",jets[0].Eta() );
    reader->setValue("jetEta1",jets[1].Eta() );
    reader->setValue("jetEta2",jets[2].Eta() );
    reader->setValue("jetEta3",jets[3].Eta() );
    reader->setValue("jetEta4",(jets.size()>4) ? jets[4].Eta() : -3. );
    reader->setValue("jetEta5",(jets.size()>5) ? jets[5].Eta() : -3. );
    reader->setValue("dPhiR"  ,dphiMinRest  );
    reader->setValue("Ap"     ,Ap           );
    // add for B0
    reader->setValue("Sp"     ,Sp           );
    reader->setValue("nJet"   ,Njets        );
    float BDT = reader->getBDT();

    if(BDT>0.92) accept("SR_BDT_GGonestep");
  }

  return;
}
