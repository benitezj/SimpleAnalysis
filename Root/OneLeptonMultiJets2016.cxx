#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(OneLeptonMultiJets2016)

void OneLeptonMultiJets2016::Init()
{
	addRegions({"preSelection","8j80_0b","8j80_3b","10j60_3b","10j80_3b"});
}

void OneLeptonMultiJets2016::ProcessEvent(AnalysisEvent *event)
{
	auto baseElectrons = event->getElectrons(10,2.47,ELooseBLLH);
	auto baseMuons     = event->getMuons(10,2.4);
	auto baseJets      = event->getJets(20.,2.8);
	auto baseLeptons   = baseElectrons + baseMuons;

	// overlap removal
	baseElectrons  	= overlapRemoval(baseElectrons, baseMuons, 0.01);
	auto Jets       = overlapRemoval(baseJets, baseElectrons, 0.2, NOT(BTag85MV2c10));
	auto muJetSpecial = [] (const AnalysisObject& jet, const AnalysisObject& muon) { 
	  if (jet.pass(NOT(BTag85MV2c10)) && (jet.pass(LessThan3Tracks) || muon.Pt()/jet.Pt()>0.5)) return 0.4;
	  else return 0.;
	};
	Jets            = overlapRemoval(Jets, baseMuons, muJetSpecial, NOT(BTag85MV2c10));
	auto Electrons  = overlapRemoval(baseElectrons, Jets, 0.4);
	auto Muons      = overlapRemoval(baseMuons, Jets, 0.4);

	// signal objects
	auto signalElectrons = filterObjects(Electrons,30,2.47,ETightLH && EIsoGradient);
	auto signalMuons     = filterObjects(Muons,30,2.4,MuMedium && MuIsoGradient);
	auto signalLeptons   = signalElectrons + signalMuons;
	auto signalJets40    = filterObjects(Jets,40.,2.4);
	auto signalJets60    = filterObjects(Jets,60.,2.4);
	auto signalJets80    = filterObjects(Jets,80.,2.4);
	auto bjets40         = filterObjects(Jets,40.,2.4,BTag77MV2c10);
	auto bjets60         = filterObjects(Jets,60.,2.4,BTag77MV2c10);
	auto bjets80         = filterObjects(Jets,80.,2.4,BTag77MV2c10);

	// Signal regions
	if ( signalLeptons.size()>=1 && signalJets40.size()>=5 ) {
		accept("preSelection");
		if ( bjets60.size()>=3 ) {
			if ( signalJets60.size()>=10 ) accept("10j60_3b");}

		if ( bjets80.size()>=3 ) {
			if ( signalJets80.size()>=8 )  accept("8j80_3b");
			if ( signalJets80.size()>=10 ) accept("10j80_3b");}
		else if ( bjets80.size()==0 ) {
			if ( signalJets80.size()>=8 )  accept("8j80_0b");
		}
	}
	return;
}
