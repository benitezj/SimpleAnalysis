#include <iostream>
#include <sstream>
#include <algorithm>

#include <TFile.h>
#include <TChain.h>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "SimpleAnalysis/SlimReader.h"
#include "SimpleAnalysis/OldSlimReader.h"
#include "SimpleAnalysis/xAODTruthReader.h"
#include "SimpleAnalysis/xAODRecoReader.h"
#include "SimpleAnalysis/AnalysisClass.h"
#include "SimpleAnalysis/OutputHandler.h"
#include "SimpleAnalysis/AnalysisRunner.h"
#include "SimpleAnalysis/TruthSmear.h"
#include "SimpleAnalysis/PDFReweight.h"

static void splitCommaString(const std::string& names,std::vector<std::string>& result) {
  std::stringstream ss(names);
  while( ss.good() ) {
    std::string substr;
    getline( ss, substr, ',' );
    result.push_back( substr );
  }
}

int main(int argc, char **argv) {
  if ( ! xAOD::Init().isSuccess() ) {
    throw std::runtime_error("Cannot initialise xAOD access !");
  }

  po::options_description desc("Run one or more truth-level analyses");
  std::string outputName;
  desc.add_options()
    ("help,h", "print usage and exit")
    ("output,o",    po::value<std::string>(&outputName), "Output name - if not supplied use analyses names")
    ("analyses,a",  po::value<std::string>(), "Comma-separated list of analyses to run")
    ("listanalyses,l", "List available analyses and exit")
    ("input-files", po::value< vector<std::string> >(), "Comma-separated list of input files")
    ("ntuple,n", "Fill ntuple")
    ("readReco,r", "Use reconstructed quantities instead of truth")
    ("pdfReweight,p",  po::value<std::string>(), "PDF reweight to '<pdfName>[,energyIn,energyOut]'")
    ("smear,s", po::value<std::string>(), "Comma-separated list smearing options (use help to see full list of options)")
    ("mcweight,w", po::value<int>()->default_value(0), "MC weight index to apply (set to -1 to ignore it, i.e. =1.)")
    ;
  po::positional_options_description p;
  p.add("input-files", -1);
  
  po::variables_map vm;
  po::store( po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);
  
  if (vm.count("listanalyses")) {
    std::cout<<"Select analyses among:"<<std::endl;
    for(auto *analysis : *getAnalysisList()) {
      std::cout<<" "<<analysis->name()<<std::endl;
    }
    return 0;
  }
  
  if (vm.count("help") || (vm.count("input-files")==0)) {
    std::cout << desc << std::endl;
    return 1;
  }
  bool doNtuple = vm.count("ntuple");
  bool mergedOutput = vm.count("output");

  TruthSmear *smearer=0;
  if (vm.count("smear")) {
    std::vector<std::string> smearingOptions;
    splitCommaString(vm["smear"].as<std::string>(),smearingOptions);
    smearer=new TruthSmear(smearingOptions);
  }

  PDFReweighter *pdfReweighter=0;
  if (vm.count("pdfReweight")) {
    std::vector<std::string> reweightOptions;
    splitCommaString(vm["pdfReweight"].as<std::string>(),reweightOptions);
    if ((reweightOptions.size()!=1) && (reweightOptions.size()!=3)) {
      std::cout<<"Please specify PDF name and optionally input and output CM energies"<<std::endl;
      return 1;
    }
    pdfReweighter=new PDFReweighter(reweightOptions);
  }

  int mcwindex = 0;
  if (vm.count("mcweight")) {
    mcwindex = vm["mcweight"].as<int>();
  }
  std::cout << "MCWeightIndex = " << mcwindex << (mcwindex>=0 ? "" : " . No MC weight will be applied.") << std::endl; 
  

  std::vector<std::string> inputFileNames;
  std::vector<std::string> analysisNames;
  
  if (vm.count("analyses")) splitCommaString(vm["analyses"].as<std::string>(),analysisNames);
  for(const auto& fileNames: vm["input-files"].as<vector<std::string> >()) {
    splitCommaString(fileNames,inputFileNames);
  }
  
  std::cout<<"Analyses to run: ";
  for(const auto& name: analysisNames) std::cout<<name<<" ";
  if (analysisNames.size()==0)  std::cout<<"all";
  std::cout<<std::endl;
  std::cout<<"Files to analyze: ";
  for(const auto& fileName: inputFileNames) std::cout<<fileName<<" ";
  std::cout<<std::endl;
  if (mergedOutput) 
    std::cout<<"Output merged into: "<<outputName<<".[txt|root]"<<std::endl;
  else
    std::cout<<"Output split per analysis"<<std::endl;

  TFile *oRoot = 0;
  if (mergedOutput) oRoot=new TFile((outputName+".root").c_str(),"RECREATE");

  std::vector<AnalysisClass*> analysisList;
  bool selectAnalysis=analysisNames.size()>0;
  for(auto *analysis : *getAnalysisList()) {
    if (selectAnalysis) {
      const auto namePtr=std::find(analysisNames.begin(),analysisNames.end(),analysis->name());
      if (namePtr==analysisNames.end()) continue;
      analysisNames.erase(namePtr);
    }
    if (!mergedOutput) {
      oRoot = new TFile((analysis->name()+".root").c_str(),"RECREATE");
    }
    OutputHandler* output=new OutputHandler(oRoot,doNtuple);
    if (outputName.size()) output->title(analysis->name());
    analysis->setOutput(output);
    analysisList.push_back(analysis);
  }
  if (analysisNames.size()!=0) {
    std::cerr<<"Unknown analysis requested: ";
    for(const auto& name: analysisNames) std::cerr<<name<<" ";
    std::cerr<<std::endl;
    return 1;
  }
  
  std::sort(analysisList.begin(), analysisList.end(), [](AnalysisClass* a, AnalysisClass* b) {
      return b->name() > a->name();   
    });
  
  TFile *fh=TFile::Open(inputFileNames[0].c_str());
  if (fh==0) {
    std::cerr<<"Failed to open the first file: "<<inputFileNames[0]<<std::endl;
    return 2;
  }
  Reader *reader=0;
  if (fh->FindKey("truth")) {
    std::cout<<"Reading (old-style) slimmed input"<<std::endl;
    reader=new OldSlimReader(analysisList);
  } else if (fh->FindKey("ntuple")) {
    std::cout<<"Reading slimmed input"<<std::endl;
    reader=new SlimReader(analysisList);
  } else if (fh->FindKey("CollectionTree")) {
    std::cout<<"Reading xAOD input"<<std::endl;
    if (vm.count("readReco")) {
      std::cout<<" using reconstructed quantities"<<std::endl;
      reader=new xAODRecoReader(analysisList);
    } else 
      reader=new xAODTruthReader(analysisList);
  } else {
    std::cerr<<"Unknown input format in: "<<inputFileNames[0]<<std::endl;
    return 2;
  }
  reader->SetSmearing(smearer);
  reader->SetReweighting(pdfReweighter);
  reader->SetMCWeightIndex(mcwindex);
  reader->processFiles(inputFileNames);
  delete reader;
  delete smearer;

  std::ofstream oFile;
  if (mergedOutput) oFile.open(outputName+".txt");
  bool first=true;
  for(const auto& analysis : analysisList) {
    if (!mergedOutput) {
      oFile.close();
      oFile.open(analysis->name()+".txt");
    }
    analysis->getOutput()->saveRegions(oFile,first);
    if (mergedOutput) first=false;
  }

  return 0;
}
