# Simplified ATLAS SUSY analysis framework

Holds a collections of SUSY analyses. These can be run over samples in different
input formats:  
*  DAOD_TRUTH      (TRUTH1 and TRUTH3 tested)  
*  xAOD            (either truth-level and or reco-level - the latter with some constraints)  
*  slimmed ntuples (Reduced ntuples produced from above input)

It provides the analysis acceptance per control and signal region as well as 
optionally histograms or ntuples with event level objects. Smearing of the truth-level objects 
can optionally be done using the latest upgrade performance smearing functions.

## Compiling:
It should compile out-of-the-box on top of any recent AnalysisBase release 
(tested in AnalysisBase,2.4.28):
```
 rcSetup Base,2.4.28
 #following line gets Ext_RestFrames RootCore package
 git clone https://github.com/crogan/Ext_RestFrames.git
 #The following line is optional and is only needed if one wants to enable smearing 
 rc checkout_pkg $SVNOFF/PhysicsAnalysis/UpgradePhys/SmearingFunctions/UpgradePerformanceFunctions/tags
 #The following line is optional and is only needed if one wants to enable b-jet weighting (for 3-bjet analysis) 
 rc checkout_pkg $SVNOFF/PhysicsAnalysis/JetTagging/JetTagPerformanceCalibration/BTaggingTruthTagging/tags
 rc find_packages
 rc compile
```

## To run:  
```
simpleAnalysis [-a listOfAnalysis] <inputFile1> [inputFile2]...
```
This will run the analyses specified with "-a" option (or all if not given)
over all of the input files and provide acceptances in a text file for each 
analysis ("analysisName.txt") and histograms in a root file ("analysisName.root").
Ntuples can be activated by adding "-n", while the different outputs can be 
merged in single text and root files by specifying "-o <name>". 
Smearing of the truth-level objects is enabled with the option "-s mu=200".
It is foreseen that later different smearing functions will be supported.
If "-r" is specified, it will try to use reco-level objects instead of truth-level.
This allows to compare directly to the standard analysis results, though no
trigger efficiencies etc. are applied and overall this is still experimental.
Events are normally weighted by their event-by-event weights. In case of samples
with multiple weights, for instance from PDF variations, the applied event weight
can be chosen by "-w <number>". The weighting can be disabled by "-w -1".
Finally the available analysis can be seen with the "-l" option.

### Submission on the grid:
Submission to the grid is straightforward:
```
lsetup panda
rc clean
prun --exec 'simpleAnalysis -a ZeroLepton2015,ThreeBjets2015 %IN' --outputs='*.txt,*.root' --useRootCore --inDS <inputDS> --outDS user.<userName>.simple.v1
```

## Slimming:
For running over large input files more than once, it can be advantageous to
first slim the files to an absolute minimum by making an ntuple with all the
input objects. This can be done trivially with `slimMaker`. Minimum object pTs
can be specified on the commandline, see `slimMaker --help`. The output can
supplied to `simpleAnalysis` in the same way as DAODs and the program will
automatically detect the type of input.

## Adding a new analysis:  
Simply make new C++ routine in `Root/MyAnalysisName.cxx` with the structure:  
```
#include "SimpleAnalysis/AnalysisClass.h"  
DefineAnalysis(MyAnalysisName)  
void MyAnalysisName::Init() {}  
void MyAnalysis::ProcessEvent(AnalysisEvent *event)  
```

Recompilation will automatically make the analysis available to run. See 
`Root/ExampleAnalysis.cxx` and the other supplied analysis examples on 
how to specify signal regions, kinematic selections and filling 
histograms/ntuples.

It is suggested to use analysis names that ends in the year of the latest used
data, i.e. "2015" or "2016", and to add "Conf" in front if it is a preliminary
analysis. Examples: `ZeroLepton2016.cxx` or `StopOneLeptonConf2016.cxx`

For more help, please contact @aagaard
