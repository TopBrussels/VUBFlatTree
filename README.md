# VUBFlatTree installation and setup

### Installation

```
cmsrel CMSSW_9_4_10
cd CMSSW_9_4_10/src
cmsenv
git cms-init

# Egamma
# Electron ID (cut-basedv2 + MVA-based v1)
git cms-merge-topic guitargeek:EgammaID_9_4_X
#EGAMMA smearing
git cms-merge-topic cms-egamma:EGM_94X_v1
git clone https://github.com/ECALELFS/ScalesSmearings.git EgammaAnalysis/ElectronTools/data/ScalesSmearings -b Run2017_17Nov2017_v1

# MET
git cms-merge-topic cms-met:METFixEE2017_949_v2

# Tools needed for AK10 jet collection
git clone https://github.com/cms-jet/JetToolbox JMEAnalysis/JetToolbox 

# Rerun DeepFlavour with newest training
git cms-addpkg RecoBTag/TensorFlow
git cherry-pick 94ceae257f846998c357fcad408986cc8a039152

# Clone this repo
git clone https://github.com/TopBrussels/VUBFlatTree.git

# Compile the monster (use -jN for multicore)
scram b -j5
```
