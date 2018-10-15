# VUBFlatTree installation and setup

### Installation

```
cmsrel CMSSW_9_4_10
cd CMSSW_9_4_10/src
cmsenv
git cms-init

# Egamma
git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP
git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3
git cms-merge-topic cms-egamma:EGM_94X_v1
git clone https://github.com/ECALELFS/ScalesSmearings.git EgammaAnalysis/ElectronTools/data/ScalesSmearings -b Run2017_17Nov2017_v1

# MET
git cms-merge-topic cms-met:METFixEE2017_949_v2

# Compile
scram b -j5

cd $CMSSW_BASE/external/slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/src
git cms-addpkg RecoEgamma/ElectronIdentification
cp ${CMSSW_BASE}/external/slc6_amd64_gcc630/data/RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_E* ${CMSSW_BASE}/src/RecoEgamma/ElectronIdentification/data/Fall17/.

# Tools needed for AK10 jet collection
git clone https://github.com/cms-jet/JetToolbox JMEAnalysis/JetToolbox 

# Clone this repo
git clone https://github.com/TopBrussels/VUBFlatTree.git

# Compile the monster (use -jN for multicore)
scram b -j5
```
