# B2GTTrees
============
INFO: Filter and Convert B2G edm ntuples to TTree ntuples

N.B (v8.0.x_v3.1_May10 only): This is a private tag for incporporating the latest MET corrections for EG issues and
switching to SUSY recommended EA values for ele/muon, otherwise use TAG: v8.0.x_v3.1_May09
I also switch to the Spring15 electron IDs privately in the beginning of Analysis/B2GTTrees/test/B2GEdmToTTreeNtupleExtra_cfg.py
(Edit them if necessary to specify for your analysis)
Everything else is the same as in B2GAnaFW tag v8.0.x_v3.1

Checkout Instructions (with B2GAnaFW)
=====================================

```Shell
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
cmsenv
git cms-init
git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTaggerV4-WithWeightFiles-v1_from-CMSSW_8_0_21
git cms-merge-topic cms-met:METRecipe_8020
git cms-merge-topic ikrav:egm_id_80X_v3_photons
git clone https://github.com/cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_80X_V3
git clone https://github.com/jkarancs/B2GAnaFW.git  Analysis/B2GAnaFW -b v8.0.x_v3.1_May10
git clone https://github.com/jkarancs/B2GTTrees.git Analysis/B2GTTrees -b v8.0.x_v3.1_May10
scram b -j 20
```

Running
=======

First copy sqlite files so you can run conveniently from the src directory
```Shell
ln -s Analysis/B2GAnaFW/test/JECs/*.db .
```

There are two (supported) ways to produce ntuples:

   * MiniAOD -> Common B2G edm ntuple -> Add ExtraVars (+Filter) -> TTree ntuple
```Shell
cmsRun Analysis/B2GTTrees/test/B2GEdmToTTreeNtupleExtra_cfg.py sample="/store/data/Run2016C/JetHT/MINIAOD/03Feb2017-v1/110000/000EC584-5BEB-E611-90BB-008CFA111230.root" outputLabel="B2GTTreeNtupleExtra_Data_80X_JetHT_Run2016C.root" DataProcessing="Data_80X_Run2016BCD_03Feb2017" wantSummary=False maxEvents=1000

cmsRun Analysis/B2GTTrees/test/B2GEdmToTTreeNtupleExtra_cfg.py sample="/store/mc/RunIISummer16MiniAODv2/TT_TuneEE5C_13TeV-powheg-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/0045D9EF-67B5-E611-B736-FA163EEF13C9.root" outputLabel="B2GTTreeNtupleExtra_MC_80X.root" DataProcessing="MC_MiniAODv2_80X_Summer16" wantSummary=False maxEvents=1000
```

CRAB3 Tools
==========
INFO: You can mass produce B2GEdmNtuples similar to multicrab

Usage - Setup:
- tcsh:
```Shell
cd Analysis/B2GTTrees/test/crab3
cmsenv
source /cvmfs/cms.cern.ch/crab3/crab.csh
voms-proxy-init --voms cms -valid 168:00
alias merged 'source crab3_B2GTTreeNtuples_Merged.csh \!*'
merged --help
```

- bash:
```Shell
cd Analysis/B2GTTrees/test/crab3
cmsenv
source /cvmfs/cms.cern.ch/crab3/crab.sh
alias merged='./crab3_B2GTTreeNtuples_Merged.csh'
merged --help
```

Examples:

   * Make TTree ntuples from MINIAOD(SIM) on grid
```Shell
cd $CMSSW_BASE/src/Analysis/B2GTTrees/test/crab3
merged create May10 B2GAnaFW_80X_V3p1_May10 MINIAODv2_80X_May10_input.txt cross_sections.txt T2_HU_Budapest /store/user/jkarancs/SusyAnalysis/B2GTTreeNtuple/May10
merged submit May10
merged status May10
merged status May10 --run #if you want automatic job handling, otherwise it prints resubmit commands
merged download May10 /data/gridout/jkarancs/SusyAnalysis/B2G/TTreeNtuple --run
merged find_missing May10
merged checkup May10
merged report May10 --run
merged get_lumi May10
```
