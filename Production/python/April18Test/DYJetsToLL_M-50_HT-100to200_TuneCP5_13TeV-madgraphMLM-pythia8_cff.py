import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D24733BB-BC41-E811-9681-0CC47AB64EF4.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/AEB933FE-BC41-E811-B454-AC1F6B23C814.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D02BDD51-BD41-E811-8272-AC1F6B23C9D8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/6E7BDB56-BD41-E811-9230-0025904CFB8A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/C0C1FEFF-BC41-E811-85D1-AC1F6B23C93A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/10880254-BD41-E811-A959-0CC47AB64AE0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/BC01225A-BD41-E811-9063-0CC47A545096.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/DAAF7F58-BD41-E811-BBA7-B4E10FD217E1.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/AA3B738F-BE41-E811-84E5-AC1F6B23C814.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7463C05B-BD41-E811-BA20-0CC47AB64AE0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/5E603ABC-BF41-E811-94CF-B4E10FA31EFB.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/067373A8-C141-E811-9B22-0CC47AE0F33A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/208105F2-C241-E811-A98A-0CC47AE0F33C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/00B6929E-C441-E811-BD5E-B4E10FD217E1.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/A401E54C-C541-E811-AFF6-AC1F6B23C970.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/EAA43DAB-C541-E811-91F8-AC1F6B23C814.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/76A0154C-C541-E811-A5E3-0025904B3092.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/8EC8C19A-C641-E811-9FA9-0CC47AE0F33A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D09981F8-C741-E811-B1A0-AC1F6B239D7A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/00440BA2-C941-E811-9C0E-0025904B3092.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/24CF38B4-CB41-E811-B49A-0025904CF972.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D643F9C6-CC41-E811-96E9-AC1F6B23C93A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/B62F202C-CF41-E811-A804-AC1F6B23C88C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/4ADB1C96-CE41-E811-819F-AC1F6B239D7A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/88DDC285-D041-E811-A738-AC1F6B239D7A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F0EC9A7B-D141-E811-8609-0025904B3092.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/96F9BA66-D241-E811-970A-B4E10FA31EE1.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/C6C6FEB6-D241-E811-8E53-0CC47AE0F33C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/4AFFD304-D541-E811-B820-AC1F6B23C94C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/60D5B8B6-D441-E811-9EEE-AC1F6B23C830.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/EC40EF13-D641-E811-A116-AC1F6B23C94C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F6914413-D641-E811-A06A-AC1F6B23C830.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/02F70E18-D741-E811-A79D-0CC47AE0F33A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/508F1287-D841-E811-90DA-AC1F6B23C86A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/FA70CBED-D741-E811-A94E-AC1F6B23C94C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/5011976F-DA41-E811-9EDB-AC1F6B23C94C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/18814571-DE41-E811-A6FD-AC1F6B23C93A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/5CBC1A2B-DF41-E811-9760-AC1F6B239D7A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/82DC03F6-E441-E811-864C-AC1F6B23C7DC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E6C767F2-E441-E811-A5DD-0025904CF976.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/02CDB7F4-E441-E811-9A9E-AC1F6B23C814.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/A8FE91F2-E441-E811-93EF-AC1F6B239D7A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/843ABBF2-E441-E811-BE14-AC1F6B23C86A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/2A6BFFFA-E441-E811-B1F6-AC1F6B23C93C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/72981CF0-E441-E811-A9C6-B4E10FA31FCF.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/2E7B0FF0-E441-E811-9EB9-0CC47AB64EF4.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/6A976EEA-E741-E811-898D-0025904B2FFC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/FE7A81E2-EC41-E811-BB95-AC1F6B23C7DE.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E22C08DA-EC41-E811-ADA2-AC1F6B23C970.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E2A356E7-EC41-E811-A784-AC1F6B23C77A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/64FA27D9-EC41-E811-8ADF-0CC47AB64E80.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/BC4A73E0-EC41-E811-9FDE-B4E10FD218CF.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D80E4BDD-EC41-E811-B3B6-AC1F6B239D7A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/72570BE5-EC41-E811-A523-0CC47AB64FB0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E4FDE8E4-EC41-E811-BDF1-AC1F6B23C9D6.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/4E67A6E0-EC41-E811-ADA1-0CC47AB58BE4.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/2AB7FDDF-EC41-E811-A767-0025904CEC20.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7A19DBE1-EC41-E811-8C15-AC1F6B23C96E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/2A863CE1-EC41-E811-B72E-0025904CEC20.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/124842DC-EC41-E811-85DD-AC1F6B23C9D6.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/5AE720DF-EC41-E811-95BE-0CC47AB64FB0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/EC606FDA-EC41-E811-83F1-0025904CEC20.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/286920DB-EC41-E811-9F12-AC1F6B23C94C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/66A214D9-EC41-E811-9260-0CC47AB58BE4.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E8C61BDB-EC41-E811-A199-AC1F6B23C94C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/4E778CDB-EC41-E811-BA53-AC1F6B23C96E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/A84ED0D9-EC41-E811-ADEF-AC1F6B23C96E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/04700ADC-EC41-E811-9750-AC1F6B23C96E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7247A9E3-EC41-E811-8E8D-0CC47AB64FB2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7E3F15D8-EC41-E811-9193-0CC47AE0F33C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/2C78B4F3-F041-E811-9B5E-AC1F6B23C93C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/DE7C09E8-F041-E811-A7A7-B4E10FD218CF.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E2766DEA-F041-E811-BAAD-AC1F6B23C7DE.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7EEB0FEC-F041-E811-8DAB-AC1F6B23C9D8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/00F3E6EA-F041-E811-A33A-0CC47AE0F33C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E681C4E8-F041-E811-8C98-0025904E239C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/CC85C7EC-F041-E811-964B-0025904B2FFC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D441CDA9-F141-E811-8683-0CC47AB65044.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/70431D5C-F241-E811-8491-AC1F6B23C9D8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/74FF76A7-F141-E811-B882-0CC47AB58BE4.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/ACF43065-F241-E811-9F4C-AC1F6B23C9D6.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/DA577B60-F341-E811-9A27-B4E10FD218CF.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/98475060-F341-E811-B71E-0CC47AB64FB2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F26D75FB-F341-E811-B9C4-AC1F6B23C77A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E0C531F9-F341-E811-85E4-0CC47AB58BE4.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/30E8A1FD-F341-E811-BFD7-AC1F6B239D78.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E84E60B1-F441-E811-8ED3-AC1F6B23C9D6.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/C65961BC-F541-E811-84D9-0CC47AB64FB2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/4C4A5975-F641-E811-95C1-AC1F6B23C93C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/6AC5607A-F641-E811-A5D7-0CC47A545096.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/08EC26E0-F641-E811-9F17-0025904CFB8A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/AED1F056-F841-E811-9816-AC1F6B23C93C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/2C016A6C-F841-E811-9BF6-0CC47AB64AE2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/C06D7ADE-F741-E811-8832-0CC47AB64FB2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/70DB7F6F-F841-E811-B280-B4E10FA326F5.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/CE3B9814-FA41-E811-9150-0CC47AB64FB2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/0E3FEE42-8242-E811-904D-0CC47A544E12.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/7A8D3108-CA41-E811-B274-AC1F6B239D7A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/48CDD289-CB41-E811-8E2A-0CC47AB64AE2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/08CD9834-8142-E811-9578-AC1F6B23C88C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/E43E4210-7742-E811-9430-AC1F6B23C96A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/426EC946-8B42-E811-AC29-0CC47AB58BE6.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/D69A2A18-9342-E811-AE6C-B4E10FA31EFB.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/1057BEEB-9942-E811-B2F3-AC1F6B239D7A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/3C8C7DE9-9942-E811-BFE5-0CC47AB64E82.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/D21E97E8-9942-E811-9F04-B4E10FD218CF.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/7C9746E9-9942-E811-9BD9-0CC47AB64E82.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/26111CEB-9942-E811-A146-B4E10FD217E1.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/085CD0EB-9942-E811-ACB7-AC1F6B23C93A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/2455B4E8-9942-E811-908C-B4E10FA32587.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/EC224AEF-9942-E811-832C-AC1F6B23C592.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/260A72EC-9942-E811-8D02-AC1F6B23C7DA.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/50D548E9-9942-E811-85CF-B4E10FD218CF.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/6287D3EA-9942-E811-99F3-AC1F6B23C96A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/528F30FF-9942-E811-8BF4-0025904E23AA.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/9EB4CFE9-9942-E811-B4DE-0CC47AB64FB0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/A075EDEF-9942-E811-B69E-0CC47AB64EF4.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/D41493EA-9942-E811-97FB-0CC47AB64AE2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/E0AFCEE9-9942-E811-BF72-0CC47AB65044.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/D4CDA5EC-9942-E811-8005-AC1F6B23C832.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/F63884EA-9942-E811-88D8-0CC47AB65044.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/22E2FAF3-9942-E811-9AD4-AC1F6B23C832.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/FCFAABEB-9942-E811-B2F1-AC1F6B23C96A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/CC0EAF3B-9A42-E811-9415-0CC47AB58E9E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/3A48613C-9A42-E811-9137-0CC47AB58BE6.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/26249342-9A42-E811-AD91-AC1F6B23C9D8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/1293F23C-9A42-E811-B5E2-0CC47AB65046.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/94AEEC90-9A42-E811-8AD4-AC1F6B23C93C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/0CA60A90-9A42-E811-81D4-AC1F6B23C592.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/8E5B793A-9A42-E811-AE67-0CC47AB6503E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/70D46890-9A42-E811-B7E7-0CC47A54529A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/98E954B1-9B42-E811-8629-B4E10FA32587.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/3EA686B8-9B42-E811-99B4-AC1F6B23C832.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/5E5DD3BC-9B42-E811-AC35-AC1F6B23C7DA.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/E265F0B8-9B42-E811-8236-B4E10FD218CF.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/9C0982C4-9B42-E811-8DD5-AC1F6B23C93A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/F27C02D5-9B42-E811-ADCD-AC1F6B23C592.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/4C60BEDB-9B42-E811-B14A-AC1F6B23C96A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/1E3DECDB-9B42-E811-B553-AC1F6B23C96A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/2A008C68-9B42-E811-B8AC-0025904E2382.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/6C5CA86C-9B42-E811-BB54-AC1F6B23C88C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/DCEC1B6C-9B42-E811-B5DA-0025904E2382.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/32B57560-9D42-E811-9947-AC1F6B23C832.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/464F0767-9D42-E811-9DBD-AC1F6B239D7A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/209BD601-9F42-E811-9B54-0CC47AB6503C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/66B489A2-A042-E811-A293-0025904E2382.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/B85FAC77-A142-E811-BF6A-AC1F6B23C96E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/D88EF69F-A342-E811-A1AF-AC1F6B23C96A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/744891D3-A242-E811-8EDF-0CC47A5450DA.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/F24036EE-A342-E811-905E-AC1F6B23C7DA.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/C6FBE950-A442-E811-A558-0CC47A54529A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/E0D1B129-A542-E811-8011-0CC47AB65046.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/FACC33BE-A642-E811-AC9D-0CC47A5450F8.root',
] )