import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/2090CAB5-83FA-E511-9FCC-0002C94CD044.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/368E42AC-83FA-E511-A3F8-24BE05C648A1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/38B94CB6-83FA-E511-A46D-0002C94DBB18.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/42F6BDB2-83FA-E511-ACF5-0002C94CDAF4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/5A4E39B3-83FA-E511-9635-0002C94CDAF4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/6EFC79AB-83FA-E511-8316-24BE05C616E1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/72039DAF-83FA-E511-A67A-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/72AEAFA5-84FA-E511-9A11-0002C94DBB18.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/7CCC549C-84FA-E511-9473-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/984B65C4-83FA-E511-8420-5065F3818291.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/9E1DAAB1-83FA-E511-968A-0002C94D5612.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/A4D5F5B6-83FA-E511-8EAB-0002C94CD9E2.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/AACFB3C0-84FA-E511-8FC9-24BE05CE1E01.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/B0A96DB4-83FA-E511-8E19-0002C94CDAF4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/B2B148B6-83FA-E511-876E-0002C94DBB18.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/CA4817AC-83FA-E511-82BB-24BE05CEFB41.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/D4F367B0-83FA-E511-934A-0002C94CD116.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/D68EAAB7-83FA-E511-B7B5-0002C94D54FA.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/20000/E804D4B5-83FA-E511-8447-0002C94CD044.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/1A7D5926-75FA-E511-A110-0002C94CD1E4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/22B3FC1C-75FA-E511-94C1-24BE05C63651.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/34FDA021-75FA-E511-AD92-24BE05CEECD1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/42714429-75FA-E511-9916-002590FD5A72.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/488D8D1D-75FA-E511-AEDB-5065F381E201.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/584E2D20-75FA-E511-8ABF-B8CA3A709648.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/62650D1F-75FA-E511-A3CB-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/6293E31F-75FA-E511-A226-B8CA3A708F98.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/70BDE328-75FA-E511-BCC6-24BE05C4D8C1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/729C8B1D-75FA-E511-AA5A-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/8816CA1C-75FA-E511-A3DA-24BE05BDBE61.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/96F44627-75FA-E511-BF1D-0002C94D5504.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/9A1D9123-75FA-E511-BD69-24BE05C3CBD1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/ACD3FC26-75FA-E511-8F0D-0002C94D5504.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/AEC9A51C-75FA-E511-8979-5065F381A2F1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/C6588F22-75FA-E511-BAE0-24BE05BDAE61.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/DA02DE1F-75FA-E511-8498-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/DA613127-75FA-E511-9D7B-24BE05C33C81.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/DA778125-75FA-E511-8DCC-0002C94CD116.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/E0EE922D-75FA-E511-B369-0002C94CD096.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/E6D0C91D-75FA-E511-A5F9-24BE05C6E711.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/F67D581D-75FA-E511-BCC7-5065F3818291.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/30000/F86E561D-75FA-E511-8124-5065F3818291.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/026042FC-2BFB-E511-8B1C-24BE05C6E711.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/04288705-37FB-E511-AE47-24BE05C6E7C1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/04A68854-13FB-E511-94AB-0002C94CDAF4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/0EFB5E74-13FB-E511-AC5E-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/1093BA31-1FFB-E511-AF74-B8CA3A709648.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/1600CA07-35FB-E511-A3CC-24BE05CEEB61.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/2E3D33BB-18FB-E511-8CB4-0002C94CDAF4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/3893BEA9-02FB-E511-A47C-24BE05C61601.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/3AB10F80-2EFB-E511-A4A7-24BE05CE6D61.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/4CC114F7-1CFB-E511-907E-B8CA3A709648.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/5491D307-3EFB-E511-90CE-24BE05C6E711.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/5E0E0636-1FFB-E511-9A2D-0002C94DBB18.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/68C9BB64-3EFB-E511-9642-24BE05C62611.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/6A5D166D-13FB-E511-A725-B083FED429D6.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/6ECC6DFD-08FB-E511-BC79-90B11C06CD59.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/761EC3F7-18FB-E511-8453-24BE05C33C61.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/78B20EF0-25FB-E511-B814-0002C94CDAF4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/7A5B83DF-25FB-E511-8D6D-0002C94CDAF4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/7CC4928E-1FFB-E511-ACB4-0002C94CD0BC.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/8291A44F-1FFB-E511-9E5A-24BE05C6C7E1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/887F0E44-1DFB-E511-9DAC-0002C94CD0BC.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/8AD0A2D6-1BFB-E511-8C8B-E839352D0E40.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/8AED3C80-2EFB-E511-8CE2-24BE05CEEB31.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/90C38742-3BFB-E511-B3DA-24BE05CEACA1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/98F61A3D-04FB-E511-9F8A-0002C94CD0D8.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/A0C8AFCF-2FFB-E511-9811-008CFAF07512.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/A23B6850-3BFB-E511-9511-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/A6C4DE56-13FB-E511-A9F5-0002C94CD120.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/A86A2B92-38FB-E511-B4AC-24BE05C48831.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/A8991C84-35FB-E511-9F84-24BE05CEDC81.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/B2695611-4EFB-E511-8E0B-5065F381E201.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/B2922230-35FB-E511-A5C2-5065F3818291.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/B64FC70F-4AFB-E511-8660-002590D9D8AE.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/B6600E44-1DFB-E511-AF57-0002C94CD0BC.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/B8AEC238-1FFB-E511-8D74-0002C94CD1E4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/BA9D6F6A-13FB-E511-913F-24BE05BD0F81.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/BC7AFB2A-1FFB-E511-961E-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/C0ECC2CB-4FFB-E511-A21C-0CC47A4D75EC.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/D0335135-1FFB-E511-ADBF-A0369F30FFD2.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/D23CBBDB-2BFB-E511-8BFD-24BE05C6D731.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/D802FE18-4EFB-E511-B82B-001E67A3F70E.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/D8D47F15-4EFB-E511-B27E-B083FED73AA1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/DC827BDC-3DFB-E511-8F43-24BE05CE1E51.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/DE3D0E44-1DFB-E511-B466-0002C94CD0BC.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/E072B93F-4EFB-E511-B08D-002590D9D8BC.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/E293A8D1-36FB-E511-9C1A-24BE05C6E711.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/E4854F37-45FB-E511-AD37-0002C94DBB20.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/E844A8D9-18FB-E511-B811-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/E8C9560B-2CFB-E511-A2C2-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/ECC69B4C-13FB-E511-9FFE-A0369F3016EC.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/FC000805-2CFB-E511-98D7-24BE05CECBD1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/40000/FCE1C4A9-02FB-E511-B5EB-24BE05C63791.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/0A2EA40C-02FB-E511-A68C-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/0C6A2007-02FB-E511-AE16-5065F3816251.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/18ACE318-02FB-E511-B24D-A0369F3102B6.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/1C9F605E-02FB-E511-8812-0025905B85BE.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/224C9008-02FB-E511-89EC-5065F381E271.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/26CD6D81-03FB-E511-B202-0002C94CDDEA.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/2EF32BD3-01FB-E511-96DA-008CFA197CE4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/368A421D-02FB-E511-97F0-0002C94CDDEA.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/3A10F107-02FB-E511-B340-24BE05C6C7F1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/3EA2D212-02FB-E511-A526-A0369F310374.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/4047AF07-02FB-E511-BD71-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/40CF6A64-02FB-E511-98D8-0CC47A4D769A.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/4687B407-02FB-E511-BDCC-5065F3818271.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/605AC919-02FB-E511-A183-0002C94CD1E4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/60B94809-02FB-E511-8629-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/66640E4B-03FB-E511-8B6F-0002C94D5504.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/70D426A1-01FB-E511-ADC1-0CC47A0AD6E0.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/7669978E-03FB-E511-A9C7-0002C94CD0D8.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/7A8B4606-02FB-E511-B4AC-5065F381B2D1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/7EF5DEA0-01FB-E511-AB64-000F530E4774.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/82254510-02FB-E511-84F2-24BE05C6E7E1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/862CF707-02FB-E511-B4E6-24BE05C6E7C1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/88DD101D-02FB-E511-9D69-0002C94CD0D8.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/A224CD94-03FB-E511-B533-0002C94CDAF6.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/A6B4FC09-02FB-E511-B910-24BE05CEEB61.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/A8BFCF06-02FB-E511-A9D6-24BE05C44BB1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/AC07F134-02FB-E511-BA74-0002C94CD12E.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/BA1571A8-01FB-E511-905A-001517E741C8.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/BE8FBC06-02FB-E511-BB1C-24BE05C33C81.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/C490FE75-03FB-E511-BA4D-A0000420FE80.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/C4FB5C0D-02FB-E511-9176-0002C94CDAF6.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/CC24BA55-03FB-E511-AEF7-0002C94CD1E4.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/D4B0FD06-02FB-E511-95B8-5065F3815221.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/DEF9681A-03FB-E511-B6C2-5065F381F1C1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/E62BF707-02FB-E511-859D-24BE05C6E7C1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/E63E8507-02FB-E511-9C0C-5065F38122A1.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/F0026937-02FB-E511-A4C9-0002C94D5504.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/F277B3C6-01FB-E511-8870-24BE05CEBD61.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/F62C6069-03FB-E511-82C5-24BE05C6D711.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/F6AD7795-02FB-E511-86E8-842B2B1815B3.root',
       '/store/mc/RunIISpring16MiniAODv1/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/50000/FEAC7864-03FB-E511-A308-0002C94CD12E.root',
] )
