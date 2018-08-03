import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/DA4C7457-AF41-E811-ADE7-0CC47A4C8EE2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/D41C0F47-BF41-E811-93E5-0CC47A4C8ECE.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6417B702-C741-E811-B3DA-0CC47A7C35C8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/587CBB04-C741-E811-BA20-0CC47A78A418.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/E6EE3675-C641-E811-A369-0CC47A4D75F2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/48FF7474-C641-E811-B01A-0CC47A4D76AA.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/AC5B5079-C641-E811-92FA-0CC47A7C35C8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/8EA85187-BF41-E811-BD4D-0025905A6090.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/54DFAF86-BF41-E811-9A89-0025905A60BE.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/D6936EE5-C541-E811-B571-0CC47A78A426.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F07DA2E4-C541-E811-9F88-0CC47A78A4A6.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/869D21E7-C541-E811-B393-0CC47A78A436.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/160C89E5-C541-E811-8BC4-0CC47A4C8EE2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/3E1EA8E6-C541-E811-959F-0CC47A7C3638.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/0C586B7A-C641-E811-8B69-0CC47A78A4B0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/9EA0C978-C641-E811-A31E-0CC47A4D7632.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/D8074BED-C541-E811-AE3A-0025905A6056.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/C6E630EF-C541-E811-9003-0025905A497A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/60C01F24-CA41-E811-9F44-0CC47A4C8EE2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/486C3025-CA41-E811-A2DD-0CC47A7C35D2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/86745325-CA41-E811-8001-0CC47A4D768E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/1055A62B-CA41-E811-BFE0-0025905B857A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/D429152E-CA41-E811-98FE-003048FFD72C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B403302C-CA41-E811-992F-0025905B85DC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/32216A04-CC41-E811-9708-0CC47A4C8E3C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/A291780B-CC41-E811-9300-0025905A60A0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/3887CE0B-CC41-E811-AE0D-0025905A605E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6825F845-CC41-E811-B67A-0CC47A4D7634.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/96BC1347-CC41-E811-BE42-0CC47A4C8E7E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/86974B44-CC41-E811-8435-0CC47A74524E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/7AADCD4E-CC41-E811-91A9-0025905A6090.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/8A49DE46-CC41-E811-8411-0CC47A7C351E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/ACC98246-CC41-E811-83B5-0CC47A4D767E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/1E762945-CC41-E811-B86A-0CC47A7C3412.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/FC29DDA7-CE41-E811-86FE-0CC47A7C357A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/62F167A7-CB41-E811-BA7F-0025905A60DA.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F887CEA0-CB41-E811-94DD-0CC47A4D7634.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/CEB75FA5-CB41-E811-ACAD-0025905A60A6.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/5E0DB5A6-CB41-E811-BE20-0025905A6056.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/AE8B209F-CB41-E811-A395-0CC47A7C357E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/E81802A1-CB41-E811-AB54-0CC47A4D7640.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/00FF479E-CB41-E811-A04B-0CC47A4C8E28.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/EEE59E9F-CB41-E811-92E6-0CC47A4D76B2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/DCBF58AC-D341-E811-A110-0CC47A7C35F8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/A07CB0C4-D341-E811-9A45-0025905A60FE.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6A28C89A-CF41-E811-980F-0CC47A7C3412.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/46DED7E2-D141-E811-9F2D-0CC47A7C357A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/AC9CEF50-D241-E811-B4D2-0CC47A4D769C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B21ED8A7-CE41-E811-91D1-0CC47A4D76D0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/7AF289AF-CE41-E811-85AF-0025905B85F6.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/5CF58DA8-CE41-E811-BFB8-0CC47A4D769C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/A22E89D3-CE41-E811-B6EB-0CC47A4D76D2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/6C1364A9-CE41-E811-8951-0CC47A4D7604.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/96FA538A-CE41-E811-959A-0CC47A4C8E98.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/1E031CAA-CE41-E811-AD44-0CC47A4D760A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/A6B7A46F-D041-E811-B326-0025905A6090.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/EC15FAA4-CE41-E811-B047-0CC47A7C35D2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/54F2E4A7-CE41-E811-A185-0CC47A7C3612.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/3E325AAE-CE41-E811-B88C-0025905B85C6.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/CC94B259-DA41-E811-B59C-0CC47A7C35C8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F465725A-DA41-E811-9353-0CC47A7C3610.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/56C45B5B-DA41-E811-AF7C-0CC47A4D7654.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/54A68A5B-DA41-E811-B9E4-0CC47A78A3EC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/E07EEB59-DA41-E811-B280-0CC47A4C8E3C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/2CB7D35B-DA41-E811-B1D5-0CC47A4C8EA8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/62665059-DA41-E811-BF02-0CC47A4C8ECE.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/CA35905F-DA41-E811-80AA-0025905B85BA.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/68B64E71-DC41-E811-9681-0CC47A7C34B0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/18D3662A-DF41-E811-9030-0CC47A7C34B0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/D0D066A9-E141-E811-99F5-0CC47A4D769C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/F490BCE5-C541-E811-8938-0CC47A4D7662.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/CEA3FF8E-CE41-E811-99C6-0025905B857A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/CC04222D-FF41-E811-AA9C-0CC47A4D768C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/36893732-0F42-E811-9CA4-0CC47A4C8EEA.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/56BC8DC4-1E42-E811-8122-0CC47A4C8EB0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/FCC44E1F-CD41-E811-A561-001E67F8F727.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/EC5591B0-7242-E811-91F3-001E6739C801.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/66D4D470-B841-E811-98F5-0CC47A7C3410.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/0499E42C-BC41-E811-957D-0025905A6088.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/920957D5-C041-E811-A05E-0CC47A4D7604.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7418369A-C741-E811-8A08-0CC47A7C347E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7C51A4D6-C741-E811-B1AE-0025905A607A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F437E5BC-C941-E811-9D36-0CC47A4D76C8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/3CD3349A-C841-E811-A733-0CC47A78A2EC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/B6BCB944-CC41-E811-B41D-0CC47A78A3EC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/C84D88D9-CA41-E811-956A-0025905B857E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/64F7264E-CC41-E811-9A62-0025905B857E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/06741C42-CE41-E811-9A60-0CC47A74525A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/CAD51A0F-CE41-E811-B711-0CC47A4C8ECE.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/A4CBAD10-CE41-E811-B4B3-0CC47A78A3EC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/24811E42-CE41-E811-B026-0CC47A74525A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/98287E43-CF41-E811-A876-0CC47A4D767E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/CEB4A4A9-D041-E811-A63A-0025905A6056.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/70893AA2-D041-E811-811B-0CC47A4C8E70.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/92B6B753-D041-E811-AEC0-003048FFCBB2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/ACDE3790-D141-E811-BFD4-0025905B85D8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/4A6AF58F-D141-E811-809D-0025905B85EC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/4C07CA8F-D141-E811-B8B3-0025905A60B2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D46D3A5D-D241-E811-9DDC-0CC47A4D767E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/64E5D02B-D441-E811-A48A-0CC47A4D7634.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/20ABDCB3-D241-E811-BAC1-0CC47A78A468.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/D08DE784-D241-E811-9527-0025905A60E4.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/AE0D70CA-D341-E811-B96D-0CC47A4C8E70.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/8E3480B9-D241-E811-B556-0025905A48FC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/02E9A9FC-D341-E811-9519-0025905B856C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F43C8D04-D541-E811-911A-0CC47A7C3610.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/28C94C2B-D441-E811-8DA0-0CC47A78A4A0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/AA5E24FE-D541-E811-99C8-0CC47A7C34A0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/4C034DFD-D541-E811-BF40-0CC47A4C8EE2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/424A606F-D741-E811-835B-0CC47A4D75F2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/86F8A772-D641-E811-9B25-0CC47A7C34EE.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/EAF20966-D741-E811-A547-0CC47A7C3610.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/EE75E804-D841-E811-B245-0025905B8612.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F03EA793-D841-E811-A272-0CC47A4C8EA8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/2CF3365D-DA41-E811-A03D-0CC47A7C3408.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/3AA27E62-DA41-E811-B06C-0025905B85D8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/40FA85A8-D941-E811-BCD8-0CC47A78A426.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/264081E5-DE41-E811-962C-0CC47A4D769C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/F22A40E3-DE41-E811-B9D4-0CC47A4D760A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/848205EE-DE41-E811-8FC6-0025905B8612.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/323DFAE5-DE41-E811-92DC-0CC47A78A4A0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/623C2DE3-DE41-E811-888D-0CC47A4D75EC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/048666E1-DE41-E811-951D-0CC47A7C3444.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/82D651EB-DE41-E811-822A-0025905A60E4.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E2281133-DF41-E811-95D9-0CC47A4D769C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/3E328A31-DF41-E811-A889-0CC47A4C8E34.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/3C73C831-DF41-E811-924B-0CC47A4D762E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/5037D83A-DF41-E811-989D-0025905B858A.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/9C973D6A-E441-E811-8F7F-0CC47A4C8F26.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/7C90CF72-E441-E811-A97D-0025905A60B2.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/4C34586C-E441-E811-AAF5-0CC47A78A414.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/5416ED6A-E441-E811-BA1A-0CC47A4D762E.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/E05AAA73-E441-E811-AB65-0025905B85B6.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/FE4C9D6B-E441-E811-8D53-0CC47A4D75EC.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/5466396A-E441-E811-BAAB-0CC47A4D75EE.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/226706E4-EC41-E811-8970-0CC47A78A3D8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/82D0CEE3-EC41-E811-B277-0CC47A78A3D8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/A81642EB-EC41-E811-88AD-0025905A6080.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/5A0CFA2C-F141-E811-AD99-003048FFD71C.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/0AC0F328-1042-E811-8D6C-0CC47A7C34EE.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/FAD5AF59-2242-E811-8D3C-0025905A6134.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/68D880DE-0642-E811-A476-0CC47A4C8E86.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/6C267326-5C42-E811-83CE-0025905A60A0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/646AE051-7342-E811-A20F-0025905A48F0.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/5AD9E4A6-CE41-E811-8FD0-001E6739A751.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/ECB75276-DA41-E811-B000-001E67FA38A8.root',
'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/540C4C26-1F42-E811-98CE-001E67FA402D.root',
] )