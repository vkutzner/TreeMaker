import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/0002DAEF-97EA-E611-8DEC-001E67E69DEC.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/00C46274-91EA-E611-9576-002590A3CA1A.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/0CAE644F-91EA-E611-8B61-001E6779285C.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/0EFD5460-8AEA-E611-ABA2-001E6779258C.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/1CFA8097-8AEA-E611-980F-001E67E6F855.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/2CB48ADC-ADEA-E611-A903-001E677926F8.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/344B3357-91EA-E611-84B2-002590A80DDA.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/3EFC0F51-91EA-E611-BB49-001E6779242E.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/42160E53-91EA-E611-A515-001E67E6F4B8.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/42792C54-91EA-E611-9F1D-001E67E6F4EF.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/4A6B1D50-91EA-E611-A16A-001E67398CA0.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/52BD2886-91EA-E611-8274-002590200AB8.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/5C72CD4F-91EA-E611-8EC3-001E6779271E.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/60F22FEF-97EA-E611-A68C-001E67E6F80F.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/62537052-91EA-E611-B24E-001E677923B6.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/74E78264-8AEA-E611-BBED-0025902008CC.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/80A5ED5A-91EA-E611-9591-002590A37114.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/84F7EF4E-9FEA-E611-A671-001E677925E6.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/8C02E450-91EA-E611-8590-001E67E63AE6.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/90ED68EF-97EA-E611-8979-001E67E6F4B8.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/923AFC4F-91EA-E611-B051-001E67792442.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/9A7CC85C-91EA-E611-B83C-002590A887FE.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/A0ABAF60-8AEA-E611-8CE3-001E67398CA0.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/AAE4C2AF-C1EA-E611-A439-001E67E6F4A9.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/BE1475F6-EDEA-E611-82A4-001E67E71A56.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/C2D956ED-97EA-E611-A2E3-002590A831B6.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/C4EB735E-91EA-E611-9EEB-001E67398C1E.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/C8E4E1EC-97EA-E611-BC0D-001E67792738.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/D2A47E61-8AEA-E611-A621-001E677925E6.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/F4D097A4-C1EA-E611-B99E-001E67398683.root',
       '/store/data/Run2016H/JetHT/MINIAOD/03Feb2017_ver3-v1/100000/F65DBFC6-93EA-E611-A4BE-001E673968BA.root',
] )
