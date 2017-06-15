import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/0013D77B-C4EC-E611-8A2E-00259073E510.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/0052E13F-0CED-E611-ADEE-0090FAA583C4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/006AE655-CAEC-E611-96B8-0025907B4FAE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/02B14398-1DED-E611-A5A2-0090FAA57F44.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/02C9DF55-CAEC-E611-B023-00259073E488.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/06182558-FDEC-E611-ADA7-0090FAA58B94.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/0630FFF1-52ED-E611-ACC1-0090FAA578F0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/08317318-47ED-E611-BA28-0090FAA581E4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/08684A45-0CED-E611-A80F-00259073E44E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/0AC6A27A-C4EC-E611-A398-0090FAA576C0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/0C57CA19-47ED-E611-B7EE-00259073E322.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/0E0433AD-F6EC-E611-82BF-0090FAA57470.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/0E11DB46-0CED-E611-8928-00259073E32A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/0E45071B-05ED-E611-B8CE-0090FAA57430.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/1035561B-05ED-E611-9600-0090FAA57ED4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/1068034D-1BED-E611-91AB-0090FAA58B94.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/10BB6C40-0CED-E611-9B50-20CF3027A5F5.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/126FC3EC-B7EC-E611-A540-00259073E488.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/1294B53F-0CED-E611-9A34-0090FAA57360.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/141658CB-13ED-E611-AB9D-0090FAA575B0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/141F1815-47ED-E611-A254-0090FAA59124.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/149BCA44-0CED-E611-939E-00259073E4CE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/181A4199-1DED-E611-8180-0090FAA58C74.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/18597442-0CED-E611-9D08-0090FAA57CD4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/1A0D8994-B8EC-E611-8766-00259073E44E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/1C873154-CAEC-E611-ACBC-20CF3027A598.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/1CF84CD9-D1EC-E611-91AD-0090FAA58864.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/22518D5A-CAEC-E611-B9A8-00259073E49A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/22571751-1BED-E611-9873-0090FAA59634.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/22B67B60-1BED-E611-ABA0-00259074AEAE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/2411364E-1BED-E611-BE6B-0090FAA58BF4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/245B93F6-52ED-E611-A4D5-00259073E51A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/28724081-50ED-E611-AC9E-00259074AE40.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/28884E91-B8EC-E611-9B9B-0090FAA57730.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/2CB5A0F3-BEEC-E611-8AFB-0090FAA576C0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/2CC9FC3F-0CED-E611-9499-0090FAA57F04.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/2E336806-BFEC-E611-B00E-00259073E35A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/2E6DF61A-05ED-E611-A2CB-0090FAA57730.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/320E0CCF-13ED-E611-AA4B-002590D0B0BA.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/32EA3ACE-13ED-E611-811B-0090FAA573B0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/34EF3CCF-13ED-E611-8375-0090FAA57E64.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/3842054E-1BED-E611-9CB7-0090FAA573E0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/3AFC195A-FDEC-E611-B362-00259073E53E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/3ED7D619-47ED-E611-A65F-0090FAA575C0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/3EE90317-47ED-E611-AE5A-20CF3019DF0C.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/4044C214-47ED-E611-A70E-0090FAA57BF0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/40AFAA60-FDEC-E611-B18F-00259073E3D2.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/40DA4F19-47ED-E611-8E0F-0090FAA57FA4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/40DA6B1C-05ED-E611-A7E1-002590D0AFD8.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/42511BCC-13ED-E611-8A67-0090FAA57A50.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/4665BAC2-B6EC-E611-BEA7-00259073E4E6.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/4815A14C-1BED-E611-83E7-0090FAA57470.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/48EA52CD-13ED-E611-A6E8-0090FAA57430.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/4A2272F1-BEEC-E611-A3E0-0090FAA583C4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/4A960551-1BED-E611-82AF-002590D0B086.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/4AA6D0CD-13ED-E611-BCBC-0090FAA57350.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/4C5A37D9-13ED-E611-82AC-00259073E384.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/4E888F4E-1BED-E611-A204-0090FAA58BF4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/521ADD91-B8EC-E611-B439-002590747E14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/5268C756-CAEC-E611-8D23-0025907277A0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/529EDA57-1BED-E611-9C6D-0090FAA57D64.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/569CDD92-B8EC-E611-82DD-20CF3027A598.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/58294D1A-05ED-E611-8CE9-0090FAA57A50.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/5C91AAD3-04ED-E611-9CD1-0090FAA59114.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/5CA9A2CD-13ED-E611-A3B6-0025907B4FAE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/6001B43F-0CED-E611-90D3-0090FAA57F34.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/6074B9F1-BEEC-E611-A3F6-0090FAA57960.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/626FA8F3-52ED-E611-BDB9-0090FAA58434.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/64F6087B-C4EC-E611-BE6A-00259073E4E2.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/6639CB3F-0CED-E611-9B14-0090FAA1AD04.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/66785D7F-C4EC-E611-9B26-00259073E49A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/68C80D3B-51ED-E611-BF9B-00259073E3FE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/6A010ACE-13ED-E611-8A64-0025907B4FC4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/6A41B51E-05ED-E611-A48A-20CF3019DF17.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/6C68A63E-0CED-E611-9D12-0090FAA57620.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/6CB084F2-BEEC-E611-A176-0090FAA569C4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/6CF63D15-47ED-E611-9B47-20CF305B0551.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/6ECCA453-CAEC-E611-B4C2-0090FAA583C4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/705604F6-BEEC-E611-9B29-00259073E3E0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/72BA2F44-0CED-E611-BE86-00259073E4BC.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/74A84AF7-52ED-E611-A10E-0090FAA57630.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/762BBA13-47ED-E611-83E3-0090FAA58C74.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/78B7283A-51ED-E611-A99C-0090FAA57F34.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/7A20A9CD-13ED-E611-BB02-0090FAA57630.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/7A910C5B-FDEC-E611-95BD-002590D0B0C8.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/7C45CB14-47ED-E611-96B6-0090FAA57C60.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/7CD56543-0CED-E611-9293-00259073E4F6.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/7EAD349D-1DED-E611-A5A3-0090FAA57630.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/8043EE17-47ED-E611-BBDE-002590747E14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/80A30857-CAEC-E611-ADFE-00259074AE7A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/82650ED0-13ED-E611-9967-00259073E452.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/826E187B-C4EC-E611-B354-00259074AE40.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/842995CE-13ED-E611-840D-0090FAA584D4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/843173CF-13ED-E611-8AA8-0090FAA582E4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/84E73014-47ED-E611-8DBA-0090FAA57B10.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/8682FBF3-BEEC-E611-87C5-00259073E4CE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/86A10619-47ED-E611-8674-00259073E4BC.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/88115121-05ED-E611-82A5-00259073E4CA.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/88497A79-C4EC-E611-8FDF-00259073E466.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/884C5C57-CAEC-E611-8A34-0025907277A0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/888D0440-0CED-E611-8AD6-0090FAA57FA4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/88907841-0CED-E611-A534-0090FAA58D04.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/88B3E83F-0CED-E611-A2E5-0090FAA585D4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/8C8B493B-51ED-E611-BF97-0025907277BE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/8CA74339-51ED-E611-9D04-0090FAA58754.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/8CEE8F40-0CED-E611-B53F-0090FAA57260.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/8EAD5E4C-1BED-E611-9F13-0090FAA587C4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/904CFF50-1BED-E611-9B63-00259073E45E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/925F8214-47ED-E611-8297-0090FAA57CE4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/92814A4C-1BED-E611-9F82-0090FAA59114.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/96BAE44C-1BED-E611-A8B8-0090FAA56F60.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/98A9D0D6-04ED-E611-9593-0090FAA59ED4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/98D18016-47ED-E611-8DFE-00259073E460.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/9A238B58-CAEC-E611-8DC9-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/9A84F243-0CED-E611-97C5-0090FAA581E4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/9AB5F9D1-13ED-E611-9FEF-0090FAA58194.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/9C31353E-0CED-E611-911C-0090FAA579F0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/9C4C6EF1-BEEC-E611-8D74-0090FAA58254.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/9CC55444-0CED-E611-B5CF-00259073E466.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/9CD7357B-C4EC-E611-942C-00259073E370.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/9E077A05-05ED-E611-8387-0090FAA58124.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/9EA778D0-13ED-E611-8FB6-00259073E32A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/9EB9B011-47ED-E611-B12C-0090FAA58294.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/9EE67F41-0CED-E611-9889-0090FAA58434.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A0C8D152-CAEC-E611-9F0C-002590D0B086.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A24B683C-51ED-E611-B714-00259074AE54.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A281345D-CAEC-E611-93D8-0090FAA58134.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A2E9BBF3-52ED-E611-AD01-0090FAA57FA4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A2F50A4E-1BED-E611-B30F-0090FAA57FA4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A428264E-1BED-E611-9799-0090FAA58484.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A4A30459-FDEC-E611-ACF4-0090FAA569C4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A4B43B3A-51ED-E611-8478-00259074AE40.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A4D0C11A-47ED-E611-8B61-00259073E45E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A69422CD-13ED-E611-8784-0090FAA58BF4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A81B05DA-13ED-E611-B023-0090FAA576A0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A8402150-1BED-E611-A92A-0090FAA579F0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A871AFB8-52ED-E611-9F60-002590D0B086.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/A8A4E113-47ED-E611-9452-0090FAA583F4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/AA5062F0-BEEC-E611-8A39-0090FAA57730.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/AAAF3F77-C4EC-E611-9D2F-0090FAA58864.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/ACF4D115-47ED-E611-BC53-00259073E4CA.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/AEF9641D-05ED-E611-9C4E-20CF305B0591.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/B05FC651-1BED-E611-88B7-0090FAA57710.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/B093B858-FDEC-E611-98BE-002590D0AFD8.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/B0D5A5F4-BEEC-E611-B60D-00259073E534.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/B0F242C0-B6EC-E611-A442-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/B203473F-0CED-E611-9F18-002590D0AF54.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/B23668D8-B5EC-E611-94B3-00259073E516.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/B4137450-1BED-E611-B497-0090FAA576C0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/B4B15338-51ED-E611-998D-0090FAA57770.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/B4B57BCC-13ED-E611-87DE-0090FAA56994.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/B6CD515C-FDEC-E611-A4B3-00259073E4CE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/B80ED450-1BED-E611-BA74-00259074AE40.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/BC55BC1A-05ED-E611-B232-0090FAA57780.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/BC98F655-51ED-E611-B658-0090FAA57F74.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/BE122114-47ED-E611-8336-0090FAA57CD4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/BE341CD0-13ED-E611-9498-00259073E398.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/BE7BA95C-D1EC-E611-B6C3-00259073E49A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/BE892240-0CED-E611-AC34-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/C0678D59-FDEC-E611-8BAC-0025907B5048.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/C436A8F0-BEEC-E611-B5E9-0090FAA58D84.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/C4B53502-BFEC-E611-83C8-00259074AE80.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/C61EEFF3-52ED-E611-B695-0090FAA57A60.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/C62B5E7A-C4EC-E611-B81B-00259073E50A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/C635154E-1BED-E611-9C6B-0090FAA58104.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/C6671FF3-52ED-E611-9B7B-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/C85239F4-BEEC-E611-8FFD-00259073E4EA.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/CA3AA416-47ED-E611-8AD3-0025907B4FC4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/CE0CCDF5-BEEC-E611-ABB3-00259074AED2.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/CED3FE4F-1BED-E611-89A2-0090FAA57910.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/D204D058-1BED-E611-9AB2-0090FAA57400.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/D2C7FD13-47ED-E611-A72E-0090FAA57A60.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/D61FAC78-C4EC-E611-B761-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/D89B46F4-52ED-E611-AFD4-0090FAA597B4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/D8BF8357-1BED-E611-8A11-0090FAA57D64.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/DA42EB52-CAEC-E611-9170-0090FAA57BF0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/DA51655E-D1EC-E611-AEC3-0090FAA57310.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/DABF73F4-BEEC-E611-82E5-0025907277A0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/DADEC4CF-13ED-E611-9946-0090FAA58864.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/DE55BA43-0CED-E611-AD74-20CF3027A62F.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/E0703615-47ED-E611-A186-0090FAA58484.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/E08B37F2-52ED-E611-99B6-0090FAA58824.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/E0E716F1-BEEC-E611-BD2C-0090FAA1ACF4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/E26FA291-B8EC-E611-BE97-00259073E4D4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/E26FC042-0CED-E611-A537-002590D0B094.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/E2856357-FDEC-E611-8C53-0090FAA58D84.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/E2956040-0CED-E611-9F2D-0090FAA58C74.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/E43B397B-C4EC-E611-A933-00259074AE40.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/E6839451-1BED-E611-8B3E-20CF3019DF0C.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/E87BA14F-1BED-E611-B886-0090FAA57430.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/EA367877-C4EC-E611-B051-0090FAA1AD04.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/EA4D4888-B8EC-E611-8E0D-00259073E51C.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/EAEA2042-51ED-E611-A2A1-002590D0B086.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/EC360877-C4EC-E611-BA36-0090FAA58C74.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/EC585313-47ED-E611-B19C-0090FAA58C74.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/EC606EF6-BEEC-E611-AA80-00259073E466.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/F4343814-05ED-E611-9F01-0090FAA57F34.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/F49B1B19-47ED-E611-828F-002590D0B0C8.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/F6E18FF6-52ED-E611-B8E3-00259073E452.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/F837CC4E-1BED-E611-A94A-0090FAA57BE0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/F86F4152-CAEC-E611-826C-0090FAA57470.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/FA6DDAF1-BEEC-E611-8FDE-0090FAA58C74.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/100000/FCF5CC1A-05ED-E611-9EAC-0090FAA573B0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/00C6085A-20EC-E611-BF3E-00259073E3FC.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/02324756-27EC-E611-BC8F-00259073E3FC.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/023F2E4D-1AEB-E611-B31C-00259073E4CE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/024F8B45-27EC-E611-8DF6-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/06A47F66-04EC-E611-9775-0090FAA58BF4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/06C74145-27EC-E611-9195-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/06E9D7A9-19EC-E611-A930-0090FAA57E34.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/0C192C5A-20EC-E611-BC61-00259073E3FC.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/0C92EAA0-2EEC-E611-897A-0090FAA597B4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/0CB28C56-20EC-E611-8204-0090FAA569C4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/106A677F-0CEB-E611-B0A1-00259073E524.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/10C93FBC-51EB-E611-8CF9-0090FAA57560.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/129EA262-71EB-E611-8638-0090FAA57ED4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/14CFBD5A-20EC-E611-BC64-00259073E3C8.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/1664A257-20EC-E611-8439-0090FAA57410.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/16AFCC55-20EC-E611-BBD4-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/18402A54-05EB-E611-89AB-0090FAA57380.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/1877820F-28EC-E611-98AF-00259073E3C8.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/1ABB41AD-2FEC-E611-85F2-20CF30561726.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/1C665660-D8EC-E611-9363-0090FAA57B10.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/20D4CA95-30EC-E611-9425-0090FAA578F0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/2209D50C-28EC-E611-A3D6-00259073E3FE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/2287335B-20EC-E611-9862-00259073E3FE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/22CB2640-36EC-E611-9A09-0090FAA57F44.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/24F13909-28EC-E611-8CF3-0090FAA58754.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/265EAF02-14EB-E611-B908-0090FAA57360.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/283C305D-20EC-E611-BE37-0090FAA57F74.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/28EE0AEB-0AEC-E611-931E-002590D0B0C8.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/2A2E477B-0CEB-E611-B291-0CC47A4DEED6.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/2ECBF59C-30EC-E611-97CB-00259074AE8C.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/30784896-30EC-E611-961A-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/38142438-1AEB-E611-8DAB-0090FAA57420.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/386A01AA-19EC-E611-8C9C-0090FAA57E34.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/386CC258-20EC-E611-9182-002590D0AFEC.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/38DCE0D4-FCEB-E611-A384-00259073E4CA.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/3A67D744-27EC-E611-B04C-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/3C9757B8-0CEB-E611-A85F-00259073E4CA.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/3E40B066-04EC-E611-9491-0090FAA572B0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/40778445-27EC-E611-96AE-002590D0AFF2.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/40ACE19A-30EC-E611-90AC-00259073E32A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/443AEB44-27EC-E611-AB7E-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/480403A4-19EC-E611-B03B-00259073E38A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/4A14409D-35EC-E611-889C-00259073E488.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/4AB7240A-28EC-E611-B2BC-0090FAA58BF4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/4C31DC0D-28EC-E611-B2FA-00259073E38A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/509E8870-20EC-E611-B38C-485B3989725C.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/5481F966-04EC-E611-8A28-0090FAA59634.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/560CB75A-20EC-E611-A1E2-00259073E4BC.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/5629035B-20EC-E611-8EA3-00259073E4BC.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/565FECB0-19EC-E611-8D2B-002590D0B0C8.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/58B2C39A-30EC-E611-B95D-485B3989725C.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/58B70003-14EB-E611-A272-0090FAA583C4.root',
] )
readFiles.extend( [
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/5CE5CDCE-FCEB-E611-9E99-002590D0B0B6.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/6055D0A9-19EC-E611-B888-0090FAA57E34.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/607D0652-05EB-E611-B901-0090FAA578E0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/60DE9158-20EC-E611-8653-002590D0AFF2.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/62821E9C-27EC-E611-BFFB-00259022277E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/64E9D45B-20EC-E611-82B2-00259073E50A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/668F7658-20EC-E611-8CFE-002590D0AFF2.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/669E3754-05EB-E611-A700-0090FAA583F4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/6C6E0A06-14EB-E611-8BB2-002590D0AF54.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/6C9A127A-0CEB-E611-B76A-0090FAA578E0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/72EE72D7-11EC-E611-B27E-0090FAA581E4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/76039377-0CEB-E611-A656-00259073E35A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/78061D56-20EC-E611-BE60-0090FAA569C4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/7868FEA5-19EC-E611-96B5-0025907277A0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/7CA3D045-36EC-E611-9B65-0090FAA57A00.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/80C114DC-11EC-E611-936C-00259073E398.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/8292F044-27EC-E611-99C0-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/82A88805-14EB-E611-B437-20CF3019DEED.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/868330A5-19EC-E611-8401-00259022277E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/86A229A6-19EC-E611-8B9B-0025907277A0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/8ACCC364-2EEC-E611-9F33-00259022277E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/8C13FECC-FCEB-E611-8AF4-20CF3019DF17.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/8E2DFA0D-35EC-E611-A620-00259022277E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/8EC5B782-0CEB-E611-9DAE-00259074AE38.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/9092F881-0CEB-E611-ADAB-00259073E35A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/90DA140B-28EC-E611-9E5E-00259073E4E6.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/92479747-27EC-E611-B7AD-002590D0AFF2.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/9408E207-14EB-E611-A3B8-00259073E4E2.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/9695A65B-20EC-E611-8240-00259073E4E6.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/96E1150C-28EC-E611-A4EC-00259073E50A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/98DFEC20-3CEC-E611-BD9F-0090FAA57C00.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/9AA2AC5B-05EB-E611-B1A5-0090FAA57D64.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/9ABA1A36-35EC-E611-884F-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/9AD3A584-0CEB-E611-9E42-00259073E51C.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/9AFC5D55-05EB-E611-87F3-0090FAA58864.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/9AFDC5A9-19EC-E611-A00A-0090FAA57E34.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/9CF75C54-05EB-E611-9B7C-0090FAA57CE4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/9E7DA807-28EC-E611-8CEA-0090FAA573E0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/9E8CF4D1-34EC-E611-978F-00259073E38A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/A2604599-35EC-E611-A911-0090FAA57F34.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/A47F3C79-0CEB-E611-8C83-0090FAA58204.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/A601FBCC-FCEB-E611-8C3B-0090FAA56994.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/AA338937-1AEB-E611-8244-0090FAA57690.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/B054B346-27EC-E611-97BF-00259073E38A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/B6B227AC-2FEC-E611-91F7-0090FAA57E34.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/B8492F57-20EC-E611-A3BD-0090FAA57410.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/BC90373B-21EB-E611-A268-0090FAA576C0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/BEC3E258-05EB-E611-8C17-00259074AE3E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/BEEF1C5C-27EC-E611-8081-00259073E466.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/C01F4052-05EB-E611-97D4-0090FAA578E0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/C0B6C053-05EB-E611-8ED8-0090FAA58224.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/C250285C-20EC-E611-8369-00259073E4E6.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/C2CF2E43-27EC-E611-91D3-0090FAA58BF4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/C4F07DA7-19EC-E611-8B28-00259073E466.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/C61D7E81-0CEB-E611-97CE-00259073E470.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/C855FA2A-2EEC-E611-8D25-0025907277A0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/CCD1540B-28EC-E611-999C-20CF30561726.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/CCE8AD62-2EEC-E611-9D8C-00259073E4B6.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/CE6240DA-11EC-E611-A983-0CC47A4D9A10.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/CE983F0C-28EC-E611-9DBC-0025907B4F2C.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/CEC3955B-20EC-E611-8E22-00259073E3FE.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/D03BAC47-35EC-E611-A6C5-00259073E4B6.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/D272BA54-05EB-E611-9165-0090FAA57630.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/D4D12547-27EC-E611-95A0-002590D0AFF2.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/D6185FCC-FCEB-E611-BD29-0090FAA59114.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/D6832E3A-F5EB-E611-A6FB-00259073E38A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/D8B273A5-27EC-E611-BAA5-00259073E4B6.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/DA4AA864-35EC-E611-BA32-0025907277A0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/DA64A77F-0CEB-E611-A0C9-002590D0B02E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/DA9CD854-05EB-E611-B3B4-0090FAA57630.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/DC452083-0CEB-E611-9F11-00259074AE80.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/DC92312A-2EEC-E611-89F3-00259073E38A.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/DE6A1D5C-05EB-E611-9211-0090FAA57BE0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/DEBCAE5C-20EC-E611-8C6B-0090FAA57F74.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/E002E144-27EC-E611-B907-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/E0CDB65A-20EC-E611-9BBE-00259073E3C8.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/E2A1D347-27EC-E611-8C12-002590D0AFF2.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/E4995F7E-0CEB-E611-A271-0090FAA57AE0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/EC8116A6-19EC-E611-84F4-0025907277A0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/F0CA8470-20EC-E611-A010-485B3989725C.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/F0EA12A6-19EC-E611-9023-0025907277A0.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/F4A43044-27EC-E611-9963-0090FAA58BF4.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/F4D2D709-28EC-E611-905A-0090FAA57F14.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/F8269413-28EC-E611-AFC5-0090FAA57E34.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/F87030A5-19EC-E611-889E-00259022277E.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/F8CE4F7C-0CEB-E611-ACC1-0090FAA57410.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/F8DD3443-1AEB-E611-AF57-002590D0AF54.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/FC01C545-27EC-E611-A9C9-002590D0AFF2.root',
       '/store/data/Run2016B/HTMHT/MINIAOD/03Feb2017_ver2-v2/50000/FE76049A-04EB-E611-8035-0CC47A4DEF00.root',
] )