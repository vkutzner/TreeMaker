import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/150/00000/EEF8F8DC-D719-E611-9F4B-02163E01421D.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/158/00000/00DD3222-261A-E611-9FD2-02163E011E34.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/158/00000/1A45407E-761A-E611-AD5E-02163E013724.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/158/00000/1AE468DD-1D1A-E611-8781-02163E014246.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/158/00000/227D5C09-211A-E611-BA9E-02163E01353C.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/158/00000/3C2DAC51-3C1A-E611-A74D-02163E013926.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/158/00000/3EF170D2-1C1A-E611-9962-02163E011F00.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/158/00000/56A3C0DA-221A-E611-B02F-02163E0137E0.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/158/00000/CA78978F-2A1A-E611-BEC8-02163E0134C5.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/158/00000/DC7DC170-761A-E611-AA0F-02163E0136FD.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/158/00000/E68DB076-761A-E611-94C1-02163E0121BD.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/276/00000/04849400-0E1A-E611-98D5-02163E014480.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/286/00000/8A59E031-121A-E611-B451-02163E01438C.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/287/00000/ECFDD632-141A-E611-BD1C-02163E014575.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/290/00000/6C00CF32-1C1A-E611-A851-02163E0142E0.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/291/00000/1C5D27DF-4B1A-E611-894A-02163E011DCB.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/292/00000/7EE41DEF-241A-E611-BD30-02163E0143D1.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/294/00000/B2B483DC-221A-E611-A19E-02163E013841.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/295/00000/7CECD1E3-251A-E611-BBB5-02163E0141A3.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/296/00000/246D1309-411A-E611-8AF4-02163E012830.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/299/00000/B01AB6C6-4C1A-E611-8F31-02163E011EEE.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/301/00000/1A526BCC-8B1A-E611-B7A6-02163E011DEB.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/301/00000/305C7A1D-931A-E611-A2FB-02163E0133C2.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/301/00000/4C02298D-8F1A-E611-A62E-02163E014170.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/301/00000/E0B43D98-9C1A-E611-985E-02163E01420B.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/302/00000/0428404B-BE1A-E611-B5BC-02163E0128E8.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/302/00000/64D5D682-A61A-E611-8DD8-02163E011C1D.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/302/00000/D040F688-9F1A-E611-8906-02163E01184A.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/394/00000/FA19D6F8-CE1A-E611-8B9C-02163E01367A.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/402/00000/6ADA63A6-491B-E611-8557-02163E013760.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/402/00000/C628C392-491B-E611-A663-02163E0139E1.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/403/00000/E8827251-4B1B-E611-BCC6-02163E011B63.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/404/00000/8C12ABD5-501B-E611-9F49-02163E011E12.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/405/00000/4E83313C-4F1B-E611-B3CF-02163E013584.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/406/00000/E4000AB2-651B-E611-9883-02163E0136EB.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/407/00000/54DA128A-3B1B-E611-9836-02163E012338.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/408/00000/74DCD7EB-471B-E611-99F2-02163E0138B2.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/409/00000/0049F859-C01B-E611-99CB-02163E01216A.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/409/00000/943D9D4F-C01B-E611-A018-02163E0129FA.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/409/00000/AE781D84-C01B-E611-A8D0-02163E0144B4.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/410/00000/8E7A2D51-7A1B-E611-A627-02163E011E2A.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/411/00000/E85796D6-811B-E611-89C8-02163E013705.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/412/00000/3ACAE4B5-761B-E611-83C0-02163E0141C1.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/424/00000/98E97099-571B-E611-8B0B-02163E014479.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/425/00000/02C32148-CD1B-E611-96F4-02163E011D20.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/425/00000/487C5ADF-C61B-E611-8AC7-02163E0146F0.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/425/00000/983BD739-D61B-E611-8521-02163E0142D3.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/425/00000/B0D0E4C6-C51B-E611-8990-02163E014183.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/425/00000/EE534358-D31B-E611-A751-02163E01262F.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/426/00000/CCD5D9A2-CE1B-E611-A041-02163E01473E.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/440/00000/3CE8D6C8-AB1B-E611-8422-02163E014795.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/442/00000/905CE546-AC1B-E611-BBFB-02163E0119CF.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/443/00000/3AE11098-B01B-E611-9984-02163E013637.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/444/00000/A47FEBCF-AE1B-E611-96E1-02163E011C18.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/445/00000/32B80269-B51B-E611-92AF-02163E01375B.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/446/00000/3C38E226-DE1B-E611-97F3-02163E01352D.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/447/00000/60E90944-1E1C-E611-A4C1-02163E011DAE.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/447/00000/62CF34E5-F61B-E611-84C2-02163E01390A.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/447/00000/A6E70E66-F91B-E611-BCA8-02163E012961.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/447/00000/BCC1A572-011C-E611-A5AE-02163E0146F0.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/448/00000/00E3FD59-2A1C-E611-93C0-02163E0145BB.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/448/00000/D41F0869-1E1C-E611-BAA2-02163E012065.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/448/00000/DEE87B17-1A1C-E611-867D-02163E0146B7.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/449/00000/8A79815C-381C-E611-962A-02163E014603.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/449/00000/9AE2C527-821C-E611-98AA-02163E014260.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/450/00000/2248308C-661C-E611-92D2-02163E014354.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/450/00000/5434B913-4C1C-E611-B8C2-02163E012042.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/450/00000/54FDC18E-541C-E611-9D08-02163E012209.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/450/00000/7C380726-511C-E611-B17B-02163E01208E.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/450/00000/B42806BF-4F1C-E611-8BDC-02163E014713.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/450/00000/D6F2708D-661C-E611-AFBB-02163E014761.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/492/00000/628014CC-7D1D-E611-83AF-02163E0133C9.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/492/00000/944F15D9-7D1D-E611-B3B4-02163E011F32.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/492/00000/EAA21BD4-7D1D-E611-AE26-02163E01469F.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/493/00000/E2FD3F7C-5E1D-E611-AE01-02163E011E3E.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/493/00000/F2028A72-5E1D-E611-ADEC-02163E014300.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/494/00000/B6276E21-6F1D-E611-B370-02163E0139B5.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/494/00000/CA5DE631-881D-E611-8678-02163E01254F.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/502/00000/386C8A09-251F-E611-9E51-02163E012A6E.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/502/00000/7491A5DD-241F-E611-9098-02163E014604.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/502/00000/78F52FE1-241F-E611-81CF-02163E012542.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/502/00000/96B1F9D8-241F-E611-9BD9-02163E0143BC.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/502/00000/AE2A90E8-241F-E611-A140-02163E013574.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/502/00000/B631E9C6-241F-E611-ABC4-02163E0118D8.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/502/00000/FA9A1CBE-241F-E611-B33E-02163E013486.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/503/00000/0ABBB108-511F-E611-A544-02163E012516.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/503/00000/46DCC811-511F-E611-AFF0-02163E012AF2.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/503/00000/4CB8360D-511F-E611-AAE7-02163E01381E.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/503/00000/E811966A-3F1F-E611-92CB-02163E0145E2.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/504/00000/EC2A7AE0-311F-E611-A988-02163E01340E.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/523/00000/F4A648F3-351F-E611-9DDC-02163E01430D.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/526/00000/E644B476-571F-E611-98BE-02163E013447.root',
       '/store/data/Run2016B/SinglePhoton/MINIAOD/PromptReco-v2/000/273/537/00000/7207BCBF-991F-E611-B4A9-02163E013918.root',
] )