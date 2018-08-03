import os
import commands

select_datasets = ["/DYJetsToLL_M-50_HT-*/*Apr2018*/MINIAOD*", "/WJetsToLNu_HT-*/*Apr2018*/MINIAOD*", "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"]

for select_dataset in select_datasets:

    datasets = commands.getstatusoutput('dasgoclient --query="dataset=%s"' % select_dataset)[1].split("\n")

    for dataset in datasets:

        print "Writing", dataset

        shortname = dataset.split('/')[1]
        if "ext1" in dataset:
            shortname += "_ext1"
        print shortname

        files = commands.getstatusoutput('dasgoclient --query="file dataset=%s"' % dataset)[1].split("\n")
        print files

        preamble = """import FWCore.ParameterSet.Config as cms

    maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
    readFiles = cms.untracked.vstring()
    secFiles = cms.untracked.vstring()
    source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
    readFiles.extend( [
    """

        ending = """] )"""

        with open(shortname + "_cff.py", 'w+') as f:

            f.write(preamble)
            for ifile in files:
                f.write("'" + ifile + "',\n")
            f.write(ending)
        
