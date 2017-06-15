import sys

# Read parameters
from TreeMaker.Utils.CommandLineParams import CommandLineParams
parameters = CommandLineParams()
scenarioName=parameters.value("scenario","")
inputFilesConfig=parameters.value("inputFilesConfig","")
dataset=parameters.value("dataset",[])##just try to print lots of things and understand what's happening before mucking. also, how do we use dbs to get filelists from multiple data tiers?

try: sidecar = dataset.replace('MINIAOD','AOD').replace('miniAOD','AOD').replace('MiniAOD','AOD')
except: sidecar = []
print 'dataset initially=', dataset
sidecar = ['/store/mc/RunIISpring16FSPremix/pMSSM_MCMC1_mH-120to130_batch1_TuneCUETP8M1_13TeV-pythia8/AODSIM/pLHE_80X_mcRun2_asymptotic_v12-v1/00000/007EA961-C39D-E611-88DE-02163E017618.root',
'/store/mc/RunIISpring16FSPremix/pMSSM_MCMC1_mH-120to130_batch1_TuneCUETP8M1_13TeV-pythia8/AODSIM/pLHE_80X_mcRun2_asymptotic_v12-v1/00000/4065A72F-C19D-E611-8349-FA163E342AFC.root',
'/store/mc/RunIISpring16FSPremix/pMSSM_MCMC1_mH-120to130_batch1_TuneCUETP8M1_13TeV-pythia8/AODSIM/pLHE_80X_mcRun2_asymptotic_v12-v1/00000/542DD08D-DB9D-E611-B962-FA163E070AA1.root',
'/store/mc/RunIISpring16FSPremix/pMSSM_MCMC1_mH-120to130_batch1_TuneCUETP8M1_13TeV-pythia8/AODSIM/pLHE_80X_mcRun2_asymptotic_v12-v1/00000/5AFF4F67-E19D-E611-B807-FA163E7707EC.root',
'/store/mc/RunIISpring16FSPremix/pMSSM_MCMC1_mH-120to130_batch1_TuneCUETP8M1_13TeV-pythia8/AODSIM/pLHE_80X_mcRun2_asymptotic_v12-v1/00000/702991AC-CF9D-E611-A3CB-02163E01765F.root',
'/store/mc/RunIISpring16FSPremix/pMSSM_MCMC1_mH-120to130_batch1_TuneCUETP8M1_13TeV-pythia8/AODSIM/pLHE_80X_mcRun2_asymptotic_v12-v1/00000/AE9E75CE-B69D-E611-8261-FA163E481135.root',
'/store/mc/RunIISpring16FSPremix/pMSSM_MCMC1_mH-120to130_batch1_TuneCUETP8M1_13TeV-pythia8/AODSIM/pLHE_80X_mcRun2_asymptotic_v12-v1/00000/C8C28777-9A9D-E611-81D7-FA163E6AAC05.root']
print 'sidecar initially=', sidecar

nstart = parameters.value("nstart",0)
nfiles = parameters.value("nfiles",-1)
numevents=parameters.value("numevents",-1)
reportfreq=parameters.value("reportfreq",1000)
outfile=parameters.value("outfile","test_run")
dump=parameters.value("dump",False)
mp=parameters.value("mp",False)

# background estimations on by default
lostlepton=parameters.value("lostlepton", True)
hadtau=parameters.value("hadtau", True)
hadtaurecluster=parameters.value("hadtaurecluster", 1)
doZinv=parameters.value("doZinv", True)

# compute the PDF weights
doPDFs=parameters.value("doPDFs", True);

# other options off by default
debugtracks=parameters.value("debugtracks", False)
applybaseline=parameters.value("applybaseline", False)
gridcontrol=parameters.value("gridcontrol", False)

# auto configuration for different scenarios
from TreeMaker.Production.scenarios import Scenario
scenario = Scenario(scenarioName)

# take command line input (w/ defaults from scenario if specified)
globaltag=parameters.value("globaltag",scenario.globaltag)
tagname=parameters.value("tagname",scenario.tagname)
geninfo=parameters.value("geninfo",scenario.geninfo)
pmssm=parameters.value("pmssm",scenario.pmssm)
fastsim=parameters.value("fastsim",scenario.fastsim)
signal=parameters.value("signal",scenario.signal)
jsonfile=parameters.value("jsonfile",scenario.jsonfile)
jecfile=parameters.value("jecfile",scenario.jecfile)
residual=parameters.value("residual",scenario.residual)
jerfile=parameters.value("jerfile",scenario.jerfile)
pufile=parameters.value("pufile",scenario.pufile)
era=parameters.value("era",scenario.era)

#temporary redirector fix
#fastsim signal is phedexed to LPC Tier3
redir=parameters.value("redir", "root://cmseos.fnal.gov/" if fastsim and signal else "root://cmsxrootd.fnal.gov/")
redir = parameters.value("redir","root://cmsxrootd.fnal.gov/")

# The process needs to be defined AFTER reading sys.argv,
# otherwise edmConfigHash fails
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process("RA2EventSelection")
if len(era)>0:
    process = cms.Process("RA2EventSelection",getattr(eras,era))

# configure geometry & conditions
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

# Load input files
readFiles = cms.untracked.vstring()
readFiles_sidecar = cms.untracked.vstring()

if inputFilesConfig!="" :
    print 'case inputFilesConfig!=""'
    if nfiles==-1:
        print 'sam checking case 1'
        process.load("TreeMaker.Production."+inputFilesConfig+"_cff")
        readFiles.extend( process.source.fileNames )
    else:
        print 'sam checking case 2'
        readFilesImport = getattr(__import__("TreeMaker.Production."+inputFilesConfig+"_cff",fromlist=["readFiles"]),"readFiles")
        readFiles.extend( readFilesImport[nstart:(nstart+nfiles)] )

if dataset!=[] :  
    readFiles.extend( [dataset] )

if sidecar!=[] :
    readFiles_sidecar.extend( sidecar )
    
for f,val in enumerate(readFiles):
    if readFiles[f][0:6]=="/store":
        readFiles[f] = redir+readFiles[f]

for f,val in enumerate(readFiles_sidecar):
    if readFiles_sidecar[f][0:6]=="/store":
        readFiles_sidecar[f] = redir+readFiles_sidecar[f]        
    
# print out settings
print "***** SETUP ************************************"
print " dataset: "+str(readFiles)
print " sidecar: "+str(readFiles_sidecar)
print " outfile: "+outfile+"_RA2AnalysisTree"
print " "
print " storing lostlepton variables: "+str(lostlepton)
print " storing hadtau variables: "+str(hadtau)+" w/ reclustering "+str(hadtaurecluster)
print " storing Zinv variables: "+str(doZinv)
print " "
print " storing PDF weights: "+str(doPDFs)
print " "
print " storing track debugging variables: "+str(debugtracks)
print " Applying baseline selection filter: "+str(applybaseline)
print " "
print " scenario: "+scenarioName
print " global tag: "+globaltag
print " Instance name of tag information: "+tagname
print " Including gen-level information: "+str(geninfo)
print " Including pMSSM-related information: "+str(pmssm)
print " Using fastsim settings: "+str(fastsim)
print " Running signal uncertainties: "+str(signal)
if len(jsonfile)>0: print " JSON file applied: "+jsonfile
if len(jecfile)>0: print " JECs applied: "+jecfile+(" (residuals)" if residual else "")
if len(jerfile)>0: print " JERs applied: "+jerfile
if len(pufile)>0: print " PU weights stored: "+pufile
print " era of this dataset: "+era
print "************************************************"

from TreeMaker.TreeMaker.makeTreeFromMiniAOD_cff import makeTreeFromMiniAOD
process = makeTreeFromMiniAOD(process,
    outfile=outfile+"_RA2AnalysisTree",
    reportfreq=reportfreq,
    dataset=readFiles,
    sidecar=readFiles_sidecar,
    globaltag=globaltag,
    numevents=numevents,
    hadtau=hadtau,
    hadtaurecluster=hadtaurecluster,
    lostlepton=lostlepton,
    applybaseline=applybaseline,
    doZinv=doZinv,
    debugtracks=debugtracks,
    geninfo=geninfo,
    pmssm=pmssm,
    tagname=tagname,
    jsonfile=jsonfile,
    jecfile=jecfile,
    residual=residual,
    jerfile=jerfile,
    pufile=pufile,
    doPDFs=doPDFs,
    fastsim=fastsim,
    signal=signal,
    scenario=scenarioName
)

# final tweaks to process
process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')
process.TFileService.closeFileFast = cms.untracked.bool(True)
if mp:
    process.IgProfService = cms.Service("IgProfService",
        reportEventInterval = cms.untracked.int32(1),
        reportFirstEvent = cms.untracked.int32(1),
        reportToFileAtPostEndJob = cms.untracked.string('| gzip -c > '+outfile+'___memory___%I_EndOfJob.gz'),
        reportToFileAtPostEvent = cms.untracked.string('| gzip -c > '+outfile+'___memory___%I.gz')
    )


# if requested, dump and exit
if dump:
    print process.dumpPython()
    sys.exit(0)
