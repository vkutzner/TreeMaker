#!/bin/env python
import os
from subprocess import check_output
import glob

################## configure here ##################
folderAOD = "root://cmsxrootd.fnal.gov//store/user/lpcsusyhad/sbein/LongLiveTheChi/aodsim/smallchunks"
outputFolder = "srm://dcache-se-cms.desy.de/pnfs/desy.de/cms/tier2/store/user/vkutzner/DisappTrksNtupleSidecar"
scenario = "Spring16Pmssm"
################## configure here ##################

filelistFile = open("filelistBM.txt",'r')
filelist = filelistFile.read().split()

completePaths = []
for fileMiniAOD in filelist:
    folderMiniAOD = folderAOD.replace("aodsim", "miniaodsim")
    completeMiniAODPath = folderMiniAOD  + "/" + fileMiniAOD
    fileAOD = fileMiniAOD.replace("step3_miniAODSIM", "step2_AODSIM")
    completeAODPath = folderAOD  + "/" + fileAOD
    completePaths.append([completeAODPath, completeMiniAODPath])
  
jdlTemplate = '''
universe = vanilla
Executable = jobExecCondorSingle.sh
+REQUIRED_OS = "rhel6"
request_disk = 10000000
request_memory = 2100
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT_OR_EVICT
Transfer_Input_Files = jobExecCondorSingle.sh, ../runMakeTreeFromMiniAOD_cfg.py, ../data, CMSSWVER.tar.gz
Output = JOBNAME_$(Cluster).stdout
Error = JOBNAME_$(Cluster).stderr
Log = JOBNAME_$(Cluster).condor
notification = Never
x509userproxy = $ENV(X509_USER_PROXY)
Arguments = CMSSWVER OUTDIR SAMPLE NPART NSTART NFILES SCENARIO JOBNAME
want_graceful_removal = true
EXTRASTUFF
on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
on_exit_hold = ( (ExitBySignal == True) || (ExitCode != 0) )
on_exit_hold_reason = strcat("Job held by ON_EXIT_HOLD due to ",\
        ifThenElse((ExitBySignal == True), "exit by signal", \
                                strcat("exit code ",ExitCode)), ".")
Queue 1
'''

shellscript = '''
#!/bin/bash

#
# variables from arguments string in jdl
#

echo "Starting job on " `date` #Only to display the starting of production date
echo "Running on " `uname -a` #Only to display the machine where the job is running
echo "System release " `cat /etc/redhat-release` #And the system release
echo "CMSSW on Condor"

# to get condor-chirp from CMSSW
PATH="/usr/libexec/condor:$PATH"
source /cvmfs/cms.cern.ch/cmsset_default.sh

CMSSWVER=$1
OUTDIR=$2
SAMPLE=$3
NPART=$4
NSTART=$5
NFILES=$6
SCENARIO=$7
JOBNAME=$8

echo ""
echo "parameter set:"
echo "CMSSWVER:   $CMSSWVER"
echo "OUTDIR:     $OUTDIR"
echo "SAMPLE:     $SAMPLE"
echo "NPART:      $NPART"
echo "NSTART:     $NSTART"
echo "NFILES:     $NFILES"
echo "SCENARIO:   $SCENARIO"
echo "JOBNAME:    $JOBNAME"

tar -xzf ${CMSSWVER}.tar.gz
cd ${CMSSWVER}
scram b ProjectRename
# cmsenv
eval `scramv1 runtime -sh`
cd -

# run CMSSW
ARGS="outfile=${JOBNAME} dataset=${SAMPLE} scenario=${SCENARIO}"
cmsRun runMakeTreeFromMiniAOD_cfg.py ${ARGS} privateSample=True 2>&1

CMSEXIT=$?

if [[ $CMSEXIT -ne 0 ]]; then
  rm *.root
  echo "exit code $CMSEXIT, skipping gfal-copy"
  exit $CMSEXIT
fi

# copy output to eos
echo "gfal-copy output for condor"
. /cvmfs/oasis.opensciencegrid.org/mis/osg-wn-client/3.3/current/el6-x86_64/setup.sh
for FILE in *RA2AnalysisTree.root
do
  echo "gfal-copy -f ${FILE} ${OUTDIR}/${FILE}"
  gfal-copy -f ${FILE} ${OUTDIR}/${FILE} 2>&1
  XRDEXIT=$?
  if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in gfal-copy"
    exit $XRDEXIT
  fi
  rm ${FILE}
done
'''

jdlTemplate = jdlTemplate.replace("$ENV(X509_USER_PROXY)", check_output(["sh", "-c", "voms-proxy-info | grep path | cut -b 13-"]))
jdlTemplate = jdlTemplate.replace("CMSSWVER", os.environ['CMSSW_BASE'].split("/")[-1])
jdlTemplate = jdlTemplate.replace("OUTDIR", outputFolder)
jdlTemplate = jdlTemplate.replace("NPART", "0")
jdlTemplate = jdlTemplate.replace("NSTART", "0")
jdlTemplate = jdlTemplate.replace("NFILES", "1")
jdlTemplate = jdlTemplate.replace("SCENARIO", scenario)
jdlTemplate = jdlTemplate.replace("EXTRASTUFF", "")


for item in completePaths:

    miniaodpath = item[1]
    jobname = miniaodpath.split('/')[-1].split('.root')[0]

    jdlFile = jdlTemplate
    jdlFile = jdlFile.replace("JOBNAME", jobname)
    jdlFile = jdlFile.replace("SAMPLE", miniaodpath)
    
    fout = open("jobExecCondor_%s.jdl" % jobname, "w")
    fout.write(jdlFile)
    fout.close()


# write out shell script
fout = open("jobExecCondorSingle.sh", "w")
fout.write(shellscript)
fout.close()
os.system("chmod +x jobExecCondorSingle.sh")

raw_input("Ready to submit jobs, press <return>")

for jdlFile in glob.glob("jobExecCondor_pMSSM12*jdl"):
    os.system("condor_submit %s" % jdlFile)
