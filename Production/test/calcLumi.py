import os,glob,re
from optparse import OptionParser

# runs in: /afs/cern.ch/work/p/pedrok/private/SUSY2015/brilcalc

parser = OptionParser()
parser.add_option("-i", "--indir", dest="indir", default="json", help="input directory for JSON files (default = %default)")
(options, args) = parser.parse_args()

output = dict()

jsonfiles = glob.glob(options.indir+"/lumiSummary_*.json")
for f in jsonfiles:
    print f
    fbase = os.path.basename(f)
    fsplit = re.split('_|\.',fbase)
    # e.g. lumiSummary_unblind_DoubleMuon_2015D.json
    if len(fsplit)==5:
        category = fsplit[1]
        pd = fsplit[2]
        section = fsplit[3]
    # e.g. lumiSummary_DoubleMuon_2015D.json
    else:
        category = "all"
        pd = fsplit[1]
        section = fsplit[2]
    key = section+"_"+category
    
    cmd = "./brilcalc lumi --normtag /afs/cern.ch/user/c/cmsbril/public/normtag_json/OfflineNormtagV1.json -u /pb -i "+f
    brilout = os.popen(cmd).read()
    lumival = brilout.split()[-3]
    
    if not key in output:
        output[key] = [[pd,lumival]]
    else:
        output[key].append([pd,lumival])
        
#organized print
for key,val in output.iteritems():
    print key+"-"*10
    for v in val:
        print v[0]+": "+v[1]
