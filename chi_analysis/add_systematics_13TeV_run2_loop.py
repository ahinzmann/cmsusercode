import os

count=0

for n in range(480):
    name="add_systematics_"+str(n)
    with open("submit/"+name+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis
cmsenv
python add_systematics_13TeV_run2.py """+str(n)+"""
""")
    with open("submit/"+name+".submit",'w+') as htc_config:
            htc_config.write("""
#HTC Submission File for GEN sample production
#requirements      =  OpSysAndVer == "SL7"
universe          = vanilla
notification      = Error
notify_user       = andreas.hinzmann@desy.de
initialdir        = /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/
output            = submit/"""+name+""".o
error             = submit/"""+name+""".e
log               = submit/"""+name+""".log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
"""+("+RequestRuntime   = 50000" if n==0 else "+RequestRuntime   = 2000")+"""
#RequestMemory     = 10G
JobBatchName      = add_systematics
#RequestDisk       = 10G
getenv            = True
executable        = /usr/bin/sh
arguments         = " submit/"""+name+""".sh"
queue 1
""")
    string="condor_submit submit/"+name+".submit"
    #string="source "+name+".sh"
    if count%5!=4:
      string+=" &"
    print string
    count+=1
