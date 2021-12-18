import os

count=0

for n in range(84):
    name="jes_"+str(n)
    with open(name+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis
cmsenv
python plotSignal_jes_13TeV_run2.py """+str(n)+"""
""")
    with open(name+".submit",'w+') as htc_config:
            htc_config.write("""
#HTC Submission File for GEN sample production
#requirements      =  OpSysAndVer == "SL7"
universe          = vanilla
notification      = Error
notify_user       = andreas.hinzmann@desy.de
initialdir        = /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/
output            = """+name+""".o
error             = """+name+""".e
log               = """+name+""".log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
+RequestRuntime   = 170000
#RequestMemory     = 10G
JobBatchName      = jes
#RequestDisk       = 10G
getenv            = True
executable        = /usr/bin/sh
arguments         = " """+name+""".sh"
queue 1
""")
    string="condor_submit "+name+".submit"
    if count%5!=4:
      string+=" &"
    print string
    count+=1
