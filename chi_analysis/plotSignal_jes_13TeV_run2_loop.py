import os

count=0

for n in range(56,88): # 84 for data
    name="jes_"+str(n)
    with open(name+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /data/dust/user/hinzmann/jetmass/CMSSW_14_1_0_pre4/src/
cmsenv
cd /data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis 
python3 plotSignal_jes_13TeV_run2.py """+str(n)+"""
""")
    with open(name+".submit",'w+') as htc_config:
            htc_config.write("""
#HTC Submission File for GEN sample production
#requirements      =  OpSysAndVer == "SL7"
universe          = vanilla
notification      = Error
notify_user       = andreas.hinzmann@desy.de
initialdir        = /data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/
output            = """+name+""".o
error             = """+name+""".e
log               = """+name+""".log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
#+RequestRuntime   = 170000
+RequestRuntime   = 500000
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
    print(string)
    count+=1
