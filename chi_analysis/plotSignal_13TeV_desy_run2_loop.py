import os

#names=["","4","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34"]
#names=["CP5","run3","18-run3","19-run3","20-run3","21-run3","22-run3","23-run3","24-run3","25-run3","26-run3","27-run3","28-run3","29-run3","30-run3","31-run3","32-run3","33-run3"]
names=["0-run3","1-run3","2-run3","3-run3","4-run3","5-run3"]

count=0

for name in names:
    samplename="pythia8_GEN"
    with open(samplename+name+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/desy.de/user/h/hinzmann/rivet/CMSSW_10_6_16/src
cmsenv
cd /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis
python plotSignal_13TeV_desy_run2.py """+name+"""
""")
    with open(samplename+name+".submit",'w+') as htc_config:
            htc_config.write("""
#HTC Submission File for GEN sample production
#requirements      =  OpSysAndVer == "SL7"
universe          = vanilla
notification      = Error
notify_user       = andreas.hinzmann@desy.de
initialdir        = /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/
output            = """+samplename+"""_"""+name+""".o
error             = """+samplename+"""_"""+name+""".e
log               = """+samplename+"""_"""+name+""".log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
+RequestRuntime   = 170000
#RequestMemory     = 10G
JobBatchName      = """+samplename+"""_"""+name+"""
RequestDisk       = 10G
getenv            = True
executable        = /usr/bin/sh
arguments         = " """+samplename+name+""".sh"
queue 1
""")
    string="condor_submit "+samplename+name+".submit"
    if count%5!=4:
      string+=" &"
    print string
    count+=1
