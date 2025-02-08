import os

CGs=["CG0p1","CG0p05","CG0p04","CG0p03","CG0p025","CG0p02","CG0p015","CG0p01","CG0p0075","CG0p005","CG0p0025","CG0p0"]

count=0

for CG in CGs:
    samplename="tripleG_QCD"
    with open(samplename+str(CGs.index(CG))+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/desy.de/user/h/hinzmann/rivet/CMSSW_10_6_16/src
cmsenv
cd /data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis
python plotSignal_13TeV_desy_run2.py tripleG_QCD """+str(CGs.index(CG))+"""
""")
    with open(samplename+str(CGs.index(CG))+".submit",'w+') as htc_config:
            htc_config.write("""
#HTC Submission File for GEN sample production
#requirements      =  OpSysAndVer == "SL7"
universe          = vanilla
notification      = Error
notify_user       = andreas.hinzmann@desy.de
initialdir        = /data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/
output            = tripleG_QCD_"""+CG+""".o
error             = tripleG_QCD_"""+CG+""".e
log               = tripleG_QCD_"""+CG+""".log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
+RequestRuntime   = 170000
#RequestMemory     = 10G
JobBatchName      = """+CG+"""
RequestDisk       = 10G
getenv            = True
executable        = /usr/bin/sh
arguments         = "tripleG_QCD"""+str(CGs.index(CG))+""".sh"
queue 1
""")
    string="condor_submit "+samplename+str(CGs.index(CG))+".submit"
    if count%5!=4:
      string+=" &"
    print string
    count+=1
