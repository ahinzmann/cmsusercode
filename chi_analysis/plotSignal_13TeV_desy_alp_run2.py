import os

fas=['fa1000','fa1500','fa2000','fa2500','fa3000','fa3500','fa4000','fa4500','fa5000','fa50000']

count=0

for fa in fas:
    samplename="alp_QCD"
    with open(samplename+str(fas.index(fa))+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/desy.de/user/h/hinzmann/rivet/CMSSW_10_6_16/src
cmsenv
cd /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis
python plotSignal_13TeV_desy_run2.py """+samplename+""" """+str(fas.index(fa))+"""
""")
    with open(samplename+str(fas.index(fa))+".submit",'w+') as htc_config:
            htc_config.write("""
#HTC Submission File for GEN sample production
#requirements      =  OpSysAndVer == "SL7"
universe          = vanilla
notification      = Error
notify_user       = andreas.hinzmann@desy.de
initialdir        = /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/
output            = datacard_shapelimit13TeV_"""+samplename+"""_"""+fa+"""-run2_chi.root
error             = """+samplename+"""_"""+fa+""".e
log               = """+samplename+"""_"""+fa+""".log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
+RequestRuntime   = 170000
#RequestMemory     = 10G
JobBatchName      = """+fa+"""
RequestDisk       = 10G
getenv            = True
executable        = /usr/bin/sh
arguments         = " """+samplename+str(fas.index(fa))+""".sh"
queue 1
""")
    string="condor_submit "+samplename+str(fas.index(fa))+".submit"
    if count%5!=4:
      string+=" &"
    print string
    count+=1
