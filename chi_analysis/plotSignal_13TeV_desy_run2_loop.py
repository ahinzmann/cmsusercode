import os

#names=["","4","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34"]
#names=["CP5","run3"]
#names=["18-run3","19-run3","20-run3","21-run3","22-run3","23-run3","24-run3","25-run3","26-run3","27-run3","28-run3","29-run3","30-run3","31-run3","32-run3","33-run3"]
#names=["0-run3","1-run3","2-run3","3-run3","4-run3","5-run3"]
#names=["18-CP5","19-CP5","20-CP5","21-CP5","22-CP5","23-CP5","24-CP5","25-CP5","26-CP5","27-CP5","28-CP5","29-CP5","30-CP5","31-CP5","32-CP5","33-CP5"] #add
names=["34-CP5","35-CP5","36-CP5","37-CP5","38-CP5","39-CP5","40-CP5","41-CP5","42-CP5","43-CP5","44-CP5","45-CP5","46-CP5","47-CP5","48-CP5","49-CP5"] #qbh
#names=["18-CP2","19-CP2","20-CP2","21-CP2","22-CP2","23-CP2","24-CP2","25-CP2","26-CP2","27-CP2","28-CP2","29-CP2","30-CP2","31-CP2","32-CP2","33-CP2"] #add

count=0

for name in names:
    samplename="pythia8_GEN"
    with open("submit/"+samplename+name+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /data/dust/user/hinzmann/jetmass/CMSSW_14_1_0_pre4/src/
cmsenv
cd /data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis
python3 plotSignal_13TeV_desy_run2.py """+name+"""
""")
    with open("submit/"+samplename+name+".submit",'w+') as htc_config:
            htc_config.write("""
#HTC Submission File for GEN sample production
#requirements      =  OpSysAndVer == "SL7"
universe          = vanilla
notification      = Error
notify_user       = andreas.hinzmann@desy.de
initialdir        = /data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/
output            = submit/"""+samplename+"""_"""+name+""".o
error             = submit/"""+samplename+"""_"""+name+""".e
#log               = submit/"""+samplename+"""_"""+name+""".log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
+RequestRuntime   = 170000
#RequestMemory     = 10G
JobBatchName      = """+samplename+"""_"""+name+"""
RequestDisk       = 10G
getenv            = True
executable        = /usr/bin/sh
arguments         = " submit/"""+samplename+name+""".sh"
queue 1
""")
    string="condor_submit submit/"+samplename+name+".submit"
    if count%5!=4:
      string+=" &"
    print(string)
    count+=1
