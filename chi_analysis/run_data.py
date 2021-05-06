import os

samples=[#("data2016",0,74,""),
         #("data2017",1,47,""),
	 #("data2018",2,120,""),
	 #("data2018",2,120,"-HEM"),
	 #("data2016",0,74,"-L1prefire"),
         #("data2017",1,47,"-L1prefire"),
	 #("QCD2016",3,7,""),
	 #("QCD2017",4,7,""),
	 #("QCD2018",5,7,""),
	 ("QCDpy2018",6,1,""),
	 ("QCDhw2018",7,1,""),
	 ]

count=0

for year,number,ns,postfix in samples:
  for n in range(ns):
    name="run_"+str(year)+"_"+str(n)+postfix
    with open(name+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis
cmsenv
python plot_data_13TeV_desy_run2.py """+str(number)+""" """+str(n)+""" """+postfix.replace("-","")+"""
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
JobBatchName      = run_"""+str(year)+postfix+"""
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
