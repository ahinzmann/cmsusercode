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
#	 ("QCDpy2018",6,1,""),
#	 ("QCDhw2018",7,1,""),
	 #("QCDpy2018",6,1,"-GEN"),
	 #("QCDhw2018",7,1,"-GEN"),
	 #("dataUL16preVFP",11,58,""),
         #("dataUL16postVFP",12,38,""),
         #("dataUL17",13,61,""),
	 #("dataUL18",14,97,""),
	 #("dataUL16preVFP",15,58,"-L1prefire"),
         #("dataUL16postVFP",16,38,"-L1prefire"),
         #("dataUL17",17,61,"-L1prefire"),
	 #("dataUL18",18,97,"-HEM"),
	 #("qcdUL16preVFP",19,105,""),
         #("qcdUL16postVFP",20,102,""),
         #("qcdUL17",21,115,""),
	 #("qcdUL18",22,159,""),
	 #("qcdUL16preVFP",19,105,"-GEN"),
         #("qcdUL16postVFP",20,102,"-GEN"),
         #("qcdUL17",21,115,"-GEN"),
	 #("qcdUL18",22,159,"-GEN"),
	 #("data2022",23,6,""),
	 #("data2023",24,12,""),
	 #("qcd2022",25,8,""),
	 #("qcd2022EE",26,8,""),
	 #("qcd2023",27,8,""),
	 #("qcd2023BPix",28,8,""),
	 #("dataNUL16preVFP",29,6,"-L1prefire"),
         #("dataNUL16postVFP",30,3,"-L1prefire"),
         #("dataNUL17",31,5,"-L1prefire"),
	 #("dataNUL18",32,4,"-HEM"),
	 #("data2023",33,12,"-28May2024"),
	 #("qcd2023",34,8,"-28May2024"),
	 #("qcd2023BPix",35,8,"-28May2024"),
	 #("qbhUL18",36,24,""),
	 ("data2024",37,18,""),
	 ]

count=0

for year,number,ns,postfix in samples:
  for n in range(ns):
    name="run_"+str(year)+"_"+str(n)+postfix
    with open(name+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /data/dust/user/hinzmann/jetmass/CMSSW_14_1_0_pre4/src/
cmsenv
cd /data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis
python3 plot_data_13TeV_desy_run2.py """+str(number)+""" """+str(n)+""" """+postfix.replace("-","")+"""
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
    print(string)
    count+=1
