import os

signalMasses=[1000,1500,1750,2000,2250,2500,3000,3500,4000,4500,5000,6000,7000,8000]
minMass=1900
samples=[]

for signalMass in signalMasses:
  if signalMass==6000:
   mDMs=[1,2990]
  elif signalMass==7000:
   mDMs=[1,4000]
  elif signalMass==8000:
   mDMs=[1,3990]
  else:
   mDMs=[1,3000]
  for mDM in mDMs:
      for nxsec in range(13):
        samples+=[('Axial_Dijet_LO_Mphi',signalMass,mDM,("1p0","1p0"),nxsec,"Mar5"),]
      for nxsec in range(13):
        samples+=[('Vector_Dijet_LO_Mphi',signalMass,mDM,("1p0","1p0"),nxsec,"Mar5"),]
      for nxsec in range(111):# 111 for PDF uncertainties
        samples+=[('Vector_Dijet_LO_Mphi',signalMass,mDM,("1p0","1p0"),nxsec,"Jun26pdf"),]

#print samples

count=0

for sample,signalMass,mDM,coupling,nxsec,version in samples:
  
    samplename=sample+"_"+str(signalMass)+"_"+str(mDM)+"_"+coupling[0]+"_"+coupling[1]+"_"+version
    with open(samplename+str(nxsec)+".sh",'w+') as wrapper_script:
            wrapper_script.write("""#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/desy.de/user/h/hinzmann/rivet/CMSSW_10_6_16/src
cmsenv
cd /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis
python plotSignal_13TeV_desy_run2.py """+samplename+""" """+str(nxsec)+"""
""")
    with open(samplename+str(nxsec)+".submit",'w+') as htc_config:
            htc_config.write("""
#HTC Submission File for GEN sample production
#requirements      =  OpSysAndVer == "SL7"
universe          = vanilla
notification      = Error
notify_user       = andreas.hinzmann@desy.de
initialdir        = /nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/
output            = """+samplename+"""_"""+str(nxsec)+""".o
error             = """+samplename+"""_"""+str(nxsec)+""".e
log               = """+samplename+"""_"""+str(nxsec)+""".log
#Requesting CPU and DISK Memory - default +RequestRuntime of 3h stays unaltered
#+RequestRuntime   = 170000
#RequestMemory     = 10G
JobBatchName      = """+samplename+"""
#RequestDisk       = 10G
getenv            = True
executable        = /usr/bin/sh
arguments         = " """+samplename+str(nxsec)+""".sh"
queue 1
""")
    string="condor_submit "+samplename+str(nxsec)+".submit"
    if count%5!=4:
      string+=" &"
    print string
    count+=1

    #samplename=sample+"_"+str(signalMass)+"_"+str(mDM)+"_"+coupling[0]+"_"+coupling[1]+"_"+version
    #string="python plotSignal_13TeV_desy_run2.py "+samplename+" "+str(nxsec)
    #if count%5!=4:
    #  string+=" &"
    #print string
    #count+=1
