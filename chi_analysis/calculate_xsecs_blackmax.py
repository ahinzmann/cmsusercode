import os
print "export SCRAM_ARCH=slc7_amd64_gcc10"
qbhlist=[]
for md in [2000,3000,4000,5000,6000,7000,8000,9000]:
  for n in [2,4,6]:
   qbhlist+=[(md,md+1000,n)]
for md,mbh,n in qbhlist:
  #print "source blackmax.sh "+str(md)+" "+str(mbh)+" "+str(n)+" 3 13000 &"
  print "cd /data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/blackmax; source runcmsgrid.sh 1 1 1 "+str(md)+" "+str(mbh)+" "+str(n)+" 3 > QBH_MD"+str(md)+"_MBH"+str(mbh)+"_n"+str(n)+".txt"

xsecs={}
for md,mbh,n in qbhlist:
  f=open("blackmax/QBH_MD"+str(md)+"_MBH"+str(mbh)+"_n"+str(n)+".txt")
  for l in f.readlines():
    if "cross section" in l:
       xsecs["QBH_MD"+str(md)+"_MBH"+str(mbh)+"_n"+str(n)]=float(l.split(" ")[-9])
       if n==6: print md, float(l.split(" ")[-9])
print xsecs
