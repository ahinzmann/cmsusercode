import os
#d="/data/dust/user/hinzmann/dijetangular/13p6TeV/"
d="/data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/"
files=os.listdir(d)
xsecs={}

for name in files:
   if ".o" in name and "Nov2022" in name:
     f=open(d+name)
     for l in f.readlines():
       if "| sum" in l:
          xsecs[name.split("Nov2022")[-2]+"Nov2022"]=float(l.split(" ")[-4])

print xsecs
