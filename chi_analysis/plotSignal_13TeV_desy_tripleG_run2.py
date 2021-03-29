import os

CGs=["CG0p1","CG0p05","CG0p04","CG0p03","CG0p025","CG0p02","CG0p015","CG0p01","CG0p005","CG0p0"]

count=0

for CG in CGs:
    samplename="tripleG_QCD"
    string="python plotSignal_13TeV_desy_run2.py "+samplename+" "+str(CGs.index(CG))
    if count%5!=4:
      string+=" &"
    print string
    count+=1
