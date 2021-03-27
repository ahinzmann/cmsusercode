import os

fas=['fa1000','fa1500','fa2000','fa2500','fa3000','fa3500','fa4000','fa4500','fa5000','fa50000']
HT="2000toInf"

count=0

for fa in fas:
    samplename="alp_QCD_HT"+str(HT)
    string="python plotSignal_13TeV_desy_run2.py "+samplename+" "+str(fas.index(fa))
    if count%5!=4:
      string+=" &"
    print string
    count+=1
