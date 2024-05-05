import os
from ROOT import *
import array

massbinssets=[[(7000,13000)],
	      [(6000,7000)],
	      [(5400,6000)],
	      [(4800,5400)],
	      [(4200,4800)],
	      [(3600,4200)],
	      [(3000,3600)],
	      [(2400,3000)],
	      [(4800,5400),(5400,6000),(6000,7000),(7000,13000)],
	      [(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)],
	      ]

signal="QCD"
signalf="GEN-QCD-run2"
signalMass=""

jesSources=27 # 1 corresponds to the single overall variation, 27 UL (23EOY) to all
jerSources=1 # 1 corresponds to the single overall variation, 1 UL (6EOY) to all
correlatedSimUncertainties=False
uncorrelatedSimUncertainties=True
separateScaleUncertainties=False
alternateScaleUncertainty=False
uncorrelatedScaleUncertainties=False
scaleUncertainty=True
theoryStatUncertainties=True
useNNLO=True # choice for QCD
useM2=True # choice of mu-scale for QCD
runs="2" # "2" or "3" or "23"
run=runs[-1]

dire="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/"
if run=="2":
  if useM2:
    prefix=dire+"versions/run2ULNNLO_m2/datacard_shapelimit13TeV"
  elif useNNLO:
    prefix=dire+"versions/run2ULNNLO_pt12/datacard_shapelimit13TeV"
  else:
    prefix=dire+"versions/run2ULNLO_pt12/datacard_shapelimit13TeV"
elif run=="3":
  prefix=dire+"versions/run3ULNNLO_pt12/datacard_shapelimit13TeV"
else:
  whatprefix

for massbins in massbinssets:
    limits={}
    name=("_".join([s[0:4] for s in str(massbins).strip("[]").split("(")])).strip("_")
    print name
    
    statUncertainties=[]
    if theoryStatUncertainties:
      chi_bins=[(1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,3,6,9,12,16),
              ]
      for massbin in massbins:
         massindex=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)].index(massbin)
         for chibin in range(len(chi_bins[massindex])-1):
             statUncertainties+=["stat"+str(massindex)+"_"+str(chibin)]
    
    cfg=open("chi_datacard_"+name+"_bestfit_run2.txt","w")
    f=TFile(prefix+"_"+str(signalf)+str(signalMass)+"_chi.root")
    cfg.writelines("""
imax """+str(len(massbins))+""" number of channels
jmax 1 number of backgrounds
kmax """+str(2+3*correlatedSimUncertainties+3*len(massbins)*uncorrelatedSimUncertainties+jesSources+jerSources+1*scaleUncertainty+1*separateScaleUncertainties+(len(massbins)-1)*uncorrelatedScaleUncertainties+len(statUncertainties))+""" number of nuisance parameters""")
    cfg.writelines("""
-----------
""")
    for i in range(len(massbins)):
        cfg.writelines("""shapes * bin"""+str(i)+""" """+prefix+"""_"""+str(signalf)+str(signalMass)+"""_chi.root $PROCESS#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1 $PROCESS#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1_$SYSTEMATIC
""")
    cfg.writelines("""-----------
""")
    text="bin "
    for i in range(len(massbins)):
       text+="bin"+str(i)+" "
    text+="\nobservation "
    for i in range(len(massbins)):
       hData=f.Get("data_obs#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
       text+=str(hData.Integral())+" "
    cfg.writelines(text+"""
-----------
""")
    text="bin "
    for i in range(len(massbins)):
       text+="bin"+str(i)+" "+"bin"+str(i)+" "
    text+="\nprocess "
    for i in range(len(massbins)):
       text+=""+str(signal)+str(signalMass)+" "+str(signal)+str(signalMass)+"_ALT "
    text+="\nprocess "
    for i in range(len(massbins)):
       text+="0 1 "
    text+="\nrate "
    for i in range(len(massbins)):
       hQCD=f.Get("QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
       hALT=f.Get(""+str(signal)+str(signalMass)+"_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
       h=f.Get(""+str(signal)+str(signalMass)+"#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
       print "QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1",hQCD.Integral()
       print ""+str(signal)+str(signalMass)+"_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1",hALT.Integral()
       print ""+str(signal)+str(signalMass)+"#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1",h.Integral()
       text+=str(h.Integral())+" "+str(hALT.Integral())+" "
    cfg.writelines(text+"""
-----------
""")
    text=""
    ones=""
    for i in range(len(massbins)):
      ones+="1 1 - "
    if uncorrelatedSimUncertainties:
     for mn in massbins:
      text+="\nmodel"+str(mn[0])+" shape "+ones
     for mn in massbins:
      text+="\nJERtail"+str(mn[0])+" shape "+ones
     for mn in massbins:
      text+="\nsim"+str(mn[0])+" shape "+ones
    if correlatedSimUncertainties:
     text+="\nmodel shape "+ones
     text+="\nJERtail shape "+ones
     text+="\nsim shape "+ones
    if jesSources>1:
     for n in range(jesSources):
      text+="\njes"+str(n+1)+" shape "
      for i in range(len(massbins)):
         text+="1 1 "
    else:
      text+="\njes shape "
      for i in range(len(massbins)):
         text+="1 1 "
    if jerSources>1:
     for n in range(jerSources):
      text+="\njer"+str(n+1)+" shape "
      for i in range(len(massbins)):
         text+="1 1 "
    else:
      text+="\njer shape "
      for i in range(len(massbins)):
         text+="1 1 "
    text+="\nprefire shape "
    for i in range(len(massbins)):
        text+="1 1 "
    text+="\npdf shape "
    for i in range(len(massbins)):
       text+="1 1 "
    if separateScaleUncertainties:
      text+="\nscaleMuR shape "
      for i in range(len(massbins)):
         text+="1 1 "
      text+="\nscaleMuF shape "
      for i in range(len(massbins)):
         text+="1 1 "
    elif uncorrelatedScaleUncertainties:
     for mn in massbins:
      text+="\nscale"+str(mn[0])+" shape "
      for i in range(len(massbins)):
         text+="1 1 "
    else:
      if alternateScaleUncertainty:
        text+="\nscaleAlt shape "
      elif scaleUncertainty:
        text+="\nscale shape "
      for i in range(len(massbins)):
         text+="1 1 "
    for su in statUncertainties:
      text+="\n"+su+" shape "
      for i in range(len(massbins)):
         text+="1 1 "
    cfg.writelines(text+"""
-----------
""")

    cfg.close()
    #os.system("combine chi_datacard_"+name+"_bestfit.txt -M MaxLikelihoodFit -n "+name+"bestfit > bestfit"+name+".txt")
    #os.system("python diffNuisances.py -a mlfit"+name+"bestfit.root > bestfit"+name+"_nuisances.txt")
    os.system("combine chi_datacard_"+name+"_bestfit_run2.txt -M GoodnessOfFit --algo saturated --fixedSignalStrength=0 -n "+name+"goodnessfit |& tee goodnessfit"+name+".txt")
    for toy in range(0,20):
       command="combine chi_datacard_"+name+"_bestfit_run2.txt -M GoodnessOfFit --algo saturated --fixedSignalStrength=0 -t 450 --saveToys -n "+name+"goodnessfittoys_toy"+str(toy)+" |& tee goodnessfittoys"+name+"_toy"+str(toy)+".txt"
       if toy!=9 and toy!=18 and toy!=19:
          command+="&"
       os.system(command)
    os.system('hadd -f higgsCombine'+name+'goodnessfittoys.GoodnessOfFit.mH120.123456.root higgsCombine'+name+'goodnessfittoys_toy*.GoodnessOfFit.mH120.123456.root')
    os.system('root -q -b "'+dire+'/extractGoodnessStats.C(\\"'+name+'\\")" |& tee goodnessfittoysmulti'+name+'.txt')
