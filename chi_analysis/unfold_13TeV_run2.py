import os,sys
from ROOT import *
import array
import ROOT
import subprocess

def system_call(command):
    print command
    p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    return p.stdout.read()

massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]

jesSources=27 # 1 corresponds to the single overall variation, 27 UL (23EOY) to all
jerSources=1 # 1 corresponds to the single overall variation, 1 UL (6EOY) to all
correlatedSimUncertainties=False
uncorrelatedSimUncertainties=True
scaleUncertainties=False
separateScaleUncertainties=False
alternateScaleUncertainty=False
theoryStatUncertainties=False

trivialClosure=False
crossClosure=False
withUncertainties=True
runs="2" # "2" or "3" or "23"
run=runs[-1]

dire="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/"
prefix="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/versions/run"+run+"ULNNLO_m2/datacard_shapelimit13TeV"

name="unfold"
if withUncertainties: name+="_withUncertainties"
if trivialClosure: name+="_trivialClosure"
if crossClosure: name+="_crossClosure"
print name
statUncertainties=[]
chi_bins=[(1,2,3,4,5,6,7,8,9,10,12,14,16),
           (1,2,3,4,5,6,7,8,9,10,12,14,16),
           (1,2,3,4,5,6,7,8,9,10,12,14,16),
           (1,2,3,4,5,6,7,8,9,10,12,14,16),
           (1,2,3,4,5,6,7,8,9,10,12,14,16),
           (1,2,3,4,5,6,7,8,9,10,12,14,16),
           (1,2,3,4,5,6,7,8,9,10,12,14,16),
           (1,3,6,9,12,16),
          ]
all_mass_bins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
mass_bin_offset=len(all_mass_bins)-len(massbins)
if theoryStatUncertainties:
  for massbin in massbins:
     massindex=all_mass_bins.index(massbin)
     for chibin in range(len(chi_bins[massindex-mass_bin_offset])-1):
         statUncertainties+=["stat"+str(massindex)+"_"+str(chibin)]

fname=prefix+"_GEN-QCD-run2_chi.root"

if True:
    f=TFile(fname)
    cfg=open("chi_datacard13TeV"+name+"_run"+run+".txt","w")
    nbins=0
    for j in range(len(massbins)):
      for chibin in range(len(chi_bins[j])-1):
        nbins+=1
    cfg.writelines("""
imax """+str(len(massbins))+""" number of channels (reco mass bins)
jmax """+str(nbins-1)+""" number of samples (gen bins)
kmax """+str(withUncertainties*(4+correlatedSimUncertainties+len(massbins)*uncorrelatedSimUncertainties+jesSources+jerSources+1*separateScaleUncertainties+1*scaleUncertainties+len(statUncertainties)))+""" number of nuisance parameters""")
    cfg.writelines("""
-----------
""")
    for i in range(len(massbins)):
      if trivialClosure: # use LO QCD instead of data and NNLO QCD
        cfg.writelines("""shapes data_obs recomass"""+str(i)+""" """+fname+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"""_$SYSTEMATIC
""")
        cfg.writelines("""shapes * recomass"""+str(i)+""" """+fname+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$PROCESS_"+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$SYSTEMATIC"+"""_$PROCESS_
""")
      elif crossClosure: # use NNLO QCD as data
        cfg.writelines("""shapes data_obs recomass"""+str(i)+""" """+fname+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"""_$SYSTEMATIC
""")
        cfg.writelines("""shapes * recomass"""+str(i)+""" """+fname+""" QCD_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$PROCESS_"+""" QCD_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$SYSTEMATIC"+"""_$PROCESS_
""")
      else:
        cfg.writelines("""shapes data_obs recomass"""+str(i)+""" """+fname+""" data_obs#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+""" data_obs#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"""_$SYSTEMATIC
""")
        cfg.writelines("""shapes * recomass"""+str(i)+""" """+fname+""" QCD_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$PROCESS_"+""" QCD_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$SYSTEMATIC"+"""_$PROCESS_
""")
    cfg.writelines("""
-----------
""")
    text="bin "
    for i in range(len(massbins)):
       text+="recomass"+str(i)+" "
    text+="\nobservation "
    for i in range(len(massbins)):
       if trivialClosure or crossClosure:
         hData=f.Get("QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
         print "QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1",hData.Integral()
       else:
         hData=f.Get("data_obs#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
         print "data_obs#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1",hData.Integral()
       text+=str(hData.Integral())+" "
    cfg.writelines(text+"""
-----------
""")
    text="bin "
    for i in range(len(massbins)):
     for j in range(len(massbins)):
         text+=("recomass"+str(i)+" ")*(len(chi_bins[j])-1)
    text+="\nprocess "
    for i in range(len(massbins)):
     for j in range(len(massbins)):
       for chibin in range(len(chi_bins[j])-1):
        text+="bin_"+str(j)+"_"+str(chibin)+" "
    text+="\nprocess "
    for i in range(len(massbins)):
     n=0
     for j in range(len(massbins)):
        text+=(str(n)+" ")*(len(chi_bins[j])-1)
        n-=1
    text+="\nrate "
    ones=""
    for i in range(len(massbins)):
     for j in range(len(massbins)):
       for chibin in range(len(chi_bins[j])-1):
        if trivialClosure:
          hQCD=f.Get("QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_")
          print "QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_",hQCD.Integral()
        else:
          hQCD=f.Get("QCD_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_")
          print "QCD_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_",hQCD.Integral()
        text+=str(hQCD.Integral())+" "
        if hQCD.Integral()>0:
          ones+="1 "
        else:
          ones+="0 "
    cfg.writelines(text+"""
-----------
""")
    text=""

    if withUncertainties:
      print "adding uncertainties"
      text+="\nmodel shape "+ones
      if runs=="23": text+="\nnuisance edit rename * * model run"+run+"_model"
      text+="\nJERtail shape "+ones
      if runs=="23": text+="\nnuisance edit rename * * JERtail run"+run+"_JERtail"
      if uncorrelatedSimUncertainties:
       for mn in massbins:
        text+="\nsim"+str(mn[0])+" shape "+ones
        if runs=="23": text+="\nnuisance edit rename * * sim"+str(mn[0])+" run"+run+"_sim"+str(mn[0])
      if correlatedSimUncertainties:
       text+="\nsim shape "+ones
       if runs=="23": text+="\nnuisance edit rename * * sim run"+run+"_sim"
      if jesSources>1:
       for n in range(jesSources):
        text+="\njes"+str(n+1)+" shape "+ones
        if runs=="23" and n>1: text+="\nnuisance edit rename * * jes"+str(n+1)+" run"+run+"_jes"+str(n+1) # if n<2 combine would do the replace twice for n>10 and >20
      else:
        text+="\njes shape "+ones
        if runs=="23": text+="\nnuisance edit rename * * jes run"+run+"_jes"
      if jerSources>1:
       for n in range(jerSources):
        text+="\njer"+str(n+1)+" shape "+ones
        if runs=="23" and n>1: text+="\nnuisance edit rename * * jer"+str(n+1)+" run"+run+"_jer"+str(n+1) # if n<2 combine would do the replace twice for n>10 and >20
      else:
        text+="\njer shape "+ones
        if runs=="23": text+="\nnuisance edit rename * * jer run"+run+"_jer"
      text+="\nprefire shape "+ones
      if runs=="23": text+="\nnuisance edit rename * * prefire run"+run+"_prefire"
      text+="\npdf shape "+ones
      if separateScaleUncertainties:
        text+="\nscaleMuR shape "+ones
        text+="\nscaleMuF shape "+ones
      else:
        if alternateScaleUncertainty:
          text+="\nscaleAlt shape "+ones
        elif scaleUncertainties:
          text+="\nscale shape "+ones
      print "adding stat uncertainties"
      for su in statUncertainties:
        text+="\n"+su+" shape "+ones
      cfg.writelines(text+"""
-----------
""")

    cfg.close()

if True:
    if runs=="23":
      out=system_call("combineCards.py run2=chi_datacard13TeV"+name+"_run2.txt run3=chi_datacard13TeV"+name+"_run3.txt > chi_datacard13TeV"+name+"_run23.txt")
    pois_map=""
    for j in range(len(massbins)):
      for chibin in range(len(chi_bins[j])-1):
        pois_map+=" --PO map='.*"+"bin_"+str(j)+"_"+str(chibin)+".*:r_Bin_"+str(j)+"_"+str(chibin)+"[1,-1,10]'"
    print "running combine"
    out=system_call("text2workspace.py -m 125 chi_datacard13TeV"+name+"_run"+runs+".txt --X-allow-no-background -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "+pois_map+" -o "+name+"_run"+runs+".root")
    print out
    out=system_call("combine -m 125 -M MultiDimFit "+name+"_run"+runs+".root --saveFitResult")
    print out
    out=system_call("mv multidimfit.root multidimfit_"+name+"_run"+runs+".root")
    print out
