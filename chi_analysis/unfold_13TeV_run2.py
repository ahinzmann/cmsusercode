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

jesSources=23 # 1 corresponds to the single overall variation, 27 UL (23EOY) to all
jerSources=1 # 1 corresponds to the single overall variation, 1 UL (6EOY) to all
correlatedSimUncertainties=False
uncorrelatedSimUncertainties=True
scaleUncertainties=False
separateScaleUncertainties=False
alternateScaleUncertainty=False
theoryStatUncertainties=False

trivialClosure=False
crossClosureNNLO=False
crossClosureHerwig=False
crossClosureMadgraph=False
withUncertainties=True
runs="2" # "2" or "3" or "23"
run=runs[-1]

dire="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/"
prefix="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/versions/run"+run+"ULNNLO_m2/datacard_shapelimit13TeV"

name="unfold"
if withUncertainties: name+="_withUncertainties"
if trivialClosure: name+="_trivialClosure"
if crossClosureNNLO: name+="_crossClosureNNLO"
if crossClosureHerwig: name+="_crossClosureHerwig"
if crossClosureMadgraph: name+="_crossClosureMadgraph"
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
    fh=TFile(fname.replace("shapelimit13TeV","shapelimit13TeVherwigpp"))
    fm=TFile(fname.replace("shapelimit13TeV","shapelimit13TeVmadgraphMLM"))
    cfg=open("chi_datacard13TeV"+name+"_run"+run+".txt","w")
    nbins=0
    for j in range(len(massbins)):
      for chibin in range(len(chi_bins[j])-1):
        nbins+=1
    cfg.writelines("""
imax """+str(len(massbins))+""" number of channels (reco mass bins)
jmax """+str(nbins-1)+""" number of samples (gen bins)
kmax """+str(withUncertainties*(2+3*correlatedSimUncertainties+3*len(massbins)*uncorrelatedSimUncertainties+jesSources+jerSources+1*separateScaleUncertainties+1*scaleUncertainties+len(statUncertainties)))+""" number of nuisance parameters""")
    cfg.writelines("""
-----------
""")
    for i in range(len(massbins)):
      if trivialClosure: # use Pythia as data and Pythia as background
        cfg.writelines("""shapes data_obs recomass"""+str(i)+""" """+fname+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"""_$SYSTEMATIC
""")
        cfg.writelines("""shapes * recomass"""+str(i)+""" """+fname+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$PROCESS_"+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$SYSTEMATIC"+"""_$PROCESS_
""")
      elif crossClosureNNLO: # use Pythia as data and NNLO as background
        cfg.writelines("""shapes data_obs recomass"""+str(i)+""" """+fname+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"""_$SYSTEMATIC
""")
        cfg.writelines("""shapes * recomass"""+str(i)+""" """+fname+""" QCD_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$PROCESS_"+""" QCD_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$SYSTEMATIC"+"""_$PROCESS_
""")
      elif crossClosureHerwig: # use Pythia as data and Herwig-smeared Pythia as background
        cfg.writelines("""shapes data_obs recomass"""+str(i)+""" """+fname+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"""_$SYSTEMATIC
""")
        cfg.writelines("""shapes * recomass"""+str(i)+""" """+fname.replace("shapelimit13TeV","shapelimit13TeVherwigpp")+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$PROCESS_"+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$SYSTEMATIC"+"""_$PROCESS_
""")
      elif crossClosureMadgraph: # use Pythia as data and Madgraph-smeared Pythia as background
        cfg.writelines("""shapes data_obs recomass"""+str(i)+""" """+fname+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"""_$SYSTEMATIC
""")
        cfg.writelines("""shapes * recomass"""+str(i)+""" """+fname.replace("shapelimit13TeV","shapelimit13TeVmadgraphMLM")+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$PROCESS_"+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$SYSTEMATIC"+"""_$PROCESS_
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
       if trivialClosure or crossClosureNNLO or crossClosureHerwig or crossClosureMadgraph:
         hData=f.Get("QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
         print "QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1",hData.Integral()
       else:
         hData=f.Get("data_obs#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
         print "data_obs#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1",hData.Integral()
       text+=str(hData.Integral())+" "
    cfg.writelines(text+"""
-----------
""")
    bintext="bin "
    process1text="\nprocess "
    process2text="\nprocess "
    ratetext="\nrate "
    ones=""
    skip={}
    for i in range(len(massbins)):
     n=0
     for j in range(len(massbins)):
       for chibin in range(len(chi_bins[j])-1):
        if trivialClosure:
          hQCD=f.Get("QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_")
          print "QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_",hQCD.Integral()
        elif crossClosureHerwig:
          hQCD=fh.Get("QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_")
          print "QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_",hQCD.Integral()
        elif crossClosureMadgraph:
          hQCD=fm.Get("QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_")
          print "QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_",hQCD.Integral()
        else:
          hQCD=f.Get("QCD_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_")
          print "QCD_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibin)+"_",hQCD.Integral()
        if hQCD.Integral()>0:
          ones+="1 "
        else:
          ones+="0 "
        bintext+=("recomass"+str(i)+" ")
        process1text+="bin_"+str(j)+"_"+str(chibin)+" "
        process2text+=(str(n)+" ")
        ratetext+=str(hQCD.Integral())+" "
       n-=1
    text=bintext+process1text+process2text+ratetext
    cfg.writelines(text+"""
-----------
""")
    text=""

    if withUncertainties:
      print "adding uncertainties"
      if uncorrelatedSimUncertainties:
       for mn in massbins:
        text+="\nmodel"+str(mn[0])+" shape "+ones
        if runs=="23": text+="\nnuisance edit rename * * model"+str(mn[0])+" run"+run+"_model"+str(mn[0])
       for mn in massbins:
        text+="\nJERtail"+str(mn[0])+" shape "+ones
        if runs=="23": text+="\nnuisance edit rename * * JERtail"+str(mn[0])+" run"+run+"_JERtail"+str(mn[0])
       for mn in massbins:
        text+="\nsim"+str(mn[0])+" shape "+ones
        if runs=="23": text+="\nnuisance edit rename * * sim"+str(mn[0])+" run"+run+"_sim"+str(mn[0])
      if correlatedSimUncertainties:
       text+="\nmodel shape "+ones
       if runs=="23": text+="\nnuisance edit rename * * model run"+run+"_model"
       text+="\nJERtail shape "+ones
       if runs=="23": text+="\nnuisance edit rename * * JERtail run"+run+"_JERtail"
       text+="\nsim shape "+ones
       if runs=="23": text+="\nnuisance edit rename * * sim run"+run+"_sim"
      if jesSources>1:
       for n in range(jesSources):
        text+="\njes"+str(n+1)+" shape "+ones
        if runs=="23" and n>1: text+="\nnuisance edit rename * * jes"+str(n+1)+" run"+run+"_jes"+str(n+1) # if n<2 combine would do the replace twice for n>10 and >20
      elif jesSources==1:
        text+="\njes shape "+ones
        if runs=="23": text+="\nnuisance edit rename * * jes run"+run+"_jes"
      if jerSources>1:
       for n in range(jerSources):
        text+="\njer"+str(n+1)+" shape "+ones
        if runs=="23" and n>1: text+="\nnuisance edit rename * * jer"+str(n+1)+" run"+run+"_jer"+str(n+1) # if n<2 combine would do the replace twice for n>10 and >20
      elif jerSources==1:
        text+="\njer shape "+ones
        if runs=="23": text+="\nnuisance edit rename * * jer run"+run+"_jer"
      text+="\nprefire shape "+ones
      if runs=="23": text+="\nnuisance edit rename * * prefire run"+run+"_prefire"
      text+="\ntrigger shape "+ones
      if runs=="23": text+="\nnuisance edit rename * * trigger run"+run+"_trigger"
      #text+="\npdf shape "+ones
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
        pois_map+=" --PO map='.*"+"bin_"+str(j)+"_"+str(chibin)+".*:r_Bin_"+str(j)+"_"+(str(chibin) if chibin>9 else "0"+str(chibin))+"[1,0,3]'"
    print "running combine"
    out=system_call("text2workspace.py -m 125 chi_datacard13TeV"+name+"_run"+runs+".txt --X-allow-no-background -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "+pois_map+" -o "+name+"_run"+runs+".root")
    print out
    out=system_call("mkdir "+name+"_run"+runs)
    print out
    out=system_call("cd "+name+"_run"+runs+";combine -m 125 -M MultiDimFit ../"+name+"_run"+runs+".root --saveFitResult")
    print out
    out=system_call("mv "+name+"_run"+runs+"/multidimfit.root multidimfit_"+name+"_run"+runs+".root")
    print out
