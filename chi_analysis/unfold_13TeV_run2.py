import os,sys
from ROOT import *
import array
import ROOT
import subprocess

def system_call(command):
    print command
    p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    return p.stdout.read()

def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text  # or whatever
 
massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]

jesSources=27 # 1 corresponds to the single overall variation, 27 UL (23EOY) to all
jerSources=1 # 1 corresponds to the single overall variation, 1 UL (6EOY) to all
correlatedSimUncertainties=False
uncorrelatedSimUncertainties=True
separateScaleUncertainties=False
alternateScaleUncertainty=False
theoryStatUncertainties=True
withUncertainties=False

runs="2" # "2" or "3" or "23"
run=runs[-1]

dire="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/"
prefix="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/versions/run"+run+"ULNNLO_pt12/datacard_shapelimit13TeV"

name="unfold"
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
all_mass_bins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
mass_bin_offset=len(all_mass_bins)-len(massbins)
for massbin in massbins:
     massindex=all_mass_bins.index(massbin)
     for chibin in range(len(chi_bins[massindex])-1):
         statUncertainties+=["stat"+str(massindex)+"_"+str(chibin)]

fname=prefix+"_GEN-QCD-run2_chi.root"

if True:
    f=TFile(fname)
    cfg=open("chi_datacard13TeVunfold_run"+run+".txt","w")
    nbins=0
    for j in range(len(massbins)):
      for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
        nbins+=1
    cfg.writelines("""
imax """+str(len(massbins))+""" number of channels
jmax """+str(nbins-1)+""" number of backgrounds
kmax """+str(withUncertainties*(5+correlatedSimUncertainties+len(massbins)*uncorrelatedSimUncertainties+jesSources+jerSources+1*separateScaleUncertainties+len(statUncertainties)))+""" number of nuisance parameters""")
    cfg.writelines("""
-----------
""")
    for i in range(len(massbins)):
        cfg.writelines("""shapes data_obs recomass"""+str(i)+""" """+fname+""" data_obs#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+""" data_obs#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"""_$SYSTEMATIC
""")
        cfg.writelines("""shapes * recomass"""+str(i)+""" """+fname+""" QCD_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$PROCESS_"+""" QCD_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1"""+"_$PROCESS_"+"""_$SYSTEMATIC
""")
    cfg.writelines("""
-----------
""")
    text="bin "
    for i in range(len(massbins)):
       text+="recomass"+str(i)+" "
    text+="\nobservation "
    print "data_obs#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"
    for i in range(len(massbins)):
       hData=f.Get("data_obs#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
       text+=str(hData.Integral())+" "
    cfg.writelines(text+"""
-----------
""")
    text="bin "
    for i in range(len(massbins)):
     for j in range(len(massbins)):
       for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
         text+="recomass"+str(i)+" "
    text+="\nprocess "
    for i in range(len(massbins)):
     for j in range(len(massbins)):
       for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
        text+="bin_"+str(j+mass_bin_offset)+"_"+str(chibin)+" "
    text+="\nprocess "
    for i in range(len(massbins)):
     n=0
     for j in range(len(massbins)):
       for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
        text+=str(n)+" "
        n-=1
    text+="\nrate "
    for i in range(len(massbins)):
     for j in range(len(massbins)):
       for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
        hQCD=f.Get("QCD_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j+mass_bin_offset)+"_"+str(chibin)+"_")
        print "QCD_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j+mass_bin_offset)+"_"+str(chibin)+"_",hQCD.Integral()
        text+=str(hQCD.Integral())+" "
    cfg.writelines(text+"""
-----------
""")
    text=""

    #for j1 in range(len(massbins)):
    #  for chibin1 in range(len(chi_bins[j1+mass_bin_offset])-1):
    #    text+="\nr_bin_"+str(j1)+"_"+str(chibin1)+" rate "
    #    print "\nr_bin_"+str(j1)+"_"+str(chibin1)+" rate "
    #    for i in range(len(massbins)):
    #     for j in range(len(massbins)):
    #       for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
    #        if j1==j and chibin1==chibin:
    #          text+="1 "
    #        else:
    #          text+="0 "
    #    if runs=="23": text+="\nnuisance edit rename * * model run"+run+"_model"

    if withUncertainties:
      text+="\nmodel shape "
      for i in range(len(massbins)):
       for j in range(len(massbins)):
         for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
          text+="1 "
      if runs=="23": text+="\nnuisance edit rename * * model run"+run+"_model"
      text+="\nJERtail shape "
      for i in range(len(massbins)):
       for j in range(len(massbins)):
         for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
          text+="1 "
      if runs=="23": text+="\nnuisance edit rename * * JERtail run"+run+"_JERtail"
      if uncorrelatedSimUncertainties:
       for mn in massbins:
        text+="\nsim"+str(mn[0])+" shape "
        for i in range(len(massbins)):
         for j in range(len(massbins)):
          for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
           text+="1 "
        if runs=="23": text+="\nnuisance edit rename * * sim"+str(mn[0])+" run"+run+"_sim"+str(mn[0])
      if correlatedSimUncertainties:
       text+="\nsim shape "
       for i in range(len(massbins)):
        for j in range(len(massbins)):
         for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
          text+="1 "
       if runs=="23": text+="\nnuisance edit rename * * sim run"+run+"_sim"
      if jesSources>1:
       for n in range(jesSources):
        text+="\njes"+str(n+1)+" shape "
        for i in range(len(massbins)):
         for j in range(len(massbins)):
          for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
           text+="1 "
        if runs=="23" and n>1: text+="\nnuisance edit rename * * jes"+str(n+1)+" run"+run+"_jes"+str(n+1) # if n<2 combine would do the replace twice for n>10 and >20
      else:
        text+="\njes shape "
        for i in range(len(massbins)):
         for j in range(len(massbins)):
          for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
           text+="1 "
        if runs=="23": text+="\nnuisance edit rename * * jes run"+run+"_jes"
      if jerSources>1:
       for n in range(jerSources):
        text+="\njer"+str(n+1)+" shape "
        for i in range(len(massbins)):
         for j in range(len(massbins)):
          for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
           text+="1 "
        if runs=="23" and n>1: text+="\nnuisance edit rename * * jer"+str(n+1)+" run"+run+"_jer"+str(n+1) # if n<2 combine would do the replace twice for n>10 and >20
      else:
        text+="\njer shape "
        for i in range(len(massbins)):
         for j in range(len(massbins)):
          for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
           text+="1 "
        if runs=="23": text+="\nnuisance edit rename * * jer run"+run+"_jer"
      text+="\nprefire shape "
      for i in range(len(massbins)):
       for j in range(len(massbins)):
         for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
          text+="1 "
      if runs=="23": text+="\nnuisance edit rename * * prefire run"+run+"_prefire"
      text+="\npdf shape "
      for i in range(len(massbins)):
       for j in range(len(massbins)):
         for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
           text+="1 "
      if separateScaleUncertainties:
        text+="\nscaleMuR shape "
        for i in range(len(massbins)):
         for j in range(len(massbins)):
          for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
            text+="1 "
        text+="\nscaleMuF shape "
        for i in range(len(massbins)):
         for j in range(len(massbins)):
          for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
            text+="1 "
      else:
        if alternateScaleUncertainty:
          text+="\nscaleAlt shape "
        else:
          text+="\nscale shape "
        for i in range(len(massbins)):
         for j in range(len(massbins)):
          for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
            text+="1 "
      for su in statUncertainties:
        text+="\n"+su+" shape "
        for i in range(len(massbins)):
         for j in range(len(massbins)):
          for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
           text+="1 "
      cfg.writelines(text+"""
-----------
""")

    cfg.close()

    if runs=="23":
      out=system_call("combineCards.py run2=chi_datacard13TeVunfold_run2.txt run3=chi_datacard13TeVunfold_run3.txt > chi_datacard13TeVunfold_run23.txt")
    pois_map=""
    for j in range(len(massbins)):
      for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
        pois_map+=" --PO map='.*"+"bin_"+str(j+mass_bin_offset)+"_"+str(chibin)+".*:r_Bin_"+str(j+mass_bin_offset)+"_"+str(chibin)+"[1,-1,10]'"
    out=system_call("text2workspace.py -m 125 chi_datacard13TeVunfold_run"+runs+".txt --X-allow-no-background -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "+pois_map+" -o unfold_run"+runs+".root")
    out=system_call("combine -m 125 -M MultiDimFit unfold_run"+runs+".root --saveFitResult")
    #set_params=" --setParameters="
    #for j in range(len(massbins)):
    #  for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
    #    set_params+="r_Bin_"+str(j+mass_bin_offset)+"_"+str(chibin)+"=1,"
    #out=system_call("combine -m 125 -M MultiDimFit "+set_params.strip(",")+" -t -1 unfold_run"+runs+".root")
    #out=system_call("mkdir "+name)
    #out=system_call("combine -m 125 -M MaxLikelihoodFit --plots --out "+name+" -n unfold unfold_run"+runs+".root")
    #out=system_call("python diffNuisances.py -p x -a "+name+"/fitDiagnosticsunfold.root -A")
    print out
