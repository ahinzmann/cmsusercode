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

dire="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/"
prefix="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/versions/run2ULNNLONov21/datacard_shapelimit13TeV"

for massbins in massbinssets:
    limits={}
    name=("_".join([s[0:4] for s in str(massbins).strip("[]").split("(")])).strip("_")
    print name
    f1=open("goodnessfit"+name+".txt")
    observed=0
    for l in f1.readlines():
        if "Best fit test statistic: " in l:
            observed=float(l.strip(" ").split(" ")[-1])
    f2=open("goodnessfittoysmulti"+name+".txt")
    expected=0
    expected68up=0
    expected68down=0
    for l in f2.readlines():
        if "Best fit test statistic: " in l or "Expected median" in l:
            expected=float(l.strip(" ").split(" ")[-1])
        if "68%" in l:
            expected68up=float(l.strip(" ").split(" ")[-1])
            expected68down=float(l.strip(" ").split(" ")[-5])
        if "Expected - 1sigma" in l:
            expected68down=float(l.strip(" ").split(" ")[-1])
        if "Expected + 1sigma" in l:
            expected68up=float(l.strip(" ").split(" ")[-1])
        if "Expected - 2sigma" in l:
            expected95down=float(l.strip(" ").split(" ")[-1])
        if "Expected + 2sigma" in l:
            expected95up=float(l.strip(" ").split(" ")[-1])
    sigma=0
    if observed>=expected:
       sigma=(observed-expected)/max((expected68up-expected),(expected-expected68down))
    else:
       sigma=(observed-expected)/(expected68down-expected)
    if sigma>=2:
      if observed>=expected:
       sigma=2.*((observed-expected)/max((expected95up-expected),(expected-expected95down)))
      else:
       sigma=2.*((observed-expected)/(expected95down-expected))
    print observed, expected, expected68down, expected68up, sigma