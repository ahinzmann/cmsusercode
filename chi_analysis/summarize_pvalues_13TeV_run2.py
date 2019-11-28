from ROOT import RooStats

import os
from ROOT import *
import array

only6000=False

massbinssets1=[[(6000,7000),(7000,13000)],
	      [(4800,5400),(5400,6000),(6000,7000),(7000,13000)],
	      [(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)],
	      ]
massbinssets2=[[(7000,13000)],
	      [(6000,7000)],
	      [(5400,6000)],
	      [(4800,5400)],
	      [(4200,4800)],
	      [(3600,4200)],
	      [(3000,3600)],
	      [(2400,3000)],
	      [(1900,2400)],
	      [(1500,1900)],
	      [(1200,1500)],
	      ]
if only6000:
  massbinssets1=[[(4800,5400),(5400,6000),(6000,13000)],
	      [(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,13000)],
	      ]
  massbinssets2=[[(6000,13000)],
	      [(5400,6000)],
	      [(4800,5400)],
	      [(4200,4800)],
	      [(3600,4200)],
	      [(3000,3600)],
	      [(2400,3000)],
	      [(1900,2400)],
	      [(1500,1900)],
	      [(1200,1500)],
	      ]

for signal,signalMass,massbinsset in [("CIplusLL","12000",massbinssets2),
                                   ("AntiCIplusLL","12000",massbinssets2),
				   ("cs_ct14nlo_","13000",massbinssets1), #nnlo
				   ("LHCa1108_DMAxial_Dijet_LO_Mphi_exp_5000_run2_v6","",[""]),
				   ]:
  for massbins in massbinsset:
    limits={}
    if "DM" in signal:
      name="pvalue"+signal
      if only6000:
        name=name.replace("v6","v5")
    else:
      name="pvalue_LHCa"+signal+"_"+("_".join([s[0:4] for s in str(massbins).strip("[]").split("(")])).strip("_")+"_exp_"+signalMass+"_run2"
    print name
    f1=open(name+".txt")
    observed=0
    for l in f1.readlines():
        if "CLb" in l:
	    pval=float(l.split("=")[1].split("+")[0].strip(" "))
            print "pvalue",pval,"significance",RooStats.PValueToSignificance((1.-pval)/2.)
        if "Significance:" in l:
	    significance=float(l.strip().split(" ")[-1].strip(")"))
            print "significance",significance
