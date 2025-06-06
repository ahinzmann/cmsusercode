import os,sys
from ROOT import *
import array
import ROOT
import subprocess

only6000=False
   
allmassbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]

models=[]
#models+=[3] #ADD
#models+=[10] #QBH ADD6
#models+=[11] #QBH RS1
#models+=[12] #QBH Blackmax n=6
#models+=[13] #QBH Blackmax n=6 RECO
models+=[60,61,62,63,64,65,66,67,68,69] #CI
#models+=[70,71,72,73,74,75,76,77] 
#models+=[78,79,80,81,82,83,84,85] 
#models+=[30,31,32,33,34,35,36,37,38,39,40] #pvalues (CI-LO)
#models+=[41,42,43,44] #pvalues (CI-NLO)
#models+=[45,46,47,48,49,50,51,52,53,54,55] #pvalues (Anti-CI-LO)
#models+=[47]
#models+=[88,89]
#models+=[60,61]
#models+=[90,91,92,93,94] #alp
#models+=[95,96,97,98,99] #tripleG
#>100 steered by loop script #DM, can be run with DMpvalue=True/False

VectorDM=False
AxialDM=True

injectSignal=False
dataWithSignal="_DMAxial_Dijet_LO_Mphi_4000_3000_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_1_chi_inject.root"

#jesSources=27 # 1 corresponds to the single overall variation, 27 UL (23EOY) to all
jesSources=30 # 1 corresponds to the single overall variation, 30 UL (decorrelated PileUpPtEC1 across years)
jerSources=1 # 1 corresponds to the single overall variation, 1 UL (6EOY) to all
uncorrelatedSimUncertainties=True
separateScaleUncertainties=False
alternateScaleUncertainty=True
uncorrelatedScaleUncertainties=False
scaleUncertainty=True
theoryStatUncertainties=True
useNNLO=True # choice for QCD
useM2=True # choice of mu-scale for QCD
runs="2" # "2" or "3" or "23"
run=runs[-1]
use_NNPDF3=False
use_CP2=False

isGen=False

isCB=False

isInjection=False

DMpvalue=False

withDiffNuisance=False

if use_NNPDF3:
  pdfset="nn31nnlo"
else:
  pdfset="ct14nnlo"
if not useNNLO:
  pdfset="ct14nlo"

signalName={}
signalExtraName={}
  
if VectorDM:
  counter=100
  #for mdm in ["1","3000"]:
  for mdm in ["1"]:
    for gv in ["0p01","0p05","0p1","0p2","0p25","0p3","0p5","0p75","1","1p5","2p0","2p5","3p0"]:
      #models+=[counter]
      signalName[counter]="DMVector_Dijet_LO_Mphi"
      signalExtraName[counter]="_"+mdm+"_1p0_1p0_Mar5_gdmv_1p0_gdma_0_gv_"+gv+"_ga_0"
      counter+=1

if AxialDM:
  counter=1100
  #for mdm in ["1","3000"]:
  for mdm in ["1"]:
    for ga in ["0p01","0p05","0p1","0p2","0p25","0p3","0p5","0p75","1","1p5","2p0","2p5","3p0"]:
      models+=[counter]
      signalName[counter]="DMAxial_Dijet_LO_Mphi"
      signalExtraName[counter]="_"+mdm+"_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_"+ga
      counter+=1    


testStat="LHC"# in 2012 and 2015 data used "LEP", checking "TEV" and "LHC" for 2016 data
asym="a" #"a" for asymptotic CLS
#testStat="LEP"
# asym=""
# The POI for LHC-style CLS is not clear, since CI models have no freedom  in signal strength or cross section.
# The LEP-style and TEV-style CLS do not fit the POI.

version="_v6b" #version number controls how many massbin to use for DM

if len(sys.argv)>1:
   models=[int(sys.argv[1])]

for model in models:

 signalExtra=""
 includeSignalTheoryUncertainties=True # Always yes, such that background scale affects both s+b and b

 if model<100:
    version=""

 massbins=[(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)] # did not calculate CI for lower mass bins yet
 if only6000:
   massbins=[(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,13000)] # did not calculate CI for lower mass bins yet
 
 if model==1:
    signal="CIplusLL"    
    signalMasses=[9000,10000,11000,12000,13000,14000,16000,18000]
 if model==2:
    signal="CIminusLL"    
    signalMasses=[9000,10000,11000,12000,13000,14000,16000,18000]
 if model==3:
    signal="ADD"
    signalMasses=[9000,10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000]
    #massbins=[(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    massbins=[(4800,5400),(5400,6000),(6000,7000),(7000,13000)] #signal MC statistics not enough in bins <4800
    if only6000:
      massbins=[(4800,5400),(5400,6000),(6000,13000)] #signal MC statistics not enough in bins <4800
    #includeSignalTheoryUncertainties=True # from QCD-only though
 if model==4:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_LL+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000]#8000,9000,
    massbins=[(4200,4800),(4800,13000)]
 if model==5:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_LL+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000]#8000,9000,
    massbins=[(4800,13000)]
 if model==6:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_LL+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000]#8000,9000,
    massbins=[(4200,4800)]
 if model==7:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_LL+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000]#8000,9000,
    massbins=[(3600,4200)]
 if model==8:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_LL-"
    signalMasses=[12000,13000,14000,15000,16000,17000,18000,19000,20000,22000,24000]
    massbins=[(4200,4800),(4800,13000)]
 if model==9:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_LL-"
    signalMasses=[12000,13000,14000,15000,16000,17000,18000,19000,20000,22000,24000]
    massbins=[(4800,13000)]
 if model==10:
    if use_NNPDF3:
      signal="QBH_ADD6_"    
      signalMasses=[7500,8000,8500,9000,9500,10000,11000,12000]
      massbins=[(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
      signalExtra=""
    else:
      signal="QBH_"    
      signalMasses=[7500,8000,8500,9000,9500,10000,10500,11000]
      massbins=[(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
      signalExtra="_6"
 if model==11:
    if use_NNPDF3:
      signal="QBH_RS1_"    
      signalMasses=[4500,5000,5500,6000,6500,7000,7500,8000]
      massbins=[(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
      signalExtra=""
    else:
      signal="QBH_"
      signalMasses=[4500,5000,5500,6000,6500,7000]
      #massbins=[(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
      massbins=[(5400,6000),(6000,7000),(7000,13000)]
      signalExtra="_RS1"
 if model==12:
    signal="QBH_MD"
    signalMasses=[6000,7000,8000,9000]
    massbins=[(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    #massbins=[(5400,6000),(6000,7000),(7000,13000)]
    signalExtra="_n6"
 if model==13:
    signal="QBH_MD"
    signalMasses=[6000,7000,8000,9000]
    massbins=[(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    #massbins=[(5400,6000),(6000,7000),(7000,13000)]
    signalExtra="_n6"

 if model==18:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000]#8000,9000,
 if model==19:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000]#8000,9000,
    includeSignalTheoryUncertainties=True

 if model==20:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_LL+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000] #dropped 8000, 9000
 if model==21:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_LL-"
    signalMasses=[12000,13000,14000,15000,16000,17000,18000,19000,20000,22000,24000]
 if model==22:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_RR+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000] # dropped 8000, 9000
 if model==23:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_RR-"
    signalMasses=[12000,13000,14000,15000,16000,17000,18000,19000,20000,22000,24000]
 if model==24:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_VV+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000]
 if model==25:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_VV-"
    signalMasses=[13000,14000,15000,16000,17000,18000,19000,20000,22000,24000,26000]
 if model==26:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_AA+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000]
 if model==27:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_AA-"
    signalMasses=[13000,14000,15000,16000,17000,18000,19000,20000,22000,24000,26000]
 if model==28:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_V-A+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000] # dropped 8000, 9000
 if model==29:
    signal="cs_"+pdfset+"_0_"
    signalExtra="_V-A-"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000]#8000,9000,

 if model>=30 and model<50:
    includeSignalTheoryUncertainties=True

 if model==30:
    signal="CIplusLL"
    signalMasses=[12000]
    massbins=[(7000,13000),]
 if model==31:
    signal="CIplusLL"
    signalMasses=[12000]
    massbins=[(6000,7000),]
    if only6000:
      massbins=[(6000,13000),]
 if model==32:
    signal="CIplusLL"
    signalMasses=[12000]
    massbins=[(5400,6000),]
 if model==33:
    signal="CIplusLL"
    signalMasses=[12000]
    massbins=[(4800,5400),]
 if model==34:
    signal="CIplusLL"
    signalMasses=[12000]
    massbins=[(4200,4800),]
 if model==35:
    signal="CIplusLL"
    signalMasses=[12000]
    massbins=[(3600,4200),]
 if model==36:
    signal="CIplusLL"
    signalMasses=[12000]
    massbins=[(3000,3600),]
 if model==37:
    signal="CIplusLL"
    signalMasses=[12000]
    massbins=[(2400,3000),]
 if model==38:
    signal="CIplusLL"
    signalMasses=[12000]
    massbins=[(1900,2400),]
 if model==39:
    signal="CIplusLL"
    signalMasses=[12000]
    massbins=[(1500,1900),]
 if model==40:
    signal="CIplusLL"
    signalMasses=[12000]
    massbins=[(1200,1500),]
 if model==41:
    signal="cs_"+pdfset+"_" #nlo/nnlo
    signalExtra="_LL+"
    signalMasses=[13000]
    massbins=[(6000,7000),(7000,13000)]
 if model==42:
    signal="cs_"+pdfset+"_" #nlo/nnlo
    signalExtra="_LL+"
    signalMasses=[13000]
    massbins=[(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    if only6000:
      massbins=[(4800,5400),(5400,6000),(6000,13000)]
 if model==43:
    signal="cs_"+pdfset+"_" #nlo/nnlo
    signalExtra="_LL+"
    signalMasses=[13000]
    massbins=[(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    if only6000:
      massbins=[(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,13000)]
 if model==44:
    signal="cs_"+pdfset+"_" #nlo/nnlo
    signalExtra="_LL+"
    signalMasses=[13000]
    massbins=[(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    if only6000:
      massbins=[(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,13000)]

 if model==45:
    signal="AntiCIplusLL"    
    signalMasses=[12000]
    massbins=[(7000,13000),]
 if model==46:
    signal="AntiCIplusLL"    
    signalMasses=[12000]
    massbins=[(6000,7000),]
    if only6000:
      massbins=[(6000,13000),]
 if model==47:
    signal="AntiCIplusLL"    
    signalMasses=[12000]
    massbins=[(5400,6000),]
 if model==48:
    signal="AntiCIplusLL"    
    signalMasses=[12000]
    massbins=[(4800,5400),]
 if model==49:
    signal="AntiCIplusLL"    
    signalMasses=[12000]
    massbins=[(4200,4800),]
 if model==50:
    signal="AntiCIplusLL"    
    signalMasses=[12000]
    massbins=[(3600,4200),]
 if model==51:
    signal="AntiCIplusLL"    
    signalMasses=[12000]
    massbins=[(3000,3600),]
 if model==52:
    signal="AntiCIplusLL"    
    signalMasses=[12000]
    massbins=[(2400,3000),]
 if model==53:
    signal="AntiCIplusLL"    
    signalMasses=[12000]
    massbins=[(1900,2400),]
 if model==54:
    signal="AntiCIplusLL"    
    signalMasses=[12000]
    massbins=[(1500,1900),]
 if model==55:
    signal="AntiCIplusLL"    
    signalMasses=[12000]
    massbins=[(1200,1500),]

 if model>=60 and model<100:
    includeSignalTheoryUncertainties=True

 if model==60 or model==88:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,22000,24000,26000,28000,30000] #21000
 if model==61 or model==89:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL-"
    signalMasses=[14000,15000,16000,17000,18000,19000,20000,22000,24000,26000,28000,30000] # 25000
 if model==62:
    signal="cs_"+pdfset+"_"
    signalExtra="_RR+"
    signalMasses=[11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,22000,24000,26000,28000,30000] # 21000
 if model==63:
    signal="cs_"+pdfset+"_"
    signalExtra="_RR-"
    signalMasses=[16000,17000,18000,19000,20000,22000,24000,26000,28000] # 25000,27000
 if model==64:
    signal="cs_"+pdfset+"_"
    signalExtra="_VV+"
    signalMasses=[14000,15000,16000,17000,18000,19000,20000,22000,24000,26000,28000,30000] #21000,23000
 if model==65:
    signal="cs_"+pdfset+"_"
    signalExtra="_VV-"
    signalMasses=[20000,22000,24000,26000,28000,30000]#,32000,34000,36000, 27000, 29000
 if model==66:
    signal="cs_"+pdfset+"_"
    signalExtra="_AA+"
    signalMasses=[14000,15000,16000,17000,18000,19000,20000,22000,24000,26000,28000,30000] # 21000, 23000
 if model==67:
    signal="cs_"+pdfset+"_"
    signalExtra="_AA-"
    signalMasses=[20000,22000,24000,26000,28000,30000]#,32000,34000,36000 , 27000,29000
 if model==68:
    signal="cs_"+pdfset+"_"
    signalExtra="_V-A+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,22000,24000,26000] # dropped 7000,8000,9000
 if model==69:
    signal="cs_"+pdfset+"_"
    signalExtra="_V-A-"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,22000,24000,26000] # dropped 7000,8000,9000

 if model==70:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[10000,11000,12000,13000,14000,15000,16000,17000,18000,19000]
    massbins=[(6000,13000)]
 if model==71:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[11000,12000,13000,14000,15000,16000,17000,18000,19000,20000]
    massbins=[(5400,6000)]
 if model==72:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[11000,12000,13000,14000,15000,16000,17000,18000,19000,20000]
    massbins=[(4800,5400)]
 if model==73:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000]#23000
    massbins=[(5400,6000),(6000,7000),(7000,13000)]
 if model==74:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL-"
    signalMasses=[13000,14000,15000,16000,17000,18000,19000,20000,22000]
    massbins=[(6000,13000)]
 if model==75:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL-"
    signalMasses=[15000,16000,17000,18000,19000,20000,22000,24000]
    massbins=[(5400,6000)]
 if model==76:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL-"
    signalMasses=[14000,15000,16000,17000,18000,19000,20000,22000,24000]
    massbins=[(4800,5400)]
 if model==77:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL-"
    signalMasses=[14000,15000,16000,17000,18000,19000,20000,22000,24000]
    massbins=[(5400,6000),(6000,7000),(7000,13000)]
 if model==78:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[12000,13000,14000,15000,16000,17000,18000,19000,20000,22000]#21000,23000
    massbins=[(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
 if model==79:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL-"
    signalMasses=[14000,15000,16000,17000,18000,19000,20000,22000,24000,26000,28000,30000]
    massbins=[(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
 if model==80:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[12000,13000,14000,15000,16000,17000,18000,19000,20000,22000]#21000,23000
    massbins=[(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
 if model==81:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL-"
    signalMasses=[14000,15000,16000,17000,18000,19000,20000,22000,24000,26000,28000,30000]
    massbins=[(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
 if model==82:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[12000,13000,14000,15000,16000,17000,18000,19000,20000,22000]#21000,23000
    massbins=[(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
 if model==83:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL-"
    signalMasses=[14000,15000,16000,17000,18000,19000,20000,22000,24000,26000,28000,30000]
    massbins=[(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
 if model==84:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[12000,13000,14000,15000,16000,17000,18000,19000,20000,22000]#21000,23000
    massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
 if model==85:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL-"
    signalMasses=[14000,15000,16000,17000,18000,19000,20000,22000,24000,26000,28000,30000]
    massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
 if model==86:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL+"
    signalMasses=[12000,13000,14000,15000,16000,17000,18000,19000,20000,22000]#21000,23000
    massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
 if model==87:
    signal="cs_"+pdfset+"_"
    signalExtra="_LL-"
    signalMasses=[14000,15000,16000,17000,18000,19000,20000,22000,24000,26000,28000,30000]
    massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]

 if model==90:
    signal="alp_QCD_fa"
    signalMasses=[1000,1500,2000,2500,3000,3500,4000,4500,4500,5000]
    massbins=[(2400,3000)]
 if model==91:
    signal="alp_QCD_fa"
    signalMasses=[1000,1500,2000,2500,3000,3500,4000,4500,4500,5000]
    massbins=[(2400,3000),(3000,3600)]
 if model==92:
    signal="alp_QCD_fa"
    signalMasses=[1000,1500,2000,2500,3000,3500,4000,4500,4500,5000]
    massbins=[(2400,3000),(3000,3600),(3600,4200)]
 if model==93:
    signal="alp_QCD_fa"
    signalMasses=[1000,1500,2000,2500,3000,3500,4000,4500,4500,5000]
    massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800)]
 if model==94:
    signal="alp_QCD_fa"
    signalMasses=[1000,1500,2000,2500,3000,3500,4000,4500,4500,5000]
    massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400)]
 if model==95:
    signal="tripleG_QCD_CG"
    signalMasses=[0.0025,0.005,0.0075,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.1]
    massbins=[(2400,3000)]
 if model==96:
    signal="tripleG_QCD_CG"
    signalMasses=[0.0025,0.005,0.0075,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.1]
    massbins=[(2400,3000),(3000,3600)]
 if model==97:
    signal="tripleG_QCD_CG"
    signalMasses=[0.0025,0.005,0.0075,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.1]
    massbins=[(2400,3000),(3000,3600),(3600,4200)]
 if model==98:
    signal="tripleG_QCD_CG"
    signalMasses=[0.0025,0.005,0.0075,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.1]
    massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800)]
 if model==99:
    signal="tripleG_QCD_CG"
    signalMasses=[0.0025,0.005,0.0075,0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.1]
    massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400)]

 dire="/data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/"
 if run=="2":
   if use_NNPDF3:
     prefix=dire+"versions/run2ULNNLO_m2_NNPDF3/datacard_shapelimit13TeV"
   elif useM2:
     prefix=dire+"versions/run2ULNNLO_m2/datacard_shapelimit13TeV"
   elif useNNLO:
     prefix=dire+"versions/run2ULNNLO_pt12/datacard_shapelimit13TeV"
   else:
     prefix=dire+"versions/run2ULNLO_pt12/datacard_shapelimit13TeV"
 elif run=="3":
   prefix=dire+"versions/run3NNLO_m2/datacard_shapelimit13TeV"
 else:
   whatprefix

 if model>=30 and model<60:
    name="pvalue_"+testStat+asym+signal+"_"+(("2400_3000" if model==44 else "")+"_".join([s[0:4] for s in str(massbins).strip("[]").split("(")])).strip("_")
 elif model>=100:  # Dark Matter
    signal=signalName[model]
    signalExtra=signalExtraName[model]
    #if float(signalExtraName[model].split("_")[2])<=2:
    #signalMasses=[1000,1500,1750,2000,2250,2500,3000,3500,4000,4500,5000,6000,7000]
    #signalMasses=[2000,2250,2500,3000,3500,4000,4500,5000,6000,7000]
    signalMasses=[4000,4500,5000,6000,7000] # for usual limits
    #signalMasses=[7000] # for postfit plots
    #signalMasses=[2000] # for control region check
    #signalMasses=[5000] # for low bins check
    includeSignalTheoryUncertainties=True # Assign QCD-only scale uncertainty to QCD+DM
    
    if isGen:
        dire="/uscms_data/d3/jingyu/ChiAnalysis/DMlimits/CMSSW_8_0_25/src/cmsusercode/chi_analysis/"
        prefix="/uscms_data/d3/jingyu/ChiAnalysis/DMlimits/CMSSW_8_0_25/src/cmsusercode/chi_analysis/invertMatrixOct20/datacard_shapelimit13TeV"
    elif isCB:
        dire="/uscms_data/d3/jingyu/ChiAnalysis/DMlimits/CMSSW_9_2_4/src/cmsusercode/chi_analysis/"
        prefix="/uscms_data/d3/jingyu/ChiAnalysis/DMlimits/CMSSW_9_2_4/src/cmsusercode/chi_analysis/crystalBallSmearedAug30/datacard_shapelimit13TeV"
    if isInjection:
        prefix=prefix.replace("/datacard_","Injection0p75/datacard_")
    print(prefix)

    name="limits"+testStat+asym+str(model)+"_"+signal
     
    if isGen:
        name=name.replace("limits","limitsGen")
    elif isCB:
        name=name.replace("limits","limitsDetCB")
    else:
        name=name.replace("limits","limitsDet")
    if DMpvalue:
        name=name.replace("limitsDet","pvalue")
        
    if isInjection:
        name=name+"Injection0p75"
 else:
    name="limits"+testStat+asym+str(model)+"_"+signal  
 print(name)

 limits={}
 for signalMass in signalMasses:
    signalWithMass=(signal+str(signalMass)+signalExtra).replace("MD"+str(signalMass),"MD"+str(signalMass)+"_MBH"+str(signalMass+1000))
    print(signalWithMass)

    if signalWithMass=="CIplusLL8000":
            fname=prefix + '_GENnp-0-run'+run+'_chi.root'
    elif signalWithMass=="CIplusLL9000":
            fname=prefix + '_GENnp-1-run'+run+'_chi.root'
    elif signalWithMass=="CIplusLL10000":
            fname=prefix + '_GENnp-2-run'+run+'_chi.root'
    elif signalWithMass=="CIplusLL11000":
            fname=prefix + '_GENnp-3-run'+run+'_chi.root'
    elif signalWithMass=="CIplusLL12000":
            fname=prefix + '_GENnp-4-run'+run+'_chi.root'
    elif signalWithMass=="CIplusLL13000":
            fname=prefix + '_GENnp-5-run'+run+'_chi.root'
    elif signalWithMass=="CIplusLL14000":
            fname=prefix + '_GENnp-6-run'+run+'_chi.root'
    elif signalWithMass=="CIplusLL16000":
            fname=prefix + '_GENnp-7-run'+run+'_chi.root'
    elif signalWithMass=="CIplusLL18000":
            fname=prefix + '_GENnp-8-run'+run+'_chi.root'
    elif signalWithMass=="CIminusLL8000":
            fname=prefix + '_GENnp-9-run'+run+'_chi.root'
    elif signalWithMass=="CIminusLL9000":
            fname=prefix + '_GENnp-10-run'+run+'_chi.root'
    elif signalWithMass=="CIminusLL10000":
            fname=prefix + '_GENnp-11-run'+run+'_chi.root'
    elif signalWithMass=="CIminusLL11000":
            fname=prefix + '_GENnp-12-run'+run+'_chi.root'
    elif signalWithMass=="CIminusLL12000":
            fname=prefix + '_GENnp-13-run'+run+'_chi.root'
    elif signalWithMass=="CIminusLL13000":
            fname=prefix + '_GENnp-14-run'+run+'_chi.root'
    elif signalWithMass=="CIminusLL14000":
            fname=prefix + '_GENnp-15-run'+run+'_chi.root'
    elif signalWithMass=="CIminusLL16000":
            fname=prefix + '_GENnp-16-run'+run+'_chi.root'
    elif signalWithMass=="CIminusLL18000":
            fname=prefix + '_GENnp-17-run'+run+'_chi.root'
    elif model==3 and run=="3":
        fname=prefix + "_GENnp-run-"+run+"-" + signalWithMass + '_chi.root'
    elif (model==3 or model==10 or model==11) and use_NNPDF3:
        fname=prefix + "_GENnp-"+("CP2" if use_CP2 else "CP5")+"-run"+run+"-" + signalWithMass + '_chi.root'
    elif signalWithMass=="ADD6000":
        fname=prefix + '_GENnp-18-run'+run+'_chi.root'
    elif signalWithMass=="ADD7000":
        fname=prefix + '_GENnp-19-run'+run+'_chi.root'
    elif signalWithMass=="ADD8000":
        fname=prefix + '_GENnp-20-run'+run+'_chi.root'
    elif signalWithMass=="ADD9000":
        fname=prefix + '_GENnp-21-run'+run+'_chi.root'
    elif signalWithMass=="ADD10000":
        fname=prefix + '_GENnp-22-run'+run+'_chi.root'
    elif signalWithMass=="ADD11000":
        fname=prefix + '_GENnp-23-run'+run+'_chi.root'
    elif signalWithMass=="ADD12000":
        fname=prefix + '_GENnp-24-run'+run+'_chi.root'
    elif signalWithMass=="ADD13000":
        fname=prefix + '_GENnp-25-run'+run+'_chi.root'
    elif signalWithMass=="ADD14000":
        fname=prefix + '_GENnp-26-run'+run+'_chi.root'
    elif signalWithMass=="ADD15000":
        fname=prefix + '_GENnp-27-run'+run+'_chi.root'
    elif signalWithMass=="ADD16000":
        fname=prefix + '_GENnp-28-run'+run+'_chi.root'
    elif signalWithMass=="ADD17000":
        fname=prefix + '_GENnp-29-run'+run+'_chi.root'
    elif signalWithMass=="ADD18000":
        fname=prefix + '_GENnp-30-run'+run+'_chi.root'
    elif signalWithMass=="ADD19000":
        fname=prefix + '_GENnp-31-run'+run+'_chi.root'
    elif signalWithMass=="ADD20000":
        fname=prefix + '_GENnp-32-run'+run+'_chi.root'
    elif signalWithMass=="ADD21000":
        fname=prefix + '_GENnp-33-run'+run+'_chi.root'
    elif signalWithMass=="ADD22000":
        fname=prefix + '_GENnp-34-run'+run+'_chi.root'
    elif signalWithMass=="AntiCIplusLL12000":
        fname=prefix + '_GENnp-antici-run'+run+'_chi.root'
    elif signalWithMass=="QBH_"+str(signalMass)+"_6":
        fname=prefix + "_" + signalWithMass + "-run"+run+"_chi.root"
    elif signalWithMass=="QBH_"+str(signalMass)+"_RS1":
        fname=prefix + "_" + signalWithMass + "-run"+run+"_chi.root"
    elif signalWithMass=="QBH_MD"+str(signalMass)+"_MBH"+str(signalMass+1000)+"_n6" and model==13:
        fname=prefix + "RECO_run2_UL18_" + signalWithMass + "-GEN_chi.root"
    elif signalWithMass=="QBH_MD"+str(signalMass)+"_MBH"+str(signalMass+1000)+"_n6":
        fname=prefix + "_run2_UL18_" + signalWithMass + "-GEN_chi.root"
    elif "cs" in signal:
        fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
    elif "alp" in signal or "tripleG" in signal:
      signalWithMass=signal+str(signalMass).replace(".","p")+signalExtra
      fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
    elif "DM" in signal and version=="_v1":
      signalWithMass=signal+'_'+str(signalMass)+signalExtra
      fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
      fname=fname.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      signalWithMass=signalWithMass.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      if signalMass<=2000:
        massbins=[(2400,3000)]
      elif signalMass==2500:
        massbins=[(2400,3000)]
      elif signalMass==3000:
        massbins=[(2400,3000),(3000,3600)]
      elif signalMass==3500:
        massbins=[(3000,3600),(3600,4200)]
      elif signalMass==4000:
        massbins=[(3600,4200),(4200,4800)]
      elif signalMass==4500:
        massbins=[(4200,4800),(4800,5400)]
      elif signalMass==5000:
        massbins=[(4800,5400),(5400,6000)]
      elif signalMass>=6000:
        massbins=[(5400,6000),(6000,7000),(7000,13000)]
    elif "DM" in signal and version=="_v2":
      signalWithMass=signal+'_'+str(signalMass)+signalExtra
      fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
      fname=fname.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000") # FIX 7000 only available with large DM mass
      signalWithMass=signalWithMass.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000") # FIX 7000 only available with large DM mass
      massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    elif "DM" in signal and version=="_v3":
      signalWithMass=signal+'_'+str(signalMass)+signalExtra
      fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
      fname=fname.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      signalWithMass=signalWithMass.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      if signalMass<=1500:
        massbins=[(2400,3000)]
      elif signalMass<=2000:
        massbins=[(2400,3000),(3000,3600)]
      elif signalMass<=2500:
        massbins=[(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=3000:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000)]
      elif signalMass<=6000:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000)]
      else:
        massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    elif "DM" in signal and version=="_v4":
      signalWithMass=signal+'_'+str(signalMass)+signalExtra
      fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
      fname=fname.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      signalWithMass=signalWithMass.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      if signalMass<=1500:
        massbins=[(2400,3000)]
      elif signalMass<=2000:
        massbins=[(2400,3000)]
      elif signalMass<=2500:
        massbins=[(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=3000:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800)]
      elif signalMass<=3500:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400)]
      elif signalMass<=4000:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000)]
      elif signalMass<=4500:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
      elif signalMass<=5000:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
      elif signalMass<=6000:
        massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000)]
      else:
        massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    elif "DM" in signal and version=="_v5":
      signalWithMass=signal+'_'+str(signalMass)+signalExtra
      fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
      fname=fname.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      signalWithMass=signalWithMass.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      if signalMass<=2500:
        massbins=[(2400,3000)]
      elif signalMass<=3000:
        massbins=[(2400,3000),(3000,3600)]
      elif signalMass<=3500:
        massbins=[(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4000:
        massbins=[(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4500:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800)]
      elif signalMass<=5000:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000)]
      elif signalMass<=6000:
        massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,13000)]
      else:
        massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,13000)]
    elif "DM" in signal and version=="_v6":
      signalWithMass=signal+'_'+str(signalMass)+signalExtra
      fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
      fname=fname.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      signalWithMass=signalWithMass.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      if signalMass<=2500:
        massbins=[(2400,3000)]
      elif signalMass<=3000:
        massbins=[(2400,3000),(3000,3600)]
      elif signalMass<=3500:
        massbins=[(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4000:
        massbins=[(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4500:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800)]
      elif signalMass<=5000:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000)]
      elif signalMass<=6000:
        massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000)]
      else:
        massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    elif "DM" in signal and version=="_v6b":
      signalWithMass=signal+'_'+str(signalMass)+signalExtra
      fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
      fname=fname.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      signalWithMass=signalWithMass.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      if signalMass<=2500:
        massbins=[(2400,3000)]
      elif signalMass<=3000:
        massbins=[(2400,3000),(3000,3600)]
      elif signalMass<=3500:
        massbins=[(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4000:
        massbins=[(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4500:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800)]
      elif signalMass<=5000:
        #massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000)] # for usual limits
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400)] # for lowestBinsCheck
      elif signalMass<=6000:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000)]
      else:
        massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    elif "DM" in signal and version=="_v7":
      signalWithMass=signal+'_'+str(signalMass)+signalExtra
      fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
      fname=fname.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      signalWithMass=signalWithMass.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      if signalMass<=2500:
        massbins=[(1900,2400),(2400,3000)]
      elif signalMass<=3000:
        massbins=[(1900,2400),(2400,3000),(3000,3600)]
      elif signalMass<=3500:
        massbins=[(1900,2400),(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4000:
        massbins=[(1900,2400),(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4500:
        massbins=[(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800)]
      elif signalMass<=5000:
        massbins=[(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000)]
      elif signalMass<=6000:
        massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000)]
      else:
        massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    elif "DM" in signal and version=="_v8":
      signalWithMass=signal+'_'+str(signalMass)+signalExtra
      fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
      fname=fname.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      signalWithMass=signalWithMass.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      if signalMass<=2500:
        massbins=[(1500,1900),(1900,2400),(2400,3000)]
      elif signalMass<=3000:
        massbins=[(1500,1900),(1900,2400),(2400,3000),(3000,3600)]
      elif signalMass<=3500:
        massbins=[(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4000:
        massbins=[(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4500:
        massbins=[(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800)]
      elif signalMass<=5000:
        massbins=[(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000)]
      elif signalMass<=6000:
        massbins=[(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000)]
      else:
        massbins=[(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    elif "DM" in signal and version=="_v9":
      signalWithMass=signal+'_'+str(signalMass)+signalExtra
      fname=prefix+"_"+str(signalWithMass)+"-run"+run+"_chi.root"
      fname=fname.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      signalWithMass=signalWithMass.replace("6000_3000","6000_2990").replace("7000_3000","7000_4000").replace("8000_3000","8000_3990").replace("7000_1","7000_4000")
      if signalMass<=2500:
        massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000)]
      elif signalMass<=3000:
        massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600)]
      elif signalMass<=3500:
        massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4000:
        massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200)]
      elif signalMass<=4500:
        massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800)]
      elif signalMass<=5000:
        massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000)]
      elif signalMass<=6000:
        massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000)]
      else:
        massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    print(fname)
    if not "DM" in signal and not "cs" in signal and not "QBH" in signal and not "alp" in signal and not "tripleG" in signal:
        signalWithMass="QCD"+signalWithMass

    # Use all lower mass bins in fit with bg-only hypothesis where no signal is present
    signalmassbins=massbins[:]
    if "limit" in name or model==44:
      for massbin in allmassbins:
        if not massbin in massbins:
          massbins.insert(allmassbins.index(massbin),massbin)
        else:
          break
    #massbins=[(5400,6000),(6000,7000),(7000,13000)] # for tests
    #signalmassbins=[(5400,6000),(6000,7000),(7000,13000)] # for tests
    print("massbins", massbins)
    print("signalmassbins", signalmassbins)
    
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

    def system_call(command):
       print(command)
       p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
       return p.stdout.read()

    def remove_prefix(text, prefix):
       if text.startswith(prefix):
           return text[len(prefix):]
       return text  # or whatever

    f=TFile(fname)
    cfg=open("chi_datacard13TeV"+str(model)+"_"+remove_prefix(signalWithMass,"QCD")+"_run"+run+".txt","w")
    cfg.writelines("""
imax """+str(len(massbins))+""" number of channels
jmax 2 number of backgrounds
kmax """+str(3+3+3*(len(massbins)-1)*uncorrelatedSimUncertainties+jesSources+jerSources+1*scaleUncertainty+1*alternateScaleUncertainty+1*separateScaleUncertainties+(len(massbins)-1)*uncorrelatedScaleUncertainties+len(statUncertainties))+""" number of nuisance parameters""")
    cfg.writelines("""
-----------
""")
    for i in range(len(massbins)):
        if injectSignal:
          cfg.writelines("""shapes data_obs bin"""+str(i)+""" """+prefix+dataWithSignal+""" data_obs#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1 data_obs#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1_$SYSTEMATIC
""")   
        else:
          cfg.writelines("""shapes data_obs bin"""+str(i)+""" """+fname+""" data_obs#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1 data_obs#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1_$SYSTEMATIC
""")
        if not massbins[i] in signalmassbins:
          cfg.writelines("""shapes """+signalWithMass+""" bin"""+str(i)+""" """+fname+""" """+signalWithMass+"""_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1 """+signalWithMass+"""_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1_$SYSTEMATIC
""")
        else:
          cfg.writelines("""shapes """+signalWithMass+""" bin"""+str(i)+""" """+fname+""" """+signalWithMass+"""#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1 """+signalWithMass+"""#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1_$SYSTEMATIC
""")
        cfg.writelines("""shapes """+signalWithMass+"""_ALT bin"""+str(i)+""" """+fname+""" """+signalWithMass+"""_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1 """+signalWithMass+"""_ALT#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1_$SYSTEMATIC
""")
        cfg.writelines("""shapes QCD bin"""+str(i)+""" """+fname+""" QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1 QCD#chi"""+str(massbins[i][0])+"""_"""+str(massbins[i][1])+"""_rebin1_$SYSTEMATIC
""")
    cfg.writelines("""
-----------
""")
    text="bin "
    for i in range(len(massbins)):
       text+="bin"+str(i)+" "
    text+="\nobservation "
    print("data_obs#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
    for i in range(len(massbins)):
       hData=f.Get("data_obs#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
       text+=str(hData.Integral())+" "
    cfg.writelines(text+"""
-----------
""")
    text="bin "
    for i in range(len(massbins)):
       text+="bin"+str(i)+" bin"+str(i)+" bin"+str(i)+" "
    text+="\nprocess "
    for i in range(len(massbins)):
       text+=signalWithMass+" "+signalWithMass+"_ALT QCD "
    text+="\nprocess "
    for i in range(len(massbins)):
       text+="-1 0 1 "
    text+="\nrate "
    for i in range(len(massbins)):
       hQCD=f.Get("QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
       hALT=f.Get(signalWithMass+"_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
       if massbins[i] in signalmassbins:
         h=f.Get(signalWithMass+"#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
       else:
         h=hALT
       print("QCD#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1",hQCD.Integral())
       print(signalWithMass+"_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1",hALT.Integral())
       print(signalWithMass+"#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1",h.Integral())
       text+=str(h.Integral())+" "+str(hALT.Integral())+" "+str(hQCD.Integral())+" "
    cfg.writelines(text+"""
-----------
""")
    text=""
    ones=""
    for i in range(len(massbins)):
      ones+="1 1 - "
    text+="\nmodel shape "+ones
    if runs=="23": text+="\nnuisance edit rename * * model run"+run+"_model"
    text+="\nJERtail shape "+ones
    if runs=="23": text+="\nnuisance edit rename * * JERtail run"+run+"_JERtail"
    text+="\nsim shape "+ones
    if runs=="23": text+="\nnuisance edit rename * * sim run"+run+"_sim"
    if uncorrelatedSimUncertainties:
     for mn in massbins[1:]:
      text+="\nmodel"+str(mn[0])+" shape "+ones
      if runs=="23": text+="\nnuisance edit rename * * model"+str(mn[0])+" run"+run+"_model"+str(mn[0])
     for mn in massbins[1:]:
      text+="\nJERtail"+str(mn[0])+" shape "+ones
      if runs=="23": text+="\nnuisance edit rename * * JERtail"+str(mn[0])+" run"+run+"_JERtail"+str(mn[0])
     for mn in massbins[1:]:
      text+="\nsim"+str(mn[0])+" shape "+ones
      if runs=="23": text+="\nnuisance edit rename * * sim"+str(mn[0])+" run"+run+"_sim"+str(mn[0])
    if jesSources>1:
     for n in range(jesSources):
      text+="\njes"+str(n+1)+" shape "
      for i in range(len(massbins)):
         text+="1 1 - "
      if runs=="23" and n>1: text+="\nnuisance edit rename * * jes"+str(n+1)+" run"+run+"_jes"+str(n+1) # if n<2 combine would do the replace twice for n>10 and >20
    else:
      text+="\njes shape "
      for i in range(len(massbins)):
         text+="1 1 - "
      if runs=="23": text+="\nnuisance edit rename * * jes run"+run+"_jes"
    if jerSources>1:
     for n in range(jerSources):
      text+="\njer"+str(n+1)+" shape "
      for i in range(len(massbins)):
         text+="1 1 - "
      if runs=="23" and n>1: text+="\nnuisance edit rename * * jer"+str(n+1)+" run"+run+"_jer"+str(n+1) # if n<2 combine would do the replace twice for n>10 and >20
    else:
      text+="\njer shape "
      for i in range(len(massbins)):
         text+="1 1 - "
      if runs=="23": text+="\nnuisance edit rename * * jer run"+run+"_jer"
    text+="\nprefire shape "
    for i in range(len(massbins)):
        text+="1 1 - "
    if runs=="23": text+="\nnuisance edit rename * * prefire run"+run+"_prefire"
    text+="\ntrigger shape "
    for i in range(len(massbins)):
        text+="1 1 - "
    if runs=="23": text+="\nnuisance edit rename * * trigger run"+run+"_trigger"
    text+="\npdf shape "
    for i in range(len(massbins)):
      if includeSignalTheoryUncertainties:
       text+="1 1 - "
      else:
       text+="- 1 - "
    if separateScaleUncertainties:
      text+="\nscaleMuR shape "
      for i in range(len(massbins)):
        if includeSignalTheoryUncertainties:
         text+="1 1 - "
        else:
         text+="- 1 - "
      text+="\nscaleMuF shape "
      for i in range(len(massbins)):
        if includeSignalTheoryUncertainties:
         text+="1 1 - "
        else:
         text+="- 1 - "
    elif uncorrelatedScaleUncertainties:
     for mn in massbins:
      text+="\nscale"+str(mn[0])+" shape "
      for i in range(len(massbins)):
        if includeSignalTheoryUncertainties:
         text+="1 1 - "
        else:
         text+="- 1 - "
      if runs=="23": text+="\nnuisance edit rename * * scale"+str(mn[0])+" run"+run+"_scale"+str(mn[0])
    else:
      if alternateScaleUncertainty:
       text+="\nscaleAlt shape "
       for i in range(len(massbins)):
        if includeSignalTheoryUncertainties:
         text+="1 1 - "
        else:
         text+="- 1 - "
      if scaleUncertainty:
       text+="\nscale shape "
       for i in range(len(massbins)):
        if includeSignalTheoryUncertainties:
         text+="1 1 - "
        else:
         text+="- 1 - "
    for su in statUncertainties:
      text+="\n"+su+" shape "
      for i in range(len(massbins)):
         text+="1 1 - "
    cfg.writelines(text+"""
-----------
""")

    cfg.close()

    massbins=signalmassbins

    def remove_prefix(text, prefix):
       if text.startswith(prefix):
           return text[len(prefix):]
       return text  # or whatever

    #out=system_call("cp "+dire+"HiggsJPCmodified.py ${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/python")
    if runs=="23":
      out=system_call("combineCards.py run2=chi_datacard13TeV"+str(model)+"_"+remove_prefix(signalWithMass,"QCD")+"_run2.txt run3=chi_datacard13TeV"+str(model)+"_"+remove_prefix(signalWithMass,"QCD")+"_run3.txt > chi_datacard13TeV"+str(model)+"_"+remove_prefix(signalWithMass,"QCD")+"_run23.txt")
    out=system_call("text2workspace.py -m "+str(signalMass)+" chi_datacard13TeV"+str(model)+"_"+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".txt -P HiggsAnalysis.CombinedLimit.HiggsJPCmodified:twoHypothesisHiggs -o fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root")
    
    if testStat=="LEP":
     poi=""
     add=""
     method=testStat+" --singlePoint 1.0"
     ntoys=30000
    if testStat=="LHC":
     poi=" --redefineSignalPOIs x" # -H ProfileLikelihood
     method=testStat+" --frequentist"
     add="LHC"
     ntoys=1000
    if testStat=="TEV":
     poi=" --redefineSignalPOIs x"
     add=""
     method=testStat+" --singlePoint 1.0"
     ntoys=30000
    
    def remove_prefix(text, prefix):
       if text.startswith(prefix):
           return text[len(prefix):]
       return text  # or whatever

    if asym:
     if "limit" in name:
      #out=system_call("combine -m "+str(signalMass)+" -M Asymptotic -n "+signal+signalExtra+" fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root")
      out=system_call("combine -m "+str(signalMass)+" -M AsymptoticLimits -n "+signal+signalExtra+" fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root")
      # -H ProfileLikelihood
      f = open(name+"_exp_"+str(signalMass)+"_run"+runs+version+".txt","w");f.write(str(out));f.close()
     else: 
      out=system_call("combine -m "+str(signalMass)+" -M Significance -n "+signal+signalExtra+" fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root")
      #out=system_call("combine --signif -m "+str(signalMass)+" -M ProfileLikelihood -n "+signal+signalExtra+" fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root")
      f = open(name+"_exp_"+str(signalMass)+"_run"+runs+version+".txt","w");f.write(str(out));f.close()

    elif testStat=="LHC":
     if "limit" in name:
      for point in [0.1,0.2,0.4,0.6,0.8,1.0,1.3,1.6,2.0,5.0,10.0]:
       out=system_call("combine -m "+str(signalMass)+" -M HybridNew --rule CLs --saveHybridResult --singlePoint "+str(point)+" -s 10000"+str(int(point*100))+" --saveToys --testStat "+method+poi+" --fork 4 -T "+str(ntoys)+" -i 2 --clsAcc 0 -n "+signal+signalExtra+" fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root")
       f = open(name+"_"+str(signalMass)+str(point)+"_run"+runs+version+".txt","w");f.write(out);f.close()
      system_call("hadd -f grid_mX"+str(signalMass)+"_"+signal+signalExtra+".root higgsCombine"+signal+signalExtra+".HybridNew.mH"+str(signalMass)+".10000*.root")

      system_call("combine -M HybridNew --frequentist --grid grid_mX"+str(signalMass)+"_"+signal+signalExtra+".root -m "+str(signalMass) + " -n "+signal+signalExtra+" fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root &> "+name+"_"+str(signalMass)+"_run"+runs+version+".txt")
      system_call("combine -M HybridNew --frequentist --grid grid_mX"+str(signalMass)+"_"+signal+signalExtra+".root -m "+str(signalMass) + " -n "+signal+signalExtra+" --expectedFromGrid 0.5 fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root &> "+name+"_exp_"+str(signalMass)+"_run"+runs+version+".txt")
      system_call("combine -M HybridNew --frequentist --grid grid_mX"+str(signalMass)+"_"+signal+signalExtra+".root -m "+str(signalMass) + " -n "+signal+signalExtra+" --expectedFromGrid 0.84 fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root &>> "+name+"_exp_"+str(signalMass)+"_run"+runs+version+".txt")
      system_call("combine -M HybridNew --frequentist --grid grid_mX"+str(signalMass)+"_"+signal+signalExtra+".root -m "+str(signalMass) + " -n "+signal+signalExtra+" --expectedFromGrid 0.16 fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root &>> "+name+"_exp_"+str(signalMass)+"_run"+runs+version+".txt")
      system_call("combine -M HybridNew --frequentist --grid grid_mX"+str(signalMass)+"_"+signal+signalExtra+".root -m "+str(signalMass) + " -n "+signal+signalExtra+" --expectedFromGrid 0.975 fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root &>> "+name+"_exp_"+str(signalMass)+"_run"+runs+version+".txt")
      system_call("combine -M HybridNew --frequentist --grid grid_mX"+str(signalMass)+"_"+signal+signalExtra+".root -m "+str(signalMass) + " -n "+signal+signalExtra+" --expectedFromGrid 0.025 fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root &>> "+name+"_exp_"+str(signalMass)+"_run"+runs+version+".txt")
     else:
      for point in [0]:
       out=system_call("combine -m "+str(signalMass)+" -M HybridNew --LHCmode LHC-significance --fullBToys --saveHybridResult -s 10000"+str(int(point*100))+" --saveToys "+poi+" --fork 4 -T "+str(ntoys)+" -i 2 -n "+signal+signalExtra+" fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root")
       f = open(name+"_"+str(signalMass)+str(point)+"_run"+runs+version+".txt","w");f.write(out);f.close()
      system_call("hadd -f toys_mX"+str(signalMass)+"_"+signal+signalExtra+".root higgsCombine"+signal+signalExtra+".HybridNew.mH"+str(signalMass)+".10000*.root")

      system_call("combine -M HybridNew --LHCmode LHC-significance --readHybridResult --toysFile=toys_mX"+str(signalMass)+"_"+signal+signalExtra+".root "+poi+" -m "+str(signalMass) + " -n "+signal+signalExtra+" fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root &> "+name+"_"+str(signalMass)+"_run"+runs+version+".txt")

    else:
    
     out=system_call("combine -m "+str(signalMass)+" -M HybridNew --rule CLs --saveHybridResult --testStat "+method+poi+" --fork 4 -T "+str(ntoys)+" --clsAcc 0.1 -n "+signal+signalExtra+" fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root")
     f = open(name+"_"+str(signalMass)+"_run"+runs+version+".txt","w");f.write(out);f.close()
    
     out=system_call('root -q -b higgsCombine'+signal+signalExtra+'.HybridNew.mH'+str(signalMass)+'.root "${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx(\\"qmu_'+signal+str(signalMass)+signalExtra+'_'+testStat+version+'.root\\",'+str(signalMass)+',1,\\"x\\")"')
    
     out=system_call('root -q -b '+dire+'"extractSignificanceStats'+add+'.C(\\"'+signal+str(signalMass)+signalExtra+'_'+testStat+version+'\\")"')
     f = open(name+'_exp_'+str(signalMass)+'_run'+runs+version+'.txt',"w");f.write(out);f.close()
    
    # diagnostics
    diagnostic=True
    if diagnostic:
      def remove_prefix(text, prefix):
         if text.startswith(prefix):
             return text[len(prefix):]
         return text  # or whatever

      out=system_call("mkdir "+name+version)
      #out=system_call("combine -m "+str(signalMass)+" -M MaxLikelihoodFit "+poi+" --plots --out "+name+version+" -n "+remove_prefix(signalWithMass,"QCD")+" fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root")
      if withDiffNuisance:
        out=system_call("combine -m "+str(signalMass)+" -M FitDiagnostics "+poi+" --plots --out "+name+version+" -n "+remove_prefix(signalWithMass,"QCD")+" fixedMu_"+str(model)+remove_prefix(signalWithMass,"QCD")+"_run"+runs+".root")
        out=system_call("python3 diffNuisances.py -p x -a "+name+version+"/fitDiagnostics"+remove_prefix(signalWithMass,"QCD")+".root -A")
        print(out)

 for signalMass in signalMasses:
    limits[signalMass]=[]
    if testStat!="LEP" and (testStat!="LHC" or asym!=""):
     tname=name+"_exp_"+str(signalMass)+"_run"+runs+version+".txt"
    else:
     tname=name+"_"+str(signalMass)+"_run"+runs+version+".txt"
    try:
      print("open",tname)
      f=open(tname)
    except:
      print("file not found", tname)
      continue
    for lineb in f.readlines():
      for line in lineb.split("\\n"):
        if "Observed Limit" in line and asym:
           limits[signalMass]=[signalMass,float(line.strip().split(" ")[-1]),0]
        if "CLs = " in line and testStat=="LEP":
           limits[signalMass]=[signalMass,float(line.strip().split(" ")[-3]),float(line.strip().split(" ")[-1])]
        if "Observed CLs = " in line and testStat!="LEP":
           limits[signalMass]=[signalMass,float(line.strip().split(" ")[-1]),0]
        if "Limit:" in line and "95% CL" in line and testStat=="LHC" and asym=="":
           limits[signalMass]=[signalMass,float(line.strip().split(" ")[3]),0]
        if "Significance:" in line and asym:
           print("observed signficance (p-value): ",line.strip().split(" ")[-1].strip(")"))
        if "CLb = " in line and testStat=="LEP":
           print("observed signficance (p-value): ",ROOT.Math.normal_quantile_c((1.-float(line.strip().split(" ")[-3]))/2.,1),"(",(1.-float(line.strip().split(" ")[-3])),")")
        if "Observed CLb = " in line and testStat!="LEP":
           print("observed signficance (p-value): ",ROOT.Math.normal_quantile_c((1.-float(line.strip().split(" ")[-1]))/2.,1),"(",(1.-float(line.strip().split(" ")[-1])),")")
    if len(limits[signalMass])==0:
         limits[signalMass]+=[signalMass,0,0]
    try:
      tname=name+"_exp_"+str(signalMass)+"_run"+runs+version+".txt"
      f=open(tname)
    except:
      print("file not found", tname)
      continue
    for lineb in f.readlines():
      for line in lineb.split("\\n"):
        #if "Observed Limit" in line and asym:
        if "Expected" in line and asym:
           limits[signalMass]+=[float(line.strip().split(" ")[-1])]
        if "Expected CLs" in line:
          try:
           limits[signalMass]+=[float(line.strip().split(" ")[-1])]
          except:
           print("didn't find one point")
        if "Limit:" in line and "95% CL" in line and testStat=="LHC" and asym=="":
           limits[signalMass]+=[float(line.strip().split(" ")[3])]
    for i in range(len(limits[signalMass]),8):
         limits[signalMass]+=[0]

 if asym: # reorder expected limits to exp,-1,+1,-2,+2
   for signalMass in signalMasses:
     limits[signalMass]=[limits[signalMass][0],limits[signalMass][1],limits[signalMass][2],limits[signalMass][5],limits[signalMass][6],limits[signalMass][4],limits[signalMass][7],limits[signalMass][3]]

 print(limits)
 name=name+"_run"+runs+version+".txt"
 f=open(name,"w")
 f.write(str([limits[signalMass] for signalMass in signalMasses]))
 f.close()
