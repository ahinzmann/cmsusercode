import os, sys
from ROOT import * 
from DataFormats.FWLite import Events,Handle
import array
from math import *

print "root"

gSystem.Load("libFWCoreFWLite.so");
FWLiteEnabler.enable();

gROOT.SetStyle("Plain")
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetTitleOffset(1.2,"Y")
gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadBottomMargin(0.15)
gStyle.SetPadTopMargin(0.03)
gStyle.SetPadRightMargin(0.05)
gStyle.SetMarkerSize(1.5)
gStyle.SetHistLineWidth(1)
gStyle.SetStatFontSize(0.020)
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetNdivisions(510, "XYZ")
gStyle.SetLegendBorderSize(0)

colors=[1,2,3,4,6,7,8,9,10,11,12,13,1,2,3,4,6,7,8,9,10,11,12,13,1,2,3,4,6,7,8,9,10,11,12,13,1,2,3,4,6,7,8,9,10,11,12,13,]

doJES=True
use_UL=True

if doJES:
  JECsources = ["AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation",
"SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta",
"RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeFSR", "RelativeSample",
"RelativeStatFSR", "RelativeStatEC", "RelativeStatHF", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
"PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF",
"Total",
  ]
  
  JESuncertainties={}
  JESparameters={}
  print "load 2016 uncertainties"
  for source in JECsources:
    JESparameters[source+"2016"] = JetCorrectorParameters("/afs/desy.de/user/h/hinzmann/uhh102/CMSSW_10_2_10/src/UHH2/JECDatabase/textFiles/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt", source)
    JESuncertainties[source+"2016"] = JetCorrectionUncertainty(JESparameters[source+"2016"])
  print "load 2017 uncertainties"
  for source in JECsources:
    JESparameters[source+"2017"] = JetCorrectorParameters("/afs/desy.de/user/h/hinzmann/uhh102/CMSSW_10_2_10/src/UHH2/JECDatabase/textFiles/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt", source)
    JESuncertainties[source+"2017"] = JetCorrectionUncertainty(JESparameters[source+"2017"])
  print "load 2018 uncertainties"
  for source in JECsources:
    JESparameters[source+"2018"] = JetCorrectorParameters("/afs/desy.de/user/h/hinzmann/uhh102/CMSSW_10_2_10/src/UHH2/JECDatabase/textFiles/Autumn18_V19_MC/Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt", source)
    JESuncertainties[source+"2018"] = JetCorrectionUncertainty(JESparameters[source+"2018"])
  print "load UL16preVFP uncertainties"
  for source in JECsources:
    JESparameters[source+"UL16preVFP"] = JetCorrectorParameters("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JECDatabase/textFiles/Summer19UL16APV_V7_MC/Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.txt", source)
    JESuncertainties[source+"UL16preVFP"] = JetCorrectionUncertainty(JESparameters[source+"UL16preVFP"])
  print "load UL16postVFP uncertainties"
  for source in JECsources:
    JESparameters[source+"UL16postVFP"] = JetCorrectorParameters("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JECDatabase/textFiles/Summer19UL16_V7_MC/Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.txt", source)
    JESuncertainties[source+"UL16postVFP"] = JetCorrectionUncertainty(JESparameters[source+"UL16postVFP"])
  print "load UL17 uncertainties"
  for source in JECsources:
    JESparameters[source+"UL17"] = JetCorrectorParameters("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JECDatabase/textFiles/Summer19UL17_V5_MC/Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.txt", source)
    JESuncertainties[source+"UL17"] = JetCorrectionUncertainty(JESparameters[source+"UL17"])
  print "load UL18 uncertainties"
  for source in JECsources:
    JESparameters[source+"UL18"] = JetCorrectorParameters("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JECDatabase/textFiles/Summer19UL18_V5_MC/Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt", source)
    JESuncertainties[source+"UL18"] = JetCorrectionUncertainty(JESparameters[source+"UL18"])

def createPlots(sample,prefix,xsec,massbins,year):
    files=[]
    if sample.endswith(".txt"):
        filelist=open(sample)
	for line in filelist.readlines():
	    if ".root" in line:
	        files+=[line.strip()]
    elif "HT" in sample:
        print "/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcd"+year+"feb2023/"
        folders=os.listdir("/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcd"+year+"feb2023/")
	for folder in folders:
	  if sample in folder:
            files+=["dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcd"+year+"feb2023/"+folder]
    elif "UL" in sample:
        print "/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/data"+year+"-12Feb2024/"
        folders=os.listdir("/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/data"+year+"-12Feb2024/")
	for folder in folders:
            files+=["dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/data"+year+"-12Feb2024/"+folder]
    else:
        folders=os.listdir("/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dijet_angular/")
	for folder in folders:
	  if sample in folder:
            files+=["/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dijet_angular/"+folder+"/GEN.root"]
	    #break

    print files
    
    prunedgenjets_handle=Handle("std::vector<reco::GenJet>")
    prunedgenjets_label="ak4GenJets"

    plots=[]
    for massbin in massbins:
      plots += [TH1F(prefix+'#chi'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';#chi;N',15,1,16)]
    plots += [TH1F(prefix+'mass',';dijet mass;N',260,0,13000)]
    if doJES:
      for source in JECsources:
       for massbin in massbins:
        plots += [TH1F(prefix+'#chi'+str(massbin).strip("()").replace(',',"_").replace(' ',"")+"_"+source+"Up",';#chi;N',15,1,16)]
       plots += [TH1F(prefix+'mass'+"_"+source+"Up",';dijet mass;N',260,0,13000)]
       for massbin in massbins:
        plots += [TH1F(prefix+'#chi'+str(massbin).strip("()").replace(',',"_").replace(' ',"")+"_"+source+"Down",';#chi;N',15,1,16)]
       plots += [TH1F(prefix+'mass'+"_"+source+"Down",';dijet mass;N',260,0,13000)]
    
    for plot in plots:
        plot.Sumw2()
	print plot.GetName()

    event_count=0
    tree='Events'
    if "HT" in sample: tree="tree"
    events=TChain(tree)
    for f in files[:]:
      events.Add(f)
    
    nevents=events.GetEntries()
    print sample,nevents,xsec

    scales=[1]
    if doJES:
      for s in JECsources:
       scales+=[s+"_Up",s]

    for event in events:
       #if event_count>10000000:break
       #if event_count>1000000:break
       #if event_count>100000:break
      
       event_count+=1
       if event_count%10000==1: print "event",event_count

       jet1=TLorentzVector()
       jet2=TLorentzVector()
       if "NANOAOD" in f:
         if event.jetAK4_pt2_nom<30: continue
         jet1.SetPtEtaPhiM(event.jetAK4_pt1_nom,event.jetAK4_eta1,event.jetAK4_phi1,event.jetAK4_mass1*event.jetAK4_pt1_nom/event.jetAK4_pt1)
         jet2.SetPtEtaPhiM(event.jetAK4_pt2_nom,event.jetAK4_eta2,event.jetAK4_phi2,event.jetAK4_mass2*event.jetAK4_pt2_nom/event.jetAK4_pt2)
       elif "HT" in sample:
         if event.jetAK4_pt2<30: continue
         jet1.SetPtEtaPhiM(event.jetAK4_pt1,event.jetAK4_eta1,event.jetAK4_phi1,event.jetAK4_mass1)
         jet2.SetPtEtaPhiM(event.jetAK4_pt2,event.jetAK4_eta2,event.jetAK4_phi2,event.jetAK4_mass2)
       else:
         if len(event.recoGenJets_ak4GenJets__GEN.product())<2: continue
         jet1.SetPtEtaPhiM(event.recoGenJets_ak4GenJets__GEN.product()[0].pt(),event.recoGenJets_ak4GenJets__GEN.product()[0].eta(),event.recoGenJets_ak4GenJets__GEN.product()[0].phi(),event.recoGenJets_ak4GenJets__GEN.product()[0].mass())
         jet2.SetPtEtaPhiM(event.recoGenJets_ak4GenJets__GEN.product()[1].pt(),event.recoGenJets_ak4GenJets__GEN.product()[1].eta(),event.recoGenJets_ak4GenJets__GEN.product()[1].phi(),event.recoGenJets_ak4GenJets__GEN.product()[1].mass())
       if jet1.Pt()<100 or jet2.Pt()<100 or abs(jet1.Eta())>3 or abs(jet2.Eta())>3: continue
       
       irec=0
       mjj=(jet1+jet2).M()
       chi=exp(abs(jet1.Rapidity()-jet2.Rapidity()))
       yboost=abs(jet1.Rapidity()+jet2.Rapidity())/2.
       if "NANOAOD" in f:
         if not event.jetAK4_TightID1 or not event.jetAK4_TightID2: continue
       for scale in scales:
         try:
	  if mjj>1000 and chi<16. and yboost<1.11 and scale!=1:
           jes=JESuncertainties[scale.replace("_Up","")+year]
           jes.setJetPt(jet1.Pt())
           jes.setJetEta(jet1.Eta())
	   factor1=1.
	   if scale[-1]=="p":
	     factor1+=jes.getUncertainty(True)
	   else:
             factor1-=jes.getUncertainty(False)
	   jes.setJetPt(jet2.Pt())
           jes.setJetEta(jet2.Eta())
	   factor2=1.
	   if scale[-1]=="p":
             factor2+=jes.getUncertainty(True)
	   else:
             factor2-=jes.getUncertainty(False)
           mjj=(jet1*factor1+jet2*factor2).M()
	   #if mjj>4800 and mjj<5400: print scale, factor1, factor2
         except:
	   print "WARNING: JEC extraction failed", jet1.Pt(), jet1.Eta(), jet2.Pt(), jet2.Eta()
	 for massbin in massbins:
           if yboost<1.11 and mjj>=massbin[0] and mjj<massbin[1]:
             plots[irec].Fill(chi)
           irec+=1
         plots[irec].Fill(mjj)
         irec+=1

    #jets="recoGenJets_ak4GenJets__GEN.obj"
    #yboost='abs('+jets+'[0].y()+'+jets+'[1].y())/2.'
    #chi='exp(abs('+jets+'[0].y()-'+jets+'[1].y()))'
    #mass='sqrt(pow('+jets+'[0].energy()+'+jets+'[1].energy(),2)-pow('+jets+'[0].px()+'+jets+'[1].px(),2)-pow('+jets+'[0].py()+'+jets+'[1].py(),2)-pow('+jets+'[0].pz()+'+jets+'[1].pz(),2))'
    #for massbin in massbins:
    #  events.Project(prefix+'#chi'+str(massbin).strip("()").replace(',',"_").replace(' ',""),chi,'('+yboost+'<1.11)*('+mass+'>='+str(massbin[0])+')*('+mass+'<'+str(massbin[1])+')')
    #  #events.Project(prefix+'y_{boost}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),yboost,'('+chi+'<16)*('+mass+'>='+str(massbin[0])+')*('+mass+'<='+str(massbin[1])+')')
    if "HT" in sample:
      event_count/=nevents
    if not "UL" in sample:
      for plot in plots:
        if event_count>0:
          plot.Scale(xsec/event_count)
    return plots

if __name__ == '__main__':
    sets=[]
    i=0
    if use_UL:
      years=["UL16preVFP","UL16postVFP","UL17","UL18"]
    else:
      years=[2016,2017,2018]
    for name in ["QCD","QCDCIplusLL10000","QCDmadgraph"]:
      for year in years:
        for bin in [1,2,3,4,5,6,7]:
	   print i,name,bin,year
	   sets+=[(name,bin,year)]
	   i+=1
    for year in years:
      print "data",1,year
      sets+=[("data",1,year)]
    print "sets",len(sets)
    wait=False
    #name="QCD"
    #name="QCDmadgraph"
    #name="QCDCIplusLL10000"
    name="data"
    bin=1 #1-6
    year="UL16preVFP"
    if len(sys.argv)>1:
       name,bin,year = sets[int(sys.argv[1])]
    #if year!=2018: print "SKIPPING 2016 and 2017 !!!!!"; STOP
    prefix="data/datacard_shapelimit13TeV_"+name+"_JESfeb2023_"+str(year)+"_"+str(bin)
 
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
              (1,2,3,4,5,6,7,8,9,10,12,14,16),
              (1,2,3,4,5,6,7,8,9,10,12,14,16),
              (1,2,3,4,5,6,7,8,9,10,12,14,16),
              (1,2,3,4,5,6,7,8,9,10,12,14,16),
              (1,2,3,4,5,6,7,8,9,10,12,14,16),
              (1,2,3,4,5,6,7,8,9,10,12,14,16),
              (1,2,3,4,5,6,7,8,9,10,12,14,16),
              (1,2,3,4,5,6,7,8,9,10,12,14,16),
              (1,2,3,4,5,6,7,8,9,10,12,14,16),
              (1,2,3,4,5,6,7,8,9,10,12,14,16),
              (1,2,3,4,5,6,7,8,9,10,12,14,16),
              ]
    massbins=[(1200,1500),
              (1500,1900),
              (1900,2400),
              (2400,3000),
              (3000,3600),
              (3600,4200),
              (4200,4800),
              (4800,5400),
              (5400,6000),
	      (6000,13000),
              (6000,6600),
	      (6600,13000),
	      (6000,7000),
              (7000,13000),
              ]

    if bin==1 and name=="QCD": 
      samples=[("QCD",[("pythia8_ci_m4300_13000_50000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
	    ]
    if bin==2 and name=="QCD":
      samples=[("QCD",[("pythia8_ci_m3800_4300_50000_1_0_0_13TeV_Nov14",5.867e-09),
		       ]),
	    ]
    if bin==3 and name=="QCD":
      samples=[("QCD",[("pythia8_ci_m3300_3800_50000_1_0_0_13TeV_Nov14",1.863e-08),
		       ]),
	    ]
    if bin==4 and name=="QCD":
      samples=[("QCD",[("pythia8_ci_m2800_3300_50000_1_0_0_13TeV_Nov14",6.446e-08),
		       ]),
	    ]
    if bin==5 and name=="QCD":
      samples=[("QCD",[("pythia8_ci_m2400_2800_50000_1_0_0_13TeV_Nov14",1.649e-07),
		       ]),
	    ]
    if bin==6 and name=="QCD":
      samples=[("QCD",[("pythia8_ci_m1900_2400_50000_1_0_0_13TeV_Nov14",8.836e-07),
		       ]),
	    ]
    if bin==7 and name=="QCD": 
      samples=[("QCD",[("pythia8_ci_m1500_1900_50000_1_0_0_13TeV_Nov14",3.307e-06),
		       ]),
	    ]

    if bin==1 and name=="QCDCIplusLL10000": 
      samples=[("QCDCIplusLL10000",[("pythia8_ci_m4300_13000_10000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
	    ]
    if bin==2 and name=="QCDCIplusLL10000":
      samples=[("QCDCIplusLL10000",[("pythia8_ci_m3800_4300_10000_1_0_0_13TeV_Nov14",5.867e-09),
		       ]),
	    ]
    if bin==3 and name=="QCDCIplusLL10000":
      samples=[("QCDCIplusLL10000",[("pythia8_ci_m3300_3800_10000_1_0_0_13TeV_Nov14",1.863e-08),
		       ]),
	    ]
    if bin==4 and name=="QCDCIplusLL10000":
      samples=[("QCDCIplusLL10000",[("pythia8_ci_m2800_3300_10000_1_0_0_13TeV_Nov14",6.446e-08),
		       ]),
	    ]
    if bin==5 and name=="QCDCIplusLL10000":
      samples=[("QCDCIplusLL10000",[("pythia8_ci_m2400_2800_10000_1_0_0_13TeV_Nov14",1.649e-07),
		       ]),
	    ]
    if bin==6 and name=="QCDCIplusLL10000":
      samples=[("QCDCIplusLL10000",[("pythia8_ci_m1900_2400_10000_1_0_0_13TeV_Nov14",8.836e-07),
		       ]),
	    ]
    if bin==7 and name=="QCDCIplusLL10000": 
      samples=[("QCDCIplusLL10000",[("pythia8_ci_m1500_1900_10000_1_0_0_13TeV_Nov14",3.307e-06),
		       ]),
	    ]

    if str(year)=="2016":
      if bin==1 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT200to300_RunII_2016v3",1712000./56709875),
  	  		 ]),
  	      ]
      if bin==2 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT300to500_RunII_2016v3",347700./53096517),
  	  		 ]),
  	      ]
      if bin==3 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT500to700_RunII_2016v3",32100./52906552),
  	  		 ]),
  	      ]
      if bin==4 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT700to1000_RunII_2016v3",6831./36741540),
  	  		 ]),
  	      ]
      if bin==5 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT1000to1500_RunII_2016v3",1207./15210939),
  	  		 ]),
  	      ]
      if bin==6 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT1500to2000_RunII_2016v3",119.9/11839357),
  	  		 ]),
  	      ]
      if bin==7 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT2000toInf_RunII_2016v3",25.24/5947849),
  	  		 ]),
  	      ]
    if str(year)=="2017":
      if bin==1 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT200to300_RunII_94X_v2",1545000./58990434),
  	  		 ]),
  	      ]
      if bin==2 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT300to500_RunII_94X_v2",323300./58748739),
  	  		 ]),
  	      ]
      if bin==3 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT500to700_RunII_94X_v2",30000./54366431),
  	  		 ]),
  	      ]
      if bin==4 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT700to1000_RunII_94X_v2",6324./46924322),
  	  		 ]),
  	      ]
      if bin==5 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT1000to1500_RunII_94X_v2",1090./16495598),
  	  		 ]),
  	      ]
      if bin==6 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT1500to2000_RunII_94X_v2",101./11196479),
  	  		 ]),
  	      ]
      if bin==7 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT2000toInf_RunII_94X_v2",20.43/5362513),
  	  		 ]),
  	      ]
    if str(year)=="2018":
      if bin==1 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT200to300_RunII_102X_v1",1461000./54289442),
  	  		 ]),
  	      ]
      if bin==2 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT300to500_RunII_102X_v1",311900./54512704),
  	  		 ]),
  	      ]
      if bin==3 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT500to700_RunII_102X_v1",29070./53919811),
  	  		 ]),
  	      ]
      if bin==4 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT700to1000_RunII_102X_v1",5962./48158738),
  	  		 ]),
  	      ]
      if bin==5 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT1000to1500_RunII_102X_v1",1005./14945819),
  	  		 ]),
  	      ]
      if bin==6 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT1500to2000_RunII_102X_v1",101.8/10707847),
  	  		 ]),
  	      ]
      if bin==7 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("dijetChiQCD_HT2000toInf_RunII_102X_v1",20.54/5329144),
  	  		 ]),
  	      ]

    if str(year)=="UL16preVFP":
      if bin==1 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT200to300_RunII_106X",1710000./17969592),
  	  		 ]),
  	      ]
      if bin==2 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT300to500_RunII_106X",347500./13586390),
  	  		 ]),
  	      ]
      if bin==3 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT500to700_RunII_106X",30363.051/55497082),
  	  		 ]),
  	      ]
      if bin==4 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT700to1000_RunII_106X",6428.869/15242034),
  	  		 ]),
  	      ]
      if bin==5 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1000to1500_RunII_106X",1122.659/13559959),
  	  		 ]),
  	      ]
      if bin==6 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1500to2000_RunII_106X",108.163/9661950),
  	  		 ]),
  	      ]
      if bin==7 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT2000toInf_RunII_106X",22.008/4827641),
  	  		 ]),
  	      ]
    if str(year)=="UL16postVFP":
      if bin==1 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT200to300_RunII_106X",1710000./42723038),
  	  		 ]),
  	      ]
      if bin==2 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT300to500_RunII_106X",347500./45502889),
  	  		 ]),
  	      ]
      if bin==3 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT500to700_RunII_106X",30363.051/15066884),
  	  		 ]),
  	      ]
      if bin==4 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT700to1000_RunII_106X",6428.869/13714842),
  	  		 ]),
  	      ]
      if bin==5 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1000to1500_RunII_106X",1122.659/12416669),
  	  		 ]),
  	      ]
      if bin==6 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1500to2000_RunII_106X",108.163/9244228),
  	  		 ]),
  	      ]
      if bin==7 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT2000toInf_RunII_106X",22.008/4843949),
  	  		 ]),
  	      ]
    if str(year)=="UL17":
      if bin==1 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT200to300_RunII_106X",1710000./42316128),
  	  		 ]),
  	      ]
      if bin==2 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT300to500_RunII_106X",347500./42914024),
  	  		 ]),
  	      ]
      if bin==3 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT500to700_RunII_106X",30363.051/35745565),
  	  		 ]),
  	      ]
      if bin==4 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT700to1000_RunII_106X",6428.869/33646855),
  	  		 ]),
  	      ]
      if bin==5 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1000to1500_RunII_106X",1122.659/10136610),
  	  		 ]),
  	      ]
      if bin==6 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1500to2000_RunII_106X",108.163/7528926),
  	  		 ]),
  	      ]
      if bin==7 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT2000toInf_RunII_106X",22.008/4089387),
  	  		 ]),
  	      ]
    if str(year)=="UL18":
      if bin==1 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT200to300_RunII_106X",1710000./56298746),
  	  		 ]),
  	      ]
      if bin==2 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT300to500_RunII_106X",347500./60991701),
  	  		 ]),
  	      ]
      if bin==3 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT500to700_RunII_106X",30363.051/48640047),
  	  		 ]),
  	      ]
      if bin==4 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT700to1000_RunII_106X",6428.869/47925782),
  	  		 ]),
  	      ]
      if bin==5 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1000to1500_RunII_106X",1122.659/14244456),
  	  		 ]),
  	      ]
      if bin==6 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1500to2000_RunII_106X",108.163/10751607),
  	  		 ]),
  	      ]
      if bin==7 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT2000toInf_RunII_106X",22.008/5278880),
  	  		 ]),
  	      ]

    if str(year)=="UL16preVFP" and name=="data": 
    	samples=[("data",[("UL16preVFP",1),
  	  		 ]),
  	      ]
    if str(year)=="UL16postVFP" and name=="data": 
    	samples=[("data",[("UL16postVFP",1),
  	  		 ]),
  	      ]
    if str(year)=="UL17" and name=="data": 
    	samples=[("data",[("UL17",1),
  	  		 ]),
  	      ]
    if str(year)=="UL18" and name=="data": 
    	samples=[("data",[("UL18",1),
  	  		 ]),
  	      ]

    #xsecs=eval(open("xsecs_13TeV.txt").readline())
    #print xsecs

    chi_binnings=[]
    for mass_bin in chi_bins:
        chi_binnings+=[array.array('d')]
        for chi_bin in mass_bin:
            chi_binnings[-1].append(chi_bin)
  
    print prefix, samples

    plots=[]
    for name,files in samples:
      plots+=[[]]
      i=0
      for filename,xsec in files:
        i+=1
        ps=createPlots(filename,name,xsec,massbins,str(year))
        if i==1:
          plots[-1]+=ps
	else:
	  for j in range(len(plots[-1])):
            plots[-1][j].Add(ps[j])

    out=TFile(prefix + '_chi.root','RECREATE')
    for j in range(len(massbins)+1):
      for i in range(len(samples)):
       if j==len(massbins): # mass histogram
         plots[i][j].Write()
       else: 
        #if plots[i][j].Integral()>0:
        #  plots[i][j].Scale(expectedevents[j]/plots[i][j].Integral())
        plots[i][j]=plots[i][j].Rebin(len(chi_binnings[j])-1,plots[i][j].GetName()+"_rebin1",chi_binnings[j])
	if samples[i][0]=="QCD":
	   # data
	   plots[i][j].Write(plots[i][j].GetName().replace("QCD","data_obs"))
	   # ALT
	   clone=plots[i][j].Clone(plots[i][j].GetName().replace("QCD",samples[-1][0]+"_ALT"))
	   clone.Write()
	   # QCD
           plots[i][j].Scale(1e-10)
           plots[i][j].Write()
	   # QCD backup
	   clonebackup=plots[i][j].Clone(plots[i][j].GetName()+"_backup")
	   clonebackup.Write()
	else:
	   # signal
	   clone=plots[i][j]
	   clone.Write()
	   # signal backup
	   clonebackup=plots[i][j].Clone(plots[i][j].GetName()+"_backup")
	   clonebackup.Write()
    if doJES:
     for j in range((len(massbins)+1),(len(massbins)+1)+2*(len(massbins)+1)*len(JECsources)):
      for i in range(len(samples)):
       if j%(len(massbins)+1)==len(massbins): # mass histogram
         plots[i][j].Write()
       else: 
        plots[i][j]=plots[i][j].Rebin(len(chi_binnings[j%len(massbins)])-1,plots[i][j].GetName()+"_rebin1",chi_binnings[j%len(massbins)])
	if samples[i][0]=="QCD":
	   # QCD
	   clone=plots[i][j]
	   clone.Write()
	   # QCD backup
	   clonebackup=plots[i][j].Clone(plots[i][j].GetName()+"_backup")
	   clonebackup.Write()
	else:
	   # signal
	   clone=plots[i][j]
	   clone.Write()
	   # signal backup
	   clonebackup=plots[i][j].Clone(plots[i][j].GetName()+"_backup")
	   clonebackup.Write()

    for j in range((len(massbins)+1)+doJES*2*(len(massbins)+1)*len(JECsources)):
      for i in range(len(samples)):
	if plots[i][j].Integral()>0:
          plots[i][j].Scale(1./plots[i][j].Integral())
	for b in range(plots[i][j].GetXaxis().GetNbins()):
          plots[i][j].SetBinContent(b+1,plots[i][j].GetBinContent(b+1)/plots[i][j].GetBinWidth(b+1))
          plots[i][j].SetBinError(b+1,plots[i][j].GetBinError(b+1)/plots[i][j].GetBinWidth(b+1))

    canvas = TCanvas("","",0,0,800,600)
    canvas.Divide(4,3)

    legends=[]
    ratios=[]
    for j in range(len(massbins)):
      canvas.cd(j+1)
      plots[0][j].Draw("he")
      ratio=plots[0][j].Clone("ratio"+str(j))
      ratios+=[ratio]
      ratio.Divide(plots[0][j],plots[0][j])
      ratio.GetYaxis().SetRangeUser(0.95,1.1)
      for b in range(ratio.GetNbinsX()):
	  ratio.SetBinError(b+1,0)
      ratio.Draw("he")
      print "number of events passed:",plots[0][j].GetEntries()
      legend1=TLegend(0.6,0.6,0.9,0.9,(str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV"))
      legends+=[legend1]
      legend1.AddEntry(plots[0][j],samples[0][0],"l")
      for i in range(1,len(samples)):
        plots[i][j].SetLineColor(i+2)
        plots[i][j].Draw("hesame")
        legend1.AddEntry(plots[i][j],samples[i][0],"l")
      if doJES:
       for k in range(len(JECsources)):
        plots[0][j+(len(massbins)+1)+2*(len(massbins)+1)*k].SetLineColor(colors[k+1])
        plots[0][j+(len(massbins)+1)+2*(len(massbins)+1)*k].SetLineStyle(2)
	plots[0][j+(len(massbins)+1)+2*(len(massbins)+1)*k].Divide(plots[0][j+(len(massbins)+1)+2*(len(massbins)+1)*k],plots[0][j])
        for b in range(plots[0][j+(len(massbins)+1)+2*(len(massbins)+1)*k].GetNbinsX()):
	    plots[0][j+(len(massbins)+1)+2*(len(massbins)+1)*k].SetBinError(b+1,0)
        plots[0][j+(len(massbins)+1)+2*(len(massbins)+1)*k].Draw("hesame")
        legend1.AddEntry(plots[0][j+(len(massbins)+1)+2*(len(massbins)+1)*k],plots[0][j+(len(massbins)+1)+2*(len(massbins)+1)*k].GetName().split("_")[-2],"l")
        plots[0][j+2*(len(massbins)+1)+2*(len(massbins)+1)*k].SetLineColor(colors[k+1])
        plots[0][j+2*(len(massbins)+1)+2*(len(massbins)+1)*k].SetLineStyle(3)
	plots[0][j+2*(len(massbins)+1)+2*(len(massbins)+1)*k].Divide(plots[0][j+2*(len(massbins)+1)+2*(len(massbins)+1)*k],plots[0][j])
        for b in range(plots[0][j+2*(len(massbins)+1)+2*(len(massbins)+1)*k].GetNbinsX()):
	    plots[0][j+2*(len(massbins)+1)+2*(len(massbins)+1)*k].SetBinError(b+1,0)
        plots[0][j+2*(len(massbins)+1)+2*(len(massbins)+1)*k].Draw("hesame")
        legend1.AddEntry(plots[0][j+2*(len(massbins)+1)+2*(len(massbins)+1)*k],plots[0][j+2*(len(massbins)+1)+2*(len(massbins)+1)*k].GetName().split("_")[-2],"l")
      legend1.SetTextSize(0.04)
      legend1.SetFillStyle(0)
      legend1.Draw("same")

    canvas.SaveAs(prefix + '_chi.pdf')
    #canvas.SaveAs(prefix + '_chi.eps')
    #if wait:
    #    os.system("ghostview "+prefix + '_chi.eps')

