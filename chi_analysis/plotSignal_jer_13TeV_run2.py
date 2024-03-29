import os, sys
from ROOT import * 
from DataFormats.FWLite import Events,Handle
import array
from math import *
import numpy

print "root"

gSystem.Load("libFWCoreFWLite.so");
FWLiteEnabler.enable();

gROOT.SetStyle("Plain")
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

doJER=True
use_UL=True

if doJER:
  JERuncertainties=["JER","JER1","JER2"]
  JER={}
  JERSF={}
  print "load 2016 JER"
  JER["2016"] = JME.JetResolutionObject("Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt")
  JERSF["2016"] = JME.JetResolutionObject("Summer16_25nsV1_MC_SF_AK4PFchs.txt")
  print "load 2017 JER"
  JER["2017"] = JME.JetResolutionObject("Fall17_V3_MC_PtResolution_AK4PFchs.txt")
  JERSF["2017"] = JME.JetResolutionObject("Fall17_V3_MC_SF_AK4PFchs.txt")
  print "load 2018 JER"
  JER["2018"] = JME.JetResolutionObject("Autumn18_V7_MC_PtResolution_AK4PFchs.txt")
  JERSF["2018"] = JME.JetResolutionObject("Autumn18_V7_MC_SF_AK4PFchs.txt")
  print "load UL16preVFP JER"
  JER["UL16preVFP"] = JME.JetResolutionObject("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JRDatabase/textFiles/Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.txt")
  JERSF["UL16preVFP"] = JME.JetResolutionObject("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JRDatabase/textFiles/Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt")
  print "load UL16postVFP JER"
  JER["UL16postVFP"] = JME.JetResolutionObject("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JRDatabase/textFiles/Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt")
  JERSF["UL16postVFP"] = JME.JetResolutionObject("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JRDatabase/textFiles/Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt")
  print "load UL17 JER"
  JER["UL17"] = JME.JetResolutionObject("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JRDatabase/textFiles/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.txt")
  JERSF["UL17"] = JME.JetResolutionObject("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JRDatabase/textFiles/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_SF_AK4PFchs.txt")
  print "load UL18 JER"
  JER["UL18"] = JME.JetResolutionObject("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt")
  JERSF["UL18"] = JME.JetResolutionObject("/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt")

def createPlots(sample,prefix,xsec,massbins,year):
    files=[]
    print sample,year
    if sample.endswith(".txt"):
        filelist=open(sample)
	for line in filelist.readlines():
	    if ".root" in line:
	        files+=[line.strip()]
    elif "HT" in sample:
        folders=os.listdir("/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcd"+year+"oct/")
	for folder in folders:
	  if sample in folder:
            files+=["dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcd"+year+"oct/"+folder]
    else:
        folders=os.listdir("/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dijet_angular/")
	for folder in folders:
	  if sample in folder:
            files+=["/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dijet_angular/"+folder+"/GEN.root"]
	    #break

    prunedgenjets_handle=Handle("std::vector<reco::GenJet>")
    prunedgenjets_label="ak4GenJets"

    plots=[]
    for massbin in massbins:
      plots += [TH1F(prefix+'#chi'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';#chi;N',15,1,16)]
    plots += [TH1F(prefix+'mass',';dijet mass;N',260,0,13000)]
    if doJER:
      for source in JERuncertainties:
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
    if doJER:
      for s in JERuncertainties:
       scales+=[s+"_Up",s]

    for event in events:
       #if event_count>10000000:break
       #if event_count>1000000:break
       #if event_count>100000:break
      
       event_count+=1
       if event_count%10000==1: print "event",event_count

       jet1=TLorentzVector()
       jet2=TLorentzVector()
       if "HT" in sample:
         if event.genJetAK4_pt2<30: continue
         jet1.SetPtEtaPhiM(event.genJetAK4_pt1,event.genJetAK4_eta1,event.genJetAK4_phi1,event.genJetAK4_mass1)
         jet2.SetPtEtaPhiM(event.genJetAK4_pt2,event.genJetAK4_eta2,event.genJetAK4_phi2,event.genJetAK4_mass2)
       else:
         if len(event.recoGenJets_ak4GenJets__GEN.product())<2: continue
         jet1.SetPtEtaPhiM(event.recoGenJets_ak4GenJets__GEN.product()[0].pt(),event.recoGenJets_ak4GenJets__GEN.product()[0].eta(),event.recoGenJets_ak4GenJets__GEN.product()[0].phi(),event.recoGenJets_ak4GenJets__GEN.product()[0].mass())
         jet2.SetPtEtaPhiM(event.recoGenJets_ak4GenJets__GEN.product()[1].pt(),event.recoGenJets_ak4GenJets__GEN.product()[1].eta(),event.recoGenJets_ak4GenJets__GEN.product()[1].phi(),event.recoGenJets_ak4GenJets__GEN.product()[1].mass())
       if jet1.Pt()<100 or jet2.Pt()<100 or abs(jet1.Eta())>3 or abs(jet2.Eta())>3: continue

       irec=0
       mjj=(jet1+jet2).M()
       chi=exp(abs(jet1.Rapidity()-jet2.Rapidity()))
       yboost=abs(jet1.Rapidity()+jet2.Rapidity())/2.
       for scale in scales:
         if True:
	 #try:
	  if mjj>1000 and chi<16. and yboost<1.11:
           pars=JME.JetParameters()
	   pars.setJetPt(jet1.Pt())
	   pars.setJetEta(jet1.Eta())
	   pars.setRho(1.) #rho
           record=JER[year].getRecord(pars)
	   factor1=JER[year].evaluateFormula(record,pars)
           record=JERSF[year].getRecord(pars)
	   if scale!=1 and (("JER1_" in scale and abs(jet1.Eta())<=1.93) or ("JER2_" in scale and abs(jet1.Eta())>1.93) or ("JER_" in scale)):
	     if scale[-1]=="p":
	       factor1*=record.getParametersValues()[2]
	     else:
               factor1*=record.getParametersValues()[1]
           else:
             factor1*=record.getParametersValues()[0]
	   factor1=numpy.random.normal(1,factor1)

           pars=JME.JetParameters()
	   pars.setJetPt(jet2.Pt())
	   pars.setJetEta(jet2.Eta())
	   pars.setRho(1.) #rho
           record=JER[year].getRecord(pars)
	   factor2=JER[year].evaluateFormula(record,pars)
           record=JERSF[year].getRecord(pars)
	   if scale!=1 and (("JER1_" in scale and abs(jet2.Eta())<=1.93) or ("JER2_" in scale and abs(jet2.Eta())>1.93) or ("JER_" in scale)):
	     if scale[-1]=="p":
	       factor2*=record.getParametersValues()[2]
	     else:
               factor2*=record.getParametersValues()[1]
           else:
             factor2*=record.getParametersValues()[0]
	   factor2=numpy.random.normal(1,factor2)

	   mjj=(jet1*factor1+jet2*factor2).M()
	   #if mjj>4800 and mjj<5400: print scale, factor1, factor2
         #except:
	 #  print "WARNING: JER extraction failed", jet1.Pt(), jet1.Eta(), jet2.Pt(), jet2.Eta()
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
    for name in ["QCD","QCDCIplusLL10000","QCDmadgraph"]: #"QCDherwig",
      for year in years:
        for bin in [1,2,3,4,5,6,7]:
	   print i,name,bin,year
	   sets+=[(name,bin,year)]
	   i+=1
    print "sets",len(sets)
    wait=False
    name="QCD"
    #name="QCDCIplusLL10000"
    bin=1 #1-6
    year="2016"
    if len(sys.argv)>1:
       name,bin,year = sets[int(sys.argv[1])]
    #if year!=2018: print "SKIPPING 2016 and 2017 !!!!!"; STOP
    prefix="datacard_shapelimit13TeV_"+name+"_JER_"+str(year)+"_"+str(bin)
 
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
              #(1,2,3,4,5,6,7,8,9,10,12,14,16),
              #(1,2,3,4,5,6,7,8,9,10,12,14,16),
              #(1,2,3,4,5,6,7,8,9,10,12,14,16),
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
	      #(6000,13000),
              #(6000,6600),
	      #(6600,13000),
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

    if bin==1 and name=="QCDherwig": 
      samples=[("QCDherwig",[("herwigpp_qcd_m4300_13000___Nov28",0.00253067e-6),
		       ]),
	    ]
    if bin==2 and name=="QCDherwig":
      samples=[("QCDherwig",[("herwigpp_qcd_m3800_4300___Nov28",0.00404967e-6),
		       ]),
	    ]
    if bin==3 and name=="QCDherwig":
      samples=[("QCDherwig",[("herwigpp_qcd_m3300_3800___Nov28",0.0126252e-6),
		       ]),
	    ]
    if bin==4 and name=="QCDherwig":
      samples=[("QCDherwig",[("herwigpp_qcd_m2800_3300___Nov28",0.0431171e-6),
		       ]),
	    ]
    if bin==5 and name=="QCDherwig":
      samples=[("QCDherwig",[("herwigpp_qcd_m2400_2800___Nov28",0.109694e-6),
		       ]),
	    ]
    if bin==6 and name=="QCDherwig":
      samples=[("QCDherwig",[("herwigpp_qcd_m1900_2400___Nov28",0.588676e-6),
		       ]),
	    ]
    if bin==7 and name=="QCDherwig": 
      samples=[("QCDherwig",[("herwigpp_qcd_m1500_1900___Nov28",2.26487e-6),
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
    	samples=[("QCDmadgraph",[("QCD_HT200to300_RunII_106X_v1",1710000./44805214),
  	  		 ]),
  	      ]
      if bin==2 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT300to500_RunII_106X_v1",347500./48404535),
  	  		 ]),
  	      ]
      if bin==3 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT500to700_RunII_106X_v1",30363.051/46063160),
  	  		 ]),
  	      ]
      if bin==4 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT700to1000_RunII_106X_v1",6428.869/37259115),
  	  		 ]),
  	      ]
      if bin==5 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1000to1500_RunII_106X_v1",1122.659/13511726),
  	  		 ]),
  	      ]
      if bin==6 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1500to2000_RunII_106X_v1",108.163/6059830),
  	  		 ]),
  	      ]
      if bin==7 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT2000toInf_RunII_106X_v1",22.008/3812684),
  	  		 ]),
  	      ]
    if str(year)=="UL16postVFP":
      if bin==1 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT200to300_RunII_106X_v1",1710000./41210455),
  	  		 ]),
  	      ]
      if bin==2 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT300to500_RunII_106X_v1",347500./47426214),
  	  		 ]),
  	      ]
      if bin==3 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT500to700_RunII_106X_v1",30363.051/49068426),
  	  		 ]),
  	      ]
      if bin==4 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT700to1000_RunII_106X_v1",6428.869/38188739),
  	  		 ]),
  	      ]
      if bin==5 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1000to1500_RunII_106X_v1",1122.659/10707004),
  	  		 ]),
  	      ]
      if bin==6 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1500to2000_RunII_106X_v1",108.163/7591790),
  	  		 ]),
  	      ]
      if bin==7 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT2000toInf_RunII_106X_v1",22.008/3620418),
  	  		 ]),
  	      ]
    if str(year)=="UL17":
      if bin==1 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT200to300_RunII_106X_v1",1710000./57721120),
  	  		 ]),
  	      ]
      if bin==2 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT300to500_RunII_106X_v1",347500./57191140),
  	  		 ]),
  	      ]
      if bin==3 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT500to700_RunII_106X_v1",30363.051/9188310),
  	  		 ]),
  	      ]
      if bin==4 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT700to1000_RunII_106X_v1",6428.869/45812757),
  	  		 ]),
  	      ]
      if bin==5 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1000to1500_RunII_106X_v1",1122.659/15346629),
  	  		 ]),
  	      ]
      if bin==6 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1500to2000_RunII_106X_v1",108.163/10598209),
  	  		 ]),
  	      ]
      if bin==7 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT2000toInf_RunII_106X_v1",22.008/5416717),
  	  		 ]),
  	      ]
    if str(year)=="UL18":
      if bin==1 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT200to300_RunII_106X_v1",1710000./22826901),
  	  		 ]),
  	      ]
      if bin==2 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT300to500_RunII_106X_v1",347500./54463611),
  	  		 ]),
  	      ]
      if bin==3 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT500to700_RunII_106X_v1",30363.051/58487165),
  	  		 ]),
  	      ]
      if bin==4 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT700to1000_RunII_106X_v1",6428.869/47703400),
  	  		 ]),
  	      ]
      if bin==5 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1000to1500_RunII_106X_v1",1122.659/15675643),
  	  		 ]),
  	      ]
      if bin==6 and name=="QCDmadgraph":
    	samples=[("QCDmadgraph",[("QCD_HT1500to2000_RunII_106X_v1",108.163/10612885),
  	  		 ]),
  	      ]
      if bin==7 and name=="QCDmadgraph": 
    	samples=[("QCDmadgraph",[("QCD_HT2000toInf_RunII_106X_v1",22.008/4504262),
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
    for j in range((len(massbins)+1)):
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
    if doJER:
     for j in range((len(massbins)+1),len(massbins)+1+2*(len(massbins)+1)*len(JERuncertainties)):
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

    for j in range((len(massbins)+1)+doJER*2*(len(massbins)+1)*len(JERuncertainties)):
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
      if doJER:
       for k in range(len(JERuncertainties)):
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

