import os, sys
from ROOT import * 
from DataFormats.FWLite import Events,Handle
import array, math

#gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
#gROOT.Reset()
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

def createPlots(sample,prefix,weightname,massbins):
    files=[]
    print "list files"
    if sample.endswith(".txt"):
        filelist=open(sample)
	for line in filelist.readlines():
	    if ".root" in line:
	        files+=[line.strip()]
    else:
      if "alp" in prefix:
        foldersM=os.listdir("/nfs/dust/cms/user/hinzmann/dijetangular")
	for folderM in foldersM:
	 if not os.path.exists("/nfs/dust/cms/user/hinzmann/dijetangular/"+folderM+"/"+sample): continue
	 if not "newsamples0" in folderM: continue
         folders=os.listdir("/nfs/dust/cms/user/hinzmann/dijetangular/"+folderM+"/"+sample)
	 for folder in folders:
	  if sample in folder and ".root" in folder:
            files+=["file:///nfs/dust/cms/user/hinzmann/dijetangular/"+folderM+"/"+sample+"/"+folder]
	    #break
      elif "tripleG" in prefix:
        foldersM=os.listdir("/nfs/dust/cms/user/hinzmann/dijetangular")
	for folderM in foldersM:
	 if not os.path.exists("/nfs/dust/cms/user/hinzmann/dijetangular/"+folderM+"/"+sample): continue
	 if not "mg" in folderM: continue
         folders=os.listdir("/nfs/dust/cms/user/hinzmann/dijetangular/"+folderM+"/"+sample)
  	 for folder in folders:
	  if sample in folder and ".root" in folder:
            files+=["file:///nfs/dust/cms/user/hinzmann/dijetangular/"+folderM+"/"+sample+"/"+folder]
	    #break
      elif "DM" in prefix:
        folders=os.listdir("/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dijet_angular/dm/")
	for folder in folders:
	  if sample in folder and ".root" in folder:
            files+=["dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dijet_angular/dm/"+folder]
	    #break
      else:
        folders=os.listdir("/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dijet_angular/")
	for folder in folders:
	  if sample in folder:
            files+=["dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dijet_angular/"+folder+"/GEN.root"]
	    #break

    print files
    prunedgenjets_handle=Handle("std::vector<reco::GenJet>")
    prunedgenjets_label="ak4GenJets"

    plots=[]
    for massbin in massbins:
      plots += [TH1F(prefix+'#chi'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';#chi;N',15,1,16)]
      #plots += [TH1F(prefix+'y_{boost}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';y_{boost};N',20,0,2)]
    
    for plot in plots:
        plot.Sumw2()

    event_count=0
    print "open chain"
    events=TChain('Events')
    for f in files[:]:
      events.Add(f)
    
    nevents=events.GetEntries()
    print sample,nevents,weightname
    event_count=0
    acceptance=0
    sumweights=0
    firstxsec=0
    for event in events:
         if "DM" in prefix or "alp" in prefix or "tripleG" in prefix:
	  xsec=event.LHEEventProduct_externalLHEProducer__GEN.product().originalXWGTUP()*1000. # convert to pb
	  if firstxsec>0 and xsec!=firstxsec:
	     print "inconsistent xsec",firstxsec,xsec
	     inconsistent
	  #weights=[(w.id,w.wgt) for w in event.LHEEventProduct_externalLHEProducer__GEN.product().weights()]
	  #print weights
	  try:
	   weight=[w.wgt for w in event.LHEEventProduct_externalLHEProducer__GEN.product().weights() if weightname==w.id][0]/event.LHEEventProduct_externalLHEProducer__GEN.product().weights()[0].wgt
	  except:
	   print "error reading weight"
 	   break
         else:
	  xsec=weightname
	  weight=1.
	 event_count+=1
	 #if event_count>1000: break
         if event_count%10000==1: print "event",event_count
	 sumweights+=weight
         jet1=TLorentzVector()
         jet2=TLorentzVector()
	 jets=event.recoGenJets_ak4GenJets__GEN.product()
	 if len(jets)<2: continue
	 j1=jets[0]
	 j2=jets[1]
         jet1.SetPtEtaPhiM(j1.pt(),j1.eta(),j1.phi(),j1.mass())
         jet2.SetPtEtaPhiM(j2.pt(),j2.eta(),j2.phi(),j2.mass())
         mjj=(jet1+jet2).M()
         chi=math.exp(abs(jet1.Rapidity()-jet2.Rapidity()))
         yboost=abs(jet1.Rapidity()+jet2.Rapidity())/2.
         if mjj<1000 or chi>16. or yboost>1.11: continue
         if mjj>2400: acceptance+=weight
         irec=0
	 for massbin in massbins:
            if yboost<1.11 and mjj>=massbin[0] and mjj<massbin[1]:
               plots[irec].Fill(chi,weight)
	    irec+=1
    print sample,weightname,"acceptance",acceptance/sumweights, "xsec",xsec*sumweights/event_count
    for plot in plots:
      if event_count>0:
        plot.Scale(xsec/event_count)
    return plots

if __name__ == '__main__':

    wait=False
 
    prefix="datacard_shapelimit13TeV_GEN-QCD-run2"
    if len(sys.argv)>1:
      if len(sys.argv)>2:
       point=sys.argv[1]
       nxsec=int(sys.argv[2])
       if "Vector" in point:
         weights=['gdmv_1p0_gdma_0_gv_0p01_ga_0', 'gdmv_1p0_gdma_0_gv_0p05_ga_0', 'gdmv_1p0_gdma_0_gv_0p1_ga_0', 'gdmv_1p0_gdma_0_gv_0p2_ga_0', 'gdmv_1p0_gdma_0_gv_0p25_ga_0', 'gdmv_1p0_gdma_0_gv_0p3_ga_0', 'gdmv_1p0_gdma_0_gv_0p5_ga_0', 'gdmv_1p0_gdma_0_gv_0p75_ga_0', 'gdmv_1p0_gdma_0_gv_1_ga_0', 'gdmv_1p0_gdma_0_gv_1p5_ga_0', 'gdmv_1p0_gdma_0_gv_2p0_ga_0', 'gdmv_1p0_gdma_0_gv_2p5_ga_0', 'gdmv_1p0_gdma_0_gv_3p0_ga_0']
         prefix="datacard_shapelimit13TeV_DM"+point+"_"+weights[nxsec]+"-run2"
       elif "Axial" in point:
         weights=['gdmv_0_gdma_1p0_gv_0_ga_0p01', 'gdmv_0_gdma_1p0_gv_0_ga_0p05', 'gdmv_0_gdma_1p0_gv_0_ga_0p1', 'gdmv_0_gdma_1p0_gv_0_ga_0p2', 'gdmv_0_gdma_1p0_gv_0_ga_0p25', 'gdmv_0_gdma_1p0_gv_0_ga_0p3', 'gdmv_0_gdma_1p0_gv_0_ga_0p5', 'gdmv_0_gdma_1p0_gv_0_ga_0p75', 'gdmv_0_gdma_1p0_gv_0_ga_1', 'gdmv_0_gdma_1p0_gv_0_ga_1p5', 'gdmv_0_gdma_1p0_gv_0_ga_2p0', 'gdmv_0_gdma_1p0_gv_0_ga_2p5', 'gdmv_0_gdma_1p0_gv_0_ga_3p0']
         prefix="datacard_shapelimit13TeV_DM"+point+"_"+weights[nxsec]+"-run2"
       elif "alp" in point:
         weights=['fa1000','fa1500','fa2000','fa2500','fa3000','fa3500','fa4000','fa4500','fa5000','fa50000']
         prefix="datacard_shapelimit13TeV_"+point+"_"+weights[nxsec]+"-run2"
       elif "tripleG" in point:
         weights=["CG0p1","CG0p05","CG0p04","CG0p03","CG0p025","CG0p02","CG0p015","CG0p01","CG0p0075","CG0p005","CG0p0025","CG0p0"]
         prefix="datacard_shapelimit13TeV_"+point+"_"+weights[nxsec]+"-run2"
      else:
       prefix="datacard_shapelimit13TeV_GENnp-"+sys.argv[1]+"-run2"
 
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
	      (6000,7000),
	      (7000,13000),
              ]
 
    samples3=[("QCD",[("pythia8_ci_m1000_1500_50000_1_0_0_13TeV_Nov14",3.769e-05),
		       ("pythia8_ci_m1500_1900_50000_1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_50000_1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_50000_1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_50000_1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_50000_1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_50000_1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_50000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
	    ]
    samples2=[("QCDNonPert",[("pythia8_ciNonPert_m1000_1500_50000_1_0_0_13TeV_Nov14",3.769e-05),
		       ("pythia8_ciNonPert_m1500_1900_50000_1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ciNonPert_m1900_2400_50000_1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ciNonPert_m2400_2800_50000_1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ciNonPert_m2800_3300_50000_1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ciNonPert_m3300_3800_50000_1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ciNonPert_m3800_4300_50000_1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ciNonPert_m4300_13000_50000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
	    ]
    samples2=[("QCDHerwig",[("herwigpp_qcd_m1000_1500___Nov28",3.769e-05),
		       ("herwigpp_qcd_m1500_1900___Nov28",3.307e-06),
		       ("herwigpp_qcd_m1900_2400___Nov28",8.836e-07),
		       ("herwigpp_qcd_m2400_2800___Nov28",1.649e-07),
		       ("herwigpp_qcd_m2800_3300___Nov28",6.446e-08),
		       ("herwigpp_qcd_m3300_3800___Nov28",1.863e-08),
		       ("herwigpp_qcd_m3800_4300___Nov28",5.867e-09),
		       ("herwigpp_qcd_m4300_13000___Nov28",3.507e-09),
		       ]),
	    ]
    samples2=[("QCDHerwigNonPert",[("herwigpp_qcdNonPert_m1000_1500___Nov28",3.769e-05),
		       ("herwigpp_qcdNonPert_m1500_1900___Nov28",3.307e-06),
		       ("herwigpp_qcdNonPert_m1900_2400___Nov28",8.836e-07),
		       ("herwigpp_qcdNonPert_m2400_2800___Nov28",1.649e-07),
		       ("herwigpp_qcdNonPert_m2800_3300___Nov28",6.446e-08),
		       ("herwigpp_qcdNonPert_m3300_3800___Nov28",1.863e-08),
		       ("herwigpp_qcdNonPert_m3800_4300___Nov28",5.867e-09),
		       ("herwigpp_qcdNonPert_m4300_13000___Nov28",3.507e-09),
		       ]),
	    ]
    samples=[("QCDCIplusLL8000",[("pythia8_ci_m1500_1900_8000_1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_8000_1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_8000_1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_8000_1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_8000_1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_8000_1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_8000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIplusLL9000",[("pythia8_ci_m1500_1900_9000_1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_9000_1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_9000_1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_9000_1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_9000_1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_9000_1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_9000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIplusLL10000",[("pythia8_ci_m1500_1900_10000_1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_10000_1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_10000_1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_10000_1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_10000_1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_10000_1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_10000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIplusLL11000",[("pythia8_ci_m1500_1900_11000_1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_11000_1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_11000_1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_11000_1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_11000_1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_11000_1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_11000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIplusLL12000",[("pythia8_ci_m1500_1900_12000_1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_12000_1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_12000_1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_12000_1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_12000_1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_12000_1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_12000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIplusLL13000",[("pythia8_ci_m1500_1900_13000_1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_13000_1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_13000_1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_13000_1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_13000_1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_13000_1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_13000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIplusLL14000",[("pythia8_ci_m1500_1900_14000_1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_14000_1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_14000_1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_14000_1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_14000_1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_14000_1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_14000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIplusLL16000",[("pythia8_ci_m1500_1900_16000_1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_16000_1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_16000_1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_16000_1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_16000_1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_16000_1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_16000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIplusLL18000",[("pythia8_ci_m1500_1900_18000_1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_18000_1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_18000_1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_18000_1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_18000_1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_18000_1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_18000_1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ]
    samples+=[("QCDCIminusLL8000",[("pythia8_ci_m1500_1900_8000_-1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_8000_-1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_8000_-1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_8000_-1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_8000_-1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_8000_-1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_8000_-1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIminusLL9000",[("pythia8_ci_m1500_1900_9000_-1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_9000_-1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_9000_-1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_9000_-1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_9000_-1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_9000_-1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_9000_-1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIminusLL10000",[("pythia8_ci_m1500_1900_10000_-1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_10000_-1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_10000_-1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_10000_-1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_10000_-1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_10000_-1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_10000_-1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIminusLL11000",[("pythia8_ci_m1500_1900_11000_-1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_11000_-1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_11000_-1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_11000_-1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_11000_-1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_11000_-1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_11000_-1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIminusLL12000",[("pythia8_ci_m1500_1900_12000_-1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_12000_-1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_12000_-1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_12000_-1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_12000_-1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_12000_-1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_12000_-1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIminusLL13000",[("pythia8_ci_m1500_1900_13000_-1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_13000_-1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_13000_-1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_13000_-1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_13000_-1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_13000_-1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_13000_-1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIminusLL14000",[("pythia8_ci_m1500_1900_14000_-1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_14000_-1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_14000_-1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_14000_-1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_14000_-1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_14000_-1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_14000_-1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIminusLL16000",[("pythia8_ci_m1500_1900_16000_-1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_16000_-1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_16000_-1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_16000_-1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_16000_-1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_16000_-1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_16000_-1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDCIminusLL18000",[("pythia8_ci_m1500_1900_18000_-1_0_0_13TeV_Nov14",3.307e-06),
		       ("pythia8_ci_m1900_2400_18000_-1_0_0_13TeV_Nov14",8.836e-07),
		       ("pythia8_ci_m2400_2800_18000_-1_0_0_13TeV_Nov14",1.649e-07),
		       ("pythia8_ci_m2800_3300_18000_-1_0_0_13TeV_Nov14",6.446e-08),
		       ("pythia8_ci_m3300_3800_18000_-1_0_0_13TeV_Nov14",1.863e-08),
		       ("pythia8_ci_m3800_4300_18000_-1_0_0_13TeV_Nov14",5.867e-09),
		       ("pythia8_ci_m4300_13000_18000_-1_0_0_13TeV_Nov14",3.507e-09),
		       ]),
             ]
    samples+=[("QCDADD6000",[("pythia8_add_m1500_1900_6000_0_0_0_1_13TeV_Nov14",3.307e-06),
		       ("pythia8_add_m1900_2400_6000_0_0_0_1_13TeV_Nov14",8.836e-07),
		       ("pythia8_add_m2400_2800_6000_0_0_0_1_13TeV_Nov14",1.649e-07),
		       ("pythia8_add_m2800_3300_6000_0_0_0_1_13TeV_Nov14",6.446e-08),
		       ("pythia8_add_m3300_3800_6000_0_0_0_1_13TeV_Nov14",1.863e-08),
		       ("pythia8_add_m3800_4300_6000_0_0_0_1_13TeV_Nov14",5.867e-09),
		       ("pythia8_add_m4300_13000_6000_0_0_0_1_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDADD7000",[("pythia8_add_m1500_1900_7000_0_0_0_1_13TeV_Nov14",3.307e-06),
		       ("pythia8_add_m1900_2400_7000_0_0_0_1_13TeV_Nov14",8.836e-07),
		       ("pythia8_add_m2400_2800_7000_0_0_0_1_13TeV_Nov14",1.649e-07),
		       ("pythia8_add_m2800_3300_7000_0_0_0_1_13TeV_Nov14",6.446e-08),
		       ("pythia8_add_m3300_3800_7000_0_0_0_1_13TeV_Nov14",1.863e-08),
		       ("pythia8_add_m3800_4300_7000_0_0_0_1_13TeV_Nov14",5.867e-09),
		       ("pythia8_add_m4300_13000_7000_0_0_0_1_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDADD8000",[("pythia8_add_m1500_1900_8000_0_0_0_1_13TeV_Nov14",3.307e-06),
		       ("pythia8_add_m1900_2400_8000_0_0_0_1_13TeV_Nov14",8.836e-07),
		       ("pythia8_add_m2400_2800_8000_0_0_0_1_13TeV_Nov14",1.649e-07),
		       ("pythia8_add_m2800_3300_8000_0_0_0_1_13TeV_Nov14",6.446e-08),
		       ("pythia8_add_m3300_3800_8000_0_0_0_1_13TeV_Nov14",1.863e-08),
		       ("pythia8_add_m3800_4300_8000_0_0_0_1_13TeV_Nov14",5.867e-09),
		       ("pythia8_add_m4300_13000_8000_0_0_0_1_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDADD9000",[("pythia8_add_m1500_1900_9000_0_0_0_1_13TeV_Nov14",3.307e-06),
		       ("pythia8_add_m1900_2400_9000_0_0_0_1_13TeV_Nov14",8.836e-07),
		       ("pythia8_add_m2400_2800_9000_0_0_0_1_13TeV_Nov14",1.649e-07),
		       ("pythia8_add_m2800_3300_9000_0_0_0_1_13TeV_Nov14",6.446e-08),
		       ("pythia8_add_m3300_3800_9000_0_0_0_1_13TeV_Nov14",1.863e-08),
		       ("pythia8_add_m3800_4300_9000_0_0_0_1_13TeV_Nov14",5.867e-09),
		       ("pythia8_add_m4300_13000_9000_0_0_0_1_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDADD10000",[("pythia8_add_m1500_1900_10000_0_0_0_1_13TeV_Nov14",3.307e-06),
		       ("pythia8_add_m1900_2400_10000_0_0_0_1_13TeV_Nov14",8.836e-07),
		       ("pythia8_add_m2400_2800_10000_0_0_0_1_13TeV_Nov14",1.649e-07),
		       ("pythia8_add_m2800_3300_10000_0_0_0_1_13TeV_Nov14",6.446e-08),
		       ("pythia8_add_m3300_3800_10000_0_0_0_1_13TeV_Nov14",1.863e-08),
		       ("pythia8_add_m3800_4300_10000_0_0_0_1_13TeV_Nov14",5.867e-09),
		       ("pythia8_add_m4300_13000_10000_0_0_0_1_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDADD11000",[("pythia8_add_m1500_1900_11000_0_0_0_1_13TeV_Nov14",3.307e-06),
		       ("pythia8_add_m1900_2400_11000_0_0_0_1_13TeV_Nov14",8.836e-07),
		       ("pythia8_add_m2400_2800_11000_0_0_0_1_13TeV_Nov14",1.649e-07),
		       ("pythia8_add_m2800_3300_11000_0_0_0_1_13TeV_Nov14",6.446e-08),
		       ("pythia8_add_m3300_3800_11000_0_0_0_1_13TeV_Nov14",1.863e-08),
		       ("pythia8_add_m3800_4300_11000_0_0_0_1_13TeV_Nov14",5.867e-09),
		       ("pythia8_add_m4300_13000_11000_0_0_0_1_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDADD12000",[("pythia8_add_m1500_1900_12000_0_0_0_1_13TeV_Nov14",3.307e-06),
		       ("pythia8_add_m1900_2400_12000_0_0_0_1_13TeV_Nov14",8.836e-07),
		       ("pythia8_add_m2400_2800_12000_0_0_0_1_13TeV_Nov14",1.649e-07),
		       ("pythia8_add_m2800_3300_12000_0_0_0_1_13TeV_Nov14",6.446e-08),
		       ("pythia8_add_m3300_3800_12000_0_0_0_1_13TeV_Nov14",1.863e-08),
		       ("pythia8_add_m3800_4300_12000_0_0_0_1_13TeV_Nov14",5.867e-09),
		       ("pythia8_add_m4300_13000_12000_0_0_0_1_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDADD13000",[("pythia8_add_m1500_1900_13000_0_0_0_1_13TeV_Nov14",3.307e-06),
		       ("pythia8_add_m1900_2400_13000_0_0_0_1_13TeV_Nov14",8.836e-07),
		       ("pythia8_add_m2400_2800_13000_0_0_0_1_13TeV_Nov14",1.649e-07),
		       ("pythia8_add_m2800_3300_13000_0_0_0_1_13TeV_Nov14",6.446e-08),
		       ("pythia8_add_m3300_3800_13000_0_0_0_1_13TeV_Nov14",1.863e-08),
		       ("pythia8_add_m3800_4300_13000_0_0_0_1_13TeV_Nov14",5.867e-09),
		       ("pythia8_add_m4300_13000_13000_0_0_0_1_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDADD14000",[("pythia8_add_m1500_1900_14000_0_0_0_1_13TeV_Nov14",3.307e-06),
		       ("pythia8_add_m1900_2400_14000_0_0_0_1_13TeV_Nov14",8.836e-07),
		       ("pythia8_add_m2400_2800_14000_0_0_0_1_13TeV_Nov14",1.649e-07),
		       ("pythia8_add_m2800_3300_14000_0_0_0_1_13TeV_Nov14",6.446e-08),
		       ("pythia8_add_m3300_3800_14000_0_0_0_1_13TeV_Nov14",1.863e-08),
		       ("pythia8_add_m3800_4300_14000_0_0_0_1_13TeV_Nov14",5.867e-09),
		       ("pythia8_add_m4300_13000_14000_0_0_0_1_13TeV_Nov14",3.507e-09),
		       ]),
             ("QCDADD15000",[("pythia8_add_m1500_1900_15000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m1900_2400_15000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2400_2800_15000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2800_3300_15000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3300_3800_15000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3800_4300_15000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m4300_5200_15000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m5200_13000_15000_0_0_0_1_13TeV_Nov14",1),
		       ]),
             ("QCDADD16000",[("pythia8_add_m1500_1900_16000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m1900_2400_16000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2400_2800_16000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2800_3300_16000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3300_3800_16000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3800_4300_16000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m4300_5200_16000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m5200_13000_16000_0_0_0_1_13TeV_Nov14",1),
		       ]),
             ("QCDADD17000",[("pythia8_add_m1500_1900_17000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m1900_2400_17000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2400_2800_17000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2800_3300_17000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3300_3800_17000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3800_4300_17000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m4300_5200_17000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m5200_13000_17000_0_0_0_1_13TeV_Nov14",1),
		       ]),
             ("QCDADD18000",[("pythia8_add_m1500_1900_18000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m1900_2400_18000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2400_2800_18000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2800_3300_18000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3300_3800_18000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3800_4300_18000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m4300_5200_18000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m5200_13000_18000_0_0_0_1_13TeV_Nov14",1),
		       ]),
             ("QCDADD19000",[("pythia8_add_m1500_1900_19000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m1900_2400_19000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2400_2800_19000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2800_3300_19000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3300_3800_19000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3800_4300_19000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m4300_5200_19000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m5200_13000_19000_0_0_0_1_13TeV_Nov14",1),
		       ]),
             ("QCDADD20000",[("pythia8_add_m1500_1900_20000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m1900_2400_20000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2400_2800_20000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2800_3300_20000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3300_3800_20000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3800_4300_20000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m4300_5200_20000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m5200_13000_20000_0_0_0_1_13TeV_Nov14",1),
		       ]),
             ("QCDADD21000",[("pythia8_add_m1500_1900_21000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m1900_2400_21000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2400_2800_21000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2800_3300_21000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3300_3800_21000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3800_4300_21000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m4300_5200_21000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m5200_13000_21000_0_0_0_1_13TeV_Nov14",1),
		       ]),
             ("QCDADD22000",[("pythia8_add_m1500_1900_22000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m1900_2400_22000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2400_2800_22000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m2800_3300_22000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3300_3800_22000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m3800_4300_22000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m4300_5200_22000_0_0_0_1_13TeV_Nov14",1),
		       ("pythia8_add_m5200_13000_22000_0_0_0_1_13TeV_Nov14",1),
		       ]),
             ]
    
    if "np" in prefix:
       samples=[samples[int(prefix.split("-")[1])]]
    elif "DM" in prefix:
       samples=[("DM"+point+"_"+weights[nxsec],[(point,weights[nxsec])])]
    elif "alp" in prefix or "tripleG" in prefix:
       samples=[(""+point+"_"+weights[nxsec],[(point+"_Mjj200to1000",weights[nxsec]),
                                              (point+"_Mjj1000to2000",weights[nxsec]),
                                              (point+"_Mjj2000to4000",weights[nxsec]),
                                              (point+"_Mjj4000toInf",weights[nxsec]),
                                              ]),
	       ]
    elif "QCD" in prefix:
       samples=samples3
    
    #xsecs=eval(open("xsecs_13TeV.txt").readline())
    #print xsecs

    chi_binnings=[]
    for mass_bin in chi_bins:
        chi_binnings+=[array.array('d')]
        for chi_bin in mass_bin:
            chi_binnings[-1].append(chi_bin)
        
#    if len(sys.argv)>1:
#        newsamples=[]
#        for sample in samples:
#          found=False
#	  for arg in sys.argv:
#	    if sample[0]==arg or sample[0]=="QCD":
#	        newsamples+=[sample]
#		break
#	samples=newsamples
#	if samples[-1][0]=="QCD":
    #        prefix+="_"+samples[-1][0]
    #    else:
    #	    prefix+="_"+samples[-1][0].replace("QCD","")
  
    print prefix, samples

    plots=[]
    for name,files in samples:
      plots+=[[]]
      i=0
      for filename,xsec in files:
        i+=1
	if "DM" in prefix or "alp" in prefix or "tripleG" in prefix:
          ps=createPlots(filename,name,xsec,massbins)
	#elif ("alp" in prefix or "tripleG" in prefix) and "HT2000" in filename:
        #  ps=createPlots(filename,name,xsec,massbins[:7]) # For mass bins < 4800 GeV
	#elif ("alp" in prefix or "tripleG" in prefix) and "HT4000" in filename:
        #  ps=createPlots(filename,name,xsec,massbins[7:]) # For mass bins >= 4800 GeV
        else:
	  ps=createPlots(filename,name,float(xsecs[filename]),massbins)
        if i==1: # or "alp" in prefix or "tripleG" in prefix
          plots[-1]+=ps
	else:
	  for i in range(len(plots[-1])):
            plots[-1][i].Add(ps[i])

    out=TFile(prefix + '_chi.root','RECREATE')
    for j in range(len(massbins)):
      for i in range(len(samples)):
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

    for j in range(len(massbins)):
      for i in range(len(samples)):
        if plots[i][j].Integral()>0:
          plots[i][j].Scale(1./plots[i][j].Integral())
        for b in range(plots[i][j].GetXaxis().GetNbins()):
          plots[i][j].SetBinContent(b+1,plots[i][j].GetBinContent(b+1)/plots[i][j].GetBinWidth(b+1))
          plots[i][j].SetBinError(b+1,plots[i][j].GetBinError(b+1)/plots[i][j].GetBinWidth(b+1))
        plots[i][j].GetYaxis().SetRangeUser(0,0.2)

    canvas = TCanvas("","",0,0,800,600)
    canvas.Divide(4,3)

    legends=[]
    for j in range(len(massbins)):
      canvas.cd(j+1)
      plots[0][j].Draw("he")
      print "number of events passed:",plots[0][j].GetEntries()
      legend1=TLegend(0.6,0.6,0.9,0.9,(str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV").replace("4200<m_{jj}<13000","m_{jj}>4200"))
      legends+=[legend1]
      legend1.AddEntry(plots[0][j],samples[0][0],"l")
      for i in range(1,len(samples)):
        plots[i][j].SetLineColor(i+2)
        plots[i][j].Draw("hesame")
        legend1.AddEntry(plots[i][j],samples[i][0],"l")
      legend1.SetTextSize(0.04)
      legend1.SetFillStyle(0)
      legend1.Draw("same")

    canvas.SaveAs(prefix + '_chi.pdf')
    #canvas.SaveAs(prefix + '_chi.eps')
    #if wait:
    #    os.system("ghostview "+prefix + '_chi.eps')

