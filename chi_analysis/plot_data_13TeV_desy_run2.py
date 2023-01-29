import os, sys
from ROOT import * 
import array
import math

#gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
#gROOT.Reset()
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
gStyle.SetNdivisions(506, "XYZ")
gStyle.SetLegendBorderSize(0)

prefire2016file=TFile.Open("L1prefiring_jetpt_2016BtoH.root")
prefire2016=prefire2016file.Get("L1prefiring_jetpt_2016BtoH")
prefire2017file=TFile.Open("L1prefiring_jetpt_2017BtoF.root")
prefire2017=prefire2017file.Get("L1prefiring_jetpt_2017BtoF")
prefireULfile=TFile.Open("L1PrefiringMaps_UL17andUL16_new.root")
prefireUL16preVFP=prefireULfile.Get("L1prefiring_jetptvseta_UL2016preVFP")
prefireUL16postVFP=prefireULfile.Get("L1prefiring_jetptvseta_UL2016postVFP")
prefireUL17=prefireULfile.Get("L1prefiring_jetptvseta_UL2017BtoF")

def deltaPhi(phi1, phi2):
  deltaphi = phi2 - phi1
  if abs(deltaphi) > math.pi:
    deltaphi = 2 * math.pi - abs(deltaphi)
  return deltaphi

def createPlots(sample,prefix,triggers,massbins,chi_bins):
    files=[]
    print "list files"
    if not ".root" in sample:
      folders=os.listdir(sample)
      for f in folders:
       if "pnfs" in sample:
    	files+=["dcap://dcache-cms-dcap.desy.de/"+sample+"/"+f]
       elif "dijetChi" in f:
        files+=[sample+"/"+f]
    else:
      files=[sample]
    #files=["/afs/desy.de/user/h/hinzmann/uhh106/CMSSW_10_6_26/src/UHH2/dijetangular/dijetChi_tree.root"]
    #files=["/nfs/dust/cms/user/hinzmann/dijetangular/workdir/dijetChiUL16preVFP_RunB_ver2_RunII_106X_v1_0_tree.root"]
    #files=["/nfs/dust/cms/user/hinzmann/dijetangular/dijetChiUL16preVFP_RunB_ver2_RunII_106X_v1-0.root"]
    if correctPrefire and "UL16preVFP" in sample:
      prefiremap=prefireUL16preVFP
    elif correctPrefire and "UL16postVFP" in sample:
      prefiremap=prefireUL16postVFP
    elif correctPrefire and "UL17" in sample:
      prefiremap=prefireUL17
    elif correctPrefire and "2016" in sample:
      prefiremap=prefire2016
    elif correctPrefire and "2017" in sample:
      prefiremap=prefire2017
    else:
      prefiremap=None  

    print files
    
    print triggers

    plots=[]
    for massbin in massbins:
      plots += [TH1F(prefix+'#chi'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';#chi;N',15,1,16)]
    for massbin in massbins:
      plots += [TH1F(prefix+'y_{boost}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';y_{boost};N',12,0,1.2)]
    for massbin in massbins:
      plots += [TH1F(prefix+'p_{T1}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';p_{T1};N',50,0,5000)]
    for massbin in massbins:
      plots += [TH1F(prefix+'p_{T2}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';p_{T2};N',50,0,5000)]
    for massbin in massbins:
      plots += [TH1F(prefix+'y_{1}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';y_{1};N',25,-2.5,2.5)]
    for massbin in massbins:
      plots += [TH1F(prefix+'y_{2}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';y_{2};N',25,-2.5,2.5)]
    for massbin in massbins:
      plots += [TH1F(prefix+'#phi_{1}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';#phi_{1};N',64,-3.1415,3.1415)]
    for massbin in massbins:
      plots += [TH1F(prefix+'#phi_{2}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';#phi_{2};N',64,-3.1415,3.1415)]
    for massbin in massbins:
      plots += [TH1F(prefix+'METsumET'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';missing E_{T} / #sum E_{T};N',50,0,1)]
    for massbin in massbins:
      plots += [TH1F(prefix+'dPtsumPt'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';(p_{T1}-p_{T2})/(p_{T1}+p_{T2});N',25,0,0.5)]
    for massbin in massbins:
      plots += [TH1F(prefix+'#Delta#phi'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';#Delta#phi;N',32,0,3.1415)]
    plots += [TH1F(prefix+'mass',';dijet mass;N',260,0,13000)]
    for c in range(len(chi_bins[0])-1):
      plots += [TH1F(prefix+'mass-reftrig-chi-'+str(chi_bins[0][c]),';dijet mass;N',260,0,13000)]
    plots += [TH1F(prefix+'mass-reftrig',';dijet mass;N',260,0,13000)]
    for t in range(len(triggers)-len(massbins)-1):
      for c in range(len(chi_bins[0])-1):
        plots += [TH1F(prefix+'mass-trig'+"or".join(triggers[len(massbins)+1+t])+'-chi-'+str(chi_bins[0][c]),';dijet mass;N',260,0,13000)]
      plots += [TH1F(prefix+'mass-trig'+"or".join(triggers[len(massbins)+1+t]),';dijet mass;N',260,0,13000)]
    
    for plot in plots:
        plot.Sumw2()

    event_count=0
    for f in files[:]:
     try:
       fil=TFile.Open(f)
       
       events=fil.Get("tree")
       nevents=events.GetEntries()
     except:
       print "error opening", f
       continue
     print f,nevents
     #if event_count>1000: break ###
     trigger_indices_len=0
     for event in events:
	 #if not int(event.EVENT_event)==971086788: continue
         event_count+=1
	 #if event_count>100000 and not "QCD" in prefix: break ###
	 #if event_count>1000: break
         if event_count%10000==1:
	   print "event",event_count
	 if len(event.HLT_isFired)!=trigger_indices_len:
	   trigger_indices={}
	   trigger_indices_len=len(event.HLT_isFired)
	   print "read trigger names", event_count
           for t in range(len(triggers)):
	     trigger_indices[t]=[]
           for a,b in event.HLT_isFired:
             for t in range(len(triggers)):
	       for trig in triggers[t]:
	         if a.startswith(trig):
	           trigger_indices[t]+=[a]
	 jet1=TLorentzVector()
         jet2=TLorentzVector()
         jet1.SetPtEtaPhiM(event.jetAK4_pt1,event.jetAK4_eta1,event.jetAK4_phi1,event.jetAK4_mass1)
         jet2.SetPtEtaPhiM(event.jetAK4_pt2,event.jetAK4_eta2,event.jetAK4_phi2,event.jetAK4_mass2)
         if jet2.Pt()>jet1.Pt():
           jetTmp=jet1
           jet1=jet2
           jet2=jetTmp
         mjj=(jet1+jet2).M()
         chi=math.exp(abs(jet1.Rapidity()-jet2.Rapidity()))
         yboost=abs(jet1.Rapidity()+jet2.Rapidity())/2.
	 weight=1.0
	 if vetoHEM and (event.EVENT_run>=319077) and ((-1.57<jet1.Phi()) and (jet1.Phi()< -0.87) or (-1.57<jet2.Phi()) and (jet2.Phi()< -0.87)): continue
	 if prefiremap:
            if abs(jet1.Eta())>2:
	      weight/=1.-prefiremap.GetBinContent(prefiremap.FindBin(jet1.Eta(),min(499,jet1.Pt())))
	    if abs(jet2.Eta())>2:
	      weight/=1.-prefiremap.GetBinContent(prefiremap.FindBin(jet2.Eta(),min(499,jet2.Pt())))
         if "qcdpy" in sample or "qcdhw" in sample: weight*=event.genWeight
         if mjj<massbins[0][0] or chi>16. or yboost>1.11: continue
	 if mjj>6000: print "found",long(event.EVENT_event), int(event.EVENT_lumiBlock), int(event.EVENT_run), mjj,chi,yboost,jet1.Pt(),jet1.Rapidity(),jet1.Phi(),jet2.Pt(),jet2.Rapidity(),jet2.Phi()
         if jet1.Pt()>13000: continue
         irec=0
	 for massbin in massbins:
	    passedHLT=len(triggers[massbins.index(massbin)])==0
	    for i in trigger_indices[massbins.index(massbin)]:
	      if event.HLT_isFired.find(i)!=event.HLT_isFired.end() and event.HLT_isFired[i]:
	        passedHLT=True
		break
	    if passedHLT and mjj>=massbin[0] and mjj<massbin[1]:
               plots[irec].Fill(chi,weight)
               plots[irec+1*len(massbins)].Fill(yboost,weight)
               plots[irec+2*len(massbins)].Fill(jet1.Pt(),weight)
               plots[irec+3*len(massbins)].Fill(jet2.Pt(),weight)
               plots[irec+4*len(massbins)].Fill(jet1.Rapidity(),weight)
               plots[irec+5*len(massbins)].Fill(jet2.Rapidity(),weight)
               plots[irec+6*len(massbins)].Fill(jet1.Phi(),weight)
               plots[irec+7*len(massbins)].Fill(jet2.Phi(),weight)
               plots[irec+8*len(massbins)].Fill(event.MET_et/event.MET_sumEt,weight)
               plots[irec+9*len(massbins)].Fill((jet1.Pt()-jet2.Pt())/(jet1.Pt()+jet2.Pt()),weight)
               plots[irec+10*len(massbins)].Fill(deltaPhi(jet1.Phi(),jet2.Phi()),weight)
	    irec+=1
         irec=11*len(massbins)
	 plots[irec].Fill(mjj,weight)
	 irec+=1
	 if not trigger_indices.has_key(len(massbins)): continue # for MC
         for i in trigger_indices[len(massbins)]:
          if event.HLT_isFired.find(i)!=event.HLT_isFired.end() and event.HLT_isFired[i]:
            for c in range(len(chi_bins[0])-1):
              if chi>=chi_bins[0][c] and chi<chi_bins[0][c+1]:
                plots[irec].Fill(mjj,weight)
              irec+=1
            plots[irec].Fill(mjj,weight)
            irec+=1
            for t in range(len(triggers)-len(massbins)-1):
              passHLT=False
	      for i in trigger_indices[len(massbins)+1+t]:
                if event.HLT_isFired.find(i)!=event.HLT_isFired.end() and event.HLT_isFired[i]:
                  passHLT=True
                  break
              if passHLT:  
                for c in range(len(chi_bins[0])-1):
                  if chi>=chi_bins[0][c] and chi<chi_bins[0][c+1]:
                    plots[irec].Fill(mjj,weight)
                  irec+=1
                plots[irec].Fill(mjj,weight)
                irec+=1
	      else:
	        irec+=len(chi_bins[0])
            break
    #for plot in plots:
    #  if event_count>0:
    #    plot.Scale(xsec/event_count)
    return plots

if __name__ == '__main__':

    wait=False
    vetoHEM=False
    correctPrefire=False
 
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
	      (6000,7000),
	      (7000,13000),
              ]
 
    samples=[#("datacard_shapelimit13TeV_run2_2016old","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dijet_angular/jobtmpFeb5_data9"),
            ("datacard_shapelimit13TeV_run2_2016","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/data2016dec"),
            ("datacard_shapelimit13TeV_run2_2017","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/data2017dec"),
            ("datacard_shapelimit13TeV_run2_2018","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/data2018dec"),
            ("datacard_shapelimit13TeV_run2_2016_QCDmadgraph","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcd2016"),
            ("datacard_shapelimit13TeV_run2_2017_QCDmadgraph","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcd2017"),
            ("datacard_shapelimit13TeV_run2_2018_QCDmadgraph","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcd2018"),
            ("datacard_shapelimit13TeV_run2_2018_QCDpythia","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcdpy2018"),
            ("datacard_shapelimit13TeV_run2_2018_QCDherwig","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcdhw2018"),
            ("datacard_shapelimit13TeV_run2_2016_SingleMuon","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dataSingleMuon2016"),
            ("datacard_shapelimit13TeV_run2_2017_SingleMuon","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dataSingleMuon2017"),
            ("datacard_shapelimit13TeV_run2_2018_SingleMuon","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dataSingleMuon2018"),
            ("datacard_shapelimit13TeV_run2_UL16preVFP","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dataUL16preVFPoct"),
            ("datacard_shapelimit13TeV_run2_UL16postVFP","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dataUL16postVFPoct"),
            ("datacard_shapelimit13TeV_run2_UL17","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dataUL17oct"),
            ("datacard_shapelimit13TeV_run2_UL18","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dataUL18oct"),
            ("datacard_shapelimit13TeV_run2_UL16preVFP","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dataUL16preVFPoct"),
            ("datacard_shapelimit13TeV_run2_UL16postVFP","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dataUL16postVFPoct"),
            ("datacard_shapelimit13TeV_run2_UL17","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dataUL17oct"),
            ("datacard_shapelimit13TeV_run2_UL18","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/dataUL18oct"),
            ("datacard_shapelimit13TeV_run2_UL16preVFP_QCDmadgraph","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcdUL16preVFPjan2023"),
            ("datacard_shapelimit13TeV_run2_UL16postVFP_QCDmadgraph","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcdUL16postVFPjan2023"),
            ("datacard_shapelimit13TeV_run2_UL17_QCDmadgraph","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcdUL17jan2023"),
            ("datacard_shapelimit13TeV_run2_UL18_QCDmadgraph","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/qcdUL18jan2023"),
            ]

    triggers=[[["HLT_PFHT475","HLT_PFJet260"], #2016
          ["HLT_PFHT475","HLT_PFJet260"],
          ["HLT_PFHT600","HLT_PFHT475","HLT_PFJet320"],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_PFHT475"],
	  ["HLT_PFHT900","HLT_PFHT800","HLT_PFJet450","HLT_PFJet500","HLT_CaloJet500_NoJetID"],
         ],
	  [["HLT_PFHT510","HLT_PFJet260"], #2017
          ["HLT_PFHT590","HLT_PFHT510","HLT_PFJet260"],
          ["HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet320"],
          [],#["HLT_PFHT890","HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet450"],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_PFHT510"],
	  ["HLT_PFHT1050","HLT_PFJet500","HLT_PFJet550","HLT_CaloJet500_NoJetID","HLT_CaloJet550_NoJetID"],
         ],
	  [["HLT_PFHT510","HLT_PFJet260"], #2018
          ["HLT_PFHT590","HLT_PFHT510","HLT_PFJet260"],
          ["HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet320"],
          [],#["HLT_PFHT890","HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet450"],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_PFHT510"],
	  ["HLT_PFHT1050","HLT_PFJet500","HLT_PFJet550","HLT_CaloJet500_NoJetID","HLT_CaloJet550_NoJetID"],
         ],
	  [[], #QCD 2016
          [],
          [],
          [],
          [], 
          [],
          [],
          [],
          [],
          [],
          [],
          [],
         ],
	  [[], #QCD 2017
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
         ],
	  [[], #QCD 2018
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
         ],
	  [[], #QCD Py 2018
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
         ],
	  [[], #QCD Hw 2018
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
         ],
	  [[], #2016 SingleMuon
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_Mu45"],
	  ["HLT_PFHT900","HLT_PFHT800","HLT_PFJet450","HLT_PFJet500","HLT_CaloJet500_NoJetID"],
	  ["HLT_PFHT900"],
	  ["HLT_PFHT800"],
	  ["HLT_PFJet450"],
	  ["HLT_PFJet500"],
	  ["HLT_CaloJet500_NoJetID"],
	  ["HLT_PFHT475"],
          ["HLT_PFHT600"],
          ["HLT_PFHT650"],
	  ["HLT_PFJet260"],
	  ["HLT_PFJet320"],
	  ["HLT_PFJet400"],
         ],
	  [[], #2017 SingleMuon
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_Mu50"],
	  ["HLT_PFHT1050","HLT_PFJet500","HLT_PFJet550","HLT_CaloJet500_NoJetID","HLT_CaloJet550_NoJetID"],
	  ["HLT_PFHT1050"],
	  ["HLT_PFJet500"],
	  ["HLT_PFJet550"],
	  ["HLT_CaloJet500_NoJetID"],
	  ["HLT_CaloJet550_NoJetID"],
	  ["HLT_PFHT510"],
	  ["HLT_PFHT590"],
	  ["HLT_PFHT680"],
	  ["HLT_PFHT780"],
	  ["HLT_PFHT890"],
	  ["HLT_PFJet260"],
	  ["HLT_PFJet320"],
	  ["HLT_PFJet400"],
	  ["HLT_PFJet450"],
         ],
	  [[], #2018 SingleMuon
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_Mu50"],#"HLT_PFHT510",
	  ["HLT_PFHT1050","HLT_PFJet500","HLT_PFJet550","HLT_CaloJet500_NoJetID","HLT_CaloJet550_NoJetID"],
	  ["HLT_PFHT1050"],
	  ["HLT_PFJet500"],
	  ["HLT_PFJet550"],
	  ["HLT_CaloJet500_NoJetID"],
	  ["HLT_CaloJet550_NoJetID"],
	  ["HLT_PFHT510"],
	  ["HLT_PFHT590"],
	  ["HLT_PFHT680"],
	  ["HLT_PFHT780"],
	  ["HLT_PFHT890"],
	  ["HLT_PFJet260"],
	  ["HLT_PFJet320"],
	  ["HLT_PFJet400"],
	  ["HLT_PFJet450"],
         ],
          [["HLT_PFHT475","HLT_PFJet260"], #2016 preVFP
          ["HLT_PFHT475","HLT_PFJet260"],
          ["HLT_PFHT600","HLT_PFHT475","HLT_PFJet320"],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_PFHT475"],
	  ["HLT_PFHT900","HLT_PFHT800","HLT_PFJet450","HLT_PFJet500","HLT_CaloJet500_NoJetID"],
         ],
          [["HLT_PFHT475","HLT_PFJet260"], #2016 postVFP
          ["HLT_PFHT475","HLT_PFJet260"],
          ["HLT_PFHT600","HLT_PFHT475","HLT_PFJet320"],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_PFHT475"],
	  ["HLT_PFHT900","HLT_PFHT800","HLT_PFJet450","HLT_PFJet500","HLT_CaloJet500_NoJetID"],
         ],
	  [["HLT_PFHT510","HLT_PFJet260"], #2017
          ["HLT_PFHT590","HLT_PFHT510","HLT_PFJet260"],
          ["HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet320"],
          [],#["HLT_PFHT890","HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet450"],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_PFHT510"],
	  ["HLT_PFHT1050","HLT_PFJet500","HLT_PFJet550","HLT_CaloJet500_NoJetID","HLT_CaloJet550_NoJetID"],
         ],
	  [["HLT_PFHT510","HLT_PFJet260"], #2018
          ["HLT_PFHT590","HLT_PFHT510","HLT_PFJet260"],
          ["HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet320"],
          [],#["HLT_PFHT890","HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet450"],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_PFHT510"],
	  ["HLT_PFHT1050","HLT_PFJet500","HLT_PFJet550","HLT_CaloJet500_NoJetID","HLT_CaloJet550_NoJetID"],
         ],
          [["HLT_PFHT475","HLT_PFJet260"], #2016 preVFP
          ["HLT_PFHT475","HLT_PFJet260"],
          ["HLT_PFHT600","HLT_PFHT475","HLT_PFJet320"],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_PFHT475"],
	  ["HLT_PFHT900","HLT_PFHT800","HLT_PFJet450","HLT_PFJet500","HLT_CaloJet500_NoJetID"],
         ],
          [["HLT_PFHT475","HLT_PFJet260"], #2016 postVFP
          ["HLT_PFHT475","HLT_PFJet260"],
          ["HLT_PFHT600","HLT_PFHT475","HLT_PFJet320"],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_PFHT475"],
	  ["HLT_PFHT900","HLT_PFHT800","HLT_PFJet450","HLT_PFJet500","HLT_CaloJet500_NoJetID"],
         ],
	  [["HLT_PFHT510","HLT_PFJet260"], #2017
          ["HLT_PFHT590","HLT_PFHT510","HLT_PFJet260"],
          ["HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet320"],
          [],#["HLT_PFHT890","HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet450"],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_PFHT510"],
	  ["HLT_PFHT1050","HLT_PFJet500","HLT_PFJet550","HLT_CaloJet500_NoJetID","HLT_CaloJet550_NoJetID"],
         ],
	  [["HLT_PFHT510","HLT_PFJet260"], #2018
          ["HLT_PFHT590","HLT_PFHT510","HLT_PFJet260"],
          ["HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet320"],
          [],#["HLT_PFHT890","HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet450"],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
	  ["HLT_PFHT510"],
	  ["HLT_PFHT1050","HLT_PFJet500","HLT_PFJet550","HLT_CaloJet500_NoJetID","HLT_CaloJet550_NoJetID"],
         ],
	  [[], #QCD 2016 preVFP
          [],
          [],
          [],
          [], 
          [],
          [],
          [],
          [],
          [],
          [],
          [],
         ],
	  [[], #QCD 2016 postVFP
          [],
          [],
          [],
          [], 
          [],
          [],
          [],
          [],
          [],
          [],
          [],
         ],
	  [[], #QCD 2017
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
         ],
	  [[], #QCD 2018
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
         ],

          ]

    if len(sys.argv)>1:
      samplenum=int(sys.argv[1])
      samples=[samples[samplenum]]
      triggers=[triggers[samplenum]]

    if len(sys.argv)>2:
      folders=os.listdir(samples[0][2])
      folders=[f for f in folders if "dijetChi" in f]
      print len(folders)
      name=sys.argv[2]
      for s in folders[int(sys.argv[2])].split("_"):
        if "HT" in s: name=s+"_"+name; break
      samples=[(samples[0][0],"-"+name,samples[0][2].replace("/pnfs/","dcap://dcache-cms-dcap.desy.de//pnfs/")+"/"+folders[int(sys.argv[2])])]
      print samples
 
    if len(sys.argv)>3:
      if "HEM" in sys.argv[3]:
        vetoHEM=True
      if "L1prefire" in sys.argv[3]:
        correctPrefire=True

    chi_binnings=[]
    for mass_bin in chi_bins:
        chi_binnings+=[array.array('d')]
        for chi_bin in mass_bin:
            chi_binnings[-1].append(chi_bin)
        
    print samples

    for prefix,postfix,files in samples:
      plots=[createPlots(files,prefix,triggers[samples.index((prefix,postfix,files))],massbins,chi_bins)]
      if vetoHEM:
        postfix+="-HEM"
      if correctPrefire:
        postfix+="-L1prefire"

      out=TFile(prefix+postfix + '_chi.root','RECREATE')
      for j in range(len(massbins)):
  	for i in range(1):
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
      for j in range(len(massbins),len(plots[0])):
  	for i in range(1):
          plots[i][j].Write()
  
      for j in range(len(massbins)):
  	for i in range(1):
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
  	legend1=TLegend(0.6,0.6,0.9,0.9,(str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV").replace("<13000",""))
  	legends+=[legend1]
  	legend1.AddEntry(plots[0][j],samples[0][0],"l")
  	legend1.SetTextSize(0.04)
  	legend1.SetFillStyle(0)
  	legend1.Draw("same")

      canvas.SaveAs(prefix+postfix + '_chi.pdf')
      #canvas.SaveAs(prefix+postfix + '_chi.eps')
      #if wait:
      #  os.system("ghostview "+prefix + '_chi.eps')
