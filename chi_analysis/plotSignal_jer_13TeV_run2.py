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

if doJER:
  JERuncertainties=["JER1","JER2"]
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

def createPlots(sample,prefix,xsec,massbins,year):
    files=[]
    if sample.endswith(".txt"):
        filelist=open(sample)
	for line in filelist.readlines():
	    if ".root" in line:
	        files+=[line.strip()]
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
      #plots += [TH1F(prefix+'y_{boost}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';y_{boost};N',20,0,2)]
    if doJER:
      for source in JERuncertainties:
       for massbin in massbins:
        plots += [TH1F(prefix+'#chi'+str(massbin).strip("()").replace(',',"_").replace(' ',"")+"_"+source+"Up",';#chi;N',15,1,16)]
       for massbin in massbins:
        plots += [TH1F(prefix+'#chi'+str(massbin).strip("()").replace(',',"_").replace(' ',"")+"_"+source+"Down",';#chi;N',15,1,16)]
    
    for plot in plots:
        plot.Sumw2()
	print plot.GetName()

    event_count=0
    events=TChain('Events')
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

       if len(event.recoGenJets_ak4GenJets__GEN.product())<2: continue
       if event.recoGenJets_ak4GenJets__GEN.product()[0].pt()<100 or event.recoGenJets_ak4GenJets__GEN.product()[1].pt()<100 or abs(event.recoGenJets_ak4GenJets__GEN.product()[0].eta())>3 or abs(event.recoGenJets_ak4GenJets__GEN.product()[1].eta())>3: continue
       
       irec=0
       jet1=TLorentzVector()
       jet2=TLorentzVector()
       jet1.SetPtEtaPhiM(event.recoGenJets_ak4GenJets__GEN.product()[0].pt(),event.recoGenJets_ak4GenJets__GEN.product()[0].eta(),event.recoGenJets_ak4GenJets__GEN.product()[0].phi(),event.recoGenJets_ak4GenJets__GEN.product()[0].mass())
       jet2.SetPtEtaPhiM(event.recoGenJets_ak4GenJets__GEN.product()[1].pt(),event.recoGenJets_ak4GenJets__GEN.product()[1].eta(),event.recoGenJets_ak4GenJets__GEN.product()[1].phi(),event.recoGenJets_ak4GenJets__GEN.product()[1].mass())
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
	   if scale!=1 and (("JER1" in scale and abs(jet1.Eta())<=1.93) or ("JER2" in scale and abs(jet1.Eta())>1.93)):
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
	   if scale!=1 and (("JER1" in scale and abs(jet2.Eta())<=1.93) or ("JER2" in scale and abs(jet2.Eta())>1.93)):
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

    #jets="recoGenJets_ak4GenJets__GEN.obj"
    #yboost='abs('+jets+'[0].y()+'+jets+'[1].y())/2.'
    #chi='exp(abs('+jets+'[0].y()-'+jets+'[1].y()))'
    #mass='sqrt(pow('+jets+'[0].energy()+'+jets+'[1].energy(),2)-pow('+jets+'[0].px()+'+jets+'[1].px(),2)-pow('+jets+'[0].py()+'+jets+'[1].py(),2)-pow('+jets+'[0].pz()+'+jets+'[1].pz(),2))'
    #for massbin in massbins:
    #  events.Project(prefix+'#chi'+str(massbin).strip("()").replace(',',"_").replace(' ',""),chi,'('+yboost+'<1.11)*('+mass+'>='+str(massbin[0])+')*('+mass+'<'+str(massbin[1])+')')
    #  #events.Project(prefix+'y_{boost}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),yboost,'('+chi+'<16)*('+mass+'>='+str(massbin[0])+')*('+mass+'<='+str(massbin[1])+')')
    for plot in plots:
      if nevents>0:
        plot.Scale(xsec/nevents)
    return plots

if __name__ == '__main__':
    sets=[]
    i=0
    for name in ["QCD","QCDCIplusLL10000"]:
      for year in [2016,2017,2018]:
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
    xsecs=eval(open("xsecs_13TeV.txt").readline())
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
        ps=createPlots(filename,name,float(xsecs[filename]),massbins,str(year))
        if i==1:
          plots[-1]+=ps
	else:
	  for j in range(len(plots[-1])):
            plots[-1][j].Add(ps[j])

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
    if doJER:
     for j in range(len(massbins),len(massbins)+2*len(massbins)*len(JERuncertainties)):
      for i in range(len(samples)):
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

    for j in range(len(massbins)+doJER*2*len(massbins)*len(JERuncertainties)):
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
        plots[0][j+len(massbins)+2*len(massbins)*k].SetLineColor(colors[k+1])
        plots[0][j+len(massbins)+2*len(massbins)*k].SetLineStyle(2)
	plots[0][j+len(massbins)+2*len(massbins)*k].Divide(plots[0][j+len(massbins)+2*len(massbins)*k],plots[0][j])
        for b in range(plots[0][j+len(massbins)+2*len(massbins)*k].GetNbinsX()):
	    plots[0][j+len(massbins)+2*len(massbins)*k].SetBinError(b+1,0)
        plots[0][j+len(massbins)+2*len(massbins)*k].Draw("hesame")
        legend1.AddEntry(plots[0][j+len(massbins)+2*len(massbins)*k],plots[0][j+len(massbins)+2*len(massbins)*k].GetName().split("_")[-2],"l")
        plots[0][j+2*len(massbins)+2*len(massbins)*k].SetLineColor(colors[k+1])
        plots[0][j+2*len(massbins)+2*len(massbins)*k].SetLineStyle(3)
	plots[0][j+2*len(massbins)+2*len(massbins)*k].Divide(plots[0][j+2*len(massbins)+2*len(massbins)*k],plots[0][j])
        for b in range(plots[0][j+2*len(massbins)+2*len(massbins)*k].GetNbinsX()):
	    plots[0][j+2*len(massbins)+2*len(massbins)*k].SetBinError(b+1,0)
        plots[0][j+2*len(massbins)+2*len(massbins)*k].Draw("hesame")
        legend1.AddEntry(plots[0][j+2*len(massbins)+2*len(massbins)*k],plots[0][j+2*len(massbins)+2*len(massbins)*k].GetName().split("_")[-2],"l")
      legend1.SetTextSize(0.04)
      legend1.SetFillStyle(0)
      legend1.Draw("same")

    canvas.SaveAs(prefix + '_chi.pdf')
    canvas.SaveAs(prefix + '_chi.eps')
    if wait:
        os.system("ghostview "+prefix + '_chi.eps')

