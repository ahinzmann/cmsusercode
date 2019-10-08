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

def createPlots(sample,prefix,triggers,massbins):
    files=[]
    print "list files"
    if not ".root" in sample:
      folders=os.listdir(sample)
      for f in folders:
    	files+=["dcap://dcache-cms-dcap.desy.de/"+sample+"/"+f]
    else:
      files=[sample]

    print files
    
    print triggers

    plots=[]
    for massbin in massbins:
      plots += [TH1F(prefix+'#chi'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';#chi;N',15,1,16)]
    
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
     #if event_count>100000: break ###
     for event in events:
         event_count+=1
	 #if event_count>100000: break ###
         if event_count%10000==1:
	   print "event",event_count
	   trigger_indices={}
           for massbin in massbins:
	     trigger_indices[massbin]=[]
           for a,b in event.HLT_isFired:
             for massbin in massbins:
	       for trig in triggers[massbins.index(massbin)]:
	         if a.startswith(trig):
	           trigger_indices[massbin]+=[a]
	 jet1=TLorentzVector()
         jet2=TLorentzVector()
         jet1.SetPtEtaPhiM(event.jetAK4_pt1,event.jetAK4_eta1,event.jetAK4_phi1,event.jetAK4_mass1)
         jet2.SetPtEtaPhiM(event.jetAK4_pt2,event.jetAK4_eta2,event.jetAK4_phi2,event.jetAK4_mass2)
         mjj=(jet1+jet2).M()
         chi=math.exp(abs(jet1.Rapidity()-jet2.Rapidity()))
         yboost=abs(jet1.Rapidity()+jet2.Rapidity())/2.
         if mjj<1000 or chi>16. or yboost>1.11: continue
         irec=0
	 for massbin in massbins:
	    passedHLT=len(triggers[massbins.index(massbin)])==0
	    for i in trigger_indices[massbin]:
	      if event.HLT_isFired.find(i)!=event.HLT_isFired.end() and event.HLT_isFired[i]:
	        passedHLT=True
		break
	    if passedHLT and yboost<1.11 and mjj>=massbin[0] and mjj<massbin[1]:
               plots[irec].Fill(chi)
	    irec+=1
    #for plot in plots:
    #  if event_count>0:
    #    plot.Scale(xsec/event_count)
    return plots

if __name__ == '__main__':

    wait=False
 
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
 
    samples=[("datacard_shapelimit13TeV_run2_2016","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/data2016"),
            ("datacard_shapelimit13TeV_run2_2017","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/data2017"),
            ("datacard_shapelimit13TeV_run2_2018","","/pnfs/desy.de/cms/tier2/store/user/hinzmann/dijetangular/data2018JECv19"),
            ]

    triggers=[[["HLT_PFHT475","HLT_PFJet260"], #2016
          ["HLT_PFHT475","HLT_PFJet260"],
          ["HLT_PFHT600","HLT_PFHT475","HLT_PFJet320"],
          [], #"HLT_PFHT900","HLT_PFJet450"
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          [],
         ],
	  [["HLT_PFHT510","HLT_PFJet260"], #2017
          ["HLT_PFHT590","HLT_PFHT510","HLT_PFJet260"],
          ["HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet320"],
          ["HLT_PFHT890","HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet450"],
          [], #"HLT_PFHT1050","HLT_PFJet500"
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          ],
	  [["HLT_PFHT510","HLT_PFJet260"], #2017
          ["HLT_PFHT590","HLT_PFHT510","HLT_PFJet260"],
          ["HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet320"],
          ["HLT_PFHT890","HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet450"],
          [], #"HLT_PFHT1050","HLT_PFJet500"
          [],
          [],
          [],
          [],
          [],
          [],
          [],
          ],]

    if len(sys.argv)>1:
      samplenum=int(sys.argv[1])
      samples=[samples[samplenum]]
      triggers=[triggers[samplenum]]

    if len(sys.argv)>2:
      folders=os.listdir(samples[0][2])
      print len(folders)
      samples=[(samples[0][0],"-"+sys.argv[2],"dcap://dcache-cms-dcap.desy.de/"+samples[0][2]+"/"+folders[int(sys.argv[2])])]
 
    chi_binnings=[]
    for mass_bin in chi_bins:
        chi_binnings+=[array.array('d')]
        for chi_bin in mass_bin:
            chi_binnings[-1].append(chi_bin)
        
    print samples

    for prefix,postfix,files in samples:
      plots=[createPlots(files,prefix,triggers[samples.index((prefix,postfix,files))],massbins)]

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
      canvas.SaveAs(prefix+postfix + '_chi.eps')
      if wait:
  	  os.system("ghostview "+prefix + '_chi.eps')
