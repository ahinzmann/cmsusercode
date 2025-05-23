import os, sys
from ROOT import * 
import array
from math import *

#gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
gROOT.Reset()
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
gStyle.SetNdivisions(506, "XYZ")
gStyle.SetLegendBorderSize(0)


def createPlots(sample,prefix,massbins,factor):
    if sample.endswith(".txt"):
        files=[]
        filelist=open(sample)
	for line in filelist.readlines():
	    if ".root" in line:
	        files+=[line.strip()]
    else:
        files=[sample]

    plots=[]
    genplots=[]
    variables=["#chi","p_{T1}","p_{T2}","y_{1}","y_{2}","y_{1}","y_{2}","y_{1}","y_{2}","y_{1}","y_{2}","y_{boost}","missing E_{T} / #sum E_{T}","#sum #vec p_{T} / #sum |p_{T}|","#Delta #phi","#Delta p_{T} / #sum p_T"]
    for massbin in massbins:
      genplots += [TH1F(prefix+'#genchi'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';#genchi;N',15,1,16)]
      plots += [TH1F(prefix+'response'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';response;N',50,0,2)]
      #plots += [TH1F(prefix+'y_{boost}'+str(massbin).strip("()").replace(',',"_").replace(' ',""),';y_{boost};N',20,0,2)]
    plots += [TH1F(prefix+'pt',';jet pt;N events',1300,0,13000)]
    genplots += [TH1F(prefix+'genmass',';dijet mass;N events',1300,0,13000)]
    
    for plot in plots:
        plot.Sumw2()
    for plot in genplots:
        plot.Sumw2()
    print len(plots),"plots to fill"

    event_count=0
    for f in files[:]:
     print f
     try:
       fil=TFile.Open(f)
       events=fil.Get("ntuplizer/tree")
       nevents=events.GetEntries()
     except:
       print "File not found"
       continue
     print sample,nevents
     for event in events:
       if event_count>200:break
       #if sample!="data_obs" and event_count>3000000:break
       
       event_count+=1
       if event_count%10000==1: print "event",event_count

       if not event.passFilter_HBHE or not event.passFilter_CSCHalo or not event.passFilter_GoodVtx or not event.passFilter_EEBadSc: continue
       if len(event.jetAK4_pt)<2 or event.jetAK4_pt[0]<100 or event.jetAK4_pt[1]<100 or abs(event.jetAK4_eta[0])>3 or abs(event.jetAK4_eta[1])>3: continue
       if not bool(event.jetAK4_IDTight[0]) or not bool(event.jetAK4_IDTight[1]): continue

       if sample!="data_obs" and len(event.genJetAK4_pt)>=2:
         genJet1=TLorentzVector()
         genJet2=TLorentzVector()
         genJet1.SetPtEtaPhiM(event.genJetAK4_pt[0],event.genJetAK4_eta[0],event.genJetAK4_phi[0],event.genJetAK4_mass[0])
         genJet2.SetPtEtaPhiM(event.genJetAK4_pt[1],event.genJetAK4_eta[1],event.genJetAK4_phi[1],event.genJetAK4_mass[1])
         genmjj=(genJet1+genJet2).M()
         genchi=exp(abs(genJet1.Rapidity()-genJet2.Rapidity()))
         genyboost=abs(genJet1.Rapidity()+genJet2.Rapidity())/2.
         igen=0
         for massbin in massbins:
           if genyboost<1.11 and genmjj>=massbin[0] and genmjj<massbin[1]:
	     genplots[igen].Fill(genchi)
  	   igen+=1
         if genyboost<1.11 and genchi<16:
           genplots[igen].Fill(genmjj)
         igen+=1

       scales=[1]
       irec=0
       for scale in scales:
         jet1=TLorentzVector()
         jet2=TLorentzVector()
         jet1.SetPtEtaPhiM(event.jetAK4_pt[0],event.jetAK4_eta[0],event.jetAK4_phi[0],event.jetAK4_mass[0])
         jet2.SetPtEtaPhiM(event.jetAK4_pt[1],event.jetAK4_eta[1],event.jetAK4_phi[1],event.jetAK4_mass[1])
         mjj=(jet1+jet2).M()
         chi=exp(abs(jet1.Rapidity()-jet2.Rapidity()))
         yboost=abs(jet1.Rapidity()+jet2.Rapidity())/2.
         #if mjj<1000 or chi>16. or yboost>1.11: continue
         for massbin in massbins:
           if yboost<1.11 and chi<16 and jet1.Pt()>600 and jet1.Pt()<800 and genJet1.DeltaR(jet1)<0.25:
             plots[irec].Fill(jet1.Pt()/genJet1.Pt())
	     #print jet1.Pt()/genJet1.Pt(), event.jetAK4_jec[0]
           irec+=1
         if scale==1:
           if yboost<1.11 and chi<16:
             plots[irec].Fill(jet1.Pt())
           irec+=1

    print "analyzed",event_count,"events"
    if event_count>0:
      for plot in plots:
        plot.Scale(factor/event_count)
      for plot in genplots:
        plot.Scale(factor/event_count)
    return (plots,genplots,nevents)

if __name__ == '__main__':

    wait=False
 
    prefix="datacard_shapelimit13TeV_25nsMC_response"
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
    massbins=[(1900,2400),
              #(2400,3000),
              #(3000,3600),
              #(3600,4200),
              #(4200,4800),
              #(4800,5400),
              #(5400,6000),
	      #(2400,13000),
	      #(3000,13000),
	      #(3600,13000),
	      #(4200,13000),
	      #(4800,13000),
	      #(5400,13000),
	      #(6000,13000),
              ]
 
    samples=[#("data_obs",[("JetHT_25ns_data5.txt",1.)]),
             #("QCD",[("QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6.txt",1.)])
             ("QCD",[#("QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8.txt",0.000165),
               #("QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8.txt",0.006830),
               #("QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8.txt",0.114943),
               #("QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8.txt",0.842650),
               #("QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8.txt",9.4183),
               ("QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8.txt",32.293),
               ("QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8.txt",186.9),
               ("QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8.txt",648.2),],)
               #("QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8.txt",7823.0)],)
            ]
 
    chi_binnings=[]
    for mass_bin in chi_bins:
        chi_binnings+=[array.array('d')]
        for chi_bin in mass_bin:
            chi_binnings[-1].append(chi_bin)
        
    print prefix, samples

    plots=[]
    genplots=[]
    nevents={}
    for name,files in samples:
      plots+=[[]]
      genplots+=[[]]
      i=0
      for filename,xsec in files:
        plot,genplot,nevent=createPlots(filename,name,massbins,xsec)
	if i>0:
	  for j in range(len(plots[-1])):
	     plots[-1][j].Add(plot[j])
	  for j in range(len(genplots[-1])):
	     genplots[-1][j].Add(genplot[j])
	else:
          plots[-1]+=plot
          genplots[-1]+=genplot
	nevents[name]=nevent
	i+=1

    out=TFile(prefix + '_chi.root','RECREATE')
    for i in range(len(samples)):
      for j in range(len(massbins)):
        print "number of "+samples[i][0]+" jet events expected in",str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV:",plots[i][j].Integral()
        print "number of "+samples[i][0]+" genjet events expected in",str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV:",genplots[i][j].Integral()
        #print "mass: number of "+samples[i][0]+" events expected in",str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV:",plots[i][len(massbins)].Integral(plots[i][len(massbins)].FindBin(massbins[j][0]),plots[i][len(massbins)].FindBin(massbins[j][1]))
	if "chi" in plots[i][j].GetName():
          plots[i][j]=plots[i][j].Rebin(len(chi_binnings[j])-1,plots[i][j].GetName()+"_rebin1",chi_binnings[j])
	  clone=plots[i][j]
	  clone.Write()
	if "chi" in genplots[i][j].GetName():
          genplots[i][j]=genplots[i][j].Rebin(len(chi_binnings[j])-1,genplots[i][j].GetName()+"_rebin1",chi_binnings[j])
	  clone=genplots[i][j]
	  clone.Write()

    for i in range(len(samples)):
      for j in range(len(massbins)):
        if plots[i][j].Integral()>0:
          plots[i][j].Scale(1./plots[i][j].Integral())
        for b in range(plots[i][j].GetXaxis().GetNbins()):
          plots[i][j].SetBinContent(b+1,plots[i][j].GetBinContent(b+1)/plots[i][j].GetBinWidth(b+1))
          plots[i][j].SetBinError(b+1,plots[i][j].GetBinError(b+1)/plots[i][j].GetBinWidth(b+1))
#        plots[i][j].GetYaxis().SetRangeUser(0,0.5)
        if genplots[i][j].Integral()>0:
          genplots[i][j].Scale(1./genplots[i][j].Integral())
        for b in range(genplots[i][j].GetXaxis().GetNbins()):
          genplots[i][j].SetBinContent(b+1,genplots[i][j].GetBinContent(b+1)/genplots[i][j].GetBinWidth(b+1))
          genplots[i][j].SetBinError(b+1,genplots[i][j].GetBinError(b+1)/genplots[i][j].GetBinWidth(b+1))
#        genplots[i][j].GetYaxis().SetRangeUser(0,0.5)

    canvas = TCanvas("","",0,0,400,200)
    canvas.SetLogy()
    if len(massbins)>1:
      canvas.Divide(2,1)
    if len(massbins)>4:
      canvas = TCanvas("","",0,0,600,400)
      canvas.Divide(3,2)
    elif len(massbins)>2:
      canvas = TCanvas("","",0,0,400,400)
      canvas.Divide(2,2)

    legends=[]
    color=[1,2,4,6,7,8,9]
    for j in range(len(massbins)):
      canvas.cd(j+1)
      plots[0][j].Draw("he")
      legend1=TLegend(0.6,0.6,0.9,0.9,(str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV"))
      legends+=[legend1]
      legend1.AddEntry(plots[0][j],samples[0][0],"l")
      for i in range(1,len(samples)):
        plots[i][j].SetLineColor(color[i])
        plots[i][j].Draw("hesame")
        legend1.AddEntry(plots[i][j],samples[i][0],"l")
      legend1.SetTextSize(0.04)
      legend1.SetFillStyle(0)
      legend1.Draw("same")

    canvas.SaveAs(prefix + '_response.pdf')
    canvas.SaveAs(prefix + '_response.eps')
    if wait:
        os.system("ghostview "+prefix + '_response.eps')

    canvas = TCanvas("","",0,0,200,200)
    canvas.SetLogy()

    legends=[]
    plots[0][len(massbins)].Draw("he")
    plots[0][len(massbins)].GetXaxis().SetRangeUser(0,5000)
    legend1=TLegend(0.6,0.6,0.9,0.9)
    legends+=[legend1]
    legend1.AddEntry(plots[0][len(massbins)],samples[0][0],"l")
    for i in range(1,len(samples)):
      plots[i][len(massbins)].Scale(plots[0][len(massbins)].Integral()/plots[i][len(massbins)].Integral())
      plots[i][len(massbins)].SetLineColor(color[i])
      plots[i][len(massbins)].Draw("hesame")
      legend1.AddEntry(plots[i][len(massbins)],samples[i][0],"l")
    legend1.SetTextSize(0.04)
    legend1.SetFillStyle(0)
    legend1.Draw("same")

    canvas.SaveAs(prefix + '_pt.pdf')
    canvas.SaveAs(prefix + '_pt.eps')
    if wait:
        os.system("ghostview "+prefix + '_pt.eps')

