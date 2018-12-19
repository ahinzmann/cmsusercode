import os, sys
import array
from ROOT import * 
from math import *

#gROOT.Reset()
gROOT.SetStyle("Plain")
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetTitleOffset(1.2,"Y")
gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadBottomMargin(0.15)
gStyle.SetPadTopMargin(0.08)
gStyle.SetPadRightMargin(0.05)
gStyle.SetMarkerSize(0.5)
gStyle.SetHistLineWidth(1)
gStyle.SetStatFontSize(0.020)
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetNdivisions(510, "XYZ")
gStyle.SetLegendBorderSize(0)

if __name__ == '__main__':

   var="chi"
   label="#chi"

   massbins=[(1900,2400),
              (2400,3000),
              (3000,3600),
              (3600,4200),
              (4200,4800),
              (4800,5400),
              (5400,6000),
              #(6000,6600),
	      #(2400,13000),
	      #(3000,13000),
	      #(3600,13000),
	      #(4200,13000),
	      #(4800,13000),
	      #(5400,13000),
	      #(6000,13000),
	      #(6600,13000),
	      (6000,7000),
              (7000,13000),
              ]

   colors=[1,2,3,4,6,7,8,9,10,11,12,13]
   styles=[1,2,3,4,5,6,7,8,9,11,12,13]
   
   sources = ["AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation",
"SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta",
"RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF",
"RelativeBal", "RelativeFSR", "RelativeSample","RelativeStatFSR",
"RelativeStatEC", "RelativeStatHF", "RelativeJEREC1", "RelativeJEREC2",
"RelativeJERHF", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", 
"PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF",
"Total",
  ]

   print "sources",len(sources)
   #sourcesets=[[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15],[16,17,18,19],[20,21,22,23],[24,25,26],[27]]
   #sourcesets=[[i] for i in range(len(sources))]
   sourcesets=[[1,2,5,6,7],[8,12,13,14,16,18],[21,22,23,24],[27]]
   #sourcesets=[[21]]

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
   chi_binnings=[]
   for mass_bin in chi_bins:
        chi_binnings+=[array.array('d')]
        for chi_bin in mass_bin:
            chi_binnings[-1].append(chi_bin)
   
   for prefix in ["QCD","QCDCIplusLL10000"]:
    for year in [2016,2017]:
   #for prefix,year in [("QCD",2016)]:

     sum_in_quadrature_up=[]
     sum_in_quadrature_down=[]
     for mass in range(len(massbins)):
       sum_in_quadrature_up+=[len(chi_binnings[mass])*[0]]
       sum_in_quadrature_down+=[len(chi_binnings[mass])*[0]]

     for sourceset in sourcesets:
      canvas = TCanvas("jes","jes",0,0,600,600)
      canvas.Divide(3,3)
      log=False
      legends=[]
      hists=[]
      files=[]
      for mass in range(min(9,len(massbins))):
        b="1"
	if mass==3: b="2"
        if mass==2: b="3"
        if mass==1: b="4"
        if mass==0: b="6"
	print "datacard_shapelimit13TeV_"+prefix+"_JES_"+str(year)+"_"+b+"_chi.root"
        f_refmc=TFile.Open("datacard_shapelimit13TeV_"+prefix+"_JES_"+str(year)+"_"+b+"_chi.root")
	files+=[f_refmc]
        f_mc=f_refmc
        canvas.cd(mass+1)
        canvas.GetPad(mass+1).SetLogy(log)
    
        hist=f_refmc.Get(prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_rebin1_backup")
	hist=hist.Clone(hist.GetName()+"main"+str(sourcesets.index(sourceset)))
	hists+=[hist]
        hist=hist.Rebin(len(chi_binnings[mass])-1,hist.GetName()+"_rebin1",chi_binnings[mass])
        #if hist.Integral()>0:
        #  hist.Scale(1./hist.Integral())
        #print hist.GetBinContent(2)
        hist.SetLineWidth(2)
      	hist.SetLineColor(1)
	#miny=0
	#if hist.Integral()>0:
        #    miny=log*0.1
	histref=hist.Clone(hist.GetName()+"ref"+str(sourcesets.index(sourceset)))
	hists+=[histref]
        for b in range(hist.GetNbinsX()):
	    hist.SetBinError(b+1,0)
	    hist.SetBinContent(b+1,1)
	hist.GetXaxis().SetTitle(label)
	hist.GetYaxis().SetTitle("N")
	hist.GetYaxis().SetRangeUser(0.95,1.1)
        hist.GetXaxis().SetTitleOffset(1.1)
        hist.GetYaxis().SetTitleOffset(1.1)
        hist.GetXaxis().SetLabelSize(0.05)
        hist.GetYaxis().SetLabelSize(0.05)
        hist.GetXaxis().SetTitleSize(0.06)
        hist.GetYaxis().SetTitleSize(0.06)
	hist.SetTitle("")
	hist.SetStats(False)
        hist.Draw("le")
        legend=TLegend(0.2,0.55,0.95,0.90,(str(massbins[mass][0])+"<m_{jj}<"+str(massbins[mass][1])+" GeV"))
	legends+=[legend]
        legend.AddEntry(hist,"central","l")
	
	for i in range(len(sourceset)):
	    #if sources[sourceset[i]]=="Total": updown="Down"
	    #else:
	    updown="Up"
	    print prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+str(sources[sourceset[i]])+updown+"_rebin1"
            hist2=f_mc.Get(prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+str(sources[sourceset[i]])+updown+"_rebin1")
            hist2=hist2.Rebin(len(chi_binnings[mass])-1,hist2.GetName()+"_rebin1",chi_binnings[mass])
	    if hist2.Integral()>0:
                hist2.Scale(histref.Integral()/hist2.Integral())
            #print hist2.GetBinContent(2)
            hist2.SetLineWidth(1)
            hist2.SetLineColor(colors[i])
            hist2.SetLineStyle(2)
            hist2.Divide(hist2,histref)
            #print hist2.GetBinContent(2)
            #if sourceset[i]==23:
            #  for chi_bin in range(len(chi_binnings[mass])):
	    #    hist2.SetBinContent(chi_bin+1,(hist2.GetBinContent(chi_bin+1)-1.)*10+1.)
            #if sourceset[i]==24:
            #  for chi_bin in range(len(chi_binnings[mass])):
	    #    hist2.SetBinContent(chi_bin+1,(hist2.GetBinContent(chi_bin+1)-1.)*100+1.)
            #miny=0
            #if hist2.Integral()>0:
            #    miny=log*0.1
	    hist2.SetTitle("")
            hist2.SetStats(False)
	    hists+=[hist2]
            hist2.Draw("histsame")
	    if sources[sourceset[i]]=="Total":
              legend.AddEntry(hist2,str(year)+" JEC Single overall variation","l")
            else:
	      legend.AddEntry(hist2,str(year)+" JEC "+sources[sourceset[i]],"l")
	    if sourceset!=[27]:
              for chi_bin in range(len(chi_binnings[mass])):
	       if (hist2.GetBinContent(chi_bin+1)-1.0)*(hist2.GetBinCenter(chi_bin+1)-8.5)>0:
	        sum_in_quadrature_up[mass][chi_bin]=sqrt(pow(sum_in_quadrature_up[mass][chi_bin],2)+pow(hist2.GetBinContent(chi_bin+1)-1.0,2))
	       else:
	        sum_in_quadrature_down[mass][chi_bin]=sqrt(pow(sum_in_quadrature_down[mass][chi_bin],2)+pow(hist2.GetBinContent(chi_bin+1)-1.0,2))
            else:
	      hist2b=hist2.Clone(hist2.GetName()+"SumInQuadrature")
              for chi_bin in range(len(chi_binnings[mass])):
	       if hist2b.GetBinCenter(chi_bin+1)-8.5>0:
	        hist2b.SetBinContent(chi_bin+1,1.0+sum_in_quadrature_up[mass][chi_bin])
	       else:
	        hist2b.SetBinContent(chi_bin+1,1.0-sum_in_quadrature_up[mass][chi_bin])
              hist2b.SetLineColor(colors[i+1])
    	      hists+=[hist2b]
	      fit=TF1(hist2b.GetName()+"smooth","pol3",1,16)
	      hist2b.Fit(fit,"NQ")
              for chi_bin in range(len(chi_binnings[mass])):
	        hist2b.SetBinContent(chi_bin+1,fit.Eval(hist2b.GetBinCenter(chi_bin+1)))
              hist2b.Draw("histsame")
              legend.AddEntry(hist2b,str(year)+" JEC Total sum in quadrature","l")

	#for i in range(len(sourceset)):
	    #if sources[sourceset[i]]=="Total": updown="Down"
	    #else:
	    updown="Down"
	    print prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+str(sources[sourceset[i]])+updown+"_rebin1"
            hist3=f_mc.Get(prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+str(sources[sourceset[i]])+updown+"_rebin1")
            hist3=hist3.Rebin(len(chi_binnings[mass])-1,hist3.GetName()+"_rebin1",chi_binnings[mass])
            if hist3.Integral()>0:
                hist3.Scale(histref.Integral()/hist3.Integral())
            #print hist3.GetBinContent(2)
            hist3.SetLineWidth(1)
            hist3.SetLineColor(colors[i])
            hist3.SetLineStyle(3)
	    hist3.Divide(hist3,histref)
            #print hist3.GetBinContent(2)
            #if sourceset[i]==23:
            #  for chi_bin in range(len(chi_binnings[mass])):
	    #    hist3.SetBinContent(chi_bin+1,(hist3.GetBinContent(chi_bin+1)-1.)*10+1.)
            #if sourceset[i]==24:
            #  for chi_bin in range(len(chi_binnings[mass])):
	    #    hist3.SetBinContent(chi_bin+1,(hist3.GetBinContent(chi_bin+1)-1.)*100+1.)
            #miny=0
            #if hist3.Integral()>0:
            #    miny=log*0.1
	    hist3.SetTitle("")
            hist3.SetStats(False)
	    hists+=[hist3]
            hist3.Draw("histsame")
	    if sourceset!=[27]:
              for chi_bin in range(len(chi_binnings[mass])):
	       if (hist3.GetBinContent(chi_bin+1)-1.0)*(hist3.GetBinCenter(chi_bin+1)-8.5)>0:
	        sum_in_quadrature_up[mass][chi_bin]=sqrt(pow(sum_in_quadrature_up[mass][chi_bin],2)+pow(hist3.GetBinContent(chi_bin+1)-1.0,2))
               else:
	        sum_in_quadrature_down[mass][chi_bin]=sqrt(pow(sum_in_quadrature_down[mass][chi_bin],2)+pow(hist3.GetBinContent(chi_bin+1)-1.0,2))
            else:
	      hist3b=hist3.Clone(hist3.GetName()+"SumInQuadrature")
              for chi_bin in range(len(chi_binnings[mass])):
	       if hist3b.GetBinCenter(chi_bin+1)-8.5>0:
	        hist3b.SetBinContent(chi_bin+1,1.0-sum_in_quadrature_down[mass][chi_bin])
	       else:
	        hist3b.SetBinContent(chi_bin+1,1.0+sum_in_quadrature_down[mass][chi_bin])
              hist3b.SetLineColor(colors[i+1])
    	      hists+=[hist3b]
	      fit=TF1(hist3b.GetName()+"smooth","pol3",1,16)
	      hist3b.Fit(fit,"NQ")
              for chi_bin in range(len(chi_binnings[mass])):
	        hist3b.SetBinContent(chi_bin+1,fit.Eval(hist3b.GetBinCenter(chi_bin+1)))
              hist3b.Draw("histsame")

        legend.SetTextSize(0.04)
        legend.SetFillStyle(0)
        legend.Draw("same")

      canvas.SaveAs("chi_systematic_plots"+var+"_"+prefix+"JES_"+str(sourcesets.index(sourceset))+"_13TeV_"+str(year)+".root")
      canvas.SaveAs("chi_systematic_plots"+var+"_"+prefix+"JES_"+str(sourcesets.index(sourceset))+"_13TeV_"+str(year)+".pdf")
