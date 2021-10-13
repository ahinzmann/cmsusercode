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
   
   lumifactor={}
   lumifactor["2016"]=36.33/137.6
   lumifactor["2017"]=41.53/137.6
   lumifactor["2018"]=59.74/137.6

   samples=[("pythia8_ci_m4300_13000_50000_1_0_0_13TeV_Nov14",3.507e-09),
     	    ("pythia8_ci_m3800_4300_50000_1_0_0_13TeV_Nov14",5.867e-09),
     	    ("pythia8_ci_m3300_3800_50000_1_0_0_13TeV_Nov14",1.863e-08),
     	    ("pythia8_ci_m2800_3300_50000_1_0_0_13TeV_Nov14",6.446e-08),
     	    ("pythia8_ci_m2400_2800_50000_1_0_0_13TeV_Nov14",1.649e-07),
     	    ("pythia8_ci_m1900_2400_50000_1_0_0_13TeV_Nov14",8.836e-07),
     	    ("pythia8_ci_m1500_1900_50000_1_0_0_13TeV_Nov14",3.307e-06),]


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

   colors=[1,2,3,4,6,7,8,9,11,12,13,14]
   styles=[1,2,3,4,5,6,7,8,9,11,12,13]
   
   sourcesets=[[("AbsoluteScale",""),
                ("AbsoluteMPFBias",""),
		("SinglePionHCAL",""),
		("FlavorQCD",""),
		("TimePtEta","2016"),("TimePtEta","2017"),("TimePtEta","2018"),],
		[("RelativePtBB",""),
		("RelativeBal",""),
		("RelativeFSR",""),],
		[("RelativeSample","2016"),("RelativeSample","2017"),("RelativeSample","2018"),
		("RelativeStatEC","2016"),("RelativeStatEC","2017"),("RelativeStatEC","2018"),
		("RelativeJEREC1","2016"),("RelativeJEREC1","2017"),("RelativeJEREC1","2018"),],
		[("PileUpDataMC",""),
		("PileUpPtRef",""),
		("PileUpPtBB",""),
		("PileUpPtEC1",""),],
		[("Sum","")]]
   
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
   chi_binnings=[]
   for mass_bin in chi_bins:
        chi_binnings+=[array.array('d')]
        for chi_bin in mass_bin:
            chi_binnings[-1].append(chi_bin)
   
   for prefix in ["QCD","QCDCIplusLL10000","QCDmadgraph"]:
     sum_in_quadrature_up=[]
     sum_in_quadrature_down=[]
     for mass in range(len(massbins)):
       sum_in_quadrature_up+=[len(chi_binnings[mass])*[0]]
       sum_in_quadrature_down+=[len(chi_binnings[mass])*[0]]

     for sourceset in sourcesets:
      canvas = TCanvas("jes","jes",0,0,800,600)
      canvas.Divide(4,3)
      log=False
      legends=[]
      hists=[]
      files=[]
      for mass in range(len(massbins)):
	#if mass>=6: b="1"
	#if mass==5: b="2"
        #if mass==4: b="3"
        #if mass==3: b="4"
        #if mass==2: b="6"
        #if mass<=1: b="7"
	f_refmc={}
	histyear={}
	scalexsec={}
	for year in ["2016","2017","2018"]:
	 for b in ["1","2","3","4","5","6","7"]:
	  print "datacard_shapelimit13TeV_"+prefix+"_JES_"+str(year)+"_"+b+"_chi.root"
          f_refmc[year+b]=TFile.Open("datacard_shapelimit13TeV_"+prefix+"_JES_"+str(year)+"_"+b+"_chi.root")
	  h=f_refmc[year+b].Get(prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_rebin1_backup")
	  scalexsec[b]=samples[int(b)-1][1]
	  files+=[f_refmc[year+b]]
	  if b=="1":
            histyear[year]=h
	    histyear[year]=histyear[year].Clone(histyear[year].GetName()+year+"main"+str(sourcesets.index(sourceset)))
            histyear[year].Scale(1./scalexsec[b])
          else:
            histyear[year].Add(h,1./scalexsec[b])
         hists+=[histyear[year]]
         histyear[year]=histyear[year].Rebin(len(chi_binnings[mass])-1,histyear[year].GetName()+"_rebin1",chi_binnings[mass])
        hist=histyear["2016"].Clone(histyear["2016"].GetName()+"combine")
	hists+=[hist]
	hist.Scale(lumifactor["2016"])
	hist.Add(histyear["2017"],lumifactor["2017"])
	hist.Add(histyear["2018"],lumifactor["2018"])
        canvas.cd(mass+1)
        canvas.GetPad(mass+1).SetLogy(log)
        hist.SetLineWidth(2)
      	hist.SetLineColor(1)
	histref=hist.Clone(hist.GetName()+"ref"+str(sourcesets.index(sourceset)))
	hists+=[histref]
        for b in range(hist.GetNbinsX()):
	    hist.SetBinError(b+1,0)
	    hist.SetBinContent(b+1,1)
	hist.GetXaxis().SetTitle(label)
	hist.GetYaxis().SetTitle("N")
	hist.GetYaxis().SetRangeUser(0.9,1.15)
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
            if sourceset[i][0]=="Sum":
	      hist2b=hist.Clone(hist.GetName()+"SumInQuadratureUp")
              for chi_bin in range(len(chi_binnings[mass])):
	       if hist2b.GetBinCenter(chi_bin+1)-8.5>0:
	        hist2b.SetBinContent(chi_bin+1,1.0+sum_in_quadrature_up[mass][chi_bin])
	       else:
	        hist2b.SetBinContent(chi_bin+1,1.0-sum_in_quadrature_up[mass][chi_bin])
               hist2b.SetBinError(chi_bin+1,sum_in_quadrature_up[mass][chi_bin]/10.)
              hist2b.SetLineColor(colors[i+1])
    	      hists+=[hist2b]
	      fit=TF1(hist2b.GetName()+"smooth","pol2" if "mad" in prefix else "pol3",1,16)
	      hist2b.Fit(fit,"NQ")
              for chi_bin in range(len(chi_binnings[mass])):
	        hist2b.SetBinContent(chi_bin+1,fit.Eval(hist2b.GetBinCenter(chi_bin+1)))
              hist2b.Draw("histsame")
              legend.AddEntry(hist2b,"JEC Total sum in quadrature","l")

	      hist3b=hist.Clone(hist.GetName()+"SumInQuadratureDown")
              for chi_bin in range(len(chi_binnings[mass])):
	       if hist3b.GetBinCenter(chi_bin+1)-8.5>0:
	        hist3b.SetBinContent(chi_bin+1,1.0-sum_in_quadrature_down[mass][chi_bin])
	       else:
	        hist3b.SetBinContent(chi_bin+1,1.0+sum_in_quadrature_down[mass][chi_bin])
               hist3b.SetBinError(chi_bin+1,sum_in_quadrature_down[mass][chi_bin]/10.)
              hist3b.SetLineColor(colors[i+1])
    	      hists+=[hist3b]
	      fit=TF1(hist3b.GetName()+"smooth","pol2" if "mad" in prefix else "pol3",1,16)
	      hist3b.Fit(fit,"NQ")
              for chi_bin in range(len(chi_binnings[mass])):
	        hist3b.SetBinContent(chi_bin+1,fit.Eval(hist3b.GetBinCenter(chi_bin+1)))
              hist3b.Draw("histsame")

            else:

	     updown="Up"
	     histyear2={}
	     for year in ["2016","2017","2018"]:
	      for b in ["1","2","3","4","5","6","7"]:
	       print prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+str(sourceset[i][0])+updown+"_rebin1"
	       h=f_refmc[year+b].Get(prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+str(sourceset[i][0])+updown+"_rebin1")
	       if b=="1":
                  histyear2[year]=h
	          histyear2[year]=histyear2[year].Clone(histyear2[year].GetName()+year+"up"+str(sourcesets.index(sourceset)))
	          histyear2[year].Scale(1./scalexsec[b])
               else:
                  histyear2[year].Add(h,1./scalexsec[b])
              histyear2[year]=histyear2[year].Rebin(len(chi_binnings[mass])-1,histyear2[year].GetName()+"_rebin1",chi_binnings[mass])
              hists+=[histyear2[year]]
	     if sourceset[i][1]=="":
	      hist2=histyear2["2016"].Clone(histyear2["2016"].GetName()+"combine")
  	      hists+=[hist2]
	      hist2.Scale(lumifactor["2016"])
	      hist2.Add(histyear2["2017"],lumifactor["2017"])
	      hist2.Add(histyear2["2018"],lumifactor["2018"])
	      if hist2.Integral()>0:
                hist2.Scale(histref.Integral()/hist2.Integral())
              hist2.Divide(hist2,histref)
	     else:
              hist2=histyear2[sourceset[i][1]]
	      hist2.SetName(hist2.GetName().replace(sourceset[i][0],sourceset[i][0]+sourceset[i][1]))
	      if hist2.Integral()>0:
                hist2.Scale(histyear[sourceset[i][1]].Integral()/hist2.Integral())
 	      hist2.Divide(hist2,histyear[sourceset[i][1]])
 	      hist2.Add(hist,-1)
	      hist2.Scale(lumifactor[sourceset[i][1]])
 	      hist2.Add(hist,1)
             hist2.SetLineWidth(1)
             hist2.SetLineColor(colors[i])
             hist2.SetLineStyle(2)
	     hist2.SetTitle("")
             hist2.SetStats(False)
	     fit=TF1(hist2.GetName()+"smooth","pol2" if "mad" in prefix else "pol3",1,16)
	     hist2.Fit(fit,"NQ")
             for chi_bin in range(len(chi_binnings[mass])):
	       hist2.SetBinContent(chi_bin+1,fit.Eval(hist2.GetBinCenter(chi_bin+1)))
             hist2.Draw("histsame")
	     legend.AddEntry(hist2,"JEC "+sourceset[i][0]+" "+sourceset[i][1],"l")
             for chi_bin in range(len(chi_binnings[mass])):
	       if (hist2.GetBinContent(chi_bin+1)-1.0)*(hist2.GetBinCenter(chi_bin+1)-8.5)>0:
	        sum_in_quadrature_up[mass][chi_bin]=sqrt(pow(sum_in_quadrature_up[mass][chi_bin],2)+pow(hist2.GetBinContent(chi_bin+1)-1.0,2))
	       else:
	        sum_in_quadrature_down[mass][chi_bin]=sqrt(pow(sum_in_quadrature_down[mass][chi_bin],2)+pow(hist2.GetBinContent(chi_bin+1)-1.0,2))

	     updown="Down"
	     histyear3={}
	     for year in ["2016","2017","2018"]:
	      for b in ["1","2","3","4","5","6","7"]:
	       print prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+str(sourceset[i][0])+updown+"_rebin1"
	       h=f_refmc[year+b].Get(prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+str(sourceset[i][0])+updown+"_rebin1")
	       if b=="1":
                  histyear3[year]=h
	          histyear3[year]=histyear3[year].Clone(histyear3[year].GetName()+year+"down"+str(sourcesets.index(sourceset)))
	          histyear3[year].Scale(1./scalexsec[b])
               else:
                  histyear3[year].Add(h,1./scalexsec[b])
	      histyear3[year]=histyear3[year].Rebin(len(chi_binnings[mass])-1,histyear3[year].GetName()+"_rebin1",chi_binnings[mass])
  	      hists+=[histyear3[year]]
             if sourceset[i][1]=="":
	      hist3=histyear3["2016"].Clone(histyear3["2016"].GetName()+"combine")
  	      hists+=[hist3]
	      hist3.Scale(lumifactor["2016"])
	      hist3.Add(histyear3["2017"],lumifactor["2017"])
	      hist3.Add(histyear3["2018"],lumifactor["2018"])
	      if hist3.Integral()>0:
                hist3.Scale(histref.Integral()/hist3.Integral())
 	      hist3.Divide(hist3,histref)
	     else:
              hist3=histyear3[sourceset[i][1]]
	      hist3.SetName(hist3.GetName().replace(sourceset[i][0],sourceset[i][0]+sourceset[i][1]))
	      if hist3.Integral()>0:
                hist3.Scale(histyear[sourceset[i][1]].Integral()/hist3.Integral())
 	      hist3.Divide(hist3,histyear[sourceset[i][1]])
 	      hist3.Add(hist,-1)
	      hist3.Scale(lumifactor[sourceset[i][1]])
 	      hist3.Add(hist,1)
             hist3.SetLineWidth(1)
             hist3.SetLineColor(colors[i])
             hist3.SetLineStyle(3)
	     hist3.SetTitle("")
             hist3.SetStats(False)
	     fit=TF1(hist3.GetName()+"smooth","pol2" if "mad" in prefix else "pol3",1,16)
	     hist3.Fit(fit,"NQ")
             for chi_bin in range(len(chi_binnings[mass])):
	       hist3.SetBinContent(chi_bin+1,fit.Eval(hist3.GetBinCenter(chi_bin+1)))
	     hist3.Draw("histsame")
             for chi_bin in range(len(chi_binnings[mass])):
	       if (hist3.GetBinContent(chi_bin+1)-1.0)*(hist3.GetBinCenter(chi_bin+1)-8.5)>0:
	        sum_in_quadrature_up[mass][chi_bin]=sqrt(pow(sum_in_quadrature_up[mass][chi_bin],2)+pow(hist3.GetBinContent(chi_bin+1)-1.0,2))
               else:
	        sum_in_quadrature_down[mass][chi_bin]=sqrt(pow(sum_in_quadrature_down[mass][chi_bin],2)+pow(hist3.GetBinContent(chi_bin+1)-1.0,2))

        legend.SetTextSize(0.04)
        legend.SetFillStyle(0)
        legend.Draw("same")

      canvas.SaveAs("chi_systematic_plots"+var+"_"+prefix+"JES_"+str(sourcesets.index(sourceset))+"_13TeV_run2.root")
      canvas.SaveAs("chi_systematic_plots"+var+"_"+prefix+"JES_"+str(sourcesets.index(sourceset))+"_13TeV_run2.pdf")
