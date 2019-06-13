import os, sys
import array
from ROOT import * 
from math import *

gROOT.Reset()
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

   masses=[1900,2400,3000,3600,4200,4800,5400,6000,13000]

   colors=[1,2,3,4,6,7,8,9,10,11,12,13]
   styles=[1,2,3,4,5,6,7,8,9,11,12,13]
   
   sets = [1,2,3,4,5,7,9] # 6 scale variation avoiding combinations of muR and muF with factor 4 difference between muR and muF

   chi_bins=[(1,2,3,4,5,6,7,8,9,10,12,14,16),
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

   for signalMass in [2000,2500,3000,3500,4000,5000,6000]:   
     canvas = TCanvas("scale","scale",0,0,600,600)
     canvas.Divide(3,3)
     legends=[]
     hists=[]
     files=[]
     for mass in range(len(masses)-1):
       print signalMass,mass
       mean={}
       for set in sets:
        name="DMVector_Dijet_LO_Mphi_"+str(signalMass)+"_1_1p0_1p0_Jun26pdf_"+str(set)
        f_refmc=TFile.Open("scaleweights/datacard_shapelimit13TeV_"+name+"_chi2016.root")
	files+=[f_refmc]
        f_mc=f_refmc
        canvas.cd(mass+1)
        legend=TLegend(0.5,0.55,0.95,0.90,(str(masses[mass])+"<m_{jj}<"+str(masses[mass+1])+" GeV").replace("6000<m_{jj}<13000","m_{jj}>6000"))
	legends+=[legend]
    
        hist=f_refmc.Get(name+"#chi"+str(masses[mass])+"_"+str(masses[mass+1])+"_rebin1")
        hist=hist.Rebin(len(chi_binnings[mass])-1,hist.GetName()+"_rebin1",chi_binnings[mass])
	hists+=[hist]
        hist.SetLineWidth(1)
      	hist.SetLineColor(1)
	hist.GetXaxis().SetTitle(label)
	hist.GetYaxis().SetTitle("N")
	#hist.GetYaxis().SetRangeUser(0.9,1.3)
        hist.GetXaxis().SetTitleOffset(1.1)
        hist.GetYaxis().SetTitleOffset(1.1)
        hist.GetXaxis().SetLabelSize(0.05)
        hist.GetYaxis().SetLabelSize(0.05)
        hist.GetXaxis().SetTitleSize(0.06)
        hist.GetYaxis().SetTitleSize(0.06)
	hist.SetTitle("")
	hist.SetStats(False)
	if set==sets[0]:
	  firsthist=hist.Clone(hist.GetName()+"first")
	  files+=[firsthist]
	  for b in range(firsthist.GetNbinsX()):
	     firsthist.SetBinContent(b+1,firsthist.GetBinContent(b+1)/firsthist.GetXaxis().GetBinWidth(b+1))
	     firsthist.SetBinError(b+1,firsthist.GetBinError(b+1)/firsthist.GetXaxis().GetBinWidth(b+1))
	  firsthist.Draw("le")
	#else:
	#  hist.Draw("lesame")
	for b in range(hist.GetNbinsX()):
	   if set==sets[0]:
	     mean[b]=hist.GetBinContent(b+1)
       histmean=hists[-1].Clone(hist.GetName()+"mean")
       for b in range(hist.GetNbinsX()):
	  histmean.SetBinContent(b+1,mean[b]/hist.GetXaxis().GetBinWidth(b+1))
       histmeang=TGraphAsymmErrors(histmean)
       up={}
       down={}
       for set in sets:
         for b in range(hist.GetNbinsX()):
	   if set==sets[0]:
	     up[b]=max(hists[-len(sets)+sets.index(set)].GetBinContent(b+1)-mean[b],0)
             down[b]=min(hists[-len(sets)+sets.index(set)].GetBinContent(b+1)-mean[b],0)
	   else:
	     up[b]=max(hists[-len(sets)+sets.index(set)].GetBinContent(b+1)-mean[b],up[b])
             down[b]=min(hists[-len(sets)+sets.index(set)].GetBinContent(b+1)-mean[b],down[b])
       for b in range(hist.GetNbinsX()):
	  histmeang.SetPointEYlow(b,-down[b]/hist.GetXaxis().GetBinWidth(b+1))
	  histmeang.SetPointEYhigh(b,up[b]/hist.GetXaxis().GetBinWidth(b+1))
	  if b==0:
	    print mean[b],up[b]/max(mean[b],1e-10),-down[b]/max(mean[b],1e-10)
       histmeang.SetLineWidth(2)
       histmeang.SetLineColor(2)
       histmeang.Draw("lesame")
       files+=[histmeang]

       legend.SetTextSize(0.04)
       legend.SetFillStyle(0)
       legend.Draw("same")
       
     canvas.SaveAs("chi_dm_scale_plots"+str(signalMass)+"_13TeV_2016.root")
     canvas.SaveAs("chi_dm_scale_plots"+str(signalMass)+"_13TeV_2016.pdf")
