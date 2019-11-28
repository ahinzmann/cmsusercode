import os, sys
import array
from ROOT import * 

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
gStyle.SetNdivisions(506, "XYZ")
gStyle.SetLegendBorderSize(0)

if __name__ == '__main__':

  colors=[1,2,4,6,7,8,9,11,40,41,42,43,44,45,46,47,48,49]
  styles=[1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,]

  prefix="datacard_shapelimit13TeV_run2_"
  postfix="_run2"

  chi_bins=[(1,2,3,4,5,6,7,8,9,10,12,14,16)]

  samples=["data_2016","data_2017","data_2018"]
  for sample in samples:
    f=TFile.Open(prefix+sample+"_chi.root")

    canvas = TCanvas("","",0,0,200,200)
    legend=TLegend(0.75,0.3,0.9,0.85,"")
    hists=[]
    i=0
    for chi_bin in [""]+["-chi-"+str(c) for c in chi_bins[0][:-1]]:
      #bins=[500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,2000,2500,5000]
      #binning=array.array('d')
      #for bin in bins:
      #    binning.append(bin)
      hist1=f.Get(prefix+sample.replace("data_","")+"mass-reftrig"+chi_bin)
      #hist1=hist1.Rebin(len(binning)-1,hist1.GetName()+"_rebin",binning)
      hist2=f.Get(prefix+sample.replace("data_","")+"mass-trig"+chi_bin)
      #hist2=hist2.Rebin(len(binning)-1,hist2.GetName()+"_rebin",binning)
      hist=TGraphAsymmErrors(hist1)
      hist.Divide(hist2,hist1,"cl=0.683 b(1,1) mode")
      hist.SetLineWidth(2)
      hist.SetTitle("")
      hist.GetXaxis().SetTitle("dijet mass [GeV]")
      hist.GetYaxis().SetTitle("trigger efficiency")
      hist.GetYaxis().SetRangeUser(0.8,1)
      hist.GetXaxis().SetRangeUser(1200,3000)
      hist.SetLineColor(colors[i])
      hist.SetLineStyle(styles[i])
      if chi_bin=="":
        legend.AddEntry(hist,"1<#chi<16","le")
        hist.Draw("ale")
      else:
        legend.AddEntry(hist,chi_bin.split("-")[-1]+"<#chi<"+str(chi_bins[0][chi_bins[0].index(int(chi_bin.split("-")[-1]))+1]),"le")
        hist.Draw("lesame")
      hists+=[hist]
      i+=1

    legend.SetTextSize(0.04)
    legend.SetFillStyle(0)
    legend.Draw("same")

    canvas.SaveAs("chi_trigger_plots_"+sample+postfix+".root")
    canvas.SaveAs("chi_trigger_plots_"+sample+postfix+".pdf")
