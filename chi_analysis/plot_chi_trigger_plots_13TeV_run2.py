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

def effErf(x, p):
  return ((TMath.Erf((x[0] - p[0]) / p[1]) + 1) / 2. *  p[2] + (1. - p[2]))*p[3]

if __name__ == '__main__':

  colors=[1,2,4,6,7,8,9,11,40,41,42,43,44,45,46,47,48,49]
  styles=[1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,]

  prefix="datacard_shapelimit13TeV_run2_"
  postfix="_run2"

  chi_bins=[(1,2,3,4,5,6,7,8,9,10,12,14,16)]

  samples=["2016","2017","2018","2016_SingleMuon","2017_SingleMuon","2018_SingleMuon"]
  
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
          ["HLT_PFHT890","HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet450"],
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
          ["HLT_PFHT890","HLT_PFHT780","HLT_PFHT680","HLT_PFHT590","HLT_PFHT510","HLT_PFJet450"],
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
          ]

  for sample in samples:
   f=TFile.Open(prefix+sample+"_chi.root")
   for t in range(13,len(triggers[samples.index(sample)])):
    print sample,t
    name="or".join(triggers[samples.index(sample)][t])
    print name
    canvas = TCanvas("","",0,0,200,200)
    legend=TLegend(0.65,0.2,0.9,0.8,"")
    hists=[]
    i=0
    for chi_bin in [""]+["-chi-"+str(c) for c in chi_bins[0][:-1]]:
      bins=[1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2800,3000]
      binning=array.array('d')
      for bin in bins:
          binning.append(bin)
      hist1=f.Get(prefix+sample+"mass-reftrig"+chi_bin)
      hist1=hist1.Rebin(len(binning)-1,hist1.GetName()+"_rebin",binning)
      hist2=f.Get(prefix+sample+"mass-trig"+name+chi_bin)
      hist2=hist2.Rebin(len(binning)-1,hist2.GetName()+"_rebin",binning)
      if chi_bin=="":
       prescalefactor=1e10
       factor=1e10
       for b in range(hist1.GetNbinsX()):
        if hist2.GetBinContent(b+1)>0:
          factor=hist1.GetBinContent(b+1)/hist2.GetBinContent(b+1)
	if factor<prescalefactor:
	  prescalefactor=factor
       print prescalefactor
      hist2.Scale(prescalefactor) # Scale to account for pre-scales
      hist=TGraphAsymmErrors(hist1)
      if prescalefactor==1:
        mode="cl=0.683 b(1,1) mode"
      else:
        mode="pois"
      hist.Divide(hist2,hist1,mode)
      hist.SetLineWidth(2)
      hist.SetTitle("")
      hist.GetXaxis().SetTitle("dijet mass [GeV]")
      hist.GetYaxis().SetTitle("trigger efficiency")
      if prescalefactor==1:
        hist.GetYaxis().SetRangeUser(0.8,1)
      else:
        hist.GetYaxis().SetRangeUser(0,1.5)
      hist.GetXaxis().SetRangeUser(1200,3000)
      hist.SetLineColor(colors[i])
      hist.SetLineStyle(styles[i])
      startfit=0
      for b in reversed(range(hist.GetN())):
        x = Double(0.)
        y = Double(0.)
        hist.GetPoint(b, x, y)
	if hist.Eval(x)>0.85:
	   startfit=x
      endfit=0
      for b in range(hist.GetN()):
        x = Double(0.)
        y = Double(0.)
        hist.GetPoint(b, x, y)
	if hist.Eval(x)>0.85:
	   endfit=x
      if startfit==0: startfit=1200
      if endfit<=startfit: endfit=3000
      fit=TF1("erf"+name+chi_bin,effErf,startfit,endfit,4)
      fit.SetParameter(0, startfit+100.)
      fit.SetParameter(1, 100.)
      fit.SetParameter(2, 1.)
      if prescalefactor==1:
        fit.FixParameter(3, 1.)
      else:
        fit.SetParameter(3, 1.)
      hist.Fit(fit,"q0")
      fit.SetLineColor(colors[i])
      fit.SetLineStyle(styles[i])
      print chi_bin, fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(2), fit.GetParameter(3)
      s999=0
      for s in reversed(range(int(startfit),int(endfit))):
        if fit.Eval(s)>0.999*fit.GetParameter(3):
	  s999=s
      if chi_bin=="":
        legend.AddEntry(hist,"1<#chi<16 ("+str(s999)+")","le")
        hist.Draw("apez")
      else:
        legend.AddEntry(hist,chi_bin.split("-")[-1]+"<#chi<"+str(chi_bins[0][chi_bins[0].index(int(chi_bin.split("-")[-1]))+1])+" ("+str(s999)+")","le")
        hist.Draw("pezsame")
      fit.Draw("lsame")
      hists+=[hist]
      hists+=[fit]
      i+=1

    legend.SetTextSize(0.04)
    legend.SetFillStyle(0)
    legend.Draw("same")

    canvas.SaveAs("chi_trigger_plots_"+sample+name+postfix+".root")
    canvas.SaveAs("chi_trigger_plots_"+sample+name+postfix+".pdf")
