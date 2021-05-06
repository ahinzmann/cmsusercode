from ROOT import *
import ROOT
import array, math
import os
from math import *

def rebin2(h1,nbins,binning):
    h1=h1.Rebin(nbins,h1.GetName()+"_rebin",binning)
    h1.Scale(1./h1.Integral())
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
    return h1

if __name__=="__main__":

    print "start ROOT"
    #gROOT.Reset()
    gROOT.SetStyle("Plain")
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)
    gStyle.SetTitleOffset(1.2,"Y")
    gStyle.SetPadLeftMargin(0.16)
    gStyle.SetPadBottomMargin(0.16)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetMarkerSize(2.5)
    gStyle.SetHistLineWidth(1)
    gStyle.SetStatFontSize(0.020)
    gStyle.SetTitleSize(0.06, "XYZ")
    gStyle.SetLabelSize(0.05, "XYZ")
    gStyle.SetLegendBorderSize(0)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetEndErrorSize(5)

    f0=TFile.Open("datacard_shapelimit13TeV_GEN-QCD-run2_chi.root")
    f1=TFile.Open("datacard_shapelimit13TeV_run2_2016_chi.root")
    f2=TFile.Open("datacard_shapelimit13TeV_run2_2017_chi.root")
    f3=TFile.Open("datacard_shapelimit13TeV_run2_2018_chi.root")

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

    for i in range(len(massbins)):
      h0=f0.Get("QCD_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
      h0a=f0.Get("QCD_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1_jesUp")
      h0b=f0.Get("QCD_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1_jesDown")
      h1=f1.Get("datacard_shapelimit13TeV_run2_2016#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
      h2=f2.Get("datacard_shapelimit13TeV_run2_2017#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
      h3=f3.Get("datacard_shapelimit13TeV_run2_2018#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1")
      h0=rebin2(h0,len(chi_binnings[i])-1,chi_binnings[i])
      h0a=rebin2(h0a,len(chi_binnings[i])-1,chi_binnings[i])
      h0b=rebin2(h0b,len(chi_binnings[i])-1,chi_binnings[i])
      h1=rebin2(h1,len(chi_binnings[i])-1,chi_binnings[i])
      h2=rebin2(h2,len(chi_binnings[i])-1,chi_binnings[i])
      h3=rebin2(h3,len(chi_binnings[i])-1,chi_binnings[i])
      h0.Scale(1./h0.Integral())
      h0a.Scale(1./h0a.Integral())
      h0b.Scale(1./h0b.Integral())
      h1.Scale(1./h1.Integral())
      h2.Scale(1./h2.Integral())
      h3.Scale(1./h3.Integral())
      
      c = TCanvas("chi", "chi", 0, 0, 300, 300)
      h0.SetLineColor(2)
      h0.GetYaxis().SetTitle("1/N dN/d#chi")
      h0.GetYaxis().SetRangeUser(min(h0.GetMinimum()/2,h2.GetMinimum()/2),max(h0.GetMaximum()*1.5,h2.GetMaximum()*1.5))
      h0.SetLineStyle(2)
      h0.Draw("hist")
      h0a.SetLineColor(6)
      h0a.SetLineStyle(3)
      h0a.Draw("histsame")
      h0b.SetLineColor(6)
      h0b.SetLineStyle(3)
      h0b.Draw("histsame")
      h1.SetLineColor(4)
      h1.Draw("hsame")
      h2.SetLineColor(1)
      h2.Draw("hsame")
      h3.SetLineColor(6)
      h3.Draw("hsame")
      l=TLegend(0.2,0.6,0.80,0.90,str(massbins[i][0])+"<m_{jj}<"+str(massbins[i][1])+" GeV")
      l.SetTextSize(0.05)
      l.AddEntry(h0,"NNLO QCD+EW+detector","l")
      l.AddEntry(h0a," + JES up/down","l")
      l.AddEntry(h1,"2016 data 36.3/fb","l")
      l.AddEntry(h2,"2017 data 41.5/fb","l")
      l.AddEntry(h3,"2018 data 59.7/fb","l")
      l.SetFillStyle(0)
      l.Draw("same")
      c.SaveAs("chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_2016_2017_2018.pdf")
      c.SaveAs("chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_2016_2017_2018.png")
