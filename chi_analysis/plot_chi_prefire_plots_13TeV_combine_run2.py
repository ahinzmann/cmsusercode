
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

   dire=""

   var="chi"
   label="#chi"

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
   
   for prefix in ["run2"]:
     sum_in_quadrature_up=[]
     sum_in_quadrature_down=[]
     for mass in range(len(massbins)):
       sum_in_quadrature_up+=[len(chi_binnings[mass])*[0]]
       sum_in_quadrature_down+=[len(chi_binnings[mass])*[0]]

     canvas = TCanvas("prefire","prefire",0,0,800,600)
     canvas.Divide(4,3)
     log=False
     legends=[]
     hists=[]
     files=[]
     for mass in range(len(massbins)):
       f_refmc={}
       histyear={}
       scalexsec={}
       for year in ["2016","2017","2018"]:
        filename=(dire+"datacard_shapelimit13TeV_"+prefix+"_"+str(year)+"_L1prefire_chi.root").replace("2018_L1prefire","2018_HEM")
	print filename
     	f_refmc[year]=TFile.Open(filename)
	histname="datacard_shapelimit13TeV_"+prefix+"_"+str(year)+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_rebin1_backup"
	print histname
        h=f_refmc[year].Get(histname)
        files+=[f_refmc[year]]
        histyear[year]=h
        histyear[year]=histyear[year].Clone(histyear[year].GetName()+year+"main")
     	hists+=[histyear[year]]
     	histyear[year]=histyear[year].Rebin(len(chi_binnings[mass])-1,histyear[year].GetName()+"_rebin1",chi_binnings[mass])
       hist=histyear["2016"].Clone(histyear["2016"].GetName()+"combine")
       hists+=[hist]
       hist.Add(histyear["2017"])
       hist.Add(histyear["2018"])
       canvas.cd(mass+1)
       canvas.GetPad(mass+1).SetLogy(log)
       hist.SetLineWidth(2)
       hist.SetLineColor(1)
       histref=hist.Clone(hist.GetName()+"ref")
       hists+=[histref]
       for b in range(hist.GetNbinsX()):
           hist.SetBinError(b+1,0)
           hist.SetBinContent(b+1,1)
       hist.GetXaxis().SetTitle(label)
       hist.GetYaxis().SetTitle("N")
       hist.GetYaxis().SetRangeUser(0.9,1.1)
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

       updown="prefire_Up"
       histyear2={}
       for year in ["2016","2017","2018"]:
        filename=(dire+"datacard_shapelimit13TeV_"+prefix+"_"+str(year)+"_chi.root")
	print filename
     	f_refmc[year]=TFile.Open(filename)
	histname="datacard_shapelimit13TeV_"+prefix+"_"+str(year)+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_rebin1_backup"
	print histname
        h=f_refmc[year].Get(histname)
        files+=[f_refmc[year]]
        histyear2[year]=h
        histyear2[year]=histyear2[year].Clone(histyear2[year].GetName()+year+"prefireUp")
     	hists+=[histyear2[year]]
     	histyear2[year]=histyear2[year].Rebin(len(chi_binnings[mass])-1,histyear2[year].GetName()+"_rebin1",chi_binnings[mass])
       hist2=histyear2["2016"].Clone(histyear2["2016"].GetName()+"combine")
       hists+=[hist2]
       hist2.Add(histyear2["2017"])
       hist2.Add(histyear2["2018"])
       if hist2.Integral()>0:
         hist2.Scale(histref.Integral()/hist2.Integral())
       hist3=hist2.Clone(hist2.GetName()+"prefireDown")
       hist2.Divide(hist2,histref)
       hist2.SetLineWidth(1)
       hist2.SetLineColor(colors[0])
       hist2.SetLineStyle(2)
       hist2.SetTitle("")
       hist2.SetStats(False)
       pol="pol2"
       if mass>6: pol="pol1"
       fit2=TF1(hist2.GetName()+"smooth",pol,1,16)
       fit2.SetLineStyle(2)
       hist2.Fit(fit2,"NQ")
       hist2.Draw("hesame")
       #for chi_bin in range(len(chi_binnings[mass])):
       #  hist2.SetBinContent(chi_bin+1,fit2.Eval(hist2.GetBinCenter(chi_bin+1)))
       fit2.Draw("lsame")
       legend.AddEntry(hist2,"prefire","l")
       for chi_bin in range(len(chi_binnings[mass])):
         if (hist2.GetBinContent(chi_bin+1)-1.0)*(hist2.GetBinCenter(chi_bin+1)-8.5)>0:
          sum_in_quadrature_up[mass][chi_bin]=sqrt(pow(sum_in_quadrature_up[mass][chi_bin],2)+pow(hist2.GetBinContent(chi_bin+1)-1.0,2))
         else:
          sum_in_quadrature_down[mass][chi_bin]=sqrt(pow(sum_in_quadrature_down[mass][chi_bin],2)+pow(hist2.GetBinContent(chi_bin+1)-1.0,2))

       hist3.Divide(histref,hist3)
       hist3.SetLineWidth(1)
       hist3.SetLineColor(colors[0])
       hist3.SetLineStyle(3)
       hist3.SetTitle("")
       hist3.SetStats(False)
       fit3=TF1(hist3.GetName()+"smooth",pol,1,16)
       fit3.SetLineStyle(3)
       hist3.Fit(fit3,"NQ")
       hist3.Draw("hesame")
       #for chi_bin in range(len(chi_binnings[mass])):
       #  hist3.SetBinContent(chi_bin+1,fit3.Eval(hist3.GetBinCenter(chi_bin+1)))
       fit3.Draw("lsame")
       for chi_bin in range(len(chi_binnings[mass])):
	 if (hist3.GetBinContent(chi_bin+1)-1.0)*(hist3.GetBinCenter(chi_bin+1)-8.5)>0:
	  sum_in_quadrature_up[mass][chi_bin]=sqrt(pow(sum_in_quadrature_up[mass][chi_bin],2)+pow(hist3.GetBinContent(chi_bin+1)-1.0,2))
	 else:
          sum_in_quadrature_down[mass][chi_bin]=sqrt(pow(sum_in_quadrature_down[mass][chi_bin],2)+pow(hist3.GetBinContent(chi_bin+1)-1.0,2))

       legend.SetTextSize(0.04)
       legend.SetFillStyle(0)
       legend.Draw("same")

     canvas.SaveAs(dire+"chi_systematic_plots"+var+"_"+prefix+"prefire_13TeV_run2.root")
     canvas.SaveAs(dire+"chi_systematic_plots"+var+"_"+prefix+"prefire_13TeV_run2.pdf")
