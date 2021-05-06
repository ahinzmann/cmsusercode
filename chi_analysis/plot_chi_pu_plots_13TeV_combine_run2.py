
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

   dire="Oct29/"

   var="chi"
   label="#chi"
   
   lumifactor={}
   lumifactor["2016"]=36.33/137.60
   lumifactor["2017"]=41.53/137.60
   lumifactor["2018"]=59.74/137.60

   mc={}
   mc["2016"]=[("2016_QCDmadgraph-HT200to300",1712000./56709875),
       ("2016_QCDmadgraph-HT300to500",347700./53096517),
       ("2016_QCDmadgraph-HT500to700",32100./52906552),
       ("2016_QCDmadgraph-HT700to1000",6831./36741540),
       ("2016_QCDmadgraph-HT1000to1500",1207./15210939),
       ("2016_QCDmadgraph-HT1500to2000",119.9/11839357),
       ("2016_QCDmadgraph-HT2000toInf",25.24/5947849),
       ]
   mc["2017"]=[("2017_QCDmadgraph-HT200to300",1545000./58990434),
       ("2017_QCDmadgraph-HT300to500",323300./58748739),
       ("2017_QCDmadgraph-HT500to700",30000./54366431),
       ("2017_QCDmadgraph-HT700to1000",6324./46924322),
       ("2017_QCDmadgraph-HT1000to1500",1090./16495598),
       ("2017_QCDmadgraph-HT1500to2000",101./11196479),
       ("2017_QCDmadgraph-HT2000toInf",20.43/5362513),
       ]
   mc["2018"]=[("2018_QCDmadgraph-HT200to300",1461000./54289442),
       ("2018_QCDmadgraph-HT300to500",311900./54512704),
       ("2018_QCDmadgraph-HT500to700",29070./53919811),
       ("2018_QCDmadgraph-HT700to1000",5962./48158738),
       ("2018_QCDmadgraph-HT1000to1500",1005./14945819),
       ("2018_QCDmadgraph-HT1500to2000",101.8/10707847),
       ("2018_QCDmadgraph-HT2000toInf",20.54/5329144),
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
   
   for prefix in ["QCDmadgraph"]:
     sum_in_quadrature_up=[]
     sum_in_quadrature_down=[]
     for mass in range(len(massbins)):
       sum_in_quadrature_up+=[len(chi_binnings[mass])*[0]]
       sum_in_quadrature_down+=[len(chi_binnings[mass])*[0]]

     canvas = TCanvas("PU","PU",0,0,800,600)
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
         print dire+"datacard_shapelimit13TeV_"+prefix+"_PU_"+str(year)+"_"+b+"_chi.root"
     	 f_refmc[year+b]=TFile.Open(dire+"datacard_shapelimit13TeV_"+prefix+"_PU_"+str(year)+"_"+b+"_chi.root")
         h=f_refmc[year+b].Get(prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_rebin1_backup")
         scalexsec[b]=1.#mc[year][int(b)-1][1]
         files+=[f_refmc[year+b]]
         if b=="1":
     	   histyear[year]=h
           histyear[year]=histyear[year].Clone(histyear[year].GetName()+year+"main")
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
       histref=hist.Clone(hist.GetName()+"ref")
       hists+=[histref]
       for b in range(hist.GetNbinsX()):
           hist.SetBinError(b+1,0)
           hist.SetBinContent(b+1,1)
       hist.GetXaxis().SetTitle(label)
       hist.GetYaxis().SetTitle("N")
       hist.GetYaxis().SetRangeUser(0.98,1.03)
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

       updown="PU_Up"
       histyear2={}
       for year in ["2016","2017","2018"]:
        for b in ["1","2","3","4","5","6","7"]:
         print prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+updown+"_rebin1"
         h=f_refmc[year+b].Get(prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+updown+"_rebin1")
         if b=="1":
     	    histyear2[year]=h
            histyear2[year]=histyear2[year].Clone(histyear2[year].GetName()+year+"up")
            histyear2[year].Scale(1./scalexsec[b])
         else:
     	    histyear2[year].Add(h,1./scalexsec[b])
        histyear2[year]=histyear2[year].Rebin(len(chi_binnings[mass])-1,histyear2[year].GetName()+"_rebin1",chi_binnings[mass])
        hists+=[histyear2[year]]
       hist2=histyear2["2016"].Clone(histyear2["2016"].GetName()+"combine")
       hists+=[hist2]
       hist2.Scale(lumifactor["2016"])
       hist2.Add(histyear2["2017"],lumifactor["2017"])
       hist2.Add(histyear2["2018"],lumifactor["2018"])
       if hist2.Integral()>0:
         hist2.Scale(histref.Integral()/hist2.Integral())
       hist2.Divide(hist2,histref)
       hist2.SetLineWidth(1)
       hist2.SetLineColor(colors[0])
       hist2.SetLineStyle(2)
       hist2.SetTitle("")
       hist2.SetStats(False)
       fit2=TF1(hist2.GetName()+"smooth","pol1",1,16)
       hist2.Fit(fit2,"NQ")
       #for chi_bin in range(len(chi_binnings[mass])):
       #  hist2.SetBinContent(chi_bin+1,fit2.Eval(hist2.GetBinCenter(chi_bin+1)))
       hist2.Draw("hesame")
       fit2.Draw("lsame")
       legend.AddEntry(hist2,"PU","l")
       for chi_bin in range(len(chi_binnings[mass])):
         if (hist2.GetBinContent(chi_bin+1)-1.0)*(hist2.GetBinCenter(chi_bin+1)-8.5)>0:
          sum_in_quadrature_up[mass][chi_bin]=sqrt(pow(sum_in_quadrature_up[mass][chi_bin],2)+pow(hist2.GetBinContent(chi_bin+1)-1.0,2))
         else:
          sum_in_quadrature_down[mass][chi_bin]=sqrt(pow(sum_in_quadrature_down[mass][chi_bin],2)+pow(hist2.GetBinContent(chi_bin+1)-1.0,2))

       updown="PU_Down"
       histyear3={}
       for year in ["2016","2017","2018"]:
        for b in ["1","2","3","4","5","6","7"]:
         print prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+updown+"_rebin1"
         h=f_refmc[year+b].Get(prefix+"#chi"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+"_"+updown+"_rebin1")
         if b=="1":
     	    histyear3[year]=h
            histyear3[year]=histyear3[year].Clone(histyear3[year].GetName()+year+"down")
            histyear3[year].Scale(1./scalexsec[b])
         else:
     	    histyear3[year].Add(h,1./scalexsec[b])
        histyear3[year]=histyear3[year].Rebin(len(chi_binnings[mass])-1,histyear3[year].GetName()+"_rebin1",chi_binnings[mass])
        hists+=[histyear3[year]]
       hist3=histyear3["2016"].Clone(histyear3["2016"].GetName()+"combine")
       hists+=[hist3]
       hist3.Scale(lumifactor["2016"])
       hist3.Add(histyear3["2017"],lumifactor["2017"])
       hist3.Add(histyear3["2018"],lumifactor["2018"])
       if hist3.Integral()>0:
         hist3.Scale(histref.Integral()/hist3.Integral())
       hist3.Divide(hist3,histref)
       hist3.SetLineWidth(1)
       hist3.SetLineColor(colors[0])
       hist3.SetLineStyle(3)
       hist3.SetTitle("")
       hist3.SetStats(False)
       fit3=TF1(hist3.GetName()+"smooth","pol1",1,16)
       hist3.Fit(fit3,"NQ")
       #for chi_bin in range(len(chi_binnings[mass])):
       #  hist3.SetBinContent(chi_bin+1,fit3.Eval(hist3.GetBinCenter(chi_bin+1)))
       hist3.Draw("hesame")
       fit3.Draw("lsame")
       for chi_bin in range(len(chi_binnings[mass])):
         if (hist3.GetBinContent(chi_bin+1)-1.0)*(hist3.GetBinCenter(chi_bin+1)-8.5)>0:
          sum_in_quadrature_up[mass][chi_bin]=sqrt(pow(sum_in_quadrature_up[mass][chi_bin],2)+pow(hist3.GetBinContent(chi_bin+1)-1.0,2))
         else:
          sum_in_quadrature_down[mass][chi_bin]=sqrt(pow(sum_in_quadrature_down[mass][chi_bin],2)+pow(hist3.GetBinContent(chi_bin+1)-1.0,2))

       legend.SetTextSize(0.04)
       legend.SetFillStyle(0)
       legend.Draw("same")

     canvas.SaveAs(dire+"chi_systematic_plots"+var+"_"+prefix+"PU_13TeV_run2.root")
     canvas.SaveAs(dire+"chi_systematic_plots"+var+"_"+prefix+"PU_13TeV_run2.pdf")
