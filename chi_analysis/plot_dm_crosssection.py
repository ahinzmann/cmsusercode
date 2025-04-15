from ROOT import *
import ROOT
import array, math
import os
from math import *

#! /usr/bin/env python         

if __name__=="__main__":

  for mDM in [1,3000]:
  
    print("start ROOT")
    gROOT.Reset()
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

    print("start CMS_lumi")

    #gROOT.LoadMacro("CMS_lumi.C");
    #iPeriod = 4;	#// 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV 
    #iPos = 11;
    
    c = TCanvas("dmwidth", "dmwidth", 0, 0, 300, 300)
    c.SetLogy()
    mg=TMultiGraph()
    colors=[1,2,3,4,6,7,8,9,15,1,2,3,4,6,7,8,9,15]
    styles=[20,21,22,23,25,26,27,28,29,30,32,33,34]
    
    if mDM==1:
      l3=TLegend(0.45,0.6,0.85,0.9,"g_{DM}=1, m_{DM}=1 GeV")
    else:
      l3=TLegend(0.45,0.6,0.85,0.9,"g_{DM}=1, m_{DM}>m_{Med}")
    l3.SetTextSize(0.05)
    
    style=0
    for scenario in range(4):
       pointsX=array.array("f",[1,1.5,2,2.5,3,4,5,6,7])
       if mDM==1:
        if scenario==0:
         pointsY=array.array("f",[117.4,23.38,6.547,2.210,0.886,0.184,0.0559,0.0229,0.0112])
        elif scenario==1:
         pointsY=array.array("f",[116.1,23.02,6.509,2.242,0.877,0.189,0.0552,0.0228,0.0112])
        elif scenario==2:
         pointsY=array.array("f",[523.4,111.2,34.16,12.96,5.856,1.606,0.590,0.271,0.145])
        else:
         pointsY=array.array("f",[509.8,109.6,33.65,13.24,5.872,1.610,0.608,0.270,0.145])
       else:
        if scenario==0:
         pointsX=array.array("f",[1,1.5,2.5,3,5,6])
         pointsY=array.array("f",[142.2,28.40,2.631,1.015,0.0596,0.0232])
        elif scenario==1:
         pointsX=array.array("f",[1.5,2,3,4,6])
         pointsY=array.array("f",[28.16,7.749,1.035,0.212,0.0232])
        elif scenario==2:
         pointsX=array.array("f",[1,2,2.5,4,5])
         pointsY=array.array("f",[550.69,36.22,13.77,1.681,0.618])
        else:
         pointsX=array.array("f",[1,1.5,2.5,3,5,6])
         pointsY=array.array("f",[537.4,115.9,13.78,6.111,0.569,0.275])
       g=TGraph(len(pointsX),pointsX,pointsY)
       #g.SetMarkerStyle(styles[style])
       #g.SetMarkerSize(1)
       g.SetLineColor(colors[style])
       g.SetLineStyle(style+1)
       g.SetLineWidth(2)
       mg.Add(g)
       if scenario==0:
         l3.AddEntry(g,"Axial g_{q}=0.5","pl")
       elif scenario==1:
         l3.AddEntry(g,"Vector g_{q}=0.5","pl")
       elif scenario==2:
         l3.AddEntry(g,"Axial g_{q}=1.0","pl")
       else:
         l3.AddEntry(g,"Vector g_{q}=1.0","pl")
       style+=1
    
    mg.Draw("al")
    mg.SetTitle("")
    mg.GetXaxis().SetTitle("m_{Med} (TeV)")
    mg.GetYaxis().SetTitle("Cross section (pb)")
    mg.GetXaxis().SetRangeUser(1,7)
    mg.GetYaxis().SetRangeUser(0.01,1000)

    l3.SetFillStyle(0)
    l3.Draw("same")
    
    #// writing the lumi information and the CMS "logo"
    #CMS_lumi( c, iPeriod, iPos );
    
    c.SaveAs("dm_crosssection"+str(mDM)+".pdf")
    c.SaveAs("dm_crosssection"+str(mDM)+".eps")
