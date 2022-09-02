import os, sys
import array
from ROOT import * 

print "start ROOT"
#gROOT.Reset()
gROOT.SetStyle("Plain")
gROOT.SetBatch(True)
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

   label="#chi"

   masses=[2400,3000,3600,4200,4800,5400,6000,6600,13000]

   colors=[1,2,3,4,6,7,8,9,10,11,12,13]
   styles=[1,2,3,4,5,6,7,8,9,11,12,13]
   
   chi_bins=[(1,2,3,4,5,6,7,8,9,10,12,14,16),
             (1,2,3,4,5,6,7,8,9,10,12,14,16),
             (1,2,3,4,5,6,7,8,9,10,12,14,16),
             (1,2,3,4,5,6,7,8,9,10,12,14,16),
             (1,2,3,4,5,6,7,8,9,10,12,14,16),
             (1,2,3,4,5,6,7,8,9,10,12,14,16),
             (1,2,3,4,5,6,7,8,9,10,12,14,16),
             (1,3,5,7,10,12,14,16),
             ]
   chi_binnings=[]
   for mass_bin in chi_bins:
        chi_binnings+=[array.array('d')]
        for chi_bin in mass_bin:
            chi_binnings[-1].append(chi_bin)

   print "open file"
            
   infile=TFile.Open('fastnlo/RunII/DijetAngularCMS13_ewk.root')
   
   print "open canvas"
            
   canvas = TCanvas("ewk","ewk",0,0,400,400)
   legends=[]
   hists=[]
   legend=TLegend(0.2,0.55,0.95,0.90)
   legends+=[legend]
   for mass in range(len(masses)-1):
     var="chi"
     hist=infile.Get("chi-"+str(masses[mass])+"-"+str(masses[mass+1]))
     hists+=[hist]
     #hist=hist.Rebin(len(chi_binnings[mass])-1,hist.GetName()+"_rebin1",chi_binnings[mass])
     hist.SetLineWidth(2)
     hist.SetLineColor(colors[mass])
     hist.SetLineStyle(styles[mass])
     for b in range(hist.GetNbinsX()):
         hist.SetBinError(b+1,0)
     hist.GetXaxis().SetTitle(label)
     hist.GetYaxis().SetTitle("EWK correction factor")
     hist.GetYaxis().SetRangeUser(0.9,1.3)
     hist.GetXaxis().SetTitleOffset(1.1)
     hist.GetYaxis().SetTitleOffset(1.1)
     hist.GetXaxis().SetLabelSize(0.05)
     hist.GetYaxis().SetLabelSize(0.05)
     hist.GetXaxis().SetTitleSize(0.06)
     hist.GetYaxis().SetTitleSize(0.06)
     hist.SetTitle("")
     hist.SetStats(False)
     if mass==0:
       hist.Draw("l")
     else:
       hist.Draw("lsame")  
     legend.AddEntry(hist,str(masses[mass])+"-"+str(masses[mass+1]),"l")

   legend.SetTextSize(0.04)
   legend.SetFillStyle(0)
   legend.Draw("same")

   canvas.SaveAs("chi_ewk_plots.root")
   canvas.SaveAs("chi_ewk_plots.pdf")
