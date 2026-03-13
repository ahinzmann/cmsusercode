import os, sys
from ROOT import * 
import array
from math import *

#gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
#gROOT.Reset()
gROOT.SetStyle("Plain")
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetTitleOffset(1.2,"Y")
gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadBottomMargin(0.15)
gStyle.SetPadTopMargin(0.075)
gStyle.SetPadRightMargin(0.05)
gStyle.SetMarkerSize(1.5)
gStyle.SetHistLineWidth(1)
gStyle.SetStatFontSize(0.020)
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetNdivisions(510, "XYZ")
gStyle.SetLegendBorderSize(0)

if __name__ == '__main__':

    trivialClosure=False
    withUncertainties=True
    run="2"
    prefix="/data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/versions/run"+run+"ULNNLO_m2_NNPDF3/datacard_shapelimit13TeV"
    preliminary=False

    name="unfold"
    if withUncertainties: name+="_withUncertainties"
    if trivialClosure: name+="_trivialClosure"
    print(name)

    massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    all_mass_bins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    mass_bin_offset=len(all_mass_bins)-len(massbins)

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

    plots=[]
    legends=[]
    new_hists=[]

    filename="/data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/multidimfit_"+name+"_run"+run+".root"
    print(filename)
    fitfile = TFile.Open(filename)
    fittree=fitfile.Get("fit_mdf")
    fitParameters=[]
    fitConstraints=[]
    
    name+="_correlationMatrix"

    pois=fittree.floatParsFinal().Clone()
    pois.removeAll()
    fitted_pars=fittree.floatParsFinal()
    print("All pars", len(fitted_pars))
    num=0
    for par in range(len(fitted_pars)):
      #print num, fitted_pars[par].GetName()
      if "r_Bin_0" in fitted_pars[par].GetName(): # skip first mass bin
        continue
      if "r_Bin" in fitted_pars[par].GetName():
        pois.add(fitted_pars[par])
        num+=1
    print("POIs", len(pois))
    reducedCovarianceMatrix=fittree.reducedCovarianceMatrix(pois)
    correlationMatrix=reducedCovarianceMatrix.Clone()
    for i in range(len(pois)):
      for j in range(len(pois)):
        if i!=j:
          correlationMatrix[i][j]/=sqrt(correlationMatrix(i,i)*correlationMatrix(j,j))
        #print(correlationMatrix[i][j])
    for i in range(len(pois)):
      correlationMatrix[i][i]=1.

    canvas = TCanvas("corr","corr",0,0,600,600)
    canvas.SetRightMargin(0.15)
    gStyle.SetNumberContours(200)
    #gStyle.SetPalette(104)
    TColor.CreateGradientColorTable(3,array.array('d', (0,0.5,1.0)),array.array('d', (1,1,0)),array.array('d', (0,1,0)),array.array('d', (0,1,1)),50)
    h=TH2D(correlationMatrix)
    h.GetZaxis().SetRangeUser(-1,1)
    h.SetName("correlation")
    h.GetXaxis().SetTitle("Bin number")
    h.GetYaxis().SetTitle("Bin number")
    h.Draw("")
    h.Draw("colz2same")
    #h.GetYaxis().SetTitle("Correlation coefficient")

    # CMS
    leg1=TLatex(0.05,78+0.05,"#bf{CMS}")
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.055)
    if preliminary:
      leg2=TLatex(13+0.05,78+0.05,"#it{Preliminary}")
    else:
      leg2=TLatex(13+0.05,78+0.03,"#it{Supplementary}")
    #cmsPos=min_x+220/1000.
    #leg2=TLatex(cmsPos,ymax-0.17,"#bf{CMS}")
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.042)
    # lumi
    lumiPos=78*0.6
    leg3=TLatex(lumiPos,78+0.03,"138 fb^{-1} (13 TeV)")
    leg3.SetTextFont(42)
    leg3.SetTextSize(0.040)
    leg1.Draw("same")
    leg2.Draw("same")
    leg3.Draw("same")

    canvas.SaveAs(prefix + "_"+name+"_run"+run+".pdf")
    canvas.SaveAs(prefix + "_"+name+"_run"+run+".root")
