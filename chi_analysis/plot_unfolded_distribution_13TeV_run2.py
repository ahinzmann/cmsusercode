import os, sys
from ROOT import * 
import array

#gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
#gROOT.Reset()
gROOT.SetStyle("Plain")
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetTitleOffset(1.2,"Y")
gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadBottomMargin(0.15)
gStyle.SetPadTopMargin(0.03)
gStyle.SetPadRightMargin(0.05)
gStyle.SetMarkerSize(1.5)
gStyle.SetHistLineWidth(1)
gStyle.SetStatFontSize(0.020)
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetNdivisions(510, "XYZ")
gStyle.SetLegendBorderSize(0)

def rebin(h1,nbins,binning):
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)*h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)*h1.GetBinWidth(b+1))
    h1=h1.Rebin(nbins,h1.GetName()+"_rebin",binning)
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
    return h1

if __name__ == '__main__':

    run="2"
    prefix="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/versions/run"+run+"ULNNLO_pt12/datacard_shapelimit13TeV"

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

    canvas = TCanvas("","",0,0,600,600)
    canvas.Divide(3,3)
    plots=[]
    legends=[]
    new_hists=[]

    filename="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/multidimfit.root"
    print filename
    fitfile = TFile.Open(filename)
    fittree=fitfile.Get("fit_mdf")
    fitParameters=[]
    fitConstraints=[]

    for i in range(len(massbins)):
      filename=prefix+"_GEN-QCD-run2_chi.root"
      print filename
      f = TFile.Open(filename)
      new_hists+=[f]
      canvas.cd(i+1)
      legend1=TLegend(0.2,0.6,0.9,0.95,(str(massbins[i][0])+"<m_{jj}<"+str(massbins[i][1])+" GeV").replace("7000<m_{jj}<13000","m_{jj}>7000"))
      legends+=[legend1]

      histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmear"
      print histname
      h1gen=f.Get(histname)
      plots+=[h1gen]
      h1genpostfit=h1gen.Clone(h1gen.GetName()+"postfit")
      plots+=[h1genpostfit]
            
      first=True
      for j in range(len(massbins)):
        for chibin in range(len(chi_bins[j+mass_bin_offset])-1):
          histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_bin_"+str(j+mass_bin_offset)+"_"+str(chibin)+"_"
          name="r_Bin_"+str(j+mass_bin_offset)+"_"+str(chibin)
          fitParameter=fittree.floatParsFinal().find(name).getVal()
          fitConstraint=fittree.floatParsFinal().find(name).getError()
          print histname,fitParameter,fitConstraint

          if i==j:
            h1genpostfit.SetBinContent(chibin+1,h1genpostfit.GetBinContent(chibin+1)*fitParameter)
            h1genpostfit.SetBinError(chibin+1,h1genpostfit.GetBinContent(chibin+1)*fitConstraint/fitParameter)

          if first:
            h1=f.Get(histname)
            plots+=[h1]
            h1postfit=h1.Clone(h1.GetName()+"postfit")
            h1postfit.Scale(fitParameter)
            h1postfit.SetBinError(chibin+1,h1postfit.GetBinContent(chibin+1)*fitConstraint/fitParameter)
            plots+=[h1postfit]
            first=False
          else:
            h1add=f.Get(histname)
            h1.Add(h1add)
            h1addpostfit=h1add.Clone(h1add.GetName()+"postfit")
            h1addpostfit.Scale(fitParameter)
            h1addpostfit.SetBinError(chibin+1,h1addpostfit.GetBinContent(chibin+1)*fitConstraint/fitParameter)
            h1postfit.Add(h1addpostfit)

      histname='data_obs#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
      print histname
      h2=f.Get(histname)
      plots+=[h2]
            
      histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
      print histname
      h3=f.Get(histname)
      plots+=[h3]
            
      h1.Scale(h3.Integral()/h1.Integral())
      h1postfit.Scale(h3.Integral()/h1postfit.Integral())
      h1gen.Scale(h3.Integral()/h1gen.Integral())
      h1genpostfit.Scale(h3.Integral()/h1genpostfit.Integral())
      h2.Scale(h3.Integral()/h2.Integral())

      for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
        h1postfit.SetBinContent(b+1,h1postfit.GetBinContent(b+1)/h1postfit.GetBinWidth(b+1))
        h1postfit.SetBinError(b+1,h1postfit.GetBinError(b+1)/h1postfit.GetBinWidth(b+1))
        h1gen.SetBinContent(b+1,h1gen.GetBinContent(b+1)/h1gen.GetBinWidth(b+1))
        h1gen.SetBinError(b+1,h1gen.GetBinError(b+1)/h1gen.GetBinWidth(b+1))
        h1genpostfit.SetBinContent(b+1,h1genpostfit.GetBinContent(b+1)/h1genpostfit.GetBinWidth(b+1))
        h1genpostfit.SetBinError(b+1,h1genpostfit.GetBinError(b+1)/h1genpostfit.GetBinWidth(b+1))
        h2.SetBinContent(b+1,h2.GetBinContent(b+1)/h2.GetBinWidth(b+1))
        h2.SetBinError(b+1,h2.GetBinError(b+1)/h2.GetBinWidth(b+1))
        h3.SetBinContent(b+1,h3.GetBinContent(b+1)/h3.GetBinWidth(b+1))
        h3.SetBinError(b+1,h3.GetBinError(b+1)/h3.GetBinWidth(b+1))

      h2.SetLineColor(1)
      h2.SetTitle("")
      h2.Draw("he")
      h2.GetYaxis().SetRangeUser(0,h2.GetMaximum()*1.5)

      h1.SetLineColor(2)
      h1.SetTitle("")
      h1.Draw("hesame")
      
      h3.SetLineColor(6)
      h3.SetLineStyle(3)
      #h3.Draw("hesame")

      h1postfit.SetLineColor(2)
      h1postfit.SetLineStyle(2)
      h1postfit.Draw("hesame")
      
      h1gen.SetLineColor(4)
      h1gen.Draw("hesame")
      
      h1genpostfit.SetLineColor(4)
      h1genpostfit.SetLineStyle(2)
      h1genpostfit.Draw("hesame")
      
      # PLOTS
      legend1.AddEntry(h2,"Raw data","le")
      legend1.AddEntry(h1,"Smeared-Prefit","le")
      #legend1.AddEntry(h3,"Smeared-Prefit","le")
      legend1.AddEntry(h1postfit,"Smeared-Postfit","le")
      legend1.AddEntry(h1gen,"Gen","le")
      legend1.AddEntry(h1genpostfit,"Unfolded data","le")

      legend1.SetTextSize(0.04)
      legend1.SetFillStyle(0)
      legend1.Draw("same")

    canvas.SaveAs(prefix + "_unfolded.pdf")
    canvas.SaveAs(prefix + "_unfolded.eps")
