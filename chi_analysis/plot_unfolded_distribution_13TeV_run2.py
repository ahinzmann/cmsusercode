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

    trivialClosure=False
    crossClosure=False
    onlyMC=False
    onlyData=True
    withUncertainties=True
    run="2"
    prefix="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/versions/run"+run+"ULNNLO_m2/datacard_shapelimit13TeV"

    name="unfold"
    if withUncertainties: name+="_withUncertainties"
    if trivialClosure: name+="_trivialClosure"
    if crossClosure: name+="_crossClosure"
    print name

    massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]

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

    canvas = TCanvas("","",0,0,600,600)
    canvas.Divide(3,3)
    plots=[]
    legends=[]
    new_hists=[]

    filename="/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/multidimfit_"+name+"_run"+run+".root"
    print filename
    fitfile = TFile.Open(filename)
    fittree=fitfile.Get("fit_mdf")
    fitParameters=[]
    fitConstraints=[]

    fout=TFile.Open("/nfs/dust/cms/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/datacards/datacard_shapelimit13TeV_"+name+"_run"+run+".root","RECREATE")

    if onlyMC: name+="_MC"
    if onlyData: name+="_data"

    for i in range(len(massbins)):
      filename=prefix+"_GEN-QCD-run2_chi.root"
      print filename
      f = TFile.Open(filename)
      new_hists+=[f]
      canvas.cd(i+1)
      legend1=TLegend(0.3,0.6,0.9,0.95,(str(massbins[i][0])+"<m_{jj}<"+str(massbins[i][1])+" GeV").replace("7000<m_{jj}<13000","m_{jj}>7000"))
      legends+=[legend1]

      if trivialClosure or crossClosure: # use LO QCD instead of NNLO QCD
        histname='QCD#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmear"
      else:
        histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmear"
      print histname
      h1gen=f.Get(histname)
      plots+=[h1gen]
      if crossClosure: # use LO QCD instead of NNLO QCD
        histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmear"
        print histname
        h1genpostfit=f.Get(histname).Clone(histname+"postfit")
      else:
        h1genpostfit=h1gen.Clone(h1gen.GetName()+"postfit")
      plots+=[h1genpostfit]
            
      first=True
      for j in range(len(massbins)):
        for chibin in range(len(chi_bins[j])-1):
          if trivialClosure or crossClosure: # use LO QCD instead of NNLO QCD
            histname='QCD#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_bin_"+str(j)+"_"+str(chibin)+"_"
          else:
            histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_bin_"+str(j)+"_"+str(chibin)+"_"
          binname="r_Bin_"+str(j)+"_"+str(chibin)
          fitParameter=fittree.floatParsFinal().find(binname).getVal()
          fitConstraint=fittree.floatParsFinal().find(binname).getError()
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

      fout.cd()
      h1genpostfit.Write()

      histname='data_obs#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
      print histname
      h2=f.Get(histname)
      plots+=[h2]
            
      histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
      print histname
      h3=f.Get(histname)
      plots+=[h3]
            
      h1.Scale(1./h1.Integral())
      h1postfit.Scale(1./h1postfit.Integral())
      h1gen.Scale(1./h1gen.Integral())
      h1genpostfit.Scale(1./h1genpostfit.Integral())
      h2.Scale(1./h2.Integral())
      h3.Scale(1./h3.Integral())

      for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,0)
        h1postfit.SetBinContent(b+1,h1postfit.GetBinContent(b+1)/h1postfit.GetBinWidth(b+1))
        h1postfit.SetBinError(b+1,0)
        h1gen.SetBinContent(b+1,h1gen.GetBinContent(b+1)/h1gen.GetBinWidth(b+1))
        h1gen.SetBinError(b+1,0)
        h1genpostfit.SetBinContent(b+1,h1genpostfit.GetBinContent(b+1)/h1genpostfit.GetBinWidth(b+1))
        h1genpostfit.SetBinError(b+1,h1genpostfit.GetBinError(b+1)/h1genpostfit.GetBinWidth(b+1))
        h2.SetBinContent(b+1,h2.GetBinContent(b+1)/h2.GetBinWidth(b+1))
        h2.SetBinError(b+1,h2.GetBinError(b+1)/h2.GetBinWidth(b+1))
        h3.SetBinContent(b+1,h3.GetBinContent(b+1)/h3.GetBinWidth(b+1))
        h3.SetBinError(b+1,0)

      
      h1gen.SetLineColor(4)
      if not onlyData:
        h1gen.Draw("he")

      if not trivialClosure and not crossClosure and not onlyMC:
        h2.SetLineColor(4)
        h2.SetTitle("")
        if onlyData:
          h2.Draw("he")
          h1gen=h2
        else:
          h2.Draw("hesame")

      if i>=7:
          h1gen.GetYaxis().SetRangeUser(0.00,0.25)
      elif i>=6:
          h1gen.GetYaxis().SetRangeUser(0.02,0.22)
      elif i>=4:
          h1gen.GetYaxis().SetRangeUser(0.045,0.12)
      else:
          h1gen.GetYaxis().SetRangeUser(0.045,0.12)
      
      if not onlyMC:
        h1genpostfit.SetLineColor(1)
        h1genpostfit.SetLineStyle(2)
        h1genpostfit.SetLineWidth(2)
        h1genpostfit.Draw("hesame")
      
      if not trivialClosure and not crossClosure and not onlyData:
        h1.SetLineColor(2)
        h1.SetTitle("")
        h1.Draw("hesame")
      
        h3.SetLineColor(6)
        h3.SetLineStyle(3)
        #h3.Draw("hesame")

      if not trivialClosure and not crossClosure and not onlyMC and not onlyData:
        h1postfit.SetLineColor(2)
        h1postfit.SetLineStyle(2)
        h1postfit.Draw("hesame")
      
      # PLOTS
      if not onlyData:
        legend1.AddEntry(h1gen,"Gen","l")
      if not trivialClosure and not crossClosure and not onlyMC:
        legend1.AddEntry(h2,"Reco data","le")
      if trivialClosure or crossClosure:
        legend1.AddEntry(h1genpostfit,"Unfolded Smeared","le")
      elif not onlyMC:
        legend1.AddEntry(h1genpostfit,"Unfolded data","le")
      if not trivialClosure and not crossClosure and not onlyData:
        legend1.AddEntry(h1,"Smeared-Prefit","l")
        #legend1.AddEntry(h3,"Smeared-Prefit","l")
      if not trivialClosure and not crossClosure and not onlyMC and not onlyData:
        legend1.AddEntry(h1postfit,"Smeared-Postfit","l")

      legend1.SetTextSize(0.04)
      legend1.SetFillStyle(0)
      legend1.Draw("same")

    canvas.SaveAs(prefix + "_"+name+"_run"+run+".pdf")
    #canvas.SaveAs(prefix + "_unfolded.eps")
    fout.Close()
