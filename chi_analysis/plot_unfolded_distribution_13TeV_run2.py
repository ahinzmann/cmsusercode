import os, sys
from ROOT import gROOT,gStyle,TCanvas,TFile,TLegend
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

    trivialClosure=False # Pythia vs Pythia templates+response
    crossClosureNNLO=False # NNLO vs Pythia templates
    crossClosureHerwig=False # Herwig vs Pythia response
    crossClosureMadgraph=False # Madgraph vs Pythia response
    modelUncertainty=False # Herwig vs Madgraph vs Pythia smearing
    modelComparison=False # NNLO vs Pythia vs Madgraph and generator level
    onlyMC=False # before/after smearing
    onlyData=False # before/after unfolding
    withUncertainties=False
    compareNoUncertainty=False
    normalize=False
    run="2"
    postfix=""
    if not normalize:
      postfix="_nonorm"
    prefix="/data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/versions/run"+run+"ULNNLO_m2_NNPDF3"+postfix+"/datacard_shapelimit13TeV"

    name="unfold"
    if withUncertainties: name+="_withUncertainties"
    if trivialClosure: name+="_trivialClosure"
    if crossClosureNNLO: name+="_crossClosureNNLO"
    if crossClosureHerwig: name+="_crossClosureHerwig"
    if crossClosureMadgraph: name+="_crossClosureMadgraph"
    print(name)

    if not trivialClosure and not crossClosureNNLO and not crossClosureHerwig and not crossClosureMadgraph and not onlyMC and not modelUncertainty and not modelComparison:
      massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    else:
      massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]

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

    filename="/data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/multidimfit_"+name+"_run"+run+".root"
    print(filename)
    fitfile = TFile.Open(filename)
    fittree=fitfile.Get("fit_mdf")
    fitParameters=[]
    fitConstraints=[]

    if compareNoUncertainty:
      filenamed=("/data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/datacards/datacard_shapelimit13TeV_"+name.replace("_withUncertainties","")+"_run"+run+".root")
      print(filenamed)
      fd = TFile.Open(filenamed)
      new_hists+=[fd]
    
    fout=TFile.Open("/data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/datacards/datacard_shapelimit13TeV_"+name+postfix+"_run"+run+".root","RECREATE")

    if onlyMC: name+="_MC"
    if onlyData: name+="_data"
    if modelUncertainty: name+="_modelUncertainty"
    if modelComparison: name+="_modelComparison"
    name+=postfix

    filename=prefix+"_GEN-QCD-run2_chi.root"
    print(filename)
    f = TFile.Open(filename)
    new_hists+=[f]
    
    filenameh=(prefix+"_GEN-QCD-run2_chi.root").replace("shapelimit13TeV","shapelimit13TeVherwigpp")
    print(filenameh)
    fh = TFile.Open(filenameh)
    new_hists+=[fh]
    
    filenamem=(prefix+"_GEN-QCD-run2_chi.root").replace("shapelimit13TeV","shapelimit13TeVmadgraphMLM")
    print(filenamem)
    fm = TFile.Open(filenamem)
    new_hists+=[fm]

    filenameu='datacards/datacard_shapelimit13TeV_unfold_withUncertainties'+postfix+'_run2.root'
    print(filenameu)
    fDataUnfolded = TFile.Open(filenameu)

    chi2detector=0
    chi2unfolded=0
    normsh1=1
    normsh1postfit=1
    normsh1gen=1
    normsh1genpostfit=1
    normsh2gen=1
    normsh3gen=1
    normsh4gen=1
    normsh2=1
    normsh3=1
    normsh4=1
   
    for i in range(len(massbins)):
      canvas.cd(i+1)
      legend1=TLegend(0.3,0.6,0.9,0.95,(str(massbins[i][0])+"<m_{jj}<"+str(massbins[i][1])+" GeV").replace("7000<m_{jj}<13000","m_{jj}>7000"))
      legends+=[legend1]

      if trivialClosure or crossClosureNNLO or crossClosureHerwig or crossClosureMadgraph: # use Pythia instead of NNLO QCD
        histname='QCD#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmear"
      elif not normalize:
        histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nonorm_nosmear"
      else:
        histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmear"
      print(histname)
      h1gen=f.Get(histname)
      plots+=[h1gen]
      if modelComparison:
        #histnameu="QCD_ALT#chi"+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmearpostfit"
        #print(histnameu)
        #h1dataUnfolded=fDataUnfolded.Get(histnameu)
        histname='QCD#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmear"
        print(histname)
        h1gen2=f.Get(histname)
        plots+=[h1gen2]
        sampleMadgraph="plots/chi_model_plots_chi_"+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_UL_run2_GEN.root"
        print(sampleMadgraph)
        inMadgraph=TFile(sampleMadgraph,'READ')
        new_hists+=[inMadgraph]
        cv=inMadgraph.Get("cv")
        new_hists+=[cv]
        h=cv.GetListOfPrimitives()[0].GetListOfPrimitives()[1]
        print(h.GetName())
        h1gen3=h.Clone("mg"+str(i))
        plots+=[h1gen3]
        h=cv.GetListOfPrimitives()[0].GetListOfPrimitives()[3]
        print(h.GetName())
        h1gen4=h.Clone("hw"+str(i))
        plots+=[h1gen4]
        for b in range(h1gen3.GetXaxis().GetNbins()):
          h1gen3.SetBinContent(b+1,h1gen3.GetBinContent(b+1)*h1gen3.GetBinWidth(b+1))
          h1gen4.SetBinContent(b+1,h1gen4.GetBinContent(b+1)*h1gen4.GetBinWidth(b+1))
        canvas.cd(i+1)
      if crossClosureNNLO: # use NNLO smeared with Pythia response
        histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmear"
        print(histname)
        h1genpostfit=f.Get(histname).Clone(histname+"postfit")
      elif crossClosureHerwig: # use Pythia smeared with Herwig-response instead of NNLO QCD
        histname='QCD#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmear"
        print(histname)
        h1genpostfit=fh.Get(histname).Clone(histname+"postfit")
      elif crossClosureMadgraph: # use Pythia smeared with Madgraph-response instead of NNLO QCD
        histname='QCD#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmear"
        print(histname)
        h1genpostfit=fm.Get(histname).Clone(histname+"postfit")
      else:
        h1genpostfit=h1gen.Clone(h1gen.GetName()+"postfit")
      plots+=[h1genpostfit]
            
      first=True
      for j in range(len(massbins)):
        for chibin in range(len(chi_bins[j])-1):
          if trivialClosure or crossClosureNNLO or crossClosureHerwig or crossClosureMadgraph: # use LO QCD instead of NNLO QCD
            histname='QCD#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_bin_"+str(j)+"_"+str(chibin)+"_"
          else:
            histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_bin_"+str(j)+"_"+str(chibin)+"_"
          binname="r_Bin_"+str(j)+"_"+(str(chibin) if chibin>9 else "0"+str(chibin))
          fitParameter=fittree.floatParsFinal().find(binname).getVal()
          fitConstraint=fittree.floatParsFinal().find(binname).getError()
          print(histname,fitParameter,fitConstraint)

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

      #if not normalize:
      #  histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nonorm"
      #  print(histname)
      #  h1=f.Get(histname)
      #  plots+=[h1]
      #  histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
      #  print(histname)
      #  h1genpostfit.Scale(f.Get(histname).Integral()/h1.Integral())

      #if modelComparison:
      #  h1genpostfit=h1dataUnfolded

      histname='data_obs#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
      print(histname)
      h2=f.Get(histname)
      plots+=[h2]
            
      histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
      print(histname)
      h3=fh.Get(histname)
      plots+=[h3]
      
      histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
      print(histname)
      h4=fm.Get(histname)
      plots+=[h4]

      if compareNoUncertainty:
        histname='QCD_ALT#chi'+str(massbins[i]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_nosmearpostfit"
        print(histname)
        h5=fd.Get(histname)
        plots+=[h5]
      
      if normalize or i==0:
        normsh1=h1.Integral()
        normsh1postfit=h1postfit.Integral()
        normsh1gen=h1gen.Integral()
        if modelComparison:
          normsh1gen2=h1gen2.Integral()
          normsh1gen3=h1gen3.Integral()
          normsh1gen4=h1gen4.Integral()
        if not normalize:
          normsh1genpostfit=normsh1gen
        normsh2=h2.Integral()
        normsh3=h3.Integral()
        normsh4=h4.Integral()
        if compareNoUncertainty:
          normsh5=h5.Integral()

      h1.Scale(normsh1genpostfit/normsh1)
      h1postfit.Scale(normsh1genpostfit/normsh1postfit)
      h1gen.Scale(normsh1genpostfit/normsh1gen)
      if modelComparison:
        h1gen2.Scale(normsh1genpostfit/normsh1gen2)
        h1gen3.Scale(normsh1genpostfit/normsh1gen3)
        h1gen4.Scale(normsh1genpostfit/normsh1gen4)
        #print(h1genpostfit.Integral(),h1gen.Integral(),h1gen3.Integral(),h1gen4.Integral())
      if normalize:
        h1genpostfit.Scale(1./h1genpostfit.Integral())
      h2.Scale(normsh1genpostfit/normsh2)
      h3.Scale(normsh1genpostfit/normsh3)
      h4.Scale(normsh1genpostfit/normsh4)
      if compareNoUncertainty:
        h5.Scale(normsh1genpostfit/normsh5)
      yscale=max(h1gen.Integral(),h1genpostfit.Integral())


      for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,0)
        h1postfit.SetBinContent(b+1,h1postfit.GetBinContent(b+1)/h1postfit.GetBinWidth(b+1))
        h1postfit.SetBinError(b+1,0)
        h1gen.SetBinContent(b+1,h1gen.GetBinContent(b+1)/h1gen.GetBinWidth(b+1))
        h1gen.SetBinError(b+1,0)
        if modelComparison:
          h1gen2.SetBinContent(b+1,h1gen2.GetBinContent(b+1)/h1gen2.GetBinWidth(b+1))
          h1gen2.SetBinError(b+1,0)
          h1gen3.SetBinContent(b+1,h1gen3.GetBinContent(b+1)/h1gen3.GetBinWidth(b+1))
          h1gen3.SetBinError(b+1,0)
          h1gen4.SetBinContent(b+1,h1gen4.GetBinContent(b+1)/h1gen4.GetBinWidth(b+1))
          h1gen4.SetBinError(b+1,0)
        h1genpostfit.SetBinContent(b+1,h1genpostfit.GetBinContent(b+1)/h1genpostfit.GetBinWidth(b+1))
        h1genpostfit.SetBinError(b+1,h1genpostfit.GetBinError(b+1)/h1genpostfit.GetBinWidth(b+1))
        h2.SetBinContent(b+1,h2.GetBinContent(b+1)/h2.GetBinWidth(b+1))
        h2.SetBinError(b+1,h2.GetBinError(b+1)/h2.GetBinWidth(b+1))
        h3.SetBinContent(b+1,h3.GetBinContent(b+1)/h3.GetBinWidth(b+1))
        h3.SetBinError(b+1,h3.GetBinError(b+1)/h3.GetBinWidth(b+1))
        h4.SetBinContent(b+1,h4.GetBinContent(b+1)/h4.GetBinWidth(b+1))
        h4.SetBinError(b+1,h4.GetBinError(b+1)/h4.GetBinWidth(b+1))
        if compareNoUncertainty:
          h5.SetBinContent(b+1,h5.GetBinContent(b+1)/h5.GetBinWidth(b+1))
          h5.SetBinError(b+1,h5.GetBinError(b+1)/h5.GetBinWidth(b+1))

      if modelUncertainty:
        h3.Divide(h3,h1)
        h4.Divide(h4,h1)
        h1.Divide(h1,h1)
        h1gen=h1
      
      h1gen.SetLineColor(4)
      if not onlyData:
        h1gen.Draw("he")
      if modelComparison:
        h1gen2.SetLineColor(2)
        h1gen2.Draw("hesame")
        h1gen3.SetLineColor(6)
        h1gen3.Draw("hesame")
        h1gen4.SetLineColor(7)
        #h1gen4.Draw("hesame")
      
      if not trivialClosure and not crossClosureNNLO and not crossClosureHerwig and not crossClosureMadgraph and not onlyMC and not modelUncertainty and not modelComparison:
        h2.SetLineColor(6)
        h2.SetTitle("")
        if onlyData:
          h2.Draw("he")
          h1gen=h2
        else:
          h2.Draw("hesame")
 
      if i>=7:
          h1gen.GetYaxis().SetRangeUser(0.00*yscale,0.25*yscale)
      elif i>=6:
          h1gen.GetYaxis().SetRangeUser(0.02*yscale,0.22*yscale)
      elif i>=4:
          h1gen.GetYaxis().SetRangeUser(0.045*yscale,0.12*yscale)
      else:
          h1gen.GetYaxis().SetRangeUser(0.045*yscale,0.12*yscale)
      if not normalize:
        h1gen.GetYaxis().SetTitle("N")
      else:
        h1gen.GetYaxis().SetTitle("1/#sigma"+("" if normalize else "")+" d#sigma/d#chi")
      h1gen.GetXaxis().SetTitle("#chi")
     
      if not onlyMC and not modelUncertainty and not compareNoUncertainty:
        h1genpostfit.SetLineColor(1)
        h1genpostfit.SetLineStyle(2)
        h1genpostfit.SetLineWidth(2)
        h1genpostfit.Draw("hesame")
      
      if compareNoUncertainty:
        h5.SetLineColor(1)
        h5.SetLineStyle(2)
        h2.SetLineWidth(2)
        h5.SetTitle("")
        h5.Draw("he1same")
        h1genpostfit.SetLineColor(2)
        h1genpostfit.SetLineStyle(1)
        h1genpostfit.Draw("hesame")
      
      if not trivialClosure and not crossClosureNNLO and not crossClosureHerwig and not crossClosureMadgraph and not onlyData and not modelUncertainty and not modelComparison:
        h1.SetLineColor(2)
        h1.SetTitle("")
        h1.Draw("hesame")
      
      if modelUncertainty:
        h1.GetYaxis().SetRangeUser(0.95,1.1)
        h1.GetYaxis().SetTitle("Alternative / Pythia")
        h1.GetXaxis().SetTitle("#chi")
        h3.SetLineColor(6)
        h3.SetLineStyle(3)
        h3.Draw("hesame")
        fit=TF1(h3.GetName()+"fithw","pol1",1,16)
        plots+=[fit]
        h3.Fit(fit,"RWN")
        fit.SetLineColor(6)
        fit.SetLineWidth(2)
        fit.SetLineStyle(3)
        fit.Draw("lsame")
        hint = h3.Clone(h3.GetName()+"bandhw")
        TVirtualFitter.GetFitter().GetConfidenceIntervals(hint)
        hint.SetStats(False)
        hint.SetFillColor(6)
        hint.SetFillStyle(3395)
        hint.Draw("e3 same")
        print(fit.Eval(1),hint.GetBinError(1))
      
        h4.SetLineColor(7)
        h4.SetLineStyle(3)
        h4.Draw("hesame")
        fit=TF1(h4.GetName()+"fitmg","pol1",1,16)
        plots+=[fit]
        h4.Fit(fit,"RWN")
        fit.SetLineColor(7)
        fit.SetLineWidth(2)
        fit.SetLineStyle(3)
        fit.Draw("lsame")
        hint = h4.Clone(h4.GetName()+"bandmg")
        TVirtualFitter.GetFitter().GetConfidenceIntervals(hint)
        hint.SetStats(False)
        hint.SetFillColor(7)
        hint.SetFillStyle(3395)
        hint.Draw("e3 same")
        print(fit.Eval(1),hint.GetBinError(1))
   
      if not trivialClosure and not crossClosureNNLO and not crossClosureHerwig and not crossClosureMadgraph and not onlyMC and not onlyData and not modelUncertainty and not modelComparison:
        h1postfit.SetLineColor(2)
        h1postfit.SetLineStyle(2)
        h1postfit.Draw("hesame")
      
      # PLOTS
      if not onlyData and (trivialClosure or crossClosureNNLO or crossClosureHerwig or crossClosureMadgraph):
        legend1.AddEntry(h1gen,"Pythia particle-level","le")
      elif not modelUncertainty and not onlyData:
        legend1.AddEntry(h1gen,"NNLO particle-level","le")
      if modelComparison:
        legend1.AddEntry(h1gen2,"Pythia particle-level","l")
        legend1.AddEntry(h1gen3,"Madgraph particle-level","l")
        #legend1.AddEntry(h1gen4,"Herwig particle-level","l")
      if not trivialClosure and not crossClosureNNLO and not crossClosureHerwig and not crossClosureMadgraph and not onlyMC and not modelUncertainty and not modelComparison:
        legend1.AddEntry(h2,"Reco data (stat)","le")
      if crossClosureHerwig:
        legend1.AddEntry(h1genpostfit,"Unfolded (Herwig response)","le") # Pythia templates
      if crossClosureMadgraph:
        legend1.AddEntry(h1genpostfit,"Unfolded (Madgraph response)","le") # Pythia templates
      elif trivialClosure:
        legend1.AddEntry(h1genpostfit,"Unfolded (Pythia response)","le") # Pythia templates
      elif crossClosureNNLO:
        legend1.AddEntry(h1genpostfit,"Unfolded (Pythia response)","le") # NNLO templates
      elif onlyData and withUncertainties:
        legend1.AddEntry(h5,"Unfolded data (stat)","le")
        legend1.AddEntry(h1genpostfit,"Unfolded data (stat+sys)","le")
      elif not onlyMC and not modelUncertainty:
        legend1.AddEntry(h1genpostfit,"Unfolded data","le")
      if not trivialClosure and not crossClosureNNLO and not crossClosureHerwig and not crossClosureMadgraph and not onlyData and not modelComparison:
        legend1.AddEntry(h1,"Smeared (Pythia response)","l")
      if modelUncertainty:
        legend1.AddEntry(h3,"Smeared (Herwig response)","l")
        legend1.AddEntry(h4,"Smeared (Madgraph response)","l")
      if not trivialClosure and not crossClosureNNLO and not crossClosureHerwig and not crossClosureMadgraph and not onlyMC and not onlyData and not modelUncertainty and not modelComparison:
        legend1.AddEntry(h1postfit,"Smeared-Postfit","l")

      legend1.SetTextSize(0.04)
      legend1.SetFillStyle(0)
      legend1.Draw("same")
      
      if compareNoUncertainty:
        for b in range(h2.GetXaxis().GetNbins()):
          chi2detector+=pow((h2.GetBinContent(b+1)-h1.GetBinContent(b+1))/h2.GetBinError(b+1),2)
          chi2unfolded+=pow((h5.GetBinContent(b+1)-h1gen.GetBinContent(b+1))/h5.GetBinError(b+1),2)
      
    print("detector-level chi2", chi2detector, "unfolded chi2", chi2unfolded)

    canvas.SaveAs(prefix + "_"+name+"_run"+run+".pdf")
    fout.Close()
